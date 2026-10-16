//! Native flow-network/solid-air/adjoint/UQ path with real geometry and cards.
#[path = "../src/json_read.rs"]
mod json_read;
#[path = "native_air_inlet_adjoint/coupled_thermal.rs"]
mod coupled_thermal;
use json_read::JsonValue as J;
use std::path::{Path, PathBuf};
use fs_project::ThermalBoundaryCondition as B;

const STUDY: &str = r#"(fsim-uncertainty-study :version 1
    :project "cooling.fsim" :samples 4 :seed 29 :wall-time 600s
    :method monte-carlo :correlation independent :qoi "temperature-max"
    :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
    :materials ("aa6061.fsmcdpk") :interfaces () :parameters (
        (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
        (uniform :name "inlet" :target air-inlet-temperature :entity "air" :low 294K :high 300K)))"#;

fn command(args: &[&str], exit: u8) -> J {
    let output = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(output.exit_code, exit, "{args:?}\n{}\n{}", output.stdout, output.stderr);
    J::parse(&output.stdout).unwrap()
}
fn bytes(ledger: &Path, identity: &str) -> Vec<u8> {
    fs_ledger::Ledger::open(ledger.to_str().unwrap()).unwrap()
        .get_artifact(&fs_blake3::ContentHash::from_hex(identity).unwrap()).unwrap().unwrap()
}
fn artifact(ledger: &Path, identity: &str) -> J {
    J::parse(std::str::from_utf8(&bytes(ledger, identity)).unwrap()).unwrap()
}
fn inlet(project: &mut fs_project::ProjectSpec, value: f64) {
    let B::AirflowConvection { inlet_temperature, .. } = &mut project.cooling.as_mut().unwrap()
        .conduction.as_mut().unwrap().boundaries[0].condition else { panic!("air inlet fixture"); };
    inlet_temperature.value = value;
}
struct Fixture { dir: PathBuf, sources: PathBuf, ledger: PathBuf, project: fs_project::ProjectSpec }
impl Fixture {
    fn new() -> Self {
        let dir = (0..10000).find_map(|i| {
            let path = std::env::temp_dir().join(format!("fs-native-inlet-adjoint-{}-{i}", std::process::id()));
            match std::fs::create_dir(&path) {
                Ok(()) => Some(path),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => None,
                Err(e) => panic!("test directory: {e}"),
            }
        }).expect("fresh fixture directory");
        let sources = dir.join("sources"); std::fs::create_dir(&sources).unwrap();
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
        let mut project = fs_project::parse_sexpr_migrating(
            &std::fs::read_to_string(root.join("cooling-reference.fsim")).unwrap()).unwrap().decoded.spec;
        for file in ["plate.stl", "aa6061.fsmcdpk"] { std::fs::copy(root.join(file), sources.join(file)).unwrap(); }
        use fs_qty::QtyAny;
        use fs_project::spec::dims;
        project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition = B::AirflowConvection {
            branch: "air".into(), order: 0, inlet_temperature: QtyAny::new(297.0, dims::TEMPERATURE),
            hydraulic_diameter: QtyAny::new(0.02, dims::LENGTH), flow_area: QtyAny::new(0.004, dims::AREA),
            channel_length: QtyAny::new(0.3, dims::LENGTH), correlation: "convection.gnielinski".into(),
        };
        project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
        project.power.as_mut().unwrap()[0].duty = 0.37;
        std::fs::write(sources.join("cooling.fsim"), fs_project::print_sexpr(&project).unwrap()).unwrap();
        Self { ledger: dir.join("study.db"), dir, sources, project }
    }
    fn solve(&self, project: &fs_project::ProjectSpec, ordinal: usize) -> (J, J) {
        let path = self.dir.join(format!("project-{ordinal}.fsim"));
        std::fs::write(&path, fs_project::print_sexpr(project).unwrap()).unwrap();
        let ledger = self.dir.join("solves.db");
        command(&["--json","import",path.to_str().unwrap(),self.sources.join("plate.stl").to_str().unwrap(),
            ledger.to_str().unwrap(),"--unit","m","--max-hole-edges","0"], fs_cli::exit::SUCCESS);
        let run = command(&["--json","solve",path.to_str().unwrap(),ledger.to_str().unwrap(),
            "--materials",self.sources.join("aa6061.fsmcdpk").to_str().unwrap()], fs_cli::exit::SUCCESS);
        let receipt = artifact(&ledger, run.str_field("run_receipt").unwrap());
        let stage = receipt.get("stages").unwrap().as_array().unwrap().iter()
            .find(|s| s.str_field("stage") == Some("conduction")).unwrap();
        let receipt = artifact(&ledger, stage.str_field("receipt").unwrap());
        let field = artifact(&ledger, receipt.str_field("solution_artifact").unwrap());
        (receipt, field)
    }
    fn study(&self, text: &str, budget: Option<&str>, exit: u8) -> J {
        let path = self.sources.join("study.fsim"); std::fs::write(&path, text).unwrap();
        let mut args = vec!["--json", "study", path.to_str().unwrap(), self.ledger.to_str().unwrap()];
        if let Some(budget) = budget { args.extend(["--budget", budget]); }
        command(&args, exit)
    }
    fn resume(&self, run: &str, budget: Option<&str>, exit: u8) -> J {
        let mut args = vec!["--json", "study", "--resume", run, self.ledger.to_str().unwrap()];
        if let Some(budget) = budget { args.extend(["--budget", budget]); }
        command(&args, exit)
    }
    fn report(&self, result: &J) -> (Vec<u8>, J) {
        let hash = result.get("receipt").unwrap().str_field("report_json").unwrap();
        let bytes = bytes(&self.ledger, hash);
        let report = J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
        (bytes, report)
    }
}
fn study_text(copula: bool, qmc: bool, controlled: bool) -> String {
    let mut text = STUDY.to_string();
    if copula { text = text.replace("independent", "(gaussian-copula :latent-correlation ((1 0.5) (0.5 1)))"); }
    if qmc { text = text.replace("monte-carlo", "quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)"); }
    if controlled { text = text.replace(":version 1", ":version 3 :mean-control (nominal-adjoint :max-solves 1)"); }
    text
}

#[test]
fn native_air_inlet_adjoint_matches_physical_shift_without_changing_the_primal() {
    let f = Fixture::new();
    let (baseline, field) = f.solve(&f.project, 0);
    let mut requested = f.project.clone();
    requested.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
        name: "temperature-max-adjoint".into(), kind: "report".into(), region: None,
    });
    let (receipt, requested_field) = f.solve(&requested, 1);
    assert_eq!(field.get("temperature"), requested_field.get("temperature"));
    assert_eq!(baseline.get("energy"), receipt.get("energy"));
    assert_eq!(baseline.get("conjugate"), receipt.get("conjugate"));
    let adjoint = receipt.get("nominal_adjoint").unwrap();
    assert_eq!(adjoint.str_field("mode"), Some("linear-solid-full-air-feedback"));
    let inlet_rows: Vec<_> = adjoint.get("parameters").unwrap().as_array().unwrap().iter()
        .filter(|r| r.str_field("target") == Some("air-inlet-temperature")).collect();
    assert_eq!(inlet_rows.len(), 1);
    assert_eq!(inlet_rows[0].str_field("entity"), Some("air"));
    let derivative = inlet_rows[0].f64_field("derivative").unwrap();
    let vertex = adjoint.f64_field("selected_vertex").unwrap() as usize;
    let mut values = Vec::new();
    for (i, t) in [296.9, 297.1].into_iter().enumerate() {
        let mut project = f.project.clone(); inlet(&mut project, t);
        let (_, shifted) = f.solve(&project, i+2);
        values.push(shifted.get("temperature").unwrap().as_array().unwrap()[vertex].as_f64().unwrap());
    }
    assert!((derivative-(values[1]-values[0])/0.2).abs() < 1e-6);
    assert!((derivative-1.0).abs() < 1e-6, "one prescribed inlet shifts this whole linear field");
    assert!(!adjoint.get("unsupported").unwrap().as_array().unwrap().iter()
        .any(|r| r.str_field("target") == Some("air-inlet-temperature")));
}

#[test]
fn native_air_inlet_mean_control_preserves_all_four_sampling_laws() {
    for (copula, qmc) in [(false,false),(true,false),(false,true),(true,true)] {
        let raw = Fixture::new();
        let result = raw.study(&study_text(copula,qmc,false), None, fs_cli::exit::SUCCESS);
        let (_, baseline) = raw.report(&result);
        let f = Fixture::new();
        let result = f.study(&study_text(copula,qmc,true), None, fs_cli::exit::SUCCESS);
        let (_, report) = f.report(&result);
        for key in ["observations", "statistics", "qmc"] { assert_eq!(report.get(key), baseline.get(key), "{key}"); }
        let control = report.get("mean_control").unwrap();
        assert_eq!(control.str_field("method"), Some("nominal-adjoint"));
        assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
        assert_eq!(control.get("estimate").unwrap().f64_field("samples_in_estimate"), Some(4.0));
        assert!(control.get("estimate").unwrap().f64_field("variance_ratio").unwrap() < 1e-4,
            "fixed-flow linear power/inlet dependence should be removed from the controlled mean");
    }
}

#[test]
fn native_air_inlet_calibration_resumes_from_retained_sources_and_coefficients() {
    let text = study_text(true,true,true);
    let whole = Fixture::new();
    let complete = whole.study(&text, None, fs_cli::exit::SUCCESS);
    let (expected, _) = whole.report(&complete);
    let f = Fixture::new();
    let prefix = f.study(&text, Some("1"), fs_cli::exit::BUDGET);
    let (_, report) = f.report(&prefix);
    assert!(report.get("observations").unwrap().as_array().unwrap().is_empty());
    let calibration = report.get("mean_control").unwrap().get("calibration").unwrap().clone();
    std::fs::rename(&f.sources, f.dir.join("relocated-sources")).unwrap();
    let one = f.resume(prefix.str_field("run").unwrap(), Some("1"), fs_cli::exit::BUDGET);
    let (_, report) = f.report(&one);
    assert_eq!(report.get("observations").unwrap().as_array().unwrap().len(), 1);
    assert_eq!(report.get("mean_control").unwrap().get("calibration"), Some(&calibration));
    let final_run = f.resume(one.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(f.report(&final_run).0, expected);
    assert_eq!(f.resume(final_run.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS), final_run);
}
