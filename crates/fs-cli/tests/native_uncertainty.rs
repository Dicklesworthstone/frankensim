//! G1/G3/G4/G5 product checks for probability studies of native `.fsim` cooling
//! projects: real imported geometry, real material cards, and retained child QoIs.

#[path = "native_uncertainty/compliance.rs"]
mod compliance;
#[path = "native_uncertainty/fan_speed.rs"]
mod fan_speed;
#[path = "native_uncertainty/qmc.rs"]
mod qmc;
#[path = "native_uncertainty/copula.rs"]
mod copula;
#[path = "native_uncertainty/nominal_adjoint.rs"]
mod nominal_adjoint;
#[path = "../src/json_read.rs"]
mod json_read;

use json_read::JsonValue;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const STUDY: &str = r#"(fsim-uncertainty-study
    :version 1 :project "cooling-reference.fsim" :samples 4 :seed 29
    :wall-time 120s :method monte-carlo :correlation independent :qoi "temperature-max"
    :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
    :materials ("aa6061.fsmcdpk") :interfaces ()
    :parameters (
        (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
        (uniform :name "ambient" :target convection-temperature :entity "air" :low 294K :high 300K)))"#;

static NEXT: AtomicU64 = AtomicU64::new(0);

fn scratch() -> PathBuf {
    loop {
        let path = std::env::temp_dir().join(format!(
            "fs-cli-native-uncertainty-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        match std::fs::create_dir(&path) {
            Ok(()) => return path,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => panic!("scratch directory: {error}"),
        }
    }
}

fn command(args: &[&str], expected_exit: u8) -> JsonValue {
    let output = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(
        output.exit_code, expected_exit,
        "{args:?}\n{}\n{}",
        output.stdout, output.stderr
    );
    JsonValue::parse(&output.stdout).expect("command result JSON")
}

fn artifact(ledger: &Path, hash: &str) -> Vec<u8> {
    let hash = fs_blake3::ContentHash::from_hex(hash).expect("content hash");
    fs_ledger::Ledger::open(ledger.to_str().unwrap())
        .unwrap()
        .get_artifact(&hash)
        .unwrap()
        .expect("artifact retained in the native study ledger")
}

fn json_artifact(ledger: &Path, hash: &str) -> JsonValue {
    JsonValue::parse(std::str::from_utf8(&artifact(ledger, hash)).unwrap()).unwrap()
}

fn rows(report: &JsonValue) -> &[JsonValue] {
    report.get("observations").unwrap().as_array().unwrap()
}

struct Fixture {
    dir: PathBuf,
    sources: PathBuf,
    ledger: PathBuf,
    project: fs_project::ProjectSpec,
}

impl Fixture {
    fn new() -> Self {
        let dir = scratch();
        let sources = dir.join("sources");
        std::fs::create_dir(&sources).unwrap();
        let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
        let mut project = fs_project::parse_sexpr_migrating(
            &std::fs::read_to_string(reference.join("cooling-reference.fsim")).unwrap(),
        )
        .unwrap()
        .decoded
        .spec;
        // Exercise the declared margin, rather than comparing samples with the
        // unadjusted temperature limit. The effective threshold is 297 K.
        let requirement = &mut project.requirements.as_mut().unwrap()[0];
        requirement.limit.value = 302.0;
        requirement.margin.value = 5.0;
        std::fs::write(
            sources.join("cooling-reference.fsim"),
            fs_project::print_sexpr(&project).unwrap(),
        )
        .unwrap();
        for name in ["plate.stl", "aa6061.fsmcdpk"] {
            std::fs::copy(reference.join(name), sources.join(name)).unwrap();
        }
        std::fs::write(sources.join("study.fsim"), STUDY).unwrap();
        Self {
            ledger: dir.join("study.db"),
            dir,
            sources,
            project,
        }
    }

    fn study(&self, budget: Option<&str>, expected_exit: u8) -> JsonValue {
        let source = self.sources.join("study.fsim");
        let mut args = vec![
            "--json",
            "study",
            source.to_str().unwrap(),
            self.ledger.to_str().unwrap(),
        ];
        if let Some(budget) = budget {
            args.extend(["--budget", budget]);
        }
        command(&args, expected_exit)
    }

    fn resume(&self, run: &str, budget: Option<&str>, expected_exit: u8) -> JsonValue {
        let mut args = vec![
            "--json",
            "study",
            "--resume",
            run,
            self.ledger.to_str().unwrap(),
        ];
        if let Some(budget) = budget {
            args.extend(["--budget", budget]);
        }
        command(&args, expected_exit)
    }

    fn report(&self, result: &JsonValue) -> (Vec<u8>, JsonValue) {
        assert_eq!(result.str_field("command"), Some("study"));
        assert_eq!(result.str_field("run_id"), result.str_field("run"));
        assert!(result.str_field("run").unwrap().starts_with("study-"));
        let receipt = result.get("receipt").unwrap();
        assert_eq!(
            receipt.str_field("schema"),
            Some("frankensim.cli.study-run-receipt.v1")
        );
        assert_eq!(
            receipt.str_field("driver"),
            Some("native-cooling-uncertainty-v1")
        );
        assert_eq!(receipt.str_field("status"), result.str_field("status"));
        assert_eq!(receipt.f64_field("samples_planned"), Some(4.0));
        let bytes = artifact(&self.ledger, receipt.str_field("report_json").unwrap());
        let report = JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
        (bytes, report)
    }

    fn independent_solve(&self, ordinal: usize, parameters: &[f64]) -> (String, String, f64) {
        assert_eq!(parameters.len(), 2);
        // Deliberately bind directly to the ordinary native project, independently
        // of the study's sample_project helper and evaluation/report code.
        let mut project = self.project.clone();
        project.power.as_mut().unwrap()[0].watts.value = parameters[0];
        let fs_project::ThermalBoundaryCondition::Convection {
            reference_temperature,
            ..
        } = &mut project
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .boundaries[0]
            .condition
        else {
            panic!("reference fixture uses a declared convection reservoir")
        };
        reference_temperature.value = parameters[1];
        self.solve_project(ordinal, &project)
    }

    fn solve_project(&self, ordinal: usize, project: &fs_project::ProjectSpec) -> (String, String, f64) {
        let source = fs_project::print_sexpr(project).unwrap();
        let project_hash = fs_project::parse_sexpr(&source).unwrap().hash().to_hex();
        let project_path = self.dir.join(format!("independent-{ordinal}.fsim"));
        std::fs::write(&project_path, source).unwrap();
        let ledger = self.dir.join("independent.db");
        command(
            &[
                "--json",
                "import",
                project_path.to_str().unwrap(),
                self.sources.join("plate.stl").to_str().unwrap(),
                ledger.to_str().unwrap(),
                "--unit",
                "m",
                "--max-hole-edges",
                "0",
            ],
            fs_cli::exit::SUCCESS,
        );
        let solved = command(
            &[
                "--json",
                "solve",
                project_path.to_str().unwrap(),
                ledger.to_str().unwrap(),
                "--materials",
                self.sources.join("aa6061.fsmcdpk").to_str().unwrap(),
            ],
            fs_cli::exit::SUCCESS,
        );
        assert_eq!(solved.str_field("status"), Some("completed"));
        let run = solved.str_field("run").unwrap().to_string();
        let receipt = json_artifact(&ledger, solved.str_field("run_receipt").unwrap());
        let qoi_hash = receipt
            .get("stages")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .find(|stage| stage.str_field("stage") == Some("qoi"))
            .unwrap()
            .str_field("receipt")
            .unwrap();
        let qoi = json_artifact(&ledger, qoi_hash);
        let value = qoi.get("qoi").unwrap().as_array().unwrap()[0]
            .f64_field("value")
            .unwrap();
        (project_hash, run, value)
    }
}

#[test]
fn g1_native_study_matches_independent_project_solves_and_retains_each_qoi() {
    let fixture = Fixture::new();
    let result = fixture.study(None, fs_cli::exit::SUCCESS);
    assert_eq!(result.str_field("status"), Some("completed"));
    let (_, report) = fixture.report(&result);
    assert!(
        report.get("compliance").is_none(),
        "version 1 retains its fixed-count report schema"
    );
    assert_eq!(report.f64_field("samples_evaluated"), Some(4.0));
    assert_eq!(rows(&report).len(), 4);
    let mut values = Vec::new();
    for (ordinal, row) in rows(&report).iter().enumerate() {
        assert_eq!(row.f64_field("ordinal"), Some(ordinal as f64));
        let parameters = row
            .get("parameters")
            .unwrap()
            .as_array()
            .unwrap()
            .iter()
            .map(|p| p.as_f64().unwrap())
            .collect::<Vec<_>>();
        assert!((4.0..=6.0).contains(&parameters[0]));
        assert!((294.0..=300.0).contains(&parameters[1]));
        let (project, run, value) = fixture.independent_solve(ordinal, &parameters);
        assert_eq!(row.str_field("project_hash"), Some(project.as_str()));
        assert_eq!(row.str_field("run"), Some(run.as_str()));
        assert_eq!(row.f64_field("value_k").unwrap().to_bits(), value.to_bits());
        let retained_qoi = json_artifact(&fixture.ledger, row.str_field("qoi_receipt").unwrap());
        assert_eq!(retained_qoi.str_field("run"), Some(run.as_str()));
        assert_eq!(
            retained_qoi.path(&["lineage", "project"]).unwrap().as_str(),
            Some(project.as_str())
        );
        let qoi = &retained_qoi.get("qoi").unwrap().as_array().unwrap()[0];
        assert_eq!(qoi.str_field("name"), Some("temperature-max"));
        assert_eq!(qoi.str_field("color"), Some("estimated"));
        assert_eq!(qoi.f64_field("value").unwrap().to_bits(), value.to_bits());
        values.push(value);
    }
    assert!(
        values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - values.iter().copied().fold(f64::INFINITY, f64::min)
            > 0.1
    );
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let std_dev = (values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64)
        .sqrt();
    let pass_rate =
        values.iter().filter(|&&value| value <= 297.0).count() as f64 / values.len() as f64;
    let statistics = report.get("statistics").unwrap();
    assert!((statistics.f64_field("mean_k").unwrap() - mean).abs() < 1.0e-11);
    assert!((statistics.f64_field("std_dev_k").unwrap() - std_dev).abs() < 1.0e-11);
    assert_eq!(
        statistics.f64_field("empirical_probability_of_compliance"),
        Some(pass_rate)
    );
}

#[test]
fn g4_g5_zero_budget_and_resume_without_source_files_match_uninterrupted_report() {
    let uninterrupted = Fixture::new();
    let complete = uninterrupted.study(None, fs_cli::exit::SUCCESS);
    let (expected_report, expected) = uninterrupted.report(&complete);
    let resumed = Fixture::new();
    let zero = resumed.study(Some("0"), fs_cli::exit::BUDGET);
    assert_eq!(zero.str_field("status"), Some("budget-exhausted"));
    let (_, report) = resumed.report(&zero);
    assert_eq!(report.f64_field("samples_evaluated"), Some(0.0));
    assert!(rows(&report).is_empty());
    assert_eq!(report.get("statistics"), Some(&JsonValue::Null));
    // Preserve every file while removing all original source paths. The ledger
    // must be sufficient for every subsequent import, material load and solve.
    std::fs::rename(&resumed.sources, resumed.dir.join("relocated-sources")).unwrap();
    let first = resumed.resume(
        zero.str_field("run").unwrap(),
        Some("1"),
        fs_cli::exit::BUDGET,
    );
    let (_, report) = resumed.report(&first);
    assert_eq!(first.str_field("status"), Some("budget-exhausted"));
    assert_eq!(report.f64_field("samples_evaluated"), Some(1.0));
    assert_eq!(rows(&report), &rows(&expected)[..1]);
    assert_eq!(report.get("statistics"), Some(&JsonValue::Null));
    let finished = resumed.resume(first.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(finished.str_field("status"), Some("completed"));
    let (actual_report, _) = resumed.report(&finished);
    assert_eq!(actual_report, expected_report);
    // Export reads the retained artifacts. It must not need the original files
    // or invoke fresh numerical work, including on the existing command routes.
    let run = finished.str_field("run").unwrap();
    let exported = command(
        &["--json", "report", run, resumed.ledger.to_str().unwrap()],
        fs_cli::exit::SUCCESS,
    );
    assert_eq!(
        std::fs::read(exported.str_field("report_json").unwrap()).unwrap(),
        expected_report
    );
    let receipt = finished.get("receipt").unwrap();
    assert_eq!(
        std::fs::read(exported.str_field("report_html").unwrap()).unwrap(),
        artifact(&resumed.ledger, receipt.str_field("report_html").unwrap())
    );
    let exported = command(
        &["--json", "package", run, resumed.ledger.to_str().unwrap()],
        fs_cli::exit::SUCCESS,
    );
    assert_eq!(
        std::fs::read(exported.str_field("package").unwrap()).unwrap(),
        artifact(&resumed.ledger, receipt.str_field("package").unwrap())
    );
    let repeated = resumed.resume(run, None, fs_cli::exit::SUCCESS);
    assert_eq!(
        repeated, finished,
        "completed resume returns the retained terminal result"
    );
}

#[test]
fn g4_native_child_refusal_cannot_be_counted_skipped_or_resumed_as_success() {
    let fixture = Fixture::new();
    std::fs::write(
        fixture.sources.join("study.fsim"),
        STUDY.replace(":materials (\"aa6061.fsmcdpk\")", ":materials ()"),
    )
    .unwrap();
    let result = fixture.study(None, fs_cli::exit::REFUSED);
    assert_eq!(result.str_field("status"), Some("refused"));
    let (_, report) = fixture.report(&result);
    assert_eq!(report.f64_field("samples_evaluated"), Some(0.0));
    assert_eq!(report.f64_field("evaluations_attempted"), Some(1.0));
    assert!(rows(&report).is_empty());
    assert_eq!(report.get("statistics"), Some(&JsonValue::Null));
    assert!(
        report
            .str_field("failure")
            .unwrap()
            .contains("project-material-card-unknown")
    );
    let repeated = fixture.resume(
        result.str_field("run").unwrap(),
        None,
        fs_cli::exit::REFUSED,
    );
    assert_eq!(
        repeated, result,
        "a failed physical sample cannot be silently retried or skipped"
    );
}
