use super::*;
use crate::json_read::JsonValue as J;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

const SOURCE: &str = r#"(fsim-uncertainty-study
 :version 1 :project "cooling-reference.fsim" :samples 4 :seed 29
 :wall-time 120s :method monte-carlo :correlation independent :qoi "temperature-max"
 :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
 :materials ("aa6061.fsmcdpk") :interfaces ()
 :parameters (
  (uniform :name "power" :target power :entity "air" :low 4W :high 6W)
  (uniform :name "ambient" :target convection-temperature :entity "air" :low 294K :high 300K)))"#;

struct Fixture { root: PathBuf, assets: PathBuf, source: PathBuf, pins: StudyPins }
impl Fixture {
    fn new(text: &str) -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let root = loop {
            let path = std::env::temp_dir().join(format!("fs-prepared-study-{}-{}",
                std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
            match std::fs::create_dir(&path) {
                Ok(()) => break path,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("fixture: {error}"),
            }
        };
        let assets = root.join("assets");
        std::fs::create_dir(&assets).unwrap();
        let reference = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project");
        for name in ["cooling-reference.fsim", "plate.stl", "aa6061.fsmcdpk"] {
            std::fs::copy(reference.join(name), assets.join(name)).unwrap();
        }
        let source = assets.join("study.fsim");
        std::fs::write(&source, text).unwrap();
        let base = crate::read_project_for_solve(&assets.join("cooling-reference.fsim"), OutputMode::Json)
            .unwrap_or_else(|output| panic!("{}", output.stderr));
        let study = fs_project::uncertainty::UncertaintyStudy::parse(text).unwrap();
        let pins = StudyPins { project: base.hash(),
            source: Some(hash_bytes(study.canonical().as_bytes())), wall_seconds: Some(120.0) };
        Self { root, assets, source, pins }
    }
    fn prepare(&self) -> PreparedStudy {
        PreparedStudy::load(&self.source, self.pins, 16 * 1024 * 1024).unwrap()
    }
    fn json(output: CommandOutput, code: u8) -> J {
        assert_eq!(output.exit_code, code, "{}\n{}", output.stdout, output.stderr);
        J::parse(&output.stdout).unwrap()
    }
}

#[test]
fn prepared_execution_uses_original_project_geometry_cards_and_sampling_law() {
    for controlled in [false, true] {
        let text = if controlled {
            SOURCE.replace(":version 1", ":version 3 :mean-control (nominal-adjoint :max-solves 1)")
                .replace("monte-carlo", "quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)")
        } else { SOURCE.to_string() };
        let f = Fixture::new(&text);
        let prepared = f.prepare();
        let direct = Fixture::json(crate::study::study_path(&f.source,
            &f.root.join("direct.db"), None, OutputMode::Json), crate::exit::SUCCESS);
        // Remove ALL original paths after binding. The snapshot must remain
        // executable, not just notice the mismatch before reopening a file.
        std::fs::rename(&f.assets, f.root.join("moved-assets")).unwrap();
        let actual = Fixture::json(prepared.run(&f.root.join("prepared.db"), None,
            OutputMode::Json), crate::exit::SUCCESS);
        let expected = direct.get("receipt").unwrap();
        let receipt = actual.get("receipt").unwrap();
        assert_eq!(receipt.str_field("study_id"), expected.str_field("study_id"));
        assert_eq!(receipt.str_field("model"), expected.str_field("model"));
        assert_eq!(receipt.str_field("report_json"), expected.str_field("report_json"));
        assert_eq!(receipt.str_field("observations"), expected.str_field("observations"));
        // Receipt pointers include elapsed time; they are not physical model IDs.
    }
}

#[test]
fn prepared_prefix_uses_the_existing_source_free_recovery_path() {
    let f = Fixture::new(SOURCE);
    let prepared = f.prepare();
    let direct = Fixture::json(prepared.run(&f.root.join("complete.db"), None,
        OutputMode::Json), crate::exit::SUCCESS);
    std::fs::rename(&f.assets, f.root.join("moved-assets")).unwrap();
    let path = f.root.join("split.db");
    let prefix = Fixture::json(prepared.run(&path, Some("1"), OutputMode::Json), crate::exit::BUDGET);
    assert_eq!(prefix.get("receipt").unwrap().f64_field("samples_completed"), Some(1.0));
    let resumed = Fixture::json(crate::study::resume_path(prefix.str_field("run").unwrap(),
        &path, None, OutputMode::Json), crate::exit::SUCCESS);
    assert_eq!(resumed.get("receipt").unwrap().str_field("report_json"),
        direct.get("receipt").unwrap().str_field("report_json"));
}

#[test]
fn wrong_pins_or_grants_refuse_before_attempting_to_read_assets() {
    let f = Fixture::new(SOURCE);
    std::fs::rename(f.assets.join("plate.stl"), f.assets.join("hidden-plate.stl")).unwrap();
    let wrong = hash_bytes(b"not the admitted model");
    for (pins, fragment) in [
        (StudyPins { source: Some(wrong), ..f.pins }, "pinned :hash"),
        (StudyPins { project: wrong, ..f.pins }, "physical project differs"),
        (StudyPins { wall_seconds: Some(60.0), ..f.pins }, "wall allowance"),
    ] {
        let error = PreparedStudy::load(&f.source, pins, 16 * 1024 * 1024).unwrap_err();
        assert!(error.contains(fragment), "{error}");
        assert!(!error.contains("cannot open"), "pins must precede geometry reads");
    }
}

#[test]
fn caller_storage_limit_is_real_and_reused_snapshots_recheck_constraints() {
    let f = Fixture::new(SOURCE);
    let prepared = f.prepare();
    assert!(prepared.input_bytes() > 0);
    assert!(PreparedStudy::load(&f.source, f.pins, prepared.input_bytes() - 1).is_err());
    assert!(PreparedStudy::load(&f.source, f.pins, 0).is_err());
    prepared.check(f.pins).unwrap();
    assert!(prepared.check(StudyPins { source: Some(hash_bytes(b"different")), ..f.pins }).is_err());
    let ledger = f.root.join("invalid-budget.db");
    assert_eq!(prepared.run(&ledger, Some("257"), OutputMode::Json).exit_code, crate::exit::REFUSED);
    assert!(!ledger.exists(), "an invalid invocation cannot open a ledger");
}
