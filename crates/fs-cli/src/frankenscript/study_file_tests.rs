//! The program must execute its admitted physics after source paths change.
use super::*;
use crate::json_read::JsonValue;
use std::sync::atomic::{AtomicU64, Ordering};

struct Fixture {
    root: PathBuf,
    source: PathBuf,
    ledger: PathBuf,
}

impl Fixture {
    fn new(source: &str) -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let root = loop {
            let path = std::env::temp_dir().join(format!("fs-script-study-file-{}-{}",
                std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
            match std::fs::create_dir(&path) {
                Ok(()) => break path,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("fixture: {error}"),
            }
        };
        let path = root.join("study.fsim");
        std::fs::write(&path, source).unwrap();
        Self { ledger: root.join("program.db"), root, source: path }
    }

    fn binder(&self) -> Binder<'_> {
        Binder {
            base: self.root.clone(), ledger: &self.ledger, seed: Some(1337),
            wall_seconds: Some(600.0),
            native_studies: native_study::SnapshotCache::new(Some(1024 * 1024 * 1024)),
            projects: Vec::new(), studies: Vec::new(),
        }
    }
}

fn thermal_source() -> String {
    include_str!("../../../../examples/marquee/thermal-2d.fsim")
        .replace(":mesh-level 4", ":mesh-level 3")
        .replace(":max-iterations 8", ":max-iterations 1")
        .replace(":steps 8", ":steps 1")
        // These words are metadata, not a request to change the PDE.
        .replace("Normalized scalar thermal cooling-hole radius study",
            "Compare against elasticity-2d and fsim-sdf3-study later")
}

fn elasticity_source() -> String {
    include_str!("../../../../examples/marquee/bracket-2d.fsim")
        .replace(":mesh-level 5", ":mesh-level 2")
        .replace(":max-iterations 32", ":max-iterations 1")
        .replace(":steps 32", ":steps 1")
        .replace(":move-cells 0.35", ":move-cells 0.05")
        .replace(":nucleation-period 4", ":nucleation-period 0")
}

fn bind(binder: &mut Binder<'_>, pin: &str, budget: &str) -> (Option<Step>, Refusals) {
    let node = fs_ir::sexpr::parse(&format!(
        "(study.run \"study.fsim\" :hash \"{pin}\" :budget {budget})",
    )).unwrap();
    let mut refusals = Refusals {
        mode: OutputMode::Json, label: "snapshot-test".into(),
        stderr: String::new(), count: 0,
    };
    let step = binder.bind(&node, &mut refusals);
    (step, refusals)
}

fn captured(binder: &mut Binder<'_>, pin: &str) -> native_study::FileStudyStep {
    let (step, refusals) = bind(binder, pin, "1");
    assert_eq!(refusals.count, 0, "{}", refusals.stderr);
    let Some(Step::StudyFile(step)) = step else { panic!("bound study") };
    step
}

fn assert_physical_snapshot(source: &str, replacement: &str) {
    let fixture = Fixture::new(source);
    let pin = fs_blake3::hash_bytes(source.as_bytes()).to_hex();
    let mut binder = fixture.binder();
    let prepared = captured(&mut binder, &pin);
    assert!(!fixture.ledger.exists(), "binding must not execute or open the ledger");

    // The ordinary command exercises the same producer and source-based routing.
    let expected = crate::study::study_path(
        &fixture.source, &fixture.root.join("direct.db"), Some("1"), OutputMode::Json,
    );
    let expected_json = JsonValue::parse(&expected.stdout)
        .unwrap_or_else(|error| panic!("{error:?}: {}", expected.stderr));
    let expected_receipt = expected_json.get("receipt")
        .unwrap_or_else(|| panic!("no physical study receipt: {}", expected.stderr));
    assert_eq!(expected_receipt.get("iterations_completed").and_then(JsonValue::number_raw),
        Some("1"), "the regression must execute actual numerical work");
    assert_ne!(source, replacement);
    std::fs::write(&fixture.source, replacement).unwrap();

    // A second clause on the same path keeps the first captured model.
    let repeated = captured(&mut binder, &pin);
    assert_eq!(prepared.source_hash(), repeated.source_hash());
    assert_eq!(binder.studies.len(), 2);
    let actual = prepared.run(&fixture.ledger, OutputMode::Json);
    assert_eq!(actual.exit_code, expected.exit_code, "{}", actual.stderr);
    let actual_json = JsonValue::parse(&actual.stdout).unwrap();
    let actual_receipt = actual_json.get("receipt").expect("physical receipt");
    for key in ["study_id", "source", "trace_hash", "iterations", "design", "report_json"] {
        let expected = expected_receipt.str_field(key).expect("retained physical artifact");
        assert_eq!(actual_receipt.str_field(key), Some(expected), "{key}");
    }
}

#[test]
fn g3_thermal_snapshot_ignores_later_file_edits_and_physics_names_in_metadata() {
    let source = thermal_source();
    let changed = source.replace(":radius 0.12", ":radius 0.08");
    assert_physical_snapshot(&source, &changed);
}

#[test]
fn g3_elasticity_snapshot_executes_the_original_load_after_file_replacement() {
    let source = elasticity_source();
    let changed = source.replace(":load-traction-pa 1000000.0", ":load-traction-pa 2000000.0");
    assert_physical_snapshot(&source, &changed);
}

#[test]
fn reused_snapshot_rechecks_source_pin_typed_seed_and_producer_budget() {
    let source = elasticity_source();
    let fixture = Fixture::new(&source);
    let pin = fs_blake3::hash_bytes(source.as_bytes()).to_hex();
    let mut binder = fixture.binder();
    captured(&mut binder, &pin);
    std::fs::rename(&fixture.source, fixture.root.join("moved.fsim")).unwrap();
    captured(&mut binder, &pin);
    let retained = binder.studies.len();
    let wrong = fs_blake3::hash_bytes(b"unreviewed study").to_hex();
    for (hash, seed, budget, code) in [
        (wrong.as_str(), 1337, "1", "frankenscript-study-hash"),
        (pin.as_str(), 1338, "1", "frankenscript-explicit-seed"),
        (pin.as_str(), 1337, "33", "frankenscript-study-budget"),
    ] {
        binder.seed = Some(seed);
        let (step, refusals) = bind(&mut binder, hash, budget);
        assert!(step.is_none());
        assert!(refusals.stderr.contains(code), "{}", refusals.stderr);
        assert_eq!(binder.studies.len(), retained);
        assert!(!fixture.ledger.exists());
    }
}

#[test]
fn g4_invalid_later_physics_or_budget_refuses_before_the_first_program_step() {
    let source = elasticity_source();
    for bad_physics in [true, false] {
        let fixture = Fixture::new(&source);
        let bad_source = if bad_physics {
            source.replace(":youngs-modulus-pa 70000000000.0", ":youngs-modulus-pa -1.0")
        } else { source.clone() };
        std::fs::write(fixture.root.join("second.fsim"), bad_source).unwrap();
        let budget = if bad_physics { 1 } else { 33 };
        let program = fixture.root.join("program.fs");
        std::fs::write(&program, format!(
            "(study \"atomic-admission\" (seed 0x539)
             (versions (constellation :lock \"2026-07\"))
             (budget (wall 600s) (mem 1024MiB))
             (capability :cores 1 :mem 1024MiB :wall 600s :ops (study.*))
             (study.run \"study.fsim\" :budget 1)
             (study.run \"second.fsim\" :budget {budget}))"
        )).unwrap();
        let output = run_program_path(&program, &fixture.ledger, OutputMode::Json);
        assert_eq!(output.exit_code, exit::REFUSED, "{}", output.stderr);
        let code = if bad_physics { "frankenscript-study-input" }
            else { "frankenscript-study-budget" };
        assert!(output.stderr.contains(code), "{}", output.stderr);
        assert!(!fixture.ledger.exists(), "a later admission refusal must prevent the first solve");
    }
}

#[cfg(feature = "sdf3-study")]
#[test]
fn sdf3_snapshot_retains_owned_geometry_loads_and_original_source_pin() {
    let source = include_str!("../../../../examples/marquee/bracket-3d-adaptive.fsim");
    let fixture = Fixture::new(source);
    let pin = fs_blake3::hash_bytes(source.as_bytes()).to_hex();
    let mut binder = fixture.binder();
    let first = captured(&mut binder, &pin);
    std::fs::rename(&fixture.source, fixture.root.join("moved-3d.fsim")).unwrap();
    let second = captured(&mut binder, &pin);
    assert_eq!(first.source_hash(), second.source_hash());
    assert!(!fixture.ledger.exists());
}
