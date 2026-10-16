use super::*;
use crate::study::elasticity::parse as parse_study;
use fs_topols::resolution::ResolutionStage;
use fs_topols::evaluated::DesignEvaluationStage;

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-mesh-checked-2d.fsim"));
fn strict() -> ElasticitySpec {
    parse_study(&EXAMPLE.replace(":absolute-compliance-tolerance-j 0.000001", ":absolute-compliance-tolerance-j 0.00000000000001")
        .replace(":relative-compliance-tolerance 0.05", ":relative-compliance-tolerance 0.0")).unwrap()
}
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn constraints(out: &Outcome) -> JsonValue { json(out).path(&["continuation", "constraints"]).unwrap().clone() }

#[test]
fn mesh_checked_native_declaration_is_canonical_explicit_and_changes_identity() {
    let spec = parse_study(EXAMPLE).unwrap();
    assert_eq!(parse_study(&spec.canonical).unwrap().id, spec.id);
    for (old, new) in [(":extra-levels 1", ":extra-levels 0"),
        (":extra-levels 1", ":extra-levels 3"),
        (":relative-compliance-tolerance 0.05", ":relative-compliance-tolerance 1.0"),
        (":relative-compliance-tolerance 0.05", ""),
        (":extra-levels 1", ":extra-levels 1 :extra-levels 1")] {
        assert!(parse_study(&EXAMPLE.replace(old, new)).is_err());
    }
    let plain = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/marquee/bracket-projected-volume-2d.fsim"));
    let plain = parse_study(plain).unwrap();
    assert!(!plain.canonical.contains("mesh-check"));
    assert_ne!(plain.id, spec.id);
}

#[test]
fn mesh_unresolved_native_baseline_is_measured_retained_and_not_a_success() {
    let spec = strict();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = drive_observed(&spec, &ledger, None, &gate(), None, |stage| {
        if let VolumeStage::Mesh(stage) = stage {
            assert!(matches!(stage, MeshCheckStage::Baseline(_)));
        }
    }).unwrap();
    assert_eq!(out.status, "mesh-unresolved");
    let value = json(&out);
    assert_eq!(integer(&value, "iterations_completed").unwrap(), 0);
    let report = document(&linked(&ledger, &value, "report_json", "study-report-json").unwrap()).unwrap();
    let check = report.path(&["constraints", "mesh_resolution", "last_check"]).unwrap();
    assert_eq!(check.str_field("outcome"), Some("baseline-unresolved"));
    assert_eq!(check.path(&["baseline", "rungs"]).and_then(JsonValue::as_array).unwrap().len(), 2);
    assert!(check.path(&["baseline", "max_compliance_change_j"]).and_then(JsonValue::as_f64).unwrap() > 1e-14);
    let old = load(&ledger, &out.pointer).unwrap();
    let again = drive_observed(&spec, &ledger, None, &gate(), Some(&old), |_| panic!("terminal recovery must not solve")).unwrap();
    assert_eq!(again.receipt, out.receipt);
    assert_eq!(render(out, OutputMode::Json).exit_code, exit::BUDGET);
}

#[test]
fn mesh_checked_cancellation_in_fine_solve_retains_seed_and_exact_retry() {
    let spec = strict();
    let ledger = Ledger::open(":memory:").unwrap();
    let stop = gate();
    let out = drive_observed(&spec, &ledger, None, &stop, None, |stage| {
        if matches!(stage, VolumeStage::Mesh(MeshCheckStage::Baseline(ResolutionStage::Evaluate {
            level: 4, stage: DesignEvaluationStage::Solve(n) })) if n > 0) { stop.request(); }
    }).unwrap();
    assert_eq!(out.status, "cancelled");
    assert_eq!(constraints(&out).path(&["mesh_resolution", "last_check"]), Some(&JsonValue::Null));
    let old = load(&ledger, &out.pointer).unwrap();
    let retried = drive(&spec, &ledger, None, &gate(), Some(&old)).unwrap();
    let other = Ledger::open(":memory:").unwrap();
    let whole = drive(&spec, &other, None, &gate(), None).unwrap();
    assert_eq!(retried.status, "mesh-unresolved");
    assert_eq!(constraints(&retried), constraints(&whole));
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &json(&retried), key, kind).unwrap(), linked(&other, &json(&whole), key, kind).unwrap());
    }
    assert_eq!(load(&ledger, &out.pointer).unwrap().bytes, out.receipt);
}

#[test]
fn mesh_checked_public_study_and_report_export_the_actual_resolution_rows() {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-mesh-check-{}-{nonce}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let source = dir.join("study.fsim");
    let db = dir.join("study.db");
    fs::write(&source, strict().canonical).unwrap();
    let invoke = |args: &[&str]| crate::run(args.iter().map(|s| s.to_string()));
    let out = invoke(&["--json", "study", source.to_str().unwrap(), db.to_str().unwrap()]);
    assert_eq!(out.exit_code, exit::BUDGET, "{}", out.stderr);
    let value = JsonValue::parse(&out.stdout).unwrap();
    assert_eq!(value.str_field("status"), Some("mesh-unresolved"));
    let pointer = value.str_field("run_id").unwrap();
    let exported = invoke(&["--json", "report", pointer, db.to_str().unwrap()]);
    assert_eq!(exported.exit_code, exit::SUCCESS, "{}", exported.stderr);
    let paths = JsonValue::parse(&exported.stdout).unwrap();
    let html = fs::read_to_string(paths.str_field("report_html").unwrap()).unwrap();
    assert!(html.contains("Observed mesh-resolution check"));
    assert!(html.contains("not a continuum error bound"));
    let repeated = invoke(&["--json", "study", "--resume", pointer, db.to_str().unwrap()]);
    assert_eq!(repeated.exit_code, out.exit_code);
    assert_eq!(repeated.stdout, out.stdout);
}
