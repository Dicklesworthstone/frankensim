//! Real native volume-only studies: actual projection, elasticity, and ledger.
use super::*;
use fs_topols::evaluated::DesignEvaluationStage;

const ORIGINAL: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/marquee/bracket-2d.fsim"));
const POLICY: &str = "    :constraint-mode projected-volume\n    :area-tolerance-m2 0.0001\n    :max-projection-shift 2.0\n    :max-area-evaluations 64\n    :max-candidates 16\n    :contraction 0.5\n    :min-relative-improvement 0.00000001\n    :cg-poll-iters 1)";
fn source() -> String {
    ORIGINAL.replace(":mesh-level 4", ":mesh-level 3")
        .replace(":youngs-modulus-pa 70000000000.0", ":youngs-modulus-pa 2.0")
        .replace(":load-traction-pa 1000000.0", ":load-traction-pa 1.0")
        .replace(":volume-fraction 0.45", ":volume-fraction 0.75")
        .replace(":max-iterations 32", ":max-iterations 2")
        .replace(":move-cells 0.35", ":move-cells 0.1")
        .replace(":nucleation-period 4", ":nucleation-period 0")
        .replace("    :steps 32)", &format!("    :steps 2\n{POLICY}"))
}
fn spec() -> ElasticitySpec { parse(&source()).expect("admitted volume-only study") }
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn native(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>, prior: Option<&Loaded>) -> Outcome {
    super::super::super::drive(spec, ledger, cap, &gate(), prior).unwrap()
}
fn read(ledger: &Ledger, out: &Outcome, spec: &ElasticitySpec) -> (GridSdf, OptimizeReport, VolumeEvidence) {
    let receipt = json(out);
    let design = document(&linked(ledger, &receipt, "design", "study-design").unwrap()).unwrap();
    let rows = document(&linked(ledger, &receipt, "iterations", "study-iterations").unwrap()).unwrap();
    let (phi, report) = decode(spec, &receipt, &design, &rows).unwrap();
    let Some(ProjectedControls::Volume(policy)) = &spec.projected else { panic!("volume policy") };
    let history = VolumeEvidence::read(receipt.path(&["continuation", "constraints"]).unwrap(), &report, policy).unwrap();
    (phi, report, history)
}

#[test]
fn volume_only_policy_is_explicit_and_canonical_without_a_stress_sentinel() {
    let spec = spec();
    assert!(matches!(&spec.projected, Some(ProjectedControls::Volume(_))));
    assert!(!spec.canonical.contains("stress-limit"));
    let again = parse(&spec.canonical).unwrap();
    assert_eq!(again.id, spec.id);
    assert_eq!(again.canonical, spec.canonical);
    assert_ne!(spec.id, parse(ORIGINAL).unwrap().id);
    for field in POLICY.lines().skip(1) {
        assert!(parse(&source().replace(field.trim_end_matches(')'), "")).is_err(), "missing {field}");
    }
    for (from, to) in [
        (":cg-poll-iters 1)", ":cg-poll-iters 1 :sampled-stress-limit-pa 1000000000000.0)"),
        (":max-candidates 16", ":max-candidates 0"),
        (":area-tolerance-m2 0.0001", ":area-tolerance-m2 0.1"),
        (":contraction 0.5", ":contraction 1.0"),
        (":cg-poll-iters 1", ":cg-poll-iters 0"),
        (":constraint-mode projected-volume", ":constraint-mode projected-volume :constraint-mode projected-volume"),
    ] { assert!(parse(&source().replace(from, to)).is_err(), "must reject {to}"); }
}

#[test]
fn same_area_baseline_and_accepted_geometry_are_real_and_resume_bit_exactly() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let first = native(&spec, &ledger, Some(1), None);
    assert_eq!(first.status, "budget-exhausted", "{}", first.receipt);
    let (phi, report, history) = read(&ledger, &first, &spec);
    assert_eq!(report.rows.len(), 1, "must perform a genuine accepted design update");
    assert!(report.compliance[0] < history.baseline.compliance);
    assert!((history.baseline.volume - 0.75).abs() <= 1e-4);
    assert!((report.volume[0] - 0.75).abs() <= 1e-4);
    assert_ne!(history.baseline.snapshot, snapshot(&initial_phi(&spec)), "project the baseline before claiming improvement");
    let oracle = fs_topols::evaluate_compliance_design(&phi, fixture(&spec), settings(&spec, spec.steps)).unwrap();
    assert!(Measured::from(oracle).same(history.current()));
    let untouched = linked(&ledger, &json(&first), "design", "study-design").unwrap();
    let old = load(&ledger, &first.pointer).unwrap();
    let resumed = native(&spec, &ledger, None, Some(&old));
    let full_ledger = Ledger::open(":memory:").unwrap();
    let full = native(&spec, &full_ledger, None, None);
    assert_eq!(resumed.status, full.status);
    let (_, _, resumed_history) = read(&ledger, &resumed, &spec);
    let (_, _, full_history) = read(&full_ledger, &full, &spec);
    assert_eq!(resumed_history.json(), full_history.json());
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &json(&resumed), key, kind).unwrap(),
            linked(&full_ledger, &json(&full), key, kind).unwrap());
    }
    assert_eq!(linked(&ledger, &json(&first), "design", "study-design").unwrap(), untouched);
    assert_eq!(json(&resumed).path(&["continuation", "legacy_prefix_updates_replayed"]).and_then(JsonValue::as_f64), Some(0.0));
    let summary = document(&linked(&ledger, &json(&resumed), "report_json", "study-report-json").unwrap()).unwrap();
    assert_eq!(summary.path(&["constraints", "stress_evaluation"]).and_then(JsonValue::as_str), Some("not-requested"));
    assert!(summary.path(&["constraints", "stress_limit_pa"]).is_none());
    assert!(summary.path(&["constraints", "baseline", "sample_count"]).is_none());
}

#[test]
fn interruption_after_candidate_solve_retains_feasible_seed_and_retries_exactly() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let cancelled = gate();
    let mut reached = false;
    let out = drive_observed(&spec, &ledger, None, &cancelled, None, |stage| {
        if matches!(stage, VolumeStage::Update(ProjectedStage::Publish(_))) {
            reached = true; cancelled.request();
        }
    }).unwrap();
    assert!(reached, "stop at a real candidate publication boundary");
    assert_eq!(out.status, "cancelled");
    let (phi, report, history) = read(&ledger, &out, &spec);
    assert!(report.rows.is_empty());
    assert_eq!(snapshot(&phi), history.baseline.snapshot);
    let old = load(&ledger, &out.pointer).unwrap();
    let resumed = native(&spec, &ledger, Some(1), Some(&old));
    let full_ledger = Ledger::open(":memory:").unwrap();
    let full = native(&spec, &full_ledger, Some(1), None);
    assert_eq!(integer(&json(&resumed), "iterations_completed").unwrap(), 1);
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &json(&resumed), key, kind).unwrap(),
            linked(&full_ledger, &json(&full), key, kind).unwrap());
    }
}

#[test]
fn setup_cancellation_cannot_publish_an_unevaluated_feasible_baseline() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let cancelled = gate();
    let mut reached = false;
    let error = drive_observed(&spec, &ledger, None, &cancelled, None, |stage| {
        if matches!(stage, VolumeStage::Setup(ProjectedSetupStage::Evaluation(DesignEvaluationStage::Solve(n))) if n > 0) {
            reached = true; cancelled.request();
        }
    }).unwrap_err();
    assert!(reached);
    assert_eq!(error.exit, exit::CANCELLED);
    assert_eq!(ledger.table_count("ops").unwrap(), 0);
}

#[test]
fn bounded_failed_search_preserves_the_feasible_design_and_terminal_receipt() {
    let spec = parse(&source().replace(":min-relative-improvement 0.00000001", ":min-relative-improvement 0.999999")
        .replace(":max-candidates 16", ":max-candidates 2")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = native(&spec, &ledger, None, None);
    assert_eq!(out.status, "no-feasible-descent");
    let (phi, report, history) = read(&ledger, &out, &spec);
    assert!(report.rows.is_empty());
    assert_eq!(snapshot(&phi), history.baseline.snapshot);
    assert_eq!(history.refusals.len(), 2);
    let old = load(&ledger, &out.pointer).unwrap();
    let again = native(&spec, &ledger, None, Some(&old));
    assert_eq!(out.pointer, again.pointer);
    assert_eq!(out.receipt, again.receipt);
    assert_eq!(render(out, OutputMode::Json).exit_code, exit::REFUSED);
}

#[test]
fn retained_area_or_baseline_mutations_refuse_before_recovery_physics() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = native(&spec, &ledger, Some(1), None);
    let (_, report, _) = read(&ledger, &out, &spec);
    let Some(ProjectedControls::Volume(policy)) = &spec.projected else { panic!("volume") };
    let mut value = json(&out).path(&["continuation", "constraints"]).unwrap().clone();
    let JsonValue::Object(ref mut fields) = value else { panic!("constraint object") };
    fields.iter_mut().find(|(key, _)| key == "area_target_m2").unwrap().1 = JsonValue::Number { value: 0.9, raw: "0.9".into() };
    assert!(VolumeEvidence::read(&value, &report, policy).is_err());
    let mut value = json(&out).path(&["continuation", "constraints"]).unwrap().clone();
    let JsonValue::Object(ref mut fields) = value else { panic!("constraints") };
    fields.iter_mut().find(|(key, _)| key == "relative_reduction").unwrap().1 = JsonValue::Number { value: 0.99, raw: "0.99".into() };
    assert!(VolumeEvidence::read(&value, &report, policy).is_err());
}

#[test]
fn public_study_resume_and_exports_use_the_volume_only_example_and_retained_results() {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};
    const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-projected-volume-2d.fsim"));
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-volume-native-{}-{nonce}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let path = dir.join("study.fsim");
    let db = dir.join("study.db");
    fs::write(&path, EXAMPLE).unwrap();
    let invoke = |args: &[&str]| crate::run(args.iter().map(|value| value.to_string()));
    let first = invoke(&["--json", "study", path.to_str().unwrap(), db.to_str().unwrap(), "--budget", "1"]);
    assert_eq!(first.exit_code, exit::BUDGET, "{}", first.stderr);
    let first_json = JsonValue::parse(&first.stdout).unwrap();
    assert_eq!(first_json.path(&["receipt", "iterations_completed"]).and_then(JsonValue::as_f64), Some(1.0));
    let pointer = first_json.str_field("run_id").unwrap();
    let resumed = invoke(&["--json", "study", "--resume", pointer, db.to_str().unwrap()]);
    assert!(matches!(resumed.exit_code, exit::SUCCESS | exit::REFUSED), "{}", resumed.stderr);
    let result = JsonValue::parse(&resumed.stdout).unwrap();
    assert!(matches!(result.str_field("status"), Some("completed" | "no-feasible-descent")));
    let pointer = result.str_field("run_id").unwrap();
    let report = invoke(&["--json", "report", pointer, db.to_str().unwrap()]);
    assert_eq!(report.exit_code, exit::SUCCESS, "{}", report.stderr);
    let exported = JsonValue::parse(&report.stdout).unwrap();
    let summary = document(&fs::read(exported.str_field("report_json").unwrap()).unwrap()).unwrap();
    assert_eq!(summary.path(&["constraints", "mode"]).and_then(JsonValue::as_str), Some("projected-volume-v1"));
    assert_eq!(summary.path(&["constraints", "stress_evaluation"]).and_then(JsonValue::as_str), Some("not-requested"));
    assert_eq!(summary.path(&["constraints", "area_constraint_satisfied"]), Some(&JsonValue::Bool(true)));
    assert!((summary.f64_field("final_material_area_m2").unwrap() - 0.75).abs() <= 1e-4);
    assert!(summary.path(&["constraints", "relative_reduction"]).and_then(JsonValue::as_f64).unwrap() > 0.0);
    let html = fs::read_to_string(exported.str_field("report_html").unwrap()).unwrap();
    assert!(html.contains("No stress limit or stress evaluation was requested"));
    let package = invoke(&["--json", "package", pointer, db.to_str().unwrap()]);
    assert_eq!(package.exit_code, exit::SUCCESS, "{}", package.stderr);
    let exported_package = JsonValue::parse(&package.stdout).unwrap();
    let package = EvidencePackage::from_json(&fs::read_to_string(exported_package.str_field("package").unwrap()).unwrap()).unwrap();
    assert!(fs_checker::check(&package).passed());
    let repeated = invoke(&["--json", "study", "--resume", pointer, db.to_str().unwrap()]);
    assert_eq!(repeated.exit_code, resumed.exit_code);
    assert_eq!(repeated.stdout, resumed.stdout, "sealed terminal has no extra solve or new receipt");
    fs::remove_dir_all(&dir).unwrap();
}
