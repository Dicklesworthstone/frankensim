//! Native weighted DWR regressions using real elasticity and ledger consumers.
use super::*;
use crate::study::elasticity::parse as study_spec;
use fs_topols::robust_resolution::goal::MultiLoadDwrStage;
use fs_topols::WeightedComplianceDwrStage;

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-weighted-dwr-2d.fsim"));
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn family(out: &Outcome) -> JsonValue {
    json(out).path(&["continuation", "constraints", "load_family"]).unwrap().clone()
}
fn native(spec: &ElasticitySpec, ledger: &Ledger, prior: Option<&Loaded>) -> Outcome {
    super::super::drive(spec, ledger, None, &gate(), prior).unwrap()
}
fn recovered(ledger: &Ledger, out: &Outcome, spec: &ElasticitySpec) -> ConstraintEvidence {
    let receipt = json(out);
    let design = document(&linked(ledger, &receipt, "design", "study-design").unwrap()).unwrap();
    let rows = document(&linked(ledger, &receipt, "iterations", "study-iterations").unwrap()).unwrap();
    let (_, report) = decode(spec, &receipt, &design, &rows).unwrap();
    ConstraintEvidence::read(receipt.path(&["continuation", "constraints"]).unwrap(),
        &report, spec.projected.as_ref().unwrap()).unwrap()
}

#[test]
fn weighted_dwr_native_policy_is_explicit_and_single_load_stress_is_supported() {
    let spec = study_spec(EXAMPLE).unwrap();
    assert_eq!(study_spec(&spec.canonical).unwrap().id, spec.id);
    for (from, to) in [
        (":max-error-fraction 0.5", ":max-error-fraction 0.0"),
        (":max-error-fraction 0.5", ":max-error-fraction 1.0"),
        (":max-error-fraction 0.5", ""),
        (":max-error-fraction 0.5", ":max-error-fraction 0.5 :max-error-fraction 0.5"),
        (":aggregate weighted-sum", ":aggregate worst-weighted-case"),
    ] { assert!(study_spec(&EXAMPLE.replace(from, to)).is_err(), "{to}"); }
    let single = EXAMPLE.replace("        (case :band (0.375 0.625) :traction-pa (0.0 1.0) :weight 1.0)\n", "");
    let single = study_spec(&single).unwrap();
    assert_eq!(stress_controls(&single).unwrap().family.as_ref().unwrap().cases(&single).unwrap().len(), 1);
    assert_eq!(study_spec(&single.canonical).unwrap().id, single.id);
    let plain = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-projected-stress-2d.fsim"));
    assert!(stress_controls(&study_spec(plain).unwrap()).unwrap().family.is_none());
}

#[test]
fn weighted_dwr_native_real_baseline_reports_each_case_and_exhausts_without_false_success() {
    let spec = study_spec(&EXAMPLE.replace(":max-solves 512", ":max-solves 6")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = native(&spec, &ledger, None);
    assert_eq!(out.status, "budget-exhausted", "{}", out.receipt);
    assert_eq!(integer(&json(&out), "iterations_completed").unwrap(), 0);
    let history = family(&out);
    assert_eq!(integer(&history, "solves_started").unwrap(), 6);
    let check = history.path(&["goal_resolution", "last_check"]).unwrap();
    assert_eq!(check.str_field("outcome"), Some("solve-budget"));
    let cases = check.get("baseline").and_then(JsonValue::as_array).unwrap();
    assert_eq!(cases.len(), 2);
    for row in cases {
        assert!(row.f64_field("coarse_compliance_j").unwrap() > 0.0);
        assert!(row.f64_field("enriched_compliance_j").unwrap() > 0.0);
        assert!(row.f64_field("coarse_relative_residual").unwrap() < 1e-12);
        assert!(row.f64_field("enriched_relative_residual").unwrap() < 1e-12);
    }
    assert!((cases[0].f64_field("coarse_compliance_j").unwrap()
        / cases[1].f64_field("coarse_compliance_j").unwrap() - 1.0).abs() < 1e-8);
    let report = document(&linked(&ledger, &json(&out), "report_json", "study-report-json").unwrap()).unwrap();
    assert_eq!(report.path(&["constraints", "load_family", "goal_resolution"]), history.get("goal_resolution"));
    let html = String::from_utf8(linked(&ledger, &json(&out), "report_html", "study-report-html").unwrap()).unwrap();
    assert!(html.contains("Independent-load DWR gain check"));
    let old = load(&ledger, &out.pointer).unwrap();
    let again = native(&spec, &ledger, Some(&old));
    assert_eq!(again.receipt, out.receipt);
    assert_eq!(again.pointer, out.pointer);
    assert_eq!(render(out, OutputMode::Json).exit_code, exit::BUDGET);
}

#[test]
fn weighted_dwr_native_memory_admission_happens_before_geometry_and_ledger_work() {
    let spec = study_spec(&EXAMPLE.replace(":memory 536870912 B", ":memory 134217728 B")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let error = super::super::drive(&spec, &ledger, None, &gate(), None).unwrap_err();
    assert!(error.message.contains("256 MiB"));
    assert_eq!(ledger.table_count("ops").unwrap(), 0);
}

#[test]
fn weighted_dwr_native_cancel_resume_keeps_geometry_and_retains_dwr_work_charges() {
    let spec = study_spec(EXAMPLE).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let stop = gate();
    let mut reached = false;
    let out = super::super::driver::drive_goal_observed(&spec, &ledger, None, &stop, None,
        |_| {}, |_| {}, |_| {}, |stage| {
            if matches!(stage, MultiLoadDwrStage::Baseline(WeightedComplianceDwrStage::BeforeCase { case: 1 })) {
                reached = true; stop.request();
            }
        }).unwrap();
    assert!(reached, "must interrupt actual DWR work after the first independent case");
    assert_eq!(out.status, "cancelled");
    assert_eq!(integer(&json(&out), "iterations_completed").unwrap(), 0);
    let history = family(&out);
    assert_eq!(integer(&history, "solves_started").unwrap(), 4);
    assert_eq!(history.path(&["goal_resolution", "last_check"]), Some(&JsonValue::Null));
    let original_design = linked(&ledger, &json(&out), "design", "study-design").unwrap();
    let old = load(&ledger, &out.pointer).unwrap();
    let resumed = native(&spec, &ledger, Some(&old));
    let clean_db = Ledger::open(":memory:").unwrap();
    let clean = native(&spec, &clean_db, None);
    assert_eq!(resumed.status, clean.status);
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &json(&resumed), key, kind).unwrap(),
            linked(&clean_db, &json(&clean), key, kind).unwrap());
    }
    let a = family(&resumed); let b = family(&clean);
    assert_eq!(a.get("goal_resolution"), b.get("goal_resolution"));
    assert_eq!(integer(&a, "solves_started").unwrap(), integer(&b, "solves_started").unwrap() + 2);
    assert_eq!(integer(&a, "recovery_solves_used").unwrap(), 4);
    assert_eq!(linked(&ledger, &json(&out), "design", "study-design").unwrap(), original_design);
}

fn replace(value: &mut JsonValue, name: &str, new: JsonValue) {
    let JsonValue::Object(fields) = value else { panic!("object") };
    fields.iter_mut().find(|(key, _)| key == name).unwrap().1 = new;
}
#[test]
fn weighted_dwr_native_recovery_checks_actual_case_bindings_and_residuals() {
    let spec = study_spec(&EXAMPLE.replace(":max-solves 512", ":max-solves 6")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = native(&spec, &ledger, None);
    let evidence = recovered(&ledger, &out, &spec);
    let retained = evidence.family.as_ref().unwrap();
    let original = family(&out);
    let valid = read(&original, &evidence.policy, &retained.cases, &retained.baseline, &retained.accepted).unwrap();
    assert!(valid.is_some());
    for key in ["coarse_compliance_j", "coarse_relative_residual", "level"] {
        let mut row = original.path(&["goal_resolution", "last_check", "baseline"])
            .and_then(JsonValue::as_array).unwrap().to_vec();
        replace(&mut row[0], key, JsonValue::Number { value: 999.0, raw: "999.0".into() });
        let mut check = original.path(&["goal_resolution", "last_check"]).unwrap().clone();
        replace(&mut check, "baseline", JsonValue::Array(row));
        let mut goal = original.get("goal_resolution").unwrap().clone();
        replace(&mut goal, "last_check", check);
        let mut changed = original.clone(); replace(&mut changed, "goal_resolution", goal);
        assert!(read(&changed, &evidence.policy, &retained.cases, &retained.baseline, &retained.accepted).is_err());
    }
}
