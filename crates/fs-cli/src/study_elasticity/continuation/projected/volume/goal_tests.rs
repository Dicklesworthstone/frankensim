use super::*;
use crate::study::elasticity::parse as parse_study;
use fs_topols::ComplianceDwrStage;

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-dwr-checked-2d.fsim"));
fn spec() -> ElasticitySpec { parse_study(EXAMPLE).unwrap() }
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn constraints(out: &Outcome) -> JsonValue {
    json(out).path(&["continuation", "constraints"]).unwrap().clone()
}
fn policy(spec: &ElasticitySpec) -> &Controls {
    let Some(ProjectedControls::Volume(policy)) = &spec.projected else { panic!("volume policy") };
    policy
}

#[test]
fn signed_dwr_native_policy_is_canonical_and_admits_the_entire_candidate_family() {
    let spec = spec();
    assert_eq!(parse_study(&spec.canonical).unwrap().id, spec.id);
    assert_eq!(policy(&spec).goal.unwrap().max_solves_per_attempt, 34);
    let plain = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-projected-volume-2d.fsim"));
    assert_ne!(parse_study(plain).unwrap().id, spec.id);
    for (old, new) in [
        (":max-solves-per-attempt 34", ":max-solves-per-attempt 33"),
        (":max-solves-per-attempt 34", ":max-solves-per-attempt 131"),
        (":max-error-fraction 0.5", ":max-error-fraction 0.0"),
        (":max-error-fraction 0.5", ":max-error-fraction 1.0"),
        (":max-error-fraction 0.5", ""),
        (":max-error-fraction 0.5", ":max-error-fraction 0.5 :max-error-fraction 0.5"),
        (":goal-check (elasticity-dwr", ":goal-check (refinement"),
        (":cg-poll-iters 1", ":cg-poll-iters 1 :mesh-check (refinement :extra-levels 1 :absolute-compliance-tolerance-j 1.0 :relative-compliance-tolerance 0.1 :area-tolerance-m2 0.01)"),
    ] {
        assert!(EXAMPLE.contains(old));
        assert!(parse_study(&EXAMPLE.replace(old, new)).is_err(), "must reject {new}");
    }
    assert!(parse_study(&EXAMPLE.replace(":max-solves-per-attempt 34", ":max-solves-per-attempt 35")).is_ok());
}

#[test]
fn signed_dwr_native_memory_admission_precedes_geometry_or_solver_work() {
    let source = EXAMPLE.replace(":memory 536870912 B", ":memory 268435455 B");
    let spec = parse_study(&source).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let error = drive_observed(&spec, &ledger, None, &gate(), None,
        |_| panic!("underfunded DWR study started geometry or physics")).unwrap_err();
    assert!(error.message.contains("256 MiB"));
    assert_eq!(ledger.table_count("ops").unwrap(), 0);
}

#[test]
fn signed_dwr_native_search_reports_actual_bound_estimates_and_replays_terminal_state() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = drive(&spec, &ledger, None, &gate(), None).unwrap();
    assert!(matches!(out.status, "completed" | "no-feasible-descent"), "{}", out.receipt);
    let value = json(&out);
    let history = constraints(&out);
    let check = history.path(&["goal_resolution", "last_check"]).unwrap();
    assert!(matches!(check.str_field("outcome"), Some("accepted" | "no-descent")));
    assert!(check.path(&["baseline", "coarse_relative_residual"]).and_then(JsonValue::as_f64).unwrap() < 1e-12);
    assert!(check.path(&["baseline", "enriched_relative_residual"]).and_then(JsonValue::as_f64).unwrap() < 1e-12);
    let summary = document(&linked(&ledger, &value, "report_json", "study-report-json").unwrap()).unwrap();
    assert_eq!(summary.path(&["constraints", "goal_resolution"]), history.get("goal_resolution"));
    let html = String::from_utf8(linked(&ledger, &value, "report_html", "study-report-html").unwrap()).unwrap();
    assert!(html.contains("Signed-DWR improvement check"));
    assert!(html.contains("not a continuum error bound"));
    let old = load(&ledger, &out.pointer).unwrap();
    let again = drive_observed(&spec, &ledger, None, &gate(), Some(&old),
        |_| panic!("terminal recovery must not start more physics")).unwrap();
    assert_eq!(again.receipt, out.receipt);
    if out.status == "no-feasible-descent" {
        assert_eq!(render(out, OutputMode::Json).exit_code, exit::REFUSED);
    }
}

#[test]
fn signed_dwr_native_cancellation_charges_work_preserves_seed_and_retries_exactly() {
    let spec = spec();
    let whole_ledger = Ledger::open(":memory:").unwrap();
    let whole = drive(&spec, &whole_ledger, Some(1), &gate(), None).unwrap();
    for target in [DwrCheckStage::Baseline(ComplianceDwrStage::BeforeEstimate),
        DwrCheckStage::Baseline(ComplianceDwrStage::BeforePublish),
        DwrCheckStage::Candidate { index: 0, stage: ComplianceDwrStage::BeforeEstimate },
        DwrCheckStage::Candidate { index: 0, stage: ComplianceDwrStage::BeforePublish }] {
        let ledger = Ledger::open(":memory:").unwrap();
        let stop = gate();
        let mut reached = false;
        let out = drive_observed(&spec, &ledger, Some(1), &stop, None, |stage| {
            let VolumeStage::Goal(stage) = stage else { return };
            let hit = stage == target || matches!((stage, target),
                (DwrCheckStage::Candidate { stage: a, .. }, DwrCheckStage::Candidate { stage: b, .. }) if a == b);
            if hit { reached = true; stop.request(); }
        }).unwrap();
        assert!(reached, "real boundary not reached: {target:?}");
        assert_eq!(out.status, "cancelled");
        assert_eq!(integer(&json(&out), "iterations_completed").unwrap(), 0);
        let retained = constraints(&out);
        let work = retained.get("goal_resolution").unwrap();
        assert_eq!(integer(work, "attempts_admitted").unwrap(), 1);
        assert_eq!(integer(work, "total_reserved_solves").unwrap(), 34);
        assert_eq!(work.get("last_check"), Some(&JsonValue::Null));
        // The worst-case debit was persisted before the actual attempt as well.
        let receipt = json(&out);
        let previous = load(&ledger, &format!("study-{}", receipt.str_field("predecessor").unwrap())).unwrap();
        assert_eq!(previous.value.path(&["continuation", "constraints", "goal_resolution", "total_reserved_solves"])
            .and_then(JsonValue::as_f64), Some(34.0));
        let old = load(&ledger, &out.pointer).unwrap();
        let retried = drive(&spec, &ledger, Some(1), &gate(), Some(&old)).unwrap();
        assert_eq!(retried.status, whole.status);
        for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
            assert_eq!(linked(&ledger, &json(&retried), key, kind).unwrap(),
                linked(&whole_ledger, &json(&whole), key, kind).unwrap());
        }
        let retried_constraints = constraints(&retried);
        assert_eq!(retried_constraints.path(&["goal_resolution", "last_check"]),
            constraints(&whole).path(&["goal_resolution", "last_check"]));
        assert_eq!(retried_constraints.path(&["goal_resolution", "attempts_admitted"])
            .and_then(JsonValue::as_f64), Some(2.0));
        assert_eq!(load(&ledger, &out.pointer).unwrap().bytes, out.receipt);
    }
}

fn replace(value: &mut JsonValue, path: &[&str], replacement: JsonValue) {
    if path.is_empty() { *value = replacement; return; }
    let JsonValue::Object(fields) = value else { panic!("object required") };
    let child = &mut fields.iter_mut().find(|(key, _)| key == path[0]).unwrap().1;
    replace(child, &path[1..], replacement);
}

#[test]
fn signed_dwr_native_recomputes_decisions_identity_residuals_and_work_on_read() {
    let spec = spec();
    let p = policy(&spec);
    let baseline = DwrMeasurement { snapshot: 1, level: 3, area: 0.75,
        coarse_compliance: 10.0, enriched_compliance: 11.0, eta_signed: 1.0,
        coarse_relative_residual: 1e-14, enriched_relative_residual: 1e-14 };
    let candidate = DwrMeasurement { snapshot: 2, coarse_compliance: 9.0,
        enriched_compliance: 10.0, ..baseline };
    let decision = compare_dwr_gains(baseline, candidate, p.goal.unwrap(), p.search.min_relative_improvement).unwrap();
    assert!(decision.refusal.is_none());
    let evidence = GoalEvidence { attempts: 1, failure: None, last: Some(LastCheck {
        outcome: "accepted", baseline, candidate: Some(candidate), decision: Some(decision), checked: 1, reserved: 4,
    }) };
    let previous = Measured { snapshot: 1, compliance: 10.0, volume: 0.75 };
    let current = Measured { snapshot: 2, compliance: 9.0, volume: 0.75 };
    let value = document(format!("{{\"goal_resolution\":{}}}", evidence.json(p.goal.unwrap())).as_bytes()).unwrap();
    assert!(read(&value, p, previous, current, 1).is_ok());
    for (path, raw) in [
        (vec!["total_reserved_solves"], "0"),
        (vec!["attempts_admitted"], "0"),
        (vec!["authority"], "\"Certified\""),
        (vec!["policy", "max_error_fraction"], "0.9"),
        (vec!["last_check", "decision", "allowance_j"], "1000.0"),
        (vec!["last_check", "reserved_solves"], "2"),
        (vec!["last_check", "baseline", "snapshot"], "\"0x0000000000000003\""),
        (vec!["last_check", "accepted_candidate", "eta_signed_j"], "999.0"),
        (vec!["last_check", "accepted_candidate", "enriched_relative_residual"], "0.001"),
        (vec!["last_check", "accepted_candidate"], "null"),
    ] {
        let mut edited = value.clone();
        let mut full_path = vec!["goal_resolution"];
        full_path.extend(path);
        replace(&mut edited, &full_path, document(raw.as_bytes()).unwrap());
        assert!(read(&edited, p, previous, current, 1).is_err(), "tamper admitted: {full_path:?}");
    }
    let mut undeclared = p.clone();
    undeclared.goal = None;
    assert!(read(&value, &undeclared, previous, current, 1).is_err());
    assert!(read(&document(b"{}").unwrap(), p, previous, current, 1).is_err());
}
