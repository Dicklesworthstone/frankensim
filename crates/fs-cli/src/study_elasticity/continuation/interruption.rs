//! Exercise the real continuation, kernel callbacks and sealed ledger together.
use super::*;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-2d.fsim"));

fn spec() -> ElasticitySpec {
    parse(&FIXTURE.replace(":mesh-level 4", ":mesh-level 3")
        .replace(":max-iterations 32", ":max-iterations 2")
        .replace(":steps 32", ":steps 2")
        .replace(":move-cells 0.35", ":move-cells 0.05")
        .replace(":hole-radius-cells 1.5", ":hole-radius-cells 0.5")
        .replace(":nucleation-period 4", ":nucleation-period 2"))
        .expect("admitted study")
}

fn interrupted_then_resumed(stop: impl Fn(CheckpointStage) -> bool) {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let prefix = drive(&spec, &ledger, Some(1), &gate, None).unwrap();
    let prefix = load(&ledger, &prefix.pointer).unwrap();
    let geometry = linked(&ledger, &prefix.value, "design", "study-design").unwrap();
    let rows = linked(&ledger, &prefix.value, "iterations", "study-iterations").unwrap();
    let mut reached = false;
    let cancelled = drive_observed(&spec, &ledger, None, &gate, Some(&prefix), |ordinal, stage| {
        assert_eq!(ordinal, 1);
        if stop(stage) { reached = true; gate.request(); }
    }).unwrap();
    assert!(reached, "a real requested kernel boundary must execute");
    assert_eq!(cancelled.status, "cancelled");
    let cancelled = load(&ledger, &cancelled.pointer).unwrap();
    assert_eq!(integer(&cancelled.value, "iterations_completed").unwrap(), 1);
    assert_eq!(linked(&ledger, &cancelled.value, "design", "study-design").unwrap(), geometry);
    assert_eq!(linked(&ledger, &cancelled.value, "iterations", "study-iterations").unwrap(), rows);
    assert_eq!(cancelled.value.path(&["continuation", "updates_this_invocation"])
        .and_then(JsonValue::as_f64), Some(0.0));
    assert!(cancelled.value.f64_field("consumed_wall_s").unwrap()
        >= prefix.value.f64_field("consumed_wall_s").unwrap());

    let resumed = drive(&spec, &ledger, None, &CancelGate::new_clock_free(), Some(&cancelled)).unwrap();
    assert_eq!(resumed.status, "constraint-unmet");
    let resumed = load(&ledger, &resumed.pointer).unwrap();
    let oracle_ledger = Ledger::open(":memory:").unwrap();
    let oracle = drive(&spec, &oracle_ledger, None, &CancelGate::new_clock_free(), None).unwrap();
    let oracle = load(&oracle_ledger, &oracle.pointer).unwrap();
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &resumed.value, key, kind).unwrap(),
            linked(&oracle_ledger, &oracle.value, key, kind).unwrap());
    }
    assert_eq!(resumed.value.str_field("status"), oracle.value.str_field("status"));
    assert_eq!(resumed.value.str_field("trace_hash"), oracle.value.str_field("trace_hash"));
    assert_eq!(linked(&ledger, &prefix.value, "design", "study-design").unwrap(), geometry);
}

#[test]
fn cancellation_inside_initial_cg_preserves_prior_receipt_and_resumes_exactly() {
    interrupted_then_resumed(|stage| matches!(stage, CheckpointStage::InitialSolve(n) if n > 0));
}

#[test]
fn cancellation_inside_candidate_cg_discards_evolved_geometry() {
    interrupted_then_resumed(|stage| matches!(stage, CheckpointStage::CandidateSolve(n) if n > 0));
}

#[test]
fn cancellation_after_solved_candidate_before_publication_retains_prior_design() {
    interrupted_then_resumed(|stage| stage == CheckpointStage::Publish);
}

#[test]
fn wall_gate_uses_accumulated_charge_and_cancellation_has_stable_precedence() {
    assert_eq!(stop_status(false, 4.0 + 5.0, 10.0), None);
    assert_eq!(stop_status(false, 4.0 + 6.0, 10.0), Some("budget-exhausted"));
    assert_eq!(stop_status(false, 4.0 + 7.0, 10.0), Some("budget-exhausted"));
    assert_eq!(stop_status(true, 11.0, 10.0), Some("cancelled"));
}
