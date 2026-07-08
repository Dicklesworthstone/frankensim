//! The epistemic type-system end-to-end battery (addendum, Layer-2
//! conformance). Asserts the whole suite passes and each stage's load-bearing
//! FAIL-CLOSED behavior holds, with structured per-stage logging.

use fs_epi_e2e::{
    run_battery, stage_evidence_roundtrip, stage_falsifier, stage_goodhart_guard, stage_laundering,
    stage_objective_epistemics,
};

#[test]
fn the_full_layer2_battery_passes() {
    let report = run_battery();
    assert!(report.passed(), "battery failed: {report:#?}");
    // all five stages present, each with logged events.
    assert_eq!(report.stages.len(), 5);
    for s in &report.stages {
        assert!(s.passed, "stage {} failed", s.stage);
        assert!(!s.events.is_empty(), "stage {} logged nothing", s.stage);
    }
}

#[test]
fn laundering_fails_closed() {
    // composition cannot upgrade a color; out-of-regime demotes; in-regime kept.
    let s = stage_laundering();
    assert!(s.passed, "{:#?}", s.events);
    assert!(s.events.iter().any(|e| e.contains("no laundering")));
    assert!(s.events.iter().any(|e| e.contains("demotion=true")));
}

#[test]
fn the_no_falsifier_no_ship_gate_blocks() {
    let s = stage_falsifier();
    assert!(s.passed, "{:#?}", s.events);
    assert!(s.events.iter().any(|e| e.contains("no-falsifier-no-ship")));
}

#[test]
fn the_goodhart_guard_refuses_exploits_but_honors_genuine_optima() {
    let s = stage_goodhart_guard();
    assert!(s.passed, "{:#?}", s.events);
    // the exploit is refused (not honored) and the smooth optimum is honored.
    assert!(
        s.events
            .iter()
            .any(|e| e.contains("exploit") && e.contains("honored=false"))
    );
    assert!(
        s.events
            .iter()
            .any(|e| e.contains("smooth") && e.contains("honored=true"))
    );
    assert!(s.events.iter().any(|e| e.contains("provisional")));
}

#[test]
fn objective_epistemics_holds_the_contract_and_weakest_input_rule() {
    let s = stage_objective_epistemics();
    assert!(s.passed, "{:#?}", s.events);
    assert!(s.events.iter().any(|e| e.contains("refused = true")));
}

#[test]
fn the_evidence_package_round_trips_and_tamper_is_caught() {
    let s = stage_evidence_roundtrip();
    assert!(s.passed, "{:#?}", s.events);
    assert!(
        s.events
            .iter()
            .any(|e| e.contains("tampered package caught = true"))
    );
    assert!(s.events.iter().any(|e| e.contains("no solver")));
}

#[test]
fn the_battery_is_deterministic() {
    assert_eq!(run_battery(), run_battery());
}
