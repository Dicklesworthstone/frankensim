//! Bead 7tv.21.12: generated ASCENT optimizer-budget trend manifest.

#[path = "support/budget_trend.rs"]
mod budget_trend;

use budget_trend::{
    BBOB_COMPONENT, BUDGET_TREND_ROWS, BUDGET_TREND_SCHEMA, GRADIENT_COMPONENT,
    audit_budget_trend_manifest, canonical_budget_trend_manifest_json,
    gate_and_emit_budget_observation,
};

#[test]
fn canonical_budget_trend_manifest_is_complete_and_deterministic() {
    assert!(
        audit_budget_trend_manifest(&BUDGET_TREND_ROWS).is_empty(),
        "canonical trend manifest must admit"
    );
    assert_eq!(
        BUDGET_TREND_ROWS
            .iter()
            .filter(|row| row.suite == BBOB_COMPONENT)
            .count(),
        8
    );
    assert_eq!(
        BUDGET_TREND_ROWS
            .iter()
            .filter(|row| row.suite == GRADIENT_COMPONENT)
            .count(),
        6
    );
    let first = canonical_budget_trend_manifest_json();
    let second = canonical_budget_trend_manifest_json();
    assert_eq!(first, second);
    assert!(first.contains(BUDGET_TREND_SCHEMA));
    assert!(first.contains("\"authority\":\"regression-gate-declaration\""));
    assert_eq!(first.matches("\"kernel\":").count(), 14);
    assert_eq!(first.matches("\"ceiling\":").count(), 14);
    println!("{first}");
}

#[test]
fn duplicate_missing_and_metadata_drift_fail_closed() {
    let mut duplicate = BUDGET_TREND_ROWS.to_vec();
    let duplicate_key = duplicate[0];
    duplicate[1].suite = duplicate_key.suite;
    duplicate[1].kernel = duplicate_key.kernel;
    let duplicate_diagnostics = audit_budget_trend_manifest(&duplicate);
    assert!(
        duplicate_diagnostics
            .iter()
            .any(|diagnostic| diagnostic.contains("duplicate budget trend row"))
    );
    assert!(
        duplicate_diagnostics
            .iter()
            .any(|diagnostic| diagnostic.contains("is missing"))
    );

    let mut missing = BUDGET_TREND_ROWS.to_vec();
    missing.pop();
    assert!(
        audit_budget_trend_manifest(&missing)
            .iter()
            .any(|diagnostic| diagnostic.contains("sqp/shared-constrained is missing"))
    );

    let mut drifted = BUDGET_TREND_ROWS.to_vec();
    let invalid_ceiling = drifted[0].sanity_floor_exclusive;
    drifted[0].ceiling = invalid_ceiling;
    let drift_diagnostics = audit_budget_trend_manifest(&drifted);
    assert!(
        drift_diagnostics
            .iter()
            .any(|diagnostic| diagnostic.contains("does not exceed sanity floor"))
    );
    assert!(
        drift_diagnostics
            .iter()
            .any(|diagnostic| diagnostic.contains("metadata drifted"))
    );

    let mut reordered = BUDGET_TREND_ROWS.to_vec();
    reordered.reverse();
    assert!(
        audit_budget_trend_manifest(&reordered)
            .iter()
            .any(|diagnostic| diagnostic.contains("not in canonical"))
    );
}

#[test]
fn observation_gate_accepts_declared_boundary_and_catches_mutants() {
    let mut admitted = fs_obs::Emitter::new(BBOB_COMPONENT, BUDGET_TREND_SCHEMA);
    gate_and_emit_budget_observation(&mut admitted, BBOB_COMPONENT, "de/rastrigin2", 2_500, 3, 5);

    let over_ceiling = std::panic::catch_unwind(|| {
        let mut emitter = fs_obs::Emitter::new(BBOB_COMPONENT, BUDGET_TREND_SCHEMA);
        gate_and_emit_budget_observation(
            &mut emitter,
            BBOB_COMPONENT,
            "de/rastrigin2",
            2_501,
            3,
            5,
        );
    });
    assert!(over_ceiling.is_err(), "ceiling mutant must be caught");

    let impossible_successes = std::panic::catch_unwind(|| {
        let mut emitter = fs_obs::Emitter::new(BBOB_COMPONENT, BUDGET_TREND_SCHEMA);
        gate_and_emit_budget_observation(
            &mut emitter,
            BBOB_COMPONENT,
            "de/rastrigin2",
            2_500,
            6,
            5,
        );
    });
    assert!(
        impossible_successes.is_err(),
        "success-count mutant must be caught"
    );
}
