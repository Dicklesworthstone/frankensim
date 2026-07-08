//! The differentiation & reality end-to-end battery (addendum, Layer-3
//! conformance). Asserts the whole suite passes and each stage's load-bearing
//! behavior: adjoint-vs-FD agreement + missing-VJP blocking, the as-built loop
//! (defect localization + misfit reduction), tolerance allocation, and the
//! honestly-gated spacetime stage.

use fs_diffreal_e2e::{
    differentiate_path, run_battery, stage_as_built_loop, stage_differentiation,
    stage_spacetime_gated, stage_tolerance_allocation,
};

#[test]
fn the_full_layer3_battery_passes() {
    let report = run_battery();
    assert!(report.passed(), "battery failed: {report:#?}");
    assert_eq!(report.stages.len(), 4);
    for s in &report.stages {
        assert!(s.passed, "stage {} failed", s.stage);
        assert!(!s.events.is_empty());
    }
}

#[test]
fn differentiation_agrees_with_fd_and_blocks_a_missing_vjp() {
    let s = stage_differentiation();
    assert!(s.passed, "{:#?}", s.events);
    assert!(s.events.iter().any(|e| e.contains("agree=true")));
    assert!(s.events.iter().any(|e| e.contains("never silent-zero")));

    // a full-coverage path differentiates; a missing VJP blocks with a message.
    let full = |op: &str| matches!(op, "sdf" | "spline" | "solve");
    assert!(differentiate_path(&["sdf", "spline", "solve"], full, 1.0).is_ok());
    let blocked = differentiate_path(&["sdf", "remesh", "solve"], full, 1.0);
    assert!(blocked.is_err());
    assert!(blocked.unwrap_err().contains("remesh"));
}

#[test]
fn the_as_built_loop_localizes_a_defect_and_reduces_misfit() {
    let s = stage_as_built_loop();
    assert!(s.passed, "{:#?}", s.events);
    // the seeded defect (0.3 at idx 1) is localized and the δ is validated.
    assert!(
        s.events
            .iter()
            .any(|e| e.contains("idx Some(1)") && e.contains("validated=true"))
    );
    // assimilation reduced the misfit.
    assert!(s.events.iter().any(|e| e.contains("assimilation misfit")));
}

#[test]
fn tolerance_allocation_tightens_high_sensitivity_and_loosens_low() {
    let s = stage_tolerance_allocation();
    assert!(s.passed, "{:#?}", s.events);
    assert!(
        s.events
            .iter()
            .any(|e| e.contains("critical -> Tighten") && e.contains("slack -> Loosen"))
    );
    assert!(
        s.events
            .iter()
            .any(|e| e.contains("robustness confirmed = true"))
    );
}

#[test]
fn the_spacetime_stage_is_honestly_gated() {
    let s = stage_spacetime_gated();
    assert!(s.passed);
    assert!(s.events.iter().any(|e| e.contains("GATED")));
}

#[test]
fn the_battery_is_deterministic() {
    assert_eq!(run_battery(), run_battery());
}
