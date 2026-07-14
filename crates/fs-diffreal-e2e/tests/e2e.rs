//! The differentiation & reality end-to-end battery (addendum, Layer-3
//! conformance). Asserts each stage's load-bearing behavior and the report's
//! fail-closed truth table: adjoint-vs-FD agreement + missing-VJP blocking, the
//! as-built loop (defect localization + misfit reduction), tolerance
//! allocation, and the honestly-gated required spacetime stage.

use fs_diffreal_e2e::{
    AS_BUILT_EVIDENCE_IDENTITY, DIFFERENTIATION_EVIDENCE_IDENTITY, DiffRealError,
    SPACETIME_EVIDENCE_IDENTITY, SPACETIME_STAGE, StageLog, StageReason, StageRequirement,
    StageStatus, TOLERANCE_EVIDENCE_IDENTITY, differentiate_path, run_battery, stage_as_built_loop,
    stage_differentiation, stage_spacetime_gated, stage_tolerance_allocation,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x6469_6666_7265_616c,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn active_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx(&CancelGate::new(), f)
}

#[test]
fn the_full_layer3_battery_reports_its_required_gate_fail_closed() {
    let report = active_cx(run_battery).expect("battery succeeds");
    assert!(!report.complete(), "a required gated stage is incomplete");
    assert!(
        !report.all_required_passed(),
        "a required gated stage did not pass"
    );
    assert!(
        !report.promotion_ready(),
        "a required gated stage cannot authorize promotion"
    );
    assert!(
        !report.passed(),
        "the compatibility predicate must also fail closed"
    );
    assert_eq!(report.stages().len(), 4);
    for stage in &report.stages()[..3] {
        assert!(stage.passed(), "stage {} failed: {stage}", stage.stage);
    }
    let spacetime = report
        .stage(SPACETIME_STAGE)
        .expect("required spacetime record exists");
    assert!(matches!(spacetime.status, StageStatus::Gated(_)));
    for stage in report.stages() {
        assert!(!stage.events.is_empty(), "stage {} has no log", stage.stage);
    }
}

#[test]
fn differentiation_agrees_with_fd_and_blocks_a_missing_vjp() {
    let s = stage_differentiation();
    assert!(s.passed(), "{:#?}", s.events);
    assert_eq!(s.status, StageStatus::Passed);
    assert_eq!(s.evidence_identity, DIFFERENTIATION_EVIDENCE_IDENTITY);
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
    let s = active_cx(stage_as_built_loop).expect("as-built stage succeeds");
    assert!(s.passed(), "{:#?}", s.events);
    assert_eq!(s.status, StageStatus::Passed);
    assert_eq!(s.evidence_identity, AS_BUILT_EVIDENCE_IDENTITY);
    // The seeded defect (0.3 at idx 1) is localized without upgrading the
    // calibration candidate beyond Estimated.
    assert!(
        s.events
            .iter()
            .any(|e| e.contains("idx Some(1)") && e.contains("estimated=true"))
    );
    // assimilation reduced the misfit.
    assert!(s.events.iter().any(|e| e.contains("assimilation misfit")));
}

#[test]
fn tolerance_allocation_tightens_high_sensitivity_and_loosens_low() {
    let s = stage_tolerance_allocation();
    assert!(s.passed(), "{:#?}", s.events);
    assert_eq!(s.status, StageStatus::Passed);
    assert_eq!(s.evidence_identity, TOLERANCE_EVIDENCE_IDENTITY);
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
    assert!(!s.passed());
    let reason = match &s.status {
        StageStatus::Gated(reason) => reason,
        other => panic!("spacetime must be gated, got {other:?}"),
    };
    assert_eq!(reason.code, "diffreal.spacetime.integration-not-activated");
    assert!(reason.detail.contains("not integrated and activated"));
    assert_eq!(s.evidence_identity, SPACETIME_EVIDENCE_IDENTITY);
    assert!(s.events.iter().any(|e| e.contains("GATED")));
    assert!(s.events.iter().any(|e| e.contains("bk0o.7 is shipped")));
}

#[test]
fn status_display_is_stable_and_distinguishes_failure_from_unavailability() {
    let failed = StageStatus::Failed(StageReason::new(
        "test.failed",
        "an evaluated assertion was false",
    ));
    let gated = StageStatus::Gated(StageReason::new(
        "test.gated",
        "the capability was unavailable",
    ));
    let refused = StageStatus::Refused(StageReason::new(
        "test.refused",
        "the input was inadmissible",
    ));

    assert_eq!(failed.code(), "failed");
    assert_eq!(gated.code(), "gated");
    assert_eq!(refused.code(), "refused");
    assert_eq!(
        failed.to_string(),
        "failed[test.failed]: an evaluated assertion was false"
    );
    assert_eq!(
        gated.to_string(),
        "gated[test.gated]: the capability was unavailable"
    );
    assert_eq!(
        refused.to_string(),
        "refused[test.refused]: the input was inadmissible"
    );
    assert_ne!(failed.to_string(), gated.to_string());

    let log = StageLog::new(
        "display-fixture",
        StageRequirement::Optional,
        gated,
        "display-fixture/v1",
        vec!["gate recorded".to_string()],
    );
    assert_eq!(
        log.to_string(),
        "stage=display-fixture requirement=optional status=gated[test.gated]: the capability was unavailable evidence_identity=display-fixture/v1"
    );
}

#[test]
fn the_battery_is_deterministic() {
    let first = active_cx(run_battery).expect("first battery succeeds");
    let second = active_cx(run_battery).expect("replay succeeds");
    assert_eq!(first, second);
}

#[test]
fn cancellation_propagates_without_a_partial_battery() {
    let gate = CancelGate::new();
    gate.request();
    let result = with_cx(&gate, run_battery);
    assert!(
        matches!(
            result,
            Err(DiffRealError::AsBuilt(
                fs_asbuilt::RegError::Cancelled { .. }
            ))
        ),
        "pre-cancelled battery must propagate the structured refusal"
    );
}
