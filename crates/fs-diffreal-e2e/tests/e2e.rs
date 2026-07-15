//! Layer-3 differentiation-and-reality conformance. These tests falsify the
//! shared tape/VJP path independently, exercise typed hostile-input refusals,
//! and keep sampled linearization evidence separate from probability claims.

use fs_adjoint::transpose::Vjp;
use fs_diffreal_e2e::{
    AS_BUILT_EVIDENCE_IDENTITY, AS_BUILT_STAGE, DIFFERENTIATION_EVIDENCE_IDENTITY,
    DIFFERENTIATION_STAGE, DiffRealError, DifferentiationError, PRODUCTION_DIFFERENTIATION_PATH,
    SPACETIME_EVIDENCE_IDENTITY, SPACETIME_STAGE, StageEvent, StageLog, StageReason,
    StageRequirement, StageStatus, TOLERANCE_EVIDENCE_IDENTITY, TOLERANCE_STAGE,
    differentiate_path, production_vjp_registry, run_battery, stage_as_built_loop,
    stage_differentiation, stage_differentiation_with_registry, stage_spacetime_gated,
    stage_tolerance_allocation, stage_tolerance_allocation_with_samples, verify_sensitivity,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_toleralloc::Action;
use std::sync::Arc;

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
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
            budget,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn active_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx(&CancelGate::new(), Budget::INFINITE, f)
}

#[test]
fn the_full_layer3_battery_reports_its_required_gate_fail_closed() {
    let report = active_cx(run_battery).expect("battery succeeds");
    assert!(!report.complete(), "a required gated stage is incomplete");
    assert!(!report.all_required_passed());
    assert!(!report.promotion_ready());
    assert!(!report.passed());
    assert_eq!(report.stages().len(), 4);
    for stage in &report.stages()[..3] {
        assert!(stage.passed(), "stage {} failed: {stage}", stage.stage);
        assert!(!stage.events.is_empty());
    }
    let spacetime = report
        .stage(SPACETIME_STAGE)
        .expect("required spacetime record exists");
    assert!(matches!(spacetime.status, StageStatus::Gated(_)));
}

#[test]
fn production_tape_returns_the_fixture_primal_and_reverse_gradient() {
    active_cx(|cx| {
        let registry = production_vjp_registry().expect("fixed registry");
        let derivative = differentiate_path(&PRODUCTION_DIFFERENTIATION_PATH, &registry, 1.5, cx)
            .expect("covered path");
        assert_eq!(derivative.value().to_bits(), 16.0_f64.to_bits());
        assert_eq!(derivative.gradient().to_bits(), 16.0_f64.to_bits());

        let blocked = differentiate_path(&["sdf", "remesh", "solve"], &registry, 1.5, cx);
        assert!(matches!(
            blocked,
            Err(DifferentiationError::MissingVjp { ref op }) if op == "remesh"
        ));

        let stage = stage_differentiation(cx).expect("stage evaluates");
        assert_eq!(stage.status, StageStatus::Passed);
        assert_eq!(stage.evidence_identity, DIFFERENTIATION_EVIDENCE_IDENTITY);
        assert!(
            stage
                .events
                .iter()
                .any(|event| matches!(event, StageEvent::GradientVerified { .. }))
        );
        assert!(stage.events.iter().any(|event| matches!(
            event,
            StageEvent::MissingVjpProbe { op, blocked: true } if op == "remesh"
        )));
    });
}

#[derive(Debug)]
struct PerturbedSquareVjp;

impl Vjp for PerturbedSquareVjp {
    fn vjp(&self, primal_inputs: &[&[f64]], out_cotangent: &[f64]) -> Vec<Vec<f64>> {
        let primal = primal_inputs[0][0];
        let cotangent = out_cotangent[0];
        vec![vec![(2.0 * primal + 0.25) * cotangent]]
    }
}

#[derive(Debug)]
struct InfiniteSquareVjp;

impl Vjp for InfiniteSquareVjp {
    fn vjp(&self, _primal_inputs: &[&[f64]], _out_cotangent: &[f64]) -> Vec<Vec<f64>> {
        vec![vec![f64::INFINITY]]
    }
}

#[derive(Debug)]
struct CancellingSquareVjp {
    gate: Arc<CancelGate>,
}

impl Vjp for CancellingSquareVjp {
    fn vjp(&self, primal_inputs: &[&[f64]], out_cotangent: &[f64]) -> Vec<Vec<f64>> {
        self.gate.request();
        vec![vec![2.0 * primal_inputs[0][0] * out_cotangent[0]]]
    }
}

#[test]
fn a_perturbed_production_vjp_fails_the_e2e_gate() {
    active_cx(|cx| {
        let mut registry = production_vjp_registry().expect("fixed registry");
        registry
            .register("spline", PerturbedSquareVjp)
            .expect("bounded replacement");
        let stage = stage_differentiation_with_registry(cx, &registry).expect("evaluated stage");
        assert!(matches!(
            stage.status,
            StageStatus::Failed(ref reason)
                if reason.code == "diffreal.differentiation.production-rejected"
        ));
        assert!(stage.events.iter().any(|event| matches!(
            event,
            StageEvent::DifferentiationRejected {
                error: DifferentiationError::OracleDisagreement { .. }
            }
        )));
    });
}

#[test]
fn sealed_sensitivity_carries_dual_fd_step_and_unit_evidence() {
    active_cx(|cx| {
        let registry = production_vjp_registry().expect("fixed registry");
        let receipt = verify_sensitivity(&PRODUCTION_DIFFERENTIATION_PATH, &registry, 1.5, cx)
            .expect("oracles agree");
        let replay = verify_sensitivity(&PRODUCTION_DIFFERENTIATION_PATH, &registry, 1.5, cx)
            .expect("replay agrees");
        assert_eq!(receipt, replay);
        assert!(receipt.verifies_integrity());
        assert!((receipt.gradient() - receipt.dual_gradient()).abs() <= receipt.fd_tolerance());
        assert!((receipt.gradient() - receipt.fd_fine()).abs() <= receipt.fd_tolerance());
        assert!(3.0 * (receipt.fd_coarse() - receipt.fd_fine()).abs() <= receipt.fd_tolerance());
        assert_eq!(
            receipt
                .gradient_in_input_units(0.001)
                .expect("positive scale")
                .to_bits(),
            (receipt.gradient() * 0.001).to_bits()
        );
        for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                receipt.gradient_in_input_units(scale),
                Err(DifferentiationError::InvalidInputScale { bits })
                    if bits == scale.to_bits()
            ));
        }
        assert!(matches!(
            receipt.gradient_in_input_units(f64::MAX),
            Err(DifferentiationError::NonFiniteRescaledGradient { .. })
        ));

        assert!(matches!(
            verify_sensitivity(&["sdf"], &registry, 1.5, cx),
            Err(DifferentiationError::OraclePathMismatch { .. })
        ));
    });
}

#[test]
fn hostile_scalars_fail_typed_and_missing_vjp_has_precedence() {
    active_cx(|cx| {
        let registry = production_vjp_registry().expect("fixed registry");
        for input in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                differentiate_path(&PRODUCTION_DIFFERENTIATION_PATH, &registry, input, cx),
                Err(DifferentiationError::NonFiniteInput { bits })
                    if bits == input.to_bits()
            ));
        }

        let overflow = differentiate_path(
            &PRODUCTION_DIFFERENTIATION_PATH,
            &registry,
            f64::MAX / 4.0,
            cx,
        );
        assert!(matches!(
            overflow,
            Err(DifferentiationError::NonFinitePrimal { ref op, bits })
                if op == "spline" && bits == f64::INFINITY.to_bits()
        ));

        let mut poisoned = production_vjp_registry().expect("fixed registry");
        poisoned
            .register("spline", InfiniteSquareVjp)
            .expect("bounded replacement");
        assert!(matches!(
            differentiate_path(
                &PRODUCTION_DIFFERENTIATION_PATH,
                &poisoned,
                1.0,
                cx,
            ),
            Err(DifferentiationError::NonFiniteGradient { bits })
                if bits == f64::INFINITY.to_bits()
        ));

        assert!(matches!(
            differentiate_path(
                &["sdf", "remesh", "solve"],
                &registry,
                f64::NAN,
                cx,
            ),
            Err(DifferentiationError::MissingVjp { ref op }) if op == "remesh"
        ));
    });
}

#[test]
fn the_as_built_loop_localizes_a_defect_and_reduces_misfit() {
    let stage = active_cx(stage_as_built_loop).expect("as-built stage succeeds");
    assert_eq!(stage.status, StageStatus::Passed);
    assert_eq!(stage.evidence_identity, AS_BUILT_EVIDENCE_IDENTITY);
    assert!(stage.events.iter().any(|event| matches!(
        event,
        StageEvent::AsBuiltDelta {
            defect_index: Some(1),
            estimated: true,
            ..
        }
    )));
    assert!(
        stage
            .events
            .iter()
            .any(|event| matches!(event, StageEvent::Assimilation { reduced: true, .. }))
    );
}

#[test]
fn tolerance_uses_sealed_sensitivities_without_a_probability_claim() {
    let stage = active_cx(stage_tolerance_allocation).expect("tolerance stage evaluates");
    assert_eq!(stage.status, StageStatus::Passed);
    assert_eq!(stage.evidence_identity, TOLERANCE_EVIDENCE_IDENTITY);
    assert_eq!(
        stage
            .events
            .iter()
            .filter(|event| matches!(event, StageEvent::GradientVerified { .. }))
            .count(),
        2
    );
    assert!(stage.events.iter().any(|event| matches!(
        event,
        StageEvent::ToleranceActions {
            critical: Some(Action::Tighten),
            slack: Some(Action::Loosen)
        }
    )));
    assert!(stage.events.iter().any(|event| matches!(
        event,
        StageEvent::GdtJustification {
            loosened: 1,
            all_verified: true
        }
    )));
    assert!(stage.events.iter().any(|event| matches!(
        event,
        StageEvent::SampledLinearization {
            samples: 3,
            confirmed: true,
            probability_claimed: false,
            ..
        }
    )));
}

#[test]
fn an_adverse_sample_fails_without_becoming_probability_evidence() {
    let stage = active_cx(|cx| stage_tolerance_allocation_with_samples(cx, &[100.0]))
        .expect("adverse samples are evaluated");
    assert!(matches!(stage.status, StageStatus::Failed(_)));
    assert!(stage.events.iter().any(|event| matches!(
        event,
        StageEvent::SampledLinearization {
            samples: 1,
            confirmed: false,
            probability_claimed: false,
            ..
        }
    )));
}

#[test]
fn the_spacetime_stage_is_honestly_gated() {
    let stage = active_cx(stage_spacetime_gated).expect("gate recording admitted");
    let reason = match &stage.status {
        StageStatus::Gated(reason) => reason,
        other => panic!("spacetime must be gated, got {other:?}"),
    };
    assert_eq!(reason.code, "diffreal.spacetime.integration-not-activated");
    assert_eq!(stage.evidence_identity, SPACETIME_EVIDENCE_IDENTITY);
    assert!(stage.events.iter().any(|event| matches!(
        event,
        StageEvent::Gate { code, detail }
            if *code == "diffreal.spacetime.integration-not-activated"
                && detail.contains("bk0o.7 is shipped")
    )));
}

#[test]
fn status_and_typed_event_displays_are_stable() {
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

    let log = StageLog::new(
        "display-fixture",
        StageRequirement::Optional,
        gated,
        "display-fixture/v1",
        vec![StageEvent::Gate {
            code: "test.gated",
            detail: "gate recorded".to_string(),
        }],
    );
    assert_eq!(
        log.to_string(),
        "stage=display-fixture requirement=optional status=gated[test.gated]: the capability was unavailable evidence_identity=display-fixture/v1"
    );

    let event = StageEvent::SampledLinearization {
        samples: 3,
        confirmed: true,
        linearized_std_bits: 0.25_f64.to_bits(),
        probability_claimed: false,
    };
    assert_eq!(
        event.to_string(),
        "event=sampled-linearization samples=3 confirmed=true linearized_std=0x3fd0000000000000 probability_claimed=false"
    );
}

#[test]
fn the_battery_is_deterministic() {
    let first = active_cx(run_battery).expect("first battery succeeds");
    let second = active_cx(run_battery).expect("replay succeeds");
    assert_eq!(first, second);
}

fn assert_cancelled(stage: &'static str, result: Result<StageLog, DiffRealError>) {
    assert!(matches!(
        result,
        Err(DiffRealError::Cancelled { stage: observed }) if observed == stage
    ));
}

fn assert_zero_cost_refused(stage: &'static str, result: Result<StageLog, DiffRealError>) {
    assert!(matches!(
        result,
        Err(DiffRealError::WorkBudgetExceeded {
            stage: observed,
            required,
            available: 0
        }) if observed == stage && required > 0
    ));
}

#[test]
fn every_stage_polls_cancellation_and_admits_its_fixed_work_budget() {
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_cancelled(
        DIFFERENTIATION_STAGE,
        with_cx(&cancelled, Budget::INFINITE, stage_differentiation),
    );
    assert_cancelled(
        AS_BUILT_STAGE,
        with_cx(&cancelled, Budget::INFINITE, stage_as_built_loop),
    );
    assert_cancelled(
        TOLERANCE_STAGE,
        with_cx(&cancelled, Budget::INFINITE, stage_tolerance_allocation),
    );
    assert_cancelled(
        SPACETIME_STAGE,
        with_cx(&cancelled, Budget::INFINITE, stage_spacetime_gated),
    );

    let zero = Budget::INFINITE.with_cost_quota(0);
    let active = CancelGate::new();
    assert_zero_cost_refused(
        DIFFERENTIATION_STAGE,
        with_cx(&active, zero, stage_differentiation),
    );
    assert_zero_cost_refused(AS_BUILT_STAGE, with_cx(&active, zero, stage_as_built_loop));
    assert_zero_cost_refused(
        TOLERANCE_STAGE,
        with_cx(&active, zero, stage_tolerance_allocation),
    );
    assert_zero_cost_refused(
        SPACETIME_STAGE,
        with_cx(&active, zero, stage_spacetime_gated),
    );

    let mid_stage_gate = Arc::new(CancelGate::new());
    let mid_stage_result = with_cx(&mid_stage_gate, Budget::INFINITE, |cx| {
        let mut registry = production_vjp_registry().expect("fixed registry");
        registry
            .register(
                "spline",
                CancellingSquareVjp {
                    gate: Arc::clone(&mid_stage_gate),
                },
            )
            .expect("bounded replacement");
        stage_differentiation_with_registry(cx, &registry)
    });
    assert_cancelled(DIFFERENTIATION_STAGE, mid_stage_result);
}

#[test]
fn cancellation_and_budget_refusal_publish_no_partial_battery() {
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        with_cx(&cancelled, Budget::INFINITE, run_battery),
        Err(DiffRealError::Cancelled {
            stage: DIFFERENTIATION_STAGE
        })
    ));
    assert!(matches!(
        with_cx(
            &CancelGate::new(),
            Budget::INFINITE.with_cost_quota(0),
            run_battery,
        ),
        Err(DiffRealError::WorkBudgetExceeded {
            stage: DIFFERENTIATION_STAGE,
            available: 0,
            ..
        })
    ));
}
