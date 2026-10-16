//! Layer-3 differentiation-and-reality conformance. These tests falsify the
//! shared tape/VJP path independently, exercise typed hostile-input refusals,
//! and keep sampled linearization evidence separate from probability claims.

use fs_adjoint::transpose::Vjp;
use fs_diffreal_e2e::{
    AS_BUILT_EVIDENCE_IDENTITY, AS_BUILT_STAGE, DIFFERENTIATION_EVIDENCE_IDENTITY,
    DIFFERENTIATION_STAGE, DiffRealError, DiffRealReport, DifferentiationError,
    NoPromotionReceiptVerifier, PRODUCTION_DIFFERENTIATION_PATH, PromotionAttestation,
    PromotionReceiptError, PromotionReceiptVerifier, PromotionVerdict,
    PromotionVerificationDecision, PromotionVerificationRequest, SPACETIME_EVIDENCE_IDENTITY,
    SPACETIME_STAGE, StageEvent, StageLog, StageReason, StageRequirement, StageStatus,
    TOLERANCE_EVIDENCE_IDENTITY, TOLERANCE_STAGE, differentiate_path, production_vjp_registry,
    promotion_policy_fingerprint, run_battery, run_battery_with_clock, stage_as_built_loop,
    stage_as_built_loop_with_clock, stage_differentiation, stage_differentiation_with_registry,
    stage_spacetime_gated, stage_tolerance_allocation, stage_tolerance_allocation_with_samples,
    verify_sensitivity,
};
use fs_exec::{
    Budget, CancelGate, Cx, ExecMode, InvocationDisposition, InvocationError, StreamKey, Time,
    VirtualClock,
};
use fs_toleralloc::Action;
use std::sync::Arc;

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_identity_cx(
        gate,
        budget,
        StreamKey {
            seed: 0x6469_6666_7265_616c,
            kernel_id: 1,
            tile: 0,
            iteration: 0,
        },
        ExecMode::Deterministic,
        f,
    )
}

fn with_identity_cx<R>(
    gate: &CancelGate,
    budget: Budget,
    stream_key: StreamKey,
    mode: ExecMode,
    f: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(gate, arena, stream_key, budget, mode);
        f(&cx)
    })
}

fn active_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx(&CancelGate::new(), Budget::INFINITE, f)
}

fn clocked_battery(budget: Budget) -> Result<DiffRealReport, DiffRealError> {
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    with_cx(&gate, budget, |cx| run_battery_with_clock(cx, &clock))
}

#[test]
fn the_full_layer3_battery_reports_its_required_gate_fail_closed() {
    let report = active_cx(run_battery).expect("battery succeeds");
    assert!(!report.complete(), "a required gated stage is incomplete");
    assert!(!report.all_required_passed());
    assert!(!report.structurally_ready());
    assert!(report.verifies_integrity());
    assert_eq!(report.receipts().len(), 4);
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

#[derive(Debug)]
struct FixturePromotionVerifier;

impl PromotionReceiptVerifier for FixturePromotionVerifier {
    fn verify(
        &self,
        request: &PromotionVerificationRequest<'_>,
        attestation: &PromotionAttestation,
    ) -> PromotionVerificationDecision {
        let verdict = match attestation.key_id() {
            "trusted-diffreal-fixture"
                if attestation.signature() == request.subject().as_bytes() =>
            {
                PromotionVerdict::Authorized
            }
            "trusted-diffreal-fixture" => PromotionVerdict::WrongSignature,
            "revoked-diffreal-fixture" => PromotionVerdict::RevokedKey,
            _ => PromotionVerdict::UnknownKey,
        };
        PromotionVerificationDecision::new(verdict, request.policy_fingerprint())
    }
}

#[derive(Debug)]
struct WrongPolicyDecisionVerifier;

impl PromotionReceiptVerifier for WrongPolicyDecisionVerifier {
    fn verify(
        &self,
        request: &PromotionVerificationRequest<'_>,
        _attestation: &PromotionAttestation,
    ) -> PromotionVerificationDecision {
        let mut wrong = request.policy_fingerprint();
        wrong.0[0] ^= 1;
        PromotionVerificationDecision::new(PromotionVerdict::Authorized, wrong)
    }
}

fn attestation_for(key_id: &str, signature: Vec<u8>) -> PromotionAttestation {
    PromotionAttestation::new(key_id, signature, promotion_policy_fingerprint())
}

#[test]
fn report_authentication_is_external_fail_closed_and_not_a_scientific_pass() {
    active_cx(|cx| {
        let report = run_battery(cx).expect("battery succeeds");
        let signature = report.promotion_verification_subject().as_bytes().to_vec();
        let authenticated = report
            .authenticate(
                cx,
                attestation_for("trusted-diffreal-fixture", signature.clone()),
                &FixturePromotionVerifier,
            )
            .expect("fixture authority authenticates the exact ordered root");
        assert_eq!(authenticated.receipt_root(), report.receipt_root());
        assert!(!authenticated.promotion_ready(), "spacetime is still gated");

        let denied = report.authenticate(
            cx,
            attestation_for("trusted-diffreal-fixture", signature.clone()),
            &NoPromotionReceiptVerifier,
        );
        assert!(matches!(
            denied,
            Err(PromotionReceiptError::Unauthorized {
                verdict: PromotionVerdict::UnknownKey
            })
        ));

        for (key_id, signature, expected) in [
            (
                "trusted-diffreal-fixture",
                b"wrong-signature".to_vec(),
                PromotionVerdict::WrongSignature,
            ),
            (
                "unknown-diffreal-fixture",
                signature.clone(),
                PromotionVerdict::UnknownKey,
            ),
            (
                "revoked-diffreal-fixture",
                signature.clone(),
                PromotionVerdict::RevokedKey,
            ),
        ] {
            assert!(matches!(
                report.authenticate(
                    cx,
                    attestation_for(key_id, signature),
                    &FixturePromotionVerifier,
                ),
                Err(PromotionReceiptError::Unauthorized { verdict }) if verdict == expected
            ));
        }

        let mut wrong_policy = promotion_policy_fingerprint();
        wrong_policy.0[0] ^= 1;
        assert!(matches!(
            report.authenticate(
                cx,
                PromotionAttestation::new(
                    "trusted-diffreal-fixture",
                    signature.clone(),
                    wrong_policy,
                ),
                &FixturePromotionVerifier,
            ),
            Err(PromotionReceiptError::AttestationPolicyMismatch)
        ));
        assert!(matches!(
            report.authenticate(
                cx,
                attestation_for("trusted-diffreal-fixture", signature),
                &WrongPolicyDecisionVerifier,
            ),
            Err(PromotionReceiptError::DecisionPolicyMismatch)
        ));
    });
}

#[test]
fn every_cx_identity_field_is_replay_bound() {
    let base_stream = StreamKey {
        seed: 0x6469_6666_7265_616c,
        kernel_id: 1,
        tile: 0,
        iteration: 0,
    };
    let base_budget = Budget::with_deadline_at_ns(100)
        .with_poll_quota(1_000)
        .with_cost_quota(20_000)
        .with_priority(7);
    let clock = VirtualClock::new();
    let report = with_identity_cx(
        &CancelGate::new(),
        base_budget,
        base_stream,
        ExecMode::Deterministic,
        |cx| run_battery_with_clock(cx, &clock),
    )
    .expect("battery succeeds");
    let execution = report.execution_identity();
    let attestation = attestation_for(
        "trusted-diffreal-fixture",
        report.promotion_verification_subject().as_bytes().to_vec(),
    );

    let stream_variants = [
        StreamKey {
            seed: base_stream.seed ^ 1,
            ..base_stream
        },
        StreamKey {
            kernel_id: base_stream.kernel_id ^ 1,
            ..base_stream
        },
        StreamKey {
            tile: base_stream.tile ^ 1,
            ..base_stream
        },
        StreamKey {
            iteration: base_stream.iteration ^ 1,
            ..base_stream
        },
    ];
    for stream in stream_variants {
        with_identity_cx(
            &CancelGate::new(),
            base_budget,
            stream,
            execution.mode(),
            |cx| {
                assert!(matches!(
                    report.authenticate(cx, attestation.clone(), &FixturePromotionVerifier),
                    Err(PromotionReceiptError::ReplayContextMismatch)
                ));
            },
        );
    }

    let budget_variants = [
        Budget {
            deadline: Budget::with_deadline_at_ns(101).deadline,
            ..base_budget
        },
        Budget {
            deadline: None,
            ..base_budget
        },
        base_budget.with_poll_quota(base_budget.poll_quota + 1),
        base_budget.with_cost_quota(20_001),
        Budget {
            cost_quota: None,
            ..base_budget
        },
        base_budget.with_priority(base_budget.priority.wrapping_add(1)),
    ];
    for budget in budget_variants {
        with_identity_cx(
            &CancelGate::new(),
            budget,
            base_stream,
            execution.mode(),
            |cx| {
                assert!(matches!(
                    report.authenticate(cx, attestation.clone(), &FixturePromotionVerifier),
                    Err(PromotionReceiptError::ReplayContextMismatch)
                ));
            },
        );
    }

    with_identity_cx(
        &CancelGate::new(),
        base_budget,
        base_stream,
        ExecMode::Fast,
        |cx| {
            assert!(matches!(
                report.authenticate(cx, attestation, &FixturePromotionVerifier),
                Err(PromotionReceiptError::ReplayContextMismatch)
            ));
        },
    );
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
fn g0_as_built_receipt_is_one_completed_affine_transaction() {
    let discovery = clocked_battery(Budget::INFINITE).expect("discover fixed invocation plan");
    let required = discovery
        .receipt(AS_BUILT_STAGE)
        .and_then(|stage| stage.invocation())
        .expect("as-built receipt carries invocation evidence")
        .required();
    let exact_budget = Budget::INFINITE
        .with_poll_quota(required.polls().get())
        .with_cost_quota(required.cost().get());
    let report = clocked_battery(exact_budget).expect("exact typed envelope admits the battery");
    let stage_receipt = report
        .receipt(AS_BUILT_STAGE)
        .expect("as-built stage receipt exists");
    let invocation = stage_receipt
        .invocation()
        .expect("as-built stage retains the full invocation receipt");

    assert_eq!(invocation.disposition(), InvocationDisposition::Completed);
    assert!(invocation.failure().is_none());
    assert!(invocation.verifies_integrity());
    assert!(report.verifies_integrity());
    assert!(
        report.receipts().iter().all(|receipt| {
            (receipt.stage() == AS_BUILT_STAGE) == receipt.invocation().is_some()
        })
    );
    assert_eq!(invocation.limits().resources(), invocation.required());
    assert_eq!(
        invocation.memory_requested_bytes(),
        invocation.memory_released_bytes(),
        "every temporary scientific allocation is released before sealing"
    );
    assert!(invocation.memory_peak_bytes() <= required.memory().get());
    assert!(invocation.output_retained_bytes() <= required.output().get());

    let children = invocation.children();
    let phases: Vec<_> = children.iter().map(|child| child.phase()).collect();
    assert_eq!(
        phases.as_slice(),
        &[
            "as-built.transaction",
            "as-built.registration",
            "as-built.comparison",
            "as-built.prior",
            "as-built.assimilation",
            "as-built.publication",
        ]
    );
    assert_eq!(children[0].parent(), None);
    assert_eq!(children[0].granted(), required);
    assert!(
        children[1..]
            .iter()
            .all(|child| child.parent() == Some(children[0].id()))
    );
    assert!(
        children
            .iter()
            .all(|child| child.disposition() == InvocationDisposition::Completed)
    );
    assert_eq!(
        children[4].consumed().polls().get(),
        46,
        "the fixed assimilation consumes its exact mixed-stride poll count, not its cap"
    );
    assert_eq!(
        children[0].consumed().polls().get(),
        69,
        "the transaction receipt retains the exact cumulative poll spend"
    );
    assert!(children[1..].iter().all(|child| {
        child.direct_memory_peak_bytes() == child.memory_peak_bytes()
            && child.memory_requested_bytes() == child.memory_released_bytes()
            && child.output_retained_bytes() == child.direct_consumed().output().get()
    }));
    assert!(
        phases.iter().all(|phase| !phase.contains("misfit")),
        "misfit evidence comes from the single assimilation result, not reissued standalone work"
    );
}

#[test]
fn g3_as_built_fast_mode_derives_its_own_exact_work_plan() {
    let deterministic = clocked_battery(Budget::INFINITE).expect("deterministic battery");
    let fast_clock = VirtualClock::new();
    let fast = with_identity_cx(
        &CancelGate::new(),
        Budget::INFINITE,
        StreamKey {
            seed: 0x6469_6666_7265_616c,
            kernel_id: 1,
            tile: 0,
            iteration: 0,
        },
        ExecMode::Fast,
        |cx| run_battery_with_clock(cx, &fast_clock),
    )
    .expect("fast-mode battery has a mode-aware preflight");
    let deterministic_invocation = deterministic
        .receipt(AS_BUILT_STAGE)
        .and_then(|receipt| receipt.invocation())
        .expect("deterministic invocation receipt");
    let fast_invocation = fast
        .receipt(AS_BUILT_STAGE)
        .and_then(|receipt| receipt.invocation())
        .expect("fast invocation receipt");

    assert_eq!(
        deterministic_invocation.required().work().get() - fast_invocation.required().work().get(),
        9,
        "candidate hashing accounts for the nine-byte mode-name difference"
    );
    assert_eq!(fast_invocation.children()[4].consumed().polls().get(), 46);
    assert!(fast.verifies_integrity());
}

#[test]
fn g3_as_built_exact_envelope_succeeds_and_one_below_refuses() {
    let discovery = clocked_battery(Budget::INFINITE).expect("discover fixed invocation plan");
    let required = discovery
        .receipt(AS_BUILT_STAGE)
        .and_then(|stage| stage.invocation())
        .expect("as-built invocation receipt")
        .required();
    assert!(required.polls().get() > 0);
    assert!(required.cost().get() > 0);

    let exact = Budget::INFINITE
        .with_poll_quota(required.polls().get())
        .with_cost_quota(required.cost().get());
    let exact_report = clocked_battery(exact).expect("exact typed resources succeed");
    let exact_stage = exact_report
        .receipt(AS_BUILT_STAGE)
        .expect("exact run has as-built receipt");
    let exact_invocation = exact_stage.invocation().expect("invocation is retained");
    assert_eq!(exact_invocation.limits().resources(), required);

    let poll_available = required.polls().get() - 1;
    let poll_clock = VirtualClock::new();
    let poll_refusal = with_cx(
        &CancelGate::new(),
        Budget::INFINITE
            .with_poll_quota(poll_available)
            .with_cost_quota(required.cost().get()),
        |cx| stage_as_built_loop_with_clock(cx, &poll_clock),
    );
    assert!(matches!(
        poll_refusal,
        Err(DiffRealError::Invocation(InvocationError::ResourceExceeded {
            resource: "polls",
            requested,
            available,
        })) if requested == u128::from(required.polls().get())
            && available == u128::from(poll_available)
    ));

    let cost_available = required.cost().get() - 1;
    let cost_clock = VirtualClock::new();
    let cost_refusal = with_cx(
        &CancelGate::new(),
        Budget::INFINITE
            .with_poll_quota(required.polls().get())
            .with_cost_quota(cost_available),
        |cx| stage_as_built_loop_with_clock(cx, &cost_clock),
    );
    assert!(matches!(
        cost_refusal,
        Err(DiffRealError::Invocation(InvocationError::ResourceExceeded {
            resource: "cost",
            requested,
            available,
        })) if requested == u128::from(required.cost().get())
            && available == u128::from(cost_available)
    ));

    let roomy = Budget::INFINITE
        .with_poll_quota(required.polls().get() + 1)
        .with_cost_quota(required.cost().get());
    let roomy_report = clocked_battery(roomy).expect("one extra poll remains admissible");
    let roomy_stage = roomy_report
        .receipt(AS_BUILT_STAGE)
        .expect("roomy run has as-built receipt");
    assert_eq!(
        exact_report.stage(AS_BUILT_STAGE),
        roomy_report.stage(AS_BUILT_STAGE),
        "resource-envelope mutation does not change the scientific diagnostic"
    );
    assert_ne!(
        exact_invocation.root(),
        roomy_stage
            .invocation()
            .expect("invocation retained")
            .root(),
        "the immutable envelope is bound into invocation identity"
    );
    let roomy_invocation = roomy_stage.invocation().expect("invocation retained");
    assert_eq!(
        exact_invocation.limits().accuracy_obligation(),
        roomy_invocation.limits().accuracy_obligation(),
        "ambient capacity cannot weaken the fixed accuracy obligation"
    );
    assert_eq!(
        exact_invocation.limits().capability_scope(),
        roomy_invocation.limits().capability_scope(),
        "ambient capacity cannot widen the fixed capability scope"
    );
    assert_ne!(
        exact_stage.root(),
        roomy_stage.root(),
        "the authenticated stage receipt binds the invocation receipt"
    );
    assert!(roomy_report.verifies_integrity());
}

#[test]
fn g4_as_built_cancellation_and_deadline_publish_no_stage_or_report() {
    let discovery = clocked_battery(Budget::INFINITE).expect("discover fixed invocation plan");
    let required = discovery
        .receipt(AS_BUILT_STAGE)
        .and_then(|stage| stage.invocation())
        .expect("as-built invocation receipt")
        .required();
    let exact = Budget::INFINITE
        .with_poll_quota(required.polls().get())
        .with_cost_quota(required.cost().get());

    let cancelled = CancelGate::new();
    cancelled.request();
    let cancellation_clock = VirtualClock::new();
    let cancellation = with_cx(&cancelled, exact, |cx| {
        stage_as_built_loop_with_clock(cx, &cancellation_clock)
    });
    match cancellation {
        Err(DiffRealError::InvocationDidNotComplete {
            disposition,
            receipt_root,
            receipt,
            cause: Some(cause),
        }) => {
            assert_eq!(disposition, InvocationDisposition::Cancelled);
            assert_eq!(receipt_root, receipt.root());
            assert_eq!(receipt.disposition(), InvocationDisposition::Cancelled);
            assert!(receipt.verifies_integrity());
            assert!(matches!(
                *cause,
                DiffRealError::Invocation(InvocationError::Cancelled {
                    phase: "as-built.setup.begin"
                })
            ));
        }
        other => panic!("expected drained cancellation receipt, got {other:?}"),
    }

    let deadline_ns = 17;
    let deadline_clock = VirtualClock::starting_at(Time::from_nanos(deadline_ns));
    let deadline_budget = Budget::with_deadline_at_ns(deadline_ns)
        .with_poll_quota(required.polls().get())
        .with_cost_quota(required.cost().get());
    let deadline = with_cx(&CancelGate::new(), deadline_budget, |cx| {
        run_battery_with_clock(cx, &deadline_clock)
    });
    assert!(matches!(
        deadline,
        Err(DiffRealError::Invocation(
            InvocationError::DeadlineExpired {
                phase: "invocation-admission",
                deadline_ns: 17,
                observed_ns: 17,
            }
        ))
    ));
}

#[test]
fn g5_as_built_invocation_receipt_replays_bit_for_bit() {
    let discovery = clocked_battery(Budget::INFINITE).expect("discover fixed invocation plan");
    let required = discovery
        .receipt(AS_BUILT_STAGE)
        .and_then(|stage| stage.invocation())
        .expect("as-built invocation receipt")
        .required();
    let exact = Budget::INFINITE
        .with_poll_quota(required.polls().get())
        .with_cost_quota(required.cost().get());
    let first = clocked_battery(exact).expect("first deterministic run");
    let replay = clocked_battery(exact).expect("deterministic replay");
    let first_stage = first
        .receipt(AS_BUILT_STAGE)
        .expect("first as-built receipt");
    let replay_stage = replay
        .receipt(AS_BUILT_STAGE)
        .expect("replay as-built receipt");

    assert_eq!(first_stage.invocation(), replay_stage.invocation());
    assert_eq!(first_stage.root(), replay_stage.root());
    assert_eq!(first.receipt_root(), replay.receipt_root());
    assert_eq!(
        first_stage
            .invocation()
            .expect("first invocation")
            .last_deadline_observation(),
        None,
        "a deadline-free replay does not bind nondeterministic wall time"
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

fn assert_invocation_cancelled(result: Result<StageLog, DiffRealError>) {
    match result {
        Err(DiffRealError::InvocationDidNotComplete {
            disposition,
            receipt_root,
            receipt,
            cause: Some(cause),
        }) => {
            assert_eq!(disposition, InvocationDisposition::Cancelled);
            assert_eq!(receipt_root, receipt.root());
            assert_eq!(receipt.disposition(), InvocationDisposition::Cancelled);
            assert!(receipt.verifies_integrity());
            assert!(matches!(
                *cause,
                DiffRealError::Invocation(InvocationError::Cancelled { .. })
            ));
        }
        other => panic!("expected drained invocation cancellation, got {other:?}"),
    }
}

fn assert_invocation_resource_refused(
    resource: &'static str,
    available: u128,
    result: Result<StageLog, DiffRealError>,
) {
    assert!(matches!(
        result,
        Err(DiffRealError::Invocation(InvocationError::ResourceExceeded {
            resource: observed,
            requested,
            available: observed_available,
        })) if observed == resource && requested > available && observed_available == available
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
    assert_invocation_cancelled(with_cx(&cancelled, Budget::INFINITE, stage_as_built_loop));
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
    assert_invocation_resource_refused("cost", 0, with_cx(&active, zero, stage_as_built_loop));
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
