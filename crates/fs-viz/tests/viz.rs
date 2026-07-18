//! Battery for scientific visualization (fs-viz). Each test checks a primitive
//! against ANALYTIC ground truth: rotation streamlines are circles, saddle
//! streamlines conserve xy, Hessian classification recovers the known Morse
//! type, and a circle-SDF isocontour lies on the circle.

use fs_blake3::DomainHasher;
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey};
use fs_viz::{
    CriticalKind, Grid2, Grid2Error, Grid3, Grid3Error, HessianSymmetryPolicy,
    ISO_CONTOUR_ARTIFACT_IDENTITY_DOMAIN, IsoContourDisposition, IsoContourError,
    IsoContourResource, IsoSurfaceError, SCALAR_FIELD3_ARTIFACT_KIND, SCALAR_FIELD3_SCHEMA_VERSION,
    STREAMLINE_ARTIFACT_IDENTITY_DOMAIN, STREAMLINE_ARTIFACT_IDENTITY_VERSION, ScalarField3,
    ScalarField3Error, ScalarFieldSemantics, ScalarLayout3, StreamlineBoundaryPolicy,
    StreamlineDisposition, StreamlineDomain2, StreamlineError, StreamlineRequest,
    StreamlineResource, StreamlineStage, StreamlineStagnationPolicy, StreamlineTermination,
    classify_hessian, classify_hessian_with_policy, required_streamline_budget, streamline,
    streamline_plan, streamline_with_cx,
};
use std::cell::Cell;
use std::mem::size_of;

fn radius(p: [f64; 2]) -> f64 {
    (p[0] * p[0] + p[1] * p[1]).sqrt()
}

fn lower_left_collapse_error(
    lower: f64,
    upper: f64,
    lower_left: f64,
    other: f64,
    iso: f64,
    crossing_limit: usize,
) -> IsoContourError {
    let grid = Grid2::from_fn(2, 2, [lower; 2], [upper; 2], 4, |point| {
        if point[0].to_bits() == lower.to_bits() && point[1].to_bits() == lower.to_bits() {
            lower_left
        } else {
            other
        }
    })
    .expect("adjacent finite endpoints form an admitted 2x2 grid");
    grid.isocontour_crossings(iso, crossing_limit)
        .expect_err("the strict real crossing must not collapse to a binary64 endpoint")
}

fn with_contour_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xF5_71_2,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

#[test]
fn a_rotation_field_streams_along_a_circle() {
    // u = (-y, x): rigid rotation, so the radius is conserved.
    let line = streamline(|p| [-p[1], p[0]], [1.0, 0.0], 0.01, 400);
    for p in &line {
        assert!(
            (radius(*p) - 1.0).abs() < 1e-3,
            "radius {} drifted",
            radius(*p)
        );
    }
    // it actually goes somewhere (not a fixed point).
    assert!((line.last().unwrap()[1]).abs() > 0.1);
}

#[test]
fn a_saddle_field_conserves_the_hyperbola_invariant() {
    // u = (x, -y): flow x·y is invariant along a streamline.
    let line = streamline(|p| [p[0], -p[1]], [1.0, 1.0], 0.01, 50);
    for p in &line {
        assert!(
            (p[0] * p[1] - 1.0).abs() < 1e-4,
            "xy = {} drifted",
            p[0] * p[1]
        );
    }
    // x grows, y shrinks (the saddle's unstable/stable manifolds).
    assert!(line.last().unwrap()[0] > 1.4 && line.last().unwrap()[1] < 0.7);
}

#[test]
fn g0_streamline_zero_step_plan_and_receipt_are_exact() {
    let request = StreamlineRequest::dimensionless_rk4([1.0, -2.0], 0.25, 0, 3);
    let budget = required_streamline_budget(request).expect("zero-step plan is representable");
    let plan = streamline_plan(request, budget).expect("exact budget admits zero-step plan");
    assert_eq!(plan.steps, 0);
    assert_eq!(plan.field_evaluations, 0);
    assert_eq!(plan.output_points, 1);
    assert_eq!(plan.output_bytes, size_of::<[f64; 2]>());
    assert_eq!(plan.polls, 3);
    let one_request = StreamlineRequest::dimensionless_rk4([1.0, -2.0], 0.25, 1, 3);
    let one_budget = required_streamline_budget(one_request).expect("one-step exact budget");
    let one_plan = streamline_plan(one_request, one_budget).expect("one-step exact plan");
    assert_eq!(one_plan.field_evaluations, 4);
    assert_eq!(one_plan.output_points, 2);

    let gate = CancelGate::new();
    let callback_calls = Cell::new(0usize);
    let output = with_contour_cx(&gate, Budget::INFINITE, |cx| {
        streamline_with_cx(
            cx,
            |_| {
                callback_calls.set(callback_calls.get() + 1);
                [0.0; 2]
            },
            request,
            budget,
        )
        .expect("zero-step request publishes its seed")
    });
    assert_eq!(callback_calls.get(), 0);
    assert_eq!(output.points(), &[[1.0, -2.0]]);
    assert_eq!(output.report().request, request);
    assert_eq!(output.report().operation_budget, Some(budget));
    assert_eq!(output.report().plan, Some(plan));
    assert_eq!(
        output.report().termination,
        Some(StreamlineTermination::StepsComplete)
    );
    assert_eq!(
        output.report().disposition,
        StreamlineDisposition::Completed
    );
    assert_eq!(output.report().field_evaluations, 0);
    assert_eq!(output.report().output_points, 1);
    assert_eq!(output.report().polls, plan.polls);
    assert_eq!(output.report().work_units, plan.work_units);
    assert_eq!(output.report().reserved_output_bytes, plan.output_bytes);
    assert_eq!(output.report().peak_live_bytes, plan.live_bytes);
    assert_eq!(output.report().identity_bytes_hashed, plan.identity_bytes);
    assert_eq!(output.report().error_estimate, None);
    assert!(output.report().terminal && output.report().published);
    assert!(output.report().artifact_identity.is_some());
    assert!(!STREAMLINE_ARTIFACT_IDENTITY_DOMAIN.is_empty());
}

#[test]
#[allow(clippy::too_many_lines)] // One admission matrix shares a callback-side-effect oracle.
fn g0_streamline_invalid_requests_and_one_short_resources_fail_before_callbacks() {
    assert!(matches!(
        required_streamline_budget(StreamlineRequest::dimensionless_rk4(
            [f64::NAN, 0.0],
            0.1,
            1,
            1,
        )),
        Err(StreamlineError::NonFiniteSeed { component: 0, .. })
    ));
    assert_eq!(
        required_streamline_budget(StreamlineRequest::dimensionless_rk4([0.0; 2], 0.0, 1, 1,)),
        Err(StreamlineError::InvalidStepSize { dt: 0.0 })
    );
    assert_eq!(
        required_streamline_budget(StreamlineRequest::dimensionless_rk4(
            [0.0; 2],
            f64::INFINITY,
            1,
            1,
        )),
        Err(StreamlineError::InvalidStepSize { dt: f64::INFINITY })
    );
    let mut unsupported_version = StreamlineRequest::dimensionless_rk4([0.0; 2], 0.1, 1, 1);
    unsupported_version.method_version += 1;
    assert_eq!(
        required_streamline_budget(unsupported_version),
        Err(StreamlineError::UnsupportedMethodVersion { version: 2 })
    );
    assert_eq!(
        required_streamline_budget(StreamlineRequest::dimensionless_rk4([0.0; 2], 0.1, 1, 0,)),
        Err(StreamlineError::InvalidPollStride { items_per_poll: 0 })
    );
    let mut invalid_domain = StreamlineRequest::dimensionless_rk4([0.0; 2], 0.1, 1, 1);
    invalid_domain.domain = Some(StreamlineDomain2 {
        lower: [-f64::MAX, -1.0],
        upper: [f64::MAX, 1.0],
    });
    assert!(matches!(
        required_streamline_budget(invalid_domain),
        Err(StreamlineError::InvalidDomain { axis: 0, .. })
    ));
    let mut outside_domain = StreamlineRequest::dimensionless_rk4([2.0, 0.0], 0.1, 1, 1);
    outside_domain.domain = Some(StreamlineDomain2 {
        lower: [-1.0; 2],
        upper: [1.0; 2],
    });
    assert!(matches!(
        required_streamline_budget(outside_domain),
        Err(StreamlineError::SeedOutsideDomain { axis: 0, .. })
    ));
    assert_eq!(
        required_streamline_budget(StreamlineRequest::dimensionless_rk4(
            [0.0; 2],
            0.1,
            usize::MAX,
            1,
        )),
        Err(StreamlineError::PlanOverflow {
            resource: StreamlineResource::OutputPoints
        })
    );
    let wrapper_calls = Cell::new(0usize);
    assert!(
        streamline(
            |_| {
                wrapper_calls.set(wrapper_calls.get() + 1);
                [1.0; 2]
            },
            [0.0; 2],
            0.1,
            usize::MAX,
        )
        .is_empty()
    );
    assert_eq!(wrapper_calls.get(), 0);
    assert!(
        streamline(|_| panic!("compatibility callback panic"), [0.0; 2], 0.1, 1,).is_empty(),
        "the no-authority wrapper contains callback unwinds"
    );

    let request = StreamlineRequest::dimensionless_rk4([0.0; 2], 0.125, 8, 2);
    let budget = required_streamline_budget(request).expect("finite exact request");
    let mut cases = Vec::new();
    macro_rules! one_short {
        ($resource:expr, $field:ident) => {{
            let mut limited = budget;
            limited.$field -= 1;
            cases.push(($resource, limited));
        }};
    }
    one_short!(StreamlineResource::Steps, step_limit);
    one_short!(StreamlineResource::FieldEvaluations, field_evaluation_limit);
    one_short!(StreamlineResource::OutputPoints, output_point_limit);
    one_short!(StreamlineResource::OutputBytes, output_byte_limit);
    one_short!(StreamlineResource::ScratchBytes, scratch_byte_limit);
    one_short!(
        StreamlineResource::DiagnosticRecords,
        diagnostic_record_limit
    );
    one_short!(StreamlineResource::DiagnosticBytes, diagnostic_byte_limit);
    one_short!(StreamlineResource::LiveBytes, live_byte_limit);
    one_short!(StreamlineResource::IdentityBytes, identity_byte_limit);
    one_short!(StreamlineResource::Polls, poll_limit);
    one_short!(StreamlineResource::WorkUnits, work_unit_limit);

    for (resource, limited) in cases {
        let calls = Cell::new(0usize);
        let gate = CancelGate::new();
        let refusal = with_contour_cx(&gate, Budget::INFINITE, |cx| {
            streamline_with_cx(
                cx,
                |_| {
                    calls.set(calls.get() + 1);
                    [1.0, 0.0]
                },
                request,
                limited,
            )
            .expect_err("one-short operation resource must refuse")
        });
        assert!(matches!(
            refusal.error,
            StreamlineError::OperationBudgetExceeded {
                resource: rejected,
                ..
            } if rejected == resource
        ));
        assert_eq!(calls.get(), 0, "{resource} admitted a callback");
        assert_eq!(refusal.report.polls, 0);
        assert_eq!(refusal.report.output_points, 0);
        assert!(!refusal.report.published);
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One G3 fixture shares exact constant-flow semantics.
fn g3_streamline_reverse_time_scaling_domain_and_stagnation_policies() {
    let forward_request = StreamlineRequest::dimensionless_rk4([0.0; 2], 0.25, 4, 2);
    let forward_budget =
        required_streamline_budget(forward_request).expect("forward constant-flow plan");
    let forward_gate = CancelGate::new();
    let forward = with_contour_cx(&forward_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, -2.0], forward_request, forward_budget)
            .expect("forward constant flow")
    });
    assert_eq!(forward.points().last().copied(), Some([1.0, -2.0]));

    let reverse_request = StreamlineRequest::dimensionless_rk4([1.0, -2.0], -0.25, 4, 2);
    let reverse_budget =
        required_streamline_budget(reverse_request).expect("reverse constant-flow plan");
    let reverse_gate = CancelGate::new();
    let reverse = with_contour_cx(&reverse_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, -2.0], reverse_request, reverse_budget)
            .expect("negative dt is explicitly supported")
    });
    assert_eq!(reverse.points().last().copied(), Some([0.0, 0.0]));

    let scaled_request = StreamlineRequest::dimensionless_rk4([0.0; 2], 0.25, 4, 2);
    let scaled_budget =
        required_streamline_budget(scaled_request).expect("scaled constant-flow plan");
    let scaled_gate = CancelGate::new();
    let scaled = with_contour_cx(&scaled_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [2.0, -4.0], scaled_request, scaled_budget)
            .expect("scaled constant flow")
    });
    for (base, doubled) in forward.points().iter().zip(scaled.points()) {
        assert_eq!([2.0 * base[0], 2.0 * base[1]], *doubled);
    }

    let mut domain_request = StreamlineRequest::dimensionless_rk4([0.0; 2], 1.0, 4, 1);
    domain_request.domain = Some(StreamlineDomain2 {
        lower: [-1.0; 2],
        upper: [1.5, 1.0],
    });
    domain_request.boundary_policy = StreamlineBoundaryPolicy::StopBeforeExit;
    let domain_budget = required_streamline_budget(domain_request).expect("bounded-domain plan");
    let domain_gate = CancelGate::new();
    let stopped = with_contour_cx(&domain_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, 0.0], domain_request, domain_budget)
            .expect("declared boundary stop publishes the in-domain prefix")
    });
    assert_eq!(stopped.points(), &[[0.0, 0.0], [1.0, 0.0]]);
    assert_eq!(
        stopped.report().termination,
        Some(StreamlineTermination::DomainExit { attempted_step: 1 })
    );
    assert_eq!(
        stopped.report().disposition,
        StreamlineDisposition::Terminated
    );

    let mut refusing_request = domain_request;
    refusing_request.boundary_policy = StreamlineBoundaryPolicy::RefuseExit;
    let refusing_budget =
        required_streamline_budget(refusing_request).expect("refusing-domain plan");
    let refusing_gate = CancelGate::new();
    let refusal = with_contour_cx(&refusing_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, 0.0], refusing_request, refusing_budget)
            .expect_err("RefuseExit cannot publish a prefix")
    });
    assert!(matches!(
        refusal.error,
        StreamlineError::DomainExit { step: 1, .. }
    ));
    assert!(!refusal.report.published);

    let mut stagnating_request = StreamlineRequest::dimensionless_rk4([3.0, -4.0], 1.0, 8, 2);
    stagnating_request.stagnation_policy = StreamlineStagnationPolicy::StopBeforeRepeat;
    let stagnating_budget =
        required_streamline_budget(stagnating_request).expect("stagnation plan");
    let stagnating_gate = CancelGate::new();
    let stagnated = with_contour_cx(&stagnating_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [0.0; 2], stagnating_request, stagnating_budget)
            .expect("declared stagnation publishes one unique point")
    });
    assert_eq!(stagnated.points(), &[[3.0, -4.0]]);
    assert_eq!(
        stagnated.report().termination,
        Some(StreamlineTermination::Stagnated { attempted_step: 0 })
    );
    assert_eq!(streamline(|_| [0.0; 2], [3.0, -4.0], 1.0, 2).len(), 3);
}

#[test]
fn g3_streamline_rk4_refinement_and_poll_chunking_preserve_claimed_semantics() {
    let exact = [1.0_f64.cos(), 1.0_f64.sin()];
    let coarse = streamline(|point| [-point[1], point[0]], [1.0, 0.0], 0.2, 5);
    let fine = streamline(|point| [-point[1], point[0]], [1.0, 0.0], 0.1, 10);
    let endpoint_error = |points: &[[f64; 2]]| {
        let endpoint = points[points.len() - 1];
        ((endpoint[0] - exact[0]).powi(2) + (endpoint[1] - exact[1]).powi(2)).sqrt()
    };
    assert!(
        endpoint_error(&fine) < endpoint_error(&coarse) / 8.0,
        "halving h should exhibit the declared fourth-order RK4 trend"
    );

    let request_one = StreamlineRequest::dimensionless_rk4([1.0, 0.0], 0.1, 10, 1);
    let request_four = StreamlineRequest::dimensionless_rk4([1.0, 0.0], 0.1, 10, 4);
    let budget_one = required_streamline_budget(request_one).expect("stride-one plan");
    let budget_four = required_streamline_budget(request_four).expect("stride-four plan");
    let gate_one = CancelGate::new();
    let one = with_contour_cx(&gate_one, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |point| [-point[1], point[0]], request_one, budget_one)
            .expect("stride-one run")
    });
    let gate_four = CancelGate::new();
    let four = with_contour_cx(&gate_four, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |point| [-point[1], point[0]], request_four, budget_four)
            .expect("stride-four run")
    });
    assert_eq!(one.points(), four.points());
    assert_ne!(one.report().polls, four.report().polls);
    assert_eq!(
        one.report().artifact_identity,
        four.report().artifact_identity,
        "artifact identity is invariant to equivalent polling chunking"
    );
}

#[test]
fn g5_streamline_identity_has_fixed_little_endian_field_order() {
    let mut request = StreamlineRequest::dimensionless_rk4([-0.0, 2.5], -0.125, 0, 7);
    request.domain = Some(StreamlineDomain2 {
        lower: [-1.0, -2.0],
        upper: [3.0, 4.0],
    });
    request.boundary_policy = StreamlineBoundaryPolicy::StopBeforeExit;
    request.stagnation_policy = StreamlineStagnationPolicy::StopBeforeRepeat;
    let budget = required_streamline_budget(request).expect("fixed identity request");
    let gate = CancelGate::new();
    let output = with_contour_cx(&gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [0.0; 2], request, budget).expect("zero-step identity fixture")
    });

    let mut expected = DomainHasher::new(STREAMLINE_ARTIFACT_IDENTITY_DOMAIN);
    expected.update(&STREAMLINE_ARTIFACT_IDENTITY_VERSION.to_le_bytes());
    expected.update(&[1]); // RK4.
    expected.update(&[1]); // Dimensionless.
    expected.update(&[2]); // StopBeforeExit.
    expected.update(&[2]); // StopBeforeRepeat.
    for coordinate in request.seed {
        expected.update(&coordinate.to_bits().to_le_bytes());
    }
    expected.update(&request.dt.to_bits().to_le_bytes());
    expected.update(&0u64.to_le_bytes()); // Requested steps.
    expected.update(&[1]); // Domain present.
    for coordinate in [-1.0_f64, -2.0, 3.0, 4.0] {
        expected.update(&coordinate.to_bits().to_le_bytes());
    }
    expected.update(&[1]); // StepsComplete.
    expected.update(&0u64.to_le_bytes()); // Termination step.
    expected.update(&0u64.to_le_bytes()); // Completed steps.
    expected.update(&0u64.to_le_bytes()); // Field evaluations.
    expected.update(&1u64.to_le_bytes()); // Published points.
    expected.update(&[0]); // No embedded error estimate.
    expected.update(&0u64.to_le_bytes());
    for coordinate in request.seed {
        expected.update(&coordinate.to_bits().to_le_bytes());
    }

    assert_eq!(output.report().artifact_identity, Some(expected.finalize()));
    assert_eq!(
        output.report().identity_bytes_hashed,
        budget.identity_byte_limit
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One G4 matrix shares the same finite admitted request.
fn g4_streamline_cancellation_budget_and_callback_faults_are_atomic() {
    let request = StreamlineRequest::dimensionless_rk4([0.0; 2], 0.125, 8, 1);
    let budget = required_streamline_budget(request).expect("finite cancellable plan");

    let cancelled_gate = CancelGate::new();
    cancelled_gate.request();
    let cancelled = with_contour_cx(&cancelled_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, 0.0], request, budget)
            .expect_err("pre-requested cancellation must precede allocation/callbacks")
    });
    assert!(matches!(
        cancelled.error,
        StreamlineError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::Cancelled { .. }
        }
    ));
    assert_eq!(
        cancelled.report.disposition,
        StreamlineDisposition::Cancelled
    );
    assert_eq!(cancelled.report.field_evaluations, 0);
    assert_eq!(cancelled.report.reserved_output_bytes, 0);
    assert!(!cancelled.report.published);

    let cost_quota = budget.work_unit_limit - 1;
    let cost_gate = CancelGate::new();
    let cost_refusal = with_contour_cx(
        &cost_gate,
        Budget::INFINITE.with_cost_quota(cost_quota),
        |cx| {
            streamline_with_cx(cx, |_| [1.0, 0.0], request, budget)
                .expect_err("ambient one-short cost plan refuses before callbacks")
        },
    );
    assert!(matches!(
        cost_refusal.error,
        StreamlineError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::CostPlanExceedsQuota { .. }
        }
    ));
    assert_eq!(cost_refusal.report.field_evaluations, 0);

    let deadline_gate = CancelGate::new();
    let missing_clock = with_contour_cx(&deadline_gate, Budget::with_deadline_at_ns(10), |cx| {
        streamline_with_cx(cx, |_| [1.0, 0.0], request, budget)
            .expect_err("deadline authority without a clock fails closed")
    });
    assert_eq!(
        missing_clock.error,
        StreamlineError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::DeadlineWithoutClock { deadline_ns: 10 }
        }
    );
    assert_eq!(missing_clock.report.polls, 0);
    assert_eq!(missing_clock.report.field_evaluations, 0);

    let poll_gate = CancelGate::new();
    let poll_refusal = with_contour_cx(&poll_gate, Budget::INFINITE.with_poll_quota(1), |cx| {
        streamline_with_cx(cx, |_| [1.0, 0.0], request, budget)
            .expect_err("second checkpoint exhausts the ambient quota")
    });
    assert!(matches!(
        poll_refusal.error,
        StreamlineError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::PollsExhausted { quota: 1, .. }
        }
    ));
    assert_eq!(poll_refusal.report.field_evaluations, 0);
    assert!(!poll_refusal.report.published);

    let panic_gate = CancelGate::new();
    let panicked = with_contour_cx(&panic_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| panic!("injected callback panic"), request, budget)
            .expect_err("callback unwind is contained")
    });
    assert_eq!(
        panicked.error,
        StreamlineError::CallbackPanicked {
            step: 0,
            stage: StreamlineStage::K1,
        }
    );
    assert_eq!(panicked.report.field_evaluations, 1);
    assert!(!panicked.report.published);

    let nonfinite_gate = CancelGate::new();
    let nonfinite = with_contour_cx(&nonfinite_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [f64::INFINITY, 0.0], request, budget)
            .expect_err("non-finite callback evidence refuses")
    });
    assert!(matches!(
        nonfinite.error,
        StreamlineError::NonFiniteFieldValue {
            step: 0,
            component: 0,
            ..
        }
    ));
    assert_eq!(nonfinite.report.field_evaluations, 1);
    assert!(!nonfinite.report.published);

    let overflow_request = StreamlineRequest::dimensionless_rk4([0.0; 2], f64::MAX, 1, 1);
    let overflow_budget =
        required_streamline_budget(overflow_request).expect("finite extreme request admits");
    let overflow_gate = CancelGate::new();
    let overflow = with_contour_cx(&overflow_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [f64::MAX, 0.0], overflow_request, overflow_budget)
            .expect_err("finite stage arithmetic overflow refuses")
    });
    assert!(matches!(
        overflow.error,
        StreamlineError::NonFiniteIntermediate {
            step: 0,
            component: 0,
            ..
        }
    ));
    assert!(!overflow.report.published);

    let replay_gate = CancelGate::new();
    let replay = with_contour_cx(&replay_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, -0.5], request, budget).expect("finite replay")
    });
    let direct_gate = CancelGate::new();
    let direct = with_contour_cx(&direct_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, -0.5], request, budget).expect("finite direct run")
    });
    assert_eq!(replay, direct);
}

#[test]
fn g4_streamline_callback_requested_cancellation_stops_at_the_next_step_boundary() {
    let request = StreamlineRequest::dimensionless_rk4([0.0; 2], 0.125, 8, 1);
    let budget = required_streamline_budget(request).expect("finite cancellable plan");
    let gate = CancelGate::new();
    let calls = Cell::new(0usize);

    let cancellation = with_contour_cx(&gate, Budget::INFINITE, |cx| {
        streamline_with_cx(
            cx,
            |_| {
                let call = calls.get() + 1;
                calls.set(call);
                if call == 2 {
                    // The request arrives during K2. One RK4 step is the
                    // documented indivisible tile, so K3/K4 finish before the
                    // next complete-step checkpoint observes cancellation.
                    gate.request();
                }
                [1.0, -0.5]
            },
            request,
            budget,
        )
        .expect_err("callback-requested cancellation must prevent publication")
    });

    assert_eq!(
        calls.get(),
        4,
        "the in-flight RK4 tile must drain exactly once"
    );
    assert!(matches!(
        cancellation.error,
        StreamlineError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::Cancelled { .. }
        }
    ));
    assert_eq!(cancellation.report.completed_steps, 1);
    assert_eq!(cancellation.report.field_evaluations, 4);
    assert_eq!(cancellation.report.output_points, 2);
    assert_eq!(
        cancellation.report.disposition,
        StreamlineDisposition::Cancelled
    );
    assert!(cancellation.report.terminal);
    assert!(!cancellation.report.published);
    assert!(cancellation.report.artifact_identity.is_none());

    let retry_gate = CancelGate::new();
    let retry = with_contour_cx(&retry_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, -0.5], request, budget)
            .expect("fresh-context retry after callback cancellation")
    });
    let direct_gate = CancelGate::new();
    let direct = with_contour_cx(&direct_gate, Budget::INFINITE, |cx| {
        streamline_with_cx(cx, |_| [1.0, -0.5], request, budget)
            .expect("independent direct execution")
    });
    assert_eq!(retry, direct);
}

#[test]
fn hessian_classification_recovers_the_morse_type() {
    let t = 1e-9;
    // f = x² + y²  -> minimum, index 0.
    assert_eq!(
        classify_hessian([[2.0, 0.0], [0.0, 2.0]], t).kind,
        CriticalKind::Minimum
    );
    // f = x² - y²  -> saddle, index 1.
    let s = classify_hessian([[2.0, 0.0], [0.0, -2.0]], t);
    assert_eq!(s.kind, CriticalKind::Saddle);
    assert_eq!(s.morse_index, 1);
    // f = -(x² + y²) -> maximum, index 2.
    assert_eq!(
        classify_hessian([[-2.0, 0.0], [0.0, -2.0]], t).morse_index,
        2
    );
    // f = xy -> saddle (off-diagonal Hessian, eigenvalues ±1).
    assert_eq!(
        classify_hessian([[0.0, 1.0], [1.0, 0.0]], t).kind,
        CriticalKind::Saddle
    );
    // a zero eigenvalue is degenerate.
    assert_eq!(
        classify_hessian([[2.0, 0.0], [0.0, 0.0]], t).kind,
        CriticalKind::Degenerate
    );
    // Positive scaling preserves inertia even when directly squaring the
    // finite off-diagonal entry would overflow.
    let normalized = classify_hessian([[1.0, 1e-108], [1e-108, 1.0]], t);
    let large = classify_hessian([[1e308, 1e200], [1e200, 1e308]], t);
    assert_eq!(large.kind, CriticalKind::Minimum);
    assert_eq!(large.morse_index, 0);
    assert_eq!(large, normalized);
    // Invalid numerics cannot manufacture a confident Morse type.
    let invalid = classify_hessian([[1.0, f64::INFINITY], [f64::INFINITY, 1.0]], t);
    assert_eq!(invalid.kind, CriticalKind::Degenerate);
    assert_eq!(invalid.morse_index, 0);
    assert_eq!(
        classify_hessian([[1.0, 0.0], [0.0, 1.0]], f64::NAN).kind,
        CriticalKind::Degenerate
    );
}

#[test]
fn hessian_symmetry_admission_is_explicit_and_fail_closed() {
    let t = 1e-9;
    let exact = [[3.0, 1.0], [1.0, 3.0]];
    assert_eq!(
        classify_hessian(exact, t),
        classify_hessian_with_policy(exact, t, HessianSymmetryPolicy::Exact)
    );
    assert_eq!(classify_hessian(exact, t).kind, CriticalKind::Minimum);
    let axis_swapped = [[2.0, 1.0], [1.0, 4.0]];
    assert_eq!(
        classify_hessian([[4.0, 1.0], [1.0, 2.0]], t),
        classify_hessian(axis_swapped, t)
    );

    // G0: a nominal Hessian with conflicting mixed partials has no symmetric
    // eigenvalue authority, even when one triangle alone looks positive.
    let asymmetric = [[1.0, 100.0], [-100.0, 1.0]];
    for refused in [asymmetric, [[1.0, -100.0], [100.0, 1.0]]] {
        assert_eq!(
            classify_hessian(refused, t),
            fs_viz::CriticalPoint {
                kind: CriticalKind::Degenerate,
                morse_index: 0,
            }
        );
    }

    // G3: the eigenvalue tolerance is not an implicit symmetry tolerance.
    // Both a near mismatch and its transpose fail closed, before and after an
    // exact power-of-two scale transformation.
    let next_after_one = f64::from_bits(1.0_f64.to_bits() + 1);
    let near = [[3.0, 1.0], [next_after_one, 3.0]];
    assert!((near[0][1] - near[1][0]).abs() < t);
    let near_transpose = [[3.0, next_after_one], [1.0, 3.0]];
    let scaled_near = [
        [near[0][0] * 2.0, near[0][1] * 2.0],
        [near[1][0] * 2.0, near[1][1] * 2.0],
    ];
    for refused in [near, near_transpose, scaled_near] {
        let classified = classify_hessian(refused, t);
        assert_eq!(classified.kind, CriticalKind::Degenerate);
        assert_eq!(classified.morse_index, 0);
    }

    // IEEE-754 signed zeros are numerically equal under Exact; changing which
    // triangle carries the sign bit cannot change the Morse evidence.
    let upper_negative_zero = classify_hessian([[2.0, -0.0], [0.0, 2.0]], t);
    let lower_negative_zero = classify_hessian([[2.0, 0.0], [-0.0, 2.0]], t);
    assert_eq!(upper_negative_zero.kind, CriticalKind::Minimum);
    assert_eq!(upper_negative_zero, lower_negative_zero);

    // The previously ignored lower triangle is also covered by numeric
    // admission; it cannot hide a non-finite value.
    for non_finite in [
        [[1.0, 0.0], [f64::INFINITY, 1.0]],
        [[1.0, 0.0], [f64::NAN, 1.0]],
    ] {
        let classified = classify_hessian(non_finite, t);
        assert_eq!(classified.kind, CriticalKind::Degenerate);
        assert_eq!(classified.morse_index, 0);
    }
    for invalid_tolerance in [f64::NAN, f64::INFINITY, -f64::EPSILON] {
        let classified = classify_hessian(exact, invalid_tolerance);
        assert_eq!(classified.kind, CriticalKind::Degenerate);
        assert_eq!(classified.morse_index, 0);
    }
}

#[test]
fn a_circle_sdf_isocontour_lies_on_the_circle() {
    // f(x,y) = sqrt(x²+y²) - 1, zero level set is the unit circle.
    let grid = Grid2::from_fn(41, 41, [-2.0, -2.0], [2.0, 2.0], 41 * 41, |p| {
        radius(p) - 1.0
    })
    .expect("finite circle grid within its exact node budget");
    let crossing_limit = 2 * 41 * 40;
    let crossings = grid
        .isocontour_crossings(0.0, crossing_limit)
        .expect("finite non-coincident circle crossings within edge budget");
    assert!(!crossings.is_empty());
    for c in &crossings {
        assert!(
            (radius(*c) - 1.0).abs() < 0.02,
            "crossing radius {}",
            radius(*c)
        );
    }
    // a level set outside the field's range has no crossings.
    assert!(
        grid.isocontour_crossings(100.0, crossing_limit)
            .expect("finite out-of-range levels are valid")
            .is_empty()
    );
}

#[test]
fn the_grid_samples_and_addresses_correctly() {
    let grid = Grid2::from_fn(3, 3, [-0.0, 0.0], [2.0, 2.0], 9, |p| p[0] + p[1])
        .expect("finite 3x3 grid within its exact node budget");
    let (p00, p22) = (grid.point(0, 0), grid.point(2, 2));
    assert_eq!(p00[0].to_bits(), (-0.0_f64).to_bits());
    assert!(p00[0].abs() < 1e-12 && p00[1].abs() < 1e-12);
    assert!((p22[0] - 2.0).abs() < 1e-12 && (p22[1] - 2.0).abs() < 1e-12);
    assert!((grid.at(1, 1) - 2.0).abs() < 1e-12); // (1,1) -> value 1+1
}

#[test]
fn grid2_layout_admission_precedes_sampling() {
    let calls = Cell::new(0usize);
    let mut sample = |_| {
        calls.set(calls.get() + 1);
        0.0
    };
    assert!(matches!(
        Grid2::from_fn(1, 2, [0.0; 2], [1.0; 2], 2, &mut sample),
        Err(Grid2Error::InvalidDimensions { dimensions: [1, 2] })
    ));
    assert!(matches!(
        Grid2::from_fn(2, 1, [0.0; 2], [1.0; 2], 2, &mut sample),
        Err(Grid2Error::InvalidDimensions { dimensions: [2, 1] })
    ));
    assert!(matches!(
        Grid2::from_fn(usize::MAX, 2, [0.0; 2], [1.0; 2], usize::MAX, &mut sample),
        Err(Grid2Error::NodeCountOverflow { .. })
    ));
    assert!(matches!(
        Grid2::from_fn(2, 2, [0.0; 2], [1.0; 2], 3, &mut sample),
        Err(Grid2Error::NodeBudgetExceeded {
            required: 4,
            limit: 3
        })
    ));

    let invalid_bounds = [
        ([f64::NAN, 0.0], [1.0, 1.0]),
        ([0.0, 0.0], [f64::INFINITY, 1.0]),
        ([0.0, 0.0], [0.0, 1.0]),
        ([1.0, 0.0], [0.0, 1.0]),
        ([-f64::MAX, 0.0], [f64::MAX, 1.0]),
    ];
    for (lo, hi) in invalid_bounds {
        assert!(matches!(
            Grid2::from_fn(2, 2, lo, hi, 4, &mut sample),
            Err(Grid2Error::InvalidBounds { axis: 0, .. })
        ));
    }
    let adjacent_to_one = 1.0_f64.next_up();
    assert!(matches!(
        Grid2::from_fn(3, 2, [1.0, 0.0], [adjacent_to_one, 1.0], 6, &mut sample),
        Err(Grid2Error::UnrepresentableCoordinates {
            axis: 0,
            first_index: 0,
            first,
            second_index: 1,
            second
        }) if first.to_bits() == 1.0_f64.to_bits()
            && second.to_bits() == 1.0_f64.to_bits()
    ));
    assert_eq!(calls.get(), 0, "invalid layouts must not invoke the field");
}

#[test]
fn grid2_rejects_the_first_nonfinite_sample() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let index = Cell::new(0usize);
        let result = Grid2::from_fn(3, 2, [0.0; 2], [1.0; 2], 6, |_| {
            let current = index.get();
            index.set(current + 1);
            if current == 4 { bad } else { current as f64 }
        });
        assert!(matches!(
            result,
            Err(Grid2Error::NonFiniteValue {
                index: 4,
                value
            }) if value.to_bits() == bad.to_bits()
        ));
        assert_eq!(index.get(), 5, "sampling stops at the first bad value");
    }
}

#[test]
fn isocontour_admission_distinguishes_invalid_from_empty() {
    let grid = Grid2::from_fn(2, 2, [-1.0; 2], [1.0; 2], 4, |p| p[0]).expect("finite affine grid");
    for iso in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            grid.isocontour_crossings(iso, 4),
            Err(IsoContourError::NonFiniteIso { iso: rejected })
                if rejected.to_bits() == iso.to_bits()
        ));
    }
    assert_eq!(
        grid.isocontour_crossings(0.0, 0),
        Err(IsoContourError::ZeroCrossingLimit)
    );
    assert_eq!(
        grid.isocontour_crossings(0.0, 1),
        Err(IsoContourError::CrossingBudgetExceeded { limit: 1 })
    );
    assert!(
        grid.isocontour_crossings(2.0, 1)
            .expect("a finite absent level is valid")
            .is_empty()
    );
}

#[test]
fn exact_level_nodes_are_unique_and_coincident_edges_are_refused() {
    let grid = Grid2::from_fn(3, 3, [-1.0; 2], [1.0; 2], 9, |p| p[0] + 2.0 * p[1])
        .expect("finite affine grid");
    let crossings = grid
        .isocontour_crossings(0.0, 16)
        .expect("isolated exact nodes have unique point intersections");
    assert_eq!(
        crossings
            .iter()
            .filter(|point| point[0] == 0.0 && point[1] == 0.0)
            .count(),
        1,
        "the exact center is shared by incident edges but emitted once"
    );

    let signed_zero = Grid2::from_fn(3, 2, [0.0; 2], [2.0, 1.0], 6, |p| {
        if p[0] == 1.0 && p[1] == 0.0 {
            -0.0
        } else {
            p[0] + p[1] - 1.0
        }
    })
    .expect("finite signed-zero grid");
    let signed_crossings = signed_zero
        .isocontour_crossings(0.0, 8)
        .expect("signed zero is one exact level node");
    assert_eq!(
        signed_crossings
            .iter()
            .filter(|point| point[0] == 1.0 && point[1] == 0.0)
            .count(),
        1
    );

    let plateau = Grid2::from_fn(3, 2, [-1.0, 0.0], [1.0, 1.0], 6, |p| p[0])
        .expect("finite grid containing a level-coincident vertical edge");
    assert_eq!(
        plateau.isocontour_crossings(0.0, 8),
        Err(IsoContourError::CoincidentLevelEdge {
            first: [1, 0],
            second: [1, 1]
        })
    );
}

#[test]
fn g0_checkerboard_exact_node_ownership_is_static_and_budget_exact() {
    for (nx, ny, exact_even) in [(65, 63, true), (63, 65, false)] {
        let sample_index = Cell::new(0usize);
        let grid = Grid2::from_fn(
            nx,
            ny,
            [0.0; 2],
            [(nx - 1) as f64, (ny - 1) as f64],
            nx * ny,
            |_| {
                let index = sample_index.get();
                sample_index.set(index + 1);
                let even = index % 2 == 0;
                if even == exact_even { 0.0 } else { 1.0 }
            },
        )
        .expect("checkerboard dimensions and samples are admitted");

        let is_exact = |i: usize, j: usize| ((j * nx + i) % 2 == 0) == exact_even;
        let mut expected = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                if i + 1 < nx && j == 0 {
                    if i == 0 && is_exact(0, 0) {
                        expected.push(grid.point(0, 0));
                    } else if is_exact(i + 1, 0) {
                        expected.push(grid.point(i + 1, 0));
                    }
                }
                if j + 1 < ny && is_exact(i, j + 1) {
                    expected.push(grid.point(i, j + 1));
                }
            }
        }

        let parity_tail = if exact_even { 1 } else { 0 };
        assert_eq!(expected.len(), (nx * ny + parity_tail) / 2);
        assert_eq!(
            grid.isocontour_crossings(0.0, expected.len())
                .expect("the exact output budget admits every canonical owner"),
            expected,
            "static ownership must retain first-incident-edge traversal order"
        );
        assert_eq!(
            grid.isocontour_crossings(0.0, expected.len())
                .expect("replay uses the same static owners"),
            expected
        );

        let scoped_budget = grid
            .required_isocontour_budget(expected.len(), 128)
            .expect("checkerboard work envelope is representable");
        let scoped_gate = CancelGate::new();
        let scoped = with_contour_cx(&scoped_gate, Budget::INFINITE, |cx| {
            grid.isocontour_crossings_with_cx(cx, 0.0, scoped_budget)
                .expect("the exact checkerboard envelope admits linear ownership")
        });
        let scoped_plan = scoped.report().plan.expect("scoped plan retained");
        assert_eq!(scoped.crossings(), expected.as_slice());
        assert_eq!(scoped.report().edge_visits, scoped_plan.edge_visits);
        assert_eq!(
            scoped.report().exact_ownership_checks,
            scoped_plan.edge_visits,
            "every checkerboard edge performs one constant-time ownership decision"
        );
        assert_eq!(scoped.report().exact_crossings, expected.len());
        assert_eq!(scoped.report().interpolated_crossings, 0);
        assert_eq!(
            scoped.report().unowned_exact_candidates,
            scoped.report().exact_ownership_checks - scoped.report().exact_crossings,
            "every non-owner exact candidate is discarded in constant time"
        );
        assert_eq!(
            scoped.report().crossings,
            scoped.report().exact_crossings + scoped.report().interpolated_crossings
        );
        assert_eq!(scoped.report().interpolations, 0);
        assert!(scoped.report().work_units <= scoped_plan.work_units);

        let one_short = expected.len() - 1;
        assert_eq!(
            grid.isocontour_crossings(0.0, one_short),
            Err(IsoContourError::CrossingBudgetExceeded { limit: one_short })
        );
        let one_short_budget = grid
            .required_isocontour_budget(one_short, 128)
            .expect("one-short checkerboard envelope is representable");
        let one_short_gate = CancelGate::new();
        let one_short_refusal = with_contour_cx(&one_short_gate, Budget::INFINITE, |cx| {
            grid.isocontour_crossings_with_cx(cx, 0.0, one_short_budget)
                .expect_err("one extra canonical owner must refuse atomically")
        });
        assert_eq!(
            one_short_refusal.error,
            IsoContourError::CrossingBudgetExceeded { limit: one_short }
        );
        assert_eq!(one_short_refusal.report.crossings, one_short);
        assert_eq!(one_short_refusal.report.exact_crossings, one_short);
        assert_eq!(one_short_refusal.report.interpolated_crossings, 0);
        assert!(!one_short_refusal.report.published);
    }
}

#[test]
fn g0_scoped_contour_plan_is_complete_and_reported() {
    let grid = Grid2::from_fn(5, 4, [0.0; 2], [4.0, 3.0], 20, |point| point[0] - 1.3)
        .expect("finite affine Grid2");
    let budget = grid
        .required_isocontour_budget(32, 3)
        .expect("checked full-grid envelope");
    let plan = grid
        .isocontour_plan(budget)
        .expect("the exact derived budget admits its plan");
    assert_eq!(plan.dimensions, [5, 4]);
    assert_eq!(plan.nodes, 20);
    assert_eq!(plan.cells, 12);
    assert_eq!(plan.edge_visits, 31);
    assert_eq!(plan.exact_ownership_checks, 31);
    assert_eq!(plan.interpolations, 31);
    assert_eq!(plan.output_bytes, 32 * size_of::<[f64; 2]>());
    assert_eq!(plan.polls, 2 + 31usize.div_ceil(3) + 32usize.div_ceil(3));

    let gate = CancelGate::new();
    let first = with_contour_cx(&gate, Budget::INFINITE, |cx| {
        grid.isocontour_crossings_with_cx(cx, 0.0, budget)
            .expect("scoped affine contour")
    });
    let report = *first.report();
    assert_eq!(report.operation_budget, Some(budget));
    assert_eq!(report.plan, Some(plan));
    assert_eq!(report.node_visits, plan.nodes);
    assert_eq!(report.cell_visits, plan.cells);
    assert_eq!(report.edge_visits, plan.edge_visits);
    assert_eq!(report.crossings, first.crossings().len());
    assert_eq!(report.exact_crossings, 0);
    assert_eq!(report.unowned_exact_candidates, 0);
    assert_eq!(report.interpolated_crossings, first.crossings().len());
    assert_eq!(
        report.crossings,
        report.exact_crossings + report.interpolated_crossings
    );
    assert!(report.polls <= plan.polls);
    assert!(report.work_units <= plan.work_units);
    assert!(report.identity_bytes_hashed <= plan.identity_bytes);
    assert_eq!(report.diagnostic_records, 1);
    assert!(report.terminal && report.published);
    assert_eq!(report.disposition, IsoContourDisposition::Completed);
    assert!(report.artifact_identity.is_some());
    assert_eq!(first.crossings().len(), 4);

    let replay_gate = CancelGate::new();
    let replay = with_contour_cx(&replay_gate, Budget::INFINITE, |cx| {
        grid.isocontour_crossings_with_cx(cx, 0.0, budget)
            .expect("scoped replay")
    });
    assert_eq!(replay, first);
    assert_eq!(
        grid.isocontour_crossings(0.0, budget.crossing_limit)
            .expect("compatibility result")
            .as_slice(),
        first.crossings()
    );
    assert!(!ISO_CONTOUR_ARTIFACT_IDENTITY_DOMAIN.is_empty());
}

#[test]
fn g0_contour_resources_refuse_before_poll_or_edge_work() {
    let grid = Grid2::from_fn(7, 6, [-1.0; 2], [1.0; 2], 42, |point| {
        point[0] + 0.3 * point[1]
    })
    .expect("finite affine Grid2");
    let budget = grid
        .required_isocontour_budget(40, 5)
        .expect("checked exact budget");
    let mut cases = Vec::new();
    macro_rules! one_short {
        ($resource:expr, $field:ident) => {{
            let mut limited = budget;
            limited.$field -= 1;
            cases.push(($resource, limited));
        }};
    }
    one_short!(IsoContourResource::Cells, cell_limit);
    one_short!(IsoContourResource::EdgeVisits, edge_visit_limit);
    one_short!(
        IsoContourResource::ExactOwnershipChecks,
        exact_ownership_limit
    );
    one_short!(IsoContourResource::Interpolations, interpolation_limit);
    one_short!(IsoContourResource::OutputBytes, output_byte_limit);
    one_short!(IsoContourResource::ScratchBytes, scratch_byte_limit);
    one_short!(
        IsoContourResource::DiagnosticRecords,
        diagnostic_record_limit
    );
    one_short!(IsoContourResource::DiagnosticBytes, diagnostic_byte_limit);
    one_short!(IsoContourResource::LiveBytes, live_byte_limit);
    one_short!(IsoContourResource::IdentityBytes, identity_byte_limit);
    one_short!(IsoContourResource::Polls, poll_limit);
    one_short!(IsoContourResource::WorkUnits, work_unit_limit);

    for (resource, limited) in cases {
        let gate = CancelGate::new();
        let refusal = with_contour_cx(&gate, Budget::INFINITE, |cx| {
            grid.isocontour_crossings_with_cx(cx, 0.0, limited)
                .expect_err("one-short complete resource must refuse")
        });
        assert!(matches!(
            refusal.error,
            IsoContourError::OperationBudgetExceeded {
                resource: rejected,
                ..
            } if rejected == resource
        ));
        assert_eq!(refusal.report.polls, 0, "{resource} refused after a poll");
        assert_eq!(
            refusal.report.edge_visits, 0,
            "{resource} refused after edge work"
        );
        assert!(!refusal.report.published);
        assert_eq!(refusal.report.diagnostic_records, 1);
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One transactional matrix shares an exact checked plan.
fn g4_scoped_contour_cancellation_and_poll_refusal_are_transactional() {
    let grid = Grid2::from_fn(17, 17, [-1.0; 2], [1.0; 2], 17 * 17, |point| {
        point[0] + point[1] - 0.03
    })
    .expect("finite affine Grid2");
    let budget = grid
        .required_isocontour_budget(64, 2)
        .expect("checked cancellable envelope");

    let cancelled_gate = CancelGate::new();
    cancelled_gate.request();
    let cancelled = with_contour_cx(&cancelled_gate, Budget::INFINITE, |cx| {
        grid.isocontour_crossings_with_cx(cx, 0.0, budget)
            .expect_err("pre-requested cancellation must refuse before allocation")
    });
    assert!(matches!(
        cancelled.error,
        IsoContourError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::Cancelled { .. }
        }
    ));
    assert_eq!(
        cancelled.report.disposition,
        IsoContourDisposition::Cancelled
    );
    assert_eq!(cancelled.report.polls, 1);
    assert_eq!(cancelled.report.edge_visits, 0);
    assert_eq!(cancelled.report.reserved_output_bytes, 0);
    assert!(!cancelled.report.published);

    let deadline_gate = CancelGate::new();
    let missing_clock = with_contour_cx(&deadline_gate, Budget::with_deadline_at_ns(10), |cx| {
        grid.isocontour_crossings_with_cx(cx, 0.0, budget)
            .expect_err("an ambient deadline without a clock must fail closed")
    });
    assert_eq!(
        missing_clock.error,
        IsoContourError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::DeadlineWithoutClock { deadline_ns: 10 }
        }
    );
    assert_eq!(missing_clock.report.polls, 0);
    assert_eq!(missing_clock.report.reserved_output_bytes, 0);

    let cost_quota = budget.work_unit_limit - 1;
    let cost_gate = CancelGate::new();
    let cost_refusal = with_contour_cx(
        &cost_gate,
        Budget::INFINITE.with_cost_quota(cost_quota),
        |cx| {
            grid.isocontour_crossings_with_cx(cx, 0.0, budget)
                .expect_err("an ambient one-short work quota must refuse at admission")
        },
    );
    assert_eq!(
        cost_refusal.error,
        IsoContourError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::CostPlanExceedsQuota {
                planned: budget.work_unit_limit,
                quota: cost_quota,
            }
        }
    );
    assert_eq!(cost_refusal.report.polls, 0);
    assert_eq!(cost_refusal.report.edge_visits, 0);
    assert_eq!(cost_refusal.report.reserved_output_bytes, 0);

    let quota_gate = CancelGate::new();
    let exhausted = with_contour_cx(&quota_gate, Budget::INFINITE.with_poll_quota(3), |cx| {
        grid.isocontour_crossings_with_cx(cx, 0.0, budget)
            .expect_err("the exact ambient poll quota must stop a later edge chunk")
    });
    assert!(matches!(
        exhausted.error,
        IsoContourError::ExecutionBudgetRefused {
            refusal: BudgetRefusal::PollsExhausted { quota: 3, .. }
        }
    ));
    assert_eq!(exhausted.report.edge_visits, 4);
    assert!(!exhausted.report.published);
    assert!(exhausted.report.artifact_identity.is_none());

    let retry_gate = CancelGate::new();
    let retry = with_contour_cx(&retry_gate, Budget::INFINITE, |cx| {
        grid.isocontour_crossings_with_cx(cx, 0.0, budget)
            .expect("retry under a fresh Cx")
    });
    let direct_gate = CancelGate::new();
    let direct = with_contour_cx(&direct_gate, Budget::INFINITE, |cx| {
        grid.isocontour_crossings_with_cx(cx, 0.0, budget)
            .expect("direct reference under a fresh Cx")
    });
    assert_eq!(retry, direct, "retry must be byte-identical to direct work");
}

#[test]
fn isocontour_interpolation_handles_extreme_finite_values() {
    for magnitude in [f64::MAX, f64::from_bits(1)] {
        let grid = Grid2::from_fn(2, 2, [0.0; 2], [1.0; 2], 4, |p| {
            if p[0] == 0.0 { -magnitude } else { magnitude }
        })
        .expect("finite extreme samples are admissible");
        let crossings = grid
            .isocontour_crossings(0.0, 2)
            .expect("scaled interpolation remains finite");
        assert_eq!(crossings, vec![[0.5, 0.0], [0.5, 1.0]]);
    }
}

#[test]
fn g0_isocontour_refuses_strict_crossings_that_round_to_an_endpoint() {
    let lower = 1.0_f64;
    let upper = lower.next_up();
    let tiny = f64::from_bits(1);
    let error = lower_left_collapse_error(lower, upper, -tiny, 1.0, 0.0, 1);
    assert_eq!(
        lower_left_collapse_error(lower, upper, -tiny, 1.0, 0.0, 2),
        error,
        "crossing budget cannot replace the earlier representability refusal"
    );

    let IsoContourError::UnrepresentableIntersection {
        first,
        second,
        first_point_bits,
        second_point_bits,
        first_value_bits,
        second_value_bits,
        iso_bits,
        first_distance_bits,
        second_distance_bits,
        interpolation_bits,
        point_bits,
        collapsed_axis,
    } = error
    else {
        panic!("strict crossing collapse must return its typed evidence: {error:?}")
    };
    assert_eq!(first, [0, 0]);
    assert_eq!(second, [1, 0]);
    assert_eq!(first_point_bits, [lower.to_bits(), lower.to_bits()]);
    assert_eq!(second_point_bits, [upper.to_bits(), lower.to_bits()]);
    assert_eq!(first_value_bits, (-tiny).to_bits());
    assert_eq!(second_value_bits, 1.0_f64.to_bits());
    assert_eq!(iso_bits, 0.0_f64.to_bits());
    assert_eq!(first_distance_bits, tiny.to_bits());
    assert_eq!(second_distance_bits, 1.0_f64.to_bits());
    assert_eq!(interpolation_bits, tiny.to_bits());
    assert_eq!(point_bits, first_point_bits);
    assert_eq!(collapsed_axis, 0);
}

#[test]
fn g3_unrepresentable_intersection_refusal_tracks_axis_sign_and_scale_neighbors() {
    let tiny = f64::from_bits(1);
    let lower = 1.0_f64;
    let upper = lower.next_up();
    let horizontal = lower_left_collapse_error(lower, upper, -tiny, 1.0, 0.0, 1);

    let vertical_grid = Grid2::from_fn(2, 2, [lower; 2], [upper; 2], 4, |point| {
        if point[1].to_bits() == lower.to_bits() {
            -tiny
        } else {
            1.0
        }
    })
    .expect("axis-permuted adjacent grid admits");
    let vertical = vertical_grid
        .isocontour_crossings(0.0, 1)
        .expect_err("axis-permuted crossing must refuse identically");
    let (
        IsoContourError::UnrepresentableIntersection {
            first: horizontal_first,
            second: horizontal_second,
            interpolation_bits: horizontal_t,
            collapsed_axis: horizontal_axis,
            ..
        },
        IsoContourError::UnrepresentableIntersection {
            first: vertical_first,
            second: vertical_second,
            interpolation_bits: vertical_t,
            collapsed_axis: vertical_axis,
            ..
        },
    ) = (horizontal, vertical)
    else {
        panic!("axis permutations must retain typed representability evidence")
    };
    assert_eq!((horizontal_first, horizontal_second), ([0, 0], [1, 0]));
    assert_eq!((vertical_first, vertical_second), ([0, 0], [0, 1]));
    assert_eq!((horizontal_axis, vertical_axis), (0, 1));
    assert_eq!(horizontal_t, vertical_t);

    for (case, lower, upper, small, iso) in [
        ("next-up/min-subnormal", 1.0, 1.0_f64.next_up(), tiny, 0.0),
        (
            "next-down/min-normal/signed-zero",
            1.0_f64.next_down(),
            1.0,
            f64::MIN_POSITIVE,
            -0.0,
        ),
        (
            "power-of-two scale neighbor",
            2.0,
            2.0_f64.next_up(),
            tiny,
            0.0,
        ),
    ] {
        assert!(
            matches!(
                lower_left_collapse_error(lower, upper, -small, 1.0, iso, 1),
                IsoContourError::UnrepresentableIntersection {
                    collapsed_axis: 0,
                    ..
                }
            ),
            "{case} must fail closed as unrepresentable"
        );
    }

    assert!(matches!(
        lower_left_collapse_error(lower, upper, tiny, -1.0, -0.0, 1),
        IsoContourError::UnrepresentableIntersection {
            first_value_bits,
            second_value_bits,
            iso_bits,
            collapsed_axis: 0,
            ..
        } if first_value_bits == tiny.to_bits()
            && second_value_bits == (-1.0_f64).to_bits()
            && iso_bits == (-0.0_f64).to_bits()
    ));
    assert!(matches!(
        lower_left_collapse_error(lower, upper, 1.0, -tiny, 0.0, 1),
        IsoContourError::UnrepresentableIntersection {
            interpolation_bits,
            point_bits,
            second_point_bits,
            collapsed_axis: 0,
            ..
        } if interpolation_bits == 1.0_f64.to_bits() && point_bits == second_point_bits
    ));
}

#[test]
fn g3_isocontour_value_transformations_preserve_crossings() {
    let base = Grid2::from_fn(4, 4, [-1.0; 2], [1.0; 2], 16, |p| p[0] + 0.5 * p[1])
        .expect("finite affine base grid");
    let inverted = Grid2::from_fn(4, 4, [-1.0; 2], [1.0; 2], 16, |p| -(p[0] + 0.5 * p[1]))
        .expect("finite sign-inverted grid");
    let scaled = Grid2::from_fn(4, 4, [-1.0; 2], [1.0; 2], 16, |p| 8.0 * (p[0] + 0.5 * p[1]))
        .expect("finite power-of-two-scaled grid");

    let expected = base
        .isocontour_crossings(0.1, 24)
        .expect("base affine contour");
    assert_eq!(
        inverted
            .isocontour_crossings(-0.1, 24)
            .expect("sign inversion contour"),
        expected
    );
    assert_eq!(
        scaled
            .isocontour_crossings(0.8, 24)
            .expect("positive scaling contour"),
        expected
    );
}

#[test]
fn visualization_is_deterministic() {
    let a = streamline(|p| [-p[1], p[0]], [1.0, 0.0], 0.01, 100);
    let b = streamline(|p| [-p[1], p[0]], [1.0, 0.0], 0.01, 100);
    assert_eq!(a.len(), b.len());
    assert_eq!(
        a.last().unwrap()[0].to_bits(),
        b.last().unwrap()[0].to_bits()
    );
}

#[test]
fn marching_tetrahedra_extracts_an_exact_oriented_plane() {
    let dimensions = [9, 10, 11];
    let node_limit = dimensions.into_iter().product();
    let grid = Grid3::from_fn(dimensions, [-1.0; 3], [1.0; 3], node_limit, |point| {
        point[0] - 0.13
    })
    .expect("bounded finite plane grid");
    assert_eq!(grid.dimensions(), dimensions);
    assert!((grid.at(0, 0, 0).expect("in bounds") + 1.13).abs() < 1e-15);
    let upper = grid.point(8, 9, 10).expect("upper node is in bounds");
    assert!(
        upper
            .into_iter()
            .all(|coordinate| (coordinate - 1.0).abs() < 1e-15)
    );
    assert_eq!(grid.point(9, 0, 0), None);

    let mesh = grid.isosurface(0.0, 10_000).expect("plane isosurface");
    assert!(!mesh.triangles().is_empty());
    assert!(mesh.vertices().len() < mesh.triangles().len() * 3);
    assert!((mesh.surface_area() - 4.0).abs() < 1e-12);
    for vertex in mesh.vertices() {
        assert!((vertex[0] - 0.13).abs() < 1e-12);
    }
    for triangle in mesh.triangles() {
        let a = mesh.vertices()[triangle[0] as usize];
        let b = mesh.vertices()[triangle[1] as usize];
        let c = mesh.vertices()[triangle[2] as usize];
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let normal_x = ab[1] * ac[2] - ab[2] * ac[1];
        assert!(
            normal_x > 0.0,
            "plane triangle must point toward increasing field"
        );
    }
    assert!(matches!(
        grid.isosurface(0.0, 1),
        Err(IsoSurfaceError::TriangleBudgetExceeded { limit: 1 })
    ));
}

#[test]
fn sphere_isosurface_area_converges_under_refinement() {
    let radius = 0.7;
    let sphere = |resolution: usize| {
        let dimensions = [resolution; 3];
        let node_limit = dimensions.into_iter().product();
        Grid3::from_fn(dimensions, [-1.2; 3], [1.2; 3], node_limit, |point| {
            point[0]
                .mul_add(point[0], point[1].mul_add(point[1], point[2] * point[2]))
                .sqrt()
                - radius
        })
        .expect("bounded finite sphere grid")
        .isosurface(0.0, 200_000)
        .expect("sphere isosurface")
    };
    let coarse = sphere(17);
    let fine = sphere(33);
    let exact_area = 4.0 * std::f64::consts::PI * radius * radius;
    let coarse_error = (coarse.surface_area() - exact_area).abs();
    let fine_error = (fine.surface_area() - exact_area).abs();
    assert!(
        fine_error < coarse_error,
        "sphere area must converge: coarse {coarse_error:.3e}, fine {fine_error:.3e}"
    );
    assert!(fine_error / exact_area < 0.03);

    // Negative values are inside the sphere, so outward winding must align
    // each nondegenerate face normal with its centroid radius.
    for triangle in fine.triangles() {
        let a = fine.vertices()[triangle[0] as usize];
        let b = fine.vertices()[triangle[1] as usize];
        let c = fine.vertices()[triangle[2] as usize];
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let normal = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        let centroid = [
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ];
        let orientation = normal[0].mul_add(
            centroid[0],
            normal[1].mul_add(centroid[1], normal[2] * centroid[2]),
        );
        assert!(orientation > 0.0);
    }
}

#[test]
fn gyroid_extraction_is_indexed_symmetric_and_deterministic() {
    let dimensions = [19; 3];
    let node_limit = dimensions.into_iter().product();
    let bound = std::f64::consts::PI;
    let grid = Grid3::from_fn(dimensions, [-bound; 3], [bound; 3], node_limit, |point| {
        point[0].sin() * point[1].cos()
            + point[1].sin() * point[2].cos()
            + point[2].sin() * point[0].cos()
    })
    .expect("bounded finite gyroid grid");
    let first = grid.isosurface(0.0, 100_000).expect("gyroid surface");
    let replay = grid.isosurface(0.0, 100_000).expect("gyroid replay");
    assert_eq!(first, replay);
    assert!(!first.triangles().is_empty());
    assert!(first.vertices().len() < first.triangles().len() * 3);

    let mut lower = [f64::INFINITY; 3];
    let mut upper = [f64::NEG_INFINITY; 3];
    for vertex in first.vertices() {
        for axis in 0..3 {
            lower[axis] = lower[axis].min(vertex[axis]);
            upper[axis] = upper[axis].max(vertex[axis]);
        }
    }
    for axis in 0..3 {
        assert!((lower[axis] + upper[axis]).abs() < 1e-12);
    }
}

#[test]
fn grid3_admission_fails_before_unbounded_or_nonfinite_work() {
    let calls = std::cell::Cell::new(0usize);
    let over_budget = Grid3::from_fn([100, 100, 100], [-1.0; 3], [1.0; 3], 1_000, |_| {
        calls.set(calls.get() + 1);
        0.0
    });
    assert!(matches!(
        over_budget,
        Err(Grid3Error::NodeBudgetExceeded {
            required: 1_000_000,
            limit: 1_000
        })
    ));
    assert_eq!(calls.get(), 0);
    assert!(matches!(
        Grid3::from_values([2, 2, 2], [-1.0; 3], [1.0; 3], 8, vec![0.0; 7]),
        Err(Grid3Error::ValueCountMismatch {
            expected: 8,
            actual: 7
        })
    ));
    assert!(matches!(
        Grid3::from_fn([2, 2, 2], [-1.0; 3], [1.0; 3], 8, |_| f64::NAN),
        Err(Grid3Error::NonFiniteValue { index: 0, .. })
    ));
    let grid = Grid3::from_fn([2, 2, 2], [-1.0; 3], [1.0; 3], 8, |point| point[0])
        .expect("small admitted grid");
    assert!(matches!(
        grid.isosurface(f64::INFINITY, 10),
        Err(IsoSurfaceError::NonFiniteIso { .. })
    ));
    assert!(matches!(
        grid.isosurface(0.0, 0),
        Err(IsoSurfaceError::ZeroTriangleLimit)
    ));
}

fn density_semantics() -> ScalarFieldSemantics {
    ScalarFieldSemantics {
        quantity: "density".to_string(),
        coordinate_unit: "m".to_string(),
        value_unit: "kg/m^3".to_string(),
    }
}

#[test]
fn scalar_field_artifact_round_trips_bit_exactly_into_node_viz() {
    assert_eq!(SCALAR_FIELD3_ARTIFACT_KIND, "frankensim.scalar-field3");
    assert_eq!(SCALAR_FIELD3_SCHEMA_VERSION, 1);
    let values: Vec<f64> = (0..12).map(|index| f64::from(index) * 0.125).collect();
    let field = ScalarField3::new(
        ScalarLayout3::NodeCentered,
        [3, 2, 2],
        [-1.0, -2.0, 0.0],
        [0.5, 2.0, 1.0],
        density_semantics(),
        12,
        values,
    )
    .expect("valid node field");
    assert_eq!(field.world_bounds(), [[-1.0, -2.0, 0.0], [0.0, 0.0, 1.0]]);
    let encoded = field.encode(4096).expect("bounded encode");
    assert_eq!(encoded, field.encode(4096).expect("replay encode"));
    let decoded = ScalarField3::decode(&encoded, 12, encoded.len()).expect("bounded decode");
    assert_eq!(decoded, field);
    assert_eq!(decoded.encode(encoded.len()).expect("re-encode"), encoded);
    let grid = decoded.into_node_grid(12).expect("node-grid conversion");
    assert_eq!(grid.dimensions(), [3, 2, 2]);
    assert_eq!(grid.bounds(), [[-1.0, -2.0, 0.0], [0.0, 0.0, 1.0]]);
    assert_eq!(grid.at(2, 1, 1), Some(11.0 * 0.125));
}

#[test]
fn scalar_field_artifact_keeps_one_cell_thick_lbm_layout_honest() {
    let values: Vec<f64> = (0..12).map(|index| f64::from(index) / 11.0).collect();
    let field = ScalarField3::new(
        ScalarLayout3::CellCentered,
        [4, 3, 1],
        [0.0; 3],
        [1.0, 1.0, 24.0],
        ScalarFieldSemantics {
            quantity: "liquid_mass_fraction".to_string(),
            coordinate_unit: "cell".to_string(),
            value_unit: "1".to_string(),
        },
        12,
        values,
    )
    .expect("valid one-cell-thick field");
    assert_eq!(field.world_bounds(), [[0.0; 3], [4.0, 3.0, 24.0]]);
    let encoded = field.encode(4096).expect("bounded encode");
    let decoded = ScalarField3::decode(&encoded, 12, 4096).expect("bounded decode");
    assert_eq!(decoded.layout(), ScalarLayout3::CellCentered);
    assert_eq!(decoded.dimensions(), [4, 3, 1]);
    assert_eq!(decoded.origin(), [0.0; 3]);
    assert_eq!(decoded.spacing(), [1.0, 1.0, 24.0]);
    assert_eq!(decoded.semantics().value_unit, "1");
    assert!(matches!(
        decoded.into_node_grid(12),
        Err(ScalarField3Error::NotNodeCentered)
    ));
}

#[test]
fn scalar_field_codec_refuses_before_unbounded_or_ambiguous_work() {
    assert!(matches!(
        ScalarField3::new(
            ScalarLayout3::NodeCentered,
            [2, 2, 2],
            [0.0; 3],
            [1.0; 3],
            density_semantics(),
            7,
            vec![0.0; 8],
        ),
        Err(ScalarField3Error::SampleBudgetExceeded {
            required: 8,
            limit: 7,
        })
    ));
    assert!(matches!(
        ScalarField3::new(
            ScalarLayout3::CellCentered,
            [1, 1, 1],
            [0.0; 3],
            [1.0; 3],
            ScalarFieldSemantics {
                quantity: String::new(),
                coordinate_unit: "m".to_string(),
                value_unit: "1".to_string(),
            },
            1,
            vec![0.0],
        ),
        Err(ScalarField3Error::InvalidSemantic { field: "quantity" })
    ));
    let field = ScalarField3::new(
        ScalarLayout3::CellCentered,
        [1, 1, 1],
        [0.0; 3],
        [1.0; 3],
        density_semantics(),
        1,
        vec![0.5],
    )
    .expect("small valid field");
    let encoded = field.encode(1024).expect("bounded encode");
    assert!(matches!(
        field.encode(encoded.len() - 1),
        Err(ScalarField3Error::ByteBudgetExceeded { .. })
    ));
    assert!(matches!(
        ScalarField3::decode(&encoded, 1, encoded.len() - 1),
        Err(ScalarField3Error::ByteBudgetExceeded { .. })
    ));
    assert!(matches!(
        ScalarField3::decode(&encoded[..encoded.len() - 1], 1, encoded.len()),
        Err(ScalarField3Error::Malformed { .. })
    ));

    let mut bad_magic = encoded.clone();
    bad_magic[0] ^= 1;
    assert!(matches!(
        ScalarField3::decode(&bad_magic, 1, bad_magic.len()),
        Err(ScalarField3Error::Malformed { what: "bad magic" })
    ));
    let mut future = encoded.clone();
    future[8..12].copy_from_slice(&(SCALAR_FIELD3_SCHEMA_VERSION + 1).to_le_bytes());
    assert!(matches!(
        ScalarField3::decode(&future, 1, future.len()),
        Err(ScalarField3Error::UnsupportedSchema { found: 2 })
    ));
    let mut nonfinite = encoded;
    let tail = nonfinite.len() - 8;
    nonfinite[tail..].copy_from_slice(&f64::NAN.to_bits().to_le_bytes());
    assert!(matches!(
        ScalarField3::decode(&nonfinite, 1, nonfinite.len()),
        Err(ScalarField3Error::NonFiniteValue { index: 0 })
    ));
}
