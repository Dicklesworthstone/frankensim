//! Restoration-path resource-contract gates (bead
//! frankensim-constraint-restoration-budget-receipts-x5sev): checked
//! [`RestorationWorkPlan`] arithmetic, budget/cancellation authority
//! enforced at deterministic tile boundaries, retained receipts with NO
//! partial evidence, honest memory boundaries, and bit-stable replay.
//!
//! Coverage map:
//! - G0: cap/zero/overflow refusals, duplicate/out-of-range skips,
//!   plan-identity binding, tampered aggregates, schema mismatches,
//!   admission-time refusals (cost quota, zero polls, deadline
//!   without clock, expired deadline).
//! - G3/G5: replay stability of reports AND receipts, skip-order
//!   equivalence sharing one plan identity.
//! - G4: pre-admission vs admitted-stop taxonomy, mid-run cancellation
//!   with bounded extra work, mid-run deadline expiry, no-partial-
//!   success receipts, memory no-lease/no-claim boundary.
//!
//! Every stop is driven by deterministic counters or deterministic
//! clocks — never wall-clock sleeps.

#![allow(clippy::float_cmp)]

use asupersync::time::TimeSource;
use asupersync::types::{Budget, Time};
use fs_constraint::{
    ConError, ConstraintKind, ConstraintSpec, DomainBox, RestorationError, RestorationWorkLimits,
    RestorationWorkPlan, RestorationWorkShape, RESTORATION_MAX_FEASIBILITY_SAMPLES,
    RESTORATION_MAX_STARTS, RESTORATION_MAX_STEPS_PER_START, RESTORATION_WORK_PLAN_SCHEMA_VERSION,
    diagnose_infeasibility, elastic_solve, elastic_solve_with_plan,
};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::{Manifold, NodeId, Problem, ProblemBuilder};
use fs_qty::Dims;
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};

const EXECUTION_SEED: u64 = 0xC07A;

/// One Rn(n) host with objective `|x|^2`; nodes are `coeff * x_k - rhs`
/// hinge constraints (`g <= 0` semantics).
fn host_with(nodes_spec: &[(f64, u32, f64)], dim: u32) -> (Problem, Vec<NodeId>) {
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim }, Dims::NONE)
        .expect("var");
    let vr = b.var_ref(v).expect("ref");
    let mut nodes = Vec::new();
    for &(coeff, component, rhs) in nodes_spec {
        let xc = b.component(vr, component).expect("component");
        let c = b.konst(coeff, Dims::NONE).expect("coeff");
        let scaled = b.mul(c, xc).expect("scaled");
        let r = b.konst(rhs, Dims::NONE).expect("rhs");
        nodes.push(b.sub(scaled, r).expect("g"));
    }
    let obj = b.norm_sq(vr).expect("obj");
    b.objective(obj, fs_opt::Sense::Minimize, 1.0).expect("o");
    (b.finish(), nodes)
}

fn hard(name: &str, node: NodeId) -> ConstraintSpec {
    ConstraintSpec {
        name: name.to_string(),
        node,
        kind: ConstraintKind::Hard,
        active_tol: 1e-9,
    }
}

fn infinite_budget() -> Budget {
    Budget {
        deadline: None,
        poll_quota: u32::MAX,
        cost_quota: None,
        priority: 0,
    }
}

/// Run `f` with a hand-assembled context: NO operation-memory lease, an
/// optional deterministic clock, and the caller's gate/budget.
fn with_cx_raw<R>(
    gate: &CancelGate,
    budget: Budget,
    clock: Option<&dyn TimeSource>,
    f: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let stream = StreamKey {
            seed: EXECUTION_SEED,
            kernel_id: 9,
            tile: 0,
            iteration: 0,
        };
        let mut cx = Cx::new(gate, arena, stream, budget, ExecMode::Deterministic);
        if let Some(clock) = clock {
            cx = cx.with_time_source(clock);
        }
        f(&cx)
    })
}

/// Deterministic clock that requests cancellation on its `after`-th
/// observation. Admission consumes one observation first.
struct CancelAfterObservations {
    gate: Arc<CancelGate>,
    observations: AtomicU32,
    after: u32,
}

impl TimeSource for CancelAfterObservations {
    fn now(&self) -> Time {
        let n = self.observations.fetch_add(1, Ordering::SeqCst);
        if n >= self.after {
            self.gate.request();
        }
        Time::from_nanos(1_000 + u64::from(n))
    }
}

/// Deterministic clock whose reading crosses the budget's fixed deadline
/// on its `crossing`-th observation.
struct DeadlineCrossingClock {
    observations: AtomicU32,
    crossing: u32,
}

impl TimeSource for DeadlineCrossingClock {
    fn now(&self) -> Time {
        let n = self.observations.fetch_add(1, Ordering::SeqCst);
        if n >= self.crossing {
            Time::from_nanos(u64::from(self.crossing) * 100 + 1)
        } else {
            Time::from_nanos(u64::from(n) * 10 + 1)
        }
    }
}

struct FixedClock(u64);

impl TimeSource for FixedClock {
    fn now(&self) -> Time {
        Time::from_nanos(self.0)
    }
}

// ---------------------------------------------------------------- G0

#[test]
fn plan_counts_above_versioned_caps_or_zero_refuse() {
    let shape = |starts, steps, samples| RestorationWorkShape {
        dimensions: 2,
        constraints_total: 2,
        skipped_count: 0,
        limits: RestorationWorkLimits {
            starts,
            steps_per_start: steps,
            feasibility_samples: samples,
        },
    };
    for (starts, steps, samples) in [
        (0, 300, 400),
        (RESTORATION_MAX_STARTS + 1, 300, 400),
        (8, 0, 400),
        (8, RESTORATION_MAX_STEPS_PER_START + 1, 400),
        (8, 300, 0),
        (8, 300, RESTORATION_MAX_FEASIBILITY_SAMPLES + 1),
    ] {
        let error = RestorationWorkPlan::plan(shape(starts, steps, samples))
            .expect_err("counts outside their versioned caps refuse");
        assert!(matches!(error, ConError::BadParam { .. }), "{error:?}");
    }
    // Impossible skip cardinality refuses too.
    let bad = RestorationWorkShape {
        dimensions: 1,
        constraints_total: 2,
        skipped_count: 3,
        limits: RestorationWorkLimits::default(),
    };
    assert!(RestorationWorkPlan::plan(bad).is_err());
    // The historical schedule itself admits.
    assert!(RestorationWorkPlan::plan(shape(8, 300, 400)).is_ok());
}

#[test]
fn overflowing_dimension_arithmetic_refuses_instead_of_wrapping() {
    // `2 * u32::MAX * u32::MAX` exceeds u64; the checked product must
    // refuse rather than wrap.
    let shape = RestorationWorkShape {
        dimensions: u32::MAX,
        constraints_total: u32::MAX,
        skipped_count: 0,
        limits: RestorationWorkLimits::default(),
    };
    let error = RestorationWorkPlan::plan(shape)
        .expect_err("u32::MAX dimensions times constraints cannot state their FD step cost");
    assert!(matches!(error, ConError::BadParam { .. }), "{error:?}");
}

#[test]
fn duplicate_and_out_of_range_skip_indices_refuse_before_any_work() {
    let (problem, nodes) = host_with(&[(1.0, 0, 0.5), (1.0, 1, 0.5)], 2);
    let specs = [hard("a", nodes[0]), hard("b", nodes[1])];
    let domain = DomainBox {
        ranges: vec![(0.0, 1.0), (0.0, 1.0)],
    };
    let gate = CancelGate::new();
    with_cx_raw(&gate, infinite_budget(), None, |cx| {
        for skip in [&[0usize, 0usize][..], &[7usize][..]] {
            let outcome = elastic_solve(&problem, &specs, &domain, skip, cx)
                .expect_err("malformed skip lists refuse");
            let RestorationError::Invalid(ConError::BadParam { .. }) = outcome else {
                panic!("expected a typed input fault, got {outcome:?}");
            };
        }
    });
}

#[test]
fn plan_identity_is_order_free_and_tamper_detectable() {
    let make = || {
        RestorationWorkPlan::plan(RestorationWorkShape {
            dimensions: 3,
            constraints_total: 5,
            skipped_count: 2,
            limits: RestorationWorkLimits::default(),
        })
        .expect("shape")
    };
    let a = make();
    let b = make();
    assert_eq!(a.identity(), b.identity(), "identical shapes share identity");

    let mut tampered = a;
    tampered.total_work_units += 1;
    assert!(
        tampered.verify_consistency().is_err(),
        "an aggregate edited behind the constructor must fail re-derivation"
    );
}

#[test]
fn active_count_above_total_is_a_typed_refusal() {
    let mut plan = RestorationWorkPlan::plan(RestorationWorkShape {
        dimensions: 1,
        constraints_total: 1,
        skipped_count: 0,
        limits: RestorationWorkLimits::default(),
    })
    .expect("plan");
    plan.active_constraints = 2;

    let verification = std::panic::catch_unwind(|| plan.verify_consistency())
        .expect("an untrusted plan must not unwind while being verified");
    assert!(
        matches!(verification, Err(ConError::BadParam { .. })),
        "the invalid count relationship must be a typed refusal"
    );

    let (problem, nodes) = host_with(&[(1.0, 0, 0.5)], 1);
    let specs = [hard("a", nodes[0])];
    let domain = DomainBox {
        ranges: vec![(0.0, 1.0)],
    };
    let gate = CancelGate::new();
    with_cx_raw(&gate, infinite_budget(), None, |cx| {
        let outcome = elastic_solve_with_plan(&problem, &specs, &domain, &[], plan, cx)
            .expect_err("an impossible caller-supplied active count must refuse");
        assert!(matches!(
            outcome,
            RestorationError::Invalid(ConError::BadParam { .. })
        ));
    });
}

#[test]
fn stale_schema_plans_refuse_admission_site_checks() {
    let (problem, nodes) = host_with(&[(1.0, 0, 0.5)], 1);
    let specs = [hard("a", nodes[0])];
    let domain = DomainBox { ranges: vec![(0.0, 1.0)] };
    let mut plan = RestorationWorkPlan::plan(RestorationWorkShape {
        dimensions: 1,
        constraints_total: 1,
        skipped_count: 0,
        limits: RestorationWorkLimits::default(),
    })
    .expect("plan");
    plan.schema_version = RESTORATION_WORK_PLAN_SCHEMA_VERSION + 1;

    let gate = CancelGate::new();
    with_cx_raw(&gate, infinite_budget(), None, |cx| {
        let outcome = elastic_solve_with_plan(&problem, &specs, &domain, &[], plan, cx)
            .expect_err("a foreign-schema plan cannot bind this host");
        assert!(matches!(
            outcome,
            RestorationError::Invalid(ConError::BadParam { .. })
        ));
    });
}

// ------------------------------------------- G4 admission-time refusals

#[test]
fn zero_cost_quota_refuses_admission_with_no_consumption() {
    let (problem, nodes) = host_with(&[(1.0, 0, 2.0)], 1);
    let specs = [hard("a", nodes[0])];
    let domain = DomainBox { ranges: vec![(0.0, 1.0)] };
    let budget = Budget {
        cost_quota: Some(0),
        ..infinite_budget()
    };
    let gate = CancelGate::new();
    with_cx_raw(&gate, budget, None, |cx| {
        let outcome = elastic_solve(&problem, &specs, &domain, &[], cx)
            .expect_err("a zero cost quota cannot admit positive planned work");
        let RestorationError::Refused { refusal, receipt } = outcome else {
            panic!("expected an admission refusal, got {outcome:?}");
        };
        assert!(matches!(
            refusal,
            fs_exec::BudgetRefusal::CostPlanExceedsQuota { .. }
        ));
        assert!(
            receipt.consumption.is_none(),
            "no budget contract existed to report"
        );
        assert_eq!(receipt.work_units_charged, 0);
    });
}

#[test]
fn zero_poll_quota_stops_at_the_first_checkpoint_with_nothing_charged() {
    let (problem, nodes) = host_with(&[(1.0, 0, 2.0)], 1);
    let specs = [hard("a", nodes[0])];
    let domain = DomainBox { ranges: vec![(0.0, 1.0)] };
    let budget = Budget {
        poll_quota: 0,
        ..infinite_budget()
    };
    let gate = CancelGate::new();
    with_cx_raw(&gate, budget, None, |cx| {
        let outcome = elastic_solve(&problem, &specs, &domain, &[], cx)
            .expect_err("zero polls cannot pass even the mask checkpoint");
        let RestorationError::Refused { refusal, receipt } = outcome else {
            panic!("expected PollsExhausted, got {outcome:?}");
        };
        assert!(matches!(refusal, fs_exec::BudgetRefusal::PollsExhausted { .. }));
        let consumption = receipt.consumption.expect("admitted runs retain consumption");
        assert_eq!(consumption.cost_charged, 0);
        assert_eq!(consumption.refusal, Some(refusal));
    });
}

#[test]
fn deadlines_fail_closed_without_a_clock_and_expire_at_admission_when_passed() {
    let (problem, nodes) = host_with(&[(1.0, 0, 2.0)], 1);
    let specs = [hard("a", nodes[0])];
    let domain = DomainBox { ranges: vec![(0.0, 1.0)] };

    let unenforceable = Budget {
        deadline: Some(Time::from_nanos(10)),
        ..infinite_budget()
    };
    let gate = CancelGate::new();
    with_cx_raw(&gate, unenforceable, None, |cx| {
        let outcome = elastic_solve(&problem, &specs, &domain, &[], cx)
            .expect_err("an enforceable-deadline budget never admits without a clock");
        assert!(matches!(
            outcome,
            RestorationError::Refused {
                refusal: fs_exec::BudgetRefusal::DeadlineWithoutClock { .. },
                ..
            }
        ));
    });

    let already_past = Budget {
        deadline: Some(Time::from_nanos(50)),
        ..infinite_budget()
    };
    let clock = FixedClock(100);
    with_cx_raw(&gate, already_past, Some(&clock), |cx| {
        let outcome = elastic_solve(&problem, &specs, &domain, &[], cx)
            .expect_err("an already-expired deadline refuses admission");
        assert!(matches!(
            outcome,
            RestorationError::Refused {
                refusal: fs_exec::BudgetRefusal::DeadlineExpiredAtAdmission { .. },
                ..
            }
        ));
    });
}

// -------------------------------------------------- G4 mid-flight stops

#[test]
fn mid_run_cancellation_stops_with_receipt_and_bounded_extra_work() {
    let (problem, nodes) = host_with(&[(1.0, 0, 2.0)], 1);
    let specs = [hard("a", nodes[0])];
    let domain = DomainBox { ranges: vec![(0.0, 1.0)] };

    let plan = RestorationWorkPlan::plan(RestorationWorkShape {
        dimensions: 1,
        constraints_total: 1,
        skipped_count: 0,
        limits: RestorationWorkLimits::default(),
    })
    .expect("plan");

    let budget = Budget {
        deadline: Some(Time::from_nanos(1 << 40)), // far future: only the clock's cancel side-effect fires
        ..infinite_budget()
    };
    let gate = Arc::new(CancelGate::new());
    let clock = CancelAfterObservations {
        gate: Arc::clone(&gate),
        observations: AtomicU32::new(0),
        after: 6, // admission + mask + seed pass + a few early checkpoints
    };
    with_cx_raw(&gate, budget, Some(&clock), |cx| {
        let outcome = elastic_solve(&problem, &specs, &domain, &[], cx)
            .expect_err("the requested cancellation must stop the run");
        let RestorationError::Refused { refusal, receipt } = outcome else {
            panic!("expected a typed budget stop, got {outcome:?}");
        };
        assert!(matches!(refusal, fs_exec::BudgetRefusal::Cancelled { .. }));
        assert!(
            receipt.work_units_charged > 0,
            "some admitted tiles completed before the request landed"
        );
        assert!(
            receipt.starts_completed < plan.limits.starts,
            "the stop must land inside the descent schedule"
        );
        let consumption = receipt.consumption.expect("admitted stops keep consumption");
        assert_eq!(consumption.refusal, Some(refusal), "latched reason is stable");
        // Bounded-extra-work proof: charged units stay well under the
        // declared worst case because the drain aborts at the next tile.
        assert!(
            receipt.work_units_charged < plan.total_work_units,
            "a cancelled run may not silently complete its whole allowance"
        );
    });
}

#[test]
fn mid_run_deadline_expiry_stops_with_the_crossing_observation() {
    let (problem, nodes) = host_with(&[(1.0, 0, 2.0)], 1);
    let specs = [hard("a", nodes[0])];
    let domain = DomainBox { ranges: vec![(0.0, 1.0)] };
    let budget = Budget {
        deadline: Some(Time::from_nanos(60)),
        ..infinite_budget()
    };
    let gate = CancelGate::new();
    let clock = DeadlineCrossingClock {
        observations: AtomicU32::new(0),
        crossing: 5,
    };
    with_cx_raw(&gate, budget, Some(&clock), |cx| {
        let outcome = elastic_solve(&problem, &specs, &domain, &[], cx)
            .expect_err("crossing the deadline must stop the run");
        assert!(matches!(
            outcome,
            RestorationError::Refused {
                refusal: fs_exec::BudgetRefusal::DeadlineExpired { .. },
                ..
            }
        ));
    });
}

// ------------------------------------------------- G3 / G5 / boundaries

#[test]
fn replays_are_bit_stable_across_reports_and_receipts() {
    let (problem, nodes) = host_with(&[(1.0, 0, 2.0), (1.0, 1, 2.0)], 2);
    let specs = [hard("a", nodes[0]), hard("b", nodes[1])];
    let domain = DomainBox {
        ranges: vec![(0.0, 1.0), (0.0, 1.0)],
    };
    let gate = CancelGate::new();

    let first = with_cx_raw(&gate, infinite_budget(), None, |cx| {
        elastic_solve(&problem, &specs, &domain, &[], cx).expect("first run")
    });
    let second = with_cx_raw(&gate, infinite_budget(), None, |cx| {
        elastic_solve(&problem, &specs, &domain, &[], cx).expect("second run")
    });

    assert_eq!(first.x, second.x, "minimizer replays bitwise");
    #[allow(clippy::float_cmp)]
    {
        assert_eq!(first.total_violation, second.total_violation);
        assert_eq!(first.violations, second.violations);
    }
    assert_eq!(first.evals, second.evals);
    assert_eq!(
        first.work, second.work,
        "resource receipts are part of the replay surface"
    );
}

#[test]
fn equivalent_skip_orderings_share_one_plan_identity_and_result() {
    let (problem, nodes) = host_with(&[(1.0, 0, 0.25), (1.0, 1, 0.25)], 2);
    let specs = [hard("keep-a", nodes[0]), hard("drop-b", nodes[1])];
    let domain = DomainBox {
        ranges: vec![(0.0, 1.0), (0.0, 1.0)],
    };
    let gate = CancelGate::new();
    let forward = with_cx_raw(&gate, infinite_budget(), None, |cx| {
        elastic_solve(&problem, &specs, &domain, &[1], cx).expect("forward")
    });
    // Same single-element set, different orderings of a multi-skip list
    // exercise the canonicalization path directly.
    let reversed_pair = with_cx_raw(&gate, infinite_budget(), None, |cx| {
        elastic_solve(&problem, &specs, &domain, &[], cx).expect("no skips")
    });
    assert_ne!(
        forward.work.plan_identity, reversed_pair.work.plan_identity,
        "different masks are different workloads"
    );
    // Order equivalence: a two-element skip list in both orders.
    let order_one = with_cx_raw(&gate, infinite_budget(), None, |cx| {
        elastic_solve(&problem, &specs, &domain, &[0, 1], cx).expect("order one")
    });
    let order_two = with_cx_raw(&gate, infinite_budget(), None, |cx| {
        elastic_solve(&problem, &specs, &domain, &[1, 0], cx).expect("order two")
    });
    assert_eq!(order_one.work.plan_identity, order_two.work.plan_identity);
    assert_eq!(order_one.x, order_two.x);
    assert_eq!(order_one.work, order_two.work);
}

#[test]
fn leaseless_contexts_record_the_no_lease_no_claim_boundary() {
    let (problem, nodes) = host_with(&[(1.0, 0, 0.25)], 1);
    let specs = [hard("a", nodes[0])];
    let domain = DomainBox { ranges: vec![(0.0, 1.0)] };
    let gate = CancelGate::new();
    with_cx_raw(&gate, infinite_budget(), None, |cx| {
        let report = elastic_solve(&problem, &specs, &domain, &[], cx).expect("solve");
        assert_eq!(
            report.work.memory,
            fs_constraint::RestorationMemoryAuthority::NoLeaseNoClaim,
            "hand-built contexts carry no lease; the receipt must say so"
        );
        let consumption = report.work.consumption.expect("successful runs keep it");
        assert!(consumption.refusal.is_none(), "success implies no latched refusal");
        assert!(report.work.starts_completed == RESTORATION_MAX_STARTS);
    });
}

#[test]
fn diagnosis_binds_one_shared_receipt_across_every_phase() {
    // Pairwise infeasible: g_a <= 0 needs x >= 0.75, g_b <= 0 needs x <= 0.25.
    let (problem, nodes) = host_with(&[(-1.0, 0, -0.75), (1.0, 0, 0.25)], 1);
    let specs = [hard("upper", nodes[0]), hard("lower", nodes[1])];
    let domain = DomainBox { ranges: vec![(0.0, 1.0)] };
    let gate = CancelGate::new();
    with_cx_raw(&gate, infinite_budget(), None, |cx| {
        let diagnosis =
            diagnose_infeasibility(&problem, &specs, &domain, cx).expect("diagnosis");
        assert!(!diagnosis.feasible);
        assert_eq!(diagnosis.core.len(), 2, "both members are necessary here");
        assert!(!diagnosis.repairs.is_empty());
        let consumption = diagnosis.work.consumption.expect("admitted run");
        assert!(consumption.refusal.is_none());
        assert!(diagnosis.work.work_units_charged > 0);
        // The shared accountant covered base solve + filter + repairs.
        assert!(
            diagnosis.work.work_units_charged > diagnosis.elastic.evals,
            "filter and repair phases charge beyond the base solve's evals"
        );
        assert_eq!(diagnosis.elastic.work, diagnosis.work);

        // Replay: the whole diagnosis (numbers AND receipts) is stable.
        let again = diagnose_infeasibility(&problem, &specs, &domain, cx).expect("replay");
        assert_eq!(diagnosis.core, again.core);
        assert_eq!(diagnosis.elastic.x, again.elastic.x);
        assert_eq!(diagnosis.work, again.work);
    });
}
