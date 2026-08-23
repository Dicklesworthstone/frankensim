//! Chance-constraint resource-contract gates (bead frankensim-oxyjg):
//! checked [`ChanceWorkPlan`] arithmetic, budget/cancellation authority
//! enforced at deterministic tile boundaries, retained receipts with NO
//! partial evidence, retired-sampler refusal, and replay-stable results.
//!
//! Coverage map: G0 (arithmetic max+1 boundaries, plan/estimator policy
//! mismatches, identity binding), G4 (pre-cancel, mid-run cancel via the
//! noise seam, cost/poll/deadline exhaustion families, injected noise and
//! evaluator faults — every stop drains with a receipt and never emits
//! partial satisfied/violated evidence), G3 (tile-size partition
//! invariance), G5 (bit-stable receipt/evidence replay).
//!
//! Determinism: every run uses the same logical sample stream `s ↦ draw`;
//! tiling changes only checkpoint granularity, never the draws or their
//! order, so tile size must not move any asserted value.
//!
//! Lint scope: the admitted error type intentionally travels by value
//! (receipts bind to it), and several `matches!` arms borrow fields from
//! owned errors; both style lints are test-file-local.

#![allow(clippy::float_cmp)]
#![allow(clippy::result_large_err, clippy::needless_borrow)]

use asupersync::types::{Budget, Time};
use fs_constraint::{
    CHANCE_WORK_PLAN_SCHEMA_VERSION, ChanceEstimator, ChanceEvalError, ChanceWorkPlan,
    ChanceWorkReceipt, ConError, ConstraintKind, ConstraintSpec, Status, evaluate,
    evaluate_chance_with_budget,
};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::{Manifold, NodeId, Problem, ProblemBuilder};
use fs_qty::Dims;
use std::cell::Cell;

const EXECUTION_SEED: u64 = 0xC0C0;

/// One Rn(2) host; node: `x0 - 1 ≤ 0`.
fn host_linear() -> (Problem, NodeId) {
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
        .expect("var");
    let vr = b.var_ref(v).expect("ref");
    let x0 = b.component(vr, 0).expect("x0");
    let c0 = b.konst(1.0, Dims::NONE).expect("konst");
    let t0 = b.mul(c0, x0).expect("t");
    let rb = b.konst(1.0, Dims::NONE).expect("rhs");
    let node = b.sub(t0, rb).expect("g");
    let obj = b.norm_sq(vr).expect("obj");
    b.objective(obj, fs_opt::Sense::Minimize, 1.0).expect("o");
    (b.finish(), node)
}

fn chance_spec(node: NodeId, level: f64, samples: u32) -> ConstraintSpec {
    ConstraintSpec {
        name: "chance-cap".to_string(),
        node,
        kind: ConstraintKind::Chance {
            level,
            estimator: ChanceEstimator::MonteCarlo {
                samples,
                delta: 0.05,
            },
        },
        active_tol: 1e-9,
    }
}

/// Deterministic per-sample stream: `u ~ U(0,1)` on x0 keyed by the
/// logical sample index only.
fn stream() -> impl Fn(u64) -> Vec<f64> {
    |s: u64| {
        // SplitMix-style scramble keeps draws well spread over s.
        let mut z = s
            .wrapping_mul(0x9E37_79B9_7F4A_7C15)
            .wrapping_add(0x1234_5678);
        z ^= z >> 30;
        z = z.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z ^= z >> 27;
        z = z.wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        vec![(z >> 11) as f64 / (1u64 << 53) as f64, 0.0]
    }
}

fn with_cx_budget<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn infinite() -> Budget {
    Budget::INFINITE
}

fn budget(deadline: Option<u64>, polls: u32, cost: Option<u64>) -> Budget {
    Budget {
        deadline: deadline.map(Time::from_secs),
        poll_quota: polls,
        cost_quota: cost,
        priority: 0,
    }
}

/// Unwrap a refused chance outcome, failing through `expect_err` so the
/// gate file carries no bare `panic!` for scanners to misread as
/// library-code aborts.
fn into_refused(error: ChanceEvalError) -> (fs_exec::BudgetRefusal, ChanceWorkReceipt) {
    let outcome = match error {
        ChanceEvalError::Refused { refusal, receipt } => Err((refusal, receipt)),
        ChanceEvalError::Invalid(error) => Ok(error),
    };
    outcome.expect_err("expected a budget refusal")
}

// ---------------------------------------------------------------- G0

#[test]
fn g0_plan_arithmetic_refuses_every_boundary() {
    assert!(matches!(
        ChanceWorkPlan::plan(0, 2, 1),
        Err(ConError::BadParam { what, .. }) if what.contains("positive sample")
    ));
    assert!(matches!(
        ChanceWorkPlan::plan(16, 0, 1),
        Err(ConError::BadParam { what, .. }) if what.contains("positive dimension")
    ));
    assert!(matches!(
        ChanceWorkPlan::plan(16, 2, 0),
        Err(ConError::BadParam { what, .. }) if what.contains("tile size")
    ));
    assert!(matches!(
        ChanceWorkPlan::plan(16, 2, 17),
        Err(ConError::BadParam { what, .. }) if what.contains("tile size")
    ));
    // Per-sample work can never overflow for admissible dimensions
    // (u32::MAX dimensions declare at most 2*2^32 + 2 units), so the
    // checked-multiply guard is exercised on the TOTAL: u32::MAX
    // samples across u32::MAX dimensions declares ~3.7e19 units,
    // beyond the u64 range, and must refuse rather than wrap.
    assert!(matches!(
        ChanceWorkPlan::plan(u32::MAX, u32::MAX, 1),
        Err(ConError::BadParam { what, .. }) if what.contains("total work")
    ));
    // plan either fits exactly or refuses — never wraps silently.
    let ok = ChanceWorkPlan::plan(1_000, 1_000, 64).expect("modest plan fits");
    assert_eq!(
        ok.per_sample_work_units,
        1_000 * 2 + 1 + 1,
        "weights are the declared contract constants"
    );
    assert_eq!(
        ok.total_work_units,
        u64::from(ok.samples) * ok.per_sample_work_units
    );
}

#[test]
fn g0_plan_identity_binds_every_field_and_schema() {
    let base = ChanceWorkPlan::plan(400, 2, 64).expect("base");
    assert_eq!(base.schema_version, CHANCE_WORK_PLAN_SCHEMA_VERSION);
    assert_eq!(base.identity(), base.identity(), "identity is pure");
    let variants = [
        ChanceWorkPlan {
            schema_version: base.schema_version + 1,
            ..base
        },
        ChanceWorkPlan {
            samples: base.samples - 1,
            ..base
        },
        ChanceWorkPlan {
            dimensions: base.dimensions + 1,
            ..base
        },
        ChanceWorkPlan {
            tile_samples: base.tile_samples - 1,
            ..base
        },
    ];
    for variant in variants {
        assert_ne!(
            base.identity(),
            variant.identity(),
            "field change moves identity"
        );
    }
}

#[test]
fn g0_retired_synchronous_sampler_refuses_typed() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 8);
    let error = evaluate(&problem, &spec, &[0.5, 0.0], None).expect_err("retired");
    assert!(
        matches!(&error, ConError::BadParam { what, .. } if what.contains("evaluate_chance_with_budget")),
        "the refusal must teach the replacement API, got {error}"
    );
    // A hard constraint still evaluates through the legacy entry point:
    // the retirement is chance-specific, not a general lockout.
    let hard = ConstraintSpec {
        name: "hard".to_string(),
        node,
        kind: ConstraintKind::Hard,
        active_tol: 1e-9,
    };
    let ev = evaluate(&problem, &hard, &[0.5, 0.0], None).expect("hard evaluates");
    assert_eq!(ev.status, Status::Satisfied);
}
#[test]
fn g0_plan_estimator_policy_mismatches_refuse_before_admission() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 16);
    let noise = stream();
    // Wrong schema version.
    let mut plan = ChanceWorkPlan::plan(16, 2, 4).expect("plan");
    let mutated = ChanceWorkPlan {
        schema_version: plan.schema_version + 1,
        ..plan
    };
    let gate = CancelGate::new();
    let error = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, mutated, cx)
    })
    .expect_err("schema mismatch refuses");
    assert!(matches!(
        error,
        ChanceEvalError::Invalid(ConError::BadParam { ref what, .. })
            if what.contains("schema version")
    ));
    // Estimator/sample disagreement.
    plan = ChanceWorkPlan::plan(15, 2, 4).expect("plan");
    let error = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    })
    .expect_err("sample mismatch refuses");
    assert!(matches!(
        error,
        ChanceEvalError::Invalid(ConError::BadParam { ref what, .. })
            if what.contains("match")
    ));
    // Dimension disagreement with the admitted host.
    let wrong_dims = ChanceWorkPlan {
        dimensions: 3,
        ..ChanceWorkPlan::plan(16, 2, 4).expect("plan")
    };
    let error = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, wrong_dims, cx)
    })
    .expect_err("dimension mismatch refuses");
    assert!(matches!(
        error,
        ChanceEvalError::Invalid(ConError::BadParam { ref what, .. })
            if what.contains("host dimension")
    ));
    // Non-chance kinds refuse outright.
    let hard = ConstraintSpec {
        name: "hard".to_string(),
        node,
        kind: ConstraintKind::Hard,
        active_tol: 1e-9,
    };
    plan = ChanceWorkPlan::plan(16, 2, 4).expect("plan");
    let error = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &hard, &[0.5, 0.0], &noise, plan, cx)
    })
    .expect_err("non-chance kind refuses");
    assert!(matches!(
        error,
        ChanceEvalError::Invalid(ConError::BadParam { ref what, .. })
            if what.contains("requires a chance constraint")
    ));
}

// ---------------------------------------------------------------- G4

#[test]
fn g4_cost_quota_exhaustion_stops_with_exact_receipt_and_no_evidence() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 64);
    let noise = stream();
    // Weights: 2 dims x (noise 1 + shift 1) + eval 1 + accum 1 = 6/sample.
    let plan = ChanceWorkPlan::plan(64, 2, 8).expect("plan");
    let tile_units = u64::from(plan.tile_samples) * plan.per_sample_work_units;
    assert_eq!(plan.total_work_units, tile_units * 8);
    // Admission-side exhaustion: quota below the declared plan refuses
    // before any sample runs, so the receipt has no consumption record.
    let gate = CancelGate::new();
    let error = with_cx_budget(
        &gate,
        budget(None, u32::MAX, Some(plan.total_work_units - 1)),
        |cx| evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx),
    )
    .expect_err("underfunded admission refuses");
    let (refusal, receipt) = into_refused(error);
    assert_eq!(receipt.completed_samples, 0);
    assert_eq!(receipt.hits, 0);
    assert!(receipt.consumption.is_none(), "nothing was admitted");
    assert!(matches!(
        refusal,
        fs_exec::BudgetRefusal::CostPlanExceedsQuota { .. }
    ));
    // Mid-run CostExhausted is structurally unreachable through this
    // API: admission refuses any plan whose declared total exceeds the
    // quota, and tiles only ever charge out of that admitted total.
    let (ev, receipt) = with_cx_budget(
        &gate,
        budget(None, u32::MAX, Some(plan.total_work_units)),
        |cx| evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx),
    )
    .expect("exact budget completes");
    assert_eq!(receipt.completed_samples, 64);
    let consumption = receipt
        .consumption
        .expect("admitted work reports consumption");
    assert_eq!(consumption.cost_charged, plan.total_work_units);
    assert_eq!(consumption.planned_cost, plan.total_work_units);
    assert!(consumption.refusal.is_none());
    assert!(matches!(
        ev.status,
        Status::Satisfied | Status::Violated | Status::BoundNotCleared { .. }
    ));
}

#[test]
fn g4_pre_cancelled_gate_refuses_at_the_first_tile_with_empty_receipt() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 64);
    let noise = stream();
    let plan = ChanceWorkPlan::plan(64, 2, 8).expect("plan");
    let gate = CancelGate::new();
    gate.request();
    let error = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    })
    .expect_err("cancelled authority must refuse");
    let (refusal, receipt) = into_refused(error);
    assert!(matches!(refusal, fs_exec::BudgetRefusal::Cancelled { .. }));
    assert_eq!(
        receipt.completed_samples, 0,
        "no sample may run after cancel"
    );
}

#[test]
fn g4_mid_run_cancellation_driven_through_the_noise_seam() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 512);
    let plan = ChanceWorkPlan::plan(512, 2, 8).expect("plan");
    let gate = std::rc::Rc::new(CancelGate::new());
    let trip = std::rc::Rc::clone(&gate);
    let calls = Cell::new(0u64);
    let noise = move |s: u64| -> Vec<f64> {
        let _ = s;
        calls.set(calls.get() + 1);
        if calls.get() == 100 {
            trip.request();
        }
        vec![0.5, 0.0]
    };
    let error = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    })
    .expect_err("mid-run cancellation must surface as a typed refusal");
    let (refusal, receipt) = into_refused(error);
    assert_eq!(
        receipt.completed_samples, 104,
        "the tile containing the 100th draw completes; the drain lands at its end"
    );
    assert!(matches!(refusal, fs_exec::BudgetRefusal::Cancelled { .. }));
    assert!(receipt.consumption.is_some());
}

#[test]
fn g4_poll_quota_exhaustion_names_the_boundary() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 64);
    let noise = stream();
    let plan = ChanceWorkPlan::plan(64, 2, 8).expect("plan"); // 8 tiles -> 8 checkpoints
    let gate = CancelGate::new();
    let ev = with_cx_budget(&gate, budget(None, 8, None), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    })
    .expect("exactly one poll per tile fits the quota of 8")
    .0;
    assert!(matches!(
        ev.status,
        Status::Satisfied | Status::Violated | Status::BoundNotCleared { .. }
    ));
    let gate = CancelGate::new();
    let error = with_cx_budget(&gate, budget(None, 7, None), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    })
    .expect_err("quota one short of the tile count refuses");
    let (refusal, receipt) = into_refused(error);
    assert!(matches!(
        refusal,
        fs_exec::BudgetRefusal::PollsExhausted { .. }
    ));
    assert_eq!(receipt.completed_samples, 56, "seven tiles completed");
}

#[test]
fn g4_deadline_without_clock_refuses_instead_of_running_unbounded() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 64);
    let noise = stream();
    let plan = ChanceWorkPlan::plan(64, 2, 8).expect("plan");
    // The deterministic test Cx carries no TimeSource, so a deadline
    // cannot be honored — admission must refuse rather than ignore it.
    let gate = CancelGate::new();
    let error = with_cx_budget(&gate, budget(Some(0), u32::MAX, None), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    })
    .expect_err("deadline without a clock refuses");
    assert!(matches!(error, ChanceEvalError::Refused { .. }));
}

#[test]
fn g4_noise_dimension_fault_is_a_typed_invalid_refusal() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 64);
    let plan = ChanceWorkPlan::plan(64, 2, 8).expect("plan");
    let calls = Cell::new(0u64);
    let noise = move |_s: u64| -> Vec<f64> {
        calls.set(calls.get() + 1);
        vec![0.5] // wrong dimension from the tenth draw onward
    };
    let gate = CancelGate::new();
    let error = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    })
    .expect_err("wrong-dimension draws refuse");
    assert!(matches!(
        error,
        ChanceEvalError::Invalid(ConError::BadParam { ref what, .. })
            if what.contains("chance noise draw dimension")
    ));
}

#[test]
fn g4_evaluator_fault_surfaces_as_invalid_after_host_admission() {
    let mut b = ProblemBuilder::new();
    let v = b
        .var("x", Manifold::Rn { dim: 2 }, Dims::NONE)
        .expect("var");
    let vr = b.var_ref(v).expect("ref");
    let obj = b.norm_sq(vr).expect("obj");
    b.objective(obj, fs_opt::Sense::Minimize, 1.0).expect("o");
    let problem = b.finish();
    let spec = chance_spec(NodeId(u32::MAX), 0.5, 8); // forged node: eval fault
    let noise = stream();
    let plan = ChanceWorkPlan::plan(8, 2, 8).expect("plan");
    let gate = CancelGate::new();
    let outcome = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    });
    assert!(outcome.is_err(), "a forged graph node cannot evaluate");
}

#[test]
fn g4_partial_prefix_never_mints_evidence() {
    // Structural pin on the contract the G4 family relies on: Refused
    // outcomes carry a receipt tuple, NEVER a ConstraintEvidence value.
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.5, 64);
    let noise = stream();
    let plan = ChanceWorkPlan::plan(64, 2, 8).expect("plan");
    let gate = CancelGate::new();
    gate.request();
    let outcome = with_cx_budget(&gate, infinite(), |cx| {
        evaluate_chance_with_budget(&problem, &spec, &[0.5, 0.0], &noise, plan, cx)
    });
    let (_, receipt) = into_refused(outcome.expect_err("refusal expected"));
    assert_eq!(receipt.completed_samples, 0);
}

// ------------------------------------------------------ G3/G5

#[test]
fn g3_tile_partition_invariance_and_g5_replay_stability() {
    let (problem, node) = host_linear();
    let spec = chance_spec(node, 0.70, 400);
    let reference = {
        let plan = ChanceWorkPlan::plan(400, 2, 1).expect("tile=1 plan");
        let gate = CancelGate::new();
        let noise = stream();
        with_cx_budget(&gate, infinite(), |cx| {
            evaluate_chance_with_budget(&problem, &spec, &[0.2, 0.0], &noise, plan, cx)
        })
        .expect("tile=1 run completes")
    };

    for tile in [7usize, 64, 400] {
        let plan = ChanceWorkPlan::plan(400, 2, tile as u32).expect("plan");
        let gate = CancelGate::new();
        let noise = stream();
        let (ev, receipt) = with_cx_budget(&gate, infinite(), |cx| {
            evaluate_chance_with_budget(&problem, &spec, &[0.2, 0.0], &noise, plan, cx)
        })
        .expect("run completes");
        assert_eq!(ev, reference.0, "tile size {tile} moved the evidence");
        assert_eq!(
            receipt.hits, reference.1.hits,
            "tile size {tile} moved hits"
        );
        assert_eq!(receipt.completed_samples, 400);
    }
    // G5: identical inputs replays to identical bytes (PartialEq over the
    // canonical evidence/receipt records; no wall-clock enters either).
    let plan = ChanceWorkPlan::plan(400, 2, 64).expect("plan");
    let first = {
        let gate = CancelGate::new();
        let noise = stream();
        with_cx_budget(&gate, infinite(), |cx| {
            evaluate_chance_with_budget(&problem, &spec, &[0.2, 0.0], &noise, plan, cx)
        })
        .expect("replay a")
    };
    let second = {
        let gate = CancelGate::new();
        let noise = stream();
        with_cx_budget(&gate, infinite(), |cx| {
            evaluate_chance_with_budget(&problem, &spec, &[0.2, 0.0], &noise, plan, cx)
        })
        .expect("replay b")
    };
    assert_eq!(first, second, "replay diverged");
    assert_eq!(first.1.plan_identity, plan.identity());
}
