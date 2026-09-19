//! End-to-end solves consume real reverse-mode IR derivatives, not test callbacks.
use fs_ascent::{ReverseSqpError, ReverseSqpStudy};
use fs_ascent::sqp::SqpStop;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, EvalLimit, Manifold, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use std::num::NonZeroU64;

const LIMITS: ReverseLimits = ReverseLimits { max_nodes: 100, max_scalar_slots: 1000 };
fn limit(n: u64) -> EvalLimit { EvalLimit::Limited(NonZeroU64::new(n).unwrap()) }

fn quadratic(budget: EvalLimit, weighted: bool) -> Problem {
    let mut b = ProblemBuilder::new();
    // Separate variable blocks exercise packing as well as constraint ordering.
    let vx = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let vy = b.var("y", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let rx = b.var_ref(vx).unwrap();
    let ry = b.var_ref(vy).unwrap();
    let x = b.component(rx, 0).unwrap();
    let y = b.component(ry, 0).unwrap();
    let one = b.konst(1.0, Dims::NONE).unwrap();
    let two = b.konst(2.0, Dims::NONE).unwrap();
    let three = b.konst(3.0, Dims::NONE).unwrap();
    let cap = b.konst(1.2, Dims::NONE).unwrap();
    let dx = b.sub(x, two).unwrap();
    let dy = b.sub(y, one).unwrap();
    let sx = b.powi(dx, 2).unwrap();
    let sy = b.powi(dy, 2).unwrap();
    b.objective(sx, Sense::Minimize, if weighted { 3.0 } else { 1.0 }).unwrap();
    if weighted {
        let negative = b.neg(sy).unwrap();
        b.objective(negative, Sense::Maximize, 2.0).unwrap();
    } else { b.objective(sy, Sense::Minimize, 1.0).unwrap(); }
    let sum = b.add(x, y).unwrap();
    let equality = b.sub(sum, two).unwrap();
    let upper_x = b.sub(x, cap).unwrap();
    let upper_y = b.sub(y, three).unwrap();
    b.constraint(upper_x, ConstraintKind::LeZero, "active-upper-x").unwrap();
    b.constraint(equality, ConstraintKind::EqZero, "sum").unwrap();
    b.constraint(upper_y, ConstraintKind::LeZero, "inactive-upper-y").unwrap();
    b.set_eval_limit(budget);
    b.finish()
}

fn scalar_problem(logarithm: bool, budget: EvalLimit) -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let x = b.component(r, 0).unwrap();
    let square = b.powi(x, 2).unwrap();
    let objective = if logarithm {
        let ln = b.ln(x).unwrap(); b.sub(square, ln).unwrap()
    } else { square };
    b.objective(objective, Sense::Minimize, 1.0).unwrap();
    let ceiling = b.konst(20.0, Dims::NONE).unwrap();
    let c = b.sub(x, ceiling).unwrap();
    b.constraint(c, ConstraintKind::LeZero, "loose-upper-bound").unwrap();
    b.set_eval_limit(budget);
    b.finish()
}

#[test]
fn solves_mixed_constraints_and_returns_ir_ordered_multipliers() {
    let problem = quadratic(limit(20), false);
    let oracle = ReverseProblem::new(&problem, LIMITS).unwrap();
    let mut study = ReverseSqpStudy::new(&oracle, &[0.0, 0.0], 16, None).unwrap();
    let report = study.run(1e-8, 20, None).unwrap();
    assert_eq!(report.stop, SqpStop::Converged);
    assert!((report.solution.x[0] - 1.2).abs() < 1e-8);
    assert!((report.solution.x[1] - 0.8).abs() < 1e-8);
    for (got, want) in report.constraint_multipliers.iter().zip([1.2, 0.4, 0.0]) {
        assert!((got - want).abs() < 1e-8);
    }
    let tape = oracle.evaluate(&[vec![report.solution.x[0]], vec![report.solution.x[1]]]).unwrap();
    let gradient = tape.lagrangian_gradient(&report.constraint_multipliers).unwrap();
    assert!(gradient.iter().flatten().all(|v| v.abs() < 1e-8));
    assert!(study.optimizer().sample().ci[1] < -2.0, "inactive residual was clipped");
    assert_eq!(study.optimizer().evaluations(), 2);
}

#[test]
fn respects_sense_weights_and_multiple_variable_blocks() {
    let problem = quadratic(limit(50), true);
    let oracle = ReverseProblem::new(&problem, LIMITS).unwrap();
    let mut study = ReverseSqpStudy::new(&oracle, &[0.0, 0.0], 16, None).unwrap();
    let report = study.run(1e-8, 30, None).unwrap();
    assert!(report.solution.converged, "{report:?}");
    assert!((report.solution.f - 2.0).abs() < 1e-8);
    for (got, want) in report.constraint_multipliers.iter().zip([4.0, 0.8, 0.0]) {
        assert!((got - want).abs() < 1e-8);
    }
}

#[test]
fn problem_budget_cannot_be_raised_by_continuation() {
    let problem = quadratic(limit(1), false);
    let oracle = ReverseProblem::new(&problem, LIMITS).unwrap();
    let mut study = ReverseSqpStudy::new(&oracle, &[0.0, 0.0], 16, None).unwrap();
    for _ in 0..3 {
        let report = study.run_with_budget(1e-8, 30, 10_000, None).unwrap();
        assert_eq!(report.stop, SqpStop::EvaluationLimit);
        assert_eq!(study.optimizer().evaluations(), 1);
        assert_eq!(study.optimizer().point(), &[0.0, 0.0]);
        assert!(!report.solution.converged);
    }
}

#[test]
fn stricter_trial_budget_preserves_state_and_counts_repeated_search_work() {
    let problem = scalar_problem(false, limit(10));
    let oracle = ReverseProblem::new(&problem, LIMITS).unwrap();
    let mut study = ReverseSqpStudy::new(&oracle, &[2.0], 8, None).unwrap();
    let report = study.run_with_budget(1e-10, 20, 2, None).unwrap();
    assert_eq!(report.stop, SqpStop::EvaluationLimit);
    assert_eq!(study.optimizer().point(), &[2.0]);
    let report = study.run(1e-10, 20, None).unwrap();
    assert!(report.solution.converged);
    assert_eq!(study.optimizer().point(), &[0.0]);
    assert_eq!(study.optimizer().evaluations(), 4);
}

#[test]
fn logarithmic_trial_errors_backtrack_but_invalid_initial_values_refuse() {
    let problem = scalar_problem(true, limit(100));
    let oracle = ReverseProblem::new(&problem, LIMITS).unwrap();
    assert!(ReverseSqpStudy::new(&oracle, &[-1.0], 8, None).is_err());
    let mut study = ReverseSqpStudy::new(&oracle, &[4.0], 8, None).unwrap();
    let report = study.run(1e-7, 60, None).unwrap();
    assert!(report.solution.converged, "{report:?}");
    assert!((study.optimizer().point()[0] - 0.5_f64.sqrt()).abs() < 1e-6);
    assert!(study.domain_rejections() > 0);
    assert!(study.last_rejection().is_some());
}

#[test]
fn cloned_checkpoints_continue_without_reinitialization() {
    let problem = scalar_problem(true, limit(100));
    let oracle = ReverseProblem::new(&problem, LIMITS).unwrap();
    let initial = ReverseSqpStudy::new(&oracle, &[4.0], 8, None).unwrap();
    let mut straight = initial.clone();
    straight.run(1e-12, 4, None).unwrap();
    assert!(straight.optimizer().iterations() > 1);
    for split in 0..=4 {
        let mut resumed = initial.clone();
        resumed.run(1e-12, split, None).unwrap();
        resumed.run(1e-12, 4 - split, None).unwrap();
        assert_eq!(straight.optimizer().point(), resumed.optimizer().point());
        assert_eq!(straight.optimizer().history(), resumed.optimizer().history());
        assert_eq!(straight.optimizer().evaluations(), resumed.optimizer().evaluations());
        assert_eq!(straight.domain_rejections(), resumed.domain_rejections());
        assert_eq!(straight.last_rejection(), resumed.last_rejection());
    }
}

#[test]
fn structural_dimension_and_manifold_refusals_are_not_objective_trials() {
    let problem = quadratic(limit(20), false);
    let oracle = ReverseProblem::new(&problem, LIMITS).unwrap();
    assert!(matches!(ReverseSqpStudy::new(&oracle, &[0.0, 0.0], 4, None), Err(ReverseSqpError::DimensionCap)));
    assert!(matches!(ReverseSqpStudy::new(&oracle, &[0.0], 16, None), Err(ReverseSqpError::PackedPointLength { .. })));
    let mut b = ProblemBuilder::new();
    let v = b.var("sphere", Manifold::Sphere { ambient: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let x = b.component(r, 0).unwrap();
    b.objective(x, Sense::Minimize, 1.0).unwrap();
    let sphere = b.finish();
    let oracle = ReverseProblem::new(&sphere, LIMITS).unwrap();
    assert!(matches!(ReverseSqpStudy::new(&oracle, &[0.0, 0.0], 16, None), Err(ReverseSqpError::NonEuclidean { variable: 0 })));
}

#[test]
fn real_cx_cancellation_retains_the_checkpoint_and_resumes() {
    let problem = quadratic(limit(20), false);
    let oracle = ReverseProblem::new(&problem, LIMITS).unwrap();
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 0, kernel_id: 1, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let mut study = ReverseSqpStudy::new(&oracle, &[0.0, 0.0], 16, Some(&cx)).unwrap();
        gate.request();
        assert_eq!(study.run(1e-8, 20, Some(&cx)).unwrap_err(), ReverseSqpError::Cancelled);
        assert_eq!(study.optimizer().evaluations(), 1);
        assert_eq!(study.optimizer().point(), &[0.0, 0.0]);
        assert!(study.run(1e-8, 20, None).unwrap().solution.converged);
    });
}
