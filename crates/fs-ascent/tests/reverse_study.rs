//! End-to-end live-IR optimization: no finite-difference callbacks in the solve.
use fs_ascent::lbfgs::LbfgsError;
use fs_ascent::{ReverseStudy, ReverseStudyError, StopReason, StopRule};
use fs_opt::reverse::{ReverseError, ReverseLimits};
use fs_opt::{ConstraintKind, EvalLimit, Manifold, OptError, Problem, ProblemBuilder, ReverseProblem, ReverseProblemError, Sense};
use fs_qty::Dims;
use std::num::NonZeroU64;

fn limits() -> ReverseLimits {
    ReverseLimits { max_nodes: 10_000, max_scalar_slots: 100_000 }
}

fn limited(n: u64) -> EvalLimit {
    EvalLimit::Limited(NonZeroU64::new(n).unwrap())
}

fn quadratic(dim: u32, limit: EvalLimit) -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let f = b.norm_sq(r).unwrap();
    b.objective(f, Sense::Minimize, 1.0).unwrap();
    b.set_eval_limit(limit);
    b.finish()
}

fn rosenbrock() -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("design", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let x = b.component(r, 0).unwrap();
    let y = b.component(r, 1).unwrap();
    let one = b.konst(1.0, Dims::NONE).unwrap();
    let hundred = b.konst(100.0, Dims::NONE).unwrap();
    let xx = b.powi(x, 2).unwrap();
    let dy = b.sub(y, xx).unwrap();
    let dy2 = b.powi(dy, 2).unwrap();
    let scaled = b.mul(hundred, dy2).unwrap();
    let dx = b.sub(one, x).unwrap();
    let dx2 = b.powi(dx, 2).unwrap();
    let f = b.add(scaled, dx2).unwrap();
    b.objective(f, Sense::Minimize, 1.0).unwrap();
    b.set_eval_limit(limited(1000));
    b.finish()
}

#[test]
fn reverse_ir_reaches_rosenbrock_optimum_without_numerical_gradients() {
    let p = rosenbrock();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseStudy::new(&oracle, &[-1.2, 1.0], 7, None).unwrap();
    let report = study.run(&StopRule::GradNorm(1e-7), 500, None).unwrap();
    assert_eq!(report.reason, StopReason::GradNorm);
    assert!(report.grad_norm <= 1e-7);
    assert!(study.optimizer().x.iter().all(|x| (x - 1.0).abs() < 1e-5));
    assert!(report.evals < 1000);
    assert_eq!(study.optimizer().history.len(), report.iters + 1);
}

#[test]
fn split_reverse_studies_reuse_initial_value_and_curvature() {
    let p = rosenbrock();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let initial = ReverseStudy::new(&oracle, &[-1.2, 1.0], 7, None).unwrap();
    let mut full = initial.clone();
    let rule = StopRule::GradNorm(0.0);
    full.run(&rule, 12, None).unwrap();
    for split in 0..=12 {
        let mut resumed = initial.clone();
        resumed.run(&rule, split, None).unwrap();
        resumed.run(&rule, 12 - split, None).unwrap();
        assert_eq!(format!("{:?}", full.optimizer()), format!("{:?}", resumed.optimizer()));
        assert_eq!(full.rejected_trials(), resumed.rejected_trials());
        assert_eq!(full.last_rejection(), resumed.last_rejection());
    }
}

#[test]
fn high_dimensional_gradient_does_not_spend_one_evaluation_per_coordinate() {
    let p = quadratic(256, limited(3));
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseStudy::new(&oracle, &[1.0; 256], 7, None).unwrap();
    let report = study.run(&StopRule::GradNorm(0.0), 100, None).unwrap();
    assert_eq!(report.reason, StopReason::Budget);
    assert_eq!(report.evals, 3);
    assert_eq!(report.iters, 1);
    assert_eq!(study.optimizer().x, vec![0.0; 256]);
    assert_eq!(report.f, 0.0);
}

#[test]
fn problem_budget_is_enforced_even_mid_line_search_and_on_resume() {
    for cap in [1, 2] {
        let p = quadratic(1, limited(cap));
        let oracle = ReverseProblem::new(&p, limits()).unwrap();
        let mut study = ReverseStudy::new(&oracle, &[1.0], 7, None).unwrap();
        let report = study.run(&StopRule::GradNorm(0.0), 100, None).unwrap();
        assert_eq!(report.reason, StopReason::Budget);
        assert_eq!(report.evals as u64, cap);
        assert_eq!(study.optimizer().x, vec![1.0]);
        assert_eq!(report.iters, 0);
        let checkpoint = format!("{:?}", study.optimizer());
        study.run(&StopRule::GradNorm(0.0), 100, None).unwrap();
        assert_eq!(format!("{:?}", study.optimizer()), checkpoint);
    }
}

fn logarithmic_objective() -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("positive-x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let x = b.component(r, 0).unwrap();
    let square = b.powi(x, 2).unwrap();
    let hundred = b.konst(100.0, Dims::NONE).unwrap();
    let penalty = b.mul(hundred, square).unwrap();
    let log = b.ln(x).unwrap();
    let f = b.sub(penalty, log).unwrap();
    b.objective(f, Sense::Minimize, 1.0).unwrap();
    b.set_eval_limit(limited(200));
    b.finish()
}

#[test]
fn invalid_logarithmic_trials_backtrack_without_poisoning_accepted_state() {
    let p = logarithmic_objective();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseStudy::new(&oracle, &[1.0], 7, None).unwrap();
    let report = study.run(&StopRule::GradNorm(1e-7), 100, None).unwrap();
    assert_eq!(report.reason, StopReason::GradNorm);
    assert!((study.optimizer().x[0] - 0.005f64.sqrt()).abs() < 1e-8);
    assert!(study.rejected_trials() > 0);
    assert!(study.last_rejection().is_some());
    assert!(study.optimizer().history.iter().all(|f| f.is_finite()));
    assert!(study.optimizer().history.windows(2).all(|f| f[1] <= f[0]));
    assert!(report.evals <= 200);
}

#[test]
fn invalid_initial_arithmetic_is_a_typed_error_not_a_rejected_success() {
    let p = logarithmic_objective();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    assert!(matches!(ReverseStudy::new(&oracle, &[-1.0], 7, None),
        Err(ReverseStudyError::Optimizer(LbfgsError::Evaluation(
            ReverseProblemError::Reverse(ReverseError::Evaluation(OptError::EvalNonFinite { .. }))
        )))
    ));
}

#[test]
fn weights_senses_and_multiple_variable_blocks_reach_their_joint_optimum() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let s = b.var("s", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let xr = b.var_ref(x).unwrap();
    let sr = b.var_ref(s).unwrap();
    let norm = b.norm_sq(xr).unwrap();
    let component = b.component(sr, 0).unwrap();
    let three = b.konst(3.0, Dims::NONE).unwrap();
    let shift = b.sub(component, three).unwrap();
    let square = b.powi(shift, 2).unwrap();
    let reward = b.neg(square).unwrap();
    b.objective(norm, Sense::Minimize, 2.0).unwrap();
    b.objective(reward, Sense::Maximize, 0.5).unwrap();
    let p = b.finish();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseStudy::new(&oracle, &[2.0, -1.0, -4.0], 7, None).unwrap();
    let report = study.run(&StopRule::GradNorm(1e-8), 100, None).unwrap();
    assert_eq!(report.reason, StopReason::GradNorm);
    assert!(study.optimizer().x[0].abs() < 1e-8);
    assert!(study.optimizer().x[1].abs() < 1e-8);
    assert!((study.optimizer().x[2] - 3.0).abs() < 1e-8);
}

#[test]
fn unsupported_constraints_manifolds_and_packed_dimensions_fail_closed() {
    for manifold in [Manifold::Rn { dim: 3 }, Manifold::Sphere { ambient: 3 }] {
        let mut b = ProblemBuilder::new();
        let v = b.var("x", manifold, Dims::NONE).unwrap();
        let r = b.var_ref(v).unwrap();
        let f = b.norm_sq(r).unwrap();
        b.objective(f, Sense::Minimize, 1.0).unwrap();
        if matches!(manifold, Manifold::Rn { .. }) {
            b.constraint(f, ConstraintKind::LeZero, "do-not-drop").unwrap();
        }
        let p = b.finish();
        let oracle = ReverseProblem::new(&p, limits()).unwrap();
        let error = ReverseStudy::new(&oracle, &[1.0, 0.0, 0.0], 7, None).unwrap_err();
        assert!(matches!(error, ReverseStudyError::ConstraintsUnsupported | ReverseStudyError::NonEuclidean { variable: 0 }));
    }
    let p = quadratic(2, limited(1));
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    assert!(matches!(ReverseStudy::new(&oracle, &[1.0], 7, None), Err(ReverseStudyError::PackedPointLength { expected: 2, actual: 1 })));
}

#[test]
fn cancelled_context_preserves_checkpoint_and_fresh_context_resumes() {
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    let p = rosenbrock();
    let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let initial = ReverseStudy::new(&oracle, &[-1.2, 1.0], 7, None).unwrap();
    let rule = StopRule::GradNorm(0.0);
    let mut straight = initial.clone();
    straight.run(&rule, 12, None).unwrap();
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey { seed: 0, kernel_id: 1, tile: 0, iteration: 0 };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        let mut interrupted = initial;
        interrupted.run(&rule, 4, Some(&cx)).unwrap();
        let before = format!("{:?}", interrupted.optimizer());
        gate.request();
        assert_eq!(interrupted.run(&rule, 8, Some(&cx)).unwrap_err(), ReverseStudyError::Cancelled);
        assert_eq!(format!("{:?}", interrupted.optimizer()), before);
        let fresh_gate = CancelGate::new();
        let fresh = Cx::new(&fresh_gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        interrupted.run(&rule, 8, Some(&fresh)).unwrap();
        assert_eq!(format!("{:?}", interrupted.optimizer()), format!("{:?}", straight.optimizer()));
    });
}
