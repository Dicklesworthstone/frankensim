//! End-to-end reverse-mode solves with heterogeneous manifold variables.
#[allow(dead_code)]
#[path = "../examples/reverse_manifold.rs"]
mod example;

use fs_ascent::{ReverseManifoldError, ReverseManifoldStudy, StopReason, StopRule};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::{
    ConstraintKind, EvalLimit, Manifold, Problem, ProblemBuilder, ReverseProblem,
    Sense, reverse::ReverseLimits,
};
use fs_qty::Dims;
use std::num::NonZeroU64;

fn limits() -> ReverseLimits { ReverseLimits { max_nodes: 4096, max_scalar_slots: 16384 } }
fn bits(v: &[f64]) -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() }
fn same(a: &ReverseManifoldStudy<'_, '_>, b: &ReverseManifoldStudy<'_, '_>) {
    assert_eq!(bits(a.point()), bits(b.point()));
    assert_eq!(bits(a.gradient()), bits(b.gradient()));
    assert_eq!(bits(a.history()), bits(b.history()));
    assert_eq!(a.evaluations(), b.evaluations());
}
fn log_sphere(cap: u64) -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("unit-direction", Manifold::Sphere { ambient: 3 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap();
    let x = b.component(r, 0).unwrap();
    let y = b.component(r, 1).unwrap();
    let l = b.ln(x).unwrap();
    let ten = b.konst(10.0, Dims::NONE).unwrap();
    let term = b.mul(ten, y).unwrap();
    let sum = b.add(l, term).unwrap();
    b.objective(sum, Sense::Maximize, 1.0).unwrap();
    b.set_eval_limit(EvalLimit::Limited(NonZeroU64::new(cap).unwrap()));
    b.finish()
}

#[test]
fn coupled_position_direction_rotation_and_frame_reach_known_minimum() {
    let (problem, point) = example::fixture();
    let oracle = ReverseProblem::new(&problem, limits()).unwrap();
    let mut s = ReverseManifoldStudy::new(&oracle, &point, 7, None).unwrap();
    assert_eq!(s.point().len(), 15);
    assert_eq!(s.gradient().len(), 14);
    let report = s.run(&StopRule::GradNorm(1e-7), 100, None).unwrap();
    assert_eq!(report.reason, StopReason::GradNorm);
    assert!(report.f < 1e-12);
    assert!(report.evals < 100);
    assert!((s.point()[1]+0.25).abs() < 1e-6);
    s.manifold().validate_point(s.point()).unwrap();
    s.manifold().validate_parameter_tangent(s.point(), s.gradient()).unwrap();
    assert!(s.history().windows(2).all(|f| f[1] <= f[0]));
}

#[test]
fn every_completed_step_split_preserves_curvature_and_work() {
    let (problem, point) = example::fixture();
    let oracle = ReverseProblem::new(&problem, limits()).unwrap();
    let initial = ReverseManifoldStudy::new(&oracle, &point, 7, None).unwrap();
    let rule = StopRule::GradNorm(0.0);
    let mut straight = initial.clone();
    straight.run(&rule, 12, None).unwrap();
    for split in 0..=12 {
        let mut s = initial.clone();
        s.run(&rule, split, None).unwrap();
        let mut checkpoint = s.clone();
        s.run(&rule, 12-split, None).unwrap();
        checkpoint.run(&rule, 12-split, None).unwrap();
        same(&s, &straight); same(&checkpoint, &straight);
        // Exercise the retained pairs after reaching the same endpoint.
        s.run(&StopRule::GradNorm(1e-7), 30, None).unwrap();
        checkpoint.run(&StopRule::GradNorm(1e-7), 30, None).unwrap();
        same(&s, &checkpoint);
    }
}

#[test]
fn primitive_domain_rejections_backtrack_without_false_zero_gradient() {
    let p = log_sphere(100);
    let o = ReverseProblem::new(&p, limits()).unwrap();
    let mut s = ReverseManifoldStudy::new(&o, &[0.6, 0.8, 0.0], 7, None).unwrap();
    let r = s.run(&StopRule::GradNorm(1e-7), 100, None).unwrap();
    assert_eq!(r.reason, StopReason::GradNorm);
    assert!(r.rejected_trials > 0);
    assert!(s.last_rejection().is_some());
    let y = (-1.0+401.0_f64.sqrt())/20.0;
    assert!((s.point()[1]-y).abs() < 1e-6);
    assert!(s.point()[0] > 0.0);
}

#[test]
fn hard_budget_inside_all_preserves_rejected_endpoint_and_allows_retry() {
    let p = log_sphere(100);
    let o = ReverseProblem::new(&p, limits()).unwrap();
    let mut s = ReverseManifoldStudy::new(&o, &[0.6, 0.8, 0.0], 7, None).unwrap();
    let before = bits(s.point());
    let r = s.run(&StopRule::All(vec![StopRule::GradNorm(0.0), StopRule::Budget(2)]), 100, None).unwrap();
    assert_eq!(r.reason, StopReason::Budget);
    assert_eq!(r.evals, 2); assert_eq!(r.iters, 0);
    assert_eq!(bits(s.point()), before);
    let r = s.run(&StopRule::GradNorm(1e-7), 100, None).unwrap();
    assert_eq!(r.reason, StopReason::GradNorm);
}

#[test]
fn problem_budget_cannot_be_raised_and_cached_runs_spend_nothing() {
    let p = log_sphere(2);
    let o = ReverseProblem::new(&p, limits()).unwrap();
    let mut s = ReverseManifoldStudy::new(&o, &[0.6, 0.8, 0.0], 7, None).unwrap();
    let before = s.clone();
    s.run(&StopRule::Budget(100), 0, None).unwrap(); same(&s, &before);
    s.run(&StopRule::GradNorm(0.0), 100, None).unwrap();
    let checkpoint = s.clone();
    assert_eq!(s.run(&StopRule::Budget(1000), 100, None).unwrap().reason, StopReason::Budget);
    same(&s, &checkpoint);
}

#[test]
fn cancellation_retains_checkpoint_and_fresh_context_resumes() {
    let (p, x) = example::fixture();
    let o = ReverseProblem::new(&p, limits()).unwrap();
    let mut s = ReverseManifoldStudy::new(&o, &x, 7, None).unwrap();
    s.run(&StopRule::GradNorm(0.0), 3, None).unwrap();
    let mut reference = s.clone();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new(); gate.request();
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 7, kernel_id: 1, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        assert!(matches!(s.run(&StopRule::GradNorm(1e-7), 100, Some(&cx)), Err(ReverseManifoldError::Cancelled)));
    });
    same(&s, &reference);
    let fresh = CancelGate::new();
    pool.scope(|arena| {
        let cx = Cx::new(&fresh, arena, StreamKey { seed: 7, kernel_id: 1, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic);
        s.run(&StopRule::GradNorm(1e-7), 100, Some(&cx)).unwrap();
    });
    reference.run(&StopRule::GradNorm(1e-7), 100, None).unwrap(); same(&s, &reference);
}

#[test]
fn quaternion_antipodes_canonicalize_before_initial_valuation() {
    let (p, x) = example::fixture();
    let o = ReverseProblem::new(&p, limits()).unwrap();
    let mut opposite = x.clone(); for q in &mut opposite[5..9] { *q = -*q; }
    let mut a = ReverseManifoldStudy::new(&o, &x, 7, None).unwrap();
    let mut b = ReverseManifoldStudy::new(&o, &opposite, 7, None).unwrap();
    same(&a, &b);
    a.run(&StopRule::GradNorm(1e-7), 100, None).unwrap();
    b.run(&StopRule::GradNorm(1e-7), 100, None).unwrap(); same(&a, &b);
}

#[test]
fn zero_memory_is_valid_retracted_steepest_descent() {
    let (p, x) = example::fixture();
    let o = ReverseProblem::new(&p, limits()).unwrap();
    let mut s = ReverseManifoldStudy::new(&o, &x, 0, None).unwrap();
    let r = s.run(&StopRule::GradNorm(1e-7), 300, None).unwrap();
    assert_eq!(r.reason, StopReason::GradNorm); assert!(r.f < 1e-12);
}

#[test]
fn constraints_bad_points_and_bad_rules_are_typed_refusals() {
    let mut b = ProblemBuilder::new();
    let v = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let root = b.var_ref(v).unwrap(); let x = b.component(root, 0).unwrap();
    b.objective(x, Sense::Minimize, 1.0).unwrap();
    b.constraint(x, ConstraintKind::LeZero, "upper-bound").unwrap();
    let p = b.finish(); let o = ReverseProblem::new(&p, limits()).unwrap();
    assert!(matches!(ReverseManifoldStudy::new(&o, &[1.0], 7, None), Err(ReverseManifoldError::ConstraintsUnsupported)));
    let p = log_sphere(100); let o = ReverseProblem::new(&p, limits()).unwrap();
    assert!(ReverseManifoldStudy::new(&o, &[1.0, 1.0, 0.0], 7, None).is_err());
    assert!(ReverseManifoldStudy::new(&o, &[0.0, 1.0, 0.0], 7, None).is_err());
    let mut s = ReverseManifoldStudy::new(&o, &[0.6, 0.8, 0.0], 7, None).unwrap();
    let before = s.clone();
    assert!(s.run(&StopRule::GradNorm(f64::NAN), 1, None).is_err()); same(&s, &before);
}
