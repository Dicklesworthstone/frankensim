//! End-to-end IR/AD/Steihaug solves and retained-work boundaries.
use fs_ascent::{ReverseNewtonError, ReverseNewtonStop, ReverseNewtonStudy, StopReason, StopRule};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, EvalLimit, Manifold, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use std::num::NonZeroU64;

fn limits() -> ReverseLimits { ReverseLimits { max_nodes: 100_000, max_scalar_slots: 1_000_000 } }
fn quadratic(n: u32, cap: Option<u64>) -> Problem {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: n }, Dims::NONE).unwrap();
    let r = b.var_ref(x).unwrap(); let z = b.norm_sq(r).unwrap();
    b.objective(z, Sense::Minimize, 1.0).unwrap();
    if let Some(cap) = cap { b.set_eval_limit(EvalLimit::Limited(NonZeroU64::new(cap).unwrap())); }
    b.finish()
}
fn rosenbrock() -> Problem {
    let mut b = ProblemBuilder::new();
    let v = b.var("xy", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let r = b.var_ref(v).unwrap(); let x = b.component(r, 0).unwrap(); let y = b.component(r, 1).unwrap();
    let one = b.konst(1.0, Dims::NONE).unwrap(); let hundred = b.konst(100.0, Dims::NONE).unwrap();
    let xx = b.powi(x, 2).unwrap(); let a = b.sub(y, xx).unwrap(); let aa = b.powi(a, 2).unwrap();
    let aa = b.mul(hundred, aa).unwrap(); let c = b.sub(one, x).unwrap(); let cc = b.powi(c, 2).unwrap();
    let root = b.add(aa, cc).unwrap(); b.objective(root, Sense::Minimize, 1.0).unwrap(); b.finish()
}
fn same(a: &ReverseNewtonStudy<'_, '_>, b: &ReverseNewtonStudy<'_, '_>) {
    let (a0,b0) = (a.snapshot(),b.snapshot());
    assert_eq!(a0.x.iter().map(|x|x.to_bits()).collect::<Vec<_>>(),b0.x.iter().map(|x|x.to_bits()).collect::<Vec<_>>());
    assert_eq!(a0.f.to_bits(), b0.f.to_bits()); assert_eq!(a0.evals,b0.evals); assert_eq!(a0.hv_evals,b0.hv_evals);
    assert_eq!(a0.iters,b0.iters); assert_eq!(a0.negative_curvature_hits,b0.negative_curvature_hits);
    assert_eq!(a.radius().to_bits(), b.radius().to_bits()); assert_eq!(a.history(),b.history()); assert_eq!(a.gradient(),b.gradient());
}

#[test]
fn high_dimensional_quadratic_needs_no_coordinate_perturbations() {
    let p = quadratic(256, None); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseNewtonStudy::new(&oracle, &vec![0.125;256], None).unwrap();
    let report = study.run(&StopRule::GradNorm(1e-12), 20, 100, None).unwrap();
    assert_eq!(report.stop, ReverseNewtonStop::Stopped(StopReason::GradNorm));
    assert!(report.solution.x.iter().all(|x|x.abs()<1e-12));
    assert_eq!(report.solution.evals, 3);
    assert_eq!(report.solution.hv_evals, 4);
}

#[test]
fn rosenbrock_runs_through_the_actual_reverse_hessian() {
    let p = rosenbrock(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let mut study = ReverseNewtonStudy::new(&oracle, &[-1.2,1.0], None).unwrap();
    let report = study.run(&StopRule::GradNorm(1e-8), 100, 1000, None).unwrap();
    assert_eq!(report.stop, ReverseNewtonStop::Stopped(StopReason::GradNorm));
    assert!(report.solution.f < 1e-16);
    assert!(report.solution.x.iter().all(|x|(x-1.0).abs()<1e-7));
    assert!(report.solution.hv_evals>report.solution.evals);
}

#[test]
fn all_thirteen_complete_iteration_splits_preserve_cache_radius_and_work() {
    let p = rosenbrock(); let oracle = ReverseProblem::new(&p, limits()).unwrap();
    let start = ReverseNewtonStudy::new(&oracle, &[-1.2,1.0], None).unwrap();
    let rule = StopRule::GradNorm(1e-12);
    let mut straight=start.clone();straight.run(&rule,12,1000,None).unwrap();
    for cut in 0..=12 {
        let mut split=start.clone();split.run(&rule,cut,1000,None).unwrap();split.run(&rule,12-cut,1000,None).unwrap();
        same(&straight,&split);
    }
}

#[test]
fn hessian_budget_preserves_the_accepted_model_and_counts_partial_krylov_work() {
    let p=quadratic(1,None);let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut study=ReverseNewtonStudy::new(&oracle,&[1.0],None).unwrap();
    let rule=StopRule::GradNorm(1e-12);
    assert_eq!(study.run(&rule,10,0,None).unwrap().stop,ReverseNewtonStop::HessianBudget);
    assert_eq!(study.snapshot().hv_evals,0);
    let report=study.run(&rule,10,1,None).unwrap();
    assert_eq!(report.stop,ReverseNewtonStop::HessianBudget);
    assert_eq!(report.solution.x,vec![1.0]);assert_eq!(report.solution.evals,1);assert_eq!(report.solution.hv_evals,1);
    assert_eq!(report.solution.iters,0);assert_eq!(study.radius(),1.0);
    let report=study.run(&rule,10,3,None).unwrap();
    assert_eq!(report.stop,ReverseNewtonStop::Stopped(StopReason::GradNorm));
    assert_eq!(report.solution.x,vec![0.0]);assert_eq!(report.solution.evals,2);assert_eq!(report.solution.hv_evals,3);
}

#[test]
fn objective_limits_include_initial_and_rejected_trials_and_nested_budget_leaves() {
    let p=quadratic(1,Some(2));let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut study=ReverseNewtonStudy::new(&oracle,&[1.0],None).unwrap();
    let report=study.run(&StopRule::GradNorm(1e-8),10,100,None).unwrap();
    assert_eq!(report.stop,ReverseNewtonStop::Stopped(StopReason::Budget));assert_eq!(report.solution.evals,2);
    assert_eq!(report.solution.x,vec![0.0]);let checkpoint=study.clone();
    study.run(&StopRule::GradNorm(1e-8),10,100,None).unwrap();same(&study,&checkpoint);
    let p=rosenbrock();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut study=ReverseNewtonStudy::new(&oracle,&[-1.2,1.0],None).unwrap();
    let rule=StopRule::All(vec![StopRule::GradNorm(0.0),StopRule::Budget(3)]);
    let report=study.run(&rule,100,1000,None).unwrap();
    assert_eq!(report.stop,ReverseNewtonStop::Stopped(StopReason::Budget));assert_eq!(report.solution.evals,3);
    let checkpoint=study.clone();study.run(&StopRule::Budget(1),100,1000,None).unwrap();same(&study,&checkpoint);
}

#[test]
fn logarithmic_domain_rejection_contracts_radius_then_recovers() {
    let mut b=ProblemBuilder::new();let v=b.var("x",Manifold::Rn{dim:1},Dims::NONE).unwrap();
    let r=b.var_ref(v).unwrap();let x=b.component(r,0).unwrap();let ln=b.ln(x).unwrap();let f=b.sub(x,ln).unwrap();
    b.objective(f,Sense::Minimize,1.0).unwrap();let p=b.finish();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut study=ReverseNewtonStudy::new(&oracle,&[3.0],None).unwrap();
    let report=study.run(&StopRule::GradNorm(1e-8),100,1000,None).unwrap();
    assert_eq!(report.stop,ReverseNewtonStop::Stopped(StopReason::GradNorm));
    assert!((report.solution.x[0]-1.0).abs()<1e-7);assert!(study.domain_rejections()>0);assert!(study.last_rejection().is_some());
}

#[test]
fn second_order_overflow_is_an_error_not_a_corrupted_or_converged_point() {
    let mut b=ProblemBuilder::new();let v=b.var("x",Manifold::Rn{dim:1},Dims::NONE).unwrap();
    let r=b.var_ref(v).unwrap();let x=b.component(r,0).unwrap();let a=b.konst(1e155,Dims::NONE).unwrap();
    let ax=b.mul(a,x).unwrap();let square=b.powi(ax,2).unwrap();let half=b.konst(0.5,Dims::NONE).unwrap();
    let f=b.mul(half,square).unwrap();b.objective(f,Sense::Minimize,1.0).unwrap();
    let p=b.finish();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut study=ReverseNewtonStudy::new(&oracle,&[1e-310],None).unwrap();let before=study.snapshot();
    assert!(matches!(study.run(&StopRule::GradNorm(1e-12),10,100,None),Err(ReverseNewtonError::Evaluation(_))));
    let after=study.snapshot();assert_eq!(after.x,before.x);assert_eq!(after.f.to_bits(),before.f.to_bits());
    assert_eq!(after.evals,1);assert_eq!(after.iters,0);assert_eq!(after.hv_evals,1);assert_eq!(study.radius(),1.0);
}

#[test]
fn cancelled_context_preserves_checkpoint_and_can_resume() {
    use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
    let p=rosenbrock();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut study=ReverseNewtonStudy::new(&oracle,&[-1.2,1.0],None).unwrap();let checkpoint=study.clone();
    let gate=CancelGate::new();gate.request();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey{seed:0,kernel_id:1,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        assert!(matches!(study.run(&StopRule::GradNorm(1e-8),100,1000,Some(&cx)),Err(ReverseNewtonError::Cancelled)));
    });
    same(&study,&checkpoint);
    assert_eq!(study.run(&StopRule::GradNorm(1e-8),100,1000,None).unwrap().stop,ReverseNewtonStop::Stopped(StopReason::GradNorm));
}

#[test]
fn structural_and_rule_refusals_do_not_silently_change_the_problem() {
    for manifold in [Manifold::So3,Manifold::Rn{dim:1}] {
        let mut b=ProblemBuilder::new();let v=b.var("x",manifold,Dims::NONE).unwrap();
        let r=b.var_ref(v).unwrap();let f=b.norm_sq(r).unwrap();b.objective(f,Sense::Minimize,1.0).unwrap();
        if matches!(manifold,Manifold::Rn{..}) { b.constraint(f,ConstraintKind::LeZero,"bound").unwrap(); }
        let p=b.finish();let oracle=ReverseProblem::new(&p,limits()).unwrap();
        assert!(matches!(ReverseNewtonStudy::new(&oracle,&[],None),Err(ReverseNewtonError::NonEuclidean{..}|ReverseNewtonError::ConstraintsUnsupported)));
    }
    let p=quadratic(2,None);let oracle=ReverseProblem::new(&p,limits()).unwrap();
    assert!(matches!(ReverseNewtonStudy::new(&oracle,&[1.0],None),Err(ReverseNewtonError::PackedPointLength{..})));
    let mut study=ReverseNewtonStudy::new(&oracle,&[1.0,2.0],None).unwrap();let before=study.clone();
    assert!(matches!(study.run(&StopRule::GradNorm(f64::NAN),10,100,None),Err(ReverseNewtonError::InvalidRule)));
    same(&study,&before);
}
