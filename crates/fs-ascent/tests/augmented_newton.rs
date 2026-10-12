//! Live IR constraints on the production matrix-free Newton kernel.
use fs_ascent::{ReverseNewtonError, ReverseNewtonStudy, StopRule};
use fs_ascent::reverse_newton::{AugmentedNewtonBudget as Budget, AugmentedNewtonOptions as Options,
    AugmentedNewtonReport as Report, AugmentedNewtonStop as Stop, AugmentedNewtonStudy as Study};
use fs_opt::{AugmentedLagrangian, ConstraintKind, Manifold, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_opt::reverse::ReverseLimits;
use fs_qty::Dims;

fn limits() -> ReverseLimits { ReverseLimits { max_nodes: 512, max_scalar_slots: 4096 } }
fn budget() -> Budget { Budget { objective_evals: 1000, hessian_products: 4000, multiplier_updates: 100 } }
fn close(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(a.len(), b.len());
    assert!(a.iter().zip(b).all(|(a,b)|(a-b).abs()<tol*(1.0+b.abs())), "{a:?} != {b:?}");
}
fn equality() -> Problem {
    let mut b=ProblemBuilder::new();let v=b.var("x",Manifold::Rn {dim:2},Dims::NONE).unwrap();
    let r=b.var_ref(v).unwrap();let f=b.norm_sq(r).unwrap();b.objective(f,Sense::Minimize,1.0).unwrap();
    let x=b.component(r,0).unwrap();let y=b.component(r,1).unwrap();let sum=b.add(x,y).unwrap();
    let one=b.konst(1.0,Dims::NONE).unwrap();let c=b.sub(sum,one).unwrap();
    b.constraint(c,ConstraintKind::EqZero,"sum").unwrap();b.finish()
}
fn inequality(nonlinear: bool) -> Problem {
    let mut b=ProblemBuilder::new();let v=b.var("x",Manifold::Rn {dim:1},Dims::NONE).unwrap();
    let r=b.var_ref(v).unwrap();let x=b.component(r,0).unwrap();let two=b.konst(2.0,Dims::NONE).unwrap();
    let dx=b.sub(x,two).unwrap();let f=b.powi(dx,2).unwrap();b.objective(f,Sense::Minimize,1.0).unwrap();
    let one=b.konst(1.0,Dims::NONE).unwrap();
    let left=if nonlinear { b.powi(x,2).unwrap() } else {x};
    let c=b.sub(left,one).unwrap();b.constraint(c,ConstraintKind::LeZero,"upper").unwrap();
    let lower=b.neg(x).unwrap();b.constraint(lower,ConstraintKind::LeZero,"lower").unwrap();b.finish()
}
fn verify_kkt(oracle:&ReverseProblem<'_>, report:&Report, tol:f64) {
    let tape=oracle.evaluate(&[report.x.clone()]).unwrap();
    let g=tape.lagrangian_gradient(&report.multipliers).unwrap();
    let norm=g.iter().flatten().map(|x|x.abs()).fold(0.0_f64,f64::max);
    close(&[report.kkt.stationarity],&[norm],1e-12);
    assert!(report.kkt.within_tolerance(tol),"{report:?}");
    assert_eq!(report.constraints,tape.constraint_values());
    assert_eq!(report.objective.to_bits(),tape.objective_value().to_bits());
}

#[test]
fn equality_reaches_analytic_primal_dual_solution_without_dense_jacobians() {
    let p=equality();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut study=Study::new(&oracle,&[3.0,-2.0],Options::default(),None).unwrap();
    let report=study.run(1e-6,200,budget(),None).unwrap();
    assert_eq!(report.stop,Stop::Converged);close(&report.x,&[0.5,0.5],1e-6);
    close(&report.multipliers,&[-1.0],1e-6);verify_kkt(&oracle,&report,1e-6);
    assert!(report.multiplier_updates>0 && report.hessian_products>0);
}

#[test]
fn active_inactive_and_nonlinear_inequalities_retain_correct_dual_signs() {
    for nonlinear in [false,true] {
        let p=inequality(nonlinear);let oracle=ReverseProblem::new(&p,limits()).unwrap();
        let mut study=Study::new(&oracle,&[0.25],Options::default(),None).unwrap();
        let report=study.run(2e-5,300,budget(),None).unwrap();
        assert_eq!(report.stop,Stop::Converged,"{report:?}");close(&report.x,&[1.0],2e-5);
        close(&report.multipliers,&[if nonlinear {1.0} else {2.0},0.0],3e-5);
        verify_kkt(&oracle,&report,2e-5);
    }
}

#[test]
fn every_complete_newton_split_replays_points_merits_multipliers_and_work() {
    let p=equality();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let original=Study::new(&oracle,&[3.0,-2.0],Options::default(),None).unwrap();
    let mut whole=original.clone();let expected=whole.run(1e-6,14,budget(),None).unwrap();
    for cut in 0..=14 {
        let mut split=original.clone();split.run(1e-6,cut,budget(),None).unwrap();
        let result=split.run(1e-6,14-cut,budget(),None).unwrap();
        assert_eq!(result,expected,"cut {cut}");assert_eq!(split.frozen_multipliers(),whole.frozen_multipliers());
        assert_eq!(split.inner().history(),whole.inner().history());assert_eq!(split.inner().radius(),whole.inner().radius());
    }
}

#[test]
fn exact_objective_hessian_and_multiplier_ceilings_resume_without_work_refunds() {
    let p=equality();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let original=Study::new(&oracle,&[0.0,0.0],Options::default(),None).unwrap();
    let mut study=original.clone();
    let objective=study.run(1e-6,5,Budget {objective_evals:1,..budget()},None).unwrap();
    assert_eq!(objective.stop,Stop::ObjectiveBudget);assert_eq!(objective.objective_evals,1);assert_eq!(objective.hessian_products,0);
    let before=study.inner().snapshot();let radius=study.inner().radius();
    let hv=study.run(1e-6,5,Budget {hessian_products:1,..budget()},None).unwrap();
    assert_eq!(hv.stop,Stop::HessianBudget);assert_eq!(hv.hessian_products,1);
    assert_eq!(study.inner().snapshot().x,before.x);assert_eq!(study.inner().radius(),radius);
    assert_eq!(study.run(1e-6,5,Budget {hessian_products:1,..budget()},None).unwrap(),hv);
    let outer=study.run(1e-6,50,Budget {multiplier_updates:0,..budget()},None).unwrap();
    assert_eq!(outer.stop,Stop::MultiplierBudget);assert_eq!(outer.multiplier_attempts,0);
    assert!(!outer.kkt.within_tolerance(1e-6));assert!(outer.hessian_products>hv.hessian_products);
    let end=study.run(1e-6,100,budget(),None).unwrap();assert_eq!(end.stop,Stop::Converged);
    assert!(end.multiplier_attempts>0 && end.hessian_products>outer.hessian_products);
}

#[test]
fn infeasible_stationary_subproblem_stops_at_policy_not_false_convergence() {
    let mut b=ProblemBuilder::new();let v=b.var("x",Manifold::Rn {dim:1},Dims::NONE).unwrap();let r=b.var_ref(v).unwrap();
    let f=b.norm_sq(r).unwrap();b.objective(f,Sense::Minimize,1.0).unwrap();
    let one=b.konst(1.0,Dims::NONE).unwrap();b.constraint(one,ConstraintKind::EqZero,"impossible").unwrap();let p=b.finish();
    let oracle=ReverseProblem::new(&p,limits()).unwrap();let mut study=Study::new(&oracle,&[0.0],Options::new(1.0,10.0,10.0,0.5).unwrap(),None).unwrap();
    let report=study.run(1e-6,20,budget(),None).unwrap();
    assert_eq!(report.stop,Stop::PenaltyLimit);assert_eq!(report.kkt.stationarity,0.0);
    assert_eq!(report.kkt.feasibility,1.0);assert_eq!(report.objective_evals,1);assert_eq!(report.hessian_products,0);
    assert_eq!(report.multiplier_attempts,2);assert_eq!(report.penalty,10.0);
}

#[test]
fn explicit_frozen_augmentation_supports_preconditioning_but_not_kkt_claims() {
    let p=equality();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    assert!(matches!(ReverseNewtonStudy::new(&oracle,&[0.0,0.0],None),Err(ReverseNewtonError::ConstraintsUnsupported)));
    let model=AugmentedLagrangian::new(&p,&[0.0],10.0).unwrap();
    let mut direct=ReverseNewtonStudy::new_augmented(&oracle,&[0.0,0.0],model,None).unwrap();let mut pre=direct.clone();
    let a=direct.run(&StopRule::GradNorm(1e-9),10,100,None).unwrap();
    let b=pre.run_preconditioned(&mut |_,v| v.to_vec(),&StopRule::GradNorm(1e-9),10,100,100,None).unwrap();
    assert_eq!(a.solution.x,b.solution.x);assert_eq!(direct.history(),pre.history());
    close(&a.solution.x,&[5.0/11.0,5.0/11.0],1e-9);
    assert!(direct.constraint_values()[0].abs()>0.01,"an inner solution is not feasible");
    assert_ne!(direct.objective_value().to_bits(),a.solution.f.to_bits());
}

#[test]
fn invalid_controls_identity_and_cancellation_preserve_current_checkpoint() {
    let p=equality();let oracle=ReverseProblem::new(&p,limits()).unwrap();
    for args in [(0.0,10.0,100.0,0.5),(1.0,1.0,100.0,0.5),(2.0,10.0,1.0,0.5),(1.0,10.0,100.0,1.0)] {
        assert!(Options::new(args.0,args.1,args.2,args.3).is_err());
    }
    let other=equality();let wrong=AugmentedLagrangian::new(&other,&[0.0],10.0).unwrap();
    assert!(ReverseNewtonStudy::new_augmented(&oracle,&[0.0,0.0],wrong,None).is_err());
    let mut study=Study::new(&oracle,&[3.0,-2.0],Options::default(),None).unwrap();
    let before=study.run(1e-6,0,budget(),None).unwrap();
    for tol in [0.0,-1.0,f64::NAN,f64::INFINITY,f64::from_bits(1)] {assert!(study.run(tol,1,budget(),None).is_err());}
    let gate=fs_exec::CancelGate::new();gate.request();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=fs_exec::Cx::new(&gate,arena,fs_exec::StreamKey {seed:0,kernel_id:1,tile:0,iteration:0},fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        assert!(matches!(study.run(1e-6,1,budget(),Some(&cx)),Err(ReverseNewtonError::Cancelled)));
    });
    assert_eq!(study.run(1e-6,0,budget(),None).unwrap(),before);
    assert_eq!(study.run(1e-6,200,budget(),None).unwrap().stop,Stop::Converged);
}
