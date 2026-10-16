//! G0/G3/G4: real IR/PHR derivatives through the existing inverse-action kernel.
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions,
    AugmentedNewtonStop, AugmentedNewtonStudy, ReverseNewtonError};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{AugmentedLagrangian, ConstraintKind, Manifold, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;

fn limits() -> ReverseLimits { ReverseLimits {max_nodes:1024,max_scalar_slots:8192} }
fn budget() -> AugmentedNewtonBudget {
    AugmentedNewtonBudget {objective_evals:1000,hessian_products:10_000,multiplier_updates:100}
}
fn quadratic(diagonal: &[f64], target: &[f64], inequality: bool) -> Problem {
    let mut b = ProblemBuilder::new();
    let var = b.var("design",Manifold::Rn {dim:diagonal.len() as u32},Dims::NONE).unwrap();
    let x = b.var_ref(var).unwrap(); let mut total = b.konst(0.0,Dims::NONE).unwrap();
    for (i,(&d,&t)) in diagonal.iter().zip(target).enumerate() {
        let xi=b.component(x,i as u32).unwrap(); let t=b.konst(t,Dims::NONE).unwrap();
        let r=b.sub(xi,t).unwrap(); let r2=b.powi(r,2).unwrap();
        let weight=b.konst(0.5*d,Dims::NONE).unwrap(); let term=b.mul(weight,r2).unwrap();
        total=b.add(total,term).unwrap();
    }
    b.objective(total,Sense::Minimize,1.0).unwrap();
    let first=b.component(x,0).unwrap();
    let cap=b.konst(if inequality {1.0} else {0.0},Dims::NONE).unwrap();
    let c=b.sub(first,cap).unwrap();
    b.constraint(c,if inequality {ConstraintKind::LeZero} else {ConstraintKind::EqZero},"cap").unwrap();
    b.finish()
}
fn inverse(r:&[f64],model:&AugmentedLagrangian<'_>,constraints:&[f64],d:&[f64],inequality:bool)->Vec<f64> {
    let active=!inequality || model.multipliers()[0]+model.penalty()*constraints[0]>0.0;
    r.iter().zip(d).enumerate().map(|(i,(&r,&d))|
        r/(d+if i==0 && active {model.penalty()} else {0.0})).collect()
}

#[test]
fn g0_stiff_constrained_quadratic_uses_one_krylov_product_with_exact_inverse() {
    let d:Vec<_>=(0..12).map(|i|4.0_f64.powi(i)).collect();
    let mut t=vec![0.005;d.len()];t[0]=0.0;
    let p=quadratic(&d,&t,false);let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut fast=AugmentedNewtonStudy::new(&oracle,&vec![0.0;d.len()],AugmentedNewtonOptions::default(),None).unwrap();
    let mut plain=fast.clone();
    let mut apply=|_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]|inverse(r,m,c,&d,false);
    let report=fast.run_preconditioned(&mut apply,1e-8,1,budget(),100,None).unwrap();
    let slow=plain.run(1e-8,1,budget(),None).unwrap();
    assert_eq!(report.stop,AugmentedNewtonStop::Converged);
    assert_eq!(report.hessian_products,2,"one Krylov and one final-model product");
    assert_eq!(report.preconditioner_evals,1);
    assert!(slow.hessian_products>report.hessian_products);
    for (x,t) in report.x.iter().zip(t) {assert!((x-t).abs()<1e-10);}
}

#[test]
fn g3_inverse_observes_new_penalty_multipliers_and_accepted_constraints_after_retarget() {
    let d=[1.0,1e4];let p=quadratic(&d,&[2.0,0.05],true);
    let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut whole=AugmentedNewtonStudy::new(&oracle,&[0.0,0.0],AugmentedNewtonOptions::default(),None).unwrap();
    let mut split=whole.clone();let mut seen=Vec::new();
    let mut apply=|x:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]| {
        assert_eq!(c.len(),1);assert_eq!(c[0].to_bits(),(x[0]-1.0).to_bits());
        seen.push((m.multipliers()[0],m.penalty()));inverse(r,m,c,&d,true)
    };
    let end=whole.run_preconditioned(&mut apply,1e-6,150,budget(),1000,None).unwrap();
    assert_eq!(end.stop,AugmentedNewtonStop::Converged,"{end:?}");
    assert!(end.kkt.within_tolerance(1e-6));assert!((end.x[0]-1.0).abs()<2e-6);
    assert!((end.x[1]-0.05).abs()<1e-10);assert!(end.multiplier_updates>0);
    assert!(seen.iter().any(|(lambda,_)|*lambda>0.0));
    assert_eq!(end.preconditioner_evals,seen.len());
    let mut replay=|_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]|inverse(r,m,c,&d,true);
    split.run_preconditioned(&mut replay,1e-6,3,budget(),1000,None).unwrap();
    assert_eq!(split.run_preconditioned(&mut replay,1e-6,147,budget(),1000,None).unwrap(),end);
}

#[test]
fn g4_malformed_nonpositive_and_panicking_actions_retain_state_and_spent_work() {
    let d=[1.0,100.0];let p=quadratic(&d,&[0.0,0.1],false);
    let oracle=ReverseProblem::new(&p,limits()).unwrap();
    for fault in 0..4 {
        let mut s=AugmentedNewtonStudy::new(&oracle,&[0.0,0.0],AugmentedNewtonOptions::default(),None).unwrap();
        let before=s.inner().snapshot();let radius=s.inner().radius();let history=s.inner().history().to_vec();
        let mut bad=|_:&[f64],r:&[f64],_:&AugmentedLagrangian<'_>,_:&[f64]|match fault {
            0=>Vec::new(),1=>vec![f64::NAN;r.len()],2=>r.iter().map(|x|-x).collect(),_=>panic!("injected inverse fault"),
        };
        assert!(matches!(s.run_preconditioned(&mut bad,1e-8,1,budget(),10,None),Err(ReverseNewtonError::Preconditioner(_))));
        assert_eq!(s.inner().snapshot().x,before.x);assert_eq!(s.inner().snapshot().f.to_bits(),before.f.to_bits());
        assert_eq!(s.inner().radius(),radius);assert_eq!(s.inner().history(),history);
        assert_eq!(s.preconditioner_evals(),1);assert_eq!(s.newton_iterations(),0);
        let mut good=|_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]|inverse(r,m,c,&d,false);
        let exhausted=s.run_preconditioned(&mut good,1e-8,1,budget(),1,None).unwrap();
        assert_eq!(exhausted.stop,AugmentedNewtonStop::PreconditionerBudget);
        assert_eq!(exhausted.preconditioner_evals,1);
        let repaired=s.run_preconditioned(&mut good,1e-8,1,budget(),10,None).unwrap();
        assert_eq!(repaired.stop,AugmentedNewtonStop::Converged);
        assert_eq!(repaired.preconditioner_evals,2);
    }
}

#[test]
fn g4_zero_iterations_and_distinct_hard_work_limits_do_not_issue_an_inverse() {
    let p=quadratic(&[1.0,100.0],&[0.0,0.1],false);
    let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let mut s=AugmentedNewtonStudy::new(&oracle,&[0.0,0.0],AugmentedNewtonOptions::default(),None).unwrap();
    let mut never=|_:&[f64],_:&[f64],_:&AugmentedLagrangian<'_>,_:&[f64]|->Vec<f64>{panic!("unfunded inverse")};
    assert_eq!(s.run_preconditioned(&mut never,1e-8,0,budget(),0,None).unwrap().stop,AugmentedNewtonStop::IterationCap);
    assert_eq!(s.run_preconditioned(&mut never,1e-8,1,budget(),0,None).unwrap().stop,AugmentedNewtonStop::PreconditionerBudget);
    let no_hessian=AugmentedNewtonBudget {hessian_products:0,..budget()};
    assert_eq!(s.run_preconditioned(&mut never,1e-8,1,no_hessian,10,None).unwrap().stop,AugmentedNewtonStop::HessianBudget);
    let no_sample=AugmentedNewtonBudget {objective_evals:1,..budget()};
    assert_eq!(s.run_preconditioned(&mut never,1e-8,1,no_sample,10,None).unwrap().stop,AugmentedNewtonStop::ObjectiveBudget);
    assert_eq!(s.preconditioner_evals(),0);assert_eq!(s.inner().snapshot().hv_evals,0);
}

#[test]
fn g4_cancellation_inside_inverse_preserves_the_checkpoint() {
    let p=quadratic(&[1.0,100.0],&[0.0,0.1],false);
    let oracle=ReverseProblem::new(&p,limits()).unwrap();
    let pool=TilePool::for_host(1,12);let gate=CancelGate::new_clock_free();
    pool.arena_pool().scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:12,kernel_id:17,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        let mut s=AugmentedNewtonStudy::new(&oracle,&[0.0,0.0],AugmentedNewtonOptions::default(),Some(&cx)).unwrap();
        let before=s.inner().snapshot();
        let mut cancel=|_:&[f64],r:&[f64],_:&AugmentedLagrangian<'_>,_:&[f64]| {gate.request();r.to_vec()};
        assert!(s.run_preconditioned(&mut cancel,1e-8,1,budget(),10,Some(&cx)).is_err());
        assert_eq!(s.preconditioner_evals(),1);assert_eq!(s.inner().snapshot().x,before.x);
        assert_eq!(s.newton_iterations(),0);assert_eq!(s.inner().radius(),1.0);
        let mut good=|_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]|inverse(r,m,c,&[1.0,100.0],false);
        assert_eq!(s.run_preconditioned(&mut good,1e-8,1,budget(),10,None).unwrap().stop,AugmentedNewtonStop::Converged);
        assert_eq!(s.preconditioner_evals(),2);
    });
}

use fs_ascent::newton_continuation::{ContinuationError, ContinuationStop, NewtonContinuation, NewtonStage};
use fs_ascent::trust::KrylovOptions;

#[test]
fn g3_continuation_translates_one_inverse_ceiling_across_stages_and_split_runs() {
    let diagonals=[[1.0,16.0,256.0],[1.0,32.0,512.0],[1.0,64.0,1024.0]];
    let problems:Vec<_>=diagonals.iter().enumerate().map(|(i,d)|
        quadratic(d,&[0.0,0.01*(i+1) as f64,0.005*(i+1) as f64],false)).collect();
    let oracles:Vec<_>=problems.iter().map(|p|ReverseProblem::new(p,limits()).unwrap()).collect();
    let stages:Vec<_>=oracles.iter().map(|oracle|NewtonStage {oracle,tolerance:1e-8,
        options:AugmentedNewtonOptions::default(),krylov:KrylovOptions::default()}).collect();
    let mut whole=NewtonContinuation::new(&stages,&[0.0;3],None).unwrap();
    let mut limited=whole.clone();let mut split=whole.clone();
    let mut seen=Vec::new();
    let mut apply=|stage:usize,_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]| {
        seen.push(stage);inverse(r,m,c,&diagonals[stage],false)
    };
    let zero=whole.run_preconditioned(&mut apply,0,budget(),0,None).unwrap();
    assert_eq!(zero.work.preconditioner_evals,0);assert_eq!(zero.work.stage_initializations,0);
    let end=whole.run_preconditioned(&mut apply,10,budget(),3,None).unwrap();
    assert_eq!(end.stop,ContinuationStop::Complete);assert_eq!(end.work.preconditioner_evals,3);
    assert_eq!(seen,vec![0,1,2]);assert_eq!(whole.completed_stages().len(),2);
    let mut again=|stage:usize,_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]|inverse(r,m,c,&diagonals[stage],false);
    let stopped=limited.run_preconditioned(&mut again,10,budget(),1,None).unwrap();
    assert_eq!(stopped.stop,ContinuationStop::StageStopped(AugmentedNewtonStop::PreconditionerBudget));
    assert_eq!(stopped.stage,1);assert_eq!(stopped.work.preconditioner_evals,1);
    assert_eq!(stopped.solution.unwrap().preconditioner_evals,0,"new stage did not inherit a fresh allowance");
    assert_eq!(limited.run_preconditioned(&mut again,10,budget(),3,None).unwrap(),end);
    split.run_preconditioned(&mut again,1,budget(),3,None).unwrap();
    assert_eq!(split.run_preconditioned(&mut again,2,budget(),3,None).unwrap(),end);
}

#[test]
fn g4_failed_second_stage_inverse_is_charged_to_the_lifetime_budget() {
    let d=[1.0,16.0,256.0];
    let p0=quadratic(&d,&[0.0,0.01,0.02],false);let p1=quadratic(&d,&[0.0,0.03,0.04],false);
    let o0=ReverseProblem::new(&p0,limits()).unwrap();let o1=ReverseProblem::new(&p1,limits()).unwrap();
    let stages:Vec<_>=[&o0,&o1].into_iter().map(|oracle|NewtonStage {oracle,tolerance:1e-8,
        options:AugmentedNewtonOptions::default(),krylov:KrylovOptions::default()}).collect();
    let mut s=NewtonContinuation::new(&stages,&[0.0;3],None).unwrap();
    let mut fail=|stage:usize,_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]| {
        if stage==1 {panic!("second-stage inverse fault");} inverse(r,m,c,&d,false)
    };
    assert!(matches!(s.run_preconditioned(&mut fail,10,budget(),10,None),
        Err(ContinuationError::Stage {stage:1,source:ReverseNewtonError::Preconditioner(_)})));
    assert_eq!(s.active_stage(),1);assert_eq!(s.work().preconditioner_evals,2);
    assert_eq!(s.completed_stages().len(),1);
    assert_eq!(s.active().unwrap().inner().snapshot().x,s.completed_stages()[0].x);
    assert_eq!(s.active().unwrap().inner().radius(),1.0);
    let mut good=|_:usize,_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]|inverse(r,m,c,&d,false);
    let denied=s.run_preconditioned(&mut good,10,budget(),2,None).unwrap();
    assert_eq!(denied.stop,ContinuationStop::StageStopped(AugmentedNewtonStop::PreconditionerBudget));
    let end=s.run_preconditioned(&mut good,10,budget(),3,None).unwrap();
    assert_eq!(end.stop,ContinuationStop::Complete);assert_eq!(end.work.preconditioner_evals,3);
}

#[test]
fn g3_unpreconditioned_interlude_cannot_reset_retired_inverse_work() {
    let d=[1.0,16.0,256.0];
    let problems:Vec<_>=(0..3).map(|i|quadratic(&d,&[0.0,0.01*(i+1) as f64,0.02*(i+1) as f64],false)).collect();
    let oracles:Vec<_>=problems.iter().map(|p|ReverseProblem::new(p,limits()).unwrap()).collect();
    let stages:Vec<_>=oracles.iter().map(|oracle|NewtonStage {oracle,tolerance:1e-8,
        options:AugmentedNewtonOptions::default(),krylov:KrylovOptions::default()}).collect();
    let mut s=NewtonContinuation::new(&stages,&[0.0;3],None).unwrap();
    let mut inverse_map=|_:usize,_:&[f64],r:&[f64],m:&AugmentedLagrangian<'_>,c:&[f64]|inverse(r,m,c,&d,false);
    s.run_preconditioned(&mut inverse_map,1,budget(),1,None).unwrap();
    s.run(1,budget(),None).unwrap();
    assert_eq!(s.work().preconditioner_evals,1);
    let stopped=s.run_preconditioned(&mut inverse_map,10,budget(),1,None).unwrap();
    assert_eq!(stopped.stage,2);
    assert_eq!(stopped.stop,ContinuationStop::StageStopped(AugmentedNewtonStop::PreconditionerBudget));
    let end=s.run_preconditioned(&mut inverse_map,10,budget(),2,None).unwrap();
    assert_eq!(end.stop,ContinuationStop::Complete);assert_eq!(end.work.preconditioner_evals,2);
}
