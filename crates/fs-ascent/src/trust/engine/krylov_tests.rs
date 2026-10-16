use super::*;
use std::cell::Cell;

fn fg(x: &[f64]) -> (f64, Vec<f64>) {
    (2.0*x[0]*x[0] + x[0]*x[1] + 1.5*x[1]*x[1],
        vec![4.0*x[0]+x[1], x[0]+3.0*x[1]])
}
fn hv(_: &[f64], v: &[f64]) -> Vec<f64> { vec![4.0*v[0]+v[1], v[0]+3.0*v[1]] }
fn pc(_: &[f64], r: &[f64]) -> Vec<f64> { vec![r[0]/4.0,r[1]/3.0] }
fn retained(a: &TrustRegionState, b: &TrustRegionState) {
    let bits=|x: &[f64]| x.iter().map(|x|x.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&a.x),bits(&b.x)); assert_eq!(a.f.to_bits(),b.f.to_bits());
    assert_eq!(bits(&a.g),bits(&b.g)); assert_eq!(bits(&a.history),bits(&b.history));
    assert_eq!(a.delta.to_bits(),b.delta.to_bits()); assert_eq!(a.iters,b.iters);
    assert_eq!(a.negative_curvature_hits,b.negative_curvature_hits);
    assert_eq!(a.stalled,b.stalled); assert_eq!(a.krylov_options(),b.krylov_options());
}

#[test]
fn local_cap_takes_a_verified_inexact_step_without_unused_inverse_work() {
    let mut state=TrustRegionState::new(&[0.2,-0.3],&mut fg);
    let before=state.report();
    state.set_krylov_options(KrylovOptions::new(1e-8,1).unwrap());
    let out=state.run_preconditioned_bounded(&mut fg,&mut hv,&mut pc,
        &StopRule::GradNorm(1e-12),1,2,1).unwrap();
    assert_eq!(out.progress,TrustRegionProgress::Stopped(StopReason::IterationCap));
    assert_eq!(out.solution.hv_evals,2); assert_eq!(out.solution.evals,2);
    assert_eq!(state.preconditioner_evals(),1); assert_eq!(out.solution.iters,1);
    // Independent first-PCG-step formula using the declared two-by-two matrix.
    let z=[-0.5/4.0,0.7/3.0]; let hz=[4.0*z[0]+z[1],z[0]+3.0*z[1]];
    let alpha=(-0.5*z[0]+0.7*z[1])/(z[0]*hz[0]+z[1]*hz[1]);
    for i in 0..2 { assert!((out.solution.x[i]-(before.x[i]+alpha*z[i])).abs()<1e-15); }
    assert!(out.solution.f<before.f); assert!(out.solution.grad_norm>1e-3);
    let (truth,gradient)=fg(&out.solution.x);
    assert_eq!(out.solution.f.to_bits(),truth.to_bits()); assert_eq!(state.gradient(),gradient.as_slice());
    assert_eq!(state.history(),&[before.f,truth]);
}

#[test]
fn local_cap_never_borrows_the_final_model_product_from_a_hard_budget() {
    let mut state=TrustRegionState::new(&[0.2,-0.3],&mut fg);
    state.set_krylov_options(KrylovOptions::new(1e-8,1).unwrap());
    let before=state.clone(); let rule=StopRule::GradNorm(1e-12);
    let out=state.run_preconditioned_bounded(&mut fg,&mut hv,&mut pc,&rule,10,1,1).unwrap();
    assert_eq!(out.progress,TrustRegionProgress::HessianBudget);
    retained(&state,&before); assert_eq!(state.report().evals,1);
    assert_eq!(state.report().hv_evals,1); assert_eq!(state.preconditioner_evals(),1);
    let out=state.run_preconditioned_bounded(&mut fg,&mut hv,&mut pc,&rule,1,3,2).unwrap();
    assert_eq!(out.progress,TrustRegionProgress::Stopped(StopReason::IterationCap));
    assert_eq!(out.solution.hv_evals,3); assert_eq!(out.solution.evals,2);
    assert_eq!(state.preconditioner_evals(),2); assert!(out.solution.f<before.f);
}

#[test]
fn configurable_residual_accuracy_changes_work_not_outer_convergence() {
    let mut evaluate=|x:&[f64]|(0.5*x[0]*x[0]+x[1]*x[1],vec![x[0],2.0*x[1]]);
    let mut product=|_:&[f64],v:&[f64]|vec![v[0],2.0*v[1]];
    let mut loose=TrustRegionState::new(&[0.25,0.25],&mut evaluate);
    let mut tight=loose.clone(); let rule=StopRule::GradNorm(1e-12);
    loose.set_krylov_options(KrylovOptions::new(0.5,usize::MAX).unwrap());
    let l=loose.run_bounded(&mut evaluate,&mut product,&rule,1,100);
    let t=tight.run_bounded(&mut evaluate,&mut product,&rule,1,100);
    assert_eq!(l.solution.hv_evals,2); assert_eq!(t.solution.hv_evals,3);
    assert_eq!(l.progress,TrustRegionProgress::Stopped(StopReason::IterationCap));
    assert_eq!(t.progress,TrustRegionProgress::Stopped(StopReason::GradNorm));
    assert!(l.solution.grad_norm>0.05); assert!(t.solution.grad_norm<1e-15);
    let done=loose.run_bounded(&mut evaluate,&mut product,&rule,100,1000);
    assert_eq!(done.progress,TrustRegionProgress::Stopped(StopReason::GradNorm));
}

#[test]
fn local_policy_is_retained_by_every_complete_iteration_split() {
    let mut start=TrustRegionState::new(&[2.0,-3.0],&mut fg);
    let options=KrylovOptions::new(0.1,1).unwrap();start.set_krylov_options(options);
    let rule=StopRule::GradNorm(1e-12);
    for preconditioned in [false,true] {
        let run=|state:&mut TrustRegionState,count| {
            if preconditioned {
                state.run_preconditioned_bounded(&mut fg,&mut hv,&mut pc,&rule,count,100,100).unwrap();
            } else {state.run_bounded(&mut fg,&mut hv,&rule,count,100);}
        };
        let mut straight=start.clone();run(&mut straight,12);
        for cut in 0..=12 {
            let mut split=start.clone();run(&mut split,cut);run(&mut split,12-cut);
            retained(&straight,&split);
            assert_eq!(straight.evals,split.evals); assert_eq!(straight.hv_evals,split.hv_evals);
            assert_eq!(straight.preconditioner_evals(),split.preconditioner_evals());
            assert_eq!(split.krylov_options(),options);
        }
    }
}

#[test]
fn local_cap_still_handles_negative_curvature_and_cancellation_atomically() {
    let mut concave=|x:&[f64]|(-0.5*(x[0]+1.0).powi(2),vec![-x[0]-1.0]);
    let mut state=TrustRegionState::new(&[0.0],&mut concave);
    state.set_krylov_options(KrylovOptions::new(0.5,1).unwrap());
    let out=state.run_bounded(&mut concave,&mut |_,v|vec![-v[0]],&StopRule::GradNorm(1e-12),1,2);
    assert_eq!(out.solution.negative_curvature_hits,1);assert_eq!(out.solution.x,vec![1.0]);
    assert_eq!(state.radius(),2.0);
    let mut state=TrustRegionState::new(&[0.2,-0.3],&mut fg);
    state.set_krylov_options(KrylovOptions::new(1e-8,1).unwrap());
    let before=state.clone();let calls=Cell::new(0usize);
    let out=state.run_with_preconditioner(&mut fg,&mut |x,v| {
        calls.set(calls.get()+1);hv(x,v)
    },Some(&mut pc),&StopRule::GradNorm(1e-12),1,100,100,&mut ||calls.get()==1).unwrap();
    assert_eq!(out.progress,TrustRegionProgress::Paused);retained(&state,&before);
    assert_eq!(state.hv_evals,1);assert_eq!(state.preconditioner_evals(),1);assert_eq!(state.evals,1);
}

#[test]
fn invalid_krylov_options_are_refused_before_installation() {
    let mut state=TrustRegionState::new(&[0.2,-0.3],&mut fg);let before=state.clone();
    for tolerance in [f64::NAN,f64::INFINITY,f64::NEG_INFINITY,-0.1,0.0,1.0,2.0] {
        assert_eq!(KrylovOptions::new(tolerance,1),Err(KrylovOptionsError::RelativeTolerance));
    }
    assert_eq!(KrylovOptions::new(1e-8,0),Err(KrylovOptionsError::ZeroProducts));
    assert_eq!(state.krylov_options().relative_tolerance(),1e-8);
    assert_eq!(state.krylov_options().max_products(),usize::MAX);
    state.set_krylov_options(KrylovOptions::default());retained(&state,&before);
    assert_eq!(state.evals,before.evals);assert_eq!(state.hv_evals,before.hv_evals);
}
