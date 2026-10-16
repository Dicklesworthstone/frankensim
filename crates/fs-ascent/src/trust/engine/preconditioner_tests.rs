use super::*;
use std::cell::Cell;

fn quadratic(x: &[f64]) -> (f64, Vec<f64>) {
    (0.5 * (4.0*x[0]*x[0] + 2.0*x[0]*x[1] + 3.0*x[1]*x[1]),
        vec![4.0*x[0] + x[1], x[0] + 3.0*x[1]])
}
fn hessian(_: &[f64], v: &[f64]) -> Vec<f64> {
    vec![4.0*v[0] + v[1], v[0] + 3.0*v[1]]
}
fn diagonal(_: &[f64], r: &[f64]) -> Vec<f64> { vec![r[0]/4.0, r[1]/3.0] }
fn unchanged(a: &TrustRegionState, b: &TrustRegionState) {
    let bits = |x: &[f64]| x.iter().map(|v|v.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&a.x), bits(&b.x)); assert_eq!(a.f.to_bits(), b.f.to_bits());
    assert_eq!(bits(&a.g), bits(&b.g)); assert_eq!(a.delta.to_bits(), b.delta.to_bits());
    assert_eq!(bits(&a.history), bits(&b.history)); assert_eq!(a.iters, b.iters);
    assert_eq!(a.negative_curvature_hits, b.negative_curvature_hits);
    assert_eq!(a.stalled, b.stalled);
}

#[test]
fn diagonal_inverse_solves_stiff_problem_with_one_krylov_product() {
    let weights: Vec<_> = (0..16).map(|i| (2.0f64).powi(i)).collect();
    let mut fg = |x: &[f64]| {
        (0.5*x.iter().zip(&weights).map(|(x,h)|h*x*x).sum::<f64>(),
            x.iter().zip(&weights).map(|(x,h)|h*x).collect())
    };
    let mut hv = |_: &[f64], v: &[f64]| v.iter().zip(&weights).map(|(v,h)|h*v).collect();
    let mut pc = |_: &[f64], r: &[f64]| r.iter().zip(&weights).map(|(r,h)|r/h).collect();
    let mut plain = TrustRegionState::new(&[0.125;16], &mut fg);
    let mut accelerated = plain.clone();
    let rule = StopRule::GradNorm(1e-10);
    plain.run_bounded(&mut fg, &mut hv, &rule, 1, 100);
    let report = accelerated.run_preconditioned_bounded(&mut fg, &mut hv,
        &mut pc, &rule, 1, 2, 1).unwrap();
    assert_eq!(report.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
    assert_eq!(report.solution.x, vec![0.0;16]);
    assert_eq!(report.solution.hv_evals, 2); // One Krylov + final model.
    assert_eq!(report.solution.evals, 2);
    assert_eq!(accelerated.preconditioner_evals(), 1);
    assert!(plain.report().hv_evals > 10);
    assert_eq!(plain.preconditioner_evals(), 0);
}

#[test]
fn noncommuting_spd_preconditioner_reaches_the_original_stationary_point() {
    let mut state = TrustRegionState::new(&[0.2,-0.3], &mut quadratic);
    let out = state.run_preconditioned_bounded(&mut quadratic, &mut hessian,
        &mut diagonal, &StopRule::GradNorm(1e-12), 2, 3, 2).unwrap();
    assert_eq!(out.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
    assert!(out.solution.x.iter().all(|v|v.abs()<1e-15));
    assert_eq!(out.solution.hv_evals, 3); assert_eq!(state.preconditioner_evals(), 2);
}

#[test]
fn inverse_budget_preserves_state_and_restart_does_not_refund_work() {
    let mut state = TrustRegionState::new(&[0.2,-0.3], &mut quadratic);
    let start = state.clone(); let rule = StopRule::GradNorm(1e-12);
    for _ in 0..2 {
        let out = state.run_preconditioned_bounded(&mut quadratic, &mut hessian,
            &mut diagonal, &rule, 2, 100, 0);
        assert!(matches!(out, Err(PreconditionerError::Budget {used:0,limit:0})));
        unchanged(&state,&start); assert_eq!(state.hv_evals,0);
    }
    assert!(matches!(state.run_preconditioned_bounded(&mut quadratic, &mut hessian,
        &mut diagonal, &rule, 2, 100, 1),
        Err(PreconditionerError::Budget {used:1,limit:1})));
    unchanged(&state,&start); assert_eq!(state.hv_evals,1); assert_eq!(state.evals,1);
    let out = state.run_preconditioned_bounded(&mut quadratic, &mut hessian,
        &mut diagonal, &rule, 2, 100, 3).unwrap();
    assert_eq!(out.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
    assert_eq!(state.preconditioner_evals(),3); assert_eq!(state.hv_evals,4);
}

#[test]
fn no_inverse_work_after_hessian_exhaustion_convergence_or_iteration_cap() {
    let rule=StopRule::GradNorm(1e-12);
    let mut state=TrustRegionState::new(&[0.2,-0.3],&mut quadratic);
    let mut forbidden=|_:&[f64],_:&[f64]| panic!("unadmitted inverse action");
    assert_eq!(state.run_preconditioned_bounded(&mut quadratic,&mut hessian,
        &mut forbidden,&rule,2,0,100).unwrap().progress,TrustRegionProgress::HessianBudget);
    state.run_preconditioned_bounded(&mut quadratic,&mut hessian,&mut forbidden,&rule,0,100,100).unwrap();
    assert_eq!(state.preconditioner_evals(),0);
    let mut done=TrustRegionState::new(&[0.0,0.0],&mut quadratic);
    assert_eq!(done.run_preconditioned_bounded(&mut quadratic,&mut hessian,
        &mut forbidden,&rule,2,100,0).unwrap().progress,
        TrustRegionProgress::Stopped(StopReason::GradNorm));
    // Only the initial inverse can run: the next HVP is already unaffordable.
    state.run_preconditioned_bounded(&mut quadratic,&mut hessian,&mut diagonal,&rule,2,1,100).unwrap();
    assert_eq!(state.preconditioner_evals(),1); assert_eq!(state.hv_evals,1);
}

#[test]
fn inverse_faults_are_typed_atomic_and_repairable() {
    for fault in 0..5 {
        let mut state=TrustRegionState::new(&[0.2,-0.3],&mut quadratic);
        let start=state.clone(); let rule=StopRule::GradNorm(1e-12);
        let out=state.run_preconditioned_bounded(&mut quadratic,&mut hessian,
            &mut |_,r| match fault {
                0=>vec![0.0], 1=>vec![f64::NAN;r.len()],
                2=>r.iter().map(|v|-v).collect(), 3=>vec![0.0;r.len()],
                _=>panic!("injected inverse fault"),
            },&rule,2,100,100);
        match (fault,out) {
            (0,Err(PreconditionerError::Dimension{expected:2,actual:1}))
            | (1,Err(PreconditionerError::NonFinite{component:0}))
            | (2|3,Err(PreconditionerError::NonPositiveProduct))
            | (4,Err(PreconditionerError::CallbackPanicked))=>{},
            (_,other)=>panic!("wrong attribution {other:?}"),
        }
        unchanged(&state,&start); assert_eq!(state.preconditioner_evals(),1);
        assert_eq!(state.hv_evals,0); assert_eq!(state.evals,1);
        let out=state.run_preconditioned_bounded(&mut quadratic,&mut hessian,
            &mut diagonal,&rule,2,100,3).unwrap();
        assert_eq!(out.progress,TrustRegionProgress::Stopped(StopReason::GradNorm));
        assert_eq!(state.preconditioner_evals(),3);
    }
}

#[test]
fn cancellation_inside_each_inverse_action_keeps_the_previous_iterate() {
    for at in 1..=2 {
        let mut state=TrustRegionState::new(&[0.2,-0.3],&mut quadratic);
        let start=state.clone(); let calls=Cell::new(0usize);
        let mut pc=|x:&[f64],r:&[f64]| {calls.set(calls.get()+1);diagonal(x,r)};
        let out=state.run_with_preconditioner(&mut quadratic,&mut hessian,Some(&mut pc),
            &StopRule::GradNorm(1e-12),2,100,100,&mut || calls.get()==at).unwrap();
        assert_eq!(out.progress,TrustRegionProgress::Paused);
        unchanged(&state,&start); assert_eq!(state.preconditioner_evals(),at);
        assert_eq!(state.hv_evals,at-1); assert_eq!(state.evals,1);
        let out=state.run_preconditioned_bounded(&mut quadratic,&mut hessian,&mut diagonal,
            &StopRule::GradNorm(1e-12),2,100,at+2).unwrap();
        assert_eq!(out.progress,TrustRegionProgress::Stopped(StopReason::GradNorm));
    }
}

fn rosenbrock(x:&[f64])->(f64,Vec<f64>) {
    let a=x[1]-x[0]*x[0];
    (100.0*a*a+(1.0-x[0]).powi(2),vec![-400.0*x[0]*a+2.0*(x[0]-1.0),200.0*a])
}
fn rosen_hv(x:&[f64],v:&[f64])->Vec<f64> {
    vec![(1200.0*x[0]*x[0]-400.0*x[1]+2.0)*v[0]-400.0*x[0]*v[1],
        -400.0*x[0]*v[0]+200.0*v[1]]
}

#[test]
fn identity_preconditioning_keeps_legacy_bits_and_splits_replay() {
    let mut plain=TrustRegionState::new(&[-1.2,1.0],&mut rosenbrock);
    let start=plain.clone();let rule=StopRule::GradNorm(1e-12);
    let mut identity=|_:&[f64],r:&[f64]|r.to_vec();
    plain.run_bounded(&mut rosenbrock,&mut rosen_hv,&rule,12,1000);
    for cut in 0..=12 {
        let mut split=start.clone();
        split.run_preconditioned_bounded(&mut rosenbrock,&mut rosen_hv,&mut identity,&rule,cut,1000,1000).unwrap();
        split.run_preconditioned_bounded(&mut rosenbrock,&mut rosen_hv,&mut identity,&rule,12-cut,1000,1000).unwrap();
        unchanged(&split,&plain); assert_eq!(split.evals,plain.evals); assert_eq!(split.hv_evals,plain.hv_evals);
    }
    let mut pc=|x:&[f64],r:&[f64]|vec![r[0]/(1200.0*x[0]*x[0]+400.0*x[1].abs()+2.0),r[1]/200.0];
    let mut straight=start.clone();
    straight.run_preconditioned_bounded(&mut rosenbrock,&mut rosen_hv,&mut pc,&rule,12,1000,1000).unwrap();
    for cut in 0..=12 {
        let mut split=start.clone();
        split.run_preconditioned_bounded(&mut rosenbrock,&mut rosen_hv,&mut pc,&rule,cut,1000,1000).unwrap();
        split.run_preconditioned_bounded(&mut rosenbrock,&mut rosen_hv,&mut pc,&rule,12-cut,1000,1000).unwrap();
        unchanged(&split,&straight); assert_eq!(split.preconditioner_evals(),straight.preconditioner_evals());
        assert_eq!(split.evals,straight.evals); assert_eq!(split.hv_evals,straight.hv_evals);
    }
}

#[test]
fn negative_curvature_and_large_inverse_units_keep_the_euclidean_ball() {
    for scale in [1e-250,1.0,1e250] {
        let mut fg=|x:&[f64]|(-0.5*(x[0]+1.0).powi(2)-(x[1]+1.0).powi(2),
            vec![-(x[0]+1.0),-2.0*(x[1]+1.0)]);
        let mut state=TrustRegionState::new(&[0.0,0.0],&mut fg);
        let out=state.run_preconditioned_bounded(&mut fg,
            &mut |_,v|vec![-v[0],-2.0*v[1]],&mut |_,r|vec![scale*r[0],0.5*scale*r[1]],
            &StopRule::GradNorm(1e-12),1,2,1).unwrap();
        assert_eq!(out.solution.negative_curvature_hits,1);
        assert!((out.solution.x.iter().map(|v|v*v).sum::<f64>()-1.0).abs()<1e-14);
        assert!((out.solution.x[0]-out.solution.x[1]).abs()<1e-14);
        assert!(out.solution.f < -1.5); assert_eq!(state.radius(),2.0);
    }
}
