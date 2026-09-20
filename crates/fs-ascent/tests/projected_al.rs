//! G0/G1/G4/G5: actual linear-storage optimization, bounds and restart behavior.
use std::ops::ControlFlow;
use fs_ascent::projected_al::{ProjectedAlError,ProjectedAlOptions,ProjectedAlSample,ProjectedAlState,ProjectedAlStop};
fn quadratic(x:&[f64])->Result<Option<ProjectedAlSample>, &'static str> {
    let target:Vec<f64>=(0..x.len()).map(|i|0.7+0.01*(i%9) as f64).collect();
    let g:Vec<f64>=x.iter().zip(target).map(|(x,t)|x-t).collect();
    Ok(Some(ProjectedAlSample{objective:0.5*g.iter().map(|g|g*g).sum::<f64>(),gradient:g,
        constraint:x.iter().sum::<f64>()/x.len() as f64-0.3,constraint_gradient:vec![1.0/x.len() as f64;x.len()]}))
}
fn start(n:usize,initial:f64,options:ProjectedAlOptions)->ProjectedAlState {
    ProjectedAlState::try_new(&vec![initial;n],&vec![0.0;n],&vec![1.0;n],options,&mut quadratic,|_|ControlFlow::Continue(())).unwrap()
}
#[test]
fn g1_4096_variables_converge_without_dense_bound_or_kkt_rows() {
    for initial in [0.1,0.9] {
        let n=4096;let mut state=start(n,initial,Default::default());
        let report=state.try_run(2000,&mut quadratic,|_|ControlFlow::Continue(())).unwrap();
        assert_eq!(report.stop,ProjectedAlStop::Converged,"{report:?}");
        assert!(report.kkt.within_tolerance(1e-6));
        let mean=(0..n).map(|i|0.7+0.01*(i%9) as f64).sum::<f64>()/n as f64;
        for (i,&x) in state.point().iter().enumerate() {
            assert!((x-(0.7+0.01*(i%9) as f64-mean+0.3)).abs()<2e-6);
        }
        assert!(report.multiplier>0.0);
        assert!(report.work.evaluations<1000);
    }
}
#[test]
fn g1_nonlinear_active_constraint_and_exact_box_faces() {
    let target=[1.0,0.8,0.4,0.0,-0.3];
    let mut evaluate=|x:&[f64]|->Result<Option<ProjectedAlSample>, &'static str> {
        assert!(x.iter().all(|x|*x>=0.0&&*x<=1.0));
        let g:Vec<f64>=x.iter().zip(target).map(|(x,t)|x-t).collect();
        Ok(Some(ProjectedAlSample{objective:0.5*g.iter().map(|g|g*g).sum::<f64>(),gradient:g,
            constraint:x.iter().map(|x|x*x).sum::<f64>()/5.0-0.1,
            constraint_gradient:x.iter().map(|x|2.0*x/5.0).collect()}))
    };
    let mut state=ProjectedAlState::try_new(&[0.9;5],&[0.0;5],&[1.0;5],Default::default(),&mut evaluate,|_|ControlFlow::Continue(())).unwrap();
    let result=state.try_run(1000,&mut evaluate,|_|ControlFlow::Continue(())).unwrap();
    assert_eq!(result.stop,ProjectedAlStop::Converged);
    assert_eq!(state.point()[4],0.0);assert_eq!(state.point()[3],0.0);
    let factor=(0.5_f64/(1.0+0.64+0.16)).sqrt();
    for (i,t) in target[..3].iter().enumerate(){assert!((state.point()[i]-factor*t).abs()<2e-6);}
    assert!(result.kkt.within_tolerance(1e-6));
}
#[test]
fn g5_split_runs_preserve_spectral_steps_multipliers_and_all_work() {
    let initial=start(17,0.9,Default::default());let mut full=initial.clone();let mut split=initial;
    let a=full.try_run(1000,&mut quadratic,|_|ControlFlow::Continue(())).unwrap();
    loop {
        let r=split.try_run(3,&mut quadratic,|_|ControlFlow::Continue(())).unwrap();
        if r.stop!=ProjectedAlStop::IterationLimit {
            assert_eq!(r.stop,a.stop);assert_eq!(r.multiplier,a.multiplier);assert_eq!(r.penalty,a.penalty);break;
        }
    }
    assert_eq!(split.point(),full.point());assert_eq!(split.sample(),full.sample());assert_eq!(split.work(),full.work());
    let spent=split.work();let report=split.try_run(0,&mut quadratic,|_|ControlFlow::Continue(())).unwrap();
    assert_eq!(report.work,spent);
}
#[test]
fn g4_cancelled_trial_and_callback_errors_do_not_publish_a_candidate() {
    let mut state=start(17,0.9,Default::default());let point=state.point().to_vec();let sample=state.sample().clone();
    let stopped=state.try_run(1,&mut quadratic,|w|if w.evaluations>1{ControlFlow::Break(())}else{ControlFlow::Continue(())});
    assert!(matches!(stopped,Err(ProjectedAlError::Cancelled)));
    assert_eq!(state.point(),point);assert_eq!(state.sample(),&sample);assert_eq!(state.work().evaluations,2);
    let stopped=state.try_run(1,&mut |_|->Result<Option<ProjectedAlSample>,&'static str>{Err("physics failed")},|_|ControlFlow::Continue(()));
    assert!(matches!(stopped,Err(ProjectedAlError::Evaluation("physics failed"))));
    assert_eq!(state.point(),point);assert_eq!(state.work().evaluations,3);
    let r=state.try_run(1000,&mut quadratic,|_|ControlFlow::Continue(())).unwrap();assert_eq!(r.stop,ProjectedAlStop::Converged);
}
#[test]
fn g4_budget_unavailable_and_invalid_trials_are_not_convergence() {
    let mut state=start(17,0.9,ProjectedAlOptions{max_evaluations:1,..Default::default()});
    let result=state.try_run(100,&mut quadratic,|_|ControlFlow::Continue(())).unwrap();
    assert_eq!(result.stop,ProjectedAlStop::EvaluationLimit);assert!(!result.kkt.within_tolerance(1e-6));
    assert!(result.constraint>0.0);assert_eq!(result.work.evaluations,1);
    let mut state=start(17,0.9,Default::default());let point=state.point().to_vec();
    let result=state.try_run(1,&mut |_|->Result<Option<ProjectedAlSample>,&'static str>{Ok(None)},|_|ControlFlow::Continue(())).unwrap();
    assert_eq!(result.stop,ProjectedAlStop::Stalled);assert_eq!(state.point(),point);assert_eq!(result.work.rejected_trials,40);
    let result=state.try_run(1,&mut |_|->Result<Option<ProjectedAlSample>,&'static str>{Ok(Some(ProjectedAlSample{objective:f64::NAN,gradient:vec![0.0;17],constraint:0.0,constraint_gradient:vec![0.0;17]}))},|_|ControlFlow::Continue(()));
    assert!(matches!(result,Err(ProjectedAlError::Invalid(_))));assert_eq!(state.point(),point);
    let mut calls=0;let result=ProjectedAlState::try_new(&[0.5;5],&[0.0;5],&[1.0;5],ProjectedAlOptions{max_dimension:4,..Default::default()},
        &mut |x|{calls+=1;quadratic(x)},|_|ControlFlow::Continue(()));
    assert!(matches!(result,Err(ProjectedAlError::Invalid(_))));assert_eq!(calls,0);
}
