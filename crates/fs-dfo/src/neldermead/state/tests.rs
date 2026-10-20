use super::*;
use std::cell::Cell;

fn config(n: usize, cap: usize) -> NelderMeadConfig {
    NelderMeadConfig { max_evaluations:cap, max_simplex_entries:n*(n+1), target:None,
        x_tolerance:1e-9, f_tolerance:1e-12, bounds:None }
}
fn run(state: &mut NelderMeadState, cap: usize, f: &mut impl FnMut(&[f64]) -> f64)
    -> NelderMeadReport
{
    state.advance_with(cap, &mut |x| Ok::<_, Infallible>(f(x)), &mut || true).unwrap()
}

#[test]
fn g0_every_budget_boundary_is_hard_including_initialization_and_shrink() {
    for cap in 0..65 {
        let mut state = NelderMeadState::new(&[0.,0.,0.], &[1.;3], config(3, cap)).unwrap();
        let mut called = 0;
        let report = run(&mut state, usize::MAX, &mut |_| { called += 1; 1.0 });
        assert_eq!(called, cap);
        assert_eq!(report.evaluations, cap);
        assert_eq!(report.stop, NelderMeadStop::EvaluationLimit);
        assert_eq!(report.best.is_some(), cap != 0);
    }
    for cap in 1..20 {
        let mut called = 0;
        let (_, value, evaluations) = super::super::nelder_mead(
            &mut |_| { called += 1; 1.0 }, &[0.,0.,0.], 1.0, cap, f64::NEG_INFINITY);
        assert_eq!((called, evaluations, value), (cap, cap, 1.0));
    }
}

#[test]
fn g3_every_pause_preserves_the_exact_search_trace_and_pending_shrink() {
    for constant in [true, false] {
        let mut uninterrupted = NelderMeadState::new(&[-1.2,1.0], &[0.25,0.25], config(2, 240)).unwrap();
        let mut reference = Vec::new();
        let objective = |x: &[f64]| if constant {1.0} else {(1.-x[0]).powi(2)+100.*(x[1]-x[0]*x[0]).powi(2)};
        let expected = run(&mut uninterrupted, 240, &mut |x| {reference.push(x.to_vec());objective(x)});
        for chunk in [1,2,3,7,19] {
            let mut resumed = NelderMeadState::new(&[-1.2,1.0], &[0.25,0.25], config(2, 0)).unwrap();
            let mut trace = Vec::new();
            let actual = loop {
                resumed.extend_evaluation_budget((resumed.evaluation_limit()+chunk).min(240)).unwrap();
                let report = run(&mut resumed, usize::MAX, &mut |x| {trace.push(x.to_vec());objective(x)});
                if report.stop != NelderMeadStop::EvaluationLimit || report.evaluations == 240 { break report; }
            };
            assert_eq!(trace, reference, "chunk={chunk}, constant={constant}");
            assert_eq!(actual, expected);
            assert_eq!(resumed.points, uninterrupted.points);
            assert_eq!(resumed.values, uninterrupted.values);
        }
    }
}

#[test]
fn g0_shrink_is_atomic_and_cancelled_completed_callback_is_not_repeated() {
    let mut state = NelderMeadState::new(&[0.,0.], &[1.,1.], config(2,20)).unwrap();
    // 3 initial + reflection + contraction + first shrunk vertex.
    let before = state.points.clone();
    let paused = run(&mut state,6,&mut |_|1.);
    assert_eq!(paused.pending_shrink_vertices,1);
    assert_eq!(state.points,before);
    let cancel = Cell::new(false);
    let mut calls=0;
    let error = state.advance_with(1,&mut |_| {
        calls+=1; cancel.set(true); Ok::<_, Infallible>(1.)
    }, &mut || !cancel.get()).unwrap_err();
    assert!(matches!(error,NelderMeadRunError::Solver(NelderMeadError::Cancelled)));
    assert_eq!(calls,1);
    assert_eq!(state.evaluations(),7);
    assert_eq!(state.points,before);
    let report=run(&mut state,0,&mut |_| panic!("completed callback must not repeat"));
    assert_eq!(report.stop,NelderMeadStop::Paused);
    assert_eq!(report.iterations,1);
    assert_eq!(report.pending_shrink_vertices,0);
    assert_eq!(state.points,vec![vec![0.,0.],vec![0.5,0.],vec![0.,0.5]]);
}

#[test]
fn g0_errors_spend_an_attempt_but_keep_the_candidate_and_original_error() {
    let mut state=NelderMeadState::new(&[1.],&[0.5],config(1,4)).unwrap();
    let error=state.advance_with(4,&mut |_| Err::<f64,_>("material validity"),&mut ||true).unwrap_err();
    assert!(matches!(error,NelderMeadRunError::Evaluation("material validity")));
    assert_eq!(state.evaluations(),1);
    assert!(state.best().is_none());
    let error=state.advance_with(1,&mut |_|Ok::<_,Infallible>(f64::NAN),&mut ||true).unwrap_err();
    assert!(matches!(error,NelderMeadRunError::Solver(NelderMeadError::NonFinite(_))));
    assert_eq!(state.evaluations(),2);
    let mut retried=Vec::new();
    run(&mut state,1,&mut |x|{retried.push(x.to_vec());x[0]*x[0]});
    assert_eq!(retried,vec![vec![1.]]);
    assert_eq!(state.evaluations(),3);
    assert_eq!(state.initialized,1);
    let error=state.advance_with(1,&mut |_|Ok::<_,Infallible>(0.),&mut ||false).unwrap_err();
    assert!(matches!(error,NelderMeadRunError::Solver(NelderMeadError::Cancelled)));
    assert_eq!(state.evaluations(),3);
}

#[test]
fn g1_box_search_activates_from_zero_and_reports_real_boundary_values() {
    for optimum in [[2.,5.],[0.,3.]] {
        let mut cfg=config(2,2000);
        cfg.bounds=Some(vec![[0.,8.],[0.,8.]]);
        let mut state=NelderMeadState::new(&[0.,0.],&[1.,1.],cfg).unwrap();
        let result=run(&mut state,2000,&mut |x| {
            assert!(x.iter().all(|&v|(0.0..=8.0).contains(&v)));
            (x[0]-optimum[0]).powi(2)+(x[1]-optimum[1]).powi(2)
        });
        assert_eq!(result.stop,NelderMeadStop::SimplexSmall);
        let best=result.best.unwrap();
        assert!(best.value<1e-15,"{best:?}");
        assert!(best.coordinates.iter().zip(optimum).all(|(&a,b)|(a-b).abs()<1e-7));
    }
    let mut cfg=config(1,20);cfg.target=Some(0.);
    let mut state=NelderMeadState::new(&[2.],&[1.],cfg).unwrap();
    let report=run(&mut state,20,&mut |x|(x[0]-2.).powi(2));
    assert_eq!(report.stop,NelderMeadStop::TargetReached);
    assert_eq!(report.evaluations,1);
}

#[test]
fn g1_rosenbrock_and_nonsmooth_absolute_value_are_solved_without_derivatives() {
    let mut state=NelderMeadState::new(&[-1.2,1.],&[0.25,0.25],config(2,2000)).unwrap();
    let report=run(&mut state,2000,&mut |x|(1.-x[0]).powi(2)+100.*(x[1]-x[0]*x[0]).powi(2));
    assert_eq!(report.stop,NelderMeadStop::SimplexSmall);
    assert!(report.best.unwrap().value<1e-15);
    let mut state=NelderMeadState::new(&[-1.,2.],&[0.5,0.5],config(2,2000)).unwrap();
    let report=run(&mut state,2000,&mut |x|(x[0]-0.3).abs()+2.*(x[1]+0.4).abs());
    assert_eq!(report.stop,NelderMeadStop::SimplexSmall);
    assert!(report.best.unwrap().value<1e-11);
}

#[test]
fn g0_malformed_simplex_and_unrepresentable_candidates_refuse() {
    assert!(NelderMeadState::new(&[],&[],config(0,10)).is_err());
    assert!(NelderMeadState::new(&[0.],&[],config(1,10)).is_err());
    for step in [0.,f64::INFINITY,f64::NAN] {
        assert!(NelderMeadState::new(&[0.],&[step],config(1,10)).is_err());
    }
    assert!(NelderMeadState::new(&[1e300],&[1.],config(1,10)).is_err());
    let mut cfg=config(1,10);cfg.max_simplex_entries=1;
    assert!(NelderMeadState::new(&[0.],&[1.],cfg).is_err());
    let mut cfg=config(1,10);cfg.bounds=Some(vec![[0.,1.]]);
    assert!(NelderMeadState::new(&[1.],&[0.5],cfg.clone()).is_err());
    assert!(NelderMeadState::new(&[1.],&[-0.5],cfg).is_ok());
    let mut cfg=config(1,10);cfg.target=Some(f64::INFINITY);
    assert!(NelderMeadState::new(&[0.],&[1.],cfg).is_err());
    let mut state=NelderMeadState::new(&[-1e308],&[1e308],config(1,10)).unwrap();
    let mut called=0;
    let err=state.advance_with(10,&mut |x| {called+=1;Ok::<_,Infallible>(-x[0]/1e308)},&mut ||true);
    assert!(matches!(err,Err(NelderMeadRunError::Solver(NelderMeadError::NonFinite(_)))));
    assert_eq!(called,3); // Expansion to 2e308 refuses before a fourth callback.
    assert_eq!(state.best().unwrap().value,-1.);
    assert!(state.extend_evaluation_budget(10).is_err());
}
