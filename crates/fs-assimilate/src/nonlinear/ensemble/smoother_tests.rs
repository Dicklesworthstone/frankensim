use super::*;

fn initial(capacity: usize) -> EnsembleSmoother {
    EnsembleSmoother::new(Ensemble::new(0.0, 1, &[-1.0, 1.0], 2).unwrap(), capacity).unwrap()
}
fn control() -> EnsembleControl { EnsembleControl::new(10000, 100000) }
fn double(_: usize, _: f64, _: f64, x: &[f64], out: &mut [f64], _: &mut dyn FnMut() -> bool)
    -> Result<(), String> { assert_eq!(x.len(), 1); out[0] = 2.0*x[0]; Ok(()) }
fn advance(state: EnsembleSmoother, time: f64, c: &mut EnsembleControl) -> EnsembleSmoother {
    let mut job = state.forecast(time, c, &mut || false).unwrap();
    job.advance(usize::MAX, &mut double, c, &mut || false).unwrap(); job.finish().unwrap()
}
fn close(a: f64, b: f64) { assert!((a-b).abs() < 2e-12*(1.0+b.abs()), "{a} != {b}"); }
fn obs(id: u64, time: f64, value: f64) -> EnsembleObservation { EnsembleObservation { id, time, value, sigma: 1.0 } }

#[test]
fn later_measurement_updates_past_and_present_with_kalman_cross_covariance() {
    let mut c = control(); let mut state = advance(initial(3), 1.0, &mut c);
    state.assimilate_scalar(obs(4, 1.0, 2.0), &mut |_, t, x| { assert_eq!(t,1.0); Ok(x[0]) }, &mut c, &mut || false).unwrap();
    let old = state.moments_at(0,&c,&mut||false).unwrap(); let new = state.moments_at(1,&c,&mut||false).unwrap();
    close(old.mean[0],8.0/9.0); close(old.std[0].powi(2),2.0/9.0);
    close(new.mean[0],16.0/9.0); close(new.std[0].powi(2),8.0/9.0);
    assert_eq!(c.model_calls(),4); assert_eq!(state.times(),[0.0,1.0]);
}

#[test]
fn delayed_measurement_updates_latest_state_without_reforecasting() {
    let mut c = control(); let mut state = advance(initial(3),1.0,&mut c);
    state.assimilate_scalar(obs(10,0.0,1.0),&mut|_,t,x| {assert_eq!(t,0.0);Ok(x[0])},&mut c,&mut||false).unwrap();
    close(state.moments_at(0,&c,&mut||false).unwrap().mean[0],2.0/3.0);
    close(state.moments_at(1,&c,&mut||false).unwrap().mean[0],4.0/3.0);
    assert_eq!(state.time(),1.0); assert_eq!(c.model_calls(),4);
}

#[test]
fn cross_time_correlated_block_matches_independent_information_solution() {
    let mut c=control();let mut state=advance(initial(3),1.0,&mut c);
    let block=ObservationBlock::new(1.0,&[3,7],&[1.0,2.4],&[1.0,0.0,0.5,0.7],2,8).unwrap();
    let report=state.assimilate_correlated(&block,&[0.0,1.0],&mut|_,states,out,_| {
        assert_eq!(states.len(),2);assert_eq!(states.time(0),Some(0.0));assert_eq!(states.time(1),Some(1.0));
        assert!(states.state(2).is_none());out[0]=states.state(0).unwrap()[0];out[1]=states.state(1).unwrap()[0];Ok(())
    },&mut c,&mut||false).unwrap();
    let h=1.5/0.7;let y=(2.4-0.5)/0.7;let variance=1.0/(0.5+1.0+h*h);let mean=variance*(1.0+h*y);
    let a=state.moments_at(0,&c,&mut||false).unwrap();let b=state.moments_at(1,&c,&mut||false).unwrap();
    close(a.mean[0],mean);close(a.std[0].powi(2),variance);close(b.mean[0],2.0*mean);close(b.std[0].powi(2),4.0*variance);
    assert_eq!(report.model_calls,2);assert_eq!(c.model_calls(),4);assert_eq!(state.last_observation(),Some(7));
}

#[test]
fn eviction_is_explicit_and_consumed_ids_survive_forecasts() {
    let mut c=control();let mut state=initial(2);
    state.assimilate_scalar(obs(12,0.0,0.5),&mut|_,_,x|Ok(x[0]),&mut c,&mut||false).unwrap();
    state=advance(state,1.0,&mut c);state=advance(state,2.0,&mut c);
    assert_eq!(state.times(),[1.0,2.0]);let before=state.clone();let spent=c.model_calls();
    for reading in [obs(13,0.0,1.0),obs(13,1.5,1.0),obs(13,3.0,1.0),obs(12,1.0,1.0)] {
        assert!(state.assimilate_scalar(reading,&mut|_,_,x|Ok(x[0]),&mut c,&mut||false).is_err());assert_eq!(state,before);
    }
    assert_eq!(c.model_calls(),spent);
    state.assimilate_scalar(obs(13,1.0,1.0),&mut|_,_,x|Ok(x[0]),&mut c,&mut||false).unwrap();
    assert_eq!(state.last_observation(),Some(13));
}

#[test]
fn partial_forecast_preserves_history_member_identity_and_charges_retry() {
    let mut c=control();let state=advance(initial(3),1.0,&mut c);let mut split=state.forecast(2.0,&c,&mut||false).unwrap();
    split.advance(1,&mut double,&mut c,&mut||false).unwrap();assert_eq!(split.completed_members(),1);
    let saved=split.clone();let before=c.model_calls();
    let error=split.advance(1,&mut|member,_,_,_,_,_| {assert_eq!(member,1);Err("failed".into())},&mut c,&mut||false);
    assert!(matches!(error,Err(EnsembleError::Model{member:1,..})));assert_eq!(split,saved);assert_eq!(c.model_calls(),before+1);
    split.advance(1,&mut double,&mut c,&mut||false).unwrap();let finished=split.finish().unwrap();
    let straight=advance(state.clone(),2.0,&mut control());assert_eq!(finished,straight);
    assert_eq!(finished.member_at(0,0),Some([-1.0].as_slice()));assert_eq!(finished.member_at(1,0),Some([-2.0].as_slice()));
    assert_eq!(finished.member_at(2,0),Some([-4.0].as_slice()));assert_eq!(state.times(),[0.0,1.0]);
    assert!(saved.finish().is_err());
}

#[test]
fn every_block_cancellation_boundary_is_atomic_and_retryable() {
    let mut c=control();let original=advance(initial(3),1.0,&mut c);
    let block=ObservationBlock::new(1.0,&[1,2],&[0.3,0.7],&[1.0,0.0,0.2,0.9],2,8).unwrap();
    let mut predict=|_:usize,states:ObservationStates<'_>,out:&mut[f64],_:&mut dyn FnMut()->bool| {
        out[0]=states.state(0).unwrap()[0];out[1]=states.state(1).unwrap()[0];Ok(())
    };
    let mut complete=original.clone();let mut polls=0;
    complete.assimilate_correlated(&block,&[0.0,1.0],&mut predict,&mut control(),&mut||{polls+=1;false}).unwrap();
    for stop in 1..=polls {
        let mut state=original.clone();let mut checks=0;let mut budget=control();
        assert_eq!(state.assimilate_correlated(&block,&[0.0,1.0],&mut predict,&mut budget,&mut||{checks+=1;checks==stop}),Err(EnsembleError::Cancelled));
        assert_eq!(state,original);let spent=budget.model_calls();
        state.assimilate_correlated(&block,&[0.0,1.0],&mut predict,&mut budget,&mut||false).unwrap();
        assert_eq!(state,complete);assert_eq!(budget.model_calls(),spent+2);
    }
}

#[test]
fn capacity_one_matches_existing_filter_and_forecast() {
    let mut c=control();let mut plain=Ensemble::new(0.0,1,&[-1.0,1.0],2).unwrap();
    let mut smooth=EnsembleSmoother::new(plain.clone(),1).unwrap();
    for (id,time) in [(1,1.0),(2,2.0)] {
        let mut job=plain.forecast(time,&c,&mut||false).unwrap();job.advance(2,&mut double,&mut c,&mut||false).unwrap();plain=job.finish().unwrap();
        smooth=advance(smooth,time,&mut c);
        plain.assimilate_scalar(obs(id,time,0.7),None,&mut|_,_,x|Ok(x[0]),&mut c,&mut||false).unwrap();
        smooth.assimilate_scalar(obs(id,time,0.7),&mut|_,_,x|Ok(x[0]),&mut c,&mut||false).unwrap();
        assert_eq!(plain.values(),smooth.joint.values());assert_eq!(smooth.times(),[time]);
    }
}

#[test]
fn exact_workspace_and_callback_limits_precede_model_work() {
    let state=initial(3); // target: 2 frames * 1 state * 2 members = 4 entries
    let required=2*4+2+2;
    assert!(matches!(state.forecast(1.0,&EnsembleControl::new(2,required-1),&mut||false),Err(EnsembleError::WorkspaceLimit{..})));
    let mut c=EnsembleControl::new(1,required);let mut job=state.forecast(1.0,&c,&mut||false).unwrap();
    assert_eq!(job.advance(2,&mut double,&mut c,&mut||false),Err(EnsembleError::ModelCallLimit));assert_eq!(c.model_calls(),0);
    c.extend(2,required).unwrap();job.advance(2,&mut double,&mut c,&mut||false).unwrap();let mut state=job.finish().unwrap();
    let block=ObservationBlock::new(1.0,&[1,2],&[0.0,0.0],&[1.0,0.0,0.0,1.0],2,8).unwrap();
    let required=state.joint.correlated_workspace(&block).unwrap()+2;
    let mut c=EnsembleControl::new(2,required-1);
    assert!(matches!(state.assimilate_correlated(&block,&[0.0,1.0],&mut|_,_,_,_|panic!("not admitted"),&mut c,&mut||false),Err(EnsembleError::WorkspaceLimit{..})));
    c.extend(2,required).unwrap();
    state.assimilate_correlated(&block,&[0.0,1.0],&mut|_,states,out,_| {out[0]=states.state(0).unwrap()[0];out[1]=states.state(1).unwrap()[0];Ok(())},&mut c,&mut||false).unwrap();
    assert_eq!(c.model_calls(),2);
}

#[test]
fn physical_forecast_receives_only_latest_multicomponent_field() {
    let initial=Ensemble::new(0.0,2,&[1.0,2.0,3.0,4.0],4).unwrap();let mut state=EnsembleSmoother::new(initial,3).unwrap();let mut c=control();
    for end in [1.0,2.0,3.0] {
        let mut job=state.forecast(end,&c,&mut||false).unwrap();
        job.advance(2,&mut|member,from,to,x,out,_| {
            assert_eq!(x.len(),2);assert_eq!(out.len(),2);assert_eq!(to-from,1.0);
            close(x[0],(2*member+1) as f64+from);out[0]=x[0]+1.0;out[1]=x[1]*2.0;Ok(())
        },&mut c,&mut||false).unwrap();state=job.finish().unwrap();
    }
    assert_eq!(state.times(),[1.0,2.0,3.0]);assert_eq!(state.member_at(2,1),Some([6.0,32.0].as_slice()));
    assert!(state.member_at(3,0).is_none());assert!(state.member_at(0,2).is_none());
}
