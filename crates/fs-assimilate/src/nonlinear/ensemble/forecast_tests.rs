use super::*;
use super::super::EnsembleObservation;

fn work() -> EnsembleControl { EnsembleControl::new(1000,10000) }
fn fixture() -> Ensemble { Ensemble::new(0.0,2,&[1.0,2.0,3.0,4.0,5.0,6.0],6).unwrap() }
fn advance(_:usize,start:f64,end:f64,x:&[f64],out:&mut[f64],check:&mut dyn FnMut()->bool) -> Result<(),String> {
    if check() { return Err("cancelled inside model".into()); }
    for (o,x) in out.iter_mut().zip(x) { *o = *x + end-start; }
    Ok(())
}
#[test]
fn split_forecasts_and_cloned_checkpoints_match_uninterrupted_results() {
    let e=fixture();let mut w=work();let mut a=e.forecast(2.0,&w,&mut || false).unwrap();
    a.advance(usize::MAX,&mut advance,&mut w,&mut || false).unwrap();let expected=a.finish().unwrap();
    let mut b=e.forecast(2.0,&w,&mut || false).unwrap();
    assert_eq!(b.advance(1,&mut advance,&mut w,&mut || false).unwrap().status,ForecastStatus::MemberLimit);
    let mut b=b.clone();b.advance(10,&mut advance,&mut w,&mut || false).unwrap();
    assert_eq!(b.finish().unwrap(),expected);assert_eq!(w.model_calls(),6);assert_eq!(e,fixture());
    assert_eq!(expected.values(),&[3.0,4.0,5.0,6.0,7.0,8.0]);assert_eq!(expected.time(),2.0);
}
#[test]
fn failed_member_keeps_completed_prefix_and_never_publishes_mixed_times() {
    let e=fixture();let mut w=work();let mut job=e.forecast(1.0,&w,&mut || false).unwrap();
    let mut bad=|id,start,end,x:&[f64],out:&mut[f64],c:&mut dyn FnMut()->bool| {
        if id==1 {out[0]=999.0;return Err("failure".into());} advance(id,start,end,x,out,c)
    };
    assert!(matches!(job.advance(3,&mut bad,&mut w,&mut || false),Err(EnsembleError::Model{member:1,..})));
    assert_eq!(job.completed_members(),1);assert_eq!(w.model_calls(),2);
    assert!(job.clone().finish().is_err());
    job.advance(3,&mut advance,&mut w,&mut || false).unwrap();assert_eq!(w.model_calls(),4);
    assert_eq!(job.finish().unwrap().values(),&[2.0,3.0,4.0,5.0,6.0,7.0]);assert_eq!(e,fixture());
}
#[test]
fn cancellation_inside_callback_is_latched_even_for_a_one_shot_request() {
    let e=fixture();let mut w=work();let mut job=e.forecast(1.0,&w,&mut || false).unwrap();let mut polls=0;
    let r=job.advance(3,&mut advance,&mut w,&mut || {polls+=1;polls==3}).unwrap();
    assert_eq!(r.status,ForecastStatus::Cancelled);assert_eq!(job.completed_members(),0);assert_eq!(w.model_calls(),1);
    job.advance(3,&mut advance,&mut w,&mut || false).unwrap();assert_eq!(w.model_calls(),4);
    assert_eq!(job.finish().unwrap().values(),&[2.0,3.0,4.0,5.0,6.0,7.0]);
}
#[test]
fn missing_output_and_exhausted_budgets_do_not_advance_the_member_cursor() {
    let e=fixture();let mut w=EnsembleControl::new(1,14);let mut job=e.forecast(1.0,&w,&mut || false).unwrap();
    assert_eq!(job.advance(2,&mut advance,&mut w,&mut || false),Err(EnsembleError::ModelCallLimit));assert_eq!(w.model_calls(),0);
    assert!(matches!(job.advance(1,&mut |_,_,_,_,out:&mut[f64],_| {out[0]=1.0;Ok(())},&mut w,&mut || false),Err(EnsembleError::NonFinite(_))));
    assert_eq!(job.completed_members(),0);assert_eq!(w.model_calls(),1);
    w.extend(4,14).unwrap();job.advance(3,&mut advance,&mut w,&mut || false).unwrap();assert!(job.is_complete());
}
#[test]
fn spread_rescaling_preserves_mean_and_cannot_reuse_observation_ids() {
    let mut e=fixture();let w=work();let before=e.moments(&w,&mut || false).unwrap();
    e.rescale_spread(4.0,&w,&mut || false).unwrap();let after=e.moments(&w,&mut || false).unwrap();
    assert_eq!(after.mean,before.mean);
    for (a,b) in after.std.iter().zip(before.std) {assert!((a-2.0*b).abs()<1e-12);}
    let snapshot=e.clone();e.rescale_spread(1.0,&w,&mut || false).unwrap();assert_eq!(e,snapshot);
    let mut w=work();let o=EnsembleObservation{id:0,time:0.0,value:2.0,sigma:1.0};
    e.assimilate_scalar(o,None,&mut |_,_,x| Ok(x[0]),&mut w,&mut || false).unwrap();
    e.rescale_spread(1.1,&w,&mut || false).unwrap();
    assert_eq!(e.assimilate_scalar(o,None,&mut |_,_,x| Ok(x[0]),&mut w,&mut || false),Err(EnsembleError::ObservationOrder));
    let mut job=e.forecast(1.0,&w,&mut || false).unwrap();job.advance(3,&mut advance,&mut w,&mut || false).unwrap();
    let mut next=job.finish().unwrap();next.assimilate_scalar(EnsembleObservation{time:1.0,..o},None,&mut |_,_,x| Ok(x[0]),&mut w,&mut || false).unwrap();
}
#[test]
fn invalid_forecast_and_cancelled_inflation_leave_source_untouched() {
    let mut e=fixture();let w=work();
    for t in [0.0,-1.0,f64::NAN,f64::INFINITY] {assert!(e.forecast(t,&w,&mut || false).is_err());}
    let before=e.clone();let mut polls=0;
    assert_eq!(e.rescale_spread(2.0,&w,&mut || {polls+=1;polls==3}),Err(EnsembleError::Cancelled));assert_eq!(e,before);
    assert!(e.rescale_spread(0.0,&w,&mut || false).is_err());assert_eq!(e,before);
}
