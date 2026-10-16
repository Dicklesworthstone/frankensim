use super::*;

fn control() -> EnsembleControl { EnsembleControl::new(100_000, 1_000_000) }
fn obs(id: u64, value: f64, sigma: f64) -> EnsembleObservation {
    EnsembleObservation { id, time: 0.0, value, sigma }
}
fn moments(e: &Ensemble) -> (Vec<f64>, Vec<f64>) {
    let n = e.dimension(); let m = e.member_count();
    let center: Vec<f64> = (0..n).map(|i| e.values.chunks_exact(n).map(|x| x[i]).sum::<f64>() / m as f64).collect();
    let mut covariance = vec![0.0; n*n];
    for x in e.values.chunks_exact(n) {
        for i in 0..n { for j in 0..n { covariance[i*n+j] += (x[i]-center[i])*(x[j]-center[j])/(m-1) as f64; } }
    }
    (center, covariance)
}
fn close(a: f64, b: f64) { assert!((a-b).abs() <= 2e-12*(1.0+a.abs()+b.abs()), "{a} != {b}"); }

#[test]
fn scalar_square_root_update_retains_the_kalman_variance() {
    let mut e = Ensemble::new(0.0, 1, &[-1.0, 1.0], 2).unwrap();
    let r = e.assimilate_scalar(obs(0, 1.0, 1.0), None, &mut |_,_,x| Ok(x[0]), &mut control(), &mut || false).unwrap();
    let (mean, p) = moments(&e);
    close(mean[0], 2.0/3.0); close(p[0], 2.0/3.0);
    close(r.square_root_factor, 1.0/(1.0+1.0/3.0_f64.sqrt()));
    assert_eq!(r.updated_components, 1);
    // Applying the full gain to each member would give variance 2/9, not 2/3.
    assert!((p[0]-2.0/9.0).abs() > 0.4);
}

#[test]
fn linear_analysis_matches_independent_dense_kalman_mean_and_covariance() {
    let mut e = Ensemble::new(0.0, 3,
        &[1.0,2.0,0.0, -1.0,1.0,3.0, 2.0,-2.0,1.0, 0.0,3.0,-2.0, 4.0,0.0,2.0], 15).unwrap();
    let h = [0.7,-0.4,1.2]; let (mu,p) = moments(&e);
    let ph: Vec<f64> = (0..3).map(|i| (0..3).map(|j| p[i*3+j]*h[j]).sum()).collect();
    let s = 0.64 + h.iter().zip(&ph).map(|(a,b)| a*b).sum::<f64>();
    let innovation = 0.3-h.iter().zip(&mu).map(|(a,b)| a*b).sum::<f64>();
    e.assimilate_scalar(obs(0,0.3,0.8), None,
        &mut |_,_,x| Ok(x.iter().zip(h).map(|(a,b)| a*b).sum()), &mut control(), &mut || false).unwrap();
    let (got,q) = moments(&e);
    for i in 0..3 {
        close(got[i], mu[i]+ph[i]*innovation/s);
        for j in 0..3 { close(q[i*3+j], p[i*3+j]-ph[i]*ph[j]/s); }
    }
}

#[test]
fn sequential_independent_readings_match_scalar_information_form() {
    let mut e = Ensemble::new(0.0,1,&[-1.0,1.0],2).unwrap(); let mut work = control();
    for (id,y) in [1.0,2.0,-0.5].into_iter().enumerate() {
        e.assimilate_scalar(obs(id as u64,y,0.5),None,&mut |_,_,x| Ok(x[0]),&mut work,&mut || false).unwrap();
    }
    let (mu,p) = moments(&e);
    close(p[0], 1.0/(0.5+12.0)); close(mu[0], (4.0*(1.0+2.0-0.5))/(0.5+12.0));
    assert_eq!(work.model_calls(),6);
}

#[test]
fn absent_spread_and_zero_localization_preserve_original_bits() {
    let mut constant = Ensemble::new(0.0,2,&[0.1,2.0,0.1,-3.0,0.1,0.4],6).unwrap();
    let before = constant.values().to_vec();
    let report = constant.assimilate_scalar(obs(4,8.0,1.0),None,&mut |_,_,x| Ok(x[0]),&mut control(),&mut || false).unwrap();
    assert_eq!(report.updated_components,0); assert_eq!(constant.values(),before); assert_eq!(constant.last_observation(),Some(4));
    let mut e = Ensemble::new(0.0,2,&[-1.0,-2.0,1.0,2.0],4).unwrap();
    e.assimilate_scalar(obs(0,1.0,1.0),Some(&[0.0,1.0]),&mut |_,_,x| Ok(x[0]),&mut control(),&mut || false).unwrap();
    assert_eq!(e.values()[0].to_bits(),(-1.0_f64).to_bits()); assert_eq!(e.values()[2].to_bits(),1.0_f64.to_bits());
    assert_ne!(e.values()[1],-2.0);
}

#[test]
fn larger_than_dense_state_cap_uses_only_member_linear_workspace() {
    let n=513; let m=8; let values: Vec<f64> = (0..m).flat_map(|j| (0..n).map(move |i| 290.0+j as f64*0.1+i as f64*0.001)).collect();
    let mut e=Ensemble::new(0.0,n,&values,n*m).unwrap(); let before=e.clone();
    let mut short=EnsembleControl::new(m,n*m+m-1);
    assert!(matches!(e.assimilate_scalar(obs(0,293.0,0.2),None,&mut |_,_,x| Ok(x[256]),&mut short,&mut || false),Err(EnsembleError::WorkspaceLimit{..})));
    assert_eq!(e,before);assert_eq!(short.model_calls(),0);
    short.extend(m,n*m+m).unwrap();
    let r=e.assimilate_scalar(obs(0,293.0,0.2),None,&mut |_,_,x| Ok(x[256]),&mut short,&mut || false).unwrap();
    assert_eq!(r.updated_components,n);assert_eq!(short.model_calls(),m);
    assert!(e.values().iter().all(|x|x.is_finite()));
}

#[test]
fn every_cancellation_boundary_rolls_back_and_spent_calls_are_not_refunded() {
    let initial=Ensemble::new(0.0,2,&[-1.0,-2.0,1.0,2.0],4).unwrap();
    let mut full=initial.clone();let mut polls=0;
    full.assimilate_scalar(obs(0,1.0,1.0),None,&mut |_,_,x| Ok(x[0]),&mut control(),&mut || {polls+=1;false}).unwrap();
    for stop in 1..=polls {
        let mut e=initial.clone();let mut work=control();let mut count=0;
        assert_eq!(e.assimilate_scalar(obs(0,1.0,1.0),None,&mut |_,_,x| Ok(x[0]),&mut work,&mut || {count+=1;count==stop}),Err(EnsembleError::Cancelled));
        assert_eq!(e,initial);let spent=work.model_calls();
        e.assimilate_scalar(obs(0,1.0,1.0),None,&mut |_,_,x| Ok(x[0]),&mut work,&mut || false).unwrap();
        assert_eq!(e,full);assert_eq!(work.model_calls(),spent+2);
    }
}

#[test]
fn failures_and_reapplied_observations_do_not_replace_members() {
    let mut e=Ensemble::new(0.0,1,&[-1.0,1.0],2).unwrap();let before=e.clone();let mut work=control();
    assert!(matches!(e.assimilate_scalar(obs(0,1.0,1.0),None,&mut |i,_,x| if i==1 {Err("failed".into())} else {Ok(x[0])},&mut work,&mut || false),Err(EnsembleError::Model{member:1,..})));
    assert_eq!(e,before);assert_eq!(work.model_calls(),2);
    assert!(matches!(e.assimilate_scalar(obs(0,1.0,1.0),None,&mut |_,_,_| Ok(f64::NAN),&mut work,&mut || false),Err(EnsembleError::NonFinite(_))));
    assert_eq!(e,before); assert_eq!(work.model_calls(),3);
    e.assimilate_scalar(obs(3,1.0,1.0),None,&mut |_,_,x| Ok(x[0]),&mut work,&mut || false).unwrap();let accepted=e.clone();
    for id in [2,3] {assert_eq!(e.assimilate_scalar(obs(id,1.0,1.0),None,&mut |_,_,_| panic!("stale observation ran"),&mut work,&mut || false),Err(EnsembleError::ObservationOrder));}
    assert_eq!(e,accepted);assert_eq!(work.model_calls(),5);
}

#[test]
fn invalid_shapes_time_noise_and_limits_refuse_before_prediction() {
    for (dim,values) in [(0,vec![1.0,2.0]),(1,vec![1.0]),(2,vec![1.0,2.0,3.0]),(1,vec![f64::NAN,1.0])] {
        assert!(Ensemble::new(0.0,dim,&values,10).is_err());
    }
    let mut e=Ensemble::new(0.0,1,&[-1.0,1.0],2).unwrap();let before=e.clone();
    for bad in [obs(0,0.0,0.0),obs(0,f64::NAN,1.0),EnsembleObservation{time:1.0,..obs(0,0.0,1.0)}] {
        assert!(matches!(e.assimilate_scalar(bad,None,&mut |_,_,_| panic!("invalid input ran"),&mut control(),&mut || false),Err(EnsembleError::Invalid(_))));
    }
    let mut small=EnsembleControl::new(1,100);
    assert_eq!(e.assimilate_scalar(obs(0,0.0,1.0),None,&mut |_,_,_| panic!("underbudget ran"),&mut small,&mut || false),Err(EnsembleError::ModelCallLimit));
    assert_eq!(e,before);assert_eq!(small.model_calls(),0);
}
