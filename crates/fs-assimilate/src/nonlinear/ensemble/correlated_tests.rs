use super::*;

fn prior() -> Ensemble {
    Ensemble::new(2.0, 2, &[-1.0,0.2, 0.0,-0.7, 0.5,1.3, 1.8,0.4, 2.2,2.0, -0.3,1.7], 12).unwrap()
}
fn block() -> ObservationBlock {
    ObservationBlock::new(2.0, &[4,9], &[1.1,-0.4], &[0.5,0.0,0.35,0.4], 2, 8).unwrap()
}
fn predict(_: usize, _: f64, x: &[f64], out: &mut [f64], _: &mut dyn FnMut()->bool) -> Result<(),String> {
    out[0]=x[0]+0.2*x[1]; out[1]=-0.3*x[0]+1.1*x[1]; Ok(())
}
fn moments(e: &Ensemble) -> ([f64;2], [[f64;2];2]) {
    let n=e.member_count() as f64; let mut mean=[0.0;2];
    for x in e.values().chunks(2) { for i in 0..2 {mean[i]+=x[i]/n;} }
    let mut p=[[0.0;2];2];
    for x in e.values().chunks(2) {for i in 0..2 {for j in 0..2 {p[i][j]+=(x[i]-mean[i])*(x[j]-mean[j])/(n-1.0);}}}
    (mean,p)
}
// Independent observation-space batch Kalman calculation from sample moments.
// The numerical method under test never constructs these dense covariances.
fn oracle(e: &Ensemble, data: &ObservationBlock, h: impl Fn(&[f64])->[f64;2]) -> ([f64;2],[[f64;2];2]) {
    let (xmean,p)=moments(e);let count=e.member_count() as f64;
    let y=e.values().chunks(2).map(h).collect::<Vec<_>>();let mut ymean=[0.0;2];
    for v in &y {for i in 0..2 {ymean[i]+=v[i]/count;}}
    let mut cross=[[0.0;2];2];let mut s=[[0.0;2];2];
    for (x,y) in e.values().chunks(2).zip(y) {for i in 0..2 {for j in 0..2 {
        cross[i][j]+=(x[i]-xmean[i])*(y[j]-ymean[j])/(count-1.0);
        s[i][j]+=(y[i]-ymean[i])*(y[j]-ymean[j])/(count-1.0);
    }}}
    for i in 0..2 {for j in 0..2 {for k in 0..2 {s[i][j]+=data.lower[2*i+k]*data.lower[2*j+k];}}}
    let det=s[0][0]*s[1][1]-s[0][1]*s[1][0];
    let inv=[[s[1][1]/det,-s[0][1]/det],[-s[1][0]/det,s[0][0]/det]];
    let mut gain=[[0.0;2];2];let mut mean=xmean;let mut cov=p;
    for i in 0..2 {for j in 0..2 {for k in 0..2 {gain[i][j]+=cross[i][k]*inv[k][j];}}}
    for i in 0..2 {for j in 0..2 {
        mean[i]+=gain[i][j]*(data.values[j]-ymean[j]);
        for k in 0..2 {cov[i][j]-=gain[i][k]*cross[j][k];}
    }}
    (mean,cov)
}
fn assert_moments(e: &Ensemble, expected: ([f64;2],[[f64;2];2])) {
    let (mean,cov)=moments(e);
    for i in 0..2 {
        assert!((mean[i]-expected.0[i]).abs()<2e-12,"mean {mean:?} != {:?}",expected.0);
        for j in 0..2 {assert!((cov[i][j]-expected.1[i][j]).abs()<2e-12,"cov {cov:?} != {:?}",expected.1);}
    }
}

#[test]
fn correlated_linear_batch_matches_dense_kalman_moments() {
    let mut e=prior();let data=block();let expected=oracle(&e,&data,|x|[x[0]+0.2*x[1],-0.3*x[0]+1.1*x[1]]);
    let mut calls=EnsembleControl::new(6,1000);
    let report=e.assimilate_correlated(&data,&mut predict,&mut calls,&mut||false).unwrap();
    assert_moments(&e,expected);assert_eq!(calls.model_calls(),6);assert_eq!(report.model_calls,6);
    assert_eq!(e.last_observation(),Some(9));assert_eq!(report.standardized_innovations.len(),2);
    assert!(report.square_root_factors.iter().all(|a|(0.5..=1.0).contains(a)));
}

#[test]
fn nonlinear_block_conditions_one_frozen_joint_ensemble() {
    let mut e=prior();let data=block();let expected=oracle(&e,&data,|x|[x[0]*x[0]+0.2*x[1],x[1]*x[1]-0.4*x[0]]);
    let mut calls=EnsembleControl::new(6,1000);
    let report=e.assimilate_correlated(&data,&mut|_,_,x,out,_| {
        out[0]=x[0]*x[0]+0.2*x[1];out[1]=x[1]*x[1]-0.4*x[0];Ok(())
    },&mut calls,&mut||false).unwrap();
    assert_moments(&e,expected);assert_eq!(report.model_calls,6);
}

#[test]
fn diagonal_block_matches_existing_scalar_linear_updates() {
    let data=ObservationBlock::new(2.0,&[4,9],&[1.1,-0.4],&[0.5,0.0,0.0,0.4],2,8).unwrap();
    let mut joint=prior();let mut scalar=prior();let mut calls=EnsembleControl::new(100,1000);
    joint.assimilate_correlated(&data,&mut predict,&mut calls,&mut||false).unwrap();
    for i in 0..2 {
        scalar.assimilate_scalar(EnsembleObservation {id:data.ids[i],time:2.0,value:data.values[i],sigma:data.lower[2*i+i]},None,
            &mut|_,_,x|Ok(if i==0 {x[0]+0.2*x[1]} else {-0.3*x[0]+1.1*x[1]}),&mut calls,&mut||false).unwrap();
    }
    assert_moments(&joint,moments(&scalar));
}

#[test]
fn consistent_sensor_units_and_permutations_preserve_physical_moments() {
    let mut original=prior();original.assimilate_correlated(&block(),&mut predict,&mut EnsembleControl::new(6,1000),&mut||false).unwrap();
    let scale=[2.0_f64.powi(300),2.0_f64.powi(-300)];
    let data=ObservationBlock::new(2.0,&[4,9],&[1.1*scale[0],-0.4*scale[1]],
        &[0.5*scale[0],0.0,0.35*scale[1],0.4*scale[1]],2,8).unwrap();
    let mut scaled=prior();scaled.assimilate_correlated(&data,&mut|i,t,x,out,c|{
        predict(i,t,x,out,c)?;for j in 0..2 {out[j]*=scale[j];}Ok(())
    },&mut EnsembleControl::new(6,1000),&mut||false).unwrap();
    assert_moments(&scaled,moments(&original));
    // Explicitly permute R and recompute its 2x2 root for the oracle fixture.
    let a=(0.35_f64*0.35+0.4*0.4).sqrt();let b=0.175/a;let c=(0.25-b*b).sqrt();
    let swapped=ObservationBlock::new(2.0,&[4,9],&[-0.4,1.1],&[a,0.0,b,c],2,8).unwrap();
    let mut e=prior();e.assimilate_correlated(&swapped,&mut|i,t,x,out,c|{predict(i,t,x,out,c)?;out.swap(0,1);Ok(())},
        &mut EnsembleControl::new(6,1000),&mut||false).unwrap();
    assert_moments(&e,moments(&original));
}

#[test]
fn every_cancellation_boundary_rolls_back_the_whole_block() {
    let data=block();let before=prior();let mut expected=before.clone();let mut polls=0;
    expected.assimilate_correlated(&data,&mut predict,&mut EnsembleControl::new(6,1000),&mut||{polls+=1;false}).unwrap();
    for stop in 1..=polls {
        let mut e=before.clone();let mut control=EnsembleControl::new(12,1000);let mut seen=0;
        assert_eq!(e.assimilate_correlated(&data,&mut predict,&mut control,&mut||{seen+=1;seen==stop}),Err(EnsembleError::Cancelled));
        assert_eq!(e,before);let spent=control.model_calls();
        e.assimilate_correlated(&data,&mut predict,&mut control,&mut||false).unwrap();
        assert_eq!(e,expected);assert_eq!(control.model_calls(),spent+6);
    }
}

#[test]
fn model_cancellation_failures_and_unwritten_outputs_are_atomic() {
    let data=block();let before=prior();
    for mode in 0..3 {
        let mut e=before.clone();let mut control=EnsembleControl::new(12,1000);
        let error=e.assimilate_correlated(&data,&mut|i,t,x,out,c|{
            if i==2 {
                if mode==0 {return Err("broken sensor".into());}
                if mode==1 {out[0]=1.0;return Ok(());}
                out.fill(f64::INFINITY);return Ok(());
            }
            predict(i,t,x,out,c)
        },&mut control,&mut||false).unwrap_err();
        assert!(matches!(error,EnsembleError::Model{member:2,..}|EnsembleError::NonFinite(_)));
        assert_eq!(e,before);assert_eq!(control.model_calls(),3);
    }
    let mut e=before.clone();let trigger=std::cell::Cell::new(false);let mut control=EnsembleControl::new(12,1000);
    assert_eq!(e.assimilate_correlated(&data,&mut|i,t,x,out,check|{
        trigger.set(true);assert!(check());trigger.set(false);predict(i,t,x,out,check)
    },&mut control,&mut||trigger.get()),Err(EnsembleError::Cancelled));
    assert_eq!(e,before);assert_eq!(control.model_calls(),1);
}

#[test]
fn admission_and_observation_order_refuse_before_model_calls() {
    for (ids,lower) in [([4,4],[0.5,0.0,0.35,0.4]),([9,4],[0.5,0.0,0.35,0.4]),
        ([4,9],[0.5,0.1,0.35,0.4]),([4,9],[0.0,0.0,0.35,0.4]),([4,9],[0.5,0.0,f64::NAN,0.4])] {
        assert!(ObservationBlock::new(2.0,&ids,&[1.0,2.0],&lower,2,8).is_err());
    }
    assert!(ObservationBlock::new(2.0,&[4,9],&[1.0,2.0],&[0.5,0.0,0.35,0.4],1,8).is_err());
    assert!(ObservationBlock::new(2.0,&[4,9],&[1.0,2.0],&[0.5,0.0,0.35,0.4],2,7).is_err());
    let data=block();let mut e=prior();let before=e.clone();let required=e.correlated_workspace(&data).unwrap();
    let mut control=EnsembleControl::new(6,required-1);
    assert_eq!(e.assimilate_correlated(&data,&mut predict,&mut control,&mut||false),
        Err(EnsembleError::WorkspaceLimit{required,limit:required-1}));
    assert_eq!(e,before);assert_eq!(control.model_calls(),0);
    assert_eq!(e.assimilate_correlated(&data,&mut predict,&mut EnsembleControl::new(5,required),&mut||false),Err(EnsembleError::ModelCallLimit));
    control.extend(12,required).unwrap();e.assimilate_correlated(&data,&mut predict,&mut control,&mut||false).unwrap();
    let accepted=e.clone();assert_eq!(e.assimilate_correlated(&data,&mut predict,&mut control,&mut||false),Err(EnsembleError::ObservationOrder));
    assert_eq!(e,accepted);assert_eq!(control.model_calls(),6);
    let wrong_time=ObservationBlock::new(3.0,&[10,11],&[1.0,2.0],&[1.0,0.0,0.0,1.0],2,8).unwrap();
    assert!(e.assimilate_correlated(&wrong_time,&mut predict,&mut control,&mut||false).is_err());
}

#[test]
fn field_sized_block_uses_no_state_covariance_or_per_sensor_model_calls() {
    let n=513;let values=(0..8*n).map(|i|300.0+(i/n) as f64*(1.0+(i%n) as f64/n as f64)).collect::<Vec<_>>();
    let mut e=Ensemble::new(2.0,n,&values,values.len()).unwrap();let data=block();
    let required=e.correlated_workspace(&data).unwrap();assert!(required<n*n);
    let mut control=EnsembleControl::new(8,required);
    e.assimilate_correlated(&data,&mut|_,_,x,out,_|{out[0]=x[0];out[1]=x[n-1];Ok(())},&mut control,&mut||false).unwrap();
    assert_eq!(control.model_calls(),8);assert!(e.values().iter().all(|v|v.is_finite()));assert_eq!(e.dimension(),n);
}

#[test]
fn shared_reference_factor_reconstructs_declared_covariance_without_jitter() {
    let sigma=[0.2_f64,1.5,0.7,0.03];let common=0.8_f64;
    let data=ObservationBlock::shared_reference(2.0,&[1,2,3,4],&[0.0;4],&sigma,common,4,40,&mut||false).unwrap();
    for i in 0..4 {for j in 0..4 {
        let r=(0..4).map(|k|data.lower[4*i+k]*data.lower[4*j+k]).sum::<f64>();
        let expected=common*common+if i==j {sigma[i]*sigma[i]} else {0.0};
        assert!((r-expected).abs()<1e-14);
    }}
    let diagonal=ObservationBlock::shared_reference(2.0,&[1,2,3,4],&[0.0;4],&sigma,0.0,4,40,&mut||false).unwrap();
    for i in 0..4 {for j in 0..4 {assert_eq!(diagonal.lower[4*i+j],if i==j {sigma[i]} else {0.0});}}
}

#[test]
fn many_common_reference_readings_retain_the_information_ceiling() {
    let m=32;let ids=(0..m as u64).collect::<Vec<_>>();
    let data=ObservationBlock::shared_reference(0.0,&ids,&vec![3.0;m],&vec![0.5;m],2.0,m,2*m*m+2*m,&mut||false).unwrap();
    let mut e=Ensemble::new(0.0,1,&[-1.0,0.0,1.0],3).unwrap();
    e.assimilate_correlated(&data,&mut|_,_,x,out,_|{out.fill(x[0]);Ok(())},&mut EnsembleControl::new(3,1000),&mut||false).unwrap();
    let moments=e.moments(&EnsembleControl::new(0,2),&mut||false).unwrap();
    let effective=4.0+0.25/m as f64;
    assert!((moments.mean[0]-3.0/(1.0+effective)).abs()<1e-12);
    assert!((moments.std[0]*moments.std[0]-effective/(1.0+effective)).abs()<1e-12);
    assert!(moments.std[0]*moments.std[0]>0.8);
}

#[test]
fn shared_reference_rejects_invalid_scales_and_cancels_construction() {
    for bad in [f64::NAN,-1.0,f64::INFINITY] {
        assert!(ObservationBlock::shared_reference(0.0,&[1,2],&[0.0;2],&[1.0;2],bad,2,12,&mut||false).is_err());
    }
    assert!(ObservationBlock::shared_reference(0.0,&[1,2],&[0.0;2],&[1.0,0.0],0.1,2,12,&mut||false).is_err());
    assert_eq!(ObservationBlock::shared_reference(0.0,&[1,2],&[0.0;2],&[1.0;2],0.1,2,11,&mut||false),
        Err(EnsembleError::WorkspaceLimit{required:12,limit:11}));
    let mut polls=0;
    assert_eq!(ObservationBlock::shared_reference(0.0,&[1,2],&[0.0;2],&[1.0;2],0.1,2,12,&mut||{polls+=1;polls==3}),Err(EnsembleError::Cancelled));
    // Factor generation handles scales whose squares overflow or underflow.
    for scale in [1e-200_f64,1e200] {
        let b=ObservationBlock::shared_reference(0.0,&[1,2],&[0.0;2],&[scale;2],scale,2,12,&mut||false).unwrap();
        assert!(b.noise_lower().iter().all(|x|x.is_finite()));
        assert!((b.noise_lower()[0]/scale-2.0_f64.sqrt()).abs()<1e-14);
    }
}
