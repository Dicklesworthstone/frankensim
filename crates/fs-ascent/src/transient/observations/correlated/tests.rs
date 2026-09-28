use super::*;
use crate::transient::{TransientConfig, TransientStudy, evaluate_transient};
use super::super::SensorReading;
use fs_time::{PiController, adaptive::adjoint::trajectory::RecordingConfig};

fn limits() -> SharedNoiseLimits {
    SharedNoiseLimits {max_components:100_000,max_work:10_000_000,relative_tolerance:1e-9}
}
fn data(values: &[f64], sigma: &[f64]) -> SensorData {
    SensorData::new(&values.iter().zip(sigma).map(|(&v,&s)|
        SensorReading::new(0.0,0,v,s,SensorLoss::Quadratic).unwrap()).collect::<Vec<_>>(),1024).unwrap()
}
fn dense(readings: &SensorData, k: usize, loadings: &[f64], predictions: &[f64]) -> (f64,Vec<f64>) {
    let m=predictions.len();let mut covariance=vec![0.0;m*m];let mut residual=Vec::new();
    for i in 0..m {
        residual.push(predictions[i]-readings.readings()[i].value());
        for j in 0..m {
            covariance[i*m+j]=if i==j {readings.readings()[i].sigma().powi(2)} else {0.0};
            for source in 0..k {covariance[i*m+j]+=loadings[i*k+source]*loadings[j*k+source];}
        }
    }
    // Different matrix and factorization from the rank-space production solve.
    let mut weights=residual.clone();fs_la::factor::lu(&covariance,m).unwrap().solve(&mut weights);
    (0.5*residual.iter().zip(&weights).map(|(a,b)|a*b).sum::<f64>(),weights)
}

#[test]
fn low_rank_score_matches_dense_covariance_and_sensor_permutation() {
    let d=data(&[0.2,-0.1,0.5,0.8],&[0.1,0.3,0.2,0.15]);
    let u=[0.4,0.0, 0.2,-0.1, 0.0,0.3, -0.2,0.1];let y=[0.6,0.3,0.1,1.1];
    let noise=SharedNoise::new(d.clone(),2,&u,limits(),&mut||false).unwrap();
    let got=noise.score(&y,&mut||false).unwrap();let expected=dense(&d,2,&u,&y);
    assert!((got.value-expected.0).abs()<1e-11);
    for (a,b) in got.prediction_bar.iter().zip(&expected.1) {assert!((a-b).abs()<1e-10);}
    assert!(got.relative_stationarity<=limits().relative_tolerance);
    let permutation=[3,1,0,2];
    let rows:Vec<_>=permutation.iter().map(|&i|d.readings()[i].clone()).collect();
    let loadings:Vec<_>=permutation.iter().flat_map(|&i|[u[2*i],u[2*i+1]]).collect();
    let predictions:Vec<_>=permutation.iter().map(|&i|y[i]).collect();
    let permuted=SharedNoise::new(SensorData::new(&rows,4).unwrap(),2,&loadings,limits(),&mut||false).unwrap();
    let other=permuted.score(&predictions,&mut||false).unwrap();
    assert!((other.value-got.value).abs()<1e-11);
    for (j,&i) in permutation.iter().enumerate() {assert!((other.prediction_bar[j]-got.prediction_bar[i]).abs()<1e-10);}
}

#[test]
fn a_shared_reference_error_does_not_average_away_with_more_readings() {
    for m in [1,8,64] {
        let d=data(&vec![0.0;m],&vec![0.1;m]);
        let noise=SharedNoise::new(d,1,&vec![0.4;m],limits(),&mut||false).unwrap();
        let result=noise.score(&vec![1.0;m],&mut||false).unwrap();
        let precision=m as f64/(0.01+m as f64*0.16);
        assert!((result.value-0.5*precision).abs()<1e-10);
        assert!((result.prediction_bar.iter().sum::<f64>()-precision).abs()<1e-9);
        assert!(precision<1.0/0.16);
        if m>1 {assert!(precision<(m as f64/0.17)*0.6);}
    }
}

#[test]
fn consistent_reading_unit_changes_preserve_cost_and_rescale_gradient() {
    let values=[0.2,-0.1,0.5];let sigmas=[0.1,0.3,0.2];let loadings=[0.4,0.2,-0.1];let y=[0.6,0.3,0.1];
    let original=SharedNoise::new(data(&values,&sigmas),1,&loadings,limits(),&mut||false).unwrap().score(&y,&mut||false).unwrap();
    let scales=[2.0f64.powi(-40),2.0f64.powi(20),2.0f64.powi(-10)];
    let v:Vec<_>=values.iter().zip(scales).map(|(v,s)|v*s).collect();
    let s:Vec<_>=sigmas.iter().zip(scales).map(|(v,s)|v*s).collect();
    let u:Vec<_>=loadings.iter().zip(scales).map(|(v,s)|v*s).collect();
    let p:Vec<_>=y.iter().zip(scales).map(|(v,s)|v*s).collect();
    let scaled=SharedNoise::new(data(&v,&s),1,&u,limits(),&mut||false).unwrap().score(&p,&mut||false).unwrap();
    assert_eq!(scaled.value.to_bits(),original.value.to_bits());
    for i in 0..3 {assert!((scaled.prediction_bar[i]*scales[i]-original.prediction_bar[i]).abs()<1e-10);}
}

#[test]
fn malformed_noise_huber_and_resource_exhaustion_refuse() {
    let d=data(&[0.0,0.0],&[0.1,0.1]);
    assert!(SharedNoise::new(d.clone(),1,&[0.2],limits(),&mut||false).is_err());
    assert!(SharedNoise::new(d.clone(),1,&[0.2,f64::NAN],limits(),&mut||false).is_err());
    assert!(SharedNoise::new(d.clone(),1,&[1e300,1e300],limits(),&mut||false).is_err());
    assert!(SharedNoise::new(d.clone(),1,&[0.2,0.2],SharedNoiseLimits {max_work:0,..limits()},&mut||false).is_err());
    assert!(SharedNoise::new(d,1,&[0.2,0.2],SharedNoiseLimits {max_components:0,..limits()},&mut||false).is_err());
    let huber=SensorData::new(&[SensorReading::new(0.0,0,0.0,0.1,SensorLoss::Huber {threshold:1.5}).unwrap()],1).unwrap();
    assert!(matches!(SharedNoise::new(huber,1,&[0.2],limits(),&mut||false),Err(SharedNoiseError::Invalid(_))));
}

struct Family;
struct Model { point:[f64;3],initial:[f64;1] }
impl SensorFamily for Family {
    type Model=Model;
    fn bounds(&self)->&[[f64;2]] {&[[0.1,2.0],[0.2,3.0],[-1.0,1.0]]}
    fn instantiate(&self,p:&[f64])->Result<Model,String> {Ok(Model {point:[p[0],p[1],p[2]],initial:[p[1]]})}
}
impl OdeVjp for Model {
    fn dimension(&self)->usize {1}fn parameter_count(&self)->usize {3}
    fn rhs(&self,_:f64,x:&[f64],out:&mut[f64]) {out[0]=-self.point[0]*x[0];}
    fn rhs_vjp(&self,_:f64,x:&[f64],b:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        xb[0]=-self.point[0]*b[0];pb[0]=-x[0]*b[0];pb[1]=0.0;pb[2]=0.0;Ok(())
    }
}
impl SensorModel for Model {
    fn initial_values(&self)->&[f64] {&self.initial}
    fn initial_vjp(&self,b:&[f64],pb:&mut[f64])->Result<(),String> {pb.fill(0.0);pb[1]=b[0];Ok(())}
    fn predict(&self,c:u64,_:f64,x:&[f64])->Result<f64,String> {Ok((c+1) as f64*x[0]+self.point[2])}
    fn prediction_vjp(&self,c:u64,_:f64,_:&[f64],b:f64,xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        xb[0]=(c+1) as f64*b;pb.fill(0.0);pb[2]=b;Ok(())
    }
}
fn readings()->SensorData {
    let rows=[0.0,0.2,0.6,0.6,1.0,2.0,3.0,4.0].into_iter().enumerate().map(|(i,t)| {
        let c=(i%2) as u64;let value=(c+1) as f64*1.2*fs_math::det::exp(-0.7*t)+0.15;
        SensorReading::new(t,c,value,0.1,SensorLoss::Quadratic).unwrap()
    }).collect::<Vec<_>>();SensorData::new(&rows,8).unwrap()
}
fn config()->TransientConfig {TransientConfig {start:0.0,initial_step:0.1,
    recording:RecordingConfig {end:4.0,rtol:1e-10,atol:1e-12,controller:PiController::default(),max_workspace_components:1024},
    max_state_components:1,max_samples:8,max_attempts:10_000,max_records:10_000,
    replay:ReplayBudget {checkpoints:64,replayed_steps:100_000},max_kkt_dimension:9}}

#[test]
fn correlated_transient_gradient_matches_independent_dense_analytic_model() {
    let d=readings();let u=[0.2,0.2,0.0,0.0,-0.1,-0.1,0.3,0.3];let p=[0.9,1.1,0.1];
    let family=CorrelatedFamily::new(ObservedFamily::new(Family,d.clone()),1,&u,limits(),&mut||false).unwrap();
    let got=evaluate_transient(&family,&config(),&p,&mut||false).unwrap().unwrap();
    let y:Vec<_>=d.readings().iter().map(|r|(r.channel()+1) as f64*p[1]*(-p[0]*r.time()).exp()+p[2]).collect();
    let (value,seeds)=dense(&d,1,&u,&y);let mut gradient=[0.0;3];
    for (r,seed) in d.readings().iter().zip(seeds) {
        let e=(-p[0]*r.time()).exp();let scale=(r.channel()+1) as f64;
        gradient[0]-=seed*scale*p[1]*r.time()*e;gradient[1]+=seed*scale*e;gradient[2]+=seed;
    }
    assert!((got.value-value).abs()<2e-7);
    for (a,b) in got.gradient.iter().zip(gradient) {assert!((a-b).abs()<2e-6,"{a} != {b}");}
    assert!(got.replayed_steps>=2*got.forward.accepted);
}

#[test]
fn zero_shared_factors_keep_the_original_fit_path_bit_for_bit() {
    let d=readings();let p=[0.9,1.1,0.1];
    let ordinary=ObservedFamily::new(Family,d.clone());
    let correlated=CorrelatedFamily::new(ObservedFamily::new(Family,d),1,&[0.0;8],limits(),&mut||false).unwrap();
    assert!(correlated.noise().is_diagonal());
    let a=evaluate_transient(&ordinary,&config(),&p,&mut||false).unwrap();
    let b=evaluate_transient(&correlated,&config(),&p,&mut||false).unwrap();assert_eq!(a,b);
}

#[test]
fn existing_sqp_fits_with_shared_cross_time_noise() {
    let family=CorrelatedFamily::new(ObservedFamily::new(Family,readings()),1,&[0.2;8],limits(),&mut||false).unwrap();
    let mut study=TransientStudy::new(&family,&[1.1,1.5,-0.1],config(),&mut||false).unwrap();
    let before=study.accepted().value;let report=study.run(1e-6,100,1500,&mut||false).unwrap();
    assert_eq!(report.stop,crate::SqpStop::Converged,"{report:?}");
    for (value,target) in study.optimizer().point().iter().zip([0.7,1.2,0.15]) {assert!((value-target).abs()<2e-4);}
    assert!(study.accepted().value<before*1e-8);assert_eq!(study.accepted().point,study.optimizer().point());
}

#[test]
fn cancelled_profile_is_retryable_without_changing_prepared_noise() {
    let noise=SharedNoise::new(data(&[0.0;4],&[0.1;4]),1,&[0.2;4],limits(),&mut||false).unwrap();
    let mut calls=0;let expected=noise.score(&[0.1,0.2,0.3,0.4],&mut||{calls+=1;false}).unwrap();
    for stop in [1,3,calls/2,calls] {
        let mut count=0;assert_eq!(noise.score(&[0.1,0.2,0.3,0.4],&mut||{count+=1;count==stop}),Err(SharedNoiseError::Cancelled));
    }
    assert_eq!(noise.score(&[0.1,0.2,0.3,0.4],&mut||false).unwrap(),expected);
}
