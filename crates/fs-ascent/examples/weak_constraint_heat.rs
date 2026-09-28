//! Weak-constraint reconstruction of a heat trajectory with omitted heating.
//!
//! A manufactured 17-interior-node rod has fixed 300 K boundaries and
//! diffusivity 0.02 m^2/s. A synthetic endpoint heat increment is deliberately
//! absent from the forecast model. Infer all knot states using sparse readings,
//! an explicit initial prior, and explicit independent model-increment scales.
//! This example supplies only the spatial RHS and its transpose; RK45 and
//! L-BFGS are the existing production implementations.
//!
//! This is not validated conduction physics, an inferred continuous heat source,
//! a confidence interval, or a certified derivative/physical-error bound.

use fs_ascent::transient::variational::{IntervalPolicy, WeakConstraintWindow,
    WindowControl, WindowError, WindowObjective};
use fs_ascent::transient::variational::study::{StudySettings, WeakConstraintStudy};
use fs_ascent::StopReason;
use fs_time::{AdaptiveState, PiController};
use fs_time::adaptive::adjoint::OdeVjp;
use fs_time::adaptive::adjoint::trajectory::{RecordedRk45,RecordingConfig,RecordingStatus,ReplayBudget};

const N: usize = 17;
const TIMES: [f64;5] = [0.0,0.1,0.3,0.6,1.0];
const NODES: [usize;3] = [3,8,13];
struct Heat;
impl OdeVjp for Heat {
    fn dimension(&self)->usize {N}
    fn parameter_count(&self)->usize {0}
    fn rhs(&self,_:f64,x:&[f64],out:&mut[f64]) {
        let c=0.02*((N+1)*(N+1)) as f64;
        for i in 0..N {
            let left=if i==0 {300.0} else {x[i-1]};let right=if i+1==N {300.0} else {x[i+1]};
            out[i]=c*((left-x[i])+(right-x[i]));
        }
    }
    fn rhs_vjp(&self,_:f64,_:&[f64],seed:&[f64],x:&mut[f64],p:&mut[f64])->Result<(),String> {
        x.fill(0.0);p.fill(0.0);let c=0.02*((N+1)*(N+1)) as f64;
        for i in 0..N {
            x[i]-=2.0*c*seed[i];if i>0 {x[i-1]+=c*seed[i];}if i+1<N {x[i+1]+=c*seed[i];}
        }
        Ok(())
    }
}
struct Sensors { truth: Vec<f64> }
impl WindowObjective for Sensors {
    fn evaluate(&self,times:&[f64],n:usize,x:&[f64],bar:&mut[f64],cancel:&mut dyn FnMut()->bool)->Result<f64,String> {
        if n!=N || times!=TIMES.as_slice() || x.len()!=self.truth.len() {return Err("sensor layout mismatch".into());}
        bar.fill(0.0);let mut cost=0.0;
        for k in 0..TIMES.len() {
            if cancel() {return Err("sensor evaluation cancelled".into());}
            for i in NODES {let j=k*N+i;let r=(x[j]-self.truth[j])/0.1;cost+=0.5*r*r;bar[j]=r/0.1;}
        }
        Ok(cost)
    }
}
fn policy()->IntervalPolicy {IntervalPolicy {
    recording:RecordingConfig {end:1.0,rtol:1e-11,atol:1e-12,controller:PiController::default(),max_workspace_components:17*N},
    initial_step:0.01,max_attempts:10000,max_records:10000,replay:ReplayBudget{checkpoints:32,replayed_steps:100000},
}}
fn reference(heated:bool)->Result<Vec<f64>,Box<dyn std::error::Error>> {
    let model=Heat;let p=policy();let mut states=vec![0.0;N*TIMES.len()];
    for (i,x) in states[..N].iter_mut().enumerate() {*x=300.0+5.0*(((i+1) as f64)*std::f64::consts::PI/(N+1) as f64).sin();}
    for k in 0..TIMES.len()-1 {
        let mut config=p.recording.clone();config.end=TIMES[k+1];
        let initial=AdaptiveState::new(TIMES[k],&states[k*N..(k+1)*N],p.initial_step);
        let mut tape=RecordedRk45::new(&model,initial,config)?;
        if tape.advance(p.max_attempts,p.max_records,&mut||false)?.status!=RecordingStatus::ReachedEnd {
            return Err("synthetic forward solve did not reach its endpoint".into());
        }
        states[(k+1)*N..(k+2)*N].copy_from_slice(&tape.state().u);
        if heated && k==1 {
            for (i,x) in states[(k+1)*N..(k+2)*N].iter_mut().enumerate() {
                *x+=1.5*(((i+1) as f64)*std::f64::consts::PI/(N+1) as f64).sin();
            }
        }
    }
    Ok(states)
}
fn window(reference:&[f64],model_sigma:f64)->Result<WeakConstraintWindow,WindowError> {
    WeakConstraintWindow::new(&TIMES,reference,&[0.25;N],&[0.15;N],&vec![model_sigma;N*(TIMES.len()-1)],1000)
}
fn settings()->StudySettings {StudySettings {
    memory:12,gradient_tolerance:1e-5,max_evaluations:2000,max_optimizer_components:100000,
}}
fn main()->Result<(),Box<dyn std::error::Error>> {
    if std::env::args().len()!=1 {return Err("usage: weak_constraint_heat (manufactured demonstration; no input file)".into());}
    let prior=reference(false)?;let sensors=Sensors{truth:reference(true)?};let model=Heat;
    let w=window(&prior,0.6)?;let mut control=WindowControl::new(2000,8000,10000);
    let mut study=WeakConstraintStudy::new(&w,&model,&sensors,&vec![0.0;w.control_dimension()],policy(),settings(),&mut control,&mut||false)?;
    let before=study.accepted().clone();let report=study.run(600,&mut control,&mut||false)?;let result=study.accepted();
    println!("source=synthetic-noiseless model=rod-with-omitted-heat-increment nodes={N} knots={} controls={} stop={:?}",TIMES.len(),w.control_dimension(),report.reason);
    println!("objective_before={:.10e} objective_after={:.10e} observation_loss={:.10e} background_penalty={:.10e} model_penalty={:.10e}",before.value,result.value,result.observation_value,result.background_value,result.model_value);
    println!("gradient_inf={:.5e} evaluations={} interval_attempts={}",report.grad_norm,control.evaluations(),control.interval_attempts());
    println!("time_s,center_temperature_k,center_model_increment_k");
    for (k,t) in TIMES.iter().enumerate() {
        let defect=if k==0 {0.0} else {result.defects[(k-1)*N+N/2]};
        println!("{t},{:.8},{defect:.8}",result.states[k*N+N/2]);
    }
    println!("scope=weak-constraint numerical state estimate; endpoint defects are not identified heat-input physics");
    if report.reason!=StopReason::GradNorm {return Err("study stopped without meeting its declared gradient tolerance".into());}
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_heat_data_are_reconciled_using_nonzero_model_increments() {
        let prior=reference(false).unwrap();let obs=Sensors{truth:reference(true).unwrap()};let model=Heat;
        let w=window(&prior,0.6).unwrap();let mut c=WindowControl::new(2000,8000,10000);
        let mut study=WeakConstraintStudy::new(&w,&model,&obs,&vec![0.0;w.control_dimension()],policy(),settings(),&mut c,&mut||false).unwrap();
        let initial=study.accepted().observation_value;
        let report=study.run(600,&mut c,&mut||false).unwrap();
        assert_eq!(report.reason,StopReason::GradNorm,"{report:?}");
        assert!(study.accepted().observation_value<initial*0.02);
        assert!(study.accepted().model_value>0.0);
        assert!(study.accepted().defects.iter().any(|x|x.abs()>0.2));
        assert_eq!(study.accepted().controls,study.optimizer().x);
    }
    #[test]
    fn tighter_declared_model_error_reduces_the_fitted_defect_norm() {
        let prior=reference(false).unwrap();let obs=Sensors{truth:reference(true).unwrap()};let model=Heat;
        let mut defect_norms=Vec::new();let mut losses=Vec::new();
        for sigma in [0.6,0.05] {
            let w=window(&prior,sigma).unwrap();let mut c=WindowControl::new(2000,8000,10000);
            let mut study=WeakConstraintStudy::new(&w,&model,&obs,&vec![0.0;w.control_dimension()],policy(),settings(),&mut c,&mut||false).unwrap();
            let report=study.run(600,&mut c,&mut||false).unwrap();assert_eq!(report.reason,StopReason::GradNorm,"{report:?}");
            defect_norms.push(study.accepted().defects.iter().map(|x|x*x).sum::<f64>());
            losses.push(study.accepted().observation_value);
        }
        assert!(defect_norms[1]<defect_norms[0]);assert!(losses[1]>losses[0]);
    }
}
