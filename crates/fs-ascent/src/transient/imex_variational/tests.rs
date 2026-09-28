use super::*;
use crate::transient::variational::{WeakConstraintWindow,WindowControl,WindowObjective};
use crate::transient::variational::study::{WeakConstraintStudy,StudySettings};
use fs_solver::LinearOp;
use fs_time::stiff::IdentityPreconditioner;
use std::cell::Cell;

// Unequal thermal capacities give a nonsymmetric stiff operator. The fixture
// models two temperatures relative to a fixed ambient reference, in kelvin.
struct Thermal { source:f64, nonlinear:bool, calls:Cell<usize>, fail:Cell<bool> }
impl Thermal {fn new(source:f64)->Self {Self{source,nonlinear:false,calls:Cell::new(0),fail:Cell::new(false)}}}
impl LinearOp for Thermal {
    fn n(&self)->usize {2}
    fn apply(&self,x:&[f64],out:&mut[f64]) {out[0]=-401.0*x[0]+400.0*x[1];out[1]=4.0*x[0]-5.0*x[1];}
    fn apply_transpose(&self,x:&[f64],out:&mut[f64]) {out[0]=-401.0*x[0]+4.0*x[1];out[1]=400.0*x[0]-5.0*x[1];}
}
impl ImexVjp for Thermal {
    fn parameter_count(&self)->usize {1}
    fn nonlinear(&self,x:&[f64],out:&mut[f64]) {
        self.calls.set(self.calls.get()+1);out[0]=0.0;
        out[1]=self.source+if self.nonlinear {0.03*x[0]*x[0]} else {0.0};
    }
    fn nonlinear_vjp(&self,x:&[f64],b:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        if self.fail.get() {return Err("injected thermal derivative refusal".into());}
        xb[0]=if self.nonlinear {0.06*x[0]*b[1]} else {0.0};xb[1]=0.0;pb[0]=b[1];Ok(())
    }
    fn linear_parameter_vjp(&self,_:&[f64],_:&[f64],out:&mut[f64])->Result<(),String> {out.fill(0.0);Ok(())}
}
struct Data(Vec<f64>);
impl WindowObjective for Data {
    fn evaluate(&self,_:&[f64],_:usize,x:&[f64],bar:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {
        let mut value=0.0;
        for ((x,y),g) in x.iter().zip(&self.0).zip(bar) {let r=(x-y)/0.1;*g=r/0.1;value+=0.5*r*r;}
        Ok(value)
    }
}
fn config()->ImexWindowConfig {ImexWindowConfig{
    step:0.125,solve:ImexSolveConfig{tolerance:1e-12,restart:2,max_cycles:8},
    max_workspace_components:4096,max_forward_steps:16,max_records:16,
    replay:ImexReplayBudget{checkpoints:16,forward_steps:1000},max_clock_steps:1024,max_grid_components:64,
}}
fn make_window(times:&[f64],sigma:f64)->WeakConstraintWindow {
    let n=times.len()*2;
    WeakConstraintWindow::new(times,&vec![1.0;n],&[1.0,2.0],&[0.5,0.5],&vec![sigma;n-2],4*n).unwrap()
}
fn control()->WindowControl {WindowControl::new(5000,20000,10000)}
fn settings()->StudySettings {StudySettings{memory:7,gradient_tolerance:1e-6,max_evaluations:2000,max_optimizer_components:100000}}

#[test]
fn imex_interval_is_the_production_recording_and_transposed_pullback() {
    let model=Thermal::new(0.4);let identity=IdentityPreconditioner;
    let p=ImexWindowPolicy::new(0.0,&[3],config(),&identity,&identity,&mut||false).unwrap();
    let tape=p.record(&model,0,p.times()[0],p.times()[1],&[1.0,0.2],&mut||false).unwrap();
    let mut direct=RecordedImex2::new(OperatorImex2::new(2,config().step,config().solve),&model,&identity,
        0.0,&[1.0,0.2],ImexRecordingConfig{steps:3,max_workspace_components:4096}).unwrap();
    direct.advance(16,16,&mut||false).unwrap();
    assert_eq!(tape.endpoint(),direct.state());assert_eq!(tape.end_time(),direct.time());
    let got=tape.pullback(&[0.3,-0.8],&[0.1],&mut||false).unwrap();
    let expected=direct.pullback(&[0.3,-0.8],&[0.1],&identity,config().replay,&mut||false).unwrap();
    assert_eq!(got.initial,expected.initial);assert_eq!(got.parameters,expected.parameters);
    assert_eq!(got.replayed_steps,expected.replayed_steps);
}

#[test]
fn nonlinear_stiff_window_gradients_match_state_and_parameter_differences() {
    let mut model=Thermal::new(0.4);model.nonlinear=true;let identity=IdentityPreconditioner;
    let p=ImexWindowPolicy::new(0.0,&[2,3],config(),&identity,&identity,&mut||false).unwrap();
    let w=make_window(p.times(),0.3);let data=Data(vec![0.7;6]);let z=[0.1,-0.03,0.05,0.01,-0.07,0.08];
    let got=w.evaluate_using(&model,&data,&z,&p,&mut control(),&mut||false).unwrap();
    for i in 0..z.len() {
        let (mut plus,mut minus)=(z,z);let h=1e-5;plus[i]+=h;minus[i]-=h;
        let a=w.evaluate_using(&model,&data,&plus,&p,&mut control(),&mut||false).unwrap();
        let b=w.evaluate_using(&model,&data,&minus,&p,&mut control(),&mut||false).unwrap();
        assert!((got.gradient[i]-(a.value-b.value)/(2.0*h)).abs()<2e-5,"coordinate {i}");
    }
    let h=1e-5;model.source+=h;
    let a=w.evaluate_using(&model,&data,&z,&p,&mut control(),&mut||false).unwrap().value;
    model.source-=2.0*h;
    let b=w.evaluate_using(&model,&data,&z,&p,&mut control(),&mut||false).unwrap().value;
    assert!((got.parameter_gradient[0]-(a-b)/(2.0*h)).abs()<2e-5);
    assert_eq!(got.accepted_steps,5);
}

#[test]
fn fixed_step_clock_is_explicit_and_never_rounded_to_observation_times() {
    let id=IdentityPreconditioner;let mut c=config();c.step=0.1;
    let p=ImexWindowPolicy::new(0.0,&[3,2],c,&id,&id,&mut||false).unwrap();
    let mut t=0.0;for _ in 0..3 {t+=0.1;}assert_eq!(p.times()[1],t);
    let w=make_window(&[0.0,0.3,0.5],0.3);let m=Thermal::new(0.0);let mut work=control();
    assert!(w.evaluate_using(&m,&Data(vec![0.0;6]),&[0.0;6],&p,&mut work,&mut||false).is_err());
    assert_eq!(work.evaluations(),0);assert_eq!(m.calls.get(),0);
    assert!(ImexWindowPolicy::new(1e100,&[1],c,&id,&id,&mut||false).is_err());
    assert!(ImexWindowPolicy::new(0.0,&[0],c,&id,&id,&mut||false).is_err());
    c.max_clock_steps=2;assert!(ImexWindowPolicy::new(0.0,&[3],c,&id,&id,&mut||false).is_err());
}

#[test]
fn forward_record_and_reverse_limits_produce_no_window_result() {
    let id=IdentityPreconditioner;let model=Thermal::new(0.4);
    for kind in 0..3 {
        let mut c=config();
        if kind==0 {c.max_forward_steps=1;} else if kind==1 {c.max_records=1;} else {c.replay.forward_steps=0;}
        let p=ImexWindowPolicy::new(0.0,&[3],c,&id,&id,&mut||false).unwrap();let w=make_window(p.times(),0.3);let mut work=control();
        assert!(matches!(w.evaluate_using(&model,&Data(vec![0.0;4]),&[0.0;4],&p,&mut work,&mut||false),Err(WindowError::Integrator{interval:0,..})));
        assert_eq!(work.evaluations(),1);assert_eq!(work.interval_attempts(),1);
    }
}

#[test]
fn adjoint_uses_the_explicit_transposed_solve_preconditioner() {
    struct Count(Cell<usize>);
    impl FlexiblePreconditioner for Count {
        fn apply(&self,_:usize,r:&[f64],out:&mut[f64]) {self.0.set(self.0.get()+1);out.copy_from_slice(r);}
    }
    let primal=Count(Cell::new(0));let adjoint=Count(Cell::new(0));let m=Thermal::new(0.4);
    let p=ImexWindowPolicy::new(0.0,&[2],config(),&primal,&adjoint,&mut||false).unwrap();
    let w=make_window(p.times(),0.3);
    w.evaluate_using(&m,&Data(vec![0.0;4]),&[0.0;4],&p,&mut control(),&mut||false).unwrap();
    assert!(primal.0.get()>0);assert!(adjoint.0.get()>0);
}

#[test]
fn cancelled_or_failed_thermal_trial_keeps_the_accepted_study_pair() {
    let id=IdentityPreconditioner;let m=Thermal::new(0.4);
    let p=ImexWindowPolicy::new(0.0,&[2,2],config(),&id,&id,&mut||false).unwrap();let w=make_window(p.times(),0.3);let data=Data(vec![0.7;6]);let mut work=control();
    let mut study=WeakConstraintStudy::new(&w,&m,&data,&[0.0;6],p,settings(),&mut work,&mut||false).unwrap();
    let before=study.accepted().clone();m.calls.set(0);
    assert!(matches!(study.run(1,&mut work,&mut||m.calls.get()>0),Err(crate::LbfgsError::Evaluation(WindowError::Cancelled))));
    assert_eq!(study.accepted(),&before);
    m.fail.set(true);
    assert!(matches!(study.run(1,&mut work,&mut||false),Err(crate::LbfgsError::Evaluation(WindowError::Integrator{phase:"IMEX adjoint",..}))));
    assert_eq!(study.accepted(),&before);m.fail.set(false);
    assert!(study.run(100,&mut work,&mut||false).unwrap().f<before.value);
}

#[test]
fn stiff_missing_heating_is_reconstructed_and_model_error_allowance_matters() {
    let id=IdentityPreconditioner;let truth=Thermal::new(2.0);let model=Thermal::new(0.0);
    let p=ImexWindowPolicy::new(0.0,&[2,2,2],config(),&id,&id,&mut||false).unwrap();
    let mut readings=vec![1.0,1.0];let mut x=vec![1.0,1.0];
    for k in 0..3 {let t=p.record(&truth,k,p.times()[k],p.times()[k+1],&x,&mut||false).unwrap();x=t.endpoint().to_vec();readings.extend_from_slice(&x);}
    let data=Data(readings);
    let mut results=Vec::new();
    for sigma in [0.6,0.01] {
        let w=make_window(p.times(),sigma);let mut work=control();
        let mut study=WeakConstraintStudy::new(&w,&model,&data,&[0.0;8],p.clone(),settings(),&mut work,&mut||false).unwrap();
        let mut split=study.clone();let mut split_work=control();let initial=study.accepted().value;
        let report=study.run(300,&mut work,&mut||false).unwrap();
        assert_eq!(report.reason,crate::StopReason::GradNorm,"{report:?}");assert!(report.f<initial);
        for _ in 0..300 {if split.run(1,&mut split_work,&mut||false).unwrap().reason!=crate::StopReason::IterationCap {break;}}
        assert_eq!(study.accepted(),split.accepted());results.push(study.accepted().clone());
    }
    assert!(results[0].observation_value<results[1].observation_value);
    let norm=|v:&[f64]|v.iter().map(|x|x*x).sum::<f64>().sqrt();
    assert!(norm(&results[1].defects)<norm(&results[0].defects));
}
