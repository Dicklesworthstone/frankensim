use super::*;
use std::cell::Cell;

struct Linear { n: usize, a: f64, calls: Cell<usize> }
impl Linear { fn new(n:usize,a:f64)->Self {Self{n,a,calls:Cell::new(0)}} }
impl WindowModel for Linear {
    type Tape<'a> = f64;
    fn dimension(&self)->usize {self.n}
    fn forecast<'a>(&'a self,_:usize,_:f64,_:f64,x:&[f64],out:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {
        self.calls.set(self.calls.get()+1);for (y,x) in out.iter_mut().zip(x) {*y=self.a*x;}Ok(self.a)
    }
    fn forecast_vjp(&self,a:f64,seed:&[f64],out:&mut[f64],_:&mut dyn FnMut()->bool)->Result<(),String> {
        self.calls.set(self.calls.get()+1);for (y,x) in out.iter_mut().zip(seed) {*y=a*x;}Ok(())
    }
    fn observe(&self,channel:u64,_:f64,x:&[f64],bar:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {
        self.calls.set(self.calls.get()+1);let j=channel as usize;if j>=x.len() {return Err("unknown channel".into());}
        bar.fill(0.0);bar[j]=1.0;Ok(x[j])
    }
}
fn limits()->WindowLimits {WindowLimits{max_frames:100,max_dimension:1024,max_observations:100,max_owned_components:100_000}}
fn control()->WindowControl {WindowControl::new(1000,10000,100000)}
fn make(n:usize,times:&[f64],reference:&[f64],scales:&[f64],obs:&[WindowObservation])->WeakConstraintWindow {
    WeakConstraintWindow::new(WindowInputs{times,background:&vec![0.8;n],background_std:&vec![0.5;n],
        model_std:&vec![0.3;n*(times.len()-1)],reference,scale:scales,observations:obs},limits(),&mut||false).unwrap()
}

#[test]
fn model_error_signs_and_control_scaling_match_the_state_form_objective() {
    let obs=[WindowObservation{id:1,frame:1,channel:0,value:1.1,sigma:0.1}];
    let w=make(1,&[0.0,0.4,1.0],&[0.7,0.9,0.6],&[2.0,0.5,3.0],&obs);let m=Linear::new(1,0.6);
    let e=w.evaluate(&m,&[0.0;3],&mut control(),&mut||false).unwrap();
    let d0=0.9-0.6*0.7;let d1=0.6-0.6*0.9;let r=0.9-1.1;
    let physical=[(0.7-0.8)/0.25-0.6*d0/0.09,d0/0.09-0.6*d1/0.09+r/0.01,d1/0.09];
    for ((a,g),s) in e.gradient().iter().zip(physical).zip([2.0,0.5,3.0]) {assert!((a-g*s).abs()<1e-12);}
    assert!((e.background_cost()-0.5*(0.7-0.8_f64).powi(2)/0.25).abs()<1e-14);
    assert!((e.model_cost()-0.5*(d0*d0+d1*d1)/0.09).abs()<1e-13);
    assert!((e.observation_cost()-0.5*r*r/0.01).abs()<1e-13);
    assert_eq!(e.states(),&[0.7,0.9,0.6]);assert_eq!(e.model_errors(),&[d0,d1]);assert_eq!(m.calls.get(),5);
}

struct Nonlinear;
impl WindowModel for Nonlinear {
    // Own the linearization made at the actual forecast state, not a new probe.
    type Tape<'a> = [f64;4];
    fn dimension(&self)->usize {2}
    fn forecast<'a>(&'a self,_:usize,a:f64,b:f64,x:&[f64],y:&mut[f64],_:&mut dyn FnMut()->bool)->Result<[f64;4],String> {
        let h=b-a;y[0]=x[0]+h*(x[1]+0.1*x[0]*x[0]);y[1]=x[1]+h*(x[0].sin()-0.3*x[1]);
        Ok([1.0+0.2*h*x[0],h,h*x[0].cos(),1.0-0.3*h])
    }
    fn forecast_vjp(&self,a:[f64;4],s:&[f64],out:&mut[f64],_:&mut dyn FnMut()->bool)->Result<(),String> {
        out[0]=a[0]*s[0]+a[2]*s[1];out[1]=a[1]*s[0]+a[3]*s[1];Ok(())
    }
    fn observe(&self,_:u64,_:f64,x:&[f64],bar:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {
        bar[0]=1.0;bar[1]=0.6*x[1];Ok(x[0]+0.3*x[1]*x[1])
    }
}
#[test]
fn nonlinear_taped_gradient_matches_independent_complete_cost_differences() {
    let obs=[WindowObservation{id:1,frame:0,channel:0,value:0.9,sigma:0.3},
        WindowObservation{id:2,frame:2,channel:0,value:1.1,sigma:0.2}];
    let w=make(2,&[0.0,0.2,0.65],&[0.5,0.1,0.6,0.2,0.8,0.25],&[1.0,0.5,2.0,0.2,0.3,1.4],&obs);
    let z=[0.03,-0.04,0.02,0.01,-0.02,0.04];let e=w.evaluate(&Nonlinear,&z,&mut control(),&mut||false).unwrap();
    let cost=|z:&[f64]| {
        let x:Vec<_>=w.reference.iter().zip(&w.scale).zip(z).map(|((r,s),z)|r+s*z).collect();
        let mut f=0.5*((x[0]-0.8).powi(2)+(x[1]-0.8).powi(2))/0.25;
        for k in 0..2 {let h=w.times[k+1]-w.times[k];let a=x[2*k];let b=x[2*k+1];
            let d=x[2*k+2]-(a+h*(b+0.1*a*a));let q=x[2*k+3]-(b+h*(a.sin()-0.3*b));f+=0.5*(d*d+q*q)/0.09;}
        for o in &obs {let a=x[o.frame*2];let b=x[o.frame*2+1];let r=(a+0.3*b*b-o.value)/o.sigma;f+=0.5*r*r;}f
    };
    assert!((e.value()-cost(&z)).abs()<1e-12);
    for j in 0..z.len() {let (mut a,mut b)=(z,z);a[j]+=1e-6;b[j]-=1e-6;
        assert!((e.gradient()[j]-(cost(&a)-cost(&b))/2e-6).abs()<2e-8,"component {j}");}
}

#[test]
fn exact_resource_limits_admit_once_and_never_refund_failed_work() {
    let obs=[WindowObservation{id:0,frame:1,channel:9,value:0.0,sigma:1.0}];
    let w=make(1,&[0.0,1.0],&[0.5,0.6],&[1.0;2],&obs);let m=Linear::new(1,0.8);
    let required=w.workspace_components().unwrap();assert_eq!(required,8);
    let mut c=WindowControl::new(1,3,required-1);
    assert!(matches!(w.evaluate(&m,&[0.0;2],&mut c,&mut||false),Err(WindowError::Limit{..})));assert_eq!(m.calls.get(),0);
    c.extend(1,3,required).unwrap();
    assert!(matches!(w.evaluate(&m,&[0.0;2],&mut c,&mut||false),Err(WindowError::Model{stage:WindowStage::Observation,..})));
    assert_eq!(c.evaluations(),1);assert_eq!(c.model_calls(),3);
    assert_eq!(w.evaluate(&m,&[0.0;2],&mut c,&mut||false),Err(WindowError::EvaluationLimit));
    assert!(c.extend(0,3,required).is_err());
    let mut too_few=WindowControl::new(1,2,required);
    assert_eq!(w.evaluate(&m,&[0.0;2],&mut too_few,&mut||false),Err(WindowError::ModelCallLimit));
    assert_eq!(too_few.evaluations(),0);
}

#[test]
fn every_cancel_boundary_returns_no_partial_evaluation_and_retry_is_identical() {
    let w=make(2,&[0.0,0.3,0.8],&[0.5;6],&[1.0;6],&[WindowObservation{id:0,frame:2,channel:0,value:1.0,sigma:0.2}]);
    let z=[0.0;6];let m=Nonlinear;let mut calls=0;
    let expected=w.evaluate(&m,&z,&mut control(),&mut||{calls+=1;false}).unwrap();
    for stop in 1..=calls {
        let mut c=control();let mut count=0;
        assert_eq!(w.evaluate(&m,&z,&mut c,&mut||{count+=1;count==stop}),Err(WindowError::Cancelled),"poll {stop}");
        let spent=c.model_calls();assert_eq!(w.evaluate(&m,&z,&mut c,&mut||false).unwrap(),expected);
        assert_eq!(c.model_calls(),spent+5);
    }
}

#[test]
fn missing_outputs_and_cancellation_inside_failed_callbacks_refuse() {
    struct Broken { stage:WindowStage, cancelled:bool }
    impl WindowModel for Broken {
        type Tape<'a> = ();
        fn dimension(&self)->usize {1}
        fn forecast<'a>(&'a self,_:usize,_:f64,_:f64,x:&[f64],out:&mut[f64],c:&mut dyn FnMut()->bool)->Result<(),String> {
            if self.cancelled {let _=c();return Err("model interrupted".into());}
            if self.stage!=WindowStage::Forecast {out.copy_from_slice(x);}Ok(())
        }
        fn forecast_vjp(&self,_:(),s:&[f64],out:&mut[f64],_:&mut dyn FnMut()->bool)->Result<(),String> {
            if self.stage!=WindowStage::Pullback {out.copy_from_slice(s);}Ok(())
        }
        fn observe(&self,_:u64,_:f64,x:&[f64],out:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {
            if self.stage!=WindowStage::Observation {out[0]=1.0;}Ok(x[0])
        }
    }
    let w=make(1,&[0.0,1.0],&[0.5;2],&[1.0;2],&[WindowObservation{id:0,frame:1,channel:0,value:0.7,sigma:1.0}]);
    for stage in [WindowStage::Forecast,WindowStage::Pullback,WindowStage::Observation] {
        assert!(matches!(w.evaluate(&Broken{stage,cancelled:false},&[0.0;2],&mut control(),&mut||false),Err(WindowError::NonFinite(_))));
    }
    // The callback sees the cancellation and returns an error; cancellation wins.
    let mut c=control();let mut polls=0;
    let result=w.evaluate(&Broken{stage:WindowStage::Forecast,cancelled:true},&[0.0;2],&mut c,&mut||{polls+=1;polls==6});
    assert_eq!(result,Err(WindowError::Cancelled));assert_eq!(c.model_calls(),1);
}

#[test]
fn malformed_time_noise_shapes_and_ids_are_rejected() {
    let build=|times:&[f64],noise:&[f64],obs:&[WindowObservation]|WeakConstraintWindow::new(WindowInputs{
        times,background:&[0.0],background_std:&[1.0],model_std:noise,reference:&[0.0;2],scale:&[1.0;2],observations:obs},limits(),&mut||false);
    for t in [[0.0,0.0],[1.0,0.0],[0.0,f64::NAN],[-f64::MAX,f64::MAX]] {assert!(build(&t,&[1.0],&[]).is_err());}
    for noise in [[0.0],[-1.0],[f64::INFINITY]] {assert!(build(&[0.0,1.0],&noise,&[]).is_err());}
    let o=WindowObservation{id:3,frame:1,channel:0,value:0.0,sigma:1.0};assert!(build(&[0.0,1.0],&[1.0],&[o,o]).is_err());
    assert!(build(&[0.0,1.0],&[1.0],&[WindowObservation{frame:2,..o}]).is_err());
    let w=build(&[0.0,1.0],&[1.0],&[]).unwrap();
    assert!(w.evaluate(&Linear::new(1,1.0),&[0.0,f64::NAN],&mut control(),&mut||false).is_err());
}

#[test]
fn dimension_above_dense_covariance_cap_and_static_window_need_no_covariance() {
    let n=513;let w=make(n,&[0.0,1.0],&vec![0.8;2*n],&vec![1.0;2*n],&[]);
    let mut c=WindowControl::new(1,2,w.workspace_components().unwrap());
    let e=w.evaluate(&Linear::new(n,1.0),&vec![0.0;2*n],&mut c,&mut||false).unwrap();
    assert_eq!(e.value(),0.0);assert!(e.gradient().iter().all(|x|*x==0.0));assert_eq!(w.workspace_components().unwrap(),8*n);
    let w=make(1,&[2.0],&[0.8],&[1.0],&[]);let m=Linear::new(1,2.0);
    assert_eq!(w.evaluate(&m,&[0.0],&mut WindowControl::new(1,0,5),&mut||false).unwrap().value(),0.0);
    assert_eq!(m.calls.get(),0);
}

#[test]
fn consistent_unit_changes_preserve_cost_and_dimensionless_gradient() {
    let evaluate=|unit:f64| {
        let obs=[WindowObservation{id:0,frame:1,channel:0,value:1.1*unit,sigma:0.1*unit}];
        let w=WeakConstraintWindow::new(WindowInputs{times:&[0.0,1.0],background:&[0.8*unit],background_std:&[0.5*unit],
            model_std:&[0.3*unit],reference:&[0.7*unit,0.9*unit],scale:&[2.0*unit,0.5*unit],observations:&obs},limits(),&mut||false).unwrap();
        w.evaluate(&Linear::new(1,0.6),&[0.0;2],&mut control(),&mut||false).unwrap()
    };
    let base=evaluate(1.0);
    for exponent in [-400,-100,100,400] {let e=evaluate(2.0_f64.powi(exponent));
        assert!((base.value()-e.value()).abs()<1e-12);for (a,b) in base.gradient().iter().zip(e.gradient()) {assert!((a-b).abs()<1e-12);}}
}
