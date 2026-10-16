use super::*;
use crate::transient::{TransientConfig, evaluate_transient};
use crate::transient::observations::{ObservedFamily, SensorData, SensorLoss, SensorReading};
use fs_time::{AdaptiveState, PiController, rk45_adaptive_checked};
use fs_time::adaptive::adjoint::trajectory::{RecordedRk45, RecordingConfig, RecordingStatus, ReplayBudget};

const BOUNDS: [[f64;2];6] = [[0.1,2.0],[0.2,3.0],[0.05,2.0],[-1.0,1.0],[-2.0,2.0],[0.5,2.0]];
struct Family;
struct Model { p: [f64;6], initial: [f64;1] }
impl SensorFamily for Family {
    type Model = Model;
    fn bounds(&self) -> &[[f64;2]] { &BOUNDS }
    fn instantiate(&self, p: &[f64]) -> Result<Model,String> {
        Ok(Model { p: p.try_into().map_err(|_| "six parameters required")?, initial: [p[1]] })
    }
}
impl OdeVjp for Model {
    fn dimension(&self) -> usize {1}
    fn parameter_count(&self) -> usize {6}
    fn rhs(&self, _:f64, x:&[f64], out:&mut[f64]) {out[0]=-self.p[0]*x[0];}
    fn rhs_vjp(&self, _:f64, x:&[f64], b:&[f64], xb:&mut[f64], pb:&mut[f64]) -> Result<(),String> {
        xb[0]=-self.p[0]*b[0];pb.fill(0.0);pb[0]=-x[0]*b[0];Ok(())
    }
}
impl SensorModel for Model {
    fn initial_values(&self) -> &[f64] {&self.initial}
    fn initial_vjp(&self, b:&[f64], pb:&mut[f64]) -> Result<(),String> {pb.fill(0.0);pb[1]=b[0];Ok(())}
    fn predict(&self, channel:u64, t:f64, x:&[f64]) -> Result<f64,String> {
        match channel {0=>Ok(x[0]),10=>Ok(self.p[5]*x[0]+self.p[3]*t),_=>Err("unknown source".into())}
    }
    fn prediction_vjp(&self, channel:u64, t:f64, x:&[f64], b:f64, xb:&mut[f64], pb:&mut[f64]) -> Result<(),String> {
        pb.fill(0.0);
        match channel {0=>{xb[0]=b;},10=>{xb[0]=self.p[5]*b;pb[5]=x[0]*b;pb[3]=t*b;},_=>return Err("unknown source".into())}
        Ok(())
    }
}
fn sensors() -> [LagSensor;2] {[
    LagSensor {channel:10,source:10,time_constant:LagValue::Parameter(2),initial:LagInitial::Equilibrium,bias:LagValue::Parameter(3)},
    LagSensor {channel:100,source:10,time_constant:LagValue::Fixed(0.8),initial:LagInitial::Value(LagValue::Parameter(4)),bias:LagValue::Fixed(0.0)},
]}
fn point() -> [f64;6] {[0.7,1.2,0.4,0.15,-0.3,1.1]}
fn dot(a:&[f64],b:&[f64]) -> f64 {a.iter().zip(b).map(|(x,y)|x*y).sum()}
fn close(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<tol,"{a} != {b}");}
fn cfg(end:f64) -> TransientConfig {TransientConfig {start:0.0,initial_step:0.04,
    recording:RecordingConfig{end,rtol:1e-10,atol:1e-12,controller:PiController::default(),max_workspace_components:4096},
    max_state_components:3,max_samples:32,max_attempts:10000,max_records:10000,
    replay:ReplayBudget{checkpoints:64,replayed_steps:100000},max_kkt_dimension:18}}

#[test]
fn augmented_rhs_transpose_matches_state_and_parameter_differences() {
    let family=LaggedFamily::new(Family,&sensors(),0.3,2,3).unwrap();
    let p=point();let model=family.instantiate(&p).unwrap();let x=[0.8,0.4,-0.1];let b=[0.2,-0.9,0.7];let t=0.9;
    let (mut xb,mut pb)=(vec![0.0;3],vec![0.0;6]);model.rhs_vjp(t,&x,&b,&mut xb,&mut pb).unwrap();
    let h=1e-6;
    let value=|p:&[f64],x:&[f64]| {let m=family.instantiate(p).unwrap();let mut f=vec![0.0;3];m.rhs(t,x,&mut f);dot(&b,&f)};
    for i in 0..3 {let (mut a,mut c)=(x,x);a[i]+=h;c[i]-=h;close(xb[i],(value(&p,&a)-value(&p,&c))/(2.0*h),2e-9);}
    for i in 0..6 {let (mut a,mut c)=(p,p);a[i]+=h;c[i]-=h;close(pb[i],(value(&a,&x)-value(&c,&x))/(2.0*h),2e-9);}
}

#[test]
fn equilibrium_initialization_propagates_direct_and_physical_dependencies() {
    let family=LaggedFamily::new(Family,&sensors(),0.3,2,3).unwrap();let p=point();let b=[0.7,-1.3,0.9];
    let model=family.instantiate(&p).unwrap();let mut pb=vec![0.0;6];model.initial_vjp(&b,&mut pb).unwrap();
    let value=|p:&[f64]| dot(family.instantiate(p).unwrap().initial_values(),&b);let h=1e-6;
    for i in 0..6 {let (mut a,mut c)=(p,p);a[i]+=h;c[i]-=h;close(pb[i],(value(&a)-value(&c))/(2.0*h),1e-9);}
    close(pb[1],0.7-1.3*p[5],1e-14);close(pb[3],-1.3*0.3,1e-14);close(pb[4],0.9,1e-14);
}

#[test]
fn channel_shadowing_passthrough_and_readout_bias_have_distinct_derivatives() {
    let family=LaggedFamily::new(Family,&sensors(),0.0,2,3).unwrap();let p=point();let model=family.instantiate(&p).unwrap();
    let x=[0.8,0.4,-0.1];let (mut xb,mut pb)=(vec![0.0;3],vec![0.0;6]);
    close(model.predict(10,0.5,&x).unwrap(),0.55,1e-15);
    model.prediction_vjp(10,0.5,&x,2.0,&mut xb,&mut pb).unwrap();assert_eq!(xb,vec![0.0,2.0,0.0]);assert_eq!(pb,vec![0.0,0.0,0.0,2.0,0.0,0.0]);
    model.prediction_vjp(0,0.5,&x,2.0,&mut xb,&mut pb).unwrap();assert_eq!(xb,vec![2.0,0.0,0.0]);assert_eq!(pb,vec![0.0;6]);
    let mut f=vec![0.0;3];model.rhs(0.5,&x,&mut f);
    close(f[1],(p[5]*x[0]+p[3]*0.5-x[1])/p[2],1e-15); // source 10 is INNER, not its shadow
}

#[test]
fn lagged_decay_matches_closed_form_including_coincident_time_constants() {
    for tau in [0.05,0.4,1.0/0.7] {
        let sensor=LagSensor{channel:20,source:0,time_constant:LagValue::Fixed(tau),initial:LagInitial::Equilibrium,bias:LagValue::Fixed(0.15)};
        let family=LaggedFamily::new(Family,&[sensor],0.0,1,2).unwrap();let m=family.instantiate(&point()).unwrap();
        let mut state=AdaptiveState::new(0.0,m.initial_values(),0.1);
        rk45_adaptive_checked(&mut state,&|t,x,f|m.rhs(t,x,f),2.0,1e-10,1e-12,&PiController::default(),10000).unwrap();
        let k=0.7_f64;let e=(-2.0*k).exp();
        let z=if (1.0-k*tau).abs()<1e-14 {1.2*e*(1.0+2.0*k)} else {1.2*(e-k*tau*(-2.0/tau).exp())/(1.0-k*tau)};
        close(m.predict(20,2.0,&state.u).unwrap(),z+0.15,2e-9);
    }
}

fn observed() -> ObservedFamily<LaggedFamily<Family>> {
    let rows=[(0.0,10,1.4),(0.17,100,0.2),(0.51,10,1.1),(0.51,0,0.8),(1.2,100,0.6),(2.0,10,0.5)]
        .map(|(t,c,y)|SensorReading::new(t,c,y,0.4,SensorLoss::Quadratic).unwrap());
    ObservedFamily::new(LaggedFamily::new(Family,&sensors(),0.0,2,3).unwrap(),SensorData::new(&rows,32).unwrap())
}

#[test]
fn full_sampled_adjoint_includes_time_constant_initial_gain_and_bias() {
    let family=observed();let p=point();let c=cfg(2.0);
    let got=evaluate_transient(&family,&c,&p,&mut||false).unwrap().unwrap();let h=1e-5;
    for i in 0..6 {let (mut a,mut b)=(p,p);a[i]+=h;b[i]-=h;
        let f=evaluate_transient(&family,&c,&a,&mut||false).unwrap().unwrap().value;
        let g=evaluate_transient(&family,&c,&b,&mut||false).unwrap().unwrap().value;
        close(got.gradient[i],(f-g)/(2.0*h),2e-6);
    }
}

#[test]
fn checkpoint_resume_retains_sensor_states_without_reinitializing_them() {
    use crate::transient::{TransientFamily,TransientModel};
    let family=observed();let m=family.instantiate(&point()).unwrap();let c=cfg(2.0);
    let create=||RecordedRk45::new_sampled(&m,AdaptiveState::new(0.0,m.initial_values(),0.04),c.recording.clone(),family.sample_times(),32).unwrap();
    let mut a=create();let mut b=create();a.advance(10000,10000,&mut||false).unwrap();
    b.advance(10000,2,&mut||false).unwrap();let mut b=b.clone();
    for _ in 0..10000 {if b.advance(1,10000,&mut||false).unwrap().status==RecordingStatus::ReachedEnd {break;}}
    assert_eq!(a.state().u,b.state().u);assert_eq!(a.state().h.to_bits(),b.state().h.to_bits());
    assert_eq!(a.pullback_samples(&m,c.replay,&mut||false).unwrap(),b.pullback_samples(&m,c.replay,&mut||false).unwrap());
    assert!(a.pullback_samples(&m,c.replay,&mut||true).is_err());
}

#[test]
fn malformed_sensor_layouts_and_unknown_sources_refuse() {
    let s=sensors();
    assert!(LaggedFamily::new(Family,&s,0.0,1,3).is_err());
    assert!(LaggedFamily::new(Family,&[s[0].clone(),s[0].clone()],0.0,2,3).is_err());
    for value in [LagValue::Fixed(0.0),LagValue::Fixed(-1.0),LagValue::Fixed(f64::NAN),LagValue::Parameter(100)] {
        let bad=LagSensor{time_constant:value,..s[0].clone()};assert!(LaggedFamily::new(Family,&[bad],0.0,1,2).is_err());
    }
    let family=LaggedFamily::new(Family,&s,0.0,2,2).err();assert!(family.is_some());
    let bad=LagSensor{source:999,..s[0].clone()};let family=LaggedFamily::new(Family,&[bad],0.0,1,2).unwrap();
    assert!(family.instantiate(&point()).is_err());
    let m=LaggedFamily::new(Family,&s,0.0,2,3).unwrap().instantiate(&point()).unwrap();
    assert!(m.predict(999,0.0,m.initial_values()).is_err());
    assert!(m.prediction_vjp(10,0.0,m.initial_values(),1.0,&mut[0.0;1],&mut[0.0;6]).is_err());
}
