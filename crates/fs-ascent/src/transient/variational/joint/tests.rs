use super::*;
use fs_time::PiController;
use fs_time::adaptive::adjoint::trajectory::{RecordingConfig, ReplayBudget};
use std::cell::Cell;

struct Family { calls: Cell<usize>, fail: Cell<bool>, missing: Cell<bool>, nonlinear: bool }
impl Family {
    fn new(nonlinear: bool) -> Self {
        Self { calls: Cell::new(0), fail: Cell::new(false), missing: Cell::new(false), nonlinear }
    }
}
struct Model { parameters: [f64; 2], missing: bool, nonlinear: bool }
impl ParameterFamily for Family {
    type Model = Model;
    fn instantiate(&self, p: &[f64], check: &mut dyn FnMut() -> bool) -> Result<Model, String> {
        self.calls.set(self.calls.get()+1);
        if check() { return Err("cancelled factory".into()); }
        if self.fail.get() { return Err("model data unavailable".into()); }
        Ok(Model { parameters: [p[0],p[1]], missing: self.missing.get(), nonlinear: self.nonlinear })
    }
}
const TIMES: [f64; 4] = [0.0, 0.2, 0.7, 1.4];
const TARGET: [f64; 4] = [1.1, 1.24, 1.51, 1.84];
impl OdeVjp for Model {
    fn dimension(&self) -> usize { 1 }
    fn parameter_count(&self) -> usize { 2 }
    fn rhs(&self, _: f64, x: &[f64], out: &mut [f64]) {
        out[0] = if self.nonlinear { -fs_math::det::exp(self.parameters[0])*x[0] } else { self.parameters[0] };
    }
    fn rhs_vjp(&self, _: f64, x: &[f64], b: &[f64], xb: &mut [f64], pb: &mut [f64]) -> Result<(), String> {
        if self.nonlinear {
            let rate = fs_math::det::exp(self.parameters[0]);
            xb[0] = -rate*b[0]; pb[0] = -rate*x[0]*b[0];
        } else { xb[0] = 0.0; pb[0] = b[0]; }
        pb[1] = 0.0; Ok(())
    }
}
impl WindowObjective for Model {
    fn evaluate(&self, times: &[f64], n: usize, states: &[f64], bar: &mut [f64], _: &mut dyn FnMut()->bool)
        -> Result<f64, String>
    {
        assert_eq!(times, TIMES); assert_eq!(n,1);
        let mut cost=0.0;
        for i in 0..4 {
            let r=(states[i]+self.parameters[1]-TARGET[i])/0.2;
            cost+=0.5*r*r;bar[i]=r/0.2;
        }
        Ok(cost)
    }
    fn parameter_partials(&self, _: &[f64], _: usize, states: &[f64], bar: &mut [f64], _: &mut dyn FnMut()->bool)
        -> Result<(), String>
    {
        if self.missing { return Ok(()); }
        bar[0]=0.0;bar[1]=(0..4).map(|i|(states[i]+self.parameters[1]-TARGET[i])/0.04).sum();Ok(())
    }
}
fn window() -> WeakConstraintWindow {
    WeakConstraintWindow::new(&TIMES, &[1.0,1.1,1.2,1.3], &[0.4], &[0.3], &[0.2,0.3,0.25], 20).unwrap()
}
fn joint<'a>(w: &'a WeakConstraintWindow, f: &'a Family) -> JointWindow<'a, Family> {
    JointWindow::new(w,f,&[0.2,0.05],&[0.3,0.2],&[0.8,0.5],6).unwrap()
}
fn policy() -> IntervalPolicy {
    IntervalPolicy { recording: RecordingConfig { end:1.4,rtol:1e-11,atol:1e-13,
        controller:PiController::default(),max_workspace_components:4096 },
        initial_step:0.1,max_attempts:10000,max_records:10000,
        replay:ReplayBudget{checkpoints:32,replayed_steps:100000} }
}
fn control() -> WindowControl { WindowControl::new(2000,6000,10000) }
fn settings() -> StudySettings {
    StudySettings { memory:8,gradient_tolerance:1e-7,max_evaluations:1000,max_optimizer_components:10000 }
}

#[test]
fn nonlinear_joint_gradient_includes_shared_dynamics_direct_sensor_and_prior_terms() {
    let w=window();let f=Family::new(true);let j=joint(&w,&f);let z=[0.1,-0.2,0.3,0.07,-0.15,0.25];
    let mut c=control();let got=j.evaluate(&z,&policy(),&mut c,&mut||false).unwrap();
    assert_eq!(f.calls.get(),1);assert_eq!(c.evaluations(),1);assert_eq!(c.interval_attempts(),3);
    assert_eq!(got.parameters,vec![0.3f64.mul_add(z[4],0.2),0.2f64.mul_add(z[5],0.05)]);
    assert_eq!(got.window.parameter_gradient.len(),2);
    let rate=got.parameters[0].exp();let mut parameter0=0.0;
    for k in 0..3 {
        let dt=TIMES[k+1]-TIMES[k];let predicted=got.window.states[k]*(-rate*dt).exp();
        let defect=got.window.states[k+1]-predicted;let sigma=[0.2,0.3,0.25][k];
        parameter0+=defect/(sigma*sigma)*rate*dt*predicted;
    }
    assert!((got.window.parameter_gradient[0]-parameter0).abs()<2e-8);
    let sensor=(0..4).map(|i|(got.window.states[i]+got.parameters[1]-TARGET[i])/0.04).sum::<f64>();
    assert!((got.window.parameter_gradient[1]-sensor).abs()<1e-12);
    let h=1e-5;
    for i in 0..6 {
        let (mut a,mut b)=(z,z);a[i]+=h;b[i]-=h;
        let fa=j.evaluate(&a,&policy(),&mut c,&mut||false).unwrap().value;
        let fb=j.evaluate(&b,&policy(),&mut c,&mut||false).unwrap().value;
        assert!((got.gradient[i]-(fa-fb)/(2.0*h)).abs()<2e-7,"coordinate {i}");
    }
}

// Independent dense normal equations for a linear state/forcing/bias problem.
fn exact_solution() -> Vec<f64> {
    let mut a=vec![vec![0.0;6];6];let mut b=vec![0.0;6];
    let mut term=|row:[f64;6],target:f64,sigma:f64| {
        for i in 0..6 {b[i]+=row[i]*target/(sigma*sigma);for k in 0..6 {a[i][k]+=row[i]*row[k]/(sigma*sigma);}}
    };
    for i in 0..4 {let mut row=[0.0;6];row[i]=1.0;row[5]=1.0;term(row,TARGET[i],0.2);}
    term([1.0,0.0,0.0,0.0,0.0,0.0],1.0,0.3);
    term([0.0,0.0,0.0,0.0,1.0,0.0],0.2,0.8);
    term([0.0,0.0,0.0,0.0,0.0,1.0],0.05,0.5);
    for i in 0..3 {let mut row=[0.0;6];row[i]=-1.0;row[i+1]=1.0;row[4]=-(TIMES[i+1]-TIMES[i]);term(row,0.0,[0.2,0.3,0.25][i]);}
    for k in 0..6 {
        for i in k+1..6 {let f=a[i][k]/a[k][k];for col in k..6 {let q=a[k][col];a[i][col]-=f*q;}let q=b[k];b[i]-=f*q;}
    }
    let mut x=vec![0.0;6];for i in (0..6).rev() {x[i]=(b[i]-(i+1..6).map(|k|a[i][k]*x[k]).sum::<f64>())/a[i][i];}x
}

#[test]
fn existing_optimizer_joint_map_matches_dense_gaussian_normal_equations() {
    let w=window();let f=Family::new(false);let j=joint(&w,&f);let mut c=control();
    let mut study=JointWindowStudy::new(&j,&[0.0;6],policy(),settings(),&mut c,&mut||false).unwrap();
    assert_eq!(study.run(250,&mut c,&mut||false).unwrap().reason,StopReason::GradNorm);
    let got=study.accepted();let exact=exact_solution();
    for (a,b) in got.window.states.iter().chain(&got.parameters).zip(exact) {assert!((a-b).abs()<2e-6,"{a} vs {b}");}
    assert_eq!(got.controls,study.optimizer().x);assert_eq!(got.gradient,study.optimizer().g);assert_eq!(got.value,study.optimizer().f);
    assert_eq!(c.evaluations(),study.optimizer().evals);assert_eq!(f.calls.get(),c.evaluations());
}

#[test]
fn joint_clone_resume_and_failed_factory_keep_accepted_pairing_and_spent_work() {
    let w=window();let f=Family::new(false);let j=joint(&w,&f);let mut c=control();
    let mut a=JointWindowStudy::new(&j,&[0.0;6],policy(),settings(),&mut c,&mut||false).unwrap();let mut b=a.clone();
    let before=a.accepted().clone();f.fail.set(true);
    assert!(matches!(a.run(1,&mut c,&mut||false),Err(LbfgsError::Evaluation(WindowError::Model(_)))));
    assert_eq!(a.accepted(),&before);assert_eq!(c.evaluations(),2);assert_eq!(c.interval_attempts(),3);
    f.fail.set(false);a.run(250,&mut c,&mut||false).unwrap();
    for _ in 0..250 {if b.run(1,&mut c,&mut||false).unwrap().reason!=StopReason::IterationCap {break;}}
    assert_eq!(a.accepted(),b.accepted());assert_eq!(a.optimizer().history,b.optimizer().history);
    assert_eq!(a.optimizer().evals,b.optimizer().evals+1);
}

#[test]
fn admission_precedes_factory_and_unwritten_parameter_partials_refuse() {
    let w=window();let f=Family::new(true);let j=joint(&w,&f);
    for (evals,intervals,work) in [(0,3,1000),(1,2,1000),(1,3,j.workspace_components().unwrap()-1)] {
        assert!(j.evaluate(&[0.0;6],&policy(),&mut WindowControl::new(evals,intervals,work),&mut||false).is_err());
        assert_eq!(f.calls.get(),0);
    }
    f.missing.set(true);let mut c=control();
    assert!(matches!(j.evaluate(&[0.0;6],&policy(),&mut c,&mut||false),Err(WindowError::NonFinite("observation parameter partials"))));
    assert_eq!(c.evaluations(),1);assert_eq!(c.interval_attempts(),0);
    f.missing.set(false);
    let got=j.evaluate(&[0.0;6],&policy(),&mut c,&mut||false).unwrap();
    assert_eq!(got.controls.len(),6);assert_eq!(c.evaluations(),2);
}

#[test]
fn every_joint_evaluation_cancellation_boundary_is_atomic_and_retryable() {
    let w=window();let f=Family::new(true);let j=joint(&w,&f);let mut polls=0;
    let expected=j.evaluate(&[0.0;6],&policy(),&mut control(),&mut||{polls+=1;false}).unwrap();
    for stop in 1..=polls {
        let mut seen=0;let mut c=control();
        assert_eq!(j.evaluate(&[0.0;6],&policy(),&mut c,&mut||{seen+=1;seen==stop}),Err(WindowError::Cancelled),"boundary {stop}");
        let spent=c.evaluations();assert_eq!(j.evaluate(&[0.0;6],&policy(),&mut c,&mut||false).unwrap(),expected);
        assert_eq!(c.evaluations(),spent+1);
    }
}

#[test]
fn invalid_parameter_priors_refuse_before_model_work() {
    let w=window();let f=Family::new(false);
    for bad in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(JointWindow::new(&w,&f,&[0.2,0.05],&[0.3,0.2],&[bad,0.5],6).is_err());
    }
    assert!(JointWindow::new(&w,&f,&[0.0],&[],&[1.0],6).is_err());
    assert!(JointWindow::new(&w,&f,&[0.0],&[1.0],&[1.0],2).is_err());assert_eq!(f.calls.get(),0);
}
