use super::*;
use fs_time::PiController;
use std::cell::Cell;

struct Decay { rate: f64, calls: Cell<usize>, fail: Cell<bool> }
impl Decay { fn new() -> Self { Self { rate: 0.7, calls: Cell::new(0), fail: Cell::new(false) } } }
impl OdeVjp for Decay {
    fn dimension(&self) -> usize { 1 }
    fn parameter_count(&self) -> usize { 0 }
    fn rhs(&self, _: f64, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get()+1); out[0] = if self.fail.get() { f64::NAN } else { -self.rate*x[0] };
    }
    fn rhs_vjp(&self, _: f64, _: &[f64], seed: &[f64], x: &mut [f64], p: &mut [f64]) -> Result<(), String> {
        x[0] = -self.rate*seed[0]; p.fill(0.0); Ok(())
    }
}
struct Observed { target: Vec<f64>, sigma: f64 }
impl WindowObjective for Observed {
    fn evaluate(&self, _: &[f64], _: usize, states: &[f64], bar: &mut [f64], _: &mut dyn FnMut()->bool) -> Result<f64, String> {
        if states.len() != self.target.len() { return Err("target length".into()); }
        let mut value = 0.0;
        for ((x, target), g) in states.iter().zip(&self.target).zip(bar) {
            let r = (x-target)/self.sigma; value += 0.5*r*r; *g = r/self.sigma;
        }
        Ok(value)
    }
}
fn policy() -> IntervalPolicy { IntervalPolicy {
    recording: RecordingConfig { end: 3.0, rtol: 1e-11, atol: 1e-13,
        controller: PiController::default(), max_workspace_components: 4096 },
    initial_step: 0.17, max_attempts: 10000, max_records: 10000,
    replay: ReplayBudget { checkpoints: 32, replayed_steps: 100000 },
} }
fn window() -> WeakConstraintWindow {
    WeakConstraintWindow::new(&[0.0, 0.13, 0.8, 2.0], &[1.2,1.0,0.7,0.3], &[2.0], &[0.4], &[0.1,0.2,0.3], 100).unwrap()
}
fn control() -> WindowControl { WindowControl::new(1000,3000,10000) }
fn analytic(w: &WeakConstraintWindow, model: &Decay, obs: &Observed, point: &[f64]) -> (f64, Vec<f64>) {
    let x = point.iter().zip(&w.reference).map(|(z,r)| r+2.0*z).collect::<Vec<_>>();
    let mut g = vec![0.0; x.len()]; let mut value = obs.evaluate(&w.times,1,&x,&mut g,&mut||false).unwrap();
    let r = (x[0]-w.reference[0])/w.background_sigma[0]; value += 0.5*r*r; g[0] += r/w.background_sigma[0];
    for k in 0..x.len()-1 {
        let a = (-model.rate*(w.times[k+1]-w.times[k])).exp();
        let r = (x[k+1]-a*x[k])/w.model_sigma[k]; let v = r/w.model_sigma[k];
        value += 0.5*r*r; g[k] -= a*v; g[k+1] += v;
    }
    for v in &mut g { *v *= 2.0; } (value,g)
}

#[test]
fn objective_and_both_endpoint_gradients_match_analytic_weak_constraint_cost() {
    let w=window(); let m=Decay::new(); let o=Observed{target:vec![1.3,1.1,0.9,0.4],sigma:0.2};
    let z=[0.07,-0.03,0.02,0.08]; let mut c=control();
    let result=w.evaluate(&m,&o,&z,&policy(),&mut c,&mut||false).unwrap();
    let (value,g)=analytic(&w,&m,&o,&z);
    assert!((result.value-value).abs()<1e-8);
    for (a,b) in result.gradient.iter().zip(&g) { assert!((a-b).abs()<1e-8,"{a} vs {b}"); }
    assert_eq!(c.evaluations(),1); assert_eq!(c.interval_attempts(),3);
    assert_eq!(result.controls,z); assert_eq!(result.defects.len(),3);
    assert!(result.model_value>0.0 && result.background_value>0.0 && result.observation_value>0.0);
    assert!(result.accepted_steps>0 && result.replayed_steps>=result.accepted_steps);
    let h=1e-5;
    for i in 0..z.len() {
        let (mut plus,mut minus)=(z,z); plus[i]+=h; minus[i]-=h;
        let a=w.evaluate(&m,&o,&plus,&policy(),&mut c,&mut||false).unwrap().value;
        let b=w.evaluate(&m,&o,&minus,&policy(),&mut c,&mut||false).unwrap().value;
        assert!((result.gradient[i]-(a-b)/(2.0*h)).abs()<1e-7);
    }
}

#[test]
fn interior_knot_has_incoming_and_outgoing_defect_contributions() {
    let w=WeakConstraintWindow::new(&[0.0,0.5,1.0],&[0.0;3],&[1.0],&[1.0],&[1.0;2],20).unwrap();
    let m=Decay{rate:0.0,calls:Cell::new(0),fail:Cell::new(false)};
    let o=Observed{target:vec![0.0;3],sigma:1.0};
    let got=w.evaluate(&m,&o,&[1.0,3.0,2.0],&policy(),&mut control(),&mut||false).unwrap();
    // Data [1,3,2], prior [1,0,0], defects [-2,3,-1].
    assert_eq!(got.gradient,vec![0.0,6.0,1.0]);
    assert_eq!(got.defects,vec![2.0,-1.0]);
    assert_eq!(got.value,10.0);
}

#[test]
fn no_dense_covariance_or_per_decision_tangent_for_large_state() {
    struct Constant(usize);
    impl OdeVjp for Constant {
        fn dimension(&self)->usize {self.0} fn parameter_count(&self)->usize {0}
        fn rhs(&self,_:f64,_:&[f64],out:&mut[f64]) {out.fill(0.0);}
        fn rhs_vjp(&self,_:f64,_:&[f64],_:&[f64],x:&mut[f64],p:&mut[f64])->Result<(),String> {x.fill(0.0);p.fill(0.0);Ok(())}
    }
    let n=513; let len=3*n;
    let w=WeakConstraintWindow::new(&[0.0,0.1,0.2],&vec![0.0;len],&vec![1.0;n],&vec![1.0;n],&vec![0.2;2*n],4*len).unwrap();
    let o=Observed{target:vec![0.1;len],sigma:0.3};
    let mut p=policy();p.recording.max_workspace_components=17*n;
    let required=w.workspace_components(0).unwrap(); let mut c=WindowControl::new(1,2,required);
    let got=w.evaluate(&Constant(n),&o,&vec![0.0;len],&p,&mut c,&mut||false).unwrap();
    assert_eq!(got.gradient.len(),1539);assert_eq!(got.defects,vec![0.0;2*n]);
    assert!(required<8*len);assert!(got.gradient.iter().all(|g|(*g+0.1/0.09).abs()<1e-14));
}

#[test]
fn shape_resource_and_interval_failures_never_return_partial_gradients() {
    let w=window();let m=Decay::new();let o=Observed{target:vec![0.0;4],sigma:0.2};let z=[0.0;4];
    for (evals,intervals,work) in [(0,3,1000),(1,2,1000),(1,3,w.workspace_components(0).unwrap()-1)] {
        let mut c=WindowControl::new(evals,intervals,work);
        assert!(w.evaluate(&m,&o,&z,&policy(),&mut c,&mut||false).is_err());
        assert_eq!(m.calls.get(),0); assert_eq!(c.evaluations(),0); assert_eq!(c.interval_attempts(),0);
    }
    let mut p=policy();p.max_attempts=0;let mut c=control();
    assert!(matches!(w.evaluate(&m,&o,&z,&p,&mut c,&mut||false),Err(WindowError::ForwardStopped{interval:0,status:RecordingStatus::AttemptLimit})));
    assert_eq!(c.evaluations(),1);assert_eq!(c.interval_attempts(),1);
    m.fail.set(true);
    assert!(matches!(w.evaluate(&m,&o,&z,&policy(),&mut c,&mut||false),Err(WindowError::Trajectory{interval:0,..})));
    assert_eq!(c.evaluations(),2);assert_eq!(c.interval_attempts(),2);
    m.fail.set(false);
    assert!(w.evaluate(&m,&o,&z,&policy(),&mut c,&mut||false).is_ok());
}

#[test]
fn observation_failures_unwritten_partials_and_latched_cancellation_refuse() {
    struct Broken(u8);
    impl WindowObjective for Broken {
        fn evaluate(&self,_:&[f64],_:usize,_:&[f64],bar:&mut[f64],cancel:&mut dyn FnMut()->bool)->Result<f64,String> {
            match self.0 {
                0 => Err("missing data".into()),
                1 => Ok(0.0),
                _ => {bar.fill(0.0);let _=cancel();Ok(0.0)}
            }
        }
    }
    let w=window();let m=Decay::new();let mut c=control();
    assert_eq!(w.evaluate(&m,&Broken(0),&[0.0;4],&policy(),&mut c,&mut||false),Err(WindowError::Observation("missing data".into())));
    assert!(matches!(w.evaluate(&m,&Broken(1),&[0.0;4],&policy(),&mut c,&mut||false),Err(WindowError::NonFinite(_))));
    let polls=Cell::new(0);
    // Entry poll + one state tile + callback poll. A one-shot true must latch.
    assert_eq!(w.evaluate(&m,&Broken(2),&[0.0;4],&policy(),&mut c,&mut||{polls.set(polls.get()+1);polls.get()==3}),Err(WindowError::Cancelled));
    assert_eq!(m.calls.get(),0);assert_eq!(c.evaluations(),3);assert_eq!(c.interval_attempts(),0);
}

#[test]
fn all_evaluation_cancellation_boundaries_refuse_and_retry_reproduces_output() {
    let w=window();let m=Decay::new();let o=Observed{target:vec![1.0;4],sigma:0.2};let z=[0.0;4];
    let mut polls=0;
    let expected=w.evaluate(&m,&o,&z,&policy(),&mut control(),&mut||{polls+=1;false}).unwrap();
    for stop in 1..=polls {
        let mut seen=0;let mut c=control();
        assert_eq!(w.evaluate(&m,&o,&z,&policy(),&mut c,&mut||{seen+=1;seen==stop}),Err(WindowError::Cancelled),"poll {stop}");
        let spent=c.evaluations();
        assert_eq!(w.evaluate(&m,&o,&z,&policy(),&mut c,&mut||false).unwrap(),expected);
        assert_eq!(c.evaluations(),spent+1);
    }
}

#[test]
fn invalid_covariance_scales_and_times_refuse_without_silent_repairs() {
    for times in [vec![0.0,0.0],vec![1.0,0.0],vec![0.0,f64::NAN]] {
        assert!(WeakConstraintWindow::new(&times,&[0.0;2],&[1.0],&[1.0],&[1.0],100).is_err());
    }
    for sigma in [0.0,-1.0,f64::INFINITY,f64::NAN] {
        assert!(WeakConstraintWindow::new(&[0.0,1.0],&[0.0;2],&[1.0],&[1.0],&[sigma],100).is_err());
    }
    let mut c=control();assert!(c.extend(999,3000,10000).is_err());
    assert!(c.extend(1000,3001,10001).is_ok());
}

#[test]
fn nonlinear_dynamics_and_joint_cross_time_observations_match_window_differences() {
    struct Nonlinear;
    impl OdeVjp for Nonlinear {
        fn dimension(&self)->usize {2} fn parameter_count(&self)->usize {0}
        fn rhs(&self,t:f64,x:&[f64],out:&mut[f64]) {
            out[0]=-0.3*x[0]+0.2*x[1].sin()+0.05*t;
            out[1]=-0.4*x[1]+0.1*x[0]*x[1];
        }
        fn rhs_vjp(&self,_:f64,x:&[f64],b:&[f64],xb:&mut[f64],p:&mut[f64])->Result<(),String> {
            xb[0]=-0.3*b[0]+0.1*x[1]*b[1];
            xb[1]=0.2*x[1].cos()*b[0]+(-0.4+0.1*x[0])*b[1];p.fill(0.0);Ok(())
        }
    }
    struct CrossTime;
    impl WindowObjective for CrossTime {
        fn evaluate(&self,_:&[f64],_:usize,x:&[f64],b:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {
            b.fill(0.0);let r=x[0]+2.0*x[3]-0.6;let s=x[1]*x[4]+0.1*x[5]-0.2;
            b[0]=r;b[3]=2.0*r;b[1]=x[4]*s;b[4]=x[1]*s;b[5]=0.1*s;Ok(0.5*(r*r+s*s))
        }
    }
    let w=WeakConstraintWindow::new(&[0.0,0.2,0.7],&[0.2,-0.3,0.15,-0.2,0.1,-0.1],
        &[0.8,1.3],&[0.4,0.6],&[0.1,0.2,0.3,0.1],100).unwrap();
    let point=[0.1,0.2,-0.1,0.05,0.15,-0.03];let mut c=control();
    let result=w.evaluate(&Nonlinear,&CrossTime,&point,&policy(),&mut c,&mut||false).unwrap();
    let h=1e-5;
    for j in 0..point.len() {
        let (mut plus,mut minus)=(point,point);plus[j]+=h;minus[j]-=h;
        let a=w.evaluate(&Nonlinear,&CrossTime,&plus,&policy(),&mut c,&mut||false).unwrap().value;
        let b=w.evaluate(&Nonlinear,&CrossTime,&minus,&policy(),&mut c,&mut||false).unwrap().value;
        assert!((result.gradient[j]-(a-b)/(2.0*h)).abs()<2e-7);
    }
}
