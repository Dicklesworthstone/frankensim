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
use fs_time::adaptive::adjoint::trajectory::{RecordingConfig, ReplayBudget};

fn settings() -> StudySettings { StudySettings {
    memory:8,gradient_tolerance:1e-7,max_evaluations:500,max_optimizer_components:10000,
} }
// Independent small dense solve, used ONLY as the linear-Gaussian MAP oracle.
fn oracle(w:&WeakConstraintWindow,m:&Decay,o:&Observed)->Vec<f64> {
    let n=w.times.len();let mut a=vec![vec![0.0;n];n];let mut b=vec![0.0;n];
    for i in 0..n {a[i][i]=1.0/(o.sigma*o.sigma);b[i]=o.target[i]*a[i][i];}
    let p=1.0/(w.background_sigma[0]*w.background_sigma[0]);a[0][0]+=p;b[0]+=p*w.reference[0];
    for i in 0..n-1 {
        let f=(-m.rate*(w.times[i+1]-w.times[i])).exp();let q=1.0/(w.model_sigma[i]*w.model_sigma[i]);
        a[i][i]+=f*f*q;a[i+1][i+1]+=q;a[i][i+1]-=f*q;a[i+1][i]-=f*q;
    }
    for k in 0..n {
        for i in k+1..n {let f=a[i][k]/a[k][k];for j in k..n {let pivot=a[k][j];a[i][j]-=f*pivot;}let pivot=b[k];b[i]-=f*pivot;}
    }
    let mut x=vec![0.0;n];
    for i in (0..n).rev() {x[i]=(b[i]-(i+1..n).map(|j|a[i][j]*x[j]).sum::<f64>())/a[i][i];}x
}

#[test]
fn existing_lbfgs_solves_the_weak_constraint_map_and_retains_exact_candidate_pairing() {
    let w=window();let m=Decay::new();let o=Observed{target:vec![1.2,1.3,1.5,0.9],sigma:0.2};let mut c=control();
    let mut study=WeakConstraintStudy::new(&w,&m,&o,&[0.0;4],policy(),settings(),&mut c,&mut||false).unwrap();
    let initial=study.accepted().value;
    let report=study.run(100,&mut c,&mut||false).unwrap();assert_eq!(report.reason,StopReason::GradNorm,"{report:?}");
    let exact=oracle(&w,&m,&o);
    for (x,y) in study.accepted().states.iter().zip(exact) {assert!((x-y).abs()<2e-7,"{x} vs {y}");}
    assert!(study.accepted().value<initial);
    assert_eq!(study.accepted().controls,study.optimizer().x);
    assert_eq!(study.accepted().gradient,study.optimizer().g);
    assert_eq!(study.accepted().value.to_bits(),study.optimizer().f.to_bits());
    assert_eq!(c.evaluations(),study.optimizer().evals);
}

#[test]
fn study_clone_and_accepted_iteration_continuation_are_bitwise_equal() {
    let w=window();let m=Decay::new();let o=Observed{target:vec![1.2,1.3,1.5,0.9],sigma:0.2};let mut c=control();
    let mut straight=WeakConstraintStudy::new(&w,&m,&o,&[0.0;4],policy(),settings(),&mut c,&mut||false).unwrap();
    let mut split=straight.clone();straight.run(100,&mut c,&mut||false).unwrap();
    for _ in 0..100 {if split.run(1,&mut c,&mut||false).unwrap().reason!=StopReason::IterationCap {break;}}
    assert_eq!(straight.accepted(),split.accepted());assert_eq!(straight.optimizer().history,split.optimizer().history);
    assert_eq!(straight.optimizer().evals,split.optimizer().evals);
}

#[test]
fn failed_trial_keeps_accepted_states_and_charges_work_before_retry() {
    let w=window();let m=Decay::new();let o=Observed{target:vec![1.2,1.3,1.5,0.9],sigma:0.2};let mut c=control();
    let mut study=WeakConstraintStudy::new(&w,&m,&o,&[0.0;4],policy(),settings(),&mut c,&mut||false).unwrap();
    let before=study.accepted().clone();let prior_calls=c.evaluations();m.fail.set(true);
    assert!(matches!(study.run(1,&mut c,&mut||false),Err(LbfgsError::Evaluation(WindowError::Trajectory{..}))));
    assert_eq!(study.accepted(),&before);assert_eq!(c.evaluations(),prior_calls+1);assert_eq!(study.optimizer().evals,2);
    m.fail.set(false);let previous=c.evaluations();
    assert!(matches!(study.run(1,&mut c,&mut||true),Err(LbfgsError::Evaluation(WindowError::Cancelled))));
    assert_eq!(c.evaluations(),previous);assert_eq!(study.accepted(),&before);
    assert_eq!(study.run(100,&mut c,&mut||false).unwrap().reason,StopReason::GradNorm);
}

#[test]
fn optimizer_budget_and_memory_admission_precede_extra_simulation() {
    let w=window();let m=Decay::new();let o=Observed{target:vec![1.2,1.3,1.5,0.9],sigma:0.2};let mut c=control();
    let mut s=settings();s.max_optimizer_components=0;
    assert!(WeakConstraintStudy::new(&w,&m,&o,&[0.0;4],policy(),s,&mut c,&mut||false).is_err());
    assert_eq!(m.calls.get(),0);assert_eq!(c.evaluations(),0);
    s=settings();s.max_evaluations=1;
    let mut study=WeakConstraintStudy::new(&w,&m,&o,&[0.0;4],policy(),s,&mut c,&mut||false).unwrap();
    let calls=m.calls.get();let before=study.accepted().clone();
    assert_eq!(study.run(100,&mut c,&mut||false).unwrap().reason,StopReason::Budget);
    assert_eq!(m.calls.get(),calls);assert_eq!(study.accepted(),&before);
}
