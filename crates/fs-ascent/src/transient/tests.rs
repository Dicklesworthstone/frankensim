use super::*;
use fs_time::PiController;
use std::cell::Cell;
use std::rc::Rc;

// Sensor model y(t) = amplitude*exp(-rate*t) + bias. The three parameters
// enter different derivative paths: dynamics, initial state, and observations.
struct Family {
    bounds: [[f64;2];3], times: Vec<f64>, targets: Rc<Vec<f64>>,
    instances: Rc<Cell<usize>>, samples: Rc<Cell<usize>>, fail: Rc<Cell<bool>>,
}
impl Family {
    fn new() -> Self {
        let times: Vec<f64> = vec![0.0,0.13,0.4,0.8,1.3,2.0,3.0,4.0];
        let targets = Rc::new(times.iter().map(|t| 1.2*(-0.7*t).exp()+0.15).collect());
        Self { bounds:[[0.1,2.0],[0.2,3.0],[-0.5,0.5]], times, targets,
            instances:Rc::new(Cell::new(0)),samples:Rc::new(Cell::new(0)),fail:Rc::new(Cell::new(false)) }
    }
}
struct Model { p:[f64;3], initial:[f64;1], targets:Rc<Vec<f64>>, samples:Rc<Cell<usize>>, fail:Rc<Cell<bool>> }
impl OdeVjp for Model {
    fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {3}
    fn rhs(&self,_:f64,u:&[f64],out:&mut[f64]) {out[0]=-self.p[0]*u[0];}
    fn rhs_vjp(&self,_:f64,u:&[f64],b:&[f64],x:&mut[f64],p:&mut[f64])->Result<(),String> {
        if self.fail.get() {return Err("injected derivative refusal".into());}
        x[0]=-self.p[0]*b[0];p[0]=-u[0]*b[0];p[1]=0.0;p[2]=0.0;Ok(())
    }
}
impl SampleObjective for Model {
    fn evaluate(&self,i:usize,_:f64,u:&[f64],x:&mut[f64],p:&mut[f64])->Result<f64,String> {
        self.samples.set(self.samples.get()+1);
        let r=u[0]+self.p[2]-self.targets[i];x[0]=r;p[0]=0.0;p[1]=0.0;p[2]=r;Ok(0.5*r*r)
    }
}
impl TransientModel for Model {
    fn initial_values(&self)->&[f64] {&self.initial}
    fn initial_vjp(&self,b:&[f64],out:&mut[f64])->Result<(),String> {out[0]=0.0;out[1]=b[0];out[2]=0.0;Ok(())}
}
impl TransientFamily for Family {
    type Model=Model;
    fn bounds(&self)->&[[f64;2]] {&self.bounds}
    fn sample_times(&self)->&[f64] {&self.times}
    fn instantiate(&self,p:&[f64])->Result<Model,String> {
        self.instances.set(self.instances.get()+1);
        Ok(Model{p:[p[0],p[1],p[2]],initial:[p[1]],targets:self.targets.clone(),samples:self.samples.clone(),fail:self.fail.clone()})
    }
}
fn config()->TransientConfig {TransientConfig {start:0.0,initial_step:0.1,
    recording:RecordingConfig{end:4.0,rtol:1e-10,atol:1e-12,controller:PiController::default(),max_workspace_components:4096},
    max_state_components:1,max_samples:16,max_attempts:10000,max_records:10000,
    replay:ReplayBudget{checkpoints:64,replayed_steps:100000},max_kkt_dimension:9,
}}
fn analytic(f:&Family,p:&[f64])->(f64,Vec<f64>) {
    let mut cost=0.0;let mut g=vec![0.0;3];
    for (&t,&target) in f.times.iter().zip(f.targets.iter()) {
        let e=(-p[0]*t).exp();let r=p[1]*e+p[2]-target;
        cost+=0.5*r*r;g[0]-=t*p[1]*e*r;g[1]+=e*r;g[2]+=r;
    }
    (cost,g)
}

#[test]
fn evaluation_includes_dynamics_initial_and_direct_sensor_partials() {
    let f=Family::new();let p=[1.1,1.5,-0.1];
    let got=evaluate_transient(&f,&config(),&p,&mut||false).unwrap().unwrap();let (value,g)=analytic(&f,&p);
    assert!((got.value-value).abs()<2e-9);
    for (a,b) in got.gradient.iter().zip(g) {assert!((a-b).abs()<2e-8,"{a} != {b}");}
    assert_eq!(got.observations,f.times.len());assert_eq!(f.samples.get(),f.times.len());assert_eq!(f.instances.get(),1);
    assert_eq!(got.forward.status,RecordingStatus::ReachedEnd);
    let h=1e-5;
    for j in 0..3 {
        let (mut plus,mut minus)=(p,p);plus[j]+=h;minus[j]-=h;
        let a=evaluate_transient(&f,&config(),&plus,&mut||false).unwrap().unwrap().value;
        let b=evaluate_transient(&f,&config(),&minus,&mut||false).unwrap().unwrap().value;
        assert!((got.gradient[j]-(a-b)/(2.0*h)).abs()<2e-7);
    }
}

#[test]
fn native_sqp_fits_three_parameters_from_irregular_observations() {
    let f=Family::new();let mut study=TransientStudy::new(&f,&[1.4,2.2,-0.2],config(),&mut||false).unwrap();
    let initial=study.accepted().value;let report=study.run(1e-7,100,1500,&mut||false).unwrap();
    assert_eq!(report.stop,SqpStop::Converged,"{report:?}");
    for (actual,expected) in study.optimizer().point().iter().zip([0.7,1.2,0.15]) {assert!((actual-expected).abs()<2e-5);}
    assert!(study.accepted().value<initial*1e-8);
    assert_eq!(study.accepted().point,study.optimizer().point());
    assert_eq!(study.accepted().value,study.optimizer().sample().f);
    assert_eq!(study.accepted().gradient,study.optimizer().sample().gradient);
}

#[test]
fn accepted_step_continuation_reuses_the_existing_optimizer_state() {
    let f=Family::new();let point=[1.4,2.2,-0.2];
    let mut straight=TransientStudy::new(&f,&point,config(),&mut||false).unwrap();
    let mut split=TransientStudy::new(&f,&point,config(),&mut||false).unwrap();
    straight.run(1e-7,100,1500,&mut||false).unwrap();
    for _ in 0..100 {if split.run(1e-7,1,1500,&mut||false).unwrap().stop!=SqpStop::IterationLimit {break;}}
    assert_eq!(straight.optimizer().point(),split.optimizer().point());
    assert_eq!(straight.optimizer().history(),split.optimizer().history());
    assert_eq!(straight.optimizer().evaluations(),split.optimizer().evaluations());
    assert_eq!(straight.accepted(),split.accepted());
}

#[test]
fn evaluation_cap_and_failed_derivative_keep_accepted_data_and_spent_work() {
    let f=Family::new();let mut study=TransientStudy::new(&f,&[1.4,2.2,-0.2],config(),&mut||false).unwrap();
    let before=study.accepted().clone();let calls=f.instances.get();
    assert_eq!(study.run(1e-7,10,1,&mut||false).unwrap().stop,SqpStop::EvaluationLimit);
    assert_eq!(f.instances.get(),calls);assert_eq!(study.accepted(),&before);
    f.fail.set(true);
    assert!(matches!(study.run(1e-7,1,100,&mut||false),Err(SqpError::Evaluation(TransientError::Trajectory(TrajectoryError::Step(AdjointError::Derivative(_)))))));
    assert_eq!(study.accepted(),&before);assert!(study.optimizer().evaluations()>1);
    f.fail.set(false);study.run(1e-7,100,1500,&mut||false).unwrap();assert!(study.accepted().value<before.value);
}

#[test]
fn cancellation_inside_observation_work_is_attributed_and_retryable() {
    let f=Family::new();let mut study=TransientStudy::new(&f,&[1.4,2.2,-0.2],config(),&mut||false).unwrap();
    let before=study.accepted().clone();f.samples.set(0);
    assert!(matches!(study.run(1e-7,1,100,&mut||f.samples.get()>0),Err(SqpError::Cancelled)));
    assert_eq!(study.accepted(),&before);assert_eq!(study.optimizer().point(),before.point);
    study.run(1e-7,100,1500,&mut||false).unwrap();assert!(study.accepted().value<before.value);
}

#[test]
fn bounds_timetables_and_work_caps_refuse_without_fake_penalties() {
    let f=Family::new();let cfg=config();
    assert!(evaluate_transient(&f,&cfg,&[2.1,1.0,0.0],&mut||false).unwrap().is_none());assert_eq!(f.instances.get(),0);
    assert!(evaluate_transient(&f,&cfg,&[f64::NAN,1.0,0.0],&mut||false).is_err());assert_eq!(f.instances.get(),0);
    let mut small=cfg.clone();small.max_samples=1;
    assert!(evaluate_transient(&f,&small,&[1.0,1.0,0.0],&mut||false).is_err());assert_eq!(f.instances.get(),0);
    small=cfg.clone();small.max_kkt_dimension=8;
    assert!(evaluate_transient(&f,&small,&[1.0,1.0,0.0],&mut||false).is_err());assert_eq!(f.instances.get(),0);
    small=cfg;small.max_attempts=0;
    assert_eq!(evaluate_transient(&f,&small,&[1.0,1.0,0.0],&mut||false),Err(TransientError::ForwardStopped(RecordingStatus::AttemptLimit)));
}
