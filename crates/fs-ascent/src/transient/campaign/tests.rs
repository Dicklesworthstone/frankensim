use super::*;
use crate::transient::{OdeVjp, SampleObjective, TransientModel};
use fs_time::PiController;
use fs_time::adaptive::adjoint::trajectory::{RecordingConfig, ReplayBudget};
use std::cell::Cell;
use std::rc::Rc;

const BOX: [[f64;2];2] = [[0.1,2.0],[-1.0,1.0]];
struct Family { amplitude: f64, times: Vec<f64>, fail: Rc<Cell<bool>>, calls: Rc<Cell<usize>> }
impl Family {
    fn new(amplitude:f64,times:Vec<f64>)->Self {
        Self { amplitude,times,fail:Rc::new(Cell::new(false)),calls:Rc::new(Cell::new(0)) }
    }
}
struct Model { initial:[f64;1], p:[f64;2], targets:Vec<f64> }
impl OdeVjp for Model {
    fn dimension(&self)->usize {1} fn parameter_count(&self)->usize {2}
    fn rhs(&self,_:f64,x:&[f64],out:&mut[f64]) {out[0]=-self.p[0]*x[0];}
    fn rhs_vjp(&self,_:f64,x:&[f64],b:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        xb[0]=-self.p[0]*b[0];pb[0]=-x[0]*b[0];pb[1]=0.0;Ok(())
    }
}
impl SampleObjective for Model {
    fn evaluate(&self,i:usize,_:f64,x:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<f64,String> {
        let r=x[0]+self.p[1]-self.targets[i];xb[0]=r;pb[0]=0.0;pb[1]=r;Ok(0.5*r*r)
    }
}
impl TransientModel for Model {
    fn initial_values(&self)->&[f64] {&self.initial}
    fn initial_vjp(&self,_:&[f64],p:&mut[f64])->Result<(),String> {p.fill(0.0);Ok(())}
}
impl TransientFamily for Family {
    type Model=Model;
    fn bounds(&self)->&[[f64;2]] {&BOX}
    fn sample_times(&self)->&[f64] {&self.times}
    fn instantiate(&self,p:&[f64])->Result<Model,String> {
        self.calls.set(self.calls.get()+1);
        if self.fail.get() {return Err("injected case failure".into());}
        Ok(Model {initial:[self.amplitude],p:[p[0],p[1]],targets:self.times.iter().map(|t|self.amplitude*(-0.7*t).exp()+0.15).collect()})
    }
}
fn case(id:u64,family:&Family,weight:f64)->Experiment<'_,Family> {
    Experiment {id,family,weight,config:TransientConfig {start:0.0,initial_step:0.1,
        recording:RecordingConfig{end:*family.times.last().unwrap(),rtol:1e-10,atol:1e-12,controller:PiController::default(),max_workspace_components:4096},
        max_state_components:1,max_samples:16,max_attempts:10000,max_records:10000,
        replay:ReplayBudget{checkpoints:64,replayed_steps:100000},max_kkt_dimension:6}}
}
fn families()->[Family;2] {
    [Family::new(1.2,vec![0.0,0.2,1.0,2.0]),Family::new(2.3,vec![0.1,0.5,1.4,3.0])]
}
fn analytic(f:&Family,p:&[f64])->(f64,[f64;2]) {
    let mut value=0.0;let mut grad=[0.0;2];
    for &t in &f.times {
        let x=f.amplitude*(-p[0]*t).exp();let r=x+p[1]-(f.amplitude*(-0.7*t).exp()+0.15);
        value+=0.5*r*r;grad[0]+=-t*x*r;grad[1]+=r;
    }
    (value,grad)
}

#[test]
fn joint_gradients_and_experiment_permutations_preserve_the_same_fit() {
    let f=families();let p=[1.1,-0.2];
    let a=TransientCampaign::new(vec![case(42,&f[1],0.3),case(7,&f[0],1.2)],2,14).unwrap();
    let b=TransientCampaign::new(vec![case(7,&f[0],1.2),case(42,&f[1],0.3)],2,14).unwrap();
    let mut work=CampaignControl::new(10,20);
    let got=a.evaluate(&p,&mut work,&mut||false).unwrap().unwrap();
    assert_eq!(got,b.evaluate(&p,&mut work,&mut||false).unwrap().unwrap());
    let (v0,g0)=analytic(&f[0],&p);let (v1,g1)=analytic(&f[1],&p);
    assert!((got.value-(1.2*v0+0.3*v1)).abs()<1e-8);
    for j in 0..2 {assert!((got.gradient[j]-(1.2*g0[j]+0.3*g1[j])).abs()<1e-8);}
    assert_eq!(got.experiments.iter().map(|e|e.id).collect::<Vec<_>>(),vec![7,42]);
    for e in &got.experiments {assert_eq!(e.result.point,got.point);}
}

#[test]
fn complete_trial_preflight_and_extension_do_not_reset_spent_work() {
    let f=families();let campaign=TransientCampaign::new(vec![case(1,&f[0],1.0),case(2,&f[1],1.0)],2,14).unwrap();
    let mut work=CampaignControl::new(1,1);
    assert_eq!(campaign.evaluate(&[1.1,0.0],&mut work,&mut||false),Err(CampaignError::ExperimentLimit));
    assert_eq!(work.work(),CampaignWork{trials:0,experiments:0});assert_eq!(f[0].calls.get(),0);
    work.extend(1,2).unwrap();campaign.evaluate(&[1.1,0.0],&mut work,&mut||false).unwrap();
    assert_eq!(work.work(),CampaignWork{trials:1,experiments:2});
    assert_eq!(campaign.evaluate(&[1.1,0.0],&mut work,&mut||false),Err(CampaignError::TrialLimit));
    work.extend(2,4).unwrap();campaign.evaluate(&[1.1,0.0],&mut work,&mut||false).unwrap();
    assert_eq!(work.work(),CampaignWork{trials:2,experiments:4});assert!(work.extend(1,4).is_err());
}

#[test]
fn failed_second_experiment_retains_prior_acceptance_and_attempt_counts() {
    let f=families();let campaign=TransientCampaign::new(vec![case(1,&f[0],1.0),case(2,&f[1],1.0)],2,14).unwrap();
    let mut work=CampaignControl::new(2000,4000);
    let mut study=CampaignStudy::new(&campaign,&[1.2,-0.3],&mut work,&mut||false).unwrap();
    let before=study.accepted().clone();let spent=study.work();f[1].fail.set(true);
    assert!(matches!(study.run(1e-7,1,100,&mut||false),Err(SqpError::Evaluation(CampaignError::Experiment{id:2,..}))));
    assert_eq!(study.accepted(),&before);assert_eq!(study.optimizer().point(),before.point);
    assert_eq!(study.work(),CampaignWork{trials:spent.trials+1,experiments:spent.experiments+2});
    f[1].fail.set(false);assert_eq!(study.run(1e-7,100,1500,&mut||false).unwrap().stop,SqpStop::Converged);
    assert!(study.accepted().value<before.value*1e-8);
}

#[test]
fn shared_parameter_sqp_resumes_without_repeating_accepted_campaigns() {
    let f=families();let campaign=TransientCampaign::new(vec![case(1,&f[0],1.0),case(2,&f[1],0.3)],2,14).unwrap();
    let (mut wa,mut wb)=(CampaignControl::new(2000,4000),CampaignControl::new(2000,4000));
    let mut a=CampaignStudy::new(&campaign,&[1.2,-0.3],&mut wa,&mut||false).unwrap();
    let mut b=CampaignStudy::new(&campaign,&[1.2,-0.3],&mut wb,&mut||false).unwrap();
    let before=b.work();assert_eq!(b.run(1e-7,10,1,&mut||false).unwrap().stop,SqpStop::EvaluationLimit);assert_eq!(b.work(),before);
    assert_eq!(a.run(1e-7,100,1500,&mut||false).unwrap().stop,SqpStop::Converged);
    for _ in 0..100 {if b.run(1e-7,1,1500,&mut||false).unwrap().stop!=SqpStop::IterationLimit {break;}}
    assert_eq!(a.accepted(),b.accepted());assert_eq!(a.work(),b.work());assert_eq!(a.optimizer().history(),b.optimizer().history());
    for (p,target) in a.optimizer().point().iter().zip([0.7,0.15]) {assert!((p-target).abs()<2e-6);}
}

#[test]
fn cancellation_charges_attempted_cases_but_publishes_no_partial_family() {
    let f=families();let campaign=TransientCampaign::new(vec![case(1,&f[0],1.0),case(2,&f[1],1.0)],2,14).unwrap();
    let mut work=CampaignControl::new(2000,4000);
    let mut study=CampaignStudy::new(&campaign,&[1.2,-0.3],&mut work,&mut||false).unwrap();
    let before=study.accepted().clone();let calls=f[1].calls.get();let spent=study.work();
    assert!(matches!(study.run(1e-7,1,100,&mut||f[1].calls.get()>calls),Err(SqpError::Cancelled)));
    assert_eq!(study.accepted(),&before);assert_eq!(study.work(),CampaignWork{trials:spent.trials+1,experiments:spent.experiments+2});
    assert_eq!(study.run(1e-7,100,1500,&mut||false).unwrap().stop,SqpStop::Converged);
}

#[test]
fn admission_limits_and_duplicate_ids_fail_before_any_factory() {
    let f=families();
    assert!(TransientCampaign::new(vec![case(1,&f[0],1.0),case(1,&f[1],1.0)],2,14).is_err());
    assert!(TransientCampaign::new(vec![case(1,&f[0],1.0),case(2,&f[1],1.0)],2,13).is_err());
    assert!(TransientCampaign::new(vec![case(1,&f[0],0.0)],1,9).is_err());
    assert!(TransientCampaign::new(vec![case(1,&f[0],1.0)],0,9).is_err());
    let campaign=TransientCampaign::new(vec![case(1,&f[0],1.0)],1,9).unwrap();let mut work=CampaignControl::new(1,1);
    assert!(campaign.evaluate(&[2.5,0.0],&mut work,&mut||false).unwrap().is_none());
    assert_eq!(work.work(),CampaignWork{trials:0,experiments:0});assert_eq!(f[0].calls.get(),0);
}
