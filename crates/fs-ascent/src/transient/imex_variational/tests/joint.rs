use super::*;
use crate::transient::variational::joint::{ParameterFamily, JointWindow, JointWindowStudy};
use fs_time::stiff::{Imex2, imex2_step};
use std::rc::Rc;

struct Family { data: Rc<Vec<f64>>, calls: Cell<usize>, fail: Cell<bool> }
struct Model { thermal: Thermal, bias: f64, data: Rc<Vec<f64>> }
impl LinearOp for Model {
    fn n(&self)->usize {2}
    fn apply(&self,x:&[f64],out:&mut[f64]) {self.thermal.apply(x,out);}
    fn apply_transpose(&self,x:&[f64],out:&mut[f64]) {self.thermal.apply_transpose(x,out);}
}
impl ImexVjp for Model {
    fn parameter_count(&self)->usize {2}
    fn nonlinear(&self,x:&[f64],out:&mut[f64]) {self.thermal.nonlinear(x,out);}
    fn nonlinear_vjp(&self,x:&[f64],b:&[f64],xb:&mut[f64],pb:&mut[f64])->Result<(),String> {
        self.thermal.nonlinear_vjp(x,b,xb,pb)?;pb[1]=0.0;Ok(())
    }
    fn linear_parameter_vjp(&self,_:&[f64],_:&[f64],pb:&mut[f64])->Result<(),String> {pb.fill(0.0);Ok(())}
}
impl WindowObjective for Model {
    fn evaluate(&self,_:&[f64],_:usize,x:&[f64],bar:&mut[f64],_:&mut dyn FnMut()->bool)->Result<f64,String> {
        let mut value=0.0;
        for ((x,y),b) in x.iter().zip(self.data.iter()).zip(bar) {
            let r=(x+self.bias-y)/0.1;*b=r/0.1;value+=0.5*r*r;
        }
        Ok(value)
    }
    fn parameter_partials(&self,_:&[f64],_:usize,x:&[f64],pb:&mut[f64],_:&mut dyn FnMut()->bool)->Result<(),String> {
        pb[0]=0.0;pb[1]=x.iter().zip(self.data.iter()).map(|(x,y)|(x+self.bias-y)/0.01).sum();Ok(())
    }
}
impl ParameterFamily for Family {
    type Model=Model;
    fn instantiate(&self,point:&[f64],_:&mut dyn FnMut()->bool)->Result<Model,String> {
        self.calls.set(self.calls.get()+1);
        if self.fail.get() {return Err("injected joint factory refusal".into());}
        Ok(Model{thermal:Thermal::new(point[0]),bias:point[1],data:self.data.clone()})
    }
}
fn dense_advance(x:&[f64],source:f64)->Vec<f64> {
    let method=Imex2::new(&[-401.0,400.0,4.0,-5.0],2,0.125);
    let mut x=x.to_vec();
    for _ in 0..2 {imex2_step(&method,&mut x,&|_,out| {out[0]=0.0;out[1]=source;});}
    x
}
fn family()->Family {
    let mut states=vec![1.0,1.0];let mut x=vec![1.0,1.0];
    for _ in 0..3 {x=dense_advance(&x,2.0);states.extend_from_slice(&x);}
    Family{data:Rc::new(states.iter().map(|x|x+0.12).collect()),calls:Cell::new(0),fail:Cell::new(false)}
}
// Independently assembled linear Gaussian least-squares oracle. Columns are
// exact linear residual actions; no adjoint or gradient-difference code is used.
fn oracle(data:&[f64])->Vec<f64> {
    let residual=|z:&[f64]| {
        let x=(0..8).map(|i|1.0+[1.0,2.0][i%2]*z[i]).collect::<Vec<_>>();
        let source=2.0*z[8];let bias=0.25*z[9];
        let mut r=x.iter().zip(data).map(|(x,y)|(x+bias-y)/0.1).collect::<Vec<_>>();
        r.extend(x[..2].iter().map(|x|(x-1.0)/0.5));
        for k in 0..3 {
            let y=dense_advance(&x[2*k..2*k+2],source);
            r.extend((0..2).map(|i|(x[2*k+2+i]-y[i])/0.03));
        }
        r.push(source/10.0);r.push(bias);r
    };
    let r0=residual(&[0.0;10]);let mut columns=Vec::new();
    for i in 0..10 {let mut z=[0.0;10];z[i]=1.0;columns.push(residual(&z).iter().zip(&r0).map(|(a,b)|a-b).collect::<Vec<_>>());}
    let mut h=vec![0.0;100];let mut rhs=vec![0.0;10];
    for i in 0..10 {
        rhs[i]=-columns[i].iter().zip(&r0).map(|(a,b)|a*b).sum::<f64>();
        for j in 0..10 {h[i*10+j]=columns[i].iter().zip(&columns[j]).map(|(a,b)|a*b).sum();}
    }
    fs_la::factor::lu(&h,10).unwrap().solve(&mut rhs);rhs
}

#[test]
fn implicit_joint_gradient_includes_heating_and_direct_sensor_parameters() {
    let id=IdentityPreconditioner;let f=family();
    let p=ImexWindowPolicy::new(0.0,&[2,2,2],config(),&id,&id,&mut||false).unwrap();
    let w=make_window(p.times(),0.03);let joint=JointWindow::new(&w,&f,&[0.0,0.0],&[2.0,0.25],&[10.0,1.0],6).unwrap();
    let z=[0.1,-0.03,0.05,0.01,-0.07,0.08,0.03,0.02,0.4,0.1];
    let got=joint.evaluate_using(&z,&p,&mut control(),&mut||false).unwrap();
    assert_eq!(f.calls.get(),1);assert_eq!(got.window.accepted_steps,6);
    for i in 0..10 {
        let h=1e-5;let (mut plus,mut minus)=(z,z);plus[i]+=h;minus[i]-=h;
        let a=joint.evaluate_using(&plus,&p,&mut control(),&mut||false).unwrap().value;
        let b=joint.evaluate_using(&minus,&p,&mut control(),&mut||false).unwrap().value;
        assert!((got.gradient[i]-(a-b)/(2.0*h)).abs()<5e-5,"coordinate {i}");
    }
}

#[test]
fn implicit_joint_study_matches_dense_MAP_and_checkpoint_continuation() {
    let id=IdentityPreconditioner;let f=family();let expected=oracle(&f.data);
    let p=ImexWindowPolicy::new(0.0,&[2,2,2],config(),&id,&id,&mut||false).unwrap();
    let w=make_window(p.times(),0.03);let joint=JointWindow::new(&w,&f,&[0.0,0.0],&[2.0,0.25],&[10.0,1.0],6).unwrap();
    let mut work=control();let mut settings=settings();settings.gradient_tolerance=1e-5;
    let mut study=JointWindowStudy::new(&joint,&[0.0;10],p,settings,&mut work,&mut||false).unwrap();
    let mut split=study.clone();let mut split_work=control();
    assert_eq!(study.run(600,&mut work,&mut||false).unwrap().reason,crate::StopReason::GradNorm);
    for _ in 0..600 {if split.run(1,&mut split_work,&mut||false).unwrap().reason!=crate::StopReason::IterationCap {break;}}
    for (a,b) in study.accepted().controls.iter().zip(expected) {assert!((a-b).abs()<2e-4,"{a} != {b}");}
    assert_eq!(study.accepted(),split.accepted());assert_eq!(study.optimizer().history,split.optimizer().history);
    assert!(study.accepted().window.observation_value<0.1);
}

#[test]
fn joint_factory_failure_and_cancellation_keep_the_accepted_implicit_result() {
    let id=IdentityPreconditioner;let f=family();
    let p=ImexWindowPolicy::new(0.0,&[2,2,2],config(),&id,&id,&mut||false).unwrap();
    let w=make_window(p.times(),0.03);let joint=JointWindow::new(&w,&f,&[0.0,0.0],&[2.0,0.25],&[10.0,1.0],6).unwrap();
    let mut work=control();let mut study=JointWindowStudy::new(&joint,&[0.0;10],p,settings(),&mut work,&mut||false).unwrap();
    let accepted=study.accepted().clone();let spent=work.evaluations();
    f.calls.set(0);
    assert!(matches!(study.run(1,&mut work,&mut||f.calls.get()>0),Err(crate::LbfgsError::Evaluation(WindowError::Cancelled))));
    assert_eq!(study.accepted(),&accepted);assert_eq!(work.evaluations(),spent+1);
    f.fail.set(true);
    assert!(matches!(study.run(1,&mut work,&mut||false),Err(crate::LbfgsError::Evaluation(WindowError::Model(_)))));
    assert_eq!(study.accepted(),&accepted);assert_eq!(work.evaluations(),spent+2);
    f.fail.set(false);assert!(study.run(50,&mut work,&mut||false).unwrap().f<accepted.value);
}
