use super::*;
use crate::{Storage, QuadraticStorage, discrete_gradient};

fn oscillator(r: f64, storage: Box<dyn Storage>) -> PortHamiltonian {
    PortHamiltonian::new(2,1,vec![0.0,1.0,-1.0,0.0],vec![0.0,0.0,0.0,r],
        vec![0.0,1.0],storage).unwrap()
}
fn quadratic() -> Box<dyn Storage> {
    Box::new(QuadraticStorage::new(vec![9.0,0.0,0.0,1.0],2).unwrap())
}
struct Quartic;
impl Storage for Quartic {
    fn hamiltonian(&self,x:&[f64])->f64 {0.5*(9.0*x[0]*x[0]+x[1]*x[1])+5.0*x[0].powi(4)}
    fn gradient(&self,x:&[f64],out:&mut[f64]) {out[0]=9.0*x[0]+20.0*x[0].powi(3);out[1]=x[1];}
}
fn hessian(x:&[f64],d:&[f64],out:&mut[f64])->bool {
    out[0]=(9.0+60.0*x[0]*x[0])*d[0];out[1]=d[1];true
}
// A smooth spatially varying nonlinear resistor, not a new contact law.
fn loss(x:&[f64],e:&[f64],out:&mut[f64])->bool {
    out[0]=0.0;out[1]=(0.4+x[0]*x[0])*(e[1]+e[1].powi(3));true
}
fn tangent(x:&[f64],e:&[f64],dx:&[f64],de:&[f64],out:&mut[f64])->bool {
    out[0]=0.0;
    out[1]=2.0*x[0]*dx[0]*(e[1]+e[1].powi(3))
        +(0.4+x[0]*x[0])*(1.0+3.0*e[1]*e[1])*de[1];true
}

#[test]
fn implicit_nonlinear_port_reproduces_the_static_resistance_equation() {
    let sys=oscillator(0.0,quadratic());let reference=oscillator(0.7,quadratic());
    let mut work=StepWorkspace::new(&sys).unwrap();let mut other=StepWorkspace::new(&reference).unwrap();
    let (mut next,mut expected,mut y,mut ey)=([0.0;2],[0.0;2],[0.0],[0.0]);
    let force=|_:&[f64],e:&[f64],out:&mut[f64]| {out[0]=0.0;out[1]=0.7*e[1];true};
    for dt in [0.0,0.003,0.03] {
        let a=work.step_into_dissipative_controlled(&sys,&[0.3,0.7],&[0.2],dt,
            &mut next,&mut y,&force,None,None,||false).unwrap();
        let b=other.step_into(&reference,&[0.3,0.7],&[0.2],dt,&mut expected,&mut ey).unwrap();
        for i in 0..2 {assert!((next[i]-expected[i]).abs()<1e-10);}
        assert!((a.dissipated-b.dissipated).abs()<1e-11);
        assert!((a.supplied-b.supplied).abs()<1e-11);assert!(a.balance_residual().abs()<1e-10);
    }
}

#[test]
fn analytic_dissipative_jacobian_includes_state_and_complete_effort_derivatives() {
    let sys=oscillator(0.2,Box::new(Quartic));let a=[0.25,0.4];let b=[0.29,0.31];let dt=0.03;
    let mut work=StepWorkspace::new(&sys).unwrap();work.flow.refresh(&sys).unwrap();
    work.x.copy_from_slice(&b);
    work.analytic_jacobian_into(&sys,&a,dt,&hessian,
        Some(Dissipation {force:&loss,tangent:Some(&tangent)}),&mut ||Ok(())).unwrap();
    let residual=|x:&[f64]| {
        let e=discrete_gradient(sys.storage.as_ref(),&a,x);
        let mid=[f64::midpoint(a[0],x[0]),f64::midpoint(a[1],x[1])];
        let mut d=[0.0;2];loss(&mid,&e,&mut d);
        [x[0]-a[0]-dt*e[1],x[1]-a[1]-dt*(-e[0]-0.2*e[1]-d[1])]
    };
    for col in 0..2 {
        let mut hi=b;let mut lo=b;let h=1e-6;hi[col]+=h;lo[col]-=h;
        let hi=residual(&hi);let lo=residual(&lo);
        for row in 0..2 {assert!((work.jacobian[row*2+col]-(hi[row]-lo[row])/(2.0*h)).abs()<2e-8);}
    }
}

#[test]
fn nonlinear_passive_port_preserves_full_work_balance_and_fd_analytic_parity() {
    let sys=oscillator(0.2,Box::new(Quartic));let mut fd=StepWorkspace::new(&sys).unwrap();
    let mut analytic=StepWorkspace::new(&sys).unwrap();let (mut a,mut b)=([0.25,0.4],[0.25,0.4]);
    let (mut next,mut expected,mut y,mut ey)=([0.0;2],[0.0;2],[0.0],[0.0]);
    let initial=sys.hamiltonian(&a);let mut net=0.0;
    for tick in 0..200 {
        let input=if tick<50 {0.3}else{0.0};
        let f=analytic.step_into_dissipative_controlled(&sys,&a,&[input],0.002,&mut next,&mut y,
            &loss,Some(&hessian),Some(&tangent),||false).unwrap();
        fd.step_into_dissipative_controlled(&sys,&b,&[input],0.002,&mut expected,&mut ey,
            &loss,None,None,||false).unwrap();
        for i in 0..2 {assert!((next[i]-expected[i]).abs()<2e-8);}
        assert!(f.dissipated>=0.0);assert!(f.balance_residual().abs()<1e-9);
        assert!((f.supplied-0.002*input*y[0]).abs()<1e-14);
        net+=f.supplied-f.dissipated;a=next;b=expected;
    }
    assert!((sys.hamiltonian(&a)-initial-net).abs()<1e-8);
}

#[test]
fn bad_dissipation_cancellation_and_missing_tangent_preserve_outputs_and_retry() {
    let sys=oscillator(0.0,Box::new(Quartic));let mut work=StepWorkspace::new(&sys).unwrap();
    let state=[0.25,0.4];let (mut next,mut y)=([123.0;2],[456.0]);
    for kind in 0..3 {
        let bad=|_:&[f64],e:&[f64],o:&mut[f64]| {
            o[1]=if kind==0 {-e[1]}else{f64::NAN};kind!=2
        };
        assert!(work.step_into_dissipative_controlled(&sys,&state,&[0.2],0.002,&mut next,&mut y,
            &bad,None,None,||false).is_err());
        assert_eq!(next,[123.0;2]);assert_eq!(y,[456.0]);
    }
    assert!(work.step_into_dissipative_controlled(&sys,&state,&[0.2],0.002,&mut next,&mut y,
        &loss,Some(&hessian),None,||false).is_err());
    let mut polls=0;
    let error=work.step_into_dissipative_controlled(&sys,&state,&[0.2],0.002,&mut next,&mut y,
        &loss,Some(&hessian),Some(&tangent),||{polls+=1;polls==4}).unwrap_err();
    assert_eq!(error,PreparedStepError::Cancelled);assert_eq!(next,[123.0;2]);assert_eq!(y,[456.0]);
    work.step_into_dissipative_controlled(&sys,&state,&[0.2],0.002,&mut next,&mut y,
        &loss,Some(&hessian),Some(&tangent),||false).unwrap();
    let mut clean=StepWorkspace::new(&sys).unwrap();let (mut expected,mut ey)=([0.0;2],[0.0]);
    clean.step_into_dissipative_controlled(&sys,&state,&[0.2],0.002,&mut expected,&mut ey,
        &loss,Some(&hessian),Some(&tangent),||false).unwrap();
    assert_eq!(next,expected);assert_eq!(y,ey);
}
