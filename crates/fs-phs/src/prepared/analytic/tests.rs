use super::*;
use crate::{Storage, QuadraticStorage, discrete_gradient};
use std::cell::Cell;
use std::rc::Rc;

struct Quartic;
impl Storage for Quartic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let s = x[0]*x[0] + 0.3*x[2]*x[2];
        0.5*(7.0*x[0]*x[0] + x[1]*x[1] + 11.0*x[2]*x[2] + x[3]*x[3]) + 20.0*s*s
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let s = x[0]*x[0] + 0.3*x[2]*x[2];
        out[0] = 7.0*x[0] + 80.0*s*x[0]; out[1] = x[1];
        out[2] = 11.0*x[2] + 24.0*s*x[2]; out[3] = x[3];
    }
}
impl Quartic {
    fn hessian_vector(&self, x: &[f64], d: &[f64], out: &mut [f64]) -> bool {
        let s = x[0]*x[0] + 0.3*x[2]*x[2];
        let ds = 2.0*x[0]*d[0] + 0.6*x[2]*d[2];
        out[0] = (7.0 + 80.0*s)*d[0] + 80.0*x[0]*ds; out[1] = d[1];
        out[2] = (11.0 + 24.0*s)*d[2] + 24.0*x[2]*ds; out[3] = d[3];
        true
    }
}
fn system(storage: Box<dyn Storage>) -> PortHamiltonian {
    let mut j = vec![0.0;16];j[1]=1.0;j[4]=-1.0;j[11]=1.0;j[14]=-1.0;
    let mut r = vec![0.0;16];r[5]=0.3;r[15]=0.5;r[7]=0.1;r[13]=0.1;
    PortHamiltonian::new(4,1,j,r,vec![0.0,1.0,0.0,0.4],storage).unwrap()
}

#[test]
fn analytic_jacobian_differentiates_the_gonzalez_correction_not_just_midpoint() {
    let sys = system(Box::new(Quartic));
    let a = [0.25,0.4,-0.13,0.1]; let b = [0.29,0.31,-0.07,0.22]; let dt = 0.03;
    let mut work = StepWorkspace::new(&sys).unwrap();work.flow.refresh(&sys).unwrap();
    work.x.copy_from_slice(&b);work.analytic_jacobian_into(&sys,&a,dt,&|x,d,o|Quartic.hessian_vector(x,d,o),None,&mut ||Ok(())).unwrap();
    let mut differs_from_midpoint = false;
    for col in 0..4 {
        let mut hi=b;let mut lo=b;let h=1e-6;hi[col]+=h;lo[col]-=h;
        let gp=discrete_gradient(sys.storage.as_ref(),&a,&hi);
        let gm=discrete_gradient(sys.storage.as_ref(),&a,&lo);
        let mid=std::array::from_fn::<f64,4,_>(|i| f64::midpoint(a[i],b[i]));
        let mut e=[0.0;4];e[col]=1.0;let mut hv=[0.0;4];Quartic.hessian_vector(&mid,&e,&mut hv);
        for row in 0..4 {
            let mut fd=if row==col {1.0}else{0.0}; let mut midpoint=fd;
            for k in 0..4 {
                let flow=sys.j[row*4+k]-sys.r[row*4+k];
                fd-=dt*flow*(gp[k]-gm[k])/(2.0*h);midpoint-=0.5*dt*flow*hv[k];
            }
            assert!((work.jacobian[row*4+col]-fd).abs()<2e-8,"{row},{col}: {fd}");
            differs_from_midpoint |= (work.jacobian[row*4+col]-midpoint).abs()>1e-5;
        }
    }
    assert!(differs_from_midpoint);
    // At exactly zero increment the owner's roundoff guard uses half H''.
    work.x.copy_from_slice(&a);work.analytic_jacobian_into(&sys,&a,dt,&|x,d,o|Quartic.hessian_vector(x,d,o),None,&mut ||Ok(())).unwrap();
    for col in 0..4 {
        let mut e=[0.0;4];e[col]=1.0;let mut hv=[0.0;4];Quartic.hessian_vector(&a,&e,&mut hv);
        for row in 0..4 {
            let expected=(if row==col {1.0}else{0.0}) - 0.5*dt*(0..4)
                .map(|k|(sys.j[row*4+k]-sys.r[row*4+k])*hv[k]).sum::<f64>();
            assert!((work.jacobian[row*4+col]-expected).abs()<1e-14);
        }
    }
}

#[test]
fn analytic_steps_preserve_nonlinear_exchange_dissipation_and_drive_work() {
    let sys=system(Box::new(Quartic));let mut fast=StepWorkspace::new(&sys).unwrap();
    let mut reference=StepWorkspace::new(&sys).unwrap();
    let mut a=[0.25,0.4,-0.13,0.1];let mut b=a;let (mut next,mut ref_next)=([0.0;4],[0.0;4]);
    let (mut y,mut ref_y)=([0.0],[0.0]);let initial=sys.hamiltonian(&a);let mut net=0.0;
    for tick in 0..200 {
        let force=if tick<50 {0.3}else{0.0};
        let record=fast.step_into_analytic(&sys,&a,&[force],0.002,&mut next,&mut y,&|x,d,o|Quartic.hessian_vector(x,d,o)).unwrap();
        reference.step_into(&sys,&b,&[force],0.002,&mut ref_next,&mut ref_y).unwrap();
        for i in 0..4 {assert!((next[i]-ref_next[i]).abs()<2e-8);}
        assert!(record.dissipated>=0.0);assert!(record.balance_residual().abs()<1e-9);
        assert!((record.supplied-0.002*force*y[0]).abs()<1e-14);
        net+=record.supplied-record.dissipated;a=next;b=ref_next;
    }
    assert!((sys.hamiltonian(&a)-initial-net).abs()<1e-8);
}

#[test]
fn analytic_cancellation_preserves_outputs_and_bit_exact_retry() {
    let sys=system(Box::new(Quartic));let a=[0.25,0.4,-0.13,0.1];
    let mut work=StepWorkspace::new(&sys).unwrap();
    let (mut next,mut y)=([123.0;4],[456.0]);let mut polls=0;
    let error=work.step_into_analytic_controlled(&sys,&a,&[0.2],0.002,&mut next,&mut y,&|x,d,o|Quartic.hessian_vector(x,d,o),
        ||{polls+=1;polls==5}).unwrap_err();
    assert_eq!(error,PreparedStepError::Cancelled);assert_eq!(next,[123.0;4]);assert_eq!(y,[456.0]);
    work.step_into_analytic(&sys,&a,&[0.2],0.002,&mut next,&mut y,&|x,d,o|Quartic.hessian_vector(x,d,o)).unwrap();
    let mut clean=StepWorkspace::new(&sys).unwrap();
    let (mut expected,mut output)=([0.0;4],[0.0]);
    clean.step_into_analytic(&sys,&a,&[0.2],0.002,&mut expected,&mut output,&|x,d,o|Quartic.hessian_vector(x,d,o)).unwrap();
    assert_eq!(next.map(f64::to_bits),expected.map(f64::to_bits));assert_eq!(y,output);
}

#[test]
fn unavailable_or_nonfinite_tangents_refuse_without_fallback_or_publication() {
    let sys=system(Box::new(Quartic));let mut work=StepWorkspace::new(&sys).unwrap();
    for available in [false,true] {
        let (mut next,mut y)=([123.0;4],[456.0]);
        let action=|_:&[f64],_:&[f64],out:&mut[f64]|{out.fill(f64::NAN);available};
        assert!(work.step_into_analytic(&sys,&[0.25,0.4,-0.13,0.1],&[0.0],0.002,
            &mut next,&mut y,&action).is_err());
        assert_eq!(next,[123.0;4]);assert_eq!(y,[456.0]);
        work.step_into(&sys,&[0.25,0.4,-0.13,0.1],&[0.0],0.002,&mut next,&mut y).unwrap();
    }
}

#[test]
fn analytic_quadratic_with_cross_terms_matches_the_original_equation() {
    let q=vec![3.0,0.5,0.5,2.0];
    let sys=PortHamiltonian::new(2,1,vec![0.0,1.0,-1.0,0.0],vec![0.0;4],vec![0.0,1.0],
        Box::new(QuadraticStorage::new(q.clone(),2).unwrap())).unwrap();
    let action=|_:&[f64],d:&[f64],out:&mut[f64]| {
        for i in 0..2 {out[i]=q[2*i]*d[0]+q[2*i+1]*d[1];}true
    };
    let mut work=StepWorkspace::new(&sys).unwrap();let (mut next,mut y)=([0.0;2],[0.0]);
    let f=work.step_into_analytic(&sys,&[0.4,-0.3],&[0.2],0.03,&mut next,&mut y,&action).unwrap();
    let expected=crate::step(&sys,&[0.4,-0.3],&[0.2],0.03).unwrap();
    for i in 0..2 {assert!((next[i]-expected.x[i]).abs()<1e-10);}
    assert!(f.balance_residual().abs()<1e-12);
    assert!((f.supplied-0.006*y[0]).abs()<1e-15);
}

struct Counted {h:Rc<Cell<usize>>,g:Rc<Cell<usize>>,hv:Rc<Cell<usize>>}
impl Storage for Counted {
    fn hamiltonian(&self,x:&[f64])->f64 {self.h.set(self.h.get()+1);Quartic.hamiltonian(x)}
    fn gradient(&self,x:&[f64],out:&mut[f64]) {self.g.set(self.g.get()+1);Quartic.gradient(x,out)}
}
impl Counted {
    fn hessian_vector(&self,x:&[f64],d:&[f64],out:&mut[f64])->bool {
        self.hv.set(self.hv.get()+1);Quartic.hessian_vector(x,d,out)
    }
}
#[test]
fn analytic_newton_eliminates_whole_energy_and_gradient_difference_probes() {
    let counts=|analytic| {
        let (h,g,hv)=(Rc::new(Cell::new(0)),Rc::new(Cell::new(0)),Rc::new(Cell::new(0)));
        let storage=Counted {h:h.clone(),g:g.clone(),hv:hv.clone()};
        let tangent=Counted {h:h.clone(),g:g.clone(),hv:hv.clone()};
        let sys=system(Box::new(storage));
        let mut work=StepWorkspace::new(&sys).unwrap();
        h.set(0);g.set(0);hv.set(0);
        let record=if analytic {
            work.step_into_analytic(&sys,&[0.25,0.4,-0.13,0.1],&[0.2],0.002,&mut [0.0;4],&mut [0.0],
                &|x,d,o|tangent.hessian_vector(x,d,o)).unwrap()
        } else {
            work.step_into(&sys,&[0.25,0.4,-0.13,0.1],&[0.2],0.002,&mut [0.0;4],&mut [0.0]).unwrap()
        };
        (h.get(),g.get(),hv.get(),record.newton_iters)
    };
    let fd=counts(false);let analytic=counts(true);
    assert_eq!(fd.2,0);assert_eq!(analytic.2,4*analytic.3);
    assert!(analytic.0<fd.0 && analytic.1<fd.1,"FD {fd:?}, analytic {analytic:?}");
    // This is a callback-work regression, deliberately not a wall-time/RT claim.
}
