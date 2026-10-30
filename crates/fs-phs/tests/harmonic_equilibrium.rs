//! Original nonlinear storage -> admitted stationary point -> shared resonance.
use fs_math::c64::C64;
use fs_phs::{DescriptorPortHamiltonian, PortHamiltonian, RelaxationBranch, Storage};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel, PoleOptions};
use fs_phs::harmonic::equilibrium::{
    EquilibriumLinearization as Linearized, EquilibriumPoint as Point,
    LinearizationBudget, LinearizationError,
};
use std::cell::Cell;

#[derive(Clone, Copy)]
struct Spring { mass:f64, stiffness:f64, cubic:f64 }
impl Storage for Spring {
    fn hamiltonian(&self,x:&[f64])->f64 {
        0.5*self.stiffness*x[0]*x[0]+0.25*self.cubic*x[0].powi(4)+0.5*x[1]*x[1]/self.mass
    }
    fn gradient(&self,x:&[f64],out:&mut[f64]) {
        out.fill(0.0);
        out[0]=self.stiffness*x[0]+self.cubic*x[0].powi(3);
        out[1]=x[1]/self.mass;
    }
}
impl Spring {
    fn hessian(&self,x:&[f64],v:&[f64],out:&mut[f64])->bool {
        out.fill(0.0);
        out[0]=(self.stiffness+3.0*self.cubic*x[0]*x[0])*v[0];
        out[1]=v[1]/self.mass;
        true
    }
    fn system(self,damping:f64)->PortHamiltonian {
        PortHamiltonian::new(2,1,vec![0.,1.,-1.,0.],vec![0.,0.,0.,damping],
            vec![0.,1.],Box::new(self)).unwrap()
    }
    fn force(&self,x:f64)->f64 {self.stiffness*x+self.cubic*x.powi(3)}
}
fn budget()->LinearizationBudget {LinearizationBudget {
    harmonic:HarmonicBudget {max_dimension:16,max_dense_bytes:2_000_000,
        max_dense_work:10_000_000,relative_residual_tolerance:1e-11},max_hessian_actions:16,
}}
fn options(scale:f64)->PoleOptions {PoleOptions {spectral_scale:scale,
    step_tolerance:1e-12,residual_tolerance:1e-12,max_iterations:32,max_dimension:16,
    max_dense_bytes:2_000_000,max_work:10_000_000}}
fn close(a:C64,b:C64,tol:f64) {assert!((a-b).abs()<tol*(1.0+b.abs()),"{a:?} vs {b:?}");}
fn point<'a>(x:&'a[f64],u:&'a[f64],scales:&'a[f64])->Point<'a> {
    Point{state:x,input:u,balance_scales:scales,relative_tolerance:1e-11}
}

#[test]
fn g1_preload_changes_frequency_through_original_nonlinear_storage() {
    let law=Spring{mass:2.5,stiffness:1600.,cubic:40_000.};
    let source=law.system(0.6);
    assert!(HarmonicModel::from_ode(&source,budget().harmonic,||false).is_err());
    for q in [0.0,0.1,0.2] {
        let x=[q,0.];let u=[law.force(q)];
        let local=Linearized::from_ode(&source,point(&x,&u,&[1.,1000.]),budget(),
            &|x,v,out|law.hessian(x,v,out),||false).unwrap();
        let k=law.stiffness+3.0*law.cubic*q*q;
        let gamma=0.6/(2.0*law.mass);
        let exact=C64::new(-gamma,-(k/law.mass-gamma*gamma).sqrt());
        let pole=local.model().resonance_pole(exact+C64::new(0.01,0.1),
            &[C64::ONE,exact.scale(law.mass)],&options(60.),||false).unwrap();
        close(pole.laplace_pole_per_s,exact,1e-10);
        assert_eq!(local.equilibrium_state(),&x);
        assert_eq!(local.equilibrium_input(),&u);
        assert_eq!(local.equilibrium_output(),&[0.]);
        assert_eq!(local.equilibrium_stored_energy_j(),law.hamiltonian(&x));
        assert!(local.maximum_scaled_residual()<1e-11);
        for w in [15.,40.,70.] {
            let f=C64::new(1.7,-0.3);
            let response=local.model().factor(w,||false).unwrap().response(&[f],||false).unwrap();
            let expected=f/C64::new(k-law.mass*w*w,-0.6*w);
            close(response.state[0],expected,1e-12);
            close(response.output[0],C64::new(0.,-w)*expected,1e-12);
        }
    }
}

#[test]
fn g3_local_linear_response_is_the_limit_of_original_nonlinear_gradient() {
    let law=Spring{mass:2.5,stiffness:1600.,cubic:40_000.};let source=law.system(0.6);
    let x=[0.15,0.];let u=[law.force(x[0])];
    let local=Linearized::from_ode(&source,point(&x,&u,&[1.,1000.]),budget(),
        &|x,v,out|law.hessian(x,v,out),||false).unwrap();
    let q=local.tangent_storage_matrix();
    let mut errors=Vec::new();
    for epsilon in [1e-3,5e-4,2.5e-4] {
        let shifted=[x[0]+epsilon,epsilon*0.2];let mut actual=[0.;2];let mut original=[0.;2];
        law.gradient(&shifted,&mut actual);law.gradient(&x,&mut original);
        let tangent=q[0]*epsilon+q[1]*epsilon*0.2;
        errors.push((actual[0]-original[0]-tangent).abs());
    }
    assert!(errors[0]/errors[1]>3.9 && errors[1]/errors[2]>3.9,"{errors:?}");
    // An exactly quadratic source retains the existing complete AC semantics.
    let exact=fs_phs::mass_spring_damper(2.5,1600.,0.6).unwrap();
    let reference=HarmonicModel::from_ode(&exact,budget().harmonic,||false).unwrap();
    let linear=Spring{cubic:0.,..law};let bx=[0.2,0.];let bu=[320.];
    let tangent=Linearized::from_ode(&exact,point(&bx,&bu,&[1.,1000.]),budget(),
        &|x,v,out|linear.hessian(x,v,out),||false).unwrap();
    let a=reference.factor(30.,||false).unwrap().response(&[C64::ONE],||false).unwrap();
    let b=tangent.model().factor(30.,||false).unwrap().response(&[C64::ONE],||false).unwrap();
    assert_eq!(a.state,b.state);assert_eq!(a.output,b.output);
}

#[test]
fn g1_relaxed_nonlinear_material_keeps_oscillation_and_internal_relaxation() {
    let law=Spring{mass:1.,stiffness:30_000.,cubic:1_000_000.};let kb=50_000.0_f64;
    let source=law.system(0.).with_relaxation_branches(vec![RelaxationBranch {
        projection:vec![1.,0.],stiffness:kb,relaxation_time_s:1.0/300.,
    }]).unwrap();
    let x=[0.1,0.,kb.sqrt()*0.1];let u=[law.force(x[0])];
    let local=Linearized::from_ode(&source,point(&x,&u,&[1.,10_000.,100.]),budget(),
        &|x,v,out| {law.hessian(x,v,out);let e=kb.sqrt()*v[0]-v[2];
            out[0]+=kb.sqrt()*e;out[2]=-e;true},||false).unwrap();
    for s in [C64::new(-50.,-87_500.0_f64.sqrt()),C64::from_re(-200.)] {
        let z=C64::from_re(kb.sqrt())/(C64::ONE+s.scale(1.0/300.));
        let p=local.model().resonance_pole(s,&[C64::ONE,s,z],&options(400.),||false).unwrap();
        close(p.laplace_pole_per_s,s,1e-10);assert_eq!(p.state.len(),3);
        assert!(p.decay_power_residual_w.abs()<1e-8*(1.0+p.dissipated_power_w.abs()));
    }
}

#[test]
fn g1_descriptor_keeps_original_algebraic_efforts_and_resistive_loading() {
    let law=Spring{mass:2.5,stiffness:1600.,cubic:40_000.};
    // Algebraic resistive port: 0=-b*(p/m)-r_a*lambda. No stored multiplier.
    let (b,ra,d)=(2.0,4.0,0.6);
    let source=DescriptorPortHamiltonian::new(3,2,1,
        vec![0.,1.,0.,-1.,0.,b,0.,-b,0.],vec![0.,0.,0.,0.,d,0.,0.,0.,ra],
        vec![0.,1.,0.],Box::new(law)).unwrap();
    let x=[0.1,0.,0.];let u=[law.force(x[0])];
    let local=Linearized::from_descriptor(&source,point(&x,&u,&[1.,1000.,1.]),budget(),
        &|x,v,out|law.hessian(x,v,out),||false).unwrap();
    let k=2800.;let gamma=(d+b*b/ra)/(2.0*law.mass);
    let s=C64::new(-gamma,-(k/law.mass-gamma*gamma).sqrt());
    let p=local.model().resonance_pole(s,&[C64::ONE,s.scale(law.mass),s.scale(-b/ra)],
        &options(40.),||false).unwrap();
    close(p.laplace_pole_per_s,s,1e-10);
    close(p.state[2],p.state[1].scale(-b/(ra*law.mass)),1e-12);
    let mut wrong=x;wrong[2]=0.001;
    assert!(matches!(Linearized::from_descriptor(&source,point(&wrong,&u,&[1.,1000.,1.]),budget(),
        &|x,v,out|law.hessian(x,v,out),||false),Err(LinearizationError::NotEquilibrium{..})));
}

#[test]
fn g1_negative_curvature_is_not_repaired_into_a_stable_resonance() {
    let law=Spring{mass:1.,stiffness:-4.,cubic:1.};let source=law.system(0.4);
    let local=Linearized::from_ode(&source,point(&[0.,0.],&[0.],&[1.,1.]),budget(),
        &|x,v,out|law.hessian(x,v,out),||false).unwrap();
    assert_eq!(local.tangent_storage_matrix()[0],-4.);
    let growth=(-0.4+(0.16_f64+16.).sqrt())*0.5;
    let p=local.model().resonance_pole(C64::from_re(growth),&[C64::ONE,C64::from_re(growth)],
        &options(3.),||false).unwrap();
    close(p.laplace_pole_per_s,C64::from_re(growth),1e-12);
    assert!(p.laplace_pole_per_s.re>0.);
}

#[test]
fn g0_g4_equilibrium_source_shape_hessian_and_cancellation_refusals_are_atomic() {
    let law=Spring{mass:2.5,stiffness:1600.,cubic:40_000.};let source=law.system(0.6);
    let x=[0.1,0.];let u=[law.force(x[0])];let scales=[1.,1000.];
    let p=point(&x,&u,&scales);let seen=Cell::new(0);
    let action=|x:&[f64],v:&[f64],out:&mut[f64]| {seen.set(seen.get()+1);law.hessian(x,v,out)};
    let mut short=budget();short.max_hessian_actions=1;
    assert!(Linearized::from_ode(&source,p,short,&action,||false).is_err());assert_eq!(seen.get(),0);
    assert!(Linearized::from_ode(&source,point(&x,&[0.],&scales),budget(),&action,||false).is_err());
    assert_eq!(seen.get(),0);
    assert!(Linearized::from_ode(&source,p,budget(),&|_,_,_|true,||false).is_err());
    assert!(Linearized::from_ode(&source,p,budget(),&|_,_,_|false,||false).is_err());
    assert!(Linearized::from_ode(&source,p,budget(),&|_,v,o| {o[0]=v[0]+v[1];o[1]=v[1];true},||false).is_err());
    let polls=Cell::new(0);
    let reference=Linearized::from_ode(&source,p,budget(),&action,||{polls.set(polls.get()+1);false}).unwrap();
    let total=polls.get();
    for cut in [1,total/2,total] {
        let count=Cell::new(0);
        assert!(Linearized::from_ode(&source,p,budget(),&action,||{count.set(count.get()+1);count.get()>=cut}).is_err());
    }
    let replay=Linearized::from_ode(&source,p,budget(),&action,||false).unwrap();
    assert_eq!(reference.tangent_storage_matrix(),replay.tangent_storage_matrix());
    assert_eq!(x,[0.1,0.]);assert_eq!(u,[law.force(0.1)]);
}
