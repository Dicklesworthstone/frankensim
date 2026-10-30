//! Original nonlinear loss, complete local equations, and independent references.
use fs_math::c64::C64;
use fs_phs::{DescriptorPortHamiltonian, PortHamiltonian, QuadraticStorage};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel, PoleOptions};
use fs_phs::harmonic::equilibrium::{EquilibriumLinearization as Local, EquilibriumPoint,
    LinearizationBudget, LinearizationError, NonlinearDissipation};
use fs_phs::harmonic::response::{HarmonicParameterDirection, HarmonicPoleResponseOptions};
use fs_la::pole_response::PoleResponseOptions;
use std::cell::Cell;

fn budget(n: usize) -> LinearizationBudget {
    LinearizationBudget { harmonic: HarmonicBudget { max_dimension: 16,
        max_dense_bytes: 2_000_000, max_dense_work: 10_000_000,
        relative_residual_tolerance: 1e-11 }, max_hessian_actions: n }
}
fn source() -> PortHamiltonian {
    PortHamiltonian::new(2,2,vec![0.,1.,-1.,0.],vec![0.,0.,0.,0.1],
        vec![1.,0.,0.,1.],Box::new(QuadraticStorage::new(vec![2.,0.,0.,0.5],2).unwrap())).unwrap()
}
fn hessian(_: &[f64], d: &[f64], out: &mut [f64]) -> bool {
    out[0]=2.*d[0]; out[1]=0.5*d[1]; true
}
fn options() -> PoleOptions {
    PoleOptions { spectral_scale: 2., step_tolerance: 1e-13, residual_tolerance: 1e-13,
        max_iterations: 32, max_dimension: 16, max_dense_bytes: 2_000_000, max_work: 10_000_000 }
}
fn point<'a>(x: &'a [f64], u: &'a [f64], scales: &'a [f64]) -> EquilibriumPoint<'a> {
    EquilibriumPoint { state:x, input:u, balance_scales:scales, relative_tolerance:1e-12 }
}
fn near(a: C64, b: C64) { assert!((a-b).abs()<1e-10*(1.+b.abs()),"{a:?} != {b:?}"); }
fn cubic(beta: f64) -> Local {
    let source=source(); let x=[0.2,0.6]; let e=[0.4,0.3];
    let flow=|x: &[f64],e: &[f64],out: &mut [f64]| {
        out[0]=(1.+beta*x[1]*x[1])*(0.1*e[0]+0.2*e[0].powi(3)); out[1]=0.; true
    };
    let tangent=|x: &[f64],e: &[f64],dx: &[f64],de: &[f64],out: &mut [f64]| {
        out[0]=2.*beta*x[1]*(0.1*e[0]+0.2*e[0].powi(3))*dx[1]
            +(1.+beta*x[1]*x[1])*(0.1+0.6*e[0]*e[0])*de[0]; out[1]=0.; true
    };
    let mut d=[0.;2]; flow(&x,&e,&mut d);
    let u=[d[0]-e[1],e[0]+0.1*e[1]];
    Local::from_dissipative_ode(&source,point(&x,&u,&[1.;2]),budget(2),&hessian,
        NonlinearDissipation { flow:&flow,tangent:&tangent,max_tangent_actions:2 },||false).unwrap()
}

#[test]
fn g1_cubic_resistor_uses_differential_not_secant_damping() {
    let local=cubic(0.); let model=local.model();
    assert!((local.equilibrium_nonlinear_dissipated_power_w()-0.02112).abs()<1e-14);
    let tangent=model.nonlinear_resisting_tangent().unwrap();
    assert!((tangent[0]-0.392).abs()<1e-14);
    assert_eq!(&tangent[1..],&[0.;3]);
    let t=0.392; let s=C64::new(-0.221,-(1.0_f64-(0.342_f64/2.).powi(2)).sqrt());
    let pole=model.resonance_pole(s,&[C64::ONE,(s+C64::from_re(t)).scale(2.)],&options(),||false).unwrap();
    near(pole.laplace_pole_per_s,s);
    assert!(pole.decay_power_residual_w.abs()<1e-11);
    let f=[C64::new(0.7,0.2),C64::new(-0.4,0.3)];
    for w in [0.2,0.7,1.,2.] {
        let response=model.factor(w,||false).unwrap().response(&f,||false).unwrap();
        let z=C64::new(0.,-w);
        let determinant=(z+C64::from_re(t))*(z+C64::from_re(0.05))+C64::ONE;
        near(response.state[0],((z+C64::from_re(0.05))*f[0]+f[1].scale(0.5))/determinant);
        near(response.state[1],((z+C64::from_re(t))*f[1]-f[0].scale(2.))/determinant);
        assert!(response.power_residual_w.abs()<1e-11);
    }
}

#[test]
fn g1_nonnegative_original_heat_does_not_hide_negative_incremental_damping() {
    let source=source(); let x=[1.,0.6]; let u=[0.1,2.03];
    let flow=|_: &[f64],e: &[f64],o: &mut [f64]| {o[0]=e[0]/(1.+e[0]*e[0]);o[1]=0.;true};
    let tangent=|_: &[f64],e: &[f64],_: &[f64],de: &[f64],o: &mut [f64]| {
        o[0]=(1.-e[0]*e[0])/(1.+e[0]*e[0]).powi(2)*de[0];o[1]=0.;true
    };
    let local=Local::from_dissipative_ode(&source,point(&x,&u,&[3.;2]),budget(2),&hessian,
        NonlinearDissipation {flow:&flow,tangent:&tangent,max_tangent_actions:2},||false).unwrap();
    assert!((local.equilibrium_nonlinear_dissipated_power_w()-0.8).abs()<1e-14);
    let s=C64::new(0.095,-(1.0_f64-0.145_f64.powi(2)).sqrt());
    let pole=local.model().resonance_pole(s,&[C64::ONE,(s-C64::from_re(0.24)).scale(2.)],&options(),||false).unwrap();
    near(pole.laplace_pole_per_s,s);
    assert!(pole.dissipated_power_w<0. && pole.laplace_pole_per_s.re>0.);
    assert!(pole.decay_power_residual_w.abs()<1e-11);
    assert!(matches!(Local::from_ode(&source,point(&x,&u,&[3.;2]),budget(2),&hessian,||false),
        Err(LinearizationError::NotEquilibrium{..})));
    let active=|_: &[f64],e: &[f64],o: &mut [f64]| {o[0]=-e[0];o[1]=0.;true};
    assert!(matches!(Local::from_dissipative_ode(&source,point(&x,&u,&[3.;2]),budget(2),&hessian,
        NonlinearDissipation {flow:&active,tangent:&tangent,max_tangent_actions:2},||false),
        Err(LinearizationError::ActiveDissipation{..})));
}

#[test]
fn g3_state_and_effort_cross_tangent_is_retained_without_symmetrizing() {
    let local=cubic(0.7); let t=local.model().nonlinear_resisting_tangent().unwrap();
    assert!((t[0]-0.490784).abs()<1e-14 && (t[1]-0.044352).abs()<1e-14);
    assert_eq!(t[2],0.);
    let s=C64::new(-0.270392,-0.9288290296583112);
    let pole=local.model().resonance_pole(s,&[C64::ONE,
        (s+C64::from_re(t[0])).scale(1./(0.5-t[1]))],&options(),||false).unwrap();
    near(pole.laplace_pole_per_s,s);
    let out=local.model().factor(1.,||false).unwrap().response(&[C64::ONE,C64::ZERO],||false).unwrap();
    let nonlinear=0.5*(out.effort[0].conj()*(out.state[0].scale(t[0])+out.state[1].scale(t[1]))).re;
    let linear=0.05*out.effort[1].norm_sq();
    assert!((out.dissipated_power_w-linear-nonlinear).abs()<1e-12);
    assert!(out.power_residual_w.abs()<1e-11);
}

#[test]
fn g1_algebraic_nonlinear_resistor_keeps_kirchhoff_reaction() {
    let j=vec![0.,1.,0., -1.,0.,-1., 0.,1.,0.];
    let mut r=vec![0.;9];r[8]=1.;
    let source=DescriptorPortHamiltonian::new(3,2,2,j,r,vec![1.,0.,0.,1.,0.,0.],
        Box::new(QuadraticStorage::new(vec![4.,0.,0.,1.],2).unwrap())).unwrap();
    let x=[0.2,0.625,0.5];let u=[-0.625,1.3];
    let h=|_: &[f64],d: &[f64],o: &mut [f64]| {o[0]=4.*d[0];o[1]=d[1];o[2]=0.;true};
    let flow=|_: &[f64],e: &[f64],o: &mut [f64]| {o.fill(0.);o[2]=e[2].powi(3);true};
    let tangent=|_: &[f64],e: &[f64],_: &[f64],de: &[f64],o: &mut [f64]| {
        o.fill(0.);o[2]=3.*e[2]*e[2]*de[2];true
    };
    let local=Local::from_dissipative_descriptor(&source,point(&x,&u,&[2.;3]),budget(3),&h,
        NonlinearDissipation {flow:&flow,tangent:&tangent,max_tangent_actions:3},||false).unwrap();
    let decay=0.5/1.75;let s=C64::new(-decay,-(4.0_f64-decay*decay).sqrt());
    let p=local.model().resonance_pole(s,&[C64::ONE,s,s.scale(1./1.75)],&options(),||false).unwrap();
    near(p.laplace_pole_per_s,s);near(p.state[2],p.state[1].scale(1./1.75));
    assert_eq!(p.state.len(),3);assert!(p.decay_power_residual_w.abs()<1e-11);
}

#[test]
fn g1_retained_loss_parameter_derivative_uses_total_constitutive_chain_rule() {
    let local=cubic(0.);let model=local.model();let s=C64::new(-0.221,-0.9852710287022551);
    let p=model.resonance_pole(s,&[C64::ONE,(s+C64::from_re(0.392)).scale(2.)],&options(),||false).unwrap();
    let prepared=model.prepare_pole_response(&p,HarmonicPoleResponseOptions {
        numerical:PoleResponseOptions {spectral_scale:2.,pole_residual_tolerance:1e-12,
            pole_step_tolerance:1e-12,linear_residual_tolerance:1e-11,max_dimension:16,
            max_vectors:4,max_dense_bytes:2_000_000,max_work:10_000_000},
        max_physical_bytes:2_000_000,max_physical_work:10_000_000,
    },||false).unwrap();
    let dt=[1.,0.,0.,0.];
    let d=prepared.parameter_derivatives(&[HarmonicParameterDirection {
        resisting_flow:Some(&dt),..Default::default()
    }],||false).unwrap().remove(0);
    near(d.laplace_pole_per_s,-(s+C64::from_re(0.05))/(s.scale(2.)+C64::from_re(0.442)));
    assert!(d.decay_power_residual_w.abs()<1e-10);
    assert!(d.relative_residual<1e-11);
}

#[test]
fn g0_g4_force_tangent_budget_and_cancellation_fail_without_mutating_the_source() {
    let source=source();let x=[0.;2];let u=[0.;2];let scales=[1.;2];
    let flow=|_: &[f64],_: &[f64],o: &mut [f64]| {o.fill(0.);true};
    let tangent=|_: &[f64],_: &[f64],_: &[f64],_: &[f64],o: &mut [f64]| {o.fill(0.);true};
    let law=NonlinearDissipation {flow:&flow,tangent:&tangent,max_tangent_actions:2};
    let polls=Cell::new(0);
    let make=|law,cancel:&mut dyn FnMut()->bool| Local::from_dissipative_ode(&source,
        point(&x,&u,&scales),budget(2),&hessian,law,cancel);
    let local=make(law,&mut ||{polls.set(polls.get()+1);false}).unwrap();let total=polls.get();
    for cut in [1,total/2,total] {
        let hit=Cell::new(0);
        assert!(make(law,&mut ||{hit.set(hit.get()+1);hit.get()>=cut}).is_err());
    }
    assert!(make(NonlinearDissipation {max_tangent_actions:1,..law},&mut ||false).is_err());
    let incomplete=|_: &[f64],_: &[f64],_: &[f64],_: &[f64],o:&mut[f64]| {o[0]=0.;true};
    assert!(matches!(make(NonlinearDissipation {tangent:&incomplete,..law},&mut ||false),
        Err(LinearizationError::DissipationRefused{column:Some(0)})));
    let ordinary=HarmonicModel::from_ode(&source,budget(2).harmonic,||false).unwrap();
    let f=[C64::ONE,C64::ZERO];
    let a=local.model().factor(1.,||false).unwrap().response(&f,||false).unwrap();
    let b=ordinary.factor(1.,||false).unwrap().response(&f,||false).unwrap();
    assert_eq!(a.state,b.state);assert_eq!(a.output,b.output);assert_eq!(a.dissipated_power_w,b.dissipated_power_w);
}
