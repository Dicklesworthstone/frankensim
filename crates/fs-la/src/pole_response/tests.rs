use super::*;
use crate::eigen_complex::lu_complex;
use std::cell::Cell;

fn options() -> PoleResponseOptions {
    PoleResponseOptions { spectral_scale:3.0,pole_residual_tolerance:1e-12,
        pole_step_tolerance:1e-12,linear_residual_tolerance:1e-12,
        max_dimension:8,max_vectors:8,max_dense_bytes:1_000_000,max_work:1_000_000 }
}
fn c(x:f64) -> C64 { C64::from_re(x) }
fn close(a:C64,b:C64) { assert!((a-b).abs()<2e-11*(1.0+b.abs()),"{a:?} != {b:?}"); }
fn nonsymmetric(s:C64,scale:f64) -> PoleEvaluation {
    let a=C64::new(-0.1,2.0);
    let b=C64::new(-3.0,1.0);
    let beta=C64::new(1.2,-0.7);
    PoleEvaluation { matrix:vec![s-a,beta,C64::ZERO,s-b].into_iter().map(|z|z.scale(scale)).collect(),
        derivative:vec![c(scale),C64::ZERO,C64::ZERO,c(scale)],residual_scale:10.0*scale }
}
fn prepared(scale:f64) -> PreparedPole {
    let a=C64::new(-0.1,2.0);
    PreparedPole::new(a,&[C64::new(0.2,-0.4),C64::ZERO],&nonsymmetric(a,scale),options(),||false).unwrap()
}

#[test]
fn g2_nonnormal_complex_force_residue_matches_exact_resolvent_and_limit() {
    let prepared=prepared(1.0);
    let s=prepared.value();
    let b=C64::new(-3.0,1.0);
    let beta=C64::new(1.2,-0.7);
    let f=[C64::new(0.7,-0.2),C64::new(-0.3,0.8)];
    let result=prepared.forcing_residues(&[&f],||false).unwrap();
    close(result[0].vector[0],f[0]-beta*f[1]/(s-b));
    close(result[0].vector[1],C64::ZERO);
    assert!(result[0].linear_residual<1e-13);
    // Original linear solves on a ladder: (s-s0)T(s)^-1 f -> R f.
    let mut errors=Vec::new();
    for h in [1e-2,5e-3,2.5e-3] {
        let shift=C64::new(h,0.4*h);
        let a=nonsymmetric(s+shift,1.0);
        let mut x=f.to_vec();lu_complex(&a.matrix,2).unwrap().solve(&mut x);
        errors.push((x[0]*shift-result[0].vector[0]).abs());
    }
    for pair in errors.windows(2) { assert!(pair[0]/pair[1]>1.98 && pair[0]/pair[1]<2.02); }
}

#[test]
fn g1_parameter_tangent_satisfies_full_nonnormal_eigenvector_equations() {
    let prepared=prepared(1.0);
    let action=[C64::ZERO,C64::ONE]; // lower-left coefficient theta
    let d=prepared.parameter_derivatives(&[&action],||false).unwrap().remove(0);
    let gap=prepared.value()-C64::new(-3.0,1.0);
    close(d.value,C64::new(1.2,-0.7)/gap);
    close(d.vector[0],C64::ZERO);
    close(d.vector[1],-gap.recip());
    let a=nonsymmetric(prepared.value(),1.0);
    for i in 0..2 {
        let sum=(0..2).fold(action[i],|sum,j|sum+a.matrix[2*i+j]*d.vector[j]
            +a.derivative[2*i+j]*prepared.right_vector()[j]*d.value);
        assert!(sum.abs()<1e-12);
    }
    // Analytic quadratic roots of (s-a)(s-b)-beta theta, followed on the same branch.
    let a=prepared.value();let b=C64::new(-3.0,1.0);let beta=C64::new(1.2,-0.7);
    let root=|theta:f64| {
        let disc=((a-b)*(a-b)+beta.scale(4.0*theta)).sqrt();
        let plus=(a+b+disc).scale(0.5);let minus=(a+b-disc).scale(0.5);
        if (plus-a).abs()<(minus-a).abs() {plus} else {minus}
    };
    let errors:Vec<_>=[1e-3,1e-4,1e-5].into_iter().map(|h|
        ((root(h)-root(-h)).scale(0.5/h)-d.value).abs()).collect();
    assert!(errors[0]/errors[1]>80.0 && errors[0]/errors[1]<120.0);
    assert!(errors[2]<1e-10);
}

#[test]
fn g2_scalar_damped_resonance_and_zero_descriptor_pole_keep_exact_derivatives() {
    // T(s)=s^2+0.4s+4, derivative with respect to stiffness and damping.
    let s=C64::new(-0.2,3.96_f64.sqrt());
    let a=PoleEvaluation {matrix:vec![s*s+s.scale(0.4)+c(4.0)],
        derivative:vec![s.scale(2.0)+c(0.4)],residual_scale:9.0};
    let p=PreparedPole::new(s,&[C64::ONE],&a,options(),||false).unwrap();
    let dk=[C64::ONE];let dc=[s];
    let d=p.parameter_derivatives(&[&dk,&dc],||false).unwrap();
    close(d[0].value,-a.derivative[0].recip());
    close(d[1].value,-s/a.derivative[0]);
    let f=[C64::new(0.5,0.8)];
    close(p.forcing_residues(&[&f],||false).unwrap()[0].vector[0],f[0]/a.derivative[0]);
    // Singular descriptor mass diag(1,0), finite root zero. No inverse mass.
    let a=PoleEvaluation {matrix:vec![C64::ZERO,c(2.0),C64::ZERO,C64::ONE],
        derivative:vec![C64::ONE,C64::ZERO,C64::ZERO,C64::ZERO],residual_scale:3.0};
    let p=PreparedPole::new(C64::ZERO,&[C64::ONE,C64::ZERO],&a,options(),||false).unwrap();
    close(p.forcing_residues(&[&[C64::ZERO,C64::ONE]],||false).unwrap()[0].vector[0],c(-2.0));
}

#[test]
fn g3_operator_scaling_gauge_and_batch_linearity_preserve_physical_residues() {
    let f=[C64::new(0.3,0.7),C64::new(-0.2,0.4)];
    let baseline=prepared(1.0).forcing_residues(&[&f],||false).unwrap();
    for scale in [1e-120,1e120] {
        let p=prepared(scale);let scaled=f.map(|z|z.scale(scale));
        let r=p.forcing_residues(&[&scaled],||false).unwrap();
        for i in 0..2 {close(r[0].vector[i],baseline[0].vector[i]);}
        let d=p.parameter_derivatives(&[&[C64::ZERO,c(scale)]],||false).unwrap();
        close(d[0].value,C64::new(1.2,-0.7)/(p.value()-C64::new(-3.0,1.0)));
    }
    let p=prepared(1.0);let g=[C64::new(0.1,-0.5),C64::new(0.8,0.4)];
    let h=[f[0]+g[0],f[1]+g[1]];let zero=[C64::ZERO;2];
    let r=p.forcing_residues(&[&f,&g,&h,&zero],||false).unwrap();
    for i in 0..2 {close(r[2].vector[i],r[0].vector[i]+r[1].vector[i]);}
    assert!(r[3].vector.iter().all(|z|z.abs()==0.0));
    assert_eq!(r[0].vector,p.forcing_residues(&[&f],||false).unwrap()[0].vector);
}

#[test]
fn g0_unresolved_roots_repeated_poles_and_resource_boundaries_refuse() {
    let s=C64::new(-0.1,2.0);
    assert!(matches!(PreparedPole::new(s+c(0.1),&[C64::ONE,C64::ZERO],
        &nonsymmetric(s+c(0.1),1.0),options(),||false),Err(PoleResponseError::NotPole{..})));
    for matrix in [vec![C64::ZERO;4],vec![C64::ZERO,C64::ONE,C64::ZERO,C64::ZERO]] {
        let a=PoleEvaluation {matrix,derivative:vec![C64::ONE,C64::ZERO,C64::ZERO,C64::ONE],residual_scale:1.0};
        assert!(matches!(PreparedPole::new(C64::ZERO,&[C64::ONE,C64::ZERO],&a,options(),||false),
            Err(PoleResponseError::Singular{..})));
    }
    for which in 0..5 {
        let mut o=options();
        match which {0=>o.max_dimension=1,1=>o.max_vectors=0,2=>o.max_dense_bytes=0,
            3=>o.max_work=0,_=>o.pole_step_tolerance=f64::NAN};
        assert!(PreparedPole::new(s,&[C64::ONE,C64::ZERO],&nonsymmetric(s,1.0),o,||false).is_err());
    }
    let p=prepared(1.0);
    assert!(p.forcing_residues(&[],||false).is_err());
    assert!(p.parameter_derivatives(&[&[C64::ONE]],||false).is_err());
    assert!(p.forcing_residues(&[&[C64::new(f64::NAN,0.0),C64::ZERO]],||false).is_err());
}

#[test]
fn g4_g5_late_batch_cancellation_never_mutates_prepared_factor() {
    let p=prepared(1.0);let f=[C64::ONE,C64::ONE];
    let polls=Cell::new(0);
    let baseline=p.forcing_residues(&[&f,&f],||{polls.set(polls.get()+1);false}).unwrap();
    for stop in [1,polls.get()/2,polls.get()] {
        let count=Cell::new(0);
        assert!(matches!(p.forcing_residues(&[&f,&f],||{count.set(count.get()+1);count.get()==stop}),
            Err(PoleResponseError::Cancelled)));
    }
    let replay=p.forcing_residues(&[&f,&f],||false).unwrap();
    assert_eq!(baseline[0].vector,replay[0].vector);
    assert_eq!(baseline[1].coefficient,replay[1].coefficient);
}
