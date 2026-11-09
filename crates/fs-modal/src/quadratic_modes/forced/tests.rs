use super::*;
use crate::quadratic_modes::{PoleOptions,PreparedQuadraticModes};

fn original<'a>(m:&'a[f64],c:&'a[f64],k:&'a[f64],n:usize)->PreparedQuadraticModes<'a> {
    PreparedQuadraticModes::new(m,c,k,n,PoleOptions {
        spectral_scale:1.0,step_tolerance:1e-12,residual_tolerance:1e-12,
        max_iterations:8,max_dimension:8,max_dense_bytes:4_000_000,max_work:10_000_000,
    },||false).unwrap()
}
fn options()->QuadraticForcedOptions {
    QuadraticForcedOptions {maximum_relative_residual:1e-11,
        maximum_payload_bytes:1_000_000,maximum_work:10_000_000}
}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"{a:.17e} != {b:.17e}");}
fn complex_near(a:C64,b:C64,tol:f64) {near(a.re,b.re,tol);near(a.im,b.im,tol);}

#[test]
fn static_compliance_and_damped_resonance_use_physical_force_and_power() {
    let (m,c,k)=([2.0],[3.0],[8.0]);
    let dc=original(&m,&c,&k,1).into_forced_response(0.0,options(),||false).unwrap();
    let x=dc.solve(&[C64::from_re(4.0)],||false).unwrap();
    complex_near(x.displacement[0],C64::from_re(0.5),1e-15);
    assert_eq!(x.velocity[0],C64::ZERO);
    assert_eq!(x.balance.average_input_power_w,0.0);
    let ac=original(&m,&c,&k,1).into_forced_response(2.0,options(),||false).unwrap();
    let x=ac.solve(&[C64::from_re(6.0)],||false).unwrap();
    complex_near(x.displacement[0],C64::new(0.0,1.0),1e-14);
    complex_near(x.velocity[0],C64::from_re(2.0),1e-14);
    near(x.balance.average_input_power_w,6.0,1e-13);
    near(x.balance.average_velocity_power_w,6.0,1e-13);
    near(x.balance.average_stiffness_power_w,0.0,1e-13);
    near(x.balance.average_inertial_power_w,0.0,1e-13);
    near(x.balance.power_balance_residual_w,0.0,1e-13);
}

#[test]
fn signed_unstable_stiffness_and_active_velocity_work_remain_signed() {
    let (m,c,k)=([2.0],[-3.0],[-4.0]);
    let prepared=original(&m,&c,&k,1).into_forced_response(2.0,options(),||false).unwrap();
    let f=C64::new(2.0,-1.0);
    let x=prepared.solve(&[f],||false).unwrap();
    complex_near(x.displacement[0],f/C64::new(-12.0,6.0),2e-15);
    assert!(x.balance.average_velocity_power_w<0.0);
    near(x.balance.average_input_power_w,x.balance.average_velocity_power_w,1e-14);
    assert!(x.balance.relative_force_residual<1e-13);
}

#[test]
fn nonsymmetric_stiffness_and_skew_gyroscopic_velocity_keep_distinct_work() {
    let m=[2.0,0.0,0.0,3.0];
    let c=[1.0,-4.0,4.0,2.0];
    let k=[8.0,3.0,-1.0,9.0];
    let omega=1.4;
    let q=[C64::new(1.0,0.2),C64::new(-0.3,0.5)];
    let s=C64::new(0.0,-omega);
    let force=[
        (s*s).scale(2.0)*q[0]+s*(q[0]-q[1].scale(4.0))+q[0].scale(8.0)+q[1].scale(3.0),
        (s*s).scale(3.0)*q[1]+s*(q[0].scale(4.0)+q[1].scale(2.0))-q[0]+q[1].scale(9.0),
    ];
    let response=original(&m,&c,&k,2).into_forced_response(omega,options(),||false).unwrap()
        .solve(&force,||false).unwrap();
    for i in 0..2 {complex_near(response.displacement[i],q[i],4e-14);}
    // G does zero mean work. The skew stiffness supplies -1.568 W.
    near(response.balance.average_velocity_power_w,1.6856,2e-13);
    near(response.balance.average_stiffness_power_w,-1.568,2e-13);
    near(response.balance.average_input_power_w,0.1176,3e-13);
    near(response.balance.average_inertial_power_w,0.0,2e-13);
    near(response.balance.power_balance_residual_w,0.0,3e-13);
}

#[test]
fn represented_narrow_detuning_is_not_replaced_by_a_false_singularity() {
    let omega=1.0+2.0_f64.powi(-27);
    let k=[1.0+2.0_f64.powi(-26)];
    let response=original(&[1.0],&[0.0],&k,1)
        .into_forced_response(omega,options(),||false).unwrap()
        .solve(&[C64::from_re(-2.0_f64.powi(-54))],||false).unwrap();
    complex_near(response.displacement[0],C64::ONE,2e-14);
    // Exact undamped resonance still has no unique finite particular response.
    assert!(matches!(original(&[2.0],&[0.0],&[8.0],1)
        .into_forced_response(2.0,options(),||false),
        Err(QuadraticForcedError::Linear(LuComplexError::Singular{..}))));
}

#[test]
fn force_amplitude_scaling_preserves_dual_units_and_physical_power() {
    let (m,c,k)=([2.0],[3.0],[8.0]);
    let factor=original(&m,&c,&k,1).into_forced_response(2.0,options(),||false).unwrap();
    for e in [-400,0,400] {
        let amplitude=2.0_f64.powi(e);
        let x=factor.solve(&[C64::from_re(6.0*amplitude)],||false).unwrap();
        complex_near(x.displacement[0].scale(1.0/amplitude),C64::new(0.0,1.0),2e-14);
        near(x.balance.average_input_power_w/(amplitude*amplitude),6.0,2e-13);
        near(x.balance.average_velocity_power_w/(amplitude*amplitude),6.0,2e-13);
    }
}

#[test]
fn independent_physical_balance_and_late_cancellation_refuse_atomically() {
    let zero=[C64::ZERO];
    assert!(evaluate_harmonic_balance(&[C64::ONE],&zero,&zero,&zero,&zero,&zero,
        1e-12,||false).is_err());
    let locked=evaluate_harmonic_balance(&[C64::ONE],&[C64::from_re(-1.0)],
        &zero,&zero,&zero,&zero,1e-12,||false).unwrap();
    assert_eq!(locked.relative_force_residual,0.0);
    assert_eq!(locked.average_reaction_power_w,0.0);
    let (m,c,k)=([2.0],[3.0],[8.0]);
    let factor=original(&m,&c,&k,1).into_forced_response(2.0,options(),||false).unwrap();
    let f=[C64::new(3.0,2.0)];
    let mut checkpoints=0;
    let before=factor.solve(&f,||{checkpoints+=1;false}).unwrap();
    let mut calls=0;
    assert!(matches!(factor.solve(&f,||{calls+=1;calls==checkpoints}),
        Err(QuadraticForcedError::Cancelled)));
    assert_eq!(before,factor.solve(&f,||false).unwrap());
    assert!(factor.solve(&[C64::new(f64::NAN,0.0)],||false).is_err());
    let mut tiny=options();tiny.maximum_payload_bytes=1;
    assert!(original(&m,&c,&k,1).into_forced_response(2.0,tiny,||false).is_err());
    assert_eq!(m,[2.0]);assert_eq!(c,[3.0]);assert_eq!(k,[8.0]);
}
