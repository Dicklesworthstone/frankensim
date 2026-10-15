//! Exercise the actual implicit effort and endpoint publication boundary.
use super::*;
#[path="../../../../../fs-contact/tests/source_point/common.rs"]
mod common;

#[test]
fn discrete_gradient_rate_and_nonfinite_closure_cannot_bypass_source_limits() {
    let c=MaterialContact::new("test".into(),common::source(1e7,Some(0.4)),vec![1.0,-1.0],0.0).unwrap();
    let (a,_)=ImpactBody::free_mass(1.0,0.0,0.0).unwrap();
    let mut s=ImpactSystem::new_with_material_contacts(vec![a.clone(),a],vec![c],vec![],vec![],
        ImpactConfig {dt_s:1e-5,max_steps:10,maximum_energy_j:1.0,
            energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-7,maximum_generalized_force:1.0}).unwrap();
    let x=[1e-5,0.0,0.0,0.0];let e=[0.0,0.3,0.0,0.0];let mut out=[17.0;4];
    assert!(s.validate_material_contacts(&x,&x).is_ok());
    assert!(!s.dissipative_flow_into(&x,&e,&mut out));assert_eq!(out,[17.0;4]);
    assert!(!s.dissipative_flow_tangent_into(&x,&e,&x,&x,&mut out));
    let mut nan=x;nan[0]=f64::NAN;assert!(s.validate_material_contacts(&nan,&x).is_err());
    let before=s.state().to_vec();let mut endpoint=vec![2e-4,0.0,0.0,0.0];
    let mut histories=vec![];let energy=s.stored_energy_j();
    let result=s.accept_step(&mut endpoint,&mut histories,energy,
        fs_phs::PreparedStepRecord {delta_h:0.0,dissipated:0.0,supplied:0.0,newton_iters:0,solver_residual:0.0},
        &fs_exec::CancelGate::new_clock_free());
    assert!(matches!(result,Err(ImpactError::MaterialContact {..})));
    assert_eq!(s.state(),before);assert_eq!(s.samples(),0);
}
