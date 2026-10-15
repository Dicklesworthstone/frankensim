//! Source-bound Hertz and Hunt--Crossley through the existing patch owner.
#[path="source_point/common.rs"]
mod common;
use common::*;
use fs_contact::normal_patch::source::*;
use fs_contact::normal_patch::NormalPatchError;
use fs_tribo::InputAuthority;

#[test]
fn source_data_changes_coefficient_and_original_pressure_together() {
    let soft=source(1e7,None); let stiff=source(3e7,None);
    let ratio=stiff.model().normal_state().elastic().reduced_modulus_pa()
        /soft.model().normal_state().elastic().reduced_modulus_pa();
    close(stiff.stiffness_n_per_m_3_2()/soft.stiffness_n_per_m_3_2(),ratio);
    let a=soft.evaluate("sample",2e-5,0.05,1e-4).unwrap();
    let b=stiff.evaluate("sample",2e-5,0.05,1e-4).unwrap();
    close(b.normal_force_n/a.normal_force_n,ratio);
    close(b.pressure.peak_pressure_pa/a.pressure.peak_pressure_pa,ratio);
    close(a.patch_radius_m,b.patch_radius_m);
    assert_ne!(soft.model().identity(),stiff.model().identity());
    assert_eq!(a.authority,InputAuthority::SyntheticFixture);
    assert_eq!(a.uncertainty,envelope().uncertainty);
}

#[test]
fn lowered_potential_force_tangent_and_loss_match_the_original_law() {
    for chi in [None,Some(0.4)] {
        let law=source(1e7,chi); let k=law.stiffness_n_per_m_3_2(); let c=law.dissipation_s_per_m();
        for d in [0.0_f64,1e-8,2e-5,1e-4] { for rate in [-0.2,0.0,0.2] {
            let r=law.evaluate("sample",d,rate,1e-4).unwrap();
            close(r.normal_force_n,k*d*d.sqrt()*(1.0+c*rate));
            close(r.tangent_n_per_m,1.5*k*d.sqrt()*(1.0+c*rate));
            close(r.reversible_energy_j,0.4*k*d*d*d.sqrt());
            close(r.dissipated_power_w,c*k*d*d.sqrt()*rate*rate);
        }}
    }
}

#[test]
fn actual_elliptic_geometry_is_not_replaced_by_a_scalar_radius() {
    let a=solid("A",1000.0,1e7); let b=solid("B",2000.0,2e7);
    let model=bind(&a,&b,Some(0.4),0.0);
    let law=SourcePointContact::new(model.clone(),PointContactGeometry::EllipticParaboloid {
        maximum_curvature_m_inverse:20.0,minimum_curvature_m_inverse:5.0 },envelope()).unwrap();
    for d in [1e-6_f64,1e-5,1e-4] {
        let r=law.evaluate("ellipse",d,0.0,1e-4).unwrap();
        let axes=r.elliptic_patch_axes.unwrap(); assert!(axes.aspect_ratio>1.0);
        close(r.normal_force_n,law.stiffness_n_per_m_3_2()*d*d.sqrt());
    }
    let circular=SourcePointContact::new(model,PointContactGeometry::EllipticParaboloid {
        maximum_curvature_m_inverse:10.0,minimum_curvature_m_inverse:10.0 },envelope()).unwrap();
    close(circular.stiffness_n_per_m_3_2(),source(1e7,Some(0.4)).stiffness_n_per_m_3_2());
}

#[test]
fn applicability_is_checked_before_lowering_including_separation_force() {
    let a=solid("A",1000.0,1e7);let b=solid("B",2000.0,2e7);
    let geom=PointContactGeometry::SpherePlane{effective_radius_m:0.1};
    assert!(matches!(SourcePointContact::new(bind(&a,&b,None,0.1),geom,envelope()),
        Err(NormalPatchError::AdhesionUnsupported{..})));
    assert!(SourcePointContact::new(bind(&a,&b,Some(6.0),0.0),geom,envelope()).is_err());
    let mut bad=envelope();bad.maximum_indentation_m=0.01;
    assert!(SourcePointContact::new(bind(&a,&b,None,0.0),geom,bad).is_err());
    bad=envelope();bad.layer_thickness_m=1e-6;
    assert!(SourcePointContact::new(bind(&a,&b,None,0.0),geom,bad).is_err());
    bad=envelope();bad.maximum_indentation_m=f64::NAN;
    assert!(SourcePointContact::new(bind(&a,&b,None,0.0),geom,bad).is_err());
}

#[test]
fn runtime_envelope_and_sample_identity_remain_explicit() {
    let law=source(1e7,Some(0.4));
    for (d,v,dt) in [(1.01e-4,0.0,1e-4),(0.0,0.201,1e-4),(0.0,-0.201,1e-4),
        (0.0,0.0,0.002),(-1e-6,0.0,1e-4),(0.0,f64::NAN,1e-4),(0.0,0.0,0.0)] {
        assert!(law.validate_kinematics(d,v,dt).is_err());
    }
    let a=law.evaluate("first",1e-5,0.02,1e-4).unwrap();
    let b=law.evaluate("second",1e-5,0.02,1e-4).unwrap();
    assert_ne!(a.request_id,b.request_id);close(a.normal_force_n,b.normal_force_n);
    assert!(law.evaluate("",1e-5,0.02,1e-4).is_err());
}

#[path="source_point/temperature.rs"]
mod temperature;
