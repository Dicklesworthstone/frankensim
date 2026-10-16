//! One source substitution updates the entire homogeneous collision model.
use super::*;
use fs_couple::render::plate::impact::material::pair::{TranslatingMaterialPair,TranslatingSpecimen};

fn compiled(rho:f64,young:f64,chi:Option<f64>) -> TranslatingMaterialPair {
    let a=solid("A",rho,young);let b=solid("B",2000.0,2e7);
    // The operating patch stays small relative to a 1 cm specimen depth.
    let mut operating=envelope();operating.maximum_indentation_m=5e-5;
    operating.half_space_depth_m=0.01;operating.layer_thickness_m=0.01;
    let source=SourcePointContact::new(bind(&a,&b,chi,0.0),
        PointContactGeometry::SpherePlane{effective_radius_m:0.01},operating).unwrap();
    TranslatingMaterialPair::new("uniform-pair".into(),source,[
        TranslatingSpecimen{volume_m3:1e-6,initial_position_m:-2e-5,initial_velocity_m_per_s:0.04},
        TranslatingSpecimen{volume_m3:1e-6,initial_position_m:0.0,initial_velocity_m_per_s:0.0},
    ],0.0).unwrap()
}

#[test]
fn one_material_pair_supplies_masses_contact_weights_and_original_sources() {
    let pair=compiled(1000.0,1e7,None);let masses=pair.masses_kg();
    close(masses[0],0.001);close(masses[1],0.002);
    let system=pair.build(config()).unwrap();
    let c=&system.material_contacts()[0];
    assert_eq!(c.source(),pair.source());close(c.weights()[0],1.0/masses[0].sqrt());
    close(c.weights()[1],-1.0/masses[1].sqrt());
    close(system.state()[1]/masses[0].sqrt(),0.04);
    close(system.stored_energy_j(),0.5*masses[0]*0.04_f64.powi(2));
    let second=pair.build(config()).unwrap();assert_eq!(system.state(),second.state());
    assert_eq!(system.material_contacts(),second.material_contacts());
}

#[test]
fn density_and_elasticity_swaps_propagate_through_the_actual_same_pair_builder() {
    let base=compiled(1000.0,1e7,Some(0.4));let heavy=compiled(2000.0,1e7,Some(0.4));
    let stiff=compiled(1000.0,3e7,Some(0.4));
    close(heavy.masses_kg()[0]/base.masses_kg()[0],2.0);
    close(heavy.source().stiffness_n_per_m_3_2(),base.source().stiffness_n_per_m_3_2());
    assert_eq!(stiff.masses_kg(),base.masses_kg());
    assert_eq!(base.specimens(),heavy.specimens());assert_eq!(base.specimens(),stiff.specimens());
    let mut a=base.build(config()).unwrap().prepare_analytic().unwrap();
    let mut b=heavy.build(config()).unwrap().prepare_analytic().unwrap();
    let mut c=stiff.build(config()).unwrap().prepare_analytic().unwrap();
    let mut peak=[0.0_f64;2];let gate=CancelGate::new_clock_free();
    for _ in 0..2000 {
        a.step(&[0.0;2],&gate).unwrap();b.step(&[0.0;2],&gate).unwrap();c.step(&[0.0;2],&gate).unwrap();
        peak[0]=peak[0].max(a.material_contact_observation(0).unwrap().normal_force_n);
        peak[1]=peak[1].max(c.material_contact_observation(0).unwrap().normal_force_n);
    }
    assert!(b.state()[3]>a.state()[3]+0.0001);
    assert!(peak[1]>1.2*peak[0]);
    assert_eq!(a.material_contacts()[0].source().model().identity(),base.source().model().identity());
}

#[test]
fn nonphysical_volume_motion_and_invalid_time_cannot_create_a_pair_runtime() {
    let pair=compiled(1000.0,1e7,None);
    for volume in [0.0,-1.0,f64::NAN,f64::INFINITY,f64::MAX] {
        let mut specimens=*pair.specimens();specimens[0].volume_m3=volume;
        assert!(TranslatingMaterialPair::new("invalid".into(),pair.source().clone(),specimens,0.0).is_err());
    }
    for (position,velocity) in [(f64::NAN,0.0),(0.0,f64::INFINITY)] {
        let mut specimens=*pair.specimens();specimens[1].initial_position_m=position;
        specimens[1].initial_velocity_m_per_s=velocity;
        assert!(TranslatingMaterialPair::new("invalid".into(),pair.source().clone(),specimens,0.0).is_err());
    }
    let mut slow=config();slow.dt_s=0.1;
    assert!(pair.build(slow).is_err());
    assert!(TranslatingMaterialPair::new(" ".into(),pair.source().clone(),*pair.specimens(),0.0).is_err());
    assert!(TranslatingMaterialPair::new("invalid".into(),pair.source().clone(),*pair.specimens(),f64::NAN).is_err());
}
