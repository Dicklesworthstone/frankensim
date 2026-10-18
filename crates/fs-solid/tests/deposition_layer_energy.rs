//! Layer heat is measured at the original material updates, not split by mass.
use std::collections::BTreeMap;
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_solid::deposition::*;

fn law(branch: f64, upper: f64) -> ThermoMaxwell {
    ThermoMaxwell::new(1e8,vec![(branch,1.0)],vec![
        ThermalKnot {temperature_k:280.0,modulus_scale:1.0,relaxation_scale:1.0,free_strain:0.0},
        ThermalKnot {temperature_k:upper,modulus_scale:1.0,relaxation_scale:1.0,free_strain:0.0},
    ],Limits {max_strain:0.04,max_strain_rate_s:10.0}).unwrap()
}
fn initial<'a>(a: &'a ThermoMaxwell, b: &'a ThermoMaxwell) -> LayeredStrip<'a> {
    let mut strip = LayeredStrip::new(0.1,Budget {max_layers:2,max_branch_states:4}).unwrap();
    for (id,z,law) in [(7,-0.0005,a),(11,0.0005,b)] {
        strip = strip.deposit(id,LayerGeometry {area_m2:1e-5,centroid_z_m:z,thickness_m:0.001},law,300.0).unwrap();
    }
    strip
}
fn config() -> Config {
    Config {force_tolerance_n:1e-8,moment_tolerance_nm:1e-10,energy_tolerance_j:1e-10,
        max_strain:0.04,max_rotation_rad:0.1}
}
fn control() -> Control {
    Control {axial:AxialControl::Strain(0.01),bending:BendingControl::CurvaturePerM(0.0)}
}
#[test]
fn local_heat_is_from_the_actual_points_and_does_not_change_the_existing_solve() {
    let a=law(3e8,320.0); let b=law(1e8,320.0); let old=initial(&a,&b);
    let t=BTreeMap::from([(7,301.0),(11,303.0)]);
    let (plain,original)=old.advance_to(0.2,&t,control(),config(),||Ok(())).unwrap();
    let (next,total,layers)=old.advance_with_layer_energy(0.2,&t,control(),config(),||Ok(())).unwrap();
    assert_eq!(plain,next);
    assert_eq!(original.axial_force_n.to_bits(),total.axial_force_n.to_bits());
    assert_eq!(original.energy_j.dissipated.to_bits(),total.energy_j.dissipated.to_bits());
    assert_eq!(total.material_point_updates,4);
    let q=0.01/1.2;
    for (id,e) in [(7,3e8),(11,1e8)] {
        let expected=1e-6*e*0.2*q*q;
        assert!((layers[&id].dissipated-expected).abs()<1e-14);
        assert!(layers[&id].numerical_dissipation>0.0);
        assert!(layers[&id].residual.abs()<1e-12);
    }
    assert!((layers.values().map(|e|e.dissipated).sum::<f64>()-total.energy_j.dissipated).abs()<1e-14);
    assert!((layers.values().map(|e|e.stored_change).sum::<f64>()-total.energy_j.stored_change).abs()<1e-14);
    assert_eq!(old.axial_strain(),0.0);
}
#[test]
fn late_source_or_cancellation_failure_publishes_no_layer_history() {
    let a=law(3e8,320.0); let b=law(1e8,310.0); let old=initial(&a,&b);
    let saved=old.clone();
    assert!(old.advance_with_layer_energy(0.2,&BTreeMap::from([(7,305.0),(11,315.0)]),
        control(),config(),||Ok(())).is_err());
    let t=BTreeMap::from([(7,305.0),(11,305.0)]);
    let mut polls=0;
    old.advance_with_layer_energy(0.2,&t,control(),config(),||{polls+=1;Ok(())}).unwrap();
    let mut current=0;
    assert!(matches!(old.advance_with_layer_energy(0.2,&t,control(),config(),||{
        current+=1; if current==polls {Err(Error::Cancelled)} else {Ok(())}
    }),Err(Error::Cancelled)));
    assert_eq!(old,saved);
}
