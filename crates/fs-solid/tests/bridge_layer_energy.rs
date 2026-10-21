//! Actual spatial material heat without changing the existing beam solution.
use std::collections::BTreeMap;
use fs_material::thermovisco::{Limits, ThermalKnot, ThermoMaxwell};
use fs_solid::{bridge::{self, Bridge, Support}, deposition};

fn law() -> ThermoMaxwell {
    ThermoMaxwell::new(1e8, vec![(3e8, 1.0)], [280.0, 320.0].map(|temperature_k|
        ThermalKnot {temperature_k, modulus_scale:1.0, relaxation_scale:1.0, free_strain:0.0}).to_vec(),
        Limits {max_strain:0.04, max_strain_rate_s:10.0}).unwrap()
}
fn config() -> bridge::Config {
    bridge::Config {section:deposition::Config {force_tolerance_n:1e-8,
        moment_tolerance_nm:1e-10, energy_tolerance_j:1e-10, max_strain:0.04,
        max_rotation_rad:0.1}, maximum_deflection_m:0.009}
}
fn initial(law: &ThermoMaxwell, support: Support) -> Bridge<'_> {
    let mut bridge = Bridge::new(&[(20,0.05),(3,0.05)], support,
        bridge::Budget {max_cells:2,max_layers_per_cell:2,max_branch_states:8}).unwrap();
    for (id,z) in [(7,-0.0005),(11,0.0005)] {
        bridge = bridge.deposit_bonded_layer(id, deposition::LayerGeometry {
            area_m2:1e-5,centroid_z_m:z,thickness_m:0.001}, 1000.0, law,
            &BTreeMap::from([(20,300.0),(3,300.0)]), &mut ||false).unwrap();
    }
    bridge
}
fn temperatures() -> BTreeMap<(u64,u64),f64> {
    BTreeMap::from([((20,7),300.0),((20,11),300.0),((3,7),300.0),((3,11),300.0)])
}
#[test]
fn local_heat_preserves_both_support_solutions_and_actual_spatial_partition() {
    let law=law();
    for support in [Support::CantileverLeft,Support::SimplySupported] {
        let old=initial(&law,support);
        let (plain,a)=old.advance_controlled(0.2,&temperatures(),9.81,&[1.0;2],config(),&mut ||false).unwrap();
        let (next,b,local)=old.advance_with_cell_energy(0.2,&temperatures(),9.81,&[1.0;2],config(),&mut ||false).unwrap();
        assert_eq!(next,plain);
        assert_eq!(a.energy_j.dissipated.to_bits(),b.energy_j.dissipated.to_bits());
        assert_eq!(a.load_work_j.to_bits(),b.load_work_j.to_bits());
        assert_eq!(a.vertical_reactions_n,b.vertical_reactions_n);
        assert_eq!(b.material_point_updates,8);
        assert_eq!(local.len(),4);
        let sum=local.values().map(|e|e.dissipated).sum::<f64>();
        assert!((sum-b.energy_j.dissipated).abs()<1e-14);
        assert!(local.values().all(|e|e.dissipated>0.0 && e.numerical_dissipation>0.0));
        if support==Support::CantileverLeft {
            // Equal cells/materials: exact average moment ratio is 7, so heat ratio is 49.
            assert!((local[&(20,7)].dissipated/local[&(3,7)].dissipated-49.0).abs()<1e-9);
        }
        assert_eq!(old.time_s(),0.0);
    }
}
#[test]
fn late_refusal_and_cancellation_return_no_partial_cell_energy_or_history() {
    let law=law(); let old=initial(&law,Support::CantileverLeft); let saved=old.clone();
    let mut invalid=temperatures(); invalid.insert((3,11),321.0);
    assert!(old.advance_with_cell_energy(0.2,&invalid,9.81,&[1.0;2],config(),&mut ||false).is_err());
    let mut polls=0;
    old.advance_with_cell_energy(0.2,&temperatures(),9.81,&[1.0;2],config(),&mut ||{polls+=1;false}).unwrap();
    let mut current=0;
    assert!(matches!(old.advance_with_cell_energy(0.2,&temperatures(),9.81,&[1.0;2],config(),
        &mut ||{current+=1;current==polls}),Err(bridge::Error::Cancelled)));
    assert_eq!(old,saved);
}
