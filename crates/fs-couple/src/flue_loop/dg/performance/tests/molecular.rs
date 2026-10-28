use super::*;

mod field;

// Deliberately synthetic constitutive rates: this test verifies the shared
// material/duct/render connection, not N2/O2 source-data accuracy.
const LAW: &str =
    "dg_molecular declared 2 N2 1 vibration 1 1000 1 power 10 293.15 0 273.15 323.15 O2 0\n";

#[test]
fn native_flue_molecular_data_changes_solved_pressure_and_preserves_block_identity() {
    let original = input();
    let changed = original.replace("\nflue ", &format!("\n{LAW}flue "));
    let gate = CancelGate::new_clock_free();
    let mut baseline = FlueDgFilePerformance::from_bytes(original.as_bytes(), 512, &gate)
        .unwrap()
        .into_renderer();
    let mut whole = FlueDgFilePerformance::from_bytes(changed.as_bytes(), 512, &gate)
        .unwrap()
        .into_renderer();
    let mut partitioned = FlueDgFilePerformance::from_bytes(changed.as_bytes(), 37, &gate)
        .unwrap()
        .into_renderer();
    let base_duct = baseline.voice().dg().unwrap().load().duct();
    let duct = whole.voice().dg().unwrap().load().duct();
    assert_eq!(base_duct.inlet_gas(), duct.inlet_gas());
    assert!(duct.molecular_relaxation().is_some());
    assert_eq!(duct.total_dofs() - base_duct.total_dofs(), 3 * 8);
    let (mut a, mut b, mut c) = (vec![0.0; 512], vec![0.0; 512], vec![0.0; 512]);
    baseline.block(&mut a).unwrap();
    whole.block(&mut b).unwrap();
    for block in c.chunks_mut(37) {
        partitioned.block(block).unwrap();
    }
    assert_eq!(b, c);
    assert!(a.iter().zip(&b).any(|(a, b)| (a - b).abs() > 1e-10));
    assert_eq!(
        whole.voice().dg().unwrap().load().state(),
        partitioned.voice().dg().unwrap().load().state()
    );
    let faster = changed.replace("power 10 ", "power 1 ");
    let mut renderer = FlueDgFilePerformance::from_bytes(faster.as_bytes(), 512, &gate)
        .unwrap()
        .into_renderer();
    renderer.block(&mut a).unwrap();
    assert!(a.iter().zip(&b).any(|(a, b)| (a - b).abs() > 1e-10));
}

#[test]
fn native_molecular_gas_refuses_missing_species_extra_states_and_incompatible_laws() {
    let changed = input().replace("\nflue ", &format!("\n{LAW}flue "));
    for refused in [
        changed.replace("N2 0.8 O2 0.2", "N2 0.8 He 0.2"),
        changed.replace("1e-11 256", "1e-11 70"),
        changed.replace(
            "gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2",
            "gas 293.15 101325 0",
        ),
        changed.replace(LAW.trim(), "dg_molecular iso9613-terrestrial-air"),
    ] {
        assert!(
            FlueDgFilePerformance::from_bytes(
                refused.as_bytes(),
                37,
                &CancelGate::new_clock_free()
            )
            .is_err(),
            "{refused}"
        );
    }
}

#[test]
fn native_conical_molecular_bore_uses_its_actual_quadrature_and_accepted_pressure() {
    let original = input().replace(
        "dg_tube 0.05 0.001 open 1 2 8 50 2000 0.01 1e-11 256",
        "dg_bore 1 open 1 2 12 50 2000 0.01 1e-11 256 4\ndg_segment cone 0.001 0.002 0.05",
    );
    let changed = original.replace("\nflue ", &format!("\n{LAW}flue "));
    let gate = CancelGate::new_clock_free();
    let mut baseline = FlueDgFilePerformance::from_bytes(original.as_bytes(), 512, &gate)
        .unwrap()
        .into_renderer();
    let mut whole = FlueDgFilePerformance::from_bytes(changed.as_bytes(), 512, &gate)
        .unwrap()
        .into_renderer();
    let mut partitioned = FlueDgFilePerformance::from_bytes(changed.as_bytes(), 37, &gate)
        .unwrap()
        .into_renderer();
    let duct = whole.voice().dg().unwrap().load().duct();
    assert!(duct.spatial_viscothermal_admission().is_some());
    assert!(duct.molecular_relaxation().is_some());
    assert_eq!(duct.total_dofs(), 2 * 3 + 2 * 4 * 12 + 4 * 12);
    assert_eq!(
        duct.inlet_gas(),
        baseline.voice().dg().unwrap().load().duct().inlet_gas()
    );
    assert!(duct.outlet_area_m2() > duct.inlet_area_m2());
    let (mut a, mut b, mut c) = (vec![0.0; 512], vec![0.0; 512], vec![0.0; 512]);
    baseline.block(&mut a).unwrap();
    whole.block(&mut b).unwrap();
    for block in c.chunks_mut(37) {
        partitioned.block(block).unwrap();
    }
    assert_eq!(b, c);
    assert!(a.iter().zip(&b).any(|(a, b)| (a - b).abs() > 1e-10));
    assert_eq!(
        whole.voice().dg().unwrap().load().state(),
        partitioned.voice().dg().unwrap().load().state()
    );
    // 150 states are required: counting nodal points instead of the four
    // thermal quadrature points would incorrectly admit this 149-state budget.
    let under_budget = changed.replace("1e-11 256 4", "1e-11 149 4");
    assert!(FlueDgFilePerformance::from_bytes(under_budget.as_bytes(), 37, &gate).is_err());
}

#[test]
fn canonical_molecular_flue_binds_internal_and_exterior_species_separately() {
    for bytes in [
        include_bytes!("../../../../../examples/flue-dg-molecular.performance").as_slice(),
        include_bytes!("../../../../../examples/flue-dg-conical-molecular.performance").as_slice(),
    ] {
        let source =
            FlueDgFilePerformance::from_bytes(bytes, 37, &CancelGate::new_clock_free()).unwrap();
        let load = source.renderer().voice().dg().unwrap().load();
        let molecular = load.duct().molecular_relaxation().unwrap();
        assert_eq!(molecular.gas(), load.duct().inlet_gas());
        assert_eq!(molecular.processes().len(), 2);
        let observer = source.renderer().observer_gas().unwrap();
        assert_eq!(
            observer.law,
            crate::air_path::GasAbsorptionModel::MolecularRelaxation
        );
        assert_eq!(
            &observer.gas,
            load.radiation_record().unwrap().gas_record().state()
        );
        assert!(
            observer
                .filters
                .iter()
                .all(|filter| filter.absorption.relaxation_omitted.is_none())
        );
    }
}
