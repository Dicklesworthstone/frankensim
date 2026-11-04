use super::*;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn law(solute: &str, saturation: f64, enthalpy: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new(solute).unwrap(),
        Temperature::new(298.15),
        AmountConcentration::new(saturation),
        enthalpy,
        [Temperature::new(280.0), Temperature::new(350.0)],
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn concentrating_and_diluting_preserve_both_phases_and_only_dissolved_osmotic_particles() {
    let initial = FiniteSolubilitySolution::new(
        Temperature::new(298.15),
        18e-6,
        Volume::new(1e-3),
        vec![
            (law("neutral-a", 100.0, 0.0), Amount::new(0.2)),
            (law("neutral-b", 30.0, 0.0), Amount::new(0.02)),
        ],
    )
    .unwrap();
    assert_eq!(initial.solution().solute_moles(), &[0.1, 0.02]);
    assert_eq!(initial.phase_partitions()[0].solid.value(), 0.1);
    assert_eq!(initial.phase_partitions()[1].solid.value(), 0.0);
    near(
        initial.solution().osmotic_pressure_pa().unwrap(),
        R * 298.15 * 120.0,
        1e-9,
    );
    near(
        initial.solution().water_activity().unwrap(),
        (-18e-6_f64 * 120.0).exp(),
        2e-15,
    );

    let concentrated = initial.with_solvent_volume_m3(5e-4).unwrap();
    assert_eq!(concentrated.solution().solute_moles(), &[0.05, 0.015]);
    assert!(std::ptr::eq(
        initial.components(),
        concentrated.components()
    ));
    near(
        concentrated.phase_partitions()[0].solid.value(),
        0.15,
        1e-16,
    );
    near(
        concentrated.phase_partitions()[1].solid.value(),
        0.005,
        1e-17,
    );
    near(
        concentrated.solution().osmotic_pressure_pa().unwrap(),
        R * 298.15 * 130.0,
        1e-9,
    );
    let diluted = concentrated.with_solvent_volume_m3(1e-2).unwrap();
    assert_eq!(diluted.solution().solute_moles(), &[0.2, 0.02]);
    assert!(
        diluted
            .phase_partitions()
            .iter()
            .all(|p| p.solid.value() == 0.0)
    );
    for state in [&initial, &concentrated, &diluted] {
        for (component, phase) in state.components().iter().zip(state.phase_partitions()) {
            near(
                phase.dissolved.value() + phase.solid.value(),
                component.total().value(),
                1e-16,
            );
        }
    }
    let returned = diluted.with_solvent_volume_m3(1e-3).unwrap();
    assert_eq!(returned, initial);
}

#[test]
fn undersaturation_and_zero_inventory_recover_original_ideal_solution() {
    for amount in [0.0, 0.01] {
        let soluble = FiniteSolubilitySolution::new(
            Temperature::new(298.15),
            18e-6,
            Volume::new(1e-3),
            vec![(law("neutral-a", 100.0, 0.0), Amount::new(amount))],
        )
        .unwrap();
        let original = IdealDiluteSolution::new(298.15, 18e-6, 1e-3, vec![amount]).unwrap();
        assert_eq!(soluble.solution(), &AqueousSolution::from(original));
        assert_eq!(soluble.phase_partitions()[0].solid.value(), 0.0);
        let diluted = soluble.with_solvent_volume_m3(2e-3).unwrap();
        assert_eq!(diluted.solution().solute_moles(), &[amount]);
    }
}

#[test]
fn prescribed_temperature_changes_solubility_and_redissolves_the_original_stock() {
    let initial = FiniteSolubilitySolution::new(
        Temperature::new(298.15),
        18e-6,
        Volume::new(1e-3),
        vec![(law("neutral-a", 100.0, 12_000.0), Amount::new(0.2))],
    )
    .unwrap();
    let heated = initial.with_temperature(Temperature::new(350.0)).unwrap();
    assert_eq!(heated.solution().solute_moles(), &[0.2]);
    assert_eq!(heated.phase_partitions()[0].solid.value(), 0.0);
    assert_eq!(heated.solution().temperature_k(), 350.0);
    near(
        heated.components()[0].saturation_concentration().value(),
        100.0 * (12_000.0 / R * (1.0 / 298.15 - 1.0 / 350.0)).exp(),
        1e-11,
    );
    let cooled = heated.with_temperature(Temperature::new(298.15)).unwrap();
    assert_eq!(cooled, initial);
    assert!(initial.with_temperature(Temperature::new(351.0)).is_err());
    assert_eq!(initial.phase_partitions()[0].solid.value(), 0.1);
}

#[test]
fn species_layout_and_failed_volume_trials_do_not_lose_original_inventories() {
    let valid = law("neutral-a", 100.0, 0.0);
    let build = |components| {
        FiniteSolubilitySolution::new(
            Temperature::new(298.15),
            18e-6,
            Volume::new(1e-3),
            components,
        )
    };
    assert!(build(Vec::new()).is_err());
    assert!(build(vec![(valid.clone(), Amount::new(1.0)); 65]).is_err());
    assert!(
        build(vec![
            (valid.clone(), Amount::new(0.1)),
            (valid.clone(), Amount::new(0.2)),
        ])
        .is_err()
    );
    let other_solvent = IdealSolubilityV1::new(
        SpeciesId::new("other-solvent").unwrap(),
        SpeciesId::new("neutral-b").unwrap(),
        Temperature::new(298.15),
        AmountConcentration::new(100.0),
        0.0,
        [Temperature::new(280.0), Temperature::new(350.0)],
    )
    .unwrap();
    assert!(
        build(vec![
            (valid.clone(), Amount::new(0.1)),
            (other_solvent, Amount::new(0.1)),
        ])
        .is_err()
    );
    assert!(build(vec![(valid.clone(), Amount::new(-0.1))]).is_err());
    let state = build(vec![(valid, Amount::new(0.2))]).unwrap();
    let before = state.clone();
    for volume in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(state.with_solvent_volume_m3(volume).is_err());
        assert_eq!(state, before);
    }
}

#[test]
fn transported_totals_repartition_in_prepared_storage_and_continue_with_updated_stocks() {
    let original = FiniteSolubilitySolution::new(
        Temperature::new(298.15),
        18e-6,
        Volume::new(1e-3),
        vec![
            (law("neutral-a", 100.0, 0.0), Amount::new(0.2)),
            (law("neutral-b", 30.0, 0.0), Amount::new(0.02)),
        ],
    )
    .unwrap();
    let mut evolved = original.independent_clone();
    let trial = original.partition_at_inventory(1, 5e-4, 0.08).unwrap();
    assert_eq!(trial.total_solute.value(), 0.08);
    assert_eq!(trial.dissolved.value(), 0.015);
    assert!(original.partition_at_inventory(2, 5e-4, 0.08).is_err());
    let component_storage = evolved.components().as_ptr();
    let phase_storage = evolved.phase_partitions().as_ptr();
    let dissolved_storage = evolved.solution().solute_moles().as_ptr();
    evolved.replace_inventories(5e-4, &[0.04, 0.08]).unwrap();
    assert_eq!(evolved.solution().solute_moles(), &[0.04, 0.015]);
    assert_eq!(evolved.phase_partitions()[0].solid.value(), 0.0);
    near(evolved.phase_partitions()[1].solid.value(), 0.065, 1e-16);
    near(
        evolved.solution().osmotic_pressure_pa().unwrap(),
        R * 298.15 * 110.0,
        1e-9,
    );
    assert_eq!(evolved.components()[0].total().value(), 0.04);
    assert_eq!(evolved.components()[1].total().value(), 0.08);

    evolved.replace_inventories(2e-3, &[0.3, 0.01]).unwrap();
    assert_eq!(evolved.solution().solute_moles(), &[0.2, 0.01]);
    near(evolved.phase_partitions()[0].solid.value(), 0.1, 1e-16);
    assert_eq!(evolved.phase_partitions()[1].solid.value(), 0.0);
    assert_eq!(evolved.components().as_ptr(), component_storage);
    assert_eq!(evolved.phase_partitions().as_ptr(), phase_storage);
    assert_eq!(
        evolved.solution().solute_moles().as_ptr(),
        dissolved_storage
    );
    let diluted = evolved.with_solvent_volume_m3(1e-2).unwrap();
    assert_eq!(diluted.solution().solute_moles(), &[0.3, 0.01]);
    assert_eq!(diluted.components()[0].total().value(), 0.3);
    assert_eq!(original.components()[0].total().value(), 0.2);
    assert_eq!(original.solution().solute_moles(), &[0.1, 0.02]);
    assert!(
        evolved.heap_bytes()
            >= evolved.solution().heap_bytes()
                + core::mem::size_of_val(evolved.components())
                + core::mem::size_of_val(evolved.phase_partitions())
    );
}

#[test]
fn failed_or_shared_inventory_replacement_leaves_every_phase_and_view_unchanged() {
    let mut state = FiniteSolubilitySolution::new(
        Temperature::new(298.15),
        18e-6,
        Volume::new(1e-3),
        vec![
            (law("neutral-a", 100.0, 0.0), Amount::new(0.2)),
            (law("neutral-b", 30.0, 0.0), Amount::new(0.02)),
        ],
    )
    .unwrap();
    let before = state.independent_clone();
    for (volume, totals) in [
        (0.0, [0.04, 0.08]),
        (1e-3, [0.04, -0.01]),
        (1e-3, [0.04, f64::NAN]),
        (1e-3, [0.04, f64::INFINITY]),
    ] {
        assert!(state.replace_inventories(volume, &totals).is_err());
        assert_eq!(state, before);
    }
    assert!(state.replace_inventories(1e-3, &[0.04]).is_err());
    assert_eq!(state, before);
    let shared = state.clone();
    assert!(state.replace_inventories(5e-4, &[0.04, 0.08]).is_err());
    assert_eq!(state, shared);
    let mut independent = shared.independent_clone();
    independent
        .replace_inventories(5e-4, &[0.04, 0.08])
        .unwrap();
    assert_eq!(shared, before);

    let mut huge_saturation = FiniteSolubilitySolution::new(
        Temperature::new(298.15),
        18e-6,
        Volume::new(1e-3),
        vec![(law("neutral-c", 1e9, 0.0), Amount::new(0.0))],
    )
    .unwrap();
    let zero = huge_saturation.independent_clone();
    assert!(huge_saturation.replace_inventories(1e-3, &[1e6]).is_err());
    assert_eq!(huge_saturation, zero);
}
