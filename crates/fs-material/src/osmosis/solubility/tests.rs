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
