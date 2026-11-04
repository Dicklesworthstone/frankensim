use super::*;
use crate::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_thermal::crystallization::CrystallizationSolve;

fn branch(phase: SolidLiquidPhase, datum: f64, rho: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([0x75; 32]),
        phase,
        datum,
        &[280.0, 380.0].map(|temperature_k| HeatCapacityKnot {
            temperature_k,
            specific_heat_capacity_j_kg_k: 1000.0,
            bulk_density_kg_m3: rho,
        }),
        1e-6,
    )
    .unwrap()
}

fn limits() -> CrystallizationStrainLimits {
    CrystallizationStrainLimits {
        maximum_absolute_total_mandel_strain: 0.02,
        maximum_absolute_elastic_mandel_strain: 0.02,
        maximum_absolute_free_isotropic_strain: 0.01,
    }
}

fn elastic() -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(1000.0, 1e10, 4e9, ContentHash([0x76; 32]), 0.02)
        .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn g1_anisotropic_release_retains_original_energy_and_fixed_reference_mass() {
    let parent = branch(SolidLiquidPhase::Liquid, 100_000.0, 1000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 0.0, 1030.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let mut tensor = [[0.0; 6]; 6];
    for (i, value) in [9e9, 8e9, 7e9, 3e9, 2e9, 1e9].into_iter().enumerate() {
        tensor[i][i] = value;
    }
    tensor[0][3] = 5e8;
    tensor[3][0] = 5e8;
    tensor[1][5] = -2e8;
    tensor[5][1] = -2e8;
    let elastic =
        ElasticAcousticState::try_new_mandel(tensor, 1000.0, ContentHash([0x77; 32]), 0.02)
            .unwrap();
    let reference_strain = [0.001, -0.0005, 0.0002, 0.0004, 0.0, -0.0001];
    let total = [0.004, -0.002, 0.001, 0.003, -0.001, 0.002];
    let lock = CrystallizationLock::new(&pair, 350.0, 0.2, reference_strain).unwrap();
    let frozen = lock.freeze(&pair, 320.0, 0.5, &elastic, limits()).unwrap();
    assert!(std::ptr::eq(frozen.source(), &elastic));
    let strain = crate::tensor::to_mandel(&frozen.elastic_strain(&total).unwrap());
    let expected: f64 = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| 0.5 * strain[i] * tensor[i][j] * strain[j])
                .sum::<f64>()
        })
        .sum();
    let free_volume = frozen.volume_state().specific_volume_m3_kg;
    let reference_volume = 0.0013;
    assert_ne!(free_volume, reference_volume);
    let released = frozen
        .unlock_at_fixed_strain(total, reference_volume, 0.375)
        .unwrap();
    close(
        released.energy.removed_elastic_energy_density_j_m3,
        expected,
        1e-9,
    );
    close(
        released.energy.removed_elastic_energy_j_kg,
        expected * reference_volume,
        2e-12,
    );
    close(
        released.energy.thermalized_energy_j_kg,
        expected * reference_volume * 0.375,
        1e-12,
    );
    close(
        released.energy.exported_energy_j_kg,
        expected * reference_volume * 0.625,
        2e-12,
    );
    close(released.energy.balance_residual_j_kg, 0.0, 1e-12);
    assert_eq!(
        released.state.reference_specific_volume_m3_kg(),
        reference_volume
    );
    assert_eq!(released.state.previous_lock(), &lock);
    assert_eq!(released.state.release_total_strain(), &total);
    assert_eq!(
        released.state.release_volume_state(),
        &pair.volume_state(320.0, 0.5).unwrap()
    );
    assert_eq!(
        released.state.elastic_state_identity(),
        elastic.material_state_identity()
    );
    assert_eq!(released.state.stress(&total, &()), [0.0; 6]);
    assert_eq!(released.state.tangent(&total, &()), [[0.0; 6]; 6]);
    assert_eq!(released.state.stress(&[0.0; 6], &()), [0.0; 6]);
    assert!(
        released
            .state
            .stress(&[0.03; 6], &())
            .iter()
            .all(|v| v.is_nan())
    );
    assert!(
        released
            .state
            .tangent(&[f64::NAN; 6], &())
            .iter()
            .flatten()
            .all(|v| v.is_nan())
    );
}

#[test]
fn g1_released_heat_enters_original_remelting_once_and_cooling_can_relock() {
    let parent = branch(SolidLiquidPhase::Liquid, 100_000.0, 1000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 0.0, 1030.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let rate = |t: f64| Ok(if t < 350.0 { 0.01 * (350.0 - t) } else { 0.0 });
    let law = pair
        .law(&rate, ContentHash([0x78; 32]), 1.0, 0.8, 350.0)
        .unwrap();
    let remelting = law.equivalent_age_remelting().unwrap();
    let initial = law.state_at_temperature(350.0, 0.5).unwrap();
    let lock =
        CrystallizationLock::new(&pair, 350.0, initial.crystalline_mass_fraction(), [0.0; 6])
            .unwrap();
    let elastic = elastic();
    let total = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
    let fixed_v0_per_mass = 0.001;
    let released = lock
        .freeze(
            &pair,
            350.0,
            initial.crystalline_mass_fraction(),
            &elastic,
            limits(),
        )
        .unwrap()
        .unlock_at_fixed_strain(total, fixed_v0_per_mass, 1.0)
        .unwrap();
    let heat = released.energy.thermalized_energy_j_kg;
    assert!(heat > 0.0);
    assert_eq!(heat, released.energy.removed_elastic_energy_j_kg);
    assert_eq!(released.energy.exported_energy_j_kg, 0.0);
    let final_h = initial.specific_enthalpy_j_kg() + heat;
    let melted = law
        .advance_specific_enthalpy_with_remelting(
            initial,
            final_h,
            0.1,
            CrystallizationSolve::default(),
            remelting,
            || false,
        )
        .unwrap();
    assert_eq!(melted.state.specific_enthalpy_j_kg(), final_h);
    assert_eq!(melted.state.temperature_k(), 350.0);
    close(
        melted.state.crystalline_mass_fraction(),
        initial.crystalline_mass_fraction() - heat / 100_000.0,
        2e-13,
    );
    close(
        initial.specific_enthalpy_j_kg() + released.energy.removed_elastic_energy_j_kg,
        melted.state.specific_enthalpy_j_kg(),
        1e-9,
    );
    let hold = law
        .advance_specific_enthalpy_with_remelting(
            melted.state,
            final_h,
            10.0,
            CrystallizationSolve::default(),
            remelting,
            || false,
        )
        .unwrap();
    assert_eq!(hold.state, melted.state);
    assert_eq!(hold.specific_heat_j_kg, 0.0);

    let cooled = law
        .advance_specific_enthalpy_with_remelting(
            hold.state,
            final_h - 20_000.0,
            0.01,
            CrystallizationSolve::default(),
            remelting,
            || false,
        )
        .unwrap()
        .state;
    assert!(cooled.temperature_k() < 350.0);
    assert!(cooled.crystalline_mass_fraction() > melted.state.crystalline_mass_fraction());
    assert_eq!(released.state.stress(&total, &()), [0.0; 6]);
    let relocked = released
        .state
        .relock_at(
            &pair,
            cooled.temperature_k(),
            cooled.crystalline_mass_fraction(),
            total,
        )
        .unwrap();
    let refrozen = relocked
        .freeze(
            &pair,
            cooled.temperature_k(),
            cooled.crystalline_mass_fraction(),
            &elastic,
            limits(),
        )
        .unwrap();
    assert_eq!(refrozen.stress(&total, &()), [0.0; 6]);
    assert_eq!(
        refrozen.tangent(&total, &()),
        elastic.tangent(&[0.0; 6], &())
    );
    // A new stress-free reference does not reset V0/m0 or make a second heat credit.
    let released_again = refrozen
        .unlock_at_fixed_strain(total, released.state.reference_specific_volume_m3_kg(), 1.0)
        .unwrap();
    assert_eq!(released_again.energy.removed_elastic_energy_j_kg, 0.0);
    assert_eq!(released_again.energy.thermalized_energy_j_kg, 0.0);
    assert_eq!(
        released_again.state.reference_specific_volume_m3_kg(),
        fixed_v0_per_mass
    );
}

#[test]
fn g0_invalid_release_and_relock_refuse_without_changing_source_or_checkpoint() {
    let parent = branch(SolidLiquidPhase::Liquid, 100_000.0, 1000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 0.0, 1030.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let elastic = elastic();
    let lock = CrystallizationLock::new(&pair, 350.0, 0.2, [0.0; 6]).unwrap();
    let snapshot = || lock.freeze(&pair, 350.0, 0.2, &elastic, limits()).unwrap();
    let total = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
    for v in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(snapshot().unlock_at_fixed_strain(total, v, 1.0).is_err());
    }
    for fraction in [-0.01, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            snapshot()
                .unlock_at_fixed_strain(total, 0.001, fraction)
                .is_err()
        );
    }
    assert!(
        snapshot()
            .unlock_at_fixed_strain([0.03; 6], 0.001, 1.0)
            .is_err()
    );
    let tighter = ElasticAcousticState::from_bulk_and_shear(
        1000.0,
        1e10,
        4e9,
        ContentHash([0x79; 32]),
        0.001,
    )
    .unwrap();
    assert!(
        lock.freeze(&pair, 350.0, 0.2, &tighter, limits())
            .unwrap()
            .unlock_at_fixed_strain(total, 0.001, 1.0)
            .is_err()
    );

    let released = snapshot()
        .unlock_at_fixed_strain(total, 0.001, 0.0)
        .unwrap();
    assert_eq!(released.energy.thermalized_energy_j_kg, 0.0);
    assert_eq!(
        released.energy.exported_energy_j_kg,
        released.energy.removed_elastic_energy_j_kg
    );
    let other_crystal = branch(SolidLiquidPhase::Solid, 0.0, 1040.0);
    let other_pair =
        CrystallizationCaloricPair::new(&parent, &other_crystal, [280.0, 380.0]).unwrap();
    assert!(matches!(
        released.state.relock_at(&other_pair, 350.0, 0.2, total),
        Err(CrystallizationMechanicsError::PairMismatch)
    ));
    assert!(released.state.relock_at(&pair, 500.0, 0.2, total).is_err());
    assert!(released.state.relock_at(&pair, 350.0, -0.1, total).is_err());
    assert!(
        released
            .state
            .relock_at(&pair, 350.0, 0.2, [0.03; 6])
            .is_err()
    );
    assert_eq!(released.state.previous_lock(), &lock);
    assert_eq!(released.state.stress(&total, &()), [0.0; 6]);
    assert_eq!(snapshot().stress(&total, &()), elastic.stress(&total, &()));
}

#[test]
fn g3_unlocking_does_not_depend_on_absolute_caloric_energy_datum() {
    let elastic = elastic();
    let release = |datum| {
        let parent = branch(SolidLiquidPhase::Liquid, datum + 100_000.0, 1000.0);
        let crystal = branch(SolidLiquidPhase::Solid, datum, 1030.0);
        let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
        CrystallizationLock::new(&pair, 350.0, 0.2, [0.0; 6])
            .unwrap()
            .freeze(&pair, 320.0, 0.5, &elastic, limits())
            .unwrap()
            .unlock_at_fixed_strain([0.001; 6], 0.001, 0.5)
            .unwrap()
            .energy
    };
    assert_eq!(release(0.0), release(1e8));
}
