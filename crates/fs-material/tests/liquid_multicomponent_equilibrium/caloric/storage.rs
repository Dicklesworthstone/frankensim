//! Same-branch phase redistribution on native mass-specific thermal storage.

use super::*;
use fs_material::liquid::solution::multicomponent::caloric::ResolvedMulticomponentLiquidEquilibriumCaloricState;
use fs_thermal::enthalpy::ThermalLaw;

#[path = "storage/conduction.rs"]
mod conduction;
#[path = "storage/stream.rs"]
mod stream;

fn initial(
    caloric: &MulticomponentLiquidCaloricLaw<'_, '_>,
    temperature: f64,
    amounts: [f64; 3],
) -> ResolvedMulticomponentLiquidEquilibriumCaloricState {
    let equilibrium = caloric
        .liquid_law()
        .resolve_equilibrium(
            candidate(caloric.model().thermodynamics(), temperature, amounts, true),
            || false,
        )
        .unwrap();
    caloric
        .equilibrium_caloric_state(equilibrium, Default::default(), || false)
        .unwrap()
}

#[test]
fn g1_equilibrium_storage_tracks_original_phase_energy_entropy_and_heat_capacity() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let original = initial(&caloric, 400.0, AMOUNTS);
    let retained = original.clone();
    let storage = caloric
        .equilibrium_thermal_storage(&original, controls(), || false)
        .unwrap();
    let mass: f64 = AMOUNTS.into_iter().zip(MASSES).map(|(n, m)| n * m).sum();
    assert_eq!(storage.sources(), caloric.sources());
    assert_eq!(storage.component_moles(), AMOUNTS);
    near(storage.mass_kg(), mass);
    assert_eq!(storage.specific_enthalpy_domain_j_kg()[0], 0.0);
    for t in [380.0, 400.0, 440.0, 450.0, 470.0, 490.0] {
        let h = storage
            .specific_enthalpy_j_kg(Temperature::new(t), || false)
            .unwrap();
        let expected_h = (enthalpy_oracle(t) - enthalpy_oracle(380.0)) / mass;
        assert!(
            (h - expected_h).abs() < 3e-6,
            "T={t}, h={h}, expected={expected_h}"
        );
        let state = storage
            .caloric_state_at_specific_enthalpy(h, || false)
            .unwrap();
        assert!((state.temperature().value() - t).abs() < 2e-7);
        assert_eq!(state.candidate().phases().len(), 3);
        assert_eq!(state.component_moles(), AMOUNTS);
        let entropy = storage
            .specific_entropy_change_j_kg_k(0.0, h)
            .unwrap()
            .unwrap();
        let expected_s = (entropy_oracle(t) - entropy_oracle(380.0)) / mass;
        assert!((entropy - expected_s).abs() < 1e-7);
        assert_eq!(
            storage.specific_entropy_change_j_kg_k(h, h).unwrap(),
            Some(0.0)
        );
        if t != 380.0 && t != 490.0 && t != 450.0 {
            let cp = storage
                .specific_heat_cp_j_kg_k(Temperature::new(t), || false)
                .unwrap();
            let derivative =
                (enthalpy_oracle(t + 1e-3) - enthalpy_oracle(t - 1e-3)) / (2e-3 * mass);
            assert!(
                (cp - derivative).abs() < 1e-5,
                "T={t}, Cp={cp}, oracle={derivative}"
            );
            assert!(state.redistribution_heat_capacity_j_per_k() > 0.0);
            assert!(1.0 / cp <= storage.maximum_temperature_slope());
        }
    }
    assert_eq!(original, retained);
}

#[test]
fn g3_equilibrium_storage_normalization_preserves_specific_physics() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let small = caloric
        .equilibrium_thermal_storage(&initial(&caloric, 400.0, AMOUNTS), controls(), || false)
        .unwrap();
    let large = caloric
        .equilibrium_thermal_storage(
            &initial(&caloric, 400.0, AMOUNTS.map(|n| 7.0 * n)),
            controls(),
            || false,
        )
        .unwrap();
    near(large.mass_kg(), 7.0 * small.mass_kg());
    let t = Temperature::new(470.0);
    let h = small.specific_enthalpy_j_kg(t, || false).unwrap();
    assert!((h - large.specific_enthalpy_j_kg(t, || false).unwrap()).abs() < 3e-6);
    assert!((small.temperature_k(h).unwrap() - large.temperature_k(h).unwrap()).abs() < 2e-7);
    let ds_small = small
        .specific_entropy_change_j_kg_k(0.0, h)
        .unwrap()
        .unwrap();
    let ds_large = large
        .specific_entropy_change_j_kg_k(0.0, h)
        .unwrap()
        .unwrap();
    assert!((ds_small - ds_large).abs() < 1e-7);
    near(
        small.maximum_temperature_slope(),
        large.maximum_temperature_slope(),
    );
}

#[test]
fn g3_phase_loss_source_domain_and_work_budget_cannot_publish_storage() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let original = initial(&caloric, 400.0, AMOUNTS);
    let retained = original.clone();
    let mut unsupported = controls();
    unsupported.thermal.temperature_bracket[1] = Temperature::new(595.0);
    assert!(
        caloric
            .equilibrium_thermal_storage(&original, unsupported, || false)
            .is_err()
    );
    let mut exhausted = controls();
    exhausted.maximum_phase_solves = 1;
    assert!(
        caloric
            .equilibrium_thermal_storage(&original, exhausted, || false)
            .is_err()
    );
    let foreign = MulticomponentLiquidCaloricLaw::try_new(
        &liquid,
        (0..3)
            .map(|i| (species(i), cp_source(i, 100_000.0, 310.0)))
            .collect(),
        "p",
    )
    .unwrap();
    assert!(matches!(
        foreign.equilibrium_thermal_storage(&original, controls(), || false),
        Err(MulticomponentLiquidCaloricError::Invalid(_))
    ));
    assert_eq!(original, retained);
}

#[test]
fn g4_final_poll_cancellation_and_extrapolation_preserve_the_original_branch() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let original = initial(&caloric, 400.0, AMOUNTS);
    let retained = original.clone();
    let mut complete_polls = 0;
    let storage = caloric
        .equilibrium_thermal_storage(&original, controls(), || {
            complete_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert_eq!(
        caloric.equilibrium_thermal_storage(&original, controls(), || {
            polls += 1;
            polls == complete_polls
        }),
        Err(MulticomponentLiquidCaloricError::Cancelled)
    );
    assert_eq!(polls, complete_polls);
    let [lo, hi] = storage.specific_enthalpy_domain_j_kg();
    for h in [lo.next_down(), hi.next_up(), f64::NAN, f64::INFINITY] {
        assert!(storage.temperature_k(h).is_err());
    }
    for t in [379.0, 491.0, f64::NAN] {
        assert!(
            storage
                .specific_enthalpy_j_kg(Temperature::new(t), || false)
                .is_err()
        );
    }
    assert_eq!(original, retained);
}
