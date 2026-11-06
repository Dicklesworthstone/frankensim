//! Actual material phases after conserved relaxation and unequal-temperature mixing.

use super::*;
use fs_material::liquid::solution::multicomponent::caloric::ResolvedMulticomponentLiquidEquilibriumCaloricState;

fn branch(
    law: &MulticomponentLiquidCaloricLaw<'_, '_>,
    temperature: f64,
    amounts: [f64; 3],
) -> ResolvedMulticomponentLiquidEquilibriumCaloricState {
    law.equilibrium_caloric_state(
        law.liquid_law()
            .resolve_equilibrium(
                candidate(
                    law.liquid_law().thermodynamics(),
                    temperature,
                    amounts,
                    true,
                ),
                || false,
            )
            .unwrap(),
        MulticomponentCaloricResponseControlV1::default(),
        || false,
    )
    .unwrap()
}

fn mixing_oracle(temperature: f64) -> [f64; 2] {
    let (r, q) = phase_composition(temperature);
    [
        GAS_CONSTANT * 1600.0 * (2.0 * r * q + q * q),
        -GAS_CONSTANT * (r * r.ln() + 2.0 * q * q.ln()),
    ]
}

fn pure_oracle(temperature: f64) -> [f64; 2] {
    let mixing = mixing_oracle(temperature);
    [
        (enthalpy_oracle(temperature) - 10.0 * mixing[0]) / 840.0,
        (entropy_oracle(temperature) - 10.0 * mixing[1]) / 840.0,
    ]
}

fn homogeneous_oracle(temperature: f64) -> [f64; 2] {
    let pure = pure_oracle(temperature);
    [
        840.0 * pure[0] + 10.0 * GAS_CONSTANT * 1600.0 * (0.2 * 0.3 + 0.2 * 0.5 + 0.3 * 0.5),
        840.0 * pure[1]
            - 10.0 * GAS_CONSTANT * (0.2 * 0.2_f64.ln() + 0.3 * 0.3_f64.ln() + 0.5 * 0.5_f64.ln()),
    ]
}

fn verify_actual_sources(
    law: &MulticomponentLiquidCaloricLaw<'_, '_>,
    state: &ResolvedMulticomponentLiquidEquilibriumCaloricState,
    card: &MaterialCard,
) {
    let resolved = state.equilibrium();
    assert_eq!(state.sources(), law.sources());
    assert_eq!(
        state.caloric_state().candidate(),
        resolved.equilibrium_candidate()
    );
    let independently_resolved = law
        .liquid_law()
        .resolve_equilibrium(state.caloric_state().candidate().clone(), || false)
        .unwrap();
    assert_eq!(resolved, &independently_resolved);
    let mut actual_volume = 0.0;
    for phase in resolved.phases() {
        near(
            phase.liquid().density_kg_m3(),
            source_value(phase.mole_fractions(), DENSITIES),
        );
        near(
            phase.liquid().thermal_conductivity_w_m_k(),
            source_value(phase.mole_fractions(), CONDUCTIVITIES),
        );
        actual_volume += phase.mass_kg() / phase.liquid().density_kg_m3();
        assert_eq!(
            phase.liquid().material().card_identity(),
            card.content_hash()
        );
        assert_eq!(
            phase.liquid().material().state_coordinate("T"),
            Some(state.temperature().value())
        );
        assert_eq!(
            phase.liquid().material().state_coordinate("p"),
            Some(100_000.0)
        );
        for property in phase.liquid().material().properties() {
            let original = law
                .liquid_law()
                .reference_state()
                .material()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
        }
    }
    near(resolved.total_volume_m3(), actual_volume);
    near(
        state.specific_heat_cp_j_kg_k(),
        state.caloric_state().isobaric_heat_capacity_j_per_k() / resolved.total_mass_kg(),
    );
    assert_ne!(
        resolved.phases()[0].liquid().density_kg_m3(),
        resolved.phases()[1].liquid().density_kg_m3()
    );
}

#[test]
fn g1_adiabatic_demixing_crosses_original_cp_knot_and_resolves_actual_phase_sources() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = caloric
        .homogeneous(Temperature::new(420.0), &AMOUNTS, || false)
        .unwrap();
    let initial_snapshot = initial.clone();
    let phase_branch = branch(&caloric, 430.0, AMOUNTS);
    let result = caloric
        .relax_homogeneous(&initial, &phase_branch, 0.0, controls(), || false)
        .unwrap();
    let input = homogeneous_oracle(420.0);
    let mut bracket = [420.0, 490.0];
    for _ in 0..100 {
        let t = 0.5 * bracket[0] + 0.5 * bracket[1];
        if enthalpy_oracle(t) > input[0] {
            bracket[1] = t;
        } else {
            bracket[0] = t;
        }
    }
    let expected = 0.5 * bracket[0] + 0.5 * bracket[1];
    assert!(
        expected > 450.0,
        "demixing must carry the input through the original Cp knot"
    );
    assert!((result.final_state().temperature().value() - expected).abs() < 2e-6);
    assert!(result.balance().residual_j().abs() <= controls().thermal.energy_tolerance_j);
    assert!(
        (result.final_state().caloric_state().enthalpy_j() - initial.caloric_state().enthalpy_j())
            .abs()
            < 2e-6
    );
    let entropy = entropy_oracle(expected) - input[1];
    assert!(entropy > 0.0);
    assert!((result.entropy_change_j_per_k() - entropy).abs() < 2e-6);
    assert_eq!(
        result.final_state().caloric_state().component_moles(),
        AMOUNTS
    );
    near(
        result.final_state().equilibrium().total_mass_kg(),
        initial.liquid().mass_kg(),
    );
    verify_actual_sources(&caloric, result.final_state(), &card);
    assert_eq!(initial, initial_snapshot);
    assert_eq!(phase_branch.temperature().value(), 430.0);
}

#[test]
fn g1_mixing_different_thermal_phase_inventories_preserves_energy_entropy_and_sources() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let left = branch(&caloric, 400.0, AMOUNTS);
    let right = caloric
        .homogeneous(Temperature::new(470.0), &[1.0; 3], || false)
        .unwrap();
    let phase_branch = branch(&caloric, 440.0, [3.0, 4.0, 6.0]);
    let result = caloric
        .mix_equilibrium(&left, &right, &phase_branch, 0.0, controls(), || false)
        .unwrap();
    let initial_h = enthalpy_oracle(400.0) + 220.0 * pure_oracle(470.0)[0] + GAS_CONSTANT * 1600.0;
    let initial_s =
        entropy_oracle(400.0) + 220.0 * pure_oracle(470.0)[1] + 3.0 * GAS_CONSTANT * 3.0_f64.ln();
    let mut interval = [400.0, 490.0];
    for _ in 0..100 {
        let t = 0.5 * interval[0] + 0.5 * interval[1];
        if 1060.0 * pure_oracle(t)[0] + 13.0 * mixing_oracle(t)[0] > initial_h {
            interval[1] = t;
        } else {
            interval[0] = t;
        }
    }
    let expected = 0.5 * interval[0] + 0.5 * interval[1];
    let entropy = 1060.0 * pure_oracle(expected)[1] + 13.0 * mixing_oracle(expected)[1] - initial_s;
    assert!((result.final_state().temperature().value() - expected).abs() < 2e-6);
    assert!((result.entropy_change_j_per_k() - entropy).abs() < 2e-6);
    assert!(entropy > 0.0);
    assert!(
        (result.final_state().caloric_state().enthalpy_j()
            - left.caloric_state().enthalpy_j()
            - right.caloric_state().enthalpy_j())
        .abs()
            < 2e-6
    );
    assert_eq!(
        result.final_state().caloric_state().component_moles(),
        &[3.0, 4.0, 6.0]
    );
    near(
        result.final_state().equilibrium().total_mass_kg(),
        left.equilibrium().total_mass_kg() + right.liquid().mass_kg(),
    );
    verify_actual_sources(&caloric, result.final_state(), &card);
}

#[test]
fn g3_source_refusal_foreign_owner_and_final_cancellation_leave_all_inputs_intact() {
    let restricted = card(0.96);
    let liquid = law(&restricted, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = caloric
        .homogeneous(Temperature::new(450.0), &AMOUNTS, || false)
        .unwrap();
    let phase_branch = branch(&caloric, 470.0, AMOUNTS);
    let originals = (initial.clone(), phase_branch.clone());
    let heat = enthalpy_oracle(400.0) - homogeneous_oracle(450.0)[0];
    let numeric = caloric
        .model()
        .relax_homogeneous(
            initial.caloric_state(),
            phase_branch.caloric_state(),
            heat,
            controls(),
            || false,
        )
        .unwrap();
    assert!((numeric.final_state().temperature().value() - 400.0).abs() < 2e-6);
    assert!(matches!(
        caloric.relax_homogeneous(&initial, &phase_branch, heat, controls(), || false),
        Err(MulticomponentLiquidCaloricError::Material(
            MulticomponentLiquidError::Source(_)
        ))
    ));
    assert_eq!((initial, phase_branch), originals);

    let full = card(1.0);
    let liquid = law(&full, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = caloric
        .homogeneous(Temperature::new(420.0), &AMOUNTS, || false)
        .unwrap();
    let phase_branch = branch(&caloric, 430.0, AMOUNTS);
    let originals = (initial.clone(), phase_branch.clone());
    let foreign = MulticomponentLiquidCaloricLaw::try_new(
        &liquid,
        (0..3)
            .map(|i| (species(i), cp_source(i, 100_000.0, 310.0)))
            .collect(),
        "p",
    )
    .unwrap();
    assert!(matches!(
        foreign.relax_homogeneous(&initial, &phase_branch, 0.0, controls(), || false),
        Err(MulticomponentLiquidCaloricError::Invalid(_))
    ));
    let mut count = 0;
    caloric
        .relax_homogeneous(&initial, &phase_branch, 0.0, controls(), || {
            count += 1;
            false
        })
        .unwrap();
    for stop in [1, count / 2, count] {
        let mut polls = 0;
        assert!(
            caloric
                .relax_homogeneous(&initial, &phase_branch, 0.0, controls(), || {
                    polls += 1;
                    polls == stop
                })
                .is_err()
        );
        assert_eq!(polls, stop);
        assert_eq!((initial.clone(), phase_branch.clone()), originals);
    }
}
