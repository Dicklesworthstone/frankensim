//! Finite source-bound heat exchange against an independent ternary oracle.

use super::*;
use fs_thermochem::solution::multicomponent::caloric::MulticomponentCaloricSolveV1;
use fs_thermochem::solution::multicomponent::caloric::equilibrium::MulticomponentEquilibriumCaloricSolveV1;

#[path = "transactions.rs"]
mod transactions;

const AMOUNTS: [f64; 3] = [2.0, 3.0, 5.0];
const GAS_CONSTANT: f64 = 8.314_462_618_153_24;

fn controls() -> MulticomponentEquilibriumCaloricSolveV1 {
    MulticomponentEquilibriumCaloricSolveV1 {
        thermal: MulticomponentCaloricSolveV1 {
            temperature_bracket: [Temperature::new(380.0), Temperature::new(490.0)],
            temperature_tolerance_k: 1e-8,
            energy_tolerance_j: 1e-6,
            entropy_tolerance_j_per_k: 1e-8,
            maximum_iterations: 100,
        },
        response: MulticomponentCaloricResponseControlV1::default(),
        maximum_temperature_step_k: 15.0,
        maximum_phase_solves: 512,
    }
}

// Each phase is a permutation of (r,q,q). This scalar common-tangent equation
// is independent of both the finite-inventory optimizer and caloric inverse.
fn phase_composition(temperature: f64) -> (f64, f64) {
    let chi = 1600.0 / temperature;
    let (mut lo, mut hi) = (1e-12, 0.1);
    for _ in 0..100 {
        let q = 0.5 * (lo + hi);
        if ((1.0_f64 - 2.0 * q) / q).ln() - chi * (1.0 - 3.0 * q) > 0.0 {
            lo = q;
        } else {
            hi = q;
        }
    }
    let q = 0.5 * (lo + hi);
    (1.0 - 2.0 * q, q)
}

fn enthalpy_oracle(temperature: f64) -> f64 {
    let sensible_factor = if temperature <= 450.0 {
        let d = temperature - 300.0;
        0.6 * d + 0.002 * d * d
    } else {
        let d = temperature - 450.0;
        135.0 + 1.2 * d + 0.003 * d * d
    };
    let sensible: f64 = AMOUNTS
        .into_iter()
        .zip(MOLAR_CP_AT_400)
        .map(|(n, cp)| n * cp * sensible_factor)
        .sum();
    let (r, q) = phase_composition(temperature);
    sensible + 10.0 * GAS_CONSTANT * 1600.0 * (2.0 * r * q + q * q)
}

fn entropy_oracle(temperature: f64) -> f64 {
    let sensible_factor = if temperature <= 450.0 {
        -0.6 * (temperature / 300.0).ln() + 0.004 * (temperature - 300.0)
    } else {
        -0.6 * 1.5_f64.ln() + 0.6 - 1.5 * (temperature / 450.0).ln() + 0.006 * (temperature - 450.0)
    };
    let sensible: f64 = AMOUNTS
        .into_iter()
        .zip(MOLAR_CP_AT_400)
        .map(|(n, cp)| n * cp * sensible_factor)
        .sum();
    let (r, q) = phase_composition(temperature);
    sensible - 10.0 * GAS_CONSTANT * (r * r.ln() + 2.0 * q * q.ln())
}

#[test]
fn g1_finite_three_phase_heat_crosses_cp_knot_and_refreshes_actual_material() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial_equilibrium = liquid
        .resolve_equilibrium(
            candidate(liquid.thermodynamics(), 400.0, AMOUNTS, true),
            || false,
        )
        .unwrap();
    let initial = caloric
        .equilibrium_caloric_state(
            initial_equilibrium.clone(),
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        )
        .unwrap();
    assert_eq!(initial.equilibrium(), &initial_equilibrium);
    assert_eq!(initial.sources(), caloric.sources());

    let heat = enthalpy_oracle(470.0) - enthalpy_oracle(400.0);
    let result = caloric
        .add_heat_equilibrium(&initial, heat, controls(), || false)
        .unwrap();
    let final_state = result.final_state();
    assert!((final_state.temperature().value() - 470.0).abs() < 2e-6);
    assert_eq!(final_state.caloric_state().component_moles(), AMOUNTS);
    assert_eq!(final_state.equilibrium().phases().len(), 3);
    assert_eq!(
        final_state.caloric_state().candidate(),
        final_state.equilibrium().equilibrium_candidate()
    );
    assert_eq!(final_state.sources(), initial.sources());
    assert!(result.balance().completed_phase_solves() > 1);
    assert!(result.balance().residual_j().abs() <= controls().thermal.energy_tolerance_j);
    assert!((result.balance().enthalpy_change_j() - heat).abs() < 2e-6);
    let expected_entropy = entropy_oracle(470.0) - entropy_oracle(400.0);
    assert!((result.entropy_change_j_per_k() - expected_entropy).abs() < 2e-6);
    near(
        caloric
            .equilibrium_enthalpy_change_j(&initial, final_state)
            .unwrap(),
        result.balance().enthalpy_change_j(),
    );
    near(
        caloric
            .equilibrium_entropy_change_j_per_k(&initial, final_state)
            .unwrap(),
        result.entropy_change_j_per_k(),
    );
    near(
        final_state.specific_heat_cp_j_kg_k(),
        final_state.caloric_state().isobaric_heat_capacity_j_per_k()
            / final_state.equilibrium().total_mass_kg(),
    );
    near(
        final_state.equilibrium().total_mass_kg(),
        initial.equilibrium().total_mass_kg(),
    );
    let (r, q) = phase_composition(470.0);
    for (i, phase) in final_state.equilibrium().phases().iter().enumerate() {
        let expected: Vec<f64> = (0..3).map(|j| if j == i { r } else { q }).collect();
        for (&actual, &reference) in phase.mole_fractions().iter().zip(&expected) {
            assert!((actual - reference).abs() < 2e-9);
        }
        near(
            phase.liquid().density_kg_m3(),
            source_value(phase.mole_fractions(), DENSITIES),
        );
        near(
            phase.liquid().thermal_conductivity_w_m_k(),
            source_value(phase.mole_fractions(), CONDUCTIVITIES),
        );
        assert_ne!(
            phase.liquid().density_kg_m3(),
            initial.equilibrium().phases()[i].liquid().density_kg_m3()
        );
    }

    let cooled = caloric
        .add_heat_equilibrium(final_state, -heat, controls(), || false)
        .unwrap();
    assert!((cooled.final_state().temperature().value() - 400.0).abs() < 2e-6);
    assert!((cooled.entropy_change_j_per_k() + expected_entropy).abs() < 2e-6);
    assert_eq!(initial.equilibrium(), &initial_equilibrium);
}

#[test]
fn g3_final_source_refusal_keeps_original_finite_equilibrium_intact() {
    // The original chart admits the warm, less-separated phases. Cooling
    // increases the A-rich composition beyond its actual source support.
    let card = card(0.96);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = caloric
        .equilibrium_caloric_state(
            liquid
                .resolve_equilibrium(
                    candidate(liquid.thermodynamics(), 470.0, AMOUNTS, true),
                    || false,
                )
                .unwrap(),
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        )
        .unwrap();
    let original = initial.clone();
    let heat = enthalpy_oracle(400.0) - enthalpy_oracle(470.0);
    assert!(matches!(
        caloric.add_heat_equilibrium(&initial, heat, controls(), || false),
        Err(MulticomponentLiquidCaloricError::Material(
            MulticomponentLiquidError::Source(_)
        ))
    ));
    assert_eq!(initial, original);
    assert_eq!(caloric.sources(), original.sources());
}

#[test]
fn g3_foreign_caloric_owner_and_last_poll_cancellation_publish_no_endpoint() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = caloric
        .equilibrium_caloric_state(
            liquid
                .resolve_equilibrium(
                    candidate(liquid.thermodynamics(), 400.0, AMOUNTS, true),
                    || false,
                )
                .unwrap(),
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        )
        .unwrap();
    let original = initial.clone();
    let foreign = MulticomponentLiquidCaloricLaw::try_new(
        &liquid,
        (0..3)
            .map(|i| (species(i), cp_source(i, 100_000.0, 310.0)))
            .collect(),
        "p",
    )
    .unwrap();
    assert!(matches!(
        foreign.add_heat_equilibrium(&initial, 100.0, controls(), || false),
        Err(MulticomponentLiquidCaloricError::Invalid(_))
    ));
    assert!(
        foreign
            .equilibrium_enthalpy_change_j(&initial, &initial)
            .is_err()
    );
    assert!(
        foreign
            .equilibrium_entropy_change_j_per_k(&initial, &initial)
            .is_err()
    );
    let heat = enthalpy_oracle(470.0) - enthalpy_oracle(400.0);
    let mut polls = 0;
    caloric
        .add_heat_equilibrium(&initial, heat, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut current = 0;
        assert!(
            caloric
                .add_heat_equilibrium(&initial, heat, controls(), || {
                    current += 1;
                    current == stop
                })
                .is_err()
        );
        assert_eq!(initial, original);
    }
}
