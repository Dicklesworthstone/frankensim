//! True ternary coexistence joins the two original sensible thermal branches.

use super::*;
use fs_material::liquid::solution::multicomponent::caloric::equilibrium::storage::{
    MulticomponentLiquidEquilibriumThermalStorage,
    coexistence::MulticomponentLiquidCoexistenceThermalStorage,
};
use fs_thermochem::solution::multicomponent::caloric::equilibrium::coexistence::MulticomponentCaloricCoexistenceControlV1;

#[path = "coexistence/conduction.rs"]
mod conduction;

const JOIN_FEED: [f64; 3] = [1.0; 3];

fn transition_temperature() -> f64 {
    1600.0 / (4.0 * 2.0_f64.ln())
}

fn endpoint(
    caloric: &MulticomponentLiquidCaloricLaw<'_, '_>,
    temperature: f64,
    split: bool,
) -> ResolvedMulticomponentLiquidEquilibriumCaloricState {
    let seed = if split {
        MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![
                vec![0.9, 0.05, 0.05],
                vec![0.05, 0.9, 0.05],
                vec![0.05, 0.05, 0.9],
            ],
            phase_weights: vec![1.0; 3],
        }
    } else {
        MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![vec![1.0 / 3.0; 3]],
            phase_weights: vec![1.0],
        }
    };
    let candidate = caloric
        .model()
        .thermodynamics()
        .equilibrium_candidate(
            Temperature::new(temperature),
            &JOIN_FEED,
            &seed,
            MulticomponentEquilibriumControlV1 {
                chemical_potential_tolerance: 2e-13,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    let resolved = caloric
        .liquid_law()
        .resolve_equilibrium(candidate, || false)
        .unwrap();
    caloric
        .equilibrium_caloric_state(resolved, Default::default(), || false)
        .unwrap()
}

fn branches(
    caloric: &MulticomponentLiquidCaloricLaw<'_, '_>,
    common_temperature: f64,
) -> [MulticomponentLiquidEquilibriumThermalStorage; 2] {
    let mut low_controls = controls();
    low_controls.thermal.temperature_bracket = [
        Temperature::new(380.0),
        Temperature::new(common_temperature),
    ];
    let mut high_controls = controls();
    high_controls.thermal.temperature_bracket = [
        Temperature::new(common_temperature),
        Temperature::new(595.0),
    ];
    [
        caloric
            .equilibrium_thermal_storage(&endpoint(caloric, 400.0, true), low_controls, || false)
            .unwrap(),
        caloric
            .equilibrium_thermal_storage(&endpoint(caloric, 590.0, false), high_controls, || false)
            .unwrap(),
    ]
}

fn joined(
    caloric: &MulticomponentLiquidCaloricLaw<'_, '_>,
) -> MulticomponentLiquidCoexistenceThermalStorage {
    let [low, high] = branches(caloric, transition_temperature());
    low.join_coexistence(&high, Default::default(), || false)
        .unwrap()
}

fn phase_fractions(temperature: f64) -> (f64, f64) {
    if temperature == transition_temperature() {
        return (2.0 / 3.0, 1.0 / 6.0);
    }
    let chi = 1600.0 / temperature;
    let (mut lo, mut hi) = (1e-14, 1.0 / 6.0);
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

fn source_enthalpy(temperature: f64) -> f64 {
    let sensible = if temperature <= 450.0 {
        let d = temperature - 300.0;
        220.0 * (0.6 * d + 0.002 * d * d)
    } else {
        let d = temperature - 450.0;
        220.0 * (135.0 + 1.2 * d + 0.003 * d * d)
    };
    let mixing = if temperature > transition_temperature() {
        GAS_CONSTANT * 1600.0
    } else {
        let (r, q) = phase_fractions(temperature);
        3.0 * GAS_CONSTANT * 1600.0 * (2.0 * r * q + q * q)
    };
    sensible + mixing
}

fn source_entropy(temperature: f64) -> f64 {
    let sensible = if temperature <= 450.0 {
        220.0 * (-0.6 * (temperature / 300.0).ln() + 0.004 * (temperature - 300.0))
    } else {
        220.0
            * (-0.6 * 1.5_f64.ln() + 0.6 - 1.5 * (temperature / 450.0).ln()
                + 0.006 * (temperature - 450.0))
    };
    let mixing = if temperature > transition_temperature() {
        3.0 * GAS_CONSTANT * 3.0_f64.ln()
    } else {
        let (r, q) = phase_fractions(temperature);
        -3.0 * GAS_CONSTANT * (r * r.ln() + 2.0 * q * q.ln())
    };
    sensible + mixing
}

#[test]
fn g1_joined_storage_has_a_true_latent_interval_with_original_species_and_entropy() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let joined = joined(&caloric);
    let mass: f64 = MASSES.into_iter().sum();
    let transition = transition_temperature();
    let [lo, hi] = joined.coexistence_specific_enthalpy_domain_j_kg();
    assert_eq!(lo, 0.0);
    assert!(joined.specific_enthalpy_domain_j_kg()[0] < 0.0);
    assert_eq!(joined.sources(), caloric.sources());
    assert_eq!(joined.component_moles(), JOIN_FEED);
    near(joined.mass_kg(), mass);
    assert!((hi - 400.0 * GAS_CONSTANT / mass).abs() < 1e-6);
    assert!(
        joined
            .specific_enthalpy_j_kg(Temperature::new(transition), || false)
            .is_err()
    );
    assert!(
        joined
            .specific_heat_cp_j_kg_k(Temperature::new(transition), || false)
            .is_err()
    );
    for fraction in [0.0, 0.25, 0.75, 1.0] {
        let h = fraction * hi;
        let state = joined
            .caloric_state_at_specific_enthalpy(h, || false)
            .unwrap();
        assert_eq!(state.temperature().value(), transition);
        assert_eq!(state.temperature_slope_kg_k_j(), 0.0);
        assert!((state.upper_branch_mass_fraction() - fraction).abs() < 4.0 * f64::EPSILON);
        let phases = state.coexistence_state().unwrap();
        assert_eq!(
            phases.phases().count(),
            if fraction == 0.0 {
                3
            } else if fraction == 1.0 {
                1
            } else {
                4
            }
        );
        let mut amounts = [0.0; 3];
        for phase in phases.phases() {
            for (sum, amount) in amounts.iter_mut().zip(phase.component_moles()) {
                *sum += amount;
            }
        }
        for amount in amounts {
            near(amount, 1.0);
        }
        let integral = joined
            .specific_entropy_change_j_kg_k(0.0, h)
            .unwrap()
            .unwrap();
        assert_eq!(integral, h / transition);
        let original = joined
            .specific_source_entropy_change_j_kg_k(0.0, h, || false)
            .unwrap();
        assert!((original - fraction * GAS_CONSTANT * 2.0_f64.ln() / mass).abs() < 1e-8);
        let expected_residual =
            fraction * joined.coexistence().latent_clausius_residual_j_per_k() / mass;
        assert!((original - integral - expected_residual).abs() < 1e-13);
    }
}

#[test]
fn g1_sensible_branches_keep_source_knots_and_conserved_heat_across_the_latent_join() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let joined = joined(&caloric);
    let mass = joined.mass_kg();
    let transition = transition_temperature();
    for t in [
        380.0,
        400.0,
        450.0,
        540.0,
        transition - 0.1,
        transition + 0.1,
        590.0,
        595.0,
    ] {
        let h = joined
            .specific_enthalpy_j_kg(Temperature::new(t), || false)
            .unwrap();
        let expected = (source_enthalpy(t) - source_enthalpy(transition)) / mass;
        assert!(
            (h - expected).abs() < 3e-6,
            "T={t}, h={h}, expected={expected}"
        );
        let state = joined
            .caloric_state_at_specific_enthalpy(h, || false)
            .unwrap();
        assert!((state.temperature().value() - t).abs() < 2e-7);
        assert_eq!(
            state.branch_state().unwrap().candidate().phases().len(),
            if t < transition { 3 } else { 1 }
        );
        assert!(state.temperature_slope_kg_k_j() > 0.0);
        assert!(state.temperature_slope_kg_k_j() <= joined.maximum_temperature_slope());
        let source = joined
            .specific_source_entropy_change_j_kg_k(0.0, h, || false)
            .unwrap();
        let expected = (source_entropy(t) - source_entropy(transition)) / mass;
        assert!(
            (source - expected).abs() < 1e-8,
            "T={t}, S={source}, expected={expected}"
        );
    }
    let cold = joined
        .specific_enthalpy_j_kg(Temperature::new(400.0), || false)
        .unwrap();
    let hot = joined
        .specific_enthalpy_j_kg(Temperature::new(590.0), || false)
        .unwrap();
    let forward = joined
        .specific_entropy_change_j_kg_k(cold, hot)
        .unwrap()
        .unwrap();
    let reverse = joined
        .specific_entropy_change_j_kg_k(hot, cold)
        .unwrap()
        .unwrap();
    assert_eq!(forward, -reverse);
    assert!(forward > 0.0);
    assert_eq!(
        joined.specific_entropy_change_j_kg_k(cold, cold).unwrap(),
        Some(0.0)
    );
}

#[test]
fn g1_coexistence_join_uses_source_increments_for_stream_datum_branches() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let [closed_low, closed_high] = branches(&caloric, transition_temperature());
    let closed = closed_low
        .join_coexistence(&closed_high, Default::default(), || false)
        .unwrap();
    let total_mass: f64 = MASSES.into_iter().sum();
    let first = MASSES[0] / total_mass;
    let second = MASSES[1] / total_mass;
    let fractions = [first, second, 1.0 - (first + second)];
    let prepare_stream = |initial_temperature, split, temperature_bracket| {
        let mut c = controls();
        c.thermal.temperature_bracket = temperature_bracket;
        caloric
            .equilibrium_stream_storage_from_mass_fractions(
                &endpoint(&caloric, initial_temperature, split),
                &fractions,
                c,
                || false,
            )
            .unwrap()
    };
    let stream_low = prepare_stream(400.0, true, closed_low.temperature_domain());
    let stream_high = prepare_stream(590.0, false, closed_high.temperature_domain());
    for branch in [&stream_low, &stream_high] {
        assert!(branch.specific_enthalpy_domain_j_kg()[0] > 0.0);
        assert!(branch.advected_composition().is_some());
        assert!(branch.advected_entropy_anchor().is_some());
    }
    for (low, high) in [
        (&stream_low, &stream_high),
        (&stream_low, &closed_high),
        (&closed_low, &stream_high),
    ] {
        let joined = low
            .join_coexistence(high, Default::default(), || false)
            .unwrap();
        assert_eq!(
            joined.specific_enthalpy_domain_j_kg(),
            closed.specific_enthalpy_domain_j_kg()
        );
        assert_eq!(
            joined.coexistence_specific_enthalpy_domain_j_kg(),
            closed.coexistence_specific_enthalpy_domain_j_kg()
        );
        assert!(joined.caloric_identity().is_none());
        assert!(joined.advected_composition().is_none());
        assert!(joined.advected_entropy_anchor().is_none());
        let [cold, hot] = joined.specific_enthalpy_domain_j_kg();
        let expected_cold =
            (source_enthalpy(380.0) - source_enthalpy(transition_temperature())) / joined.mass_kg();
        assert!((cold - expected_cold).abs() < 3e-6);
        assert_eq!(joined.temperature_k(cold).unwrap(), 380.0);
        assert_eq!(joined.temperature_k(hot).unwrap(), 595.0);
        for t in [400.0, 450.0, 590.0] {
            let h = joined
                .specific_enthalpy_j_kg(Temperature::new(t), || false)
                .unwrap();
            assert_eq!(
                h,
                closed
                    .specific_enthalpy_j_kg(Temperature::new(t), || false)
                    .unwrap()
            );
            assert!((joined.temperature_k(h).unwrap() - t).abs() < 2e-7);
        }
        assert_eq!(
            joined
                .specific_source_entropy_change_j_kg_k(cold, hot, || false)
                .unwrap(),
            closed
                .specific_source_entropy_change_j_kg_k(cold, hot, || false)
                .unwrap()
        );
    }
}

#[test]
fn g3_metastable_noncoexistence_and_foreign_source_joins_refuse() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let [low, high] = branches(&caloric, transition_temperature());
    let [wrong_low, wrong_high] = branches(&caloric, transition_temperature() - 2.0);
    assert!(
        wrong_low
            .join_coexistence(&wrong_high, Default::default(), || false)
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
    let [_, foreign_high] = branches(&foreign, transition_temperature());
    assert!(matches!(
        low.join_coexistence(&foreign_high, Default::default(), || false),
        Err(MulticomponentLiquidCaloricError::Invalid(_))
    ));
    let controls = MulticomponentCaloricCoexistenceControlV1 {
        maximum_retained_phases: 3,
        ..Default::default()
    };
    assert!(low.join_coexistence(&high, controls, || false).is_err());
}

#[test]
fn g4_join_cancellation_and_unrepresentable_queries_keep_original_branches() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let [low, high] = branches(&caloric, transition_temperature());
    let retained = [low.clone(), high.clone()];
    let mut complete_polls = 0;
    let joined = low
        .join_coexistence(&high, Default::default(), || {
            complete_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert_eq!(
        low.join_coexistence(&high, Default::default(), || {
            polls += 1;
            polls == complete_polls
        }),
        Err(MulticomponentLiquidCaloricError::Cancelled)
    );
    assert_eq!(polls, complete_polls);
    let [a, b] = joined.specific_enthalpy_domain_j_kg();
    for h in [a.next_down(), b.next_up(), f64::NAN] {
        assert!(joined.temperature_k(h).is_err());
    }
    assert!(joined.mass_kg() < 0.5);
    assert!(joined.temperature_k(f64::from_bits(1)).is_err());
    assert_eq!([low, high], retained);
}
