//! Independent ternary common-tangent and actual finite coexistence inventories.

use super::*;
use crate::solution::multicomponent::caloric::equilibrium::coexistence::MulticomponentCaloricCoexistenceControlV1;

fn pair(
    model: &MulticomponentLiquidCaloricModelV1,
    temperature: f64,
) -> [MulticomponentEquilibriumCaloricStateV1; 2] {
    let low = state(model, temperature, &[1.0; 3], &seed(3, 2.0 / 3.0));
    let high = state(
        model,
        temperature,
        &[1.0; 3],
        &MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![vec![1.0 / 3.0; 3]],
            phase_weights: vec![1.0],
        },
    );
    [low, high]
}

#[test]
fn g1_true_ternary_latent_heat_retains_original_entropy_and_phase_inventories() {
    let temperature = 1000.0 / (4.0 * 2.0_f64.ln() - 0.5);
    let expected_heat = R * 1000.0 / 4.0;
    let expected_entropy = R * (2.0_f64.ln() - 0.125);
    for datum in [100.0, 1e20] {
        let model = model(3, 0.5, 1000.0, datum, false);
        let [low, high] = pair(&model, temperature);
        if datum == 1e20 {
            assert_eq!(low.enthalpy_j(), high.enthalpy_j());
        }
        let plateau = model
            .coexistence_plateau(&low, &high, Default::default(), || false)
            .unwrap();
        assert_eq!(plateau.lower(), &low);
        assert_eq!(plateau.upper(), &high);
        near(plateau.latent_enthalpy_j(), expected_heat, 1e-8, 0.0);
        near(
            plateau.latent_entropy_j_per_k(),
            expected_entropy,
            1e-10,
            0.0,
        );
        assert!(plateau.latent_clausius_residual_j_per_k().abs() < 1e-10);
        assert!(plateau.maximum_chemical_potential_residual() < 1e-11);
        for fraction in [0.0, 0.2, 0.5, 0.8, 1.0] {
            let heat = fraction * plateau.latent_enthalpy_j();
            let mixed = plateau.state_at_enthalpy_increment(heat, || false).unwrap();
            assert_eq!(mixed.temperature().value(), temperature);
            near(
                mixed.upper_ensemble_fraction(),
                fraction,
                2.0 * f64::EPSILON,
                0.0,
            );
            near(
                mixed.enthalpy_increment_j(),
                fraction * expected_heat,
                1e-8,
                0.0,
            );
            near(
                mixed.entropy_increment_j_per_k(),
                fraction * expected_entropy,
                1e-10,
                0.0,
            );
            assert!(mixed.enthalpy_residual_j().abs() < 1e-12);
            assert_eq!(
                mixed.phases().count(),
                if fraction == 0.0 {
                    3
                } else if fraction == 1.0 {
                    1
                } else {
                    4
                }
            );
            let mut amounts = [0.0; 3];
            for phase in mixed.phases() {
                assert!(phase.total_moles() > 0.0);
                for (total, n) in amounts.iter_mut().zip(phase.component_moles()) {
                    *total += n;
                }
            }
            for amount in amounts {
                near(amount, 1.0, 64.0 * f64::EPSILON, 0.0);
            }
        }
    }
}

#[test]
fn g3_noncoexisting_metastable_branches_foreign_sources_and_reversed_heat_refuse() {
    let model = model(3, 0.5, 1000.0, 100.0, false);
    let temperature = 1000.0 / (4.0 * 2.0_f64.ln() - 0.5);
    let [low, high] = pair(&model, temperature);
    let original = [low.clone(), high.clone()];
    let [off_low, off_high] = pair(&model, temperature - 2.0);
    assert!(
        model
            .coexistence_plateau(&off_low, &off_high, Default::default(), || false)
            .is_err()
    );
    assert!(
        model
            .coexistence_plateau(&high, &low, Default::default(), || false)
            .is_err()
    );
    assert!(
        model
            .coexistence_plateau(&low, &low, Default::default(), || false)
            .is_err()
    );
    let foreign = super::model(3, 0.5, 1000.0, 120.0, false);
    assert!(
        foreign
            .coexistence_plateau(&low, &high, Default::default(), || false)
            .is_err()
    );
    let mut controls = MulticomponentCaloricCoexistenceControlV1::default();
    controls.maximum_retained_phases = 3;
    assert!(
        model
            .coexistence_plateau(&low, &high, controls, || false)
            .is_err()
    );
    assert_eq!([low, high], original);
}

#[test]
fn g4_coexistence_final_poll_and_latent_domain_refuse_without_mutation() {
    let model = model(3, 0.5, 1000.0, 100.0, false);
    let temperature = 1000.0 / (4.0 * 2.0_f64.ln() - 0.5);
    let [low, high] = pair(&model, temperature);
    let mut complete_polls = 0;
    let plateau = model
        .coexistence_plateau(&low, &high, Default::default(), || {
            complete_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert_eq!(
        model.coexistence_plateau(&low, &high, Default::default(), || {
            polls += 1;
            polls == complete_polls
        }),
        Err(MulticomponentEquilibriumCaloricErrorV1::Cancelled)
    );
    assert_eq!(polls, complete_polls);
    let retained = plateau.clone();
    for heat in [-1.0, plateau.latent_enthalpy_j().next_up(), f64::NAN] {
        assert!(plateau.state_at_enthalpy_increment(heat, || false).is_err());
    }
    assert_eq!(plateau, retained);
}
