//! Native primary resistance recovery: analytic return, permanent history,
//! discrete tangent, independently shifted temperature kinetics and work.
use fs_material::creep::recovery::{
    CreepRecoveryParameters, RecoveringCreepState, RecoveringNortonCreep,
};
use fs_material::creep::thermal::expansion::CreepHistoryStrain;
use fs_material::creep::{
    CreepLaw, CreepState, CreepTolerances, GAS_CONSTANT_J_MOL_K, NortonCreep,
    NortonCreepParameters, StrainHardeningCreep,
};
use fs_material::tensor::von_mises;
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt};

const E: f64 = 1.2e9;
const STRESS: f64 = 3e6;
const RATE: f64 = 2e-3;
const A: f64 = 8e-4;

fn norton(n: f64, dt: f64, temperature: f64) -> NortonCreep {
    NortonCreep::new(
        IsotropicElastic::new(E, 0.2, 0.05).unwrap(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: n,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 500.0,
            temperature_range_k: [300.0, 900.0],
        },
        temperature,
        dt,
    )
    .unwrap()
}

fn primary(n: f64, dt: f64, temperature: f64) -> NortonCreep {
    norton(n, dt, temperature)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: A,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap()
}

fn parameters(rate: f64) -> CreepRecoveryParameters {
    CreepRecoveryParameters {
        reference_rate_s: rate,
        activation_energy_j_mol: 0.0,
        reference_temperature_k: 500.0,
    }
}

fn recovering(n: f64, dt: f64, rate: f64) -> RecoveringNortonCreep {
    RecoveringNortonCreep::new(primary(n, dt, 500.0), parameters(rate)).unwrap()
}

fn tolerance() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 3e-6,
        energy_j_m3: 2e-7,
    }
}

fn dev(a: f64) -> Voigt {
    [a, -0.5 * a, -0.5 * a, 0.0, 0.0, 0.0]
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b} (tolerance {tolerance})"
    );
}

#[test]
fn zero_recovery_reduces_exactly_to_original_primary_history_tangent_and_energy() {
    let native = primary(3.2, 0.4, 500.0);
    let law = RecoveringNortonCreep::new(native.clone(), parameters(0.0)).unwrap();
    assert_eq!(law.norton(), &native);
    assert_eq!(law.kinetics(), &native);
    assert_eq!(law.recovery_parameters(), &parameters(0.0));
    assert_eq!(law.recovery_rate_s(), 0.0);
    let mut state = law.initial_state();
    let mut old_strain = [0.0; 6];
    for strain in [
        [0.002, -0.0005, 0.0001, 0.001, -0.0002, 0.0006],
        [-0.001, 0.0018, -0.0004, 0.0002, 0.0012, -0.0008],
        [0.001, -0.0001, -0.0007, -0.001, 0.0002, 0.0003],
    ] {
        let expected = native
            .advance_with_energy(&old_strain, &state.native, &strain, tolerance())
            .unwrap();
        let actual = law
            .advance_with_energy(&old_strain, &state, &strain, tolerance())
            .unwrap();
        assert_eq!(actual.state.native, expected.state);
        assert_eq!(
            actual.state.resistance_strain,
            expected.state.equivalent_creep_strain
        );
        assert_eq!(actual.stress_pa, expected.stress_pa);
        assert_eq!(actual.tangent_pa, expected.tangent_pa);
        assert_eq!(actual.energy, expected.energy);
        assert_eq!(actual.state.inelastic_strain(), expected.state.creep_strain);
        state = actual.state;
        old_strain = strain;
    }
}

#[test]
fn recovering_primary_return_matches_independent_quadratic_and_native_heat() {
    let dt = 0.7;
    let recovery_rate = 1.3;
    let law = recovering(1.0, dt, recovery_rate);
    let old = RecoveringCreepState {
        native: CreepState {
            creep_strain: dev(2e-4),
            equivalent_creep_strain: 9e-4,
        },
        resistance_strain: 6e-4,
    };
    let strain = core::array::from_fn(|i| old.native.creep_strain[i] + dev(0.004)[i]);
    let modulus = 3.0 * law.norton().elastic().lame().1;
    let q_trial = modulus * 0.004;
    let divisor = 1.0 + recovery_rate * dt;
    let effective_a = A * divisor;
    let coefficient = dt * RATE * effective_a / STRESS;
    let b = effective_a + old.resistance_strain + coefficient * modulus;
    let c = coefficient * q_trial;
    let increment = 2.0 * c / (b + (b * b + 4.0 * c).sqrt());
    let q = q_trial - modulus * increment;
    let actual = law
        .advance_with_energy(&old.native.creep_strain, &old, &strain, tolerance())
        .unwrap();
    close(
        actual.state.native.equivalent_creep_strain,
        9e-4 + increment,
        3e-17,
    );
    close(
        actual.state.resistance_strain,
        (6e-4 + increment) / divisor,
        3e-17,
    );
    close(von_mises(&actual.stress_pa), q, 3e-8);
    close(actual.energy.creep_dissipation_j_m3, q * increment, 1e-8);
    assert_eq!(
        actual.energy.physical_dissipation_j_m3,
        actual.energy.creep_dissipation_j_m3
    );
    assert_eq!(actual.energy.plastic_dissipation_j_m3, 0.0);
    close(actual.energy.residual_j_m3, 0.0, tolerance().energy_j_m3);
    let radial = 1.0
        / (1.0 + modulus * coefficient / (effective_a + old.resistance_strain + 2.0 * increment));
    let tangent_action: f64 = (0..6).map(|j| actual.tangent_pa[0][j] * dev(1.0)[j]).sum();
    close(
        tangent_action,
        2.0 * law.norton().elastic().lame().1 * radial,
        3e-5,
    );
    // Material identity always exposes the authored scale, never a* from dt.
    assert_eq!(law.norton().strain_hardening().unwrap().reference_strain, A);
}

#[test]
fn unloaded_recovery_preserves_permanent_strain_and_converges_to_exponential() {
    let initial = RecoveringCreepState {
        native: CreepState {
            creep_strain: dev(2e-4),
            equivalent_creep_strain: 1e-3,
        },
        resistance_strain: 1e-3,
    };
    let time = 4.0;
    let rate = 0.8;
    let exact = initial.resistance_strain * f64::exp(-rate * time);
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let dt = time / f64::from(steps);
        let law = recovering(1.0, dt, rate);
        let mut state = initial.clone();
        for _ in 0..steps {
            let frozen = state.clone();
            let strain = state.native.creep_strain;
            let actual = law
                .creep_increment(&strain, &state, &strain, tolerance())
                .unwrap();
            assert_eq!(state, frozen);
            assert_eq!(actual.state.native, initial.native);
            assert_eq!(actual.stress_pa, [0.0; 6]);
            assert_eq!(actual.energy.physical_dissipation_j_m3, 0.0);
            assert_eq!(actual.energy.mechanical_work_j_m3, 0.0);
            assert_eq!(actual.energy.numerical_dissipation_j_m3, 0.0);
            assert_eq!(actual.energy.residual_j_m3, 0.0);
            assert_eq!(
                actual.state.resistance_strain,
                state.resistance_strain / (1.0 + rate * dt)
            );
            assert!(actual.state.resistance_strain < state.resistance_strain);
            state = actual.state;
        }
        errors.push((state.resistance_strain - exact).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.2,
            "{errors:?}"
        );
    }
}

#[test]
fn dwell_recovery_renews_faster_creep_on_reload_without_resetting_accumulated_path() {
    let initial = RecoveringCreepState {
        native: CreepState {
            creep_strain: dev(1e-3),
            equivalent_creep_strain: 3e-3,
        },
        resistance_strain: 3e-3,
    };
    let strain = initial.native.creep_strain;
    let recovered = recovering(3.0, 2.0, 1.0)
        .advance_with_energy(&strain, &initial, &strain, tolerance())
        .unwrap()
        .state;
    assert_eq!(recovered.native, initial.native);
    let reload = recovering(3.0, 0.2, 0.0);
    let end = core::array::from_fn(|i| strain[i] + dev(0.003)[i]);
    let without_dwell = reload
        .advance_with_energy(&strain, &initial, &end, tolerance())
        .unwrap();
    let with_dwell = reload
        .advance_with_energy(&strain, &recovered, &end, tolerance())
        .unwrap();
    assert!(
        with_dwell.state.native.equivalent_creep_strain
            > without_dwell.state.native.equivalent_creep_strain
    );
    assert!(von_mises(&with_dwell.stress_pa) < von_mises(&without_dwell.stress_pa));
    assert!(
        with_dwell.state.native.equivalent_creep_strain > initial.native.equivalent_creep_strain
    );
    assert!(with_dwell.state.resistance_strain < with_dwell.state.native.equivalent_creep_strain);
}

#[test]
fn nonproportional_recovery_tangent_differentiates_the_coupled_endpoint_and_closes_work() {
    let native = primary(3.2, 0.4, 500.0)
        .with_strain_hardening(StrainHardeningCreep {
            reference_strain: A,
            exponent: 0.8,
            residual_rate_fraction: 0.06,
        })
        .unwrap();
    let law = RecoveringNortonCreep::new(native, parameters(1.7)).unwrap();
    let first_strain = [0.002, -0.0005, 0.0001, 0.001, -0.0002, 0.0006];
    let old = law.update_state(&first_strain, &law.initial_state());
    let frozen = old.clone();
    let strain = [-0.001, 0.0018, -0.0004, 0.0002, 0.0012, -0.0008];
    let response = law
        .advance_with_energy(&first_strain, &old, &strain, tolerance())
        .unwrap();
    for column in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        let step = 1e-8;
        plus[column] += step;
        minus[column] -= step;
        let p = law.stress(&plus, &old);
        let m = law.stress(&minus, &old);
        for row in 0..6 {
            close(
                (p[row] - m[row]) / (2.0 * step),
                response.tangent_pa[row][column],
                0.5,
            );
        }
    }
    assert_eq!(old, frozen);
    assert!(response.energy.creep_dissipation_j_m3 > 0.0);
    close(response.energy.residual_j_m3, 0.0, tolerance().energy_j_m3);
    assert_eq!(law.update_state(&strain, &old), response.state);
}

#[test]
fn solved_increment_survives_roundoff_of_the_permanent_accumulated_path() {
    let dt = 1e-17;
    let rate = 1e17;
    let law = recovering(1.0, dt, rate);
    let old = RecoveringCreepState {
        native: CreepState {
            creep_strain: [0.0; 6],
            equivalent_creep_strain: 1.0,
        },
        resistance_strain: 0.0,
    };
    let strain = dev(1e6 / (3.0 * law.norton().elastic().lame().1));
    let response = law
        .advance_with_energy(&strain, &old, &strain, tolerance())
        .unwrap();
    assert_eq!(response.state.native.equivalent_creep_strain, 1.0);
    assert!(response.state.resistance_strain > 0.0);
    assert!(response.energy.creep_dissipation_j_m3 > 0.0);
    let increment = response.state.resistance_strain * (1.0 + rate * dt);
    let q = von_mises(&response.stress_pa);
    let rate = law
        .norton()
        .equivalent_creep_rate_at_strain_s(q, response.state.resistance_strain)
        .unwrap();
    close(increment, dt * rate, 1e-32);
    close(response.energy.creep_dissipation_j_m3, q * increment, 1e-26);
}

#[test]
fn recovery_has_independent_arrhenius_kinetics_and_zero_is_an_exact_selection() {
    let recovery = CreepRecoveryParameters {
        reference_rate_s: 0.7,
        activation_energy_j_mol: 40_000.0,
        reference_temperature_k: 500.0,
    };
    let law = RecoveringNortonCreep::new(primary(2.0, 0.1, 700.0), recovery).unwrap();
    let expected = 0.7 * f64::exp(40_000.0 / GAS_CONSTANT_J_MOL_K * (1.0 / 500.0 - 1.0 / 700.0));
    close(law.recovery_rate_s(), expected, 3e-13 * expected);
    assert_eq!(law.recovery_parameters(), &recovery);
    assert_eq!(law.norton().parameters().activation_energy_j_mol, 0.0);
    let reference = RecoveringNortonCreep::new(primary(2.0, 0.1, 500.0), recovery).unwrap();
    assert_eq!(reference.recovery_rate_s(), 0.7);
    let zero = RecoveringNortonCreep::new(
        primary(2.0, 0.1, 900.0),
        CreepRecoveryParameters {
            reference_rate_s: 0.0,
            activation_energy_j_mol: f64::MAX,
            reference_temperature_k: 300.0,
        },
    )
    .unwrap();
    assert_eq!(zero.recovery_rate_s(), 0.0);
}

#[test]
fn invalid_primary_selection_parameters_and_resistance_histories_refuse() {
    assert!(RecoveringNortonCreep::new(norton(2.0, 0.1, 500.0), parameters(1.0)).is_err());
    for recovery in [
        parameters(-1.0),
        parameters(f64::NAN),
        parameters(f64::INFINITY),
        CreepRecoveryParameters {
            activation_energy_j_mol: -1.0,
            ..parameters(1.0)
        },
        CreepRecoveryParameters {
            activation_energy_j_mol: f64::INFINITY,
            ..parameters(1.0)
        },
        CreepRecoveryParameters {
            reference_temperature_k: 0.0,
            ..parameters(1.0)
        },
        CreepRecoveryParameters {
            reference_temperature_k: 901.0,
            ..parameters(1.0)
        },
    ] {
        assert!(RecoveringNortonCreep::new(primary(2.0, 0.1, 500.0), recovery).is_err());
    }
    assert!(RecoveringNortonCreep::new(primary(2.0, 1e308, 500.0), parameters(2.0)).is_err());
    for (temperature, reference) in [(900.0, 300.0), (300.0, 900.0)] {
        assert!(
            RecoveringNortonCreep::new(
                primary(2.0, 0.1, temperature),
                CreepRecoveryParameters {
                    activation_energy_j_mol: 1e8,
                    reference_temperature_k: reference,
                    ..parameters(1.0)
                },
            )
            .is_err()
        );
    }
    let law = recovering(2.0, 0.1, 1.0);
    for resistance_strain in [-1e-6, f64::NAN, f64::INFINITY, 1.01e-3] {
        let bad = RecoveringCreepState {
            native: CreepState {
                creep_strain: [0.0; 6],
                equivalent_creep_strain: 1e-3,
            },
            resistance_strain,
        };
        assert!(law.validate_history(&[0.0; 6], &bad).is_err());
        assert!(
            law.advance_with_energy(&[0.0; 6], &bad, &[0.0; 6], tolerance())
                .is_err()
        );
        assert!(law.stress(&[0.0; 6], &bad).iter().all(|x| x.is_nan()));
    }
    let bad_native = RecoveringCreepState {
        native: CreepState {
            creep_strain: dev(2e-3),
            equivalent_creep_strain: 1e-3,
        },
        resistance_strain: 0.0,
    };
    assert!(
        law.validate_history(&bad_native.native.creep_strain, &bad_native)
            .is_err()
    );
    let virgin_resistance = RecoveringCreepState {
        resistance_strain: f64::from_bits(1),
        ..law.initial_state()
    };
    assert!(law.validate_history(&[0.0; 6], &virgin_resistance).is_err());
}
