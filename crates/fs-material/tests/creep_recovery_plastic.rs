//! Coupled primary resistance recovery/J2: independent branch solutions,
//! permanent histories, consistent tangents and the original work/heat split.
use fs_material::creep::plastic::{NortonJ2, NortonJ2State};
use fs_material::creep::recovery::plastic::{RecoveringNortonJ2, RecoveringNortonJ2State};
use fs_material::creep::recovery::{
    CreepRecoveryParameters, RecoveringCreepState, RecoveringNortonCreep,
};
use fs_material::creep::thermal::expansion::CreepHistoryStrain;
use fs_material::creep::{
    CreepLaw, CreepState, CreepTolerances, NortonCreep, NortonCreepParameters, StrainHardeningCreep,
};
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::tensor::von_mises;
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt};

const E: f64 = 1.2e9;
const STRESS: f64 = 3e6;
const RATE: f64 = 2e-3;
const A: f64 = 8e-4;
const YIELD: f64 = 2e6;

fn elastic() -> IsotropicElastic {
    IsotropicElastic::new(E, 0.2, 0.05).unwrap()
}

fn creep(n: f64, dt: f64) -> NortonCreep {
    NortonCreep::new(
        elastic(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: n,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 500.0,
            temperature_range_k: [300.0, 900.0],
        },
        500.0,
        dt,
    )
    .unwrap()
    .with_strain_hardening(StrainHardeningCreep {
        reference_strain: A,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    })
    .unwrap()
}

fn native(n: f64, dt: f64, hardening: f64, yield_pa: f64) -> NortonJ2 {
    NortonJ2::new(
        creep(n, dt),
        J2Plasticity::new(elastic(), yield_pa, hardening).unwrap(),
    )
    .unwrap()
}

fn recovery(rate: f64) -> CreepRecoveryParameters {
    CreepRecoveryParameters {
        reference_rate_s: rate,
        activation_energy_j_mol: 0.0,
        reference_temperature_k: 500.0,
    }
}

fn law(n: f64, dt: f64, hardening: f64, rate: f64) -> RecoveringNortonJ2 {
    RecoveringNortonJ2::new(native(n, dt, hardening, YIELD), recovery(rate)).unwrap()
}

fn tolerances() -> CreepTolerances {
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

// Positive root of x*x + b*x - c = 0 without cancellation.
fn quadratic(b: f64, c: f64) -> f64 {
    2.0 * c / (b + (b * b + 4.0 * c).sqrt())
}

#[test]
fn zero_recovery_preserves_original_joint_response_history_and_heat_exactly() {
    for hardening in [0.0, 1e8] {
        let native = native(3.2, 0.25, hardening, YIELD);
        let law = RecoveringNortonJ2::new(native.clone(), recovery(0.0)).unwrap();
        assert_eq!(law.norton_j2(), &native);
        assert_eq!(law.kinetics(), native.creep());
        assert_eq!(law.recovery_parameters(), &recovery(0.0));
        assert_eq!(law.recovery_rate_s(), 0.0);
        assert_eq!(RecoveringNortonJ2State::STRAIN_FIELDS, 2);
        let mut state = law.initial_state();
        let mut old_strain = [0.0; 6];
        for step in 0..4 {
            let strain = match step {
                0 | 1 => dev(0.006),
                2 => state.inelastic_strain(),
                _ => dev(-0.006),
            };
            let expected = native
                .advance_with_energy(&old_strain, &state.native, &strain, tolerances())
                .unwrap();
            let actual = law
                .advance_with_energy(&old_strain, &state, &strain, tolerances())
                .unwrap();
            assert_eq!(actual.state.native, expected.state);
            assert_eq!(
                actual.state.resistance_strain,
                expected.state.creep.equivalent_creep_strain
            );
            assert_eq!(actual.stress_pa, expected.stress_pa);
            assert_eq!(actual.tangent_pa, expected.tangent_pa);
            assert_eq!(actual.energy, expected.energy);
            assert_eq!(
                actual.state.inelastic_strain(),
                expected.state.inelastic_strain()
            );
            state = actual.state;
            old_strain = strain;
        }
    }
}

#[test]
fn recovering_joint_return_matches_quadratic_for_perfect_and_hardening_plasticity() {
    let dt = 0.25;
    let recovery_rate = 1.3;
    let divisor = 1.0 + recovery_rate * dt;
    let effective_a = A * divisor;
    let coefficient = dt * RATE * effective_a / STRESS;
    let modulus = 3.0 * elastic().lame().1;
    let q_trial = modulus * 0.006;
    for hardening in [0.0, 1e8] {
        let law = law(1.0, dt, hardening, recovery_rate);
        let old = RecoveringNortonJ2State {
            native: NortonJ2State {
                creep: CreepState {
                    creep_strain: dev(2e-4),
                    equivalent_creep_strain: 9e-4,
                },
                plastic: J2State {
                    plastic_strain: dev(6e-4),
                    alpha: 6e-4,
                },
            },
            resistance_strain: 6e-4,
        };
        let permanent = old.inelastic_strain();
        let strain = core::array::from_fn(|i| permanent[i] + dev(0.006)[i]);
        let old_yield = YIELD + hardening * old.native.plastic.alpha;
        let weight = hardening / (modulus + hardening);
        let effective_trial = old_yield + weight * (q_trial - old_yield);
        let effective_modulus = modulus * weight;
        let dc = quadratic(
            effective_a + old.resistance_strain + coefficient * effective_modulus,
            coefficient * effective_trial,
        );
        let q = effective_trial - effective_modulus * dc;
        let dp = (q_trial - old_yield - modulus * dc) / (modulus + hardening);
        assert!(dc > 0.0 && dp > 0.0);
        let next = law
            .creep_increment(&permanent, &old, &strain, tolerances())
            .unwrap();
        close(von_mises(&next.stress_pa), q, 5e-8);
        close(
            next.state.native.creep.equivalent_creep_strain,
            9e-4 + dc,
            3e-17,
        );
        close(next.state.native.plastic.alpha, 6e-4 + dp, 3e-17);
        close(next.state.resistance_strain, (6e-4 + dc) / divisor, 3e-17);
        for i in 0..6 {
            close(
                next.state.native.creep.creep_strain[i],
                old.native.creep.creep_strain[i] + dev(dc)[i],
                3e-17,
            );
            close(
                next.state.native.plastic.plastic_strain[i],
                old.native.plastic.plastic_strain[i] + dev(dp)[i],
                3e-17,
            );
        }
        let stored = q * q / (2.0 * modulus) + 0.5 * hardening * (6e-4 + dp).powi(2);
        close(
            next.energy.initial_stored_j_m3,
            0.5 * hardening * (6e-4_f64).powi(2),
            2e-8,
        );
        close(next.energy.final_stored_j_m3, stored, 2e-8);
        close(
            next.energy.numerical_dissipation_j_m3,
            q * q / (2.0 * modulus) + 0.5 * hardening * dp * dp,
            2e-8,
        );
        close(next.energy.creep_dissipation_j_m3, q * dc, 2e-8);
        close(next.energy.plastic_dissipation_j_m3, YIELD * dp, 2e-8);
        close(
            next.energy.physical_dissipation_j_m3,
            q * dc + YIELD * dp,
            2e-8,
        );
        close(next.energy.mechanical_work_j_m3, q * 0.006, 3e-8);
        close(next.energy.residual_j_m3, 0.0, tolerances().energy_j_m3);
        let radial = weight
            / (1.0
                + effective_modulus * coefficient
                    / (effective_a + old.resistance_strain + 2.0 * dc));
        let action: f64 = (0..6).map(|j| next.tangent_pa[0][j] * dev(1.0)[j]).sum();
        close(action, 2.0 * elastic().lame().1 * radial, 3e-5);
        assert_eq!(
            law.kinetics().strain_hardening().unwrap().reference_strain,
            A
        );
    }
}

#[test]
fn yield_crossing_uses_recovered_resistance_and_endpoint_scale_not_permanent_alpha() {
    let baseline = law(1.0, 0.25, 1e8, 0.0);
    let old = RecoveringNortonJ2State {
        native: NortonJ2State {
            creep: CreepState {
                creep_strain: [0.0; 6],
                equivalent_creep_strain: 3e-3,
            },
            plastic: baseline.norton_j2().plastic().initial_state(),
        },
        resistance_strain: 3e-3,
    };
    let strain = dev(2.2e6 / (3.0 * elastic().lame().1));
    let unchanged = baseline
        .advance_with_energy(&[0.0; 6], &old, &strain, tolerances())
        .unwrap();
    assert!(unchanged.state.native.plastic.alpha > 0.0);
    // Either recovery inside this interval, or a retained recovered resistance,
    // allows creep to relax the above-yield predictor without plastic flow.
    for (rate, resistance_strain) in [(40.0, 3e-3), (0.0, 1e-4)] {
        let law = law(1.0, 0.25, 1e8, rate);
        let state = RecoveringNortonJ2State {
            resistance_strain,
            ..old.clone()
        };
        let next = law
            .advance_with_energy(&[0.0; 6], &state, &strain, tolerances())
            .unwrap();
        assert_eq!(next.state.native.plastic, state.native.plastic);
        assert!(von_mises(&next.stress_pa) < YIELD);
        let pure = RecoveringNortonCreep::new(law.kinetics().clone(), recovery(rate)).unwrap();
        let expected = pure
            .advance_with_energy(
                &[0.0; 6],
                &RecoveringCreepState {
                    native: state.native.creep.clone(),
                    resistance_strain,
                },
                &strain,
                tolerances(),
            )
            .unwrap();
        assert_eq!(next.state.native.creep, expected.state.native);
        assert_eq!(
            next.state.resistance_strain,
            expected.state.resistance_strain
        );
        assert_eq!(next.stress_pa, expected.stress_pa);
        assert_eq!(next.tangent_pa, expected.tangent_pa);
        close(
            next.energy.creep_dissipation_j_m3,
            expected.energy.creep_dissipation_j_m3,
            1e-8,
        );
        assert_eq!(next.energy.plastic_dissipation_j_m3, 0.0);
    }
}

#[test]
fn nonproportional_joint_recovery_has_consistent_tangent_and_physical_heat_split() {
    for hardening in [0.0, 2e8] {
        let creep = creep(3.2, 0.125)
            .with_strain_hardening(StrainHardeningCreep {
                reference_strain: A,
                exponent: 0.8,
                residual_rate_fraction: 0.06,
            })
            .unwrap();
        let law = RecoveringNortonJ2::new(
            NortonJ2::new(
                creep,
                J2Plasticity::new(elastic(), YIELD, hardening).unwrap(),
            )
            .unwrap(),
            recovery(1.7),
        )
        .unwrap();
        let first_strain = dev(0.004);
        let old = law
            .advance_with_energy(&[0.0; 6], &law.initial_state(), &first_strain, tolerances())
            .unwrap()
            .state;
        let frozen = old.clone();
        let direction = [0.005, -0.001, -0.002, 0.0003, -0.0005, 0.0002];
        let permanent = old.inelastic_strain();
        let strain = core::array::from_fn(|i| permanent[i] + direction[i]);
        let next = law
            .advance_with_energy(&first_strain, &old, &strain, tolerances())
            .unwrap();
        let dc = next.state.native.creep.equivalent_creep_strain
            - old.native.creep.equivalent_creep_strain;
        let dp = next.state.native.plastic.alpha - old.native.plastic.alpha;
        assert!(dc > 0.0 && dp > 0.0);
        let q = von_mises(&next.stress_pa);
        close(
            dc,
            0.125
                * law
                    .kinetics()
                    .equivalent_creep_rate_at_strain_s(q, next.state.resistance_strain)
                    .unwrap(),
            3e-17,
        );
        close(q, YIELD + hardening * next.state.native.plastic.alpha, 1e-7);
        close(next.energy.creep_dissipation_j_m3, q * dc, 2e-8);
        close(next.energy.plastic_dissipation_j_m3, YIELD * dp, 2e-8);
        close(next.energy.residual_j_m3, 0.0, tolerances().energy_j_m3);
        let step = 2e-9;
        for column in 0..6 {
            let mut plus = strain;
            let mut minus = strain;
            plus[column] += step;
            minus[column] -= step;
            let p = law.response(&plus, &old).unwrap().stress_pa;
            let m = law.response(&minus, &old).unwrap().stress_pa;
            for row in 0..6 {
                close(
                    (p[row] - m[row]) / (2.0 * step),
                    next.tangent_pa[row][column],
                    2e-6 * next.tangent_pa[row][column].abs().max(1e6),
                );
                let wi = if row < 3 { 1.0 } else { 2.0 };
                let wj = if column < 3 { 1.0 } else { 2.0 };
                close(
                    wi * next.tangent_pa[row][column],
                    wj * next.tangent_pa[column][row],
                    2e-6,
                );
            }
        }
        assert_eq!(old, frozen);
        assert_eq!(law.response(&strain, &old).unwrap().state, next.state);
        assert_eq!(law.update_state(&strain, &old), next.state);
    }
}

#[test]
fn unloaded_dwell_recovers_only_resistance_and_redistributes_reload_flow() {
    let old = RecoveringNortonJ2State {
        native: NortonJ2State {
            creep: CreepState {
                creep_strain: dev(1.0 / 4096.0),
                equivalent_creep_strain: 1e-3,
            },
            plastic: J2State {
                plastic_strain: dev(1.0 / 2048.0),
                alpha: 1.0 / 2048.0,
            },
        },
        resistance_strain: 8e-4,
    };
    let strain = old.inelastic_strain();
    let dwell = law(1.0, 2.0, 1e8, 1.5)
        .advance_with_energy(&strain, &old, &strain, tolerances())
        .unwrap();
    assert_eq!(dwell.state.native, old.native);
    assert_eq!(dwell.state.resistance_strain, old.resistance_strain / 4.0);
    assert_eq!(dwell.stress_pa, [0.0; 6]);
    assert_eq!(dwell.energy.creep_dissipation_j_m3, 0.0);
    assert_eq!(dwell.energy.plastic_dissipation_j_m3, 0.0);
    assert_eq!(dwell.energy.physical_dissipation_j_m3, 0.0);
    assert_eq!(
        dwell.energy.initial_stored_j_m3,
        dwell.energy.final_stored_j_m3
    );
    assert_eq!(dwell.energy.numerical_dissipation_j_m3, 0.0);
    let reload = law(1.0, 0.25, 1e8, 0.0);
    let end = core::array::from_fn(|i| strain[i] + dev(0.006)[i]);
    let before = reload
        .advance_with_energy(&strain, &old, &end, tolerances())
        .unwrap();
    let after = reload
        .advance_with_energy(&strain, &dwell.state, &end, tolerances())
        .unwrap();
    assert!(
        after.state.native.creep.equivalent_creep_strain
            > before.state.native.creep.equivalent_creep_strain
    );
    assert!(after.state.native.plastic.alpha < before.state.native.plastic.alpha);
    assert!(after.state.native.plastic.alpha > old.native.plastic.alpha);
    assert!(von_mises(&after.stress_pa) < von_mises(&before.stress_pa));
    assert!(after.energy.creep_dissipation_j_m3 > 0.0);
    assert!(after.energy.plastic_dissipation_j_m3 > 0.0);
}

#[test]
fn actual_coupled_creep_increment_survives_accumulated_path_roundoff() {
    let dt = 1e-17;
    let recovery_rate = 1e17;
    for hardening in [0.0, 1e8] {
        let law = law(1.0, dt, hardening, recovery_rate);
        let old = RecoveringNortonJ2State {
            native: NortonJ2State {
                creep: CreepState {
                    creep_strain: [0.0; 6],
                    equivalent_creep_strain: 1.0,
                },
                plastic: law.norton_j2().plastic().initial_state(),
            },
            resistance_strain: 0.0,
        };
        let next = law
            .advance_with_energy(&[0.0; 6], &old, &dev(0.006), tolerances())
            .unwrap();
        assert_eq!(next.state.native.creep.equivalent_creep_strain, 1.0);
        assert!(next.state.resistance_strain > 0.0);
        assert!(next.energy.creep_dissipation_j_m3 > 0.0);
        assert!(next.energy.plastic_dissipation_j_m3 > 0.0);
        let dc = next.state.resistance_strain * (1.0 + recovery_rate * dt);
        let q = von_mises(&next.stress_pa);
        let rate = law
            .kinetics()
            .equivalent_creep_rate_at_strain_s(q, next.state.resistance_strain)
            .unwrap();
        close(dc, dt * rate, 1e-31);
        close(next.energy.creep_dissipation_j_m3, q * dc, 1e-25);
        close(
            next.energy.plastic_dissipation_j_m3,
            YIELD * next.state.native.plastic.alpha,
            2e-8,
        );
    }
}

#[test]
fn changed_yield_at_held_strain_keeps_native_storage_and_both_heat_channels() {
    let strong = RecoveringNortonJ2::new(native(1.0, 0.25, 1e8, 20e6), recovery(0.2)).unwrap();
    let softened = law(1.0, 0.25, 1e8, 1.5);
    let strain = dev(0.006);
    let first = strong
        .advance_with_energy(&[0.0; 6], &strong.initial_state(), &strain, tolerances())
        .unwrap();
    assert_eq!(first.state.native.plastic.alpha, 0.0);
    let next = softened
        .advance_with_energy(&strain, &first.state, &strain, tolerances())
        .unwrap();
    assert!(next.state.native.plastic.alpha > 0.0);
    assert!(
        next.state.native.creep.equivalent_creep_strain
            > first.state.native.creep.equivalent_creep_strain
    );
    assert!(next.energy.creep_dissipation_j_m3 > 0.0);
    assert!(next.energy.plastic_dissipation_j_m3 > 0.0);
    assert_eq!(next.energy.mechanical_work_j_m3, 0.0);
    close(
        next.energy.initial_stored_j_m3,
        first.energy.final_stored_j_m3,
        2e-8,
    );
    close(
        next.energy.initial_stored_j_m3 - next.energy.final_stored_j_m3,
        next.energy.physical_dissipation_j_m3 + next.energy.numerical_dissipation_j_m3,
        2e-7,
    );
}

#[test]
fn resistance_is_bounded_by_creep_alone_and_both_native_histories_remain_checked() {
    let law = law(2.0, 0.1, 1e8, 1.0);
    let state = RecoveringNortonJ2State {
        native: NortonJ2State {
            creep: CreepState {
                creep_strain: [0.0; 6],
                equivalent_creep_strain: 1e-3,
            },
            plastic: J2State {
                plastic_strain: [0.0; 6],
                alpha: 0.01,
            },
        },
        resistance_strain: 2e-3,
    };
    let frozen = state.clone();
    assert!(law.validate_history(&[0.0; 6], &state).is_err());
    assert!(
        law.advance_with_energy(&[0.0; 6], &state, &[0.0; 6], tolerances())
            .is_err()
    );
    assert!(law.stress(&[0.0; 6], &state).iter().all(|x| x.is_nan()));
    assert_eq!(state, frozen);
    for invalid_creep in [true, false] {
        let mut bad = law.initial_state();
        if invalid_creep {
            bad.native.creep.creep_strain = dev(1e-3);
        } else {
            bad.native.plastic.plastic_strain = dev(1e-3);
        }
        assert!(law.response(&[0.0; 6], &bad).is_err());
    }
    let original = law.kinetics();
    let secondary = NortonCreep::new(
        elastic(),
        original.parameters().clone(),
        original.temperature_k(),
        original.duration_s(),
    )
    .unwrap();
    assert!(
        RecoveringNortonJ2::new(
            NortonJ2::new(secondary, law.norton_j2().plastic().clone()).unwrap(),
            recovery(1.0),
        )
        .is_err()
    );
    assert!(RecoveringNortonJ2::new(native(2.0, 0.1, 1e8, YIELD), recovery(-1.0)).is_err());
}
