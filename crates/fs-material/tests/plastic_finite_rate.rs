//! Synthetic finite Perzyna oracles; these are not empirical calibrations.

use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Error, FiniteJ2Limits, FiniteJ2Parameters};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const Y: f64 = 2e6;
const H: f64 = 8e6;
const A: f64 = 3.0 * MU + H;

fn native() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 3.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}
fn kinetics(eta: f64, n: f64, dt: f64) -> FiniteOverstressJ2 {
    FiniteOverstressJ2::new(
        native(),
        FiniteOverstressParameters {
            viscosity_pa_s: eta,
            stress_scale_pa: 2e6,
            exponent: n,
        },
        dt,
    )
    .unwrap()
}
fn coaxial(mean: f64, axial: f64) -> [f64; 9] {
    [
        (mean + axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
        0.0,
        0.0,
        0.0,
        (mean - 0.5 * axial).exp(),
    ]
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual {actual:.16e}, expected {expected:.16e}, error {:.8e}",
        (actual - expected).abs()
    );
}
fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}
fn equivalent(tau: &[f64; 9]) -> f64 {
    let mean = (tau[0] + tau[4] + tau[8]) / 3.0;
    let dev: [f64; 9] =
        core::array::from_fn(|i| tau[i] - if [0, 4, 8].contains(&i) { mean } else { 0.0 });
    (1.5 * dev.iter().map(|x| x * x).sum::<f64>()).sqrt()
}

#[test]
fn linear_return_has_independent_finite_hencky_stress_history_and_distinct_rate_potential() {
    let (eta, dt, mean, axial) = (12e6, 0.25, 0.03, 0.21);
    let law = kinetics(eta, 1.0, dt);
    let initial = law.initial_state();
    let f = coaxial(mean, axial);
    let result = law.response(&f, &initial).unwrap();
    let gamma = (3.0 * MU * axial - Y) / (A + eta / dt);
    let viscous = eta * gamma / dt;
    let elastic = axial - gamma;
    let stored =
        0.5 * K * (3.0 * mean).powi(2) + 1.5 * MU * elastic * elastic + 0.5 * H * gamma * gamma;
    near(result.plastic_increment, gamma, 1e-14, 1e-13);
    near(result.state.alpha, gamma, 1e-14, 1e-13);
    near(
        result.state.inverse_plastic_cauchy_green[0],
        (-2.0 * gamma).exp(),
        1e-14,
        1e-13,
    );
    for i in [1, 2] {
        near(
            result.state.inverse_plastic_cauchy_green[i],
            gamma.exp(),
            1e-14,
            1e-13,
        );
    }
    near(
        result.equivalent_kirchhoff_pa,
        Y + H * gamma + viscous,
        1e-6,
        1e-13,
    );
    near(result.jacobian, (3.0 * mean).exp(), 1e-14, 1e-14);
    for (i, tau) in [
        (0, 3.0 * K * mean + 2.0 * MU * elastic),
        (4, 3.0 * K * mean - MU * elastic),
        (8, 3.0 * K * mean - MU * elastic),
    ] {
        near(result.kirchhoff_pa[i], tau, 1e-6, 1e-13);
        near(result.first_piola_pa[i], tau / f[i], 1e-6, 1e-13);
    }
    near(result.stored_energy_j_m3, stored, 1e-7, 1e-13);
    near(
        result.plastic_dissipation_j_m3,
        (Y + viscous) * gamma,
        1e-7,
        1e-13,
    );
    near(
        result.incremental_potential_j_m3,
        stored + Y * gamma + 0.5 * viscous * gamma,
        1e-7,
        1e-13,
    );
    assert!(
        result.stored_energy_j_m3 + result.plastic_dissipation_j_m3
            - result.incremental_potential_j_m3
            > 1e5
    );
    assert!(
        native()
            .validate_committed(&f, &result.state, 1e-4)
            .is_err()
    );
    law.validate_committed(&f, &result.state, 1e-4).unwrap();
    assert_eq!(
        law.stress_and_potential(&f, &initial).unwrap(),
        (result.first_piola_pa, result.incremental_potential_j_m3)
    );
    assert_eq!(initial, law.initial_state());
}

#[test]
fn held_finite_deformation_relaxes_in_physical_time_and_refines_toward_the_continuous_solution() {
    let eta = 400e6;
    let f = coaxial(0.0, 0.18);
    let preload = kinetics(eta, 1.0, 0.02)
        .response(&f, &native().initial_state())
        .unwrap();
    let initial_excess = preload.equivalent_kirchhoff_pa - Y - H * preload.state.alpha;
    let duration = 2.0;
    let exact = initial_excess * (-A * duration / eta).exp();
    let mut previous_error = f64::INFINITY;
    for intervals in [10, 20, 40] {
        let dt = duration / f64::from(intervals);
        let law = kinetics(eta, 1.0, dt);
        let mut state = preload.state.clone();
        let mut energy = law.stored_energy(&f, &state).unwrap();
        let mut excess = initial_excess;
        for _ in 0..intervals {
            let next = law.response(&f, &state).unwrap();
            let actual = next.equivalent_kirchhoff_pa - Y - H * next.state.alpha;
            let expected = excess * (eta / dt) / (A + eta / dt);
            near(actual, expected, 2e-5, 2e-12);
            assert!(next.state.alpha > state.alpha);
            assert!(next.stored_energy_j_m3 < energy);
            assert!(next.plastic_dissipation_j_m3 > 0.0);
            // A held total F has exactly zero mechanical work. Stored energy
            // supplies real heat and a separate positive endpoint error.
            let numerical = energy - next.stored_energy_j_m3 - next.plastic_dissipation_j_m3;
            assert!(numerical > 0.0);
            excess = actual;
            energy = next.stored_energy_j_m3;
            state = next.state;
        }
        let error = (excess - exact).abs();
        assert!(error < previous_error);
        if previous_error.is_finite() {
            assert!(error < 0.55 * previous_error);
        }
        previous_error = error;
    }
}

#[test]
fn held_kirchhoff_load_accumulates_finite_plastic_creep_with_hardening_storage() {
    let (eta, dt, q) = (12e6, 0.1, 5e6);
    let law = kinetics(eta, 1.0, dt);
    let mut state = law.initial_state();
    let mut expected_alpha = 0.0;
    let mut final_f = [0.0; 9];
    for _ in 0..20 {
        let old = expected_alpha;
        expected_alpha = (old + dt * (q - Y) / eta) / (1.0 + dt * H / eta);
        let f = coaxial(0.02, expected_alpha + q / (3.0 * MU));
        let result = law.response(&f, &state).unwrap();
        near(result.state.alpha, expected_alpha, 2e-13, 2e-12);
        near(equivalent(&result.kirchhoff_pa), q, 2e-5, 2e-12);
        let gamma = expected_alpha - old;
        near(
            result.plastic_dissipation_j_m3,
            (q - H * expected_alpha) * gamma,
            2e-6,
            2e-11,
        );
        let elastic = 0.5 * K * 0.06_f64.powi(2) + q * q / (6.0 * MU);
        near(
            result.stored_energy_j_m3,
            elastic + 0.5 * H * expected_alpha.powi(2),
            2e-6,
            2e-11,
        );
        assert!(result.state.alpha > state.alpha);
        state = result.state;
        final_f = f;
    }
    assert!(state.alpha > 0.25 && final_f[0] > 1.3);
    near(
        state.inverse_plastic_cauchy_green[0],
        (-2.0 * expected_alpha).exp(),
        1e-12,
        1e-12,
    );
}

#[test]
fn nonlinear_flow_and_all_total_f_tangent_entries_differentiate_the_same_reduced_potential() {
    let (eta, dt, n) = (9e6, 0.3, 3.2);
    let law = kinetics(eta, n, dt);
    let initial = law.initial_state();
    let retained = law.response(&coaxial(0.01, 0.08), &initial).unwrap().state;
    for (f, old) in [
        (coaxial(0.02, 0.20), &initial),
        (
            [1.25, 0.13, 0.03, 0.02, 0.92, 0.07, 0.01, -0.02, 0.96],
            &retained,
        ),
    ] {
        let response = law.response(&f, old).unwrap();
        let gamma = response.plastic_increment;
        assert!(gamma > 0.0);
        let viscous = 2e6 * (eta * gamma / (dt * 2e6)).powf(1.0 / n);
        let actual = equivalent(&response.kirchhoff_pa) - Y - H * response.state.alpha;
        near(actual, viscous, 2e-6, 2e-12);
        near(
            response.plastic_dissipation_j_m3,
            (Y + viscous) * gamma,
            1e-7,
            1e-12,
        );
        near(
            response.incremental_potential_j_m3,
            response.stored_energy_j_m3 + Y * gamma + n / (n + 1.0) * viscous * gamma,
            1e-7,
            1e-12,
        );
        for column in 0..9 {
            let h = 2e-6;
            let mut plus = f;
            plus[column] += h;
            let mut minus = f;
            minus[column] -= h;
            let (p_plus, w_plus) = law.stress_and_potential(&plus, old).unwrap();
            let (p_minus, w_minus) = law.stress_and_potential(&minus, old).unwrap();
            near(
                (w_plus - w_minus) / (2.0 * h),
                response.first_piola_pa[column],
                0.25,
                3e-6,
            );
            for row in 0..9 {
                near(
                    (p_plus[row] - p_minus[row]) / (2.0 * h),
                    response.tangent_pa[row][column],
                    2.0,
                    3e-6,
                );
                near(
                    response.tangent_pa[row][column],
                    response.tangent_pa[column][row],
                    3e-5,
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn viscous_relaxation_is_objective_and_free_plastic_shape_preserves_complete_history() {
    let law = kinetics(25e6, 2.5, 0.2);
    let f = coaxial(0.01, 0.2);
    let old = law.response(&f, &law.initial_state()).unwrap().state;
    let saved = old.clone();
    let response = law.response(&f, &old).unwrap();
    assert!(response.plastic_increment > 0.0);
    let (s, c) = 1.15_f64.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let rotated = law.response(&multiply(&rotation, &f), &old).unwrap();
    near(rotated.state.alpha, response.state.alpha, 2e-13, 2e-12);
    for i in 0..6 {
        near(
            rotated.state.inverse_plastic_cauchy_green[i],
            response.state.inverse_plastic_cauchy_green[i],
            2e-13,
            2e-12,
        );
    }
    let expected = multiply(&rotation, &response.first_piola_pa);
    for i in 0..9 {
        near(rotated.first_piola_pa[i], expected[i], 2e-5, 2e-12);
    }
    near(
        rotated.plastic_dissipation_j_m3,
        response.plastic_dissipation_j_m3,
        2e-7,
        2e-12,
    );
    assert_eq!(old, saved);

    let free = coaxial(0.0, old.alpha);
    let unloaded = law.response(&free, &old).unwrap();
    assert_eq!(unloaded.state, old);
    assert_eq!(unloaded.plastic_increment, 0.0);
    assert_eq!(unloaded.plastic_dissipation_j_m3, 0.0);
    assert!(unloaded.first_piola_pa.iter().all(|x| x.abs() < 1e-6));
    near(
        unloaded.stored_energy_j_m3,
        0.5 * H * old.alpha.powi(2),
        1e-7,
        1e-12,
    );
}

#[test]
fn kinetics_and_retained_elastic_domains_refuse_without_equilibrium_projection_or_mutation() {
    let valid = FiniteOverstressParameters {
        viscosity_pa_s: 12e6,
        stress_scale_pa: 2e6,
        exponent: 2.0,
    };
    for p in [
        FiniteOverstressParameters {
            viscosity_pa_s: 0.0,
            ..valid
        },
        FiniteOverstressParameters {
            stress_scale_pa: -1.0,
            ..valid
        },
        FiniteOverstressParameters {
            exponent: 0.8,
            ..valid
        },
        FiniteOverstressParameters {
            exponent: f64::INFINITY,
            ..valid
        },
    ] {
        assert!(FiniteOverstressJ2::new(native(), p, 0.1).is_err());
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(FiniteOverstressJ2::new(native(), valid, dt).is_err());
    }
    assert!(
        FiniteOverstressJ2::new(
            native(),
            FiniteOverstressParameters {
                viscosity_pa_s: f64::MAX,
                exponent: 1.0,
                ..valid
            },
            f64::MIN_POSITIVE
        )
        .is_err()
    );
    let backbone = FiniteJ2::new(
        *native().parameters(),
        FiniteJ2Limits {
            maximum_elastic_log_strain: 0.05,
            ..*native().limits()
        },
    )
    .unwrap();
    let f = coaxial(0.0, 0.2);
    let state = backbone.initial_state();
    let saved = state.clone();
    assert!(backbone.response(&f, &state).is_ok());
    let slow = FiniteOverstressJ2::new(
        backbone,
        FiniteOverstressParameters {
            viscosity_pa_s: 1e10,
            exponent: 1.0,
            ..valid
        },
        1e-4,
    )
    .unwrap();
    assert!(matches!(
        slow.response(&f, &state),
        Err(FiniteJ2Error::Domain { .. })
    ));
    assert_eq!(state, saved);
    assert!(slow.response(&[0.0; 9], &state).is_err());
    assert_eq!(state, saved);
}
