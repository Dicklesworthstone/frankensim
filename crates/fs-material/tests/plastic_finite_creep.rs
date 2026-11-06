//! Independent finite Norton references; these synthetic laws are not measured
//! material calibration or a simultaneous creep-and-yield model.

use fs_material::creep::StrainHardeningCreep;
use fs_material::plastic::finite::creep::FiniteNortonCreep;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{
    FiniteJ2, FiniteJ2Error, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2State,
};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const Y: f64 = 2e6;
const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.1,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 3.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}
fn parameters(eta: f64, n: f64) -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: eta,
        stress_scale_pa: 2e6,
        exponent: n,
    }
}
fn law(eta: f64, n: f64, dt: f64) -> FiniteNortonCreep {
    FiniteNortonCreep::new(MU, K, limits(), parameters(eta, n), dt).unwrap()
}
fn plastic() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: 0.0,
        },
        limits(),
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
        "actual={actual:.16e}, expected={expected:.16e}, difference={:.8e}",
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

#[test]
fn independent_linear_finite_return_has_only_elastic_storage_and_physical_creep_heat() {
    let (eta, dt, mean, axial) = (60e6, 0.4, 0.03, 0.15);
    let law = law(eta, 1.0, dt);
    let old = law.initial_state();
    let f = coaxial(mean, axial);
    let r = law.response(&f, &old).unwrap();
    let gamma = 3.0 * MU * axial / (3.0 * MU + eta / dt);
    let elastic = axial - gamma;
    let q = 3.0 * MU * elastic;
    let stored = 0.5 * K * (3.0 * mean).powi(2) + 1.5 * MU * elastic * elastic;
    near(r.plastic_increment, gamma, 1e-14, 2e-13);
    near(r.state.alpha, gamma, 1e-14, 2e-13);
    near(
        r.state.inverse_plastic_cauchy_green[0],
        (-2.0 * gamma).exp(),
        1e-14,
        2e-13,
    );
    for i in [1, 2] {
        near(
            r.state.inverse_plastic_cauchy_green[i],
            gamma.exp(),
            1e-14,
            2e-13,
        );
    }
    for (i, tau) in [
        (0, 3.0 * K * mean + 2.0 * MU * elastic),
        (4, 3.0 * K * mean - MU * elastic),
        (8, 3.0 * K * mean - MU * elastic),
    ] {
        near(r.kirchhoff_pa[i], tau, 2e-6, 2e-13);
        near(r.first_piola_pa[i], tau / f[i], 2e-6, 2e-13);
    }
    near(r.equivalent_kirchhoff_pa, q, 2e-6, 2e-13);
    near(
        law.endpoint_equivalent_stress_pa(gamma).unwrap(),
        q,
        2e-6,
        2e-13,
    );
    near(r.stored_energy_j_m3, stored, 2e-7, 2e-13);
    near(r.plastic_dissipation_j_m3, q * gamma, 2e-7, 2e-13);
    near(
        r.incremental_potential_j_m3,
        stored + 0.5 * q * gamma,
        2e-7,
        2e-13,
    );
    assert!(r.incremental_potential_j_m3 < stored + r.plastic_dissipation_j_m3);
    law.validate_committed(&f, &r.state).unwrap();
    assert_eq!(
        law.stress_and_potential(&f, &old).unwrap(),
        (r.first_piola_pa, r.incremental_potential_j_m3)
    );
    assert_eq!(old, law.initial_state());
}

#[test]
fn primary_finite_return_matches_retained_history_quadratic_and_cubic_flow_potential() {
    let (eta, dt, mean, axial, reference) = (60e6, 0.4, 0.02, 0.18, 0.025);
    let hardening = StrainHardeningCreep {
        reference_strain: reference,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    };
    let material = law(eta, 1.0, dt).with_strain_hardening(hardening).unwrap();
    let old_alpha = 0.04_f64;
    let old = FiniteJ2State {
        inverse_plastic_cauchy_green: [
            (-2.0 * old_alpha).exp(),
            old_alpha.exp(),
            old_alpha.exp(),
            0.0,
            0.0,
            0.0,
        ],
        alpha: old_alpha,
        creep_alpha: 0.0,
    };
    material.validate_history(&old).unwrap();
    let f = coaxial(mean, axial);
    let r = material.response(&f, &old).unwrap();
    let modulus = eta / dt;
    let trial = 3.0 * MU * (axial - old_alpha);
    // dc^2 + (a+alpha_old+a*3mu*dt/eta)*dc - a*dt*q_trial/eta = 0.
    let b = reference + old_alpha + reference * 3.0 * MU / modulus;
    let drive = reference * trial / modulus;
    let dc = 2.0 * drive / (b + (b * b + 4.0 * drive).sqrt());
    let q = trial - 3.0 * MU * dc;
    let alpha = old_alpha + dc;
    let stored = 0.5 * K * (3.0 * mean).powi(2) + q * q / (6.0 * MU);
    let potential =
        modulus * (0.5 * (1.0 + old_alpha / reference) * dc * dc + dc.powi(3) / (3.0 * reference));
    near(r.plastic_increment, dc, 2e-13, 2e-11);
    near(r.state.alpha, alpha, 2e-13, 2e-11);
    near(r.equivalent_kirchhoff_pa, q, 3e-5, 2e-11);
    near(r.stored_energy_j_m3, stored, 2e-5, 2e-11);
    near(r.plastic_dissipation_j_m3, q * dc, 2e-5, 2e-11);
    near(r.flow_potential_j_m3, potential, 2e-5, 2e-10);
    near(
        r.incremental_potential_j_m3,
        stored + potential,
        2e-5,
        2e-10,
    );
    assert!(r.flow_potential_j_m3 < 0.5 * r.plastic_dissipation_j_m3);
    assert_eq!(r.state.creep_alpha, 0.0);
    assert_eq!(r.creep_increment, 0.0);
    assert_eq!(r.creep_dissipation_j_m3, 0.0);
    for (actual, expected) in r.state.inverse_plastic_cauchy_green.iter().zip([
        (-2.0 * alpha).exp(),
        alpha.exp(),
        alpha.exp(),
        0.0,
        0.0,
        0.0,
    ]) {
        near(*actual, expected, 2e-13, 2e-11);
    }
    assert_eq!(material.strain_hardening(), Some(&hardening));
    assert_eq!(
        material.stress_and_potential(&f, &old).unwrap(),
        (r.first_piola_pa, r.incremental_potential_j_m3)
    );
    material.validate_committed(&f, &r.state).unwrap();
}

#[test]
fn nonlinear_primary_creep_has_objective_retained_history_and_all_potential_derivatives() {
    let hardening = StrainHardeningCreep {
        reference_strain: 0.015,
        exponent: 0.8,
        residual_rate_fraction: 0.17,
    };
    let native = law(20e6, 2.5, 0.3);
    let material = native.clone().with_strain_hardening(hardening).unwrap();
    let old = material
        .response(&coaxial(0.01, 0.07), &material.initial_state())
        .unwrap()
        .state;
    let f = [1.21, 0.11, 0.02, 0.02, 0.94, 0.06, 0.01, -0.02, 0.96];
    let r = material.response(&f, &old).unwrap();
    let multiplier = 0.17 + 0.83 * (1.0 + r.state.alpha / 0.015).powf(-0.8);
    let expected_increment =
        0.3 * 2e6 / 20e6 * (r.equivalent_kirchhoff_pa / 2e6).powf(2.5) * multiplier;
    near(r.plastic_increment, expected_increment, 2e-13, 2e-11);
    assert!(old.alpha > 0.0 && r.state.alpha > old.alpha);
    for column in 0..9 {
        let h = 2e-6;
        let (mut plus, mut minus) = (f, f);
        plus[column] += h;
        minus[column] -= h;
        let (pp, wp) = material.stress_and_potential(&plus, &old).unwrap();
        let (pm, wm) = material.stress_and_potential(&minus, &old).unwrap();
        near((wp - wm) / (2.0 * h), r.first_piola_pa[column], 0.25, 3e-6);
        for row in 0..9 {
            near(
                (pp[row] - pm[row]) / (2.0 * h),
                r.tangent_pa[row][column],
                2.0,
                3e-6,
            );
            near(
                r.tangent_pa[row][column],
                r.tangent_pa[column][row],
                3e-5,
                2e-12,
            );
        }
    }
    let (s, c) = 1.15_f64.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let rotated = material.response(&multiply(&rotation, &f), &old).unwrap();
    near(rotated.state.alpha, r.state.alpha, 2e-13, 2e-11);
    near(
        rotated.flow_potential_j_m3,
        r.flow_potential_j_m3,
        2e-5,
        2e-10,
    );
    near(
        rotated.inelastic_dissipation_j_m3(),
        r.inelastic_dissipation_j_m3(),
        2e-5,
        2e-11,
    );
    for i in 0..6 {
        near(
            rotated.state.inverse_plastic_cauchy_green[i],
            r.state.inverse_plastic_cauchy_green[i],
            2e-13,
            2e-11,
        );
    }
    let expected = multiply(&rotation, &r.first_piola_pa);
    for i in 0..9 {
        near(rotated.first_piola_pa[i], expected[i], 3e-5, 3e-11);
    }
    for disabled in [
        StrainHardeningCreep {
            exponent: 0.0,
            ..hardening
        },
        StrainHardeningCreep {
            residual_rate_fraction: 1.0,
            ..hardening
        },
    ] {
        let reduced = native.clone().with_strain_hardening(disabled).unwrap();
        assert_eq!(
            reduced.response(&f, &old).unwrap(),
            native.response(&f, &old).unwrap()
        );
    }
}

#[test]
fn sustained_stress_below_plastic_yield_accumulates_finite_creep_without_artificial_hardening() {
    let (q, eta, dt) = (0.8e6, 8e6, 0.2);
    assert!(q < Y);
    let law = law(eta, 1.0, dt);
    let mut state = law.initial_state();
    let gamma = dt * q / eta;
    let mut final_f = ID;
    for step in 1..=20 {
        let alpha = f64::from(step) * gamma;
        let f = coaxial(0.01, alpha + q / (3.0 * MU));
        let r = law.response(&f, &state).unwrap();
        near(r.state.alpha, alpha, 2e-13, 2e-12);
        near(r.plastic_increment, gamma, 2e-13, 2e-12);
        near(r.equivalent_kirchhoff_pa, q, 3e-5, 2e-12);
        near(r.plastic_dissipation_j_m3, q * gamma, 2e-6, 2e-11);
        near(
            r.stored_energy_j_m3,
            0.5 * K * 0.03_f64.powi(2) + q * q / (6.0 * MU),
            2e-6,
            2e-11,
        );
        assert!(r.state.alpha > state.alpha);
        state = r.state;
        final_f = f;
    }
    assert!(state.alpha > 0.39 && final_f[0] > 1.5);
    near(
        state.inverse_plastic_cauchy_green[0],
        (-2.0 * state.alpha).exp(),
        2e-12,
        2e-12,
    );
    // The independently selected threshold law does not flow under this stress.
    let plastic = plastic();
    let elastic_f = coaxial(0.01, q / (3.0 * MU));
    let r = plastic
        .response(&elastic_f, &plastic.initial_state())
        .unwrap();
    assert_eq!(r.plastic_increment, 0.0);
    assert_eq!(r.state.alpha, 0.0);
}

#[test]
fn held_finite_shape_relaxes_below_yield_and_refines_to_the_continuous_norton_solution() {
    let (eta, axial) = (400e6, 0.08);
    let f = coaxial(0.0, axial);
    let initial = law(eta, 1.0, 0.01)
        .response(&f, &law(eta, 1.0, 0.01).initial_state())
        .unwrap();
    let exact = initial.equivalent_kirchhoff_pa * (-3.0 * MU * 2.0 / eta).exp();
    let mut previous_error = f64::INFINITY;
    for count in [10, 20, 40] {
        let dt = 2.0 / f64::from(count);
        let law = law(eta, 1.0, dt);
        let mut state = initial.state.clone();
        let mut q = initial.equivalent_kirchhoff_pa;
        let mut energy = law.stored_energy(&f, &state).unwrap();
        for _ in 0..count {
            let next = law.response(&f, &state).unwrap();
            near(
                next.equivalent_kirchhoff_pa,
                q / (1.0 + 3.0 * MU * dt / eta),
                3e-5,
                2e-12,
            );
            let defect = energy - next.stored_energy_j_m3 - next.plastic_dissipation_j_m3;
            near(
                defect,
                1.5 * MU * next.plastic_increment.powi(2),
                2e-6,
                2e-10,
            );
            assert!(defect > 0.0 && next.plastic_dissipation_j_m3 > 0.0);
            q = next.equivalent_kirchhoff_pa;
            energy = next.stored_energy_j_m3;
            state = next.state;
        }
        let error = (q - exact).abs();
        if previous_error.is_finite() {
            assert!(error < 0.55 * previous_error);
        }
        previous_error = error;
    }
    let creep = law(eta, 1.0, 0.5);
    let threshold = FiniteOverstressJ2::new(plastic(), parameters(eta, 1.0), 0.5).unwrap();
    let mut c = creep.initial_state();
    let mut p = threshold.initial_state();
    let (mut cq, mut pq) = (0.0, 0.0);
    for _ in 0..24 {
        let cr = creep.response(&f, &c).unwrap();
        let pr = threshold.response(&f, &p).unwrap();
        cq = cr.equivalent_kirchhoff_pa;
        pq = pr.equivalent_kirchhoff_pa;
        c = cr.state;
        p = pr.state;
    }
    assert!(cq < 0.1 * Y && pq > Y);
}

#[test]
fn zero_and_resolvable_near_zero_stress_have_the_smooth_linear_creep_tangent_and_potential() {
    let (eta, dt) = (60e6, 0.4);
    let law = law(eta, 1.0, dt);
    let r = law.response(&ID, &law.initial_state()).unwrap();
    let ratio = (eta / dt) / (3.0 * MU + eta / dt);
    assert_eq!(r.state, law.initial_state());
    assert_eq!(r.first_piola_pa, [0.0; 9]);
    assert_eq!(r.incremental_potential_j_m3, 0.0);
    for ij in 0..9 {
        let (i, j) = (ij / 3, ij % 3);
        for kl in 0..9 {
            let (k, l) = (kl / 3, kl % 3);
            let delta = |a: usize, b: usize| if a == b { 1.0 } else { 0.0 };
            let expected = K * delta(i, j) * delta(k, l)
                + MU * ratio
                    * (delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k)
                        - (2.0 / 3.0) * delta(i, j) * delta(k, l));
            near(r.tangent_pa[ij][kl], expected, 2e-6, 2e-14);
        }
    }
    let tiny = coaxial(0.0, 1e-15);
    let logs = [tiny[0].ln(), tiny[4].ln(), tiny[8].ln()];
    let mean = logs.iter().sum::<f64>() / 3.0;
    let dev = logs.map(|x| x - mean);
    let q_trial = (6.0 * dev.iter().map(|x| x * x).sum::<f64>()).sqrt() * MU;
    assert!(q_trial > 0.0 && q_trial < 64.0 * f64::EPSILON * (2.0 * MU));
    let tiny_response = law.response(&tiny, &law.initial_state()).unwrap();
    assert!(tiny_response.plastic_increment > 0.0);
    near(
        tiny_response.plastic_increment,
        q_trial / (3.0 * MU + eta / dt),
        1e-29,
        1e-11,
    );
    let trial_energy = MU * dev.iter().map(|x| x * x).sum::<f64>();
    let volume_energy = 0.5 * K * logs.iter().sum::<f64>().powi(2);
    near(
        tiny_response.incremental_potential_j_m3,
        volume_energy + ratio * trial_energy,
        1e-34,
        1e-10,
    );
    near(
        tiny_response.plastic_dissipation_j_m3,
        2.0 * ratio * (1.0 - ratio) * trial_energy,
        1e-34,
        1e-10,
    );
    let nonlinear = FiniteNortonCreep::new(MU, K, limits(), parameters(eta, 3.0), dt).unwrap();
    let n = nonlinear.response(&ID, &nonlinear.initial_state()).unwrap();
    near(n.tangent_pa[1][1], MU, 1e-7, 1e-14);
}

#[test]
fn nonlinear_finite_creep_has_objective_history_reversal_and_the_full_potential_tangent() {
    let law = law(20e6, 2.7, 0.3);
    let initial = law.initial_state();
    let old = law.response(&coaxial(0.01, 0.07), &initial).unwrap().state;
    let f = [1.21, 0.11, 0.02, 0.02, 0.94, 0.06, 0.01, -0.02, 0.96];
    let r = law.response(&f, &old).unwrap();
    assert!(r.plastic_increment > 0.0);
    let q = 2e6 * (20e6 * r.plastic_increment / (0.3 * 2e6)).powf(1.0 / 2.7);
    near(r.equivalent_kirchhoff_pa, q, 2e-6, 2e-12);
    near(
        r.plastic_dissipation_j_m3,
        q * r.plastic_increment,
        2e-7,
        2e-12,
    );
    near(
        r.incremental_potential_j_m3,
        r.stored_energy_j_m3 + (2.7 / 3.7) * q * r.plastic_increment,
        2e-7,
        2e-12,
    );
    for column in 0..9 {
        let h = 2e-6;
        let mut plus = f;
        plus[column] += h;
        let mut minus = f;
        minus[column] -= h;
        let (pp, wp) = law.stress_and_potential(&plus, &old).unwrap();
        let (pm, wm) = law.stress_and_potential(&minus, &old).unwrap();
        near((wp - wm) / (2.0 * h), r.first_piola_pa[column], 0.25, 3e-6);
        for row in 0..9 {
            near(
                (pp[row] - pm[row]) / (2.0 * h),
                r.tangent_pa[row][column],
                2.0,
                3e-6,
            );
            near(
                r.tangent_pa[row][column],
                r.tangent_pa[column][row],
                3e-5,
                2e-12,
            );
        }
    }
    let (s, c) = 1.15_f64.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let rotated = law.response(&multiply(&rotation, &f), &old).unwrap();
    near(rotated.state.alpha, r.state.alpha, 2e-13, 2e-12);
    for i in 0..6 {
        near(
            rotated.state.inverse_plastic_cauchy_green[i],
            r.state.inverse_plastic_cauchy_green[i],
            2e-13,
            2e-12,
        );
    }
    let expected = multiply(&rotation, &r.first_piola_pa);
    for i in 0..9 {
        near(rotated.first_piola_pa[i], expected[i], 3e-5, 3e-12);
    }
    let reverse = law.response(&coaxial(0.0, -0.10), &old).unwrap();
    assert!(reverse.state.alpha > old.alpha);
    assert!(reverse.state.inverse_plastic_cauchy_green[0] > old.inverse_plastic_cauchy_green[0]);
    assert!(reverse.kirchhoff_pa[0] < 0.0 && reverse.plastic_dissipation_j_m3 > 0.0);
}

#[test]
fn creep_is_explicit_and_invalid_or_unresolvable_steps_leave_retained_history_unchanged() {
    assert!(
        FiniteJ2::new(
            FiniteJ2Parameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: K,
                yield_stress_pa: 0.0,
                hardening_modulus_pa: 0.0,
            },
            limits()
        )
        .is_err()
    );
    for (mu, k) in [(0.0, K), (MU, -1.0), (f64::INFINITY, K)] {
        assert!(FiniteNortonCreep::new(mu, k, limits(), parameters(60e6, 1.0), 1.0).is_err());
    }
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(FiniteNortonCreep::new(MU, K, limits(), parameters(60e6, 1.0), dt).is_err());
    }
    for p in [
        parameters(0.0, 1.0),
        parameters(60e6, 0.5),
        FiniteOverstressParameters {
            stress_scale_pa: 0.0,
            ..parameters(60e6, 1.0)
        },
    ] {
        assert!(FiniteNortonCreep::new(MU, K, limits(), p, 1.0).is_err());
    }
    let ordinary = law(60e6, 1.0, 0.4);
    let f = coaxial(0.0, 0.12);
    let old = ordinary
        .response(&f, &ordinary.initial_state())
        .unwrap()
        .state;
    let saved = old.clone();
    let unresolved = law(1e30, 1.0, 1.0);
    assert!(matches!(
        unresolved.response(&f, &old),
        Err(FiniteJ2Error::Numerical(_))
    ));
    let bounded = FiniteNortonCreep::new(
        MU,
        K,
        FiniteJ2Limits {
            maximum_accumulated_plastic_strain: old.alpha + 1e-5,
            ..limits()
        },
        parameters(60e6, 1.0),
        0.4,
    )
    .unwrap();
    assert!(matches!(
        bounded.response(&f, &old),
        Err(FiniteJ2Error::Domain { .. })
    ));
    assert!(ordinary.response(&[0.0; 9], &old).is_err());
    assert_eq!(old, saved);
    let mut invalid = old.clone();
    invalid.inverse_plastic_cauchy_green[0] = -1.0;
    assert!(ordinary.validate_committed(&f, &invalid).is_err());
    ordinary.validate_committed(&f, &old).unwrap();
}
