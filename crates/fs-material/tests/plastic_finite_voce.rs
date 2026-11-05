//! Independent scalar oracles and full tensor derivatives for finite saturation.
//! Synthetic coefficients verify the implementation, not a fitted material.

use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2Response};

const MU: f64 = 80e6;
const K: f64 = 160e6;
const Y: f64 = 2e6;
const H: f64 = 8e6;
const Q: f64 = 15e6;
const B: f64 = 12.0;
const ID: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn linear() -> FiniteJ2 {
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
            maximum_accumulated_plastic_strain: 5.0,
        },
    )
    .unwrap()
}
fn law() -> FiniteJ2 {
    linear().with_voce_hardening(Q, B).unwrap()
}
fn kinetics(exponent: f64) -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: 16e6,
        stress_scale_pa: 2e6,
        exponent,
    }
}
fn coaxial(axial: f64, mean: f64) -> [f64; 9] {
    [
        (mean + axial).exp(), 0.0, 0.0,
        0.0, (mean - 0.5 * axial).exp(), 0.0,
        0.0, 0.0, (mean - 0.5 * axial).exp(),
    ]
}
fn radius(alpha: f64) -> f64 {
    Y + H * alpha + Q * (1.0 - (-B * alpha).exp())
}
fn storage(alpha: f64) -> f64 {
    0.5 * H * alpha * alpha + Q * (alpha - (1.0 - (-B * alpha).exp()) / B)
}
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "actual={a:.16e}, expected={b:.16e}, difference={:.8e}",
        (a - b).abs()
    );
}
fn bisect(mut low: f64, mut high: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(low) <= 0.0 && residual(high) >= 0.0);
    for _ in 0..160 {
        let middle = 0.5 * (low + high);
        if residual(middle) < 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}
fn derivatives(
    r: &FiniteJ2Response,
    f: &[f64; 9],
    update: impl Fn(&[f64; 9]) -> ([f64; 9], f64),
) {
    for column in 0..9 {
        let step = 2e-6;
        let mut plus = *f;
        let mut minus = *f;
        plus[column] += step;
        minus[column] -= step;
        let (pp, wp) = update(&plus);
        let (pm, wm) = update(&minus);
        near((wp - wm) / (2.0 * step), r.first_piola_pa[column], 0.3, 4e-6);
        for row in 0..9 {
            near(
                (pp[row] - pm[row]) / (2.0 * step),
                r.tangent_pa[row][column],
                3.0,
                4e-6,
            );
            near(r.tangent_pa[row][column], r.tangent_pa[column][row], 3e-5, 3e-12);
        }
    }
}
fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| (0..3).map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3]).sum())
}

#[test]
fn finite_saturation_matches_independent_coaxial_flow_and_recoverable_energy() {
    let material = law();
    let old = material.initial_state();
    let axial = 0.24;
    let mean = 0.025;
    let f = coaxial(axial, mean);
    let gamma = bisect(0.0, axial, |g| 3.0 * MU * g + radius(g) - 3.0 * MU * axial);
    let r = material.response(&f, &old).unwrap();
    near(r.plastic_increment, gamma, 2e-15, 2e-13);
    near(r.equivalent_kirchhoff_pa, radius(gamma), 2e-6, 2e-13);
    near(r.state.inverse_plastic_cauchy_green[0], (-2.0 * gamma).exp(), 2e-14, 2e-14);
    let expected_storage = 0.5 * K * (3.0 * mean).powi(2)
        + 1.5 * MU * (axial - gamma).powi(2) + storage(gamma);
    near(r.stored_energy_j_m3, expected_storage, 2e-7, 2e-13);
    near(r.plastic_dissipation_j_m3, Y * gamma, 2e-8, 2e-13);
    near(r.incremental_potential_j_m3, expected_storage + Y * gamma, 2e-7, 2e-13);
    near(material.stored_energy(&f, &r.state).unwrap(), expected_storage, 2e-7, 2e-13);
    material.validate_committed(&f, &r.state, 1e-4).unwrap();
    assert_eq!(old, material.initial_state());
}

#[test]
fn finite_saturation_retains_history_through_unloading_and_reversal() {
    let material = law();
    let mut state = material.initial_state();
    let mut plastic = 0.0;
    let mut alpha = 0.0;
    for axial in [0.1, 0.15, 0.13, 0.0, -0.1, -0.18, 0.0, 0.12] {
        let trial = 3.0 * MU * (axial - plastic);
        let increment = if trial.abs() <= radius(alpha) {
            0.0
        } else {
            bisect(0.0, trial.abs() / (3.0 * MU), |g| {
                3.0 * MU * g + radius(alpha + g) - trial.abs()
            })
        };
        plastic += trial.signum() * increment;
        alpha += increment;
        let f = coaxial(axial, 0.0);
        let r = material.response(&f, &state).unwrap();
        near(r.state.alpha, alpha, 2e-14, 3e-13);
        near(r.kirchhoff_pa[0], 2.0 * MU * (axial - plastic), 3e-6, 3e-13);
        near(r.state.inverse_plastic_cauchy_green[0], (-2.0 * plastic).exp(), 3e-14, 3e-13);
        near(r.plastic_dissipation_j_m3, Y * increment, 2e-7, 3e-13);
        let held = material.response(&f, &r.state).unwrap();
        assert_eq!(held.state, r.state);
        assert_eq!(held.plastic_dissipation_j_m3, 0.0);
        state = r.state;
    }
}

#[test]
fn finite_saturation_full_tangent_and_potential_are_objective_on_noncoaxial_history() {
    let material = law();
    let old = material.response(&coaxial(0.1, 0.0), &material.initial_state()).unwrap().state;
    let f = [1.21, 0.13, 0.03, 0.02, 0.92, 0.07, 0.01, -0.02, 0.96];
    let r = material.response(&f, &old).unwrap();
    derivatives(&r, &f, |trial| material.stress_and_potential(trial, &old).unwrap());
    let (s, c) = 0.63_f64.sin_cos();
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let rotated = material.response(&multiply(&rotation, &f), &old).unwrap();
    for (a, b) in rotated.first_piola_pa.into_iter().zip(multiply(&rotation, &r.first_piola_pa)) {
        near(a, b, 5e-6, 3e-12);
    }
    for i in 0..6 {
        near(rotated.state.inverse_plastic_cauchy_green[i], r.state.inverse_plastic_cauchy_green[i], 5e-14, 3e-13);
    }
    near(rotated.state.alpha, r.state.alpha, 5e-15, 3e-13);
    near(rotated.stored_energy_j_m3, r.stored_energy_j_m3, 2e-7, 3e-13);
}

#[test]
fn finite_voce_perzyna_satisfies_joint_flow_and_all_potential_derivatives() {
    let f = [1.22, 0.13, 0.03, 0.02, 0.92, 0.07, 0.01, -0.02, 0.96];
    for n in [1.0, 3.2] {
        let material = FiniteOverstressJ2::new(law(), kinetics(n), 0.3).unwrap();
        let old = material.initial_state();
        let r = material.response(&f, &old).unwrap();
        let viscous = 2e6 * (16e6 * r.plastic_increment / (0.3 * 2e6)).powf(1.0 / n);
        near(r.equivalent_kirchhoff_pa, radius(r.state.alpha) + viscous, 3e-6, 3e-13);
        near(r.plastic_dissipation_j_m3, (Y + viscous) * r.plastic_increment, 2e-7, 3e-13);
        near(
            r.incremental_potential_j_m3 - r.stored_energy_j_m3,
            (Y + n / (n + 1.0) * viscous) * r.plastic_increment,
            2e-7,
            3e-13,
        );
        derivatives(&r, &f, |trial| material.stress_and_potential(trial, &old).unwrap());
        let held = material.response(&f, &r.state).unwrap();
        assert!(held.state.alpha > r.state.alpha);
        assert!(held.equivalent_kirchhoff_pa < r.equivalent_kirchhoff_pa);
        assert!(held.stored_energy_j_m3 + held.plastic_dissipation_j_m3 < r.stored_energy_j_m3);
    }
}

#[test]
fn simultaneous_norton_and_voce_match_complementarity_and_independent_scalar_equations() {
    for n in [1.0, 2.5] {
        let material = FiniteNortonJ2::new(law(), kinetics(n), 0.1, 5.0).unwrap();
        let old = material.initial_state();
        for axial in [0.003, 0.2] {
            let f = coaxial(axial, 0.02);
            let creep_at = |q: f64| 0.1 * (2e6 / 16e6) * (q / 2e6).powf(n);
            let trial = 3.0 * MU * axial;
            let active = trial > Y + 3.0 * MU * creep_at(Y);
            let (q, dp) = if active {
                let p = bisect(0.0, axial, |p| {
                    radius(p) + 3.0 * MU * (p + creep_at(radius(p))) - trial
                });
                (radius(p), p)
            } else {
                (bisect(0.0, trial, |q| q + 3.0 * MU * creep_at(q) - trial), 0.0)
            };
            let dc = creep_at(q);
            let r = material.response(&f, &old).unwrap();
            near(r.plastic_increment, dp, 2e-14, 3e-12);
            near(r.creep_increment, dc, 2e-14, 3e-12);
            near(r.equivalent_kirchhoff_pa, q, 3e-6, 3e-12);
            near(r.state.inverse_plastic_cauchy_green[0], (-2.0 * (dp + dc)).exp(), 3e-14, 3e-12);
            near(r.plastic_dissipation_j_m3, Y * dp, 2e-7, 3e-12);
            near(r.creep_dissipation_j_m3, q * dc, 2e-7, 3e-12);
            material.validate_committed(&f, &r.state, 1e-4).unwrap();
            derivatives(&r, &f, |trial| material.stress_and_potential(trial, &old).unwrap());
        }
    }
}

#[test]
fn saturation_energy_resolves_tiny_accumulations_and_zero_q_preserves_the_original_law() {
    let native = linear();
    let zero = native.clone().with_voce_hardening(0.0, B).unwrap();
    assert_eq!(native, zero);
    let f = coaxial(0.19, 0.01);
    assert_eq!(
        native.response(&f, &native.initial_state()).unwrap(),
        zero.response(&f, &zero.initial_state()).unwrap(),
    );
    let material = law();
    let mut state = material.initial_state();
    state.alpha = 1e-12;
    let expected = 0.5 * (H + Q * B) * state.alpha * state.alpha;
    near(material.stored_energy(&ID, &state).unwrap(), expected, 0.0, 1e-11);
    let mut limits = *material.limits();
    limits.maximum_accumulated_plastic_strain = 0.01;
    let bounded = FiniteJ2::new(*material.parameters(), limits).unwrap()
        .with_voce_hardening(Q, B).unwrap();
    let old = bounded.initial_state();
    assert!(bounded.response(&coaxial(0.2, 0.0), &old).is_err());
    assert_eq!(old, bounded.initial_state());
    for (q, b) in [
        (-1.0, B), (Q, 0.0), (Q, -1.0), (f64::NAN, B),
        (Q, f64::INFINITY), (f64::MAX, 2.0),
    ] {
        assert!(linear().with_voce_hardening(q, b).is_err());
    }
    assert!(material.with_voce_hardening(Q * 2.0, B).is_err());
}
