//! G0/G1: nonlinear crush hardening through the original finite porous return.
//! Manufactured constitutive data establish numerical behavior, not calibration.

use fs_material::plastic::finite::porous::{
    FinitePorous, FinitePorousLimits, FinitePorousParameters,
};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};

const MU: f64 = 80.0;
const K: f64 = 160.0;
const Y: f64 = 2.0;
const H: f64 = 8.0;
const A: f64 = 1.2;
const Q: f64 = 12.0;
const B: f64 = 3.0;

fn base(sensitivity: f64, scale: f64) -> FinitePorous {
    FinitePorous::new(
        FinitePorousParameters {
            shear_modulus_pa: MU * scale,
            bulk_modulus_pa: K * scale,
            yield_stress_pa: Y * scale,
            hardening_modulus_pa: H * scale,
            pressure_sensitivity: sensitivity,
        },
        FinitePorousLimits {
            minimum_jacobian: 0.2,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn law() -> FinitePorous {
    base(A, 1.0).with_voce_hardening(Q, B).unwrap()
}

fn diagonal(value: f64) -> [f64; 9] {
    [value, 0.0, 0.0, 0.0, value, 0.0, 0.0, 0.0, value]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, allowance {tolerance:.3e}"
    );
}

fn radius(alpha: f64) -> f64 {
    Y + H * alpha + Q * -(-B * alpha).exp_m1()
}

fn hardening_energy(alpha: f64) -> f64 {
    0.5 * H * alpha * alpha + Q * (alpha + (-B * alpha).exp_m1() / B)
}

#[test]
fn hydrostatic_saturation_retains_plastic_volume_and_separate_hardening_storage() {
    let law = law();
    let old = law.initial_state();
    let stretch = 0.8_f64;
    let beta = A / (1.0 + A * A / 9.0).sqrt();
    // Independent one-dimensional hydrostatic equation, without the native
    // normalized pressure/deviatoric return or its Newton derivative.
    let mut lower = 0.0;
    let mut upper = -3.0 * stretch.ln() / beta;
    for _ in 0..100 {
        let middle = 0.5 * lower + 0.5 * upper;
        if K * beta * (-3.0 * stretch.ln() - beta * middle) > radius(middle) {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    let alpha = 0.5 * lower + 0.5 * upper;
    let tau = -radius(alpha) / beta;
    let jp = (-beta * alpha).exp();
    let response = law.response(&diagonal(stretch), &old).unwrap();
    near(response.state.alpha, alpha, 2e-12);
    near(response.plastic_increment, alpha, 2e-12);
    near(response.plastic_jacobian, jp, 2e-12);
    near(response.jacobian, stretch.powi(3), 2e-14);
    near(response.kirchhoff_pa[0], tau, 2e-10);
    near(response.first_piola_pa[0], tau / stretch, 2e-10);
    near(response.yield_residual_pa, 0.0, 2e-10);
    let stored = 0.5 * tau * tau / K + hardening_energy(alpha);
    near(response.stored_energy_j_m3, stored, 2e-10);
    near(response.plastic_dissipation_j_m3, Y * alpha, 2e-12);
    near(
        response.incremental_potential_j_m3,
        stored + Y * alpha,
        2e-10,
    );
    assert_eq!(old, law.initial_state());
    law.validate_committed(&diagonal(stretch), &response.state, 2e-9)
        .unwrap();

    let released_stretch = (-beta * alpha / 3.0).exp();
    let released = law
        .response(&diagonal(released_stretch), &response.state)
        .unwrap();
    assert_eq!(released.state, response.state);
    assert_eq!(released.plastic_dissipation_j_m3, 0.0);
    near(released.plastic_jacobian, jp, 2e-12);
    near(released.jacobian, jp, 2e-12);
    near(released.stored_energy_j_m3, hardening_energy(alpha), 2e-10);
    for stress in released.first_piola_pa {
        near(stress, 0.0, 2e-10);
    }
    let reload = law
        .response(
            &diagonal((stretch * released_stretch).sqrt()),
            &released.state,
        )
        .unwrap();
    assert_eq!(reload.state, released.state);
    near(reload.kirchhoff_pa[0], 0.5 * tau, 2e-10);

    let linear = base(A, 1.0).response(&diagonal(stretch), &old).unwrap();
    assert!(response.state.alpha < linear.state.alpha);
    assert!(response.plastic_jacobian > linear.plastic_jacobian);
    assert!(response.kirchhoff_pa[0].abs() > linear.kirchhoff_pa[0].abs());
}

#[test]
fn nonlinear_return_has_all_nine_consistent_columns_and_incremental_energy_gradient() {
    let law = law();
    let compacted = law
        .response(&diagonal(0.8), &law.initial_state())
        .unwrap()
        .state;
    let previous_f = [0.89, 0.11, -0.02, 0.03, 0.93, 0.06, -0.01, 0.04, 0.84];
    let old = law.response(&previous_f, &compacted).unwrap().state;
    let saved = old.clone();
    let f = [1.03, -0.06, 0.04, 0.09, 0.91, 0.03, 0.02, -0.01, 0.87];
    let response = law.response(&f, &old).unwrap();
    assert!(response.plastic_increment > 0.01);
    let scale = response
        .tangent_pa
        .iter()
        .flatten()
        .fold(1.0_f64, |s, v| s.max(v.abs()));
    for column in 0..9 {
        let h = 2e-6;
        let mut plus = f;
        let mut minus = f;
        plus[column] += h;
        minus[column] -= h;
        let (p_plus, w_plus) = law.stress_and_potential(&plus, &old).unwrap();
        let (p_minus, w_minus) = law.stress_and_potential(&minus, &old).unwrap();
        for row in 0..9 {
            near(
                response.tangent_pa[row][column],
                (p_plus[row] - p_minus[row]) / (2.0 * h),
                2e-6 * scale,
            );
            near(
                response.tangent_pa[row][column],
                response.tangent_pa[column][row],
                2e-8 * scale,
            );
        }
        near(
            (w_plus - w_minus) / (2.0 * h),
            response.first_piola_pa[column],
            2e-6,
        );
    }
    assert_eq!(old, saved);
}

#[test]
fn zero_saturation_preserves_native_response_and_zero_sensitivity_matches_voce_j2() {
    let linear = base(A, 1.0);
    let explicit_zero = linear.clone().with_voce_hardening(0.0, B).unwrap();
    assert_eq!(linear, explicit_zero);
    let f = [0.91, 0.08, -0.02, 0.01, 0.94, 0.04, 0.0, 0.03, 0.87];
    assert_eq!(
        linear.response(&f, &linear.initial_state()).unwrap(),
        explicit_zero.response(&f, &linear.initial_state()).unwrap()
    );
    let porous = base(0.0, 1.0).with_voce_hardening(Q, B).unwrap();
    let j2 = FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.2,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
    .with_voce_hardening(Q, B)
    .unwrap();
    let mut porous_state = porous.initial_state();
    let mut j2_state = j2.initial_state();
    for f in [f, [1.02, -0.06, 0.0, 0.01, 0.91, 0.0, 0.0, 0.0, 0.96]] {
        let p = porous.response(&f, &porous_state).unwrap();
        let j = j2.response(&f, &j2_state).unwrap();
        near(p.state.alpha, j.state.alpha, 2e-10);
        near(p.plastic_jacobian, 1.0, 2e-10);
        near(p.stored_energy_j_m3, j.stored_energy_j_m3, 2e-9);
        near(p.plastic_dissipation_j_m3, j.plastic_dissipation_j_m3, 2e-9);
        for i in 0..9 {
            near(p.first_piola_pa[i], j.first_piola_pa[i], 2e-8);
            for k in 0..9 {
                near(p.tangent_pa[i][k], j.tangent_pa[i][k], 2e-7);
            }
        }
        porous_state = p.state;
        j2_state = j.state;
    }
}

#[test]
fn saturation_keeps_stress_scaling_and_refuses_invalid_or_replaced_coefficients() {
    let native = law();
    let scaled = base(A, 1e6).with_voce_hardening(Q * 1e6, B).unwrap();
    let f = diagonal(0.8);
    let a = native.response(&f, &native.initial_state()).unwrap();
    let b = scaled.response(&f, &scaled.initial_state()).unwrap();
    near(a.state.alpha, b.state.alpha, 2e-12);
    near(a.plastic_jacobian, b.plastic_jacobian, 2e-12);
    near(a.stored_energy_j_m3, b.stored_energy_j_m3 / 1e6, 2e-10);
    near(
        a.plastic_dissipation_j_m3,
        b.plastic_dissipation_j_m3 / 1e6,
        2e-10,
    );
    for i in 0..9 {
        near(a.first_piola_pa[i], b.first_piola_pa[i] / 1e6, 2e-9);
        for j in 0..9 {
            near(a.tangent_pa[i][j], b.tangent_pa[i][j] / 1e6, 2e-8);
        }
    }
    for (q, rate) in [
        (-1.0, B),
        (Q, 0.0),
        (f64::NAN, B),
        (Q, f64::INFINITY),
        (f64::MAX, B),
    ] {
        assert!(base(A, 1.0).with_voce_hardening(q, rate).is_err());
    }
    assert_eq!(native.clone().with_voce_hardening(Q, B).unwrap(), native);
    assert!(native.clone().with_voce_hardening(Q + 1.0, B).is_err());
    assert!(native.with_voce_hardening(0.0, B).is_err());
}
