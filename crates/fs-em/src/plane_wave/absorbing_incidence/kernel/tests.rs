//! G1/G3: independent complex Fresnel, six-component Maxwell and Poynting checks.
use super::*;
use crate::plane_wave::{Polarization, interface_response};

fn data() -> Input {
    Input {
        omega: 1.7,
        kx: 0.45,
        epsilon_i: C64::new(2.5, -0.8),
        epsilon_t: C64::new(0.7, -0.2),
        mu_i: 1.4,
        mu_t: 0.9,
        electric: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn complex_near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:?} != {b:?}, tolerance {tolerance:e}"
    );
}
fn cross(a: [C64; 3], b: [C64; 3]) -> [C64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn power(field: InterfaceBoundaryField) -> f64 {
    0.5 * cross(field.electric_v_m, field.magnetic_a_m.map(C64::conj))[2].re
}
fn combined(a: InterfaceBoundaryField, b: InterfaceBoundaryField) -> InterfaceBoundaryField {
    InterfaceBoundaryField {
        electric_v_m: core::array::from_fn(|c| a.electric_v_m[c] + b.electric_v_m[c]),
        magnetic_a_m: core::array::from_fn(|c| a.magnetic_a_m[c] + b.magnetic_a_m[c]),
    }
}
fn dimensional_root(omega: f64, kx: f64, epsilon: C64, mu: f64) -> C64 {
    let q = (epsilon.scale(omega * omega * mu) - C64::from_re(kx * kx)).sqrt();
    if q.im > 0.0 { -q } else { q }
}
fn maxwell(field: InterfaceBoundaryField, kx: f64, q: C64, omega: f64, epsilon: C64, mu: f64) {
    let k = [C64::from_re(kx), C64::ZERO, q];
    let ke = cross(k, field.electric_v_m);
    let kh = cross(k, field.magnetic_a_m);
    for c in 0..3 {
        complex_near(ke[c], field.magnetic_a_m[c].scale(omega * mu), 3e-13);
        complex_near(
            kh[c],
            -(epsilon * field.electric_v_m[c]).scale(omega),
            3e-13,
        );
    }
}
fn boundary_and_power(input: Input, result: AbsorbingInterfaceReflection) {
    let left = combined(result.incident, result.reflected);
    for c in 0..2 {
        complex_near(
            left.electric_v_m[c],
            result.transmitted.electric_v_m[c],
            3e-13,
        );
        complex_near(
            left.magnetic_a_m[c],
            result.transmitted.magnetic_a_m[c],
            3e-13,
        );
    }
    complex_near(
        input.epsilon_i * left.electric_v_m[2],
        input.epsilon_t * result.transmitted.electric_v_m[2],
        3e-13,
    );
    complex_near(
        left.magnetic_a_m[2].scale(input.mu_i),
        result.transmitted.magnetic_a_m[2].scale(input.mu_t),
        3e-13,
    );
    near(
        power(result.incident),
        result.incident_normal_power_w_m2,
        2e-14,
    );
    near(
        power(result.reflected),
        result.reflected_normal_power_w_m2,
        2e-14,
    );
    near(
        power(result.transmitted),
        result.transmitted_normal_power_w_m2,
        2e-14,
    );
    let cross_power = 0.5
        * (cross(
            result.incident.electric_v_m,
            result.reflected.magnetic_a_m.map(C64::conj),
        )[2]
        .re + cross(
            result.reflected.electric_v_m,
            result.incident.magnetic_a_m.map(C64::conj),
        )[2]
        .re);
    near(cross_power, result.interference_normal_power_w_m2, 2e-14);
    near(power(left), power(result.transmitted), 3e-14);
    near(
        result.incident_normal_power_w_m2
            + result.reflected_normal_power_w_m2
            + result.interference_normal_power_w_m2,
        result.transmitted_normal_power_w_m2,
        3e-14,
    );
}

#[test]
fn lossy_oblique_ports_match_complex_fresnel_and_both_maxwell_curls() {
    for kx in [-1.4, -0.45, 0.0, 0.45, 1.4] {
        let input = Input { kx, ..data() };
        let result = solve(input).unwrap();
        let qi = dimensional_root(input.omega, kx, input.epsilon_i, input.mu_i);
        let qt = dimensional_root(input.omega, kx, input.epsilon_t, input.mu_t);
        complex_near(result.incident_normal_wavenumber_rad_m, qi, 2e-14);
        complex_near(result.transmitted_normal_wavenumber_rad_m, qt, 2e-14);
        let yi = [
            input.epsilon_i.scale(input.omega) / qi,
            qi.scale(1.0 / (input.omega * input.mu_i)),
        ];
        let yt = [
            input.epsilon_t.scale(input.omega) / qt,
            qt.scale(1.0 / (input.omega * input.mu_t)),
        ];
        for c in 0..2 {
            complex_near(
                result.reflected_tangential_e_gain[c],
                (yi[c] - yt[c]) / (yi[c] + yt[c]),
                2e-14,
            );
            complex_near(
                result.transmitted_tangential_e_gain[c],
                yi[c].scale(2.0) / (yi[c] + yt[c]),
                2e-14,
            );
        }
        maxwell(
            result.incident,
            kx,
            qi,
            input.omega,
            input.epsilon_i,
            input.mu_i,
        );
        maxwell(
            result.reflected,
            kx,
            -qi,
            input.omega,
            input.epsilon_i,
            input.mu_i,
        );
        maxwell(
            result.transmitted,
            kx,
            qt,
            input.omega,
            input.epsilon_t,
            input.mu_t,
        );
        boundary_and_power(input, result);
    }
}

#[test]
fn interference_can_make_transmitted_power_exceed_the_separate_incident_power_without_gain() {
    // Exact dimensional oracle: k_i=1-2i, k_t=1, r=(1-i)/2,
    // t=(3-i)/2. Fluxes are +1/2, -1/4, +1, +5/4 W/m².
    for electric in [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]] {
        let input = Input {
            omega: 1.0,
            kx: 0.0,
            epsilon_i: C64::new(-3.0, -4.0),
            epsilon_t: C64::ONE,
            mu_i: 1.0,
            mu_t: 1.0,
            electric,
        };
        let result = solve(input).unwrap();
        complex_near(
            result.incident_normal_wavenumber_rad_m,
            C64::new(1.0, -2.0),
            2e-15,
        );
        for r in result.reflected_tangential_e_gain {
            complex_near(r, C64::new(0.5, -0.5), 2e-15);
        }
        near(result.incident_normal_power_w_m2, 0.5, 2e-15);
        near(result.reflected_normal_power_w_m2, -0.25, 2e-15);
        near(result.interference_normal_power_w_m2, 1.0, 2e-15);
        near(result.transmitted_normal_power_w_m2, 1.25, 3e-15);
        assert!(result.transmitted_normal_power_w_m2 > result.incident_normal_power_w_m2);
        boundary_and_power(input, result);
    }
}

#[test]
fn identical_absorbing_sources_are_reflection_free() {
    let input = Input {
        epsilon_t: data().epsilon_i,
        mu_t: data().mu_i,
        ..data()
    };
    let result = solve(input).unwrap();
    assert_eq!(result.reflected_tangential_e_gain, [C64::ZERO; 2]);
    assert_eq!(result.transmitted_tangential_e_gain, [C64::ONE; 2]);
    assert_eq!(result.reflected_normal_power_w_m2, 0.0);
    assert_eq!(result.interference_normal_power_w_m2, 0.0);
    boundary_and_power(input, result);
}

#[test]
fn lossless_reduction_preserves_the_original_interface_response() {
    for epsilon_t in [0.3, 2.25] {
        for mu_t in [0.8, 1.0] {
            for cosine in [0.3, 0.8, 1.0] {
                for (channel, polarization) in [(0, Polarization::Tm), (1, Polarization::Te)] {
                    let mut electric = [C64::ZERO; 2];
                    electric[channel] = C64::ONE;
                    let kx = ((1.0_f64 - cosine) * (1.0 + cosine)).sqrt();
                    let result = solve(Input {
                        omega: 1.0,
                        kx,
                        epsilon_i: C64::ONE,
                        epsilon_t: C64::from_re(epsilon_t),
                        mu_i: 1.0,
                        mu_t,
                        electric,
                    })
                    .unwrap();
                    let original = interface_response(
                        WaveMedium::new(C64::from_re(epsilon_t), mu_t, 1.0, 1.0).unwrap(),
                        1.0,
                        cosine,
                        polarization,
                    )
                    .unwrap();
                    complex_near(
                        result.reflected_tangential_e_gain[channel],
                        original.reflected_tangential_e,
                        2e-14,
                    );
                    complex_near(
                        result.transmitted_tangential_e_gain[channel],
                        original.transmitted_tangential_e,
                        2e-14,
                    );
                    assert_eq!(result.interference_normal_power_w_m2, 0.0);
                    near(
                        -result.reflected_normal_power_w_m2 / result.incident_normal_power_w_m2,
                        original.reflected_power_fraction,
                        2e-14,
                    );
                    near(
                        result.transmitted_normal_power_w_m2 / result.incident_normal_power_w_m2,
                        original.transmitted_power_fraction,
                        2e-14,
                    );
                }
            }
        }
    }
}

#[test]
fn critical_exit_retains_longitudinal_tm_and_te_fields_without_dividing_by_zero() {
    for electric in [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]] {
        let input = Input {
            omega: 1.0,
            kx: 1.0,
            epsilon_i: C64::new(2.0, -0.5),
            epsilon_t: C64::ONE,
            mu_i: 1.0,
            mu_t: 1.0,
            electric,
        };
        let result = solve(input).unwrap();
        assert_eq!(result.transmitted_normal_wavenumber_rad_m, C64::ZERO);
        assert_eq!(result.transmitted_normal_power_w_m2, 0.0);
        assert_eq!(result.reflected_tangential_e_gain, [-C64::ONE, C64::ONE]);
        if electric[0] != C64::ZERO {
            assert_eq!(result.transmitted.electric_v_m[0], C64::ZERO);
            assert!(result.transmitted.electric_v_m[2].abs() > 1.0);
            assert!(result.transmitted.magnetic_a_m[1].abs() > 1.0);
        } else {
            assert_eq!(result.transmitted.magnetic_a_m[0], C64::ZERO);
            complex_near(result.transmitted.electric_v_m[1], C64::from_re(2.0), 2e-15);
            complex_near(result.transmitted.magnetic_a_m[2], C64::from_re(2.0), 2e-15);
        }
        maxwell(result.transmitted, 1.0, C64::ZERO, 1.0, C64::ONE, 1.0);
        boundary_and_power(input, result);
    }
}

#[test]
fn evanescent_exit_keeps_zero_normal_flux_and_the_complex_reflection() {
    let input = Input {
        omega: 1.0,
        kx: 1.0,
        epsilon_i: C64::new(2.0, -0.1),
        epsilon_t: C64::from_re(0.25),
        mu_i: 1.0,
        mu_t: 1.0,
        ..data()
    };
    let result = solve(input).unwrap();
    assert_eq!(result.transmitted_normal_power_w_m2, 0.0);
    assert_eq!(result.transmitted_normal_wavenumber_rad_m.re, 0.0);
    assert!(result.transmitted_normal_wavenumber_rad_m.im < 0.0);
    assert!(result.interference_normal_power_w_m2.abs() > 1e-4);
    boundary_and_power(input, result);
}

#[test]
fn original_product_residual_survives_absolute_wave_number_admission() {
    // .625*1.6 rounds to 1, but its exact binary value is 1+2^-54.
    let input = Input {
        omega: 1.0,
        kx: 1.0,
        epsilon_t: C64::from_re(0.625),
        mu_t: 1.6,
        ..data()
    };
    assert_eq!(0.625 * 1.6, 1.0);
    let result = solve(input).unwrap();
    near(
        result.transmitted_normal_wavenumber_rad_m.re,
        2.0_f64.powi(-27),
        1e-22,
    );
    assert_eq!(result.transmitted_normal_wavenumber_rad_m.im, 0.0);
    assert!(result.transmitted_normal_power_w_m2 > 0.0);
}

#[test]
fn actual_amplitude_and_global_phase_scale_fields_and_power() {
    let input = data();
    let base = solve(input).unwrap();
    let factor = C64::new(0.0, 2.0);
    let scaled = solve(Input {
        electric: input.electric.map(|e| e * factor),
        ..input
    })
    .unwrap();
    assert_eq!(
        base.reflected_tangential_e_gain,
        scaled.reflected_tangential_e_gain
    );
    for (a, b) in [
        (base.incident, scaled.incident),
        (base.reflected, scaled.reflected),
        (base.transmitted, scaled.transmitted),
    ] {
        for c in 0..3 {
            complex_near(b.electric_v_m[c], a.electric_v_m[c] * factor, 2e-14);
            complex_near(b.magnetic_a_m[c], a.magnetic_a_m[c] * factor, 2e-14);
        }
    }
    for (a, b) in [
        (
            base.incident_normal_power_w_m2,
            scaled.incident_normal_power_w_m2,
        ),
        (
            base.reflected_normal_power_w_m2,
            scaled.reflected_normal_power_w_m2,
        ),
        (
            base.interference_normal_power_w_m2,
            scaled.interference_normal_power_w_m2,
        ),
        (
            base.transmitted_normal_power_w_m2,
            scaled.transmitted_normal_power_w_m2,
        ),
    ] {
        near(b, 4.0 * a, 2e-14);
    }
}

#[test]
fn frequency_and_transverse_scaling_preserve_the_boundary_solution() {
    let input = data();
    let base = solve(input).unwrap();
    for exponent in [-100, 100] {
        let scale = 2.0_f64.powi(exponent);
        let result = solve(Input {
            omega: input.omega * scale,
            kx: input.kx * scale,
            ..input
        })
        .unwrap();
        for c in 0..2 {
            complex_near(
                result.reflected_tangential_e_gain[c],
                base.reflected_tangential_e_gain[c],
                3e-14,
            );
        }
        complex_near(
            result.incident_normal_wavenumber_rad_m.scale(1.0 / scale),
            base.incident_normal_wavenumber_rad_m,
            3e-14,
        );
        near(
            result.incident_normal_power_w_m2,
            base.incident_normal_power_w_m2,
            3e-14,
        );
        near(
            result.transmitted_normal_power_w_m2,
            base.transmitted_normal_power_w_m2,
            3e-14,
        );
    }
}

#[test]
fn binary_references_avoid_overflow_and_underflow_of_the_unscaled_bulk_product() {
    // epsilon*mu is outside binary64, but k_z and the actual E/H and power
    // are representable. The independent root is sqrt(1-i) times the scale.
    let root = C64::new(1.0, -1.0).sqrt();
    for scale in [1e-300, 1e300] {
        let epsilon = C64::new(scale, -scale);
        let result = solve(Input {
            omega: 1.0,
            kx: 0.0,
            epsilon_i: epsilon,
            epsilon_t: epsilon,
            mu_i: scale,
            mu_t: scale,
            electric: [C64::ZERO, C64::ONE],
        })
        .unwrap();
        complex_near(
            divide_real(result.incident_normal_wavenumber_rad_m, scale),
            root,
            4e-15,
        );
        assert_eq!(result.reflected_tangential_e_gain, [C64::ZERO; 2]);
        near(result.incident_normal_power_w_m2, 0.5 * root.re, 4e-15);
        near(result.transmitted_normal_power_w_m2, 0.5 * root.re, 4e-15);
    }
}

#[test]
fn zero_illumination_and_refusals_keep_the_admitted_domain_explicit() {
    let zero = solve(Input {
        electric: [C64::ZERO; 2],
        ..data()
    })
    .unwrap();
    assert_eq!(zero.incident_normal_power_w_m2, 0.0);
    assert_eq!(zero.reflected_normal_power_w_m2, 0.0);
    assert_eq!(zero.interference_normal_power_w_m2, 0.0);
    assert_eq!(zero.transmitted_normal_power_w_m2, 0.0);
    for input in [
        Input {
            omega: 0.0,
            ..data()
        },
        Input {
            omega: f64::NAN,
            ..data()
        },
        Input {
            kx: f64::INFINITY,
            ..data()
        },
        Input {
            mu_i: 0.0,
            ..data()
        },
        Input {
            epsilon_i: C64::new(1.0, 0.1),
            ..data()
        },
        Input {
            epsilon_t: C64::ZERO,
            ..data()
        },
        Input {
            electric: [C64::new(f64::NAN, 0.0), C64::ZERO],
            ..data()
        },
        Input {
            omega: 1.0,
            kx: 1.0,
            epsilon_i: C64::ONE,
            mu_i: 1.0,
            ..data()
        },
        Input {
            omega: 1.0,
            kx: 2.0,
            epsilon_i: C64::ONE,
            mu_i: 1.0,
            ..data()
        },
        Input {
            electric: [C64::from_re(1e-300), C64::ZERO],
            ..data()
        },
    ] {
        assert!(solve(input).is_err());
    }
}
