use super::*;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}");
}

#[test]
fn original_mode_jets_match_independent_complex_rational_derivatives() {
    for force in [3e8_f64, -3e8] {
        for rate in [2e5_f64, 2e7] {
            for omega in [1e5_f64, 4e7, 1e10] {
                let sigma = C64::from_re(force * force) / C64::new(rate, omega);
                let chi = C64::new(sigma.im / omega, -sigma.re / omega);
                let r = C64::new(0.0, omega) / C64::new(rate, omega);
                let first = -(chi * (C64::ONE + r));
                let second = (chi * (C64::ONE + r + r * r)).scale(2.0);
                let actual = mode_derivatives(force, rate, omega).unwrap();
                let expected = [first.re, first.im, second.re, second.im,
                    -(chi * r).im, (chi * r * r).scale(2.0).im];
                for i in 0..6 {
                    near(actual[i], expected[i], 3e-10);
                }
            }
        }
    }
}

#[test]
fn weak_loss_moments_and_inertial_real_slopes_keep_representable_tails() {
    let low = mode_derivatives(1.0, 1.0, 1e-100).unwrap();
    near(low[0], 2e-200, 4e-15);
    near(low[2], 2e-200, 4e-15);
    near(low[4], 2e-100, 4e-15);
    near(low[5], 2e-100, 4e-15);
    let high = mode_derivatives(1.0, 1.0, 1e100).unwrap();
    near(high[0], 2e-200, 4e-15);
    near(high[1], 3e-300, 4e-15);
    near(high[2], -6e-200, 4e-15);
    near(high[3], -12e-300, 4e-15);
    near(high[4], 2e-300, 4e-15);
    near(high[5], -6e-300, 4e-15);
}

#[test]
fn dark_modes_contribute_nothing_and_unrepresentable_modes_refuse() {
    assert_eq!(mode_derivatives(0.0, 1.0, 1.0).unwrap(), [0.0; 6]);
    assert!(mode_derivatives(1e308, 1.0, 1.0).is_err());
}
