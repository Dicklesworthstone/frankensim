//! Rest-safe evaluation of Churchill's explicitly selected pipe-friction law.

use fs_math::det;

/// `(f_D Re / 64) - 1`, retaining small excesses above Poiseuille resistance.
/// Inputs have already been admitted to Re in [0, 1e8] and epsilon/D in [0, .05].
///
/// Churchill (1977), equations 5 and 18: the paper's friction factor is f_D/8.
/// <https://files.engineering.com/files/85c0f3a6-a102-4a22-9d35-f15858c0dd2b/CEM_-_Friction-factor_equation_%281977%29.pdf>
/// With A=a^16, B=b^16, and q=Re/[8(A+B)^(1/8)], the original relation is
/// f_D Re = 64 (1+q^12)^(1/12). This avoids the divergent f_D at rest and the
/// catastrophic subtraction f_D - 64/Re when forming quadratic resistance.
pub(crate) fn excess_over_laminar(reynolds: f64, relative_roughness: f64) -> f64 {
    if reynolds <= 1.0 {
        // q <= Re^3/(8*37530^2), so the relative excess here is below 2e-122.
        // Its omission cannot change represented total resistance. This also
        // covers subnormal Re without evaluating an overflowing 37530/Re.
        return 0.0;
    }
    let argument = det::pow(7.0 / reynolds, 0.9) + 0.27 * relative_roughness;
    let a = -2.457 * det::ln(argument);
    let b = 37_530.0 / reynolds;
    let sum = pow16(a) + pow16(b);
    let q = reynolds / (8.0 * det::sqrt(det::sqrt(det::sqrt(sum))));
    let q2 = q * q;
    let q4 = q2 * q2;
    let q12 = (q4 * q4) * q4;
    // Correct ln(1+x) for rounding in 1+x; expm1 then retains the small
    // difference from unity. Over the admitted range all powers are finite.
    let one_plus = 1.0 + q12;
    let log1p = if one_plus == 1.0 {
        q12
    } else {
        det::ln(one_plus) * (q12 / (one_plus - 1.0))
    };
    det::expm1(log1p / 12.0)
}

fn pow16(x: f64) -> f64 {
    let x2 = x * x;
    let x4 = x2 * x2;
    let x8 = x4 * x4;
    x8 * x8
}

#[cfg(test)]
mod tests {
    use super::excess_over_laminar;

    #[test]
    fn churchill_g1_original_equation_and_darcy_normalization() {
        // Independent unscaled published equation, using platform functions.
        for roughness in [0.0, 1e-8, 0.001, 0.01, 0.05] {
            for re in [
                1.0_f64, 100.0, 1000.0, 2100.0, 2300.0, 3000.0, 4000.0, 1e4, 1e5, 1e6, 1e8,
            ] {
                let a = (-2.457 * ((7.0 / re).powf(0.9) + 0.27 * roughness).ln()).powi(16);
                let b = (37_530.0 / re).powi(16);
                let expected =
                    8.0 * re * ((8.0 / re).powi(12) + (a + b).powf(-1.5)).powf(1.0 / 12.0);
                let actual = 64.0 * (1.0 + excess_over_laminar(re, roughness));
                assert!(
                    (actual / expected - 1.0).abs() < 3e-13,
                    "Re={re} epsilon/D={roughness}: {actual} != {expected}"
                );
                assert!(actual >= 64.0);
            }
        }
    }

    #[test]
    fn churchill_g1_rest_limit_and_small_positive_quadratic_resistance() {
        for re in [0.0, f64::from_bits(1), 1e-300, 1e-100, 0.5, 1.0] {
            for roughness in [0.0, 0.05] {
                assert_eq!(64.0 * (1.0 + excess_over_laminar(re, roughness)), 64.0);
            }
        }
        // These small corrections would disappear if computed as f_D - 64/Re.
        for re in [10.0, 100.0, 500.0] {
            let excess = excess_over_laminar(re, 0.001);
            assert!(excess > 0.0 && excess < f64::EPSILON);
            assert!((32.0 * excess / re).is_finite());
        }
    }

    #[test]
    fn churchill_g1_continuous_transition_and_fully_rough_limit() {
        for re in [2100.0_f64, 2300.0, 3000.0, 4000.0] {
            let left = 1.0 + excess_over_laminar(re * (1.0 - 1e-8), 0.001);
            let right = 1.0 + excess_over_laminar(re * (1.0 + 1e-8), 0.001);
            assert!((right / left - 1.0).abs() < 1e-6);
        }
        let mut previous = 0.0;
        for roughness in [0.0, 0.001, 0.01, 0.05] {
            let factor = 64.0 * (1.0 + excess_over_laminar(1e8, roughness)) / 1e8;
            assert!(factor > previous);
            previous = factor;
        }
        let fully_rough = 8.0 / (2.457 * (1.0_f64 / (0.27 * 0.05)).ln()).powi(2);
        assert!((previous / fully_rough - 1.0).abs() < 2e-5);
    }
}
