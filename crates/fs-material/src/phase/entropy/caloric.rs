//! Positive, dimensionless moments of an original affine Cp(T) interval.

/// All kernels use v=(T-Ta)/(Tb-Ta), with y=(h(T)-ha)/(hb-ha).
/// Keeping the kernels positive avoids entropy-quotient cancellation in the
/// finite thermal conjugate, including enthalpy steps with equal rounded T.
pub(super) struct CaloricEntropySegment {
    pub entropy: f64,
    pub loss: f64,
    pub reciprocal_mean_capacity: f64,
    pub zeroth_kernel: f64,
    pub lower_kernel: f64,
    pub upper_kernel: f64,
}

/// A clipped interval measured in enthalpy. Recover its small temperature
/// increment from dh/Cp_mean instead of subtracting rounded temperatures.
pub(super) fn segment(dh: f64, ta: f64, capacities: [f64; 2]) -> Option<CaloricEntropySegment> {
    let mean_capacity = 0.5 * capacities[0] + 0.5 * capacities[1];
    let x = (dh / mean_capacity) / ta;
    let (ratio, loss, zeroth, lower, upper) = moments(x, capacities)?;
    let entropy = (dh / ta) * ratio;
    let reciprocal_mean_capacity = mean_capacity.recip();
    if !entropy.is_finite() || entropy <= 0.0 {
        return None;
    }
    Some(CaloricEntropySegment {
        entropy,
        loss,
        reciprocal_mean_capacity,
        zeroth_kernel: zeroth,
        lower_kernel: lower,
        upper_kernel: upper,
    })
}

/// Exact integral of Cp(T)/T over ascending, finite positive temperatures.
/// Endpoint capacities are the original affine law, not an apparent Cp.
pub(super) fn affine_cp_entropy(ta: f64, tb: f64, capacities: [f64; 2]) -> Option<f64> {
    if !ta.is_finite() || !tb.is_finite() || ta <= 0.0 || tb < ta {
        return None;
    }
    let x = (tb - ta) / ta;
    let (ratio, _, _, _, _) = moments(x, capacities)?;
    if ta == tb {
        return Some(0.0);
    }
    let mean_capacity = 0.5 * capacities[0] + 0.5 * capacities[1];
    let entropy = (mean_capacity * ratio) * x;
    (entropy.is_finite() && entropy > 0.0).then_some(entropy)
}

/// Entropy divided by dh/Ta; its cancellation-free complement; and integrals
/// of 1, 1-y, y against 1/(1+x*v)^2. Cp endpoint weights are positive, even
/// when its slope is negative or its intercept at zero temperature is negative.
fn moments(x: f64, capacities: [f64; 2]) -> Option<(f64, f64, f64, f64, f64)> {
    if !x.is_finite() || x < 0.0 || capacities.iter().any(|cp| !cp.is_finite() || *cp <= 0.0) {
        return None;
    }
    let scale = capacities[0].max(capacities[1]);
    let scaled_a = capacities[0] / scale;
    let scaled_b = capacities[1] / scale;
    let denominator = scaled_a + scaled_b;
    let weight_a = scaled_a / denominator;
    let weight_b = scaled_b / denominator;
    let (basis_a, basis_b, loss, first, second) = if x <= 0.125 {
        // Every term is an exact monomial integral. Twenty terms at x<=1/8
        // leave less than 1e-18 absolute truncation error without tiny
        // differences of logarithms, temperatures, or entropy quotients.
        let mut power = 1.0;
        let mut basis_a = 0.0;
        let mut basis_b = 0.0;
        let mut loss_over_x = 0.0;
        let mut first = 0.0;
        let mut second = 0.0;
        for n in 0..20 {
            let n1 = f64::from(n + 1);
            let n2 = f64::from(n + 2);
            let n3 = f64::from(n + 3);
            basis_a += power / (n1 * n2);
            basis_b += power / n2;
            loss_over_x += power * (weight_a / (n2 * n3) + weight_b / n3);
            first += power * n1 / n2;
            second += power * n1 / n3;
            power *= -x;
        }
        (basis_a, basis_b, 2.0 * x * loss_over_x, first, second)
    } else {
        let rounded = 1.0 + x;
        let logarithm = fs_math::det::ln(rounded) * (x / (rounded - 1.0));
        let average = logarithm / x;
        let basis_b = (1.0 - average) / x;
        let basis_a = average - basis_b;
        let first = ((logarithm - x / (1.0 + x)) / x) / x;
        let second = (basis_b - first) / x;
        let ratio = 2.0 * (weight_a * basis_a + weight_b * basis_b);
        (basis_a, basis_b, 1.0 - ratio, first, second)
    };
    let ratio = 2.0 * (weight_a * basis_a + weight_b * basis_b);
    let zeroth = (1.0 + x).recip();
    let upper = weight_a * (2.0 * first - second) + weight_b * second;
    // For x>=0 this difference retains at least a third of zeroth, so it
    // cannot suffer the cancellation that the enthalpy quotient would.
    let lower = zeroth - upper;
    // Scalar entropy need not refuse solely because a derivative kernel is
    // unrepresentably small. The conjugate consumer validates its required
    // positive kernels and scaled derivatives independently.
    if !ratio.is_finite() || !loss.is_finite() || ratio <= 0.0 || loss < 0.0 {
        return None;
    }
    Some((ratio, loss, zeroth, lower, upper))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn affine_capacity_entropy_matches_original_primitive_in_both_slope_directions() {
        for capacities in [[1000.0, 4000.0], [4000.0, 1000.0], [2000.0, 2000.0]] {
            let ta = 300.0_f64;
            let tb = 900.0_f64;
            let slope = (capacities[1] - capacities[0]) / (tb - ta);
            let expected = (capacities[0] - slope * ta) * (tb / ta).ln() + slope * (tb - ta);
            let actual = affine_cp_entropy(ta, tb, capacities).unwrap();
            assert!((actual - expected).abs() < 2.0e-12);
        }
    }

    #[test]
    fn kernels_retain_small_step_heat_and_nonzero_derivatives() {
        let result = segment(1e-300, 300.0, [2000.0, 2000.0]).unwrap();
        assert_eq!(result.entropy, 1e-300 / 300.0);
        assert_eq!(result.reciprocal_mean_capacity, 0.0005);
        assert_eq!(result.lower_kernel, 0.5);
        assert_eq!(result.upper_kernel, 0.5);
        let ta = 300.0_f64;
        let tb = ta.next_up();
        let entropy = affine_cp_entropy(ta, tb, [1000.0, 4000.0]).unwrap();
        let heat = 2500.0 * (tb - ta);
        assert!(entropy >= heat / tb && entropy <= heat / ta);
    }

    #[test]
    fn capacity_proportional_to_temperature_has_exact_arithmetic_mean_and_jacobians() {
        // Cp(T)=a*T makes ds=a*dT and theta=(Ta+Tb)/2. This independently
        // checks both nonconstant-Cp kernels without differencing a quotient.
        for x in [1e-15, 0.001, 0.125, 0.126, 1.0, 10.0, 1e5] {
            let ta = 300.0;
            let delta_temperature = ta * x;
            let capacities = [1000.0, 1000.0 * (1.0 + x)];
            let heat = (0.5 * capacities[0] + 0.5 * capacities[1]) * delta_temperature;
            let result = segment(heat, ta, capacities).unwrap();
            let theta = ta + heat * result.loss / result.entropy;
            let scale = theta / ta;
            let lower = result.reciprocal_mean_capacity * result.lower_kernel * scale * scale;
            let upper = result.reciprocal_mean_capacity * result.upper_kernel * scale * scale;
            assert!((theta / (ta + 0.5 * delta_temperature) - 1.0).abs() < 1e-13);
            assert!(
                (lower / (0.5 / capacities[0]) - 1.0).abs() < 1e-13,
                "x={x}, lower={lower}"
            );
            assert!(
                (upper / (0.5 / capacities[1]) - 1.0).abs() < 1e-13,
                "x={x}, upper={upper}"
            );
            assert!(result.zeroth_kernel > 0.0);
        }
    }

    #[test]
    fn finite_scalar_entropy_does_not_require_representable_jacobian_moments() {
        let entropy = affine_cp_entropy(300.0, 1e200, [1000.0, 1000.0]).unwrap();
        let expected = 1000.0 * (1e200_f64 / 300.0).ln();
        assert!((entropy / expected - 1.0).abs() < 1e-14);
    }

    #[test]
    fn positive_source_weights_admit_extreme_slopes_and_refuse_invalid_support() {
        for x in [0.0, 1e-15, 0.124, 0.126, 1.0, 1000.0, 1e100] {
            for capacities in [[1e-100, 1e100], [1e100, 1e-100]] {
                let (ratio, loss, zero, lower, upper) = moments(x, capacities).unwrap();
                assert!(ratio > 0.0 && ratio <= 1.0);
                assert!(loss >= 0.0 && loss < 1.0 || loss == 1.0);
                assert!((ratio + loss - 1.0).abs() <= 5.0 * f64::EPSILON);
                assert!((lower + upper - zero).abs() <= f64::EPSILON * zero);
            }
        }
        for (ta, tb, capacities) in [
            (0.0, 100.0, [1000.0, 1000.0]),
            (300.0, 200.0, [1000.0, 1000.0]),
            (300.0, f64::INFINITY, [1000.0, 1000.0]),
            (300.0, 400.0, [0.0, 1000.0]),
        ] {
            assert!(affine_cp_entropy(ta, tb, capacities).is_none());
        }
    }
}
