//! Analytic averages for one original ideal-gas species in constant gravity.
//!
//! The caller supplies its source-owned datum density and `beta=M*g/(R*T)`.
//! This leaf neither changes composition nor supplies molecular dynamics.

use fs_math::det;

/// Average `datum_density * exp(beta*x)` over an offset interval. A zero-width
/// interval evaluates a physical face. Authored zeros stay zero; an active
/// species that cannot be represented positively is refused, never removed.
pub(crate) fn exponential_average(
    datum_density: f64,
    beta: f64,
    left: f64,
    right: f64,
) -> Option<f64> {
    if !datum_density.is_finite()
        || datum_density < 0.0
        || !beta.is_finite()
        || !left.is_finite()
        || !right.is_finite()
        || left > right
    {
        return None;
    }
    if datum_density == 0.0 || beta == 0.0 {
        return Some(datum_density);
    }
    let span = beta.abs() * (right - left);
    let exponent = beta * if beta > 0.0 { right } else { left };
    if !span.is_finite() || !exponent.is_finite() {
        return None;
    }
    // Anchoring at the denser end avoids sinh overflow. expm1 retains the
    // correct small-gravity limit without subtracting nearly equal values.
    let factor = if span == 0.0 {
        1.0
    } else {
        -det::expm1(-span) / span
    };
    let direct = datum_density * det::exp(exponent) * factor;
    let value = if direct.is_finite() && direct > 0.0 {
        direct
    } else {
        // A large exponent and a small datum can still have a finite product.
        // This path also avoids an overflowing endpoint when its average fits.
        det::exp(det::ln(datum_density) + exponent + det::ln(factor))
    };
    (value.is_finite() && value > 0.0).then_some(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64, relative: f64) {
        assert!(
            (actual - expected).abs() <= relative * expected.abs(),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn analytic_species_averages_obey_gravity_balance_and_partition_additivity() {
        for beta in [-40.0, -2.0, 0.2, 3.0, 40.0] {
            let (a, split, b, rho0) = (-0.2, 0.07, 0.4, 1.7);
            let whole = exponential_average(rho0, beta, a, b).unwrap();
            // Independent antiderivative using the platform exponential.
            let exact = rho0 * ((beta * b).exp() - (beta * a).exp()) / (beta * (b - a));
            near(whole, exact, 8e-15);
            let left = exponential_average(rho0, beta, a, split).unwrap();
            let right = exponential_average(rho0, beta, split, b).unwrap();
            near(
                whole * (b - a),
                left * (split - a) + right * (b - split),
                8e-15,
            );
            let lower = exponential_average(rho0, beta, a, a).unwrap();
            let upper = exponential_average(rho0, beta, b, b).unwrap();
            near(beta * whole * (b - a), upper - lower, 8e-15);
        }
    }

    #[test]
    fn small_gravity_and_exact_zeros_have_no_cancellation_floor() {
        for beta in [-1e-14, -1e-8, 0.0, 1e-8, 1e-14] {
            let average = exponential_average(1.7, beta, -0.25, 0.25).unwrap();
            let z = 0.25 * beta;
            near(average, 1.7 * (1.0 + z * z / 6.0), 3e-16);
        }
        assert_eq!(exponential_average(1.7, 0.0, -2.0, 8.0), Some(1.7));
        assert_eq!(exponential_average(0.0, 1e300, -1e300, 1e300), Some(0.0));
        assert_eq!(exponential_average(1.0, -1000.0, 1.0, 1.0), None);
        assert_eq!(exponential_average(1.0, 1000.0, 1.0, 1.0), None);
        assert_eq!(exponential_average(-1.0, 1.0, 0.0, 1.0), None);
    }

    #[test]
    fn representable_products_survive_intermediate_exponential_overflow() {
        near(
            exponential_average(1e-300, 720.0, 1.0, 1.0).unwrap(),
            (720.0 + (1e-300_f64).ln()).exp(),
            2e-13,
        );
        near(
            exponential_average(1e300, -750.0, 1.0, 1.0).unwrap(),
            (-750.0 + (1e300_f64).ln()).exp(),
            2e-13,
        );
        let avg = exponential_average(1.0, 714.0, 0.0, 1.0).unwrap();
        near(avg, (714.0 - 714.0_f64.ln()).exp(), 2e-13);
    }
}
