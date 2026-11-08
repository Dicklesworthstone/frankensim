//! Analytic blackbody bounds outside a declared wavelength band.
//!
//! A passive equilibrium surface has emissivity in [0, 1]. Its unknown
//! out-of-band power is therefore bounded by blackbody power there. These
//! bounds do not continue a measured spectrum or estimate in-band quadrature
//! error. Their logarithmic evaluation keeps emitting area in the calculation,
//! including when exitance alone would overflow.

use fs_math::det;

pub(super) const PLANCK_J_S: f64 = 6.626_070_15e-34;
pub(super) const LIGHT_M_S: f64 = 299_792_458.0;
pub(super) const BOLTZMANN_J_K: f64 = 1.380_649e-23;
pub(super) const C1_NM: f64 = 2.0 * PLANCK_J_S * LIGHT_M_S * LIGHT_M_S * 1.0e36;
pub(super) const C2_NM_K: f64 = PLANCK_J_S * LIGHT_M_S / BOLTZMANN_J_K * 1.0e9;

/// Conservative omitted hemispherical power [W] for one emitting face.
///
/// With x = hc/(lambda k T), the dimensionless Planck integrand is
/// x^3/(exp(x)-1). Above a lower x=a, replacing its denominator by
/// exp(x)*(1-exp(-a)) gives the integrated Wien majorant
/// exp(-a)*(a^3+3a^2+6a+6)/(1-exp(-a)). At long wavelengths,
/// exp(x)-1 >= x gives the integrated Rayleigh--Jeans majorant x^3/3.
/// Their sum is capped by the complete Planck integral pi^4/15.
///
/// All constants are the SAME exact-SI inputs as the pointwise source. The
/// logarithmic operation tree is padded outward beyond the strict elementary
/// functions' declared ULP budgets, including their large-argument arithmetic.
/// This is an analytic model bound with floating-point padding, not an interval
/// certificate for a measured material. Nonzero underflow retains one positive
/// subnormal; invalid inputs or an unrepresentable bound return None.
pub(crate) fn omitted_blackbody_power_upper_bound_w(
    temperature_k: f64,
    wavelength_bounds_nm: (f64, f64),
    area_m2: f64,
) -> Option<f64> {
    let (lower, upper) = wavelength_bounds_nm;
    if !temperature_k.is_finite()
        || temperature_k <= 0.0
        || !lower.is_finite()
        || lower <= 0.0
        || !upper.is_finite()
        || upper <= lower
        || !area_m2.is_finite()
        || area_m2 < 0.0
    {
        return None;
    }
    if area_m2 == 0.0 {
        return Some(0.0);
    }
    let log_t = det::ln(temperature_k);
    let log_area = det::ln(area_m2);
    let log_lower = det::ln(lower);
    let log_upper = det::ln(upper);
    let log_c1 = det::ln(C1_NM);
    let log_c2 = det::ln(C2_NM_K);
    let log_pi = det::ln(core::f64::consts::PI);

    // A*pi*C1*(T/C2)^4 multiplies the dimensionless Planck integral.
    let log_scale = log_area + log_pi + log_c1 + 4.0 * (log_t - log_c2);
    let scale_magnitude =
        log_area.abs() + log_pi.abs() + log_c1.abs() + 4.0 * (log_t.abs() + log_c2.abs());
    let total = padded_log(
        log_scale + 4.0 * log_pi - det::ln(15.0),
        scale_magnitude + 4.0 * log_pi.abs() + det::ln(15.0),
    );

    // Integral from lambda=upper to infinity of the Rayleigh--Jeans law.
    let long = padded_log(
        log_area + log_pi + log_c1 - log_c2 + log_t - det::ln(3.0) - 3.0 * log_upper,
        log_area.abs()
            + log_pi.abs()
            + log_c1.abs()
            + log_c2.abs()
            + log_t.abs()
            + det::ln(3.0)
            + 3.0 * log_upper.abs(),
    );

    // The short-wave majorant decreases with x. Round its input DOWN so the
    // bound cannot shrink because the dimensionless endpoint rounded upward.
    let log_x = log_c2 - log_lower - log_t;
    let log_x_lower = log_x - padding(log_c2.abs() + log_lower.abs() + log_t.abs());
    let mut x = det::exp(log_x_lower);
    for _ in 0..=det::EXP_ULP_BUDGET {
        x = x.next_down().max(0.0);
    }
    let short = if x < 1.0 {
        // The complete power is a sharper bound than this Wien majorant here.
        total
    } else if x.is_infinite() {
        f64::NEG_INFINITY
    } else {
        let inverse = 1.0 / x;
        let log_polynomial =
            3.0 * det::ln(x) + det::ln(1.0 + inverse * (3.0 + inverse * (6.0 + 6.0 * inverse)));
        let log_denominator = det::ln(-det::expm1(-x));
        padded_log(
            log_scale - x + log_polynomial - log_denominator,
            scale_magnitude + x + log_polynomial.abs() + log_denominator.abs(),
        )
    };
    let combined = log_add_upper(short, long);
    let log_bound = total.min(combined);
    let mut bound = det::exp(log_bound);
    for _ in 0..=det::EXP_ULP_BUDGET {
        bound = bound.next_up();
    }
    bound.is_finite().then_some(bound)
}

fn padding(magnitude: f64) -> f64 {
    64.0 * f64::EPSILON * (1.0 + magnitude)
}

fn padded_log(value: f64, magnitude: f64) -> f64 {
    value + padding(magnitude)
}

fn log_add_upper(a: f64, b: f64) -> f64 {
    let high = a.max(b);
    let low = a.min(b);
    if low == f64::NEG_INFINITY {
        high
    } else {
        padded_log(
            high + det::ln(1.0 + det::exp(low - high)),
            high.abs() + (low - high).abs(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Independent Simpson integration with the platform elementary-function
    // oracle, only in tests. x=80 leaves less than 1e-27 dimensionless power.
    fn integral(a: f64, b: f64) -> f64 {
        let n = 20_000;
        let dx = (b - a) / f64::from(n);
        let value = |x: f64| {
            if x == 0.0 {
                0.0
            } else {
                x * x * x / x.exp_m1()
            }
        };
        let mut sum = value(a) + value(b);
        for i in 1..n {
            sum += if i % 2 == 0 { 2.0 } else { 4.0 } * value(a + f64::from(i) * dx);
        }
        sum * dx / 3.0
    }

    #[test]
    fn analytic_bounds_cover_independent_planck_tail_integrals() {
        for temperature in [300.0, 2500.0, 6000.0] {
            for band in [(380.0, 780.0), (100.0, 30_000.0), (10.0, 1.0e7)] {
                let prefactor = core::f64::consts::PI * C1_NM * (temperature / C2_NM_K).powi(4);
                let long_x = (C2_NM_K / (band.1 * temperature)).min(80.0);
                let short_x = (C2_NM_K / (band.0 * temperature)).min(80.0);
                let reference = prefactor * (integral(0.0, long_x) + integral(short_x, 80.0));
                let bound = omitted_blackbody_power_upper_bound_w(temperature, band, 1.0).unwrap();
                assert!(
                    bound >= reference * (1.0 - 2.0e-12),
                    "{temperature} {band:?}: {bound} < {reference}"
                );
                let complete = prefactor * core::f64::consts::PI.powi(4) / 15.0;
                assert!(bound <= complete * (1.0 + 1.0e-11));
            }
        }
    }

    #[test]
    fn spectral_support_area_and_temperature_scaling_remain_physical() {
        let narrow = omitted_blackbody_power_upper_bound_w(2600.0, (380.0, 780.0), 0.02).unwrap();
        let wide = omitted_blackbody_power_upper_bound_w(2600.0, (100.0, 100_000.0), 0.02).unwrap();
        assert!(wide < 0.0001 * narrow);
        let doubled_area =
            omitted_blackbody_power_upper_bound_w(2600.0, (100.0, 100_000.0), 0.04).unwrap();
        assert!((doubled_area / (2.0 * wide) - 1.0).abs() < 1.0e-12);
        let scaled =
            omitted_blackbody_power_upper_bound_w(5200.0, (50.0, 50_000.0), 0.02 / 16.0).unwrap();
        assert!((scaled / wide - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn area_safe_extremes_and_invalid_inputs_have_explicit_results() {
        assert!(
            omitted_blackbody_power_upper_bound_w(1.0e100, (1.0e-95, 1.0e-91), 1.0e-200).unwrap()
                > 0.0
        );
        assert!(
            omitted_blackbody_power_upper_bound_w(1.0e-100, (1.0, 2.0), 1.0e-100).unwrap() > 0.0
        );
        assert_eq!(
            omitted_blackbody_power_upper_bound_w(2600.0, (1.0, 2.0), 0.0),
            Some(0.0)
        );
        for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(omitted_blackbody_power_upper_bound_w(temperature, (1.0, 2.0), 1.0).is_none());
        }
        for band in [(0.0, 2.0), (2.0, 1.0), (1.0, f64::INFINITY)] {
            assert!(omitted_blackbody_power_upper_bound_w(2600.0, band, 1.0).is_none());
        }
        assert!(omitted_blackbody_power_upper_bound_w(2600.0, (1.0, 2.0), -1.0).is_none());
        assert!(omitted_blackbody_power_upper_bound_w(1.0e100, (1.0, 2.0), 1.0).is_none());
    }
}
