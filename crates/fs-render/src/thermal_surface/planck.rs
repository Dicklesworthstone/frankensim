//! Stable differences and temperature slopes of the original SI Planck law.
//!
//! Thermal exchange needs B(Ta)-B(Tb), not two separately rounded gross powers.
//! The operations below share the radiance owner's constants/log expression.
//! Logarithmic final conversion also retains representable slopes when a tiny
//! radiance difference underflows, and finite differences of overflowing gross
//! radiances. These are scalar source evaluations, not spectral quadrature or
//! emissivity assumptions; optical support stays with the caller's material.

use super::{
    BOLTZMANN_J_K, C2_NM_K, LIGHT_M_S, ThermalSurfaceError, checked_planck_exp, planck_log_radiance,
};
use fs_math::det;

/// Signed B(first)-B(second), in W/(m² sr nm), evaluated without cancellation at
/// nearly equal temperatures. Equal supported temperatures give exactly zero;
/// swapping unequal inputs exactly negates the result, including signed zero.
/// The admitted domain is the original positive finite kelvin/wavelength law.
pub(crate) fn planck_radiance_difference_w_m2_sr_nm(
    first_k: f64,
    second_k: f64,
    wavelength_nm: f64,
) -> Result<f64, ThermalSurfaceError> {
    validate(first_k, second_k, wavelength_nm)?;
    if first_k == second_k {
        return Ok(0.0);
    }
    let (hot, cold, sign) = ordered(first_k, second_k);
    let (difference, _) = difference_logs(hot, cold, wavelength_nm);
    Ok(sign * checked_planck_exp(difference)?)
}

/// Positive secant [B(first)-B(second)]/(first-second), in W/(m² sr nm K).
/// Equality returns the analytic tangent. No rounded radiance difference is
/// divided by a tiny ΔT; an underflowed difference can have a finite secant.
pub(crate) fn planck_radiance_temperature_secant_w_m2_sr_nm_k(
    first_k: f64,
    second_k: f64,
    wavelength_nm: f64,
) -> Result<f64, ThermalSurfaceError> {
    validate(first_k, second_k, wavelength_nm)?;
    if first_k == second_k {
        return planck_radiance_temperature_derivative_w_m2_sr_nm_k(first_k, wavelength_nm);
    }
    let (hot, cold, _) = ordered(first_k, second_k);
    let (_, secant) = difference_logs(hot, cold, wavelength_nm);
    checked_planck_exp(secant)
}

/// Net radiative power and its positive temperature secant, after applying a
/// nonnegative exchange coefficient (m² sr) and wavelength quadrature weight
/// (nm). `None` is an explicit vacuum bath with B=0; every supplied temperature
/// remains in the original positive finite domain. Vacuum secant is B(T)/T.
///
/// Both rates are evaluated together, adding the two log scale factors BEFORE
/// rounding to watts. Neither K*weight nor unscaled ΔB needs to be representable.
/// This prevents an underflowed intermediate from silently erasing real heat.
/// Equal temperatures have zero power and the scaled analytic tangent. Zero
/// coefficient or quadrature weight gives zero after domain validation.
pub(crate) fn planck_radiative_exchange_w_and_w_k(
    first_k: f64,
    second_k: Option<f64>,
    wavelength_nm: f64,
    coefficient_m2_sr: f64,
    quadrature_weight_nm: f64,
) -> Result<(f64, f64), ThermalSurfaceError> {
    validate(first_k, second_k.unwrap_or(first_k), wavelength_nm)?;
    if [coefficient_m2_sr, quadrature_weight_nm]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
    {
        return Err(ThermalSurfaceError::ExchangeCoefficient);
    }
    if coefficient_m2_sr == 0.0 || quadrature_weight_nm == 0.0 {
        return Ok((0.0, 0.0));
    }
    let scale = det::ln(coefficient_m2_sr) + det::ln(quadrature_weight_nm);
    let Some(second) = second_k else {
        let source = planck_log_radiance(first_k, wavelength_nm, exponent(first_k, wavelength_nm));
        return check_rates(
            checked_planck_exp(source + scale)?,
            checked_planck_exp(source - det::ln(first_k) + scale)?,
        );
    };
    if first_k == second {
        return Ok((
            0.0,
            checked_planck_exp(tangent_log(first_k, wavelength_nm) + scale)?,
        ));
    }
    let (hot, cold, sign) = ordered(first_k, second);
    let (difference, secant) = difference_logs(hot, cold, wavelength_nm);
    check_rates(
        sign * checked_planck_exp(difference + scale)?,
        checked_planck_exp(secant + scale)?,
    )
}

fn check_rates(power: f64, secant: f64) -> Result<(f64, f64), ThermalSurfaceError> {
    // A native caloric slope or small cell mass can amplify a subnormal
    // conductance. Publishing nonzero heat with a rounded-zero secant would
    // silently remove the only corresponding thermal stability restriction.
    if power != 0.0 && secant == 0.0 {
        Err(ThermalSurfaceError::Unrepresentable)
    } else {
        Ok((power, secant))
    }
}

/// Analytic dB/dT, in W/(m² sr nm K), from the same Planck law. The dimensionless
/// correction to the Rayleigh–Jeans slope is x² exp(-x)/(1-exp(-x))². Its
/// small-x log series avoids cancellation of the two first-order terms.
pub(crate) fn planck_radiance_temperature_derivative_w_m2_sr_nm_k(
    temperature_k: f64,
    wavelength_nm: f64,
) -> Result<f64, ThermalSurfaceError> {
    validate(temperature_k, temperature_k, wavelength_nm)?;
    checked_planck_exp(tangent_log(temperature_k, wavelength_nm))
}

fn tangent_log(temperature_k: f64, wavelength_nm: f64) -> f64 {
    let x = exponent(temperature_k, wavelength_nm);
    let correction = if x < 1.0e-3 {
        let x2 = x * x;
        // log[(x/(2 sinh(x/2)))²]; next term is O(x^8), below binary64
        // resolution throughout this branch, including x underflow to zero.
        x2 * (-1.0 / 12.0 + x2 * (1.0 / 1440.0 - x2 / 90_720.0))
    } else if x.is_infinite() {
        f64::NEG_INFINITY
    } else {
        2.0 * det::ln(x) - x - 2.0 * log_removed(x)
    };
    log_rayleigh_jeans_slope(wavelength_nm) + correction
}

fn validate(first: f64, second: f64, wavelength: f64) -> Result<(), ThermalSurfaceError> {
    if !first.is_finite() || first <= 0.0 || !second.is_finite() || second <= 0.0 {
        return Err(ThermalSurfaceError::Temperature);
    }
    if !wavelength.is_finite() || wavelength <= 0.0 {
        return Err(ThermalSurfaceError::Wavelength);
    }
    Ok(())
}

fn ordered(first: f64, second: f64) -> (f64, f64, f64) {
    if first > second {
        (first, second, 1.0)
    } else {
        (second, first, -1.0)
    }
}

fn exponent(temperature: f64, wavelength: f64) -> f64 {
    let numerator = C2_NM_K / wavelength;
    if numerator.is_finite() {
        numerator / temperature
    } else {
        // Avoid an overflowing intermediate C2/lambda when the final exponent
        // is finite. This does not change the established radiance entry point.
        (C2_NM_K / temperature) / wavelength
    }
}

fn log_rayleigh_jeans_slope(wavelength: f64) -> f64 {
    det::ln(2.0 * LIGHT_M_S * BOLTZMANN_J_K * 1.0e27) - 4.0 * det::ln(wavelength)
}

/// Inputs are ordered distinct positive temperatures. For x_c=C2/(lambda*T_c),
/// d=x_c*(T_h-T_c)/T_h, the exact identity is
/// B_h-B_c = B_h * (1-exp(-d))/(1-exp(-x_c)). No subtraction of rounded
/// exponents or radiances is involved. Computing d from actual ΔT is essential
/// when the temperatures are adjacent binary64 numbers.
fn difference_logs(hot: f64, cold: f64, wavelength: f64) -> (f64, f64) {
    let delta = hot - cold;
    let relative = delta / hot;
    let x_hot = exponent(hot, wavelength);
    let x_cold = exponent(cold, wavelength);
    let d = x_cold * relative;
    if x_cold < 1.0 {
        // (1-e^-d)/(1-e^-xc) = (ΔT/Th) phi(d)/phi(xc), with
        // phi(x)=(1-e^-x)/x. This stays valid when d or xc underflows.
        let correction = log_phi(d) - log_phi(x_cold);
        let log_b_over_t = if x_hot < 1.0e-3 {
            // B/T = RJ_slope * exp(-x)/phi(x). Keep log(T) out of the
            // secant so extreme temperatures cannot cancel useful digits.
            log_rayleigh_jeans_slope(wavelength) - x_hot - log_phi(x_hot)
        } else {
            planck_log_radiance(hot, wavelength, x_hot) - det::ln(hot)
        };
        let secant = log_b_over_t + correction;
        (secant + det::ln(delta), secant)
    } else {
        let log_fraction = log_removed(d) - log_removed(x_cold);
        let difference = planck_log_radiance(hot, wavelength, x_hot) + log_fraction;
        (difference, difference - det::ln(delta))
    }
}

fn log_removed(x: f64) -> f64 {
    det::ln(-det::expm1(-x))
}

/// log phi(x), with the analytic phi(0)=1 limit and a small-x series whose
/// omitted x^8 term is far below binary64 resolution on this branch.
fn log_phi(x: f64) -> f64 {
    if x < 1.0e-3 {
        let x2 = x * x;
        -0.5 * x + x2 * (1.0 / 24.0 + x2 * (-1.0 / 2880.0 + x2 / 181_440.0))
    } else {
        log_removed(x) - det::ln(x)
    }
}

#[cfg(test)]
mod tests;
