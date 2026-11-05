//! Hat-function integrals of the normalized causal semicircle kernel.

use crate::det;

/// Integrate `4/pi * sqrt(1-s*s)`, `0 <= s <= 1`, against the linear
/// interpolation hats centered on `j * step_over_support`.
///
/// This is the regular part of the uniform baffled piston's radiation impulse
/// response, with time normalized by the aperture transit time `2a/c`.
/// Physical impedance scales and stepping belong to the consuming solvers.
/// The nonnegative weights sum to one up to roundoff. Their first moment is
/// `4/(3*pi)`; linear interpolation reproduces that moment exactly.
///
/// Refuses invalid spacing, an exceeded caller tap budget, failed allocation,
/// or unrepresentable/nonpositive quadrature. No allocation occurs after setup.
pub fn semicircle_hat_weights(
    step_over_support: f64,
    maximum_taps: usize,
) -> Result<Vec<f64>, &'static str> {
    let h = step_over_support;
    if !h.is_finite() || h <= 0.0 {
        return Err("semicircle kernel needs finite positive normalized spacing");
    }
    let taps = (1.0 / h).ceil() + 1.0;
    if !taps.is_finite() || taps < 2.0 || taps > maximum_taps as f64 {
        return Err("semicircle kernel exceeds its tap budget");
    }
    let mut weights = Vec::new();
    weights
        .try_reserve_exact(taps as usize)
        .map_err(|_| "semicircle kernel allocation failed")?;
    let f0 = |s: f64| 0.5 * (s * det::sqrt(1.0 - s * s) + det::asin(s));
    let f1 = |s: f64| -(1.0 - s * s) * det::sqrt(1.0 - s * s) / 3.0;
    for j in 0..taps as usize {
        let centre = j as f64 * h;
        let mut weight = 0.0;
        if j > 0 {
            let (a, b) = ((centre - h).max(0.0), centre.min(1.0));
            if a < b {
                weight += (f1(b) - f1(a) - a * (f0(b) - f0(a))) / h;
            }
        }
        let (a, b) = (centre, (centre + h).min(1.0));
        if a < b {
            weight += ((centre + h) * (f0(b) - f0(a)) - (f1(b) - f1(a))) / h;
        }
        weight *= 4.0 / core::f64::consts::PI;
        if !weight.is_finite() || weight < 0.0 {
            return Err("semicircle kernel weight is not finite and nonnegative");
        }
        weights.push(weight);
    }
    let total: f64 = weights.iter().sum();
    if (total - 1.0).abs() > 1e-9 {
        return Err("semicircle kernel lost its unit integral");
    }
    Ok(weights)
}
