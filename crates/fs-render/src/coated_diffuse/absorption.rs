//! Bounded preparation of incoherent gray Beer escape and return losses.

use super::{DielectricError, GAUSS_8, det, fresnel_dielectric};

/// Eight-point panels are recursively bisected only where their two-panel
/// estimate has not converged. The original lower-index hemisphere removes
/// the Fresnel critical cusp; a separate mapped interval retains trapped rays.
pub(super) fn probabilities(
    exterior: f64,
    coating: f64,
    depth: f64,
    clear_escape: f64,
) -> Result<[f64; 2], DielectricError> {
    let normal = det::exp(-depth);
    if normal == 0.0 {
        // Every oblique ray has at least this optical depth. This is the
        // native Beer law's floating-point opaque limit, not a gain clamp.
        return Ok([0.0, 1.0]);
    }
    let high_coat = coating > exterior;
    let ratio = if high_coat { exterior / coating } else { 1.0 };
    let cone = ratio * ratio;
    let trapped = if high_coat {
        ((1.0 - ratio) * (1.0 + ratio)).sqrt()
    } else {
        0.0
    };
    let integrand = |mu: f64| -> Result<[f64; 2], DielectricError> {
        let interface = if high_coat {
            fresnel_dielectric(mu, exterior, coating)?
        } else {
            fresnel_dielectric(mu, coating, exterior)?
        };
        let internal = if high_coat {
            interface.transmitted_cosine
        } else {
            mu
        };
        let path_depth = depth / internal;
        let a = det::exp(-path_depth);
        let loss = -det::expm1(-path_depth);
        // Factoring out exp(-depth) resolves escape near the opaque limit;
        // the remaining integrand never exceeds the original clear escape.
        let relative = det::exp(-depth * ((1.0 - internal) / internal));
        let escape = 2.0 * cone * mu * relative * interface.transmittance;
        let layer = 2.0 * cone * mu * loss * (1.0 + a * interface.reflectance);
        let trapped_loss = if trapped == 0.0 {
            0.0
        } else {
            2.0 * trapped * trapped * mu * -det::expm1(-2.0 * (depth / (trapped * mu)))
        };
        Ok([escape, layer + trapped_loss])
    };
    // C <= min(1,4*tau) follows from 1-exp(-tau/mu) <= tau/mu.
    // It gives weak absorption a relative-scale tolerance without erasing
    // it in a subtraction of two nearly unit probabilities.
    let absorption_bound = if depth >= 0.25 { 1.0 } else { 4.0 * depth };
    let tolerance = [clear_escape, absorption_bound]
        .map(|scale| (2.0e-13 * scale).max(4.0 * f64::from_bits(1)));
    let initial = panel(&integrand, 0.0, 1.0)?;
    let mut remaining_panels = 2048;
    let [escape, layer] = refine(
        &integrand,
        0.0,
        1.0,
        initial,
        tolerance,
        64,
        &mut remaining_panels,
    )?;
    let result = [normal * escape, layer];
    if result.iter().any(|x| !x.is_finite() || *x < 0.0) || result[0] + result[1] > 1.0 + 8.0e-13 {
        return Err(DielectricError::InvalidAbsorption);
    }
    Ok(result)
}

fn panel(
    f: &impl Fn(f64) -> Result<[f64; 2], DielectricError>,
    lo: f64,
    hi: f64,
) -> Result<[f64; 2], DielectricError> {
    let midpoint = 0.5 * (lo + hi);
    let half_width = 0.5 * (hi - lo);
    let mut sum = [0.0; 2];
    let mut correction = [0.0; 2];
    for (node, weight) in GAUSS_8 {
        for sign in [-1.0, 1.0] {
            let value = f(midpoint + sign * half_width * node)?;
            for i in 0..2 {
                let y = weight * value[i] - correction[i];
                let next = sum[i] + y;
                correction[i] = (next - sum[i]) - y;
                sum[i] = next;
            }
        }
    }
    Ok(sum.map(|value| half_width * value))
}

fn refine(
    f: &impl Fn(f64) -> Result<[f64; 2], DielectricError>,
    lo: f64,
    hi: f64,
    coarse: [f64; 2],
    tolerance: [f64; 2],
    depth: usize,
    remaining: &mut usize,
) -> Result<[f64; 2], DielectricError> {
    if *remaining < 2 || depth == 0 {
        return Err(DielectricError::InvalidAbsorption);
    }
    *remaining -= 2;
    let mid = 0.5 * (lo + hi);
    let left = panel(f, lo, mid)?;
    let right = panel(f, mid, hi)?;
    let fine = [left[0] + right[0], left[1] + right[1]];
    if (0..2).all(|i| (fine[i] - coarse[i]).abs() <= tolerance[i]) {
        return Ok(fine);
    }
    let half = tolerance.map(|t| 0.5 * t);
    let left = refine(f, lo, mid, left, half, depth - 1, remaining)?;
    let right = refine(f, mid, hi, right, half, depth - 1, remaining)?;
    Ok([left[0] + right[0], left[1] + right[1]])
}
