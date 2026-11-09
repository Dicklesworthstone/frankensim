//! Finite-cutoff Maxwellian (LF7) and evaporation (LF9) source moments.
//!
//! The density is E^(s-1) exp(-E/theta), normalized on the declared [0,R],
//! with s=3/2 or 2. Analytic normalization defines the supplied law; no source
//! probabilities are repaired. Small R/theta uses E/R coordinates so a tiny
//! normalization is never formed. Narrow bins integrate a local power series
//! directly. Broad tail bins use upper gamma values, never two CDFs near one.
//!
//! The lower positive series and half-integer recurrence follow DLMF
//! <https://dlmf.nist.gov/8.7.E1>, <https://dlmf.nist.gov/8.4.E1> and
//! <https://dlmf.nist.gov/8.8.E2>. The strict elementary functions and scaled
//! positive products retain their original fs-math and ENDF owners.

use fs_math::{dd::Dd, decay::attenuate, det};

use crate::endf::decay::model::emissions::{DecayEmissionLimits, visit};
use crate::endf::spectrum::PositiveIntegral;
use crate::endf::{EndfError, budget};
use crate::{Accumulator, FissionError, nonnegative, positive};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(in crate::endf) enum AnalyticFissionSpectrum {
    Maxwellian,
    Evaporation,
}

/// Add unit-normalized count and first energy moments on a grid containing
/// the complete positive support. A lost positive or subnormal prepared bin
/// refuses: an arbitrarily large future stock must not amplify its rounding.
pub(in crate::endf) fn integrate(
    shape: AnalyticFissionSpectrum,
    theta_j: f64,
    cutoff_j: f64,
    edges: &[f64],
    counts: &mut [PositiveIntegral],
    energies: &mut [PositiveIntegral],
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), EndfError> {
    checkpoint()?;
    if !positive(theta_j)
        || !positive(cutoff_j)
        || counts.is_empty()
        || counts.len() != energies.len()
        || counts.len().checked_add(1) != Some(edges.len())
        || edges[0] != 0.0
        || edges[edges.len() - 1] < cutoff_j
    {
        return Err(invalid(
            "positive analytic parameters and full source support on matching groups required",
        ));
    }
    if counts.len() > limits.max_groups
        || counts
            .len()
            .checked_mul(2)
            .is_none_or(|n| n > limits.max_group_values)
    {
        return Err(budget("analytic fission spectrum groups"));
    }
    for pair in edges.windows(2) {
        visit(work, limits, checkpoint)?;
        if !nonnegative(pair[0]) || !positive(pair[1]) || pair[0] >= pair[1] {
            return Err(invalid("finite ascending analytic spectrum grid required"));
        }
    }
    let exposure = cutoff_j / theta_j;
    if !positive(exposure) {
        return Err(FissionError::Arithmetic.into());
    }
    let s = match shape {
        AnalyticFissionSpectrum::Maxwellian => 1.5,
        AnalyticFissionSpectrum::Evaporation => 2.0,
    };
    let (scale, denominator) = if exposure <= 1.0 {
        (
            cutoff_j,
            det::exp(-exposure) * lower_series(s, exposure, work, limits, checkpoint)?,
        )
    } else {
        (theta_j, lower_gamma(s, exposure, work, limits, checkpoint)?)
    };
    if !positive(denominator) {
        return Err(FissionError::Arithmetic.into());
    }
    for (group, bin) in edges.windows(2).enumerate() {
        visit(work, limits, checkpoint)?;
        let left = bin[0];
        let right = bin[1].min(cutoff_j);
        if left >= right {
            continue;
        }
        let width = right - left;
        let middle = left + 0.5 * width;
        let q = 0.5 * (width / middle);
        let z = 0.5 * (width / theta_j);
        let mut count = PositiveIntegral::new();
        let mut energy = PositiveIntegral::new();
        if left > 0.0 && q <= 0.5 && z <= 0.5 {
            // (1+q*t)^p exp(-z*t), -1<=t<=1. Keeping width separate
            // preserves adjacent energy knots even if dividing both by theta
            // rounds them to the same floating-point coordinate.
            let scaled_width = width / scale;
            let scaled_middle = middle / scale;
            // A midpoint between adjacent floats can round to an endpoint.
            // Its tiny residual becomes a substantial relative tail error
            // inside exp(-E/theta). Preserve the exposure, including division
            // roundoff, with the original double-double arithmetic owner.
            let theta = Dd::from_f64(theta_j);
            let midpoint_exposure =
                Dd::from_f64(left) / theta + (Dd::from_f64(width) / theta) * Dd::from_f64(0.5);
            if !nonnegative(midpoint_exposure.hi) || !midpoint_exposure.lo.is_finite() {
                return Err(FissionError::Arithmetic.into());
            }
            let half_decay =
                det::exp(-0.5 * midpoint_exposure.hi) * det::exp(-0.5 * midpoint_exposure.lo);
            let power = if s == 1.5 {
                scaled_middle.sqrt()
            } else {
                scaled_middle
            };
            let area = midpoint_series(s - 1.0, q, z, work, limits, checkpoint)?;
            let first = midpoint_series(s, q, z, work, limits, checkpoint)?;
            if !positive(scaled_width) || !positive(power) || !positive(half_decay) {
                return Err(FissionError::Arithmetic.into());
            }
            count.add_triple_product(
                scaled_width,
                power * area,
                half_decay,
                half_decay,
                denominator,
            )?;
            energy.add_four_product(
                scaled_width,
                power * first,
                half_decay,
                half_decay,
                middle,
                denominator,
            )?;
        } else if right / theta_j <= 1.0 {
            // Broad low-energy interval. Relative primitives are bounded away
            // from cancellation; powers of the dimensional theta never occur.
            let xr = right / theta_j;
            let xl = left / theta_j;
            let ratio = left / right;
            let relative_power = if s == 1.5 {
                ratio * ratio.sqrt()
            } else {
                ratio * ratio
            };
            let correction = det::exp(xr - xl);
            let mut area = lower_series(s, xr, work, limits, checkpoint)?;
            let mut first = lower_series(s + 1.0, xr, work, limits, checkpoint)?;
            if left > 0.0 {
                area -=
                    relative_power * correction * lower_series(s, xl, work, limits, checkpoint)?;
                first -= relative_power
                    * ratio
                    * correction
                    * lower_series(s + 1.0, xl, work, limits, checkpoint)?;
            }
            let coordinate = right / scale;
            let power = if s == 1.5 {
                coordinate.sqrt()
            } else {
                coordinate
            };
            let half_decay = det::exp(-0.5 * xr);
            if !positive(area) || !positive(first) || !positive(coordinate) || !positive(power) {
                return Err(FissionError::Arithmetic.into());
            }
            count.add_triple_product(
                coordinate,
                power * area,
                half_decay,
                half_decay,
                denominator,
            )?;
            energy.add_four_product(
                coordinate,
                power * first,
                half_decay,
                half_decay,
                right,
                denominator,
            )?;
        } else {
            let xl = left / theta_j;
            let xr = right / theta_j;
            // If both endpoints exceed theta, every remaining bin is wider
            // than theta. Upper-tail values are then well separated; no narrow
            // interval or probabilities rounded to one are subtracted.
            let (area, first) = if xl >= 1.0 {
                (
                    upper_gamma(s, xl) - upper_gamma(s, xr),
                    upper_gamma(s + 1.0, xl) - upper_gamma(s + 1.0, xr),
                )
            } else {
                (
                    lower_gamma(s, xr, work, limits, checkpoint)?
                        - lower_gamma(s, xl, work, limits, checkpoint)?,
                    lower_gamma(s + 1.0, xr, work, limits, checkpoint)?
                        - lower_gamma(s + 1.0, xl, work, limits, checkpoint)?,
                )
            };
            if !positive(area) || !positive(first) {
                return Err(FissionError::Arithmetic.into());
            }
            count.add_product(1.0, area, 1.0, denominator)?;
            energy.add_product(theta_j, first, 1.0, denominator)?;
        }
        let count = count.value()?;
        let energy = energy.value()?;
        if !count.is_normal() || !energy.is_normal() {
            return Err(FissionError::Arithmetic.into());
        }
        counts[group].add_product(1.0, count, 1.0, 1.0)?;
        energies[group].add_product(1.0, energy, 1.0, 1.0)?;
    }
    checkpoint()?;
    Ok(())
}

// gamma(s,x)=x^s exp(-x) S(s,x), where S=sum x^k/(s)_(k+1).
// Only 0<=x<=1 reaches this positive, geometrically bounded series.
fn lower_series(
    s: f64,
    x: f64,
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<f64, EndfError> {
    let mut term = 1.0 / s;
    let mut total = Accumulator::new(term);
    for k in 1..=64 {
        visit(work, limits, checkpoint)?;
        term *= x / (s + f64::from(k));
        total.add(term)?;
        let tail = term * (x / (s + f64::from(k + 1))) / (1.0 - x / (s + f64::from(k + 2)));
        if tail <= 0.125 * f64::EPSILON * total.value() {
            return Ok(total.value());
        }
    }
    Err(budget("analytic lower-gamma series"))
}

fn lower_gamma(
    s: f64,
    x: f64,
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<f64, EndfError> {
    visit(work, limits, checkpoint)?;
    if x == 0.0 {
        return Ok(0.0);
    }
    if x <= 1.0 {
        // This raw primitive is used only beside a well-separated O(1)
        // endpoint. Small-only bins use the separately scaled path above.
        let power = match s {
            1.5 => x * x.sqrt(),
            2.0 => x * x,
            2.5 => x * x * x.sqrt(),
            3.0 => x * x * x,
            _ => unreachable!("only LF7/LF9 count and energy shapes"),
        };
        return Ok(power * det::exp(-x) * lower_series(s, x, work, limits, checkpoint)?);
    }
    let complete = match s {
        1.5 => 0.5 * core::f64::consts::PI.sqrt(),
        2.0 => 1.0,
        2.5 => 0.75 * core::f64::consts::PI.sqrt(),
        3.0 => 2.0,
        _ => unreachable!("only LF7/LF9 count and energy shapes"),
    };
    Ok(complete - upper_gamma(s, x))
}

fn upper_gamma(s: f64, x: f64) -> f64 {
    // Even x^2 exp(-x) is below binary64 at this cutoff. Returning the
    // rounded primitive is safe; an isolated lost positive bin still refuses.
    if x >= 1024.0 {
        return 0.0;
    }
    match s {
        1.5 | 2.5 => {
            let root = x.sqrt();
            let lower = attenuate(root, x) + 0.5 * core::f64::consts::PI.sqrt() * det::erfc(root);
            if s == 1.5 {
                lower
            } else {
                1.5 * lower + attenuate(x * root, x)
            }
        }
        2.0 => attenuate(x + 1.0, x),
        3.0 => attenuate(x * x + 2.0 * x + 2.0, x),
        _ => unreachable!("only LF7/LF9 count and energy shapes"),
    }
}

// The differential identity (1+q*t)F'=(p*q-z-z*q*t)F generates
// F=(1+q*t)^p exp(-z*t). Integrating even terms over [-1,1] avoids
// neighboring incomplete-gamma subtraction. On the complex circle of radius
// rho=(7/8)/max(q,z), |F| <= (15/8)^2 exp(7/8) < 9. Cauchy's bound therefore
// bounds the entire omitted even tail, independently of coefficient zeros.
fn midpoint_series(
    p: f64,
    q: f64,
    z: f64,
    work: &mut usize,
    limits: DecayEmissionLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<f64, EndfError> {
    let mut previous = 1.0;
    let mut current = p * q - z;
    let mut total = Accumulator::new(1.0);
    let inverse_radius = q.max(z) / 0.875;
    let inverse_square = inverse_radius * inverse_radius;
    let mut radius_power = inverse_radius;
    for n in 1..128 {
        visit(work, limits, checkpoint)?;
        radius_power *= inverse_radius;
        let next = ((p * q - z - f64::from(n) * q) * current - z * q * previous) / f64::from(n + 1);
        if (n + 1) % 2 == 0 {
            let term = next / f64::from(n + 2);
            total.add(term)?;
            let tail =
                9.0 * radius_power * inverse_square / (f64::from(n + 4) * (1.0 - inverse_square));
            if tail <= 0.125 * f64::EPSILON * total.value().abs() {
                let value = total.value();
                return if positive(value) {
                    Ok(value)
                } else {
                    Err(FissionError::Arithmetic.into())
                };
            }
        }
        previous = current;
        current = next;
    }
    Err(budget("analytic fission local-moment series"))
}

fn invalid(message: &'static str) -> EndfError {
    FissionError::InvalidInput(message).into()
}

#[cfg(test)]
mod tests;
