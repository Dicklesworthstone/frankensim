//! Analytic positive moments, including integrable singular yield surfaces.

use super::{SlurryError, require};
use fs_math::det;

/// Product/quotient with a logarithmic fallback for intermediate range loss.
/// A positive mathematical response that is not representable is refused.
pub(super) fn product(numerator: &[f64], denominator: &[f64]) -> Result<f64, SlurryError> {
    require(
        numerator.iter().all(|x| x.is_finite() && *x >= 0.0)
            && denominator.iter().all(|x| x.is_finite() && *x > 0.0),
        "invalid slit product factors",
    )?;
    if numerator.contains(&0.0) {
        return Ok(0.0);
    }
    let direct = numerator.iter().product::<f64>() / denominator.iter().product::<f64>();
    let value = if direct.is_normal() {
        direct
    } else {
        det::exp(
            numerator.iter().map(|x| det::ln(*x)).sum::<f64>()
                - denominator.iter().map(|x| det::ln(*x)).sum::<f64>(),
        )
    };
    require(
        value.is_finite() && value > 0.0,
        "positive slit response is unrepresentable",
    )?;
    Ok(value)
}

/// Normalized moments of `(1-d*t)^(a-1)` on `0 <= t <= 1`.
/// `span` preserves the known stress difference when endpoint subtraction
/// loses a small pressure increment. `low=0` is an actual yield surface.
#[derive(Clone, Copy, Debug)]
pub(super) struct Distribution {
    pub average: f64,
    pub mean: f64,
    pub variance: f64,
    pub wall_moment: f64,
}

pub(super) fn distribution(
    a: f64,
    low: f64,
    high: f64,
    span: f64,
) -> Result<Distribution, SlurryError> {
    require(
        a.is_finite()
            && a > 0.0
            && low.is_finite()
            && low >= 0.0
            && high.is_finite()
            && high > 0.0
            && low <= high
            && span.is_finite()
            && span > 0.0
            && span <= high,
        "invalid slit power-moment interval",
    )?;
    let d = span / high;
    require(
        d > 0.0,
        "slit stress interval is below numerical resolution",
    )?;
    let result = if low == 0.0 {
        // Beta(1,a). Keeping the positive variance and E[t(1-t)] separately
        // avoids subtracting nearly equal moments when a is very small.
        require(span == high, "yield-surface interval must reach zero")?;
        let inverse1 = 1.0 / (a + 1.0);
        let inverse2 = 1.0 / (a + 2.0);
        let wall_moment = product(&[a, inverse1, inverse2], &[])?;
        Distribution {
            average: 1.0 / a,
            mean: inverse1,
            variance: product(&[wall_moment, inverse1], &[])?,
            wall_moment,
        }
    } else if d < 0.25 && (a - 1.0).abs() * d <= 0.5 {
        // Convergent divided-difference series, not spatial quadrature.
        let (mut f0, mut f1, mut f2, mut fw) = (1.0, 0.5, 1.0 / 3.0, 1.0 / 6.0);
        let mut coefficient = 1.0;
        for k in 1..=48 {
            let kf = f64::from(k);
            coefficient *= (kf - a) * d / kf;
            f0 += coefficient / (kf + 1.0);
            f1 += coefficient / (kf + 2.0);
            f2 += coefficient / (kf + 3.0);
            fw += coefficient / ((kf + 2.0) * (kf + 3.0));
            if coefficient == 0.0 || (k >= 4 && coefficient.abs() < f64::EPSILON * f0.abs()) {
                break;
            }
        }
        let mean = f1 / f0;
        Distribution {
            average: f0,
            mean,
            variance: f2 / f0 - mean * mean,
            wall_moment: fw / f0,
        }
    } else {
        let r = 1.0 - d;
        let logarithm = if r > 0.0 && r < 1.0 {
            det::ln(r) * (d / (1.0 - r))
        } else if r == 1.0 {
            -d
        } else {
            det::ln(low) - det::ln(high)
        };
        let exponent = a * logarithm;
        let denominator = -det::expm1(exponent);
        let tail = det::exp(exponent) / denominator;
        let inverse1 = (1.0 / (a + 1.0)) / d;
        let inverse2 = (1.0 / (a + 2.0)) / d;
        let correction = if tail == 0.0 { 0.0 } else { (a * d) * tail };
        let mean = (1.0 - correction) * inverse1;
        let second = (2.0 - correction * (2.0 + (a + 1.0) * d)) * inverse1 * inverse2;
        Distribution {
            average: (denominator / a) / d,
            mean,
            variance: second - mean * mean,
            wall_moment: mean - second,
        }
    };
    require(
        result.average.is_finite()
            && result.average > 0.0
            && result.mean.is_finite()
            && (0.0..=1.0).contains(&result.mean)
            && result.variance.is_finite()
            && result.variance > 0.0
            && result.wall_moment.is_finite()
            && result.wall_moment > 0.0,
        "slit power moments are unrepresentable",
    )?;
    Ok(result)
}

/// Moments of inverse differential viscosity about z=1/2-y/h.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Moments {
    pub weight: f64,
    pub mean: f64,
    pub variance: f64,
    pub wall_weight: f64,
}

impl Moments {
    pub fn add(
        &mut self,
        weight: f64,
        mean: f64,
        variance: f64,
        wall_moment: f64,
    ) -> Result<(), SlurryError> {
        let total = self.weight + weight;
        require(
            total.is_finite() && total > 0.0,
            "slit moment weight overflow",
        )?;
        let delta = mean - self.mean;
        self.variance += product(&[weight, variance], &[])?
            + product(&[self.weight, weight, delta.abs(), delta.abs()], &[total])?;
        self.mean += (weight / total) * delta;
        self.weight = total;
        self.wall_weight += product(&[weight, wall_moment], &[])?;
        require(
            self.variance.is_finite()
                && self.variance > 0.0
                && self.wall_weight.is_finite()
                && self.wall_weight > 0.0,
            "slit positive moment overflow",
        )
    }
}
