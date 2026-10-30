//! Normalized piecewise-linear angular laws and matched inverse-CDF sampling.

use std::sync::Arc;

use super::{DomainHasher, VolumeTransportError, det};

/// Maximum knot count in one admitted angular phase table.
pub const MAX_PHASE_KNOTS: usize = 4096;

#[derive(Debug, PartialEq)]
struct AngularTable {
    // Cosine of propagation deflection and density per steradian.
    knots: Box<[(f64, f64)]>,
    cumulative: Box<[f64]>,
    mean_cosine: f64,
}

/// Reciprocal, azimuthally symmetric scattering supplied as angular data.
///
/// Between knots the density is linear in `mu = cos(theta)`, not in degrees.
/// Normalization uses `2 pi integral p(mu) dmu = 1`, including both forward
/// and backward hemispheres. The same interpolant supplies evaluation and
/// inverse-CDF sampling. The immutable storage is shared by cheap clones.
///
/// This is a declared wavelength-independent angular model. Normalizing a
/// table does not establish measurement provenance, angular resolution, an
/// interpolation-error bound, polarization, or a particle-size model.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedPhase {
    data: Arc<AngularTable>,
}

impl TabulatedPhase {
    /// Admit 2–[`MAX_PHASE_KNOTS`] strictly increasing cosine knots with
    /// endpoints exactly -1 and +1, and finite nonnegative relative values.
    /// At least one value must be positive. Common positive scaling is removed.
    ///
    /// Zero-density intervals retain zero density and receive no samples.
    /// Positive mass lost to underflow or to cumulative-probability rounding,
    /// unrepresentable PDFs, and active intervals without a floating-point
    /// interior refuse with [`VolumeTransportError::Unrepresentable`].
    pub fn try_new(knots: &[(f64, f64)]) -> Result<Self, VolumeTransportError> {
        if !(2..=MAX_PHASE_KNOTS).contains(&knots.len())
            || knots.first().map(|knot| knot.0) != Some(-1.0)
            || knots.last().map(|knot| knot.0) != Some(1.0)
        {
            return Err(VolumeTransportError::Phase);
        }
        let mut maximum = 0.0_f64;
        for (index, &(cosine, value)) in knots.iter().enumerate() {
            if !cosine.is_finite()
                || !value.is_finite()
                || value < 0.0
                || (index > 0 && cosine <= knots[index - 1].0)
            {
                return Err(VolumeTransportError::Phase);
            }
            maximum = maximum.max(value);
        }
        if maximum == 0.0 {
            return Err(VolumeTransportError::Phase);
        }
        let mut normalized = Vec::new();
        let mut cumulative = Vec::new();
        normalized
            .try_reserve_exact(knots.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        cumulative
            .try_reserve_exact(knots.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        for &(cosine, value) in knots {
            let scaled = value / maximum;
            if value > 0.0 && scaled == 0.0 {
                return Err(VolumeTransportError::Unrepresentable);
            }
            normalized.push((if cosine == 0.0 { 0.0 } else { cosine }, scaled));
        }
        let (mut total, mut correction) = (0.0, 0.0);
        for pair in normalized.windows(2) {
            let width = pair[1].0 - pair[0].0;
            let active = pair[0].1 > 0.0 || pair[1].1 > 0.0;
            let mass = width * (0.5 * pair[0].1 + 0.5 * pair[1].1);
            if active && (mass == 0.0 || pair[0].0.next_up() >= pair[1].0) {
                return Err(VolumeTransportError::Unrepresentable);
            }
            add(&mut total, &mut correction, mass);
        }
        let total = total + correction;
        let solid_angle_integral = core::f64::consts::TAU * total;
        if !(solid_angle_integral > 0.0 && solid_angle_integral.is_finite()) {
            return Err(VolumeTransportError::Unrepresentable);
        }
        let (mut sum, mut sum_correction) = (0.0, 0.0);
        let (mut mean, mut mean_correction) = (0.0, 0.0);
        cumulative.push(0.0);
        for (index, pair) in normalized.windows(2).enumerate() {
            let width = pair[1].0 - pair[0].0;
            let left = pair[0].1;
            let right = pair[1].1;
            let mass = width * (0.5 * left + 0.5 * right);
            add(&mut sum, &mut sum_correction, mass);
            let probability = if index + 2 == normalized.len() {
                1.0
            } else {
                (sum + sum_correction) / total
            };
            let previous = cumulative[index];
            if !(previous..=1.0).contains(&probability) || (mass > 0.0 && probability <= previous) {
                return Err(VolumeTransportError::Unrepresentable);
            }
            cumulative.push(probability);
            if mass > 0.0 {
                // The centroid of a linear density is analytic, including a
                // triangular interval with one zero endpoint.
                let center = pair[0].0 + width * (left + 2.0 * right) / (3.0 * (left + right));
                add(&mut mean, &mut mean_correction, (mass / total) * center);
            }
        }
        for (_, value) in &mut normalized {
            let density = *value / solid_angle_integral;
            if !density.is_finite() || (*value > 0.0 && density == 0.0) {
                return Err(VolumeTransportError::Unrepresentable);
            }
            *value = if density == 0.0 { 0.0 } else { density };
        }
        Ok(Self {
            data: Arc::new(AngularTable {
                knots: normalized.into_boxed_slice(),
                cumulative: cumulative.into_boxed_slice(),
                mean_cosine: mean + mean_correction,
            }),
        })
    }

    /// The admitted cosine knots and their normalized density in sr⁻¹.
    #[must_use]
    pub fn knots_per_sr(&self) -> &[(f64, f64)] {
        &self.data.knots
    }

    /// Mean propagation deflection cosine of the piecewise-linear model.
    #[must_use]
    pub fn mean_cosine(&self) -> f64 {
        self.data.mean_cosine
    }

    /// Evaluate the normalized table in sr⁻¹ without angular extrapolation.
    /// Tiny unit-vector dot-product roundoff clamps to the endpoint.
    pub fn eval_cos(&self, cosine: f64) -> Result<f64, VolumeTransportError> {
        if !cosine.is_finite() || cosine.abs() > 1.0 + 2.0e-10 {
            return Err(VolumeTransportError::Phase);
        }
        Ok(self.eval_clamped(cosine.clamp(-1.0, 1.0)))
    }

    pub(super) fn eval_clamped(&self, cosine: f64) -> f64 {
        let knots = &self.data.knots;
        let upper = knots.partition_point(|&(mu, _)| mu <= cosine);
        let index = upper.saturating_sub(1).min(knots.len() - 2);
        let (lo, left) = knots[index];
        let (hi, right) = knots[index + 1];
        let fraction = ((cosine - lo) / (hi - lo)).clamp(0.0, 1.0);
        left * (1.0 - fraction) + right * fraction
    }

    fn interval(&self, uniform: f64) -> usize {
        // Right-biased lookup skips all zero-mass intervals, including a
        // leading dark cap and a dark gap at an exact CDF knot.
        self.data
            .cumulative
            .partition_point(|&value| value <= uniform)
            - 1
    }

    pub(super) fn sample_cosine(&self, uniform: f64) -> Result<f64, VolumeTransportError> {
        let index = self.interval(uniform);
        let (lo, left) = self.data.knots[index];
        let (hi, right) = self.data.knots[index + 1];
        let start = self.data.cumulative[index];
        let end = self.data.cumulative[index + 1];
        let local = ((uniform - start) / (end - start)).clamp(0.0, 1.0);
        // Scale each interval before solving its quadratic CDF: squaring a
        // physical PDF can overflow, while these endpoints lie in [0,1].
        let scale = left.max(right);
        let a = left / scale;
        let b = right / scale;
        let mass = 0.5 * local * (a + b);
        let root = det::sqrt((2.0 * (b - a)).mul_add(mass, a * a).max(0.0));
        let fraction = if mass == 0.0 {
            0.0
        } else {
            (2.0 * mass / (a + root)).clamp(0.0, 1.0)
        };
        let mut cosine = (hi - lo).mul_add(fraction, lo).clamp(lo, hi);
        if self.eval_clamped(cosine) == 0.0 {
            // CDF endpoints with zero density have measure zero. Resolve them
            // inside the same active interval instead of creating a PDF floor.
            cosine = if cosine == lo {
                lo.next_up()
            } else {
                hi.next_down()
            };
            if self.eval_clamped(cosine) == 0.0 {
                cosine = self.interior_cosine(uniform);
            }
        }
        if !cosine.is_finite() || self.eval_clamped(cosine) == 0.0 {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(cosine)
    }

    pub(super) fn interior_cosine(&self, uniform: f64) -> f64 {
        let index = self.interval(uniform);
        let lo = self.data.knots[index].0;
        let hi = self.data.knots[index + 1].0;
        (hi - lo).mul_add(0.5, lo)
    }

    pub(super) fn hash_into(&self, hash: &mut DomainHasher) {
        hash.update(&(self.data.knots.len() as u64).to_le_bytes());
        for (&(cosine, density), &cumulative) in
            self.data.knots.iter().zip(self.data.cumulative.iter())
        {
            hash.update(&cosine.to_bits().to_le_bytes());
            hash.update(&density.to_bits().to_le_bytes());
            // Retain the represented proposal too: its cumulative rounding
            // affects exact replay even when density knots round identically.
            hash.update(&cumulative.to_bits().to_le_bytes());
        }
        hash.update(&self.data.mean_cosine.to_bits().to_le_bytes());
    }
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let next = *sum + value;
    *correction += if sum.abs() >= value.abs() {
        (*sum - next) + value
    } else {
        (value - next) + *sum
    };
    *sum = next;
}

#[cfg(test)]
mod tests;
