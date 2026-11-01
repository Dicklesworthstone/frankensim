//! Importance sampling of the source's represented angular intensity.
//!
//! Tables are piecewise bilinear in cosine and azimuth. Their polar marginal
//! and conditional azimuth densities are therefore piecewise linear, with
//! analytic inverse CDFs. Preparation takes linear work and bounded storage;
//! each sample uses two binary searches and no allocation. Analytic profiles
//! need neither a table nor a search. This samples the declared profile, not
//! an inferred manufacturer's beam or a spatially extended emitter.

use super::{AngularDistribution, AngularLaw, AngularTable, PhotometricSource, TAU, unit};
use crate::instances::RigidTransform;
use fs_exec::{Cancelled, Cx};
use fs_geom::Vec3;
use fs_math::det;

/// Admission, cancellation, or floating-point representability refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SourceSamplingError {
    /// Samples must be finite and strictly inside the unit square.
    InvalidSample,
    /// A positive cell mass, direction, or density is not representable.
    NumericalRange,
    /// Preparation was cancelled before publishing a sampler.
    Cancelled(Cancelled),
}

impl From<Cancelled> for SourceSamplingError {
    fn from(error: Cancelled) -> Self {
        Self::Cancelled(error)
    }
}
impl core::fmt::Display for SourceSamplingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "angular source sampling: {self:?}")
    }
}
impl core::error::Error for SourceSamplingError {}

/// A sampled outgoing direction and its probability per solid angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DirectionSample {
    /// Unit local direction for an angular sampler, world direction for a source.
    pub direction: Vec3,
    /// Density per steradian, including the prepared CDF's normalization.
    pub pdf_per_sr: f64,
    angular_weight: f64,
}

impl DirectionSample {
    /// Original represented angular intensity divided by the sampling PDF.
    ///
    /// This restores the declared source distribution, including its floating
    /// normalization residue. Multiply by source lumens and W/(lm nm) to form
    /// a physical W/nm packet. The density and weight are evaluated in local
    /// sampling coordinates; no inverse-rotation or narrow-cone boundary
    /// reevaluation is needed. It is not a measured-source accuracy bound.
    #[must_use]
    pub const fn angular_weight(self) -> f64 {
        self.angular_weight
    }
}

#[derive(Debug)]
struct TableCdf {
    azimuth: Vec<f64>,
    polar: Vec<f64>,
}

/// A reusable immutable angular proposal borrowing its original profile.
#[derive(Debug)]
pub struct PreparedAngularSampler<'a> {
    distribution: &'a AngularDistribution,
    table: Option<TableCdf>,
    cosine_power_deficit: f64,
    integral: f64,
}

/// The same proposal placed by the source's proper-rigid world transform.
#[derive(Debug)]
pub struct PreparedSourceSampler<'a> {
    transform: RigidTransform,
    angular: PreparedAngularSampler<'a>,
}

impl AngularDistribution {
    /// Prepare the exact represented profile, up to floating-point rounding.
    ///
    /// Table storage is one CDF value per admitted nodal entry plus one per
    /// polar row. Cancellation is polled every 1024 table entries. Positive
    /// masses lost to CDF representability are refused instead of discarded.
    pub fn prepare_direction_sampler(
        &self,
        cx: &Cx<'_>,
    ) -> Result<PreparedAngularSampler<'_>, SourceSamplingError> {
        cx.checkpoint()?;
        let mut table = None;
        let mut cosine_power_deficit = 0.0;
        let integral = match &self.law {
            AngularLaw::Isotropic => self.peak_per_sr * (2.0 * TAU),
            AngularLaw::Cone { cosine_cutoff } => self.peak_per_sr * (TAU * (1.0 - cosine_cutoff)),
            AngularLaw::CosinePower {
                exponent,
                cosine_cutoff,
                ..
            } => {
                let power = exponent + 1.0;
                cosine_power_deficit = if *cosine_cutoff == 0.0 {
                    1.0
                } else {
                    -det::expm1(power * det::ln(*cosine_cutoff))
                };
                self.peak_per_sr * ((TAU * cosine_power_deficit) / power)
            }
            AngularLaw::Table(profile) => {
                let prepared = TableCdf::new(cx, profile)?;
                let integral = *prepared.polar.last().expect("admitted nonempty table");
                table = Some(prepared);
                integral
            }
        };
        if !integral.is_finite() || integral <= 0.0 {
            return Err(SourceSamplingError::NumericalRange);
        }
        cx.checkpoint()?;
        Ok(PreparedAngularSampler {
            distribution: self,
            table,
            cosine_power_deficit,
            integral,
        })
    }
}

impl PhotometricSource {
    /// Prepare once for repeated world-space emission sampling.
    pub fn prepare_direction_sampler(
        &self,
        cx: &Cx<'_>,
    ) -> Result<PreparedSourceSampler<'_>, SourceSamplingError> {
        Ok(PreparedSourceSampler {
            transform: self.transform,
            angular: self.distribution.prepare_direction_sampler(cx)?,
        })
    }
}

impl PreparedSourceSampler<'_> {
    /// Sample the full authored support in world space without allocation.
    pub fn sample(&self, sample: [f64; 2]) -> Result<DirectionSample, SourceSamplingError> {
        let mut result = self.angular.sample(sample)?;
        result.direction = unit(self.transform.transform_vector(result.direction))
            .map_err(|_| SourceSamplingError::NumericalRange)?;
        Ok(result)
    }
}

impl PreparedAngularSampler<'_> {
    /// Sample cosine and azimuth from two independent open-unit inputs.
    ///
    /// A sample that rounds onto a zero-density boundary returns a numerical
    /// refusal. No rejection loop, density floor, or support clipping is used.
    pub fn sample(&self, [u, v]: [f64; 2]) -> Result<DirectionSample, SourceSamplingError> {
        if [u, v]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0 || *x >= 1.0)
        {
            return Err(SourceSamplingError::InvalidSample);
        }
        let (mu, phi, density) = match &self.distribution.law {
            AngularLaw::Isotropic => (1.0 - 2.0 * u, TAU * v, self.distribution.peak_per_sr),
            AngularLaw::Cone { cosine_cutoff } => (
                1.0 - u * (1.0 - cosine_cutoff),
                TAU * v,
                self.distribution.peak_per_sr,
            ),
            AngularLaw::CosinePower {
                exponent,
                cosine_cutoff,
                ..
            } => {
                let mu = if *exponent == 0.0 {
                    1.0 - u * (1.0 - cosine_cutoff)
                } else {
                    let t = u * self.cosine_power_deficit;
                    let complement = 1.0 - t;
                    // Correct the rounded subtraction before taking a very
                    // small logarithm (narrow profiles otherwise collapse).
                    let log = det::ln(complement) + ((1.0 - complement) - t) / complement;
                    det::exp(log / (exponent + 1.0))
                };
                (
                    mu,
                    TAU * v,
                    self.distribution.peak_per_sr * det::pow(mu, *exponent),
                )
            }
            AngularLaw::Table(profile) => self
                .table
                .as_ref()
                .expect("prepared table matches immutable profile")
                .sample(profile, self.integral, u, v)?,
        };
        let radial = det::sqrt((1.0 - mu) * (1.0 + mu));
        let direction = Vec3::new(radial * det::cos(phi), radial * det::sin(phi), mu);
        let pdf_per_sr = density / self.integral;
        if !pdf_per_sr.is_finite()
            || pdf_per_sr <= 0.0
            || ![direction.x, direction.y, direction.z]
                .iter()
                .all(|x| x.is_finite())
        {
            return Err(SourceSamplingError::NumericalRange);
        }
        Ok(DirectionSample {
            direction,
            pdf_per_sr,
            angular_weight: self.integral,
        })
    }
}

impl TableCdf {
    fn new(cx: &Cx<'_>, profile: &AngularTable) -> Result<Self, SourceSamplingError> {
        let columns = profile.phi.len();
        let rows = profile.cosine.len();
        let mut azimuth = vec![0.0; profile.values.len()];
        let mut work = 0_usize;
        for row in 0..rows {
            let mut sum = Cumulative::default();
            for column in 1..columns {
                if work % 1024 == 0 {
                    cx.checkpoint()?;
                }
                work += 1;
                let i = row * columns + column;
                let mass = (0.5 * profile.values[i - 1] + 0.5 * profile.values[i])
                    * (profile.phi[column] - profile.phi[column - 1]);
                if mass == 0.0 && (profile.values[i - 1] > 0.0 || profile.values[i] > 0.0) {
                    return Err(SourceSamplingError::NumericalRange);
                }
                azimuth[i] = sum.add(mass)?;
            }
        }
        let mut polar = vec![0.0; rows];
        let mut sum = Cumulative::default();
        for row in 1..rows {
            if row % 1024 == 0 {
                cx.checkpoint()?;
            }
            let upper = azimuth[row * columns - 1];
            let lower = azimuth[(row + 1) * columns - 1];
            let mass =
                (0.5 * upper + 0.5 * lower) * (profile.cosine[row - 1] - profile.cosine[row]);
            if mass == 0.0 && (upper > 0.0 || lower > 0.0) {
                return Err(SourceSamplingError::NumericalRange);
            }
            polar[row] = sum.add(mass)?;
        }
        Ok(Self { azimuth, polar })
    }

    fn sample(
        &self,
        profile: &AngularTable,
        integral: f64,
        u: f64,
        v: f64,
    ) -> Result<(f64, f64, f64), SourceSamplingError> {
        let columns = profile.phi.len();
        let target = u * integral;
        let row = self
            .polar
            .partition_point(|x| *x <= target)
            .saturating_sub(1);
        if row + 1 >= profile.cosine.len() {
            return Err(SourceSamplingError::NumericalRange);
        }
        let upper = &self.azimuth[row * columns..(row + 1) * columns];
        let lower = &self.azimuth[(row + 1) * columns..(row + 2) * columns];
        let t = inverse_linear(
            upper[columns - 1],
            lower[columns - 1],
            (target - self.polar[row]) / (self.polar[row + 1] - self.polar[row]),
        );
        if !t.is_finite() || !(0.0..=1.0).contains(&t) {
            return Err(SourceSamplingError::NumericalRange);
        }
        let conditional = |column: usize| (1.0 - t) * upper[column] + t * lower[column];
        let target = v * conditional(columns - 1);
        // Binary search the interpolated row CDF without allocating a row.
        let mut low = 0;
        let mut high = columns;
        while low < high {
            let mid = low + (high - low) / 2;
            if conditional(mid) <= target {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        let column = low.saturating_sub(1);
        if column + 1 >= columns {
            return Err(SourceSamplingError::NumericalRange);
        }
        let i = row * columns + column;
        let a = (1.0 - t) * profile.values[i] + t * profile.values[i + columns];
        let b = (1.0 - t) * profile.values[i + 1] + t * profile.values[i + columns + 1];
        let s = inverse_linear(
            a,
            b,
            (target - conditional(column)) / (conditional(column + 1) - conditional(column)),
        );
        if !s.is_finite() || !(0.0..=1.0).contains(&s) {
            return Err(SourceSamplingError::NumericalRange);
        }
        let mu = (1.0 - t) * profile.cosine[row] + t * profile.cosine[row + 1];
        let phi = (1.0 - s) * profile.phi[column] + s * profile.phi[column + 1];
        Ok((mu, phi, (1.0 - s) * a + s * b))
    }
}

/// Invert q = (a*t + (b-a)*t²/2) / ((a+b)/2) without cancellation.
fn inverse_linear(a: f64, b: f64, q: f64) -> f64 {
    if q == 0.0 {
        return 0.0;
    }
    if q == 1.0 {
        return 1.0;
    }
    let scale = a.max(b);
    let a = a / scale;
    let b = b / scale;
    q * (a + b) / (a + det::sqrt((1.0 - q) * a * a + q * b * b))
}

#[derive(Default)]
struct Cumulative {
    sum: f64,
    correction: f64,
}

impl Cumulative {
    fn add(&mut self, mass: f64) -> Result<f64, SourceSamplingError> {
        if mass == 0.0 {
            return Ok(self.sum);
        }
        let y = mass - self.correction;
        let next = self.sum + y;
        if !mass.is_finite()
            || mass < 0.0
            || !next.is_finite()
            || (mass > 0.0 && next <= self.sum)
            || next < self.sum
        {
            return Err(SourceSamplingError::NumericalRange);
        }
        self.correction = (next - self.sum) - y;
        self.sum = next;
        Ok(next)
    }
}

#[cfg(test)]
mod tests;
