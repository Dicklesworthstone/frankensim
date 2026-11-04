//! Physical spectral detector channels, independent of a color observer.
//!
//! The caller supplies external quantum efficiency in collected electrons per
//! incident photon, without gain, on an original vacuum-wavelength grid. The
//! piecewise-linear response is retained exactly; it is never extrapolated.
//! A common proposal proportional to the sum of all channel responses lets
//! every traced wavelength contribute to every channel. Preparation integrates
//! that proposal analytically; queries allocate nothing.
//!
//! Photon conversion uses E = h*c/lambda_vac and the exact SI definitions of
//! h and c: <https://www.bipm.org/en/measurement-units/si-defining-constants>.
//! External QE and spectral power responsivity follow NIST SP 250-41 (2008),
//! <https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication250-41e2008.pdf>.
//! These authored response curves do not infer detector calibration, optical
//! etendue, exposure, gain, noise, saturation or a measured material response.

use fs_math::{dd::Dd, det};

pub const MAX_SENSOR_CHANNELS: usize = 16;
pub const MAX_SENSOR_KNOTS: usize = 1024;
const PHOTONS_PER_JOULE_PER_NM: f64 = 1.0e-9 / (6.626_070_15e-34 * 299_792_458.0);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectralSensorError {
    /// Between 1 and 16 channels and 2 to 1024 ordered, positive finite knots.
    InvalidGrid,
    /// A response is outside [0,1], nonfinite, or all channels are dark.
    InvalidEfficiency,
    /// A query is outside the original closed source band.
    InvalidWavelength,
    /// The inverse-CDF variate must be finite and in [0,1).
    InvalidVariate,
    /// A nonzero probability or weighted response cannot be represented.
    NumericalRange,
    AllocationRefused,
}

impl core::fmt::Display for SpectralSensorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spectral sensor: {self:?}")
    }
}
impl core::error::Error for SpectralSensorError {}

/// Immutable external quantum-efficiency curves on their original shared grid.
/// Zero-response intervals are permitted, including gaps between disjoint
/// channels. A completely dark channel is retained with exact zero response;
/// the complete sensor must have at least one nonzero channel.
#[derive(Clone, Debug)]
pub struct QuantumEfficiencySensor<const CHANNELS: usize> {
    knots: Vec<(f64, [f64; CHANNELS])>,
    total_efficiency: Vec<f64>,
    cumulative_nm: Vec<f64>,
    scale: f64,
    area_nm: f64,
}

/// One physical wavelength with the exact channel/QE importance weights.
/// Divide each accumulated contribution by the number of completed draws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DetectorWavelength<const CHANNELS: usize> {
    wavelength_nm: f64,
    radiance_weights_nm: [f64; CHANNELS],
    electron_weights_nm_per_joule: [f64; CHANNELS],
}

impl<const C: usize> DetectorWavelength<C> {
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }

    /// QE_i(lambda)/p(lambda), in nm. Multiplied by spectral radiance, this
    /// estimates the channel's QE-weighted radiance integral in W/(m² sr).
    pub const fn radiance_weights_nm(&self) -> &[f64; C] {
        &self.radiance_weights_nm
    }

    /// QE_i(lambda)*lambda_vac/(h*c*p(lambda)), in nm/J. Multiplying spectral
    /// radiance estimates collected-electron rate per m² sr, before etendue
    /// and exposure. This is not electrical gain or radiant heat deposition.
    pub const fn electron_weights_nm_per_joule(&self) -> &[f64; C] {
        &self.electron_weights_nm_per_joule
    }
}

impl<const C: usize> QuantumEfficiencySensor<C> {
    /// Retain `(vacuum_wavelength_nm, [electrons_per_photon; C])` knots.
    /// Preparation costs O(C*K), stores O(C*K) numbers, and admits allocation
    /// only after checking the explicit channel/knot caps and source values.
    pub fn try_new(knots: Vec<(f64, [f64; C])>) -> Result<Self, SpectralSensorError> {
        if C == 0
            || C > MAX_SENSOR_CHANNELS
            || !(2..=MAX_SENSOR_KNOTS).contains(&knots.len())
            || knots.iter().any(|(w, _)| !w.is_finite() || *w <= 0.0)
            || knots.windows(2).any(|pair| pair[0].0 >= pair[1].0)
        {
            return Err(SpectralSensorError::InvalidGrid);
        }
        if knots
            .iter()
            .flat_map(|(_, q)| q)
            .any(|q| !q.is_finite() || !(0.0..=1.0).contains(q))
        {
            return Err(SpectralSensorError::InvalidEfficiency);
        }
        let mut total_efficiency = Vec::new();
        let mut cumulative_nm = Vec::new();
        total_efficiency
            .try_reserve_exact(knots.len())
            .map_err(|_| SpectralSensorError::AllocationRefused)?;
        cumulative_nm
            .try_reserve_exact(knots.len())
            .map_err(|_| SpectralSensorError::AllocationRefused)?;
        for (_, q) in &knots {
            total_efficiency.push(
                q.iter()
                    .fold(Dd::ZERO, |sum, &q| sum + Dd::from_f64(q))
                    .to_f64(),
            );
        }
        let scale = total_efficiency.iter().copied().fold(0.0, f64::max);
        if scale == 0.0 {
            return Err(SpectralSensorError::InvalidEfficiency);
        }
        // Scale before integration so weak but finite detectors retain their
        // wavelength distribution without underflowing its normalization.
        let mut area = Dd::ZERO;
        cumulative_nm.push(0.0);
        for i in 1..knots.len() {
            let a = total_efficiency[i - 1] / scale;
            let b = total_efficiency[i] / scale;
            if (total_efficiency[i - 1] > 0.0 && a == 0.0)
                || (total_efficiency[i] > 0.0 && b == 0.0)
            {
                return Err(SpectralSensorError::NumericalRange);
            }
            let segment = (0.5 * a + 0.5 * b) * (knots[i].0 - knots[i - 1].0);
            if !segment.is_finite() || ((a > 0.0 || b > 0.0) && segment == 0.0) {
                return Err(SpectralSensorError::NumericalRange);
            }
            area = area + Dd::from_f64(segment);
            let end = area.to_f64();
            if !end.is_finite() || (segment > 0.0 && end <= cumulative_nm[i - 1]) {
                return Err(SpectralSensorError::NumericalRange);
            }
            cumulative_nm.push(end);
        }
        let area_nm = area.to_f64();
        if area_nm <= 0.0 {
            return Err(SpectralSensorError::NumericalRange);
        }
        Ok(Self {
            knots,
            total_efficiency,
            cumulative_nm,
            scale,
            area_nm,
        })
    }

    pub fn knots(&self) -> &[(f64, [f64; C])] {
        &self.knots
    }

    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        (self.knots[0].0, self.knots[self.knots.len() - 1].0)
    }

    pub fn quantum_efficiency(&self, wavelength_nm: f64) -> Result<[f64; C], SpectralSensorError> {
        let (lo, hi) = self.wavelength_bounds_nm();
        if !wavelength_nm.is_finite() || !(lo..=hi).contains(&wavelength_nm) {
            return Err(SpectralSensorError::InvalidWavelength);
        }
        let upper = self.knots.partition_point(|(w, _)| *w < wavelength_nm);
        if upper == 0 || self.knots[upper].0 == wavelength_nm {
            return Ok(self.knots[upper].1);
        }
        let a = self.knots[upper - 1];
        let b = self.knots[upper];
        let t = (wavelength_nm - a.0) / (b.0 - a.0);
        let efficiency: [f64; C] =
            core::array::from_fn(|channel| (1.0 - t) * a.1[channel] + t * b.1[channel]);
        // This query is strictly inside a source interval. A positive endpoint
        // therefore gives a positive response, even if both weighted terms
        // underflow when evaluated separately. Exact source knots above keep
        // their original zeros; a lost positive interpolation must refuse.
        if efficiency
            .iter()
            .enumerate()
            .any(|(channel, &q)| q == 0.0 && (a.1[channel] > 0.0 || b.1[channel] > 0.0))
        {
            return Err(SpectralSensorError::NumericalRange);
        }
        Ok(efficiency)
    }

    /// Invert the represented sum-of-QE CDF in O(log K + C), without an FFT,
    /// spectral quadrature, rejection loop or query allocation. Zero-response
    /// gaps have no probability. At an exactly zero-density endpoint the
    /// channel/proposal ratio uses its analytic one-sided limit.
    pub fn sample_wavelength(&self, u: f64) -> Result<DetectorWavelength<C>, SpectralSensorError> {
        if !u.is_finite() || !(0.0..1.0).contains(&u) {
            return Err(SpectralSensorError::InvalidVariate);
        }
        let target = u * self.area_nm;
        let upper = self.cumulative_nm.partition_point(|&area| area <= target);
        if upper == 0 || upper >= self.knots.len() {
            return Err(SpectralSensorError::NumericalRange);
        }
        let a = self.knots[upper - 1];
        let b = self.knots[upper];
        let local_u = (target - self.cumulative_nm[upper - 1])
            / (self.cumulative_nm[upper] - self.cumulative_nm[upper - 1]);
        let local_scale = self.total_efficiency[upper - 1].max(self.total_efficiency[upper]);
        let qa = self.total_efficiency[upper - 1] / local_scale;
        let qb = self.total_efficiency[upper] / local_scale;
        let t = if local_u == 0.0 {
            0.0
        } else {
            // Rationalized inverse of a*t + (b-a)*t²/2. Interpolate squared
            // endpoints to avoid subtracting nearly equal terms on flat QE.
            local_u * (qa + qb) / (qa + det::sqrt((1.0 - local_u) * qa * qa + local_u * qb * qb))
        };
        if !t.is_finite() || !(0.0..=1.0).contains(&t) {
            return Err(SpectralSensorError::NumericalRange);
        }
        let wavelength_nm = a.0 + t * (b.0 - a.0);
        let density = (1.0 - t) * qa + t * qb;
        let ratio: [f64; C] = core::array::from_fn(|channel| {
            if density == 0.0 {
                if t == 0.0 {
                    b.1[channel] / self.total_efficiency[upper]
                } else {
                    a.1[channel] / self.total_efficiency[upper - 1]
                }
            } else {
                ((1.0 - t) * (a.1[channel] / local_scale) + t * (b.1[channel] / local_scale))
                    / density
            }
        });
        // Dividing a positive subnormal channel by the proposal scale can
        // round it to zero before the band width restores a representable
        // weight. Retain authored dark channels, but never report this loss
        // of a positive source response as darkness. At a zero-density
        // endpoint, positivity belongs to the same one-sided limit as above.
        for channel in 0..C {
            let positive_response = if density == 0.0 {
                if t == 0.0 {
                    b.1[channel] > 0.0
                } else {
                    a.1[channel] > 0.0
                }
            } else {
                (t < 1.0 && a.1[channel] > 0.0) || (t > 0.0 && b.1[channel] > 0.0)
            };
            if positive_response && ratio[channel] == 0.0 {
                return Err(SpectralSensorError::NumericalRange);
            }
        }
        let radiance_weights_nm = ratio.map(|ratio| (ratio * self.scale) * self.area_nm);
        let photon_conversion = wavelength_nm * PHOTONS_PER_JOULE_PER_NM;
        let electron_weights_nm_per_joule =
            radiance_weights_nm.map(|weight| weight * photon_conversion);
        if !wavelength_nm.is_finite()
            || !photon_conversion.is_finite()
            || radiance_weights_nm
                .iter()
                .chain(&electron_weights_nm_per_joule)
                .any(|v| !v.is_finite() || *v < 0.0)
            || ratio
                .iter()
                .zip(&radiance_weights_nm)
                .any(|(&q, &r)| q > 0.0 && r == 0.0)
            || radiance_weights_nm
                .iter()
                .zip(&electron_weights_nm_per_joule)
                .any(|(&r, &e)| r > 0.0 && e == 0.0)
        {
            return Err(SpectralSensorError::NumericalRange);
        }
        Ok(DetectorWavelength {
            wavelength_nm,
            radiance_weights_nm,
            electron_weights_nm_per_joule,
        })
    }
}

#[cfg(test)]
mod tests;
