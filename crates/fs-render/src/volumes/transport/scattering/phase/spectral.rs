//! Spectral angular laws and one matched proposal for every packet lane.

use super::{PhaseFunction, PhaseSample, Vec3, VolumeTransportError, unit_variate};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use fs_blake3::DomainHasher;
use std::sync::Arc;

/// Maximum retained wavelength samples in one angular law.
pub const MAX_PHASE_WAVELENGTHS: usize = 16;

/// Explicit wavelength discretization of normalized, reciprocal phase laws.
///
/// Each interval is the convex mixture of its two endpoint laws. This exactly
/// normalizes the retained interpolant, not an underlying continuous Mie or
/// measured law. Callers must refine wavelength and angular data independently.
/// Nested spectral tables are refused; storage is shared outside path loops.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralPhaseTable {
    knots: Arc<[(f64, PhaseFunction)]>,
}

/// Wavelength-specific physical response and the complete sampling density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhasePacketEvaluation<const N: usize> {
    /// Physical phase probability in sr⁻¹ for each retained wavelength.
    pub value: [f64; N],
    /// Uniform mixture of every lane's phase law, in sr⁻¹.
    pub pdf: f64,
}

/// One direction shared by all wavelengths, with unbiased spectral weights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhasePacketSample<const N: usize> {
    /// Unit outgoing travel direction.
    pub direction: Vec3,
    /// Complete packet-mixture PDF, also used by direct-light MIS.
    pub pdf: f64,
    /// Per-lane `phase(lambda, theta) / pdf`; multiply continuation throughput.
    pub weight: [f64; N],
}

impl SpectralPhaseTable {
    /// Admit 2–16 increasing vacuum-wavelength knots covering 380–780 nm.
    pub fn try_new(knots: &[(f64, PhaseFunction)]) -> Result<Self, VolumeTransportError> {
        if !(2..=MAX_PHASE_WAVELENGTHS).contains(&knots.len())
            || knots[0].0 != LAMBDA_MIN
            || knots[knots.len() - 1].0 != LAMBDA_MAX
            || knots.iter().any(|(wavelength, phase)| {
                !wavelength.is_finite() || matches!(phase, PhaseFunction::SpectralTable(_))
            })
            || knots.windows(2).any(|pair| pair[0].0 >= pair[1].0)
        {
            return Err(VolumeTransportError::Phase);
        }
        for (_, phase) in knots {
            phase.validate()?;
        }
        Ok(Self {
            knots: Arc::from(knots),
        })
    }

    /// Original wavelength knots and admitted achromatic laws.
    #[must_use]
    pub fn knots(&self) -> &[(f64, PhaseFunction)] {
        &self.knots
    }

    fn interval(&self, wavelength: f64) -> (&PhaseFunction, &PhaseFunction, f64) {
        let upper = self
            .knots
            .partition_point(|(lambda, _)| *lambda < wavelength)
            .clamp(1, self.knots.len() - 1);
        let (left, a) = &self.knots[upper - 1];
        let (right, b) = &self.knots[upper];
        (a, b, (wavelength - left) / (right - left))
    }

    pub(super) fn hash_into(&self, hash: &mut DomainHasher) {
        hash.update(&(self.knots.len() as u64).to_le_bytes());
        for (wavelength, phase) in self.knots.iter() {
            hash.update(&wavelength.to_bits().to_le_bytes());
            phase.hash_into(hash);
        }
    }
}

impl PhaseFunction {
    /// Retain the declared wavelength/angular discretization without fitting HG.
    pub fn spectral_table(knots: &[(f64, Self)]) -> Result<Self, VolumeTransportError> {
        SpectralPhaseTable::try_new(knots).map(Self::SpectralTable)
    }

    /// Physical phase response at a vacuum wavelength in the visible band.
    pub fn eval_at(&self, wavelength_nm: f64, cos_theta: f64) -> Result<f64, VolumeTransportError> {
        validate_wavelength(wavelength_nm)?;
        match self {
            Self::SpectralTable(table) => {
                let (a, b, t) = table.interval(wavelength_nm);
                // Exact endpoints preserve the original angular law's bits.
                if t == 0.0 {
                    return a.eval_cos(cos_theta);
                }
                if t == 1.0 {
                    return b.eval_cos(cos_theta);
                }
                Ok((1.0 - t) * a.eval_cos(cos_theta)? + t * b.eval_cos(cos_theta)?)
            }
            _ => self.eval_cos(cos_theta),
        }
    }

    /// Evaluate physical lane responses and their complete mixture PDF.
    /// Achromatic laws retain their original scalar arithmetic exactly.
    pub fn eval_packet<const N: usize>(
        &self,
        wavelengths_nm: &[f64; N],
        cos_theta: f64,
    ) -> Result<PhasePacketEvaluation<N>, VolumeTransportError> {
        validate_packet(wavelengths_nm)?;
        if !matches!(self, Self::SpectralTable(_)) {
            let pdf = self.eval_cos(cos_theta)?;
            return Ok(PhasePacketEvaluation {
                value: [pdf; N],
                pdf,
            });
        }
        let mut value = [0.0; N];
        let mut pdf = 0.0;
        for (lane, &wavelength) in wavelengths_nm.iter().enumerate() {
            value[lane] = self.eval_at(wavelength, cos_theta)?;
            pdf += value[lane] / N as f64;
        }
        if !pdf.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(PhasePacketEvaluation { value, pdf })
    }

    /// Sample the uniform mixture of the complete wavelength packet.
    ///
    /// Conditional CDF remapping uses one variate for lane selection, spectral
    /// endpoint selection, and angular sampling. The second supplies azimuth.
    /// No allocation or cloned phase table occurs at a scattering vertex.
    pub fn sample_packet<const N: usize>(
        &self,
        wavelengths_nm: &[f64; N],
        incoming_travel: Vec3,
        u_cos: f64,
        u_azimuth: f64,
    ) -> Result<PhasePacketSample<N>, VolumeTransportError> {
        validate_packet(wavelengths_nm)?;
        if !unit_variate(u_cos) || !unit_variate(u_azimuth) {
            return Err(VolumeTransportError::Sample);
        }
        let Self::SpectralTable(table) = self else {
            let PhaseSample { direction, pdf } = self.sample(incoming_travel, u_cos, u_azimuth)?;
            return Ok(PhasePacketSample {
                direction,
                pdf,
                weight: [1.0; N],
            });
        };
        let scaled = u_cos * N as f64;
        let lane = (scaled as usize).min(N - 1);
        let u = unit_below_one(scaled - lane as f64);
        let (a, b, t) = table.interval(wavelengths_nm[lane]);
        let (phase, u) = if t == 0.0 {
            (a, u)
        } else if t == 1.0 {
            (b, u)
        } else if u < 1.0 - t {
            (a, unit_below_one(u / (1.0 - t)))
        } else {
            (b, unit_below_one((u - (1.0 - t)) / t))
        };
        let sampled = phase.sample(incoming_travel, u, u_azimuth)?;
        let evaluation =
            self.eval_packet(wavelengths_nm, incoming_travel.dot(sampled.direction))?;
        if !(evaluation.pdf > 0.0) {
            return Err(VolumeTransportError::Unrepresentable);
        }
        let weight = evaluation.value.map(|value| value / evaluation.pdf);
        if weight.iter().any(|value| !value.is_finite()) {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(PhasePacketSample {
            direction: sampled.direction,
            pdf: evaluation.pdf,
            weight,
        })
    }
}

fn unit_below_one(value: f64) -> f64 {
    // Only conditional-probability endpoint roundoff is adjusted, never a PDF.
    value.clamp(0.0, fs_math::next_down(1.0))
}

fn validate_wavelength(wavelength: f64) -> Result<(), VolumeTransportError> {
    if !wavelength.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength) {
        return Err(VolumeTransportError::Wavelength);
    }
    Ok(())
}

fn validate_packet<const N: usize>(wavelengths: &[f64; N]) -> Result<(), VolumeTransportError> {
    if N == 0 {
        return Err(VolumeTransportError::Sample);
    }
    for &wavelength in wavelengths {
        validate_wavelength(wavelength)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
