//! Spectral angular laws and one matched proposal for every packet lane.

use super::{PhaseFunction, PhaseSample, Vec3, VolumeTransportError, unit_variate};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use fs_blake3::{ContentHash, DomainHasher};
use std::sync::Arc;

/// Maximum retained wavelength samples in one angular law.
pub const MAX_PHASE_WAVELENGTHS: usize = 16;

/// Explicit wavelength discretization of normalized, reciprocal phase laws.
///
/// Each interval is the convex mixture of its two endpoint laws. This exactly
/// normalizes the retained interpolant, not an underlying continuous Mie or
/// measured law. Callers must refine wavelength and angular data independently.
/// Spectral tables and scattering mixtures cannot be endpoint laws; storage is
/// shared outside path loops.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralPhaseTable {
    knots: Arc<[(f64, PhaseFunction)]>,
    source_identity: Option<ContentHash>,
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
        {
            return Err(VolumeTransportError::Phase);
        }
        Self::admit(knots, None)
    }

    /// Admit an identified source's 2–16 angular laws over its exact positive
    /// wavelength domain. No phase response is extrapolated beyond the knots;
    /// the source identity does not assert experimental qualification.
    pub fn source_table(
        knots: &[(f64, PhaseFunction)],
        source_identity: ContentHash,
    ) -> Result<Self, VolumeTransportError> {
        if source_identity == ContentHash([0; 32]) {
            return Err(VolumeTransportError::SourceIdentity);
        }
        Self::admit(knots, Some(source_identity))
    }

    fn admit(
        knots: &[(f64, PhaseFunction)],
        source_identity: Option<ContentHash>,
    ) -> Result<Self, VolumeTransportError> {
        if !(2..=MAX_PHASE_WAVELENGTHS).contains(&knots.len())
            || knots.iter().any(|(wavelength, phase)| {
                !wavelength.is_finite()
                    || *wavelength <= 0.0
                    || matches!(
                        phase,
                        PhaseFunction::SpectralTable(_) | PhaseFunction::SpectralMixture(_)
                    )
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
            source_identity,
        })
    }

    /// Original wavelength knots and admitted achromatic laws.
    #[must_use]
    pub fn knots(&self) -> &[(f64, PhaseFunction)] {
        &self.knots
    }

    /// Inclusive vacuum-wavelength support, in nm.
    #[must_use]
    pub fn wavelength_bounds_nm(&self) -> (f64, f64) {
        (self.knots[0].0, self.knots[self.knots.len() - 1].0)
    }

    /// Original source identity, when this law was explicitly source-bound.
    #[must_use]
    pub const fn source_identity(&self) -> Option<ContentHash> {
        self.source_identity
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
        if let Some(source) = self.source_identity {
            hash.update(b"source-table");
            hash.update(source.as_bytes());
        }
    }
}

impl PhaseFunction {
    /// Retain the declared wavelength/angular discretization without fitting HG.
    pub fn spectral_table(knots: &[(f64, Self)]) -> Result<Self, VolumeTransportError> {
        SpectralPhaseTable::try_new(knots).map(Self::SpectralTable)
    }

    /// Retain explicitly sourced angular laws over their exact UV/visible/IR
    /// support, with the same normalized mixture used for response and sampling.
    pub fn source_spectral_table(
        knots: &[(f64, Self)],
        source_identity: ContentHash,
    ) -> Result<Self, VolumeTransportError> {
        SpectralPhaseTable::source_table(knots, source_identity).map(Self::SpectralTable)
    }

    pub(in crate::volumes::transport) fn validate_wavelength(
        &self,
        wavelength: f64,
    ) -> Result<(), VolumeTransportError> {
        if !wavelength.is_finite() || wavelength <= 0.0 {
            return Err(VolumeTransportError::Wavelength);
        }
        let bounds = match self {
            Self::SpectralTable(table) => Some(table.wavelength_bounds_nm()),
            Self::SpectralMixture(mixture) => Some(mixture.wavelength_bounds_nm()),
            _ => None,
        };
        if let Some((lower, upper)) = bounds {
            if !(lower..=upper).contains(&wavelength) {
                return Err(VolumeTransportError::Wavelength);
            }
        }
        Ok(())
    }

    /// Physical phase response at an admitted positive vacuum wavelength.
    pub fn eval_at(&self, wavelength_nm: f64, cos_theta: f64) -> Result<f64, VolumeTransportError> {
        self.validate_wavelength(wavelength_nm)?;
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
            Self::SpectralMixture(mixture) => mixture.eval_at(wavelength_nm, cos_theta),
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
        validate_packet(self, wavelengths_nm)?;
        if !matches!(self, Self::SpectralTable(_) | Self::SpectralMixture(_)) {
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
    /// endpoint/component selection, and angular sampling. The second supplies
    /// azimuth. Mixture components retain their own complete spectral laws.
    /// No allocation or cloned phase table occurs at a scattering vertex.
    pub fn sample_packet<const N: usize>(
        &self,
        wavelengths_nm: &[f64; N],
        incoming_travel: Vec3,
        u_cos: f64,
        u_azimuth: f64,
    ) -> Result<PhasePacketSample<N>, VolumeTransportError> {
        validate_packet(self, wavelengths_nm)?;
        if !unit_variate(u_cos) || !unit_variate(u_azimuth) {
            return Err(VolumeTransportError::Sample);
        }
        if !matches!(self, Self::SpectralTable(_) | Self::SpectralMixture(_)) {
            let PhaseSample { direction, pdf } = self.sample(incoming_travel, u_cos, u_azimuth)?;
            return Ok(PhasePacketSample {
                direction,
                pdf,
                weight: [1.0; N],
            });
        }
        let scaled = u_cos * N as f64;
        let lane = (scaled as usize).min(N - 1);
        let u = unit_below_one(scaled - lane as f64);
        let direction =
            self.sample_direction_at(wavelengths_nm[lane], incoming_travel, u, u_azimuth)?;
        let evaluation = self.eval_packet(wavelengths_nm, incoming_travel.dot(direction))?;
        if !(evaluation.pdf > 0.0) {
            return Err(VolumeTransportError::Unrepresentable);
        }
        let weight = evaluation.value.map(|value| value / evaluation.pdf);
        if weight.iter().any(|value| !value.is_finite()) {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(PhasePacketSample {
            direction,
            pdf: evaluation.pdf,
            weight,
        })
    }

    // Inputs are admitted by sample_packet. Mixtures are flat, and spectral
    // table endpoints are achromatic, so this dispatch has bounded depth.
    pub(super) fn sample_direction_at(
        &self,
        wavelength_nm: f64,
        incoming_travel: Vec3,
        u_cos: f64,
        u_azimuth: f64,
    ) -> Result<Vec3, VolumeTransportError> {
        match self {
            Self::SpectralTable(table) => {
                let (a, b, t) = table.interval(wavelength_nm);
                let (phase, u) = if t == 0.0 {
                    (a, u_cos)
                } else if t == 1.0 {
                    (b, u_cos)
                } else if u_cos < 1.0 - t {
                    (a, unit_below_one(u_cos / (1.0 - t)))
                } else {
                    (b, unit_below_one((u_cos - (1.0 - t)) / t))
                };
                Ok(phase.sample(incoming_travel, u, u_azimuth)?.direction)
            }
            Self::SpectralMixture(mixture) => {
                mixture.sample_direction_at(wavelength_nm, incoming_travel, u_cos, u_azimuth)
            }
            _ => Ok(self.sample(incoming_travel, u_cos, u_azimuth)?.direction),
        }
    }
}

pub(super) fn unit_below_one(value: f64) -> f64 {
    // Only conditional-probability endpoint roundoff is adjusted, never a PDF.
    value.clamp(0.0, fs_math::next_down(1.0))
}

fn validate_packet<const N: usize>(
    phase: &PhaseFunction,
    wavelengths: &[f64; N],
) -> Result<(), VolumeTransportError> {
    if N == 0 {
        return Err(VolumeTransportError::Sample);
    }
    for &wavelength in wavelengths {
        phase.validate_wavelength(wavelength)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
