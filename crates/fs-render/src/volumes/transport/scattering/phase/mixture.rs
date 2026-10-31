//! Exact composition of retained independent-particle angular/spectral laws.

use super::super::MassScatteringSpectrum;
use super::spectral::unit_below_one;
use super::{DomainHasher, PhaseFunction, Vec3, VolumeTransportError};
use crate::volumes::transport::{MAX_OPTICAL_MIXTURE_COMPONENTS, OpticalMixtureBudget, checkpoint};
use fs_exec::Cx;
use std::sync::Arc;

/// One independent scattering population relative to a declared density field.
///
/// If that field is total mixture mass density, `mass_weight` is its component
/// mass fraction. With another explicitly declared reference density it is
/// kg of component per kg of reference material. Weights are never silently
/// normalized: the same values must form the medium's scattering coefficient.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringMixtureComponent {
    /// Finite nonnegative component mass loading, including explicit zero.
    pub mass_weight: f64,
    /// Original coefficient per kg of this component, in m²/kg.
    pub scattering: MassScatteringSpectrum,
    /// Original normalized angular law, possibly a wavelength-resolved table.
    pub phase: PhaseFunction,
}

/// Cross-section-weighted phase of coexisting independent scatterers.
///
/// At each wavelength, the normalized weight of component `i` is
/// `mass_weight[i] * kappa_s[i](lambda) / sum_j(mass_weight[j] * kappa_s[j](lambda))`.
/// Original analytic spectra and angular laws are retained without fitting or
/// wavelength resampling. This composes the admitted independent-scattering
/// models; it does not introduce dependent scattering between nearby particles.
///
/// A wavelength with zero total scattering uses the normalized isotropic law.
/// There is no physical collision in that lane, but its normalized proposal
/// remains defined when other lanes of a shared packet do scatter. Mixtures
/// cannot contain mixtures, including inside a spectral-table endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralPhaseMixture {
    components: Arc<[ScatteringMixtureComponent]>,
    wavelength_bounds_nm: (f64, f64),
}

impl SpectralPhaseMixture {
    /// Retain a nonempty, bounded list with identical coefficient domains.
    /// Each phase must cover that domain; no missing response is extrapolated.
    /// All-zero weights are admitted. Cancellation is checked per component,
    /// and immutable storage is shared by clones outside the path loop.
    pub fn try_new(
        cx: &Cx<'_>,
        components: &[ScatteringMixtureComponent],
        budget: OpticalMixtureBudget,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if components.len() > budget.max_components
            || components.len() > MAX_OPTICAL_MIXTURE_COMPONENTS
        {
            return Err(VolumeTransportError::Budget);
        }
        let Some(first) = components.first() else {
            return Err(VolumeTransportError::Mixture);
        };
        let bounds = first.scattering.wavelength_bounds_nm();
        let mut retained = Vec::new();
        retained
            .try_reserve_exact(components.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        for component in components {
            checkpoint(cx)?;
            if !component.mass_weight.is_finite() || component.mass_weight < 0.0 {
                return Err(VolumeTransportError::Mixture);
            }
            if component.scattering.wavelength_bounds_nm() != bounds {
                return Err(VolumeTransportError::Wavelength);
            }
            component.phase.validate()?;
            if matches!(component.phase, PhaseFunction::SpectralMixture(_))
                || matches!(&component.phase, PhaseFunction::SpectralTable(table)
                    if table.knots().iter().any(|(_, phase)| matches!(phase,
                        PhaseFunction::SpectralMixture(_) | PhaseFunction::SpectralTable(_))))
            {
                return Err(VolumeTransportError::Phase);
            }
            component.phase.validate_wavelength(bounds.0)?;
            component.phase.validate_wavelength(bounds.1)?;
            retained.push(component.clone());
        }
        checkpoint(cx)?;
        Ok(Self {
            components: retained.into(),
            wavelength_bounds_nm: bounds,
        })
    }

    /// Original declaration order, component loadings, spectra, and phases.
    #[must_use]
    pub fn components(&self) -> &[ScatteringMixtureComponent] {
        &self.components
    }

    /// Common inclusive coefficient support in vacuum nm.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }

    pub(super) fn hash_into(&self, hash: &mut DomainHasher) {
        hash.update(&(self.components.len() as u64).to_le_bytes());
        for component in self.components.iter() {
            hash.update(&component.mass_weight.to_bits().to_le_bytes());
            component.scattering.hash_into(hash);
            component.phase.hash_into(hash);
        }
    }

    pub(super) fn eval_at(
        &self,
        wavelength_nm: f64,
        cos_theta: f64,
    ) -> Result<f64, VolumeTransportError> {
        let weights = self.weights_at(wavelength_nm)?;
        if weights.active == 0 {
            return PhaseFunction::Isotropic.eval_cos(cos_theta);
        }
        if weights.active == 1 {
            return self.components[weights.last]
                .phase
                .eval_at(wavelength_nm, cos_theta);
        }
        let (mut value, mut correction) = (0.0, 0.0);
        for (component, &probability) in self.components.iter().zip(&weights.probability) {
            if probability == 0.0 {
                continue;
            }
            let phase = component.phase.eval_at(wavelength_nm, cos_theta)?;
            let contribution = probability * phase;
            if phase > 0.0 && contribution == 0.0 {
                return Err(VolumeTransportError::Unrepresentable);
            }
            add(&mut value, &mut correction, contribution);
        }
        let value = value + correction;
        if !value.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(value)
    }

    pub(super) fn sample_direction_at(
        &self,
        wavelength_nm: f64,
        incoming_travel: Vec3,
        u_cos: f64,
        u_azimuth: f64,
    ) -> Result<Vec3, VolumeTransportError> {
        let weights = self.weights_at(wavelength_nm)?;
        if weights.active == 0 {
            return Ok(PhaseFunction::Isotropic
                .sample(incoming_travel, u_cos, u_azimuth)?
                .direction);
        }
        let (mut sum, mut correction, mut previous) = (0.0, 0.0, 0.0);
        let mut selected = None;
        for (index, &probability) in weights.probability[..self.components.len()]
            .iter()
            .enumerate()
        {
            if probability == 0.0 {
                continue;
            }
            add(&mut sum, &mut correction, probability);
            let end = if index == weights.last {
                1.0
            } else {
                sum + correction
            };
            // A positive population must retain sampling support. Never hide
            // cumulative-rounding loss by replacing its physical PDF or mass.
            if !(end > previous && end <= 1.0) {
                return Err(VolumeTransportError::Unrepresentable);
            }
            if selected.is_none() && u_cos < end {
                selected = Some((index, unit_below_one((u_cos - previous) / (end - previous))));
            }
            previous = end;
        }
        let (index, conditional) = selected.ok_or(VolumeTransportError::Unrepresentable)?;
        self.components[index].phase.sample_direction_at(
            wavelength_nm,
            incoming_travel,
            conditional,
            u_azimuth,
        )
    }

    fn weights_at(&self, wavelength_nm: f64) -> Result<ScatteringWeights, VolumeTransportError> {
        let mut weights = ScatteringWeights {
            probability: [0.0; MAX_OPTICAL_MIXTURE_COMPONENTS],
            active: 0,
            last: 0,
        };
        let mut maximum = 0.0_f64;
        for (index, component) in self.components.iter().enumerate() {
            if component.mass_weight == 0.0 {
                continue;
            }
            let scattering = component.scattering.at(wavelength_nm)?;
            if scattering == 0.0 {
                continue;
            }
            let weight = component.mass_weight * scattering;
            if !(weight > 0.0 && weight.is_finite()) {
                return Err(VolumeTransportError::Unrepresentable);
            }
            weights.probability[index] = weight;
            maximum = maximum.max(weight);
            weights.active += 1;
            weights.last = index;
        }
        if weights.active == 0 {
            return Ok(weights);
        }
        // Scaling bounds the sum by the hard component cap, even when an
        // unscaled sum would overflow. Positive components never underflow
        // silently into a different scattering model.
        let (mut total, mut correction) = (0.0, 0.0);
        for probability in &mut weights.probability[..self.components.len()] {
            if *probability == 0.0 {
                continue;
            }
            *probability /= maximum;
            if *probability == 0.0 {
                return Err(VolumeTransportError::Unrepresentable);
            }
            add(&mut total, &mut correction, *probability);
        }
        let total = total + correction;
        for probability in &mut weights.probability[..self.components.len()] {
            if *probability == 0.0 {
                continue;
            }
            *probability /= total;
            if *probability == 0.0 {
                return Err(VolumeTransportError::Unrepresentable);
            }
        }
        Ok(weights)
    }
}

impl PhaseFunction {
    /// Compose independent populations using their wavelength-specific
    /// scattering contributions. The supplied mass loadings must also form
    /// the transport coefficient; they are not number fractions or angular
    /// interpolation weights.
    pub fn scattering_mixture(
        cx: &Cx<'_>,
        components: &[ScatteringMixtureComponent],
        budget: OpticalMixtureBudget,
    ) -> Result<Self, VolumeTransportError> {
        SpectralPhaseMixture::try_new(cx, components, budget).map(Self::SpectralMixture)
    }
}

struct ScatteringWeights {
    probability: [f64; MAX_OPTICAL_MIXTURE_COMPONENTS],
    active: usize,
    last: usize,
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
