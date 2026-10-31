//! Absorption-linked emission from independently thermalized constituents.

use super::mixture::{
    PositiveSum, hash_absorption, validate_count, validate_weight, weighted_product,
};
use super::{MassAbsorptionSpectrum, OpticalMixtureBudget, VolumeTransportError, checkpoint};
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_math::det;
use std::sync::Arc;

/// One original constituent's absorption and accepted thermodynamic state.
///
/// `mass_weight` is kg of this constituent per kg of the density field's
/// reference mass, exactly as in [`MassAbsorptionSpectrum::weighted_sum`].
/// Each absorber obeys Kirchhoff's law at its own temperature. Constituents
/// need not be in thermal equilibrium with one another. Scattering creates
/// neither thermal emission nor heat in this source model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalEmissionComponent {
    /// Nonnegative mass weight; admission never normalizes it.
    pub mass_weight: f64,
    /// Original mass-absorption spectrum in m²/kg of constituent.
    pub absorption: MassAbsorptionSpectrum,
    /// Accepted kelvin temperature. May be absent only for a zero-weight or
    /// identically nonabsorbing constituent; a supplied value must be positive.
    pub temperature_k: Option<f64>,
}

/// Bounded immutable constituent sources paired with their exact optical sum.
///
/// The emission coefficient is `rho * sum(w_i * kappa_ai * B_lambda(T_i))`.
/// The native cell integrator uses the absorption-weighted source radiance,
/// rather than assigning an average temperature to the mixture. Component
/// spectra, weights, and temperatures retain their original declaration order.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalEmissionModel {
    components: Arc<[ThermalEmissionComponent]>,
    absorption: MassAbsorptionSpectrum,
    maximum_temperature_k: Option<f64>,
    common_temperature_k: Option<f64>,
    maximum_absorption_m2_kg: f64,
    identity: ContentHash,
}

impl ThermalEmissionModel {
    /// Prepare a bounded source and its absorption together. Every source must
    /// have the same wavelength endpoints, and the exact linear-knot union
    /// must fit the native spectrum. No temperature or spectral extrapolation
    /// is introduced. Zero-weight entries remain in the source identity.
    pub fn try_new(
        cx: &Cx<'_>,
        components: &[ThermalEmissionComponent],
        budget: OpticalMixtureBudget,
    ) -> Result<Self, VolumeTransportError> {
        validate_count(cx, components.len(), budget)?;
        let mut hash = DomainHasher::new("org.frankensim.render.constituent-thermal-source.v1");
        hash.update(&(components.len() as u64).to_le_bytes());
        let mut optical_components = Vec::with_capacity(components.len());
        let mut maximum_temperature: Option<f64> = None;
        let mut first_temperature = None;
        let mut common_temperature = true;
        for component in components {
            checkpoint(cx)?;
            validate_weight(component.mass_weight)?;
            if component
                .temperature_k
                .is_some_and(|temperature| !temperature.is_finite() || temperature <= 0.0)
            {
                return Err(VolumeTransportError::Temperature);
            }
            let absorbs = component.mass_weight > 0.0
                && component.absorption.knots[..component.absorption.count]
                    .iter()
                    .any(|&(_, coefficient)| coefficient > 0.0);
            if absorbs {
                let temperature = component
                    .temperature_k
                    .ok_or(VolumeTransportError::Temperature)?;
                maximum_temperature = Some(
                    maximum_temperature.map_or(temperature, |previous| previous.max(temperature)),
                );
                if let Some(first) = first_temperature {
                    common_temperature &= first == temperature;
                } else {
                    first_temperature = Some(temperature);
                }
            }
            optical_components.push((component.mass_weight, component.absorption));
            hash.update(&component.mass_weight.to_bits().to_le_bytes());
            hash_absorption(&component.absorption, &mut hash);
            hash.update(&[u8::from(component.temperature_k.is_some())]);
            if let Some(temperature) = component.temperature_k {
                hash.update(&temperature.to_bits().to_le_bytes());
            }
        }
        let absorption = MassAbsorptionSpectrum::weighted_sum(cx, &optical_components, budget)?;
        let maximum_absorption_m2_kg = absorption.knots[..absorption.count]
            .iter()
            .fold(0.0_f64, |maximum, &(_, value)| maximum.max(value));
        checkpoint(cx)?;
        Ok(Self {
            components: Arc::from(components),
            absorption,
            maximum_temperature_k: maximum_temperature,
            common_temperature_k: common_temperature.then_some(first_temperature).flatten(),
            maximum_absorption_m2_kg,
            identity: hash.finalize(),
        })
    }

    /// Sum of the same original absorbers that create this thermal emission.
    #[must_use]
    pub const fn absorption_spectrum(&self) -> MassAbsorptionSpectrum {
        self.absorption
    }

    /// Original constituents, including zero-weight and nonabsorbing entries.
    #[must_use]
    pub fn components(&self) -> &[ThermalEmissionComponent] {
        &self.components
    }

    /// Exact shared wavelength support of all original constituents, vacuum nm.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.absorption.wavelength_bounds_nm()
    }

    /// Hottest active absorber, used only to bound an emission proposal.
    /// Transport always evaluates each original constituent temperature.
    #[must_use]
    pub const fn maximum_temperature_k(&self) -> Option<f64> {
        self.maximum_temperature_k
    }

    /// Maximum combined absorption over the complete piecewise-linear support.
    #[must_use]
    pub const fn maximum_absorption_m2_kg(&self) -> f64 {
        self.maximum_absorption_m2_kg
    }

    /// Original source spectra, component weights, and accepted temperatures.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Absorption-weighted source radiance, W/(m² sr nm), at one wavelength.
    ///
    /// A transparent wavelength has exactly zero source. Equal-temperature
    /// absorbers reduce to the original Planck evaluation without a second
    /// weighted sum. Otherwise division precedes multiplication to avoid an
    /// overflowing `kappa * B`; a logarithmic fallback retains representable
    /// contributions from very small absorption fractions.
    pub fn source_radiance_at(&self, wavelength_nm: f64) -> Result<f64, VolumeTransportError> {
        let absorption = self.absorption.at(wavelength_nm)?;
        if absorption == 0.0 {
            return Ok(0.0);
        }
        if let Some(temperature) = self.common_temperature_k {
            return planck_radiance_w_m2_sr_nm(temperature, wavelength_nm)
                .map_err(|_| VolumeTransportError::Unrepresentable);
        }
        let mut radiance = PositiveSum::default();
        for component in self.components.iter() {
            if component.mass_weight == 0.0 {
                continue;
            }
            let coefficient = component.absorption.at(wavelength_nm)?;
            if coefficient == 0.0 {
                continue;
            }
            let weighted_absorption = weighted_product(component.mass_weight, coefficient)?;
            let source = planck_radiance_w_m2_sr_nm(
                component
                    .temperature_k
                    .ok_or(VolumeTransportError::Temperature)?,
                wavelength_nm,
            )
            .map_err(|_| VolumeTransportError::Unrepresentable)?;
            let mut contribution = (weighted_absorption / absorption) * source;
            if contribution == 0.0 && source > 0.0 {
                contribution =
                    det::exp(det::ln(weighted_absorption) - det::ln(absorption) + det::ln(source));
            }
            radiance.add(contribution)?;
        }
        radiance.finish()
    }
}

#[cfg(test)]
mod tests;
