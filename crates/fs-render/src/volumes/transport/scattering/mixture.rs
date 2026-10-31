//! Prepared independent-scatterer sums with an exact molecular contribution.

use super::{
    MassAbsorptionSpectrum, MassScatteringSpectrum, PhaseFunction, ScatteringMixtureComponent,
    ScatteringOptics, VolumeTransportError,
};
use crate::volumes::transport::OpticalMixtureBudget;
use crate::volumes::transport::mixture::{
    PositiveSum, validate_count, validate_weight, weighted_product,
};
use fs_blake3::DomainHasher;
use fs_exec::Cx;

impl MassScatteringSpectrum {
    /// Compile an exact weighted sum of the retained constituent spectra.
    ///
    /// Weights are explicit kg-constituent/kg-reference-mass, as for
    /// [`MassAbsorptionSpectrum::weighted_sum`]. All domains must match. Linear
    /// terms use their exact knot union (at most 16); all molecular inverse-fourth
    /// terms combine analytically. Preparation never samples that molecular
    /// term into a linear table. The result remains fixed-size and `Copy`.
    pub fn weighted_sum(
        cx: &Cx<'_>,
        components: &[(f64, Self)],
        budget: OpticalMixtureBudget,
    ) -> Result<Self, VolumeTransportError> {
        validate_count(cx, components.len(), budget)?;
        let bounds = components[0].1.wavelength_bounds_nm();
        let reference = bounds.0;
        let mut linear = Vec::new();
        linear
            .try_reserve_exact(components.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        let mut rayleigh = PositiveSum::default();
        let mut hash = DomainHasher::new("org.frankensim.render.scattering-mixture.v1");
        hash.update(&(components.len() as u64).to_le_bytes());
        for &(weight, spectrum) in components {
            super::checkpoint(cx)?;
            validate_weight(weight)?;
            if spectrum.wavelength_bounds_nm() != bounds {
                return Err(VolumeTransportError::Wavelength);
            }
            hash.update(&weight.to_bits().to_le_bytes());
            spectrum.hash_into(&mut hash);
            linear.push((
                if spectrum.rayleigh_reference_nm.is_none() {
                    weight
                } else {
                    0.0
                },
                spectrum.spectrum,
            ));
            if weight == 0.0 {
                continue;
            }
            if let Some(original_reference) = spectrum.rayleigh_reference_nm {
                let ratio = original_reference / reference;
                let square = ratio * ratio;
                let coefficient = spectrum.spectrum.at(reference)? * square * square;
                rayleigh.add(weighted_product(weight, coefficient)?)?;
            }
            if let Some((original_reference, coefficient)) = spectrum.rayleigh_addend {
                let ratio = original_reference / reference;
                let square = ratio * ratio;
                rayleigh.add(weighted_product(weight, coefficient * square * square)?)?;
            }
        }
        if components.len() == 1 && components[0].0 == 1.0 {
            return Ok(components[0].1);
        }
        let mut spectrum = MassAbsorptionSpectrum::weighted_sum(cx, &linear, budget)?;
        // Bind the ORIGINAL scattering laws, including their power-law references,
        // rather than only the zeroed linear placeholders used during compilation.
        spectrum.source_identity = Some(hash.finalize());
        let coefficient = rayleigh.finish()?;
        let result = Self {
            spectrum,
            rayleigh_reference_nm: None,
            rayleigh_addend: (coefficient > 0.0).then_some((reference, coefficient)),
        };
        result.at(bounds.0)?;
        result.at(bounds.1)?;
        Ok(result)
    }
}

impl ScatteringOptics {
    /// Prepare extinction and phase from the SAME constituent contributions.
    ///
    /// The physical phase is `sum(w_i*kappa_s_i(lambda)*p_i(lambda,mu)) /
    /// sum(w_i*kappa_s_i(lambda))`; it is neither mass-weighted alone nor fitted
    /// to an average asymmetry. Original angular/spectral models remain retained.
    /// Independent-scatterer superposition does not model dependent scattering,
    /// coherent effective indices, particle correlations, or interfaces.
    pub fn weighted_sum(
        cx: &Cx<'_>,
        components: &[ScatteringMixtureComponent],
        budget: OpticalMixtureBudget,
    ) -> Result<Self, VolumeTransportError> {
        validate_count(cx, components.len(), budget)?;
        let phase = PhaseFunction::scattering_mixture(cx, components, budget)?;
        let mut spectra = Vec::new();
        spectra
            .try_reserve_exact(components.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        spectra.extend(
            components
                .iter()
                .map(|component| (component.mass_weight, component.scattering)),
        );
        let scattering = MassScatteringSpectrum::weighted_sum(cx, &spectra, budget)?;
        Ok(Self { scattering, phase })
    }
}

#[cfg(test)]
mod tests;
