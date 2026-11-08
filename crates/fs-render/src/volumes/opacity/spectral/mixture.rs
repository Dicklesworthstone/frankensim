//! Prepare point-sampled mass mixtures once for the indexed opacity kernel.
//!
//! This is the same independent-constituent mass-weighted absorption law used
//! by native volume transport. No refractive-index, phase or temperature
//! averaging is implied; no mass basis or material data is inferred.

use super::{DensityAbsorptionSpectrum, DvrError, Sum};

/// Hard bounds on preparation of a reusable mixture absorption palette.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AbsorptionMixtureSettings {
    /// Maximum original constituent spectra, including zero-weight entries.
    pub max_components: usize,
    /// Maximum prepared mixture rows.
    pub max_mixtures: usize,
    /// Maximum mixtures times constituents times supplied wavelengths.
    /// Zero weights and zero absorption still count toward this work limit.
    pub max_spectral_evaluations: u64,
}

/// Prepare mass-weighted absorption spectra for the original indexed renderer.
///
/// Constituent coefficients are component-major, wavelength-fastest, in m²/kg
/// of that constituent at the supplied common vacuum wavelengths. Weights are
/// mixture-major, component-fastest, in kg of constituent per kg of the density
/// field's declared reference mass. The returned table is mixture-major,
/// wavelength-fastest, in m²/kg of that same reference mass. It can be borrowed
/// directly as IndexedDensityAbsorptionSpectrum::mass_absorption_m2_kg.
///
/// For total mixture density, the caller supplies mass fractions that sum to
/// one. For carrier-relative density, supplied loadings need not sum to one.
/// This operation never normalizes weights or changes their density basis.
/// Every coefficient and weight must be finite and nonnegative, including
/// unused components. An all-zero recipe explicitly prepares zero absorption.
///
/// At each supplied wavelength the model is exactly sum(weight * coefficient).
/// It assumes independently additive absorption at the declared state; it does
/// not account for concentration-dependent optical interactions. Homogeneous
/// composition within a cell is distinct from unresolved spatial coverage:
/// the existing renderer still averages transmitted intensities across rays.
///
/// Preparation allocates only the final coefficient table. No constituent
/// evaluation or mixture allocation is added to the rendering hot loop. All
/// preparation work and output size are admitted before unbounded input scans
/// or allocation. Positive products that round to zero and overflowing
/// coefficients are refused rather than losing absorption that a subsequent
/// density/path scale could make observable. Input wavelengths are not
/// interpolated, resampled or assigned measurement authority.
pub fn prepare_density_absorption_mixtures(
    wavelengths_nm: &[f64],
    component_absorption_m2_kg: &[f64],
    mixture_mass_weights: &[f64],
    settings: AbsorptionMixtureSettings,
) -> Result<Vec<f64>, DvrError> {
    let lanes = wavelengths_nm.len();
    DensityAbsorptionSpectrum {
        wavelengths_nm,
        mass_absorption_m2_kg: component_absorption_m2_kg
            .get(..lanes)
            .ok_or(DvrError::InvalidAbsorptionSpectrum)?,
        density_unit_kg_m3: 1.0,
        length_unit_m: 1.0,
    }
    .validate()?;
    if component_absorption_m2_kg.len() % lanes != 0 {
        return Err(DvrError::InvalidAbsorptionSpectrum);
    }
    let components = component_absorption_m2_kg.len() / lanes;
    if components > settings.max_components {
        return Err(DvrError::AbsorptionComponentBudgetExceeded {
            requested: components,
            limit: settings.max_components,
        });
    }
    if mixture_mass_weights.is_empty() || mixture_mass_weights.len() % components != 0 {
        return Err(DvrError::InvalidAbsorptionMixtureWeights);
    }
    let mixtures = mixture_mass_weights.len() / components;
    if mixtures > settings.max_mixtures {
        return Err(DvrError::AbsorptionMaterialBudgetExceeded {
            requested: mixtures,
            limit: settings.max_mixtures,
        });
    }
    let evaluations = (mixtures as u128)
        .checked_mul(components as u128)
        .and_then(|value| value.checked_mul(lanes as u128))
        .and_then(|value| u64::try_from(value).ok())
        .ok_or(DvrError::SpectralEvaluationCountOverflow)?;
    if evaluations > settings.max_spectral_evaluations {
        return Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: evaluations,
            limit: settings.max_spectral_evaluations,
        });
    }
    let samples = mixtures
        .checked_mul(lanes)
        .ok_or(DvrError::SpectralEvaluationCountOverflow)?;
    if component_absorption_m2_kg
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(DvrError::InvalidAbsorptionSpectrum);
    }
    if mixture_mass_weights
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(DvrError::InvalidAbsorptionMixtureWeights);
    }
    let mut palette = Vec::new();
    palette.try_reserve_exact(samples).map_err(|_| {
        DvrError::AbsorptionMixtureAllocationRefused {
            coefficients: samples,
        }
    })?;
    for recipe in mixture_mass_weights.chunks_exact(components) {
        for lane in 0..lanes {
            let mut total = Sum::default();
            for (component, &weight) in recipe.iter().enumerate() {
                let coefficient = component_absorption_m2_kg[component * lanes + lane];
                let contribution = weight * coefficient;
                if !contribution.is_finite()
                    || (weight > 0.0 && coefficient > 0.0 && contribution == 0.0)
                {
                    return Err(DvrError::UnrepresentableAbsorptionMixture);
                }
                total.add(contribution);
            }
            if !total.value.is_finite() {
                return Err(DvrError::UnrepresentableAbsorptionMixture);
            }
            palette.push(total.value);
        }
    }
    Ok(palette)
}

#[cfg(test)]
mod tests;
