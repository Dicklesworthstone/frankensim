//! Source-bound bulk absorption at the original accepted material condition.
//!
//! This bridge evaluates the existing optical law at requested wavelengths.
//! It does not replace a complex-index law by interpolating converted mass
//! coefficients, infer absorption from RGB, or evolve the material state.

use super::bulk_optics::{BulkOpticalLaw, accepted_density_property};
use super::{MaterialOpticalBinding, MaterialOpticalError};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use crate::volumes::{DensityAbsorptionSpectrum, MAX_OPACITY_WAVELENGTHS};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::PropertyKey;
use fs_material::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, ResolvedScalarProperty,
    SpectralComplexIndexStatePoint, VisibleDielectricAbsorption, VisibleOpticalStatePoint,
};
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq)]
enum AbsorptionLaw<'a> {
    Visible(BulkOpticalLaw),
    Spectral(&'a SpectralComplexIndexStatePoint),
}

/// Bulk intensity absorption and intrinsic density at one accepted condition.
///
/// The original accepted and optical bundles and selected density receipt are
/// borrowed. Mass absorption is `kappa = alpha / rho_material` in m²/kg of this
/// material, with `alpha = 4*pi*k/lambda_vacuum` for complex-index sources.
/// This density is the material's intrinsic density, not its dispersed mass per
/// bulk volume or a lattice-unit conversion.
///
/// State changes require refreshing the original source claims and rebinding.
/// No linear density/temperature extrapolation, material mixture, scattering,
/// interface reflection or thermal emission is inferred by this conversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialAbsorptionBinding<'a> {
    accepted: &'a ResolvedMaterialStatePoint,
    optical: &'a ResolvedMaterialStatePoint,
    density: &'a ResolvedScalarProperty,
    density_kg_m3: f64,
    law: AbsorptionLaw<'a>,
    source_identity: ContentHash,
    identity: ContentHash,
}

impl<'a> MaterialAbsorptionBinding<'a> {
    /// Bind the original visible conductor or spectral dielectric absorption.
    ///
    /// The complete accepted card, condition and quantity conventions must
    /// match. Density is selected by its full property key and the same strict
    /// scalar kg/m³ admission used by material particle optics. Legacy
    /// RGB-derived absorption is refused; surface finish is not required.
    pub fn try_visible(
        accepted: &'a ResolvedMaterialStatePoint,
        optical: &'a VisibleOpticalStatePoint,
        density_property: &PropertyKey,
    ) -> Result<Self, MaterialAbsorptionError> {
        MaterialOpticalBinding::validate_bulk_state(accepted, optical.resolved())?;
        let absorption = match optical {
            VisibleOpticalStatePoint::Dielectric(value) => Some(value.absorption()),
            VisibleOpticalStatePoint::DielectricSellmeier(value) => Some(value.absorption()),
            _ => None,
        };
        if matches!(
            absorption,
            Some(VisibleDielectricAbsorption::ReferenceRgb { .. })
        ) {
            return Err(MaterialAbsorptionError::SpectralAbsorptionRequired);
        }
        Self::bind(
            accepted,
            optical.resolved(),
            optical.resolved().identity(),
            AbsorptionLaw::Visible(BulkOpticalLaw::from_state(optical)?),
            density_property,
        )
    }

    /// Bind an original finite-band complex-index source, including infrared.
    ///
    /// Every wavelength is evaluated through that source's declared n/k law.
    /// The original closed support and complete wavelength/property-key
    /// interpretation remain intact; no visible-to-infrared extension occurs.
    pub fn try_spectral(
        accepted: &'a ResolvedMaterialStatePoint,
        optical: &'a SpectralComplexIndexStatePoint,
        density_property: &PropertyKey,
    ) -> Result<Self, MaterialAbsorptionError> {
        MaterialOpticalBinding::validate_bulk_state(accepted, optical.resolved())?;
        Self::bind(
            accepted,
            optical.resolved(),
            optical.identity(),
            AbsorptionLaw::Spectral(optical),
            density_property,
        )
    }

    fn bind(
        accepted: &'a ResolvedMaterialStatePoint,
        optical: &'a ResolvedMaterialStatePoint,
        source_identity: ContentHash,
        law: AbsorptionLaw<'a>,
        density_property: &PropertyKey,
    ) -> Result<Self, MaterialAbsorptionError> {
        let density = accepted_density_property(accepted, density_property)
            .ok_or(MaterialAbsorptionError::Density)?;
        let density_kg_m3 = density.value_si();
        let mut hash = DomainHasher::new("org.frankensim.render.material-absorption.v1");
        hash.update(accepted.identity().as_bytes());
        hash.update(source_identity.as_bytes());
        hash.update(&[match law {
            AbsorptionLaw::Visible(BulkOpticalLaw::Conductor(_)) => 0,
            AbsorptionLaw::Visible(BulkOpticalLaw::Dielectric(_)) => 1,
            AbsorptionLaw::Spectral(_) => 2,
        }]);
        hash.update(density.answer().receipt.content_hash().as_bytes());
        hash.update(&density_kg_m3.to_bits().to_le_bytes());
        Ok(Self {
            accepted,
            optical,
            density,
            density_kg_m3,
            law,
            source_identity,
            identity: hash.finalize(),
        })
    }

    /// Original accepted bundle, including its complete state coordinates.
    #[must_use]
    pub const fn accepted_state(&self) -> &'a ResolvedMaterialStatePoint {
        self.accepted
    }

    /// Original optical bundle, without reselecting or discarding source claims.
    #[must_use]
    pub const fn optical_state(&self) -> &'a ResolvedMaterialStatePoint {
        self.optical
    }

    /// Exact original density requirement, answer and usage receipt.
    #[must_use]
    pub const fn density_property(&self) -> &'a ResolvedScalarProperty {
        self.density
    }

    /// Intrinsic material density used to convert inverse metres to m²/kg.
    #[must_use]
    pub const fn material_density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Original optical state and wavelength/key interpretation identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Accepted state, optical source and original density receipt identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Closed original vacuum-wavelength support [nm], including lossless states.
    #[must_use]
    pub fn wavelength_bounds_nm(&self) -> [f64; 2] {
        match self.law {
            AbsorptionLaw::Visible(_) => [LAMBDA_MIN, LAMBDA_MAX],
            AbsorptionLaw::Spectral(source) => source.wavelength_bounds_nm(),
        }
    }

    /// Evaluate the original law at this exact wavelength [m²/kg material].
    ///
    /// Dielectrics use their native intensity absorption directly, without an
    /// alpha-to-k-to-alpha roundtrip. Complex-index sources retain their original
    /// interpolation before the k/lambda conversion. Positive coefficients that
    /// cannot be represented, including underflow to zero, are refused.
    pub fn mass_absorption_at(&self, wavelength_nm: f64) -> Result<f64, MaterialAbsorptionError> {
        let [minimum_nm, maximum_nm] = self.wavelength_bounds_nm();
        if !wavelength_nm.is_finite() || wavelength_nm < minimum_nm || wavelength_nm > maximum_nm {
            return Err(MaterialAbsorptionError::WavelengthOutsideSource {
                wavelength_nm,
                minimum_nm,
                maximum_nm,
            });
        }
        match self.law {
            AbsorptionLaw::Visible(BulkOpticalLaw::Dielectric(glass)) => {
                let alpha = glass
                    .absorption()
                    .extinction_per_m(wavelength_nm)
                    .map_err(MaterialOpticalError::from)?;
                checked_coefficient(alpha / self.density_kg_m3, alpha == 0.0)
            }
            AbsorptionLaw::Visible(BulkOpticalLaw::Conductor(optics)) => {
                let index = optics
                    .complex_ior(wavelength_nm)
                    .map_err(MaterialOpticalError::from)?;
                complex_mass_absorption(index.k(), wavelength_nm, self.density_kg_m3)
            }
            AbsorptionLaw::Spectral(source) => {
                let index = source
                    .at_wavelength(wavelength_nm)
                    .map_err(MaterialAbsorptionError::Source)?;
                complex_mass_absorption(index.k, wavelength_nm, self.density_kg_m3)
            }
        }
    }

    /// Prepare one to sixteen strictly increasing monochromatic point samples.
    ///
    /// This fixed-size preparation allocates no table or field. The supplied
    /// wavelengths remain borrowed; the resulting coefficients imply no law
    /// between those points and no wider spectral support.
    pub fn sample<'w>(
        &self,
        wavelengths_nm: &'w [f64],
    ) -> Result<MaterialAbsorptionSamples<'w>, MaterialAbsorptionError> {
        if !(1..=MAX_OPACITY_WAVELENGTHS).contains(&wavelengths_nm.len())
            || wavelengths_nm
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            || wavelengths_nm.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(MaterialAbsorptionError::InvalidWavelengthSamples);
        }
        let mut coefficients = [0.0; MAX_OPACITY_WAVELENGTHS];
        for (index, &wavelength) in wavelengths_nm.iter().enumerate() {
            coefficients[index] = self.mass_absorption_at(wavelength)?;
        }
        Ok(MaterialAbsorptionSamples {
            wavelengths_nm,
            coefficients,
            accepted_state_identity: self.accepted.identity(),
            source_identity: self.source_identity,
            binding_identity: self.identity,
            source_bounds_nm: self.wavelength_bounds_nm(),
        })
    }
}

/// Exact evaluations of one frozen material law at the supplied wavelengths.
///
/// These are point values, not an interpolated mass-absorption law. They can
/// feed a uniform density observation or an explicitly prepared mixture row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialAbsorptionSamples<'a> {
    wavelengths_nm: &'a [f64],
    coefficients: [f64; MAX_OPACITY_WAVELENGTHS],
    accepted_state_identity: ContentHash,
    source_identity: ContentHash,
    binding_identity: ContentHash,
    source_bounds_nm: [f64; 2],
}

impl<'a> MaterialAbsorptionSamples<'a> {
    /// Original requested wavelengths [nm], without resampling.
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &'a [f64] {
        self.wavelengths_nm
    }

    /// Corresponding intensity mass absorption [m²/kg of the bound material].
    #[must_use]
    pub fn mass_absorption_m2_kg(&self) -> &[f64] {
        &self.coefficients[..self.wavelengths_nm.len()]
    }

    /// Accepted material snapshot used for every prepared value.
    #[must_use]
    pub const fn accepted_state_identity(&self) -> ContentHash {
        self.accepted_state_identity
    }

    /// Original optical-source interpretation used for every prepared value.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// The source, accepted condition and original density receipt.
    #[must_use]
    pub const fn binding_identity(&self) -> ContentHash {
        self.binding_identity
    }

    /// Original closed source support [nm]; points do not widen this interval.
    #[must_use]
    pub const fn source_wavelength_bounds_nm(&self) -> [f64; 2] {
        self.source_bounds_nm
    }

    /// Borrow these samples directly into the existing lean opacity renderer.
    ///
    /// Unit scales describe the caller's original grid. They do not replace
    /// the intrinsic material density already used in kappa; the renderer
    /// validates both scales under its existing SI admission.
    #[must_use]
    pub fn as_density_spectrum(
        &self,
        density_unit_kg_m3: f64,
        length_unit_m: f64,
    ) -> DensityAbsorptionSpectrum<'_> {
        DensityAbsorptionSpectrum {
            wavelengths_nm: self.wavelengths_nm,
            mass_absorption_m2_kg: self.mass_absorption_m2_kg(),
            density_unit_kg_m3,
            length_unit_m,
        }
    }
}

/// Refusal before a source-bound absorption sample can be published.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialAbsorptionError {
    /// Accepted material and optical source disagree, or native optics refused.
    Optical(MaterialOpticalError),
    /// The original source could not evaluate the requested wavelength.
    Source(MaterialStatePointError),
    /// Missing, nonpositive or incorrectly typed scalar bulk density.
    Density,
    /// RGB-derived attenuation cannot supply a physical absorption spectrum.
    SpectralAbsorptionRequired,
    /// Point sampling requires 1--16 finite positive increasing wavelengths.
    InvalidWavelengthSamples,
    /// A requested wavelength falls outside the original closed source support.
    WavelengthOutsideSource {
        wavelength_nm: f64,
        minimum_nm: f64,
        maximum_nm: f64,
    },
    /// Positive mass absorption overflowed or underflowed to zero.
    UnrepresentableCoefficient,
}

impl From<MaterialOpticalError> for MaterialAbsorptionError {
    fn from(error: MaterialOpticalError) -> Self {
        Self::Optical(error)
    }
}

impl core::fmt::Display for MaterialAbsorptionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Optical(error) => error.fmt(f),
            Self::Source(error) => error.fmt(f),
            Self::Density => {
                f.write_str("bulk absorption requires the accepted scalar density in kg/m3")
            }
            Self::SpectralAbsorptionRequired => {
                f.write_str("bulk absorption requires spectral source data, not RGB attenuation")
            }
            Self::InvalidWavelengthSamples => {
                f.write_str("bulk absorption requires 1--16 positive increasing wavelength samples")
            }
            Self::WavelengthOutsideSource {
                wavelength_nm,
                minimum_nm,
                maximum_nm,
            } => write!(
                f,
                "bulk absorption wavelength {wavelength_nm} nm is outside [{minimum_nm}, {maximum_nm}] nm"
            ),
            Self::UnrepresentableCoefficient => {
                f.write_str("positive mass absorption cannot be represented")
            }
        }
    }
}

impl core::error::Error for MaterialAbsorptionError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Optical(error) => Some(error),
            Self::Source(error) => Some(error),
            _ => None,
        }
    }
}

fn checked_coefficient(value: f64, lossless: bool) -> Result<f64, MaterialAbsorptionError> {
    if lossless {
        return Ok(0.0);
    }
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MaterialAbsorptionError::UnrepresentableCoefficient)
    }
}

fn complex_mass_absorption(
    k: f64,
    wavelength_nm: f64,
    density: f64,
) -> Result<f64, MaterialAbsorptionError> {
    if k == 0.0 {
        return Ok(0.0);
    }
    // Avoid forming alpha first: it may overflow although alpha/rho is finite.
    let by_density = k / density;
    let by_wavelength = by_density / wavelength_nm;
    let ordinary = by_wavelength * (4.0 * core::f64::consts::PI * 1e9);
    let value = if by_density.is_finite()
        && by_density >= f64::MIN_POSITIVE
        && by_wavelength.is_finite()
        && by_wavelength >= f64::MIN_POSITIVE
        && ordinary.is_finite()
        && ordinary >= f64::MIN_POSITIVE
    {
        ordinary
    } else {
        det::exp(
            det::ln(k) - det::ln(density) - det::ln(wavelength_nm)
                + det::ln(4.0 * core::f64::consts::PI * 1e9),
        )
    };
    checked_coefficient(value, false)
}

#[cfg(test)]
mod tests;
