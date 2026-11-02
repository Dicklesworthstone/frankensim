//! The original composition reconstruction supplies curved molecular collisions.
//!
//! Nonnegative molecular cross sections make the original ordered sum monotone
//! in every partial density. Evaluating the declared componentwise upper vector
//! at the actual wavelength therefore bounds the local scattering rate without
//! fitting a gray coefficient, changing the mixture phase, or inflating density.

use super::{
    ContentHash, DomainHasher, GasCompositionGradientField, GasPartialDensitySample,
    GasScatteringError, GradientIndexError, MAX_GAS_MIXTURE_COMPONENTS,
    MaterialGasCompositionRefractiveField, Point3, UNWRITTEN,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighField, CurvedRayleighSample,
};
use crate::volumes::transport::VolumeTransportError;
use fs_material::gas_scattering::composition::CompositionRayleighScattering;

/// One immutable scattering view of the exact original refractive field.
///
/// Bounds use that field's original component order and density basis. They
/// must hold everywhere the medium is reachable, including outside any sampled
/// camera frustum. They are caller declarations, checked at every candidate;
/// this binding does not certify a reconstruction's global range. Original
/// provider/domain failures are never interpreted as outside-domain vacuum.
/// In particular, an all-zero envelope cannot validate an unsampled provider.
pub struct MaterialGasCompositionScatteringField<'s, 'a, F: GasCompositionGradientField + ?Sized> {
    refractive: &'s MaterialGasCompositionRefractiveField<'a, F>,
    upper: [f64; MAX_GAS_MIXTURE_COMPONENTS],
    identity: ContentHash,
}

impl<'a, F: GasCompositionGradientField + ?Sized> MaterialGasCompositionRefractiveField<'a, F> {
    /// Share this exact source, provider, density basis and physical snapshot
    /// with molecular delta tracking. At most 128 upper bounds are copied;
    /// construction and local queries allocate no memory.
    pub fn bind_scattering_field<'s>(
        &'s self,
        component_upper_bounds: &[f64],
    ) -> Result<MaterialGasCompositionScatteringField<'s, 'a, F>, GasCompositionScatteringError>
    {
        MaterialGasCompositionScatteringField::try_new(self, component_upper_bounds)
    }
}

impl<'s, 'a, F: GasCompositionGradientField + ?Sized>
    MaterialGasCompositionScatteringField<'s, 'a, F>
{
    /// Admit finite nonnegative bounds in the original field's declared units.
    /// The complete upper vector must itself satisfy the material's dilute
    /// law at both source endpoints. Nonnegative Cauchy coefficients put the
    /// largest refractivity at the shorter wavelength; runtime queries still
    /// use the original material admission at their exact wavelength.
    pub fn try_new(
        refractive: &'s MaterialGasCompositionRefractiveField<'a, F>,
        component_upper_bounds: &[f64],
    ) -> Result<Self, GasCompositionScatteringError> {
        if component_upper_bounds.len() != refractive.model.components().len() {
            return Err(GasCompositionScatteringError::Components);
        }
        let mut upper = [0.0; MAX_GAS_MIXTURE_COMPONENTS];
        let mut hash = DomainHasher::new("org.frankensim.render.gas-composition-scattering.v1");
        hash.update(refractive.content_identity().as_bytes());
        hash.update(&(component_upper_bounds.len() as u64).to_le_bytes());
        for (component, &value) in component_upper_bounds.iter().enumerate() {
            if !value.is_finite() || value < 0.0 {
                return Err(GasCompositionScatteringError::InvalidBound { component });
            }
            upper[component] = if value == 0.0 { 0.0 } else { value };
            hash.update(&upper[component].to_bits().to_le_bytes());
        }
        let result = Self {
            refractive,
            upper,
            identity: hash.finalize(),
        };
        for wavelength in refractive.model.wavelength_range_nm() {
            result.scattering_majorant_optics(wavelength)?;
        }
        Ok(result)
    }

    /// The actual borrowed GRIN owner, not a reconstructed optical copy.
    pub const fn refractive_field(&self) -> &'s MaterialGasCompositionRefractiveField<'a, F> {
        self.refractive
    }

    /// Original snapshot identity plus the exact component bounds. Native
    /// scene attachment uses `source_field_identity`, which remains the GRIN
    /// owner's identity; this separate identity also records tracking policy.
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    pub fn component_upper_bounds(&self) -> &[f64] {
        &self.upper[..self.refractive.model.components().len()]
    }

    /// Diagnostic majorant evaluation with the material's exact refusal intact.
    pub fn scattering_majorant_optics(
        &self,
        wavelength_nm: f64,
    ) -> Result<CompositionRayleighScattering, GasCompositionScatteringError> {
        self.validate_source()?;
        self.refractive
            .model
            .at_wavelength(wavelength_nm, self.refractive.basis)
            .map_err(GasCompositionScatteringError::Gas)?
            .scattering(self.component_upper_bounds())
            .map_err(GasCompositionScatteringError::Gas)
    }

    /// Evaluate the same local partial densities that bend the ray. The
    /// molecular owner supplies both the wavelength-specific rate and the
    /// scattering-weighted depolarization; no gray phase or gas law is fitted.
    pub fn sample_scattering_optics(
        &self,
        position_m: Point3,
        wavelength_nm: f64,
    ) -> Result<CompositionRayleighScattering, GasCompositionScatteringError> {
        self.validate_source()?;
        if [position_m.x, position_m.y, position_m.z]
            .into_iter()
            .any(|value| !value.is_finite())
        {
            return Err(GasCompositionScatteringError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let spectral = self
            .refractive
            .model
            .at_wavelength(wavelength_nm, self.refractive.basis)
            .map_err(GasCompositionScatteringError::Gas)?;
        let count = self.refractive.model.components().len();
        let mut samples = [UNWRITTEN; MAX_GAS_MIXTURE_COMPONENTS];
        self.refractive
            .field
            .sample_partial_densities(position_m, &mut samples[..count])
            .map_err(GasCompositionScatteringError::Field)?;
        let mut densities = [0.0; MAX_GAS_MIXTURE_COMPONENTS];
        for (component, sample) in samples[..count].iter().enumerate() {
            let GasPartialDensitySample {
                density,
                gradient_per_m: g,
            } = *sample;
            if !density.is_finite()
                || density < 0.0
                || [g.x, g.y, g.z].into_iter().any(|value| !value.is_finite())
            {
                return Err(GasCompositionScatteringError::InvalidSample { component });
            }
            if density > self.upper[component] {
                return Err(GasCompositionScatteringError::BoundExceeded { component });
            }
            densities[component] = density;
        }
        spectral
            .scattering(&densities[..count])
            .map_err(GasCompositionScatteringError::Gas)
    }

    fn validate_source(&self) -> Result<(), GasCompositionScatteringError> {
        if self.refractive.field.component_count() != self.refractive.model.components().len() {
            return Err(GasCompositionScatteringError::Components);
        }
        if self.refractive.field.model_identity() != self.refractive.model.content_identity()
            || self.refractive.field.density_basis() != self.refractive.basis
        {
            return Err(GasCompositionScatteringError::SourceMismatch);
        }
        Ok(())
    }
}

impl<F: GasCompositionGradientField + Sync + ?Sized> CurvedRayleighField
    for MaterialGasCompositionScatteringField<'_, '_, F>
{
    fn source_field_identity(&self) -> ContentHash {
        self.refractive.content_identity()
    }

    fn scattering_majorant_per_m(&self, wavelength_nm: f64) -> Result<f64, VolumeTransportError> {
        self.scattering_majorant_optics(wavelength_nm)
            .map(|sample| sample.scattering_per_m)
            .map_err(GasCompositionScatteringError::into_volume_error)
    }

    fn sample_scattering(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError> {
        self.sample_scattering_optics(position, wavelength_nm)
            .map(|sample| CurvedRayleighSample {
                scattering_per_m: sample.scattering_per_m,
                depolarization_ratio: sample.depolarization_ratio,
            })
            .map_err(GasCompositionScatteringError::into_volume_error)
    }
}

/// Original material/provider failures and explicit reconstruction-bound errors.
#[derive(Clone, Debug, PartialEq)]
pub enum GasCompositionScatteringError {
    /// Row count differs from the original molecular model.
    Components,
    /// Provider source identity or density units changed after GRIN binding.
    SourceMismatch,
    InvalidBound {
        component: usize,
    },
    /// A row was unwritten, its density negative/nonfinite, or its gradient nonfinite.
    InvalidSample {
        component: usize,
    },
    BoundExceeded {
        component: usize,
    },
    Field(GradientIndexError),
    Gas(GasScatteringError),
}

impl GasCompositionScatteringError {
    fn into_volume_error(self) -> VolumeTransportError {
        match self {
            Self::SourceMismatch | Self::Gas(GasScatteringError::StateMismatch) => {
                VolumeTransportError::SourceIdentity
            }
            Self::Components | Self::Gas(GasScatteringError::Components) => {
                VolumeTransportError::Mixture
            }
            Self::InvalidBound { component }
            | Self::InvalidSample { component }
            | Self::BoundExceeded { component } => {
                VolumeTransportError::Density { index: component }
            }
            Self::Field(GradientIndexError::Cancelled) => VolumeTransportError::Cancelled,
            Self::Field(GradientIndexError::NumericalRange)
            | Self::Gas(GasScatteringError::Unrepresentable) => {
                VolumeTransportError::Unrepresentable
            }
            Self::Gas(GasScatteringError::Wavelength) => VolumeTransportError::Wavelength,
            Self::Field(_) => VolumeTransportError::Grid,
            Self::Gas(_) => VolumeTransportError::ScatteringSpectrum,
        }
    }
}

impl core::fmt::Display for GasCompositionScatteringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material gas composition scattering refused: {self:?}")
    }
}
impl core::error::Error for GasCompositionScatteringError {}

#[cfg(test)]
mod tests;
