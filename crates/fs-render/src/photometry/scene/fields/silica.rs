//! Original temperature-resolved silica in the native scene's medium topology.
//!
//! The reference is the released absolute CHARMS law at an explicit reference
//! temperature, with the field's original independent absorption. It names the
//! original scene medium; every optical stage evaluates the actual frozen
//! temperature reconstruction instead of adding a correction to that reference.

use super::*;
use crate::material_state::silica_refraction::{
    MaterialSilicaRefractiveField, PreparedMaterialSilicaRefractiveField, corning7980_glass,
};

impl<'a> NativeMaterialFieldBinding<'a> {
    /// Observe an original Corning 7980 temperature field through the native
    /// scene, camera and beam owners. The complete reference glass must equal
    /// the source at `reference_temperature_k`, including its original band
    /// and provenance. Matching only a sampled index is insufficient.
    ///
    /// The reference temperature is not the field temperature. Actual spatial
    /// temperatures, source limits and Beer absorption stay with the field.
    /// Prepared and unprepared queries name the same immutable optical state.
    pub fn from_silica_field<F: TemperatureGradientField + Sync + ?Sized>(
        medium: NativeOpticalMedium,
        reference: &'a DielectricGlass,
        reference_temperature_k: f64,
        field: &'a MaterialSilicaRefractiveField<'_, F>,
    ) -> Result<Self, NativeOpticalSceneError> {
        if reference.absorption_ref() != field.absorption_source() {
            return Err(NativeOpticalSceneError::SilicaAbsorptionMismatch { medium });
        }
        let original = corning7980_glass(reference_temperature_k, *field.absorption_source())
            .map_err(|_| NativeOpticalSceneError::SilicaReferenceRequired { medium })?;
        if *reference != original {
            return Err(NativeOpticalSceneError::SilicaReferenceRequired { medium });
        }
        let mut hash = DomainHasher::new("org.frankensim.render.native-silica-medium.v1");
        let reference_identity = Material::Dielectric {
            glass: original,
            surface: DielectricSurface::SMOOTH,
        }
        .content_identity();
        hash.update(reference_identity.as_bytes());
        hash.update(&reference_temperature_k.to_bits().to_le_bytes());
        hash.update(field.law_identity().as_bytes());
        hash.update(field.snapshot_identity().as_bytes());
        Ok(Self {
            medium,
            reference,
            field: Field::Gradient(field),
            snapshot_identity: field.snapshot_identity(),
            field_identity: hash.finalize(),
        })
    }

    /// Retain the source's analytic thermal curvature for ray bundles and
    /// native point/area connections as well as ordinary camera transport.
    pub fn from_differentiable_silica_field<F: TemperatureSecondOrderField + Sync + ?Sized>(
        medium: NativeOpticalMedium,
        reference: &'a DielectricGlass,
        reference_temperature_k: f64,
        field: &'a MaterialSilicaRefractiveField<'_, F>,
    ) -> Result<Self, NativeOpticalSceneError> {
        Ok(Self {
            field: Field::Differentiable(field),
            ..Self::from_silica_field(medium, reference, reference_temperature_k, field)?
        })
    }

    /// Reuse one original wavelength preparation across a batch of native
    /// camera/beam observations. A different queried wavelength still refuses.
    pub fn from_prepared_silica_field<F: TemperatureGradientField + Sync + ?Sized>(
        medium: NativeOpticalMedium,
        reference: &'a DielectricGlass,
        reference_temperature_k: f64,
        prepared: &'a PreparedMaterialSilicaRefractiveField<'_, '_, F>,
    ) -> Result<Self, NativeOpticalSceneError> {
        Ok(Self {
            field: Field::Gradient(prepared),
            ..Self::from_silica_field(
                medium,
                reference,
                reference_temperature_k,
                prepared.source_field(),
            )?
        })
    }

    /// Reuse the same wavelength preparation with original analytic Hessians.
    pub fn from_prepared_differentiable_silica_field<
        F: TemperatureSecondOrderField + Sync + ?Sized,
    >(
        medium: NativeOpticalMedium,
        reference: &'a DielectricGlass,
        reference_temperature_k: f64,
        prepared: &'a PreparedMaterialSilicaRefractiveField<'_, '_, F>,
    ) -> Result<Self, NativeOpticalSceneError> {
        Ok(Self {
            field: Field::Differentiable(prepared),
            ..Self::from_prepared_silica_field(
                medium,
                reference,
                reference_temperature_k,
                prepared,
            )?
        })
    }
}

#[cfg(test)]
mod tests;
