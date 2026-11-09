//! Original tabulated thermal optics in the native scene's medium topology.
//!
//! The field's complete admitted reference glass names the original scene
//! medium. Actual optical stages retain the same source path and accepted
//! temperature snapshot. Prepared and unprepared forms identify the same
//! physical field, independent of which derivative capability is requested.

use super::*;
use crate::material_state::spectral_refraction::{
    MaterialSpectralRefractiveField, PreparedMaterialSpectralRefractiveField,
};

impl<'a> NativeMaterialFieldBinding<'a> {
    /// Bind an admitted original spectral/temperature source to a native
    /// medium. The complete reference glass still undergoes the native scene's
    /// medium-consistency check before observation; a single matching index
    /// is insufficient. The complete phase/loss/provenance and gap policy
    /// must agree.
    pub fn from_spectral_field<F: TemperatureGradientField + Sync + ?Sized>(
        medium: NativeOpticalMedium,
        field: &'a MaterialSpectralRefractiveField<'_, F>,
    ) -> Self {
        Self {
            medium,
            reference: field.reference_glass(),
            field: Field::Gradient(field),
            snapshot_identity: field.snapshot_identity(),
            field_identity: spectral_identity(field),
        }
    }

    /// Retain original analytic thermal curvature for native bundles and
    /// point/area connections as well as ordinary camera transport.
    pub fn from_differentiable_spectral_field<
        F: TemperatureSecondOrderField + Sync + ?Sized,
    >(
        medium: NativeOpticalMedium,
        field: &'a MaterialSpectralRefractiveField<'_, F>,
    ) -> Self {
        Self {
            field: Field::Differentiable(field),
            ..Self::from_spectral_field(medium, field)
        }
    }

    /// Reuse a single wavelength preparation throughout native observation.
    /// Querying another wavelength retains the prepared owner's refusal.
    pub fn from_prepared_spectral_field<F: TemperatureGradientField + Sync + ?Sized>(
        medium: NativeOpticalMedium,
        prepared: &'a PreparedMaterialSpectralRefractiveField<'_, '_, F>,
    ) -> Self {
        Self {
            field: Field::Gradient(prepared),
            ..Self::from_spectral_field(medium, prepared.source_field())
        }
    }

    /// The same prepared original source with full analytic index Hessians.
    pub fn from_prepared_differentiable_spectral_field<
        F: TemperatureSecondOrderField + Sync + ?Sized,
    >(
        medium: NativeOpticalMedium,
        prepared: &'a PreparedMaterialSpectralRefractiveField<'_, '_, F>,
    ) -> Self {
        Self {
            field: Field::Differentiable(prepared),
            ..Self::from_prepared_spectral_field(medium, prepared)
        }
    }
}

fn spectral_identity<F: TemperatureGradientField + ?Sized>(
    field: &MaterialSpectralRefractiveField<'_, F>,
) -> ContentHash {
    let mut hash = DomainHasher::new("org.frankensim.render.native-tabulated-thermal-medium.v1");
    // Canonical complete bulk reference identity; this fixed surface tag
    // does not scatter a ray or replace any original scene surface law.
    let reference = Material::Dielectric {
        glass: *field.reference_glass(),
        surface: DielectricSurface::SMOOTH,
    }
    .content_identity();
    hash.update(reference.as_bytes());
    hash.update(field.source().identity().as_bytes());
    hash.update(field.snapshot_identity().as_bytes());
    hash.finalize()
}
