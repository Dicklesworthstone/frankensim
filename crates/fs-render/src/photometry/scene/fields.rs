//! Material-driven thermal refraction in the original scene's medium topology.
//!
//! Each binding borrows a complete accepted thermo-optic field. Its original
//! reference glass must equal the selected original scene medium in every
//! phase, absorption and provenance value, not merely at one wavelength. The
//! field supplies the actual local state; nominal glass is never added to it
//! a second time. Geometry, original surface laws and medium indices remain
//! unchanged, including equal-index bodies with distinct physical ownership.

use super::{NativeOpticalMedium, NativeOpticalScene, NativeOpticalSceneError};
use crate::dielectric::{DielectricGlass, DielectricSurface};
use crate::gradient_index::GradientIndexField;
use crate::gradient_index::bundle::SecondOrderIndexField;
use crate::material_state::thermo_optic::{
    MaterialThermoOpticField, TemperatureGradientField, TemperatureSecondOrderField,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::tracer::Material;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;

#[derive(Clone, Copy)]
enum Field<'a> {
    Gradient(&'a (dyn GradientIndexField + Sync)),
    Differentiable(&'a (dyn SecondOrderIndexField + Sync)),
}

/// An explicitly selected original medium and one immutable thermal/material
/// field. Positions are world metres, temperature is kelvin, absorption is per
/// metre, and wavelength is vacuum nanometres, as in the original field owner.
/// The snapshot must remain fixed for the entire observation batch; no time
/// interpolation, heat step, moving geometry or unsupported constitutive law
/// is inferred. The ordinary `render_scene()` remains the reference-glass view;
/// thermal-field observations use the explicitly selected native entry points.
#[derive(Clone, Copy)]
pub struct NativeThermoOpticBinding<'a> {
    medium: NativeOpticalMedium,
    reference: &'a DielectricGlass,
    field: Field<'a>,
    snapshot_identity: ContentHash,
    field_identity: ContentHash,
}

impl<'a> NativeThermoOpticBinding<'a> {
    /// Retain an original first-order field for native ray/camera transport.
    /// Differentiated bundle queries deliberately refuse this capability.
    pub fn from_gradient_field<F: TemperatureGradientField + Sync + ?Sized>(
        medium: NativeOpticalMedium,
        field: &'a MaterialThermoOpticField<'_, F>,
    ) -> Self {
        Self {
            medium,
            reference: field.reference_glass(),
            field: Field::Gradient(field),
            snapshot_identity: field.snapshot_identity(),
            field_identity: identity(field),
        }
    }

    /// Retain the original analytic Hessian for both native rays and bundles.
    /// No finite differences or zero-curvature substitution are introduced.
    pub fn from_differentiable_field<F: TemperatureSecondOrderField + Sync + ?Sized>(
        medium: NativeOpticalMedium,
        field: &'a MaterialThermoOpticField<'_, F>,
    ) -> Self {
        Self {
            medium,
            reference: field.reference_glass(),
            field: Field::Differentiable(field),
            snapshot_identity: field.snapshot_identity(),
            field_identity: identity(field),
        }
    }

    pub const fn medium(&self) -> NativeOpticalMedium {
        self.medium
    }

    /// Original accepted temperature/reconstruction snapshot supplied to the
    /// material field, without promoting that identity to measured authority.
    pub const fn source_snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }

    /// Native medium identity binds the reference optics, complete thermo-optic
    /// constitutive law and original frozen temperature/reconstruction identity.
    /// Choosing first- or second-order capability does not change the physics.
    pub const fn field_identity(&self) -> ContentHash {
        self.field_identity
    }

    fn native_medium(self) -> SceneMedium<'a> {
        match self.field {
            Field::Gradient(field) => SceneMedium::GradientIndex {
                field,
                identity: self.field_identity,
            },
            Field::Differentiable(field) => SceneMedium::Differentiable {
                field,
                identity: self.field_identity,
            },
        }
    }
}

fn identity<F: TemperatureGradientField + ?Sized>(
    field: &MaterialThermoOpticField<'_, F>,
) -> ContentHash {
    let mut hash = DomainHasher::new("org.frankensim.render.native-thermo-optic-medium.v1");
    // Reuse the renderer's canonical encoding of the complete original glass.
    // The fixed surface tag identifies bulk data only; it never scatters a ray.
    let reference = Material::Dielectric {
        glass: *field.reference_glass(),
        surface: DielectricSurface::SMOOTH,
    }
    .content_identity();
    hash.update(reference.as_bytes());
    hash.update(field.law().identity().as_bytes());
    hash.update(field.snapshot_identity().as_bytes());
    hash.finalize()
}

impl NativeOpticalScene {
    /// Observe accepted thermo-optic states over original geometry and media.
    /// Each selected medium's complete reference glass must match the original
    /// scene. Unlisted media keep their original homogeneous law. At most one
    /// field may bind a medium, including the ambient; primitive references
    /// resolve through the same topology as ordinary optical observations.
    ///
    /// The field's physical/domain errors remain native propagation errors.
    /// It must support the original contact-resolution collar on both sides of
    /// a boundary; no extrapolation or nominal-glass fallback is supplied.
    /// Fields with unresolved internal material/contact jumps require their
    /// existing separate-region event treatment. Surface geometry is not
    /// thermally deformed and source temperature is not inferred from a field.
    ///
    /// Preparation is O(media + bindings), with one bounded medium table and
    /// no field copies. Native acceleration is built once for the callback;
    /// each path/bundle query retains its own limits and cancellation points.
    pub fn with_thermo_optic_scene<R>(
        &self,
        cx: &Cx<'_>,
        max_build_work: u64,
        bindings: &[NativeThermoOpticBinding<'_>],
        observe: impl FnOnce(&CurvedPolarizedScene<'_>) -> R,
    ) -> Result<R, NativeOpticalSceneError> {
        cx.checkpoint()?;
        if bindings.len() > self.media.len() {
            return Err(NativeOpticalSceneError::TooManyFieldBindings {
                bindings: bindings.len(),
                media: self.media.len(),
            });
        }
        let mut media = self.native_media(cx)?;
        for binding in bindings {
            cx.checkpoint()?;
            let index = self.medium_index(binding.medium)?;
            if !matches!(media[index], SceneMedium::Homogeneous(_)) {
                return Err(NativeOpticalSceneError::DuplicateFieldBinding {
                    medium: binding.medium,
                });
            }
            if self.media[index] != *binding.reference {
                return Err(NativeOpticalSceneError::ReferenceGlassMismatch {
                    medium: binding.medium,
                });
            }
            media[index] = binding.native_medium();
        }
        self.with_native_media(cx, max_build_work, &media, observe)
    }
}

#[cfg(test)]
mod tests;
