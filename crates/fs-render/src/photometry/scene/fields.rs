//! Source-backed material refraction in the original scene's medium topology.
//!
//! Each binding borrows a complete accepted thermo-optic field. Its original
//! reference glass must equal the selected original scene medium in every
//! phase, absorption and provenance value, not merely at one wavelength. The
//! field supplies the actual local state; nominal glass is never added to it
//! a second time. Geometry, original surface laws and medium indices remain
//! unchanged, including equal-index bodies with distinct physical ownership.

use super::{
    NativeOpticalBinding, NativeOpticalMedium, NativeOpticalScene, NativeOpticalSceneError,
    NativeReferenceMedium, reserved,
};
use crate::dielectric::{DielectricGlass, DielectricSurface};
use crate::gradient_index::GradientIndexField;
use crate::gradient_index::bundle::SecondOrderIndexField;
use crate::instances::InstanceScene;
use crate::material_state::gas_scattering::{
    GasCompositionGradientField, GasCompositionSecondOrderField, GasDensityGradientField,
    GasDensitySecondOrderField, MaterialGasCompositionRefractiveField, MaterialGasRefractiveField,
    PreparedMaterialGasCompositionRefractiveField,
};
use crate::material_state::thermo_optic::{
    MaterialThermoOpticField, TemperatureGradientField, TemperatureSecondOrderField,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::tracer::{Material, Scene};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;

#[derive(Clone, Copy)]
enum Field<'a> {
    Gradient(&'a (dyn GradientIndexField + Sync)),
    Differentiable(&'a (dyn SecondOrderIndexField + Sync)),
}

impl<'a> Field<'a> {
    fn native_medium(self, identity: ContentHash) -> SceneMedium<'a> {
        match self {
            Self::Gradient(field) => SceneMedium::GradientIndex { field, identity },
            Self::Differentiable(field) => SceneMedium::Differentiable { field, identity },
        }
    }
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
        self.field.native_medium(self.field_identity)
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
        let mut media = self.native_media(cx)?;
        self.bind_thermal_media(cx, &mut media, bindings)?;
        self.with_native_media(cx, max_build_work, &media, observe)
    }

    fn bind_thermal_media<'a>(
        &self,
        cx: &Cx<'_>,
        media: &mut [SceneMedium<'a>],
        bindings: &[NativeThermoOpticBinding<'a>],
    ) -> Result<(), NativeOpticalSceneError> {
        cx.checkpoint()?;
        if bindings.len() > self.media.len() {
            return Err(NativeOpticalSceneError::TooManyFieldBindings {
                bindings: bindings.len(),
                media: self.media.len(),
            });
        }
        for binding in bindings {
            cx.checkpoint()?;
            let index = self.medium_index(binding.medium)?;
            if !matches!(media[index], SceneMedium::Homogeneous(_)) {
                return Err(NativeOpticalSceneError::DuplicateFieldBinding {
                    medium: binding.medium,
                });
            }
            if self.media[index].glass() != Some(binding.reference) {
                return Err(NativeOpticalSceneError::ReferenceGlassMismatch {
                    medium: binding.medium,
                });
            }
            media[index] = binding.native_medium();
        }
        Ok(())
    }
}

/// One accepted molecular material and frozen density/composition snapshot.
/// Only the original material gas owners can construct this ambient; there is
/// no arbitrary GRIN injection, nominal glass fit or inferred species/EOS.
/// The source model, density basis and supported vacuum-wavelength band remain
/// authoritative. Gas true absorption is exactly zero under these owners;
/// molecular scattering has its separate source-coupled native attachment.
#[derive(Clone, Copy)]
pub struct NativeGasAmbient<'a> {
    material_identity: ContentHash,
    field_identity: ContentHash,
    wavelength_range_nm: [f64; 2],
    field: Field<'a>,
}

impl<'a> NativeGasAmbient<'a> {
    /// Original fixed-composition density field, including
    /// `IsobaricGasField::bind_refraction`. Its declared constant-polarizability
    /// visible model has equal phase/group indices; no dispersion is invented.
    pub fn from_density_field<F: GasDensityGradientField + Sync + ?Sized>(
        field: &'a MaterialGasRefractiveField<'_, F>,
    ) -> Self {
        Self {
            material_identity: field.state().model_identity(),
            field_identity: field.content_identity(),
            // This is the original RayleighGasState refraction domain.
            wavelength_range_nm: [380.0, 780.0],
            field: Field::Gradient(field),
        }
    }

    /// Retain both nonlinear Lorentz--Lorenz density-jet terms for native
    /// bundle focusing; a locally zero gradient never substitutes for a Hessian.
    pub fn from_differentiable_density_field<F: GasDensitySecondOrderField + Sync + ?Sized>(
        field: &'a MaterialGasRefractiveField<'_, F>,
    ) -> Self {
        Self {
            field: Field::Differentiable(field),
            ..Self::from_density_field(field)
        }
    }

    /// Original partial-density field with its exact component order, units,
    /// optical sources and, where supplied, molecular spectral dispersion.
    pub fn from_composition_field<F: GasCompositionGradientField + Sync + ?Sized>(
        field: &'a MaterialGasCompositionRefractiveField<'_, F>,
    ) -> Self {
        Self {
            material_identity: field.model().content_identity(),
            field_identity: field.content_identity(),
            wavelength_range_nm: field.model().wavelength_range_nm(),
            field: Field::Gradient(field),
        }
    }

    /// Original component Hessians and full molecular chain rule, including
    /// focusing when total mass density is constant but composition varies.
    pub fn from_differentiable_composition_field<
        F: GasCompositionSecondOrderField + Sync + ?Sized,
    >(
        field: &'a MaterialGasCompositionRefractiveField<'_, F>,
    ) -> Self {
        Self {
            field: Field::Differentiable(field),
            ..Self::from_composition_field(field)
        }
    }

    /// Reuse the original stack-owned wavelength preparation across every ray
    /// stage and connection trial. Its exact physical snapshot must match the
    /// supplied source field; another queried wavelength remains an error.
    /// Source-band metadata is retained, not narrowed to a fictitious material.
    pub fn from_prepared_composition_field<F: GasCompositionGradientField + Sync + ?Sized>(
        source: &MaterialGasCompositionRefractiveField<'_, F>,
        prepared: &'a PreparedMaterialGasCompositionRefractiveField<'_, '_, F>,
    ) -> Result<Self, NativeOpticalSceneError> {
        if source.content_identity() != prepared.content_identity() {
            return Err(NativeOpticalSceneError::GasSnapshotMismatch);
        }
        Ok(Self {
            material_identity: source.model().content_identity(),
            field_identity: source.content_identity(),
            wavelength_range_nm: source.model().wavelength_range_nm(),
            field: Field::Gradient(prepared),
        })
    }

    /// Keep that same spectral preparation and the source's analytic Hessian.
    pub fn from_prepared_differentiable_composition_field<
        F: GasCompositionSecondOrderField + Sync + ?Sized,
    >(
        source: &MaterialGasCompositionRefractiveField<'_, F>,
        prepared: &'a PreparedMaterialGasCompositionRefractiveField<'_, '_, F>,
    ) -> Result<Self, NativeOpticalSceneError> {
        Ok(Self {
            field: Field::Differentiable(prepared),
            ..Self::from_prepared_composition_field(source, prepared)?
        })
    }

    pub const fn material_identity(&self) -> ContentHash {
        self.material_identity
    }

    /// Exactly the original field identity, also used by the corresponding
    /// `MaterialGasCompositionScatteringField` attachment. No rehash disconnects
    /// molecular collisions from their actual refractive snapshot.
    pub const fn field_identity(&self) -> ContentHash {
        self.field_identity
    }

    /// Original material source band. A prepared binding admits only its one
    /// selected wavelength even though it retains this complete source band.
    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        self.wavelength_range_nm
    }
}

/// Original static geometry/materials with a source-backed molecular ambient.
///
/// This type deliberately has no ordinary `render_scene` or `into_scene` view:
/// that renderer cannot carry the gas field, and substituting vacuum or fitted
/// nominal glass would change the physics. Every native observation explicitly
/// supplies an accepted snapshot of the material admitted at construction.
/// Geometry is shared exactly as in `NativeOpticalScene`; no gas surface or
/// replacement volume geometry is created. Positions and density derivatives
/// retain their original SI units and the snapshot stays fixed for the batch.
pub struct NativeGasOpticalScene {
    pub(super) scene: NativeOpticalScene,
}

impl NativeGasOpticalScene {
    pub fn try_from_scene(
        cx: &Cx<'_>,
        scene: Scene,
        ambient: NativeGasAmbient<'_>,
        bindings: &[NativeOpticalBinding],
        max_objects: usize,
    ) -> Result<Self, NativeOpticalSceneError> {
        Ok(Self {
            scene: NativeOpticalScene::try_from_scene_with_ambient(
                cx,
                scene,
                NativeReferenceMedium::Molecular {
                    material_identity: ambient.material_identity,
                    wavelength_range_nm: ambient.wavelength_range_nm,
                },
                bindings,
                max_objects,
            )?,
        })
    }

    /// The unchanged original geometry, including original physical blockers.
    pub const fn geometry(&self) -> &InstanceScene {
        self.scene.geometry()
    }

    pub fn medium_index(
        &self,
        medium: NativeOpticalMedium,
    ) -> Result<usize, NativeOpticalSceneError> {
        self.scene.medium_index(medium)
    }

    /// Original interior reference glass; the gas ambient returns `None` and
    /// never fabricates a homogeneous optical law. Indices use `medium_index`.
    pub fn medium(&self, index: usize) -> Option<&DielectricGlass> {
        self.scene.medium(index)
    }

    /// Observe gas refraction and optional original interior thermo-optic
    /// fields through the same native geometry and local boundary optics.
    /// The ambient argument is mandatory; an ambient thermal binding would
    /// duplicate its ownership and refuses. Original molecular material/band
    /// mismatches refuse even at equal index or zero density. A newly accepted
    /// frozen snapshot of the same material may be supplied for another batch.
    ///
    /// Source spatial support, dilute limits, wavelength domains and contact
    /// collars remain authoritative. The gas field is not extrapolated into
    /// vacuum, and internal discontinuities need their original event treatment.
    /// With a composition field, its existing source-bound scattering adapter
    /// can attach to `CurvedRayleighScene` using the unchanged field identity;
    /// this coherent refraction view alone does not add scattered radiance.
    pub fn with_curved_scene<R>(
        &self,
        cx: &Cx<'_>,
        max_build_work: u64,
        ambient: NativeGasAmbient<'_>,
        thermal: &[NativeThermoOpticBinding<'_>],
        observe: impl FnOnce(&CurvedPolarizedScene<'_>) -> R,
    ) -> Result<R, NativeOpticalSceneError> {
        cx.checkpoint()?;
        match self.scene.media.first() {
            Some(NativeReferenceMedium::Molecular {
                material_identity,
                wavelength_range_nm,
            }) if *material_identity == ambient.material_identity
                && *wavelength_range_nm == ambient.wavelength_range_nm => {}
            _ => return Err(NativeOpticalSceneError::GasMaterialMismatch),
        }
        let mut media = reserved(self.scene.media.len())?;
        media.push(ambient.field.native_medium(ambient.field_identity));
        for medium in self.scene.media.iter().skip(1) {
            cx.checkpoint()?;
            media.push(SceneMedium::Homogeneous(
                medium
                    .glass()
                    .ok_or(NativeOpticalSceneError::MissingGasAmbient)?,
            ));
        }
        self.scene.bind_thermal_media(cx, &mut media, thermal)?;
        self.scene
            .with_native_media(cx, max_build_work, &media, observe)
    }
}

#[cfg(test)]
mod gas_tests;
#[cfg(test)]
mod tests;
