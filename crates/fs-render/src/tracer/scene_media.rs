//! Physical cell fields owned by declared closed dielectric interiors.
//!
//! Grids remain in world metres and borrow accepted density/temperature arrays.
//! An explicit binding replaces the glass owner's homogeneous extinction,
//! including when the bound grid is vacuum. The dielectric still supplies its
//! interface law and phase index. Properly nested boundaries select the top
//! owner; leaving an inner body restores the enclosing field.

use super::*;
use crate::volumes::transport::timeline::PhysicalMedium;

mod execution;
pub use execution::{
    SceneMediaFilm, render_scene_media, render_scene_media_range,
    render_scene_media_range_with_execution, render_scene_media_with_execution,
};

/// A physical field replacing the homogeneous extinction of one glass body.
///
/// The primitive must declare an outward-oriented, closed dielectric boundary.
/// The ordinary path walk checks boundary nesting and closure; this binding
/// does not certify mesh watertightness. The grid's finite support is clipped by
/// the next surface, and is vacuum outside its own world-space bounds.
#[derive(Clone, Copy, Debug)]
pub struct InteriorMedium<'view, 'field> {
    /// Index in the admitted scene's primitive array.
    pub boundary_primitive: usize,
    /// Borrowed SI absorption/emission or scattering field.
    pub medium: PhysicalMedium<'view, 'field>,
}

/// Bounds on admission scans and once-prepared light metadata.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SceneMediaBudget {
    /// Maximum number of supplied interior bindings.
    pub max_interior_bindings: usize,
    /// Maximum scene primitive count inspected during admission.
    pub max_primitives: usize,
    /// Maximum thermal source triangles plus two triangles per rectangle.
    /// Borrowed geometry and field storage is not duplicated by this budget.
    pub max_light_triangles: usize,
}

/// Invalid physical ownership is refused before rendering starts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SceneMediaError {
    /// Admission metadata exceeds the caller's limit or integer domain.
    Budget,
    /// Bindings must have strictly increasing primitive IDs, without duplicates.
    UnsortedBindings {
        /// Index of the first non-increasing binding.
        index: usize,
    },
    /// A binding names no scene primitive.
    MissingBoundary {
        /// Invalid primitive index.
        boundary_primitive: usize,
    },
    /// A physical interior requires a dielectric boundary owner.
    NotDielectric {
        /// Primitive whose material is not dielectric.
        boundary_primitive: usize,
    },
    /// This static scene entry does not move world-space field support.
    AnimatedGeometry {
        /// Primitive requiring an exposure-aware scene/field binding.
        primitive: usize,
    },
}

impl core::fmt::Display for SceneMediaError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "physical scene media refused: {self:?}")
    }
}

impl core::error::Error for SceneMediaError {}

impl From<SceneMediaError> for TracerError {
    fn from(error: SceneMediaError) -> Self {
        Self::SceneMedia(error)
    }
}

/// Immutable scene, explicit medium ownership and prepared light proposals.
///
/// Absorption/scattering coefficients use vacuum wavelength, density is kg/m³,
/// distances are metres, and temperatures are kelvin. The host dielectric's
/// spatially constant phase index governs refraction and the equilibrium
/// interior source factor `n(lambda)^2`; this does not model gradient-index
/// media. Camera rays start outside every declared dielectric, as in the
/// ordinary tracer. Arbitrary overlapping bodies and moving fields are not
/// admitted by this entry.
///
/// Numeric continuity binds the scene and all supplied medium states once.
/// Direct chart backends retain the existing process-local replacement guard;
/// it is not a persisted content-authority digest for chart internals.
pub struct SceneMedia<'scene, 'field> {
    pub(super) scene: &'scene Scene,
    pub(super) lighting: SceneLighting<'scene>,
    ambient: Option<PhysicalMedium<'scene, 'field>>,
    interiors: &'scene [InteriorMedium<'scene, 'field>],
    identity: ContentHash,
    has_thermal_emission: bool,
}

impl<'scene, 'field> SceneMedia<'scene, 'field> {
    /// Admit sorted bindings and prepare ordinary source proposals once.
    ///
    /// Existing grid constructors own field validation and per-segment budgets.
    /// This constructor scans only scene/binding metadata and geometry identity;
    /// it never copies or rehashes density/temperature buffers. A zero-density
    /// interior remains an explicit replacement owner, while an ambient vacuum
    /// can use the ordinary surface fast path.
    pub fn try_new(
        cx: &Cx<'_>,
        scene: &'scene Scene,
        ambient: Option<PhysicalMedium<'scene, 'field>>,
        interiors: &'scene [InteriorMedium<'scene, 'field>],
        budget: SceneMediaBudget,
    ) -> Result<Self, TracerError> {
        cx.checkpoint()?;
        if interiors.len() > budget.max_interior_bindings
            || scene.primitives.len() > budget.max_primitives
        {
            return Err(SceneMediaError::Budget.into());
        }
        let mut triangles = scene
            .lights
            .len()
            .checked_mul(2)
            .filter(|&value| value <= budget.max_light_triangles)
            .ok_or(SceneMediaError::Budget)?;
        for (index, primitive) in scene.primitives.iter().enumerate() {
            cx.checkpoint()?;
            if matches!(&primitive.shape, Shape::AnimatedInstance(_)) {
                return Err(SceneMediaError::AnimatedGeometry { primitive: index }.into());
            }
            if primitive
                .material
                .thermal_surface()
                .is_some_and(|surface| surface.is_emitting())
            {
                let count = match &primitive.shape {
                    Shape::Mesh(mesh) => mesh.triangles.len(),
                    Shape::Instance(instance) => match instance.geometry() {
                        SharedGeometry::Mesh(mesh) => mesh.triangles.len(),
                        SharedGeometry::Chart(_) => 0,
                    },
                    Shape::Chart(_) | Shape::AnimatedInstance(_) => 0,
                };
                triangles = triangles
                    .checked_add(count)
                    .filter(|&value| value <= budget.max_light_triangles)
                    .ok_or(SceneMediaError::Budget)?;
            }
        }

        let mut hash = DomainHasher::new("org.frankensim.render.scene-media.v1");
        hash.update(scene_continuity_fingerprint(scene, cx)?.as_bytes());
        match ambient {
            Some(medium) => {
                hash.update(&[1]);
                update_medium_identity(&mut hash, medium);
            }
            None => hash.update(&[0]),
        }
        hash.update(&(interiors.len() as u64).to_le_bytes());
        let mut has_thermal_emission = ambient.is_some_and(PhysicalMedium::has_thermal_emission);
        for (index, binding) in interiors.iter().enumerate() {
            cx.checkpoint()?;
            if index > 0 && interiors[index - 1].boundary_primitive >= binding.boundary_primitive {
                return Err(SceneMediaError::UnsortedBindings { index }.into());
            }
            let primitive = scene.primitives.get(binding.boundary_primitive).ok_or(
                SceneMediaError::MissingBoundary {
                    boundary_primitive: binding.boundary_primitive,
                },
            )?;
            if !matches!(primitive.material, Material::Dielectric { .. }) {
                return Err(SceneMediaError::NotDielectric {
                    boundary_primitive: binding.boundary_primitive,
                }
                .into());
            }
            hash.update(&(binding.boundary_primitive as u64).to_le_bytes());
            update_medium_identity(&mut hash, binding.medium);
            has_thermal_emission |= binding.medium.has_thermal_emission();
        }
        // A dark accepted physical state is meaningful. This flag permits
        // an empty source set; it does not create emission or change fields.
        let lighting = validate_scene_with_medium_emission(scene, cx, None, true)?;
        cx.checkpoint()?;
        Ok(Self {
            scene,
            lighting,
            ambient,
            interiors,
            identity: hash.finalize(),
            has_thermal_emission,
        })
    }

    /// Immutable geometry, optical interfaces and light-source owners.
    #[must_use]
    pub const fn scene(&self) -> &'scene Scene {
        self.scene
    }

    /// Optional index-one field outside declared dielectric interiors.
    #[must_use]
    pub const fn ambient(&self) -> Option<PhysicalMedium<'scene, 'field>> {
        self.ambient
    }

    /// Sorted explicit interior owners, borrowing the caller's binding array.
    #[must_use]
    pub const fn interiors(&self) -> &'scene [InteriorMedium<'scene, 'field>] {
        self.interiors
    }

    /// Numeric scene/field continuity identity used by progressive films.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Whether any supplied physical field has an LTE source.
    #[must_use]
    pub const fn has_thermal_emission(&self) -> bool {
        self.has_thermal_emission
    }

    pub(super) fn ambient_path(&self) -> Option<PathMedium<'scene, 'field>> {
        self.ambient
            .filter(|medium| !medium.is_vacuum())
            .map(path_medium)
    }

    pub(super) fn interior_path(
        &self,
        boundary_primitive: usize,
    ) -> Option<PathMedium<'scene, 'field>> {
        self.interiors
            .binary_search_by_key(&boundary_primitive, |binding| binding.boundary_primitive)
            .ok()
            .map(|index| path_medium(self.interiors[index].medium))
    }
}

fn path_medium<'view, 'field>(medium: PhysicalMedium<'view, 'field>) -> PathMedium<'view, 'field> {
    match medium {
        PhysicalMedium::Absorbing(medium) => PathMedium::Absorbing(medium),
        PhysicalMedium::Scattering(medium) if medium.has_scattering() => {
            PathMedium::Scattering(medium)
        }
        PhysicalMedium::Scattering(medium) => PathMedium::Absorbing(medium.absorption_grid()),
    }
}

fn update_medium_identity(hash: &mut DomainHasher, medium: PhysicalMedium<'_, '_>) {
    hash.update(&[match medium {
        PhysicalMedium::Absorbing(_) => 0,
        PhysicalMedium::Scattering(_) => 1,
    }]);
    hash.update(medium.content_identity().as_bytes());
}

#[cfg(test)]
mod tests;
