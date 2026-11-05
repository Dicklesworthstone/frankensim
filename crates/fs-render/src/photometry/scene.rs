//! Original tracer scenes admitted to calibrated curved optical transport.
//!
//! [`NativeOpticalScene`] moves the original static meshes, charts and placed
//! instances into shared geometry. The ordinary rendering view retains every
//! material, emitter, light, environment and camera; the curved view uses those
//! same geometry allocations and exact supported optical laws. It does not
//! replace a lens by planes, fit a spectral table, or infer medium containment.
//!
//! Callers declare each boundary's surrounding medium and each detector. The
//! normal of a dielectric boundary must point from its interior to that declared
//! exterior. Parent references are checked for cycles, but geometric nesting is
//! the caller's responsibility. Native traversal checks the arriving medium at
//! every interaction, including equal-index media with different absorption.
//!
//! Use [`NativeOpticalScene::with_curved_scene`] to prepare the native boundary
//! accelerator once for a batch of point connections, image discovery or finite
//! sensor observations. Their existing budgets, branch policies, polarization,
//! Jacobians, visibility and convergence semantics remain in force. A selected
//! connection is not an all-image measurement. Source observations measure the
//! explicitly supplied source. Use [`NativeOpticalScene::with_curved_camera`]
//! with an explicit [`NativeRadianceCalibration`] to observe the original scene
//! emitters and environment through its own pinhole camera. The reciprocal
//! camera and area-light strategies retain the same geometry and source laws.
//! With `material-state`, `with_thermo_optic_scene` and
//! `with_thermo_optic_camera` bind accepted thermal fields to these same medium
//! references. Their native rays and bundles then see the actual local index,
//! its analytic derivatives, group delay and original absorption.
//! `NativeGasOpticalScene` admits an original molecular ambient separately;
//! its observations require that source's accepted gas field and never expose
//! an ordinary rendering view that would silently drop the ambient physics.

use core::fmt;
use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::{Cancelled, Cx};
use fs_geom::Vec3;

use crate::conductor::ConductorOptics;
use crate::dielectric::DielectricGlass;
use crate::instances::{
    GeometryInstance, InstanceError, InstanceScene, RigidTransform, SharedGeometry,
};
use crate::polarization::dielectric::scene::curved::camera::DiffuseThermalEmitter;
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedSceneError, SceneMedium,
};
use crate::polarization::dielectric::scene::{
    Interaction, MAX_POLARIZED_SCENE_ITEMS, SurfaceBinding,
};
use crate::spectral::SpectralReflectance;
use crate::tracer::{Material, Primitive, Scene, Shape};

mod camera;
pub use camera::{NativeCameraError, NativeRadianceCalibration};
#[cfg(feature = "material-state")]
mod fields;
#[cfg(feature = "material-state")]
pub use fields::{
    NativeGasAmbient, NativeGasOpticalScene, NativeMaterialFieldBinding, NativeThermoOpticBinding,
};

/// Stable identity for geometry that an ordinary tracer scene has not placed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NativeGeometryBinding {
    /// Retain the existing instance's object ID, immutable identity and placement.
    OriginalInstance,
    /// Move a raw mesh or chart into an identity-placed shared instance. The
    /// caller supplies its real immutable geometry identity; it is never guessed.
    Unplaced {
        object_id: u64,
        geometry_identity: ContentHash,
    },
}

/// Physical medium references use original primitive indices, not material
/// equality. Distinct dielectric bodies retain distinct native medium indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NativeOpticalMedium {
    Ambient,
    /// Interior of a dielectric primitive whose role is `Material`. A detector
    /// intercepts rays and therefore cannot also define a traversable interior.
    InteriorOf(usize),
}

/// The role of an original surface in calibrated transport.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NativeOpticalRole {
    /// Use the original material in this explicitly declared surrounding medium.
    Material { outside: NativeOpticalMedium },
    /// Observe the original geometry instead of scattering from its material.
    /// The world reference axis is projected transverse to each arriving ray.
    Detector { reference_axis: Vec3 },
}

/// Exactly one entry, in original primitive order, is required for every shape.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NativeOpticalBinding {
    pub geometry: NativeGeometryBinding,
    pub role: NativeOpticalRole,
}

/// Explicit admission failures; unsupported physics never becomes a dark sample.
#[derive(Clone, Debug, PartialEq)]
pub enum NativeOpticalSceneError {
    InvalidObjectLimit,
    ObjectBudgetExceeded {
        objects: usize,
        max_objects: usize,
    },
    BindingCount {
        primitives: usize,
        bindings: usize,
    },
    GeometryBindingMismatch {
        primitive: usize,
    },
    UnsupportedGeometry {
        primitive: usize,
        reason: &'static str,
    },
    UnsupportedMaterial {
        primitive: usize,
        reason: &'static str,
    },
    InvalidMaterial {
        primitive: usize,
    },
    ConflictingThermalEmission {
        primitive: usize,
    },
    InvalidDetectorAxis {
        primitive: usize,
    },
    InvalidMediumReference {
        medium: NativeOpticalMedium,
    },
    CyclicMedia {
        primitive: usize,
    },
    MediumBudgetExceeded,
    TooManyFieldBindings {
        bindings: usize,
        media: usize,
    },
    DuplicateFieldBinding {
        medium: NativeOpticalMedium,
    },
    ReferenceGlassMismatch {
        medium: NativeOpticalMedium,
    },
    /// Water fields require the original IAPWS water law, not an equal index fit.
    WaterReferenceRequired {
        medium: NativeOpticalMedium,
    },
    /// The field must retain the original independent water absorption law.
    WaterAbsorptionMismatch {
        medium: NativeOpticalMedium,
    },
    /// A source-backed gas ambient requires its original material field.
    MissingGasAmbient,
    /// A different molecular material cannot replace the admitted ambient.
    GasMaterialMismatch,
    /// Spectral preparation must belong to the same frozen physical field.
    GasSnapshotMismatch,
    Geometry(InstanceError),
    Curved(CurvedSceneError),
    AllocationRefused,
    Cancelled,
}

impl From<Cancelled> for NativeOpticalSceneError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}
impl From<InstanceError> for NativeOpticalSceneError {
    fn from(error: InstanceError) -> Self {
        match error {
            InstanceError::Cancelled => Self::Cancelled,
            error => Self::Geometry(error),
        }
    }
}
impl From<CurvedSceneError> for NativeOpticalSceneError {
    fn from(error: CurvedSceneError) -> Self {
        Self::Curved(error)
    }
}
impl fmt::Display for NativeOpticalSceneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "native optical scene: {self:?}")
    }
}
impl std::error::Error for NativeOpticalSceneError {}

#[derive(Clone, Copy)]
enum NativeReferenceMedium {
    Glass(DielectricGlass),
    #[cfg(feature = "material-state")]
    Molecular {
        material_identity: ContentHash,
        wavelength_range_nm: [f64; 2],
    },
}

impl NativeReferenceMedium {
    fn glass(&self) -> Option<&DielectricGlass> {
        match self {
            Self::Glass(glass) => Some(glass),
            #[cfg(feature = "material-state")]
            Self::Molecular { .. } => None,
        }
    }

    fn wavelength_bounds_nm(&self) -> (f64, f64) {
        match self {
            Self::Glass(glass) => glass.ior().wavelength_bounds_nm(),
            #[cfg(feature = "material-state")]
            Self::Molecular {
                wavelength_range_nm: [lo, hi],
                ..
            } => (*lo, *hi),
        }
    }
}

/// One immutable scene shared by ordinary rendering and calibrated optics.
/// Repeated observations borrow this owner; no mesh, chart or placement is rebuilt.
pub struct NativeOpticalScene {
    render: Scene,
    geometry: InstanceScene,
    media: Vec<NativeReferenceMedium>,
    interior_media: Vec<Option<usize>>,
    conductors: Vec<ConductorOptics>,
    reflectances: Vec<SpectralReflectance>,
    thermal_emitters: Vec<DiffuseThermalEmitter>,
    // Original primitive order; native admission sorts its bounded working copy.
    bindings: Vec<SurfaceBinding>,
}

impl NativeOpticalScene {
    /// Own an ordinary static scene without replacing its geometry or optics.
    /// `max_objects` is a caller allocation bound in `1..=65_536`. Native
    /// boundary preparation is separately charged by `with_curved_scene`.
    ///
    /// Dielectric, smooth/rough/brushed conductor, spectral diffuse and coated
    /// diffuse laws are retained exactly. Lifted Lambertian spectra use their
    /// original polynomial throughout the native visible band without a table
    /// fit. Original gray/spectral diffuse thermal materials retain their
    /// accepted temperature and exact reciprocal reflection; the camera adds
    /// their original SI Planck emission without authored-radiance calibration.
    /// Legacy Schlick GGX refuses until the curved owner supports its law.
    /// Animated geometry requires an explicit caller-selected static snapshot.
    pub fn try_from_scene(
        cx: &Cx<'_>,
        scene: Scene,
        ambient: DielectricGlass,
        bindings: &[NativeOpticalBinding],
        max_objects: usize,
    ) -> Result<Self, NativeOpticalSceneError> {
        Self::try_from_scene_with_ambient(
            cx,
            scene,
            NativeReferenceMedium::Glass(ambient),
            bindings,
            max_objects,
        )
    }

    fn try_from_scene_with_ambient(
        cx: &Cx<'_>,
        scene: Scene,
        ambient: NativeReferenceMedium,
        bindings: &[NativeOpticalBinding],
        max_objects: usize,
    ) -> Result<Self, NativeOpticalSceneError> {
        cx.checkpoint()?;
        if max_objects == 0 || max_objects > MAX_POLARIZED_SCENE_ITEMS {
            return Err(NativeOpticalSceneError::InvalidObjectLimit);
        }
        let count = scene.primitives.len();
        if count > max_objects {
            return Err(NativeOpticalSceneError::ObjectBudgetExceeded {
                objects: count,
                max_objects,
            });
        }
        if count != bindings.len() {
            return Err(NativeOpticalSceneError::BindingCount {
                primitives: count,
                bindings: bindings.len(),
            });
        }
        let (mut dielectric_count, mut conductor_count, mut reflectance_count, mut thermal_count) =
            (0, 0, 0, 0);
        for (index, (primitive, binding)) in scene.primitives.iter().zip(bindings).enumerate() {
            cx.checkpoint()?;
            if primitive.material.thermal_surface().is_some() && primitive.emission.is_some() {
                return Err(NativeOpticalSceneError::ConflictingThermalEmission {
                    primitive: index,
                });
            }
            if matches!(binding.role, NativeOpticalRole::Material { .. }) {
                match primitive.material {
                    Material::Dielectric { .. } => dielectric_count += 1,
                    Material::Conductor { .. } | Material::BrushedConductor { .. } => {
                        conductor_count += 1
                    }
                    Material::SpectralDiffuse { .. } | Material::CoatedDiffuse { .. } => {
                        reflectance_count += 1
                    }
                    Material::GrayDiffuseThermal { .. }
                    | Material::SpectralDiffuseThermal { .. } => {
                        reflectance_count += 1;
                        thermal_count += 1;
                    }
                    _ => {}
                }
            }
        }
        if dielectric_count >= MAX_POLARIZED_SCENE_ITEMS {
            return Err(NativeOpticalSceneError::MediumBudgetExceeded);
        }
        let mut media = reserved(dielectric_count + 1)?;
        media.push(ambient);
        let mut interior_media = reserved(count)?;
        for (primitive, binding) in scene.primitives.iter().zip(bindings) {
            cx.checkpoint()?;
            let interior =
                if let (Material::Dielectric { glass, .. }, NativeOpticalRole::Material { .. }) =
                    (primitive.material, binding.role)
                {
                    if media.len() == MAX_POLARIZED_SCENE_ITEMS {
                        return Err(NativeOpticalSceneError::MediumBudgetExceeded);
                    }
                    let index = media.len();
                    media.push(NativeReferenceMedium::Glass(glass));
                    Some(index)
                } else {
                    None
                };
            interior_media.push(interior);
        }
        validate_topology(cx, bindings, &interior_media)?;

        let mut instances = reserved(count)?;
        let mut primitives = reserved(count)?;
        let mut surfaces = reserved(count)?;
        let mut conductors = reserved(conductor_count)?;
        let mut reflectances = reserved(reflectance_count)?;
        let mut thermal_emitters = reserved(thermal_count)?;
        for (index, (primitive, binding)) in scene.primitives.into_iter().zip(bindings).enumerate()
        {
            cx.checkpoint()?;
            let instance = match (primitive.shape, binding.geometry) {
                (Shape::Instance(instance), NativeGeometryBinding::OriginalInstance) => instance,
                (
                    Shape::Mesh(mesh),
                    NativeGeometryBinding::Unplaced {
                        object_id,
                        geometry_identity,
                    },
                ) => GeometryInstance::try_new(
                    object_id,
                    geometry_identity,
                    SharedGeometry::mesh(mesh),
                    RigidTransform::identity(),
                )?,
                (
                    Shape::Chart(chart),
                    NativeGeometryBinding::Unplaced {
                        object_id,
                        geometry_identity,
                    },
                ) => GeometryInstance::try_new(
                    object_id,
                    geometry_identity,
                    SharedGeometry::Chart(Arc::from(chart)),
                    RigidTransform::identity(),
                )?,
                (Shape::AnimatedInstance(_), _) => {
                    return Err(NativeOpticalSceneError::UnsupportedGeometry {
                        primitive: index,
                        reason: "animated geometry requires an explicit static snapshot",
                    });
                }
                _ => {
                    return Err(NativeOpticalSceneError::GeometryBindingMismatch {
                        primitive: index,
                    });
                }
            };
            let interaction = match binding.role {
                NativeOpticalRole::Detector { reference_axis } => {
                    if ![reference_axis.x, reference_axis.y, reference_axis.z]
                        .iter()
                        .all(|v| v.is_finite())
                        || reference_axis
                            .x
                            .abs()
                            .max(reference_axis.y.abs())
                            .max(reference_axis.z.abs())
                            == 0.0
                    {
                        return Err(NativeOpticalSceneError::InvalidDetectorAxis {
                            primitive: index,
                        });
                    }
                    Interaction::Detector { reference_axis }
                }
                NativeOpticalRole::Material { outside } => {
                    let medium = resolve_medium(outside, &interior_media)?;
                    match primitive.material {
                        Material::Dielectric { surface, .. } => Interaction::Dielectric {
                            outside: medium,
                            // The first pass allocated every material dielectric's interior.
                            inside: interior_media[index].ok_or(
                                NativeOpticalSceneError::InvalidMaterial { primitive: index },
                            )?,
                            surface,
                        },
                        Material::Conductor { optics, surface } => {
                            let index = conductors.len();
                            conductors.push(optics);
                            if surface.is_delta() {
                                Interaction::SmoothConductor {
                                    optics: index,
                                    medium,
                                }
                            } else {
                                Interaction::RoughConductor {
                                    optics: index,
                                    medium,
                                    surface,
                                }
                            }
                        }
                        Material::BrushedConductor { optics, detail } => {
                            let index = conductors.len();
                            conductors.push(optics);
                            Interaction::BrushedConductor {
                                optics: index,
                                medium,
                                detail,
                            }
                        }
                        Material::SpectralDiffuse { reflectance } => {
                            let index = reflectances.len();
                            reflectances.push(reflectance);
                            Interaction::DepolarizingDiffuse {
                                medium,
                                reflectance: index,
                            }
                        }
                        Material::CoatedDiffuse {
                            coating,
                            reflectance,
                        } => {
                            let index = reflectances.len();
                            reflectances.push(reflectance);
                            Interaction::CoatedDiffuse {
                                medium,
                                reflectance: index,
                                coating,
                            }
                        }
                        Material::Lambertian { reflectance } => {
                            if reflectance
                                .c
                                .iter()
                                .any(|coefficient| !coefficient.is_finite())
                            {
                                return Err(NativeOpticalSceneError::InvalidMaterial {
                                    primitive: index,
                                });
                            }
                            Interaction::DepolarizingLiftedDiffuse {
                                medium,
                                reflectance,
                            }
                        }
                        Material::Ggx { .. } => {
                            return Err(NativeOpticalSceneError::UnsupportedMaterial {
                                primitive: index,
                                reason: "legacy Schlick GGX has no exact curved material owner",
                            });
                        }
                        Material::GrayDiffuseThermal { surface } => {
                            let emitter = DiffuseThermalEmitter::try_gray_in_band(
                                instance.object_id(),
                                surface,
                                media[medium].wavelength_bounds_nm(),
                            )
                            .map_err(|_| {
                                NativeOpticalSceneError::InvalidMaterial { primitive: index }
                            })?;
                            bind_thermal_surface(
                                index,
                                medium,
                                emitter,
                                &mut reflectances,
                                &mut thermal_emitters,
                            )?
                        }
                        Material::SpectralDiffuseThermal { surface } => bind_thermal_surface(
                            index,
                            medium,
                            DiffuseThermalEmitter::new(instance.object_id(), surface),
                            &mut reflectances,
                            &mut thermal_emitters,
                        )?,
                    }
                }
            };
            surfaces.push(SurfaceBinding {
                object_id: instance.object_id(),
                interaction,
            });
            primitives.push(Primitive {
                shape: Shape::Instance(instance.clone()),
                material: primitive.material,
                emission: primitive.emission,
            });
            instances.push(instance);
        }
        cx.checkpoint()?;
        let geometry = InstanceScene::try_new(instances)?;
        cx.checkpoint()?;
        Ok(Self {
            render: Scene {
                primitives,
                lights: scene.lights,
                environment: scene.environment,
                camera: scene.camera,
            },
            geometry,
            media,
            interior_media,
            conductors,
            reflectances,
            thermal_emitters,
            bindings: surfaces,
        })
    }

    /// Original renderer-facing state, with static geometry now shared through
    /// identity-preserving instances. Mutation is excluded while optics borrows it.
    #[must_use]
    pub const fn render_scene(&self) -> &Scene {
        &self.render
    }

    /// Recover the ordinary rendering scene and release the optical view.
    #[must_use]
    pub fn into_scene(self) -> Scene {
        self.render
    }

    /// The original shapes and placements, sorted by their stable object IDs.
    #[must_use]
    pub const fn geometry(&self) -> &InstanceScene {
        &self.geometry
    }

    /// Native interaction for an original primitive index, before ID sorting.
    #[must_use]
    pub fn surface_binding(&self, primitive: usize) -> Option<&SurfaceBinding> {
        self.bindings.get(primitive)
    }

    /// Resolve a source's initial medium using the same explicit topology as
    /// the boundaries. A source inside glass must supply its interior reference.
    pub fn medium_index(
        &self,
        medium: NativeOpticalMedium,
    ) -> Result<usize, NativeOpticalSceneError> {
        resolve_medium(medium, &self.interior_media)
    }

    /// Original physical medium at a resolved native index, including absorption.
    #[must_use]
    pub fn medium(&self, index: usize) -> Option<&DielectricGlass> {
        self.media.get(index).and_then(NativeReferenceMedium::glass)
    }

    /// Prepare the existing curved transport owner once for a batch of native
    /// source connections, image discovery, camera paths or finite detectors.
    /// `max_build_work` bounds the native geometry accelerator. Query, solve and
    /// path budgets remain explicit in each observation's existing request.
    ///
    /// The scoped view avoids a self-referential scene and never copies geometry.
    /// The callback may return its normal `Result`; use two `?` operators to
    /// propagate admission and observation errors separately.
    pub fn with_curved_scene<R>(
        &self,
        cx: &Cx<'_>,
        max_build_work: u64,
        observe: impl FnOnce(&CurvedPolarizedScene<'_>) -> R,
    ) -> Result<R, NativeOpticalSceneError> {
        let media = self.native_media(cx)?;
        self.with_native_media(cx, max_build_work, &media, observe)
    }

    fn native_media(&self, cx: &Cx<'_>) -> Result<Vec<SceneMedium<'_>>, NativeOpticalSceneError> {
        cx.checkpoint()?;
        let mut media = reserved(self.media.len())?;
        for medium in &self.media {
            cx.checkpoint()?;
            media.push(SceneMedium::Homogeneous(
                medium
                    .glass()
                    .ok_or(NativeOpticalSceneError::MissingGasAmbient)?,
            ));
        }
        Ok(media)
    }

    fn with_native_media<R>(
        &self,
        cx: &Cx<'_>,
        max_build_work: u64,
        media: &[SceneMedium<'_>],
        observe: impl FnOnce(&CurvedPolarizedScene<'_>) -> R,
    ) -> Result<R, NativeOpticalSceneError> {
        cx.checkpoint()?;
        let mut bindings = reserved(self.bindings.len())?;
        for binding in &self.bindings {
            cx.checkpoint()?;
            bindings.push(*binding);
        }
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &self.geometry,
            media,
            &self.conductors,
            &[],
            &self.reflectances,
            bindings,
            max_build_work,
        )?;
        cx.checkpoint()?;
        Ok(observe(&scene))
    }
}

// Material reflection and camera emission retain the same accepted Kirchhoff
// source. The native material owns 1-epsilon; emission is attached separately
// only when observing the original scene camera, never fitted into a lamp.
fn bind_thermal_surface(
    primitive: usize,
    medium: usize,
    emitter: DiffuseThermalEmitter,
    reflectances: &mut Vec<SpectralReflectance>,
    emitters: &mut Vec<DiffuseThermalEmitter>,
) -> Result<Interaction, NativeOpticalSceneError> {
    let reflectance = reflectances.len();
    reflectances.push(
        emitter
            .reflectance_spectrum()
            .map_err(|_| NativeOpticalSceneError::InvalidMaterial { primitive })?,
    );
    emitters.push(emitter);
    Ok(Interaction::DepolarizingDiffuse {
        medium,
        reflectance,
    })
}

fn reserved<T>(count: usize) -> Result<Vec<T>, NativeOpticalSceneError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| NativeOpticalSceneError::AllocationRefused)?;
    Ok(values)
}

fn resolve_medium(
    reference: NativeOpticalMedium,
    interiors: &[Option<usize>],
) -> Result<usize, NativeOpticalSceneError> {
    match reference {
        NativeOpticalMedium::Ambient => Ok(0),
        NativeOpticalMedium::InteriorOf(primitive) => interiors
            .get(primitive)
            .copied()
            .flatten()
            .ok_or(NativeOpticalSceneError::InvalidMediumReference { medium: reference }),
    }
}

// Each dielectric has one declared parent. Three-color traversal visits each
// node at most twice, so nested media admission remains O(objects), not O(depth²).
fn validate_topology(
    cx: &Cx<'_>,
    bindings: &[NativeOpticalBinding],
    interiors: &[Option<usize>],
) -> Result<(), NativeOpticalSceneError> {
    for binding in bindings {
        cx.checkpoint()?;
        if let NativeOpticalRole::Material { outside } = binding.role {
            resolve_medium(outside, interiors)?;
        }
    }
    let mut colors = reserved(bindings.len())?;
    colors.resize(bindings.len(), 0_u8);
    let parent = |index: usize| match bindings[index].role {
        NativeOpticalRole::Material {
            outside: NativeOpticalMedium::InteriorOf(parent),
        } => Some(parent),
        _ => None,
    };
    for start in 0..bindings.len() {
        cx.checkpoint()?;
        if interiors[start].is_none() || colors[start] != 0 {
            continue;
        }
        let mut current = Some(start);
        while let Some(index) = current {
            cx.checkpoint()?;
            match colors[index] {
                1 => return Err(NativeOpticalSceneError::CyclicMedia { primitive: index }),
                2 => break,
                _ => {
                    colors[index] = 1;
                    current = parent(index);
                }
            }
        }
        current = Some(start);
        while let Some(index) = current {
            cx.checkpoint()?;
            if colors[index] == 2 {
                break;
            }
            colors[index] = 2;
            current = parent(index);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
