//! Calibrated reciprocal camera observations of an original tracer scene.
//!
//! Original primitive emission is two-sided and independent of the scattering
//! material. Rectangular lights are proposal metadata naming that same source;
//! they never register a second emitter. Native area proposals can be prepared
//! from the returned camera owner with their existing geometric and image budgets.

#[cfg(feature = "material-state")]
use super::{NativeGasAmbient, NativeGasOpticalScene, NativeThermoOpticBinding};
use super::{NativeOpticalScene, NativeOpticalSceneError, reserved};
use crate::camera::{Aperture, CameraError, PhysicalCamera};
use crate::lighting::{AdmittedLighting, LightingError};
use crate::polarization::dielectric::scene::curved::CurvedPolarizedScene;
use crate::polarization::dielectric::scene::curved::camera::{
    CameraRadianceError, CurvedCameraScene, EmissionSides, EnvironmentRadiance, RadianceSpectrum,
    SurfaceEmitter,
};
use crate::tracer::rectangle_light_matches_primitive;
use fs_exec::{Cancelled, Cx};

/// Explicit conversion of original scene radiance units to local SI spectral
/// radiance. Surface and environment units need not share a calibration. Both
/// scales apply to the original spectrum times its original authored strength;
/// they do not normalize total flux, alter color or infer a measured spectrum.
/// Original thermal materials already emit SI radiance from their accepted
/// temperature and emissivity; neither calibration scale applies to them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NativeRadianceCalibration {
    surface_w_m2_sr_nm_per_unit: f64,
    environment_w_m2_sr_nm_per_unit: f64,
}

impl NativeRadianceCalibration {
    pub fn try_new(
        surface_w_m2_sr_nm_per_unit: f64,
        environment_w_m2_sr_nm_per_unit: f64,
    ) -> Result<Self, NativeCameraError> {
        if [surface_w_m2_sr_nm_per_unit, environment_w_m2_sr_nm_per_unit]
            .into_iter()
            .any(|scale| !scale.is_finite() || scale < 0.0)
        {
            return Err(NativeCameraError::InvalidCalibration);
        }
        Ok(Self {
            surface_w_m2_sr_nm_per_unit,
            environment_w_m2_sr_nm_per_unit,
        })
    }

    pub const fn surface_w_m2_sr_nm_per_unit(self) -> f64 {
        self.surface_w_m2_sr_nm_per_unit
    }

    pub const fn environment_w_m2_sr_nm_per_unit(self) -> f64 {
        self.environment_w_m2_sr_nm_per_unit
    }
}

/// Original-camera and source admission failures retain their physical owner.
#[derive(Clone, Debug, PartialEq)]
pub enum NativeCameraError {
    InvalidCalibration,
    InvalidEmission { primitive: usize },
    LightBudgetExceeded { lights: usize, objects: usize },
    LightPrimitiveMismatch { light: usize, primitive: usize },
    OpticalScene(NativeOpticalSceneError),
    Camera(CameraError),
    Radiance(CameraRadianceError),
    Lighting(LightingError),
}

impl From<NativeOpticalSceneError> for NativeCameraError {
    fn from(error: NativeOpticalSceneError) -> Self {
        Self::OpticalScene(error)
    }
}
impl From<CameraError> for NativeCameraError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}
impl From<CameraRadianceError> for NativeCameraError {
    fn from(error: CameraRadianceError) -> Self {
        Self::Radiance(error)
    }
}
impl From<LightingError> for NativeCameraError {
    fn from(error: LightingError) -> Self {
        Self::Lighting(error)
    }
}
impl From<Cancelled> for NativeCameraError {
    fn from(_: Cancelled) -> Self {
        Self::OpticalScene(NativeOpticalSceneError::Cancelled)
    }
}
impl core::fmt::Display for NativeCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original scene camera: {self:?}")
    }
}
impl core::error::Error for NativeCameraError {}

struct CameraSources<'a> {
    camera: PhysicalCamera,
    emitters: Vec<SurfaceEmitter>,
    environment: Option<EnvironmentRadiance<'a>>,
}

#[cfg(feature = "material-state")]
impl NativeGasOpticalScene {
    /// Original camera/source radiance through the explicitly selected accepted
    /// molecular ambient, with optional original interior thermal-index fields.
    /// Source calibration, geometry and reflection are shared with the ordinary
    /// native camera; the gas field supplies the actual wavelength-dependent
    /// phase/group law and spatial derivatives. It never supplies an absorption
    /// or thermal-emission law that the molecular material does not own.
    ///
    /// Every query retains the source spectrum's own support and emitting-host
    /// restrictions. `CurvedRayleighCameraScene` can compose this exact camera
    /// with the original composition scattering field inside the callback;
    /// the ambient's field identity already matches that native attachment.
    /// Scattering then remains elastic transport, separate from true absorption.
    #[allow(clippy::too_many_arguments)]
    pub fn with_curved_camera<R>(
        &self,
        cx: &Cx<'_>,
        max_build_work: u64,
        calibration: NativeRadianceCalibration,
        ambient: NativeGasAmbient<'_>,
        thermal: &[NativeThermoOpticBinding<'_>],
        observe: impl FnOnce(&CurvedCameraScene<'_, '_>, &PhysicalCamera) -> R,
    ) -> Result<R, NativeCameraError> {
        let sources = self.scene.camera_sources(cx, calibration)?;
        self.with_curved_scene(cx, max_build_work, ambient, thermal, |scene| {
            self.scene.observe_camera(cx, scene, &sources, observe)
        })?
    }
}

impl NativeOpticalScene {
    /// Prepare original surface/environment emission and the original pinhole
    /// camera over the shared native optical scene for one observation batch.
    /// Geometry coordinates must already be metres, as in the optical owner.
    /// The explicit calibration converts the renderer's authored radiance units
    /// to W/(m² sr nm), independently for surfaces and the environment.
    ///
    /// Every original emissive primitive becomes one two-sided source, with its
    /// exact lifted spectrum and geometry. Rectangles remain proposals naming
    /// these same primitives; their identity, emission and original geometry are
    /// checked using the renderer's own admission rules. They add no emission.
    /// `CurvedAreaEmitter::from_camera` can prepare that original surface for
    /// native area-light MIS inside the callback, with explicit image seeds.
    /// The ordinary environment uses its exact native RGB spectral reconstruction.
    /// Original gray/spectral diffuse thermal materials instead retain their
    /// own Planck/Kirchhoff source in SI units, independent of this calibration.
    /// Their source identity, two-sided emission and exact reciprocal reflection
    /// also serve native area proposals. They retain the original index-one
    /// emitting-medium restriction, even when the camera lies inside glass.
    ///
    /// The existing reciprocal camera retains n_camera²/n_source², actual Beer
    /// paths, emission before scattering, native material laws and unresolved
    /// path errors. `CurvedCameraRaster` can sample the supplied physical camera;
    /// its four wavelengths each trace their original refracted path. The caller
    /// declares the initial medium through `medium_index`; no containment is
    /// inferred from the camera position. Detector surfaces remain unsupported
    /// by the reciprocal camera owner. No automatic image completeness or NEE
    /// strategy is inferred from legacy rectangle metadata.
    ///
    /// Source storage is bounded by the original object count; native geometry
    /// preparation uses `max_build_work`. Each query retains its own existing
    /// budgets. The callback may return its own `Result`, as `with_curved_scene`
    /// does; admission and observation errors remain separate.
    pub fn with_curved_camera<R>(
        &self,
        cx: &Cx<'_>,
        max_build_work: u64,
        calibration: NativeRadianceCalibration,
        observe: impl FnOnce(&CurvedCameraScene<'_, '_>, &PhysicalCamera) -> R,
    ) -> Result<R, NativeCameraError> {
        let sources = self.camera_sources(cx, calibration)?;
        self.with_curved_scene(cx, max_build_work, |scene| {
            self.observe_camera(cx, scene, &sources, observe)
        })?
    }

    /// Observe the original camera and sources through accepted thermal-index
    /// fields bound to original media. Every wavelength retains its actual
    /// local phase/group law, gradient and absorption. Differentiable bindings
    /// also serve native point/area-light bundle connections in this callback.
    ///
    /// Reference-glass, snapshot, capability and contact-domain rules are those
    /// of `with_thermo_optic_scene`. Source preparation and radiance calibration
    /// are shared with `with_curved_camera`: field temperature never invents
    /// surface emission or changes the original accepted thermal emitters.
    #[cfg(feature = "material-state")]
    pub fn with_thermo_optic_camera<R>(
        &self,
        cx: &Cx<'_>,
        max_build_work: u64,
        calibration: NativeRadianceCalibration,
        bindings: &[NativeThermoOpticBinding<'_>],
        observe: impl FnOnce(&CurvedCameraScene<'_, '_>, &PhysicalCamera) -> R,
    ) -> Result<R, NativeCameraError> {
        let sources = self.camera_sources(cx, calibration)?;
        self.with_thermo_optic_scene(cx, max_build_work, bindings, |scene| {
            self.observe_camera(cx, scene, &sources, observe)
        })?
    }

    fn camera_sources(
        &self,
        cx: &Cx<'_>,
        calibration: NativeRadianceCalibration,
    ) -> Result<CameraSources<'_>, NativeCameraError> {
        cx.checkpoint()?;
        let original = &self.render.camera;
        let camera = PhysicalCamera::try_legacy_compatible(
            original.eye,
            original.forward,
            original.up,
            original.half_tan,
            // Focus is unused by the exact pinhole branch; this positive
            // storage value introduces no inferred optical focus or aperture.
            1.0,
            Aperture::try_circular(0.0)?,
        )?;
        self.admit_camera_rectangles(cx)?;
        let mut count = 0;
        for primitive in &self.render.primitives {
            cx.checkpoint()?;
            count += usize::from(primitive.emission.is_some());
        }
        let mut emitters = reserved(count)?;
        for (primitive, (original, binding)) in self
            .render
            .primitives
            .iter()
            .zip(&self.bindings)
            .enumerate()
        {
            cx.checkpoint()?;
            let Some((spectrum, strength)) = original.emission else {
                continue;
            };
            if !strength.is_finite()
                || strength < 0.0
                || spectrum
                    .c
                    .iter()
                    .any(|coefficient| !coefficient.is_finite())
            {
                return Err(NativeCameraError::InvalidEmission { primitive });
            }
            let scale = strength * calibration.surface_w_m2_sr_nm_per_unit;
            if !scale.is_finite()
                || (scale == 0.0 && strength > 0.0 && calibration.surface_w_m2_sr_nm_per_unit > 0.0)
            {
                return Err(NativeCameraError::InvalidEmission { primitive });
            }
            emitters.push(SurfaceEmitter::new(
                binding.object_id,
                RadianceSpectrum::try_lifted(spectrum, scale)?,
                EmissionSides::Both,
            ));
        }
        let environment = self
            .render
            .environment
            .as_ref()
            .map(|map| {
                EnvironmentRadiance::try_new(map, calibration.environment_w_m2_sr_nm_per_unit)
            })
            .transpose()?;
        cx.checkpoint()?;
        Ok(CameraSources {
            camera,
            emitters,
            environment,
        })
    }

    fn observe_camera<R>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        sources: &CameraSources<'_>,
        observe: impl FnOnce(&CurvedCameraScene<'_, '_>, &PhysicalCamera) -> R,
    ) -> Result<R, NativeCameraError> {
        let radiance =
            CurvedCameraScene::try_new(cx, scene, &sources.emitters, sources.environment)?
                .with_diffuse_thermal_emitters(cx, &self.thermal_emitters)?;
        cx.checkpoint()?;
        Ok(observe(&radiance, &sources.camera))
    }

    fn admit_camera_rectangles(&self, cx: &Cx<'_>) -> Result<(), NativeCameraError> {
        let lights = &self.render.lights;
        if lights.len() > self.render.primitives.len() {
            return Err(NativeCameraError::LightBudgetExceeded {
                lights: lights.len(),
                objects: self.render.primitives.len(),
            });
        }
        if !lights.is_empty() {
            // The existing bounded source list remains only proposal metadata.
            // Admission checks original shape/emission and duplicate ownership;
            // no second copy enters the native camera's source accumulator.
            AdmittedLighting::try_new_cancellable(cx, lights, None)?;
            for (light, proposal) in lights.iter().enumerate() {
                cx.checkpoint()?;
                let matches = self
                    .render
                    .primitives
                    .get(proposal.prim)
                    .is_some_and(|primitive| {
                        rectangle_light_matches_primitive(proposal, primitive)
                    });
                if !matches {
                    return Err(NativeCameraError::LightPrimitiveMismatch {
                        light,
                        primitive: proposal.prim,
                    });
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
