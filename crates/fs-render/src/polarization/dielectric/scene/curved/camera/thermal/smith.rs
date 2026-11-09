//! Native thermal camera with matching complete Smith conductor reflection.
//!
//! Emission and continuation use the same admitted optical/roughness law and
//! separate deterministic random walks. The thermal estimate is contracted at
//! the original geometric arrival, before the independent continuation acts.
//! Existing calibrated surface, environment and volume/LTE sources keep their
//! original camera accumulation. This explicit owner cannot silently install a
//! multiple-scattering emissivity on a single-scattering camera or NEE path.

use super::super::microsurface::SmithCameraScene;
use super::super::{
    CameraAnalyzer, CameraRadianceError, CameraRadianceSample, CurvedCameraScene, CurvedPathResult,
    CurvedSceneError, CurvedTraceSettings, Interaction, phase_index,
};
use crate::camera::SequentialCameraBeam;
use crate::charts::Ray;
use crate::conductor::ConductorSurface;
use crate::gradient_index::event::BoundarySide;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{SceneError, smith_draw};
use crate::polarization::dielectric::transport::PolarizedBeam;
use crate::polarization::dielectric::transport::connection::normalized_camera;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use crate::thermal_surface::smith::{SmithCoatingThermalEmission, SmithThermalError};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// One geometric side of the original rough, smooth or untilted brushed metal.
/// The kernel retains the original substrate/film optics and separate physical
/// temperatures. Roughness and body-locked orientation come from the scene.
#[derive(Clone, Debug)]
pub struct SmithThermalEmitter {
    object_id: u64,
    side: BoundarySide,
    source: SmithCoatingThermalEmission,
}

impl SmithThermalEmitter {
    pub fn try_new(
        object_id: u64,
        side: BoundarySide,
        source: SmithCoatingThermalEmission,
    ) -> Result<Self, CameraRadianceError> {
        if side == BoundarySide::OnSurface {
            return Err(CameraRadianceError::InvalidEmitter { object_id });
        }
        Ok(Self {
            object_id,
            side,
            source,
        })
    }

    pub const fn object_id(&self) -> u64 {
        self.object_id
    }
    pub const fn side(&self) -> BoundarySide {
        self.side
    }
    pub const fn source(&self) -> &SmithCoatingThermalEmission {
        &self.source
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum SmithThermalCameraError {
    Camera(CameraRadianceError),
    Thermal {
        object_id: u64,
        error: SmithThermalError,
    },
}
impl From<CameraRadianceError> for SmithThermalCameraError {
    fn from(error: CameraRadianceError) -> Self {
        Self::Camera(error)
    }
}
impl From<fs_exec::Cancelled> for SmithThermalCameraError {
    fn from(error: fs_exec::Cancelled) -> Self {
        Self::Camera(error.into())
    }
}
impl core::fmt::Display for SmithThermalCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Smith thermal camera: {self:?}")
    }
}
impl core::error::Error for SmithThermalCameraError {}

/// An explicit full-Smith camera and a bounded table of thermal face bindings.
/// Preparation allocates at most two indices per original surface; sampling
/// uses stack state. Sources and native scene are borrowed without duplication.
pub struct SmithThermalCameraScene<'a, 's, 'm> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    sources: &'a [SmithThermalEmitter],
    slots: Vec<[Option<usize>; 2]>,
    microsurface: SmithWalkSettings,
}

impl<'a, 's, 'm> SmithThermalCameraScene<'a, 's, 'm> {
    /// Match the exact original optics, refuse duplicate or competing surface
    /// emitters, and reject authored normal tilt. Geometric anisotropy, rigid
    /// body rotation, native roughness marks and the radial-pole fallback stay
    /// with the original finish. Thermal sources currently emit into local
    /// phase index one; wavelength/medium admission occurs at each encounter.
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedCameraScene<'s, 'm>,
        sources: &'a [SmithThermalEmitter],
        microsurface: SmithWalkSettings,
    ) -> Result<Self, SmithThermalCameraError> {
        cx.checkpoint()?;
        let mut slots = Vec::new();
        if !sources.is_empty() {
            let maximum = camera.scene.bindings.len().saturating_mul(2);
            if sources.len() > maximum {
                return Err(CameraRadianceError::InvalidEmitter {
                    object_id: sources[maximum].object_id,
                }
                .into());
            }
            slots
                .try_reserve_exact(camera.scene.bindings.len())
                .map_err(|_| CameraRadianceError::AllocationRefused)?;
            slots.resize(camera.scene.bindings.len(), [None, None]);
            for (source_index, source) in sources.iter().enumerate() {
                cx.checkpoint()?;
                let invalid = CameraRadianceError::InvalidEmitter {
                    object_id: source.object_id,
                };
                let boundary = camera
                    .scene
                    .bindings
                    .binary_search_by_key(&source.object_id, |binding| binding.object_id)
                    .map_err(|_| invalid.clone())?;
                let optics = match camera.scene.bindings[boundary].interaction {
                    Interaction::SmoothConductor { optics, .. }
                    | Interaction::RoughConductor { optics, .. } => optics,
                    Interaction::BrushedConductor { optics, detail, .. }
                        if detail.normal_perturbation_rad() == 0.0 =>
                    {
                        optics
                    }
                    _ => return Err(invalid.into()),
                };
                let side = usize::from(source.side == BoundarySide::Negative);
                if camera.scene.conductors.get(optics) != Some(source.source.optics())
                    || slots[boundary][side].is_some()
                    || camera.emitters[boundary].is_some()
                    || camera
                        .diffuse_thermal_emitters
                        .get(boundary)
                        .is_some_and(Option::is_some)
                    || camera
                        .coating_thermal_emitters
                        .get(boundary)
                        .is_some_and(|faces| faces.iter().any(Option::is_some))
                {
                    return Err(invalid.into());
                }
                slots[boundary][side] = Some(source_index);
            }
        }
        cx.checkpoint()?;
        Ok(Self {
            camera,
            sources,
            slots,
            microsurface,
        })
    }

    pub const fn source_camera(&self) -> &'a CurvedCameraScene<'s, 'm> {
        self.camera
    }
    pub const fn thermal_emitters(&self) -> &'a [SmithThermalEmitter] {
        self.sources
    }
    pub const fn microsurface_settings(&self) -> SmithWalkSettings {
        self.microsurface
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, SmithThermalCameraError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, SmithThermalCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, SmithThermalCameraError> {
        cx.checkpoint()?;
        let beam = self.camera.launch_sequential(optical, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CameraRadianceSample, SmithThermalCameraError> {
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            |event| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c32],
                    [seed as u32, (seed >> 32) as u32],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            |event, facet| smith_draw(seed, path_id, event, facet),
            |event, facet| {
                philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [seed as u32 ^ 0x6469_656c, (seed >> 32) as u32 ^ 0x6563_7472],
                )
                .map(|v| (f64::from(v) + 0.5) / 4_294_967_296.0)
            },
            // Independent effective Philox key from both continuation owners;
            // a variable thermal facet count cannot shift the physical path.
            |event, facet| smith_draw(seed ^ 0x7468_6572_6d61_6c31, path_id, event, facet),
        )
    }

    /// Independent zero-based surface, metal, glass and thermal facet streams.
    /// The thermal coordinates have the native opaque `[flight,u,v]` meaning.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        thermal: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, SmithThermalCameraError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            dielectric,
            thermal,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        thermal: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, SmithThermalCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            dielectric,
            thermal,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
        thermal: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, SmithThermalCameraError> {
        cx.checkpoint()?;
        let beam = self.camera.launch_sequential(optical, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            opaque,
            dielectric,
            thermal,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        mut dielectric: impl FnMut(u32, u32) -> [f64; 4],
        mut thermal: impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<CameraRadianceSample, SmithThermalCameraError> {
        let mut source_error = None;
        let sample = SmithCameraScene::new(self.camera, self.microsurface)
            .sample_launched_with_emission(
                cx,
                beam,
                initial_medium,
                settings,
                surfaces,
                opaque,
                Some(&mut dielectric),
                |path, initial_index| {
                    self.emitted(cx, path, initial_index, &mut thermal)
                        .map_err(|error| {
                            source_error = Some(error);
                            CameraRadianceError::NumericalRange
                        })
                },
            );
        match source_error {
            Some(error) => Err(error),
            None => sample.map_err(Into::into),
        }
    }

    fn emitted(
        &self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        initial_index: f64,
        thermal: &mut impl FnMut(u32, u32) -> [f64; 3],
    ) -> Result<f64, SmithThermalCameraError> {
        let contact = path
            .last_contact()
            .ok_or(CameraRadianceError::NumericalRange)?;
        let boundary = contact.witness.boundary_id;
        let beam = *path.beam();
        let side = usize::from(contact.witness.normal.dot(beam.ray().dir) > 0.0);
        let Some(source) = self.slots.get(boundary).and_then(|faces| faces[side]) else {
            return Ok(0.0);
        };
        let source = &self.sources[source];
        let native = self
            .camera
            .scene
            .bindings
            .get(boundary)
            .ok_or(CameraRadianceError::NumericalRange)?;
        if native.object_id != source.object_id || native.object_id != contact.object_id {
            return Err(CameraRadianceError::NumericalRange.into());
        }
        let unit = normalized_camera(beam, Stokes::unpolarized(1.0))
            .map_err(CurvedSceneError::from)
            .map_err(CameraRadianceError::from)?;
        let ordinal = path
            .surface_events()
            .checked_sub(1)
            .ok_or(CameraRadianceError::NumericalRange)?;
        let draw = |facet| thermal(ordinal, facet);
        let observed = match native.interaction {
            Interaction::SmoothConductor { medium, .. }
            | Interaction::RoughConductor { medium, .. }
                if medium != path.medium_index() =>
            {
                return Err(CameraRadianceError::NumericalRange.into());
            }
            Interaction::SmoothConductor { .. } => source.source.sample(
                cx,
                unit,
                contact.witness.normal,
                ConductorSurface::SMOOTH,
                self.microsurface,
                draw,
            ),
            Interaction::RoughConductor { surface, .. } => source.source.sample(
                cx,
                unit,
                contact.witness.normal,
                surface,
                self.microsurface,
                draw,
            ),
            Interaction::BrushedConductor { medium, detail, .. } => {
                if medium != path.medium_index() {
                    return Err(CameraRadianceError::NumericalRange.into());
                }
                let instance = self
                    .camera
                    .scene
                    .geometry
                    .instances()
                    .get(boundary)
                    .ok_or(CameraRadianceError::NumericalRange)?;
                let transform = instance.transform();
                let local = transform
                    .inverse()
                    .map_err(SceneError::from)
                    .map_err(CurvedSceneError::from)
                    .map_err(CameraRadianceError::from)?
                    .transform_point(contact.witness.surface_point);
                let frame = unit
                    .smith_instanced_brushed_frame(
                        &detail,
                        transform,
                        local,
                        contact.witness.normal,
                    )
                    .map_err(SceneError::from)
                    .map_err(CurvedSceneError::from)
                    .map_err(CameraRadianceError::from)?;
                if let Some(frame) = frame {
                    let (alpha_x, alpha_y) = detail.alpha_pair_for_mark(frame.mark_covered);
                    source.source.sample_anisotropic(
                        cx,
                        unit,
                        frame.shading_normal,
                        frame.tangent,
                        alpha_x,
                        alpha_y,
                        self.microsurface,
                        draw,
                    )
                } else {
                    let surface = ConductorSurface::try_rough(detail.alpha_pair().0)
                        .map_err(|_| CameraRadianceError::NumericalRange)?;
                    source.source.sample(
                        cx,
                        unit,
                        contact.witness.normal,
                        surface,
                        self.microsurface,
                        draw,
                    )
                }
            }
            _ => return Err(CameraRadianceError::NumericalRange.into()),
        }
        .map_err(|error| SmithThermalCameraError::Thermal {
            object_id: source.object_id,
            error,
        })?;
        // The source is already contracted with UNIT reciprocal polarization.
        // Only the original physical prefix and camera/local radiance ratio
        // remain. No cosine, Pi, additional emissivity or second analyzer.
        let local_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        super::weighted_source(
            observed.radiance_w_m2_sr_nm(),
            initial_index / local_index,
            beam.power().i,
            1.0,
        )
        .map_err(Into::into)
    }
}

#[cfg(test)]
mod tests;
