//! Reciprocal complete Smith surfaces interleaved with physical volume events.
//!
//! This path-only owner observes original surface, environment and LTE sources
//! on one native participating path. The same rough-surface Mueller walk acts
//! before or after any accepted molecular/particle collision. Every microscopic
//! escape refreshes the thermal source prefix and uses its actual final medium
//! for subsequent absorption, free flight and endpoint radiance conversion.
//! No hidden microsurface density enters MIS and no scattering survival factor
//! is applied twice. Point/area/environment NEE remains a separate capability.

use super::super::CameraAnalyzer;
use super::{
    CameraRadianceError, CurvedRayleighCameraScene, CurvedRayleighTraceSettings, PolarizedBeam,
    Ray, RayleighCameraSample, SurfaceEmission, phase_index,
};
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighTransportEvent;
use crate::polarization::dielectric::scene::curved::participating::microsurface::SmithDielectricRayleighScene;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use fs_exec::Cx;

/// Opt-in view retaining the exact admitted source camera and volume scene.
/// Construction adds no heap workspace and changes no source/optical bindings.
pub struct SmithDielectricRayleighCameraScene<'v, 'a, 's, 'm, 'o> {
    camera: &'v CurvedRayleighCameraScene<'a, 's, 'm, 'o>,
    microsurface: SmithWalkSettings,
}

impl<'v, 'a, 's, 'm, 'o> SmithDielectricRayleighCameraScene<'v, 'a, 's, 'm, 'o> {
    pub const fn new(
        camera: &'v CurvedRayleighCameraScene<'a, 's, 'm, 'o>,
        microsurface: SmithWalkSettings,
    ) -> Self {
        Self {
            camera,
            microsurface,
        }
    }

    pub const fn source_camera(&self) -> &'v CurvedRayleighCameraScene<'a, 's, 'm, 'o> {
        self.camera
    }

    pub const fn microsurface_settings(&self) -> SmithWalkSettings {
        self.microsurface
    }

    /// Counter-addressed seeded path. Zero scattering retains the exact source
    /// surface and Smith streams of `SmithDielectricCameraScene`.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self
            .camera
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// Install the complete passive detector row before all volume and facet
    /// Mueller events; no scalar post-processing analyzer is substituted.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.camera.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// Preserve the actual wavelength-specific resolved lens beam and its
    /// reciprocal Stokes seed. Pupil geometry stays with the raster owner.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self
            .camera
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera.camera, initial_index);
        let scene = self.camera.scattering.for_camera();
        let smith = SmithDielectricRayleighScene::new(&scene, self.microsurface);
        let result = smith.trace_events_observed(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            &self.camera.camera.spatial_volume_emitters,
            |event| match event {
                CurvedRayleighTransportEvent::Native(event) => {
                    emission.observe_rayleigh_event(event)
                }
                CurvedRayleighTransportEvent::MicrosurfaceDeparture(path) => {
                    emission.depart(path);
                    Ok(())
                }
            },
        );
        let (path, emitted) = emission.finish_scattering(result)?;
        self.camera.finish(cx, path, initial_index, emitted)
    }

    /// Zero-based independent surface, physical-flight, opaque-facet and
    /// dielectric-facet streams. Physical draws retain all four native tracking
    /// coordinates, so heterogeneous and local spectral laws stay supported.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        flights: impl FnMut(u32) -> [f64; 4],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self
            .camera
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_tracking_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            flights,
            opaque,
            dielectric,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        flights: impl FnMut(u32) -> [f64; 4],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.camera.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_tracking_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            flights,
            opaque,
            dielectric,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        flights: impl FnMut(u32) -> [f64; 4],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        cx.checkpoint()?;
        let beam = self
            .camera
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched_with_tracking_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            flights,
            opaque,
            dielectric,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        surfaces: impl FnMut(u32) -> [f64; 3],
        flights: impl FnMut(u32) -> [f64; 4],
        opaque: impl FnMut(u32, u32) -> [f64; 3],
        dielectric: impl FnMut(u32, u32) -> [f64; 4],
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera.camera, initial_index);
        let scene = self.camera.scattering.for_camera();
        let smith = SmithDielectricRayleighScene::new(&scene, self.microsurface);
        let result = smith.trace_with_tracking_samples_events_observed(
            cx,
            beam,
            initial_medium,
            settings,
            surfaces,
            flights,
            opaque,
            dielectric,
            &self.camera.camera.spatial_volume_emitters,
            |event| match event {
                CurvedRayleighTransportEvent::Native(event) => {
                    emission.observe_rayleigh_event(event)
                }
                CurvedRayleighTransportEvent::MicrosurfaceDeparture(path) => {
                    emission.depart(path);
                    Ok(())
                }
            },
        );
        let (path, emitted) = emission.finish_scattering(result)?;
        self.camera.finish(cx, path, initial_index, emitted)
    }
}
