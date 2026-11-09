//! Calibrated point-light connections on complete Smith camera paths.
//!
//! Point fixtures have zero probability of a random continuation hit, so their
//! directional Smith estimates add to ordinary source endpoints without MIS.
//! Both receiver evaluation and continuation use the SAME multiple-reflection
//! law, with independent counter domains. A stochastic BRDF estimate is never
//! treated as an outgoing proposal PDF or inserted into a nonlinear weight.
//!
//! Reuses native source photometry, reciprocal source polarization, bundle
//! shooting, Beer/GRIN attenuation, geometric receiver irradiance and n²
//! transport. Surface/volume LTE emission remains on the original camera path.
//! Only the original declared source branches are covered: no new image search,
//! area/environment NEE, molecular scattering, or coherent microgeometry.

use super::super::CameraAnalyzer;
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::{CurvedTransportEvent, SmithSurfaceSampler};
use crate::polarization::dielectric::scene::smith_draw;
use crate::polarization::dielectric::transport::PolarizedBeam;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use fs_rand::philox::philox4x32_10;

pub mod dielectric;
pub use dielectric::SmithDielectricPointLightCameraScene;

/// Borrowed, allocation-free opt-in view of an already admitted point-light
/// camera. The per-surface facet limit also bounds each fixed-direction solve.
pub struct SmithPointLightCameraScene<'v, 'a, 's, 'm, 'l> {
    source: &'v CurvedPointLightCameraScene<'a, 's, 'm, 'l>,
    microsurface: SmithWalkSettings,
}
impl<'v, 'a, 's, 'm, 'l> SmithPointLightCameraScene<'v, 'a, 's, 'm, 'l> {
    pub const fn new(
        source: &'v CurvedPointLightCameraScene<'a, 's, 'm, 'l>,
        microsurface: SmithWalkSettings,
    ) -> Self {
        Self {
            source,
            microsurface,
        }
    }

    pub const fn source_point_camera(&self) -> &'v CurvedPointLightCameraScene<'a, 's, 'm, 'l> {
        self.source
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
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
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.source.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    /// Retain the resolved lens's reciprocal Stokes seed; its raster owner
    /// still owns the pupil/geometric weight.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, seed, path_id)
    }

    fn sample_launched(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
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
            |event, binding, facet| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [
                        (seed as u32) ^ 0x736d_6e65,
                        ((seed >> 32) as u32) ^ binding as u32,
                    ],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
        )
    }

    /// Caller-supplied independent streams: ordinary surfaces, continuing
    /// (surface, facet), and light evaluation (surface, binding, facet).
    /// A failed/budget-exhausted estimate refuses the complete camera sample;
    /// callers must not discard failures and average only successful paths.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, usize, u32) -> [f64; 3],
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            facet_samples,
            connection_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        connection_samples: impl FnMut(u32, usize, u32) -> [f64; 3],
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        self.sample_launched_with_microsurfaces(
            cx,
            beam,
            initial_medium,
            settings,
            surface_samples,
            facet_samples,
            connection_samples,
            None,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_microsurfaces(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        mut connection_samples: impl FnMut(u32, usize, u32) -> [f64; 3],
        mut dielectric_samples: Option<&mut dyn FnMut(u32, u32) -> [f64; 4]>,
        mut dielectric_connections: Option<&mut dyn FnMut(u32, usize, u32) -> [f64; 4]>,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        validate_settings(settings)?;
        let camera = self.source.camera;
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(camera, initial_index, beam);
        let mut fixtures = FixtureContributions::new(self.source, settings, initial_index);
        let mut connector = ClearFixtureConnector {
            scene: camera.scene,
        };
        let result = camera
            .scene
            .trace_with_flights_events_and_smith(
                cx,
                beam,
                initial_medium,
                settings.surfaces,
                surface_samples,
                true,
                &camera.spatial_volume_emitters,
                |event| match event {
                    CurvedTransportEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                        emission.observe(path)?;
                        fixtures
                            .add_contact_with_response(
                                cx,
                                path,
                                &mut connector,
                                |receiver, connection, binding| {
                                    if let Some(samples) = dielectric_connections.as_mut() {
                                        if receiver.uses_smith_dielectric(camera)? {
                                            return receiver
                                                .smith_dielectric_response(
                                                    cx,
                                                    camera,
                                                    connection,
                                                    binding,
                                                    self.microsurface,
                                                    |facet| {
                                                        samples(
                                                            path.surface_events() - 1,
                                                            binding,
                                                            facet,
                                                        )
                                                    },
                                                )
                                                .map_err(Into::into);
                                        }
                                    }
                                    receiver
                                        .smith_response(
                                            cx,
                                            camera,
                                            connection,
                                            binding,
                                            self.microsurface,
                                            |facet| {
                                                connection_samples(
                                                    path.surface_events(),
                                                    binding,
                                                    facet,
                                                )
                                            },
                                        )
                                        .map_err(Into::into)
                                },
                            )
                            .map_err(|error| {
                                fixtures.error = Some(error);
                                CurvedSceneError::NumericalRange
                            })
                    }
                    CurvedTransportEvent::Surface(CurvedTraceEvent::Departure { path, .. })
                    | CurvedTransportEvent::MicrosurfaceDeparture(path) => {
                        emission.depart(path);
                        Ok(())
                    }
                },
                0,
                |_, _, _| Ok(None),
                |_, _, _, _| Err(CurvedSceneError::NumericalRange),
                |_, _, _| Ok(()),
                Some(SmithSurfaceSampler {
                    settings: self.microsurface,
                    samples: &mut facet_samples,
                    dielectric_samples: dielectric_samples
                        .as_mut()
                        .map(|samples| &mut **samples as &mut dyn FnMut(u32, u32) -> [f64; 4]),
                }),
            )
            .map(|(path, _)| path);
        self.source
            .finish(cx, initial_index, emission, fixtures, result)
    }
}
