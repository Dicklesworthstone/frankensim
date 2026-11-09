//! Environment NEE on the original polarized participating camera walk.
//!
//! Accepted molecular and retained Mie collisions use their original spectral
//! phase response and angular PDF. Null candidates keep the previous continuous
//! anchor. Original grid cell traversals integrate source-leg scattering and
//! separately bound true absorption exactly; their finite support admits escape.
//! Arbitrary active fields and GRIN tails keep the original camera estimator.
//!
//! MIS uses the same deterministic angular partition on both techniques.
//! Source survival multiplies radiance, never a nonlinear weight. Camera free
//! flights already sample their own survival; endpoint radiance does not receive
//! another attenuation factor. Delta interactions remove the continuous anchor.

use super::super::super::participating::CurvedRayleighCameraScene;
use super::source::EnvironmentSourceSample;
use super::weights::positive_product_with_log_survival;
use super::{
    CameraAnalyzer, CameraRadianceError, CameraRadianceSample, CurvedEnvironmentCameraScene,
    EnvironmentCameraError, EnvironmentCameraSample, EnvironmentCameraSettings, EnvironmentRadiance,
    RadianceEndpoint, SurfaceEmission, SurfaceReceiver, Work, clear_to_environment, counter_samples,
    finite_nonnegative, phase_index, positive_product, power_weight, supports_segment,
};
use crate::charts::Ray;
use crate::polarization::dielectric::scene::SceneError;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighEvent, CurvedRayleighTraceSettings,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedContact, CurvedPathResult, CurvedSceneError, CurvedTraceEvent,
};
use crate::polarization::dielectric::transport::PolarizedBeam;
use crate::volumes::transport::scattering::polarized::PolarizedScatteringError;
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// Limits remain cumulative across all surface and physical volume vertices.
/// Zero source-cell work still admits absent/zero laws and known uniform tails.
#[derive(Clone, Copy, Debug)]
pub struct ParticipatingEnvironmentCameraSettings {
    pub environment: EnvironmentCameraSettings,
    pub max_volume_events: u32,
    pub max_source_cells: u32,
}

impl Default for ParticipatingEnvironmentCameraSettings {
    fn default() -> Self {
        Self {
            environment: EnvironmentCameraSettings::default(),
            max_volume_events: 64,
            max_source_cells: 65_536,
        }
    }
}

/// Original path and collision counters, direct environment radiance, and the
/// actual work used by deterministic source-grid integration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticipatingEnvironmentCameraSample {
    environment: EnvironmentCameraSample,
    volume_events: u32,
    free_flight_draws: u32,
    source_cell_visits: u32,
}

impl ParticipatingEnvironmentCameraSample {
    pub const fn environment_sample(&self) -> &EnvironmentCameraSample {
        &self.environment
    }
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        self.environment.radiance_sample()
    }
    pub const fn environment_radiance_w_m2_sr_nm(&self) -> f64 {
        self.environment.environment_radiance_w_m2_sr_nm()
    }
    pub const fn environment_endpoint_weight(&self) -> f64 {
        self.environment.environment_endpoint_weight()
    }
    pub const fn volume_events(&self) -> u32 {
        self.volume_events
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
    pub const fn source_cell_visits(&self) -> u32 {
        self.source_cell_visits
    }
}

/// Borrow the already identity-checked surface/scattering composition. No grid,
/// phase table, environment CDF, or mutable per-scene workspace is copied.
pub struct CurvedParticipatingEnvironmentCameraScene<'a, 'c, 's, 'm, 'o> {
    camera: &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o>,
    environment: EnvironmentRadiance<'s>,
}

impl<'a, 'c, 's, 'm, 'o> CurvedParticipatingEnvironmentCameraScene<'a, 'c, 's, 'm, 'o> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o>,
    ) -> Result<Self, EnvironmentCameraError> {
        let environment =
            CurvedEnvironmentCameraScene::try_new(cx, camera.surface_camera())?.environment;
        Ok(Self { camera, environment })
    }

    pub const fn source_camera(&self) -> &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o> {
        self.camera
    }

    /// Preserve native surface/flight streams. Surface environment directions
    /// keep the surface-only domain; physical volume vertices use a new domain.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: ParticipatingEnvironmentCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<ParticipatingEnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_ray_with_tracking_samples(
            cx,
            ray,
            wavelength_nm,
            initial_medium,
            settings,
            |event| counter_samples(seed, path_id, event, 0x706f_6c32),
            |flight| {
                philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, flight, 0x706f_6c34],
                    [seed as u32, (seed >> 32) as u32],
                )
                .map(|lane| (f64::from(lane) + 0.5) / 4_294_967_296.0)
            },
            |event| {
                let sample = counter_samples(seed, path_id, event, 0x706f_6c39);
                [sample[0], sample[1]]
            },
            |event| {
                let sample = counter_samples(seed, path_id, event, 0x706f_6c3a);
                [sample[0], sample[1]]
            },
        )
    }

    /// Flight variates retain native order: distance, phase cosine, azimuth,
    /// independent delta acceptance. The two environment streams are indexed
    /// by original surface events and accepted physical volume events.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: ParticipatingEnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_environment_samples: impl FnMut(u32) -> [f64; 2],
        volume_environment_samples: impl FnMut(u32) -> [f64; 2],
    ) -> Result<ParticipatingEnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_with_launch(
            cx,
            || self.camera.surface_camera().launch(ray, wavelength_nm, initial_medium),
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            surface_environment_samples,
            volume_environment_samples,
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
        settings: ParticipatingEnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_environment_samples: impl FnMut(u32) -> [f64; 2],
        volume_environment_samples: impl FnMut(u32) -> [f64; 2],
    ) -> Result<ParticipatingEnvironmentCameraSample, EnvironmentCameraError> {
        self.sample_with_launch(
            cx,
            || analyzer.launch(self.camera.surface_camera(), ray, wavelength_nm, initial_medium),
            initial_medium,
            settings,
            surface_samples,
            flight_samples,
            surface_environment_samples,
            volume_environment_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_with_launch(
        &self,
        cx: &Cx<'_>,
        launch: impl FnOnce() -> Result<PolarizedBeam, CameraRadianceError>,
        initial_medium: usize,
        settings: ParticipatingEnvironmentCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        mut surface_environment_samples: impl FnMut(u32) -> [f64; 2],
        mut volume_environment_samples: impl FnMut(u32) -> [f64; 2],
    ) -> Result<ParticipatingEnvironmentCameraSample, EnvironmentCameraError> {
        cx.checkpoint()?;
        settings.environment.validate()?;
        if settings.max_volume_events == 0 {
            return Err(EnvironmentCameraError::InvalidSettings);
        }
        let camera = self.camera.surface_camera();
        let beam = launch()?;
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(camera, initial_index, beam);
        let mut work = Work::default();
        let mut source_cells = 0;
        let mut volume_ordinal = 0_u32;
        let mut direct = 0.0;
        let mut continuous = false;
        let mut continuation_pdf = None;
        let mut source_error = None;
        let active = self.environment.scale_w_m2_sr_nm > 0.0 && !self.environment.map.is_black();
        let result = self.camera.participating_scene().for_camera()
            .trace_with_tracking_samples_events_observed(
                cx,
                beam,
                initial_medium,
                CurvedRayleighTraceSettings {
                    surfaces: settings.environment.surfaces,
                    max_volume_events: settings.max_volume_events,
                },
                surface_samples,
                flight_samples,
                &camera.spatial_volume_emitters,
                |event| {
                    // Preserve the emission owner's detailed error; its finish
                    // method resolves the stored failure after traversal stops.
                    match &event {
                        CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                            emission.observe(path)?;
                        }
                        CurvedRayleighEvent::VolumeArrival { path, .. } => {
                            emission.observe_volume(path)?;
                        }
                        CurvedRayleighEvent::Surface(CurvedTraceEvent::Departure { path, .. })
                        | CurvedRayleighEvent::VolumeDeparture { path, .. } => emission.depart(path),
                    }
                    let observed: Result<(), EnvironmentCameraError> = (|| {
                        match event {
                            CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                                continuation_pdf = None;
                                let receiver = SurfaceReceiver::at_arrival(camera, path)?;
                                continuous = active && receiver.is_some();
                                if let Some(receiver) = receiver.filter(|_| active) {
                                    let ordinal = path.surface_events().checked_sub(1)
                                        .ok_or(EnvironmentCameraError::NumericalRange)?;
                                    let source = self.sample_source(path, settings, &mut work,
                                        || surface_environment_samples(ordinal))?;
                                    if let Some(source) = source {
                                        let (response, medium) = receiver
                                            .directional_response(camera, source.frame, source.stokes)?;
                                        let factor = if receiver.is_rough() {
                                            response.factor
                                        } else {
                                            response.factor / core::f64::consts::PI
                                        };
                                        let contact = *path.last_contact()
                                            .ok_or(EnvironmentCameraError::NumericalRange)?;
                                        let value = self.add_source(
                                            cx, path, initial_index, source, Some(contact), medium,
                                            factor, response.pdf,
                                            response.normal.dot(source.direction).abs().min(1.0),
                                            settings, &mut work, &mut source_cells,
                                        )?;
                                        direct = finite_nonnegative(direct + value)?;
                                    }
                                }
                            }
                            CurvedRayleighEvent::Surface(CurvedTraceEvent::Departure {
                                path, pdf, delta, ..
                            }) => {
                                continuation_pdf = self.continuation_pdf(
                                    path, continuous && !delta, pdf,
                                )?;
                            }
                            CurvedRayleighEvent::VolumeArrival { path, phase } => {
                                continuation_pdf = None;
                                continuous = active && self.supports_tail(path.medium_index(),
                                    path.beam().wavelength_nm())?;
                                let ordinal = volume_ordinal;
                                volume_ordinal = volume_ordinal.checked_add(1)
                                    .ok_or(EnvironmentCameraError::NumericalRange)?;
                                if continuous {
                                    let source = self.sample_source(path, settings, &mut work,
                                        || volume_environment_samples(ordinal))?;
                                    if let Some(source) = source {
                                        let factor = phase.polarized_connection_factor(
                                            *path.beam(), source.frame, source.stokes,
                                        ).map_err(scattering_error)?;
                                        let pdf = phase.eval_at(path.beam().wavelength_nm(),
                                            path.beam().ray().dir.dot(source.direction))
                                            .map_err(|error| scattering_error(error.into()))?;
                                        let value = self.add_source(
                                            cx, path, initial_index, source, None, path.medium_index(),
                                            factor, pdf, 1.0, settings, &mut work, &mut source_cells,
                                        )?;
                                        direct = finite_nonnegative(direct + value)?;
                                    }
                                }
                            }
                            CurvedRayleighEvent::VolumeDeparture { path, pdf, .. } => {
                                continuation_pdf = self.continuation_pdf(path, continuous, pdf)?;
                            }
                        }
                        Ok(())
                    })();
                    observed.map_err(|error| {
                        source_error = Some(error);
                        CurvedSceneError::NumericalRange
                    })
                },
            );
        if let Some(error) = source_error {
            return Err(error);
        }
        let (path, emitted) = emission.finish_scattering(result)?;
        let native_path = *path.path();
        let mut radiance = camera.finish(cx, native_path, initial_index, 0.0)?;
        let endpoint_weight = if radiance.endpoint() == RadianceEndpoint::Environment {
            if let Some(pdf) = continuation_pdf {
                let source = self.environment.map.evaluate(native_path.beam().ray().dir)
                    .ok_or(EnvironmentCameraError::NumericalRange)?;
                power_weight(pdf, source.pdf_solid_angle)?
            } else {
                1.0
            }
        } else {
            1.0
        };
        let endpoint = positive_product(&[radiance.radiance_w_m2_sr_nm(), endpoint_weight], 1.0)?;
        radiance.radiance_w_m2_sr_nm = finite_nonnegative(emitted + direct + endpoint)?;
        cx.checkpoint()?;
        Ok(ParticipatingEnvironmentCameraSample {
            environment: EnvironmentCameraSample {
                radiance,
                direct_environment: direct,
                endpoint_weight,
                work,
            },
            volume_events: path.volume_events(),
            free_flight_draws: path.free_flight_draws(),
            source_cell_visits: source_cells,
        })
    }

    fn supports_tail(&self, medium: usize, wavelength: f64)
        -> Result<bool, EnvironmentCameraError>
    {
        if !supports_segment(self.camera.surface_camera().scene, medium, wavelength)? {
            return Ok(false);
        }
        self.camera.participating_scene().supports_environment_tail(medium, wavelength)
            .map_err(CameraRadianceError::from)
            .map_err(EnvironmentCameraError::from)
    }

    fn continuation_pdf(&self, path: &CurvedPathResult, continuous: bool, pdf: f64)
        -> Result<Option<f64>, EnvironmentCameraError>
    {
        if !continuous || !self.supports_tail(path.medium_index(), path.beam().wavelength_nm())? {
            return Ok(None);
        }
        if !pdf.is_finite() || pdf <= 0.0 {
            return Err(EnvironmentCameraError::NumericalRange);
        }
        Ok(Some(pdf))
    }

    fn sample_source(
        &self,
        path: &CurvedPathResult,
        settings: ParticipatingEnvironmentCameraSettings,
        work: &mut Work,
        samples: impl FnOnce() -> [f64; 2],
    ) -> Result<Option<EnvironmentSourceSample>, EnvironmentCameraError> {
        if work.samples >= settings.environment.max_environment_samples {
            return Err(EnvironmentCameraError::SampleLimit);
        }
        work.samples += 1;
        self.environment.sample_direction(path.beam().wavelength_nm(), samples())
    }

    #[allow(clippy::too_many_arguments)]
    fn add_source(
        &self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        initial_index: f64,
        source: EnvironmentSourceSample,
        contact: Option<CurvedContact>,
        medium: usize,
        factor: f64,
        continuation_pdf: f64,
        receiving_cosine: f64,
        settings: ParticipatingEnvironmentCameraSettings,
        work: &mut Work,
        source_cells: &mut u32,
    ) -> Result<f64, EnvironmentCameraError> {
        if factor == 0.0 || receiving_cosine == 0.0
            || !self.supports_tail(medium, path.beam().wavelength_nm())?
        {
            return Ok(0.0);
        }
        let ray = Ray { origin: path.beam().ray().origin, dir: source.direction };
        work.shadows = work.shadows.checked_add(1)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        let remaining = settings.environment.max_query_work.checked_sub(work.geometry)
            .ok_or(EnvironmentCameraError::GeometryLimit)?;
        let visibility = clear_to_environment(
            cx, self.camera.surface_camera().scene, contact, ray, source.frame.horizontal(),
            path.beam().wavelength_nm(), medium, settings.environment, remaining,
        )?;
        work.fields = work.fields.checked_add(visibility.fields)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        work.geometry = work.geometry.checked_add(visibility.geometry)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        if !visibility.visible {
            return Ok(0.0);
        }
        let remaining_cells = settings.max_source_cells.checked_sub(*source_cells)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        let survival = self.camera.participating_scene().environment_scattering_survival(
            cx, medium, ray, path.beam().wavelength_nm(), remaining_cells,
        ).map_err(CameraRadianceError::from)?
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        *source_cells = source_cells.checked_add(survival.cell_visits)
            .ok_or(EnvironmentCameraError::NumericalRange)?;
        if survival.log_survival == f64::NEG_INFINITY {
            // This marker is the admitted uniform infinite tail's exact zero,
            // not underflow of a finite optical depth or an unresolved path.
            return Ok(0.0);
        }
        let weight = power_weight(source.pdf, continuation_pdf)?;
        let local_index = phase_index(path.beam().medium(), path.beam().wavelength_nm())?;
        let ratio = initial_index / local_index;
        positive_product_with_log_survival(
            &[
                path.beam().power().i, ratio, ratio, factor,
                source.radiance, receiving_cosine, weight,
            ],
            source.pdf,
            survival.log_survival,
        )
    }
}

fn scattering_error(error: PolarizedScatteringError) -> EnvironmentCameraError {
    CameraRadianceError::from(CurvedSceneError::from(SceneError::from(error))).into()
}
