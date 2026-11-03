//! Calibrated point lighting on the original molecular camera walk.
//!
//! Camera free-flight probabilities already include scattering survival. This
//! owner adds no extinction or collision-rate multiplier to that path. At each
//! selected diffuse or molecular vertex it connects the original point fixture and
//! independently estimates source-leg scattering survival on the native curved
//! bundle. Absorption, Fresnel transfer and the source solid-angle Jacobian stay
//! in the original solved bundle; the tracked replay supplies survival only.
//! It must revalidate that same optical branch and receiver before a sample is
//! accepted. The retained numerical trajectory is not an exact-path certificate.
//!
//! Uniform source scattering uses its actual curved arc length. Heterogeneous
//! source fields use bounded ratio tracking, including exact stochastic zero
//! factors. Those zeros are valid estimates, not missing or failed connections.
//! Positive survival remains logarithmic until it is combined with calibrated
//! irradiance and the camera prefix. Source draws have their own stream and do
//! not consume or change camera surface/volume variates.
//!
//! The original emitter observation and diffuse continuation happen on the same
//! camera path, including after molecular collisions. Optional volume branches
//! add the full reciprocal Stokes contraction at accepted molecular vertices,
//! before sampling their original continuation. Null candidates add no term.

use super::super::participating::CurvedRayleighCameraScene;
use super::{
    CameraRadianceError, CameraRadianceSample, ConnectionWork, CurvedCameraScene,
    CurvedPointLightCameraScene, CurvedSceneError, FixtureConnection, FixtureConnectionRequest,
    FixtureConnector, FixtureContributions, PointLightCameraError, PointLightCameraSample,
    PointLightCameraSettings, SurfaceEmission, phase_index, validate_settings,
};
use crate::charts::Ray;
use crate::polarization::dielectric::scene::SceneError;
use crate::polarization::dielectric::scene::curved::CurvedTraceEvent;
use crate::polarization::dielectric::scene::curved::participating::transmittance::{
    CurvedTransmittanceError, CurvedTransmittanceOptions,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighEvent, CurvedRayleighPath, CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

pub mod images;
mod vertices;
pub use vertices::MolecularPointLightBranch;

/// One canonical camera-surface configuration, the original volume-event cap,
/// and cumulative source tracking limits across all selected surface and volume vertices.
/// `point_lights.max_trace_evaluations` counts both geometric source solves and
/// their tracked replays. Zero tracking limits still admit clear/uniform source
/// legs that need no random flight or candidate; limits are charged on work.
#[derive(Clone, Copy, Debug)]
pub struct RayleighPointLightCameraSettings {
    pub point_lights: PointLightCameraSettings,
    pub max_volume_events: u32,
    pub max_source_free_flight_draws: u32,
    pub max_source_candidates: u32,
}

impl Default for RayleighPointLightCameraSettings {
    fn default() -> Self {
        Self {
            point_lights: PointLightCameraSettings::default(),
            max_volume_events: 64,
            max_source_free_flight_draws: 4096,
            max_source_candidates: 4096,
        }
    }
}

/// Complete camera sample, its original real/null molecular-flight counters,
/// and independently counted source-connection work. No partial path or pixel
/// is published on a camera, connection, tracking, or cancellation error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighPointLightCameraSample {
    point_lights: PointLightCameraSample,
    volume_events: u32,
    free_flight_draws: u32,
}

impl RayleighPointLightCameraSample {
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        self.point_lights.radiance_sample()
    }
    pub const fn fixture_radiance_w_m2_sr_nm(&self) -> f64 {
        self.point_lights.fixture_radiance_w_m2_sr_nm()
    }
    pub const fn connections(&self) -> u32 {
        self.point_lights.connections()
    }
    pub const fn trace_evaluations(&self) -> u32 {
        self.point_lights.trace_evaluations()
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.point_lights.field_evaluations()
    }
    pub const fn query_work(&self) -> u64 {
        self.point_lights.query_work()
    }
    pub const fn volume_events(&self) -> u32 {
        self.volume_events
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
    pub const fn source_free_flight_draws(&self) -> u32 {
        self.point_lights.work.source_draws
    }
    pub const fn source_candidates(&self) -> u32 {
        self.point_lights.work.source_candidates
    }
}

/// Borrow the two original admitted owners, with no copied scene or optical
/// tables. Their base camera must be the exact same original object, including
/// emitter and environment bindings, not merely a camera over equal geometry.
pub struct CurvedRayleighPointLightCameraScene<'a, 'c, 's, 'm, 'o, 'l> {
    camera: &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o>,
    point_lights: &'a CurvedPointLightCameraScene<'c, 's, 'm, 'l>,
    volume_branches: &'a [MolecularPointLightBranch<'l>],
}

impl<'a, 'c, 's, 'm, 'o, 'l> CurvedRayleighPointLightCameraScene<'a, 'c, 's, 'm, 'o, 'l> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o>,
        point_lights: &'a CurvedPointLightCameraScene<'c, 's, 'm, 'l>,
    ) -> Result<Self, PointLightCameraError> {
        Self::try_new_with_volume_branches(cx, camera, point_lights, &[])
    }

    /// Add one declared optical branch per original fixture/receiving medium.
    /// The branch's shooting seed must resolve each sampled volume vertex; an
    /// unresolved or blocked connection refuses the entire sample. No global
    /// multiple-image completeness or area-emitter MIS is inferred.
    pub fn try_new_with_volume_branches(
        cx: &Cx<'_>,
        camera: &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o>,
        point_lights: &'a CurvedPointLightCameraScene<'c, 's, 'm, 'l>,
        volume_branches: &'a [MolecularPointLightBranch<'l>],
    ) -> Result<Self, PointLightCameraError> {
        cx.checkpoint()?;
        if !core::ptr::eq(camera.surface_camera(), point_lights.source_camera()) {
            return Err(CameraRadianceError::Scene(SceneError::InvalidBindings.into()).into());
        }
        vertices::validate_branches(cx, point_lights, volume_branches)?;
        Ok(Self {
            camera,
            point_lights,
            volume_branches,
        })
    }

    pub const fn volume_branches(&self) -> &'a [MolecularPointLightBranch<'l>] {
        self.volume_branches
    }

    pub const fn source_camera(&self) -> &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o> {
        self.camera
    }
    pub const fn point_lights(&self) -> &'a CurvedPointLightCameraScene<'c, 's, 'm, 'l> {
        self.point_lights
    }
    pub const fn surface_camera(&self) -> &'c CurvedCameraScene<'s, 'm> {
        self.camera.surface_camera()
    }

    /// Retain both original camera streams. Source tracking uses a third,
    /// independent domain keyed by connection ordinal and replay draw ordinal.
    /// A geometric solve consumes no source draws; they begin on its validated
    /// branch replay, and zero/uniform scattering consumes none.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<RayleighPointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        validate_settings(settings.point_lights)?;
        let beam = self
            .surface_camera()
            .launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let mut emission = SurfaceEmission::new(self.surface_camera(), initial_index);
        let mut fixtures =
            FixtureContributions::new(self.point_lights, settings.point_lights, initial_index);
        let mut connector = self.connector(settings, |connection, draw| {
            let block = philox4x32_10(
                [path_id as u32, (path_id >> 32) as u32, connection, draw],
                [seed as u32, (seed >> 32) as u32 ^ 0x706f_6c36],
            );
            (f64::from(block[0]) + 0.5) / 4_294_967_296.0
        });
        let result = self
            .camera
            .participating_scene()
            .for_camera()
            .trace_events_observed(
                cx,
                beam,
                initial_medium,
                camera_settings(settings),
                seed,
                path_id,
                |event| {
                    self.observe_event(cx, event, &mut emission, &mut fixtures, &mut connector)
                },
            );
        self.finish(cx, initial_index, emission, fixtures, result)
    }

    /// Uniform-density camera flights retain [distance, phase cosine, phase
    /// azimuth]. Source tracking may still be heterogeneous, using its separate
    /// scalar sampler `(connection_ordinal, replay_flight_draw_ordinal)`.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<RayleighPointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        validate_settings(settings.point_lights)?;
        let beam = self
            .surface_camera()
            .launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let mut emission = SurfaceEmission::new(self.surface_camera(), initial_index);
        let mut fixtures =
            FixtureContributions::new(self.point_lights, settings.point_lights, initial_index);
        let mut connector = self.connector(settings, source_samples);
        let result = self
            .camera
            .participating_scene()
            .for_camera()
            .trace_with_samples_events_observed(
                cx,
                beam,
                initial_medium,
                camera_settings(settings),
                surface_samples,
                flight_samples,
                |event| {
                    self.observe_event(cx, event, &mut emission, &mut fixtures, &mut connector)
                },
            );
        self.finish(cx, initial_index, emission, fixtures, result)
    }

    /// Uniform or heterogeneous camera flights retain the native fourth,
    /// independent delta-acceptance coordinate. The source sampler's first
    /// argument distinguishes every matching fixture/vertex encounter; its
    /// second counts only that source replay's scalar free-flight draws.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<RayleighPointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        validate_settings(settings.point_lights)?;
        let beam = self
            .surface_camera()
            .launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let mut emission = SurfaceEmission::new(self.surface_camera(), initial_index);
        let mut fixtures =
            FixtureContributions::new(self.point_lights, settings.point_lights, initial_index);
        let mut connector = self.connector(settings, source_samples);
        let result = self
            .camera
            .participating_scene()
            .for_camera()
            .trace_with_tracking_samples_events_observed(
                cx,
                beam,
                initial_medium,
                camera_settings(settings),
                surface_samples,
                flight_samples,
                |event| {
                    self.observe_event(cx, event, &mut emission, &mut fixtures, &mut connector)
                },
            );
        self.finish(cx, initial_index, emission, fixtures, result)
    }

    fn observe_event(
        &self,
        cx: &Cx<'_>,
        event: CurvedRayleighEvent<'_>,
        emission: &mut SurfaceEmission<'_, 's, 'm>,
        fixtures: &mut FixtureContributions<'_, 'c, 's, 'm, 'l>,
        connector: &mut impl FixtureConnector,
    ) -> Result<(), CurvedSceneError> {
        match event {
            CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                emission.observe(path)?;
                fixtures.observe_with_connector(cx, path, connector)
            }
            CurvedRayleighEvent::VolumeArrival { path, phase } => {
                emission.observe_volume(path)?;
                fixtures.observe_volume(cx, path, phase, self.volume_branches, connector)
            }
            CurvedRayleighEvent::Surface(CurvedTraceEvent::Departure { path, .. })
            | CurvedRayleighEvent::VolumeDeparture { path, .. } => {
                emission.depart(path);
                Ok(())
            }
        }
    }

    fn connector<F: FnMut(u32, u32) -> f64>(
        &self,
        settings: RayleighPointLightCameraSettings,
        samples: F,
    ) -> ParticipatingFixtureConnector<'c, 's, 'm, 'o, F> {
        ParticipatingFixtureConnector {
            scattering: self.camera.participating_scene(),
            settings,
            samples,
        }
    }

    fn finish(
        &self,
        cx: &Cx<'_>,
        initial_index: f64,
        emission: SurfaceEmission<'_, 's, 'm>,
        fixtures: FixtureContributions<'_, 'c, 's, 'm, 'l>,
        result: Result<CurvedRayleighPath, CurvedSceneError>,
    ) -> Result<RayleighPointLightCameraSample, PointLightCameraError> {
        let counters = result
            .as_ref()
            .ok()
            .map(|path| (path.volume_events(), path.free_flight_draws()));
        let point_lights = self.point_lights.finish(
            cx,
            initial_index,
            emission,
            fixtures,
            result.map(|path| *path.path()),
        )?;
        let (volume_events, free_flight_draws) =
            counters.ok_or(PointLightCameraError::NumericalRange)?;
        Ok(RayleighPointLightCameraSample {
            point_lights,
            volume_events,
            free_flight_draws,
        })
    }
}

struct ParticipatingFixtureConnector<'q, 's, 'm, 'o, F> {
    scattering: &'q CurvedRayleighScene<'s, 'm, 'o>,
    settings: RayleighPointLightCameraSettings,
    samples: F,
}

impl<F: FnMut(u32, u32) -> f64> FixtureConnector
    for ParticipatingFixtureConnector<'_, '_, '_, '_, F>
{
    fn connect(
        &mut self,
        cx: &Cx<'_>,
        request: FixtureConnectionRequest<'_>,
        previous_work: ConnectionWork,
    ) -> Result<FixtureConnection, PointLightCameraError> {
        let remaining_draws = self
            .settings
            .max_source_free_flight_draws
            .checked_sub(previous_work.source_draws)
            .ok_or(PointLightCameraError::SourceDrawLimit)?;
        let remaining_candidates = self
            .settings
            .max_source_candidates
            .checked_sub(previous_work.source_candidates)
            .ok_or(PointLightCameraError::SourceCandidateLimit)?;
        let result = self
            .scattering
            .connect_bundle_with_transmittance_samples(
                cx,
                request.launch,
                request.target,
                request.options,
                request.policy,
                CurvedTransmittanceOptions {
                    max_free_flight_draws: remaining_draws,
                    max_tracking_candidates: remaining_candidates,
                },
                |draw| (self.samples)(request.ordinal, draw),
            )
            .map_err(|error| match error {
                CurvedTransmittanceError::DrawLimit => PointLightCameraError::SourceDrawLimit,
                CurvedTransmittanceError::CandidateLimit => {
                    PointLightCameraError::SourceCandidateLimit
                }
                error => PointLightCameraError::SourceConnection {
                    binding: request.binding,
                    error,
                },
            })?;
        Ok(FixtureConnection {
            original: *result.connection(),
            work: ConnectionWork {
                traces: result.trace_evaluations(),
                fields: result.field_evaluations(),
                geometry: result.query_work(),
                source_draws: result.free_flight_draws(),
                source_candidates: result.tracking_candidates(),
                ..ConnectionWork::default()
            },
            log_scattering_survival: result.log_scattering_survival(),
            ratio_tracking_zero: result.ratio_tracking_zero(),
        })
    }
}

fn camera_settings(settings: RayleighPointLightCameraSettings) -> CurvedRayleighTraceSettings {
    CurvedRayleighTraceSettings {
        surfaces: settings.point_lights.surfaces,
        max_volume_events: settings.max_volume_events,
    }
}

#[cfg(test)]
mod tests;
