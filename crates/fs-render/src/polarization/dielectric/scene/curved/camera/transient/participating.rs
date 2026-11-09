//! Finite-source histories on the native refractive molecular camera walk.
//!
//! Each surface history reads the group delay accumulated along the actual
//! curved, multiply scattered path before that surface response. Free-flight
//! survival, density acceptance, Mueller/phase weights, absorption and LTE
//! emission remain owned by the original participating-medium transport.
//! Histories consume no variates and add no collision or extinction weight.
//!
//! The underlying transient camera's frozen-medium, incoherent-envelope and
//! finite-surface-source assumptions apply. Molecular events are instantaneous;
//! no fluorescent lifetime or time-dependent volume source is inferred.

use super::{
    CameraRadianceError, CameraRadianceSample, CurvedSceneError, SurfaceEmission,
    TransientCameraError, TransientCameraSample, TransientCurvedCameraScene, admit_time,
    phase_index,
};
use crate::charts::Ray;
use crate::polarization::dielectric::scene::SceneError;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighEvent, CurvedRayleighPath, CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use crate::polarization::dielectric::transport::PolarizedBeam;
use fs_exec::Cx;
use fs_math::dd::Dd;

pub mod exposure;
mod sequential;

/// A completed retarded source observation and its original volume-work counts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientRayleighCameraSample {
    sample: TransientCameraSample,
    volume_events: u32,
    free_flight_draws: u32,
}

impl TransientRayleighCameraSample {
    pub const fn transient_sample(&self) -> &TransientCameraSample {
        &self.sample
    }
    pub const fn observation_time_s(&self) -> Dd {
        self.sample.observation_time_s()
    }
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        self.sample.radiance_sample()
    }
    pub const fn volume_events(&self) -> u32 {
        self.volume_events
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
}

/// Allocation-free composition of admitted histories and native molecular media.
/// Both must borrow the exact same original surface scene.
pub struct TransientCurvedRayleighCameraScene<'a, 'h, 's, 'm, 'o> {
    camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
    scattering: &'a CurvedRayleighScene<'s, 'm, 'o>,
}

impl<'a, 'h, 's, 'm, 'o> TransientCurvedRayleighCameraScene<'a, 'h, 's, 'm, 'o> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
        scattering: &'a CurvedRayleighScene<'s, 'm, 'o>,
    ) -> Result<Self, TransientCameraError> {
        cx.checkpoint()?;
        if !core::ptr::eq(
            camera.source_camera().source_scene(),
            scattering.surface_scene(),
        ) {
            return Err(CameraRadianceError::Scene(SceneError::InvalidBindings.into()).into());
        }
        Ok(Self { camera, scattering })
    }

    pub const fn transient_camera(&self) -> &'a TransientCurvedCameraScene<'h, 's, 'm> {
        self.camera
    }
    pub const fn participating_scene(&self) -> &'a CurvedRayleighScene<'s, 'm, 'o> {
        self.scattering
    }

    /// Preserve the original independent surface and free-flight Philox streams.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_at_time(
            cx,
            beam,
            initial_medium,
            time,
            Dd::ZERO,
            settings,
            seed,
            path_id,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_ray_at_time(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        time: Dd,
        launch_delay_s: Dd,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        let camera = self.camera.source_camera();
        let index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(camera, index, beam);
        let mut error = None;
        let result = self.scattering.for_camera().trace_events_observed(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            &camera.spatial_volume_emitters,
            |event| self.observe(cx, time, launch_delay_s, event, &mut emission, &mut error),
        );
        self.finish(cx, time, index, emission, result, error)
    }

    /// Uniform density retains [distance CDF, phase cosine, phase azimuth].
    /// Positive heterogeneous majorants require the independent acceptance
    /// variate provided by `sample_ray_with_tracking_samples_at_time`.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_samples_at_time(
            cx,
            beam,
            initial_medium,
            time,
            Dd::ZERO,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        time: Dd,
        launch_delay_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        let camera = self.camera.source_camera();
        let index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(camera, index, beam);
        let mut error = None;
        let result = self
            .scattering
            .for_camera()
            .trace_with_samples_events_observed(
                cx,
                beam,
                initial_medium,
                settings,
                surface_samples,
                flight_samples,
                &camera.spatial_volume_emitters,
                |event| self.observe(cx, time, launch_delay_s, event, &mut emission, &mut error),
            );
        self.finish(cx, time, index, emission, result, error)
    }

    /// Uniform or heterogeneous density retains the native fourth delta-tracking
    /// acceptance variate. Null candidates neither query histories nor reset the
    /// accumulated source delay, polarization or thermal-emission state.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let time = admit_time(observation_time_s)?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_tracking_samples_at_time(
            cx,
            beam,
            initial_medium,
            time,
            Dd::ZERO,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        time: Dd,
        launch_delay_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        let camera = self.camera.source_camera();
        let index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(camera, index, beam);
        let mut error = None;
        let result = self
            .scattering
            .for_camera()
            .trace_with_tracking_samples_events_observed(
                cx,
                beam,
                initial_medium,
                settings,
                surface_samples,
                flight_samples,
                &camera.spatial_volume_emitters,
                |event| self.observe(cx, time, launch_delay_s, event, &mut emission, &mut error),
            );
        self.finish(cx, time, index, emission, result, error)
    }

    fn observe(
        &self,
        cx: &Cx<'_>,
        time: Dd,
        launch_delay_s: Dd,
        event: CurvedRayleighEvent<'_>,
        emission: &mut SurfaceEmission<'_, '_, '_>,
        error: &mut Option<TransientCameraError>,
    ) -> Result<(), CurvedSceneError> {
        match event {
            CurvedRayleighEvent::Surface(event) => {
                self.camera
                    .observe_with_delay(cx, time, launch_delay_s, event, emission, error)
            }
            event => emission.observe_rayleigh_event(event),
        }
    }

    fn finish(
        &self,
        cx: &Cx<'_>,
        time: Dd,
        initial_index: f64,
        emission: SurfaceEmission<'_, '_, '_>,
        result: Result<CurvedRayleighPath, CurvedSceneError>,
        error: Option<TransientCameraError>,
    ) -> Result<TransientRayleighCameraSample, TransientCameraError> {
        if let Some(error) = error {
            return Err(error);
        }
        let (path, emitted) = emission.finish_scattering(result)?;
        let radiance =
            self.camera
                .source_camera()
                .finish(cx, *path.path(), initial_index, emitted)?;
        Ok(TransientRayleighCameraSample {
            sample: TransientCameraSample {
                observation_time_s: time,
                radiance,
            },
            volume_events: path.volume_events(),
            free_flight_draws: path.free_flight_draws(),
        })
    }
}

#[cfg(test)]
mod tests;
