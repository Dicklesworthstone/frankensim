//! Finite shutters on the original refracted and molecularly scattered paths.
//!
//! Frozen transport is linear in source radiance. Each actual surface arrival
//! requests one source-owned mean over its group-delay-shifted shutter; no
//! temporal quadrature, repeated paths or extra collision variates are added.
//! Volume emission and every native absorption/Mueller/tracking weight remain
//! unchanged. The base exposure module owns interval admission and units.

use super::super::exposure::{ExposureCameraSample, ExposureWindow, finish as finish_mean};
use super::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureRayleighCameraSample {
    sample: ExposureCameraSample,
    volume_events: u32,
    free_flight_draws: u32,
}

impl ExposureRayleighCameraSample {
    pub const fn exposure_sample(&self) -> &ExposureCameraSample {
        &self.sample
    }
    pub const fn window(&self) -> ExposureWindow {
        self.sample.window()
    }
    pub const fn mean_radiance_sample(&self) -> &CameraRadianceSample {
        self.sample.mean_radiance_sample()
    }
    pub const fn integrated_radiance_j_m2_sr_nm(&self) -> f64 {
        self.sample.integrated_radiance_j_m2_sr_nm()
    }
    pub const fn volume_events(&self) -> u32 {
        self.volume_events
    }
    pub const fn free_flight_draws(&self) -> u32 {
        self.free_flight_draws
    }
}

impl TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_> {
    /// Integrate one shutter with the original seeded surface/flight streams.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            Dd::ZERO,
            settings,
            seed,
            path_id,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn sample_launched_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        window: ExposureWindow,
        launch_delay_s: Dd,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
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
            |event| {
                self.observe_exposure(cx, window, launch_delay_s, event, &mut emission, &mut error)
            },
        );
        self.finish_exposure(cx, window, index, emission, result, error)
    }

    /// Uniform-density three-variate flights retain their original admission.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_samples_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            Dd::ZERO,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn sample_launched_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        window: ExposureWindow,
        launch_delay_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
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
                |event| {
                    self.observe_exposure(
                        cx,
                        window,
                        launch_delay_s,
                        event,
                        &mut emission,
                        &mut error,
                    )
                },
            );
        self.finish_exposure(cx, window, index, emission, result, error)
    }

    /// Original uniform or heterogeneous tracking; source integration never
    /// consumes a distance, angular or independent acceptance variate.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_tracking_samples_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            Dd::ZERO,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn sample_launched_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        window: ExposureWindow,
        launch_delay_s: Dd,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
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
                |event| {
                    self.observe_exposure(
                        cx,
                        window,
                        launch_delay_s,
                        event,
                        &mut emission,
                        &mut error,
                    )
                },
            );
        self.finish_exposure(cx, window, index, emission, result, error)
    }

    fn observe_exposure(
        &self,
        cx: &Cx<'_>,
        window: ExposureWindow,
        launch_delay_s: Dd,
        event: CurvedRayleighEvent<'_>,
        emission: &mut SurfaceEmission<'_, '_, '_>,
        error: &mut Option<TransientCameraError>,
    ) -> Result<(), CurvedSceneError> {
        match event {
            CurvedRayleighEvent::Surface(event) => self.camera.observe_exposure_with_delay(
                cx,
                window,
                launch_delay_s,
                event,
                emission,
                error,
            ),
            event => emission.observe_rayleigh_event(event),
        }
    }

    fn finish_exposure(
        &self,
        cx: &Cx<'_>,
        window: ExposureWindow,
        initial_index: f64,
        emission: SurfaceEmission<'_, '_, '_>,
        result: Result<CurvedRayleighPath, CurvedSceneError>,
        error: Option<TransientCameraError>,
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
        if let Some(error) = error {
            return Err(error);
        }
        let (path, emitted) = emission.finish_scattering(result)?;
        let mean = self
            .camera
            .source_camera()
            .finish(cx, *path.path(), initial_index, emitted)?;
        Ok(ExposureRayleighCameraSample {
            sample: finish_mean(window, mean)?,
            volume_events: path.volume_events(),
            free_flight_draws: path.free_flight_draws(),
        })
    }
}
