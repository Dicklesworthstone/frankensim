//! Source-owned shutter integrals on original real-lens molecular paths.

use super::super::super::exposure::ExposureWindow;
use super::super::exposure::ExposureRayleighCameraSample;
use super::*;

impl TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_> {
    /// Integrate one physical sensor shutter on a timed lens beam with the
    /// original seeded molecular transport. Every source receives its own
    /// lens-plus-scene retarded interval; volume events remain instantaneous.
    /// The beam owns its sensor analyzer and a separately configured analyzer
    /// refuses. Frozen media, absorption and LTE emission remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let (beam, delay) = self.launch_timed_sequential(optical, initial_medium)?;
        self.sample_launched_ray_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            delay,
            settings,
            seed,
            path_id,
        )
    }

    /// Uniform-density flights retain the original three variates and their
    /// admission. Source interval averages consume no transport variates.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let (beam, delay) = self.launch_timed_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_samples_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            delay,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    /// Original uniform or heterogeneous tracking with independent acceptance.
    /// Null candidates retain the original state/counts and query no history.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let (beam, delay) = self.launch_timed_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_tracking_samples_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            delay,
            settings,
            surface_samples,
            flight_samples,
        )
    }

    // The raster already placed the configured analyzer at the sensor while
    // generating this beam. Apply its full transported state only once.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super::super) fn sample_prepared_sequential_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<ExposureRayleighCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let delay = lens_delay(optical)?;
        let beam = self
            .camera
            .source_camera()
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched_ray_with_tracking_samples_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            delay,
            settings,
            surface_samples,
            flight_samples,
        )
    }
}
