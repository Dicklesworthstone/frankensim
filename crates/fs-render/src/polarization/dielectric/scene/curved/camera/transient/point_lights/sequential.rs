//! Sensor-clock point fixtures through an original resolved full-Stokes lens.

use super::super::sequential::lens_delay;
use super::*;
use crate::camera::SequentialCameraBeam;

impl TransientCurvedPointLightCameraScene<'_, '_, '_, '_, '_> {
    /// A complete beam already contains any sensor analyzer. A separately
    /// configured row must be installed by the physical raster before the lens.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientPointLightCameraSample, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if self.source_camera().analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens.into());
        }
        self.sample_prepared_sequential_ray_with_samples_at_time(
            cx,
            optical,
            initial_medium,
            observation_time_s,
            settings,
            surface_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposurePointLightCameraSample, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if self.source_camera().analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens.into());
        }
        self.sample_prepared_sequential_ray_with_samples_during_exposure(
            cx,
            optical,
            initial_medium,
            window,
            settings,
            surface_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientPointLightCameraSample, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
        let delay = lens_delay(optical)?;
        let point = self.observe_launched(
            cx,
            || {
                self.source_camera()
                    .source_camera()
                    .launch_sequential(optical, initial_medium)
            },
            initial_medium,
            Observation::Instant(observation_time_s),
            delay,
            settings,
            surface_samples,
        )?;
        Ok(TransientPointLightCameraSample {
            observation_time_s,
            point,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposurePointLightCameraSample, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let delay = lens_delay(optical)?;
        let point = self.observe_launched(
            cx,
            || {
                self.source_camera()
                    .source_camera()
                    .launch_sequential(optical, initial_medium)
            },
            initial_medium,
            Observation::Exposure(window),
            delay,
            settings,
            surface_samples,
        )?;
        finish_exposure(window, point)
    }
}
