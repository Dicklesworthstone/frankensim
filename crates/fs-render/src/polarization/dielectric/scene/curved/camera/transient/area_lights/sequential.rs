//! Finite-area NEE/MIS continued from a timed, fully polarized real-lens beam.

use super::super::sequential::lens_delay;
use super::*;
use crate::camera::SequentialCameraBeam;

impl TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_> {
    /// Continue the original timed lens response through finite-area next-event
    /// estimation and complementary camera emitter hits. Each source image
    /// receives the sensor clock and its own lens-plus-scene group delay.
    /// Lens Fresnel/Beer weights and endpoint normalization stay in the beam;
    /// the caller applies its separate geometric sampling weight once.
    ///
    /// The beam already owns any analyzer installed at the physical sensor.
    /// A separately configured transient-scene analyzer is refused here; the
    /// sequential raster installs that analyzer before the lens instead.
    /// Untimed beams are refused, including when no history is bound.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<TransientAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let time = admit_time(observation_time_s)?;
        if self.camera.analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens.into());
        }
        self.sample_prepared_sequential_ray_with_samples_at_time(
            cx,
            optical,
            initial_medium,
            time,
            settings,
            workspace,
            surface_samples,
            area_samples,
        )
    }

    /// Integrate each original source history over its individually retarded
    /// sensor shutter on one frozen lens/camera/connection walk. No temporal
    /// samples, extra transport, or additional random draws are introduced.
    /// Analyzer and qualified-delay admission match the instantaneous method.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<ExposureAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if self.camera.analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens.into());
        }
        self.sample_prepared_sequential_ray_with_samples_during_exposure(
            cx,
            optical,
            initial_medium,
            window,
            settings,
            workspace,
            surface_samples,
            area_samples,
        )
    }

    // Only the original raster may use these prepared seams: it installs the
    // configured analyzer at the sensor before generating the timed beam.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<TransientAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
        let delay = lens_delay(optical)?;
        let area = self.observe_launched(
            cx,
            || {
                self.camera
                    .source_camera()
                    .launch_sequential(optical, initial_medium)
            },
            initial_medium,
            Observation::Instant(observation_time_s),
            delay,
            settings,
            workspace,
            surface_samples,
            area_samples,
        )?;
        Ok(TransientAreaLightCameraSample {
            observation_time_s,
            area,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<ExposureAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let delay = lens_delay(optical)?;
        let area = self.observe_launched(
            cx,
            || {
                self.camera
                    .source_camera()
                    .launch_sequential(optical, initial_medium)
            },
            initial_medium,
            Observation::Exposure(window),
            delay,
            settings,
            workspace,
            surface_samples,
            area_samples,
        )?;
        let integrated_radiance_j_m2_sr_nm = integrate_rate(
            area.radiance_sample().radiance_w_m2_sr_nm(),
            window.duration_s(),
        )
        .map_err(TransientCameraError::from)?;
        Ok(ExposureAreaLightCameraSample {
            window,
            area,
            integrated_radiance_j_m2_sr_nm,
        })
    }
}
