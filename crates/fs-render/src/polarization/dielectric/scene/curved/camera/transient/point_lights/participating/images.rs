//! Each discovered molecular fixture image retains its own source clock.
//! Duplicate seeds do not duplicate histories or source-survival streams.

use super::super::super::sequential::lens_delay;
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::images::CurvedRayleighImageLightCameraScene;
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;

/// The exact original camera/fixture owner plus borrowed immutable histories.
/// Unbound fixtures remain steady and require no source-connection group data.
pub struct TransientCurvedRayleighImageLightCameraScene<'a, 'g, 'c, 'h, 's, 'm, 'o, 'l> {
    observer: PointHistoryObserver<'a, 'h, 's, 'm, 'l>,
    lights: &'a CurvedRayleighImageLightCameraScene<'g, 'c, 'h, 's, 'm, 'o, 'l>,
}

impl<'a, 'g, 'c, 'h, 's, 'm, 'o, 'l>
    TransientCurvedRayleighImageLightCameraScene<'a, 'g, 'c, 'h, 's, 'm, 'o, 'l>
{
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
        lights: &'a CurvedRayleighImageLightCameraScene<'g, 'c, 'h, 's, 'm, 'o, 'l>,
        histories: &'a [PointEmissionHistoryBinding<'l>],
    ) -> Result<Self, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if !core::ptr::eq(
            camera.source_camera(),
            lights.source_camera().surface_camera(),
        ) {
            return Err(PointLightCameraError::Camera(CameraRadianceError::Scene(
                crate::polarization::dielectric::scene::SceneError::InvalidBindings.into(),
            ))
            .into());
        }
        let observer = PointHistoryObserver::try_new(cx, camera, histories, |source| {
            lights.contains_source(source)
        })?;
        Ok(Self { observer, lights })
    }
    pub const fn source_camera(&self) -> &'a TransientCurvedCameraScene<'h, 's, 'm> {
        self.observer.camera
    }
    pub const fn source_lights(
        &self,
    ) -> &'a CurvedRayleighImageLightCameraScene<'g, 'c, 'h, 's, 'm, 'o, 'l> {
        self.lights
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<TransientRayleighPointLightCameraSample, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
        let point = self.observe_launched(
            cx,
            || {
                self.observer
                    .camera
                    .launch(ray, wavelength_nm, initial_medium)
            },
            initial_medium,
            Observation::Instant(observation_time_s),
            Dd::ZERO,
            settings,
            workspace,
            surface_samples,
            flight_samples,
            source_samples,
        )?;
        Ok(TransientRayleighPointLightCameraSample {
            observation_time_s,
            point,
        })
    }

    /// Each source integrates its own retarded shutter; no temporal draws,
    /// retraces or midpoint substitutions alter the native point estimator.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<ExposureRayleighPointLightCameraSample, TransientPointLightCameraError> {
        let point = self.observe_launched(
            cx,
            || {
                self.observer
                    .camera
                    .launch(ray, wavelength_nm, initial_medium)
            },
            initial_medium,
            Observation::Exposure(window),
            Dd::ZERO,
            settings,
            workspace,
            surface_samples,
            flight_samples,
            source_samples,
        )?;
        finish_molecular_exposure(window, point)
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_launched(
        &self,
        cx: &Cx<'_>,
        launch: impl FnOnce() -> Result<PolarizedBeam, CameraRadianceError>,
        initial_medium: usize,
        observation: Observation,
        launch_delay_s: Dd,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<RayleighPointLightCameraSample, TransientPointLightCameraError> {
        let mut source_error = None;
        let result = self.lights.sample_with_modulation(
            cx,
            launch,
            initial_medium,
            settings,
            workspace,
            surface_samples,
            flight_samples,
            source_samples,
            |path| {
                self.observer
                    .scale(cx, observation, launch_delay_s, path)
                    .map_err(|error| {
                        source_error = Some(error);
                        PointLightCameraError::NumericalRange
                    })
            },
        );
        if let Some(error) = source_error {
            return Err(error);
        }
        Ok(result?)
    }
}

impl TransientCurvedRayleighImageLightCameraScene<'_, '_, '_, '_, '_, '_, '_, '_> {
    /// A complete beam already contains any sensor analyzer. A separately
    /// configured row must be installed by the physical raster before the lens.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<TransientRayleighPointLightCameraSample, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if self.source_camera().analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens.into());
        }
        self.sample_prepared_sequential_ray_with_tracking_samples_at_time(
            cx,
            optical,
            initial_medium,
            observation_time_s,
            settings,
            workspace,
            surface_samples,
            flight_samples,
            source_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<ExposureRayleighPointLightCameraSample, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if self.source_camera().analyzer_response().is_some() {
            return Err(TransientCameraError::AnalyzerMustPrecedeLens.into());
        }
        self.sample_prepared_sequential_ray_with_tracking_samples_during_exposure(
            cx,
            optical,
            initial_medium,
            window,
            settings,
            workspace,
            surface_samples,
            flight_samples,
            source_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super::super) fn sample_prepared_sequential_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<TransientRayleighPointLightCameraSample, TransientPointLightCameraError> {
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
            workspace,
            surface_samples,
            flight_samples,
            source_samples,
        )?;
        Ok(TransientRayleighPointLightCameraSample {
            observation_time_s,
            point,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super::super) fn sample_prepared_sequential_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<ExposureRayleighPointLightCameraSample, TransientPointLightCameraError> {
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
            workspace,
            surface_samples,
            flight_samples,
            source_samples,
        )?;
        finish_molecular_exposure(window, point)
    }
}
