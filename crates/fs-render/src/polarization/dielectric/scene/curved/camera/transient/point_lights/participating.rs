//! Retarded point fixtures at native molecular surface and volume vertices.
//! Original collision/survival laws and source Stokes remain unchanged; only
//! each exact source's intensity is observed on its own qualified group clock.

use super::super::super::point_lights::participating::{
    CurvedRayleighPointLightCameraScene, RayleighPointLightCameraSample,
    RayleighPointLightCameraSettings,
};
use super::*;

pub mod images;
mod sequential;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientRayleighPointLightCameraSample {
    observation_time_s: Dd,
    point: RayleighPointLightCameraSample,
}
impl TransientRayleighPointLightCameraSample {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn point_light_sample(&self) -> &RayleighPointLightCameraSample {
        &self.point
    }
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        self.point.radiance_sample()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureRayleighPointLightCameraSample {
    window: ExposureWindow,
    point: RayleighPointLightCameraSample,
    integrated_radiance_j_m2_sr_nm: f64,
}
impl ExposureRayleighPointLightCameraSample {
    pub const fn window(&self) -> ExposureWindow {
        self.window
    }
    pub const fn mean_point_light_sample(&self) -> &RayleighPointLightCameraSample {
        &self.point
    }
    pub const fn mean_radiance_sample(&self) -> &CameraRadianceSample {
        self.point.radiance_sample()
    }
    pub const fn integrated_radiance_j_m2_sr_nm(&self) -> f64 {
        self.integrated_radiance_j_m2_sr_nm
    }
}

/// The exact original camera/fixture owner plus borrowed immutable histories.
/// Unbound fixtures remain steady and require no source-connection group data.
pub struct TransientCurvedRayleighPointLightCameraScene<'a, 'c, 'h, 's, 'm, 'o, 'l> {
    observer: PointHistoryObserver<'a, 'h, 's, 'm, 'l>,
    lights: &'a CurvedRayleighPointLightCameraScene<'c, 'h, 's, 'm, 'o, 'l>,
}

impl<'a, 'c, 'h, 's, 'm, 'o, 'l>
    TransientCurvedRayleighPointLightCameraScene<'a, 'c, 'h, 's, 'm, 'o, 'l>
{
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
        lights: &'a CurvedRayleighPointLightCameraScene<'c, 'h, 's, 'm, 'o, 'l>,
        histories: &'a [PointEmissionHistoryBinding<'l>],
    ) -> Result<Self, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if !core::ptr::eq(camera.source_camera(), lights.surface_camera()) {
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
    ) -> &'a CurvedRayleighPointLightCameraScene<'c, 'h, 's, 'm, 'o, 'l> {
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

fn finish_molecular_exposure(
    window: ExposureWindow,
    point: RayleighPointLightCameraSample,
) -> Result<ExposureRayleighPointLightCameraSample, TransientPointLightCameraError> {
    let integrated_radiance_j_m2_sr_nm = integrate_rate(
        point.radiance_sample().radiance_w_m2_sr_nm(),
        window.duration_s(),
    )
    .map_err(TransientCameraError::from)?;
    Ok(ExposureRayleighPointLightCameraSample {
        window,
        point,
        integrated_radiance_j_m2_sr_nm,
    })
}
