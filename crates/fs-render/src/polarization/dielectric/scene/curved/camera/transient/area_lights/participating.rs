//! Retarded finite-area emission on the original molecular NEE/MIS estimator.
//!
//! Camera collisions, source scattering survival, native refracted images and
//! polarization remain in their original owners. A finite surface history only
//! scales that source's radiance at its own camera-plus-connection group delay.
//! Spatial volume/LTE emission remains steady and molecular events instantaneous.

use super::super::sequential::lens_delay;
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::camera::area_lights::participating::{
    CurvedRayleighAreaLightCameraScene, RayleighAreaLightCameraSample,
    RayleighAreaLightCameraSettings, RayleighAreaLightCameraWorkspace,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientRayleighAreaLightCameraSample {
    observation_time_s: Dd,
    area: RayleighAreaLightCameraSample,
}

impl TransientRayleighAreaLightCameraSample {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    /// Original molecular camera, source-survival and native image work.
    pub const fn area_sample(&self) -> &RayleighAreaLightCameraSample {
        &self.area
    }
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        self.area.radiance_sample()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureRayleighAreaLightCameraSample {
    window: ExposureWindow,
    area: RayleighAreaLightCameraSample,
    integrated_radiance_j_m2_sr_nm: f64,
}

impl ExposureRayleighAreaLightCameraSample {
    pub const fn window(&self) -> ExposureWindow {
        self.window
    }
    pub const fn mean_area_sample(&self) -> &RayleighAreaLightCameraSample {
        &self.area
    }
    pub const fn mean_radiance_sample(&self) -> &CameraRadianceSample {
        self.area.radiance_sample()
    }
    pub const fn integrated_radiance_j_m2_sr_nm(&self) -> f64 {
        self.integrated_radiance_j_m2_sr_nm
    }
}

/// Borrow the original admitted histories and molecular area-source owner.
/// The clear observer supplies the SAME source-clock law; it performs no clear
/// transport here. Both native source owners must share the original camera.
pub struct TransientCurvedRayleighAreaLightCameraScene<'a, 'c, 'h, 's, 'm, 'o, 'l> {
    observer: TransientCurvedAreaLightCameraScene<'a, 'h, 's, 'm, 'l>,
    lights: &'a CurvedRayleighAreaLightCameraScene<'c, 'h, 's, 'm, 'o, 'l>,
}

impl<'a, 'c, 'h, 's, 'm, 'o, 'l>
    TransientCurvedRayleighAreaLightCameraScene<'a, 'c, 'h, 's, 'm, 'o, 'l>
{
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
        lights: &'a CurvedRayleighAreaLightCameraScene<'c, 'h, 's, 'm, 'o, 'l>,
    ) -> Result<Self, TransientAreaLightCameraError> {
        let observer =
            TransientCurvedAreaLightCameraScene::try_new(cx, camera, lights.area_lights())?;
        Ok(Self { observer, lights })
    }

    pub const fn source_camera(&self) -> &'a TransientCurvedCameraScene<'h, 's, 'm> {
        self.observer.source_camera()
    }
    pub const fn source_lights(
        &self,
    ) -> &'a CurvedRayleighAreaLightCameraScene<'c, 'h, 's, 'm, 'o, 'l> {
        self.lights
    }

    /// Preserve the native surface, camera-flight, surface-area, molecular-area
    /// and source-survival proposal streams. Histories consume no variates.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        volume_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        source_samples: impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<TransientRayleighAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
        let area = self.observe_launched(
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
            surface_area_samples,
            volume_area_samples,
            source_samples,
        )?;
        Ok(TransientRayleighAreaLightCameraSample {
            observation_time_s,
            area,
        })
    }

    /// Integrate every source's retarded shutter interval on one frozen
    /// camera/collision/image realization, without time samples or retracing.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        volume_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        source_samples: impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<ExposureRayleighAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let area = self.observe_launched(
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
            surface_area_samples,
            volume_area_samples,
            source_samples,
        )?;
        finish_exposure(window, area)
    }

    /// The complete timed lens beam retains Fresnel/Beer and Stokes response.
    /// A configured analyzer must already precede the lens at the sensor;
    /// raw-beam calls refuse a separate observer analyzer. The geometric lens
    /// sampling weight remains separate, to be applied once by the raster.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        volume_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        source_samples: impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<TransientRayleighAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
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
            surface_area_samples,
            volume_area_samples,
            source_samples,
        )
    }

    /// Finite exposure with the same qualified lens prefix and original
    /// molecular/source-survival realization as the instantaneous query.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        volume_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        source_samples: impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<ExposureRayleighAreaLightCameraSample, TransientAreaLightCameraError> {
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
            surface_area_samples,
            volume_area_samples,
            source_samples,
        )
    }

    // The native raster alone installs the observer analyzer at the physical
    // sensor before generating these timed full-Stokes beams.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_with_tracking_samples_at_time(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        volume_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        source_samples: impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<TransientRayleighAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
        let delay = lens_delay(optical)?;
        let area = self.observe_launched(
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
            surface_area_samples,
            volume_area_samples,
            source_samples,
        )?;
        Ok(TransientRayleighAreaLightCameraSample {
            observation_time_s,
            area,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(in super::super::super) fn sample_prepared_sequential_ray_with_tracking_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        volume_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        source_samples: impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<ExposureRayleighAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let delay = lens_delay(optical)?;
        let area = self.observe_launched(
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
            surface_area_samples,
            volume_area_samples,
            source_samples,
        )?;
        finish_exposure(window, area)
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_launched(
        &self,
        cx: &Cx<'_>,
        launch: impl FnOnce() -> Result<PolarizedBeam, CameraRadianceError>,
        initial_medium: usize,
        observation: Observation,
        launch_delay_s: Dd,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        surface_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        volume_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        source_samples: impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<RayleighAreaLightCameraSample, TransientAreaLightCameraError> {
        let mut source_error = None;
        let result = self.lights.sample_with_modulation(
            cx,
            launch,
            initial_medium,
            settings,
            workspace,
            surface_samples,
            flight_samples,
            surface_area_samples,
            volume_area_samples,
            source_samples,
            |path| {
                self.observer
                    .scale(cx, observation, launch_delay_s, path)
                    .map_err(|error| {
                        source_error = Some(error);
                        AreaLightCameraError::NumericalRange
                    })
            },
        );
        if let Some(error) = source_error {
            return Err(error.into());
        }
        Ok(result?)
    }
}

fn finish_exposure(
    window: ExposureWindow,
    area: RayleighAreaLightCameraSample,
) -> Result<ExposureRayleighAreaLightCameraSample, TransientAreaLightCameraError> {
    let integrated_radiance_j_m2_sr_nm = integrate_rate(
        area.radiance_sample().radiance_w_m2_sr_nm(),
        window.duration_s(),
    )
    .map_err(TransientCameraError::from)?;
    Ok(ExposureRayleighAreaLightCameraSample {
        window,
        area,
        integrated_radiance_j_m2_sr_nm,
    })
}
