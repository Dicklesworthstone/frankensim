//! Retarded intensity histories on original refractive point-fixture paths.
//!
//! Each exact calibrated source keeps its spectrum, polarization and native
//! connection. Histories see its actual outgoing direction and camera-plus-
//! connection group delay; a resolved lens contributes its sensor prefix once.
//! Point fixtures are delta endpoints and retain their original additive NEE
//! without an invented competing endpoint PDF or MIS factor. Geometry and
//! optics are frozen, and spatial volume/LTE emission remains steady.

use super::super::point_lights::{
    CurvedPointLightCameraScene, PointLightCameraError, PointLightCameraSample,
    PointLightCameraSettings, PointSourcePath,
};
use super::super::{CameraRadianceError, CameraRadianceSample};
use super::exposure::{ExposureWindow, integrate_rate};
use super::{SurfaceHistoryError, TransientCameraError, TransientCurvedCameraScene, admit_time};
use crate::charts::Ray;
use crate::polarization::dielectric::transport::PolarizedBeam;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_math::dd::Dd;

mod history;
mod sequential;
use history::PointHistoryObserver;
pub use history::{
    PointEmissionHistory, PointEmissionHistoryBinding, PointEmissionQuery, PointExposureQuery,
};

#[derive(Clone, Debug, PartialEq)]
pub enum TransientPointLightCameraError {
    Point(PointLightCameraError),
    /// Original finite-surface history, analyzer, lens or observation refusal.
    Source(TransientCameraError),
    InvalidHistoryBinding {
        index: usize,
    },
    MissingHistoryIdentity {
        index: usize,
    },
    MissingGroupDelay {
        source_identity: ContentHash,
    },
    UnresolvedClock {
        source_identity: ContentHash,
    },
    InvalidHistoryScale {
        source_identity: ContentHash,
    },
    History {
        source_identity: ContentHash,
        error: SurfaceHistoryError,
    },
}
impl From<PointLightCameraError> for TransientPointLightCameraError {
    fn from(error: PointLightCameraError) -> Self {
        Self::Point(error)
    }
}
impl From<TransientCameraError> for TransientPointLightCameraError {
    fn from(error: TransientCameraError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for TransientPointLightCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient native point-fixture camera: {self:?}")
    }
}
impl core::error::Error for TransientPointLightCameraError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientPointLightCameraSample {
    observation_time_s: Dd,
    point: PointLightCameraSample,
}
impl TransientPointLightCameraSample {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn point_light_sample(&self) -> &PointLightCameraSample {
        &self.point
    }
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        self.point.radiance_sample()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposurePointLightCameraSample {
    window: ExposureWindow,
    point: PointLightCameraSample,
    integrated_radiance_j_m2_sr_nm: f64,
}
impl ExposurePointLightCameraSample {
    pub const fn window(&self) -> ExposureWindow {
        self.window
    }
    pub const fn mean_point_light_sample(&self) -> &PointLightCameraSample {
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
pub struct TransientCurvedPointLightCameraScene<'a, 'h, 's, 'm, 'l> {
    observer: PointHistoryObserver<'a, 'h, 's, 'm, 'l>,
    lights: &'a CurvedPointLightCameraScene<'h, 's, 'm, 'l>,
}

#[derive(Clone, Copy)]
enum Observation {
    Instant(Dd),
    Exposure(ExposureWindow),
}
impl Observation {
    fn open_s(self) -> Dd {
        match self {
            Self::Instant(time) => time,
            Self::Exposure(window) => window.open_s(),
        }
    }
}

impl<'a, 'h, 's, 'm, 'l> TransientCurvedPointLightCameraScene<'a, 'h, 's, 'm, 'l> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
        lights: &'a CurvedPointLightCameraScene<'h, 's, 'm, 'l>,
        histories: &'a [PointEmissionHistoryBinding<'l>],
    ) -> Result<Self, TransientPointLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if !core::ptr::eq(camera.source_camera(), lights.source_camera()) {
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
    pub const fn source_lights(&self) -> &'a CurvedPointLightCameraScene<'h, 's, 'm, 'l> {
        self.lights
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientPointLightCameraSample, TransientPointLightCameraError> {
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
        )?;
        Ok(TransientPointLightCameraSample {
            observation_time_s,
            point,
        })
    }

    /// Each source integrates its own retarded shutter; no temporal draws,
    /// retraces or midpoint substitutions alter the native point estimator.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<ExposurePointLightCameraSample, TransientPointLightCameraError> {
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
        )?;
        finish_exposure(window, point)
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_launched(
        &self,
        cx: &Cx<'_>,
        launch: impl FnOnce() -> Result<PolarizedBeam, CameraRadianceError>,
        initial_medium: usize,
        observation: Observation,
        launch_delay_s: Dd,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<PointLightCameraSample, TransientPointLightCameraError> {
        let mut source_error = None;
        let result = self.lights.sample_with_modulation(
            cx,
            launch,
            initial_medium,
            settings,
            surface_samples,
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

fn finish_exposure(
    window: ExposureWindow,
    point: PointLightCameraSample,
) -> Result<ExposurePointLightCameraSample, TransientPointLightCameraError> {
    let integrated_radiance_j_m2_sr_nm = integrate_rate(
        point.radiance_sample().radiance_w_m2_sr_nm(),
        window.duration_s(),
    )
    .map_err(TransientCameraError::from)?;
    Ok(ExposurePointLightCameraSample {
        window,
        point,
        integrated_radiance_j_m2_sr_nm,
    })
}
