//! Retarded finite-area emission on the original clear camera/MIS estimator.
//!
//! Every native connection image evaluates its original source at the camera
//! clock minus the sum of the camera-prefix and source-connection **group**
//! delays. The sum is formed in DD before division by c. Camera emitter hits
//! use their own accumulated native delay and the original complementary MIS
//! weight. Source positions, image discovery, Jacobians, branch probabilities,
//! source-area samples and camera continuations remain in their native owners.
//!
//! Geometry and optics are frozen. Finite shutters request each source's mean
//! over its individually retarded interval; they do not resample transport or
//! replace a source integral with a time sample. Unbound sources remain steady.
//! This adapter does not include molecular scattering along either subpath.

use super::super::area_lights::camera::{
    AreaLightCameraError, AreaLightCameraSample, AreaLightCameraSettings, AreaLightCameraWorkspace,
    AreaSourcePath, CurvedAreaLightCameraScene,
};
use super::super::{CameraRadianceError, CameraRadianceSample};
use super::exposure::{ExposureWindow, average_history, integrate_rate};
use super::{
    SurfaceEmissionQuery, TransientCameraError, TransientCurvedCameraScene, admit_time,
    retarded_time, sample_history,
};
use crate::charts::Ray;
use crate::polarization::dielectric::transport::PolarizedBeam;
use fs_exec::Cx;
use fs_math::dd::Dd;

mod sequential;

/// Original molecular area-source transport observed at retarded source times.
#[cfg(feature = "volumes")]
pub mod participating;

#[derive(Clone, Debug, PartialEq)]
pub enum TransientAreaLightCameraError {
    /// Both prepared owners must borrow the same original camera instance.
    CameraIdentityMismatch,
    Area(AreaLightCameraError),
    Source(TransientCameraError),
}

impl From<AreaLightCameraError> for TransientAreaLightCameraError {
    fn from(error: AreaLightCameraError) -> Self {
        Self::Area(error)
    }
}
impl From<TransientCameraError> for TransientAreaLightCameraError {
    fn from(error: TransientCameraError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for TransientAreaLightCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient native area-light camera: {self:?}")
    }
}
impl core::error::Error for TransientAreaLightCameraError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientAreaLightCameraSample {
    observation_time_s: Dd,
    area: AreaLightCameraSample,
}

impl TransientAreaLightCameraSample {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    /// The complete original sample, including per-image source/MIS work.
    pub const fn area_sample(&self) -> &AreaLightCameraSample {
        &self.area
    }
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        self.area.radiance_sample()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureAreaLightCameraSample {
    window: ExposureWindow,
    area: AreaLightCameraSample,
    integrated_radiance_j_m2_sr_nm: f64,
}

impl ExposureAreaLightCameraSample {
    pub const fn window(&self) -> ExposureWindow {
        self.window
    }
    /// Mean source radiance on the one original camera and connection walk.
    pub const fn mean_area_sample(&self) -> &AreaLightCameraSample {
        &self.area
    }
    pub const fn mean_radiance_sample(&self) -> &CameraRadianceSample {
        self.area.radiance_sample()
    }
    pub const fn integrated_radiance_j_m2_sr_nm(&self) -> f64 {
        self.integrated_radiance_j_m2_sr_nm
    }
}

/// Borrow prepared original sources and their histories. Preparation and each
/// observation add no allocation; callers reuse the native area workspace.
pub struct TransientCurvedAreaLightCameraScene<'a, 'h, 's, 'm, 'l> {
    camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
    lights: &'a CurvedAreaLightCameraScene<'h, 's, 'm, 'l>,
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

impl<'a, 'h, 's, 'm, 'l> TransientCurvedAreaLightCameraScene<'a, 'h, 's, 'm, 'l> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a TransientCurvedCameraScene<'h, 's, 'm>,
        lights: &'a CurvedAreaLightCameraScene<'h, 's, 'm, 'l>,
    ) -> Result<Self, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        if !core::ptr::eq(camera.source_camera(), lights.source_camera()) {
            return Err(TransientAreaLightCameraError::CameraIdentityMismatch);
        }
        Ok(Self { camera, lights })
    }

    pub const fn source_camera(&self) -> &'a TransientCurvedCameraScene<'h, 's, 'm> {
        self.camera
    }
    pub const fn source_lights(&self) -> &'a CurvedAreaLightCameraScene<'h, 's, 'm, 'l> {
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
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<TransientAreaLightCameraSample, TransientAreaLightCameraError> {
        cx.checkpoint().map_err(TransientCameraError::from)?;
        let observation_time_s = admit_time(observation_time_s)?;
        let area = self.observe(
            cx,
            ray,
            wavelength_nm,
            initial_medium,
            Observation::Instant(observation_time_s),
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
    pub fn sample_ray_with_samples_during_exposure(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<ExposureAreaLightCameraSample, TransientAreaLightCameraError> {
        let area = self.observe(
            cx,
            ray,
            wavelength_nm,
            initial_medium,
            Observation::Exposure(window),
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

    #[allow(clippy::too_many_arguments)]
    fn observe(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation: Observation,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<AreaLightCameraSample, TransientAreaLightCameraError> {
        self.observe_launched(
            cx,
            || self.camera.launch(ray, wavelength_nm, initial_medium),
            initial_medium,
            observation,
            Dd::ZERO,
            settings,
            workspace,
            surface_samples,
            area_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn observe_launched(
        &self,
        cx: &Cx<'_>,
        launch: impl FnOnce() -> Result<PolarizedBeam, CameraRadianceError>,
        initial_medium: usize,
        observation: Observation,
        launch_delay_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<AreaLightCameraSample, TransientAreaLightCameraError> {
        let mut source_error = None;
        let result = self.lights.sample_with_modulation(
            cx,
            launch,
            initial_medium,
            settings,
            workspace,
            surface_samples,
            area_samples,
            |path| {
                self.scale(cx, observation, launch_delay_s, path).map_err(|error| {
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

    fn scale(
        &self,
        cx: &Cx<'_>,
        observation: Observation,
        launch_delay_s: Dd,
        source: AreaSourcePath<'_>,
    ) -> Result<f64, TransientCameraError> {
        let (object_id, position_m, wavelength_nm, camera_group_m, connection_group_m) =
            match source {
                AreaSourcePath::Surface(path) => {
                    return match observation {
                        Observation::Instant(time) => {
                            self.camera.source_scale_with_delay(cx, time, launch_delay_s, path)
                        }
                        Observation::Exposure(window) => {
                            self.camera.exposure_scale_with_delay(cx, window, launch_delay_s, path)
                        }
                    };
                }
                AreaSourcePath::Connected {
                    object_id,
                    position_m,
                    wavelength_nm,
                    camera_group_m,
                    connection_group_m,
                } => (
                    object_id,
                    position_m,
                    wavelength_nm,
                    camera_group_m,
                    connection_group_m,
                ),
            };
        let slot = self
            .camera
            .camera
            .scene
            .bindings
            .binary_search_by_key(&object_id, |binding| binding.object_id)
            .map_err(|_| CameraRadianceError::NumericalRange)?;
        let Some(binding) = self.camera.histories[slot] else {
            return Ok(1.0);
        };
        cx.checkpoint()?;
        let group_delay_s = super::connected_group_delay(
            camera_group_m, connection_group_m, launch_delay_s,
        ).map_err(|error| match error {
            super::ConnectionClockError::MissingGroupDelay => TransientCameraError::MissingGroupDelay { object_id },
            super::ConnectionClockError::UnresolvedClock => TransientCameraError::UnresolvedClock { object_id },
        })?;
        let observation_time_s = observation.open_s();
        let query = SurfaceEmissionQuery {
            object_id,
            position_m,
            wavelength_nm,
            observation_time_s,
            emission_time_s: retarded_time(observation_time_s, group_delay_s, object_id)?,
            group_delay_s,
        };
        match observation {
            Observation::Instant(_) => sample_history(cx, binding, query),
            Observation::Exposure(window) => average_history(cx, binding, query, window),
        }
    }
}
