//! Original CIE/QE observations of retarded finite-area source connections.
//!
//! Each wavelength retains its own refracted camera path, source images and
//! source-clock intervals. Native area sampling, continuation counter domains,
//! MIS, pupil weights and work accounting are shared with the steady adapter.
//! All detector channels share those same four paths. One failed source/image
//! or spectral lane refuses the complete observation without renormalization.

use super::super::area_lights::{AreaLightPixelSample, AreaLightRasterError};
use super::super::sensor::CurvedSensorPixelSample;
use super::exposure::{
    ExposurePixelSample, ExposureSensorPixelSample, finish_pixel, finish_sensor,
};
use super::{CurvedCameraRaster, CurvedPixelSample, CurvedRasterError, TransientRasterError};
use crate::charts::Ray;
use crate::polarization::dielectric::scene::curved::camera::area_lights::camera::{
    AreaLightCameraSettings, AreaLightCameraWorkspace,
};
use crate::polarization::dielectric::scene::curved::camera::transient::TransientCameraError;
use crate::polarization::dielectric::scene::curved::camera::transient::area_lights::{
    TransientAreaLightCameraError, TransientCurvedAreaLightCameraScene,
};
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;
use crate::polarization::dielectric::scene::curved::camera::{
    CameraRadianceError, CameraRadianceSample,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_exec::Cx;
use fs_math::dd::Dd;

mod sequential;

/// Molecular finite-area CIE/QE observations through a timed real lens.
#[cfg(feature = "volumes")]
pub mod participating;

#[derive(Debug)]
pub enum TransientAreaLightRasterError {
    Raster(CurvedRasterError),
    /// A final detector conversion refused the complete observation.
    Source(TransientCameraError),
    Transport {
        lane: usize,
        error: TransientAreaLightCameraError,
    },
}

impl From<AreaLightRasterError<TransientAreaLightCameraError>> for TransientAreaLightRasterError {
    fn from(error: AreaLightRasterError<TransientAreaLightCameraError>) -> Self {
        match error {
            AreaLightRasterError::Raster(error) => Self::Raster(error),
            AreaLightRasterError::Transport { lane, error } => Self::Transport { lane, error },
        }
    }
}
impl From<TransientRasterError> for TransientAreaLightRasterError {
    fn from(error: TransientRasterError) -> Self {
        match error {
            TransientRasterError::Raster(error) => Self::Raster(error),
            TransientRasterError::Source(error) => Self::Source(error),
        }
    }
}
impl core::fmt::Display for TransientAreaLightRasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient native area-light pixel: {self:?}")
    }
}
impl core::error::Error for TransientAreaLightRasterError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientAreaLightPixelSample<P = CurvedPixelSample> {
    observation_time_s: Dd,
    area: AreaLightPixelSample<P>,
}

impl<P> TransientAreaLightPixelSample<P> {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn pixel_sample(&self) -> &P {
        self.area.pixel_sample()
    }
    /// Original per-lane area radiance, image and native source-work accounts.
    pub const fn area_sample(&self) -> &AreaLightPixelSample<P> {
        &self.area
    }
}

#[derive(Clone, Copy)]
enum Observation {
    Instant(Dd),
    Exposure(ExposureWindow),
}

impl CurvedCameraRaster<'_> {
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_area_lights_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<TransientAreaLightPixelSample, TransientAreaLightRasterError> {
        super::admit_analyzer_band(scene.source_camera(), super::LAMBDA_MIN, super::LAMBDA_MAX)?;
        let (observation_time_s, area) = self.sample_timed_area_transport(
            scene,
            initial_medium,
            Observation::Instant(observation_time_s),
            settings,
            workspace,
            seed,
            pixel,
            sample,
            |trace| self.sample_pixel_with_transport(cx, seed, pixel, sample, trace),
        )?;
        Ok(TransientAreaLightPixelSample {
            observation_time_s,
            area,
        })
    }

    /// All QE channels share native wavelength-specific source images and
    /// retarded histories; no CIE/display conversion or extra transport occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_area_lights_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientAreaLightPixelSample<CurvedSensorPixelSample<C>>,
        TransientAreaLightRasterError,
    > {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        super::admit_analyzer_band(scene.source_camera(), lo, hi)?;
        let (observation_time_s, area) = self.sample_timed_area_transport(
            scene,
            initial_medium,
            Observation::Instant(observation_time_s),
            settings,
            workspace,
            seed,
            pixel,
            sample,
            |trace| self.sample_sensor_with_transport(cx, sensor, seed, pixel, sample, trace),
        )?;
        Ok(TransientAreaLightPixelSample {
            observation_time_s,
            area,
        })
    }

    /// One frozen native camera/source trace per lane supplies every source's
    /// shutter mean before the original pupil and CIE integration.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_area_lights_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample<ExposurePixelSample>, TransientAreaLightRasterError> {
        super::admit_analyzer_band(scene.source_camera(), super::LAMBDA_MIN, super::LAMBDA_MAX)?;
        let (_, area) = self.sample_timed_area_transport(
            scene,
            initial_medium,
            Observation::Exposure(window),
            settings,
            workspace,
            seed,
            pixel,
            sample,
            |trace| self.sample_pixel_with_transport(cx, seed, pixel, sample, trace),
        )?;
        Ok(area.try_map_pixel(|mean| Ok::<_, TransientRasterError>(finish_pixel(window, mean)))?)
    }

    /// Original QE mean rates and integrated energy/electrons per m² sr,
    /// before pixel etendue, gain, saturation or noise. Work is channel-independent.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_area_lights_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample<ExposureSensorPixelSample<C>>, TransientAreaLightRasterError>
    {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        super::admit_analyzer_band(scene.source_camera(), lo, hi)?;
        let (_, area) = self.sample_timed_area_transport(
            scene,
            initial_medium,
            Observation::Exposure(window),
            settings,
            workspace,
            seed,
            pixel,
            sample,
            |trace| self.sample_sensor_with_transport(cx, sensor, seed, pixel, sample, trace),
        )?;
        Ok(area.try_map_pixel(|mean| finish_sensor(window, mean))?)
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_timed_area_transport<P>(
        &self,
        scene: &TransientCurvedAreaLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation: Observation,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
        integrate: impl FnOnce(
            &mut dyn FnMut(
                &Cx<'_>,
                Ray,
                f64,
                usize,
            ) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<(Dd, AreaLightPixelSample<P>), TransientAreaLightRasterError> {
        let mut admitted_time = None;
        let area = self.sample_area_with_transport(
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, samples| {
                let (clock, traced) = match observation {
                    Observation::Instant(time) => {
                        let result = scene.sample_ray_with_samples_at_time(
                            cx,
                            ray,
                            wavelength,
                            initial_medium,
                            time,
                            settings,
                            workspace,
                            |event| samples.surface(event),
                            |event, binding| samples.area(event, binding),
                        )?;
                        (result.observation_time_s(), *result.area_sample())
                    }
                    Observation::Exposure(window) => {
                        let result = scene.sample_ray_with_samples_during_exposure(
                            cx,
                            ray,
                            wavelength,
                            initial_medium,
                            window,
                            settings,
                            workspace,
                            |event| samples.surface(event),
                            |event, binding| samples.area(event, binding),
                        )?;
                        (window.open_s(), *result.mean_area_sample())
                    }
                };
                admitted_time = Some(clock);
                Ok::<_, TransientAreaLightCameraError>(traced)
            },
            integrate,
        )?;
        let clock = admitted_time.ok_or(TransientAreaLightRasterError::Raster(
            CurvedRasterError::NumericalRange,
        ))?;
        Ok((clock, area))
    }
}
