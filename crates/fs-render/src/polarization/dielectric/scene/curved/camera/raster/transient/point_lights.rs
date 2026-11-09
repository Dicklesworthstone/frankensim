//! Original CIE/QE point-fixture observations at retarded source clocks.
//! Native source branches/images, counter domains and work are shared with
//! steady rendering. Failed lanes refuse the complete observation; clipped
//! real-lens lanes retain their original probability mass and zero work.

use super::super::point_lights::{PointLightPixelSample, PointLightRasterError};
use super::super::sensor::{CurvedSensorPixelSample, CurvedSequentialSensorPixelSample};
use super::super::{
    CurvedCameraRaster, CurvedPixelSample, CurvedRasterError, CurvedSequentialPixelSample,
};
use super::exposure::{
    ExposurePixelSample, ExposureSensorPixelSample, finish_pixel, finish_sensor_rates,
};
use super::{LAMBDA_MAX, LAMBDA_MIN, TransientRasterError, admit_analyzer_band};
use crate::camera::SequentialCameraBeam;
use crate::charts::Ray;
use crate::polarization::dielectric::scene::curved::camera::point_lights::PointLightCameraSettings;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::{
    TransientCurvedPointLightCameraScene, TransientPointLightCameraError,
};
use crate::polarization::dielectric::scene::curved::camera::transient::{
    TransientCameraError, admit_time,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CameraRadianceError, CameraRadianceSample,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_exec::Cx;
use fs_math::dd::Dd;

#[cfg(feature = "volumes")]
pub mod participating;
mod sequential;

#[derive(Debug)]
pub enum TransientPointLightRasterError {
    Raster(CurvedRasterError),
    Source(TransientCameraError),
    Transport {
        lane: usize,
        error: TransientPointLightCameraError,
    },
}
impl From<PointLightRasterError<TransientPointLightCameraError>>
    for TransientPointLightRasterError
{
    fn from(error: PointLightRasterError<TransientPointLightCameraError>) -> Self {
        match error {
            PointLightRasterError::Raster(error) => Self::Raster(error),
            PointLightRasterError::Transport { lane, error } => Self::Transport { lane, error },
        }
    }
}
impl From<TransientRasterError> for TransientPointLightRasterError {
    fn from(error: TransientRasterError) -> Self {
        match error {
            TransientRasterError::Raster(error) => Self::Raster(error),
            TransientRasterError::Source(error) => Self::Source(error),
        }
    }
}
impl core::fmt::Display for TransientPointLightRasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient native point-fixture pixel: {self:?}")
    }
}
impl core::error::Error for TransientPointLightRasterError {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientPointLightPixelSample<P = CurvedPixelSample> {
    observation_time_s: Dd,
    point: PointLightPixelSample<P>,
}
impl<P> TransientPointLightPixelSample<P> {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn pixel_sample(&self) -> &P {
        self.point.pixel_sample()
    }
    pub const fn point_light_sample(&self) -> &PointLightPixelSample<P> {
        &self.point
    }
}

#[derive(Clone, Copy)]
enum Observation {
    Instant(Dd),
    Exposure(ExposureWindow),
}
fn admit_observation(
    cx: &Cx<'_>,
    observation: Observation,
) -> Result<(Dd, Observation), TransientPointLightRasterError> {
    cx.checkpoint()
        .map_err(|error| TransientPointLightRasterError::Source(error.into()))?;
    match observation {
        Observation::Instant(time) => {
            let time = admit_time(time).map_err(TransientPointLightRasterError::Source)?;
            Ok((time, Observation::Instant(time)))
        }
        Observation::Exposure(window) => Ok((window.open_s(), Observation::Exposure(window))),
    }
}

impl CurvedCameraRaster<'_> {
    /// Observe each exact source at its own native group delay.
    /// Original spectral/pupil integration and source work remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_point_lights_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<TransientPointLightPixelSample<CurvedPixelSample>, TransientPointLightRasterError>
    {
        admit_analyzer_band(scene.source_camera(), LAMBDA_MIN, LAMBDA_MAX)?;
        let (observation_time_s, point) = self.sample_timed_point_lights_transport(
            cx,
            scene,
            initial_medium,
            Observation::Instant(observation_time_s),
            settings,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_pixel_with_transport(
                    cx,
                    seed,
                    pixel,
                    sample,
                    |cx, ray, wavelength, lane| trace(cx, (ray, wavelength), lane),
                )
            },
        )?;
        Ok(TransientPointLightPixelSample {
            observation_time_s,
            point,
        })
    }

    /// Integrate every original source over its own retarded shutter.
    /// Original spectral/pupil integration and source work remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_point_lights_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        PointLightPixelSample<ExposurePixelSample<CurvedPixelSample>>,
        TransientPointLightRasterError,
    > {
        admit_analyzer_band(scene.source_camera(), LAMBDA_MIN, LAMBDA_MAX)?;
        let (_, point) = self.sample_timed_point_lights_transport(
            cx,
            scene,
            initial_medium,
            Observation::Exposure(window),
            settings,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_pixel_with_transport(
                    cx,
                    seed,
                    pixel,
                    sample,
                    |cx, ray, wavelength, lane| trace(cx, (ray, wavelength), lane),
                )
            },
        )?;
        Ok(
            point
                .try_map_pixel(|mean| Ok::<_, TransientRasterError>(finish_pixel(window, mean)))?,
        )
    }

    /// Observe each exact source at its own native group delay.
    /// Original spectral/pupil integration and source work remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_point_lights_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedPointLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientPointLightPixelSample<CurvedSensorPixelSample<C>>,
        TransientPointLightRasterError,
    > {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        admit_analyzer_band(scene.source_camera(), lo, hi)?;
        let (observation_time_s, point) = self.sample_timed_point_lights_transport(
            cx,
            scene,
            initial_medium,
            Observation::Instant(observation_time_s),
            settings,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_sensor_with_transport(
                    cx,
                    sensor,
                    seed,
                    pixel,
                    sample,
                    |cx, ray, wavelength, lane| trace(cx, (ray, wavelength), lane),
                )
            },
        )?;
        Ok(TransientPointLightPixelSample {
            observation_time_s,
            point,
        })
    }

    /// Integrate every original source over its own retarded shutter.
    /// Original spectral/pupil integration and source work remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_point_lights_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedPointLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        PointLightPixelSample<ExposureSensorPixelSample<C, CurvedSensorPixelSample<C>>>,
        TransientPointLightRasterError,
    > {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        admit_analyzer_band(scene.source_camera(), lo, hi)?;
        let (_, point) = self.sample_timed_point_lights_transport(
            cx,
            scene,
            initial_medium,
            Observation::Exposure(window),
            settings,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_sensor_with_transport(
                    cx,
                    sensor,
                    seed,
                    pixel,
                    sample,
                    |cx, ray, wavelength, lane| trace(cx, (ray, wavelength), lane),
                )
            },
        )?;
        Ok(point.try_map_pixel(|mean| {
            let energy = *mean.qe_weighted_radiance_w_m2_sr();
            let electrons = *mean.photoelectron_rate_per_m2_sr();
            finish_sensor_rates(window, mean, energy, electrons)
        })?)
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_timed_point_lights_transport<P>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation: Observation,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        integrate: impl FnOnce(
            &mut dyn FnMut(
                &Cx<'_>,
                (Ray, f64),
                usize,
            ) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<(Dd, PointLightPixelSample<P>), TransientPointLightRasterError> {
        // Admission precedes optional lens lanes, including completely clipped packets.
        let (clock, observation) = admit_observation(cx, observation)?;
        let point = self.sample_point_with_launch(
            seed,
            pixel,
            sample,
            |cx, (ray, wavelength): (Ray, f64), draws| match observation {
                Observation::Instant(time) => scene
                    .sample_ray_with_samples_at_time(
                        cx,
                        ray,
                        wavelength,
                        initial_medium,
                        time,
                        settings,
                        |event| draws.surface(event),
                    )
                    .map(|result| *result.point_light_sample()),
                Observation::Exposure(window) => scene
                    .sample_ray_with_samples_during_exposure(
                        cx,
                        ray,
                        wavelength,
                        initial_medium,
                        window,
                        settings,
                        |event| draws.surface(event),
                    )
                    .map(|result| *result.mean_point_light_sample()),
            },
            integrate,
        )?;
        Ok((clock, point))
    }
}
