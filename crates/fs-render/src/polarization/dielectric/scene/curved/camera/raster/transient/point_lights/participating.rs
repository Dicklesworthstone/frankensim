//! CIE/QE point histories on the original molecular camera and source streams.

use super::*;
use super::super::super::point_lights::participating::RayleighPointLightPixelSample;
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::RayleighPointLightCameraSettings;
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::participating::TransientCurvedRayleighPointLightCameraScene;

mod images;
mod sequential;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientRayleighPointLightPixelSample<P = CurvedPixelSample> {
    observation_time_s: Dd,
    point: RayleighPointLightPixelSample<P>,
}
impl<P> TransientRayleighPointLightPixelSample<P> {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn pixel_sample(&self) -> &P {
        self.point.pixel_sample()
    }
    pub const fn molecular_sample(&self) -> &RayleighPointLightPixelSample<P> {
        &self.point
    }
}

impl CurvedCameraRaster<'_> {
    /// Observe each exact source at its own native group delay.
    /// Original spectral/pupil integration and source work remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_point_lights_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientRayleighPointLightPixelSample<CurvedPixelSample>,
        TransientPointLightRasterError,
    > {
        admit_analyzer_band(scene.source_camera(), LAMBDA_MIN, LAMBDA_MAX)?;
        let (observation_time_s, point) = self.sample_timed_rayleigh_point_lights_transport(
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
        Ok(TransientRayleighPointLightPixelSample {
            observation_time_s,
            point,
        })
    }

    /// Integrate every original source over its own retarded shutter.
    /// Original spectral/pupil integration and source work remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_point_lights_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighPointLightPixelSample<ExposurePixelSample<CurvedPixelSample>>,
        TransientPointLightRasterError,
    > {
        admit_analyzer_band(scene.source_camera(), LAMBDA_MIN, LAMBDA_MAX)?;
        let (_, point) = self.sample_timed_rayleigh_point_lights_transport(
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
    pub fn sample_pixel_sensor_rayleigh_point_lights_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientRayleighPointLightPixelSample<CurvedSensorPixelSample<C>>,
        TransientPointLightRasterError,
    > {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        admit_analyzer_band(scene.source_camera(), lo, hi)?;
        let (observation_time_s, point) = self.sample_timed_rayleigh_point_lights_transport(
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
        Ok(TransientRayleighPointLightPixelSample {
            observation_time_s,
            point,
        })
    }

    /// Integrate every original source over its own retarded shutter.
    /// Original spectral/pupil integration and source work remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_rayleigh_point_lights_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        RayleighPointLightPixelSample<ExposureSensorPixelSample<C, CurvedSensorPixelSample<C>>>,
        TransientPointLightRasterError,
    > {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        admit_analyzer_band(scene.source_camera(), lo, hi)?;
        let (_, point) = self.sample_timed_rayleigh_point_lights_transport(
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
    fn sample_timed_rayleigh_point_lights_transport<P>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation: Observation,
        settings: RayleighPointLightCameraSettings,
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
    ) -> Result<(Dd, RayleighPointLightPixelSample<P>), TransientPointLightRasterError> {
        // Admission precedes optional lens lanes, including completely clipped packets.
        let (clock, observation) = admit_observation(cx, observation)?;
        let point = self.sample_rayleigh_point_with_launch(
            seed,
            pixel,
            sample,
            |cx, (ray, wavelength): (Ray, f64), draws| match observation {
                Observation::Instant(time) => scene
                    .sample_ray_with_tracking_samples_at_time(
                        cx,
                        ray,
                        wavelength,
                        initial_medium,
                        time,
                        settings,
                        |event| draws.surface(event),
                        |flight| draws.flight(flight),
                        |connection, flight| draws.source(connection, flight),
                    )
                    .map(|result| *result.point_light_sample()),
                Observation::Exposure(window) => scene
                    .sample_ray_with_tracking_samples_during_exposure(
                        cx,
                        ray,
                        wavelength,
                        initial_medium,
                        window,
                        settings,
                        |event| draws.surface(event),
                        |flight| draws.flight(flight),
                        |connection, flight| draws.source(connection, flight),
                    )
                    .map(|result| *result.mean_point_light_sample()),
            },
            integrate,
        )?;
        Ok((clock, point))
    }
}
