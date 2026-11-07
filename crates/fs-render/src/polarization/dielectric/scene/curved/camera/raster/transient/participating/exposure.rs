//! Finite-shutter color and QE observations through native molecular media.
//!
//! Four wavelength paths are each traced once. Their source-owned means use
//! each actual refracted/scattered group delay before the existing pupil and
//! detector integration. All channels share those same paths and means.

use super::super::exposure::{
    ExposurePixelSample, ExposureSensorPixelSample, finish_pixel, finish_sensor,
};
use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureRayleighPixelSample<P = CurvedPixelSample> {
    sample: ExposurePixelSample<P>,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
}

impl<P> ExposureRayleighPixelSample<P> {
    pub const fn exposure_sample(&self) -> &ExposurePixelSample<P> {
        &self.sample
    }
    pub const fn window(&self) -> ExposureWindow {
        self.sample.window()
    }
    pub const fn mean_pixel_sample(&self) -> &P {
        self.sample.mean_pixel_sample()
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }

    pub(super) fn with_work(
        sample: ExposurePixelSample<P>,
        volume_events: [u32; CAMERA_WAVELENGTHS],
        free_flight_draws: [u32; CAMERA_WAVELENGTHS],
    ) -> Self {
        Self {
            sample,
            volume_events,
            free_flight_draws,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureRayleighSensorPixelSample<
    const CHANNELS: usize,
    P = CurvedSensorPixelSample<CHANNELS>,
> {
    sample: ExposureSensorPixelSample<CHANNELS, P>,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
}

impl<const C: usize, P> ExposureRayleighSensorPixelSample<C, P> {
    pub const fn exposure_sample(&self) -> &ExposureSensorPixelSample<C, P> {
        &self.sample
    }
    pub const fn window(&self) -> ExposureWindow {
        self.sample.window()
    }
    pub const fn mean_pixel_sample(&self) -> &P {
        self.sample.mean_pixel_sample()
    }
    pub const fn qe_weighted_energy_j_m2_sr(&self) -> &[f64; C] {
        self.sample.qe_weighted_energy_j_m2_sr()
    }
    pub const fn photoelectrons_per_m2_sr(&self) -> &[f64; C] {
        self.sample.photoelectrons_per_m2_sr()
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }

    pub(super) fn with_work(
        sample: ExposureSensorPixelSample<C, P>,
        volume_events: [u32; CAMERA_WAVELENGTHS],
        free_flight_draws: [u32; CAMERA_WAVELENGTHS],
    ) -> Self {
        Self {
            sample,
            volume_events,
            free_flight_draws,
        }
    }
}

impl CurvedCameraRaster<'_> {
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<ExposureRayleighPixelSample, TransientRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_error = None;
        let result = self.sample_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let result = scene.sample_ray_with_tracking_samples_during_exposure(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    window,
                    settings,
                    |event| {
                        let block =
                            philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |flight| {
                        philox4x32_10([pixel, sample, flight, 0x6372_7350 + lane as u32], key)
                            .map(unit)
                    },
                );
                match result {
                    Ok(exposure) => {
                        volume_events[lane] = exposure.volume_events();
                        free_flight_draws[lane] = exposure.free_flight_draws();
                        Ok(*exposure.mean_radiance_sample())
                    }
                    Err(error) => {
                        source_error = Some(error);
                        Err(CameraRadianceError::NumericalRange)
                    }
                }
            },
        );
        if let Some(error) = source_error {
            return Err(TransientRasterError::Source(error));
        }
        Ok(ExposureRayleighPixelSample::with_work(
            finish_pixel(window, result?),
            volume_events,
            free_flight_draws,
        ))
    }

    /// Original source-bounded QE means and integrated energy/electron counts,
    /// before pixel etendue and sensor electronics. Failed interval histories,
    /// unresolved paths and cancellation refuse the complete observation.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_rayleigh_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<ExposureRayleighSensorPixelSample<C>, TransientRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_error = None;
        let result = self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let result = scene.sample_ray_with_tracking_samples_during_exposure(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    window,
                    settings,
                    |event| {
                        let block =
                            philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |flight| {
                        philox4x32_10([pixel, sample, flight, 0x6372_7350 + lane as u32], key)
                            .map(unit)
                    },
                );
                match result {
                    Ok(exposure) => {
                        volume_events[lane] = exposure.volume_events();
                        free_flight_draws[lane] = exposure.free_flight_draws();
                        Ok(*exposure.mean_radiance_sample())
                    }
                    Err(error) => {
                        source_error = Some(error);
                        Err(CameraRadianceError::NumericalRange)
                    }
                }
            },
        );
        if let Some(error) = source_error {
            return Err(TransientRasterError::Source(error));
        }
        Ok(ExposureRayleighSensorPixelSample::with_work(
            finish_sensor(window, result?)?,
            volume_events,
            free_flight_draws,
        ))
    }
}
