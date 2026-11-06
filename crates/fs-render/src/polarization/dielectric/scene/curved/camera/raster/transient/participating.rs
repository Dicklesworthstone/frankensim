//! Time-resolved CIE pixels and QE rates on actual refracted/scattered paths.
//!
//! Each channel shares the same four original spectral paths. The native surface
//! and free-flight counter domains, heterogeneous acceptance variate and optical
//! pupil weights are preserved; history sampling adds no transport resampling.
//! These are instantaneous observations at the launch clock under the parent
//! transient camera's source and frozen-medium assumptions.

use super::super::CAMERA_WAVELENGTHS;
use super::{
    CameraRadianceError, CurvedCameraRaster, CurvedPixelSample, CurvedRasterError,
    CurvedSensorPixelSample, TransientPixelSample, TransientRasterError,
    TransientSensorPixelSample, philox4x32_10, unit,
};
use crate::polarization::dielectric::scene::curved::camera::transient::participating::TransientCurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighTraceSettings;
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_exec::Cx;
use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientRayleighPixelSample {
    sample: TransientPixelSample,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
}

impl TransientRayleighPixelSample {
    pub const fn transient_sample(&self) -> &TransientPixelSample {
        &self.sample
    }
    pub const fn observation_time_s(&self) -> Dd {
        self.sample.observation_time_s()
    }
    pub const fn pixel_sample(&self) -> &CurvedPixelSample {
        self.sample.pixel_sample()
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientRayleighSensorPixelSample<const CHANNELS: usize> {
    sample: TransientSensorPixelSample<CHANNELS>,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
}

impl<const C: usize> TransientRayleighSensorPixelSample<C> {
    pub const fn transient_sample(&self) -> &TransientSensorPixelSample<C> {
        &self.sample
    }
    pub const fn observation_time_s(&self) -> Dd {
        self.sample.observation_time_s()
    }
    pub const fn pixel_sample(&self) -> &CurvedSensorPixelSample<C> {
        self.sample.pixel_sample()
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }
}

impl CurvedCameraRaster<'_> {
    /// Native CIE integration with source histories evaluated at the accumulated
    /// group delay of each actual molecularly scattered wavelength path.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_transient_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<TransientRayleighPixelSample, TransientRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_error = None;
        let mut admitted_time = None;
        let result = self.sample_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let result = scene.sample_ray_with_tracking_samples_at_time(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    observation_time_s,
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
                    Ok(timed) => {
                        admitted_time = Some(timed.observation_time_s());
                        volume_events[lane] = timed.volume_events();
                        free_flight_draws[lane] = timed.free_flight_draws();
                        Ok(*timed.radiance_sample())
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
        let pixel = result?;
        Ok(TransientRayleighPixelSample {
            sample: TransientPixelSample {
                observation_time_s: admitted_time.ok_or(CurvedRasterError::NumericalRange)?,
                pixel,
            },
            volume_events,
            free_flight_draws,
        })
    }

    /// Source-bounded QE and photon integration; all channels share four paths,
    /// their history queries and their original collision work. A failed lane
    /// refuses the complete observation without renormalizing remaining lanes.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_rayleigh_transient_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<TransientRayleighSensorPixelSample<C>, TransientRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_error = None;
        let mut admitted_time = None;
        let result = self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let result = scene.sample_ray_with_tracking_samples_at_time(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    observation_time_s,
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
                    Ok(timed) => {
                        admitted_time = Some(timed.observation_time_s());
                        volume_events[lane] = timed.volume_events();
                        free_flight_draws[lane] = timed.free_flight_draws();
                        Ok(*timed.radiance_sample())
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
        let pixel = result?;
        Ok(TransientRayleighSensorPixelSample {
            sample: TransientSensorPixelSample {
                observation_time_s: admitted_time.ok_or(CurvedRasterError::NumericalRange)?,
                pixel,
            },
            volume_events,
            free_flight_draws,
        })
    }
}

#[cfg(test)]
mod tests;
