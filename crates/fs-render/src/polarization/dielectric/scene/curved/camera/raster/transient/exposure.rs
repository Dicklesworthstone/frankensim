//! Native color and QE observations integrated over a finite shutter.
//!
//! Each of the original four spectral paths is traced once. Its individual
//! retarded source intervals are averaged before the unchanged pupil and
//! spectral weights. The detector then converts the mean rates to integrated
//! energy and electrons per m² sr using the explicit shutter duration. These
//! are aperture averages before pixel etendue, gain, saturation and noise.

use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, integrate_rate,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposurePixelSample {
    window: ExposureWindow,
    mean_pixel: CurvedPixelSample,
}

impl ExposurePixelSample {
    pub const fn window(&self) -> ExposureWindow {
        self.window
    }
    /// Original CIE/XYZ observation of the shutter-averaged radiance.
    pub const fn mean_pixel_sample(&self) -> &CurvedPixelSample {
        &self.mean_pixel
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureSensorPixelSample<const CHANNELS: usize> {
    window: ExposureWindow,
    mean_pixel: CurvedSensorPixelSample<CHANNELS>,
    qe_weighted_energy_j_m2_sr: [f64; CHANNELS],
    photoelectrons_per_m2_sr: [f64; CHANNELS],
}

impl<const C: usize> ExposureSensorPixelSample<C> {
    pub const fn window(&self) -> ExposureWindow {
        self.window
    }
    pub const fn mean_pixel_sample(&self) -> &CurvedSensorPixelSample<C> {
        &self.mean_pixel
    }
    /// Time-integrated QE-weighted radiance, before pixel etendue.
    pub const fn qe_weighted_energy_j_m2_sr(&self) -> &[f64; C] {
        &self.qe_weighted_energy_j_m2_sr
    }
    /// Time-integrated expected electron rate, before pixel etendue/electronics.
    pub const fn photoelectrons_per_m2_sr(&self) -> &[f64; C] {
        &self.photoelectrons_per_m2_sr
    }
}

impl CurvedCameraRaster<'_> {
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedCameraScene<'_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<ExposurePixelSample, TransientRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut source_error = None;
        let result = self.sample_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let observed = scene.sample_ray_with_samples_during_exposure(
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
                );
                match observed {
                    Ok(exposure) => Ok(*exposure.mean_radiance_sample()),
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
        Ok(finish_pixel(window, result?))
    }

    /// All channels share the original four paths, source integrals and pupil
    /// draws. Unsupported interval histories or any failed lane refuse the
    /// complete observation, including already accumulated detector channels.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedCameraScene<'_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<ExposureSensorPixelSample<C>, TransientRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut source_error = None;
        let result = self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let observed = scene.sample_ray_with_samples_during_exposure(
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
                );
                match observed {
                    Ok(exposure) => Ok(*exposure.mean_radiance_sample()),
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
        finish_sensor(window, result?)
    }
}

pub(super) fn finish_pixel(
    window: ExposureWindow,
    mean_pixel: CurvedPixelSample,
) -> ExposurePixelSample {
    ExposurePixelSample { window, mean_pixel }
}

pub(super) fn finish_sensor<const C: usize>(
    window: ExposureWindow,
    mean_pixel: CurvedSensorPixelSample<C>,
) -> Result<ExposureSensorPixelSample<C>, TransientRasterError> {
    let mut qe_weighted_energy_j_m2_sr = [0.0; C];
    let mut photoelectrons_per_m2_sr = [0.0; C];
    for channel in 0..C {
        qe_weighted_energy_j_m2_sr[channel] = integrate_rate(
            mean_pixel.qe_weighted_radiance_w_m2_sr()[channel],
            window.duration_s(),
        )
        .map_err(|error| TransientRasterError::Source(error.into()))?;
        photoelectrons_per_m2_sr[channel] = integrate_rate(
            mean_pixel.photoelectron_rate_per_m2_sr()[channel],
            window.duration_s(),
        )
        .map_err(|error| TransientRasterError::Source(error.into()))?;
    }
    Ok(ExposureSensorPixelSample {
        window,
        mean_pixel,
        qe_weighted_energy_j_m2_sr,
        photoelectrons_per_m2_sr,
    })
}
