//! Spectral pixels and physical QE channels at one native observation clock.
//!
//! Each wavelength independently follows the original curved camera path and
//! reads finite-source histories at that path's own accumulated group delay.
//! Existing raster, pupil, wavelength and surface random streams are unchanged.
//! Every detector channel shares the same four original spectral paths; history
//! and transport work do not scale with the number of QE channels.
//!
//! The DD clock is at the ray's original launch/pupil. Results are instantaneous
//! radiance and electron-rate observations under the native aperture-average
//! convention, before exposure/etendue/electronics. No camera-lens transit time,
//! exposure integration, temporal interpolation or extra spectral broadening is
//! inferred. The underlying transient scene owns its finite-source and frozen-
//! medium assumptions; unsupported time-dependent source types remain outside it.
//! Explicit sequential-lens entry points include authored lens group transit
//! and place any configured analyzer at the physical sensor before the lens.

use super::sensor::CurvedSensorPixelSample;
use super::{
    CameraRadianceError, CurvedCameraRaster, CurvedPixelSample, CurvedRasterError,
    CurvedTraceSettings, philox4x32_10, unit,
};
use crate::polarization::dielectric::scene::curved::camera::transient::{
    TransientCameraError, TransientCurvedCameraScene,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use fs_exec::Cx;
use fs_math::dd::Dd;

pub mod area_lights;
pub mod exposure;
#[cfg(feature = "volumes")]
pub mod participating;
/// Retarded calibrated point fixtures on native CIE/QE paths.
#[cfg(feature = "photometry")]
pub mod point_lights;
mod sequential;
pub use sequential::TransientSequentialPixelSample;

#[derive(Clone, Debug, PartialEq)]
pub enum TransientRasterError {
    Raster(CurvedRasterError),
    Source(TransientCameraError),
}

impl From<CurvedRasterError> for TransientRasterError {
    fn from(error: CurvedRasterError) -> Self {
        Self::Raster(error)
    }
}
impl core::fmt::Display for TransientRasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient refractive pixel: {self:?}")
    }
}
impl core::error::Error for TransientRasterError {}

fn admit_analyzer_band(
    scene: &TransientCurvedCameraScene<'_, '_, '_>,
    min_nm: f64,
    max_nm: f64,
) -> Result<(), TransientRasterError> {
    if let Some(analyzer) = scene.analyzer_response() {
        analyzer
            .admit_band(min_nm, max_nm)
            .map_err(|error| TransientRasterError::Source(error.into()))?;
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientPixelSample {
    observation_time_s: Dd,
    pixel: CurvedPixelSample,
}

impl TransientPixelSample {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn pixel_sample(&self) -> &CurvedPixelSample {
        &self.pixel
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientSensorPixelSample<const CHANNELS: usize> {
    observation_time_s: Dd,
    pixel: CurvedSensorPixelSample<CHANNELS>,
}

impl<const C: usize> TransientSensorPixelSample<C> {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn pixel_sample(&self) -> &CurvedSensorPixelSample<C> {
        &self.pixel
    }
}

impl CurvedCameraRaster<'_> {
    /// Original CIE/XYZ pixel integration with independently retarded finite
    /// emission at each of the four native wavelength-specific paths.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_transient_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedCameraScene<'_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<TransientPixelSample, TransientRasterError> {
        admit_analyzer_band(scene, LAMBDA_MIN, LAMBDA_MAX)?;
        let key = [seed as u32, (seed >> 32) as u32];
        let mut source_error = None;
        let mut admitted_time = None;
        let result = self.sample_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let result = scene.sample_ray_with_samples_at_time(
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
                );
                match result {
                    Ok(timed) => {
                        admitted_time = Some(timed.observation_time_s());
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
        Ok(TransientPixelSample {
            observation_time_s: admitted_time.ok_or(CurvedRasterError::NumericalRange)?,
            pixel,
        })
    }

    /// Original source-bounded QE/energy/photon integration at one launch
    /// clock. Per-channel weights retain their actual vacuum wavelengths.
    /// Any failed lane refuses the entire observation without renormalizing
    /// surviving lanes or publishing partially accumulated detector channels.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_transient_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedCameraScene<'_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<TransientSensorPixelSample<C>, TransientRasterError> {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        admit_analyzer_band(scene, lo, hi)?;
        let key = [seed as u32, (seed >> 32) as u32];
        let mut source_error = None;
        let mut admitted_time = None;
        let result = self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let result = scene.sample_ray_with_samples_at_time(
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
                );
                match result {
                    Ok(timed) => {
                        admitted_time = Some(timed.observation_time_s());
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
        Ok(TransientSensorPixelSample {
            observation_time_s: admitted_time.ok_or(CurvedRasterError::NumericalRange)?,
            pixel,
        })
    }
}

#[cfg(test)]
mod tests;
