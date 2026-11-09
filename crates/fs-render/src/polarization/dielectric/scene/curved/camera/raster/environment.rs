//! Full rough-glass environment transport on the original spectral raster.
use super::super::environment::sampling::{
    EnvironmentCameraError, EnvironmentCameraSample, EnvironmentCameraSettings,
    SmithDielectricEnvironmentCameraScene,
};
use super::microsurface::dielectric::{lighting_path, lighting_seed};
use super::*;

mod detector;

/// Preserve the complete native environment receipt in each live wavelength.
/// Lens-clipped lanes are `None`, with zero contribution and original mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmithEnvironmentPixelSample<P = CurvedPixelSample> {
    pixel: P,
    lanes: [Option<EnvironmentCameraSample>; CAMERA_WAVELENGTHS],
}
impl<P> SmithEnvironmentPixelSample<P> {
    pub const fn pixel_sample(&self) -> &P {
        &self.pixel
    }
    pub const fn environment_samples(
        &self,
    ) -> &[Option<EnvironmentCameraSample>; CAMERA_WAVELENGTHS] {
        &self.lanes
    }
}

/// Retain the original source failure and exact spectral lane. No successful
/// prefix or clipped lane replaces an unresolved numerical/source failure.
#[derive(Debug)]
pub enum SmithEnvironmentRasterError {
    Raster(CurvedRasterError),
    Transport {
        lane: usize,
        error: EnvironmentCameraError,
    },
}
impl core::fmt::Display for SmithEnvironmentRasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "full Smith environment pixel: {self:?}")
    }
}
impl core::error::Error for SmithEnvironmentRasterError {}

impl CurvedCameraRaster<'_> {
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_dielectric_environment(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: EnvironmentCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<SmithEnvironmentPixelSample, SmithEnvironmentRasterError> {
        environment_pixel(
            |cx, (ray, wavelength), lane| {
                let seed = lighting_seed(seed, lane);
                let path = lighting_path(pixel, sample);
                match analyzer {
                    Some(a) => scene.sample_ray_analyzed(
                        cx,
                        ray,
                        wavelength,
                        initial_medium,
                        a,
                        settings,
                        seed,
                        path,
                    ),
                    None => {
                        scene.sample_ray(cx, ray, wavelength, initial_medium, settings, seed, path)
                    }
                }
            },
            |trace| {
                self.sample_pixel_with_transport(
                    cx,
                    seed,
                    pixel,
                    sample,
                    |cx, ray, wavelength, lane| trace(cx, (ray, wavelength), lane),
                )
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_dielectric_environment(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: EnvironmentCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<SmithEnvironmentPixelSample<CurvedSequentialPixelSample>, SmithEnvironmentRasterError>
    {
        environment_pixel(
            |cx, beam: SequentialCameraBeam, lane| {
                scene.sample_sequential_ray(
                    cx,
                    &beam,
                    initial_medium,
                    settings,
                    lighting_seed(seed, lane),
                    lighting_path(pixel, sample),
                )
            },
            |trace| {
                self.sample_sequential_pixel_with_timing(
                    cx,
                    seed,
                    pixel,
                    sample,
                    analyzer.map(AnalyzerResponse::Constant),
                    false,
                    |cx, beam, lane| trace(cx, *beam, lane),
                )
            },
        )
    }
}

#[allow(clippy::type_complexity)]
fn environment_pixel<P, L>(
    mut trace: impl FnMut(&Cx<'_>, L, usize) -> Result<EnvironmentCameraSample, EnvironmentCameraError>,
    integrate: impl FnOnce(
        &mut dyn FnMut(&Cx<'_>, L, usize) -> Result<CameraRadianceSample, CameraRadianceError>,
    ) -> Result<P, CurvedRasterError>,
) -> Result<SmithEnvironmentPixelSample<P>, SmithEnvironmentRasterError> {
    let mut lanes = [None; CAMERA_WAVELENGTHS];
    let mut failure = None;
    let result = integrate(&mut |cx, launch, lane| {
        let traced = match trace(cx, launch, lane) {
            Ok(value) => value,
            Err(error) => {
                failure = Some((lane, error));
                return Err(CameraRadianceError::NumericalRange);
            }
        };
        let radiance = *traced.radiance_sample();
        lanes[lane] = Some(traced);
        Ok(radiance)
    });
    if let Some((lane, error)) = failure {
        return Err(SmithEnvironmentRasterError::Transport { lane, error });
    }
    Ok(SmithEnvironmentPixelSample {
        pixel: result.map_err(SmithEnvironmentRasterError::Raster)?,
        lanes,
    })
}

#[cfg(test)]
mod tests;
