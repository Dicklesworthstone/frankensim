//! Native visible and calibrated detector pixels from complete Smith thermal
//! emission/reflection. The original raster owns wavelength and pupil measure;
//! the source camera owns independent thermal and continuation facet streams.

use super::super::thermal::smith::{SmithThermalCameraError, SmithThermalCameraScene};
use super::microsurface::dielectric::{lighting_path, lighting_seed};
use super::sensor::{CurvedSensorPixelSample, CurvedSequentialSensorPixelSample};
use super::*;
use crate::spectral::detector::QuantumEfficiencySensor;

/// Retain the original source error and its spectral lane. A failed live lane
/// aborts the entire pixel; it cannot become black or a partial detector result.
#[derive(Clone, Debug, PartialEq)]
pub enum SmithThermalRasterError {
    Raster(CurvedRasterError),
    Transport {
        lane: usize,
        error: SmithThermalCameraError,
    },
}
impl core::fmt::Display for SmithThermalRasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Smith thermal pixel: {self:?}")
    }
}
impl core::error::Error for SmithThermalRasterError {}

impl CurvedCameraRaster<'_> {
    /// Four independent native thermal/continuation paths with original CIE
    /// integration. The optional passive analyzer precedes scene transport.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_thermal(
        &self,
        cx: &Cx<'_>,
        scene: &SmithThermalCameraScene<'_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, SmithThermalRasterError> {
        thermal_pixel(
            |cx, (ray, wavelength), lane| {
                trace_ray(
                    cx,
                    scene,
                    ray,
                    wavelength,
                    initial_medium,
                    analyzer,
                    settings,
                    lighting_seed(seed, lane),
                    lighting_path(pixel, sample),
                )
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

    /// All QE energy and photon channels share the original four spectral
    /// paths. Changing the detector channel count never retraces a thermal walk.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_smith_thermal<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithThermalCameraScene<'_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSensorPixelSample<C>, SmithThermalRasterError> {
        thermal_pixel(
            |cx, (ray, wavelength), lane| {
                trace_ray(
                    cx,
                    scene,
                    ray,
                    wavelength,
                    initial_medium,
                    analyzer,
                    settings,
                    lighting_seed(seed, lane),
                    lighting_path(pixel, sample),
                )
            },
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
        )
    }

    /// Put the analyzer at the sensor before the resolved lens. Every live
    /// wavelength retains its complete reciprocal Stokes seed; clipped lanes
    /// keep zero contribution and their original proposal mass. Steady thermal
    /// transport never requests a lens group delay or invents a clock.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_thermal(
        &self,
        cx: &Cx<'_>,
        scene: &SmithThermalCameraScene<'_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialPixelSample, SmithThermalRasterError> {
        thermal_pixel(
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

    /// Original multichannel QE and photon measures through the physical lens,
    /// with thermal emission and reflection composed before detector weighting.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_smith_thermal<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithThermalCameraScene<'_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialSensorPixelSample<C>, SmithThermalRasterError> {
        thermal_pixel(
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
                self.sample_sequential_sensor_with_timing(
                    cx,
                    sensor,
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

#[allow(clippy::too_many_arguments)]
fn trace_ray(
    cx: &Cx<'_>,
    scene: &SmithThermalCameraScene<'_, '_, '_>,
    ray: Ray,
    wavelength: f64,
    initial_medium: usize,
    analyzer: Option<CameraAnalyzer>,
    settings: CurvedTraceSettings,
    seed: u64,
    path: u64,
) -> Result<CameraRadianceSample, SmithThermalCameraError> {
    match analyzer {
        Some(analyzer) => scene.sample_ray_analyzed(
            cx,
            ray,
            wavelength,
            initial_medium,
            analyzer,
            settings,
            seed,
            path,
        ),
        None => scene.sample_ray(cx, ray, wavelength, initial_medium, settings, seed, path),
    }
}

// The established raster/lens owners return CameraRadianceError. Preserve the
// richer source failure around that seam, then restore it before any result
// can leave this function. Pupil misses never invoke this callback.
#[allow(clippy::type_complexity)]
fn thermal_pixel<P, L>(
    mut trace: impl FnMut(&Cx<'_>, L, usize) -> Result<CameraRadianceSample, SmithThermalCameraError>,
    integrate: impl FnOnce(
        &mut dyn FnMut(&Cx<'_>, L, usize) -> Result<CameraRadianceSample, CameraRadianceError>,
    ) -> Result<P, CurvedRasterError>,
) -> Result<P, SmithThermalRasterError> {
    let mut failure = None;
    let result = integrate(&mut |cx, launch, lane| {
        trace(cx, launch, lane).map_err(|error| {
            failure = Some((lane, error));
            CameraRadianceError::NumericalRange
        })
    });
    match failure {
        Some((lane, error)) => Err(SmithThermalRasterError::Transport { lane, error }),
        None => result.map_err(SmithThermalRasterError::Raster),
    }
}

#[cfg(test)]
mod tests;
