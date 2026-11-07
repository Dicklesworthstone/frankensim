//! Passive analyzer measurements on the original geometric spectral raster.
//!
//! Every wavelength retains its own ray and full analyzed transport. The
//! original packet, counter streams, CIE or QE observer, and lens/pupil weights
//! are shared. The authored analyzer axis is in world coordinates at the
//! geometric camera launch, projected separately for each sampled ray.

use super::super::CameraAnalyzer;
use super::sensor::CurvedSensorPixelSample;
use super::*;
use crate::spectral::detector::QuantumEfficiencySensor;

#[cfg(feature = "volumes")]
mod molecular;

impl CurvedCameraRaster<'_> {
    /// Observe a passive analyzer through the original four refractive paths.
    /// XYZ is the existing equal-energy integral of analyzed radiance, before
    /// exposure/display conversion. A failing spectral lane publishes no pixel.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        self.sample_pixel_with_analyzer(
            cx,
            scene,
            initial_medium,
            AnalyzerResponse::Constant(analyzer),
            settings,
            seed,
            pixel,
            sample,
        )
    }

    /// Integrate a source-owned chromatic polarization response on the same
    /// visible spectral paths. The complete visible band is admitted before
    /// tracing; each absolute row acts before scene transport, not on XYZ.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        analyzer.admit_band(LAMBDA_MIN, LAMBDA_MAX)?;
        self.sample_pixel_with_analyzer(
            cx,
            scene,
            initial_medium,
            AnalyzerResponse::Spectral(analyzer),
            settings,
            seed,
            pixel,
            sample,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_pixel_with_analyzer(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        initial_medium: usize,
        analyzer: AnalyzerResponse<'_>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_pixel_with_transport(cx, seed, pixel, sample, |cx, ray, wavelength, lane| {
            scene.sample_ray_analyzed_with_samples(
                cx,
                ray,
                wavelength,
                initial_medium,
                analyzer.at(wavelength)?,
                settings,
                |event| {
                    let block =
                        philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                    core::array::from_fn(|i| unit(block[i]))
                },
            )
        })
    }

    /// Integrate analyzed radiance against the original source-bounded QE and
    /// QE*lambda/(h*c) responses, including admitted infrared bands. All channels
    /// share the same four physical paths; none applies an extra analyzer half,
    /// CIE normalization, etendue or electronic gain after propagation.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSensorPixelSample<C>, CurvedRasterError> {
        self.sample_sensor_with_analyzer(
            cx,
            scene,
            sensor,
            initial_medium,
            AnalyzerResponse::Constant(analyzer),
            settings,
            seed,
            pixel,
            sample,
        )
    }

    /// A chromatic polarization response followed by the original QE channels.
    /// The analyzer must cover the entire authored sensor band, including any
    /// infrared support; its absolute row and QE are each applied exactly once.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_spectrally_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSensorPixelSample<C>, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        let (lo, hi) = sensor.wavelength_bounds_nm();
        analyzer.admit_band(lo, hi)?;
        self.sample_sensor_with_analyzer(
            cx,
            scene,
            sensor,
            initial_medium,
            AnalyzerResponse::Spectral(analyzer),
            settings,
            seed,
            pixel,
            sample,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_sensor_with_analyzer<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: AnalyzerResponse<'_>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSensorPixelSample<C>, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                scene.sample_ray_analyzed_with_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    analyzer.at(wavelength)?,
                    settings,
                    |event| {
                        let block =
                            philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                )
            },
        )
    }
}
