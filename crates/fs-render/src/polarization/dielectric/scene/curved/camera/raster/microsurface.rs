//! Native spectral raster integration for path-only Smith conductor transport.

use super::super::microsurface::SmithCameraScene;
use super::*;

/// Explicit two-sided rough-glass camera transport on the native raster.
pub mod dielectric;

/// Each wavelength and geometric surface has an independent microscopic stream.
/// The variable walk length cannot shift the original raster, pupil, wavelength
/// or later surface draws. The seed/path counter map is one-to-one per lane.
pub(super) fn facet_draw(
    seed: u64,
    pixel: u32,
    sample: u32,
    lane: usize,
    event: u32,
    facet: u32,
) -> [f64; 3] {
    let block = philox4x32_10(
        [pixel, sample, event, facet],
        [
            seed as u32 ^ 0x736d_7261,
            (seed >> 32) as u32 ^ (0x7374_7200 + lane as u32),
        ],
    );
    core::array::from_fn(|i| unit(block[i]))
}

impl CurvedCameraRaster<'_> {
    /// Integrate four independent wavelengths with complete conductor
    /// microsurface paths. Uses the original raster, aperture and CIE owner.
    /// No source-connection MIS or hidden-path-to-solid-angle PDF conversion.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith(
        &self,
        cx: &Cx<'_>,
        scene: &SmithCameraScene<'_, '_, '_>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        self.sample_pixel_smith_with_analyzer(
            cx,
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            None,
        )
    }

    /// Apply the complete passive analyzer before each wavelength's original
    /// camera path, including every complex microscopic conductor response.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &SmithCameraScene<'_, '_, '_>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        self.sample_pixel_smith_with_analyzer(
            cx,
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            Some(AnalyzerResponse::Constant(analyzer)),
        )
    }

    /// A source-bound chromatic analyzer retains its absolute wavelength row;
    /// its complete visible support is admitted before tracing any pixel lane.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &SmithCameraScene<'_, '_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        analyzer.admit_band(LAMBDA_MIN, LAMBDA_MAX)?;
        self.sample_pixel_smith_with_analyzer(
            cx,
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            Some(AnalyzerResponse::Spectral(analyzer)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_pixel_smith_with_analyzer(
        &self,
        cx: &Cx<'_>,
        scene: &SmithCameraScene<'_, '_, '_>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<AnalyzerResponse<'_>>,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_pixel_with_transport(cx, seed, pixel, sample, |cx, ray, wavelength, lane| {
            let surface = |event| {
                let block = philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                core::array::from_fn(|i| unit(block[i]))
            };
            let facets = |event, facet| facet_draw(seed, pixel, sample, lane, event, facet);
            if let Some(analyzer) = analyzer {
                scene.sample_ray_analyzed_with_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    analyzer.at(wavelength)?,
                    settings,
                    surface,
                    facets,
                )
            } else {
                scene.sample_ray_with_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    settings,
                    surface,
                    facets,
                )
            }
        })
    }

    /// Retain the resolved lens's full reciprocal Stokes seed and geometry
    /// weight on every live wavelength. Physical pupil/TIR misses remain zero
    /// draws, while any failed scene/facet budget refuses the entire pixel.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith(
        &self,
        cx: &Cx<'_>,
        scene: &SmithCameraScene<'_, '_, '_>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialPixelSample, CurvedRasterError> {
        self.sample_sequential_smith_with_analyzer(
            cx,
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &SmithCameraScene<'_, '_, '_>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialPixelSample, CurvedRasterError> {
        self.sample_sequential_smith_with_analyzer(
            cx,
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            Some(AnalyzerResponse::Constant(analyzer)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_sequential_smith_with_analyzer(
        &self,
        cx: &Cx<'_>,
        scene: &SmithCameraScene<'_, '_, '_>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<AnalyzerResponse<'_>>,
    ) -> Result<CurvedSequentialPixelSample, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_sequential_pixel_with_timing(
            cx,
            seed,
            pixel,
            sample,
            analyzer,
            false,
            |cx, beam, lane| {
                scene.sample_sequential_ray_with_samples(
                    cx,
                    beam,
                    initial_medium,
                    settings,
                    |event| {
                        let block =
                            philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |event, facet| facet_draw(seed, pixel, sample, lane, event, facet),
                )
            },
        )
    }
}
