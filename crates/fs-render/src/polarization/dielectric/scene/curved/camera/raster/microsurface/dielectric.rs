//! Four-wavelength/pupil raster integration for full rough dielectric paths.
use super::super::super::microsurface::SmithDielectricCameraScene;
use super::*;

fn glass_draw(seed: u64, pixel: u32, sample: u32, lane: usize, event: u32, facet: u32) -> [f64; 4] {
    let block = philox4x32_10(
        [pixel, sample, event, facet],
        [
            seed as u32 ^ 0x676c_6173,
            (seed >> 32) as u32 ^ (0x7372_6100 + lane as u32),
        ],
    );
    core::array::from_fn(|i| unit(block[i]))
}

impl CurvedCameraRaster<'_> {
    /// Integrate the native four-wavelength raster with full two-sided glass
    /// returns. An optional passive analyzer is installed before the walk.
    /// No ray is accepted after exhausting its microscopic facet budget.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_dielectric(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricCameraScene<'_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_pixel_with_transport(cx, seed, pixel, sample, |cx, ray, wavelength, lane| {
            let surfaces = |event| {
                let block = philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                core::array::from_fn(|i| unit(block[i]))
            };
            let opaque = |event, facet| facet_draw(seed, pixel, sample, lane, event, facet);
            let glass = |event, facet| glass_draw(seed, pixel, sample, lane, event, facet);
            match analyzer {
                Some(analyzer) => scene.sample_ray_analyzed_with_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    analyzer,
                    settings,
                    surfaces,
                    opaque,
                    glass,
                ),
                None => scene.sample_ray_with_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    settings,
                    surfaces,
                    opaque,
                    glass,
                ),
            }
        })
    }

    /// Preserve the original resolved lens's reciprocal Stokes seed and pupil
    /// weight. Pupil/TIR misses remain zero draws; scene failures refuse the pixel.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_dielectric(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricCameraScene<'_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialPixelSample, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_sequential_pixel_with_timing(
            cx,
            seed,
            pixel,
            sample,
            analyzer.map(AnalyzerResponse::Constant),
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
                    |event, facet| glass_draw(seed, pixel, sample, lane, event, facet),
                )
            },
        )
    }
}
