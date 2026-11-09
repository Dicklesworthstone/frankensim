//! Four-wavelength/pupil raster integration for full rough dielectric paths.
use super::super::super::microsurface::SmithDielectricCameraScene;
use super::*;

pub(in super::super) fn glass_draw(
    seed: u64,
    pixel: u32,
    sample: u32,
    lane: usize,
    event: u32,
    facet: u32,
) -> [f64; 4] {
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

// Light wrappers retain their seeded optical/source domains. The four lanes
// occupy the low key's top two bits: source binding/channel addresses use the
// high key, and the fixed low-key domains differ below these two bits. Placing
// the lane in the high key would alias (lane 0, binding 1) with (lane 1, binding 0).
// Pixel/sample retain the complete path counter; pupil/wavelength draws use the
// original seed. The caller supplies one of the four native spectral lanes.
pub(in super::super) fn lighting_seed(seed: u64, lane: usize) -> u64 {
    seed ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30)
}
pub(in super::super) fn lighting_path(pixel: u32, sample: u32) -> u64 {
    (u64::from(pixel) << 32) | u64::from(sample)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spectral_lanes_do_not_alias_point_bindings_or_area_channels() {
        let counter = [19, 23, 7, 11];
        // Effective source keys from the original full-glass point/area owners.
        // For area sources, address 1 is binding 0's second thermal channel.
        for seed in [0, u64::MAX, 0x7a93_51c0_24ef_08d6] {
            for (domain, high_domain) in [
                (0x706d_6e65, 0),
                (0x7067_6e65, 0),
                (0x7361_7265, 0),
                (0x6761_7265, 0x6167_6c73),
            ] {
                let key = |lane, address| {
                    let lane_seed = lighting_seed(seed, lane);
                    [
                        lane_seed as u32 ^ domain,
                        (lane_seed >> 32) as u32 ^ high_domain ^ address,
                    ]
                };
                for first in 0..CAMERA_WAVELENGTHS {
                    for second in first + 1..CAMERA_WAVELENGTHS {
                        let a = key(first, 1);
                        let b = key(second, 0);
                        assert_ne!(a, b, "lane/source addresses must be distinct");
                        assert_ne!(philox4x32_10(counter, a), philox4x32_10(counter, b));
                    }
                }
            }
        }
    }
}
