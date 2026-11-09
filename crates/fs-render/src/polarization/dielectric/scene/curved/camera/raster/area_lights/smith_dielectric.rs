//! Native finite-source receipts and spectral pixels for full rough glass.
use super::super::super::analyzer::AnalyzerResponse;
use super::super::super::area_lights::camera::SmithDielectricAreaLightCameraScene;
use super::super::microsurface::dielectric::{lighting_path, lighting_seed};
use super::super::{CameraAnalyzer, CurvedSequentialPixelSample};
use super::*;
use crate::camera::SequentialCameraBeam;

impl CurvedCameraRaster<'_> {
    /// Reuse native source/lane receipts, four hero wavelengths, pupil sampling
    /// and XYZ integration. The caller's one workspace is reused per live lane.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_dielectric_area_lights(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricAreaLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample, AreaLightRasterError> {
        self.sample_area_with_transport(
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, samples| {
                let seed = lighting_seed(seed, samples.lane);
                let path = lighting_path(pixel, sample);
                match analyzer {
                    Some(a) => scene.sample_ray_analyzed(
                        cx,
                        ray,
                        wavelength,
                        initial_medium,
                        a,
                        settings,
                        workspace,
                        seed,
                        path,
                    ),
                    None => scene.sample_ray(
                        cx,
                        ray,
                        wavelength,
                        initial_medium,
                        settings,
                        workspace,
                        seed,
                        path,
                    ),
                }
            },
            |trace| self.sample_pixel_with_transport(cx, seed, pixel, sample, trace),
        )
    }

    /// Preserve the resolved lens's native optical response and projected pupil
    /// weight. A clipped wavelength keeps its zero draw without tracing a scene.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_dielectric_area_lights(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricAreaLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample<CurvedSequentialPixelSample>, AreaLightRasterError> {
        self.sample_area_with_launch(
            seed,
            pixel,
            sample,
            |cx, beam: SequentialCameraBeam, samples| {
                scene.sample_sequential_ray(
                    cx,
                    &beam,
                    initial_medium,
                    settings,
                    workspace,
                    lighting_seed(seed, samples.lane),
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
