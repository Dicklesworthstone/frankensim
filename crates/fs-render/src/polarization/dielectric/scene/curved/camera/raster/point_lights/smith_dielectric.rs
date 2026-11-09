//! Original spectral/lens pixels with complete Smith glass point lighting.
use super::super::super::analyzer::AnalyzerResponse;
use super::super::super::point_lights::SmithDielectricPointLightCameraScene;
use super::super::microsurface::dielectric::{lighting_path, lighting_seed};
use super::super::{CameraAnalyzer, CurvedSequentialPixelSample};
use super::*;
use crate::camera::SequentialCameraBeam;

impl CurvedCameraRaster<'_> {
    /// Four independent full-glass source paths in the original hero-wavelength,
    /// pupil and XYZ integration owner. Every lane retains its source work.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_dielectric_point_lights(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<PointLightPixelSample, PointLightRasterError> {
        self.sample_point_with_launch(
            seed,
            pixel,
            sample,
            |cx, (ray, wavelength), samples| {
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

    /// Carry the resolved lens's reciprocal Stokes seed once into each live
    /// lane. Pupil/TIR misses keep zero source work and their original mass.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_dielectric_point_lights(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<PointLightPixelSample<CurvedSequentialPixelSample>, PointLightRasterError> {
        self.sample_point_with_launch(
            seed,
            pixel,
            sample,
            |cx, beam: SequentialCameraBeam, samples| {
                scene.sample_sequential_ray(
                    cx,
                    &beam,
                    initial_medium,
                    settings,
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
