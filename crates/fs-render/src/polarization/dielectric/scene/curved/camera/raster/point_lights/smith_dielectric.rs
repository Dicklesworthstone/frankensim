//! Original CIE/QE lens pixels with complete Smith glass point lighting.
use super::super::super::analyzer::AnalyzerResponse;
use super::super::super::point_lights::SmithDielectricPointLightCameraScene;
use super::super::microsurface::dielectric::{lighting_path, lighting_seed};
use super::super::{
    CameraAnalyzer, CurvedSequentialPixelSample, LAMBDA_MAX, LAMBDA_MIN, SpectralCameraAnalyzer,
};
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::charts::Ray;

mod sensor;

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
        self.sample_smith_dielectric_point_transport(
            scene,
            initial_medium,
            analyzer.map(AnalyzerResponse::Constant),
            settings,
            seed,
            pixel,
            sample,
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

    /// Admit the complete visible analyzer band, then retain each wavelength's
    /// absolute reciprocal Stokes row through every glass facet and source.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_dielectric_point_lights_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<PointLightPixelSample, PointLightRasterError> {
        admit_analyzer(cx, analyzer, LAMBDA_MIN, LAMBDA_MAX)?;
        self.sample_smith_dielectric_point_transport(
            scene,
            initial_medium,
            Some(AnalyzerResponse::Spectral(analyzer)),
            settings,
            seed,
            pixel,
            sample,
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
        self.sample_smith_dielectric_point_lens_transport(
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
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

    /// Install the source-bounded chromatic detector row before the full lens.
    /// Active steady coatings require no group clock or timing approximation.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_dielectric_point_lights_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<PointLightPixelSample<CurvedSequentialPixelSample>, PointLightRasterError> {
        self.sample_smith_dielectric_point_lens_transport(
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_sequential_pixel_with_timing(
                    cx,
                    seed,
                    pixel,
                    sample,
                    Some(AnalyzerResponse::Spectral(analyzer)),
                    false,
                    |cx, beam, lane| trace(cx, *beam, lane),
                )
            },
        )
    }

    // CIE and QE consume the same seeded source owner. Wavelength selection and
    // detector weighting remain in their original raster integrators; neither
    // detector channel count nor stochastic facet length changes another lane.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_smith_dielectric_point_transport<P>(
        &self,
        scene: &SmithDielectricPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<AnalyzerResponse<'_>>,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        integrate: impl FnOnce(
            &mut dyn FnMut(
                &Cx<'_>,
                (Ray, f64),
                usize,
            ) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<PointLightPixelSample<P>, PointLightRasterError> {
        self.sample_point_with_launch(
            seed,
            pixel,
            sample,
            |cx, (ray, wavelength), samples| {
                let seed = lighting_seed(seed, samples.lane);
                let path = lighting_path(pixel, sample);
                match analyzer {
                    Some(response) => scene.sample_ray_analyzed(
                        cx,
                        ray,
                        wavelength,
                        initial_medium,
                        response.at(wavelength)?,
                        settings,
                        seed,
                        path,
                    ),
                    None => {
                        scene.sample_ray(cx, ray, wavelength, initial_medium, settings, seed, path)
                    }
                }
            },
            integrate,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_smith_dielectric_point_lens_transport<P>(
        &self,
        scene: &SmithDielectricPointLightCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        integrate: impl FnOnce(
            &mut dyn FnMut(
                &Cx<'_>,
                SequentialCameraBeam,
                usize,
            ) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<PointLightPixelSample<P>, PointLightRasterError> {
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
            integrate,
        )
    }
}

// Admission precedes the first ideal-camera lane. The original resolved-lens
// packet performs this same full-band check before any pupil-clipped lane.
fn admit_analyzer(
    cx: &Cx<'_>,
    analyzer: &SpectralCameraAnalyzer,
    min_nm: f64,
    max_nm: f64,
) -> Result<(), PointLightRasterError> {
    cx.checkpoint()
        .map_err(CameraRadianceError::from)
        .map_err(CurvedRasterError::from)
        .map_err(PointLightRasterError::Raster)?;
    analyzer
        .admit_band(min_nm, max_nm)
        .map_err(CurvedRasterError::from)
        .map_err(PointLightRasterError::Raster)
}

#[cfg(test)]
mod tests;
