//! Calibrated detector and chromatic analyzer views of the full glass environment.

use super::super::sensor::{CurvedSensorPixelSample, CurvedSequentialSensorPixelSample};
use super::*;
use crate::spectral::detector::QuantumEfficiencySensor;

impl CurvedCameraRaster<'_> {
    /// Keep the complete visible analyzer band and each original environment
    /// receipt, including its deterministic complementary endpoint partition.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_dielectric_environment_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: EnvironmentCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<SmithEnvironmentPixelSample, SmithEnvironmentRasterError> {
        admit_analyzer(cx, analyzer, LAMBDA_MIN, LAMBDA_MAX)?;
        self.sample_smith_environment_detector_transport(
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

    /// Place the wavelength-dependent passive row at the original sensor,
    /// before the complete coated lens and full rough-glass source transport.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_dielectric_environment_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: EnvironmentCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<SmithEnvironmentPixelSample<CurvedSequentialPixelSample>, SmithEnvironmentRasterError>
    {
        self.sample_smith_environment_lens_detector_transport(
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

    /// Observe four native full-glass environment paths under the original
    /// sum-of-QE proposal; all channels reuse their source and geometry work.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_smith_dielectric_environment<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: EnvironmentCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<SmithEnvironmentPixelSample<CurvedSensorPixelSample<C>>, SmithEnvironmentRasterError>
    {
        self.sample_smith_environment_detector_transport(
            scene,
            initial_medium,
            analyzer.map(AnalyzerResponse::Constant),
            settings,
            seed,
            pixel,
            sample,
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

    /// Preserve source polarization and the authored detector row at every
    /// wavelength before original QE energy and lambda/(h*c) photon weighting.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_smith_dielectric_environment_spectrally_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: EnvironmentCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<SmithEnvironmentPixelSample<CurvedSensorPixelSample<C>>, SmithEnvironmentRasterError>
    {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        admit_analyzer(cx, analyzer, lo, hi)?;
        self.sample_smith_environment_detector_transport(
            scene,
            initial_medium,
            Some(AnalyzerResponse::Spectral(analyzer)),
            settings,
            seed,
            pixel,
            sample,
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

    /// Carry the resolved lens's complete reciprocal state into every live
    /// native environment lane. Clipped lanes retain zero work and /4 mass.
    /// Steady coating response requires no group-delay approximation.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_smith_dielectric_environment<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: EnvironmentCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        SmithEnvironmentPixelSample<CurvedSequentialSensorPixelSample<C>>,
        SmithEnvironmentRasterError,
    > {
        self.sample_smith_environment_lens_detector_transport(
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
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

    /// Admit the full sensor/analyzer band before optional pupil lanes, then
    /// retain chromatic lens and glass polarization in the original detectors.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_smith_dielectric_environment_spectrally_analyzed<
        const C: usize,
    >(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: EnvironmentCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        SmithEnvironmentPixelSample<CurvedSequentialSensorPixelSample<C>>,
        SmithEnvironmentRasterError,
    > {
        self.sample_smith_environment_lens_detector_transport(
            scene,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            |trace| {
                self.sample_sequential_sensor_with_timing(
                    cx,
                    sensor,
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

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_smith_environment_detector_transport<P>(
        &self,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<AnalyzerResponse<'_>>,
        settings: EnvironmentCameraSettings,
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
    ) -> Result<SmithEnvironmentPixelSample<P>, SmithEnvironmentRasterError> {
        environment_pixel(
            |cx, (ray, wavelength), lane| {
                let seed = lighting_seed(seed, lane);
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
    fn sample_smith_environment_lens_detector_transport<P>(
        &self,
        scene: &SmithDielectricEnvironmentCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        settings: EnvironmentCameraSettings,
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
    ) -> Result<SmithEnvironmentPixelSample<P>, SmithEnvironmentRasterError> {
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
            integrate,
        )
    }
}

fn admit_analyzer(
    cx: &Cx<'_>,
    analyzer: &SpectralCameraAnalyzer,
    min_nm: f64,
    max_nm: f64,
) -> Result<(), SmithEnvironmentRasterError> {
    cx.checkpoint()
        .map_err(CameraRadianceError::from)
        .map_err(CurvedRasterError::from)
        .map_err(SmithEnvironmentRasterError::Raster)?;
    analyzer
        .admit_band(min_nm, max_nm)
        .map_err(CurvedRasterError::from)
        .map_err(SmithEnvironmentRasterError::Raster)
}
