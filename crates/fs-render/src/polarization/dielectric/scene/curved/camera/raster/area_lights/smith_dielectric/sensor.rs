//! Calibrated photon/energy observations of complete rough-glass finite-area paths.

use super::super::super::sensor::{CurvedSensorPixelSample, CurvedSequentialSensorPixelSample};
use super::*;
use crate::spectral::detector::QuantumEfficiencySensor;

impl CurvedCameraRaster<'_> {
    /// Sample the original sum-of-QE proposal and share its four full Smith
    /// paths across all detector channels. No extra source or facet sampling.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_smith_dielectric_area_lights<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricAreaLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample<CurvedSensorPixelSample<C>>, AreaLightRasterError> {
        self.sample_smith_dielectric_area_transport(
            scene,
            initial_medium,
            analyzer.map(AnalyzerResponse::Constant),
            settings,
            workspace,
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

    /// Preserve the authored spectral analyzer and source Stokes before QE and
    /// photon weighting. Admit the entire detector band before tracing a lane.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_smith_dielectric_area_lights_spectrally_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricAreaLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample<CurvedSensorPixelSample<C>>, AreaLightRasterError> {
        let (lo, hi) = sensor.wavelength_bounds_nm();
        admit_analyzer(cx, analyzer, lo, hi)?;
        self.sample_smith_dielectric_area_transport(
            scene,
            initial_medium,
            Some(AnalyzerResponse::Spectral(analyzer)),
            settings,
            workspace,
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

    /// Carry each complete resolved-lens beam through the original full Smith
    /// finite-area estimator; apply its pupil and QE/photon weights once afterward.
    /// Clipped wavelengths retain one-quarter proposal mass and zero work.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_smith_dielectric_area_lights<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricAreaLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample<CurvedSequentialSensorPixelSample<C>>, AreaLightRasterError>
    {
        self.sample_smith_dielectric_area_lens_transport(
            scene,
            initial_medium,
            settings,
            workspace,
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

    /// Apply the wavelength-dependent analyzer at the sensor before coated-lens
    /// and rough-glass Mueller transport. The original packet admits the full
    /// band even when every pupil lane is clipped; no timing model is required.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_smith_dielectric_area_lights_spectrally_analyzed<
        const C: usize,
    >(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricAreaLightCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample<CurvedSequentialSensorPixelSample<C>>, AreaLightRasterError>
    {
        self.sample_smith_dielectric_area_lens_transport(
            scene,
            initial_medium,
            settings,
            workspace,
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
}
