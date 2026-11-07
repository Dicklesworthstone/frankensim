//! Original spectral molecular pixels with a passive polarization analyzer.

use super::super::participating::RayleighPixelSample;
use super::super::sensor::CurvedRayleighSensorPixelSample;
use super::*;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighTraceSettings;

impl CurvedCameraRaster<'_> {
    /// Integrate four original analyzed molecular paths against the original
    /// XYZ observer. Each lane retains its own free-flight ordinal, independent
    /// acceptance variate and complete native path. Failure publishes no pixel.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPixelSample, CurvedRasterError> {
        self.sample_pixel_rayleigh_with_analyzer(
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

    /// Resolve a passive spectral detector before each native molecular walk.
    /// The entire visible band must be supported. Free-flight and scattering
    /// proposals retain their original counters and complete Mueller response.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_spectrally_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPixelSample, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        analyzer.admit_band(LAMBDA_MIN, LAMBDA_MAX)?;
        self.sample_pixel_rayleigh_with_analyzer(
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
    fn sample_pixel_rayleigh_with_analyzer(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: AnalyzerResponse<'_>,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPixelSample, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let result = self.sample_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let traced = scene.sample_ray_analyzed_with_tracking_samples(
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
                    |flight| {
                        philox4x32_10([pixel, sample, flight, 0x6372_7350 + lane as u32], key)
                            .map(unit)
                    },
                )?;
                volume_events[lane] = traced.volume_events();
                free_flight_draws[lane] = traced.free_flight_draws();
                Ok(*traced.radiance_sample())
            },
        )?;
        Ok(RayleighPixelSample::from_transport(
            result,
            volume_events,
            free_flight_draws,
        ))
    }

    /// Original source-bounded QE and QE*lambda/(h*c) integration, using the
    /// same four analyzed molecular paths for every detector channel. No
    /// second analyzer, scattering survival, index or exposure factor is added.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_rayleigh_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedRayleighSensorPixelSample<C>, CurvedRasterError> {
        self.sample_sensor_rayleigh_with_analyzer(
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

    /// Integrate the wavelength-dependent analyzed molecular radiance through
    /// the original QE/photon observer. All channels share the same four paths;
    /// the complete sensor band must lie within the authored detector table.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_rayleigh_spectrally_analyzed<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: &SpectralCameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedRayleighSensorPixelSample<C>, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        let (lo, hi) = sensor.wavelength_bounds_nm();
        analyzer.admit_band(lo, hi)?;
        self.sample_sensor_rayleigh_with_analyzer(
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
    fn sample_sensor_rayleigh_with_analyzer<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: AnalyzerResponse<'_>,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedRayleighSensorPixelSample<C>, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let result = self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let traced = scene.sample_ray_analyzed_with_tracking_samples(
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
                    |flight| {
                        philox4x32_10([pixel, sample, flight, 0x6372_7350 + lane as u32], key)
                            .map(unit)
                    },
                )?;
                volume_events[lane] = traced.volume_events();
                free_flight_draws[lane] = traced.free_flight_draws();
                Ok(*traced.radiance_sample())
            },
        )?;
        Ok(CurvedRayleighSensorPixelSample::from_transport(
            result,
            volume_events,
            free_flight_draws,
        ))
    }
}
