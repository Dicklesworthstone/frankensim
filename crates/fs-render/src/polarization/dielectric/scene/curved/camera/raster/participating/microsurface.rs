//! Full Smith surface/volume paths through the original spectral observers.

use super::super::microsurface::{dielectric::glass_draw, facet_draw};
use super::super::sensor::{CurvedRayleighSensorPixelSample, CurvedSequentialSensorPixelSample};
use super::super::{
    AnalyzerResponse, CameraAnalyzer, CameraRadianceError, CameraRadianceSample,
    CurvedSequentialPixelSample, CurvedSequentialRayleighPixelSample,
};
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::charts::Ray;
use crate::polarization::dielectric::scene::curved::camera::participating::{
    RayleighCameraSample, microsurface::SmithDielectricRayleighCameraScene,
};
use crate::spectral::detector::QuantumEfficiencySensor;

/// Original surface-only Smith streams plus the existing independent physical
/// free-flight stream. A variable microscopic walk cannot consume a volume
/// variate or shift a later surface, pixel, wavelength or pupil sample.
#[derive(Clone, Copy)]
struct LaneSamples {
    seed: u64,
    pixel: u32,
    sample: u32,
    lane: usize,
}

impl LaneSamples {
    fn surfaces(self, event: u32) -> [f64; 3] {
        let block = philox4x32_10(
            [
                self.pixel,
                self.sample,
                event,
                0x6372_7340 + self.lane as u32,
            ],
            [self.seed as u32, (self.seed >> 32) as u32],
        );
        core::array::from_fn(|i| unit(block[i]))
    }
    fn flights(self, flight: u32) -> [f64; 4] {
        philox4x32_10(
            [
                self.pixel,
                self.sample,
                flight,
                0x6372_7350 + self.lane as u32,
            ],
            [self.seed as u32, (self.seed >> 32) as u32],
        )
        .map(unit)
    }
    fn opaque(self, event: u32, facet: u32) -> [f64; 3] {
        facet_draw(self.seed, self.pixel, self.sample, self.lane, event, facet)
    }
    fn glass(self, event: u32, facet: u32) -> [f64; 4] {
        glass_draw(self.seed, self.pixel, self.sample, self.lane, event, facet)
    }

    #[allow(clippy::too_many_arguments)]
    fn trace_ray(
        self,
        cx: &Cx<'_>,
        scene: &SmithDielectricRayleighCameraScene<'_, '_, '_, '_, '_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedRayleighTraceSettings,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        match analyzer {
            Some(analyzer) => scene.sample_ray_analyzed_with_tracking_samples(
                cx,
                ray,
                wavelength_nm,
                initial_medium,
                analyzer,
                settings,
                |event| self.surfaces(event),
                |flight| self.flights(flight),
                |event, facet| self.opaque(event, facet),
                |event, facet| self.glass(event, facet),
            ),
            None => scene.sample_ray_with_tracking_samples(
                cx,
                ray,
                wavelength_nm,
                initial_medium,
                settings,
                |event| self.surfaces(event),
                |flight| self.flights(flight),
                |event, facet| self.opaque(event, facet),
                |event, facet| self.glass(event, facet),
            ),
        }
    }

    fn trace_lens(
        self,
        cx: &Cx<'_>,
        scene: &SmithDielectricRayleighCameraScene<'_, '_, '_, '_, '_>,
        beam: &SequentialCameraBeam,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
    ) -> Result<RayleighCameraSample, CameraRadianceError> {
        scene.sample_sequential_ray_with_tracking_samples(
            cx,
            beam,
            initial_medium,
            settings,
            |event| self.surfaces(event),
            |flight| self.flights(flight),
            |event, facet| self.opaque(event, facet),
            |event, facet| self.glass(event, facet),
        )
    }
}

#[derive(Default)]
struct VolumeWork {
    events: [u32; CAMERA_WAVELENGTHS],
    flights: [u32; CAMERA_WAVELENGTHS],
}

impl VolumeWork {
    fn record(&mut self, lane: usize, traced: RayleighCameraSample) -> CameraRadianceSample {
        self.events[lane] = traced.volume_events();
        self.flights[lane] = traced.free_flight_draws();
        *traced.radiance_sample()
    }
}

impl CurvedCameraRaster<'_> {
    /// Observe full rough-metal/glass returns before and after physical volume
    /// collisions with the original four-wavelength visible packet. Optional
    /// analyzer polarization enters before every scene Mueller event.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_smith_dielectric_rayleigh(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricRayleighCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPixelSample, CurvedRasterError> {
        let mut work = VolumeWork::default();
        let result = self.sample_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let draws = LaneSamples {
                    seed,
                    pixel,
                    sample,
                    lane,
                };
                let traced = draws.trace_ray(
                    cx,
                    scene,
                    ray,
                    wavelength,
                    initial_medium,
                    analyzer,
                    settings,
                )?;
                Ok(work.record(lane, traced))
            },
        )?;
        Ok(RayleighPixelSample::from_transport(
            result,
            work.events,
            work.flights,
        ))
    }

    /// Integrate the same full paths against the authored sum-of-QE proposal,
    /// including source-supported infrared bands. Detector channels share the
    /// four paths; QE and photon-energy conversion keep their native owners.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_smith_dielectric_rayleigh<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricRayleighCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedRayleighSensorPixelSample<C>, CurvedRasterError> {
        let mut work = VolumeWork::default();
        let result = self.sample_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let draws = LaneSamples {
                    seed,
                    pixel,
                    sample,
                    lane,
                };
                let traced = draws.trace_ray(
                    cx,
                    scene,
                    ray,
                    wavelength,
                    initial_medium,
                    analyzer,
                    settings,
                )?;
                Ok(work.record(lane, traced))
            },
        )?;
        Ok(CurvedRayleighSensorPixelSample::from_transport(
            result,
            work.events,
            work.flights,
        ))
    }

    /// Carry each resolved wavelength-specific lens beam through the complete
    /// surface/volume path. Analyzer, lens Stokes and pupil weights each enter
    /// once. Clipped lanes retain their zero contribution and original mass.
    /// This steady observation does not require or invent group delay.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_smith_dielectric_rayleigh(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricRayleighCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialRayleighPixelSample<CurvedSequentialPixelSample>, CurvedRasterError>
    {
        let mut work = VolumeWork::default();
        let result = self.sample_sequential_pixel_with_timing(
            cx,
            seed,
            pixel,
            sample,
            analyzer.map(AnalyzerResponse::Constant),
            false,
            |cx, beam, lane| {
                let traced = LaneSamples {
                    seed,
                    pixel,
                    sample,
                    lane,
                }
                .trace_lens(cx, scene, beam, initial_medium, settings)?;
                Ok(work.record(lane, traced))
            },
        )?;
        Ok(CurvedSequentialRayleighPixelSample::from_transport(
            result,
            work.events,
            work.flights,
        ))
    }

    /// Physical multichannel energy/photon observation through the actual
    /// lens and the same full Smith/volume transport, with per-lane native
    /// volume-work receipts. Any live-path error invalidates the whole pixel.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_smith_dielectric_rayleigh<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &SmithDielectricRayleighCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        analyzer: Option<CameraAnalyzer>,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        CurvedSequentialRayleighPixelSample<CurvedSequentialSensorPixelSample<C>>,
        CurvedRasterError,
    > {
        let mut work = VolumeWork::default();
        let result = self.sample_sequential_sensor_with_timing(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            analyzer.map(AnalyzerResponse::Constant),
            false,
            |cx, beam, lane| {
                let traced = LaneSamples {
                    seed,
                    pixel,
                    sample,
                    lane,
                }
                .trace_lens(cx, scene, beam, initial_medium, settings)?;
                Ok(work.record(lane, traced))
            },
        )?;
        Ok(CurvedSequentialRayleighPixelSample::from_transport(
            result,
            work.events,
            work.flights,
        ))
    }
}

#[cfg(test)]
mod tests;
