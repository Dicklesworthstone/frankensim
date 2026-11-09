//! Resolved sensor/lens polarization through native molecular scattering.

use super::super::sensor::CurvedSequentialSensorPixelSample;
use super::*;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighTraceSettings;
use crate::spectral::detector::QuantumEfficiencySensor;

/// A visible or QE sequential-lens sample with original per-lane volume work.
/// Rejected optical lanes remain absent in `pixel_sample()` and have zero
/// volume events and free-flight draws; surviving lanes keep their own counts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedSequentialRayleighPixelSample<P> {
    pixel: P,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
}

impl<P> CurvedSequentialRayleighPixelSample<P> {
    pub(in super::super) const fn from_transport(
        pixel: P,
        volume_events: [u32; CAMERA_WAVELENGTHS],
        free_flight_draws: [u32; CAMERA_WAVELENGTHS],
    ) -> Self {
        Self {
            pixel,
            volume_events,
            free_flight_draws,
        }
    }

    pub const fn pixel_sample(&self) -> &P {
        &self.pixel
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    /// Includes rejected tracking candidates and surface-interrupted flights.
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }
}

impl CurvedCameraRaster<'_> {
    /// Carry the resolved lens's complete reciprocal response through the
    /// native molecular Mueller walk, retaining the original visible packet.
    /// Zero scattering preserves the entire surface-only sequential sample;
    /// unresolved volume paths abort the pixel, including its other wavelengths.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialRayleighPixelSample<CurvedSequentialPixelSample>, CurvedRasterError>
    {
        self.sample_sequential_rayleigh_pixel(
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

    /// The same molecular observation with an authored passive analyzer at
    /// the actual sensor, before reciprocal lens and scene propagation.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_analyzed(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedSequentialRayleighPixelSample<CurvedSequentialPixelSample>, CurvedRasterError>
    {
        self.sample_sequential_rayleigh_pixel(
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
    fn sample_sequential_rayleigh_pixel(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<AnalyzerResponse<'_>>,
    ) -> Result<CurvedSequentialRayleighPixelSample<CurvedSequentialPixelSample>, CurvedRasterError>
    {
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let pixel = self.sample_sequential_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            analyzer,
            tracking_transport(
                scene,
                initial_medium,
                settings,
                seed,
                pixel,
                sample,
                &mut volume_events,
                &mut free_flight_draws,
            ),
        )?;
        Ok(CurvedSequentialRayleighPixelSample {
            pixel,
            volume_events,
            free_flight_draws,
        })
    }

    /// Observe the same resolved lens and native scattering paths with the
    /// authored sum-of-QE proposal, including supported infrared bands.
    /// Lens transmission, conditional Beer loss, Mueller/proposal weights and
    /// geometry each retain their original owner and are integrated once.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        CurvedSequentialRayleighPixelSample<CurvedSequentialSensorPixelSample<C>>,
        CurvedRasterError,
    > {
        self.sample_sequential_rayleigh_sensor(
            cx,
            scene,
            sensor,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            None,
        )
    }

    /// The QE observation with its passive analyzer at the sensor. Molecular
    /// scattering retains the complete lens/analyzer Stokes correlations.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_analyzed<const C: usize>(
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
    ) -> Result<
        CurvedSequentialRayleighPixelSample<CurvedSequentialSensorPixelSample<C>>,
        CurvedRasterError,
    > {
        self.sample_sequential_rayleigh_sensor(
            cx,
            scene,
            sensor,
            initial_medium,
            settings,
            seed,
            pixel,
            sample,
            Some(AnalyzerResponse::Constant(analyzer)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_sequential_rayleigh_sensor<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<AnalyzerResponse<'_>>,
    ) -> Result<
        CurvedSequentialRayleighPixelSample<CurvedSequentialSensorPixelSample<C>>,
        CurvedRasterError,
    > {
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let pixel = self.sample_sequential_sensor_with_transport(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            analyzer,
            tracking_transport(
                scene,
                initial_medium,
                settings,
                seed,
                pixel,
                sample,
                &mut volume_events,
                &mut free_flight_draws,
            ),
        )?;
        Ok(CurvedSequentialRayleighPixelSample {
            pixel,
            volume_events,
            free_flight_draws,
        })
    }
}

// Both observers use the original independent surface and free-flight counter
// domains, including the fourth variate required by heterogeneous tracking.
// Returning the closure retains static dispatch and allocates no per-pixel data.
#[allow(clippy::too_many_arguments)]
fn tracking_transport<'a>(
    scene: &'a CurvedRayleighCameraScene<'_, '_, '_, '_>,
    initial_medium: usize,
    settings: CurvedRayleighTraceSettings,
    seed: u64,
    pixel: u32,
    sample: u32,
    volume_events: &'a mut [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: &'a mut [u32; CAMERA_WAVELENGTHS],
) -> impl FnMut(
    &Cx<'_>,
    &SequentialCameraBeam,
    usize,
) -> Result<CameraRadianceSample, CameraRadianceError>
+ 'a {
    let key = [seed as u32, (seed >> 32) as u32];
    move |cx, beam, lane| {
        let result = scene.sample_sequential_ray_with_tracking_samples(
            cx,
            beam,
            initial_medium,
            settings,
            |event| {
                let block = philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                core::array::from_fn(|i| unit(block[i]))
            },
            |flight| {
                philox4x32_10([pixel, sample, flight, 0x6372_7350 + lane as u32], key).map(unit)
            },
        )?;
        volume_events[lane] = result.volume_events();
        free_flight_draws[lane] = result.free_flight_draws();
        Ok(*result.radiance_sample())
    }
}
