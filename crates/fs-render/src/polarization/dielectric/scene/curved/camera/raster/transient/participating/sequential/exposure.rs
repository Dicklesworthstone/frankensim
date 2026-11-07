//! Physical sensor shutters through the resolved lens and molecular medium.

use super::super::super::exposure::{finish_pixel, finish_sensor_rates};
use super::super::exposure::{ExposureRayleighPixelSample, ExposureRayleighSensorPixelSample};
use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;

impl CurvedCameraRaster<'_> {
    /// Integrate each source's original interval history along one actual
    /// refracted/scattered path per wavelength. The physical sensor shutter is
    /// retarded by lens transit and that source's own native scene delay.
    /// A configured analyzer is installed before the lens. Clipped wavelengths
    /// retain their original proposal mass, zero contribution and zero work.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<ExposureRayleighPixelSample<CurvedSequentialPixelSample>, TransientRasterError>
    {
        cx.checkpoint()
            .map_err(|error| TransientRasterError::Source(error.into()))?;
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_error = None;
        let result = self.sample_sequential_pixel_with_timing(
            cx,
            seed,
            pixel,
            sample,
            scene.transient_camera().analyzer(),
            true,
            exposure_tracking_transport(
                scene,
                initial_medium,
                window,
                settings,
                seed,
                pixel,
                sample,
                &mut volume_events,
                &mut free_flight_draws,
                &mut source_error,
            ),
        );
        if let Some(error) = source_error {
            return Err(TransientRasterError::Source(error));
        }
        Ok(ExposureRayleighPixelSample::with_work(
            finish_pixel(window, result?),
            volume_events,
            free_flight_draws,
        ))
    }

    /// The same four original molecular paths and source integrals feed every
    /// QE channel, including supported infrared bands. Integrated energy and
    /// expected photoelectrons use the admitted DD sensor duration once.
    /// Detector etendue/electronics and time-dependent volume laws are separate.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        ExposureRayleighSensorPixelSample<C, CurvedSequentialSensorPixelSample<C>>,
        TransientRasterError,
    > {
        cx.checkpoint()
            .map_err(|error| TransientRasterError::Source(error.into()))?;
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_error = None;
        let result = self.sample_sequential_sensor_with_timing(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            scene.transient_camera().analyzer(),
            true,
            exposure_tracking_transport(
                scene,
                initial_medium,
                window,
                settings,
                seed,
                pixel,
                sample,
                &mut volume_events,
                &mut free_flight_draws,
                &mut source_error,
            ),
        );
        if let Some(error) = source_error {
            return Err(TransientRasterError::Source(error));
        }
        let mean = result?;
        let qe_rate = *mean.qe_weighted_radiance_w_m2_sr();
        let electron_rate = *mean.photoelectron_rate_per_m2_sr();
        Ok(ExposureRayleighSensorPixelSample::with_work(
            finish_sensor_rates(window, mean, qe_rate, electron_rate)?,
            volume_events,
            free_flight_draws,
        ))
    }
}

#[allow(clippy::too_many_arguments)]
fn exposure_tracking_transport<'a>(
    scene: &'a TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
    initial_medium: usize,
    window: ExposureWindow,
    settings: CurvedRayleighTraceSettings,
    seed: u64,
    pixel: u32,
    sample: u32,
    volume_events: &'a mut [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: &'a mut [u32; CAMERA_WAVELENGTHS],
    source_error: &'a mut Option<TransientCameraError>,
) -> impl FnMut(
    &Cx<'_>,
    &SequentialCameraBeam,
    usize,
) -> Result<CameraRadianceSample, CameraRadianceError>
+ 'a {
    let key = [seed as u32, (seed >> 32) as u32];
    move |cx, beam, lane| {
        let result = scene.sample_prepared_sequential_ray_with_tracking_samples_during_exposure(
            cx,
            beam,
            initial_medium,
            window,
            settings,
            |event| {
                let block = philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                core::array::from_fn(|i| unit(block[i]))
            },
            |flight| {
                philox4x32_10([pixel, sample, flight, 0x6372_7350 + lane as u32], key).map(unit)
            },
        );
        match result {
            Ok(exposure) => {
                volume_events[lane] = exposure.volume_events();
                free_flight_draws[lane] = exposure.free_flight_draws();
                Ok(*exposure.mean_radiance_sample())
            }
            Err(error) => {
                *source_error = Some(error);
                Err(CameraRadianceError::NumericalRange)
            }
        }
    }
}
