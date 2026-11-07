//! Finite shutter measurements through the original resolved lens.

use super::super::exposure::{
    ExposurePixelSample, ExposureSensorPixelSample, finish_pixel, finish_sensor_rates,
};
use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;

impl CurvedCameraRaster<'_> {
    /// Integrate original source histories over a rectangular shutter at the
    /// physical sensor. Each wavelength uses one timed lens traversal and one
    /// native scene path; each emitting hit supplies its own retarded interval
    /// average. No time quadrature, source-time sampling or extra RNG is used.
    ///
    /// A configured analyzer is installed at the sensor before the lens.
    /// Clipped wavelengths retain zero contribution and their original proposal
    /// mass. Supplying a different window per pixel explicitly models row or
    /// pixel exposure timing; no readout timing or moving geometry is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_during_exposure(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedCameraScene<'_, '_, '_>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<ExposurePixelSample<CurvedSequentialPixelSample>, TransientRasterError> {
        cx.checkpoint()
            .map_err(|error| TransientRasterError::Source(error.into()))?;
        let mut source_error = None;
        let result = self.sample_sequential_pixel_with_timing(
            cx,
            seed,
            pixel,
            sample,
            scene.analyzer(),
            true,
            exposure_transport(
                scene,
                initial_medium,
                window,
                settings,
                seed,
                pixel,
                sample,
                &mut source_error,
            ),
        );
        if let Some(error) = source_error {
            return Err(TransientRasterError::Source(error));
        }
        Ok(finish_pixel(window, result?))
    }

    /// The same real-lens exposure integrated against the original QE and
    /// QE*lambda/(h*c) responses, including supported infrared bands. Returned
    /// energy and expected photoelectrons use the explicit DD sensor duration;
    /// detector etendue, gain, saturation and noise remain separate.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_during_exposure<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedCameraScene<'_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        window: ExposureWindow,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        ExposureSensorPixelSample<C, CurvedSequentialSensorPixelSample<C>>,
        TransientRasterError,
    > {
        cx.checkpoint()
            .map_err(|error| TransientRasterError::Source(error.into()))?;
        let mut source_error = None;
        let result = self.sample_sequential_sensor_with_timing(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            scene.analyzer(),
            true,
            exposure_transport(
                scene,
                initial_medium,
                window,
                settings,
                seed,
                pixel,
                sample,
                &mut source_error,
            ),
        );
        if let Some(error) = source_error {
            return Err(TransientRasterError::Source(error));
        }
        let mean = result?;
        let qe_rate = *mean.qe_weighted_radiance_w_m2_sr();
        let electron_rate = *mean.photoelectron_rate_per_m2_sr();
        finish_sensor_rates(window, mean, qe_rate, electron_rate)
    }
}

#[allow(clippy::too_many_arguments)]
fn exposure_transport<'a>(
    scene: &'a TransientCurvedCameraScene<'_, '_, '_>,
    initial_medium: usize,
    window: ExposureWindow,
    settings: CurvedTraceSettings,
    seed: u64,
    pixel: u32,
    sample: u32,
    source_error: &'a mut Option<TransientCameraError>,
) -> impl FnMut(
    &Cx<'_>,
    &SequentialCameraBeam,
    usize,
) -> Result<CameraRadianceSample, CameraRadianceError>
+ 'a {
    let key = [seed as u32, (seed >> 32) as u32];
    move |cx, beam, lane| match scene.sample_prepared_sequential_ray_with_samples_during_exposure(
        cx,
        beam,
        initial_medium,
        window,
        settings,
        |event| {
            let block = philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
            core::array::from_fn(|i| unit(block[i]))
        },
    ) {
        Ok(exposed) => Ok(*exposed.mean_radiance_sample()),
        Err(error) => {
            *source_error = Some(error);
            Err(CameraRadianceError::NumericalRange)
        }
    }
}
