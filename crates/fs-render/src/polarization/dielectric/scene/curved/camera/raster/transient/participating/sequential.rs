//! Sensor-clock visible/QE observations through a real lens and molecular media.

use super::super::super::{CurvedSequentialPixelSample, sensor::CurvedSequentialSensorPixelSample};
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::camera::CameraRadianceSample;
use crate::polarization::dielectric::scene::curved::camera::transient::{
    TransientCameraError, admit_time,
};

mod exposure;

/// A complete sensor-clock sequential-lens observation with original work.
/// Optional optical lanes retain their proposal mass. Clipped lanes have no
/// native path, no history queries, and zero volume/flight counts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientSequentialRayleighPixelSample<P> {
    observation_time_s: Dd,
    pixel: P,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
}

impl<P> TransientSequentialRayleighPixelSample<P> {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn pixel_sample(&self) -> &P {
        &self.pixel
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }
}

impl CurvedCameraRaster<'_> {
    /// Integrate source histories into XYZ after native refraction and molecular
    /// Mueller scattering. Each source uses its own lens-plus-scene group delay
    /// from the physical sensor clock. A configured analyzer is installed at
    /// the sensor before the lens; geometric weights are integrated once.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_rayleigh_transient_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientSequentialRayleighPixelSample<CurvedSequentialPixelSample>,
        TransientRasterError,
    > {
        cx.checkpoint()
            .map_err(|error| TransientRasterError::Source(error.into()))?;
        let observation_time_s =
            admit_time(observation_time_s).map_err(TransientRasterError::Source)?;
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
            timed_tracking_transport(
                scene,
                initial_medium,
                observation_time_s,
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
        Ok(TransientSequentialRayleighPixelSample {
            observation_time_s,
            pixel: result?,
            volume_events,
            free_flight_draws,
        })
    }

    /// The same four native paths integrated against authored QE and photon
    /// response, including supported infrared bands. Channels share histories,
    /// collisions and all random addresses. Any source/transport failure aborts
    /// the complete observation; surviving wavelengths are never renormalized.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_rayleigh_transient_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedRayleighTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientSequentialRayleighPixelSample<CurvedSequentialSensorPixelSample<C>>,
        TransientRasterError,
    > {
        cx.checkpoint()
            .map_err(|error| TransientRasterError::Source(error.into()))?;
        let observation_time_s =
            admit_time(observation_time_s).map_err(TransientRasterError::Source)?;
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
            timed_tracking_transport(
                scene,
                initial_medium,
                observation_time_s,
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
        Ok(TransientSequentialRayleighPixelSample {
            observation_time_s,
            pixel: result?,
            volume_events,
            free_flight_draws,
        })
    }
}

// Both detector integrations use the existing independent surface/free-flight
// domains, including the original heterogeneous acceptance draw. Static
// dispatch retains the native worker without per-pixel allocation or retracing.
#[allow(clippy::too_many_arguments)]
fn timed_tracking_transport<'a>(
    scene: &'a TransientCurvedRayleighCameraScene<'_, '_, '_, '_, '_>,
    initial_medium: usize,
    observation_time_s: Dd,
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
        let result = scene.sample_prepared_sequential_ray_with_tracking_samples_at_time(
            cx,
            beam,
            initial_medium,
            observation_time_s,
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
            Ok(timed) => {
                volume_events[lane] = timed.volume_events();
                free_flight_draws[lane] = timed.free_flight_draws();
                Ok(*timed.radiance_sample())
            }
            Err(error) => {
                *source_error = Some(error);
                Err(CameraRadianceError::NumericalRange)
            }
        }
    }
}
