//! Instantaneous visible/QE observations at the physical real-lens sensor.

use super::super::sensor::CurvedSequentialSensorPixelSample;
use super::super::{CameraRadianceSample, CurvedSequentialPixelSample};
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::camera::transient::admit_time;

mod exposure;

/// A resolved-lens observation with its original admitted sensor clock.
/// The contained visible or QE result retains optional wavelength paths and
/// the complete proposal measure, including physically clipped draws.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientSequentialPixelSample<P> {
    observation_time_s: Dd,
    pixel: P,
}

impl<P> TransientSequentialPixelSample<P> {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn pixel_sample(&self) -> &P {
        &self.pixel
    }
}

impl CurvedCameraRaster<'_> {
    /// Observe finite source histories at one physical sensor clock through
    /// the resolved lens. Each wavelength adds its real lens group delay to
    /// each original source encounter's native scene delay exactly once.
    /// A configured transient-scene analyzer is installed at the sensor before
    /// tracing the lens. Histories retain that sensor observation clock.
    ///
    /// Missing lens group laws refuse; ordinary steady entry points remain
    /// independent of them. Optical misses have no native path or history
    /// callback and retain their original zero-contribution draw. Even a
    /// completely clipped packet still admits the clock and cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sequential_transient_at_time(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedCameraScene<'_, '_, '_>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<TransientSequentialPixelSample<CurvedSequentialPixelSample>, TransientRasterError>
    {
        cx.checkpoint()
            .map_err(|error| TransientRasterError::Source(error.into()))?;
        let observation_time_s =
            admit_time(observation_time_s).map_err(TransientRasterError::Source)?;
        let mut source_error = None;
        let result = self.sample_sequential_pixel_with_timing(
            cx,
            seed,
            pixel,
            sample,
            scene.analyzer(),
            true,
            timed_transport(
                scene,
                initial_medium,
                observation_time_s,
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
        Ok(TransientSequentialPixelSample {
            observation_time_s,
            pixel: result?,
        })
    }

    /// The same sensor-clock observation integrated against authored QE and
    /// QE*lambda/(h*c), including supported infrared bands. Every channel
    /// shares the original four wavelength paths and source callbacks; lens
    /// geometry is applied once by the same steady detector integration owner.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_sequential_transient_at_time<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &TransientCurvedCameraScene<'_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<
        TransientSequentialPixelSample<CurvedSequentialSensorPixelSample<C>>,
        TransientRasterError,
    > {
        cx.checkpoint()
            .map_err(|error| TransientRasterError::Source(error.into()))?;
        let observation_time_s =
            admit_time(observation_time_s).map_err(TransientRasterError::Source)?;
        let mut source_error = None;
        let result = self.sample_sequential_sensor_with_timing(
            cx,
            sensor,
            seed,
            pixel,
            sample,
            scene.analyzer(),
            true,
            timed_transport(
                scene,
                initial_medium,
                observation_time_s,
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
        Ok(TransientSequentialPixelSample {
            observation_time_s,
            pixel: result?,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn timed_transport<'a>(
    scene: &'a TransientCurvedCameraScene<'_, '_, '_>,
    initial_medium: usize,
    observation_time_s: Dd,
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
    move |cx, beam, lane| match scene.sample_prepared_sequential_ray_with_samples_at_time(
        cx,
        beam,
        initial_medium,
        observation_time_s,
        settings,
        |event| {
            let block = philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
            core::array::from_fn(|i| unit(block[i]))
        },
    ) {
        Ok(timed) => Ok(*timed.radiance_sample()),
        Err(error) => {
            *source_error = Some(error);
            Err(CameraRadianceError::NumericalRange)
        }
    }
}

#[cfg(test)]
mod tests;
