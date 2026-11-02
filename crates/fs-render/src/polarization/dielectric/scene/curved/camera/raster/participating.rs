//! Physical camera pixels through refractive, molecularly scattering media.
//!
//! The same native camera ray, four independent wavelengths and CIE integration
//! are used as for surface-only pixels. Surface random-number domains remain
//! unchanged. Each lane has a separate free-flight counter domain, with an
//! independent fourth variate for heterogeneous delta-tracking acceptance.

use super::{CAMERA_WAVELENGTHS, CurvedCameraRaster, CurvedPixelSample, CurvedRasterError, unit};
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighTraceSettings;
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// One complete spectral pixel sample with per-wavelength native volume work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighPixelSample {
    pixel: CurvedPixelSample,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
}

impl RayleighPixelSample {
    pub const fn pixel_sample(&self) -> &CurvedPixelSample {
        &self.pixel
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    /// Counts include rejected delta-tracking candidates and flights interrupted
    /// by a surface; these are distinct from real scattering collisions.
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }
}

impl CurvedCameraRaster<'_> {
    /// Sample the original physical camera through the admitted molecular
    /// scene. Every lane retains the original surface and volume budgets;
    /// unresolved paths and cancellation abort the complete pixel sample.
    /// Zero scattering preserves the surface-only pixel's complete result.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
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
                let traced = scene.sample_ray_with_tracking_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
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
        Ok(RayleighPixelSample {
            pixel: result,
            volume_events,
            free_flight_draws,
        })
    }
}

#[cfg(test)]
mod tests;
