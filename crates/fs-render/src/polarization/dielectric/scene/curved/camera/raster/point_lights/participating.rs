//! Calibrated refractive point fixtures in molecularly scattering pixels.
//!
//! Camera flights retain the original molecular path's probabilities and random
//! streams. Each source connection independently estimates scattering survival
//! on its selected curved optical branch, at diffuse contacts and explicitly
//! bound molecular vertices. Its samples cannot consume or replace camera
//! draws. All four wavelengths retain separate source work and budgets.

use super::{
    CAMERA_WAVELENGTHS, CameraRadianceError, CurvedCameraRaster, CurvedPixelSample,
    PointLightPixelSample, PointLightRasterError, unit,
};
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::{
    CurvedRayleighPointLightCameraScene, RayleighPointLightCameraSettings,
};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// Complete spectral pixel with distinct camera and source scattering work.
/// The source counters include candidates with an exactly zero ratio weight;
/// these are legitimate samples, not failed or retried connections.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighPointLightPixelSample {
    point_lights: PointLightPixelSample,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
    source_free_flight_draws: [u32; CAMERA_WAVELENGTHS],
    source_candidates: [u32; CAMERA_WAVELENGTHS],
}

impl RayleighPointLightPixelSample {
    pub const fn point_light_sample(&self) -> &PointLightPixelSample {
        &self.point_lights
    }
    pub const fn pixel_sample(&self) -> &CurvedPixelSample {
        self.point_lights.pixel_sample()
    }
    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }
    pub const fn free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.free_flight_draws
    }
    pub const fn source_free_flight_draws(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.source_free_flight_draws
    }
    pub const fn source_candidates(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.source_candidates
    }
}

impl CurvedCameraRaster<'_> {
    /// Sample physical lens/raster coordinates and native spectral quadrature,
    /// adding selected point fixtures with their source scattering survival.
    /// The camera's own scattering survival is already in its flight sampling;
    /// no additional camera extinction factor is applied. Any unresolved lane
    /// refuses the complete pixel, retaining its original transport error.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_point_lights(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighPointLightCameraScene<'_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighPointLightPixelSample, PointLightRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut fixture_radiance_w_m2_sr_nm = [0.0; CAMERA_WAVELENGTHS];
        let mut connections = [0; CAMERA_WAVELENGTHS];
        let mut trace_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut field_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut query_work = [0; CAMERA_WAVELENGTHS];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_candidates = [0; CAMERA_WAVELENGTHS];
        let mut failure = None;
        let pixel_result = self.sample_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let traced = match scene.sample_ray_with_tracking_samples(
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
                    |connection, flight| {
                        // Separate key domain leaves all four counter words
                        // available for the full pixel/sample/connection/draw
                        // identity. No truncation or ordinal packing is used.
                        let source_key = [key[0], key[1] ^ (0x6372_7360 + lane as u32)];
                        unit(philox4x32_10([pixel, sample, connection, flight], source_key)[0])
                    },
                ) {
                    Ok(traced) => traced,
                    Err(error) => {
                        failure = Some((lane, error));
                        return Err(CameraRadianceError::NumericalRange);
                    }
                };
                fixture_radiance_w_m2_sr_nm[lane] = traced.fixture_radiance_w_m2_sr_nm();
                connections[lane] = u64::from(traced.connections());
                trace_evaluations[lane] = u64::from(traced.trace_evaluations());
                field_evaluations[lane] = traced.field_evaluations();
                query_work[lane] = traced.query_work();
                volume_events[lane] = traced.volume_events();
                free_flight_draws[lane] = traced.free_flight_draws();
                source_free_flight_draws[lane] = traced.source_free_flight_draws();
                source_candidates[lane] = traced.source_candidates();
                Ok(*traced.radiance_sample())
            },
        );
        if let Some((lane, error)) = failure {
            return Err(PointLightRasterError::Transport { lane, error });
        }
        Ok(RayleighPointLightPixelSample {
            point_lights: PointLightPixelSample {
                pixel: pixel_result.map_err(PointLightRasterError::Raster)?,
                fixture_radiance_w_m2_sr_nm,
                connections,
                trace_evaluations,
                field_evaluations,
                query_work,
            },
            volume_events,
            free_flight_draws,
            source_free_flight_draws,
            source_candidates,
        })
    }
}

#[cfg(test)]
mod tests;
