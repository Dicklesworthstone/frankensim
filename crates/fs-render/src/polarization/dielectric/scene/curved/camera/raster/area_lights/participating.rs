//! Spectral pixels with finite-area lighting along native molecular paths.
//!
//! Camera, lens, wavelength, surface and free-flight streams retain their
//! original identities. Surface-area, molecular-area and source-survival draws
//! occupy separate counter domains. One prepared workspace is reused across
//! the four wavelength lanes; a failed lane refuses the complete pixel.

use super::{AreaLightPixelSample, AreaLightRasterError};
use super::super::{
    CAMERA_WAVELENGTHS, CameraRadianceError, CurvedCameraRaster, CurvedPixelSample, unit,
};
use crate::polarization::dielectric::scene::curved::camera::area_lights::camera::MAX_AREA_LIGHT_BINDINGS;
use crate::polarization::dielectric::scene::curved::camera::area_lights::participating::{
    CurvedRayleighAreaLightCameraScene, RayleighAreaLightCameraSettings,
    RayleighAreaLightCameraWorkspace,
};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// A complete original physical-camera pixel and its per-wavelength source
/// and molecular work. The contained area result includes both weighted
/// next-event terms and complementary emitter hits in its total radiance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighAreaLightPixelSample {
    area: AreaLightPixelSample,
    volume_events: [u32; CAMERA_WAVELENGTHS],
    free_flight_draws: [u32; CAMERA_WAVELENGTHS],
    source_free_flight_draws: [u32; CAMERA_WAVELENGTHS],
    source_candidates: [u32; CAMERA_WAVELENGTHS],
}

impl RayleighAreaLightPixelSample {
    pub const fn pixel_sample(&self) -> &CurvedPixelSample {
        self.area.pixel_sample()
    }

    pub const fn area_pixel_sample(&self) -> &AreaLightPixelSample {
        &self.area
    }

    pub const fn volume_events(&self) -> &[u32; CAMERA_WAVELENGTHS] {
        &self.volume_events
    }

    /// Includes native null candidates and flights interrupted by surfaces.
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
    /// Render each wavelength through its original refractive molecular
    /// scene, including finite-area connections at diffuse and real molecular
    /// vertices. Surface-area streams agree with the surface-only adapter.
    /// Source tracking is independent of both native camera flight acceptance
    /// and source-area selection. Each lane enforces its own cumulative bounds.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_rayleigh_area_lights(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedRayleighAreaLightCameraScene<'_, '_, '_, '_, '_, '_>,
        initial_medium: usize,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<RayleighAreaLightPixelSample, AreaLightRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut area_radiance_w_m2_sr_nm = [0.0; CAMERA_WAVELENGTHS];
        let mut area_samples = [0; CAMERA_WAVELENGTHS];
        let mut images = [0; CAMERA_WAVELENGTHS];
        let mut trace_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut field_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut query_work = [0; CAMERA_WAVELENGTHS];
        let mut volume_events = [0; CAMERA_WAVELENGTHS];
        let mut free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_free_flight_draws = [0; CAMERA_WAVELENGTHS];
        let mut source_candidates = [0; CAMERA_WAVELENGTHS];
        let mut source_failure = None;
        let result = self.sample_pixel_with_transport(
            cx, seed, pixel, sample,
            |cx, ray, wavelength, lane| {
                let traced = match scene.sample_ray_with_tracking_samples(
                    cx, ray, wavelength, initial_medium, settings, workspace,
                    |event| {
                        let block = philox4x32_10(
                            [pixel, sample, event, 0x6372_7340 + lane as u32], key,
                        );
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |flight| {
                        philox4x32_10(
                            [pixel, sample, flight, 0x6372_7350 + lane as u32], key,
                        ).map(unit)
                    },
                    |event, binding| {
                        let domain = 0x6372_6100
                            + (lane * MAX_AREA_LIGHT_BINDINGS + binding) as u32;
                        let block = philox4x32_10([pixel, sample, event, domain], key);
                        [unit(block[0]), unit(block[1])]
                    },
                    |collision, binding| {
                        let domain = 0x6372_6500
                            + (lane * MAX_AREA_LIGHT_BINDINGS + binding) as u32;
                        let block = philox4x32_10([pixel, sample, collision, domain], key);
                        [unit(block[0]), unit(block[1])]
                    },
                    |area, image, flight| {
                        let source_key = [
                            key[0] ^ image,
                            key[1] ^ (0x6372_6200 + lane as u32),
                        ];
                        unit(philox4x32_10([pixel, sample, area, flight], source_key)[0])
                    },
                ) {
                    Ok(traced) => traced,
                    Err(error) => {
                        source_failure = Some((lane, error));
                        return Err(CameraRadianceError::NumericalRange);
                    }
                };
                area_radiance_w_m2_sr_nm[lane] = traced.area_radiance_w_m2_sr_nm();
                area_samples[lane] = u64::from(traced.area_samples());
                images[lane] = u64::from(traced.images());
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
        if let Some((lane, error)) = source_failure {
            return Err(AreaLightRasterError::Transport { lane, error });
        }
        Ok(RayleighAreaLightPixelSample {
            area: AreaLightPixelSample {
                pixel: result.map_err(AreaLightRasterError::Raster)?,
                area_radiance_w_m2_sr_nm,
                area_samples,
                images,
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
