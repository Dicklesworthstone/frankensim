//! Spectral pixels with finite-area next-event and emitter-hit MIS.
//!
//! Each wavelength follows its own camera walk and refractive source images.
//! Raster, lens, wavelength and continuation draws retain the original streams;
//! area sampling has a separate event/binding/lane Philox domain. One prepared
//! worker workspace is reused sequentially, and any failed lane refuses the
//! complete pixel. XYZ integration remains owned by the original raster.

use super::super::area_lights::camera::{
    AreaLightCameraError, AreaLightCameraSettings, AreaLightCameraWorkspace,
    CurvedAreaLightCameraScene, MAX_AREA_LIGHT_BINDINGS,
};
use super::{
    CAMERA_WAVELENGTHS, CameraRadianceError, CurvedCameraRaster, CurvedPixelSample,
    CurvedRasterError, unit,
};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// Keep the original physical-camera or source-transport failure with its lane.
#[derive(Debug)]
pub enum AreaLightRasterError {
    Raster(CurvedRasterError),
    Transport {
        lane: usize,
        error: AreaLightCameraError,
    },
}

impl core::fmt::Display for AreaLightRasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native area-light pixel: {self:?}")
    }
}
impl core::error::Error for AreaLightRasterError {}

/// A complete spectral pixel with per-lane NEE contributions and source work.
/// Total radiance also contains the complementary emitter-hit contribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaLightPixelSample {
    pixel: CurvedPixelSample,
    area_radiance_w_m2_sr_nm: [f64; CAMERA_WAVELENGTHS],
    area_samples: [u64; CAMERA_WAVELENGTHS],
    images: [u64; CAMERA_WAVELENGTHS],
    trace_evaluations: [u64; CAMERA_WAVELENGTHS],
    field_evaluations: [u64; CAMERA_WAVELENGTHS],
    query_work: [u64; CAMERA_WAVELENGTHS],
}

impl AreaLightPixelSample {
    pub const fn pixel_sample(&self) -> &CurvedPixelSample {
        &self.pixel
    }

    pub const fn area_radiance_w_m2_sr_nm(&self) -> &[f64; CAMERA_WAVELENGTHS] {
        &self.area_radiance_w_m2_sr_nm
    }

    pub const fn area_samples(&self) -> &[u64; CAMERA_WAVELENGTHS] {
        &self.area_samples
    }

    pub const fn images(&self) -> &[u64; CAMERA_WAVELENGTHS] {
        &self.images
    }

    pub const fn trace_evaluations(&self) -> &[u64; CAMERA_WAVELENGTHS] {
        &self.trace_evaluations
    }

    pub const fn field_evaluations(&self) -> &[u64; CAMERA_WAVELENGTHS] {
        &self.field_evaluations
    }

    pub const fn query_work(&self) -> &[u64; CAMERA_WAVELENGTHS] {
        &self.query_work
    }
}

impl CurvedCameraRaster<'_> {
    /// Reuse original lens/raster sampling and spectral color integration.
    /// Each lane retains its complete camera and cumulative source budgets;
    /// total work is bounded by four such paths. The caller can reuse one
    /// workspace for successive samples on this worker.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_area_lights(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedAreaLightCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample, AreaLightRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut area_radiance_w_m2_sr_nm = [0.0; CAMERA_WAVELENGTHS];
        let mut area_samples = [0; CAMERA_WAVELENGTHS];
        let mut images = [0; CAMERA_WAVELENGTHS];
        let mut trace_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut field_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut query_work = [0; CAMERA_WAVELENGTHS];
        let mut source_failure = None;
        let result = self.sample_pixel_with_transport(
            cx, seed, pixel, sample,
            |cx, ray, wavelength, lane| {
                let traced = match scene.sample_ray_with_samples(
                    cx, ray, wavelength, initial_medium, settings, workspace,
                    |event| {
                        let block = philox4x32_10(
                            [pixel, sample, event, 0x6372_7340 + lane as u32], key,
                        );
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |event, binding| {
                        let domain = 0x6372_6100
                            + (lane * MAX_AREA_LIGHT_BINDINGS + binding) as u32;
                        let block = philox4x32_10([pixel, sample, event, domain], key);
                        [unit(block[0]), unit(block[1])]
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
                Ok(*traced.radiance_sample())
            },
        );
        if let Some((lane, error)) = source_failure {
            return Err(AreaLightRasterError::Transport { lane, error });
        }
        Ok(AreaLightPixelSample {
            pixel: result.map_err(AreaLightRasterError::Raster)?,
            area_radiance_w_m2_sr_nm,
            area_samples,
            images,
            trace_evaluations,
            field_evaluations,
            query_work,
        })
    }
}

#[cfg(test)]
mod tests;
