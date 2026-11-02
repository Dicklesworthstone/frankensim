//! Physical-camera pixels with calibrated point-fixture connections.
//!
//! Every wavelength independently connects the declared optical branch at each
//! matching diffuse encounter. Raster, lens, wavelength and surface variates
//! retain the original streams. Native color integration consumes the complete
//! spectral radiance; failed source connections abort the whole pixel sample.
//! The selected branches do not establish completeness over other optical images.

use super::super::point_lights::{
    CurvedPointLightCameraScene, PointLightCameraError, PointLightCameraSettings,
};
use super::{
    CAMERA_WAVELENGTHS, CameraRadianceError, CurvedCameraRaster, CurvedPixelSample,
    CurvedRasterError, unit,
};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// A physical-camera failure or the original source-connection failure in one
/// spectral lane. No successfully traced prefix is published as a partial pixel.
#[derive(Debug)]
pub enum PointLightRasterError {
    Raster(CurvedRasterError),
    Transport {
        lane: usize,
        error: PointLightCameraError,
    },
}

impl core::fmt::Display for PointLightRasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved point-fixture pixel: {self:?}")
    }
}
impl core::error::Error for PointLightRasterError {}

/// Complete native pixel plus the selected fixtures' contributions and actual
/// connection work in each independently traced wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointLightPixelSample {
    pixel: CurvedPixelSample,
    fixture_radiance_w_m2_sr_nm: [f64; CAMERA_WAVELENGTHS],
    connections: [u64; CAMERA_WAVELENGTHS],
    trace_evaluations: [u64; CAMERA_WAVELENGTHS],
    field_evaluations: [u64; CAMERA_WAVELENGTHS],
    query_work: [u64; CAMERA_WAVELENGTHS],
}

impl PointLightPixelSample {
    pub const fn pixel_sample(&self) -> &CurvedPixelSample {
        &self.pixel
    }
    pub const fn fixture_radiance_w_m2_sr_nm(&self) -> &[f64; CAMERA_WAVELENGTHS] {
        &self.fixture_radiance_w_m2_sr_nm
    }
    pub const fn connections(&self) -> &[u64; CAMERA_WAVELENGTHS] {
        &self.connections
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
    /// Reuse the original physical camera, spectral stratification and XYZ
    /// quadrature. Each lane retains its cumulative camera/connection limits;
    /// total work is bounded by four such paths. Empty fixture bindings preserve
    /// the original surface-only pixel, including its random draws.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_point_lights(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPointLightCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<PointLightPixelSample, PointLightRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut fixture_radiance_w_m2_sr_nm = [0.0; CAMERA_WAVELENGTHS];
        let mut connections = [0; CAMERA_WAVELENGTHS];
        let mut trace_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut field_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut query_work = [0; CAMERA_WAVELENGTHS];
        let mut connection_failure = None;
        let pixel_result = self.sample_pixel_with_transport(
            cx,
            seed,
            pixel,
            sample,
            |cx, ray, wavelength, lane| {
                let traced = match scene.sample_ray_with_samples(
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
                ) {
                    Ok(traced) => traced,
                    Err(error) => {
                        // Keep the original error without adding a recursive
                        // camera-error variant or a boxed per-ray allocation.
                        connection_failure = Some((lane, error));
                        return Err(CameraRadianceError::NumericalRange);
                    }
                };
                fixture_radiance_w_m2_sr_nm[lane] = traced.fixture_radiance_w_m2_sr_nm();
                connections[lane] = u64::from(traced.connections());
                trace_evaluations[lane] = u64::from(traced.trace_evaluations());
                field_evaluations[lane] = traced.field_evaluations();
                query_work[lane] = traced.query_work();
                Ok(*traced.radiance_sample())
            },
        );
        if let Some((lane, error)) = connection_failure {
            return Err(PointLightRasterError::Transport { lane, error });
        }
        Ok(PointLightPixelSample {
            pixel: pixel_result.map_err(PointLightRasterError::Raster)?,
            fixture_radiance_w_m2_sr_nm,
            connections,
            trace_evaluations,
            field_evaluations,
            query_work,
        })
    }
}

#[cfg(test)]
mod tests;
