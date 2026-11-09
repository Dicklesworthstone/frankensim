//! Physical-camera pixels with calibrated point-fixture connections.
//!
//! Every wavelength independently connects the declared optical branch at each
//! matching diffuse encounter. Raster, lens, wavelength and surface variates
//! retain the original streams. Native color integration consumes the complete
//! spectral radiance; failed source connections abort the whole pixel sample.
//! The selected branches do not establish completeness over other optical images.

use super::super::point_lights::{
    CurvedPointLightCameraScene, PointLightCameraError, PointLightCameraSample,
    PointLightCameraSettings,
};
use super::{
    CAMERA_WAVELENGTHS, CameraRadianceError, CameraRadianceSample, CurvedCameraRaster,
    CurvedPixelSample, CurvedRasterError, unit,
};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// Molecularly scattering camera and source paths with the same spectral pixels.
#[cfg(feature = "volumes")]
pub mod participating;

/// A physical-camera failure or the original source-connection failure in one
/// spectral lane. No successfully traced prefix is published as a partial pixel.
#[derive(Debug)]
pub enum PointLightRasterError<E = PointLightCameraError> {
    Raster(CurvedRasterError),
    Transport { lane: usize, error: E },
}

impl<E: core::fmt::Debug> core::fmt::Display for PointLightRasterError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved point-fixture pixel: {self:?}")
    }
}
impl<E: core::fmt::Debug> core::error::Error for PointLightRasterError<E> {}

/// Complete native pixel plus the selected fixtures' contributions and actual
/// connection work in each independently traced wavelength.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointLightPixelSample<P = CurvedPixelSample> {
    pixel: P,
    fixture_radiance_w_m2_sr_nm: [f64; CAMERA_WAVELENGTHS],
    connections: [u64; CAMERA_WAVELENGTHS],
    trace_evaluations: [u64; CAMERA_WAVELENGTHS],
    field_evaluations: [u64; CAMERA_WAVELENGTHS],
    query_work: [u64; CAMERA_WAVELENGTHS],
}

impl<P> PointLightPixelSample<P> {
    pub const fn pixel_sample(&self) -> &P {
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

    pub(super) fn try_map_pixel<Q, E>(
        self,
        map: impl FnOnce(P) -> Result<Q, E>,
    ) -> Result<PointLightPixelSample<Q>, E> {
        Ok(PointLightPixelSample {
            pixel: map(self.pixel)?,
            fixture_radiance_w_m2_sr_nm: self.fixture_radiance_w_m2_sr_nm,
            connections: self.connections,
            trace_evaluations: self.trace_evaluations,
            field_evaluations: self.field_evaluations,
            query_work: self.query_work,
        })
    }
}

/// The original camera and independent point-connection counter addresses.
#[derive(Clone, Copy)]
pub(super) struct PointLaneSamples {
    key: [u32; 2],
    pixel: u32,
    sample: u32,
    lane: usize,
}
impl PointLaneSamples {
    pub(super) fn surface(self, event: u32) -> [f64; 3] {
        let block = philox4x32_10(
            [
                self.pixel,
                self.sample,
                event,
                0x6372_7340 + self.lane as u32,
            ],
            self.key,
        );
        core::array::from_fn(|i| unit(block[i]))
    }
    #[cfg(feature = "volumes")]
    pub(super) fn flight(self, flight: u32) -> [f64; 4] {
        philox4x32_10(
            [
                self.pixel,
                self.sample,
                flight,
                0x6372_7350 + self.lane as u32,
            ],
            self.key,
        )
        .map(unit)
    }
    #[cfg(feature = "volumes")]
    pub(super) fn source(self, connection: u32, flight: u32) -> f64 {
        let key = [self.key[0], self.key[1] ^ (0x6372_7360 + self.lane as u32)];
        unit(philox4x32_10([self.pixel, self.sample, connection, flight], key)[0])
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
        self.sample_point_with_launch(
            seed,
            pixel,
            sample,
            |cx, (ray, wavelength), samples| {
                scene.sample_ray_with_samples(
                    cx,
                    ray,
                    wavelength,
                    initial_medium,
                    settings,
                    |event| samples.surface(event),
                )
            },
            |trace| {
                self.sample_pixel_with_transport(
                    cx,
                    seed,
                    pixel,
                    sample,
                    |cx, ray, wavelength, lane| trace(cx, (ray, wavelength), lane),
                )
            },
        )
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn sample_point_with_launch<P, E, L>(
        &self,
        seed: u64,
        pixel: u32,
        sample: u32,
        mut trace: impl FnMut(&Cx<'_>, L, PointLaneSamples) -> Result<PointLightCameraSample, E>,
        integrate: impl FnOnce(
            &mut dyn FnMut(&Cx<'_>, L, usize) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<PointLightPixelSample<P>, PointLightRasterError<E>> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut fixture_radiance_w_m2_sr_nm = [0.0; CAMERA_WAVELENGTHS];
        let mut connections = [0; CAMERA_WAVELENGTHS];
        let mut trace_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut field_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut query_work = [0; CAMERA_WAVELENGTHS];
        let mut connection_failure = None;
        let pixel_result = integrate(&mut |cx, launch, lane| {
            let traced = match trace(
                cx,
                launch,
                PointLaneSamples {
                    key,
                    pixel,
                    sample,
                    lane,
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
        });
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
