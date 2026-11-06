//! Spectral pixels with finite-area next-event and emitter-hit MIS.
//!
//! Each wavelength follows its own camera walk and refractive source images.
//! Raster, lens, wavelength and continuation draws retain the original streams;
//! area sampling has a separate event/binding/lane Philox domain. One prepared
//! worker workspace is reused sequentially, and any failed lane refuses the
//! complete pixel. XYZ integration remains owned by the original raster.

use super::super::area_lights::camera::{
    AreaLightCameraError, AreaLightCameraSample, AreaLightCameraSettings, AreaLightCameraWorkspace,
    CurvedAreaLightCameraScene, MAX_AREA_LIGHT_BINDINGS,
};
use super::sensor::CurvedSensorPixelSample;
use super::{
    CAMERA_WAVELENGTHS, CameraRadianceError, CameraRadianceSample, CurvedCameraRaster,
    CurvedPixelSample, CurvedRasterError, unit,
};
use crate::charts::Ray;
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// Spectral pixels with native molecular flights and finite-area lighting.
#[cfg(feature = "volumes")]
pub mod participating;

/// Keep the original physical-camera or source-transport failure with its lane.
#[derive(Debug)]
pub enum AreaLightRasterError<E = AreaLightCameraError> {
    Raster(CurvedRasterError),
    Transport {
        lane: usize,
        error: E,
    },
}

impl<E: core::fmt::Debug> core::fmt::Display for AreaLightRasterError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native area-light pixel: {self:?}")
    }
}
impl<E: core::fmt::Debug> core::error::Error for AreaLightRasterError<E> {}

/// A complete spectral pixel with per-lane NEE contributions and source work.
/// Total radiance also contains the complementary emitter-hit contribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaLightPixelSample<P = CurvedPixelSample> {
    pixel: P,
    area_radiance_w_m2_sr_nm: [f64; CAMERA_WAVELENGTHS],
    area_samples: [u64; CAMERA_WAVELENGTHS],
    images: [u64; CAMERA_WAVELENGTHS],
    trace_evaluations: [u64; CAMERA_WAVELENGTHS],
    field_evaluations: [u64; CAMERA_WAVELENGTHS],
    query_work: [u64; CAMERA_WAVELENGTHS],
}

impl<P> AreaLightPixelSample<P> {
    pub const fn pixel_sample(&self) -> &P {
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

    /// Change only the detector observation while retaining every original
    /// lane's source radiance and work; used by finite-shutter unit conversion.
    pub(super) fn try_map_pixel<Q, E>(
        self,
        map: impl FnOnce(P) -> Result<Q, E>,
    ) -> Result<AreaLightPixelSample<Q>, E> {
        Ok(AreaLightPixelSample {
            pixel: map(self.pixel)?,
            area_radiance_w_m2_sr_nm: self.area_radiance_w_m2_sr_nm,
            area_samples: self.area_samples,
            images: self.images,
            trace_evaluations: self.trace_evaluations,
            field_evaluations: self.field_evaluations,
            query_work: self.query_work,
        })
    }
}

/// Typed counter address shared by steady and transient consumers. Keeping
/// draws as methods permits the original inner-loop closures to inline.
#[derive(Clone, Copy)]
pub(super) struct AreaLaneSamples {
    key: [u32; 2],
    pixel: u32,
    sample: u32,
    lane: usize,
}

impl AreaLaneSamples {
    pub(super) fn surface(self, event: u32) -> [f64; 3] {
        let block = philox4x32_10(
            [self.pixel, self.sample, event, 0x6372_7340 + self.lane as u32], self.key);
        core::array::from_fn(|i| unit(block[i]))
    }
    pub(super) fn area(self, event: u32, binding: usize) -> [f64; 2] {
        let domain = 0x6372_6100 + (self.lane * MAX_AREA_LIGHT_BINDINGS + binding) as u32;
        let block = philox4x32_10([self.pixel, self.sample, event, domain], self.key);
        [unit(block[0]), unit(block[1])]
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
        self.sample_area_transport(
            scene,
            initial_medium,
            settings,
            workspace,
            seed,
            pixel,
            sample,
            |trace| self.sample_pixel_with_transport(cx, seed, pixel, sample, trace),
        )
    }

    /// Observe finite-area next-event and complementary emitter-hit transport
    /// with the detector's original wavelength proposal and physical channel
    /// weights. All channels share four native paths, including the lens's
    /// wavelength-dependent launch rays. No visible-band or display conversion
    /// is introduced. The returned rates precede etendue and exposure.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel_sensor_area_lights<const C: usize>(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedAreaLightCameraScene<'_, '_, '_, '_>,
        sensor: &QuantumEfficiencySensor<C>,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<AreaLightPixelSample<CurvedSensorPixelSample<C>>, AreaLightRasterError> {
        self.sample_area_transport(
            scene,
            initial_medium,
            settings,
            workspace,
            seed,
            pixel,
            sample,
            |trace| self.sample_sensor_with_transport(cx, sensor, seed, pixel, sample, trace),
        )
    }

    // Both observers use exactly the same native source/continuation streams,
    // MIS and source-work accounting. The callback runs once per spectral lane;
    // it does not enter any geometry, image-solve or tracking inner loop.
    #[allow(clippy::too_many_arguments)]
    fn sample_area_transport<P>(
        &self,
        scene: &CurvedAreaLightCameraScene<'_, '_, '_, '_>,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        pixel: u32,
        sample: u32,
        integrate: impl FnOnce(
            &mut dyn FnMut(
                &Cx<'_>,
                Ray,
                f64,
                usize,
            ) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<AreaLightPixelSample<P>, AreaLightRasterError> {
        self.sample_area_with_transport(seed, pixel, sample,
            |cx, ray, wavelength, samples| {
                scene.sample_ray_with_samples(cx, ray, wavelength, initial_medium, settings,
                    workspace, |event| samples.surface(event), |event, binding| samples.area(event, binding))
            }, integrate)
    }

    /// One owner for the native continuation/area counter domains and complete
    /// per-lane work. A temporal adapter supplies only its source observation.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn sample_area_with_transport<P, E>(
        &self,
        seed: u64,
        pixel: u32,
        sample: u32,
        mut trace: impl FnMut(&Cx<'_>, Ray, f64, AreaLaneSamples) -> Result<AreaLightCameraSample, E>,
        integrate: impl FnOnce(
            &mut dyn FnMut(
                &Cx<'_>, Ray, f64, usize,
            ) -> Result<CameraRadianceSample, CameraRadianceError>,
        ) -> Result<P, CurvedRasterError>,
    ) -> Result<AreaLightPixelSample<P>, AreaLightRasterError<E>> {
        let key = [seed as u32, (seed >> 32) as u32];
        let mut area_radiance_w_m2_sr_nm = [0.0; CAMERA_WAVELENGTHS];
        let mut area_samples = [0; CAMERA_WAVELENGTHS];
        let mut images = [0; CAMERA_WAVELENGTHS];
        let mut trace_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut field_evaluations = [0; CAMERA_WAVELENGTHS];
        let mut query_work = [0; CAMERA_WAVELENGTHS];
        let mut source_failure = None;
        let result = integrate(&mut |cx, ray, wavelength, lane| {
            let traced = match trace(
                cx,
                ray,
                wavelength,
                AreaLaneSamples { key, pixel, sample, lane },
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
        });
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
