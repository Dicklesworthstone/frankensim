//! Physical-camera pixel samples with independent spectral curved paths.
//!
//! The existing pinhole/thin-lens camera supplies the geometric ray. Four
//! rotated, stratified hero wavelengths each take their own native refractive
//! path, so dispersion changes both visibility and radiance. No wavelength
//! borrows the hero's geometry. Raster, lens, wavelength and boundary variates
//! have separate counter-based streams; scheduling and progressive batches do
//! not change an absolute `(seed, pixel, sample)` result.
//!
//! XYZ uses the renderer's equal-energy normalization, not sensor irradiance,
//! exposure, a photometric luminance conversion or a display transform. The
//! original SI spectral radiances and full native paths are retained alongside
//! it. The ideal lens estimates an aperture average, as in the scalar tracer.
//! Diffraction requires a source-depth model and is not admitted here.

use super::{CameraRadianceError, CameraRadianceSample, CurvedCameraScene};
use crate::camera::{CameraError, LensSample, PhysicalCamera};
use crate::charts::Ray;
use crate::hero_wavelength_packet;
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, cie_x, cie_y, cie_z, y_integral};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

/// Spectral pixel integration with native curved molecular scattering.
#[cfg(feature = "volumes")]
pub mod participating;

/// Spectral pixels with selected refractive point-fixture connections.
#[cfg(feature = "photometry")]
pub mod point_lights;

/// Spectral pixels with finite-area refractive light and emitter-hit MIS.
pub mod area_lights;

/// Each pixel sample independently traces this many wavelengths.
pub const CAMERA_WAVELENGTHS: usize = 4;

#[derive(Clone, Debug, PartialEq)]
pub enum CurvedRasterError {
    Camera(CameraError),
    Transport(CameraRadianceError),
    /// Zero/overflowing dimensions or an out-of-range row-major pixel.
    InvalidRaster,
    NumericalRange,
}

impl From<CameraError> for CurvedRasterError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}
impl From<CameraRadianceError> for CurvedRasterError {
    fn from(error: CameraRadianceError) -> Self {
        Self::Transport(error)
    }
}
impl core::fmt::Display for CurvedRasterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved camera raster: {self:?}")
    }
}
impl core::error::Error for CurvedRasterError {}

/// One complete sample, with no partially accepted spectral lanes on failure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedPixelSample {
    xyz: [f64; 3],
    ray: Ray,
    wavelengths_nm: [f64; CAMERA_WAVELENGTHS],
    spectral: [CameraRadianceSample; CAMERA_WAVELENGTHS],
}

impl CurvedPixelSample {
    /// Add this normalized XYZ contribution to a pixel sum, then divide by
    /// that pixel's completed sample count. No exposure or display transform.
    pub const fn xyz(&self) -> [f64; 3] {
        self.xyz
    }
    pub const fn camera_ray(&self) -> Ray {
        self.ray
    }
    pub const fn wavelengths_nm(&self) -> &[f64; CAMERA_WAVELENGTHS] {
        &self.wavelengths_nm
    }
    pub const fn spectral_samples(&self) -> &[CameraRadianceSample; CAMERA_WAVELENGTHS] {
        &self.spectral
    }
}

/// Immutable camera/raster preparation. Projection and color normalization are
/// evaluated once, and each pixel sample uses fixed stack storage. A caller
/// may distribute independent absolute sample identities over its executor.
pub struct CurvedCameraRaster<'a> {
    camera: &'a PhysicalCamera,
    width: u32,
    height: u32,
    pixels: u32,
    half_tan: f64,
    aspect: f64,
    xyz_scale: f64,
}

impl<'a> CurvedCameraRaster<'a> {
    pub fn try_new(
        camera: &'a PhysicalCamera,
        width: u32,
        height: u32,
    ) -> Result<Self, CurvedRasterError> {
        let pixels = width
            .checked_mul(height)
            .filter(|&count| count > 0)
            .ok_or(CurvedRasterError::InvalidRaster)?;
        if camera.has_diffraction() {
            return Err(CameraError::InvalidDiffraction.into());
        }
        let half_tan = camera.raster_vertical_half_tan();
        let aspect = f64::from(width) / f64::from(height);
        if !(half_tan > 0.0 && (half_tan * aspect).is_finite()) {
            return Err(CameraError::InvalidProjection.into());
        }
        Ok(Self {
            camera,
            width,
            height,
            pixels,
            half_tan,
            aspect,
            xyz_scale: (LAMBDA_MAX - LAMBDA_MIN) / CAMERA_WAVELENGTHS as f64 / y_integral(),
        })
    }

    pub const fn width(&self) -> u32 {
        self.width
    }
    pub const fn height(&self) -> u32 {
        self.height
    }

    /// Trace one absolute row-major pixel/sample identity through the original
    /// source scene. `initial_medium` applies at the sampled lens position;
    /// the caller must keep the whole aperture in that admitted medium.
    ///
    /// Each wavelength retains the supplied native path budgets, so total
    /// work is bounded by four such paths. Every numerical/cancellation/source
    /// error propagates: a limited path is never counted as a black sample.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_pixel(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedCameraScene<'_, '_>,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        pixel: u32,
        sample: u32,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        let key = [seed as u32, (seed >> 32) as u32];
        self.sample_pixel_with_transport(cx, seed, pixel, sample, |cx, ray, wavelength, lane| {
            scene.sample_ray_with_samples(cx, ray, wavelength, initial_medium, settings, |event| {
                let block = philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                core::array::from_fn(|i| unit(block[i]))
            })
        })
    }

    // Both surface-only and participating estimators share raster/lens draws,
    // wavelength stratification and color integration without changing bits.
    fn sample_pixel_with_transport(
        &self,
        cx: &Cx<'_>,
        seed: u64,
        pixel: u32,
        sample: u32,
        mut trace: impl FnMut(
            &Cx<'_>,
            Ray,
            f64,
            usize,
        ) -> Result<CameraRadianceSample, CameraRadianceError>,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        cx.checkpoint().map_err(CameraError::from)?;
        if pixel >= self.pixels {
            return Err(CurvedRasterError::InvalidRaster);
        }
        let render_cx = cx.with_stream_seed(seed);
        let key = [seed as u32, (seed >> 32) as u32];
        let raster = philox4x32_10([pixel, sample, 0, 0x6372_7331], key);
        let lens = philox4x32_10([pixel, sample, 0, 0x6372_7332], key);
        let x_tan =
            (2.0 * (f64::from(pixel % self.width) + unit(raster[0])) / f64::from(self.width) - 1.0)
                * self.aspect
                * self.half_tan;
        let y_tan = (1.0
            - 2.0 * (f64::from(pixel / self.width) + unit(raster[1])) / f64::from(self.height))
            * self.half_tan;
        let ray = self.camera.generate_ray_from_tangent_offsets(
            &render_cx,
            x_tan,
            y_tan,
            LensSample::try_new(unit(lens[0]), unit(lens[1]))?,
        )?;
        let wavelengths_nm = hero_wavelength_packet::<CAMERA_WAVELENGTHS>(
            LAMBDA_MIN + unit(raster[2]) * (LAMBDA_MAX - LAMBDA_MIN),
            LAMBDA_MIN,
            LAMBDA_MAX,
        );
        let mut trace_lane = |lane: usize| trace(&render_cx, ray, wavelengths_nm[lane], lane);
        // Fixed storage and explicit propagation of the first failed lane.
        let spectral = [
            trace_lane(0)?,
            trace_lane(1)?,
            trace_lane(2)?,
            trace_lane(3)?,
        ];
        let mut xyz = [0.0; 3];
        for (lane, &wavelength) in spectral.iter().zip(&wavelengths_nm) {
            let weight = lane.radiance_w_m2_sr_nm() * self.xyz_scale;
            for (sum, observer) in
                xyz.iter_mut()
                    .zip([cie_x(wavelength), cie_y(wavelength), cie_z(wavelength)])
            {
                *sum += weight * observer;
            }
        }
        if xyz.iter().any(|value| !value.is_finite()) {
            return Err(CurvedRasterError::NumericalRange);
        }
        render_cx.checkpoint().map_err(CameraError::from)?;
        Ok(CurvedPixelSample {
            xyz,
            ray,
            wavelengths_nm,
            spectral,
        })
    }
}

fn unit(word: u32) -> f64 {
    (f64::from(word) + 0.5) / 4_294_967_296.0
}

#[cfg(test)]
mod tests;
