//! Physical-camera pixel samples with independent spectral curved paths.
//!
//! The existing pinhole/thin-lens camera supplies geometric rays. Four
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
//! Diffraction requires the prepared homogeneous source-plane proof in
//! [`wave::CurvedWaveCameraRaster`]; the generic geometric raster rejects it.
//! An explicit geometric dispersive lens retains its native pupil and phase
//! index laws, generating a separate camera ray for every sampled wavelength.
//! Resolved sequential lenses use [`CurvedCameraRaster::sample_pixel_sequential`]
//! and its detector variants: each live lane carries the full reciprocal lens
//! polarization into the scene. Aperture/TIR misses have no ray and retain
//! the original zero-contribution wavelength draw.

use super::analyzer::AnalyzerResponse;
use super::{
    CameraAnalyzer, CameraRadianceError, CameraRadianceSample, CurvedCameraScene,
    SpectralCameraAnalyzer,
};
use crate::camera::{CameraError, LensSample, PhysicalCamera, SequentialCameraBeam};
use crate::charts::Ray;
use crate::hero_wavelength_packet;
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, cie_x, cie_y, cie_z, y_integral};
use fs_exec::Cx;
use fs_rand::philox::philox4x32_10;

// Passive analyzers share the original packet/raster and detector owners.
mod analyzer;
// Path-only Smith transport shares these same geometric and optical samples.
mod microsurface;

/// Spectral pixel integration with native curved molecular scattering.
#[cfg(feature = "volumes")]
pub mod participating;

/// Spectral pixels with selected refractive point-fixture connections.
#[cfg(feature = "photometry")]
pub mod point_lights;

/// Spectral pixels with finite-area refractive light and emitter-hit MIS.
pub mod area_lights;

/// Full Smith glass pixels with original environment source/work receipts.
pub mod environment;

/// Complete Smith thermal emission and reflection in native CIE/QE pixels.
pub mod thermal;

/// Source-bounded multichannel detector signals, including infrared bands.
pub mod sensor;
mod sequential;
pub use sequential::CurvedSequentialPixelSample;
#[cfg(feature = "volumes")]
pub use sequential::CurvedSequentialRayleighPixelSample;
/// Native finite-source histories at the original observation clock.
pub mod transient;

/// Native wave imaging through a proven clear homogeneous source scene.
pub mod wave;
pub use wave::CurvedWaveCameraRaster;

/// Each pixel sample independently traces this many wavelengths.
pub const CAMERA_WAVELENGTHS: usize = 4;

#[derive(Clone, Debug, PartialEq)]
pub enum CurvedRasterError {
    Camera(CameraError),
    Transport(CameraRadianceError),
    Sensor(crate::spectral::detector::SpectralSensorError),
    WaveScene(wave::WaveSceneError),
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
impl From<crate::spectral::detector::SpectralSensorError> for CurvedRasterError {
    fn from(error: crate::spectral::detector::SpectralSensorError) -> Self {
        Self::Sensor(error)
    }
}
impl From<wave::WaveSceneError> for CurvedRasterError {
    fn from(error: wave::WaveSceneError) -> Self {
        Self::WaveScene(error)
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
    rays: [Ray; CAMERA_WAVELENGTHS],
    wavelengths_nm: [f64; CAMERA_WAVELENGTHS],
    spectral: [CameraRadianceSample; CAMERA_WAVELENGTHS],
    diffraction_weights: [f64; CAMERA_WAVELENGTHS],
}

impl CurvedPixelSample {
    /// Add this normalized XYZ contribution to a pixel sum, then divide by
    /// that pixel's completed sample count. No exposure or display transform.
    pub const fn xyz(&self) -> [f64; 3] {
        self.xyz
    }
    /// The first spectral lane's camera ray; use [`Self::camera_rays`] when
    /// an explicit geometric dispersive lens is selected.
    pub const fn camera_ray(&self) -> Ray {
        self.rays[0]
    }
    /// Camera rays in wavelength order. They coincide for an achromatic
    /// camera; a geometric dispersive lens gives each lane its own ray.
    /// Wave imaging shares a proposal ray, retaining each native spectral path.
    pub const fn camera_rays(&self) -> &[Ray; CAMERA_WAVELENGTHS] {
        &self.rays
    }
    pub const fn wavelengths_nm(&self) -> &[f64; CAMERA_WAVELENGTHS] {
        &self.wavelengths_nm
    }
    pub const fn spectral_samples(&self) -> &[CameraRadianceSample; CAMERA_WAVELENGTHS] {
        &self.spectral
    }
    /// Optical pupil/proposal weights already applied to [`Self::xyz`].
    /// The retained spectral samples are original transport radiances before
    /// these weights. Geometric cameras retain exact unit weights.
    pub const fn diffraction_weights(&self) -> &[f64; CAMERA_WAVELENGTHS] {
        &self.diffraction_weights
    }
}

// A clipped wavelength has neither an outgoing ray nor a fabricated path.
// The original mandatory-ray API is adapted only from complete packets.
#[derive(Clone, Copy, Debug, PartialEq)]
struct OptionalCameraPacket {
    rays: [Option<Ray>; CAMERA_WAVELENGTHS],
    wavelengths_nm: [f64; CAMERA_WAVELENGTHS],
    spectral: [Option<CameraRadianceSample>; CAMERA_WAVELENGTHS],
    optical_weights: [f64; CAMERA_WAVELENGTHS],
}

#[derive(Clone, Copy)]
enum CameraLaunch {
    Geometric(Ray),
    Sequential(SequentialCameraBeam),
}

impl CameraLaunch {
    fn ray(self) -> Ray {
        match self {
            Self::Geometric(ray) => ray,
            Self::Sequential(beam) => beam.ray(),
        }
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
        if camera.has_diffraction() {
            return Err(CameraError::InvalidDiffraction.into());
        }
        Self::prepare(camera, width, height)
    }

    // Only the geometric constructor and the separately proven wave wrapper
    // can prepare this owner. The wrapper never exposes its internal raster.
    fn prepare(
        camera: &'a PhysicalCamera,
        width: u32,
        height: u32,
    ) -> Result<Self, CurvedRasterError> {
        let pixels = width
            .checked_mul(height)
            .filter(|&count| count > 0)
            .ok_or(CurvedRasterError::InvalidRaster)?;
        let half_tan = camera.raster_vertical_half_tan();
        let aspect = f64::from(width) / f64::from(height);
        if !(half_tan > 0.0 && (half_tan * aspect).is_finite()) {
            return Err(CameraError::InvalidProjection.into());
        }
        camera.admit_dispersive_raster(half_tan * aspect, half_tan)?;
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
        trace: impl FnMut(&Cx<'_>, Ray, f64, usize) -> Result<CameraRadianceSample, CameraRadianceError>,
    ) -> Result<CurvedPixelSample, CurvedRasterError> {
        let (rays, wavelengths_nm, spectral, diffraction_weights) = self
            .sample_packet_with_transport(
                cx,
                seed,
                pixel,
                sample,
                |u| {
                    Ok(hero_wavelength_packet::<CAMERA_WAVELENGTHS>(
                        LAMBDA_MIN + u * (LAMBDA_MAX - LAMBDA_MIN),
                        LAMBDA_MIN,
                        LAMBDA_MAX,
                    ))
                },
                trace,
            )?;
        let mut xyz = [0.0; 3];
        for ((lane, &wavelength), &optical_weight) in spectral
            .iter()
            .zip(&wavelengths_nm)
            .zip(&diffraction_weights)
        {
            let weight =
                weighted_radiance(lane.radiance_w_m2_sr_nm(), optical_weight)? * self.xyz_scale;
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
        cx.checkpoint().map_err(CameraError::from)?;
        Ok(CurvedPixelSample {
            xyz,
            rays,
            wavelengths_nm,
            spectral,
            diffraction_weights,
        })
    }

    // One raster/lens/transport owner serves visible XYZ and source-bounded
    // detectors. The wavelength proposal never lends one lane another's path.
    fn sample_packet_with_transport(
        &self,
        cx: &Cx<'_>,
        seed: u64,
        pixel: u32,
        sample: u32,
        wavelengths: impl FnOnce(f64) -> Result<[f64; CAMERA_WAVELENGTHS], CurvedRasterError>,
        mut trace: impl FnMut(
            &Cx<'_>,
            Ray,
            f64,
            usize,
        ) -> Result<CameraRadianceSample, CameraRadianceError>,
    ) -> Result<
        (
            [Ray; CAMERA_WAVELENGTHS],
            [f64; CAMERA_WAVELENGTHS],
            [CameraRadianceSample; CAMERA_WAVELENGTHS],
            [f64; CAMERA_WAVELENGTHS],
        ),
        CurvedRasterError,
    > {
        if self.camera.sequential_lens().is_some() {
            cx.checkpoint().map_err(CameraError::from)?;
            if pixel >= self.pixels {
                return Err(CurvedRasterError::InvalidRaster);
            }
            return Err(CameraError::SequentialLensRequiresWeightedRays.into());
        }
        let packet = self.sample_optional_packet_with_transport(
            cx,
            seed,
            pixel,
            sample,
            None,
            false,
            wavelengths,
            |cx, launch, wavelength, lane| match launch {
                CameraLaunch::Geometric(ray) => Ok(trace(cx, ray, wavelength, lane)?),
                CameraLaunch::Sequential(_) => {
                    Err(CameraError::SequentialLensRequiresWeightedRays.into())
                }
            },
        )?;
        Ok((
            complete_lanes(packet.rays)?,
            packet.wavelengths_nm,
            complete_lanes(packet.spectral)?,
            packet.optical_weights,
        ))
    }

    // All camera models retain the same original raster/pupil/proposal draws
    // and transport callback. Only physically absent lanes skip that callback.
    #[allow(clippy::too_many_arguments)]
    fn sample_optional_packet_with_transport(
        &self,
        cx: &Cx<'_>,
        seed: u64,
        pixel: u32,
        sample: u32,
        analyzer: Option<AnalyzerResponse<'_>>,
        timed_lens: bool,
        wavelengths: impl FnOnce(f64) -> Result<[f64; CAMERA_WAVELENGTHS], CurvedRasterError>,
        mut trace: impl FnMut(
            &Cx<'_>,
            CameraLaunch,
            f64,
            usize,
        ) -> Result<CameraRadianceSample, CurvedRasterError>,
    ) -> Result<OptionalCameraPacket, CurvedRasterError> {
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
        let wavelengths_nm = wavelengths(unit(raster[2]))?;
        let lens_sample = LensSample::try_new(unit(lens[0]), unit(lens[1]))?;
        let diffraction_sample = if self.camera.has_diffraction() {
            let block = philox4x32_10([pixel, sample, 0, 0x6372_7333], key);
            LensSample::try_new(unit(block[0]), unit(block[1]))?
        } else {
            LensSample::CENTER
        };
        let layered = if self.camera.layered_diffraction().is_some() {
            Some(self.camera.sample_layered_diffraction_ray(
                &render_cx,
                x_tan,
                y_tan,
                diffraction_sample,
                &wavelengths_nm,
            )?)
        } else {
            None
        };
        let (launches, mut diffraction_weights) = if self.camera.sequential_lens().is_some() {
            let beams = if let Some(AnalyzerResponse::Spectral(spectrum)) = analyzer {
                // Resolve the complete source-owned packet before tracing ANY
                // lens path, including lanes that will be aperture-clipped.
                let mut rows = [crate::polarization::Stokes::unpolarized(1.0); CAMERA_WAVELENGTHS];
                for (row, &wavelength) in rows.iter_mut().zip(&wavelengths_nm) {
                    *row = spectrum.at(wavelength)?.row();
                }
                if timed_lens {
                    self.camera.generate_sequential_spectral_timed_beams_with_detector_rows_from_tangent_offsets(
                        &render_cx, x_tan, y_tan, lens_sample, &wavelengths_nm,
                        spectrum.reference_axis(), &rows,
                    )?
                } else {
                    self.camera
                        .generate_sequential_spectral_beams_with_detector_rows_from_tangent_offsets(
                            &render_cx,
                            x_tan,
                            y_tan,
                            lens_sample,
                            &wavelengths_nm,
                            spectrum.reference_axis(),
                            &rows,
                        )?
                }
            } else if timed_lens {
                if let Some(AnalyzerResponse::Constant(analyzer)) = analyzer {
                    self.camera
                        .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
                            &render_cx,
                            x_tan,
                            y_tan,
                            lens_sample,
                            &wavelengths_nm,
                            analyzer.reference_axis(),
                            analyzer.row(),
                        )?
                } else {
                    self.camera
                        .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                            &render_cx,
                            x_tan,
                            y_tan,
                            lens_sample,
                            &wavelengths_nm,
                        )?
                }
            } else if let Some(AnalyzerResponse::Constant(analyzer)) = analyzer {
                self.camera
                    .generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                        &render_cx,
                        x_tan,
                        y_tan,
                        lens_sample,
                        &wavelengths_nm,
                        analyzer.reference_axis(),
                        analyzer.row(),
                    )?
            } else {
                self.camera
                    .generate_sequential_spectral_beams_from_tangent_offsets(
                        &render_cx,
                        x_tan,
                        y_tan,
                        lens_sample,
                        &wavelengths_nm,
                    )?
            };
            (
                beams.map(|beam| beam.map(CameraLaunch::Sequential)),
                beams.map(|beam| beam.map_or(0.0, |beam| beam.geometry_weight())),
            )
        } else if let Some(token) = &layered {
            (
                [Some(CameraLaunch::Geometric(token.ray())); CAMERA_WAVELENGTHS],
                [0.0; CAMERA_WAVELENGTHS],
            )
        } else if self.camera.has_diffraction() {
            let optical = self.camera.generate_spectral_ray_from_tangent_offsets(
                &render_cx,
                x_tan,
                y_tan,
                lens_sample,
                diffraction_sample,
                &wavelengths_nm,
            )?;
            (
                [Some(CameraLaunch::Geometric(optical.ray)); CAMERA_WAVELENGTHS],
                optical.weights,
            )
        } else {
            (
                self.camera
                    .generate_geometric_spectral_rays_from_tangent_offsets(
                        &render_cx,
                        x_tan,
                        y_tan,
                        lens_sample,
                        &wavelengths_nm,
                    )?
                    .map(|ray| Some(CameraLaunch::Geometric(ray))),
                [1.0; CAMERA_WAVELENGTHS],
            )
        };
        let mut trace_lane =
            |lane: usize| -> Result<Option<CameraRadianceSample>, CurvedRasterError> {
                let Some(launch) = launches[lane] else {
                    return Ok(None);
                };
                let result = trace(&render_cx, launch, wavelengths_nm[lane], lane)?;
                if let Some(lens) = self.camera.dispersive_lens() {
                    lens.admit_exterior_phase_index(
                        wavelengths_nm[lane],
                        result.launch_phase_index(),
                    )?;
                }
                Ok(Some(result))
            };
        // Fixed storage and explicit propagation of the first failed lane.
        let spectral = [
            trace_lane(0)?,
            trace_lane(1)?,
            trace_lane(2)?,
            trace_lane(3)?,
        ];
        if self.camera.has_diffraction() {
            let resolved = complete_lanes(spectral)?;
            let points = [
                wave::source_point(&resolved[0])?,
                wave::source_point(&resolved[1])?,
                wave::source_point(&resolved[2])?,
                wave::source_point(&resolved[3])?,
            ];
            if let Some(token) = layered {
                diffraction_weights = token.resolve_spectral_sources(&render_cx, points)?.weights;
            } else {
                for ((point, &wavelength), weight) in points
                    .into_iter()
                    .zip(&wavelengths_nm)
                    .zip(&mut diffraction_weights)
                {
                    if let Some(point) = point {
                        self.camera
                            .admit_diffraction_depth(Some(point), &[wavelength])?;
                    } else {
                        *weight = 0.0;
                    }
                }
            }
        }
        Ok(OptionalCameraPacket {
            rays: launches.map(|launch| launch.map(CameraLaunch::ray)),
            wavelengths_nm,
            spectral,
            optical_weights: diffraction_weights,
        })
    }
}

fn complete_lanes<T>(
    [a, b, c, d]: [Option<T>; CAMERA_WAVELENGTHS],
) -> Result<[T; CAMERA_WAVELENGTHS], CurvedRasterError> {
    Ok([
        a.ok_or(CurvedRasterError::NumericalRange)?,
        b.ok_or(CurvedRasterError::NumericalRange)?,
        c.ok_or(CurvedRasterError::NumericalRange)?,
        d.ok_or(CurvedRasterError::NumericalRange)?,
    ])
}

fn unit(word: u32) -> f64 {
    (f64::from(word) + 0.5) / 4_294_967_296.0
}

fn weighted_radiance(radiance: f64, weight: f64) -> Result<f64, CurvedRasterError> {
    let value = radiance * weight;
    if !value.is_finite() || value < 0.0 || (radiance > 0.0 && weight > 0.0 && value == 0.0) {
        return Err(CurvedRasterError::NumericalRange);
    }
    Ok(value)
}

#[cfg(test)]
mod tests;
