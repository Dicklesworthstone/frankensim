//! Physical-camera observation of original discrete luminescence powers.
//!
//! Pixel and pupil samples average radiance and transmission independently at
//! each original spectral node, using the native absorption/emission integral.

use super::super::camera::{
    CameraSample, approximate_xyz, composite, render_camera, render_camera_radiance,
};
use super::super::checkpoint;
use super::{LuminescenceGrid, LuminescenceScatteringBudget, ScatteringLuminescenceGrid};
use crate::camera::PhysicalCamera;
use crate::volumes::DvrError;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// Per-ray illumination packet shared by physical volume cameras.
pub use super::super::camera::VolumeCameraBackground;

/// Original camera, image-admission and physical-transport refusals.
pub use super::super::camera::VolumeCameraError as LuminescenceCameraError;
/// Shared volume-camera sampling and whole-image limits. The luminescence
/// field retains its independent per-ray traversal limit.
pub use super::super::camera::VolumeCameraSettings as LuminescenceCameraSettings;

/// Nonthermal radiance and opacity averaged over one pixel/pupil footprint.
/// All arrays retain the parent image's original wavelength order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LuminescenceCameraPixel<const M: usize> {
    emission: [f64; M],
    transmission: [f64; M],
    opacity: [f64; M],
    effective_depth: [f64; M],
}

impl<const M: usize> LuminescenceCameraPixel<M> {
    /// Emitted line/bin-integrated radiance in W/(m² sr), not per nanometre.
    /// Foreground and self-absorption have already attenuated these values.
    #[must_use]
    pub const fn emitted_radiance_w_m2_sr(&self) -> &[f64; M] {
        &self.emission
    }

    /// Mean transmitted fraction at each original wavelength.
    #[must_use]
    pub const fn transmittance(&self) -> &[f64; M] {
        &self.transmission
    }

    /// Mean absorbed fraction, evaluated independently so optically thin
    /// attenuation survives when transmission rounds to one.
    #[must_use]
    pub const fn opacity(&self) -> &[f64; M] {
        &self.opacity
    }

    /// Effective negative log of mean transmission, not mean ray optical depth.
    /// Finite opaque columns retain depth even if unit-radiance survival underflows.
    #[must_use]
    pub const fn effective_optical_depth(&self) -> &[f64; M] {
        &self.effective_depth
    }

    /// Add background radiance specified in the same integrated spectral bins.
    /// Each lane uses emitted radiance plus transmitted background, with no
    /// second alpha multiplication. The background must be constant over this
    /// pixel/pupil footprint; varying backgrounds need per-ray transport to
    /// retain their correlation with intervening opacity.
    pub fn composite(
        &self,
        background_w_m2_sr: [f64; M],
    ) -> Result<[f64; M], LuminescenceCameraError> {
        composite(
            &self.emission,
            &self.transmission,
            &self.effective_depth,
            background_w_m2_sr,
        )
    }
}

/// Complete spectral image bound to one original luminescence snapshot.
/// Wavelengths and source identity are retained once, independently of pixels.
#[derive(Clone, Debug, PartialEq)]
pub struct LuminescenceCameraImage<const M: usize> {
    resolution: [usize; 2],
    wavelengths_nm: [f64; M],
    source_identity: ContentHash,
    field_identity: ContentHash,
    pixels: Vec<LuminescenceCameraPixel<M>>,
}

impl<const M: usize> LuminescenceCameraImage<M> {
    /// Width and height; row-major pixels start with the top image row.
    #[must_use]
    pub const fn resolution(&self) -> [usize; 2] {
        self.resolution
    }

    /// Original discrete spectral-node/bin wavelengths in vacuum nanometres.
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; M] {
        &self.wavelengths_nm
    }

    /// Original declared photophysical source and condition identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Complete luminescence field identity, including absorption and powers.
    #[must_use]
    pub const fn field_identity(&self) -> ContentHash {
        self.field_identity
    }

    /// Complete row-major image with no display or exposure transformation.
    #[must_use]
    pub fn pixels(&self) -> &[LuminescenceCameraPixel<M>] {
        &self.pixels
    }

    /// Project visible nodes into the ordinary renderer's approximate relative
    /// XYZ convention, using its analytic CIE fits and luminance normalization.
    /// These are not calibrated cd/m² values.
    ///
    /// Each node contributes L_bin times its matching function, divided by the
    /// renderer's y_integral(). Powers were already integrated over their bins,
    /// so no wavelength spacing or width is multiplied in. A coarse authored
    /// bin is represented by its declared node; neither this operation nor the
    /// camera bounds the unknown spectrum inside that bin. UV/IR nodes refuse
    /// this visible-only projection instead of acquiring an extrapolated
    /// observer. The original spectral image remains intact.
    ///
    /// Output allocation is bounded by the already admitted image size. No
    /// partial projection is published on cancellation or numerical refusal.
    pub fn approximate_emission_xyz(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<[f64; 3]>, LuminescenceCameraError> {
        self.approximate_composited_xyz(cx, [0.0; M])
    }

    /// Approximate relative XYZ including one uniform background spectrum.
    /// Each background lane is already line/bin-integrated in W/(m² sr), in
    /// the image's original wavelength order. No bin width or node spacing is
    /// multiplied in. Observer support and the limitations of representing a
    /// finite bin by its declared node are the same as
    /// [Self::approximate_emission_xyz].
    ///
    /// Background values must be finite and nonnegative and constant over every
    /// pixel/pupil footprint. Foreground emission and transmitted illumination
    /// are projected separately, without requiring their unweighted sum to fit in
    /// f64. Effective log transmission preserves bright-background attenuation
    /// through finite opaque columns whose stored transmission is zero.
    pub fn approximate_composited_xyz(
        &self,
        cx: &Cx<'_>,
        background_w_m2_sr: [f64; M],
    ) -> Result<Vec<[f64; 3]>, LuminescenceCameraError> {
        approximate_xyz(
            cx,
            &self.wavelengths_nm,
            &[1.0; M],
            self.pixels.len(),
            |pixel, lane| {
                let pixel = &self.pixels[pixel];
                (
                    pixel.emission[lane],
                    pixel.transmission[lane],
                    pixel.effective_depth[lane],
                    background_w_m2_sr[lane],
                )
            },
        )
    }

    /// Approximate relative XYZ over one bin-integrated background per pixel.
    /// Backgrounds must exactly match [Self::pixels] in top-row-first
    /// row-major order, with each lane in W/(m² sr) at its original node.
    ///
    /// Each supplied spectrum must be constant over its pixel/pupil footprint.
    /// Subpixel or pupil variation requires per-ray background transport to
    /// preserve its correlation with opacity. Projection uses the same separate
    /// contributions and observer convention as
    /// [Self::approximate_composited_xyz], with no intermediate spectral image
    /// and no partial result on cancellation or refusal.
    pub fn approximate_composited_xyz_over_image(
        &self,
        cx: &Cx<'_>,
        backgrounds_w_m2_sr: &[[f64; M]],
    ) -> Result<Vec<[f64; 3]>, LuminescenceCameraError> {
        checkpoint(cx)?;
        if backgrounds_w_m2_sr.len() != self.pixels.len() {
            return Err(DvrError::InvalidBackgroundSize.into());
        }
        approximate_xyz(
            cx,
            &self.wavelengths_nm,
            &[1.0; M],
            self.pixels.len(),
            |pixel, lane| {
                (
                    self.pixels[pixel].emission[lane],
                    self.pixels[pixel].transmission[lane],
                    self.pixels[pixel].effective_depth[lane],
                    backgrounds_w_m2_sr[pixel][lane],
                )
            },
        )
    }
}

/// Observed bin-integrated radiance after per-ray transport.
///
/// Original source bins retain their declared nodes and integrated units.
/// Filtering averages complete ray observations, including correlated
/// backgrounds or multiply scattered emission when explicitly selected.
#[derive(Clone, Debug, PartialEq)]
pub struct LuminescenceCameraRadianceImage<const M: usize> {
    resolution: [usize; 2],
    wavelengths_nm: [f64; M],
    source_identity: ContentHash,
    field_identity: ContentHash,
    radiance: Vec<[f64; M]>,
}

impl<const M: usize> LuminescenceCameraRadianceImage<M> {
    /// Width and height; row-major pixels start with the top image row.
    #[must_use]
    pub const fn resolution(&self) -> [usize; 2] {
        self.resolution
    }

    /// Original discrete spectral-node/bin wavelengths in vacuum nanometres.
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; M] {
        &self.wavelengths_nm
    }

    /// Original photophysical source and condition identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Accepted emission/transport field identity, independent of any
    /// callback-owned illumination source. Scattered observations retain the
    /// full original scattering and phase snapshot in this identity.
    #[must_use]
    pub const fn field_identity(&self) -> ContentHash {
        self.field_identity
    }

    /// Observed row-major bin-integrated radiance in W/(m² sr), without a
    /// wavelength-width factor, exposure or display transformation.
    #[must_use]
    pub fn radiance_w_m2_sr(&self) -> &[[f64; M]] {
        &self.radiance
    }

    /// Approximate relative XYZ of the already observed bin-integrated
    /// radiance. No wavelength width or second opacity factor is applied.
    ///
    /// The renderer's analytic CIE fits represent each bin at its declared
    /// node, with its ordinary relative normalization. Unknown within-bin
    /// spectra are not reconstructed, and these are not calibrated luminance
    /// values. Nodes outside 380–780 nm refuse visible projection while
    /// retaining the original spectral image. Cancellation/refusal publishes
    /// no partial result.
    pub fn approximate_xyz(&self, cx: &Cx<'_>) -> Result<Vec<[f64; 3]>, LuminescenceCameraError> {
        approximate_xyz(
            cx,
            &self.wavelengths_nm,
            &[1.0; M],
            self.radiance.len(),
            |pixel, lane| (self.radiance[pixel][lane], 1.0, 0.0, 0.0),
        )
    }
}

/// Render accepted nonthermal cell powers through a pinhole or ideal thin lens.
///
/// The field and camera geometry are in metres. Every ray uses the same exact
/// native cell integral; only pixel area and pupil area are sampled. One sample
/// uses the pixel centre and pupil variate (0.5, 0.5). Output retains the complete
/// original spectrum, self-absorption, foreground absorption, and
/// wavelength-dependent partial coverage. Radiance and transmission are averaged
/// independently, without an opacity threshold. Missed rays contribute zero
/// foreground and unit transmission. Actual cell visits consume one whole-image
/// budget alongside the field's per-ray limit.
///
/// These are normalized radiance averages, not aperture gain, exposure, sensor
/// irradiance or display RGB. Dispersive, diffractive and sequential lenses keep
/// their native refusal because this adapter uses achromatic geometric rays.
/// Thermal emission remains separate; no scattering, surfaces, refraction,
/// secondary fluorescence or time-of-flight is inferred. Cancellation or any
/// work/physical refusal returns no partial image.
pub fn render_luminescence_camera<const M: usize>(
    cx: &Cx<'_>,
    field: &LuminescenceGrid<'_, '_, M>,
    camera: &PhysicalCamera,
    settings: LuminescenceCameraSettings,
) -> Result<LuminescenceCameraImage<M>, LuminescenceCameraError> {
    let pixels = render_camera(
        cx,
        camera,
        settings,
        0x6c75_6d63,
        |origin, direction, remaining_cells| {
            let transported =
                field.transport_counted(cx, origin, direction, f64::INFINITY, remaining_cells)?;
            Ok(CameraSample {
                emission: transported.transport.emitted_radiance_w_m2_sr,
                transmittance: transported.transport.transmittance,
                optical_depth: transported.optical_depth,
                cell_visits: transported.cell_visits,
            })
        },
        |sample| LuminescenceCameraPixel {
            emission: sample.emission,
            transmission: sample.transmission,
            opacity: sample.opacity,
            effective_depth: sample.effective_depth,
        },
    )?;
    Ok(LuminescenceCameraImage {
        resolution: settings.resolution,
        wavelengths_nm: *field.wavelengths_nm(),
        source_identity: field.source_identity(),
        field_identity: field.content_identity(),
        pixels,
    })
}

/// Render correlated illumination through original luminescence and opacity.
///
/// Once per actual image sample, the callback receives the world-space ray
/// origin, unit direction, original wavelength nodes and cancellable Cx. It
/// supplies bin-integrated W/(m² sr) at that ray's endpoint, plus a nonnegative
/// distance in metres or positive infinity for a distant background. Native
/// cell transport stops at that endpoint, and emitted plus transmitted
/// radiance is composed before pixel/pupil filtering. No spectral width or RGB
/// interpretation is inferred.
///
/// This shares the original luminescence camera's exact raster/pupil stream
/// and work limits. The callback owns its illumination evaluation; the grid
/// is traversed once per sample. Each composed ray must be representable.
/// Cancellation, callback failure and native refusals publish no partial
/// image. Aperture/exposure gain, scattering, refraction and thermal sources
/// remain outside this nonthermal, normalized-radiance observation.
pub fn render_luminescence_camera_with_background<const M: usize>(
    cx: &Cx<'_>,
    field: &LuminescenceGrid<'_, '_, M>,
    camera: &PhysicalCamera,
    settings: LuminescenceCameraSettings,
    mut background: impl FnMut(
        &Cx<'_>,
        Point3,
        Vec3,
        &[f64; M],
    ) -> Result<VolumeCameraBackground<M>, LuminescenceCameraError>,
) -> Result<LuminescenceCameraRadianceImage<M>, LuminescenceCameraError> {
    let radiance = render_camera_radiance(
        cx,
        camera,
        settings,
        0x6c75_6d63,
        |origin, direction, remaining_cells| {
            let background = background(cx, origin, direction, field.wavelengths_nm())?;
            let sample = field.transport_counted(
                cx,
                origin,
                direction,
                background.distance_m,
                remaining_cells,
            )?;
            Ok(CameraSample {
                emission: composite(
                    &sample.transport.emitted_radiance_w_m2_sr,
                    &sample.transport.transmittance,
                    &sample.optical_depth,
                    background.radiance,
                )?,
                transmittance: sample.transport.transmittance,
                optical_depth: sample.optical_depth,
                cell_visits: sample.cell_visits,
            })
        },
    )?;
    Ok(LuminescenceCameraRadianceImage {
        resolution: settings.resolution,
        wavelengths_nm: *field.wavelengths_nm(),
        source_identity: field.source_identity(),
        field_identity: field.content_identity(),
        radiance,
    })
}

/// Render accepted nonthermal emission through elastic multiple scattering.
///
/// Each original raster/pupil sample launches one independent spectral path,
/// using the same pinhole/thin-lens geometry and the same pixel stream as the
/// absorption-only luminescence camera. Native exact extinction, the phase law
/// of each actual collision cell, and the original sparse emitting-cell map
/// determine the result. The output remains bin-integrated W/(m² sr): no LTE
/// source, wavelength width, exposure gain or display transformation is added.
///
/// Pixel samples average complete path radiances. These estimates carry Monte
/// Carlo error; this function does not infer a confidence bound from sample
/// count. Both per-history limits and the actual cumulative whole-image cell
/// allowance are binding. A history that reaches a surviving collision cap
/// refuses the whole image rather than silently dropping unresolved light.
/// Retry the entire image with the same seed and a larger work allowance.
///
/// This is an open, stationary, index-one volume observation. Background
/// illumination, surfaces, refraction, wavelength-changing scattering and
/// photon flight delays remain separate. Cancellation/refusal returns no
/// partial image. With no scattering, the original exact source integral and
/// raster/pupil draw order are retained bit-for-bit.
pub fn render_scattered_luminescence_camera<const M: usize>(
    cx: &Cx<'_>,
    field: &ScatteringLuminescenceGrid<'_, '_, '_, M>,
    camera: &PhysicalCamera,
    settings: LuminescenceCameraSettings,
    budget: LuminescenceScatteringBudget,
) -> Result<LuminescenceCameraRadianceImage<M>, LuminescenceCameraError> {
    let mut path_id = 0_u64;
    let radiance = render_camera_radiance(
        cx,
        camera,
        settings,
        0x6c75_6d63,
        |origin, direction, remaining_cells| {
            let sample = field.sample_ray(
                cx,
                origin,
                direction,
                LuminescenceScatteringBudget {
                    max_cell_visits: budget
                        .max_cell_visits
                        .min(u32::try_from(remaining_cells).unwrap_or(u32::MAX)),
                    ..budget
                },
                settings.seed,
                path_id,
            )?;
            path_id = path_id
                .checked_add(1)
                .ok_or(super::VolumeTransportError::Budget)?;
            Ok(CameraSample {
                emission: sample.radiance_w_m2_sr,
                // The shared radiance-only accumulator reads no opacity.
                // Scattered light cannot be represented by a scalar screen.
                transmittance: [1.0; M],
                optical_depth: [0.0; M],
                cell_visits: sample.cell_visits,
            })
        },
    )?;
    Ok(LuminescenceCameraRadianceImage {
        resolution: settings.resolution,
        wavelengths_nm: *field.wavelengths_nm(),
        source_identity: field.source_identity(),
        field_identity: field.content_identity(),
        radiance,
    })
}

#[cfg(test)]
mod tests;
