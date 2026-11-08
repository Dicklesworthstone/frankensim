//! Physical-camera observation of original discrete luminescence powers.
//!
//! Pixel and pupil samples average radiance and transmission independently at
//! each original spectral node, using the native absorption/emission integral.

use super::super::camera::{CameraSample, approximate_xyz, composite, render_camera};
use super::LuminescenceGrid;
use crate::camera::PhysicalCamera;
use fs_blake3::ContentHash;
use fs_exec::Cx;

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
        approximate_xyz(
            cx,
            &self.wavelengths_nm,
            &[1.0; M],
            self.pixels.iter().map(|pixel| &pixel.emission),
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

#[cfg(test)]
mod tests;
