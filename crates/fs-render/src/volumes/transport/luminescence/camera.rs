//! Physical-camera observation of original discrete luminescence powers.
//!
//! Pixel and pupil samples average radiance and transmission independently at
//! each original spectral node, using the native absorption/emission integral.

use super::{LuminescenceGrid, LuminescenceRaySample, VolumeTransportError, checkpoint};
use crate::camera::{CameraError, LensSample, PhysicalCamera};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, cie_x, cie_y, cie_z, y_integral};
use crate::volumes::DvrError;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_math::det;
use fs_rand::StreamKey;

/// Shared volume-camera sampling and whole-image limits. The luminescence
/// field retains its independent per-ray traversal limit.
pub use crate::volumes::opacity::TransferCameraSettings as LuminescenceCameraSettings;

/// Original camera, image-admission and physical-transport refusals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LuminescenceCameraError {
    /// Invalid image size, sampling limit, or output allocation.
    Image(DvrError),
    /// The camera cannot generate these achromatic geometric rays.
    Camera(CameraError),
    /// Native emission/absorption, cancellation, or cell-visit budget refusal.
    Transport(VolumeTransportError),
    /// Approximate renderer XYZ is defined only over 380–780 nm. The original
    /// spectral image remains usable on its source-supported UV/IR nodes.
    ObserverWavelength { lane: usize },
}

impl From<DvrError> for LuminescenceCameraError {
    fn from(error: DvrError) -> Self {
        Self::Image(error)
    }
}
impl From<CameraError> for LuminescenceCameraError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}
impl From<VolumeTransportError> for LuminescenceCameraError {
    fn from(error: VolumeTransportError) -> Self {
        Self::Transport(error)
    }
}
impl core::fmt::Display for LuminescenceCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "luminescence camera render refused: {self:?}")
    }
}
impl core::error::Error for LuminescenceCameraError {}

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
        let mut result = self.emission;
        for lane in 0..M {
            let background = background_w_m2_sr[lane];
            if !background.is_finite() || background < 0.0 {
                return Err(DvrError::InvalidBackgroundRadiance.into());
            }
            result[lane] = if self.transmission[lane] < f64::MIN_POSITIVE && background > 0.0 {
                self.emission[lane] + det::exp(det::ln(background) - self.effective_depth[lane])
            } else {
                self.transmission[lane].mul_add(background, self.emission[lane])
            };
            if !result[lane].is_finite() {
                return Err(VolumeTransportError::Unrepresentable.into());
            }
        }
        Ok(result)
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
        checkpoint(cx)?;
        let normalization = 1.0 / y_integral();
        let mut observer = [[0.0; 3]; M];
        for (lane, &wavelength) in self.wavelengths_nm.iter().enumerate() {
            if !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength) {
                return Err(LuminescenceCameraError::ObserverWavelength { lane });
            }
            observer[lane] = [cie_x(wavelength), cie_y(wavelength), cie_z(wavelength)]
                .map(|value| value * normalization);
        }
        let mut image = Vec::new();
        image.try_reserve_exact(self.pixels.len()).map_err(|_| {
            DvrError::ImageAllocationRefused {
                pixels: self.pixels.len(),
            }
        })?;
        for pixel in &self.pixels {
            checkpoint(cx)?;
            let mut xyz = [Compensated::default(); 3];
            for (radiance, weights) in pixel.emission.iter().zip(&observer) {
                for channel in 0..3 {
                    xyz[channel].add(radiance * weights[channel]);
                }
            }
            let xyz = xyz.map(|sum| sum.value);
            if xyz.iter().any(|value| !value.is_finite() || *value < 0.0) {
                return Err(VolumeTransportError::Unrepresentable.into());
            }
            image.push(xyz);
        }
        checkpoint(cx)?;
        Ok(image)
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
    checkpoint(cx)?;
    let [width, height] = settings.resolution;
    let pixels = width
        .checked_mul(height)
        .filter(|_| width > 0 && height > 0)
        .ok_or(DvrError::InvalidImageSize)?;
    if settings.samples_per_pixel == 0 {
        return Err(DvrError::InvalidSamplingBudget.into());
    }
    if pixels > settings.max_pixels {
        return Err(DvrError::PixelBudgetExceeded {
            requested: pixels,
            limit: settings.max_pixels,
        }
        .into());
    }
    let samples = u64::try_from(pixels)
        .ok()
        .and_then(|pixels| pixels.checked_mul(u64::from(settings.samples_per_pixel)))
        .ok_or(DvrError::SampleCountOverflow)?;
    if samples > settings.max_samples {
        return Err(DvrError::SampleBudgetExceeded {
            requested: samples,
            limit: settings.max_samples,
        }
        .into());
    }
    if u32::try_from(pixels - 1).is_err() {
        return Err(DvrError::PixelIdentityOverflow.into());
    }
    let half_y = camera.raster_vertical_half_tan();
    let half_x = half_y * (width as f64 / height as f64);
    if !half_x.is_finite() {
        return Err(CameraError::InvalidProjection.into());
    }
    let centre_lens = LensSample::try_new(0.5, 0.5)?;
    camera.generate_ray_from_tangent_offsets(cx, 0.0, 0.0, centre_lens)?;
    let mut image = Vec::new();
    image
        .try_reserve_exact(pixels)
        .map_err(|_| DvrError::ImageAllocationRefused { pixels })?;
    let mut remaining_cells = settings.max_cell_visits;
    for py in 0..height {
        for px in 0..width {
            checkpoint(cx)?;
            let pixel = py * width + px;
            let mut stream = StreamKey {
                seed: settings.seed,
                kernel: 0x6c75_6d63,
                tile: u32::try_from(pixel).map_err(|_| DvrError::PixelIdentityOverflow)?,
            }
            .stream();
            let mut average = SpectralAverage::default();
            for sample in 0..settings.samples_per_pixel {
                checkpoint(cx)?;
                let (jx, jy, lens) = if settings.samples_per_pixel == 1 {
                    (0.5, 0.5, centre_lens)
                } else {
                    (
                        stream.next_f64(),
                        stream.next_f64(),
                        LensSample::try_new(stream.next_f64(), stream.next_f64())?,
                    )
                };
                let x_tan = (2.0 * (px as f64 + jx) / width as f64 - 1.0) * half_x;
                let y_tan = (1.0 - 2.0 * (py as f64 + jy) / height as f64) * half_y;
                let primary = camera.generate_ray_from_tangent_offsets(cx, x_tan, y_tan, lens)?;
                let transported = field.transport_counted(
                    cx,
                    primary.origin,
                    primary.dir,
                    f64::INFINITY,
                    remaining_cells,
                )?;
                remaining_cells -= u64::from(transported.cell_visits);
                average.add(transported, sample + 1);
            }
            image.push(average.finish(settings.samples_per_pixel)?);
        }
    }
    checkpoint(cx)?;
    Ok(LuminescenceCameraImage {
        resolution: settings.resolution,
        wavelengths_nm: *field.wavelengths_nm(),
        source_identity: field.source_identity(),
        field_identity: field.content_identity(),
        pixels: image,
    })
}

// Online compensated means avoid overflowing sums of finite radiances.
// Separate opacity and log-sum-exp survival preserve thin layers, partial pupil
// coverage, and finite bright backgrounds behind opaque columns.
#[derive(Clone, Copy, Default)]
struct Compensated {
    value: f64,
    correction: f64,
}

impl Compensated {
    fn add(&mut self, value: f64) {
        let corrected = value - self.correction;
        let next = self.value + corrected;
        self.correction = (next - self.value) - corrected;
        self.value = next;
    }

    fn update_mean(&mut self, sample: f64, count: u32) {
        self.add((sample - self.value) / f64::from(count));
    }
}

struct SpectralAverage<const M: usize> {
    emission: [Compensated; M],
    transmission: [Compensated; M],
    opacity: [Compensated; M],
    minimum_depth: [f64; M],
    scaled_transmission: [Compensated; M],
}

impl<const M: usize> Default for SpectralAverage<M> {
    fn default() -> Self {
        Self {
            emission: [Compensated::default(); M],
            transmission: [Compensated::default(); M],
            opacity: [Compensated::default(); M],
            minimum_depth: [f64::INFINITY; M],
            scaled_transmission: [Compensated::default(); M],
        }
    }
}

impl<const M: usize> SpectralAverage<M> {
    fn add(&mut self, sample: LuminescenceRaySample<M>, count: u32) {
        for lane in 0..M {
            let depth = sample.optical_depth[lane];
            self.emission[lane]
                .update_mean(sample.transport.emitted_radiance_w_m2_sr[lane], count);
            self.transmission[lane].update_mean(sample.transport.transmittance[lane], count);
            self.opacity[lane].update_mean(-det::expm1(-depth), count);
            if depth.is_finite() {
                if depth < self.minimum_depth[lane] {
                    let scale = det::exp(depth - self.minimum_depth[lane]);
                    self.scaled_transmission[lane].value *= scale;
                    self.scaled_transmission[lane].correction *= scale;
                    self.minimum_depth[lane] = depth;
                }
                self.scaled_transmission[lane].add(det::exp(self.minimum_depth[lane] - depth));
            }
        }
    }

    fn finish(self, count: u32) -> Result<LuminescenceCameraPixel<M>, VolumeTransportError> {
        let mut pixel = LuminescenceCameraPixel {
            emission: self.emission.map(|sum| sum.value),
            transmission: self.transmission.map(|sum| sum.value),
            opacity: self.opacity.map(|sum| sum.value),
            effective_depth: [0.0; M],
        };
        for lane in 0..M {
            if !pixel.emission[lane].is_finite()
                || pixel.emission[lane] < 0.0
                || !(0.0..=1.0).contains(&pixel.transmission[lane])
                || !(0.0..=1.0).contains(&pixel.opacity[lane])
            {
                return Err(VolumeTransportError::Unrepresentable);
            }
            pixel.effective_depth[lane] = if pixel.opacity[lane] <= 0.5 {
                let transmission = 1.0 - pixel.opacity[lane];
                -det::ln(transmission)
                    - ((1.0 - transmission) - pixel.opacity[lane]) / transmission
            } else if self.scaled_transmission[lane].value > 0.0 {
                self.minimum_depth[lane]
                    - det::ln(self.scaled_transmission[lane].value / f64::from(count))
            } else {
                f64::INFINITY
            };
        }
        Ok(pixel)
    }
}

#[cfg(test)]
mod tests;
