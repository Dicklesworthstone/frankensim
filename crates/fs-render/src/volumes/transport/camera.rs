//! Physical-camera images of native spectral absorption and thermal emission.
//!
//! A shared raster/pupil sampler and stable spectral average also serve the
//! discrete luminescence camera. Units remain explicit at the public owners:
//! thermal radiance is per nanometre; luminescence powers are already integrated.

use super::{AbsorbingGrid, VolumeTransportError, checkpoint};
use crate::camera::{CameraError, LensSample, PhysicalCamera};
use crate::spectral::observer::{XyzAccumulator, XyzObserverError, XyzObserverLane};
use crate::spectral::y_integral;
use crate::volumes::DvrError;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;
use fs_rand::StreamKey;

/// Shared volume-camera sampling and whole-image limits. Each admitted physical
/// field retains its independent per-ray traversal limit.
pub use crate::volumes::opacity::TransferCameraSettings as VolumeCameraSettings;

/// Original camera, image-admission and physical-transport refusals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VolumeCameraError {
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

impl From<DvrError> for VolumeCameraError {
    fn from(error: DvrError) -> Self {
        Self::Image(error)
    }
}
impl From<CameraError> for VolumeCameraError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}
impl From<VolumeTransportError> for VolumeCameraError {
    fn from(error: VolumeTransportError) -> Self {
        Self::Transport(error)
    }
}
impl core::fmt::Display for VolumeCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spectral volume camera render refused: {self:?}")
    }
}
impl core::error::Error for VolumeCameraError {}

/// Spectral radiance and opacity averaged over a physical pixel/pupil footprint.
/// Arrays retain the parent image's original wavelength order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingCameraPixel<const M: usize> {
    emission: [f64; M],
    transmission: [f64; M],
    opacity: [f64; M],
    effective_depth: [f64; M],
}

impl<const M: usize> AbsorbingCameraPixel<M> {
    /// Native thermal spectral radiance in W/(m² sr nm).
    /// Foreground and self-absorption have already attenuated these values.
    #[must_use]
    pub const fn emitted_radiance_w_m2_sr_nm(&self) -> &[f64; M] {
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

    /// Add background spectral radiance in W/(m² sr nm) at the same nodes.
    /// Composition precedes any wavelength quadrature: emitted radiance plus
    /// transmitted background, without a second multiplication by opacity.
    /// The background must be constant over this pixel/pupil footprint; a
    /// varying background requires per-ray transport to retain correlation.
    pub fn composite(
        &self,
        background_w_m2_sr_nm: [f64; M],
    ) -> Result<[f64; M], VolumeCameraError> {
        composite(
            &self.emission,
            &self.transmission,
            &self.effective_depth,
            background_w_m2_sr_nm,
        )
    }
}

/// Complete spectral-density image bound to one accepted absorption/thermal
/// snapshot. Wavelengths and field identity are retained once per image.
#[derive(Clone, Debug, PartialEq)]
pub struct AbsorbingCameraImage<const M: usize> {
    resolution: [usize; 2],
    wavelengths_nm: [f64; M],
    field_identity: ContentHash,
    pixels: Vec<AbsorbingCameraPixel<M>>,
}

impl<const M: usize> AbsorbingCameraImage<M> {
    /// Width and height; row-major pixels start with the top image row.
    #[must_use]
    pub const fn resolution(&self) -> [usize; 2] {
        self.resolution
    }

    /// Original queried wavelengths in vacuum nanometres.
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; M] {
        &self.wavelengths_nm
    }

    /// Original density, absorption, and accepted thermal-source identity.
    #[must_use]
    pub const fn field_identity(&self) -> ContentHash {
        self.field_identity
    }

    /// Complete row-major image without display or exposure transformation.
    #[must_use]
    pub fn pixels(&self) -> &[AbsorbingCameraPixel<M>] {
        &self.pixels
    }

    /// Approximate relative renderer XYZ from explicit wavelength quadrature.
    ///
    /// Each spectral-density node contributes L_lambda times its supplied
    /// weight in nm and the renderer's analytic CIE fit, divided by y_integral().
    /// Weights must be finite and nonnegative; invalid weights use the native
    /// Sample diagnostic. No bin width, spectral interpolation or coverage is
    /// inferred. The caller owns wavelength resolution, integration support and
    /// quadrature accuracy; these values are not calibrated luminance.
    ///
    /// All nodes must be within the visible observer's 380–780 nm support.
    /// The original UV/IR spectral image remains available when this projection
    /// refuses. Output size is bounded by the already admitted pixel count,
    /// with cancellation and no partial result on refusal.
    pub fn approximate_emission_xyz(
        &self,
        cx: &Cx<'_>,
        quadrature_weights_nm: [f64; M],
    ) -> Result<Vec<[f64; 3]>, VolumeCameraError> {
        self.approximate_composited_xyz(cx, quadrature_weights_nm, [0.0; M])
    }

    /// Approximate relative XYZ including background spectral radiance at the
    /// original wavelength nodes, in W/(m² sr nm). The explicit quadrature
    /// weights are in nm, with the same observer and accuracy limits as
    /// [Self::approximate_emission_xyz].
    ///
    /// Background values must be finite and nonnegative. The uniform background
    /// must be constant over every pixel/pupil footprint. Foreground emission
    /// and transmitted illumination are projected separately, so their
    /// unweighted spectral sum need not be representable.
    /// Finite opaque columns retain attenuation of bright backgrounds even when
    /// stored unit-radiance transmission has underflowed to zero.
    pub fn approximate_composited_xyz(
        &self,
        cx: &Cx<'_>,
        quadrature_weights_nm: [f64; M],
        background_w_m2_sr_nm: [f64; M],
    ) -> Result<Vec<[f64; 3]>, VolumeCameraError> {
        approximate_xyz(
            cx,
            &self.wavelengths_nm,
            &quadrature_weights_nm,
            self.pixels.len(),
            |pixel, lane| {
                let pixel = &self.pixels[pixel];
                (
                    pixel.emission[lane],
                    pixel.transmission[lane],
                    pixel.effective_depth[lane],
                    background_w_m2_sr_nm[lane],
                )
            },
        )
    }

    /// Approximate relative XYZ over one background spectrum per pixel.
    /// Spectra are W/(m² sr nm); quadrature weights are nm. Backgrounds must
    /// exactly match [Self::pixels] in top-row-first row-major order.
    ///
    /// Each supplied spectrum must be constant over its pixel/pupil footprint.
    /// Subpixel or pupil variation requires per-ray background transport to
    /// preserve its correlation with opacity. As in
    /// [Self::approximate_composited_xyz], projection uses separate emission
    /// and transmitted-background contributions and returns no partial image
    /// on cancellation, an invalid background, or an unrepresentable XYZ value.
    pub fn approximate_composited_xyz_over_image(
        &self,
        cx: &Cx<'_>,
        quadrature_weights_nm: [f64; M],
        backgrounds_w_m2_sr_nm: &[[f64; M]],
    ) -> Result<Vec<[f64; 3]>, VolumeCameraError> {
        checkpoint(cx)?;
        if backgrounds_w_m2_sr_nm.len() != self.pixels.len() {
            return Err(DvrError::InvalidBackgroundSize.into());
        }
        approximate_xyz(
            cx,
            &self.wavelengths_nm,
            &quadrature_weights_nm,
            self.pixels.len(),
            |pixel, lane| {
                (
                    self.pixels[pixel].emission[lane],
                    self.pixels[pixel].transmission[lane],
                    self.pixels[pixel].effective_depth[lane],
                    backgrounds_w_m2_sr_nm[pixel][lane],
                )
            },
        )
    }
}

/// Observe an accepted spectral absorber and its original thermal sources
/// through a pinhole or ideal thin lens.
///
/// Wavelengths must be nonempty, finite, strictly increasing and inside the
/// medium's original source support. Every primary ray uses the same native
/// Beer–Lambert/constituent source integral as AbsorbingGrid::transport. Density
/// is kg/m³, mass absorption is m²/kg and geometry is metres. A grid without
/// thermal sources supplies attenuation and zero foreground; no temperature,
/// emissivity or RGB absorption is inferred.
///
/// Pixel/pupil samples average spectral radiance and transmission separately.
/// Output emission and composed backgrounds are W/(m² sr nm), unlike integrated
/// luminescence bins. Thin opacity and finite effective depths survive rounded
/// transmission. The public integrator's opaque limits remain available without
/// weakening the stricter observed-column API used by other consumers.
///
/// These are normalized radiance averages, not aperture gain, exposure, sensor
/// irradiance or display RGB. Dispersive, diffractive and sequential lenses keep
/// their native refusal. Surfaces, refraction, scattering and evolving fields
/// require their original transport owners. Each ray retains the medium's
/// crossing limit and consumes the remaining whole-image cell budget; failure
/// publishes no partial image.
pub fn render_absorbing_camera<const M: usize>(
    cx: &Cx<'_>,
    medium: &AbsorbingGrid<'_>,
    camera: &PhysicalCamera,
    wavelengths_nm: [f64; M],
    settings: VolumeCameraSettings,
) -> Result<AbsorbingCameraImage<M>, VolumeCameraError> {
    checkpoint(cx)?;
    if M == 0 {
        return Err(VolumeTransportError::Sample.into());
    }
    let (minimum, maximum) = medium.wavelength_bounds_nm();
    if wavelengths_nm
        .iter()
        .any(|wavelength| !wavelength.is_finite() || !(minimum..=maximum).contains(wavelength))
        || wavelengths_nm.windows(2).any(|pair| pair[0] >= pair[1])
    {
        return Err(VolumeTransportError::Wavelength.into());
    }
    let pixels = render_camera(
        cx,
        camera,
        settings,
        0x6162_7363,
        |origin, direction, remaining_cells| {
            let sample = medium.transport_counted(
                cx,
                origin,
                direction,
                f64::INFINITY,
                &wavelengths_nm,
                u32::try_from(remaining_cells).unwrap_or(u32::MAX),
            )?;
            Ok(CameraSample {
                emission: sample.transport.emission_w_m2_sr_nm,
                transmittance: sample.transport.transmittance,
                optical_depth: sample.optical_depth,
                cell_visits: sample.cell_visits,
            })
        },
        |sample| AbsorbingCameraPixel {
            emission: sample.emission,
            transmission: sample.transmission,
            opacity: sample.opacity,
            effective_depth: sample.effective_depth,
        },
    )?;
    Ok(AbsorbingCameraImage {
        resolution: settings.resolution,
        wavelengths_nm,
        field_identity: medium.content_identity(),
        pixels,
    })
}

// Internal samples are in their caller's radiance units. Public wrappers keep
// spectral densities separate from powers already integrated over source bins.
pub(super) struct CameraSample<const M: usize> {
    pub(super) emission: [f64; M],
    pub(super) transmittance: [f64; M],
    pub(super) optical_depth: [f64; M],
    pub(super) cell_visits: u32,
}

pub(super) struct CameraPixel<const M: usize> {
    pub(super) emission: [f64; M],
    pub(super) transmission: [f64; M],
    pub(super) opacity: [f64; M],
    pub(super) effective_depth: [f64; M],
}

pub(super) fn composite<const M: usize>(
    emission: &[f64; M],
    transmission: &[f64; M],
    effective_depth: &[f64; M],
    background: [f64; M],
) -> Result<[f64; M], VolumeCameraError> {
    let mut result = *emission;
    for lane in 0..M {
        let background = background[lane];
        if !background.is_finite() || background < 0.0 {
            return Err(DvrError::InvalidBackgroundRadiance.into());
        }
        result[lane] = if transmission[lane] < f64::MIN_POSITIVE && background > 0.0 {
            emission[lane] + det::exp(det::ln(background) - effective_depth[lane])
        } else {
            transmission[lane].mul_add(background, emission[lane])
        };
        if !result[lane].is_finite() {
            return Err(VolumeTransportError::Unrepresentable.into());
        }
    }
    Ok(result)
}

pub(super) fn approximate_xyz<const M: usize>(
    cx: &Cx<'_>,
    wavelengths_nm: &[f64; M],
    spectral_weights: &[f64; M],
    pixels: usize,
    sample: impl Fn(usize, usize) -> (f64, f64, f64, f64),
) -> Result<Vec<[f64; 3]>, VolumeCameraError> {
    checkpoint(cx)?;
    let normalization = 1.0 / y_integral();
    let mut observer = [XyzObserverLane::ZERO; M];
    for (lane, &wavelength) in wavelengths_nm.iter().enumerate() {
        observer[lane] = XyzObserverLane::new(wavelength, spectral_weights[lane], normalization)
            .map_err(|error| observer_error(error, lane))?;
    }
    let mut image = Vec::new();
    image
        .try_reserve_exact(pixels)
        .map_err(|_| DvrError::ImageAllocationRefused { pixels })?;
    for pixel in 0..pixels {
        checkpoint(cx)?;
        let mut xyz = XyzAccumulator::default();
        for (lane, observer) in observer.iter().enumerate() {
            let (emission, transmission, depth, background) = sample(pixel, lane);
            xyz.add(observer, emission, transmission, depth, background)
                .map_err(|error| observer_error(error, lane))?;
        }
        image.push(xyz.finish().map_err(|error| observer_error(error, 0))?);
    }
    checkpoint(cx)?;
    Ok(image)
}

fn observer_error(error: XyzObserverError, lane: usize) -> VolumeCameraError {
    match error {
        XyzObserverError::Wavelength => VolumeCameraError::ObserverWavelength { lane },
        XyzObserverError::Weight => VolumeTransportError::Sample.into(),
        XyzObserverError::Radiance => DvrError::InvalidBackgroundRadiance.into(),
        XyzObserverError::Unrepresentable => VolumeTransportError::Unrepresentable.into(),
    }
}

// One original raster/pupil stream and one stable average serve each native
// source owner. Luminescence retains its original kernel, draw order and math.
pub(super) fn render_camera<const M: usize, P>(
    cx: &Cx<'_>,
    camera: &PhysicalCamera,
    settings: VolumeCameraSettings,
    kernel: u32,
    mut trace: impl FnMut(Point3, Vec3, u64) -> Result<CameraSample<M>, VolumeCameraError>,
    finish: impl Fn(CameraPixel<M>) -> P,
) -> Result<Vec<P>, VolumeCameraError> {
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
                kernel,
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
                let transported = trace(primary.origin, primary.dir, remaining_cells)?;
                remaining_cells -= u64::from(transported.cell_visits);
                average.add(transported, sample + 1);
            }
            image.push(finish(average.finish(settings.samples_per_pixel)?));
        }
    }
    checkpoint(cx)?;
    Ok(image)
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
    fn add(&mut self, sample: CameraSample<M>, count: u32) {
        for lane in 0..M {
            let depth = sample.optical_depth[lane];
            self.emission[lane].update_mean(sample.emission[lane], count);
            self.transmission[lane].update_mean(sample.transmittance[lane], count);
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

    fn finish(self, count: u32) -> Result<CameraPixel<M>, VolumeTransportError> {
        let mut pixel = CameraPixel {
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
                -det::ln(transmission) - ((1.0 - transmission) - pixel.opacity[lane]) / transmission
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
