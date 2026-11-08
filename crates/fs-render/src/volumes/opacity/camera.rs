//! Camera images of exact scalar-volume transfer, with radiance/transmission
//! filtering across the original physical camera's pixel and pupil samples.

use super::{DvrError, Sum, TransferRayError, TransferRayGrid, TransferTransportPixel};
use crate::camera::{CameraError, LensSample, PhysicalCamera};
use crate::volumes::Ray;
use crate::volumes::transport::VolumeTransportError;
use fs_exec::Cx;
use fs_math::det;
use fs_rand::StreamKey;

/// Camera sampling and whole-image work limits. The admitted field retains
/// its own independent per-ray crossing limit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TransferCameraSettings {
    /// `[width, height]`, stored left-to-right with the top image row first.
    pub resolution: [usize; 2],
    /// Positive samples per pixel. One uses the pixel centre and lens variate
    /// (0.5, 0.5); larger counts sample pixel area and pupil independently.
    pub samples_per_pixel: u32,
    /// Stable per-pixel Philox seed, independent of execution scheduling.
    pub seed: u64,
    /// Maximum allocated output pixels.
    pub max_pixels: usize,
    /// Maximum primary samples across the complete image.
    pub max_samples: u64,
    /// Maximum actual DDA cell visits across all pixels and samples, including
    /// zero-length crossings. Transparent missed rays consume no cell visits.
    pub max_cell_visits: u64,
}

/// A failed render returns no partial image or modified physical field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransferCameraError {
    /// Image allocation, sampling admission or source-radiance diagnostic.
    Image(DvrError),
    /// The original camera cannot generate this achromatic geometric ray.
    Camera(CameraError),
    /// The original volume ray refused, including cancellation and either the
    /// per-ray or remaining whole-image cell-crossing budget.
    Ray(TransferRayError),
}

impl From<DvrError> for TransferCameraError {
    fn from(error: DvrError) -> Self {
        Self::Image(error)
    }
}
impl From<CameraError> for TransferCameraError {
    fn from(error: CameraError) -> Self {
        Self::Camera(error)
    }
}
impl From<TransferRayError> for TransferCameraError {
    fn from(error: TransferRayError) -> Self {
        Self::Ray(error)
    }
}
impl core::fmt::Display for TransferCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scalar-transfer camera render refused: {self:?}")
    }
}
impl core::error::Error for TransferCameraError {}

/// Render a perspective or ideal finite-aperture camera through a borrowed
/// scalar-transfer field. Every primary ray uses exact native cell transport;
/// only pixel-area and pupil integration use sampling.
///
/// Samples average emitted linear RGB and transmission separately. Opacity
/// retains the mean extinguished fraction, including optically thin values;
/// effective depth is `-ln(mean transmission)`, never mean optical depth.
/// This preserves partial coverage and permits ordinary background composition.
/// The background supplied to a returned pixel must be constant over that
/// pixel/pupil footprint; a varying background needs per-ray transport to retain
/// its correlation with the intervening volume.
///
/// The camera's existing projection, focus and uniform pupil-area law determine
/// rays. These are normalized radiance/transmission averages, with no aperture
/// gain, exposure, tone mapping or sensor-irradiance conversion. Dispersive,
/// diffractive and sequential lenses retain their original camera refusal;
/// an RGB transfer function supplies no missing spectral optical information.
/// Physical sources remain the prescribed cellwise RGB transfer model, without
/// inferred LTE, in-scattering or surface/refraction transport.
pub fn render_transfer_camera(
    cx: &Cx<'_>,
    medium: &TransferRayGrid<'_, '_>,
    camera: &PhysicalCamera,
    settings: TransferCameraSettings,
) -> Result<Vec<TransferTransportPixel>, TransferCameraError> {
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
    // Refuse a camera requiring spectral/throughput-weighted rays before
    // allocating the private image. No sampling state is consumed here.
    camera.generate_ray_from_tangent_offsets(cx, 0.0, 0.0, centre_lens)?;
    let mut image = Vec::new();
    image
        .try_reserve_exact(pixels)
        .map_err(|_| DvrError::ImageAllocationRefused { pixels })?;
    let mut remaining_cells = settings.max_cell_visits;
    for py in 0..height {
        for px in 0..width {
            checkpoint(cx)?;
            let pixel_index = py * width + px;
            let mut stream = StreamKey {
                seed: settings.seed,
                kernel: 0x0303,
                tile: u32::try_from(pixel_index).map_err(|_| DvrError::PixelIdentityOverflow)?,
            }
            .stream();
            let mut average = PixelAverage::default();
            for _ in 0..settings.samples_per_pixel {
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
                let (sample, visits) = medium.transport_counted(
                    cx,
                    Ray {
                        origin: [primary.origin.x, primary.origin.y, primary.origin.z],
                        dir: [primary.dir.x, primary.dir.y, primary.dir.z],
                        t0: 0.0,
                        t1: f64::INFINITY,
                    },
                    remaining_cells,
                )?;
                remaining_cells -= u64::from(visits);
                average.add(sample);
            }
            image.push(average.finish()?);
        }
    }
    checkpoint(cx)?;
    Ok(image)
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), TransferCameraError> {
    cx.checkpoint()
        .map_err(|_| TransferRayError::Traversal(VolumeTransportError::Cancelled).into())
}

struct PixelAverage {
    samples: u32,
    emission: [Sum; 3],
    opacity: Sum,
    transmission: Sum,
    minimum_depth: f64,
    scaled_transmission: Sum,
}

impl Default for PixelAverage {
    fn default() -> Self {
        Self {
            samples: 0,
            emission: [Sum::default(); 3],
            opacity: Sum::default(),
            transmission: Sum::default(),
            minimum_depth: f64::INFINITY,
            scaled_transmission: Sum::default(),
        }
    }
}

impl PixelAverage {
    fn add(&mut self, sample: TransferTransportPixel) {
        self.samples += 1;
        let count = f64::from(self.samples);
        // Compensated online means cannot overflow merely because many finite
        // bright samples were added, and retain constant subnormal sources.
        for channel in 0..3 {
            let mean = &mut self.emission[channel];
            let increment = (sample.emission_rgb[channel] - mean.value) / count;
            mean.add(increment);
        }
        self.opacity
            .add((sample.opacity - self.opacity.value) / count);
        self.transmission
            .add((sample.transmittance - self.transmission.value) / count);
        let depth = sample.optical_depth;
        if depth.is_finite() {
            if depth < self.minimum_depth {
                let scale = det::exp(depth - self.minimum_depth);
                self.scaled_transmission.value *= scale;
                self.scaled_transmission.correction *= scale;
                self.minimum_depth = depth;
            }
            self.scaled_transmission
                .add(det::exp(self.minimum_depth - depth));
        }
    }

    fn finish(self) -> Result<TransferTransportPixel, DvrError> {
        let opacity = self.opacity.value;
        let transmittance = self.transmission.value;
        let optical_depth = if opacity <= 0.5 {
            // Compensate the rounded subtraction in log(1-opacity). This
            // retains thin samples even when transmission rounds to one.
            let transmission = 1.0 - opacity;
            -det::ln(transmission) - ((1.0 - transmission) - opacity) / transmission
        } else if self.scaled_transmission.value > 0.0 {
            // Log-mean-exp retains a finite depth when all unit-radiance
            // transmissions underflow, so bright backgrounds can still survive.
            self.minimum_depth - det::ln(self.scaled_transmission.value / f64::from(self.samples))
        } else {
            f64::INFINITY
        };
        let emission_rgb = self.emission.map(|sum| sum.value);
        if emission_rgb.iter().any(|v| !v.is_finite() || *v < 0.0)
            || !opacity.is_finite()
            || !(0.0..=1.0).contains(&opacity)
            || !transmittance.is_finite()
            || !(0.0..=1.0).contains(&transmittance)
        {
            return Err(DvrError::UnrepresentableRadiance);
        }
        Ok(TransferTransportPixel {
            emission_rgb,
            transmittance,
            opacity,
            optical_depth,
        })
    }
}

#[cfg(test)]
mod tests;
