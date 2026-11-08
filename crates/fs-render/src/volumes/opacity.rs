//! Exact opacity and source transport through borrowed scalar fields.
//!
//! Coefficients are constant within each cell. The solution is exact for that
//! representation (up to floating-point arithmetic), not for an unresolved
//! continuous field. Extinction removes light from the direct beam; a transfer
//! source is an explicitly supplied radiance, not a multiple-scattering solve.

use super::{
    AbsorptionImageSettings, DvrError, TransferFunction, TransferSample, VolumeGrid,
    pixel_centre_cell, validate_transfer_grid,
};
use fs_math::det;

mod area;
pub use area::{
    render_density_transmittance_box, render_density_transport_box, render_transfer_transport_box,
};

mod spectral;
pub use spectral::{
    DensityAbsorptionSpectrum, IndexedDensityAbsorptionSpectrum,
    IndexedSpectralAbsorptionImageSettings, MAX_OPACITY_WAVELENGTHS,
    SpectralAbsorptionImageSettings, SpectralOpacityImage, SpectralOpacitySample,
    render_density_spectral_box, render_density_spectral_indexed_box,
    render_density_spectral_indexed_transmittance_box, render_density_spectral_transmittance_box,
};

#[cfg(feature = "tracer")]
mod camera;
#[cfg(feature = "tracer")]
pub use camera::{
    PreintegratedCameraPixel, TransferCameraError, TransferCameraSettings,
    render_preintegrated_camera, render_transfer_camera,
};
#[cfg(feature = "tracer")]
mod ray;
#[cfg(feature = "tracer")]
pub(in crate::volumes) use ray::physical_ray;
#[cfg(feature = "tracer")]
pub use ray::{TransferRayBudget, TransferRayError, TransferRayGrid};

/// SI calibration for a borrowed mass-density field.
///
/// A stored value v represents v * density_unit_kg_m3 kg/m³. One grid length
/// represents length_unit_m metres. The supplied mass absorption coefficient
/// belongs to the caller's wavelength and material state; no optical data or
/// conversion from lattice units is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityAbsorption {
    /// Strictly positive density represented by one stored field unit, kg/m³.
    pub density_unit_kg_m3: f64,
    /// Strictly positive metres represented by one grid coordinate unit.
    pub length_unit_m: f64,
    /// Nonnegative mass absorption coefficient at the chosen wavelength, m²/kg.
    pub mass_absorption_m2_kg: f64,
}

impl DensityAbsorption {
    fn validate(self) -> Result<(), DvrError> {
        if !self.density_unit_kg_m3.is_finite()
            || self.density_unit_kg_m3 <= 0.0
            || !self.length_unit_m.is_finite()
            || self.length_unit_m <= 0.0
            || !self.mass_absorption_m2_kg.is_finite()
            || self.mass_absorption_m2_kg < 0.0
        {
            return Err(DvrError::InvalidDensityAbsorption);
        }
        Ok(())
    }
}

/// Direct transfer along a volume ray, column, or sampled camera pixel.
///
/// Emission is premultiplied by its physical attenuation. Opacity is evaluated
/// independently with expm1, so a thin nonzero layer remains measurable even
/// when transmittance rounds to one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransferTransportPixel {
    emission_rgb: [f64; 3],
    transmittance: f64,
    opacity: f64,
    optical_depth: f64,
}

impl TransferTransportPixel {
    /// Accumulated linear source radiance, before adding a background.
    #[must_use]
    pub const fn emission_rgb(self) -> [f64; 3] {
        self.emission_rgb
    }

    /// Fraction of a background beam transmitted through the complete segment.
    #[must_use]
    pub const fn transmittance(self) -> f64 {
        self.transmittance
    }

    /// Fraction extinguished from a direct beam, including thin-layer values.
    #[must_use]
    pub const fn opacity(self) -> f64 {
        self.opacity
    }

    /// Integrated dimensionless extinction along a single ray. For a sampled
    /// camera pixel this is the effective `-ln(mean transmission)`, not the mean
    /// ray depth. Infinite means complete extinction in every retained sample.
    #[must_use]
    pub const fn optical_depth(self) -> f64 {
        self.optical_depth
    }

    /// Place this volume in front of nonnegative linear background radiance.
    ///
    /// This is L_emitted + T * L_background, without display gamma, a second
    /// multiplication by alpha, or a conversion of scattering into heat.
    pub fn composite(self, background_rgb: [f64; 3]) -> Result<[f64; 3], DvrError> {
        if background_rgb.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(DvrError::InvalidBackgroundRadiance);
        }
        let result = std::array::from_fn(|channel| {
            let background = background_rgb[channel];
            if self.transmittance < f64::MIN_POSITIVE && background > 0.0 {
                // Retain a bright, but still finite, background behind a depth
                // whose unit-radiance transmittance underflows.
                det::exp(det::ln(background) - self.optical_depth) + self.emission_rgb[channel]
            } else {
                self.transmittance
                    .mul_add(background, self.emission_rgb[channel])
            }
        });
        if result.iter().any(|v| !v.is_finite()) {
            return Err(DvrError::UnrepresentableRadiance);
        }
        Ok(result)
    }
}

/// Exact monochromatic absorption of a mass-density grid without copying it.
///
/// Optical depth is rho * kappa_a * ds in SI units. World translations and
/// pixel sampling do not alter the chosen grid column. Zero absorption is
/// transparent even at large density; scattering and refraction are separate
/// physical models. Work is bounded by the admitted grid cells and cell visits.
pub fn render_density_transmittance(
    grid: &VolumeGrid<'_>,
    settings: AbsorptionImageSettings,
    calibration: DensityAbsorption,
) -> Result<Vec<f64>, DvrError> {
    calibration.validate()?;
    let pixels = admit_image(grid, settings)?;
    if let Some(index) = grid.data.iter().position(|value| *value < 0.0) {
        return Err(DvrError::NegativeGridAbsorption { index });
    }
    let mut image = reserve_image(pixels)?;
    let [width, height] = settings.resolution;
    let [nx, ny, _] = grid.dims;
    for py in 0..height {
        let iy = pixel_centre_cell(py, height, ny);
        for px in 0..width {
            let column = pixel_centre_cell(px, width, nx) + nx * iy;
            image.push(det::exp(-density_column_depth(grid, column, calibration)));
        }
    }
    Ok(image)
}

/// Exact orthographic transfer-function emission, opacity and transmittance.
///
/// Rays look along -z. Every scalar cell is mapped through the supplied
/// transfer function and integrated analytically over its actual length:
/// dL = T_front * source * (1 - exp(-extinction * length)).
/// There is no collision sampling, global majorant, opacity cutoff or seed.
/// Cost is proportional to intersected cells, independent of optical depth.
/// Pixel centres select columns in integer grid coordinates, so translating
/// the field cannot change face ownership. Lateral pixel-area filtering and
/// continuous-field preintegration are outside this piecewise-constant model.
pub fn render_transfer_transport(
    grid: &VolumeGrid<'_>,
    transfer: &TransferFunction,
    settings: AbsorptionImageSettings,
) -> Result<Vec<TransferTransportPixel>, DvrError> {
    let pixels = admit_image(grid, settings)?;
    let mut image = reserve_image(pixels)?;
    let [width, height] = settings.resolution;
    let [nx, ny, _] = grid.dims;
    for py in 0..height {
        let iy = pixel_centre_cell(py, height, ny);
        for px in 0..width {
            let column = pixel_centre_cell(px, width, nx) + nx * iy;
            image.push(transfer_column(grid, column, transfer)?);
        }
    }
    Ok(image)
}

fn density_column(
    grid: &VolumeGrid<'_>,
    column: usize,
    calibration: DensityAbsorption,
) -> TransferTransportPixel {
    let depth = density_column_depth(grid, column, calibration);
    TransferTransportPixel {
        emission_rgb: [0.0; 3],
        transmittance: det::exp(-depth),
        opacity: -det::expm1(-depth),
        optical_depth: depth,
    }
}

fn density_column_depth(
    grid: &VolumeGrid<'_>,
    column: usize,
    calibration: DensityAbsorption,
) -> f64 {
    let plane = grid.dims[0] * grid.dims[1];
    let mut depth = Sum::default();
    for iz in (0..grid.dims[2]).rev() {
        depth.add(optical_product(&[
            grid.data[column + plane * iz],
            grid.cell[2],
            calibration.density_unit_kg_m3,
            calibration.length_unit_m,
            calibration.mass_absorption_m2_kg,
        ]));
    }
    depth.value
}

fn transfer_column(
    grid: &VolumeGrid<'_>,
    column: usize,
    transfer: &TransferFunction,
) -> Result<TransferTransportPixel, DvrError> {
    let plane = grid.dims[0] * grid.dims[1];
    let mut integral = TransferIntegral::default();
    for iz in (0..grid.dims[2]).rev() {
        let optical = transfer.sample_finite(grid.data[column + plane * iz]);
        integral.add(optical, grid.cell[2])?;
    }
    Ok(integral.finish())
}

// The orthographic specialization and arbitrary-ray DDA use exactly the same
// ordered cell solution. No source/opacity buffer is copied or pretabulated.
#[derive(Default)]
struct TransferIntegral {
    depth: Sum,
    emission: [Sum; 3],
}

impl TransferIntegral {
    fn add(&mut self, optical: TransferSample, length: f64) -> Result<(), DvrError> {
        let tau = optical_product(&[optical.extinction, length]);
        let alpha = -det::expm1(-tau);
        let transmission = det::exp(-self.depth.value);
        for channel in 0..3 {
            let source = optical.source_rgb[channel];
            let mut contribution = (source * alpha) * transmission;
            if source > 0.0 && optical.extinction > 0.0 && length > 0.0 {
                if tau < f64::MIN_POSITIVE {
                    // In this limit 1-exp(-tau) = tau to binary64 accuracy.
                    // A bright source can survive even when sigma*ds itself
                    // underflows; retain its original positive factors.
                    contribution = det::exp(
                        det::ln(source) + det::ln(optical.extinction) + det::ln(length)
                            - self.depth.value,
                    );
                } else if transmission < f64::MIN_POSITIVE || contribution == 0.0 {
                    // A bright source can survive an underflowed foreground
                    // transmission without changing ordinary-range arithmetic.
                    contribution = det::exp(det::ln(source) + det::ln(alpha) - self.depth.value);
                }
            }
            self.emission[channel].add(contribution);
            if !self.emission[channel].value.is_finite() {
                return Err(DvrError::UnrepresentableRadiance);
            }
        }
        self.depth.add(tau);
        Ok(())
    }

    fn finish(self) -> TransferTransportPixel {
        TransferTransportPixel {
            emission_rgb: self.emission.map(|sum| sum.value),
            transmittance: det::exp(-self.depth.value),
            opacity: -det::expm1(-self.depth.value),
            optical_depth: self.depth.value,
        }
    }
}

// Average transported radiance and survival, never density or optical depth.
// Integer weights cover both camera samples and exact grid/pixel overlap areas.
// Callers admit the total weight before using this allocation-free accumulator.
struct TransportAverage {
    weight: u64,
    emission: [Sum; 3],
    opacity: Sum,
    transmission: Sum,
    minimum_depth: f64,
    scaled_transmission: Sum,
}

impl Default for TransportAverage {
    fn default() -> Self {
        Self {
            weight: 0,
            emission: [Sum::default(); 3],
            opacity: Sum::default(),
            transmission: Sum::default(),
            minimum_depth: f64::INFINITY,
            scaled_transmission: Sum::default(),
        }
    }
}

impl TransportAverage {
    fn add(&mut self, sample: TransferTransportPixel, weight: u64) {
        self.weight += weight;
        let count = self.weight as f64;
        let fraction = weight as f64 / count;
        let increment = |difference: f64| {
            // Retain the existing camera's ordinary-range arithmetic.
            if weight == 1 {
                difference / count
            } else {
                difference * fraction
            }
        };
        for channel in 0..3 {
            let mean = &mut self.emission[channel];
            mean.add(increment(sample.emission_rgb[channel] - mean.value));
        }
        self.add_survival(sample, weight);
    }

    // Spectral absorption has no emission channels to average. Reuse the same
    // thin-opacity and logarithmic survival arithmetic without that RGB work.
    fn add_absorption(&mut self, depth: f64, weight: u64) {
        self.weight += weight;
        self.add_survival(
            TransferTransportPixel {
                emission_rgb: [0.0; 3],
                transmittance: det::exp(-depth),
                opacity: -det::expm1(-depth),
                optical_depth: depth,
            },
            weight,
        );
    }

    fn add_survival(&mut self, sample: TransferTransportPixel, weight: u64) {
        let count = self.weight as f64;
        let fraction = weight as f64 / count;
        let increment = |difference: f64| {
            if weight == 1 {
                difference / count
            } else {
                difference * fraction
            }
        };
        self.opacity
            .add(increment(sample.opacity - self.opacity.value));
        self.transmission
            .add(increment(sample.transmittance - self.transmission.value));
        let depth = sample.optical_depth;
        if depth.is_finite() {
            if depth < self.minimum_depth {
                let scale = det::exp(depth - self.minimum_depth);
                self.scaled_transmission.value *= scale;
                self.scaled_transmission.correction *= scale;
                self.minimum_depth = depth;
            }
            self.scaled_transmission
                .add(weight as f64 * det::exp(self.minimum_depth - depth));
        }
    }

    fn finish(self) -> Result<TransferTransportPixel, DvrError> {
        let opacity = self.opacity.value;
        let transmittance = self.transmission.value;
        let optical_depth = if opacity <= 0.5 {
            let transmission = 1.0 - opacity;
            -det::ln(transmission) - ((1.0 - transmission) - opacity) / transmission
        } else if self.scaled_transmission.value > 0.0 {
            self.minimum_depth - det::ln(self.scaled_transmission.value / self.weight as f64)
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

fn reserve_image<T>(pixels: usize) -> Result<Vec<T>, DvrError> {
    let mut image = Vec::new();
    image
        .try_reserve_exact(pixels)
        .map_err(|_| DvrError::ImageAllocationRefused { pixels })?;
    Ok(image)
}

fn admit_image(
    grid: &VolumeGrid<'_>,
    settings: AbsorptionImageSettings,
) -> Result<usize, DvrError> {
    let [width, height] = settings.resolution;
    let pixels = width
        .checked_mul(height)
        .filter(|_| width > 0 && height > 0)
        .ok_or(DvrError::InvalidImageSize)?;
    if pixels > settings.max_pixels {
        return Err(DvrError::PixelBudgetExceeded {
            requested: pixels,
            limit: settings.max_pixels,
        });
    }
    if grid.data.len() > settings.max_grid_cells {
        return Err(DvrError::GridCellBudgetExceeded {
            requested: grid.data.len(),
            limit: settings.max_grid_cells,
        });
    }
    let visits = u64::try_from(pixels)
        .ok()
        .zip(u64::try_from(grid.dims[2]).ok())
        .and_then(|(pixels, depth)| pixels.checked_mul(depth))
        .ok_or(DvrError::CellVisitCountOverflow)?;
    if visits > settings.max_cell_visits {
        return Err(DvrError::CellVisitBudgetExceeded {
            requested: visits,
            limit: settings.max_cell_visits,
        });
    }
    validate_transfer_grid(grid)?;
    Ok(pixels)
}

// Compensated sums preserve many optically thin cells. Overflow of a sum of
// nonnegative optical depths means complete extinction, never NaN/gain.
#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, value: f64) {
        if self.value.is_infinite() {
            return;
        }
        let corrected = value - self.correction;
        let next = self.value + corrected;
        self.correction = if next.is_finite() {
            (next - self.value) - corrected
        } else {
            0.0
        };
        self.value = next;
    }
}

// The common path uses ordinary products. On range loss, combine logarithms
// before exponentiating: a large density and a small mass coefficient may
// have a perfectly finite optical depth. All factors were admitted nonnegative.
fn optical_product(factors: &[f64]) -> f64 {
    if factors.contains(&0.0) {
        return 0.0;
    }
    let mut product = 1.0;
    let mut ordinary = true;
    for &factor in factors {
        product *= factor;
        ordinary &= product.is_finite() && product >= f64::MIN_POSITIVE;
    }
    if ordinary {
        product
    } else {
        det::exp(factors.iter().map(|value| det::ln(*value)).sum())
    }
}

#[cfg(test)]
mod tests;
