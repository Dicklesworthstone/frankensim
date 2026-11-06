//! Exact opacity and source transport through borrowed scalar-field columns.
//!
//! Coefficients are constant within each cell. The solution is exact for that
//! representation (up to floating-point arithmetic), not for an unresolved
//! continuous field. Extinction removes light from the direct beam; a transfer
//! source is an explicitly supplied radiance, not a multiple-scattering solve.

use super::{
    AbsorptionImageSettings, DvrError, TransferFunction, VolumeGrid, pixel_centre_cell,
    validate_transfer_grid,
};
use fs_math::det;

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

/// Exact direct transfer for one orthographic volume pixel.
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

    /// Fraction of a background beam transmitted through the complete column.
    #[must_use]
    pub const fn transmittance(self) -> f64 {
        self.transmittance
    }

    /// Fraction extinguished from a direct beam, including thin-layer values.
    #[must_use]
    pub const fn opacity(self) -> f64 {
        self.opacity
    }

    /// Integrated dimensionless extinction; infinite means complete extinction.
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
                self.transmittance.mul_add(background, self.emission_rgb[channel])
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
    if !calibration.density_unit_kg_m3.is_finite()
        || calibration.density_unit_kg_m3 <= 0.0
        || !calibration.length_unit_m.is_finite()
        || calibration.length_unit_m <= 0.0
        || !calibration.mass_absorption_m2_kg.is_finite()
        || calibration.mass_absorption_m2_kg < 0.0
    {
        return Err(DvrError::InvalidDensityAbsorption);
    }
    let pixels = admit_image(grid, settings)?;
    if let Some(index) = grid.data.iter().position(|value| *value < 0.0) {
        return Err(DvrError::NegativeGridAbsorption { index });
    }
    let mut image = reserve_image(pixels)?;
    let [width, height] = settings.resolution;
    let [nx, ny, nz] = grid.dims;
    let plane = nx * ny;
    for py in 0..height {
        let iy = pixel_centre_cell(py, height, ny);
        for px in 0..width {
            let column = pixel_centre_cell(px, width, nx) + nx * iy;
            let mut depth = Sum::default();
            for iz in (0..nz).rev() {
                depth.add(optical_product(&[
                    grid.data[column + plane * iz],
                    grid.cell[2],
                    calibration.density_unit_kg_m3,
                    calibration.length_unit_m,
                    calibration.mass_absorption_m2_kg,
                ]));
            }
            image.push(det::exp(-depth.value));
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
    let [nx, ny, nz] = grid.dims;
    let plane = nx * ny;
    for py in 0..height {
        let iy = pixel_centre_cell(py, height, ny);
        for px in 0..width {
            let column = pixel_centre_cell(px, width, nx) + nx * iy;
            let mut depth = Sum::default();
            let mut emission = [Sum::default(); 3];
            for iz in (0..nz).rev() {
                let optical = transfer.sample_finite(grid.data[column + plane * iz]);
                let tau = optical_product(&[optical.extinction, grid.cell[2]]);
                let alpha = -det::expm1(-tau);
                let transmission = det::exp(-depth.value);
                for channel in 0..3 {
                    let source = optical.source_rgb[channel];
                    let mut contribution = (source * alpha) * transmission;
                    // A bright source behind a thick foreground can remain
                    // representable after exp(-depth) itself has underflowed.
                    if (transmission < f64::MIN_POSITIVE || contribution == 0.0)
                        && source > 0.0
                        && alpha > 0.0
                    {
                        contribution =
                            det::exp(det::ln(source) + det::ln(alpha) - depth.value);
                    }
                    emission[channel].add(contribution);
                    if !emission[channel].value.is_finite() {
                        return Err(DvrError::UnrepresentableRadiance);
                    }
                }
                depth.add(tau);
            }
            image.push(TransferTransportPixel {
                emission_rgb: emission.map(|sum| sum.value),
                transmittance: det::exp(-depth.value),
                opacity: -det::expm1(-depth.value),
                optical_depth: depth.value,
            });
        }
    }
    Ok(image)
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
