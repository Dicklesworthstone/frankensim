//! Exact orthographic box footprints on the original piecewise-constant field.

use super::{
    AbsorptionImageSettings, DensityAbsorption, DvrError, TransferFunction,
    TransferTransportPixel, TransportAverage, VolumeGrid, density_column, reserve_image,
    transfer_column, validate_transfer_grid,
};

/// Area-averaged SI absorption, with stable opacity and effective optical depth.
///
/// Each pixel covers its complete rectangular footprint in the grid's XY
/// extent; parallel rays face -Z. Original columns are integrated before their
/// transmissions are averaged. Thus a half-clear, half-opaque pixel transmits
/// half its background, rather than becoming opaque at the mean density.
/// Integer overlap areas avoid coordinate rounding and require no supersampling.
///
/// The original cell-constant model is retained: subcell interfaces, scattering,
/// refraction, and unknown material coefficients are not inferred. Composition
/// assumes a background constant over each footprint. Output is x-fastest,
/// starting at the minimum Y face, as in the existing orthographic functions.
///
/// No converted field or full-resolution image is allocated. Cell work is
/// exactly `(nx + width - gcd(nx, width)) * (ny + height - gcd(ny, height)) * nz`;
/// it must fit `max_cell_visits` before scanning data or allocating output.
pub fn render_density_transport_box(
    grid: &VolumeGrid<'_>,
    settings: AbsorptionImageSettings,
    calibration: DensityAbsorption,
) -> Result<Vec<TransferTransportPixel>, DvrError> {
    render_density_box(grid, settings, calibration, |pixel| pixel)
}

/// The same physical box filter, retaining only direct transmittance.
///
/// This writes the scalar output directly, without allocating a temporary
/// radiance image. Use [`render_density_transport_box`] when thin opacity or
/// finite log attenuation behind an underflowed transmission is needed.
pub fn render_density_transmittance_box(
    grid: &VolumeGrid<'_>,
    settings: AbsorptionImageSettings,
    calibration: DensityAbsorption,
) -> Result<Vec<f64>, DvrError> {
    render_density_box(
        grid, settings, calibration, TransferTransportPixel::transmittance,
    )
}

fn render_density_box<P>(
    grid: &VolumeGrid<'_>,
    settings: AbsorptionImageSettings,
    calibration: DensityAbsorption,
    finish: impl Fn(TransferTransportPixel) -> P,
) -> Result<Vec<P>, DvrError> {
    calibration.validate()?;
    let pixels = admit(grid, settings)?;
    if let Some(index) = grid.data.iter().position(|value| *value < 0.0) {
        return Err(DvrError::NegativeGridAbsorption { index });
    }
    render_box(
        grid, settings, pixels,
        |column| Ok(density_column(grid, column, calibration)), finish,
    )
}

/// Exact box-filtered transfer emission and direct-beam opacity.
///
/// Reuses the original ordered cell solution on each column. Premultiplied
/// emitted radiance and transmission are averaged independently, preserving
/// their correlation with the original layers. Classifying an averaged scalar
/// or averaging optical depth would solve a different physical problem.
/// Footprint, work, layout and background assumptions are identical to
/// [`render_density_transport_box`]. No source or opacity cutoff is used.
pub fn render_transfer_transport_box(
    grid: &VolumeGrid<'_>,
    transfer: &TransferFunction,
    settings: AbsorptionImageSettings,
) -> Result<Vec<TransferTransportPixel>, DvrError> {
    let pixels = admit(grid, settings)?;
    render_box(
        grid, settings, pixels,
        |column| transfer_column(grid, column, transfer), |pixel| pixel,
    )
}

fn render_box<P>(
    grid: &VolumeGrid<'_>,
    settings: AbsorptionImageSettings,
    pixels: usize,
    mut column: impl FnMut(usize) -> Result<TransferTransportPixel, DvrError>,
    finish: impl Fn(TransferTransportPixel) -> P,
) -> Result<Vec<P>, DvrError> {
    let [width, height] = settings.resolution;
    let [nx, ny, _] = grid.dims;
    let mut image = reserve_image(pixels)?;
    for py in 0..height {
        for px in 0..width {
            let mut average = TransportAverage::default();
            for (iy, wy) in overlaps(py, height, ny) {
                for (ix, wx) in overlaps(px, width, nx) {
                    // Grid admission proves nx*ny fits usize and u64. These
                    // disjoint integer areas sum exactly to nx*ny per pixel.
                    let weight = (wx * wy) as u64;
                    average.add(column(ix + nx * iy)?, weight);
                }
            }
            image.push(finish(average.finish()?));
        }
    }
    Ok(image)
}

fn admit(grid: &VolumeGrid<'_>, settings: AbsorptionImageSettings) -> Result<usize, DvrError> {
    let [width, height] = settings.resolution;
    let pixels = width
        .checked_mul(height)
        .filter(|_| width > 0 && height > 0)
        .ok_or(DvrError::InvalidImageSize)?;
    if pixels > settings.max_pixels {
        return Err(DvrError::PixelBudgetExceeded {
            requested: pixels, limit: settings.max_pixels,
        });
    }
    if grid.data.len() > settings.max_grid_cells {
        return Err(DvrError::GridCellBudgetExceeded {
            requested: grid.data.len(), limit: settings.max_grid_cells,
        });
    }
    let [nx, ny, nz] = grid.dims;
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(DvrError::InvalidGrid);
    }
    let visits = overlap_count(nx, width)
        .checked_mul(overlap_count(ny, height))
        .and_then(|columns| columns.checked_mul(nz as u128))
        .and_then(|visits| u64::try_from(visits).ok())
        .ok_or(DvrError::CellVisitCountOverflow)?;
    if visits > settings.max_cell_visits {
        return Err(DvrError::CellVisitBudgetExceeded {
            requested: visits, limit: settings.max_cell_visits,
        });
    }
    validate_transfer_grid(grid)?;
    Ok(pixels)
}

fn overlap_count(cells: usize, pixels: usize) -> u128 {
    let (mut a, mut b) = (cells, pixels);
    while b != 0 {
        (a, b) = (b, a % b);
    }
    cells as u128 + pixels as u128 - a as u128
}

// Scale normalized coordinates by cells*pixels. Cell faces are multiples of
// pixels; pixel faces are multiples of cells. Every overlap is then an exact
// integer, including nondivisible resolutions and very large world origins.
fn overlaps(pixel: usize, pixels: usize, cells: usize) -> impl Iterator<Item = (usize, usize)> {
    let left = pixel as u128 * cells as u128;
    let right = (pixel as u128 + 1) * cells as u128;
    let stride = pixels as u128;
    let first = (left / stride) as usize;
    let last = ((right - 1) / stride) as usize;
    (first..=last).map(move |cell| {
        let start = (cell as u128 * stride).max(left);
        let end = ((cell as u128 + 1) * stride).min(right);
        (cell, (end - start) as usize)
    })
}

#[cfg(test)]
mod tests;
