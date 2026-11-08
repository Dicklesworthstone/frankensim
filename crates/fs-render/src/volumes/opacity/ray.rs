//! Exact direct transport for perspective, oblique and clipped volume rays.

use super::{DvrError, TransferFunction, TransferIntegral, TransferTransportPixel, VolumeGrid};
use crate::volumes::Ray;
use crate::volumes::transport::{CellRay, VolumeTransportError, validate_ray};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

/// Explicit admission and per-ray traversal limits. No optical-depth cutoff
/// or stochastic collision budget changes the physical transfer solution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TransferRayBudget {
    /// Maximum borrowed cells scanned once when the immutable view is admitted.
    pub max_grid_cells: usize,
    /// Maximum native cell crossings per ray, including zero-length crossings.
    pub max_steps_per_ray: u32,
}

/// Invalid input, unrepresentable geometry/radiance or cancelled/budgeted work
/// returns no partial ray result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransferRayError {
    /// Original scalar-grid admission or source-radiance diagnostic.
    Transfer(DvrError),
    /// Original native DDA diagnostic, including cancellation and traversal cap.
    Traversal(VolumeTransportError),
}

impl From<DvrError> for TransferRayError {
    fn from(error: DvrError) -> Self {
        Self::Transfer(error)
    }
}

impl From<VolumeTransportError> for TransferRayError {
    fn from(error: VolumeTransportError) -> Self {
        Self::Traversal(error)
    }
}

impl core::fmt::Display for TransferRayError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "scalar-transfer ray refused: {self:?}")
    }
}
impl core::error::Error for TransferRayError {}

/// Immutable, zero-copy scalar-transfer view for repeated arbitrary-ray queries.
///
/// Admission scans the borrowed field once; each ray allocates nothing and
/// visits only intersected cells. Scalars may be signed and zero need not mean
/// vacuum. Extinction is per world length, so the ray direction need not be unit.
///
/// The exact solution is for piecewise-constant cells and the supplied linear
/// RGB source radiance: dL/ds = -sigma_t L + sigma_t S. It includes direct-beam
/// extinction and that prescribed source only. Scattering redistribution, LTE
/// temperatures, refraction, continuous-field interpolation and pixel-area
/// filtering are not inferred from a visualization transfer function.
#[derive(Clone, Copy)]
pub struct TransferRayGrid<'g, 'a> {
    grid: &'g VolumeGrid<'a>,
    transfer: &'g TransferFunction,
    max_steps_per_ray: u32,
}

impl<'g, 'a> TransferRayGrid<'g, 'a> {
    /// Validate one immutable scalar snapshot and retain its original buffers.
    /// The single field scan checkpoints at least every 1024 scalar values.
    pub fn try_new(
        cx: &Cx<'_>,
        grid: &'g VolumeGrid<'a>,
        transfer: &'g TransferFunction,
        budget: TransferRayBudget,
    ) -> Result<Self, TransferRayError> {
        checkpoint(cx)?;
        if grid.data.len() > budget.max_grid_cells {
            return Err(DvrError::GridCellBudgetExceeded {
                requested: grid.data.len(),
                limit: budget.max_grid_cells,
            }
            .into());
        }
        if budget.max_steps_per_ray == 0 {
            return Err(VolumeTransportError::Budget.into());
        }
        crate::volumes::validate_transfer_grid_geometry(grid)?;
        let (lo, hi) = grid.bounds();
        // Native crossings need representable cell faces at both ends, not
        // just a representable whole-grid extent.
        if (0..3).any(|a| lo[a] + grid.cell[a] <= lo[a] || hi[a] - grid.cell[a] >= hi[a]) {
            return Err(DvrError::InvalidGrid.into());
        }
        for (chunk_index, chunk) in grid.data.chunks(1024).enumerate() {
            checkpoint(cx)?;
            if let Some(offset) = chunk.iter().position(|value| !value.is_finite()) {
                return Err(DvrError::NonFiniteGridValue {
                    index: chunk_index * 1024 + offset,
                }
                .into());
            }
        }
        checkpoint(cx)?;
        Ok(Self {
            grid,
            transfer,
            max_steps_per_ray: budget.max_steps_per_ray,
        })
    }

    /// Integrate `origin + t * dir` in increasing t, from finite `t0` to `t1`.
    /// `t1` may be positive infinity. The original bounds clip the segment;
    /// a miss or a zero-length segment returns transparent, zero-emission transfer.
    ///
    /// Cell contributions use their actual world length and camera order. No
    /// majorant, seed or stepping resolution is needed; work depends on cells
    /// crossed rather than opacity. An opaque foreground does not bypass the
    /// traversal limit or cancellation checks for the remaining segment.
    pub fn transport(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
    ) -> Result<TransferTransportPixel, TransferRayError> {
        checkpoint(cx)?;
        let (origin, direction, distance) = physical_ray(ray)?;
        let mut integral = TransferIntegral::default();
        if let Some(mut cells) = CellRay::clip(
            self.grid,
            origin,
            direction,
            distance,
            self.max_steps_per_ray,
        )? {
            while let Some(segment) = cells.next(cx)? {
                if segment.length_m > 0.0 {
                    let optical = self.transfer.sample_finite(self.grid.data[segment.cell]);
                    integral.add(optical, segment.length_m)?;
                }
            }
        }
        checkpoint(cx)?;
        Ok(integral.finish())
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), VolumeTransportError> {
    cx.checkpoint().map_err(|_| VolumeTransportError::Cancelled)
}

fn physical_ray(ray: Ray) -> Result<(Point3, Vec3, f64), VolumeTransportError> {
    if ray
        .origin
        .iter()
        .chain(ray.dir.iter())
        .any(|v| !v.is_finite())
        || !ray.t0.is_finite()
        || ray.t1.is_nan()
        || ray.t1 < ray.t0
    {
        return Err(VolumeTransportError::Ray);
    }
    // Scale before squaring so a finite nonzero direction remains usable even
    // when its unscaled squared norm would underflow or overflow.
    let scale = ray.dir.into_iter().fold(0.0_f64, |m, d| m.max(d.abs()));
    if scale == 0.0 {
        return Err(VolumeTransportError::Ray);
    }
    let scaled = ray.dir.map(|d| d / scale);
    let norm = scaled.iter().map(|d| d * d).sum::<f64>().sqrt();
    let direction = Vec3::new(scaled[0] / norm, scaled[1] / norm, scaled[2] / norm);
    let start = ray.at(ray.t0);
    let origin = Point3::new(start[0], start[1], start[2]);
    let span = ray.t1 - ray.t0;
    if ray.t1.is_finite() && !span.is_finite() {
        return Err(VolumeTransportError::Unrepresentable);
    }
    let distance = if span == f64::INFINITY {
        f64::INFINITY
    } else {
        super::optical_product(&[scale, norm, span])
    };
    if span > 0.0 && distance == 0.0 {
        return Err(VolumeTransportError::Unrepresentable);
    }
    validate_ray(origin, direction, distance)?;
    Ok((origin, direction, distance))
}

#[cfg(test)]
mod tests;
