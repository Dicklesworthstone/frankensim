//! Controlled continuous-field transport along arbitrary physical rays.
//!
//! A trilinear cell restricts to a cubic Bézier scalar curve along a ray.
//! Convex-hull and transfer-Lipschitz bounds admit linear chords before the
//! existing knot-resolving preintegrator handles them. This is an adaptive
//! approximation to that cubic, not exact linear stepping through nodal data.

use super::{
    DvrError, NodeVolumeGrid, PreintegratedPixel, PreintegrationError, PreintegrationSettings,
    RayIntegral, Sum, TransferFunction, Work, integrate_segment, lerp,
};
use crate::volumes::Ray;
use crate::volumes::opacity::physical_ray;
use crate::volumes::transport::{CellRay, CellRayGeometry, VolumeTransportError};
use fs_exec::Cx;
use fs_geom::Point3;
use fs_math::det;

const MAX_CURVE_DEPTH: usize = 40;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NodeRayBudget {
    /// Maximum original nodes scanned at admission; no scalar buffer is copied.
    pub max_grid_nodes: usize,
    /// Maximum native cell crossings per ray, including zero-length crossings.
    pub max_steps_per_ray: u32,
}

/// Absolute reconstruction targets and hard work limits for one complete ray.
/// Zero reconstruction tolerance admits only chords with zero model bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreintegratedRaySettings {
    /// Maximum dimensionless extinction-integral error from cubic-to-chord
    /// approximation, excluding floating-point and source-data uncertainty.
    pub optical_depth_tolerance: f64,
    /// Maximum emitted-radiance reconstruction error in each source-RGB unit.
    /// Includes changed foreground attenuation of downstream sources; excludes
    /// the separately returned quadrature estimates and floating-point error.
    pub source_reconstruction_tolerance_rgb: [f64; 3],
    /// Maximum cubic bisections across the complete ray. Stack depth is also
    /// limited to 40; exceeding either limit refuses rather than relaxes accuracy.
    pub max_curve_subdivisions: usize,
    /// Existing source quadrature tolerance and complete-ray work limits.
    /// Transfer-interval work includes coefficient-bound scans as well as
    /// accepted-chord preintegration; neither budget resets at a cell crossing.
    pub integration: PreintegrationSettings,
}

/// Transport through the admitted chords, with the two error sources separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreintegratedRayPixel {
    /// Integrated radiance/depth and existing Gauss source-error estimates.
    pub transport: PreintegratedPixel,
    /// Analytic reconstruction bound on the extinction integral. This excludes
    /// floating-point geometry/arithmetic and uncertainty in the original data.
    pub optical_depth_error: f64,
    /// Reconstruction bound on final emitted radiance, including optical-depth
    /// uncertainty in front of every source. Add the transport's quadrature
    /// estimates separately; their sum is not an interval certificate.
    pub source_reconstruction_error_rgb: [f64; 3],
}

#[derive(Default)]
pub(in crate::volumes) struct RayWork {
    pub cells: u32,
    pub intervals: usize,
    pub quadrature_evaluations: u64,
    pub subdivisions: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PreintegratedRayError {
    /// Original transfer/linear-preintegration refusal.
    Integration(PreintegrationError),
    /// Original geometry, ray, cancellation or cell-crossing-budget refusal.
    Traversal(VolumeTransportError),
    /// Reconstruction tolerances must be finite and nonnegative.
    InvalidReconstructionTolerance,
    /// Another cubic bisection would exceed the complete-ray limit.
    SubdivisionBudgetExceeded,
    /// The bounded stack, physical subdivision or representable error target
    /// was exhausted without meeting the requested reconstruction accuracy.
    SubdivisionResolutionLimit,
}

impl From<PreintegrationError> for PreintegratedRayError {
    fn from(error: PreintegrationError) -> Self {
        match error {
            PreintegrationError::Cancelled => Self::Traversal(VolumeTransportError::Cancelled),
            error => Self::Integration(error),
        }
    }
}
impl From<DvrError> for PreintegratedRayError {
    fn from(error: DvrError) -> Self {
        Self::Integration(PreintegrationError::Transfer(error))
    }
}
impl From<VolumeTransportError> for PreintegratedRayError {
    fn from(error: VolumeTransportError) -> Self {
        Self::Traversal(error)
    }
}
impl core::fmt::Display for PreintegratedRayError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous volume ray refused: {self:?}")
    }
}
impl core::error::Error for PreintegratedRayError {}

/// Borrowed node-centered field with its original trilinear reconstruction.
/// Admission scans nodes once; ray traversal and cubic subdivision allocate
/// nothing. Length and extinction units are the original grid/transfer units.
pub struct PreintegratedRayGrid<'g, 'a> {
    grid: &'g NodeVolumeGrid<'a>,
    transfer: &'g TransferFunction,
    geometry: CellRayGeometry,
    maximum_source_rgb: [f64; 3],
    max_steps_per_ray: u32,
}

impl<'g, 'a> PreintegratedRayGrid<'g, 'a> {
    pub fn try_new(
        cx: &Cx<'_>,
        grid: &'g NodeVolumeGrid<'a>,
        transfer: &'g TransferFunction,
        budget: NodeRayBudget,
    ) -> Result<Self, PreintegratedRayError> {
        checkpoint(cx)?;
        if grid.data.len() > budget.max_grid_nodes {
            return Err(DvrError::GridCellBudgetExceeded {
                requested: grid.data.len(),
                limit: budget.max_grid_nodes,
            }
            .into());
        }
        if budget.max_steps_per_ray == 0 {
            return Err(VolumeTransportError::Budget.into());
        }
        if grid.dims.iter().any(|n| *n < 2)
            || grid.dims.iter().try_fold(1usize, |n, &d| n.checked_mul(d)) != Some(grid.data.len())
        {
            return Err(DvrError::InvalidGrid.into());
        }
        let geometry = CellRayGeometry {
            dims: grid.dims.map(|n| n - 1),
            origin: grid.origin,
            cell: grid.spacing,
        };
        for a in 0..3 {
            let upper = (geometry.dims[a] as f64).mul_add(grid.spacing[a], grid.origin[a]);
            if !grid.origin[a].is_finite()
                || !grid.spacing[a].is_finite()
                || grid.spacing[a] <= 0.0
                || !upper.is_finite()
                || !(upper - grid.origin[a]).is_finite()
                || upper <= grid.origin[a]
                || grid.origin[a] + grid.spacing[a] <= grid.origin[a]
                || upper - grid.spacing[a] >= upper
            {
                return Err(DvrError::InvalidGrid.into());
            }
        }
        let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
        for (chunk_index, chunk) in grid.data.chunks(1024).enumerate() {
            checkpoint(cx)?;
            for (offset, &scalar) in chunk.iter().enumerate() {
                if !scalar.is_finite() {
                    return Err(DvrError::NonFiniteGridValue {
                        index: chunk_index * 1024 + offset,
                    }
                    .into());
                }
                low = low.min(scalar);
                high = high.max(scalar);
            }
        }
        let front = transfer.sample_finite(low).source_rgb;
        let back = transfer.sample_finite(high).source_rgb;
        let mut maximum_source_rgb = std::array::from_fn(|c| front[c].max(back[c]));
        let points = transfer.points();
        let first = points.partition_point(|p| p.scalar < low);
        let last = points.partition_point(|p| p.scalar <= high);
        for chunk in points[first..last].chunks(1024) {
            checkpoint(cx)?;
            for point in chunk {
                for c in 0..3 {
                    maximum_source_rgb[c] = maximum_source_rgb[c].max(point.source_rgb[c]);
                }
            }
        }
        checkpoint(cx)?;
        Ok(Self {
            grid,
            transfer,
            geometry,
            maximum_source_rgb,
            max_steps_per_ray: budget.max_steps_per_ray,
        })
    }

    /// Trace the original node reconstruction in camera order. Non-unit ray
    /// parameters and clipping retain the same physical meaning as cell rays.
    /// Each cubic's convex hull bounds all crossed transfer intervals, including
    /// narrow opacity peaks invisible at its two endpoints. Accuracy is never
    /// inferred from midpoint agreement or from an already-classified color.
    pub fn transport(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        settings: PreintegratedRaySettings,
    ) -> Result<PreintegratedRayPixel, PreintegratedRayError> {
        self.transport_counted(cx, ray, settings, u64::MAX)
            .map(|(pixel, _)| pixel)
    }

    pub(in crate::volumes) fn transport_counted(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        settings: PreintegratedRaySettings,
        remaining_cells: u64,
    ) -> Result<(PreintegratedRayPixel, RayWork), PreintegratedRayError> {
        checkpoint(cx)?;
        if !settings.optical_depth_tolerance.is_finite()
            || settings.optical_depth_tolerance < 0.0
            || settings
                .source_reconstruction_tolerance_rgb
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(PreintegratedRayError::InvalidReconstructionTolerance);
        }
        let mut poll = || cx.checkpoint().map_err(|_| PreintegrationError::Cancelled);
        let mut work = Work::with_checkpoint(settings.integration, &mut poll)?;
        let (origin, direction, distance) = physical_ray(ray)?;
        let mut integral = RayIntegral::default();
        let max_steps = self
            .max_steps_per_ray
            .min(u32::try_from(remaining_cells).unwrap_or(u32::MAX));
        let Some(mut cells) =
            CellRay::clip_geometry(self.geometry, origin, direction, distance, max_steps)?
        else {
            return Ok((
                PreintegratedRayPixel {
                    transport: integral.finish(),
                    optical_depth_error: 0.0,
                    source_reconstruction_error_rgb: [0.0; 3],
                },
                RayWork::default(),
            ));
        };
        let mut depth_tolerance = settings.optical_depth_tolerance;
        for c in 0..3 {
            if self.maximum_source_rgb[c] > 0.0 {
                depth_tolerance = depth_tolerance.min(
                    0.5 * (settings.source_reconstruction_tolerance_rgb[c]
                        / self.maximum_source_rgb[c]),
                );
            }
        }
        let coefficient_tolerance = depth_tolerance / cells.clipped_length();
        let mut depth_error = Sum::default();
        let mut maximum_source_deviation = [0.0_f64; 3];
        let mut subdivisions = 0;
        let mut cell_visits = 0;
        let mut stack = [Curve {
            control: [0.0; 4],
            length: 0.0,
            depth: 0,
        }; MAX_CURVE_DEPTH + 1];
        while let Some(segment) = cells.next(cx)? {
            cell_visits += 1;
            if segment.length_m == 0.0 {
                continue;
            }
            let control = self.cell_curve(
                segment.cell,
                cells.position(segment.start_m),
                cells.position(segment.start_m + segment.length_m),
            );
            stack[0] = Curve {
                control,
                length: segment.length_m,
                depth: 0,
            };
            let mut pending = 1;
            while pending > 0 {
                checkpoint(cx)?;
                pending -= 1;
                let piece = stack[pending];
                let bounds = coefficient_bounds(self.transfer, piece.control, &mut work)?;
                if bounds.extinction <= coefficient_tolerance
                    && (0..3).all(|c| {
                        bounds.source[c] <= 0.5 * settings.source_reconstruction_tolerance_rgb[c]
                    })
                {
                    integrate_segment(
                        self.transfer,
                        piece.control[0],
                        piece.control[3],
                        piece.length,
                        &mut work,
                        &mut integral,
                    )?;
                    depth_error.add(positive_bound_product(bounds.extinction, piece.length));
                    for c in 0..3 {
                        maximum_source_deviation[c] =
                            maximum_source_deviation[c].max(bounds.source[c]);
                    }
                } else {
                    if subdivisions >= settings.max_curve_subdivisions {
                        return Err(PreintegratedRayError::SubdivisionBudgetExceeded);
                    }
                    if piece.depth == MAX_CURVE_DEPTH || 0.5 * piece.length == 0.0 {
                        return Err(PreintegratedRayError::SubdivisionResolutionLimit);
                    }
                    subdivisions += 1;
                    let (front, back) = piece.split();
                    stack[pending] = back;
                    stack[pending + 1] = front;
                    pending += 2;
                }
            }
        }
        checkpoint(cx)?;
        // Variation of dL/ds = sigma*(S-L): |delta L| <= max|delta S|
        // + max(S)*integral|delta sigma|. Both S and downstream L lie in
        // [0,max(S)], so this includes foreground attenuation of every later
        // source, rather than merely adding local unattenuated source errors.
        let source_reconstruction_error_rgb = std::array::from_fn(|c| {
            (maximum_source_deviation[c]
                + positive_bound_product(self.maximum_source_rgb[c], depth_error.value))
            .min(self.maximum_source_rgb[c])
        });
        if !depth_error.value.is_finite()
            || depth_error.value > settings.optical_depth_tolerance
            || (0..3).any(|c| {
                source_reconstruction_error_rgb[c] > settings.source_reconstruction_tolerance_rgb[c]
            })
        {
            return Err(PreintegratedRayError::SubdivisionResolutionLimit);
        }
        Ok((
            PreintegratedRayPixel {
                transport: integral.finish(),
                optical_depth_error: depth_error.value,
                source_reconstruction_error_rgb,
            },
            RayWork {
                cells: cell_visits,
                intervals: work.intervals,
                quadrature_evaluations: work.evaluations,
                subdivisions,
            },
        ))
    }

    fn cell_curve(&self, cell: usize, front: Point3, back: Point3) -> [f64; 4] {
        let nx = self.geometry.dims[0];
        let ny = self.geometry.dims[1];
        let index = [cell % nx, (cell / nx) % ny, cell / (nx * ny)];
        let local = |point: Point3| {
            let p = [point.x, point.y, point.z];
            std::array::from_fn::<_, 3, _>(|a| {
                let lower = (index[a] as f64).mul_add(self.grid.spacing[a], self.grid.origin[a]);
                ((p[a] - lower) / self.grid.spacing[a]).clamp(0.0, 1.0)
            })
        };
        let a = local(front);
        let b = local(back);
        let scalar = |p: [f64; 3]| {
            lerp(
                self.grid
                    .plane_scalar(index[0], index[1], index[2], p[0], p[1]),
                self.grid
                    .plane_scalar(index[0], index[1], index[2] + 1, p[0], p[1]),
                p[2],
            )
        };
        // Multiaffine blossom: each control is a convex combination of node
        // values. Derivative differences and a power-basis conversion are
        // unnecessary, including on translated or anisotropic grids.
        [
            scalar(a),
            mean_three(
                scalar([b[0], a[1], a[2]]),
                scalar([a[0], b[1], a[2]]),
                scalar([a[0], a[1], b[2]]),
            ),
            mean_three(
                scalar([a[0], b[1], b[2]]),
                scalar([b[0], a[1], b[2]]),
                scalar([b[0], b[1], a[2]]),
            ),
            scalar(b),
        ]
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), PreintegratedRayError> {
    cx.checkpoint()
        .map_err(|_| VolumeTransportError::Cancelled.into())
}

fn mean_three(a: f64, b: f64, c: f64) -> f64 {
    lerp(lerp(a, b, 0.5), c, 1.0 / 3.0)
}

#[derive(Clone, Copy)]
struct Curve {
    control: [f64; 4],
    length: f64,
    depth: usize,
}

impl Curve {
    fn split(self) -> (Self, Self) {
        let [a, b, c, d] = self.control;
        let ab = lerp(a, b, 0.5);
        let bc = lerp(b, c, 0.5);
        let cd = lerp(c, d, 0.5);
        let abc = lerp(ab, bc, 0.5);
        let bcd = lerp(bc, cd, 0.5);
        let middle = lerp(abc, bcd, 0.5);
        let half = 0.5 * self.length;
        (
            Self {
                control: [a, ab, abc, middle],
                length: half,
                depth: self.depth + 1,
            },
            Self {
                control: [middle, bcd, cd, d],
                length: self.length - half,
                depth: self.depth + 1,
            },
        )
    }
}

struct CoefficientBounds {
    extinction: f64,
    source: [f64; 3],
}

fn coefficient_bounds(
    transfer: &TransferFunction,
    curve: [f64; 4],
    work: &mut Work<'_>,
) -> Result<CoefficientBounds, PreintegrationError> {
    let deviation = (curve[1] - lerp(curve[0], curve[3], 1.0 / 3.0))
        .abs()
        .max((curve[2] - lerp(curve[0], curve[3], 2.0 / 3.0)).abs());
    if deviation == 0.0 {
        return Ok(CoefficientBounds {
            extinction: 0.0,
            source: [0.0; 3],
        });
    }
    let low = curve.into_iter().fold(f64::INFINITY, f64::min);
    let high = curve.into_iter().fold(f64::NEG_INFINITY, f64::max);
    let points = transfer.points();
    let mut knot = points.partition_point(|p| p.scalar <= low);
    let mut scalar = low;
    let mut optical = transfer.sample_finite(scalar);
    let (mut min_extinction, mut max_extinction) = (optical.extinction, optical.extinction);
    let (mut min_source, mut max_source) = (optical.source_rgb, optical.source_rgb);
    let mut bounds = CoefficientBounds {
        extinction: 0.0,
        source: [0.0; 3],
    };
    loop {
        work.interval()?;
        let next = if knot < points.len() && points[knot].scalar < high {
            let scalar = points[knot].scalar;
            knot += 1;
            scalar
        } else {
            high
        };
        let next_optical = transfer.sample_finite(next);
        let width = next - scalar;
        bounds.extinction = bounds.extinction.max(lipschitz_bound(
            (next_optical.extinction - optical.extinction).abs(),
            deviation,
            width,
        ));
        min_extinction = min_extinction.min(next_optical.extinction);
        max_extinction = max_extinction.max(next_optical.extinction);
        for c in 0..3 {
            bounds.source[c] = bounds.source[c].max(lipschitz_bound(
                (next_optical.source_rgb[c] - optical.source_rgb[c]).abs(),
                deviation,
                width,
            ));
            min_source[c] = min_source[c].min(next_optical.source_rgb[c]);
            max_source[c] = max_source[c].max(next_optical.source_rgb[c]);
        }
        if next == high {
            break;
        }
        scalar = next;
        optical = next_optical;
    }
    bounds.extinction = bounds.extinction.min(max_extinction - min_extinction);
    for c in 0..3 {
        bounds.source[c] = if max_extinction == 0.0 {
            0.0
        } else {
            bounds.source[c].min(max_source[c] - min_source[c])
        };
    }
    Ok(bounds)
}

fn lipschitz_bound(change: f64, deviation: f64, width: f64) -> f64 {
    if change == 0.0 || deviation == 0.0 {
        return 0.0;
    }
    if !deviation.is_finite() || !width.is_finite() {
        return f64::INFINITY;
    }
    let direct = change * (deviation / width);
    if direct > 0.0 && direct.is_finite() {
        return direct;
    }
    det::exp(det::ln(change) + det::ln(deviation) - det::ln(width)).max(f64::from_bits(1))
}

fn positive_bound_product(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        0.0
    } else {
        (a * b).max(f64::from_bits(1))
    }
}

#[cfg(test)]
mod tests;
