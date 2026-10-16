//! Certified converter SDF → mesh via DUAL CONTOURING with QEF vertex
//! placement (plan §7.3 edge 2): sharp-feature capable, with THE
//! certificate — an exact-distance enclosure verification that the extracted
//! surface lies within tolerance of the zero set everywhere, not just at
//! samples. This is the certificate's one-sided output-surface proximity
//! claim; it does not prove that every zero-set point is near the mesh.
//!
//! The certificate's honesty: v1 accepts only a chart that advertises the
//! GLOBAL [`TraceStepClaim::ExactDistance`] theorem and supplies a rigorous
//! finite enclosure of every centroid evaluation. Then
//! `distance(x, zero_set) <= max_abs(enclosure(c)) + |x-c|` over the whole
//! triangle. Subdividing until that bound closes below tolerance proves the
//! claim; local `ChartSample::lipschitz` values never mint global authority.

use crate::winding::Soup;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{
    Aabb, Chart, ClippedChart, Point3, SamplingDomain, SamplingDomainError, TraceStepClaim, Vec3,
};
use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Vertex-placement strategy: QEF preserves sharp features; `MassPoint`
/// is the marching-cubes-class baseline the acceptance criteria compare
/// against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Placement {
    /// Quadratic-error-function minimization (sharp features).
    Qef,
    /// Edge-crossing centroid (feature-blurring baseline).
    MassPoint,
}

/// Dual-contouring options.
#[derive(Debug, Clone, Copy)]
pub struct DcOptions {
    /// Cell edge length.
    pub h: f64,
    /// Vertex placement strategy.
    pub placement: Placement,
    /// QEF regularization toward the mass point (Schaefer-style; keeps
    /// near-planar systems well-posed without SVD).
    pub regularization: f64,
}

impl DcOptions {
    /// Sharp-feature defaults at cell size `h`.
    #[must_use]
    pub fn sharp(h: f64) -> Self {
        DcOptions {
            h,
            placement: Placement::Qef,
            regularization: 0.05,
        }
    }
}

/// Contouring phase that consumed a chart sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContourSampleStage {
    /// Uniform corner-lattice sampling.
    CornerLattice,
    /// Sign-changing edge refinement.
    EdgeRefinement,
    /// Gradient query at an extracted crossing.
    CrossingGradient,
    /// Finite-difference fallback for a missing chart gradient.
    GradientProbe,
}

/// Finite arithmetic stage that could not be represented.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContourArithmeticStage {
    /// Lookup or construction of a uniform lattice point.
    LatticePlacement,
    /// Convex interpolation of a sign-changing edge.
    EdgeInterpolation,
    /// Reuse of an extracted crossing for a gradient query.
    CrossingGradient,
    /// Construction of a finite-difference probe point.
    GradientProbe,
    /// Normalization of a chart or finite-difference gradient.
    GradientNormalization,
    /// Stable averaging of Hermite crossing points.
    MassPoint,
    /// Assembly or solution of the regularized QEF.
    Qef,
}

/// Structured contouring failure.
#[derive(Debug, Clone, PartialEq)]
pub enum ContourError {
    /// The source support or explicit clip is not an admissible finite
    /// three-dimensional sampling domain.
    SamplingDomain(SamplingDomainError),
    /// The requested cell spacing was not finite and strictly positive.
    InvalidSpacing {
        /// Which constructor argument was invalid.
        field: &'static str,
        /// The offending value.
        value: f64,
    },
    /// QEF regularization must be finite and non-negative.
    InvalidRegularization {
        /// Offending value.
        value: f64,
    },
    /// The grid would exceed the per-axis cell cap.
    ResolutionTooFine {
        /// Cells/axis needed.
        need: u64,
        /// The cap.
        cap: u64,
    },
    /// Checked multiplication of the corner-lattice dimensions overflowed
    /// the addressable allocation domain.
    GridSizeOverflow {
        /// Corner-lattice dimensions that could not be multiplied safely.
        dims: [u64; 3],
    },
    /// A normalized lattice coordinate could not be represented even though
    /// its sampling domain had passed admission.
    CoordinatePlacementOverflow {
        /// Coordinate axis.
        axis: &'static str,
        /// Zero-based lattice-node index.
        index: u32,
        /// Number of cells along the axis.
        cells: u32,
    },
    /// A finite chart query point could not be constructed.
    NonRepresentablePoint {
        /// Arithmetic phase that constructed the point.
        stage: ContourArithmeticStage,
        /// Offending point.
        point: Point3,
    },
    /// A chart returned a NaN or infinite signed-field value.
    InvalidSample {
        /// Sampling phase.
        stage: ContourSampleStage,
        /// Query point.
        point: Point3,
        /// Offending nominal value.
        value: f64,
    },
    /// A chart or finite-difference fallback produced an unusable gradient.
    InvalidGradient {
        /// Crossing point whose normal was requested.
        point: Point3,
        /// Offending gradient.
        gradient: Vec3,
    },
    /// Finite inputs overflowed during a derived contour calculation.
    NonRepresentableArithmetic {
        /// Arithmetic phase that refused.
        stage: ContourArithmeticStage,
    },
    /// Contouring observed cancellation at a bounded polling point.
    Cancelled,
    /// The zero set never crossed the sampled grid.
    EmptySurface,
}

impl core::fmt::Display for ContourError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ContourError::SamplingDomain(error) => write!(f, "{error}"),
            ContourError::InvalidSpacing { field, value } => write!(
                f,
                "dual contouring refused: `{field}` must be finite and strictly positive, got {value}"
            ),
            ContourError::ResolutionTooFine { need, cap } => write!(
                f,
                "dual contouring refused: {need} cells/axis exceed the {cap} cap; coarsen h \
                 or shrink the region"
            ),
            ContourError::InvalidRegularization { value } => write!(
                f,
                "dual contouring refused: QEF regularization must be finite and non-negative, got {value}"
            ),
            ContourError::GridSizeOverflow { dims } => write!(
                f,
                "dual contouring refused: corner-lattice dimensions {dims:?} overflow the addressable grid"
            ),
            ContourError::CoordinatePlacementOverflow { axis, index, cells } => write!(
                f,
                "dual contouring refused: normalized {axis}-axis lattice coordinate {index}/{cells} is not representable"
            ),
            ContourError::NonRepresentablePoint { stage, point } => write!(
                f,
                "dual contouring refused: {stage:?} produced a non-representable point {point:?}"
            ),
            ContourError::InvalidSample {
                stage,
                point,
                value,
            } => write!(
                f,
                "dual contouring refused: chart returned non-finite value {value} during {stage:?} at {point:?}"
            ),
            ContourError::InvalidGradient { point, gradient } => write!(
                f,
                "dual contouring refused: unusable gradient {gradient:?} at {point:?}"
            ),
            ContourError::NonRepresentableArithmetic { stage } => write!(
                f,
                "dual contouring refused: finite arithmetic overflowed during {stage:?}"
            ),
            ContourError::Cancelled => {
                write!(f, "dual contouring cancelled at a bounded polling point")
            }
            ContourError::EmptySurface => write!(
                f,
                "dual contouring found no zero crossings: the chart's zero set does not \
                 intersect the sampled support (empty or out-of-band geometry)"
            ),
        }
    }
}

impl core::error::Error for ContourError {}

impl From<SamplingDomainError> for ContourError {
    fn from(error: SamplingDomainError) -> Self {
        Self::SamplingDomain(error)
    }
}

/// Per-axis cell cap.
pub const DC_MAX_CELLS_PER_AXIS: u64 = 256;

fn finite_positive(value: f64, field: &'static str) -> Result<f64, ContourError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(ContourError::InvalidSpacing { field, value })
    }
}

fn finite_nonnegative_regularization(value: f64) -> Result<f64, ContourError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(ContourError::InvalidRegularization { value })
    }
}

fn finite_point(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

fn sample_signed_distance(
    chart: &dyn Chart,
    point: Point3,
    stage: ContourSampleStage,
    cx: &Cx<'_>,
) -> Result<f64, ContourError> {
    if !finite_point(point) {
        return Err(ContourError::NonRepresentablePoint {
            stage: match stage {
                ContourSampleStage::CornerLattice => ContourArithmeticStage::LatticePlacement,
                ContourSampleStage::EdgeRefinement => ContourArithmeticStage::EdgeInterpolation,
                ContourSampleStage::CrossingGradient => ContourArithmeticStage::CrossingGradient,
                ContourSampleStage::GradientProbe => ContourArithmeticStage::GradientProbe,
            },
            point,
        });
    }
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    let value = chart.eval(point, cx).signed_distance;
    // The chart itself need not poll the gate. Observe producer-requested
    // cancellation before consuming or publishing its sample.
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ContourError::InvalidSample {
            stage,
            point,
            value,
        })
    }
}

fn checked_axis_nodes(span: f64, h: f64) -> Result<u32, ContourError> {
    let rounded_cells = (span / h).ceil();
    if !rounded_cells.is_finite()
        || rounded_cells < 0.0
        || rounded_cells > DC_MAX_CELLS_PER_AXIS as f64
    {
        let need = if rounded_cells.is_finite() && rounded_cells >= 0.0 {
            rounded_cells as u64
        } else {
            u64::MAX
        };
        return Err(ContourError::ResolutionTooFine {
            need,
            cap: DC_MAX_CELLS_PER_AXIS,
        });
    }
    // A positive representable span is a real cell even when `span / h`
    // underflows to zero. Never admit a one-node/zero-cell axis.
    let mut cells = rounded_cells.max(1.0) as u64;
    loop {
        let width_cells = u32::try_from(cells).map_err(|_| ContourError::ResolutionTooFine {
            need: cells,
            cap: DC_MAX_CELLS_PER_AXIS,
        })?;
        let realized_width = span / f64::from(width_cells);
        if !realized_width.is_finite() {
            return Err(ContourError::NonRepresentableArithmetic {
                stage: ContourArithmeticStage::LatticePlacement,
            });
        }
        if realized_width <= h {
            break;
        }
        cells = cells
            .checked_add(1)
            .ok_or(ContourError::ResolutionTooFine {
                need: u64::MAX,
                cap: DC_MAX_CELLS_PER_AXIS,
            })?;
        if cells > DC_MAX_CELLS_PER_AXIS {
            return Err(ContourError::ResolutionTooFine {
                need: cells,
                cap: DC_MAX_CELLS_PER_AXIS,
            });
        }
    }
    let nodes = cells.checked_add(1).ok_or(ContourError::GridSizeOverflow {
        dims: [u64::MAX; 3],
    })?;
    u32::try_from(nodes).map_err(|_| ContourError::GridSizeOverflow { dims: [nodes; 3] })
}

fn checked_grid_size(dims: [u32; 3]) -> Result<usize, ContourError> {
    let dims_u64 = dims.map(u64::from);
    let total = dims_u64
        .iter()
        .try_fold(1u128, |product, &dim| product.checked_mul(u128::from(dim)))
        .ok_or(ContourError::GridSizeOverflow { dims: dims_u64 })?;
    usize::try_from(total).map_err(|_| ContourError::GridSizeOverflow { dims: dims_u64 })
}

fn axis_coordinates(
    min: f64,
    max: f64,
    span: f64,
    nodes: u32,
    axis: &'static str,
) -> Result<Vec<f64>, ContourError> {
    let cells = nodes
        .checked_sub(1)
        .ok_or(ContourError::CoordinatePlacementOverflow {
            axis,
            index: 0,
            cells: 0,
        })?;
    let capacity = usize::try_from(nodes).map_err(|_| ContourError::GridSizeOverflow {
        dims: [u64::from(nodes); 3],
    })?;
    let mut coordinates = Vec::with_capacity(capacity);
    for index in 0..nodes {
        let coordinate = if index == 0 {
            min
        } else if index == cells {
            max
        } else {
            let t = f64::from(index) / f64::from(cells);
            (min + span * t).clamp(min, max)
        };
        if !coordinate.is_finite()
            || coordinates
                .last()
                .is_some_and(|previous| coordinate <= *previous)
        {
            return Err(ContourError::CoordinatePlacementOverflow { axis, index, cells });
        }
        coordinates.push(coordinate);
    }
    Ok(coordinates)
}

/// Contouring statistics (ledgered evidence).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DcStats {
    /// Cells carrying a vertex.
    pub active_cells: u64,
    /// Output triangles.
    pub triangles: u64,
    /// Hermite edge crossings found.
    pub crossings: u64,
}

impl DcStats {
    /// Canonical JSON.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(72);
        let _ = write!(
            s,
            "{{\"active_cells\":{},\"triangles\":{},\"crossings\":{}}}",
            self.active_cells, self.triangles, self.crossings
        );
        s
    }
}

/// Dual-contour any signed-distance chart on a uniform grid.
///
/// # Errors
/// [`ContourError`] refusals, including direct cancellation checkpoints in
/// every major contouring phase.
// One coherent pass (sample -> hermite -> place -> stitch); splitting
// would scatter the winding/orientation invariants rmesh-008 audits.
#[allow(clippy::too_many_lines)]
pub fn dual_contour(
    chart: &dyn Chart,
    opts: DcOptions,
    cx: &Cx<'_>,
) -> Result<(Soup, DcStats), ContourError> {
    let h = finite_positive(opts.h, "h")?;
    let regularization = finite_nonnegative_regularization(opts.regularization)?;
    let padding = finite_positive(2.0 * h, "2 * h")?;
    let raw_support = SamplingDomain::admit(chart.support(), None)?.bounds();
    let domain = SamplingDomain::admit(raw_support.inflate(padding), None)?;
    let support = domain.bounds();
    let spans = domain.spans();
    let n = [
        checked_axis_nodes(spans.x, h)?,
        checked_axis_nodes(spans.y, h)?,
        checked_axis_nodes(spans.z, h)?,
    ];
    let grid_size = checked_grid_size(n)?;
    let x = axis_coordinates(support.min.x, support.max.x, spans.x, n[0], "x")?;
    let y = axis_coordinates(support.min.y, support.max.y, spans.y, n[1], "y")?;
    let z = axis_coordinates(support.min.z, support.max.z, spans.z, n[2], "z")?;
    let pos = |i: u32, j: u32, k: u32| {
        // The public cap is 256 cells/axis, so these conversions are total on
        // every supported Rust target. `map_or` keeps the invariant panic-free.
        let i = usize::try_from(i).unwrap_or_default();
        let j = usize::try_from(j).unwrap_or_default();
        let k = usize::try_from(k).unwrap_or_default();
        Point3::new(x[i], y[j], z[k])
    };
    // Sample the corner lattice once.
    let idx = |i: u32, j: u32, k: u32| {
        (usize::try_from(k).expect("u32 fits usize")
            * usize::try_from(n[1]).expect("u32 fits usize")
            + usize::try_from(j).expect("u32 fits usize"))
            * usize::try_from(n[0]).expect("u32 fits usize")
            + usize::try_from(i).expect("u32 fits usize")
    };
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    let mut phi = vec![0.0f64; grid_size];
    for k in 0..n[2] {
        for j in 0..n[1] {
            cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
            for i in 0..n[0] {
                let point = pos(i, j, k);
                phi[idx(i, j, k)] =
                    sample_signed_distance(chart, point, ContourSampleStage::CornerLattice, cx)?;
            }
        }
    }
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    // Hermite data per sign-changing lattice edge; cell -> crossings map.
    let mut cell_hermite: BTreeMap<[u32; 3], Vec<(Point3, Vec3)>> = BTreeMap::new();
    let mut crossings = 0u64;
    let axes: [[u32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for k in 0..n[2] {
        for j in 0..n[1] {
            cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
            for i in 0..n[0] {
                for d in axes {
                    let (i2, j2, k2) = (i + d[0], j + d[1], k + d[2]);
                    if i2 >= n[0] || j2 >= n[1] || k2 >= n[2] {
                        continue;
                    }
                    let (fa, fb) = (phi[idx(i, j, k)], phi[idx(i2, j2, k2)]);
                    if (fa < 0.0) == (fb < 0.0) {
                        continue;
                    }
                    crossings += 1;
                    let (pa, pb) = (pos(i, j, k), pos(i2, j2, k2));
                    let crossing = secant_crossing(chart, pa, fa, pb, fb, cx)?;
                    let normal = gradient_at(chart, crossing, h, cx)?;
                    // Every VALID adjacent cell gets the Hermite pair
                    // (boundary edges feed fewer cells; their cells still
                    // place vertices even though no quad is emitted).
                    for (du, dv) in [(0u32, 0u32), (1, 0), (1, 1), (0, 1)] {
                        let (u, v) = match d {
                            [1, 0, 0] => ([0u32, 1, 0], [0u32, 0, 1]),
                            [0, 1, 0] => ([0u32, 0, 1], [1u32, 0, 0]),
                            _ => ([1u32, 0, 0], [0u32, 1, 0]),
                        };
                        let cell = [
                            i.wrapping_sub(du * u[0]).wrapping_sub(dv * v[0]),
                            j.wrapping_sub(du * u[1]).wrapping_sub(dv * v[1]),
                            k.wrapping_sub(du * u[2]).wrapping_sub(dv * v[2]),
                        ];
                        if cell[0] < n[0] - 1 && cell[1] < n[1] - 1 && cell[2] < n[2] - 1 {
                            cell_hermite
                                .entry(cell)
                                .or_default()
                                .push((crossing, normal));
                        }
                    }
                }
            }
        }
    }
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    if cell_hermite.is_empty() {
        return Err(ContourError::EmptySurface);
    }
    // One vertex per active cell.
    let mut vertex_of: BTreeMap<[u32; 3], u32> = BTreeMap::new();
    let mut positions: Vec<Point3> = Vec::with_capacity(cell_hermite.len());
    for (&cell, hermite) in &cell_hermite {
        cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
        let cell_min = pos(cell[0], cell[1], cell[2]);
        let cell_max = pos(cell[0] + 1, cell[1] + 1, cell[2] + 1);
        let v = match opts.placement {
            Placement::MassPoint => mass_point(hermite)?,
            Placement::Qef => solve_qef(hermite, regularization, cell_min, cell_max)?,
        };
        vertex_of.insert(cell, positions.len() as u32);
        positions.push(v);
    }
    // Quads per interior sign-changing edge, oriented negative -> positive.
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
            for i in 0..n[0] {
                for d in axes {
                    let (i2, j2, k2) = (i + d[0], j + d[1], k + d[2]);
                    if i2 >= n[0] || j2 >= n[1] || k2 >= n[2] {
                        continue;
                    }
                    let (fa, fb) = (phi[idx(i, j, k)], phi[idx(i2, j2, k2)]);
                    if (fa < 0.0) == (fb < 0.0) {
                        continue;
                    }
                    let Some(ring) = edge_cells([i, j, k], d, n) else {
                        continue; // boundary edge: no full quad (open rim)
                    };
                    let q: Vec<u32> = ring.iter().map(|c| vertex_of[c]).collect();
                    // The ring circulates with normal +d; outward normals
                    // point from negative to positive phi, so fa < 0
                    // (inside at base) keeps the ring, else reverse.
                    if fa < 0.0 {
                        triangles.push([q[0], q[1], q[2]]);
                        triangles.push([q[0], q[2], q[3]]);
                    } else {
                        triangles.push([q[0], q[2], q[1]]);
                        triangles.push([q[0], q[3], q[2]]);
                    }
                }
            }
        }
    }
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    let stats = DcStats {
        active_cells: positions.len() as u64,
        triangles: triangles.len() as u64,
        crossings,
    };
    Ok((
        Soup {
            positions,
            triangles,
        },
        stats,
    ))
}

/// Dual-contour the geometric intersection `chart ∩ clip` on a uniform grid.
/// The explicit clip is part of the sampled field, not merely a replacement
/// sampling extent.
///
/// # Errors
/// [`ContourError`] when the clip, spacing, or resulting grid is inadmissible,
/// or when the clipped zero set has no sampled crossing.
pub fn dual_contour_clipped(
    chart: &dyn Chart,
    opts: DcOptions,
    clip: Aabb,
    cx: &Cx<'_>,
) -> Result<(Soup, DcStats), ContourError> {
    finite_positive(opts.h, "h")?;
    let clipped = ClippedChart::new(chart, clip)?;
    dual_contour(&clipped, opts, cx)
}

/// The four cells sharing lattice edge `(base, base+d)` in RING order
/// whose circulation normal points along +d (cyclic (u, v) axes make the
/// orientation axis-uniform); `None` when the edge sits on the lattice
/// boundary (no full quad — the open rim).
fn edge_cells(base: [u32; 3], d: [u32; 3], n: [u32; 3]) -> Option<[[u32; 3]; 4]> {
    // Cyclic successors keep u x v = +d for every axis.
    let (u, v) = match d {
        [1, 0, 0] => ([0u32, 1, 0], [0u32, 0, 1]), // (y, z)
        [0, 1, 0] => ([0u32, 0, 1], [1u32, 0, 0]), // (z, x)
        _ => ([1u32, 0, 0], [0u32, 1, 0]),         // (x, y)
    };
    // CCW ring viewed from +d: offsets (0,0), (1,0), (1,1), (0,1).
    let mut ring = [[0u32; 3]; 4];
    for (slot, (du, dv)) in [(0u32, 0u32), (1, 0), (1, 1), (0, 1)]
        .into_iter()
        .enumerate()
    {
        let c = [
            base[0].wrapping_sub(du * u[0]).wrapping_sub(dv * v[0]),
            base[1].wrapping_sub(du * u[1]).wrapping_sub(dv * v[1]),
            base[2].wrapping_sub(du * u[2]).wrapping_sub(dv * v[2]),
        ];
        // Cells are corner-indexed: valid when strictly inside the lattice.
        if !(c[0] < n[0] - 1 && c[1] < n[1] - 1 && c[2] < n[2] - 1) {
            return None;
        }
        ring[slot] = c;
    }
    Some(ring)
}

fn secant_crossing(
    chart: &dyn Chart,
    mut a: Point3,
    mut fa: f64,
    mut b: Point3,
    mut fb: f64,
    cx: &Cx<'_>,
) -> Result<Point3, ContourError> {
    for _ in 0..8 {
        let abs_a = fa.abs();
        let abs_b = fb.abs();
        let scale = abs_a.max(abs_b);
        let t = if scale == 0.0 {
            0.5
        } else {
            let scaled_a = abs_a / scale;
            let scaled_b = abs_b / scale;
            scaled_a / (scaled_a + scaled_b)
        };
        let m = convex_point(a, b, t).ok_or(ContourError::NonRepresentablePoint {
            stage: ContourArithmeticStage::EdgeInterpolation,
            point: Point3::new(f64::NAN, f64::NAN, f64::NAN),
        })?;
        let fm = sample_signed_distance(chart, m, ContourSampleStage::EdgeRefinement, cx)?;
        if fm.abs() < 1e-12 {
            return Ok(m);
        }
        if (fm < 0.0) == (fa < 0.0) {
            a = m;
            fa = fm;
        } else {
            b = m;
            fb = fm;
        }
    }
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    stable_midpoint(a, b).ok_or(ContourError::NonRepresentablePoint {
        stage: ContourArithmeticStage::EdgeInterpolation,
        point: Point3::new(
            f64::midpoint(a.x, b.x),
            f64::midpoint(a.y, b.y),
            f64::midpoint(a.z, b.z),
        ),
    })
}

fn gradient_at(chart: &dyn Chart, p: Point3, h: f64, cx: &Cx<'_>) -> Result<Vec3, ContourError> {
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    let sample = chart.eval(p, cx);
    cx.checkpoint().map_err(|_| ContourError::Cancelled)?;
    if !sample.signed_distance.is_finite() {
        return Err(ContourError::InvalidSample {
            stage: ContourSampleStage::CrossingGradient,
            point: p,
            value: sample.signed_distance,
        });
    }
    if let Some(gradient) = sample.gradient {
        return normalize_gradient(gradient, p);
    }
    let e = 0.1 * h;
    if !e.is_finite() || e <= 0.0 {
        return Err(ContourError::NonRepresentableArithmetic {
            stage: ContourArithmeticStage::GradientProbe,
        });
    }
    let probes = [
        Point3::new(p.x + e, p.y, p.z),
        Point3::new(p.x - e, p.y, p.z),
        Point3::new(p.x, p.y + e, p.z),
        Point3::new(p.x, p.y - e, p.z),
        Point3::new(p.x, p.y, p.z + e),
        Point3::new(p.x, p.y, p.z - e),
    ];
    let mut values = [0.0; 6];
    for (slot, point) in values.iter_mut().zip(probes) {
        if !finite_point(point) {
            return Err(ContourError::NonRepresentablePoint {
                stage: ContourArithmeticStage::GradientProbe,
                point,
            });
        }
        *slot = sample_signed_distance(chart, point, ContourSampleStage::GradientProbe, cx)?;
    }
    let scale = values
        .iter()
        .fold(0.0f64, |acc, value| acc.max(value.abs()));
    if !scale.is_finite() || scale == 0.0 {
        return Err(ContourError::InvalidGradient {
            point: p,
            gradient: Vec3::new(0.0, 0.0, 0.0),
        });
    }
    // The common `2e` denominator does not affect direction. Scaling before
    // subtraction prevents opposite finite samples from overflowing.
    let gradient = Vec3::new(
        values[0] / scale - values[1] / scale,
        values[2] / scale - values[3] / scale,
        values[4] / scale - values[5] / scale,
    );
    normalize_gradient(gradient, p)
}

fn normalize_gradient(gradient: Vec3, point: Point3) -> Result<Vec3, ContourError> {
    let scale = gradient.x.abs().max(gradient.y.abs()).max(gradient.z.abs());
    if !scale.is_finite() || scale == 0.0 {
        return Err(ContourError::InvalidGradient { point, gradient });
    }
    let scaled = gradient.scale(1.0 / scale);
    let norm = scaled.norm();
    if !norm.is_finite() || norm == 0.0 {
        return Err(ContourError::InvalidGradient { point, gradient });
    }
    let normalized = scaled.scale(1.0 / norm);
    if normalized.x.is_finite() && normalized.y.is_finite() && normalized.z.is_finite() {
        Ok(normalized)
    } else {
        Err(ContourError::NonRepresentableArithmetic {
            stage: ContourArithmeticStage::GradientNormalization,
        })
    }
}

fn convex_point(a: Point3, b: Point3, t: f64) -> Option<Point3> {
    let interpolate = |lhs: f64, rhs: f64| lhs.mul_add(1.0 - t, rhs * t);
    let point = Point3::new(
        interpolate(a.x, b.x),
        interpolate(a.y, b.y),
        interpolate(a.z, b.z),
    );
    finite_point(point).then_some(point)
}

fn mass_point(hermite: &[(Point3, Vec3)]) -> Result<Point3, ContourError> {
    let Some(&(first, _)) = hermite.first() else {
        return Err(ContourError::NonRepresentableArithmetic {
            stage: ContourArithmeticStage::MassPoint,
        });
    };
    let mut mean = first;
    for (index, &(point, _)) in hermite.iter().enumerate().skip(1) {
        #[allow(clippy::cast_precision_loss)]
        let weight = 1.0 / (index + 1) as f64;
        mean =
            convex_point(mean, point, weight).ok_or(ContourError::NonRepresentableArithmetic {
                stage: ContourArithmeticStage::MassPoint,
            })?;
    }
    Ok(mean)
}

/// Regularized 3×3 QEF solve: minimize Σ(nᵢ·(x−pᵢ))² + λ|x−m|² where m is
/// the mass point (Schaefer-style regularization keeps near-planar
/// systems well-posed); the result clamps into the cell.
fn solve_qef(
    hermite: &[(Point3, Vec3)],
    lambda: f64,
    cell_min: Point3,
    cell_max: Point3,
) -> Result<Point3, ContourError> {
    let m = mass_point(hermite)?;
    // Solve in cell-local coordinates `x = m + y`. This avoids multiplying
    // normals by huge world-space offsets: A y = Σ n(n·(p-m)), while the
    // regularizer λ|x-m|² contributes λI and a zero local right-hand side.
    let mut a = [[0.0f64; 3]; 3];
    let mut b = [0.0f64; 3];
    for &(p, nrm) in hermite {
        let nv = [nrm.x, nrm.y, nrm.z];
        let local = p.delta_from(m);
        if !(local.x.is_finite() && local.y.is_finite() && local.z.is_finite()) {
            return Err(ContourError::NonRepresentableArithmetic {
                stage: ContourArithmeticStage::Qef,
            });
        }
        let nd = nrm.dot(local);
        if !nd.is_finite() {
            return Err(ContourError::NonRepresentableArithmetic {
                stage: ContourArithmeticStage::Qef,
            });
        }
        for (r, (row, rhs)) in a.iter_mut().zip(&mut b).enumerate() {
            for (entry, nc) in row.iter_mut().zip(nv) {
                *entry += nv[r] * nc;
                if !entry.is_finite() {
                    return Err(ContourError::NonRepresentableArithmetic {
                        stage: ContourArithmeticStage::Qef,
                    });
                }
            }
            *rhs += nv[r] * nd;
            if !rhs.is_finite() {
                return Err(ContourError::NonRepresentableArithmetic {
                    stage: ContourArithmeticStage::Qef,
                });
            }
        }
    }
    for (r, row) in a.iter_mut().enumerate() {
        row[r] += lambda;
        if !row[r].is_finite() {
            return Err(ContourError::NonRepresentableArithmetic {
                stage: ContourArithmeticStage::Qef,
            });
        }
    }
    // Cramer's rule (3×3; positive λ usually keeps the determinant away from
    // zero, while λ = 0 falls back to the mass point when singular).
    let det = det3(&a);
    if !det.is_finite() {
        return Err(ContourError::NonRepresentableArithmetic {
            stage: ContourArithmeticStage::Qef,
        });
    }
    if det.abs() < 1e-30 {
        return Ok(m);
    }
    let solve_col = |col: usize| -> Option<f64> {
        let mut ac = a;
        for r in 0..3 {
            ac[r][col] = b[r];
        }
        let value = det3(&ac) / det;
        value.is_finite().then_some(value)
    };
    let local = Point3::new(
        solve_col(0).ok_or(ContourError::NonRepresentableArithmetic {
            stage: ContourArithmeticStage::Qef,
        })?,
        solve_col(1).ok_or(ContourError::NonRepresentableArithmetic {
            stage: ContourArithmeticStage::Qef,
        })?,
        solve_col(2).ok_or(ContourError::NonRepresentableArithmetic {
            stage: ContourArithmeticStage::Qef,
        })?,
    );
    let x = Point3::new(m.x + local.x, m.y + local.y, m.z + local.z);
    if !finite_point(x) {
        return Err(ContourError::NonRepresentableArithmetic {
            stage: ContourArithmeticStage::Qef,
        });
    }
    Ok(Point3::new(
        x.x.clamp(cell_min.x, cell_max.x),
        x.y.clamp(cell_min.y, cell_max.y),
        x.z.clamp(cell_min.z, cell_max.z),
    ))
}

fn det3(a: &[[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// A triangle that failed the bracket certificate, with its margin.
#[derive(Debug, Clone, PartialEq)]
pub struct BracketFailure {
    /// Output triangle index.
    pub triangle: usize,
    /// The best (smallest) proven upper bound on |φ| over the triangle.
    pub proven_bound: f64,
    /// The tolerance it had to close under.
    pub tolerance: f64,
}

/// The bracket certificate report.
#[derive(Debug, Clone, PartialEq)]
pub struct BracketReport {
    /// Triangles proven within tolerance.
    pub proven: u64,
    /// Worst proven bound across all passing triangles.
    pub worst_margin: f64,
    /// Evaluations spent.
    pub evals: u64,
}

/// Why rigorous trace evidence at a triangle centroid was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BracketEvidenceIssue {
    /// The chart returned an estimate or an explicit no-claim.
    NonRigorous {
        /// Actual evidence strength.
        kind: NumericalKind,
    },
    /// The nominal chart value was NaN or infinite.
    NonFiniteNominal,
    /// One or both enclosure endpoints were NaN or infinite.
    NonFiniteBounds,
    /// Public certificate fields contained `lo > hi`.
    InvertedBounds,
    /// The purported enclosure did not contain the nominal evaluation.
    NominalOutsideBounds,
    /// `Exact` did not name one bit-identical finite nominal value.
    MalformedExact,
}

/// Stable stage at which finite triangle arithmetic became unrepresentable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BracketGeometryStage {
    /// A coordinate span between input vertices overflowed.
    TriangleSpan,
    /// The stable three-point centroid was not representable.
    Centroid,
    /// A conservative upper bound on centroid-to-vertex distance overflowed.
    Radius,
    /// Adding the trace enclosure magnitude to the radius overflowed.
    Bound,
    /// A recursive subdivision point was not representable.
    Subdivision,
}

/// Why a bracket certificate could not produce a pass/fail verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum BracketCertificateError {
    /// The requested tolerance was not finite and non-negative.
    InvalidTolerance {
        /// Offending tolerance.
        value: f64,
    },
    /// v1 has no global triangle-distance proof for this trace theorem.
    UnsupportedTraceClaim {
        /// The chart's advertised theorem.
        actual: TraceStepClaim,
    },
    /// A soup triangle named a vertex outside `positions`.
    InvalidTriangleIndex {
        /// Output triangle index.
        triangle: usize,
        /// Corner within the triangle (`0..3`).
        corner: usize,
        /// Offending vertex index.
        vertex: u32,
        /// Available position count.
        vertex_count: usize,
    },
    /// A referenced triangle vertex contained NaN or infinity.
    NonFiniteVertex {
        /// Output triangle index.
        triangle: usize,
        /// Corner within the triangle (`0..3`).
        corner: usize,
        /// Offending point.
        point: Point3,
    },
    /// Finite inputs could not support a finite conservative calculation.
    NonRepresentableGeometry {
        /// Output triangle index.
        triangle: usize,
        /// Arithmetic stage that refused.
        stage: BracketGeometryStage,
    },
    /// A centroid evaluation lacked usable rigorous trace evidence.
    InvalidTraceEvidence {
        /// Output triangle index.
        triangle: usize,
        /// Completed chart evaluations, including this centroid.
        completed_evaluations: u64,
        /// Nominal chart evaluation.
        nominal: f64,
        /// Evidence returned by `trace_value_enclosure`.
        evidence: NumericalCertificate,
        /// Deterministic refusal reason.
        issue: BracketEvidenceIssue,
    },
    /// Cancellation was observed before publishing a verdict.
    Cancelled {
        /// Triangles whose pass/fail classification was completed.
        completed_triangles: usize,
        /// Chart evaluations completed before cancellation was observed.
        completed_evaluations: u64,
    },
}

impl core::fmt::Display for BracketCertificateError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidTolerance { value } => write!(
                f,
                "bracket certificate refused: tolerance must be finite and non-negative, got {value}"
            ),
            Self::UnsupportedTraceClaim { actual } => write!(
                f,
                "bracket certificate refused: v1 requires the global ExactDistance theorem, got {actual:?}; a local Lipschitz sample is insufficient"
            ),
            Self::InvalidTriangleIndex {
                triangle,
                corner,
                vertex,
                vertex_count,
            } => write!(
                f,
                "bracket certificate refused: triangle {triangle} corner {corner} references vertex {vertex}, but only {vertex_count} positions exist"
            ),
            Self::NonFiniteVertex {
                triangle,
                corner,
                point,
            } => write!(
                f,
                "bracket certificate refused: triangle {triangle} corner {corner} is non-finite ({}, {}, {})",
                point.x, point.y, point.z
            ),
            Self::NonRepresentableGeometry { triangle, stage } => write!(
                f,
                "bracket certificate refused: triangle {triangle} is not representable during {stage:?}"
            ),
            Self::InvalidTraceEvidence {
                triangle,
                completed_evaluations,
                nominal,
                evidence,
                issue,
            } => write!(
                f,
                "bracket certificate refused: triangle {triangle} evaluation {completed_evaluations} returned invalid trace evidence {evidence:?} for nominal {nominal}: {issue:?}"
            ),
            Self::Cancelled {
                completed_triangles,
                completed_evaluations,
            } => write!(
                f,
                "bracket certificate cancelled after {completed_triangles} completed triangles and {completed_evaluations} completed chart evaluations; no verdict was published"
            ),
        }
    }
}

impl core::error::Error for BracketCertificateError {}

/// Compatibility name for the former unit refusal. New callers should match
/// [`BracketCertificateError`] so cancellation, evidence, and input refusals
/// remain distinguishable.
pub type NoLipschitz = BracketCertificateError;

/// Verify that every output triangle lies within `tol` of the zero set.
///
/// v1 accepts only [`TraceStepClaim::ExactDistance`], which supplies the
/// global 1-Lipschitz theorem needed to extend a rigorous centroid enclosure
/// over an entire triangle. Each recursive evaluation validates
/// [`Chart::trace_value_enclosure`]; estimates, no-claims, malformed evidence,
/// and non-representable geometry are refusals, never pass/fail verdicts.
///
/// # Errors
/// The outer [`BracketCertificateError`] reports authority, input, numerical,
/// or cancellation refusals. `Ok(Err(failures))` is a completed rigorous
/// verdict whose localized triangle bounds exceeded `tol`.
#[allow(clippy::type_complexity)] // refusal outside a localized proof verdict
pub fn bracket_certificate(
    chart: &dyn Chart,
    soup: &Soup,
    tol: f64,
    cx: &Cx<'_>,
) -> Result<Result<BracketReport, Vec<BracketFailure>>, BracketCertificateError> {
    if !tol.is_finite() || tol < 0.0 {
        return Err(BracketCertificateError::InvalidTolerance { value: tol });
    }
    let mut evals = 0u64;
    checkpoint_bracket(cx, 0, evals)?;
    let claim = chart.trace_step_claim();
    if claim != TraceStepClaim::ExactDistance {
        return Err(BracketCertificateError::UnsupportedTraceClaim { actual: claim });
    }

    // Preflight every referenced triangle before chart evaluation. This makes
    // malformed soup a deterministic refusal instead of an indexing panic or
    // a partial evidence-dependent verdict.
    for (triangle, indices) in soup.triangles.iter().copied().enumerate() {
        if triangle.is_multiple_of(256) {
            checkpoint_bracket(cx, 0, evals)?;
        }
        let _ = validated_triangle(soup, triangle, indices)?;
    }

    let mut failures = Vec::new();
    let mut worst = 0.0f64;
    for (triangle, indices) in soup.triangles.iter().copied().enumerate() {
        let [a, b, c] = validated_triangle(soup, triangle, indices)?;
        let mut vctx = VerifyCtx {
            chart,
            tol,
            triangle,
            completed_triangles: triangle,
            evals: &mut evals,
            cx,
        };
        let bound = verify_triangle(&mut vctx, a, b, c, 5)?;
        if bound <= tol {
            worst = worst.max(bound);
        } else {
            failures.push(BracketFailure {
                triangle,
                proven_bound: bound,
                tolerance: tol,
            });
        }
    }
    checkpoint_bracket(cx, soup.triangles.len(), evals)?;
    if failures.is_empty() {
        Ok(Ok(BracketReport {
            proven: u64::try_from(soup.triangles.len()).unwrap_or(u64::MAX),
            worst_margin: worst,
            evals,
        }))
    } else {
        Ok(Err(failures))
    }
}

fn checkpoint_bracket(
    cx: &Cx<'_>,
    completed_triangles: usize,
    completed_evaluations: u64,
) -> Result<(), BracketCertificateError> {
    cx.checkpoint()
        .map_err(|_| BracketCertificateError::Cancelled {
            completed_triangles,
            completed_evaluations,
        })
}

fn validated_triangle(
    soup: &Soup,
    triangle: usize,
    indices: [u32; 3],
) -> Result<[Point3; 3], BracketCertificateError> {
    let mut points = [Point3::new(0.0, 0.0, 0.0); 3];
    for (corner, vertex) in indices.into_iter().enumerate() {
        let index =
            usize::try_from(vertex).map_err(|_| BracketCertificateError::InvalidTriangleIndex {
                triangle,
                corner,
                vertex,
                vertex_count: soup.positions.len(),
            })?;
        let point = soup.positions.get(index).copied().ok_or(
            BracketCertificateError::InvalidTriangleIndex {
                triangle,
                corner,
                vertex,
                vertex_count: soup.positions.len(),
            },
        )?;
        if !(point.x.is_finite() && point.y.is_finite() && point.z.is_finite()) {
            return Err(BracketCertificateError::NonFiniteVertex {
                triangle,
                corner,
                point,
            });
        }
        points[corner] = point;
    }
    for coordinates in [
        [points[0].x, points[1].x, points[2].x],
        [points[0].y, points[1].y, points[2].y],
        [points[0].z, points[1].z, points[2].z],
    ] {
        let min = coordinates.into_iter().fold(f64::INFINITY, f64::min);
        let max = coordinates.into_iter().fold(f64::NEG_INFINITY, f64::max);
        if !(max - min).is_finite() {
            return Err(BracketCertificateError::NonRepresentableGeometry {
                triangle,
                stage: BracketGeometryStage::TriangleSpan,
            });
        }
    }
    Ok(points)
}

/// Best rigorous upper bound on distance to the zero set over one triangle,
/// refined by subdivision while the parent bound does not close.
struct VerifyCtx<'a, 'c> {
    chart: &'a dyn Chart,
    tol: f64,
    triangle: usize,
    completed_triangles: usize,
    evals: &'a mut u64,
    cx: &'a Cx<'c>,
}

impl VerifyCtx<'_, '_> {
    fn checkpoint(&self) -> Result<(), BracketCertificateError> {
        checkpoint_bracket(self.cx, self.completed_triangles, *self.evals)
    }

    fn centroid_trace_bound(&mut self, centroid: Point3) -> Result<f64, BracketCertificateError> {
        self.checkpoint()?;
        let sample = self.chart.eval(centroid, self.cx);
        *self.evals = self.evals.saturating_add(1);
        // The consumer, not the chart, owns cancellation observability. A
        // chart may request cancellation without polling it itself.
        self.checkpoint()?;
        let evidence = self.chart.trace_value_enclosure(centroid, &sample, self.cx);
        self.checkpoint()?;
        validate_trace_evidence(sample.signed_distance, evidence).map_err(|issue| {
            BracketCertificateError::InvalidTraceEvidence {
                triangle: self.triangle,
                completed_evaluations: *self.evals,
                nominal: sample.signed_distance,
                evidence,
                issue,
            }
        })
    }
}

fn validate_trace_evidence(
    nominal: f64,
    evidence: NumericalCertificate,
) -> Result<f64, BracketEvidenceIssue> {
    if !matches!(
        evidence.kind,
        NumericalKind::Exact | NumericalKind::Enclosure
    ) {
        return Err(BracketEvidenceIssue::NonRigorous {
            kind: evidence.kind,
        });
    }
    if !nominal.is_finite() {
        return Err(BracketEvidenceIssue::NonFiniteNominal);
    }
    if !(evidence.lo.is_finite() && evidence.hi.is_finite()) {
        return Err(BracketEvidenceIssue::NonFiniteBounds);
    }
    if evidence.lo > evidence.hi {
        return Err(BracketEvidenceIssue::InvertedBounds);
    }
    if nominal < evidence.lo || nominal > evidence.hi {
        return Err(BracketEvidenceIssue::NominalOutsideBounds);
    }
    if evidence.kind == NumericalKind::Exact
        && (evidence.lo.to_bits() != evidence.hi.to_bits()
            || evidence.lo.to_bits() != nominal.to_bits())
    {
        return Err(BracketEvidenceIssue::MalformedExact);
    }
    Ok(evidence.lo.abs().max(evidence.hi.abs()))
}

fn verify_triangle(
    v: &mut VerifyCtx<'_, '_>,
    a: Point3,
    b: Point3,
    c: Point3,
    depth: u32,
) -> Result<f64, BracketCertificateError> {
    let centroid =
        stable_centroid(a, b, c).ok_or(BracketCertificateError::NonRepresentableGeometry {
            triangle: v.triangle,
            stage: BracketGeometryStage::Centroid,
        })?;
    // An outward-rounded L1 radius is conservative for Euclidean distance and
    // remains stable where naive squaring/norm arithmetic would overflow.
    let mut radius = 0.0f64;
    for point in [a, b, c] {
        radius = radius.max(upper_l1_distance(point, centroid).ok_or(
            BracketCertificateError::NonRepresentableGeometry {
                triangle: v.triangle,
                stage: BracketGeometryStage::Radius,
            },
        )?);
    }
    let centroid_bound = v.centroid_trace_bound(centroid)?;
    let bound = upper_sum(centroid_bound, radius).ok_or(
        BracketCertificateError::NonRepresentableGeometry {
            triangle: v.triangle,
            stage: BracketGeometryStage::Bound,
        },
    )?;
    if bound <= v.tol || depth == 0 {
        return Ok(bound);
    }

    // 4-way midpoint subdivision. Each recursive evaluation remains fallible,
    // so cancellation and malformed evidence never collapse into a failure.
    let mab = stable_midpoint(a, b).ok_or(BracketCertificateError::NonRepresentableGeometry {
        triangle: v.triangle,
        stage: BracketGeometryStage::Subdivision,
    })?;
    let mbc = stable_midpoint(b, c).ok_or(BracketCertificateError::NonRepresentableGeometry {
        triangle: v.triangle,
        stage: BracketGeometryStage::Subdivision,
    })?;
    let mca = stable_midpoint(c, a).ok_or(BracketCertificateError::NonRepresentableGeometry {
        triangle: v.triangle,
        stage: BracketGeometryStage::Subdivision,
    })?;
    let sub = [(a, mab, mca), (mab, b, mbc), (mca, mbc, c), (mab, mbc, mca)];
    let mut worst = 0.0f64;
    for (x, y, z) in sub {
        worst = worst.max(verify_triangle(v, x, y, z, depth - 1)?);
    }
    Ok(worst)
}

fn stable_centroid(a: Point3, b: Point3, c: Point3) -> Option<Point3> {
    fn mean3(a: f64, b: f64, c: f64) -> Option<f64> {
        let midpoint = f64::midpoint(a, b);
        let delta = c - midpoint;
        let mean = midpoint + delta / 3.0;
        mean.is_finite().then_some(mean)
    }
    Some(Point3::new(
        mean3(a.x, b.x, c.x)?,
        mean3(a.y, b.y, c.y)?,
        mean3(a.z, b.z, c.z)?,
    ))
}

fn upper_l1_distance(a: Point3, b: Point3) -> Option<f64> {
    let dx = upper_abs_difference(a.x, b.x)?;
    let dy = upper_abs_difference(a.y, b.y)?;
    let dz = upper_abs_difference(a.z, b.z)?;
    upper_sum(upper_sum(dx, dy)?, dz)
}

fn upper_abs_difference(a: f64, b: f64) -> Option<f64> {
    let difference = a - b;
    if !difference.is_finite() {
        return None;
    }
    if is_zero(difference) {
        return Some(0.0);
    }
    let upper = difference.abs().next_up();
    upper.is_finite().then_some(upper)
}

fn upper_sum(a: f64, b: f64) -> Option<f64> {
    if is_zero(a) {
        return Some(b);
    }
    if is_zero(b) {
        return Some(a);
    }
    let upper = (a + b).next_up();
    upper.is_finite().then_some(upper)
}

fn is_zero(value: f64) -> bool {
    value.to_bits() << 1 == 0
}

fn stable_midpoint(a: Point3, b: Point3) -> Option<Point3> {
    let midpoint = Point3::new(
        f64::midpoint(a.x, b.x),
        f64::midpoint(a.y, b.y),
        f64::midpoint(a.z, b.z),
    );
    (midpoint.x.is_finite() && midpoint.y.is_finite() && midpoint.z.is_finite()).then_some(midpoint)
}

#[cfg(test)]
mod tests {
    use super::*;
    use asupersync::types::Budget;
    use fs_evidence::NumericalCertificate;
    use fs_exec::{CancelGate, ExecMode, StreamKey};
    use fs_geom::ChartSample;
    use std::sync::atomic::{AtomicU64, Ordering};

    struct CancellingPlane<'a> {
        gate: &'a CancelGate,
        evals: AtomicU64,
    }

    struct InvalidContourChart {
        value: f64,
        gradient: Option<Vec3>,
    }

    impl Chart for InvalidContourChart {
        fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
            ChartSample {
                signed_distance: if self.value.is_nan() { self.value } else { x.x },
                gradient: self.gradient,
                lipschitz: Some(1.0),
                error: NumericalCertificate::no_claim(),
            }
        }

        fn support(&self) -> Aabb {
            Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
        }

        fn name(&self) -> &'static str {
            "test/invalid-contour"
        }
    }

    impl Chart for CancellingPlane<'_> {
        fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
            if self.evals.fetch_add(1, Ordering::Relaxed) == 0 {
                self.gate.request();
            }
            ChartSample {
                signed_distance: x.x,
                gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
                lipschitz: Some(1.0),
                error: NumericalCertificate::enclosure(x.x, x.x),
            }
        }

        fn support(&self) -> Aabb {
            Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
        }

        fn trace_step_claim(&self) -> TraceStepClaim {
            TraceStepClaim::ExactDistance
        }

        fn name(&self) -> &'static str {
            "test/cancelling-plane"
        }
    }

    fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 0xDC,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    #[test]
    fn dual_contour_observes_cancellation_without_chart_cooperation() {
        let gate = CancelGate::new();
        let chart = CancellingPlane {
            gate: &gate,
            evals: AtomicU64::new(0),
        };
        with_cx(&gate, |cx| {
            let error = dual_contour(&chart, DcOptions::sharp(0.5), cx)
                .expect_err("consumer checkpoint must observe chart-requested cancellation");
            assert_eq!(error, ContourError::Cancelled);
        });
        assert!(chart.evals.load(Ordering::Relaxed) > 0);
    }

    #[test]
    fn dual_contour_refuses_nonfinite_field_and_gradient_samples() {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let invalid_field = InvalidContourChart {
                value: f64::NAN,
                gradient: None,
            };
            assert!(matches!(
                dual_contour(&invalid_field, DcOptions::sharp(0.5), cx),
                Err(ContourError::InvalidSample {
                    stage: ContourSampleStage::CornerLattice,
                    ..
                })
            ));

            let invalid_gradient = InvalidContourChart {
                value: 0.0,
                gradient: Some(Vec3::new(f64::NAN, 0.0, 0.0)),
            };
            assert!(matches!(
                dual_contour(&invalid_gradient, DcOptions::sharp(0.5), cx),
                Err(ContourError::InvalidGradient { .. })
            ));
        });
    }

    #[test]
    fn dual_contour_refuses_invalid_regularization_before_sampling() {
        let gate = CancelGate::new();
        let chart = CancellingPlane {
            gate: &gate,
            evals: AtomicU64::new(0),
        };
        with_cx(&gate, |cx| {
            let options = DcOptions {
                regularization: f64::NAN,
                ..DcOptions::sharp(0.5)
            };
            assert!(matches!(
                dual_contour(&chart, options, cx),
                Err(ContourError::InvalidRegularization { .. })
            ));
        });
        assert_eq!(chart.evals.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn normalized_axis_coordinates_hit_admitted_endpoints_without_overshoot() {
        let coordinates = axis_coordinates(-0.3, 1.0, 1.3, 5, "x").expect("finite axis");
        assert_eq!(coordinates.first().copied(), Some(-0.3));
        assert_eq!(coordinates.last().copied(), Some(1.0));
        assert!(
            coordinates
                .windows(2)
                .all(|pair| pair[0] <= pair[1] && pair[1] <= 1.0)
        );
    }

    #[test]
    fn positive_subnormal_span_always_admits_at_least_one_cell() {
        let span = f64::from_bits(1);
        let nodes = checked_axis_nodes(span, f64::MAX).expect("positive axis is admissible");
        assert_eq!(nodes, 2, "one positive span requires two endpoint nodes");
        assert_eq!(
            axis_coordinates(0.0, span, span, nodes, "x").expect("endpoint-only axis"),
            vec![0.0, span]
        );
    }

    #[test]
    fn near_integer_quotient_preserves_requested_maximum_cell_width() {
        let span = 0.900_000_000_000_000_1;
        let h = 0.1;
        let nodes = checked_axis_nodes(span, h).expect("near-integer quotient is admissible");
        assert_eq!(nodes, 11);
        let coordinates = axis_coordinates(0.0, span, span, nodes, "x").expect("finite axis");
        assert!(
            coordinates.windows(2).all(|pair| pair[1] - pair[0] <= h),
            "realized contour cells must not exceed requested h"
        );
    }

    #[test]
    fn bracket_consumer_observes_chart_requested_cancellation_with_progress() {
        let gate = CancelGate::new();
        let chart = CancellingPlane {
            gate: &gate,
            evals: AtomicU64::new(0),
        };
        let soup = Soup {
            positions: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 0.25, 0.0),
                Point3::new(0.0, 0.0, 0.25),
            ],
            triangles: vec![[0, 1, 2]],
        };
        with_cx(&gate, |cx| {
            let error = bracket_certificate(&chart, &soup, 0.5, cx)
                .expect_err("consumer checkpoint must observe chart-requested cancellation");
            assert_eq!(
                error,
                BracketCertificateError::Cancelled {
                    completed_triangles: 0,
                    completed_evaluations: 1,
                }
            );
        });
        assert_eq!(chart.evals.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn bracket_preflight_refuses_invalid_tolerance_and_geometry_without_evaluation() {
        let gate = CancelGate::new();
        let chart = CancellingPlane {
            gate: &gate,
            evals: AtomicU64::new(0),
        };
        let triangle = Soup {
            positions: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 0.25, 0.0),
                Point3::new(0.0, 0.0, 0.25),
            ],
            triangles: vec![[0, 1, 2]],
        };
        with_cx(&gate, |cx| {
            assert!(matches!(
                bracket_certificate(&chart, &triangle, f64::NAN, cx),
                Err(BracketCertificateError::InvalidTolerance { value }) if value.is_nan()
            ));
            assert_eq!(
                bracket_certificate(&chart, &triangle, -1.0, cx),
                Err(BracketCertificateError::InvalidTolerance { value: -1.0 })
            );

            let invalid_index = Soup {
                positions: triangle.positions.clone(),
                triangles: vec![[0, 1, 99]],
            };
            assert!(matches!(
                bracket_certificate(&chart, &invalid_index, 0.5, cx),
                Err(BracketCertificateError::InvalidTriangleIndex {
                    triangle: 0,
                    corner: 2,
                    vertex: 99,
                    ..
                })
            ));

            let nonrepresentable = Soup {
                positions: vec![
                    Point3::new(-f64::MAX, 0.0, 0.0),
                    Point3::new(f64::MAX, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                ],
                triangles: vec![[0, 1, 2]],
            };
            assert_eq!(
                bracket_certificate(&chart, &nonrepresentable, 0.5, cx),
                Err(BracketCertificateError::NonRepresentableGeometry {
                    triangle: 0,
                    stage: BracketGeometryStage::TriangleSpan,
                })
            );
        });
        assert_eq!(
            chart.evals.load(Ordering::Relaxed),
            0,
            "all input refusals precede chart evaluation"
        );
    }
}
