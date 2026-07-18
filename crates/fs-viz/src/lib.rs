//! fs-viz — scientific visualization primitives. Layer: L5.
//!
//! A 10⁸-cell field is illegible as raw numbers; the job of visualization is to
//! compress it into a picture an agent can read in one glance. This v0 is the
//! verifiable core of that job — the pieces with ANALYTIC ground truth, so a
//! rendered topology is a certificate rather than a guess:
//!
//! - [`streamline`] — RK4 integration of a flow field; on a rotation field the
//!   streamline is a circle (radius conserved), on a saddle field it is a
//!   hyperbola (`xy` conserved);
//! - [`classify_hessian`] — the Morse critical-point type (min / max / saddle)
//!   and index from a scalar field's Hessian — the atom of a Morse–Smale
//!   complex;
//! - [`Grid2::isocontour_crossings`] — bounded, fail-closed isocontour edge
//!   intersections of an admitted scalar grid, which on a circle SDF all lie
//!   on the circle.
//! - [`Grid3::isosurface`] — bounded deterministic marching tetrahedra from an
//!   owned x-fastest scalar grid to a renderer-ready indexed triangle mesh.
//! - [`ScalarField3`] — a bounded versioned artifact codec with explicit
//!   node/cell centering, quantity, and units for ledger composition.
//!
//! Deterministic; no dependencies.

mod isosurface;
mod scalar_field;

pub use isosurface::{Grid3, Grid3Error, IsoMesh3, IsoSurfaceError, Vec3};
pub use scalar_field::{
    SCALAR_FIELD3_ARTIFACT_KIND, SCALAR_FIELD3_SCHEMA_VERSION, ScalarField3, ScalarField3Error,
    ScalarFieldSemantics, ScalarLayout3,
};

/// A 2-D point / vector.
pub type Vec2 = [f64; 2];

/// Integrate a streamline of `field` from `seed` with `steps` RK4 steps of size
/// `dt`. Returns the ordered polyline (including the seed).
#[must_use]
pub fn streamline(field: impl Fn(Vec2) -> Vec2, seed: Vec2, dt: f64, steps: usize) -> Vec<Vec2> {
    let mut pts = Vec::with_capacity(steps + 1);
    pts.push(seed);
    let mut p = seed;
    for _ in 0..steps {
        let k1 = field(p);
        let k2 = field([p[0] + 0.5 * dt * k1[0], p[1] + 0.5 * dt * k1[1]]);
        let k3 = field([p[0] + 0.5 * dt * k2[0], p[1] + 0.5 * dt * k2[1]]);
        let k4 = field([p[0] + dt * k3[0], p[1] + dt * k3[1]]);
        p = [
            p[0] + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
            p[1] + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        ];
        pts.push(p);
    }
    pts
}

/// The Morse type of a critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriticalKind {
    /// A local minimum (Morse index 0).
    Minimum,
    /// A saddle (Morse index 1 in 2-D).
    Saddle,
    /// A local maximum (Morse index 2 in 2-D).
    Maximum,
    /// Degenerate or unclassifiable (a zero eigenvalue or invalid numeric
    /// input) — not a Morse critical point.
    Degenerate,
}

/// A classified critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CriticalPoint {
    /// The Morse type.
    pub kind: CriticalKind,
    /// The Morse index (number of negative Hessian eigenvalues).
    pub morse_index: usize,
}

/// Classify a critical point from its (symmetric) `2×2` Hessian: the Morse index
/// is the number of negative eigenvalues; a near-zero eigenvalue is degenerate.
/// Non-finite entries and non-finite or negative tolerances fail closed as
/// degenerate with no claimed negative eigenvalues.
#[must_use]
pub fn classify_hessian(hessian: [[f64; 2]; 2], tol: f64) -> CriticalPoint {
    if !tol.is_finite() || tol < 0.0 || hessian.iter().flatten().any(|value| !value.is_finite()) {
        return CriticalPoint {
            kind: CriticalKind::Degenerate,
            morse_index: 0,
        };
    }
    let (a, b, c) = (hessian[0][0], hessian[0][1], hessian[1][1]);
    let mean = f64::midpoint(a, c);
    let half_diff = f64::midpoint(a, -c);
    let scale = half_diff.abs().max(b.abs());
    let disc = if scale == 0.0 {
        0.0
    } else {
        let x = half_diff / scale;
        let y = b / scale;
        scale * (x * x + y * y).sqrt()
    };
    let (l1, l2) = (mean - disc, mean + disc);
    if l1.abs() <= tol || l2.abs() <= tol {
        return CriticalPoint {
            kind: CriticalKind::Degenerate,
            morse_index: usize::from(l1 < -tol) + usize::from(l2 < -tol),
        };
    }
    let morse_index = usize::from(l1 < 0.0) + usize::from(l2 < 0.0);
    let kind = match morse_index {
        0 => CriticalKind::Minimum,
        1 => CriticalKind::Saddle,
        _ => CriticalKind::Maximum,
    };
    CriticalPoint { kind, morse_index }
}

/// A scalar field sampled on a regular 2-D grid.
#[derive(Debug, Clone)]
pub struct Grid2 {
    nx: usize,
    ny: usize,
    lo: Vec2,
    hi: Vec2,
    values: Vec<f64>,
}

/// Admission failures for a regular 2-D scalar grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Grid2Error {
    /// Both dimensions must contain at least two sample nodes.
    InvalidDimensions {
        /// Rejected `[nx, ny]` dimensions.
        dimensions: [usize; 2],
    },
    /// The dimension product does not fit in `usize`.
    NodeCountOverflow {
        /// Rejected `[nx, ny]` dimensions.
        dimensions: [usize; 2],
    },
    /// The explicit sampling budget is smaller than the grid.
    NodeBudgetExceeded {
        /// Required sample count.
        required: usize,
        /// Caller-provided limit.
        limit: usize,
    },
    /// A world-space interval is non-finite, non-increasing, or has a
    /// non-finite extent.
    InvalidBounds {
        /// Cartesian axis, `0..2`.
        axis: usize,
        /// Rejected lower bound.
        lower: f64,
        /// Rejected upper bound.
        upper: f64,
    },
    /// Distinct logical nodes collapse onto the same or a reversed floating
    /// coordinate at the requested resolution.
    UnrepresentableCoordinates {
        /// Cartesian axis, `0..2`.
        axis: usize,
        /// Earlier logical node index.
        first_index: usize,
        /// Earlier generated coordinate.
        first: f64,
        /// Later logical node index.
        second_index: usize,
        /// Later generated coordinate.
        second: f64,
    },
    /// A sampled scalar is non-finite.
    NonFiniteValue {
        /// Linear x-fastest node index.
        index: usize,
        /// Rejected value.
        value: f64,
    },
    /// The requested sample allocation could not be reserved.
    AllocationFailed {
        /// Requested node count.
        nodes: usize,
    },
}

impl core::fmt::Display for Grid2Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDimensions { dimensions } => write!(
                f,
                "Grid2 dimensions must each be at least two (got {dimensions:?})"
            ),
            Self::NodeCountOverflow { dimensions } => write!(
                f,
                "Grid2 node count overflows for dimensions {dimensions:?}"
            ),
            Self::NodeBudgetExceeded { required, limit } => write!(
                f,
                "Grid2 requires {required} nodes, exceeding the explicit limit {limit}"
            ),
            Self::InvalidBounds { axis, lower, upper } => write!(
                f,
                "Grid2 axis {axis} bounds must be finite and increasing with finite extent (got {lower}..{upper})"
            ),
            Self::UnrepresentableCoordinates {
                axis,
                first_index,
                first,
                second_index,
                second,
            } => write!(
                f,
                "Grid2 axis {axis} nodes {first_index} and {second_index} are not strictly ordered ({first}, {second})"
            ),
            Self::NonFiniteValue { index, value } => {
                write!(f, "Grid2 value {index} is non-finite ({value})")
            }
            Self::AllocationFailed { nodes } => {
                write!(f, "Grid2 could not reserve storage for {nodes} nodes")
            }
        }
    }
}

impl std::error::Error for Grid2Error {}

/// Failures from bounded 2-D isocontour edge extraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IsoContourError {
    /// The requested isovalue is non-finite.
    NonFiniteIso {
        /// Rejected isovalue.
        iso: f64,
    },
    /// At least one output crossing must be admitted by the explicit budget.
    ZeroCrossingLimit,
    /// Extraction found more crossings than the caller admitted.
    CrossingBudgetExceeded {
        /// Caller-provided maximum crossing count.
        limit: usize,
    },
    /// A whole grid edge lies on the requested level, so its intersection is a
    /// segment that cannot be represented by a point-only result.
    CoincidentLevelEdge {
        /// First endpoint as `[i, j]`.
        first: [usize; 2],
        /// Second endpoint as `[i, j]`.
        second: [usize; 2],
    },
    /// Storage for the next crossing could not be reserved.
    AllocationFailed {
        /// Number of crossings, including the one that could not be reserved.
        required: usize,
    },
    /// Interpolation produced a non-finite point or invalid interpolation
    /// fraction from otherwise finite inputs.
    NonFiniteGeometry,
}

impl core::fmt::Display for IsoContourError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFiniteIso { iso } => {
                write!(f, "isocontour level must be finite (got {iso})")
            }
            Self::ZeroCrossingLimit => {
                write!(f, "isocontour crossing limit must be positive")
            }
            Self::CrossingBudgetExceeded { limit } => {
                write!(f, "isocontour exceeded the explicit {limit}-crossing limit")
            }
            Self::CoincidentLevelEdge { first, second } => write!(
                f,
                "isocontour edge {first:?}..{second:?} lies wholly on the requested level"
            ),
            Self::AllocationFailed { required } => write!(
                f,
                "isocontour could not reserve storage for {required} crossings"
            ),
            Self::NonFiniteGeometry => {
                write!(f, "isocontour interpolation produced non-finite geometry")
            }
        }
    }
}

impl std::error::Error for IsoContourError {}

impl Grid2 {
    /// Sample `f` on an admitted `nx × ny` regular grid spanning `[lo, hi]`.
    ///
    /// No sample is evaluated until dimensions, the explicit node budget, and
    /// world bounds have been admitted, sample storage has been reserved, and
    /// every generated axis coordinate is strictly increasing. Sampling order
    /// is row-major with x fastest.
    ///
    /// # Errors
    /// [`Grid2Error`] for invalid dimensions/bounds, count overflow,
    /// coordinate collapse, budget or allocation refusal, or the first
    /// non-finite sample.
    pub fn from_fn(
        nx: usize,
        ny: usize,
        lo: Vec2,
        hi: Vec2,
        node_limit: usize,
        mut f: impl FnMut(Vec2) -> f64,
    ) -> Result<Grid2, Grid2Error> {
        let node_count = validate_grid2_layout(nx, ny, lo, hi, node_limit)?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(node_count)
            .map_err(|_| Grid2Error::AllocationFailed { nodes: node_count })?;
        validate_grid2_coordinates(nx, ny, lo, hi)?;
        for j in 0..ny {
            for i in 0..nx {
                let p = grid_point(nx, ny, lo, hi, i, j);
                let value = f(p);
                if !value.is_finite() {
                    return Err(Grid2Error::NonFiniteValue {
                        index: values.len(),
                        value,
                    });
                }
                values.push(value);
            }
        }
        Ok(Grid2 {
            nx,
            ny,
            lo,
            hi,
            values,
        })
    }

    /// The scalar value at grid node `(i, j)`.
    ///
    /// # Panics
    /// Panics when either index is outside the admitted grid.
    #[must_use]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        assert!(i < self.nx && j < self.ny, "Grid2 node index out of bounds");
        self.values[j * self.nx + i]
    }

    /// The world coordinate of grid node `(i, j)`.
    ///
    /// # Panics
    /// Callers must keep both indices inside the admitted grid. Out-of-range
    /// coordinates are not part of the public grid domain.
    #[must_use]
    pub fn point(&self, i: usize, j: usize) -> Vec2 {
        assert!(i < self.nx && j < self.ny, "Grid2 node index out of bounds");
        grid_point(self.nx, self.ny, self.lo, self.hi, i, j)
    }

    /// The bounded isocontour edge crossings at `iso`.
    ///
    /// Strict sign changes yield one scaled-interpolation point. An exact-level
    /// endpoint is emitted once at its original bits even when several incident
    /// edges meet there. A wholly coincident edge is refused because its
    /// intersection is a segment, not a unique point. Traversal and first-node
    /// ownership are row-major, with the positive-x edge before positive-y.
    ///
    /// # Errors
    /// [`IsoContourError`] for a non-finite level, zero or exceeded output
    /// budget, a coincident level edge, allocation refusal, or non-finite
    /// interpolated geometry. Every failure returns no partial crossing vector.
    pub fn isocontour_crossings(
        &self,
        iso: f64,
        crossing_limit: usize,
    ) -> Result<Vec<Vec2>, IsoContourError> {
        if !iso.is_finite() {
            return Err(IsoContourError::NonFiniteIso { iso });
        }
        if crossing_limit == 0 {
            return Err(IsoContourError::ZeroCrossingLimit);
        }
        let mut out = Vec::new();
        let mut exact_nodes = Vec::new();
        for j in 0..self.ny {
            for i in 0..self.nx {
                let (p, v) = (self.point(i, j), self.at(i, j));
                let node = j * self.nx + i;
                if i + 1 < self.nx {
                    if let Some(crossing) = edge_crossing(
                        [i, j],
                        node,
                        p,
                        v,
                        [i + 1, j],
                        node + 1,
                        self.point(i + 1, j),
                        self.at(i + 1, j),
                        iso,
                    )? {
                        push_crossing(&mut out, &mut exact_nodes, crossing, crossing_limit)?;
                    }
                }
                if j + 1 < self.ny {
                    if let Some(crossing) = edge_crossing(
                        [i, j],
                        node,
                        p,
                        v,
                        [i, j + 1],
                        node + self.nx,
                        self.point(i, j + 1),
                        self.at(i, j + 1),
                        iso,
                    )? {
                        push_crossing(&mut out, &mut exact_nodes, crossing, crossing_limit)?;
                    }
                }
            }
        }
        Ok(out)
    }
}

fn validate_grid2_layout(
    nx: usize,
    ny: usize,
    lo: Vec2,
    hi: Vec2,
    node_limit: usize,
) -> Result<usize, Grid2Error> {
    let dimensions = [nx, ny];
    if dimensions.into_iter().any(|dimension| dimension < 2) {
        return Err(Grid2Error::InvalidDimensions { dimensions });
    }
    let node_count = nx
        .checked_mul(ny)
        .ok_or(Grid2Error::NodeCountOverflow { dimensions })?;
    if node_count > node_limit {
        return Err(Grid2Error::NodeBudgetExceeded {
            required: node_count,
            limit: node_limit,
        });
    }
    for axis in 0..2 {
        let extent = hi[axis] - lo[axis];
        if !(lo[axis].is_finite()
            && hi[axis].is_finite()
            && hi[axis] > lo[axis]
            && extent.is_finite())
        {
            return Err(Grid2Error::InvalidBounds {
                axis,
                lower: lo[axis],
                upper: hi[axis],
            });
        }
    }
    Ok(node_count)
}

fn validate_grid2_coordinates(nx: usize, ny: usize, lo: Vec2, hi: Vec2) -> Result<(), Grid2Error> {
    for (axis, nodes) in [nx, ny].into_iter().enumerate() {
        let mut previous = grid_axis_point(nodes, lo[axis], hi[axis], 0);
        for index in 1..nodes {
            let current = grid_axis_point(nodes, lo[axis], hi[axis], index);
            if !current.is_finite() || current <= previous {
                return Err(Grid2Error::UnrepresentableCoordinates {
                    axis,
                    first_index: index - 1,
                    first: previous,
                    second_index: index,
                    second: current,
                });
            }
            previous = current;
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
enum EdgeCrossing2 {
    Exact { node: usize, point: Vec2 },
    Interpolated(Vec2),
}

fn edge_crossing(
    a_index: [usize; 2],
    a_node: usize,
    a: Vec2,
    va: f64,
    b_index: [usize; 2],
    b_node: usize,
    b: Vec2,
    vb: f64,
    iso: f64,
) -> Result<Option<EdgeCrossing2>, IsoContourError> {
    let a_exact = va == iso;
    let b_exact = vb == iso;
    if a_exact && b_exact {
        return Err(IsoContourError::CoincidentLevelEdge {
            first: a_index,
            second: b_index,
        });
    }
    if a_exact {
        return Ok(Some(EdgeCrossing2::Exact {
            node: a_node,
            point: a,
        }));
    }
    if b_exact {
        return Ok(Some(EdgeCrossing2::Exact {
            node: b_node,
            point: b,
        }));
    }
    if (va < iso) == (vb < iso) {
        return Ok(None);
    }
    let scale = va.abs().max(vb.abs()).max(iso.abs());
    let scaled_iso = iso / scale;
    let a_distance = (va / scale - scaled_iso).abs();
    let b_distance = (vb / scale - scaled_iso).abs();
    let t = a_distance / (a_distance + b_distance);
    if !t.is_finite() || !(0.0 < t && t < 1.0) {
        return Err(IsoContourError::NonFiniteGeometry);
    }
    let point = [
        (b[0] - a[0]).mul_add(t, a[0]),
        (b[1] - a[1]).mul_add(t, a[1]),
    ];
    if point.into_iter().any(|coordinate| !coordinate.is_finite()) {
        return Err(IsoContourError::NonFiniteGeometry);
    }
    Ok(Some(EdgeCrossing2::Interpolated(point)))
}

fn push_crossing(
    crossings: &mut Vec<Vec2>,
    exact_nodes: &mut Vec<usize>,
    crossing: EdgeCrossing2,
    crossing_limit: usize,
) -> Result<(), IsoContourError> {
    let (point, exact_node) = match crossing {
        EdgeCrossing2::Exact { node, point } => {
            if exact_nodes.contains(&node) {
                return Ok(());
            }
            (point, Some(node))
        }
        EdgeCrossing2::Interpolated(point) => (point, None),
    };
    if crossings.len() == crossing_limit {
        return Err(IsoContourError::CrossingBudgetExceeded {
            limit: crossing_limit,
        });
    }
    let required = crossings
        .len()
        .checked_add(1)
        .ok_or(IsoContourError::AllocationFailed {
            required: usize::MAX,
        })?;
    crossings
        .try_reserve(1)
        .map_err(|_| IsoContourError::AllocationFailed { required })?;
    if let Some(node) = exact_node {
        exact_nodes
            .try_reserve(1)
            .map_err(|_| IsoContourError::AllocationFailed { required })?;
        exact_nodes.push(node);
    }
    crossings.push(point);
    Ok(())
}

fn grid_point(nx: usize, ny: usize, lo: Vec2, hi: Vec2, i: usize, j: usize) -> Vec2 {
    [
        grid_axis_point(nx, lo[0], hi[0], i),
        grid_axis_point(ny, lo[1], hi[1], j),
    ]
}

fn grid_axis_point(nodes: usize, lower: f64, upper: f64, index: usize) -> f64 {
    if index == 0 {
        lower
    } else if index + 1 == nodes {
        upper
    } else {
        let t = index as f64 / (nodes - 1) as f64;
        (upper - lower).mul_add(t, lower)
    }
}
