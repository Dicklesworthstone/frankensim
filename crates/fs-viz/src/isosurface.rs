//! Bounded deterministic isosurface extraction on regular 3-D scalar grids.

use std::collections::BTreeMap;

/// A 3-D point or vector.
pub type Vec3 = [f64; 3];

/// Validation failures for an owned regular scalar grid.
#[derive(Debug, Clone, PartialEq)]
pub enum Grid3Error {
    /// Every dimension must contain at least two sample nodes.
    InvalidDimensions {
        /// Rejected `[nx, ny, nz]` dimensions.
        dimensions: [usize; 3],
    },
    /// The dimension product does not fit in `usize`.
    NodeCountOverflow {
        /// Rejected dimensions.
        dimensions: [usize; 3],
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
        /// Cartesian axis, `0..3`.
        axis: usize,
        /// Rejected lower bound.
        lower: f64,
        /// Rejected upper bound.
        upper: f64,
    },
    /// Supplied storage does not match the dimension product.
    ValueCountMismatch {
        /// Required value count.
        expected: usize,
        /// Supplied value count.
        actual: usize,
    },
    /// A sampled or supplied scalar is non-finite.
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

impl core::fmt::Display for Grid3Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDimensions { dimensions } => {
                write!(
                    f,
                    "Grid3 dimensions must each be at least two (got {dimensions:?})"
                )
            }
            Self::NodeCountOverflow { dimensions } => {
                write!(
                    f,
                    "Grid3 node count overflows for dimensions {dimensions:?}"
                )
            }
            Self::NodeBudgetExceeded { required, limit } => write!(
                f,
                "Grid3 requires {required} nodes, exceeding the explicit limit {limit}"
            ),
            Self::InvalidBounds { axis, lower, upper } => write!(
                f,
                "Grid3 axis {axis} bounds must be finite and increasing with finite extent (got {lower}..{upper})"
            ),
            Self::ValueCountMismatch { expected, actual } => write!(
                f,
                "Grid3 requires {expected} x-fastest values but received {actual}"
            ),
            Self::NonFiniteValue { index, value } => {
                write!(f, "Grid3 value {index} is non-finite ({value})")
            }
            Self::AllocationFailed { nodes } => {
                write!(f, "Grid3 could not reserve storage for {nodes} nodes")
            }
        }
    }
}

impl std::error::Error for Grid3Error {}

/// Failures from bounded marching-tetrahedra extraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IsoSurfaceError {
    /// The requested isovalue is non-finite.
    NonFiniteIso {
        /// Rejected isovalue.
        iso: f64,
    },
    /// At least one triangle must be admitted by the explicit budget.
    ZeroTriangleLimit,
    /// Extraction reached the caller-provided triangle budget.
    TriangleBudgetExceeded {
        /// Caller-provided maximum triangle count.
        limit: usize,
    },
    /// The indexed vertex count cannot be represented by `u32` triangles.
    VertexIndexOverflow,
    /// Interpolation or orientation produced non-finite geometry.
    NonFiniteGeometry,
}

impl core::fmt::Display for IsoSurfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFiniteIso { iso } => write!(f, "isosurface level must be finite (got {iso})"),
            Self::ZeroTriangleLimit => write!(f, "isosurface triangle limit must be positive"),
            Self::TriangleBudgetExceeded { limit } => {
                write!(f, "isosurface exceeded the explicit {limit}-triangle limit")
            }
            Self::VertexIndexOverflow => write!(f, "isosurface vertex count exceeds u32 indexing"),
            Self::NonFiniteGeometry => {
                write!(f, "isosurface interpolation produced non-finite geometry")
            }
        }
    }
}

impl std::error::Error for IsoSurfaceError {}

/// Scalar samples on a regular x-fastest 3-D grid.
#[derive(Debug, Clone, PartialEq)]
pub struct Grid3 {
    dimensions: [usize; 3],
    lower: Vec3,
    upper: Vec3,
    spacing: Vec3,
    values: Vec<f64>,
}

impl Grid3 {
    /// Validate and adopt x-fastest scalar samples.
    ///
    /// `values[(k * ny + j) * nx + i]` is the node at `(i, j, k)`. The
    /// explicit `node_limit` is checked even though the caller already owns
    /// the vector, so admission policy is identical to [`Grid3::from_fn`].
    ///
    /// # Errors
    /// [`Grid3Error`] for invalid dimensions/bounds, budget or length
    /// mismatch, overflow, or a non-finite sample.
    pub fn from_values(
        dimensions: [usize; 3],
        lower: Vec3,
        upper: Vec3,
        node_limit: usize,
        values: Vec<f64>,
    ) -> Result<Self, Grid3Error> {
        let (node_count, spacing) = validate_grid3_layout(dimensions, lower, upper, node_limit)?;
        if values.len() != node_count {
            return Err(Grid3Error::ValueCountMismatch {
                expected: node_count,
                actual: values.len(),
            });
        }
        for (index, value) in values.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(Grid3Error::NonFiniteValue { index, value });
            }
        }
        Ok(Self {
            dimensions,
            lower,
            upper,
            spacing,
            values,
        })
    }

    /// Sample a scalar function on a bounded regular grid.
    ///
    /// Sampling order is fixed at z/y/x with x fastest in storage. No sample
    /// is evaluated until dimensions, bounds, overflow, and `node_limit` have
    /// been admitted.
    ///
    /// # Errors
    /// [`Grid3Error`] for invalid geometry/budget, allocation refusal, or a
    /// non-finite function value.
    pub fn from_fn(
        dimensions: [usize; 3],
        lower: Vec3,
        upper: Vec3,
        node_limit: usize,
        mut field: impl FnMut(Vec3) -> f64,
    ) -> Result<Self, Grid3Error> {
        let (node_count, spacing) = validate_grid3_layout(dimensions, lower, upper, node_limit)?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(node_count)
            .map_err(|_| Grid3Error::AllocationFailed { nodes: node_count })?;
        for k in 0..dimensions[2] {
            for j in 0..dimensions[1] {
                for i in 0..dimensions[0] {
                    let point = grid3_point(lower, spacing, i, j, k);
                    let value = field(point);
                    if !value.is_finite() {
                        return Err(Grid3Error::NonFiniteValue {
                            index: values.len(),
                            value,
                        });
                    }
                    values.push(value);
                }
            }
        }
        Ok(Self {
            dimensions,
            lower,
            upper,
            spacing,
            values,
        })
    }

    /// `[nx, ny, nz]` node dimensions.
    #[must_use]
    pub const fn dimensions(&self) -> [usize; 3] {
        self.dimensions
    }

    /// Inclusive world-space lower and upper node coordinates.
    #[must_use]
    pub const fn bounds(&self) -> [Vec3; 2] {
        [self.lower, self.upper]
    }

    /// Read-only x-fastest scalar storage.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Scalar value at a checked node coordinate.
    #[must_use]
    pub fn at(&self, i: usize, j: usize, k: usize) -> Option<f64> {
        self.linear_index(i, j, k).map(|index| self.values[index])
    }

    /// World coordinate of a checked node.
    #[must_use]
    pub fn point(&self, i: usize, j: usize, k: usize) -> Option<Vec3> {
        self.linear_index(i, j, k)
            .map(|_| grid3_point(self.lower, self.spacing, i, j, k))
    }

    /// Extract an indexed isosurface with deterministic marching tetrahedra.
    ///
    /// Every cube is decomposed into the same six tetrahedra around its `0-7`
    /// diagonal. Vertices are shared by a canonical global edge/node key,
    /// triangles are oriented from the `< iso` side toward the `>= iso` side,
    /// and no partial mesh is returned when `triangle_limit` is exceeded.
    ///
    /// # Errors
    /// [`IsoSurfaceError`] for an invalid level/budget, index exhaustion, or
    /// non-finite interpolated geometry.
    pub fn isosurface(&self, iso: f64, triangle_limit: usize) -> Result<IsoMesh3, IsoSurfaceError> {
        if !iso.is_finite() {
            return Err(IsoSurfaceError::NonFiniteIso { iso });
        }
        if triangle_limit == 0 {
            return Err(IsoSurfaceError::ZeroTriangleLimit);
        }

        let mut extractor = Extractor3 {
            grid: self,
            iso,
            triangle_limit,
            vertices: Vec::new(),
            triangles: Vec::new(),
            vertex_by_edge: BTreeMap::new(),
        };
        let [nx, ny, nz] = self.dimensions;
        for k in 0..(nz - 1) {
            for j in 0..(ny - 1) {
                for i in 0..(nx - 1) {
                    let cube = [
                        self.index_unchecked(i, j, k),
                        self.index_unchecked(i + 1, j, k),
                        self.index_unchecked(i, j + 1, k),
                        self.index_unchecked(i + 1, j + 1, k),
                        self.index_unchecked(i, j, k + 1),
                        self.index_unchecked(i + 1, j, k + 1),
                        self.index_unchecked(i, j + 1, k + 1),
                        self.index_unchecked(i + 1, j + 1, k + 1),
                    ];
                    for tetrahedron in TETRAHEDRA3 {
                        extractor.polygonize([
                            cube[tetrahedron[0]],
                            cube[tetrahedron[1]],
                            cube[tetrahedron[2]],
                            cube[tetrahedron[3]],
                        ])?;
                    }
                }
            }
        }
        Ok(IsoMesh3 {
            vertices: extractor.vertices,
            triangles: extractor.triangles,
        })
    }

    fn linear_index(&self, i: usize, j: usize, k: usize) -> Option<usize> {
        (i < self.dimensions[0] && j < self.dimensions[1] && k < self.dimensions[2])
            .then(|| self.index_unchecked(i, j, k))
    }

    fn index_unchecked(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.dimensions[1] + j) * self.dimensions[0] + i
    }

    fn point_from_index(&self, index: usize) -> Vec3 {
        let i = index % self.dimensions[0];
        let rest = index / self.dimensions[0];
        let j = rest % self.dimensions[1];
        let k = rest / self.dimensions[1];
        grid3_point(self.lower, self.spacing, i, j, k)
    }
}

/// Indexed triangle mesh produced by [`Grid3::isosurface`].
#[derive(Debug, Clone, PartialEq)]
pub struct IsoMesh3 {
    vertices: Vec<Vec3>,
    triangles: Vec<[u32; 3]>,
}

impl IsoMesh3 {
    /// Indexed vertices in deterministic first-crossing order.
    #[must_use]
    pub fn vertices(&self) -> &[Vec3] {
        &self.vertices
    }

    /// Outward-oriented triangles in grid/tetrahedron traversal order.
    #[must_use]
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// Consume the mesh into renderer-ready vertex and triangle arrays.
    #[must_use]
    pub fn into_parts(self) -> (Vec<Vec3>, Vec<[u32; 3]>) {
        (self.vertices, self.triangles)
    }

    /// Piecewise-linear surface area.
    #[must_use]
    pub fn surface_area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|triangle| {
                let a = self.vertices[triangle[0] as usize];
                let b = self.vertices[triangle[1] as usize];
                let c = self.vertices[triangle[2] as usize];
                let ab = subtract3(b, a);
                let ac = subtract3(c, a);
                0.5 * norm_squared3(cross3(ab, ac)).sqrt()
            })
            .sum()
    }
}

fn validate_grid3_layout(
    dimensions: [usize; 3],
    lower: Vec3,
    upper: Vec3,
    node_limit: usize,
) -> Result<(usize, Vec3), Grid3Error> {
    if dimensions.into_iter().any(|dimension| dimension < 2) {
        return Err(Grid3Error::InvalidDimensions { dimensions });
    }
    let node_count = dimensions
        .into_iter()
        .try_fold(1usize, |count, dimension| count.checked_mul(dimension))
        .ok_or(Grid3Error::NodeCountOverflow { dimensions })?;
    if node_count > node_limit {
        return Err(Grid3Error::NodeBudgetExceeded {
            required: node_count,
            limit: node_limit,
        });
    }
    let mut spacing = [0.0; 3];
    for axis in 0..3 {
        let extent = upper[axis] - lower[axis];
        if !(lower[axis].is_finite()
            && upper[axis].is_finite()
            && upper[axis] > lower[axis]
            && extent.is_finite())
        {
            return Err(Grid3Error::InvalidBounds {
                axis,
                lower: lower[axis],
                upper: upper[axis],
            });
        }
        spacing[axis] = extent / (dimensions[axis] - 1) as f64;
    }
    Ok((node_count, spacing))
}

fn grid3_point(lower: Vec3, spacing: Vec3, i: usize, j: usize, k: usize) -> Vec3 {
    [
        spacing[0].mul_add(i as f64, lower[0]),
        spacing[1].mul_add(j as f64, lower[1]),
        spacing[2].mul_add(k as f64, lower[2]),
    ]
}

const TETRAHEDRA3: [[usize; 4]; 6] = [
    [0, 7, 1, 3],
    [0, 7, 3, 2],
    [0, 7, 2, 6],
    [0, 7, 6, 4],
    [0, 7, 4, 5],
    [0, 7, 5, 1],
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum VertexKey3 {
    Node(usize),
    Edge(usize, usize),
}

struct Extractor3<'a> {
    grid: &'a Grid3,
    iso: f64,
    triangle_limit: usize,
    vertices: Vec<Vec3>,
    triangles: Vec<[u32; 3]>,
    vertex_by_edge: BTreeMap<VertexKey3, u32>,
}

impl Extractor3<'_> {
    fn polygonize(&mut self, nodes: [usize; 4]) -> Result<(), IsoSurfaceError> {
        let mut inside = [0usize; 4];
        let mut inside_count = 0usize;
        let mut outside = [0usize; 4];
        let mut outside_count = 0usize;
        for (local, node) in nodes.into_iter().enumerate() {
            if self.grid.values[node] < self.iso {
                inside[inside_count] = local;
                inside_count += 1;
            } else {
                outside[outside_count] = local;
                outside_count += 1;
            }
        }
        if inside_count == 0 || inside_count == 4 {
            return Ok(());
        }
        let outward =
            self.outward_direction(nodes, &inside[..inside_count], &outside[..outside_count]);
        match inside_count {
            1 => {
                let a = inside[0];
                let triangle = [
                    self.edge_vertex(nodes[a], nodes[outside[0]])?,
                    self.edge_vertex(nodes[a], nodes[outside[1]])?,
                    self.edge_vertex(nodes[a], nodes[outside[2]])?,
                ];
                self.push_triangle(triangle, outward)
            }
            3 => {
                let a = outside[0];
                let triangle = [
                    self.edge_vertex(nodes[a], nodes[inside[0]])?,
                    self.edge_vertex(nodes[a], nodes[inside[1]])?,
                    self.edge_vertex(nodes[a], nodes[inside[2]])?,
                ];
                self.push_triangle(triangle, outward)
            }
            2 => {
                let (a, b) = (inside[0], inside[1]);
                let (c, d) = (outside[0], outside[1]);
                let ac = self.edge_vertex(nodes[a], nodes[c])?;
                let ad = self.edge_vertex(nodes[a], nodes[d])?;
                let bd = self.edge_vertex(nodes[b], nodes[d])?;
                let bc = self.edge_vertex(nodes[b], nodes[c])?;
                self.push_triangle([ac, ad, bd], outward)?;
                self.push_triangle([ac, bd, bc], outward)
            }
            _ => unreachable!("tetrahedron has at most four vertices"),
        }
    }

    fn edge_vertex(&mut self, first: usize, second: usize) -> Result<u32, IsoSurfaceError> {
        let first_value = self.grid.values[first];
        let second_value = self.grid.values[second];
        let first_on_level = same_level3(first_value, self.iso);
        let second_on_level = same_level3(second_value, self.iso);
        let key = if first_on_level {
            VertexKey3::Node(first)
        } else if second_on_level {
            VertexKey3::Node(second)
        } else {
            VertexKey3::Edge(first.min(second), first.max(second))
        };
        if let Some(index) = self.vertex_by_edge.get(&key) {
            return Ok(*index);
        }

        let point = if first_on_level {
            self.grid.point_from_index(first)
        } else if second_on_level {
            self.grid.point_from_index(second)
        } else {
            let scale = first_value
                .abs()
                .max(second_value.abs())
                .max(self.iso.abs())
                .max(1.0);
            let first_distance = first_value / scale - self.iso / scale;
            let second_distance = second_value / scale - self.iso / scale;
            let denominator = first_distance - second_distance;
            if !denominator.is_finite() || denominator.abs() <= f64::MIN_POSITIVE {
                return Err(IsoSurfaceError::NonFiniteGeometry);
            }
            let t = (first_distance / denominator).clamp(0.0, 1.0);
            let a = self.grid.point_from_index(first);
            let b = self.grid.point_from_index(second);
            [
                t.mul_add(b[0] - a[0], a[0]),
                t.mul_add(b[1] - a[1], a[1]),
                t.mul_add(b[2] - a[2], a[2]),
            ]
        };
        if point.into_iter().any(|coordinate| !coordinate.is_finite()) {
            return Err(IsoSurfaceError::NonFiniteGeometry);
        }
        let index =
            u32::try_from(self.vertices.len()).map_err(|_| IsoSurfaceError::VertexIndexOverflow)?;
        self.vertices.push(point);
        self.vertex_by_edge.insert(key, index);
        Ok(index)
    }

    fn push_triangle(
        &mut self,
        mut triangle: [u32; 3],
        outward: Vec3,
    ) -> Result<(), IsoSurfaceError> {
        if triangle[0] == triangle[1] || triangle[1] == triangle[2] || triangle[0] == triangle[2] {
            return Ok(());
        }
        let a = self.vertices[triangle[0] as usize];
        let b = self.vertices[triangle[1] as usize];
        let c = self.vertices[triangle[2] as usize];
        let normal = cross3(subtract3(b, a), subtract3(c, a));
        let normal_squared = norm_squared3(normal);
        if !normal_squared.is_finite() {
            return Err(IsoSurfaceError::NonFiniteGeometry);
        }
        if normal_squared <= f64::MIN_POSITIVE {
            return Ok(());
        }
        if dot3(normal, outward) < 0.0 {
            triangle.swap(1, 2);
        }
        if self.triangles.len() == self.triangle_limit {
            return Err(IsoSurfaceError::TriangleBudgetExceeded {
                limit: self.triangle_limit,
            });
        }
        self.triangles.push(triangle);
        Ok(())
    }

    fn outward_direction(&self, nodes: [usize; 4], inside: &[usize], outside: &[usize]) -> Vec3 {
        let centroid = |locals: &[usize]| {
            let mut sum = [0.0; 3];
            for local in locals {
                let point = self.grid.point_from_index(nodes[*local]);
                for axis in 0..3 {
                    sum[axis] += point[axis];
                }
            }
            let inverse_count = 1.0 / locals.len() as f64;
            sum.map(|value| value * inverse_count)
        };
        subtract3(centroid(outside), centroid(inside))
    }
}

fn same_level3(value: f64, iso: f64) -> bool {
    value.to_bits() == iso.to_bits()
        || (matches!(value.classify(), core::num::FpCategory::Zero)
            && matches!(iso.classify(), core::num::FpCategory::Zero))
}

fn subtract3(left: Vec3, right: Vec3) -> Vec3 {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn cross3(left: Vec3, right: Vec3) -> Vec3 {
    [
        left[1].mul_add(right[2], -left[2] * right[1]),
        left[2].mul_add(right[0], -left[0] * right[2]),
        left[0].mul_add(right[1], -left[1] * right[0]),
    ]
}

fn dot3(left: Vec3, right: Vec3) -> f64 {
    left[0].mul_add(right[0], left[1].mul_add(right[1], left[2] * right[2]))
}

fn norm_squared3(vector: Vec3) -> f64 {
    dot3(vector, vector)
}
