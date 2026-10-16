//! Affine-triangle arbitrary-Lagrangian-Eulerian geometry audits.
//!
//! This module measures the discrete geometric-conservation-law (GCL)
//! identity for a fixed-connectivity [`TriMesh`] whose vertices move
//! linearly during one time step.  It does not advance a flow solution or
//! remap a field; it only returns the geometric evidence a later ALE
//! operator can bind into its own receipt.

use crate::TriMesh;
use core::fmt;

const MIN_TRIANGLE_AREA: f64 = 1.0e-14;

/// Which endpoint of the admitted mesh trajectory produced a geometry error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AleMeshState2 {
    /// The `TriMesh` coordinates at the start of the time step.
    Initial,
    /// The caller-supplied coordinates at the end of the time step.
    Final,
}

/// Structural defect found in the public `TriMesh` tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AleGclMeshIssue2 {
    /// A per-triangle table has a different length than `tris`.
    TriangleTableLength,
    /// A triangle references a vertex outside `verts`.
    TriangleVertex,
    /// A triangle repeats a vertex.
    RepeatedTriangleVertex,
    /// A stored triangle area differs from the canonical coordinate value.
    StoredArea,
    /// A stored triangle centroid differs from the canonical coordinate value.
    StoredCentroid,
    /// A local triangle edge references an edge outside `edges`.
    LocalEdgeIndex,
    /// A triangle repeats a global edge in two local slots.
    RepeatedLocalEdge,
    /// A local edge sign is not exactly `+1` or `-1`.
    LocalEdgeSign,
    /// A local opposite edge and the global edge have different vertices.
    LocalEdgeVertices,
    /// A local sign and the global owner/neighbor table disagree.
    LocalEdgeAdjacency,
    /// An interior neighbor traverses the shared edge in the owner direction.
    InteriorEdgeOrientation,
    /// A global edge has invalid or non-canonical endpoint indices.
    EdgeVertices,
    /// Stored length, midpoint, or owner-normal differs from the coordinates.
    StoredEdgeGeometry,
    /// Two global edge rows describe the same endpoint pair.
    DuplicateGlobalEdge,
    /// A global edge has an invalid owner or neighbor triangle.
    EdgeAdjacency,
}

/// Scratch collection whose fallible reservation failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AleGclCollection2 {
    /// Per-cell GCL evidence.
    Cells,
    /// Per-edge swept-area evidence.
    Edges,
}

/// Stage at which finite inputs overflowed finite arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AleGclArithmetic2 {
    /// Initial or final signed-area evaluation.
    CellArea,
    /// Canonical centroid evaluation for the initial mesh.
    CellCentroid,
    /// Continuous-in-time minimum-area evaluation.
    CellTrajectory,
    /// Oriented edge swept-area evaluation.
    EdgeSweep,
    /// Per-cell accumulation of signed edge sweeps.
    CellBalance,
    /// Mesh-wide reduction.
    GlobalBalance,
}

/// Typed refusal returned by [`audit_affine_triangle_gcl2`].
#[derive(Debug, Clone, PartialEq)]
pub enum AleGclError2 {
    /// No vertices, triangles, or edges were supplied.
    EmptyMesh,
    /// The time step is non-finite or not strictly positive.
    InvalidTimeStep {
        /// Caller-supplied time step in coherent SI seconds.
        seconds: f64,
    },
    /// The final coordinate table does not match the initial vertex count.
    VertexCountMismatch {
        /// Number of initial vertices.
        expected: usize,
        /// Number of final vertices.
        actual: usize,
    },
    /// A coordinate is not finite.
    NonFiniteCoordinate {
        /// Initial or final coordinate table.
        state: AleMeshState2,
        /// Vertex index.
        vertex: usize,
        /// Cartesian axis, zero for x and one for y.
        axis: usize,
    },
    /// The public `TriMesh` tables are internally inconsistent.
    InvalidMesh {
        /// Kind of structural inconsistency.
        issue: AleGclMeshIssue2,
        /// Triangle or edge index associated with the inconsistency.
        index: usize,
    },
    /// A cell is inverted or at/below the canonical area floor.
    DegenerateCell {
        /// Initial or final endpoint of the motion.
        state: AleMeshState2,
        /// Triangle index.
        cell: usize,
        /// Signed area at that endpoint.
        signed_area: f64,
    },
    /// A cell collapses or inverts between two admissible endpoints.
    TrajectoryCollapse {
        /// Triangle index.
        cell: usize,
        /// Normalized time in the closed interval `[0, 1]`.
        normalized_time: f64,
        /// Minimum signed area on the linear trajectory.
        signed_area: f64,
    },
    /// Finite coordinates overflowed during a derived calculation.
    NonFiniteArithmetic {
        /// Calculation that overflowed.
        stage: AleGclArithmetic2,
        /// Cell or edge index associated with the calculation.
        index: usize,
    },
    /// A receipt collection could not reserve its exact required capacity.
    AllocationFailed {
        /// Collection whose reservation failed.
        collection: AleGclCollection2,
        /// Exact requested element count.
        requested: usize,
    },
}

impl fmt::Display for AleGclError2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMesh => f.write_str("ALE GCL audit requires a non-empty triangle mesh"),
            Self::InvalidTimeStep { seconds } => write!(
                f,
                "ALE GCL time_step_seconds must be finite and positive, got {seconds}"
            ),
            Self::VertexCountMismatch { expected, actual } => write!(
                f,
                "ALE GCL final vertex count {actual} does not match initial count {expected}"
            ),
            Self::NonFiniteCoordinate {
                state,
                vertex,
                axis,
            } => write!(
                f,
                "ALE GCL {state:?} vertex {vertex} axis {axis} is not finite"
            ),
            Self::InvalidMesh { issue, index } => {
                write!(f, "ALE GCL mesh issue {issue:?} at index {index}")
            }
            Self::DegenerateCell {
                state,
                cell,
                signed_area,
            } => write!(
                f,
                "ALE GCL {state:?} cell {cell} has inadmissible signed area {signed_area}"
            ),
            Self::TrajectoryCollapse {
                cell,
                normalized_time,
                signed_area,
            } => write!(
                f,
                "ALE GCL cell {cell} collapses at normalized time {normalized_time} with signed area {signed_area}"
            ),
            Self::NonFiniteArithmetic { stage, index } => write!(
                f,
                "ALE GCL {stage:?} arithmetic became non-finite at index {index}"
            ),
            Self::AllocationFailed {
                collection,
                requested,
            } => write!(
                f,
                "ALE GCL could not reserve {requested} {collection:?} receipt entries"
            ),
        }
    }
}

impl std::error::Error for AleGclError2 {}

/// Per-triangle integrated geometric-conservation evidence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellGcl2 {
    area_before: f64,
    area_after: f64,
    minimum_area: f64,
    outward_swept_area: f64,
    area_defect: f64,
}

impl CellGcl2 {
    /// Initial signed area in squared mesh-coordinate units.
    #[must_use]
    pub const fn area_before(&self) -> f64 {
        self.area_before
    }

    /// Final signed area in squared mesh-coordinate units.
    #[must_use]
    pub const fn area_after(&self) -> f64 {
        self.area_after
    }

    /// Minimum signed area over the trajectory, in squared coordinate units.
    #[must_use]
    pub const fn minimum_area(&self) -> f64 {
        self.minimum_area
    }

    /// Time-integrated outward mesh flux, in squared coordinate units.
    #[must_use]
    pub const fn outward_swept_area(&self) -> f64 {
        self.outward_swept_area
    }

    /// `(area_after - area_before) - outward_swept_area`.
    #[must_use]
    pub const fn area_defect(&self) -> f64 {
        self.area_defect
    }
}

/// Deterministic integrated GCL evidence for one affine-triangle mesh step.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct AleGclReceipt2 {
    time_step_seconds: f64,
    cells: Vec<CellGcl2>,
    owner_edge_swept_areas: Vec<f64>,
    boundary_edges: usize,
    interior_edges: usize,
    total_area_change: f64,
    boundary_swept_area: f64,
    total_area_defect: f64,
    worst_abs_cell_defect: f64,
    worst_abs_interior_cancellation: f64,
}

impl AleGclReceipt2 {
    /// Admitted motion duration in coherent SI seconds.
    #[must_use]
    pub const fn time_step_seconds(&self) -> f64 {
        self.time_step_seconds
    }

    /// Per-cell evidence in triangle-index order.
    #[must_use]
    pub fn cells(&self) -> &[CellGcl2] {
        &self.cells
    }

    /// Per-edge swept areas in squared coordinate units, oriented from the owner.
    #[must_use]
    pub fn owner_edge_swept_areas(&self) -> &[f64] {
        &self.owner_edge_swept_areas
    }

    /// Number of boundary edges admitted by the audit.
    #[must_use]
    pub const fn boundary_edges(&self) -> usize {
        self.boundary_edges
    }

    /// Number of two-cell interior edges admitted by the audit.
    #[must_use]
    pub const fn interior_edges(&self) -> usize {
        self.interior_edges
    }

    /// Sum of all final cell areas minus all initial cell areas.
    #[must_use]
    pub const fn total_area_change(&self) -> f64 {
        self.total_area_change
    }

    /// Sum of the owner-oriented swept areas on boundary edges.
    #[must_use]
    pub const fn boundary_swept_area(&self) -> f64 {
        self.boundary_swept_area
    }

    /// `total_area_change - boundary_swept_area`.
    #[must_use]
    pub const fn total_area_defect(&self) -> f64 {
        self.total_area_defect
    }

    /// Maximum absolute per-cell GCL defect.
    #[must_use]
    pub const fn worst_abs_cell_defect(&self) -> f64 {
        self.worst_abs_cell_defect
    }

    /// Maximum absolute sum of the two signed contributions of one interior edge.
    #[must_use]
    pub const fn worst_abs_interior_cancellation(&self) -> f64 {
        self.worst_abs_interior_cancellation
    }
}

/// Audit the integrated GCL identity for linearly moving triangle vertices.
///
/// The input topology is fixed.  For every owner-oriented edge, the audit
/// exactly integrates the normal mesh velocity of its linearly moving
/// endpoints.  Each cell then compares the signed sum of those swept areas
/// with its measured area change.  Finite defects are evidence and are
/// returned rather than rejected, so this function does not certify itself.
///
/// `time_step_seconds` is recorded explicitly; the returned swept quantities
/// are integrated over that step rather than divided into rates.
///
/// # Errors
///
/// Returns a typed refusal for malformed public mesh tables, non-finite data,
/// vertex-count mismatch, endpoint or mid-step collapse, arithmetic overflow,
/// or receipt-allocation failure.
#[allow(clippy::too_many_lines)] // One ordered validate-measure-reduce receipt transaction.
pub fn audit_affine_triangle_gcl2(
    mesh: &TriMesh,
    next_vertices: &[[f64; 2]],
    time_step_seconds: f64,
) -> Result<AleGclReceipt2, AleGclError2> {
    if !time_step_seconds.is_finite() || time_step_seconds <= 0.0 {
        return Err(AleGclError2::InvalidTimeStep {
            seconds: time_step_seconds,
        });
    }
    if mesh.verts.is_empty() || mesh.tris.is_empty() || mesh.edges.is_empty() {
        return Err(AleGclError2::EmptyMesh);
    }
    if next_vertices.len() != mesh.verts.len() {
        return Err(AleGclError2::VertexCountMismatch {
            expected: mesh.verts.len(),
            actual: next_vertices.len(),
        });
    }
    validate_coordinates(&mesh.verts, AleMeshState2::Initial)?;
    validate_coordinates(next_vertices, AleMeshState2::Final)?;
    validate_mesh_tables(mesh)?;

    let mut cells = Vec::new();
    cells
        .try_reserve_exact(mesh.tris.len())
        .map_err(|_| AleGclError2::AllocationFailed {
            collection: AleGclCollection2::Cells,
            requested: mesh.tris.len(),
        })?;

    for (cell, tri) in mesh.tris.iter().enumerate() {
        let before = tri_points(&mesh.verts, *tri);
        let after = tri_points(next_vertices, *tri);
        let area_before = checked_area(before, AleMeshState2::Initial, cell)?;
        let area_after = checked_area(after, AleMeshState2::Final, cell)?;
        let minimum_area = minimum_trajectory_area(before, after, cell)?;
        cells.push(CellGcl2 {
            area_before,
            area_after,
            minimum_area,
            outward_swept_area: 0.0,
            area_defect: 0.0,
        });
    }

    let mut edge_sweeps = Vec::new();
    edge_sweeps
        .try_reserve_exact(mesh.edges.len())
        .map_err(|_| AleGclError2::AllocationFailed {
            collection: AleGclCollection2::Edges,
            requested: mesh.edges.len(),
        })?;
    for (edge_id, edge) in mesh.edges.iter().enumerate() {
        let owner = edge.tris.0;
        let Some(local) = local_edge(mesh, owner, edge_id) else {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::EdgeAdjacency,
                index: edge_id,
            });
        };
        let tri = mesh.tris[owner];
        let a = tri[(local + 1) % 3];
        let b = tri[(local + 2) % 3];
        let swept = swept_area(
            mesh.verts[a],
            mesh.verts[b],
            next_vertices[a],
            next_vertices[b],
        )
        .ok_or(AleGclError2::NonFiniteArithmetic {
            stage: AleGclArithmetic2::EdgeSweep,
            index: edge_id,
        })?;
        edge_sweeps.push(canonical_zero(swept));
    }

    let mut total_area_change = 0.0;
    let mut worst_abs_cell_defect = 0.0f64;
    for (cell, evidence) in cells.iter_mut().enumerate() {
        let mut swept = 0.0;
        for (edge, sign) in mesh.tri_edges[cell] {
            swept += sign * edge_sweeps[edge];
            if !swept.is_finite() {
                return Err(AleGclError2::NonFiniteArithmetic {
                    stage: AleGclArithmetic2::CellBalance,
                    index: cell,
                });
            }
        }
        let change = evidence.area_after - evidence.area_before;
        let defect = change - swept;
        if !change.is_finite() || !defect.is_finite() {
            return Err(AleGclError2::NonFiniteArithmetic {
                stage: AleGclArithmetic2::CellBalance,
                index: cell,
            });
        }
        evidence.outward_swept_area = canonical_zero(swept);
        evidence.area_defect = canonical_zero(defect);
        worst_abs_cell_defect = worst_abs_cell_defect.max(defect.abs());
        total_area_change += change;
        if !total_area_change.is_finite() {
            return Err(AleGclError2::NonFiniteArithmetic {
                stage: AleGclArithmetic2::GlobalBalance,
                index: cell,
            });
        }
    }

    let mut boundary_swept_area = 0.0;
    let mut boundary_edges = 0usize;
    let mut interior_edges = 0usize;
    let mut worst_abs_interior_cancellation = 0.0f64;
    for (edge_id, edge) in mesh.edges.iter().enumerate() {
        let swept = edge_sweeps[edge_id];
        if edge.tris.1 == usize::MAX {
            boundary_edges += 1;
            boundary_swept_area += swept;
            if !boundary_swept_area.is_finite() {
                return Err(AleGclError2::NonFiniteArithmetic {
                    stage: AleGclArithmetic2::GlobalBalance,
                    index: edge_id,
                });
            }
        } else {
            interior_edges += 1;
            let Some(neighbor_local) = local_edge(mesh, edge.tris.1, edge_id) else {
                return Err(AleGclError2::InvalidMesh {
                    issue: AleGclMeshIssue2::EdgeAdjacency,
                    index: edge_id,
                });
            };
            let neighbor_tri = mesh.tris[edge.tris.1];
            let neighbor_a = neighbor_tri[(neighbor_local + 1) % 3];
            let neighbor_b = neighbor_tri[(neighbor_local + 2) % 3];
            let neighbor_swept = swept_area(
                mesh.verts[neighbor_a],
                mesh.verts[neighbor_b],
                next_vertices[neighbor_a],
                next_vertices[neighbor_b],
            )
            .ok_or(AleGclError2::NonFiniteArithmetic {
                stage: AleGclArithmetic2::EdgeSweep,
                index: edge_id,
            })?;
            let cancellation = swept + neighbor_swept;
            if !cancellation.is_finite() {
                return Err(AleGclError2::NonFiniteArithmetic {
                    stage: AleGclArithmetic2::EdgeSweep,
                    index: edge_id,
                });
            }
            worst_abs_interior_cancellation =
                worst_abs_interior_cancellation.max(cancellation.abs());
        }
    }
    let total_area_defect = total_area_change - boundary_swept_area;
    if !total_area_defect.is_finite() {
        return Err(AleGclError2::NonFiniteArithmetic {
            stage: AleGclArithmetic2::GlobalBalance,
            index: mesh.tris.len(),
        });
    }

    Ok(AleGclReceipt2 {
        time_step_seconds,
        cells,
        owner_edge_swept_areas: edge_sweeps,
        boundary_edges,
        interior_edges,
        total_area_change: canonical_zero(total_area_change),
        boundary_swept_area: canonical_zero(boundary_swept_area),
        total_area_defect: canonical_zero(total_area_defect),
        worst_abs_cell_defect: canonical_zero(worst_abs_cell_defect),
        worst_abs_interior_cancellation: canonical_zero(worst_abs_interior_cancellation),
    })
}

fn validate_coordinates(vertices: &[[f64; 2]], state: AleMeshState2) -> Result<(), AleGclError2> {
    for (vertex, point) in vertices.iter().enumerate() {
        for (axis, value) in point.iter().enumerate() {
            if !value.is_finite() {
                return Err(AleGclError2::NonFiniteCoordinate {
                    state,
                    vertex,
                    axis,
                });
            }
        }
    }
    Ok(())
}

fn validate_mesh_tables(mesh: &TriMesh) -> Result<(), AleGclError2> {
    if mesh.tri_edges.len() != mesh.tris.len()
        || mesh.areas.len() != mesh.tris.len()
        || mesh.centroids.len() != mesh.tris.len()
    {
        return Err(AleGclError2::InvalidMesh {
            issue: AleGclMeshIssue2::TriangleTableLength,
            index: mesh.tris.len(),
        });
    }
    validate_triangle_tables(mesh)?;
    validate_edge_tables(mesh)
}

fn validate_triangle_tables(mesh: &TriMesh) -> Result<(), AleGclError2> {
    for (cell, tri) in mesh.tris.iter().enumerate() {
        if tri.iter().any(|&vertex| vertex >= mesh.verts.len()) {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::TriangleVertex,
                index: cell,
            });
        }
        if tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0] {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::RepeatedTriangleVertex,
                index: cell,
            });
        }
        let points = tri_points(&mesh.verts, *tri);
        let area = signed_area(points);
        if !area.is_finite() {
            return Err(AleGclError2::NonFiniteArithmetic {
                stage: AleGclArithmetic2::CellArea,
                index: cell,
            });
        }
        if area <= MIN_TRIANGLE_AREA {
            return Err(AleGclError2::DegenerateCell {
                state: AleMeshState2::Initial,
                cell,
                signed_area: area,
            });
        }
        if mesh.areas[cell].to_bits() != area.to_bits() {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::StoredArea,
                index: cell,
            });
        }
        let centroid = [
            (points[0][0] + points[1][0] + points[2][0]) / 3.0,
            (points[0][1] + points[1][1] + points[2][1]) / 3.0,
        ];
        if centroid.iter().any(|value| !value.is_finite()) {
            return Err(AleGclError2::NonFiniteArithmetic {
                stage: AleGclArithmetic2::CellCentroid,
                index: cell,
            });
        }
        if mesh.centroids[cell].iter().any(|value| !value.is_finite())
            || mesh.centroids[cell][0].to_bits() != centroid[0].to_bits()
            || mesh.centroids[cell][1].to_bits() != centroid[1].to_bits()
        {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::StoredCentroid,
                index: cell,
            });
        }

        let local_edges = mesh.tri_edges[cell];
        if local_edges[0].0 == local_edges[1].0
            || local_edges[1].0 == local_edges[2].0
            || local_edges[2].0 == local_edges[0].0
        {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::RepeatedLocalEdge,
                index: cell,
            });
        }
        for (local, (edge_id, sign)) in local_edges.into_iter().enumerate() {
            let Some(edge) = mesh.edges.get(edge_id) else {
                return Err(AleGclError2::InvalidMesh {
                    issue: AleGclMeshIssue2::LocalEdgeIndex,
                    index: cell,
                });
            };
            if sign.to_bits() != 1.0f64.to_bits() && sign.to_bits() != (-1.0f64).to_bits() {
                return Err(AleGclError2::InvalidMesh {
                    issue: AleGclMeshIssue2::LocalEdgeSign,
                    index: cell,
                });
            }
            let a = tri[(local + 1) % 3];
            let b = tri[(local + 2) % 3];
            if edge.verts != (a.min(b), a.max(b)) {
                return Err(AleGclError2::InvalidMesh {
                    issue: AleGclMeshIssue2::LocalEdgeVertices,
                    index: cell,
                });
            }
            let adjacency_matches = if sign.is_sign_positive() {
                edge.tris.0 == cell
            } else {
                edge.tris.1 == cell
            };
            if !adjacency_matches {
                return Err(AleGclError2::InvalidMesh {
                    issue: AleGclMeshIssue2::LocalEdgeAdjacency,
                    index: cell,
                });
            }
        }
    }
    Ok(())
}

fn validate_edge_tables(mesh: &TriMesh) -> Result<(), AleGclError2> {
    let mut canonical_edges = Vec::new();
    canonical_edges
        .try_reserve_exact(mesh.edges.len())
        .map_err(|_| AleGclError2::AllocationFailed {
            collection: AleGclCollection2::Edges,
            requested: mesh.edges.len(),
        })?;
    for (edge_id, edge) in mesh.edges.iter().enumerate() {
        if edge.verts.0 >= edge.verts.1 || edge.verts.1 >= mesh.verts.len() {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::EdgeVertices,
                index: edge_id,
            });
        }
        canonical_edges.push((edge.verts, edge_id));
        if edge.tris.0 >= mesh.tris.len()
            || (edge.tris.1 != usize::MAX
                && (edge.tris.1 >= mesh.tris.len() || edge.tris.0 >= edge.tris.1))
        {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::EdgeAdjacency,
                index: edge_id,
            });
        }
        let Some(owner) = local_edge(mesh, edge.tris.0, edge_id) else {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::EdgeAdjacency,
                index: edge_id,
            });
        };
        if mesh.tri_edges[edge.tris.0][owner].1.to_bits() != 1.0f64.to_bits() {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::EdgeAdjacency,
                index: edge_id,
            });
        }
        let owner_directed = validate_stored_edge_geometry(mesh, edge_id, owner)?;
        if edge.tris.1 != usize::MAX {
            let Some(neighbor) = local_edge(mesh, edge.tris.1, edge_id) else {
                return Err(AleGclError2::InvalidMesh {
                    issue: AleGclMeshIssue2::EdgeAdjacency,
                    index: edge_id,
                });
            };
            if mesh.tri_edges[edge.tris.1][neighbor].1.to_bits() != (-1.0f64).to_bits() {
                return Err(AleGclError2::InvalidMesh {
                    issue: AleGclMeshIssue2::EdgeAdjacency,
                    index: edge_id,
                });
            }
            let neighbor_tri = mesh.tris[edge.tris.1];
            let neighbor_directed = (
                neighbor_tri[(neighbor + 1) % 3],
                neighbor_tri[(neighbor + 2) % 3],
            );
            if neighbor_directed != (owner_directed.1, owner_directed.0) {
                return Err(AleGclError2::InvalidMesh {
                    issue: AleGclMeshIssue2::InteriorEdgeOrientation,
                    index: edge_id,
                });
            }
        }
    }
    canonical_edges.sort_unstable_by_key(|&(vertices, _)| vertices);
    for pair in canonical_edges.windows(2) {
        if pair[0].0 == pair[1].0 {
            return Err(AleGclError2::InvalidMesh {
                issue: AleGclMeshIssue2::DuplicateGlobalEdge,
                index: pair[1].1,
            });
        }
    }
    Ok(())
}

fn validate_stored_edge_geometry(
    mesh: &TriMesh,
    edge_id: usize,
    owner_local: usize,
) -> Result<(usize, usize), AleGclError2> {
    let edge = &mesh.edges[edge_id];
    let owner_tri = mesh.tris[edge.tris.0];
    let owner_directed = (
        owner_tri[(owner_local + 1) % 3],
        owner_tri[(owner_local + 2) % 3],
    );
    let sorted_a = mesh.verts[edge.verts.0];
    let sorted_b = mesh.verts[edge.verts.1];
    let dx = sorted_b[0] - sorted_a[0];
    let dy = sorted_b[1] - sorted_a[1];
    let canonical_len = dx.hypot(dy);
    let owner_a = mesh.verts[owner_directed.0];
    let owner_b = mesh.verts[owner_directed.1];
    let owner_dx = owner_b[0] - owner_a[0];
    let owner_dy = owner_b[1] - owner_a[1];
    let canonical_normal = [owner_dy / canonical_len, -owner_dx / canonical_len];
    let canonical_mid = [
        f64::midpoint(sorted_a[0], sorted_b[0]),
        f64::midpoint(sorted_a[1], sorted_b[1]),
    ];
    if !canonical_len.is_finite()
        || canonical_normal.iter().any(|value| !value.is_finite())
        || canonical_mid.iter().any(|value| !value.is_finite())
    {
        return Err(AleGclError2::NonFiniteArithmetic {
            stage: AleGclArithmetic2::EdgeSweep,
            index: edge_id,
        });
    }
    if edge.len.to_bits() != canonical_len.to_bits()
        || edge.normal[0].to_bits() != canonical_normal[0].to_bits()
        || edge.normal[1].to_bits() != canonical_normal[1].to_bits()
        || edge.mid[0].to_bits() != canonical_mid[0].to_bits()
        || edge.mid[1].to_bits() != canonical_mid[1].to_bits()
    {
        return Err(AleGclError2::InvalidMesh {
            issue: AleGclMeshIssue2::StoredEdgeGeometry,
            index: edge_id,
        });
    }
    Ok(owner_directed)
}

fn local_edge(mesh: &TriMesh, cell: usize, edge: usize) -> Option<usize> {
    mesh.tri_edges[cell]
        .iter()
        .position(|&(candidate, _)| candidate == edge)
}

fn tri_points(vertices: &[[f64; 2]], tri: [usize; 3]) -> [[f64; 2]; 3] {
    [vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]]
}

fn signed_area(points: [[f64; 2]; 3]) -> f64 {
    0.5 * ((points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
        - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1]))
}

fn checked_area(
    points: [[f64; 2]; 3],
    state: AleMeshState2,
    cell: usize,
) -> Result<f64, AleGclError2> {
    let area = signed_area(points);
    if !area.is_finite() {
        return Err(AleGclError2::NonFiniteArithmetic {
            stage: AleGclArithmetic2::CellArea,
            index: cell,
        });
    }
    if area <= MIN_TRIANGLE_AREA {
        return Err(AleGclError2::DegenerateCell {
            state,
            cell,
            signed_area: area,
        });
    }
    Ok(area)
}

fn minimum_trajectory_area(
    before: [[f64; 2]; 3],
    after: [[f64; 2]; 3],
    cell: usize,
) -> Result<f64, AleGclError2> {
    let displacement = [
        sub(after[0], before[0]),
        sub(after[1], before[1]),
        sub(after[2], before[2]),
    ];
    let u0 = sub(before[1], before[0]);
    let v0 = sub(before[2], before[0]);
    let du = sub(displacement[1], displacement[0]);
    let dv = sub(displacement[2], displacement[0]);
    let c0 = cross(u0, v0);
    let c1 = cross(du, v0) + cross(u0, dv);
    let c2 = cross(du, dv);
    if [c0, c1, c2].iter().any(|value| !value.is_finite()) {
        return Err(AleGclError2::NonFiniteArithmetic {
            stage: AleGclArithmetic2::CellTrajectory,
            index: cell,
        });
    }

    let mut minimum_double_area = c0.min(2.0 * signed_area(after));
    let mut minimum_time = if c0 <= 2.0 * signed_area(after) {
        0.0
    } else {
        1.0
    };
    if c2 > 0.0 {
        let stationary_time = -0.5 * (c1 / c2);
        if stationary_time > 0.0 && stationary_time < 1.0 {
            let candidate = (c2 * stationary_time + c1) * stationary_time + c0;
            if !candidate.is_finite() {
                return Err(AleGclError2::NonFiniteArithmetic {
                    stage: AleGclArithmetic2::CellTrajectory,
                    index: cell,
                });
            }
            if candidate < minimum_double_area {
                minimum_double_area = candidate;
                minimum_time = stationary_time;
            }
        }
    }
    let minimum_area = 0.5 * minimum_double_area;
    if !minimum_area.is_finite() {
        return Err(AleGclError2::NonFiniteArithmetic {
            stage: AleGclArithmetic2::CellTrajectory,
            index: cell,
        });
    }
    if minimum_area <= MIN_TRIANGLE_AREA {
        return Err(AleGclError2::TrajectoryCollapse {
            cell,
            normalized_time: minimum_time,
            signed_area: minimum_area,
        });
    }
    Ok(minimum_area)
}

fn swept_area(a0: [f64; 2], b0: [f64; 2], a1: [f64; 2], b1: [f64; 2]) -> Option<f64> {
    let da = sub(a1, a0);
    let db = sub(b1, b0);
    let e0 = sub(b0, a0);
    let e1 = sub(b1, a1);
    let velocity_sum = [da[0] + db[0], da[1] + db[1]];
    let edge_sum = [e0[0] + e1[0], e0[1] + e1[1]];
    let swept = 0.25 * cross(velocity_sum, edge_sum);
    swept.is_finite().then_some(swept)
}

fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn canonical_zero(value: f64) -> f64 {
    if value.to_bits() == (-0.0f64).to_bits() {
        0.0
    } else {
        value
    }
}
