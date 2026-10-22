//! Render the accepted material body's actual P1 tetrahedral boundary.
//!
//! The material/solid owners supply the positions and topology. The existing
//! solid load geometry extracts exterior faces; the existing renderer builds
//! its usual triangle BVH. No display displacement, second deformation solve,
//! heat state, vertex welding or surface remeshing is introduced.
//!
//! Coordinates remain in the body's metre frame. A scene placement can use the
//! existing instance transform. Supply appearance separately through the
//! original material-state bindings and replace a scene shape only after this
//! function succeeds. Regional appearances and phase-changing topology are not
//! inferred from this fixed-topology surface.

use fs_exec::Cx;
use fs_mbd_material::deformation::DeformedBody;
use fs_mbd_material::deformation::thermoelastic::ThermalBody;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::{LoadError3, TetLoadGeometry3};

use crate::charts::TriMesh;

/// A refusal before a complete rendered shape is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialGeometryError {
    /// The original solid boundary owner rejected the mesh or its budget.
    Boundary(LoadError3),
    /// Vertex, triangle or BVH indices exceed the original renderer's capacity.
    RendererCapacity,
    /// The original triangle backend cannot represent this face's unit normal.
    RendererNormal {
        /// Original shared vertex IDs of the refused exterior face.
        nodes: [usize; 3],
    },
    /// Cancellation requested before the completed mesh could be returned.
    Cancelled,
}

impl core::fmt::Display for MaterialGeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Boundary(error) => error.fmt(f),
            Self::RendererCapacity => {
                f.write_str("material surface exceeds renderer index capacity")
            }
            Self::RendererNormal { nodes } => {
                write!(
                    f,
                    "material face {nodes:?} is outside the triangle normal domain"
                )
            }
            Self::Cancelled => f.write_str("material surface construction cancelled"),
        }
    }
}

impl core::error::Error for MaterialGeometryError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Boundary(error) => Some(error),
            _ => None,
        }
    }
}

impl From<LoadError3> for MaterialGeometryError {
    fn from(error: LoadError3) -> Self {
        match error {
            LoadError3::Cancelled => Self::Cancelled,
            error => Self::Boundary(error),
        }
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), MaterialGeometryError> {
    cx.checkpoint()
        .map_err(|_| MaterialGeometryError::Cancelled)
}

/// Build the exact exterior of an accepted, mass-conserving deformation.
///
/// Vertices retain the original shared node IDs and bitwise positions. Internal
/// faces are omitted and exterior triangles wind outward, independently of tet
/// connectivity winding. Normals, bounds and visibility are recomputed by the
/// original triangle backend from these actual positions.
///
/// `budget` and `maximum_face_entries` have the same meanings as in
/// [`TetLoadGeometry3::new`]; no stiffness matrix is assembled. The caller must
/// supply a conforming, non-overlapping P1 volume mesh: the solid face-incidence
/// checks do not prove global injectivity or detect T-junctions.
///
/// Cancellation is checked per copied vertex and exterior face, in the solid
/// boundary extractor, and before/after the original BVH build. That last build
/// is a bounded synchronous section with no internal cancellation checkpoints.
/// Neither the body nor any scene is mutated, including on cancellation.
///
/// # Errors
/// Refuses the original solid owner's invalid geometry and budget cases, index
/// overflow, and faces outside the current `TriMesh` normal domain (finite
/// cross-product norm strictly greater than `1e-12` square metres). No rescaling
/// or substitute geometry is used to hide an unsupported scale.
pub fn deformed_body_mesh(
    body: &DeformedBody<'_>,
    budget: TetAssemblyBudget,
    maximum_face_entries: usize,
    cx: &Cx<'_>,
) -> Result<TriMesh, MaterialGeometryError> {
    poll(cx)?;
    u32::try_from(body.nodes_m().len()).map_err(|_| MaterialGeometryError::RendererCapacity)?;
    let boundary = TetLoadGeometry3::new(
        body.nodes_m(),
        body.tetrahedra(),
        budget,
        maximum_face_entries,
        cx,
    )?;
    // A binary BVH has fewer than twice as many nodes as triangles. Its child
    // indices and the triangle permutation are both stored as u32.
    boundary
        .faces()
        .len()
        .checked_mul(2)
        .and_then(|count| u32::try_from(count).ok())
        .ok_or(MaterialGeometryError::RendererCapacity)?;
    let mut triangles = Vec::with_capacity(boundary.faces().len());
    for face in boundary.faces() {
        poll(cx)?;
        let [a, b, c] = face.nodes.map(|node| body.nodes_m()[node]);
        let ab: [f64; 3] = core::array::from_fn(|i| b[i] - a[i]);
        let ac: [f64; 3] = core::array::from_fn(|i| c[i] - a[i]);
        let normal = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        let length = normal.iter().map(|n| n * n).sum::<f64>().sqrt();
        if !length.is_finite() || length <= 1.0e-12 {
            return Err(MaterialGeometryError::RendererNormal { nodes: face.nodes });
        }
        let alignment = (0..3)
            .map(|i| (normal[i] / length) * face.outward_normal[i])
            .sum::<f64>();
        if !alignment.is_finite() || alignment == 0.0 {
            return Err(MaterialGeometryError::RendererNormal { nodes: face.nodes });
        }
        // The boundary owner returns sorted triples, not oriented triangles.
        // All IDs fit because its connectivity check and the vertex count check
        // above precede this conversion.
        let mut triangle = face.nodes.map(|node| node as u32);
        if alignment < 0.0 {
            triangle.swap(1, 2);
        }
        triangles.push(triangle);
    }
    let mut vertices = Vec::with_capacity(body.nodes_m().len());
    for &node in body.nodes_m() {
        poll(cx)?;
        vertices.push(node);
    }
    poll(cx)?;
    let mesh = TriMesh::new(vertices, triangles);
    poll(cx)?;
    Ok(mesh)
}

/// Observe the geometry already accepted by the original thermoelastic solve.
/// This does not query a second temperature or integrate heat. Appearance can
/// bind to the same result's `sources()` through `MaterialThermalBinding`.
pub fn thermal_body_mesh(
    body: &ThermalBody<'_>,
    budget: TetAssemblyBudget,
    maximum_face_entries: usize,
    cx: &Cx<'_>,
) -> Result<TriMesh, MaterialGeometryError> {
    deformed_body_mesh(body.body(), budget, maximum_face_entries, cx)
}
