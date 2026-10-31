//! Preserve material-region ownership on the actual exterior of a deformed body.
//!
//! Boundary extraction is global: bonded internal faces never become visible
//! just because their cells have different materials or temperatures. Only then
//! is the exterior partitioned by the original region-to-element assignment.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_mbd_material::deformation::DeformedBody;
use fs_solid::linear3::TetAssemblyBudget;

use super::{MaterialGeometryError, boundary_triangles, poll};
use crate::charts::TriMesh;

#[cfg(feature = "material-state")]
mod appearance;
#[cfg(feature = "material-state")]
pub use appearance::{
    MaterialRegionAppearance, MaterialRegionBinding, MaterialRegionBindingError,
    thermal_body_primitives,
};
#[cfg(feature = "material-state")]
pub use appearance::continuation::{
    MaterialRegionContinuation, MaterialRegionContinuationError, MaterialRegionSource,
    continue_thermal_body_primitives,
};

/// A nonempty exterior patch belonging to one original material cell.
pub struct MaterialRegionMesh {
    /// Index into the reference body's canonical region inventory. This is
    /// deliberately independent of the original tetrahedron's element index.
    pub region_index: usize,
    /// Original global node IDs corresponding to the four local mesh vertices.
    /// Coordinates are copied bitwise from the accepted body, without welding.
    pub source_nodes: [usize; 4],
    /// Exterior triangles and the original renderer's deterministic BVH.
    pub mesh: TriMesh,
}

/// Partition the globally admitted exterior by its original material regions.
///
/// Returns patches in canonical region order, omitting wholly internal cells.
/// Each patch stores only its four actual cell vertices and at most four faces;
/// the full body vertex array is not duplicated per region. The existing solid
/// node/tet/face-entry budgets bound extraction and the number of primitives.
/// The same conforming, non-overlapping, fixed-topology and triangle-scale
/// requirements as [`super::deformed_body_mesh`] apply.
///
/// # Errors
/// Propagates the original geometry, budget, scale and cancellation refusals.
/// No scene, body or material state is mutated on success or failure.
pub fn deformed_body_regions(
    body: &DeformedBody<'_>,
    budget: TetAssemblyBudget,
    maximum_face_entries: usize,
    cx: &Cx<'_>,
) -> Result<Vec<MaterialRegionMesh>, MaterialGeometryError> {
    let triangles = boundary_triangles(body, budget, maximum_face_entries, cx)?;
    let mut exterior = BTreeMap::new();
    for triangle in triangles {
        poll(cx)?;
        let mut key = triangle;
        key.sort_unstable();
        exterior.insert(key, triangle);
    }
    let mut regions = Vec::new();
    for (region_index, &element) in body.region_elements().iter().enumerate() {
        poll(cx)?;
        let source_nodes = body.tetrahedra()[element];
        let mut faces = Vec::with_capacity(4);
        for local in [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]] {
            poll(cx)?;
            let mut key = local.map(|i| source_nodes[i] as u32);
            key.sort_unstable();
            if let Some(triangle) = exterior.remove(&key) {
                // Every matched ID comes from this exact tet, admitted by the
                // original solid boundary owner. Preserve its outward winding.
                faces.push(triangle.map(|node| {
                    source_nodes
                        .iter()
                        .position(|&source| source == node as usize)
                        .expect("matched exterior face belongs to this cell")
                        as u32
                }));
            }
        }
        if faces.is_empty() {
            continue;
        }
        let vertices = source_nodes.map(|node| body.nodes_m()[node]).to_vec();
        poll(cx)?;
        let mesh = TriMesh::new(vertices, faces);
        poll(cx)?;
        regions.push(MaterialRegionMesh {
            region_index,
            source_nodes,
            mesh,
        });
    }
    // DeformedBody retains a complete bijective region-to-element assignment
    // from TetrahedralMass. No exterior face can lack its original owner.
    debug_assert!(exterior.is_empty());
    poll(cx)?;
    Ok(regions)
}
