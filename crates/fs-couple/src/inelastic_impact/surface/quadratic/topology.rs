//! Bind original P2 incidence to the existing closed-boundary topology owner.

use std::collections::BTreeMap;

use super::*;
use fs_exec::Cx;
use fs_solid::{
    hyper3::pressure::ClosedTetPressureSurface,
    linear3::{
        TetAssemblyBudget,
        tet10::{TET10_EDGES, tet10_point, validate_tet10_topology},
    },
};

impl QuadraticSolidContactSurface {
    /// Prepare the full immutable original mesh and its actual six-node exterior.
    /// The existing TET10 owner rejects mismatched shared midsides, hanging roles
    /// and unused nodes before internal faces are cancelled. Corner geometry is
    /// used only for original incidence and winding, never for a curved query.
    pub fn from_reference(
        nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 10]],
        budget: TetAssemblyBudget,
        maximum_face_entries: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, QuadraticSolidSurfaceError> {
        checkpoint(cx)?;
        if nodes_m.is_empty()
            || tetrahedra.is_empty()
            || nodes_m.len() > budget.maximum_nodes
            || tetrahedra.len() > budget.maximum_tetrahedra
            || tetrahedra
                .len()
                .checked_mul(4)
                .is_none_or(|count| count > maximum_face_entries)
            || tetrahedra.len().checked_mul(6).is_none()
            || nodes_m.len().checked_mul(3).is_none()
        {
            return Err(invalid("quadratic reference mesh exceeds topology budget"));
        }
        for (index, node) in nodes_m.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            vector(*node)?;
        }
        let vertices = validate_tet10_topology(nodes_m.len(), tetrahedra, cx)?;
        // Compact only corner indices because the original P1 incidence owner
        // correctly rejects unused nodes. The ordering preserves original IDs.
        let corner_ids: Vec<_> = vertices.into_iter().collect();
        let corner_nodes: Vec<_> = corner_ids.iter().map(|&node| nodes_m[node]).collect();
        let mut to_corner = vec![usize::MAX; nodes_m.len()];
        for (compact, &node) in corner_ids.iter().enumerate() {
            if compact % 256 == 0 {
                checkpoint(cx)?;
            }
            to_corner[node] = compact;
        }
        let mut corner_tets = Vec::with_capacity(tetrahedra.len());
        let mut owners = BTreeMap::new();
        for (tetrahedron, tet) in tetrahedra.iter().enumerate() {
            checkpoint(cx)?;
            // This is the native positive isoparametric map, not a surrogate
            // linear element. Global injectivity is not inferred from one sample.
            tet10_point(
                tet.map(|node| nodes_m[node]),
                [0.25; 4],
                budget.minimum_scaled_jacobian,
                tetrahedron,
                0,
            )?;
            corner_tets.push(core::array::from_fn(|local| to_corner[tet[local]]));
            for opposite in 0..4 {
                let mut key = [0; 3];
                let mut slot = 0;
                for (local, &node) in tet[..4].iter().enumerate() {
                    if local != opposite {
                        key[slot] = node;
                        slot += 1;
                    }
                }
                key.sort_unstable();
                owners.entry(key).or_insert(tetrahedron);
            }
        }
        let exterior = ClosedTetPressureSurface::new(
            &corner_nodes,
            &corner_tets,
            budget,
            maximum_face_entries,
            cx,
        )?;
        let mut exterior_edges = BTreeMap::new();
        for compact in exterior.faces() {
            checkpoint(cx)?;
            let corners = compact.map(|node| corner_ids[node]);
            for [a, b] in [[0, 1], [0, 2], [1, 2]] {
                exterior_edges.entry(edge(corners[a], corners[b])).or_insert(0);
            }
        }
        for (id, entry) in exterior_edges.values_mut().enumerate() {
            if id % 256 == 0 {
                checkpoint(cx)?;
            }
            *entry = id;
        }
        let mut triangles = Vec::with_capacity(exterior.faces().len());
        let mut material_faces = Vec::with_capacity(exterior.faces().len());
        let mut edge_orientation = vec![0_i8; exterior_edges.len()];
        for compact in exterior.faces() {
            checkpoint(cx)?;
            let mut corners = compact.map(|node| corner_ids[node]);
            let mut key = corners;
            key.sort_unstable();
            let tetrahedron = *owners
                .get(&key)
                .ok_or_else(|| invalid("quadratic exterior owner missing"))?;
            let tet = tetrahedra[tetrahedron];
            let mut local_corners = [0; 3];
            for (slot, &node) in corners.iter().enumerate() {
                local_corners[slot] = tet[..4]
                    .iter()
                    .position(|&candidate| candidate == node)
                    .ok_or_else(|| invalid("quadratic exterior corner missing"))?;
            }
            // Closed incidence is reused above; physical winding follows the
            // ORIGINAL positive TET10 chart. A curved face's normal need not
            // agree with its corner plane. Check the actual reference map here
            // and at each selected current witness, never a corner normal.
            let opposite = (0..4)
                .find(|local| !local_corners.contains(local))
                .ok_or_else(|| invalid("quadratic face has no opposite corner"))?;
            let outward = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]][opposite];
            if !(0..3).any(|start| {
                (0..3).all(|slot| local_corners[slot] == outward[(start + slot) % 3])
            }) {
                corners.swap(1, 2);
                local_corners.swap(1, 2);
            }
            let mut lambda = [1.0 / 3.0; 4];
            lambda[opposite] = 0.0;
            tet10_point(
                tet.map(|node| nodes_m[node]),
                lambda,
                budget.minimum_scaled_jacobian,
                tetrahedron,
                opposite + 1,
            )?;
            let mut triangle = [0; 6];
            triangle[..3].copy_from_slice(&corners);
            for (slot, [a, b]) in [[0, 1], [1, 2], [2, 0]].into_iter().enumerate() {
                let local_edge = edge(local_corners[a], local_corners[b]);
                let edge_index = TET10_EDGES
                    .iter()
                    .position(|&candidate| edge(candidate[0], candidate[1]) == local_edge)
                    .ok_or_else(|| invalid("quadratic exterior midside missing"))?;
                triangle[3 + slot] = tet[4 + edge_index];
            }
            let mut query_triangle = [0; 6];
            for (to, from) in query_triangle.iter_mut().zip(triangle) {
                *to = u32::try_from(from)
                    .map_err(|_| invalid("quadratic query node index exceeds u32"))?;
            }
            let mut edge_ids = [0; 3];
            for (slot, [a, b]) in [[0, 1], [0, 2], [1, 2]].into_iter().enumerate() {
                edge_ids[slot] = *exterior_edges
                    .get(&edge(corners[a], corners[b]))
                    .ok_or_else(|| invalid("quadratic exterior edge missing"))?;
            }
            for [a, b] in [[0, 1], [1, 2], [2, 0]] {
                let id = exterior_edges[&edge(corners[a], corners[b])];
                edge_orientation[id] += if corners[a] < corners[b] { 1 } else { -1 };
            }
            triangles.push(query_triangle);
            material_faces.push(MaterialFace {
                tetrahedron,
                local_corners,
                edge_ids,
            });
        }
        if edge_orientation.into_iter().any(|orientation| orientation != 0) {
            return Err(invalid("native quadratic exterior charts are inconsistently oriented"));
        }
        let feature_count = nodes_m
            .len()
            .checked_add(exterior_edges.len())
            .and_then(|count| count.checked_add(triangles.len()))
            .and_then(|count| u64::try_from(count).ok())
            .ok_or_else(|| invalid("quadratic feature ID count overflow"))?;
        let triangle_surface = PreparedQuadraticTriangleSurface::new(nodes_m, &triangles, cx)?;
        checkpoint(cx)?;
        Ok(Self {
            reference_nodes_m: nodes_m.to_vec(),
            tetrahedra: tetrahedra.to_vec(),
            triangle_surface,
            material_faces,
            edge_count: exterior_edges.len(),
            feature_count,
        })
    }

    /// Canonical outward six-node faces: corners, then edge nodes `(01,12,20)`.
    pub fn faces(&self) -> &[[u32; 6]] {
        self.triangle_surface.triangles()
    }

    /// Match every original coordinate bit and the complete ordered connectivity,
    /// including all midsides and interior nodes, rather than just the boundary.
    pub fn matches_reference(&self, nodes_m: &[[f64; 3]], tetrahedra: &[[usize; 10]]) -> bool {
        self.tetrahedra == tetrahedra
            && self.reference_nodes_m.len() == nodes_m.len()
            && self
                .reference_nodes_m
                .iter()
                .flatten()
                .zip(nodes_m.iter().flatten())
                .all(|(a, b)| a.to_bits() == b.to_bits())
    }

    /// Borrow native accepted displacement without copying constitutive history.
    pub(crate) fn current_positions_from_reference(
        &self,
        nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 10]],
        displacement_m: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<[f64; 3]>, QuadraticSolidSurfaceError> {
        checkpoint(cx)?;
        if !self.matches_reference(nodes_m, tetrahedra) {
            return Err(FiniteSolidSurfaceError::ReferenceMismatch.into());
        }
        if displacement_m.len() != nodes_m.len() {
            return Err(invalid("complete quadratic displacement required"));
        }
        let mut positions = Vec::with_capacity(nodes_m.len());
        for (index, (node, displacement)) in nodes_m.iter().zip(displacement_m).enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            vector(*displacement)?;
            positions.push(vector(core::array::from_fn(|axis| node[axis] + displacement[axis]))?);
        }
        checkpoint(cx)?;
        Ok(positions)
    }
}
