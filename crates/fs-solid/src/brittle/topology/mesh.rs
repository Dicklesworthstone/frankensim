//! Topology and transfer kernels; public construction requires admitted geometry.

use std::collections::{BTreeMap, BTreeSet};

/// An original interior triangle and its explicitly ordered receiving bodies.
/// Node order is arbitrary. `side_a_element` selects A; the other owner is B.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InteriorSeam {
    pub id: u64,
    pub nodes: [usize; 3],
    pub side_a_element: usize,
}

/// Explicit storage/work caps, checked before any mesh-sized allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FractureTopologyBudget {
    pub maximum_output_nodes: usize,
    /// Bounds the four original face entries and four corners per tetrahedron.
    pub maximum_face_entries: usize,
    pub maximum_seams: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FractureTopologyError {
    Invalid(&'static str),
    Budget(&'static str),
    Face {
        nodes: [usize; 3],
        what: &'static str,
    },
    /// A point/edge-only original connection has no unambiguous volume star.
    DisconnectedNodeStar {
        node: usize,
    },
    /// The current P1 interface owner does not admit tied crack-front vertices.
    TiedCrackFront {
        seam: u64,
        node: usize,
    },
    Conservation,
    Cancelled,
}

impl core::fmt::Display for FractureTopologyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral fracture topology refused: {self:?}")
    }
}
impl std::error::Error for FractureTopologyError {}

/// Two exterior traces in corresponding original-node order. Seam order and
/// IDs equal the caller's order; elements and A/B roles are never inferred.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SplitFace {
    pub id: u64,
    pub original_nodes: [usize; 3],
    pub elements: [usize; 2],
    pub side_a: [usize; 3],
    pub side_b: [usize; 3],
}

/// Exact represented P1 integrals, evaluated with constant element densities.
/// Angular momentum uses the current positions `reference + displacement`
/// about the reference frame's origin. Kinetic energy uses the consistent
/// tetrahedral mass, not row-sum lumping or a rigid-fragment approximation.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct P1PhysicalInventory {
    pub mass_kg: f64,
    pub first_moment_kg_m: [f64; 3],
    pub linear_momentum_kg_m_s: [f64; 3],
    pub angular_momentum_kg_m2_s: [f64; 3],
    pub kinetic_energy_j: f64,
}

/// Atomic result. Duplicated motion has the original values bit for bit.
#[derive(Debug, Clone, PartialEq)]
pub struct FractureMotionTransfer {
    pub displacement_m: Vec<f64>,
    pub velocity_m_s: Vec<f64>,
    /// Scalar consistent-mass row sums per OUTPUT node, in kg. They partition
    /// nodal mass across the new stars and must not replace the mass operator
    /// when computing kinetic energy or subsequent deformable dynamics.
    pub row_sum_mass_kg: Vec<f64>,
    pub before: P1PhysicalInventory,
    pub after: P1PhysicalInventory,
    /// Inventories in `element_fragments()` order. Reducing these by a different
    /// summation order can differ from `after` by floating-point rounding.
    pub fragments: Vec<P1PhysicalInventory>,
}

/// Immutable explicit seam topology borrowing its original admitted mesh.
/// Only node IDs change; element indices, local corner order and volume stay.
#[derive(Debug)]
pub struct TetFractureTopology<'g> {
    original_nodes_m: &'g [[f64; 3]],
    original_tetrahedra: &'g [[usize; 4]],
    volumes_m3: &'g [f64],
    nodes_m: Vec<[f64; 3]>,
    tetrahedra: Vec<[usize; 4]>,
    original_nodes: Vec<usize>,
    seams: Vec<SplitFace>,
    element_fragments: Vec<usize>,
    fragment_count: usize,
}

#[derive(Clone, Copy)]
struct Owners {
    first: usize,
    second: Option<usize>,
}

impl<'g> TetFractureTopology<'g> {
    // Only the admitted geometry facade can invoke this input boundary.
    pub(super) fn prepare<C: FnMut() -> bool>(
        nodes: &'g [[f64; 3]],
        tetrahedra: &'g [[usize; 4]],
        volumes: &'g [f64],
        seams: &[InteriorSeam],
        budget: FractureTopologyBudget,
        cancelled: &mut C,
    ) -> Result<Self, FractureTopologyError> {
        poll(cancelled)?;
        let corners = tetrahedra
            .len()
            .checked_mul(4)
            .ok_or(FractureTopologyError::Budget("tetrahedral corner overflow"))?;
        if nodes.is_empty() || tetrahedra.is_empty() || volumes.len() != tetrahedra.len() {
            return Err(FractureTopologyError::Invalid(
                "empty or inconsistent original mesh",
            ));
        }
        if corners > budget.maximum_face_entries
            || nodes.len() > budget.maximum_output_nodes
            || seams.len() > budget.maximum_seams
            || budget.maximum_output_nodes.checked_mul(3).is_none()
        {
            return Err(FractureTopologyError::Budget(
                "original mesh, output nodes or seam count",
            ));
        }
        let mut selected = BTreeMap::new();
        let mut ids = BTreeSet::new();
        for (index, seam) in seams.iter().enumerate() {
            poll(cancelled)?;
            let mut face = seam.nodes;
            face.sort_unstable();
            if face[2] >= nodes.len()
                || face[0] == face[1]
                || face[1] == face[2]
                || seam.side_a_element >= tetrahedra.len()
            {
                return Err(FractureTopologyError::Invalid(
                    "invalid selected face or A-side element",
                ));
            }
            if selected.insert(face, index).is_some()
                || !ids.insert(seam.id)
                || seam
                    .id
                    .checked_mul(3)
                    .and_then(|id| id.checked_add(2))
                    .is_none()
            {
                return Err(FractureTopologyError::Invalid(
                    "duplicate face, duplicate ID or overflowing patch ID",
                ));
            }
        }
        // All original faces are indexed once. Original geometry already owns
        // orientation, quality, volume and manifold-face admission.
        let mut faces = BTreeMap::<[usize; 3], Owners>::new();
        for (element, tet) in tetrahedra.iter().enumerate() {
            poll(cancelled)?;
            for opposite in 0..4 {
                let key = face_without(*tet, opposite);
                match faces.entry(key) {
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        entry.insert(Owners {
                            first: element,
                            second: None,
                        });
                    }
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        if entry.get().second.is_some() {
                            return Err(FractureTopologyError::Face {
                                nodes: key,
                                what: "more than two incident elements",
                            });
                        }
                        entry.get_mut().second = Some(element);
                    }
                }
            }
        }
        for (face, &index) in &selected {
            poll(cancelled)?;
            let owners = faces.get(face).ok_or(FractureTopologyError::Face {
                nodes: *face,
                what: "selected triangle is not an original face",
            })?;
            if owners.second.is_none() {
                return Err(FractureTopologyError::Face {
                    nodes: *face,
                    what: "selected face is already exterior",
                });
            }
            let side = seams[index].side_a_element;
            if owners.first != side && owners.second != Some(side) {
                return Err(FractureTopologyError::Face {
                    nodes: *face,
                    what: "selected A-side does not own the face",
                });
            }
        }
        // A corner is a (tetrahedron, original node) pair. Union only equal
        // original nodes across surviving faces: no cut propagates elsewhere.
        let mut original_stars = UnionFind::new(corners);
        let mut split_stars = UnionFind::new(corners);
        let mut components = UnionFind::new(tetrahedra.len());
        for (face, owners) in &faces {
            poll(cancelled)?;
            let Some(b) = owners.second else { continue };
            let cut = selected.contains_key(face);
            if !cut {
                components.join(owners.first, b);
            }
            for &node in face {
                let a_corner = corner(tetrahedra, owners.first, node);
                let b_corner = corner(tetrahedra, b, node);
                original_stars.join(a_corner, b_corner);
                if !cut {
                    split_stars.join(a_corner, b_corner);
                }
            }
        }
        let mut original_root = vec![None; nodes.len()];
        let mut assigned = vec![false; nodes.len()];
        let mut output_of_root = vec![usize::MAX; corners];
        let mut output_nodes = nodes.to_vec();
        let mut original_nodes: Vec<usize> = (0..nodes.len()).collect();
        let mut output_tets = Vec::with_capacity(tetrahedra.len());
        for (element, tet) in tetrahedra.iter().enumerate() {
            poll(cancelled)?;
            let mut output = [0; 4];
            for (local, &node) in tet.iter().enumerate() {
                let index = 4 * element + local;
                let original = original_stars.find(index);
                match original_root[node] {
                    Some(previous) if previous != original => {
                        return Err(FractureTopologyError::DisconnectedNodeStar { node });
                    }
                    None => original_root[node] = Some(original),
                    _ => {}
                }
                let root = split_stars.find(index);
                if output_of_root[root] == usize::MAX {
                    let target = if assigned[node] {
                        if output_nodes.len() == budget.maximum_output_nodes {
                            return Err(FractureTopologyError::Budget("duplicated node count"));
                        }
                        let target = output_nodes.len();
                        output_nodes.push(nodes[node]);
                        original_nodes.push(node);
                        target
                    } else {
                        assigned[node] = true;
                        node
                    };
                    output_of_root[root] = target;
                }
                output[local] = output_of_root[root];
            }
            output_tets.push(output);
        }
        let mut split_faces = Vec::with_capacity(seams.len());
        for seam in seams {
            poll(cancelled)?;
            let mut original_nodes = seam.nodes;
            original_nodes.sort_unstable();
            let owners = faces[&original_nodes];
            let a = seam.side_a_element;
            let b = if owners.first == a {
                owners.second.unwrap()
            } else {
                owners.first
            };
            let side_a = original_nodes.map(|node| output_tets[a][corner(tetrahedra, a, node) % 4]);
            let side_b = original_nodes.map(|node| output_tets[b][corner(tetrahedra, b, node) % 4]);
            for i in 0..3 {
                if side_a[i] == side_b[i] {
                    return Err(FractureTopologyError::TiedCrackFront {
                        seam: seam.id,
                        node: original_nodes[i],
                    });
                }
            }
            split_faces.push(SplitFace {
                id: seam.id,
                original_nodes,
                elements: [a, b],
                side_a,
                side_b,
            });
        }
        let mut fragment_of_root = vec![usize::MAX; tetrahedra.len()];
        let mut element_fragments = Vec::with_capacity(tetrahedra.len());
        let mut fragment_count = 0;
        for element in 0..tetrahedra.len() {
            poll(cancelled)?;
            let root = components.find(element);
            if fragment_of_root[root] == usize::MAX {
                fragment_of_root[root] = fragment_count;
                fragment_count += 1;
            }
            element_fragments.push(fragment_of_root[root]);
        }
        poll(cancelled)?;
        Ok(Self {
            original_nodes_m: nodes,
            original_tetrahedra: tetrahedra,
            volumes_m3: volumes,
            nodes_m: output_nodes,
            tetrahedra: output_tets,
            original_nodes,
            seams: split_faces,
            element_fragments,
            fragment_count,
        })
    }

    pub fn nodes_m(&self) -> &[[f64; 3]] {
        &self.nodes_m
    }
    pub fn tetrahedra(&self) -> &[[usize; 4]] {
        &self.tetrahedra
    }
    /// Output-to-input node map. Original IDs are retained for the first star;
    /// additional copies follow deterministic original element/corner order.
    pub fn original_nodes(&self) -> &[usize] {
        &self.original_nodes
    }
    pub fn seams(&self) -> &[SplitFace] {
        &self.seams
    }
    /// Components of the volume mesh with selected faces removed. IDs are dense
    /// in first-original-element order, independent of the supplied seam order.
    pub fn element_fragments(&self) -> &[usize] {
        &self.element_fragments
    }
    pub const fn fragment_count(&self) -> usize {
        self.fragment_count
    }

    /// Material/source bindings and one-state-per-element histories retain their
    /// actual ownership and exact order. No law is re-created, reset or averaged.
    /// Multiple quadrature states per element keep their existing block layout
    /// (for example, pass a slice of fixed-size per-element state arrays).
    pub fn retain_element_values<'s, S>(
        &self,
        values: &'s [S],
    ) -> Result<&'s [S], FractureTopologyError> {
        if values.len() != self.tetrahedra.len() {
            return Err(FractureTopologyError::Invalid(
                "element state/binding count",
            ));
        }
        Ok(values)
    }

    /// Predicate form of the same atomic motion/inertia transfer. Element
    /// integration and local node order are unchanged, so both independently
    /// reconstructed global inventories must compare exactly. No tolerance can
    /// hide lost mass, artificial fracture impulses or added kinetic energy.
    pub fn transfer_motion_controlled<C: FnMut() -> bool>(
        &self,
        displacement_m: &[f64],
        velocity_m_s: &[f64],
        density_kg_m3: &[f64],
        cancelled: &mut C,
    ) -> Result<FractureMotionTransfer, FractureTopologyError> {
        poll(cancelled)?;
        if displacement_m.len() != 3 * self.original_nodes_m.len()
            || velocity_m_s.len() != displacement_m.len()
            || density_kg_m3.len() != self.tetrahedra.len()
        {
            return Err(FractureTopologyError::Invalid(
                "motion or element-density dimensions",
            ));
        }
        for values in [displacement_m, velocity_m_s, density_kg_m3] {
            for chunk in values.chunks(256) {
                poll(cancelled)?;
                for &value in chunk {
                    finite(value)?;
                }
            }
        }
        if density_kg_m3.iter().any(|&density| density <= 0.0) {
            return Err(FractureTopologyError::Invalid(
                "positive reference mass density required",
            ));
        }
        let before = inventory(
            self.original_nodes_m,
            self.original_tetrahedra,
            self.volumes_m3,
            density_kg_m3,
            displacement_m,
            velocity_m_s,
            None,
            cancelled,
        )?
        .0;
        let mut displacement = Vec::with_capacity(3 * self.nodes_m.len());
        let mut velocity = Vec::with_capacity(3 * self.nodes_m.len());
        for (index, &node) in self.original_nodes.iter().enumerate() {
            if index % 256 == 0 {
                poll(cancelled)?;
            }
            displacement.extend_from_slice(&displacement_m[3 * node..3 * node + 3]);
            velocity.extend_from_slice(&velocity_m_s[3 * node..3 * node + 3]);
        }
        let (after, fragments, row_sum_mass_kg) = inventory(
            &self.nodes_m,
            &self.tetrahedra,
            self.volumes_m3,
            density_kg_m3,
            &displacement,
            &velocity,
            Some((&self.element_fragments, self.fragment_count)),
            cancelled,
        )?;
        if before != after {
            return Err(FractureTopologyError::Conservation);
        }
        poll(cancelled)?;
        Ok(FractureMotionTransfer {
            displacement_m: displacement,
            velocity_m_s: velocity,
            row_sum_mass_kg,
            before,
            after,
            fragments,
        })
    }
}

fn inventory<C: FnMut() -> bool>(
    nodes: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
    volumes: &[f64],
    densities: &[f64],
    displacement: &[f64],
    velocity: &[f64],
    fragments: Option<(&[usize], usize)>,
    cancelled: &mut C,
) -> Result<(P1PhysicalInventory, Vec<P1PhysicalInventory>, Vec<f64>), FractureTopologyError> {
    let mut total = InventorySum::default();
    let mut pieces = vec![InventorySum::default(); fragments.map_or(0, |(_, count)| count)];
    let mut masses = if fragments.is_some() {
        vec![0.0; nodes.len()]
    } else {
        Vec::new()
    };
    for (element, tet) in tetrahedra.iter().enumerate() {
        poll(cancelled)?;
        let mass = finite(densities[element] * volumes[element])?;
        let scale = mass / 20.0;
        if mass <= 0.0 || scale <= 0.0 {
            return Err(FractureTopologyError::Invalid(
                "unrepresentable positive consistent P1 mass",
            ));
        }
        let mut x = [[0.0; 3]; 4];
        let mut v = [[0.0; 3]; 4];
        let mut sum_v = [0.0; 3];
        for (local, &node) in tet.iter().enumerate() {
            for c in 0..3 {
                x[local][c] = finite(nodes[node][c] + displacement[3 * node + c])?;
                v[local][c] = velocity[3 * node + c];
                sum_v[c] = finite(sum_v[c] + v[local][c])?;
            }
            if !masses.is_empty() {
                masses[node] = finite(masses[node] + mass / 4.0)?;
            }
        }
        let mut one = InventorySum::default();
        one.values[0].add(mass)?;
        // Integral N_i N_j = V (1 + delta_ij) / 20. The momentum at
        // corner i is m/20 * (sum_j v_j + v_i), including all cross terms.
        for i in 0..4 {
            let mut p = [0.0; 3];
            for c in 0..3 {
                p[c] = finite(scale * finite(sum_v[c] + v[i][c])?)?;
                one.values[1 + c].add(finite((mass / 4.0) * x[i][c])?)?;
                one.values[4 + c].add(p[c])?;
                // Positive quadratic form avoids cancellation between signed
                // nodal work terms. Multiplying the small mass weight first
                // avoids unnecessarily squaring a large velocity.
                one.values[10].add(finite(v[i][c] * ((0.5 * scale) * v[i][c]))?)?;
            }
            let moment = cross(x[i], p);
            for c in 0..3 {
                one.values[7 + c].add(finite(moment[c])?)?;
            }
        }
        for c in 0..3 {
            one.values[10].add(finite(sum_v[c] * ((0.5 * scale) * sum_v[c]))?)?;
        }
        total.add(one)?;
        if let Some((ids, _)) = fragments {
            pieces[ids[element]].add(one)?;
        }
    }
    poll(cancelled)?;
    Ok((
        total.finish(),
        pieces.into_iter().map(InventorySum::finish).collect(),
        masses,
    ))
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), FractureTopologyError> {
        let adjusted = finite(value - self.correction)?;
        let next = finite(self.value + adjusted)?;
        self.correction = finite((next - self.value) - adjusted)?;
        self.value = next;
        Ok(())
    }
}
#[derive(Clone, Copy, Default)]
struct InventorySum {
    values: [Sum; 11],
}
impl InventorySum {
    fn add(&mut self, other: Self) -> Result<(), FractureTopologyError> {
        for (target, value) in self.values.iter_mut().zip(other.values) {
            target.add(value.value)?;
        }
        Ok(())
    }
    fn finish(self) -> P1PhysicalInventory {
        let v = self.values.map(|sum| sum.value);
        P1PhysicalInventory {
            mass_kg: v[0],
            first_moment_kg_m: [v[1], v[2], v[3]],
            linear_momentum_kg_m_s: [v[4], v[5], v[6]],
            angular_momentum_kg_m2_s: [v[7], v[8], v[9]],
            kinetic_energy_j: v[10],
        }
    }
}

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}
impl UnionFind {
    fn new(count: usize) -> Self {
        Self {
            parent: (0..count).collect(),
            rank: vec![0; count],
        }
    }
    fn find(&mut self, mut value: usize) -> usize {
        while self.parent[value] != value {
            self.parent[value] = self.parent[self.parent[value]];
            value = self.parent[value];
        }
        value
    }
    fn join(&mut self, a: usize, b: usize) {
        let mut a = self.find(a);
        let mut b = self.find(b);
        if a == b {
            return;
        }
        if self.rank[a] < self.rank[b] || (self.rank[a] == self.rank[b] && a > b) {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[b] = a;
        if self.rank[a] == self.rank[b] {
            self.rank[a] += 1;
        }
    }
}
fn corner(tetrahedra: &[[usize; 4]], element: usize, node: usize) -> usize {
    // The face index contains only nodes of each incident original tetrahedron.
    4 * element
        + tetrahedra[element]
            .iter()
            .position(|&id| id == node)
            .unwrap()
}
fn face_without(tet: [usize; 4], opposite: usize) -> [usize; 3] {
    let mut face = [0; 3];
    let mut index = 0;
    for (local, node) in tet.into_iter().enumerate() {
        if local != opposite {
            face[index] = node;
            index += 1;
        }
    }
    face.sort_unstable();
    face
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn finite(value: f64) -> Result<f64, FractureTopologyError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FractureTopologyError::Invalid(
            "nonfinite physical transfer arithmetic",
        ))
    }
}
fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), FractureTopologyError> {
    if cancelled() {
        Err(FractureTopologyError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const NODES: [[f64; 3]; 5] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0],
    ];
    const TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [0, 2, 1, 4]];
    const VOLUMES: [f64; 2] = [1.0 / 6.0; 2];
    const SEAM: InteriorSeam = InteriorSeam {
        id: 41,
        nodes: [2, 0, 1],
        side_a_element: 0,
    };
    fn budget() -> FractureTopologyBudget {
        FractureTopologyBudget {
            maximum_output_nodes: 64,
            maximum_face_entries: 64,
            maximum_seams: 16,
        }
    }
    fn split() -> TetFractureTopology<'static> {
        TetFractureTopology::prepare(&NODES, &TETS, &VOLUMES, &[SEAM], budget(), &mut || false)
            .unwrap()
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3.0e-13 * actual.abs().max(expected.abs()).max(1.0),
            "actual={actual:e}, expected={expected:e}"
        );
    }

    // Independent positive four-point degree-two tetrahedral quadrature.
    // This does not reuse the production consistent-mass moment identities.
    fn quadrature(
        nodes: &[[f64; 3]],
        tets: &[[usize; 4]],
        volumes: &[f64],
        density: &[f64],
        u: &[f64],
        v: &[f64],
    ) -> P1PhysicalInventory {
        let a = 0.5854101966249685;
        let b = 0.1381966011250105;
        let mut result = P1PhysicalInventory::default();
        for (element, tet) in tets.iter().enumerate() {
            for q in 0..4 {
                let mut x = [0.0; 3];
                let mut velocity = [0.0; 3];
                for i in 0..4 {
                    let shape = if i == q { a } else { b };
                    for c in 0..3 {
                        x[c] += shape * (nodes[tet[i]][c] + u[3 * tet[i] + c]);
                        velocity[c] += shape * v[3 * tet[i] + c];
                    }
                }
                let mass = density[element] * volumes[element] / 4.0;
                result.mass_kg += mass;
                let moment = cross(x, velocity);
                for c in 0..3 {
                    result.first_moment_kg_m[c] += mass * x[c];
                    result.linear_momentum_kg_m_s[c] += mass * velocity[c];
                    result.angular_momentum_kg_m2_s[c] += mass * moment[c];
                    result.kinetic_energy_j += 0.5 * mass * velocity[c] * velocity[c];
                }
            }
        }
        result
    }
    fn compare(a: P1PhysicalInventory, b: P1PhysicalInventory) {
        close(a.mass_kg, b.mass_kg);
        close(a.kinetic_energy_j, b.kinetic_energy_j);
        for c in 0..3 {
            close(a.first_moment_kg_m[c], b.first_moment_kg_m[c]);
            close(a.linear_momentum_kg_m_s[c], b.linear_momentum_kg_m_s[c]);
            close(a.angular_momentum_kg_m2_s[c], b.angular_momentum_kg_m2_s[c]);
        }
    }

    #[test]
    fn through_seam_duplicates_only_needed_stars_and_preserves_original_element_state() {
        let topology = split();
        assert_eq!(topology.nodes_m().len(), 8);
        assert_eq!(topology.tetrahedra()[0], TETS[0]);
        assert_eq!(topology.element_fragments(), &[0, 1]);
        assert_eq!(topology.fragment_count(), 2);
        assert_eq!(topology.seams()[0].id, SEAM.id);
        assert_eq!(topology.seams()[0].elements, [0, 1]);
        assert_eq!(topology.seams()[0].original_nodes, [0, 1, 2]);
        for (element, tet) in topology.tetrahedra().iter().enumerate() {
            for i in 0..4 {
                assert_eq!(topology.original_nodes()[tet[i]], TETS[element][i]);
                assert_eq!(topology.nodes_m()[tet[i]], NODES[TETS[element][i]]);
            }
        }
        for i in 0..3 {
            let face = topology.seams()[0];
            assert_ne!(face.side_a[i], face.side_b[i]);
            assert_eq!(
                topology.nodes_m()[face.side_a[i]],
                topology.nodes_m()[face.side_b[i]]
            );
        }
        // Histories need not be Clone and are not reset or rebound to new laws.
        struct History {
            material: &'static str,
            dissipated_j: f64,
        }
        let history = [
            History {
                material: "A",
                dissipated_j: 7.0,
            },
            History {
                material: "B",
                dissipated_j: 11.0,
            },
        ];
        let retained = topology.retain_element_values(&history).unwrap();
        assert!(std::ptr::eq(retained.as_ptr(), history.as_ptr()));
        assert_eq!(retained[0].material, "A");
        assert_eq!(retained[1].dissipated_j, 11.0);
        assert!(topology.retain_element_values(&history[..1]).is_err());
    }

    #[test]
    fn split_motion_preserves_consistent_mass_momentum_and_energy_against_independent_quadrature() {
        let topology = split();
        let u: Vec<_> = (0..15)
            .map(|i| 0.001 * ((i * 7 % 13) as f64 - 6.0))
            .collect();
        let v: Vec<_> = (0..15).map(|i| 0.2 * ((i * 3 % 11) as f64 - 5.0)).collect();
        let density = [6.0, 12.0];
        let transfer = topology
            .transfer_motion_controlled(&u, &v, &density, &mut || false)
            .unwrap();
        assert_eq!(transfer.before, transfer.after);
        compare(
            transfer.before,
            quadrature(&NODES, &TETS, &VOLUMES, &density, &u, &v),
        );
        compare(
            transfer.after,
            quadrature(
                topology.nodes_m(),
                topology.tetrahedra(),
                &VOLUMES,
                &density,
                &transfer.displacement_m,
                &transfer.velocity_m_s,
            ),
        );
        close(transfer.after.mass_kg, 3.0);
        close(transfer.row_sum_mass_kg.iter().sum(), 3.0);
        assert_eq!(transfer.fragments.len(), 2);
        close(transfer.fragments[0].mass_kg, 1.0);
        close(transfer.fragments[1].mass_kg, 2.0);
        for (node, &original) in topology.original_nodes().iter().enumerate() {
            for c in 0..3 {
                assert_eq!(
                    transfer.displacement_m[3 * node + c].to_bits(),
                    u[3 * original + c].to_bits()
                );
                assert_eq!(
                    transfer.velocity_m_s[3 * node + c].to_bits(),
                    v[3 * original + c].to_bits()
                );
            }
        }
        let mut one_corner = vec![0.0; 15];
        one_corner[9] = 2.0;
        let nonrigid = topology
            .transfer_motion_controlled(&[0.0; 15], &one_corner, &density, &mut || false)
            .unwrap();
        // Only one corner of the unit-mass first tet moves. The consistent
        // kinetic energy is m*v^2/20 = 0.2 J; row-sum lumping would give 0.5 J.
        close(nonrigid.after.kinetic_energy_j, 0.2);
        assert!((nonrigid.after.kinetic_energy_j - 0.5).abs() > 0.2);
    }

    #[test]
    fn inventories_transform_as_physical_vectors_and_moments() {
        let rotate = |v: [f64; 3]| [-v[1], v[0], v[2]];
        let shift = [4.0, -2.0, 7.0];
        let rotated_nodes = NODES.map(|x| {
            let r = rotate(x);
            std::array::from_fn(|c| r[c] + shift[c])
        });
        let original = split();
        let transformed = TetFractureTopology::prepare(
            &rotated_nodes,
            &TETS,
            &VOLUMES,
            &[SEAM],
            budget(),
            &mut || false,
        )
        .unwrap();
        let u: Vec<_> = (0..15).map(|i| 0.003 * (i as f64 - 7.0)).collect();
        let v: Vec<_> = (0..15).map(|i| 0.1 * ((i * 7 % 17) as f64 - 8.0)).collect();
        let rotate_field = |field: &[f64]| {
            field
                .chunks_exact(3)
                .flat_map(|a| rotate([a[0], a[1], a[2]]))
                .collect::<Vec<_>>()
        };
        let a = original
            .transfer_motion_controlled(&u, &v, &[6.0, 12.0], &mut || false)
            .unwrap()
            .after;
        let b = transformed
            .transfer_motion_controlled(
                &rotate_field(&u),
                &rotate_field(&v),
                &[6.0, 12.0],
                &mut || false,
            )
            .unwrap()
            .after;
        close(a.mass_kg, b.mass_kg);
        close(a.kinetic_energy_j, b.kinetic_energy_j);
        let p = rotate(a.linear_momentum_kg_m_s);
        let first = rotate(a.first_moment_kg_m);
        let angular = rotate(a.angular_momentum_kg_m2_s);
        let moved_origin = cross(shift, p);
        for c in 0..3 {
            close(b.linear_momentum_kg_m_s[c], p[c]);
            close(b.first_moment_kg_m[c], first[c] + shift[c] * a.mass_kg);
            close(b.angular_momentum_kg_m2_s[c], angular[c] + moved_origin[c]);
        }
    }

    #[test]
    fn unselected_faces_stay_shared_and_seam_order_never_selects_extra_cuts() {
        let nodes = [
            NODES[0],
            NODES[1],
            NODES[2],
            NODES[3],
            NODES[4],
            [0.0, -1.0, 0.0],
        ];
        let tets = [TETS[0], TETS[1], [0, 1, 4, 5]];
        let volumes = [1.0 / 6.0; 3];
        let topology =
            TetFractureTopology::prepare(&nodes, &tets, &volumes, &[SEAM], budget(), &mut || false)
                .unwrap();
        assert_eq!(topology.element_fragments(), &[0, 1, 1]);
        for node in [0, 1, 4] {
            let left = topology.tetrahedra()[1][corner(&tets, 1, node) % 4];
            let right = topology.tetrahedra()[2][corner(&tets, 2, node) % 4];
            assert_eq!(left, right, "uncut face node {node} was split");
        }
        let second = InteriorSeam {
            id: 99,
            nodes: [4, 1, 0],
            side_a_element: 2,
        };
        let a = TetFractureTopology::prepare(
            &nodes,
            &tets,
            &volumes,
            &[SEAM, second],
            budget(),
            &mut || false,
        )
        .unwrap();
        let b = TetFractureTopology::prepare(
            &nodes,
            &tets,
            &volumes,
            &[second, SEAM],
            budget(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(a.tetrahedra(), b.tetrahedra());
        assert_eq!(a.original_nodes(), b.original_nodes());
        assert_eq!(a.element_fragments(), &[0, 1, 2]);
        assert_eq!(b.seams()[0].id, 99);
        assert_eq!(b.seams()[0].elements, [2, 1]);
        assert_eq!(a.seams()[0], b.seams()[1]);
    }

    #[test]
    fn tied_crack_front_and_original_point_connections_refuse_instead_of_adding_cuts() {
        let nodes = [
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [-0.5, 0.9, 0.0],
            [-0.5, -0.9, 0.0],
        ];
        let tets = [[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 4, 2]];
        let volumes = [0.3; 3];
        let seam = InteriorSeam {
            id: 3,
            nodes: [0, 1, 3],
            side_a_element: 0,
        };
        assert!(matches!(
            TetFractureTopology::prepare(&nodes, &tets, &volumes, &[seam], budget(), &mut || false),
            Err(FractureTopologyError::TiedCrackFront { seam: 3, node: 0 })
        ));
        let all = [
            seam,
            InteriorSeam {
                id: 4,
                nodes: [0, 1, 4],
                side_a_element: 1,
            },
            InteriorSeam {
                id: 5,
                nodes: [0, 1, 2],
                side_a_element: 2,
            },
        ];
        let separated =
            TetFractureTopology::prepare(&nodes, &tets, &volumes, &all, budget(), &mut || false)
                .unwrap();
        assert_eq!(separated.element_fragments(), &[0, 1, 2]);
        let touching_nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
        ];
        let touching_tets = [[0, 1, 2, 3], [0, 4, 5, 6]];
        assert!(matches!(
            TetFractureTopology::prepare(
                &touching_nodes,
                &touching_tets,
                &VOLUMES,
                &[],
                budget(),
                &mut || false
            ),
            Err(FractureTopologyError::DisconnectedNodeStar { node: 0 })
        ));
    }

    #[test]
    fn admission_caps_and_cancellation_leave_original_state_reusable() {
        let mut cap = budget();
        cap.maximum_output_nodes = 8;
        assert!(
            TetFractureTopology::prepare(&NODES, &TETS, &VOLUMES, &[SEAM], cap, &mut || false)
                .is_ok()
        );
        cap.maximum_output_nodes = 7;
        assert!(matches!(
            TetFractureTopology::prepare(&NODES, &TETS, &VOLUMES, &[SEAM], cap, &mut || false),
            Err(FractureTopologyError::Budget(_))
        ));
        for bad in [
            InteriorSeam {
                nodes: [0, 1, 3],
                ..SEAM
            },
            InteriorSeam {
                nodes: [0, 3, 4],
                ..SEAM
            },
            InteriorSeam {
                id: u64::MAX,
                ..SEAM
            },
        ] {
            assert!(
                TetFractureTopology::prepare(
                    &NODES,
                    &TETS,
                    &VOLUMES,
                    &[bad],
                    budget(),
                    &mut || false
                )
                .is_err()
            );
        }
        assert!(
            TetFractureTopology::prepare(
                &NODES,
                &TETS,
                &VOLUMES,
                &[SEAM, SEAM],
                budget(),
                &mut || false
            )
            .is_err()
        );
        for stop in [1, 4, 10, 16] {
            let mut polls = 0;
            assert!(matches!(
                TetFractureTopology::prepare(
                    &NODES,
                    &TETS,
                    &VOLUMES,
                    &[SEAM],
                    budget(),
                    &mut || {
                        polls += 1;
                        polls == stop
                    }
                ),
                Err(FractureTopologyError::Cancelled)
            ));
        }
        let topology = split();
        let u = [0.0; 15];
        let v = [1.0; 15];
        let before = topology
            .transfer_motion_controlled(&u, &v, &[6.0, 12.0], &mut || false)
            .unwrap();
        for stop in [1, 3, 6, 9] {
            let mut polls = 0;
            assert!(matches!(
                topology.transfer_motion_controlled(&u, &v, &[6.0, 12.0], &mut || {
                    polls += 1;
                    polls == stop
                }),
                Err(FractureTopologyError::Cancelled)
            ));
        }
        let mut invalid = u;
        invalid[0] = f64::NAN;
        assert!(
            topology
                .transfer_motion_controlled(&invalid, &v, &[6.0, 12.0], &mut || false)
                .is_err()
        );
        assert!(
            topology
                .transfer_motion_controlled(&u, &v, &[6.0, 0.0], &mut || false)
                .is_err()
        );
        assert!(
            topology
                .transfer_motion_controlled(&u[..12], &v, &[6.0, 12.0], &mut || false)
                .is_err()
        );
        assert_eq!(
            before,
            topology
                .transfer_motion_controlled(&u, &v, &[6.0, 12.0], &mut || false)
                .unwrap()
        );
    }
}
