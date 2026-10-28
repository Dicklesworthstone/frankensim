//! Matching P1 tetrahedral traces lowered to the existing cohesive operator.
//!
//! Three positive area-weighted points integrate an undamaged linear interface
//! stiffness exactly. Nonlinear softening is quadrature-discrete and requires
//! independent surface refinement. Frames stay fixed in the reference mesh.
//! Pairing is explicit and coordinates must compare equal: no nearest
//! neighbour bond discovery, projection, snapping, or remeshing is performed.

use std::collections::BTreeMap;
use super::{BilinearCohesive, CohesivePatch, Cx, Error, JumpTerm, dot, finite};
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_matdb::MaterialStateId;

/// Work-conserving overlap quadrature for explicitly bonded nonmatching traces.
pub mod nonmatching;

/// Explicitly bonded exterior faces. Node order is arbitrary on either side.
/// The two traces must use distinct DOFs even though their coordinates coincide.
#[derive(Debug, Clone, Copy)]
pub struct FacePair<'a> {
    pub id: u64,
    pub side_a: [usize; 3],
    pub side_b: [usize; 3],
    pub law: &'a BilinearCohesive,
}

/// An explicit bonded face pair with the original ordered interface material.
/// Surface material IDs describe the actual assignments on A and B; their order
/// follows the geometric jump, independently of vertex permutation on a face.
#[derive(Debug, Clone, Copy)]
pub struct MaterialFacePair<'a> {
    pub id: u64,
    pub side_a: [usize; 3],
    pub side_b: [usize; 3],
    pub source: &'a ResolvedBilinearCohesive,
    pub surface_materials: [&'a MaterialStateId; 2],
}

/// Bound allocations and admit a minimum dimensionless doubled-area/edge^2.
#[derive(Debug, Clone, Copy)]
pub struct TraceBudget {
    pub maximum_nodes: usize,
    pub maximum_tetrahedra: usize,
    pub maximum_interfaces: usize,
    pub minimum_face_quality: f64,
}

/// Quadrature location and fixed work-conjugate local frame.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfacePoint {
    pub interface_id: u64,
    /// Quadrature index within the face, or within the bonded region for
    /// nonmatching traces with overlap subdivision.
    pub point_in_face: usize,
    pub location_m: [f64; 3],
    pub area_weight_m2: f64,
    /// Unit axes as ROWS: outward normal from A, first tangent, second tangent.
    /// Positive normal jump is separation of B from A, not penetration.
    pub axes: [[f64; 3]; 3],
    pub shape_values: [f64; 3],
    /// Shape values on B; equal to `shape_values` for matching traces.
    pub side_b_shape_values: [f64; 3],
    /// Nodes on A and B in their respective shape-value order. For matching
    /// traces, B is permuted to the exact positional order of A.
    pub paired_nodes: [[usize; 3]; 2],
}

/// One immutable topology/map; reduced and prescribed motions use the SAME map.
#[derive(Debug, Clone)]
pub struct P1Interfaces<'a> {
    full_dimension: usize,
    free_dofs: Vec<usize>,
    reduced_of: Vec<Option<usize>>,
    full_patches: Vec<CohesivePatch<'a>>,
    patches: Vec<CohesivePatch<'a>>,
    points: Vec<SurfacePoint>,
}

impl<'a> P1Interfaces<'a> {
    /// Lower source-bound bonds through the same geometric pairing, quadrature
    /// and work-conjugate jump maps as [`Self::new`]. No second assembly exists.
    pub fn from_interfaces(nodes_m: &[[f64; 3]], tetrahedra: &[[usize; 4]],
        free_dofs: &[usize], pairs: &[MaterialFacePair<'a>], budget: TraceBudget,
        cx: &Cx<'_>) -> Result<Self, Error> {
        Self::from_interfaces_controlled(nodes_m, tetrahedra, free_dofs, pairs, budget,
            &mut || cx.checkpoint().is_err())
    }

    /// Predicate form with identical source, geometry and cancellation checks.
    /// Every full/reduced quadrature patch retains the original source bundle.
    pub fn from_interfaces_controlled<C: FnMut() -> bool>(nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 4]], free_dofs: &[usize], pairs: &[MaterialFacePair<'a>],
        budget: TraceBudget, cancelled: &mut C) -> Result<Self, Error> {
        poll(cancelled)?;
        if pairs.is_empty() || pairs.len() > budget.maximum_interfaces {
            return Err(Error::Budget("P1 source interface count"));
        }
        let mut geometric = Vec::with_capacity(pairs.len());
        for pair in pairs {
            poll(cancelled)?;
            pair.source.require_surface_materials(pair.surface_materials[0], pair.surface_materials[1])
                .map_err(|error| Error::Source { patch: pair.id, error })?;
            geometric.push(FacePair { id: pair.id, side_a: pair.side_a,
                side_b: pair.side_b, law: pair.source.law() });
        }
        let mut result = Self::new_controlled(nodes_m, tetrahedra, free_dofs,
            &geometric, budget, cancelled)?;
        for (index, pair) in pairs.iter().enumerate() {
            poll(cancelled)?;
            for patch in &mut result.patches[3*index..3*index+3] {
                patch.source = Some(pair.source);
            }
            for patch in &mut result.full_patches[3*index..3*index+3] {
                patch.source = Some(pair.source);
            }
        }
        poll(cancelled)?;
        Ok(result)
    }

    /// Validate that each requested face has exactly one incident tetrahedron,
    /// each partner lies on the opposite side, and each face is used once.
    /// Free DOFs must be strictly increasing full node-major xyz indices.
    /// Construction polls per tetrahedron/interface and every 256 nodes/DOFs.
    /// The volume owner still admits element quality and material properties.
    pub fn new(nodes_m: &[[f64; 3]], tetrahedra: &[[usize; 4]], free_dofs: &[usize],
        pairs: &[FacePair<'a>], budget: TraceBudget, cx: &Cx<'_>) -> Result<Self, Error> {
        Self::new_controlled(nodes_m, tetrahedra, free_dofs, pairs, budget,
            &mut || cx.checkpoint().is_err())
    }

    pub fn new_controlled<C: FnMut() -> bool>(nodes_m: &[[f64; 3]], tetrahedra: &[[usize; 4]],
        free_dofs: &[usize], pairs: &[FacePair<'a>], budget: TraceBudget,
        cancelled: &mut C) -> Result<Self, Error> {
        poll(cancelled)?;
        positive(budget.minimum_face_quality)?;
        if budget.minimum_face_quality >= 1.0 { return Err(Error::Invalid("face quality must be below one")); }
        if nodes_m.is_empty() || tetrahedra.is_empty() || pairs.is_empty()
            || nodes_m.len() > budget.maximum_nodes || tetrahedra.len() > budget.maximum_tetrahedra
            || pairs.len() > budget.maximum_interfaces { return Err(Error::Budget("P1 mesh or interface count")); }
        let full_dimension = nodes_m.len().checked_mul(3).ok_or(Error::Budget("P1 DOF overflow"))?;
        let count = pairs.len().checked_mul(3).ok_or(Error::Budget("P1 quadrature overflow"))?;
        count.checked_mul(18).ok_or(Error::Budget("P1 jump map overflow"))?;
        for chunk in nodes_m.chunks(256) { poll(cancelled)?; for node in chunk { for &x in node { finite(x)?; } } }
        if free_dofs.len() > full_dimension { return Err(Error::Invalid("too many free P1 DOFs")); }
        let mut reduced_of = vec![None; full_dimension];
        for (index, &dof) in free_dofs.iter().enumerate() {
            if index % 256 == 0 { poll(cancelled)?; }
            if dof >= full_dimension || (index > 0 && dof <= free_dofs[index-1]) {
                return Err(Error::Invalid("P1 free DOFs must be distinct, sorted and in range"));
            }
            reduced_of[dof] = Some(index);
        }
        // Index only requested faces: O(interface count) storage, not all 4E faces.
        let mut faces = BTreeMap::<[usize; 3], (usize, usize)>::new();
        let mut ids = std::collections::BTreeSet::new();
        for pair in pairs {
            poll(cancelled)?;
            if !ids.insert(pair.id) || pair.id.checked_mul(3).and_then(|n| n.checked_add(2)).is_none() {
                return Err(Error::Invalid("duplicate or overflowing interface identity"));
            }
            let mut six = Vec::with_capacity(6);
            for mut face in [pair.side_a, pair.side_b] {
                face.sort_unstable();
                if face[2] >= nodes_m.len() || face[0] == face[1] || face[1] == face[2]
                    || faces.insert(face, (0, 0)).is_some() { return Err(Error::Invalid("invalid or reused interface face")); }
                six.extend(face);
            }
            six.sort_unstable();
            if six.windows(2).any(|w| w[0] == w[1]) { return Err(Error::Invalid("cohesive traces must have distinct vertices")); }
        }
        for tet in tetrahedra {
            poll(cancelled)?;
            let mut sorted = *tet;
            sorted.sort_unstable();
            if sorted[3] >= nodes_m.len() || sorted.windows(2).any(|w| w[0] == w[1]) {
                return Err(Error::Invalid("invalid tetrahedral connectivity"));
            }
            for opposite in 0..4 {
                let face = face_without(*tet, opposite);
                if let Some((incidence, vertex)) = faces.get_mut(&face) {
                    *incidence = incidence.saturating_add(1);
                    *vertex = tet[opposite];
                }
            }
        }
        if faces.values().any(|&(count, _)| count != 1) {
            return Err(Error::Invalid("cohesive face must be an exterior trace with one incident tet"));
        }
        let mut full_patches = Vec::with_capacity(count);
        let mut patches = Vec::with_capacity(count);
        let mut points = Vec::with_capacity(count);
        for pair in pairs {
            poll(cancelled)?;
            let mut a = pair.side_a;
            a.sort_unstable();
            let mut b = [0; 3];
            for i in 0..3 {
                let matches: Vec<_> = pair.side_b.iter().copied().filter(|&j| nodes_m[a[i]] == nodes_m[j]).collect();
                if matches.len() != 1 { return Err(Error::Invalid("P1 traces must match exactly without ambiguous coordinates")); }
                b[i] = matches[0];
            }
            let mut key_b = pair.side_b;
            key_b.sort_unstable();
            let mut matched = b;
            matched.sort_unstable();
            if matched != key_b { return Err(Error::Invalid("P1 coordinate correspondence is not bijective")); }
            let (axes, area) = frame(nodes_m, a, faces[&a].1, faces[&key_b].1, budget.minimum_face_quality)?;
            for q in 0..3 {
                let shape = core::array::from_fn(|i| if i == q { 2.0/3.0 } else { 1.0/6.0 });
                let mut terms = Vec::with_capacity(18);
                for (side, sign) in [(a, -1.0), (b, 1.0)] {
                    for i in 0..3 { for c in 0..3 {
                        terms.push(JumpTerm { dof: 3*side[i]+c,
                            local: core::array::from_fn(|k| sign*shape[i]*axes[k][c]) });
                    } }
                }
                let reduced = terms.iter().filter_map(|term| reduced_of[term.dof].map(|dof| JumpTerm { dof, local: term.local })).collect();
                let area_weight_m2 = positive(area/3.0)?;
                full_patches.push(CohesivePatch::new(3*pair.id+q as u64, area_weight_m2, pair.law, terms)?);
                patches.push(CohesivePatch::new(3*pair.id+q as u64, area_weight_m2, pair.law, reduced)?);
                let mut location_m = [0.0; 3];
                for c in 0..3 { for i in 0..3 { location_m[c] = finite(location_m[c]+shape[i]*nodes_m[a[i]][c])?; } }
                points.push(SurfacePoint { interface_id: pair.id, point_in_face: q, location_m,
                    area_weight_m2, axes, shape_values: shape, side_b_shape_values: shape,
                    paired_nodes: [a, b] });
            }
        }
        poll(cancelled)?;
        Ok(Self { full_dimension, free_dofs: free_dofs.to_vec(), reduced_of, full_patches, patches, points })
    }

    pub fn free_dofs(&self) -> &[usize] { &self.free_dofs }
    pub fn points(&self) -> &[SurfacePoint] { &self.points }
    pub fn patches(&self) -> &[CohesivePatch<'a>] { &self.patches }
    pub fn full_dimension(&self) -> usize { self.full_dimension }

    /// Lower constrained-only nodal motion to the same point order as patches.
    /// Nonzero entries on free DOFs refuse, preventing double application.
    pub fn prescribed_jumps<C: FnMut() -> bool>(&self, prescribed_m: &[f64], cancelled: &mut C)
        -> Result<Vec<[f64; 3]>, Error> {
        self.check_prescribed(prescribed_m, cancelled)?;
        self.jumps(prescribed_m, cancelled)
    }

    pub fn expand<C: FnMut() -> bool>(&self, reduced_m: &[f64], prescribed_m: &[f64], cancelled: &mut C)
        -> Result<Vec<f64>, Error> {
        self.check_prescribed(prescribed_m, cancelled)?;
        if reduced_m.len() != self.free_dofs.len() { return Err(Error::Invalid("reduced displacement dimension")); }
        let mut full = prescribed_m.to_vec();
        for (i, (&dof, &value)) in self.free_dofs.iter().zip(reduced_m).enumerate() {
            if i % 256 == 0 { poll(cancelled)?; }
            full[dof] = finite(value)?;
        }
        poll(cancelled)?;
        Ok(full)
    }

    /// Pure full-space jump action. This does not advance damage history.
    pub fn jumps<C: FnMut() -> bool>(&self, displacement_m: &[f64], cancelled: &mut C)
        -> Result<Vec<[f64; 3]>, Error> {
        check_vector(displacement_m, self.full_dimension, cancelled)?;
        let mut result = Vec::with_capacity(self.patches.len());
        for patch in &self.full_patches { poll(cancelled)?; result.push(patch.jump(displacement_m, [0.0; 3])?); }
        poll(cancelled)?;
        Ok(result)
    }

    /// Area-weighted transpose: internal nodal forces B^T A t in newtons.
    /// Supply tractions in points() order. These are algebraic inputs, not
    /// implicitly accepted material histories or a complete support reaction.
    pub fn full_forces<C: FnMut() -> bool>(&self, tractions_pa: &[[f64; 3]], cancelled: &mut C)
        -> Result<Vec<f64>, Error> {
        poll(cancelled)?;
        if tractions_pa.len() != self.patches.len() { return Err(Error::Invalid("P1 traction point count")); }
        let mut force = vec![0.0; self.full_dimension];
        for (patch, traction) in self.full_patches.iter().zip(tractions_pa) {
            poll(cancelled)?;
            for &value in traction { finite(value)?; }
            patch.scatter(*traction, &mut force)?;
        }
        poll(cancelled)?;
        Ok(force)
    }

    fn check_prescribed<C: FnMut() -> bool>(&self, values: &[f64], cancelled: &mut C) -> Result<(), Error> {
        check_vector(values, self.full_dimension, cancelled)?;
        for (dof, &value) in values.iter().enumerate() {
            if dof % 256 == 0 { poll(cancelled)?; }
            if self.reduced_of[dof].is_some() && value != 0.0 {
                return Err(Error::Invalid("prescribed motion on a free DOF"));
            }
        }
        Ok(())
    }
}

fn positive(value: f64) -> Result<f64, Error> {
    super::positive(value)?;
    Ok(value)
}
fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn check_vector<C: FnMut() -> bool>(values: &[f64], n: usize, cancelled: &mut C) -> Result<(), Error> {
    poll(cancelled)?;
    if values.len() != n { return Err(Error::Invalid("full P1 vector dimension")); }
    for chunk in values.chunks(256) { poll(cancelled)?; for &v in chunk { finite(v)?; } }
    Ok(())
}
fn face_without(tet: [usize; 4], opposite: usize) -> [usize; 3] {
    let mut face = [0; 3];
    let mut i = 0;
    for (j, node) in tet.into_iter().enumerate() { if j != opposite { face[i] = node; i += 1; } }
    face.sort_unstable();
    face
}
fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { core::array::from_fn(|i| a[i]-b[i]) }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn length(v: [f64; 3]) -> Result<f64, Error> {
    let scale = v.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    positive(scale)?;
    let u = v.map(|x| x/scale);
    positive(scale*fs_math::det::sqrt(dot(u, u)))
}
fn frame(nodes: &[[f64; 3]], a: [usize; 3], opposite_a: usize, opposite_b: usize, quality: f64)
    -> Result<([[f64; 3]; 3], f64), Error> {
    let e1 = subtract(nodes[a[1]], nodes[a[0]]);
    let e2 = subtract(nodes[a[2]], nodes[a[0]]);
    let l1 = length(e1)?;
    let scale = l1.max(length(e2)?).max(length(subtract(nodes[a[2]], nodes[a[1]]))?);
    let normal = cross(e1.map(|x| x/scale), e2.map(|x| x/scale));
    let norm = length(normal)?;
    if norm < quality { return Err(Error::Invalid("degenerate or under-quality P1 face")); }
    let mut normal = normal.map(|x| x/norm);
    let da = finite(dot(normal, subtract(nodes[opposite_a], nodes[a[0]])))?;
    if da == 0.0 { return Err(Error::Invalid("coplanar opposite tetrahedron vertex")); }
    if da > 0.0 { normal = normal.map(|x| -x); }
    if finite(dot(normal, subtract(nodes[opposite_b], nodes[a[0]])))? <= 0.0 {
        return Err(Error::Invalid("paired tetrahedra must be on opposite sides"));
    }
    let tangent = e1.map(|x| x/l1);
    let area = positive((0.5*norm*scale)*scale)?;
    Ok(([normal, tangent, cross(normal, tangent)], area))
}

#[cfg(test)]
mod tests;
