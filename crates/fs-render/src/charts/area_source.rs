//! Native mesh support for area emission, sharing the original triangle BVH.
//!
//! Preparation rejects positive-area coplanar overlap against every original
//! face, including faces outside a selected source subset. Noncoplanar crossings
//! have zero area and retain the native boundary oracle's contact refusals.
//! Queries identify a unique original face in three dimensions; a projection
//! onto another face's plane is not membership. No mesh conversion is performed.

use super::TriMesh;
use fs_exec::Cx;
use fs_geom::Point3;
use fs_ivl::predicates::{Sign, orient2d, orient3d};
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MeshAreaError {
    InvalidSupport,
    InvalidGeometry,
    WorkLimit { limit: u64 },
    OverlappingFaces { first: usize, second: usize },
    AmbiguousPoint,
    NumericalRange,
    Cancelled,
}

impl core::fmt::Display for MeshAreaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native area support: {self:?}")
    }
}
impl core::error::Error for MeshAreaError {}

struct Work { used: u64, limit: u64 }
impl Work {
    fn charge(&mut self, cx: &Cx<'_>) -> Result<(), MeshAreaError> {
        cx.checkpoint().map_err(|_| MeshAreaError::Cancelled)?;
        if self.used == self.limit { return Err(MeshAreaError::WorkLimit { limit: self.limit }); }
        self.used += 1;
        Ok(())
    }
}

// Native TriMesh construction is a balanced median tree over u32 faces.
// Fixed scratch also refuses a malformed/deeper tree without hidden growth.
struct Stack { nodes: [u32; 64], len: usize }
impl Stack {
    fn root() -> Self { Self { nodes: [0; 64], len: 1 } }
    fn pop(&mut self) -> Option<u32> {
        if self.len == 0 { None } else { self.len -= 1; Some(self.nodes[self.len]) }
    }
    fn push(&mut self, node: u32) -> Result<(), MeshAreaError> {
        let slot = self.nodes.get_mut(self.len).ok_or(MeshAreaError::InvalidGeometry)?;
        *slot = node;
        self.len += 1;
        Ok(())
    }
}

impl TriMesh {
    /// The caller already borrows the immutable mesh admitted by the native
    /// boundary owner. This supplements that owner's BVH/triangle validation
    /// with the disjoint-area requirement of an emitting face subset.
    pub(crate) fn validate_area_source_faces(
        &self, cx: &Cx<'_>, faces: &[usize], max_work: u64,
    ) -> Result<(), MeshAreaError> {
        let mut work = Work { used: 0, limit: max_work };
        work.charge(cx)?;
        if faces.is_empty() || faces.windows(2).any(|w| w[0] >= w[1])
            || faces.last().is_some_and(|&i| i >= self.triangles.len()) || self.nodes.is_empty() {
            return Err(MeshAreaError::InvalidSupport);
        }
        for &first in faces {
            work.charge(cx)?;
            let a = triangle(self, first)?;
            let lower = core::array::from_fn::<_, 3, _>(|i| a[0][i].min(a[1][i]).min(a[2][i]));
            let upper = core::array::from_fn::<_, 3, _>(|i| a[0][i].max(a[1][i]).max(a[2][i]));
            let mut stack = Stack::root();
            while let Some(index) = stack.pop() {
                work.charge(cx)?;
                let node = self.nodes.get(index as usize).ok_or(MeshAreaError::InvalidGeometry)?;
                if (0..3).any(|i| upper[i] < node.lo[i] || lower[i] > node.hi[i]) { continue; }
                if !node.leaf {
                    stack.push(node.b)?;
                    stack.push(node.a)?;
                    continue;
                }
                let end = node.a.checked_add(node.b).ok_or(MeshAreaError::InvalidGeometry)?;
                for slot in node.a..end {
                    work.charge(cx)?;
                    let second = *self.order.get(slot as usize).ok_or(MeshAreaError::InvalidGeometry)? as usize;
                    if second == first || (second < first && faces.binary_search(&second).is_ok()) { continue; }
                    let b = triangle(self, second)?;
                    if (0..3).any(|i| {
                        let lo = b[0][i].min(b[1][i]).min(b[2][i]);
                        let hi = b[0][i].max(b[1][i]).max(b[2][i]);
                        upper[i] < lo || lower[i] > hi
                    }) { continue; }
                    if coplanar_overlap(a, b)? {
                        return Err(MeshAreaError::OverlappingFaces { first, second });
                    }
                }
            }
        }
        cx.checkpoint().map_err(|_| MeshAreaError::Cancelled)
    }

    /// Locate the original face at a numerical native contact. The supplied
    /// local-distance tolerance covers the original contact and rigid inverse
    /// evaluation; it cannot promote edge or multiple-face membership to a hit.
    /// Returned None means no face has supported three-dimensional membership.
    pub(crate) fn area_source_face_at(
        &self, cx: &Cx<'_>, local_point: Point3, tolerance_m: f64, max_work: u64,
    ) -> Result<Option<usize>, MeshAreaError> {
        let p = [local_point.x, local_point.y, local_point.z];
        if !tolerance_m.is_finite() || tolerance_m < 0.0 || p.iter().any(|v| !v.is_finite()) {
            return Err(MeshAreaError::NumericalRange);
        }
        let mut work = Work { used: 0, limit: max_work };
        work.charge(cx)?;
        if self.nodes.is_empty() { return Ok(None); }
        let mut stack = Stack::root();
        let mut found = None;
        while let Some(index) = stack.pop() {
            work.charge(cx)?;
            let node = self.nodes.get(index as usize).ok_or(MeshAreaError::InvalidGeometry)?;
            if (0..3).any(|i| (p[i] + tolerance_m).next_up() < node.lo[i]
                || (p[i] - tolerance_m).next_down() > node.hi[i]) { continue; }
            if !node.leaf {
                stack.push(node.b)?;
                stack.push(node.a)?;
                continue;
            }
            let end = node.a.checked_add(node.b).ok_or(MeshAreaError::InvalidGeometry)?;
            for slot in node.a..end {
                work.charge(cx)?;
                let face = *self.order.get(slot as usize).ok_or(MeshAreaError::InvalidGeometry)? as usize;
                match membership(p, triangle(self, face)?, tolerance_m)? {
                    Membership::Outside => {}
                    Membership::Edge => return Err(MeshAreaError::AmbiguousPoint),
                    Membership::Interior => {
                        if found.replace(face).is_some() { return Err(MeshAreaError::AmbiguousPoint); }
                    }
                }
            }
        }
        cx.checkpoint().map_err(|_| MeshAreaError::Cancelled)?;
        Ok(found)
    }
}

fn triangle(mesh: &TriMesh, face: usize) -> Result<[[f64; 3]; 3], MeshAreaError> {
    let indices = mesh.triangles.get(face).ok_or(MeshAreaError::InvalidGeometry)?;
    let mut vertices = [[0.0; 3]; 3];
    for (out, &index) in vertices.iter_mut().zip(indices) {
        *out = *mesh.vertices.get(index as usize).ok_or(MeshAreaError::InvalidGeometry)?;
    }
    if vertices.iter().flatten().any(|v| !v.is_finite()) { return Err(MeshAreaError::InvalidGeometry); }
    Ok(vertices)
}

// Match the native proximity predicate owner's conservative exponent envelope.
// Exact predicates decide topology only after overflow/underflow are excluded.
fn predicate_domain(points: &[[f64; 3]]) -> Result<(), MeshAreaError> {
    if points.iter().flatten().any(|v| !v.is_finite()
        || (*v != 0.0 && !(1e-60..=1e60).contains(&v.abs()))) {
        return Err(MeshAreaError::NumericalRange);
    }
    for p in points {
        for q in points {
            if (0..3).any(|i| { let d = (p[i] - q[i]).abs(); d != 0.0 && !(1e-60..=1e60).contains(&d) }) {
                return Err(MeshAreaError::NumericalRange);
            }
        }
    }
    Ok(())
}

fn coplanar_overlap(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> Result<bool, MeshAreaError> {
    predicate_domain(&[a[0], a[1], a[2], b[0], b[1], b[2]])?;
    if b.iter().any(|&p| orient3d(a[0], a[1], a[2], p) != Sign::Zero) { return Ok(false); }
    let (u, v) = projection(a)?;
    let a = a.map(|p| [p[u], p[v]]);
    let b = b.map(|p| [p[u], p[v]]);
    // Convex triangle interiors overlap iff none of their six edge lines
    // separates them. A boundary-only intersection has zero area and passes.
    for (own, other) in [(a, b), (b, a)] {
        let inside = orient2d(own[0], own[1], own[2]);
        if inside == Sign::Zero { return Err(MeshAreaError::InvalidGeometry); }
        for edge in 0..3 {
            if other.iter().all(|&p| orient2d(own[edge], own[(edge + 1) % 3], p) != inside) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn projection(t: [[f64; 3]; 3]) -> Result<(usize, usize), MeshAreaError> {
    let e = core::array::from_fn::<_, 3, _>(|i| t[1][i] - t[0][i]);
    let f = core::array::from_fn::<_, 3, _>(|i| t[2][i] - t[0][i]);
    let n = [e[1] * f[2] - e[2] * f[1], e[2] * f[0] - e[0] * f[2], e[0] * f[1] - e[1] * f[0]];
    if n.iter().any(|v| !v.is_finite()) || n.iter().all(|&v| v == 0.0) {
        return Err(MeshAreaError::InvalidGeometry);
    }
    Ok(if n[0].abs() >= n[1].abs() && n[0].abs() >= n[2].abs() { (1, 2) }
        else if n[1].abs() >= n[2].abs() { (0, 2) } else { (0, 1) })
}

enum Membership { Outside, Edge, Interior }
fn membership(p: [f64; 3], t: [[f64; 3]; 3], tolerance: f64) -> Result<Membership, MeshAreaError> {
    let (u, v) = projection(t)?;
    let e = core::array::from_fn::<_, 3, _>(|i| t[1][i] - t[0][i]);
    let f = core::array::from_fn::<_, 3, _>(|i| t[2][i] - t[0][i]);
    let offset = core::array::from_fn::<_, 3, _>(|i| p[i] - t[0][i]);
    let scale = e.iter().chain(f.iter()).fold(0.0_f64, |s, x| s.max(x.abs()));
    if !scale.is_finite() || scale == 0.0 { return Err(MeshAreaError::NumericalRange); }
    let e = e.map(|x| x / scale);
    let f = f.map(|x| x / scale);
    let n = [e[1] * f[2] - e[2] * f[1], e[2] * f[0] - e[0] * f[2], e[0] * f[1] - e[1] * f[0]];
    let length = det::hypot(det::hypot(n[0], n[1]), n[2]);
    if !length.is_finite() || length == 0.0 { return Err(MeshAreaError::NumericalRange); }
    let normal = n.map(|x| x / length);
    let distance = normal.iter().zip(offset).map(|(n, x)| n * x).sum::<f64>().abs();
    let roundoff = 32.0 * f64::EPSILON * (0..3)
        .map(|i| normal[i].abs() * (p[i].abs() + t[0][i].abs() + scale)).sum::<f64>();
    if !distance.is_finite() || !roundoff.is_finite() { return Err(MeshAreaError::NumericalRange); }
    if distance > tolerance + roundoff { return Ok(Membership::Outside); }
    predicate_domain(&[t[0], t[1], t[2], p])?;
    let t = t.map(|x| [x[u], x[v]]);
    let p = [p[u], p[v]];
    let inside = orient2d(t[0], t[1], t[2]);
    if inside == Sign::Zero { return Err(MeshAreaError::InvalidGeometry); }
    let mut unresolved = false;
    for i in 0..3 {
        let a = t[i];
        let b = t[(i + 1) % 3];
        let side = orient2d(a, b, p);
        if !resolved_edge_clearance(a, b, p, tolerance)? || side == Sign::Zero {
            unresolved = true;
        } else if side != inside {
            return Ok(Membership::Outside);
        }
    }
    Ok(if unresolved { Membership::Edge } else { Membership::Interior })
}

// The exact predicate supplies the side, while this conservative clearance
// test establishes that the uncertain inverse/contact point cannot cross the
// edge. A nearby point outside another face must remain ambiguous as well.
// Orthogonal projection cannot enlarge the supplied three-dimensional error.
fn resolved_edge_clearance(
    a: [f64; 2], b: [f64; 2], p: [f64; 2], tolerance: f64,
) -> Result<bool, MeshAreaError> {
    let edge = [b[0] - a[0], b[1] - a[1]];
    let scale = edge[0].abs().max(edge[1].abs());
    if !scale.is_finite() || scale == 0.0 { return Err(MeshAreaError::InvalidGeometry); }
    let edge = edge.map(|x| x / scale);
    let length = det::hypot(edge[0], edge[1]);
    let normal = [-edge[1] / length, edge[0] / length];
    let clearance = (normal[0] * (p[0] - a[0]) + normal[1] * (p[1] - a[1])).abs();
    // Include subtraction, edge normalization and dot-product roundoff. The
    // predicate envelope excludes subnormal edges and overflowing sums.
    let roundoff = 64.0 * f64::EPSILON
        * (p[0].abs() + p[1].abs() + a[0].abs() + a[1].abs() + scale);
    let allowance = (tolerance + roundoff).next_up();
    if !clearance.is_finite() || !allowance.is_finite() {
        return Err(MeshAreaError::NumericalRange);
    }
    Ok(clearance > allowance)
}

#[cfg(test)]
mod tests;
