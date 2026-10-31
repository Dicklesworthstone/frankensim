//! Conservative visibility of the existing opaque native mesh snapshot.
//!
//! This is an optional numerical witness, not a second transport integrator.
//! Bounds and triangle IDs are prepared once; vertices, winding, and rigid
//! placements remain owned by the original scene. Its native nominal BVH is
//! deliberately not evidence for a clear ray. This stackless metadata hierarchy
//! uses outward world bounds and interval slab/triangle predicates instead.

use super::{EnclosureUnavailable, cross, dot, finite, hash_float};
use crate::charts::TriMesh;
use crate::instances::{RigidTransform, SharedGeometry};
use crate::photometry::PhotometryError;
use crate::tracer::{Primitive, Shape};
use fs_blake3::DomainHasher;
use fs_exec::Cx;
use fs_ivl::Interval;

/// Admission ceiling for optional visibility metadata (triangles and objects).
/// No ray allocates. A larger scene still has its original nominal result, but
/// no enclosure from this implementation. Only referenced vertices are visited.
pub const MAX_VISIBILITY_TRIANGLES: usize = 262_144;
const LEAF_SIZE: usize = 8;
type Vector = [Interval; 3];

#[derive(Clone, Copy)]
struct Bounds {
    lo: [f64; 3],
    hi: [f64; 3],
}

impl Bounds {
    fn triangle(vertices: [Vector; 3]) -> Option<Self> {
        if vertices.iter().flatten().any(|v| !finite(*v)) {
            return None;
        }
        Some(Self {
            lo: std::array::from_fn(|i| vertices.iter().map(|v| v[i].lo()).fold(f64::INFINITY, f64::min)),
            hi: std::array::from_fn(|i| vertices.iter().map(|v| v[i].hi()).fold(f64::NEG_INFINITY, f64::max)),
        })
    }

    fn union(self, other: Self) -> Self {
        Self {
            lo: std::array::from_fn(|i| self.lo[i].min(other.lo[i])),
            hi: std::array::from_fn(|i| self.hi[i].max(other.hi[i])),
        }
    }

    fn may_intersect(self, origin: Vector, direction: Vector, end: Interval) -> bool {
        let mut lower = 0.0_f64;
        let mut upper = end.hi().max(0.0);
        for axis in 0..3 {
            if direction[axis].contains_zero() {
                // A direction interval straddling zero cannot safely use a
                // nominal reciprocal. Its entire swept coordinate can still
                // exclude this slab; otherwise retain the candidate.
                let swept = origin[axis] + direction[axis] * Interval::new(0.0, upper);
                if swept.hi() < self.lo[axis] || swept.lo() > self.hi[axis] {
                    return false;
                }
            } else {
                let slab = (Interval::new(self.lo[axis], self.hi[axis]) - origin[axis])
                    / direction[axis];
                lower = lower.max(slab.lo());
                upper = upper.min(slab.hi());
                if lower > upper {
                    return false;
                }
            }
        }
        true
    }
}

struct Entry {
    primitive: u32,
    triangle: u32,
    bounds: Bounds,
}

struct Node {
    bounds: Bounds,
    // Preorder traversal: first child is the next node; escape skips this
    // entire subtree. A positive count marks an entry-range leaf.
    escape: usize,
    start: usize,
    count: usize,
}

pub(super) struct VisibilityIndex {
    entries: Vec<Entry>,
    nodes: Vec<Node>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Visibility {
    Clear,
    Blocked,
    Uncertain,
}

impl VisibilityIndex {
    pub(super) fn prepare(
        primitives: &[Primitive],
        hash: &mut DomainHasher,
        cx: &Cx<'_>,
    ) -> Result<Result<Self, EnclosureUnavailable>, PhotometryError> {
        cx.checkpoint()?;
        if primitives.len() > MAX_VISIBILITY_TRIANGLES {
            return Ok(Err(EnclosureUnavailable::GeometryBudgetExceeded));
        }
        let mut count = 0_usize;
        for primitive in primitives {
            cx.checkpoint()?;
            let Some((mesh, _)) = native_mesh(primitive) else {
                return Ok(Err(EnclosureUnavailable::GeometryVisibility));
            };
            let Some(next) = count.checked_add(mesh.triangles.len()) else {
                return Ok(Err(EnclosureUnavailable::GeometryBudgetExceeded));
            };
            if next > MAX_VISIBILITY_TRIANGLES {
                return Ok(Err(EnclosureUnavailable::GeometryBudgetExceeded));
            }
            count = next;
        }
        let mut entries = Vec::new();
        entries.try_reserve_exact(count).map_err(|_| PhotometryError::AllocationRefused)?;
        hash.update(&(primitives.len() as u64).to_le_bytes());
        for (primitive_index, primitive) in primitives.iter().enumerate() {
            let Some((mesh, placement)) = native_mesh(primitive) else {
                return Ok(Err(EnclosureUnavailable::GeometryVisibility));
            };
            hash.update(&[u8::from(placement.is_some())]);
            if let Some(placement) = placement {
                for value in placement.rotation_xyzw().into_iter().chain(placement.translation_m()) {
                    hash_float(hash, value);
                }
            }
            hash.update(&(mesh.triangles.len() as u64).to_le_bytes());
            for (triangle_index, indices) in mesh.triangles.iter().enumerate() {
                cx.checkpoint()?;
                // Bind actual referenced coordinates, not the caller-supplied
                // geometry identity or a cached nominal BVH fingerprint.
                for &index in indices {
                    hash.update(&index.to_le_bytes());
                    let Some(vertex) = mesh.vertices.get(index as usize) else {
                        return Ok(Err(EnclosureUnavailable::GeometryVisibility));
                    };
                    for &value in vertex {
                        hash_float(hash, value);
                    }
                }
                let Some(bounds) = world_triangle(primitive, triangle_index).and_then(Bounds::triangle) else {
                    return Ok(Err(EnclosureUnavailable::GeometryVisibility));
                };
                entries.push(Entry { primitive: primitive_index as u32, triangle: triangle_index as u32, bounds });
            }
        }
        let mut nodes = Vec::new();
        if count > 0 {
            // Except for a single root, median-split terminal leaves contain
            // at least LEAF_SIZE/2 entries. Reserve a proved node ceiling.
            let capacity = 2 * count.div_ceil(LEAF_SIZE / 2) - 1;
            nodes.try_reserve_exact(capacity).map_err(|_| PhotometryError::AllocationRefused)?;
            build(&mut entries, 0, &mut nodes, cx)?;
        }
        cx.checkpoint()?;
        Ok(Ok(Self { entries, nodes }))
    }

    /// Test the represented source-to-receiver segment, 0 < t <= end.
    /// Unknown edge/parallel/overflow cases retain [0,1] transmission. Continue
    /// after an uncertain candidate: another definite blocker still proves zero.
    pub(super) fn visibility(
        &self,
        primitives: &[Primitive],
        origin: Vector,
        direction: Vector,
        end: Interval,
        cx: &Cx<'_>,
    ) -> Result<Visibility, PhotometryError> {
        let mut result = Visibility::Clear;
        let mut node_index = 0;
        while node_index < self.nodes.len() {
            cx.checkpoint()?;
            let node = &self.nodes[node_index];
            if !node.bounds.may_intersect(origin, direction, end) {
                node_index = node.escape;
                continue;
            }
            for entry in &self.entries[node.start..node.start + node.count] {
                cx.checkpoint()?;
                if !entry.bounds.may_intersect(origin, direction, end) {
                    continue;
                }
                let state = world_triangle(&primitives[entry.primitive as usize], entry.triangle as usize)
                    .map_or(Visibility::Uncertain, |triangle| triangle_visibility(origin, direction, end, triangle));
                match state {
                    Visibility::Blocked => return Ok(Visibility::Blocked),
                    Visibility::Uncertain => result = Visibility::Uncertain,
                    Visibility::Clear => {}
                }
            }
            node_index += 1;
        }
        cx.checkpoint()?;
        Ok(result)
    }
}

fn build(entries: &mut [Entry], start: usize, nodes: &mut Vec<Node>, cx: &Cx<'_>) -> Result<(), PhotometryError> {
    cx.checkpoint()?;
    let mut bounds = entries[0].bounds;
    for entry in &entries[1..] {
        cx.checkpoint()?;
        bounds = bounds.union(entry.bounds);
    }
    let id = nodes.len();
    nodes.push(Node { bounds, escape: 0, start, count: 0 });
    if entries.len() <= LEAF_SIZE {
        nodes[id].count = entries.len();
    } else {
        let axis = (0..3).max_by(|&a, &b| {
            (bounds.hi[a] - bounds.lo[a]).total_cmp(&(bounds.hi[b] - bounds.lo[b]))
        }).unwrap_or(0);
        let middle = entries.len() / 2;
        entries.select_nth_unstable_by(middle, |a, b| {
            a.bounds.lo[axis].midpoint(a.bounds.hi[axis])
                .total_cmp(&b.bounds.lo[axis].midpoint(b.bounds.hi[axis]))
                .then_with(|| a.primitive.cmp(&b.primitive))
                .then_with(|| a.triangle.cmp(&b.triangle))
        });
        cx.checkpoint()?;
        let (left, right) = entries.split_at_mut(middle);
        build(left, start, nodes, cx)?;
        build(right, start + middle, nodes, cx)?;
    }
    nodes[id].escape = nodes.len();
    Ok(())
}

fn native_mesh(primitive: &Primitive) -> Option<(&TriMesh, Option<RigidTransform>)> {
    match &primitive.shape {
        Shape::Mesh(mesh) => Some((mesh, None)),
        Shape::Instance(instance) => match instance.geometry() {
            SharedGeometry::Mesh(mesh) => Some((mesh.as_ref(), Some(instance.transform()))),
            SharedGeometry::Chart(_) => None,
        },
        _ => None,
    }
}

fn world_triangle(primitive: &Primitive, triangle: usize) -> Option<[Vector; 3]> {
    let (mesh, placement) = native_mesh(primitive)?;
    let indices = mesh.triangles.get(triangle)?;
    let mut result = [[Interval::point(0.0); 3]; 3];
    for (output, &index) in result.iter_mut().zip(indices) {
        let vertex = *mesh.vertices.get(index as usize)?;
        if vertex.iter().any(|value| !value.is_finite()) { return None; }
        let point = vertex.map(Interval::point);
        *output = match placement {
            None => point,
            Some(transform) => {
                let [x, y, z, w] = transform.rotation_xyzw().map(Interval::point);
                let twice = cross([x, y, z], point).map(|v| Interval::point(2.0) * v);
                let double = cross([x, y, z], twice);
                let translation = transform.translation_m();
                std::array::from_fn(|i| point[i] + w * twice[i] + double[i] + Interval::point(translation[i]))
            }
        };
    }
    Some(result)
}

// Moller-Trumbore determinants evaluated outward, not a nominal hit plus a
// tolerance. Strict tests prove misses/hits; inconclusive predicates widen
// visibility. Winding is irrelevant for an opaque two-sided blocker.
fn triangle_visibility(origin: Vector, direction: Vector, end: Interval, [a, b, c]: [Vector; 3]) -> Visibility {
    let subtract = |left: Vector, right: Vector| std::array::from_fn(|i| left[i] - right[i]);
    let vectors = [subtract(b, a), subtract(c, a), direction, subtract(origin, a)];
    if vectors.iter().flatten().any(|v| !finite(*v)) { return Visibility::Uncertain; }
    let scale = vectors.iter().flatten().map(|v| v.lo().abs().max(v.hi().abs())).fold(0.0_f64, f64::max);
    if scale == 0.0 { return Visibility::Uncertain; }
    // Common scaling leaves all three determinant ratios unchanged and avoids
    // cubic overflow on finite, otherwise well-separated large stage geometry.
    let [e1, e2, d, offset] = vectors.map(|v| v.map(|x| x / Interval::point(scale)));
    let p = cross(d, e2);
    let determinant = dot(e1, p);
    if determinant.contains_zero() || !finite(determinant) { return Visibility::Uncertain; }
    let q = cross(offset, e1);
    let u = dot(offset, p) / determinant;
    let v = dot(d, q) / determinant;
    let t = dot(e2, q) / determinant;
    let uv = u + v;
    if u.hi() < 0.0 || v.hi() < 0.0 || uv.lo() > 1.0 || t.hi() <= 0.0 || t.lo() > end.hi() {
        Visibility::Clear
    } else if u.lo() >= 0.0 && v.lo() >= 0.0 && uv.hi() <= 1.0 && t.lo() > 0.0 && t.hi() <= end.lo() {
        Visibility::Blocked
    } else {
        Visibility::Uncertain
    }
}

#[cfg(test)]
mod tests;
