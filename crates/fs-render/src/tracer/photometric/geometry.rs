//! Prepared broad phase over the original calibrated scene primitives.
//!
//! A spatial hierarchy identifies possible leaves. A bounded bitmap then
//! visits those leaves in their original order, retaining the native closest
//! hit, mixed mesh/instance tie rules, BLAS bounds and hit diagnostics. The
//! index owns metadata only; original geometry and material witnesses remain
//! borrowed by the scene and shared path jobs.

use super::*;
use crate::motion_bounds::{
    FiniteLocalAabb, RigidMotionSegment, RotationSweep, conservative_world_swept_aabb,
};
use crate::tlas::MAX_TLAS_INSTANCES;
use fs_geom::Aabb;

const MIN_INDEXED_PRIMITIVES: usize = 8;
const LEAF_SIZE: usize = 4;
const BITMAP_WORDS: usize = MAX_TLAS_INSTANCES.div_ceil(64);
const STACK_CAPACITY: usize = 32;
// For near-unit rays this loose broad-phase margin covers the native quaternion
// rotate/add chain, inverse-ray evaluation, dominant-axis parameter scaling and
// return to world coordinates. It is not a physical/numerical error certificate;
// underflow/overflow regimes and rays outside that norm envelope use the leaf
// path directly instead of relying on this filter.
const ROUNDOFF_MARGIN: f64 = 1024.0 * f64::EPSILON;
const SAFE_COORDINATE: f64 = f64::MAX / 64.0;
const UNIT_NORM_TOLERANCE: f64 = 32.0 * f64::EPSILON;
const MIN_ROUNDOFF_SCALE: f64 = f64::MIN_POSITIVE / ROUNDOFF_MARGIN;

#[derive(Clone, Copy)]
struct Entry {
    primitive: u32,
    bounds: Aabb,
}

#[derive(Clone, Copy)]
struct Node {
    bounds: Aabb,
    // Nonempty leaves use start/count into entries. Interior nodes use left
    // and right; count zero distinguishes them without a pointer or allocation.
    start: u32,
    count: u32,
    left: u32,
    right: u32,
}

struct Hierarchy {
    entries: Vec<Entry>,
    nodes: Vec<Node>,
    always: Vec<u32>,
    primitive_count: usize,
    maximum_coordinate: f64,
}

/// Immutable reusable scene broad phase, including an ordinary linear fallback.
pub(crate) struct PreparedGeometry {
    hierarchy: Option<Hierarchy>,
}

impl PreparedGeometry {
    pub(crate) fn try_new(primitives: &[Primitive], cx: &Cx<'_>) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if !(MIN_INDEXED_PRIMITIVES..=MAX_TLAS_INSTANCES).contains(&primitives.len()) {
            return Ok(Self { hierarchy: None });
        }
        let mut entries = Vec::new();
        let mut always = Vec::new();
        entries
            .try_reserve_exact(primitives.len())
            .map_err(|_| PhotometryError::AllocationRefused)?;
        always
            .try_reserve_exact(primitives.len())
            .map_err(|_| PhotometryError::AllocationRefused)?;
        let mut maximum_coordinate = 0.0_f64;
        for (primitive, object) in primitives.iter().enumerate() {
            cx.checkpoint()?;
            if let Some(bounds) = primitive_bounds(object, cx)? {
                maximum_coordinate = maximum_coordinate.max(coordinate_scale(bounds));
                entries.push(Entry {
                    primitive: primitive as u32,
                    bounds,
                });
            } else {
                // An uncertain/unrepresentable bound must never erase a native
                // hit or backend refusal. Such a primitive is always visited.
                always.push(primitive as u32);
            }
        }
        if entries.len() < MIN_INDEXED_PRIMITIVES || maximum_coordinate > SAFE_COORDINATE {
            return Ok(Self { hierarchy: None });
        }
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(2 * entries.len() - 1)
            .map_err(|_| PhotometryError::AllocationRefused)?;
        let count = entries.len();
        build_node(&mut entries, &mut nodes, 0, count, cx)?;
        cx.checkpoint()?;
        Ok(Self {
            hierarchy: Some(Hierarchy {
                entries,
                nodes,
                always,
                primitive_count: primitives.len(),
                maximum_coordinate,
            }),
        })
    }

    pub(crate) fn intersect(
        &self,
        primitives: &[Primitive],
        cx: &Cx<'_>,
        ray: &Ray,
        distance: f64,
    ) -> Result<Option<(usize, Hit)>, TracerError> {
        self.intersect_scene(primitives, cx, ray, distance)
            .map(|hit| hit.map(|hit| (hit.primitive_index, hit.hit)))
    }

    pub(in crate::tracer) fn intersect_scene(
        &self,
        primitives: &[Primitive],
        cx: &Cx<'_>,
        ray: &Ray,
        distance: f64,
    ) -> Result<Option<SceneIntersection>, TracerError> {
        cx.checkpoint()?;
        let linear = || intersect_primitives(primitives, cx, ray, None, Some(distance));
        let Some(hierarchy) = &self.hierarchy else {
            return linear();
        };
        let origin_scale = ray
            .origin
            .x
            .abs()
            .max(ray.origin.y.abs())
            .max(ray.origin.z.abs());
        let direction_norm_squared = ray.dir.dot(ray.dir);
        if primitives.len() != hierarchy.primitive_count
            || !distance.is_finite()
            || distance < f64::MIN_POSITIVE
            || [
                ray.origin.x,
                ray.origin.y,
                ray.origin.z,
                ray.dir.x,
                ray.dir.y,
                ray.dir.z,
            ]
            .into_iter()
            .any(|value| !value.is_finite())
            || origin_scale > SAFE_COORDINATE
            || !direction_norm_squared.is_finite()
            || (direction_norm_squared - 1.0).abs() > UNIT_NORM_TOLERANCE
        {
            // Preserve the leaf's existing refusal/parameterization behavior
            // when finite rigid arithmetic is not safely bounded by this filter.
            return linear();
        }
        let roundoff_scale = origin_scale.max(hierarchy.maximum_coordinate);
        if roundoff_scale < MIN_ROUNDOFF_SCALE {
            return linear();
        }
        let padding = roundoff_scale * ROUNDOFF_MARGIN;
        let mut candidates = [0_u64; BITMAP_WORDS];
        for &index in &hierarchy.always {
            cx.checkpoint()?;
            candidates[index as usize / 64] |= 1_u64 << (index % 64);
        }
        let mut stack = [0_u32; STACK_CAPACITY];
        let mut pending = 1;
        while pending > 0 {
            cx.checkpoint()?;
            pending -= 1;
            let node = hierarchy.nodes[stack[pending] as usize];
            if !may_intersect(node.bounds, ray, distance, padding) {
                continue;
            }
            if node.count > 0 {
                for entry in
                    &hierarchy.entries[node.start as usize..(node.start + node.count) as usize]
                {
                    candidates[entry.primitive as usize / 64] |= 1_u64 << (entry.primitive % 64);
                }
            } else {
                // Median splitting and the admitted count bound need at most
                // 17 pending nodes. Retain a full native fallback if that
                // invariant ever changes; a truncated traversal is not a miss.
                if pending + 2 > STACK_CAPACITY {
                    return linear();
                }
                stack[pending] = node.right;
                stack[pending + 1] = node.left;
                pending += 2;
            }
        }
        let mut best: Option<SceneIntersection> = None;
        for (word_index, &mask) in candidates[..hierarchy.primitive_count.div_ceil(64)]
            .iter()
            .enumerate()
        {
            cx.checkpoint()?;
            let mut remaining = mask;
            while remaining != 0 {
                let index = word_index * 64 + remaining.trailing_zeros() as usize;
                remaining &= remaining - 1;
                let current_best_t = best.as_ref().map(|hit| hit.hit.t);
                if let Some(candidate) = intersect_primitive(
                    index,
                    &primitives[index],
                    cx,
                    ray,
                    None,
                    Some(distance),
                    current_best_t,
                )? {
                    if scene_hit_precedes(primitives, &candidate, best.as_ref()) {
                        best = Some(candidate);
                    }
                }
            }
        }
        cx.checkpoint()?;
        Ok(best)
    }

    #[cfg(test)]
    pub(crate) const fn linear_for_test() -> Self {
        Self { hierarchy: None }
    }

    #[cfg(test)]
    pub(in crate::tracer) fn is_accelerated_for_test(&self) -> bool {
        self.hierarchy.is_some()
    }
}

fn primitive_bounds(primitive: &Primitive, cx: &Cx<'_>) -> Result<Option<Aabb>, PhotometryError> {
    let (mesh, placement) = match &primitive.shape {
        Shape::Mesh(mesh) => (mesh, None),
        Shape::Instance(instance) => match instance.geometry() {
            SharedGeometry::Mesh(mesh) => (mesh.as_ref(), Some(instance.transform())),
            SharedGeometry::Chart(_) => return Ok(None),
        },
        Shape::Chart(_) | Shape::AnimatedInstance(_) => return Ok(None),
    };
    let Some(first) = mesh.vertices.first() else {
        return Ok(None);
    };
    let mut low = *first;
    let mut high = *first;
    for vertex in &mesh.vertices {
        cx.checkpoint()?;
        if vertex.iter().any(|coordinate| !coordinate.is_finite()) {
            return Ok(None);
        }
        for axis in 0..3 {
            low[axis] = low[axis].min(vertex[axis]);
            high[axis] = high[axis].max(vertex[axis]);
        }
    }
    let local = Aabb {
        min: Point3::new(low[0], low[1], low[2]),
        max: Point3::new(high[0], high[1], high[2]),
    };
    let Some(placement) = placement else {
        return Ok(Some(local));
    };
    // Reuse the established outward-rounded rigid-motion bound at a constant
    // pose. Its origin sphere covers quaternion evaluation error and arbitrary
    // orientation; identity placement uses the tighter translated local box.
    // No animated geometry or synthetic surface is created.
    let Ok(local) = FiniteLocalAabb::try_new(local) else {
        return Ok(None);
    };
    let Ok(constant) =
        RigidMotionSegment::try_new(placement, placement, 0.0, 1.0, RotationSweep::Constant)
    else {
        return Ok(None);
    };
    Ok(conservative_world_swept_aabb(local, constant).ok())
}

fn union(left: Aabb, right: Aabb) -> Aabb {
    Aabb {
        min: Point3::new(
            left.min.x.min(right.min.x),
            left.min.y.min(right.min.y),
            left.min.z.min(right.min.z),
        ),
        max: Point3::new(
            left.max.x.max(right.max.x),
            left.max.y.max(right.max.y),
            left.max.z.max(right.max.z),
        ),
    }
}

fn coordinate_scale(bounds: Aabb) -> f64 {
    [
        bounds.min.x,
        bounds.min.y,
        bounds.min.z,
        bounds.max.x,
        bounds.max.y,
        bounds.max.z,
    ]
    .into_iter()
    .map(f64::abs)
    .fold(0.0, f64::max)
}

fn center(bounds: Aabb, axis: usize) -> f64 {
    let low = [bounds.min.x, bounds.min.y, bounds.min.z][axis];
    let high = [bounds.max.x, bounds.max.y, bounds.max.z][axis];
    0.5 * low + 0.5 * high
}

fn build_node(
    entries: &mut [Entry],
    nodes: &mut Vec<Node>,
    start: usize,
    end: usize,
    cx: &Cx<'_>,
) -> Result<u32, PhotometryError> {
    cx.checkpoint()?;
    let mut bounds = entries[start].bounds;
    for entry in &entries[start + 1..end] {
        cx.checkpoint()?;
        bounds = union(bounds, entry.bounds);
    }
    let index = nodes.len() as u32;
    nodes.push(Node {
        bounds,
        start: start as u32,
        count: (end - start) as u32,
        left: 0,
        right: 0,
    });
    if end - start <= LEAF_SIZE {
        return Ok(index);
    }
    let extent = [
        bounds.max.x - bounds.min.x,
        bounds.max.y - bounds.min.y,
        bounds.max.z - bounds.min.z,
    ];
    let axis = (0..3)
        .max_by(|&a, &b| extent[a].total_cmp(&extent[b]))
        .unwrap_or(0);
    // The original primitive index makes the comparator total, including equal
    // centroids. Sorting uses no heap allocation and has an explicit 65,536
    // entry ceiling; cancellation is checked around each bounded sort/node.
    entries[start..end].sort_unstable_by(|left, right| {
        center(left.bounds, axis)
            .total_cmp(&center(right.bounds, axis))
            .then(left.primitive.cmp(&right.primitive))
    });
    cx.checkpoint()?;
    let middle = start + (end - start) / 2;
    let left = build_node(entries, nodes, start, middle, cx)?;
    let right = build_node(entries, nodes, middle, end, cx)?;
    nodes[index as usize].count = 0;
    nodes[index as usize].left = left;
    nodes[index as usize].right = right;
    Ok(index)
}

fn may_intersect(bounds: Aabb, ray: &Ray, distance: f64, padding: f64) -> bool {
    let low = [bounds.min.x, bounds.min.y, bounds.min.z].map(|value| (value - padding).next_down());
    let high = [bounds.max.x, bounds.max.y, bounds.max.z].map(|value| (value + padding).next_up());
    if low.into_iter().chain(high).any(|value| !value.is_finite()) {
        return true;
    }
    // The mesh backend's outward difference/division slabs preserve tangencies
    // and inclusive finite endpoints; the broader padding covers the additional
    // world/local arithmetic used by native rigid instances.
    crate::charts::slab_hit(ray, low, high, distance)
}
