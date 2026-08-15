//! Deterministic top-level acceleration (TLAS) over animated instances
//! (bead frankensim-h7xu5.3.4).
//!
//! Hundreds of frames reuse immutable per-geometry BLAS/chart backends;
//! what changes per frame is instance MOTION. This module builds a
//! deterministic BVH over conservative shutter-swept instance bounds:
//!
//! - Construction order is a pure function of the inputs: instances sort
//!   by swept-centroid along the split axis with object-id tie-breaks,
//!   median split, no randomness, no pointer order.
//! - Node slabs are outward-rounded ([`f64::next_down`]/[`f64::next_up`])
//!   so a bound can never cull by a rounding ulp; non-finite bounds refuse
//!   at admission.
//! - Refit keeps the topology and recomputes slabs bottom-up; it is
//!   admitted only for the same instance set in the same order (anything
//!   else is a rebuild, by typed refusal rather than silent drift).
//! - Traversal uses an explicit stack (never recursion), polls
//!   cancellation at bounded node intervals, and reports counters so
//!   pruning is measurable instead of assumed.
//! - A brute-force oracle with the exact same tie-break doctrine
//!   (closest `t`, then smaller object id) is part of the public surface;
//!   the batteries hold TLAS results equal to it.
//!
//! No-claims: the TLAS accelerates geometry traversal only. Light
//! selection, shading, and sampling are outside it, and a fingerprint
//! equality proves structural identity, not image equality.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::Aabb;

use crate::animated_instances::{
    AnimatedGeometryInstance, AnimatedInstanceError, RigidTransformTrajectory,
};
use crate::instances::InstanceHit;
use crate::motion::{ShutterInterval, TimedRay};
use crate::motion_bounds::{FiniteLocalAabb, conservative_trajectory_swept_aabb};
use crate::tracer::Ray;

/// Versioned domain for TLAS fingerprints.
pub const TLAS_FINGERPRINT_DOMAIN: &str = "org.frankensim.fs-render.animated-tlas.v1";
/// Cancellation poll cadence during construction and traversal.
const CANCEL_POLL_NODES: usize = 64;
/// Leaf size: instances per leaf before splitting stops.
const LEAF_INSTANCES: usize = 2;

/// Typed refusals of the TLAS boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TlasError {
    /// No instances were supplied.
    EmptyScene,
    /// `local_bounds` length does not match the instance slice.
    BoundsArityMismatch,
    /// A swept instance bound was non-finite or inverted.
    NonFiniteSweptBounds {
        /// Offending instance position in the input slice.
        instance_index: u32,
    },
    /// Refit requires the same instance set in the same order.
    RefitShapeMismatch,
    /// Motion-bounds or trajectory admission refused.
    Animated(AnimatedInstanceError),
    /// Execution was cancelled at a bounded boundary.
    Cancelled,
}

impl core::fmt::Display for TlasError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "animated TLAS refusal: {self:?}")
    }
}

impl std::error::Error for TlasError {}

impl From<fs_exec::Cancelled> for TlasError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}

impl From<AnimatedInstanceError> for TlasError {
    fn from(error: AnimatedInstanceError) -> Self {
        match error {
            AnimatedInstanceError::Cancelled => Self::Cancelled,
            other => Self::Animated(other),
        }
    }
}

/// Traversal counters: pruning is measured, never assumed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TlasTraversalCounters {
    /// Interior + leaf nodes visited.
    pub nodes_visited: u64,
    /// Slab tests performed.
    pub aabb_tests: u64,
    /// Instance intersections actually executed.
    pub instance_tests: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum TlasNode {
    Interior {
        aabb: Aabb,
        left: u32,
        right: u32,
    },
    Leaf {
        aabb: Aabb,
        /// Range into `order`.
        start: u32,
        count: u32,
    },
}

impl TlasNode {
    const fn aabb(&self) -> &Aabb {
        match self {
            Self::Interior { aabb, .. } | Self::Leaf { aabb, .. } => aabb,
        }
    }
}

/// Deterministic TLAS over one shutter interval.
#[derive(Clone, Debug, PartialEq)]
pub struct AnimatedTlas {
    nodes: Vec<TlasNode>,
    order: Vec<u32>,
    swept: Vec<Aabb>,
    shutter: ShutterInterval,
    fingerprint: ContentHash,
}

fn outward(aabb: Aabb) -> Aabb {
    Aabb {
        min: fs_geom::Point3::new(
            aabb.min.x.next_down(),
            aabb.min.y.next_down(),
            aabb.min.z.next_down(),
        ),
        max: fs_geom::Point3::new(
            aabb.max.x.next_up(),
            aabb.max.y.next_up(),
            aabb.max.z.next_up(),
        ),
    }
}

fn union(a: &Aabb, b: &Aabb) -> Aabb {
    Aabb {
        min: fs_geom::Point3::new(
            a.min.x.min(b.min.x),
            a.min.y.min(b.min.y),
            a.min.z.min(b.min.z),
        ),
        max: fs_geom::Point3::new(
            a.max.x.max(b.max.x),
            a.max.y.max(b.max.y),
            a.max.z.max(b.max.z),
        ),
    }
}

fn finite(aabb: &Aabb) -> bool {
    aabb.min.x.is_finite()
        && aabb.min.y.is_finite()
        && aabb.min.z.is_finite()
        && aabb.max.x.is_finite()
        && aabb.max.y.is_finite()
        && aabb.max.z.is_finite()
        && aabb.min.x <= aabb.max.x
        && aabb.min.y <= aabb.max.y
        && aabb.min.z <= aabb.max.z
}

fn centroid(aabb: &Aabb) -> [f64; 3] {
    [
        0.5 * (aabb.min.x + aabb.max.x),
        0.5 * (aabb.min.y + aabb.max.y),
        0.5 * (aabb.min.z + aabb.max.z),
    ]
}

/// Conservative slab test over `[0, t_max]`.
fn slab_hit(aabb: &Aabb, ray: &Ray, t_max: f64) -> bool {
    let mut t_enter = 0.0_f64;
    let mut t_exit = t_max;
    let origin = [ray.origin.x, ray.origin.y, ray.origin.z];
    let direction = [ray.dir.x, ray.dir.y, ray.dir.z];
    let minimum = [aabb.min.x, aabb.min.y, aabb.min.z];
    let maximum = [aabb.max.x, aabb.max.y, aabb.max.z];
    for axis in 0..3 {
        if direction[axis] == 0.0 {
            if origin[axis] < minimum[axis] || origin[axis] > maximum[axis] {
                return false;
            }
            continue;
        }
        let inverse = direction[axis].recip();
        let mut near = (minimum[axis] - origin[axis]) * inverse;
        let mut far = (maximum[axis] - origin[axis]) * inverse;
        if near > far {
            core::mem::swap(&mut near, &mut far);
        }
        t_enter = t_enter.max(near);
        t_exit = t_exit.min(far);
        if t_enter > t_exit {
            return false;
        }
    }
    true
}

fn swept_instance_bounds(
    instances: &[AnimatedGeometryInstance],
    local_bounds: &[FiniteLocalAabb],
    shutter: ShutterInterval,
) -> Result<Vec<Aabb>, TlasError> {
    if instances.is_empty() {
        return Err(TlasError::EmptyScene);
    }
    if instances.len() != local_bounds.len() {
        return Err(TlasError::BoundsArityMismatch);
    }
    let mut swept = Vec::with_capacity(instances.len());
    for (index, (instance, local)) in instances.iter().zip(local_bounds).enumerate() {
        let trajectory: &RigidTransformTrajectory = instance.trajectory();
        trajectory.admit_shutter(shutter)?;
        let bounds =
            conservative_trajectory_swept_aabb(*local, trajectory, shutter).map_err(|_| {
                TlasError::NonFiniteSweptBounds {
                    instance_index: index as u32,
                }
            })?;
        let bounds = outward(bounds);
        if !finite(&bounds) {
            return Err(TlasError::NonFiniteSweptBounds {
                instance_index: index as u32,
            });
        }
        swept.push(bounds);
    }
    Ok(swept)
}

impl AnimatedTlas {
    /// Build deterministically over conservative shutter-swept bounds.
    ///
    /// # Errors
    /// Empty scenes, bounds-arity mismatch, non-finite swept bounds,
    /// shutter/trajectory admission refusals, and cancellation.
    pub fn build(
        cx: &Cx<'_>,
        instances: &[AnimatedGeometryInstance],
        local_bounds: &[FiniteLocalAabb],
        shutter: ShutterInterval,
    ) -> Result<Self, TlasError> {
        cx.checkpoint()?;
        let swept = swept_instance_bounds(instances, local_bounds, shutter)?;
        let mut order: Vec<u32> = (0..instances.len() as u32).collect();
        let mut nodes = Vec::new();
        build_node(
            cx,
            &swept,
            instances,
            &mut order,
            &mut nodes,
            0,
            instances.len(),
        )?;
        let fingerprint = fingerprint_of(&nodes, &order, &swept, instances, shutter);
        Ok(Self {
            nodes,
            order,
            swept,
            shutter,
            fingerprint,
        })
    }

    /// Structural fingerprint: shutter, ordered instance identities, swept
    /// slab bits, and topology. Equal fingerprints mean an identical
    /// acceleration structure, nothing about images.
    #[must_use]
    pub const fn fingerprint(&self) -> ContentHash {
        self.fingerprint
    }

    /// Node count (structure size, for refit accounting).
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Refit: recompute swept slabs for the SAME instances in the SAME
    /// order, keeping topology. Returns the maximum per-instance slab
    /// inflation ratio (new half-extent sum over old) so callers can
    /// decide when accumulated drift warrants a rebuild.
    ///
    /// # Errors
    /// Shape mismatch (different instance count is a rebuild, not a
    /// refit), plus every build-side admission refusal.
    pub fn refit(
        &mut self,
        cx: &Cx<'_>,
        instances: &[AnimatedGeometryInstance],
        local_bounds: &[FiniteLocalAabb],
        shutter: ShutterInterval,
    ) -> Result<f64, TlasError> {
        cx.checkpoint()?;
        if instances.len() != self.swept.len() {
            return Err(TlasError::RefitShapeMismatch);
        }
        let swept = swept_instance_bounds(instances, local_bounds, shutter)?;
        let mut max_inflation = 0.0_f64;
        for (new, old) in swept.iter().zip(&self.swept) {
            let extent = |aabb: &Aabb| {
                (aabb.max.x - aabb.min.x) + (aabb.max.y - aabb.min.y) + (aabb.max.z - aabb.min.z)
            };
            let old_extent = extent(old).max(f64::MIN_POSITIVE);
            max_inflation = max_inflation.max(extent(new) / old_extent);
        }
        self.swept = swept;
        // Bottom-up slab refresh: children precede parents in `nodes` by
        // construction order? They do not (parents precede children), so
        // walk indices in reverse: every child index is greater than its
        // parent's, making reverse order a valid bottom-up pass.
        for index in (0..self.nodes.len()).rev() {
            if index % CANCEL_POLL_NODES == 0 {
                cx.checkpoint()?;
            }
            let refreshed = match self.nodes[index] {
                TlasNode::Leaf { start, count, .. } => {
                    let mut aabb = self.swept[self.order[start as usize] as usize];
                    for slot in (start + 1)..(start + count) {
                        aabb = union(&aabb, &self.swept[self.order[slot as usize] as usize]);
                    }
                    TlasNode::Leaf { aabb, start, count }
                }
                TlasNode::Interior { left, right, .. } => {
                    let aabb = union(
                        self.nodes[left as usize].aabb(),
                        self.nodes[right as usize].aabb(),
                    );
                    TlasNode::Interior { aabb, left, right }
                }
            };
            self.nodes[index] = refreshed;
        }
        self.shutter = shutter;
        self.fingerprint =
            fingerprint_of(&self.nodes, &self.order, &self.swept, instances, shutter);
        Ok(max_inflation)
    }

    /// Closest hit through the TLAS with the shared tie-break doctrine
    /// (closest `t`; equal `t` resolves to the smaller object id).
    ///
    /// # Errors
    /// Instance admission/intersection refusals and cancellation.
    pub fn intersect(
        &self,
        cx: &Cx<'_>,
        instances: &[AnimatedGeometryInstance],
        timed_ray: &TimedRay<Ray>,
        t_max: f64,
        eps: f64,
        counters: &mut TlasTraversalCounters,
    ) -> Result<Option<InstanceHit>, TlasError> {
        cx.checkpoint()?;
        let ray = timed_ray.spatial();
        let mut best: Option<InstanceHit> = None;
        let mut best_t = t_max;
        // Explicit stack: arena walkers must iterate, never recurse.
        let mut stack: Vec<u32> = vec![0];
        let mut visited = 0usize;
        while let Some(index) = stack.pop() {
            visited += 1;
            if visited % CANCEL_POLL_NODES == 0 {
                cx.checkpoint()?;
            }
            counters.nodes_visited += 1;
            counters.aabb_tests += 1;
            let node = &self.nodes[index as usize];
            if !slab_hit(node.aabb(), ray, best_t) {
                continue;
            }
            match node {
                TlasNode::Interior { left, right, .. } => {
                    stack.push(*right);
                    stack.push(*left);
                }
                TlasNode::Leaf { start, count, .. } => {
                    for slot in *start..(*start + *count) {
                        let instance_index = self.order[slot as usize] as usize;
                        counters.instance_tests += 1;
                        let candidate =
                            instances[instance_index].intersect(cx, timed_ray, best_t, eps)?;
                        if let Some(candidate) = candidate {
                            let wins = match &best {
                                None => true,
                                Some(current) => {
                                    candidate.hit.t < current.hit.t
                                        || (candidate.hit.t == current.hit.t
                                            && candidate.object_id < current.object_id)
                                }
                            };
                            if wins {
                                best_t = candidate.hit.t;
                                best = Some(candidate);
                            }
                        }
                    }
                }
            }
        }
        Ok(best)
    }
}

fn build_node(
    cx: &Cx<'_>,
    swept: &[Aabb],
    instances: &[AnimatedGeometryInstance],
    order: &mut [u32],
    nodes: &mut Vec<TlasNode>,
    start: usize,
    end: usize,
) -> Result<u32, TlasError> {
    if nodes.len() % CANCEL_POLL_NODES == 0 {
        cx.checkpoint()?;
    }
    let mut aabb = swept[order[start] as usize];
    for &slot in &order[start + 1..end] {
        aabb = union(&aabb, &swept[slot as usize]);
    }
    let index = nodes.len() as u32;
    if end - start <= LEAF_INSTANCES {
        nodes.push(TlasNode::Leaf {
            aabb,
            start: start as u32,
            count: (end - start) as u32,
        });
        return Ok(index);
    }
    // Deterministic split: widest centroid axis; sort by (centroid,
    // object id) so equal centroids cannot depend on input order.
    let mut centroid_min = [f64::INFINITY; 3];
    let mut centroid_max = [f64::NEG_INFINITY; 3];
    for &slot in &order[start..end] {
        let center = centroid(&swept[slot as usize]);
        for axis in 0..3 {
            centroid_min[axis] = centroid_min[axis].min(center[axis]);
            centroid_max[axis] = centroid_max[axis].max(center[axis]);
        }
    }
    let mut axis = 0;
    let mut widest = centroid_max[0] - centroid_min[0];
    for candidate in 1..3 {
        let width = centroid_max[candidate] - centroid_min[candidate];
        if width > widest {
            widest = width;
            axis = candidate;
        }
    }
    order[start..end].sort_by(|&a, &b| {
        let ca = centroid(&swept[a as usize])[axis];
        let cb = centroid(&swept[b as usize])[axis];
        ca.total_cmp(&cb).then_with(|| {
            instances[a as usize]
                .object_id()
                .cmp(&instances[b as usize].object_id())
        })
    });
    let middle = start + (end - start) / 2;
    // Placeholder; children patch it after allocation.
    nodes.push(TlasNode::Interior {
        aabb,
        left: 0,
        right: 0,
    });
    let left = build_node(cx, swept, instances, order, nodes, start, middle)?;
    let right = build_node(cx, swept, instances, order, nodes, middle, end)?;
    nodes[index as usize] = TlasNode::Interior { aabb, left, right };
    Ok(index)
}

fn fingerprint_of(
    nodes: &[TlasNode],
    order: &[u32],
    swept: &[Aabb],
    instances: &[AnimatedGeometryInstance],
    shutter: ShutterInterval,
) -> ContentHash {
    let mut payload = Vec::new();
    payload.extend_from_slice(&shutter.open_s().to_le_bytes());
    payload.extend_from_slice(&shutter.close_s().to_le_bytes());
    payload.extend_from_slice(&(instances.len() as u64).to_le_bytes());
    for (instance, aabb) in instances.iter().zip(swept) {
        payload.extend_from_slice(&instance.object_id().to_le_bytes());
        payload.extend_from_slice(instance.geometry_identity().as_bytes());
        for value in [
            aabb.min.x, aabb.min.y, aabb.min.z, aabb.max.x, aabb.max.y, aabb.max.z,
        ] {
            payload.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }
    for &slot in order {
        payload.extend_from_slice(&slot.to_le_bytes());
    }
    for node in nodes {
        match node {
            TlasNode::Interior { left, right, .. } => {
                payload.push(0);
                payload.extend_from_slice(&left.to_le_bytes());
                payload.extend_from_slice(&right.to_le_bytes());
            }
            TlasNode::Leaf { start, count, .. } => {
                payload.push(1);
                payload.extend_from_slice(&start.to_le_bytes());
                payload.extend_from_slice(&count.to_le_bytes());
            }
        }
    }
    fs_blake3::hash_domain(TLAS_FINGERPRINT_DOMAIN, &payload)
}

/// Brute-force oracle: linear scan with the exact TLAS tie-break doctrine.
///
/// # Errors
/// Instance admission/intersection refusals and cancellation.
pub fn brute_force_intersect(
    cx: &Cx<'_>,
    instances: &[AnimatedGeometryInstance],
    timed_ray: &TimedRay<Ray>,
    t_max: f64,
    eps: f64,
) -> Result<Option<InstanceHit>, TlasError> {
    let mut best: Option<InstanceHit> = None;
    for instance in instances {
        let candidate = instance.intersect(cx, timed_ray, t_max, eps)?;
        if let Some(candidate) = candidate {
            let wins = match &best {
                None => true,
                Some(current) => {
                    candidate.hit.t < current.hit.t
                        || (candidate.hit.t == current.hit.t
                            && candidate.object_id < current.object_id)
                }
            };
            if wins {
                best = Some(candidate);
            }
        }
    }
    Ok(best)
}
