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
use crate::charts::Ray;
use crate::instances::InstanceHit;
use crate::motion::{ShutterInterval, TimedRay};
use crate::motion_bounds::{FiniteLocalAabb, conservative_trajectory_swept_aabb};

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

/// Conservatively derive one instance's LOCAL bounds from its geometry
/// (bead clause: "for analytic charts, obtain or conservatively derive
/// local bounds").
///
/// - Charts: [`fs_geom::Chart::support`], which the trait contract already
///   makes a conservative containment; an unbounded support refuses typed
///   (an instance that fills space cannot be TLAS-bounded).
/// - Meshes: the exact fold of the vertex set.
///
/// # Errors
/// [`TlasError::NonFiniteSweptBounds`]-class refusal via the caller;
/// this function refuses unbounded chart supports and empty meshes.
pub fn derived_local_bounds(
    geometry: &crate::instances::SharedGeometry,
) -> Result<FiniteLocalAabb, TlasError> {
    let aabb = match geometry {
        crate::instances::SharedGeometry::Chart(chart) => chart.support(),
        crate::instances::SharedGeometry::Mesh(mesh) => {
            let vertices = &mesh.vertices;
            let Some(first) = vertices.first() else {
                return Err(TlasError::EmptyScene);
            };
            let mut minimum = *first;
            let mut maximum = *first;
            for vertex in &vertices[1..] {
                for axis in 0..3 {
                    minimum[axis] = minimum[axis].min(vertex[axis]);
                    maximum[axis] = maximum[axis].max(vertex[axis]);
                }
            }
            Aabb {
                min: fs_geom::Point3::new(minimum[0], minimum[1], minimum[2]),
                max: fs_geom::Point3::new(maximum[0], maximum[1], maximum[2]),
            }
        }
    };
    FiniteLocalAabb::try_new(aabb).map_err(|_| TlasError::NonFiniteSweptBounds {
        instance_index: u32::MAX,
    })
}

impl AnimatedTlas {
    /// Build with per-instance local bounds DERIVED from each instance's
    /// own geometry (charts via `support()`, meshes via the vertex fold).
    ///
    /// # Errors
    /// Derivation refusals plus every [`AnimatedTlas::build`] refusal.
    pub fn build_derived(
        cx: &Cx<'_>,
        instances: &[AnimatedGeometryInstance],
        shutter: ShutterInterval,
    ) -> Result<Self, TlasError> {
        let mut local_bounds = Vec::with_capacity(instances.len());
        for (index, instance) in instances.iter().enumerate() {
            let bounds = derived_local_bounds(instance.geometry()).map_err(|error| {
                if let TlasError::NonFiniteSweptBounds { .. } = error {
                    TlasError::NonFiniteSweptBounds {
                        instance_index: index as u32,
                    }
                } else {
                    error
                }
            })?;
            local_bounds.push(bounds);
        }
        Self::build(cx, instances, &local_bounds, shutter)
    }
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

#[cfg(test)]
mod tests {
    use asupersync::types::Budget;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_blake3::hash_domain;
    use fs_exec::{CancelGate, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};

    use super::*;
    use crate::animated_instances::TransformKeyframe;
    use crate::charts::TriMesh;
    use crate::instances::RigidTransform;
    use crate::instances::SharedGeometry;
    use crate::motion::{
        NormalizedShutterTime, ShotTimeBounds, ShutterConvention, ShutterDistribution,
    };

    fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 71,
                    kernel_id: 5,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&gate, &cx)
        })
    }

    fn shutter(open_s: f64, duration_s: f64) -> ShutterInterval {
        ShutterInterval::resolve(
            open_s,
            duration_s,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::UniformCounterV1,
            ShotTimeBounds::try_new(-10.0, 10.0).unwrap(),
        )
        .unwrap()
    }

    fn keyframe(time_s: f64, translation_m: [f64; 3], velocity: [f64; 3]) -> TransformKeyframe {
        TransformKeyframe::try_new(
            time_s,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], translation_m).unwrap(),
            velocity,
        )
        .unwrap()
    }

    /// A unit triangle in the local z=0 plane, moved along a per-instance
    /// x trajectory. `speed = 0` gives a static instance.
    fn instance(object_id: u64, base_x: f64, speed: f64) -> AnimatedGeometryInstance {
        let mesh = TriMesh::new(
            vec![[-0.4, -0.4, 0.0], [0.4, -0.4, 0.0], [0.0, 0.4, 0.0]],
            vec![[0, 1, 2]],
        );
        let identity = hash_domain(
            "org.frankensim.test.tlas-instance",
            &object_id.to_le_bytes(),
        );
        let trajectory = RigidTransformTrajectory::try_new(vec![
            keyframe(0.0, [base_x, 0.0, 0.0], [speed, 0.0, 0.0]),
            keyframe(2.0, [base_x + 2.0 * speed, 0.0, 0.0], [speed, 0.0, 0.0]),
        ])
        .unwrap();
        AnimatedGeometryInstance::try_new(
            object_id,
            identity,
            SharedGeometry::mesh(mesh),
            trajectory,
        )
        .unwrap()
    }

    fn local_bounds(count: usize) -> Vec<FiniteLocalAabb> {
        (0..count)
            .map(|_| {
                FiniteLocalAabb::try_new(Aabb {
                    min: Point3::new(-0.4, -0.4, 0.0),
                    max: Point3::new(0.4, 0.4, 0.0),
                })
                .unwrap()
            })
            .collect()
    }

    fn ray_at(x: f64, normalized_time: f64, exposure: ShutterInterval) -> TimedRay<Ray> {
        TimedRay::at_normalized(
            Ray {
                origin: Point3::new(x, 0.0, 5.0),
                dir: Vec3::new(0.0, 0.0, -1.0),
            },
            exposure,
            NormalizedShutterTime::try_new(normalized_time).unwrap(),
        )
    }

    fn scene(speeds: &[f64]) -> Vec<AnimatedGeometryInstance> {
        speeds
            .iter()
            .enumerate()
            .map(|(index, &speed)| instance(index as u64 + 1, 3.0 * index as f64, speed))
            .collect()
    }

    #[test]
    fn static_and_moving_hits_match_the_brute_force_oracle() {
        with_cx(|_, cx| {
            let exposure = shutter(0.0, 2.0);
            // Mixed static and moving instances spread along x.
            let instances = scene(&[0.0, 0.5, 0.0, 1.0, 0.25, 0.0, 0.75, 0.0]);
            let bounds = local_bounds(instances.len());
            let tlas = AnimatedTlas::build(cx, &instances, &bounds, exposure).unwrap();
            let mut counters = TlasTraversalCounters::default();
            // A grid of rays across space and shutter time: every TLAS
            // answer must equal the oracle bit-for-bit (motion bounds may
            // never cull a sampled hit).
            for step in 0..48 {
                let x = -1.0 + 0.5 * f64::from(step);
                for time_step in 0..5 {
                    let normalized = f64::from(time_step) / 4.0;
                    let ray = ray_at(x, normalized, exposure);
                    let fast = tlas
                        .intersect(cx, &instances, &ray, 100.0, 1.0e-9, &mut counters)
                        .unwrap();
                    let oracle =
                        brute_force_intersect(cx, &instances, &ray, 100.0, 1.0e-9).unwrap();
                    assert_eq!(fast, oracle, "x={x} normalized={normalized}");
                }
            }
            assert!(counters.instance_tests > 0);
        });
    }

    #[test]
    fn traversal_prunes_instances_measurably() {
        with_cx(|_, cx| {
            let exposure = shutter(0.0, 2.0);
            let instances = scene(&[0.0; 8]);
            let bounds = local_bounds(instances.len());
            let tlas = AnimatedTlas::build(cx, &instances, &bounds, exposure).unwrap();
            let mut counters = TlasTraversalCounters::default();
            // One ray straight down onto instance 1's triangle: far
            // instances must be pruned by slabs, not intersected.
            let ray = ray_at(0.0, 0.0, exposure);
            let hit = tlas
                .intersect(cx, &instances, &ray, 100.0, 1.0e-9, &mut counters)
                .unwrap()
                .unwrap();
            assert_eq!(hit.object_id, 1);
            assert!(
                counters.instance_tests <= 2,
                "8 separated instances must not all be tested: {counters:?}"
            );
        });
    }

    #[test]
    fn fingerprints_are_deterministic_and_motion_sensitive() {
        with_cx(|_, cx| {
            let exposure = shutter(0.0, 2.0);
            let instances = scene(&[0.0, 0.5, 0.0]);
            let bounds = local_bounds(instances.len());
            let first = AnimatedTlas::build(cx, &instances, &bounds, exposure).unwrap();
            let again = AnimatedTlas::build(cx, &instances, &bounds, exposure).unwrap();
            assert_eq!(first.fingerprint(), again.fingerprint(), "identical builds");
            // A different trajectory moves the fingerprint.
            let moved = scene(&[0.0, 0.6, 0.0]);
            let other = AnimatedTlas::build(cx, &moved, &bounds, exposure).unwrap();
            assert_ne!(first.fingerprint(), other.fingerprint());
        });
    }

    #[test]
    fn refit_matches_a_fresh_build_and_reports_inflation() {
        with_cx(|_, cx| {
            let exposure = shutter(0.0, 2.0);
            let instances = scene(&[0.0, 0.5, 0.0, 1.0]);
            let bounds = local_bounds(instances.len());
            let mut tlas = AnimatedTlas::build(cx, &instances, &bounds, exposure).unwrap();
            // Same inputs refit: inflation ~1, identical fingerprint.
            let inflation = tlas.refit(cx, &instances, &bounds, exposure).unwrap();
            assert!(
                (inflation - 1.0).abs() < 1.0e-9,
                "same-input inflation {inflation}"
            );
            let fresh = AnimatedTlas::build(cx, &instances, &bounds, exposure).unwrap();
            assert_eq!(tlas.fingerprint(), fresh.fingerprint());
            // Faster trajectories refit: results still match the oracle
            // (topology kept, slabs conservative for the NEW motion).
            let faster = scene(&[0.0, 0.9, 0.0, 1.4]);
            let inflation = tlas.refit(cx, &faster, &bounds, exposure).unwrap();
            assert!(inflation > 1.0, "faster sweep must inflate: {inflation}");
            let mut counters = TlasTraversalCounters::default();
            for step in 0..24 {
                let x = -1.0 + 0.6 * f64::from(step);
                let ray = ray_at(x, 0.75, exposure);
                let fast = tlas
                    .intersect(cx, &faster, &ray, 100.0, 1.0e-9, &mut counters)
                    .unwrap();
                let oracle = brute_force_intersect(cx, &faster, &ray, 100.0, 1.0e-9).unwrap();
                assert_eq!(fast, oracle, "refit x={x}");
            }
            // Different instance count is a rebuild, never a refit.
            let fewer = scene(&[0.0]);
            assert_eq!(
                tlas.refit(cx, &fewer, &local_bounds(1), exposure)
                    .unwrap_err(),
                TlasError::RefitShapeMismatch
            );
        });
    }

    #[test]
    fn admission_refuses_degenerate_scenes_and_cancellation_holds() {
        with_cx(|gate, cx| {
            let exposure = shutter(0.0, 2.0);
            assert_eq!(
                AnimatedTlas::build(cx, &[], &[], exposure).unwrap_err(),
                TlasError::EmptyScene
            );
            let instances = scene(&[0.0, 0.5]);
            assert_eq!(
                AnimatedTlas::build(cx, &instances, &local_bounds(1), exposure).unwrap_err(),
                TlasError::BoundsArityMismatch
            );
            // A shutter outside the trajectories refuses through admission.
            assert!(matches!(
                AnimatedTlas::build(cx, &instances, &local_bounds(2), shutter(5.0, 1.0))
                    .unwrap_err(),
                TlasError::Animated(_)
            ));
            // Cancellation observes at the bounded boundary.
            gate.request();
            assert_eq!(
                AnimatedTlas::build(cx, &instances, &local_bounds(2), exposure).unwrap_err(),
                TlasError::Cancelled
            );
        });
    }

    #[test]
    fn derived_bounds_match_manual_bounds_for_meshes_and_admit_charts() {
        with_cx(|_, cx| {
            let exposure = shutter(0.0, 2.0);
            let instances = scene(&[0.0, 0.5, 0.0]);
            // Mesh derivation equals the manual local bounds: identical
            // fingerprints all the way through.
            let manual = AnimatedTlas::build(cx, &instances, &local_bounds(3), exposure).unwrap();
            let derived = AnimatedTlas::build_derived(cx, &instances, exposure).unwrap();
            assert_eq!(derived.fingerprint(), manual.fingerprint());

            // A finite-support chart instance derives from Chart::support
            // and never culls its own hits (oracle equality).
            #[derive(Clone, Copy)]
            struct SphereChart;
            impl fs_geom::Chart for SphereChart {
                fn eval(&self, x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
                    let distance = (x.x * x.x + x.y * x.y + x.z * x.z).sqrt() - 0.5;
                    fs_geom::ChartSample {
                        signed_distance: distance,
                        gradient: None,
                        lipschitz: None,
                        error: fs_evidence::NumericalCertificate::exact(0.0),
                    }
                }
                fn support(&self) -> Aabb {
                    Aabb {
                        min: Point3::new(-0.5, -0.5, -0.5),
                        max: Point3::new(0.5, 0.5, 0.5),
                    }
                }
                fn trace_step_claim(&self) -> fs_geom::TraceStepClaim {
                    fs_geom::TraceStepClaim::ExactDistance
                }
                fn name(&self) -> &'static str {
                    "tlas-test-sphere"
                }
            }
            let identity = hash_domain("org.frankensim.test.tlas-chart", b"sphere");
            let trajectory = RigidTransformTrajectory::try_new(vec![
                keyframe(0.0, [0.0, 0.0, 0.0], [0.0; 3]),
                keyframe(2.0, [0.0, 0.0, 0.0], [0.0; 3]),
            ])
            .unwrap();
            let chart_instance = AnimatedGeometryInstance::try_new(
                9,
                identity,
                SharedGeometry::chart(SphereChart),
                trajectory,
            )
            .unwrap();
            let chart_scene = vec![chart_instance];
            let tlas = AnimatedTlas::build_derived(cx, &chart_scene, exposure).unwrap();
            let mut counters = TlasTraversalCounters::default();
            let ray = ray_at(0.0, 0.5, exposure);
            let fast = tlas
                .intersect(cx, &chart_scene, &ray, 100.0, 1.0e-9, &mut counters)
                .unwrap();
            let oracle = brute_force_intersect(cx, &chart_scene, &ray, 100.0, 1.0e-9).unwrap();
            assert_eq!(fast, oracle, "chart TLAS equals the oracle");
            assert!(fast.is_some(), "the derived chart bound admits its own hit");

            // An unbounded chart support refuses typed.
            #[derive(Clone, Copy)]
            struct UnboundedChart;
            impl fs_geom::Chart for UnboundedChart {
                fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
                    fs_geom::ChartSample {
                        signed_distance: -1.0,
                        gradient: None,
                        lipschitz: None,
                        error: fs_evidence::NumericalCertificate::exact(0.0),
                    }
                }
                fn support(&self) -> Aabb {
                    Aabb {
                        min: Point3::new(f64::NEG_INFINITY, -1.0, -1.0),
                        max: Point3::new(f64::INFINITY, 1.0, 1.0),
                    }
                }
                fn name(&self) -> &'static str {
                    "tlas-test-unbounded"
                }
            }
            assert!(matches!(
                derived_local_bounds(&SharedGeometry::chart(UnboundedChart)),
                Err(TlasError::NonFiniteSweptBounds { .. })
            ));
        });
    }

    #[test]
    fn single_instance_scene_works_end_to_end() {
        with_cx(|_, cx| {
            let exposure = shutter(0.0, 2.0);
            let instances = scene(&[0.5]);
            let bounds = local_bounds(1);
            let tlas = AnimatedTlas::build(cx, &instances, &bounds, exposure).unwrap();
            assert_eq!(tlas.node_count(), 1);
            let mut counters = TlasTraversalCounters::default();
            let ray = ray_at(0.5, 0.5, exposure);
            let fast = tlas
                .intersect(cx, &instances, &ray, 100.0, 1.0e-9, &mut counters)
                .unwrap();
            let oracle = brute_force_intersect(cx, &instances, &ray, 100.0, 1.0e-9).unwrap();
            assert_eq!(fast, oracle);
        });
    }
}
