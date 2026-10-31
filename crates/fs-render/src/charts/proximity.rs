//! Bounded native surface clearance for curved geometric-optics propagation.
//!
//! Approximate closest-point arithmetic selects a useful separating direction;
//! it never certifies a distance by itself. For any direction `v`,
//! `min_i (p-a_i).v / |v|` is a lower bound on distance to a triangle. Outward
//! interval arithmetic evaluates that dual bound. A convex combination of the
//! three vertices supplies an independent upper bound. The same bounds drive
//! BVH pruning, including at edges and on very thin separated faces.

use super::*;
use fs_geom::Aabb;
use fs_ivl::Interval;

/// Native point-to-surface query refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProximityError {
    InvalidInput,
    InvalidGeometry,
    UnsupportedChart,
    AmbiguousContact,
    NumericalRange,
    WorkBudgetExceeded,
    Allocation,
    Cancelled,
}

impl core::fmt::Display for ProximityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native surface proximity: {self:?}")
    }
}
impl core::error::Error for ProximityError {}
impl From<Cancelled> for ProximityError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

/// One cumulative work envelope, shared by scene/instance/mesh queries.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ProximityBudget {
    remaining: u64,
    used: u64,
}
impl ProximityBudget {
    #[must_use]
    pub const fn new(maximum_work: u64) -> Self {
        Self {
            remaining: maximum_work,
            used: 0,
        }
    }
    #[must_use]
    pub const fn used(self) -> u64 {
        self.used
    }
    #[must_use]
    pub const fn remaining(self) -> u64 {
        self.remaining
    }
    pub(crate) fn charge(&mut self, cx: &Cx<'_>) -> Result<(), ProximityError> {
        if self.remaining == 0 {
            return Err(ProximityError::WorkBudgetExceeded);
        }
        self.remaining -= 1;
        self.used += 1;
        if self.used % 64 == 1 {
            cx.checkpoint()?;
        }
        Ok(())
    }
}

/// Sign relative to the outward supporting plane at the selected witness.
/// `Unresolved` includes exact contact and overlapping numerical intervals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProximitySide {
    Positive,
    Negative,
    Unresolved,
}

/// A numerical surface representative, not a straight-ray intersection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceWitness {
    /// Rounded representative of a point on the mesh, or an exact-distance
    /// chart query point whose surface residual is bounded below.
    pub surface_point: Point3,
    /// Radius enclosing a true surface point around `surface_point`.
    pub surface_point_error_m: f64,
    /// Upper bound on distance from the original query to a true surface point.
    pub distance_upper_m: f64,
    /// Unit normal when this backend supplies a regular oriented surface.
    pub normal: Option<Vec3>,
    pub side: ProximitySide,
    /// Original native mesh triangle, absent for a chart.
    pub triangle_index: Option<u32>,
    /// Approximate barycentric coordinates of the retained convex witness.
    pub barycentric: Option<[f64; 3]>,
}

/// Lower clearance to every retained surface and a nearest upper witness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceProximity {
    /// Nonnegative lower distance; +infinity means an empty surface set.
    pub clearance_m: f64,
    pub witness: Option<SurfaceWitness>,
}

pub(crate) fn point_array(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
pub(crate) fn array_point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn iv(a: [f64; 3]) -> [Interval; 3] {
    a.map(Interval::point)
}
fn iv_sub(a: [Interval; 3], b: [Interval; 3]) -> [Interval; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn iv_dot(a: [Interval; 3], b: [Interval; 3]) -> Interval {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn iv_norm(a: [Interval; 3]) -> Interval {
    let a = a.map(Interval::abs);
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Conservative point-to-box lower bound shared by BLAS and scene pruning.
pub(crate) fn box_clearance(point: Point3, bounds: Aabb) -> Result<f64, ProximityError> {
    let p = point_array(point);
    let lo = point_array(bounds.min);
    let hi = point_array(bounds.max);
    if p.iter()
        .chain(lo.iter())
        .chain(hi.iter())
        .any(|v| !v.is_finite())
        || (0..3).any(|i| lo[i] > hi[i])
    {
        return Err(ProximityError::InvalidGeometry);
    }
    let delta = std::array::from_fn(|i| {
        let nearest = p[i].clamp(lo[i], hi[i]);
        Interval::point(p[i]) - Interval::point(nearest)
    });
    Ok(iv_norm(delta).lo().max(0.0))
}

/// Native cancellable nearest-surface traversal over the existing mesh BVH.
impl TriMesh {
    /// Validate the public authoring arrays before using retained BVH bounds.
    /// Repeated queries should use [`MeshPlanarPatches`], whose immutable mesh
    /// borrow amortizes this validation across the complete optical path.
    pub fn nearest_surface_with_cx(
        &self,
        cx: &Cx<'_>,
        point: Point3,
        budget: &mut ProximityBudget,
    ) -> Result<SurfaceProximity, ProximityError> {
        validate_mesh_bvh(self, cx, budget)?;
        self.nearest_surface_filtered(cx, point, budget, None)
    }

    pub(crate) fn nearest_surface_filtered(
        &self,
        cx: &Cx<'_>,
        point: Point3,
        budget: &mut ProximityBudget,
        excluded: Option<(&MeshPlanarPatches<'_>, u32)>,
    ) -> Result<SurfaceProximity, ProximityError> {
        cx.checkpoint()?;
        if point_array(point).iter().any(|v| !v.is_finite()) {
            return Err(ProximityError::InvalidInput);
        }
        let mut result = SurfaceProximity {
            clearance_m: f64::INFINITY,
            witness: None,
        };
        if self.nodes.is_empty() {
            return Ok(result);
        }
        let mut stack = BvhTraversalStack::<BVH_INLINE_STACK_CAPACITY>::new(0)
            .map_err(|_| ProximityError::Allocation)?;
        let mut best_upper = f64::INFINITY;
        while let Some(index) = stack.pop() {
            budget.charge(cx)?;
            if excluded
                .is_some_and(|(patches, patch)| patches.node_patch[index as usize] == Some(patch))
            {
                continue;
            }
            let node = self
                .nodes
                .get(index as usize)
                .ok_or(ProximityError::InvalidGeometry)?;
            let lower = box_clearance(
                point,
                Aabb {
                    min: array_point(node.lo),
                    max: array_point(node.hi),
                },
            )?;
            if lower > best_upper {
                continue;
            }
            if node.leaf {
                for slot in node.a..node.a + node.b {
                    let triangle_index = *self
                        .order
                        .get(slot as usize)
                        .ok_or(ProximityError::InvalidGeometry)?;
                    if excluded.is_some_and(|(patches, patch)| {
                        patches.triangle_patch[triangle_index as usize] == patch
                    }) {
                        continue;
                    }
                    budget.charge(cx)?;
                    let indices = *self
                        .triangles
                        .get(triangle_index as usize)
                        .ok_or(ProximityError::InvalidGeometry)?;
                    let mut vertices = [[0.0; 3]; 3];
                    for i in 0..3 {
                        vertices[i] = *self
                            .vertices
                            .get(indices[i] as usize)
                            .ok_or(ProximityError::InvalidGeometry)?;
                    }
                    let (lower, witness) = triangle_proximity(point, vertices, triangle_index)?;
                    result.clearance_m = result.clearance_m.min(lower);
                    let wins = witness.distance_upper_m < best_upper
                        || (witness.distance_upper_m == best_upper
                            && result
                                .witness
                                .is_none_or(|old| old.triangle_index > witness.triangle_index));
                    if wins {
                        best_upper = witness.distance_upper_m;
                        result.witness = Some(witness);
                    }
                }
            } else {
                // Near-first traversal discovers a tight upper bound early.
                let left = &self.nodes[node.a as usize];
                let right = &self.nodes[node.b as usize];
                let dl = box_clearance(
                    point,
                    Aabb {
                        min: array_point(left.lo),
                        max: array_point(left.hi),
                    },
                )?;
                let dr = box_clearance(
                    point,
                    Aabb {
                        min: array_point(right.lo),
                        max: array_point(right.hi),
                    },
                )?;
                let (near, far) = if dl <= dr {
                    (node.a, node.b)
                } else {
                    (node.b, node.a)
                };
                stack
                    .try_push(far)
                    .map_err(|_| ProximityError::Allocation)?;
                stack
                    .try_push(near)
                    .map_err(|_| ProximityError::Allocation)?;
            }
        }
        cx.checkpoint()?;
        Ok(result)
    }
}

/// A candidate is selected numerically; the returned interval bounds are
/// independent of that candidate being the exact closest point.
fn closest_barycentric(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return [1.0, 0.0, 0.0];
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return [0.0, 1.0, 0.0];
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return [1.0 - v, v, 0.0];
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return [0.0, 0.0, 1.0];
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return [1.0 - w, 0.0, w];
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return [0.0, 1.0 - w, w];
    }
    let denominator = va + vb + vc;
    let v = vb / denominator;
    let w = vc / denominator;
    [1.0 - v - w, v, w]
}

fn convex_bounds(
    p: [f64; 3],
    vertices: [[f64; 3]; 3],
    barycentric: [f64; 3],
) -> Result<(f64, f64, Point3, f64, [f64; 3]), ProximityError> {
    if vertices
        .iter()
        .flatten()
        .chain(p.iter())
        .chain(barycentric.iter())
        .any(|v| !v.is_finite())
    {
        return Err(ProximityError::NumericalRange);
    }
    // Nested convex blends keep the exact real witness inside the triangle
    // even when three rounded barycentric coefficients would not sum to one.
    let v = barycentric[2].clamp(0.0, 1.0);
    let denominator = barycentric[0].max(0.0) + barycentric[1].max(0.0);
    let u = if denominator > 0.0 {
        (barycentric[1].max(0.0) / denominator).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let [a, b, c] = vertices.map(iv);
    let one = Interval::point(1.0);
    let iu = Interval::point(u);
    let jv = Interval::point(v);
    let q: [Interval; 3] =
        std::array::from_fn(|i| ((one - iu) * a[i] + iu * b[i]) * (one - jv) + jv * c[i]);
    let center = q.map(Interval::midpoint);
    let upper = iv_norm(iv_sub(iv(p), q)).hi();
    let error = iv_norm(iv_sub(iv(center), q)).hi();
    let direction = sub(p, center);
    let length = iv_norm(iv(direction)).hi();
    let mut lower = 0.0;
    if length.is_finite() && length > 0.0 {
        let support = vertices
            .into_iter()
            .map(|v| iv_dot(iv_sub(iv(p), iv(v)), iv(direction)).lo())
            .fold(f64::INFINITY, f64::min);
        if support > 0.0 {
            lower = (Interval::point(support) / Interval::point(length))
                .lo()
                .max(0.0);
        }
    }
    if !upper.is_finite() || !error.is_finite() {
        return Err(ProximityError::NumericalRange);
    }
    Ok((
        lower,
        upper,
        array_point(center),
        error,
        [(1.0 - v) * (1.0 - u), (1.0 - v) * u, v],
    ))
}

fn triangle_proximity(
    point: Point3,
    vertices: [[f64; 3]; 3],
    triangle: u32,
) -> Result<(f64, SurfaceWitness), ProximityError> {
    let p = point_array(point);
    let barycentric = closest_barycentric(p, vertices[0], vertices[1], vertices[2]);
    let (lower, upper, surface_point, error, barycentric) =
        convex_bounds(p, vertices, barycentric)?;
    let ab = iv_sub(iv(vertices[1]), iv(vertices[0]));
    let ac = iv_sub(iv(vertices[2]), iv(vertices[0]));
    let normal = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    if iv_norm(normal).lo() <= 0.0 {
        return Err(ProximityError::InvalidGeometry);
    }
    let components = normal.map(Interval::midpoint);
    let length = dot(components, components).sqrt();
    if !length.is_finite() || length <= 0.0 {
        return Err(ProximityError::NumericalRange);
    }
    let signed = iv_dot(iv_sub(iv(p), iv(vertices[0])), normal);
    // Near a face interior, the rounded closest point can carry a tiny
    // tangential error. Its dual direction remains sound but becomes weak
    // when that error dominates the squared normal gap. The supporting plane
    // is a second independent separating direction, linear in that gap.
    let plane_numerator = signed.lo().max(-signed.hi()).max(0.0);
    let plane_lower = (Interval::point(plane_numerator) / Interval::point(iv_norm(normal).hi()))
        .lo()
        .max(0.0);
    let lower = lower.max(plane_lower);
    let side = if signed.lo() > 0.0 {
        ProximitySide::Positive
    } else if signed.hi() < 0.0 {
        ProximitySide::Negative
    } else {
        ProximitySide::Unresolved
    };
    Ok((
        lower,
        SurfaceWitness {
            surface_point,
            surface_point_error_m: error,
            distance_upper_m: upper,
            normal: Some(Vec3::new(
                components[0] / length,
                components[1] / length,
                components[2] / length,
            )),
            side,
            triangle_index: Some(triangle),
            barycentric: Some(barycentric),
        },
    ))
}

pub(crate) fn segment_clearance(
    point: Point3,
    a: [f64; 3],
    b: [f64; 3],
) -> Result<f64, ProximityError> {
    let delta = sub(b, a);
    let length = dot(delta, delta);
    if !length.is_finite() || length <= 0.0 {
        return Err(ProximityError::InvalidGeometry);
    }
    let t = (dot(sub(point_array(point), a), delta) / length).clamp(0.0, 1.0);
    Ok(convex_bounds(point_array(point), [a, b, a], [1.0 - t, t, 0.0])?.0)
}

/// A chart's safe ball is independent from a geometric contact witness.
pub fn chart_proximity(
    chart: &dyn Chart,
    cx: &Cx<'_>,
    point: Point3,
    budget: &mut ProximityBudget,
) -> Result<SurfaceProximity, ProximityError> {
    budget.charge(cx)?;
    if point_array(point).iter().any(|v| !v.is_finite()) {
        return Err(ProximityError::InvalidInput);
    }
    let claim = chart.trace_step_claim();
    if claim == TraceStepClaim::NoClaim {
        return Err(ProximityError::UnsupportedChart);
    }
    let sample = chart.eval(point, cx);
    let certificate = chart.trace_value_enclosure(point, &sample, cx);
    let validated = validate_trace_sample(&sample, claim, certificate)
        .ok_or(ProximityError::InvalidGeometry)?;
    cx.checkpoint()?;
    let witness = if claim == TraceStepClaim::ExactDistance {
        let normal = sample.gradient.and_then(normalize_gradient);
        Some(SurfaceWitness {
            surface_point: point,
            surface_point_error_m: validated.geometric_hit_upper,
            distance_upper_m: validated.geometric_hit_upper,
            normal,
            side: if certificate.lo > 0.0 {
                ProximitySide::Positive
            } else if certificate.hi < 0.0 {
                ProximitySide::Negative
            } else {
                ProximitySide::Unresolved
            },
            triangle_index: None,
            barycentric: None,
        })
    } else {
        None
    };
    Ok(SurfaceProximity {
        clearance_m: validated.safe_radius,
        witness,
    })
}

/// Prepared exactly coplanar, consistently oriented face patches. Interior
/// triangulation edges disappear; geometric patch boundaries remain barriers.
/// This is a local planar continuation witness, not a watertightness claim.
pub struct MeshPlanarPatches<'mesh> {
    mesh: &'mesh TriMesh,
    triangle_patch: Vec<u32>,
    node_patch: Vec<Option<u32>>,
    boundary_edges: Vec<(u32, u32, u32)>,
}

impl<'mesh> MeshPlanarPatches<'mesh> {
    pub fn try_new(
        mesh: &'mesh TriMesh,
        cx: &Cx<'_>,
        budget: &mut ProximityBudget,
    ) -> Result<Self, ProximityError> {
        validate_mesh_bvh(mesh, cx, budget)?;
        let count = mesh.triangles.len();
        // Indexed meshes commonly duplicate vertices at material/UV seams.
        // Equal finite coordinates identify the same geometric edge even
        // when authoring indices differ. Signed zero has one canonical key.
        let mut coordinate_ids = std::collections::BTreeMap::new();
        let mut canonical = Vec::new();
        canonical
            .try_reserve_exact(mesh.vertices.len())
            .map_err(|_| ProximityError::Allocation)?;
        for (index, vertex) in mesh.vertices.iter().enumerate() {
            budget.charge(cx)?;
            if vertex.iter().any(|v| !v.is_finite()) {
                return Err(ProximityError::InvalidGeometry);
            }
            let key = vertex.map(|v| if v == 0.0 { 0 } else { v.to_bits() });
            let id = u32::try_from(index).map_err(|_| ProximityError::InvalidGeometry)?;
            canonical.push(*coordinate_ids.entry(key).or_insert(id));
        }
        let edge_count = count.checked_mul(3).ok_or(ProximityError::Allocation)?;
        let mut edges = Vec::new();
        edges
            .try_reserve_exact(edge_count)
            .map_err(|_| ProximityError::Allocation)?;
        let mut parent = Vec::new();
        parent
            .try_reserve_exact(count)
            .map_err(|_| ProximityError::Allocation)?;
        for (index, triangle) in mesh.triangles.iter().enumerate() {
            budget.charge(cx)?;
            let id = u32::try_from(index).map_err(|_| ProximityError::InvalidGeometry)?;
            parent.push(id);
            if triangle.iter().any(|&v| v as usize >= mesh.vertices.len()) {
                return Err(ProximityError::InvalidGeometry);
            }
            let vertices = triangle.map(|v| mesh.vertices[v as usize]);
            // Also validates finite, nondegenerate geometry before any query.
            let _ = triangle_proximity(array_point(vertices[0]), vertices, id)?;
            let [a, b, c] = triangle.map(|v| canonical[v as usize]);
            for (a, b) in [(a, b), (b, c), (c, a)] {
                edges.push((a.min(b), a.max(b), id, a < b));
            }
        }
        edges.sort_unstable();
        for group in edges.chunk_by(|a, b| a.0 == b.0 && a.1 == b.1) {
            budget.charge(cx)?;
            if group.len() != 2 || group[0].3 == group[1].3 {
                continue;
            }
            let a = group[0].2;
            let b = group[1].2;
            if coplanar_neighbors(mesh, a, b) {
                let a = patch_root(&mut parent, a);
                let b = patch_root(&mut parent, b);
                parent[a.max(b) as usize] = a.min(b);
            }
        }
        for i in 0..count {
            budget.charge(cx)?;
            parent[i] = patch_root(&mut parent, i as u32);
        }
        let mut boundary_edges = Vec::new();
        boundary_edges
            .try_reserve_exact(edge_count)
            .map_err(|_| ProximityError::Allocation)?;
        for group in edges.chunk_by(|a, b| a.0 == b.0 && a.1 == b.1) {
            budget.charge(cx)?;
            if group.len() == 2
                && group[0].3 != group[1].3
                && parent[group[0].2 as usize] == parent[group[1].2 as usize]
            {
                continue;
            }
            for &(a, b, triangle, _) in group {
                boundary_edges.push((parent[triangle as usize], a, b));
            }
        }
        boundary_edges.sort_unstable();
        boundary_edges.dedup();
        let mut node_patch = Vec::new();
        node_patch
            .try_reserve_exact(mesh.nodes.len())
            .map_err(|_| ProximityError::Allocation)?;
        node_patch.resize(mesh.nodes.len(), None);
        for index in (0..mesh.nodes.len()).rev() {
            budget.charge(cx)?;
            let node = &mesh.nodes[index];
            node_patch[index] = if node.leaf {
                let mut common = None;
                for slot in node.a..node.a + node.b {
                    let patch = parent[mesh.order[slot as usize] as usize];
                    if common.is_some_and(|old| old != patch) {
                        common = None;
                        break;
                    }
                    common = Some(patch);
                }
                common
            } else {
                let a = node_patch[node.a as usize];
                let b = node_patch[node.b as usize];
                if a == b { a } else { None }
            };
        }
        cx.checkpoint()?;
        Ok(Self {
            mesh,
            triangle_patch: parent,
            node_patch,
            boundary_edges,
        })
    }

    #[must_use]
    pub fn patch_id(&self, triangle: u32) -> Option<u32> {
        self.triangle_patch.get(triangle as usize).copied()
    }

    /// Retained BLAS bounds are sound after preparation validates every leaf.
    /// Unreferenced authoring vertices do not enlarge the represented surface.
    /// An empty triangle set with vertices retains a degenerate empty bound;
    /// a completely empty authoring mesh has no instance placement bound.
    pub(crate) fn bounds(&self) -> Option<Aabb> {
        self.mesh.nodes.first().map_or_else(
            || {
                self.mesh.vertices.first().map(|&point| Aabb {
                    min: array_point(point),
                    max: array_point(point),
                })
            },
            |node| {
                Some(Aabb {
                    min: array_point(node.lo),
                    max: array_point(node.hi),
                })
            },
        )
    }

    /// Accelerated query of the immutable mesh validated during preparation.
    pub fn nearest_surface_with_cx(
        &self,
        cx: &Cx<'_>,
        point: Point3,
        budget: &mut ProximityBudget,
    ) -> Result<SurfaceProximity, ProximityError> {
        self.mesh.nearest_surface_filtered(cx, point, budget, None)
    }

    /// Ball around the rounded surface representative containing only this
    /// planar patch. The caller must additionally exclude other scene objects.
    pub(crate) fn departure_radius(
        &self,
        cx: &Cx<'_>,
        witness: &SurfaceWitness,
        budget: &mut ProximityBudget,
    ) -> Result<Option<f64>, ProximityError> {
        let mesh = self.mesh;
        let triangle = witness.triangle_index.ok_or(ProximityError::InvalidInput)?;
        let patch = self
            .patch_id(triangle)
            .ok_or(ProximityError::InvalidInput)?;
        let start = self.boundary_edges.partition_point(|e| e.0 < patch);
        let end = self.boundary_edges.partition_point(|e| e.0 <= patch);
        let mut radius = f64::INFINITY;
        for &(_, a, b) in &self.boundary_edges[start..end] {
            budget.charge(cx)?;
            radius = radius.min(segment_clearance(
                witness.surface_point,
                mesh.vertices[a as usize],
                mesh.vertices[b as usize],
            )?);
        }
        let other =
            mesh.nearest_surface_filtered(cx, witness.surface_point, budget, Some((self, patch)))?;
        radius = radius.min(other.clearance_m);
        // A finite patch and positive open ball are needed to leave contact.
        let radius = conservative_difference_lower(radius, witness.surface_point_error_m);
        Ok((radius.is_finite() && radius > 0.0).then_some(radius))
    }
}

fn patch_root(parent: &mut [u32], mut index: u32) -> u32 {
    while parent[index as usize] != index {
        parent[index as usize] = parent[parent[index as usize] as usize];
        index = parent[index as usize];
    }
    index
}

fn validate_mesh_bvh(
    mesh: &TriMesh,
    cx: &Cx<'_>,
    budget: &mut ProximityBudget,
) -> Result<(), ProximityError> {
    cx.checkpoint()?;
    if mesh.order.len() != mesh.triangles.len()
        || (!mesh.triangles.is_empty() && mesh.nodes.is_empty())
    {
        return Err(ProximityError::InvalidGeometry);
    }
    // The BVH/order are private and immutable. Public authoring arrays can
    // change, so every leaf must still enclose every referenced vertex.
    for node in &mesh.nodes {
        budget.charge(cx)?;
        if !node.leaf {
            continue;
        }
        for slot in node.a..node.a + node.b {
            budget.charge(cx)?;
            let triangle = *mesh
                .order
                .get(slot as usize)
                .ok_or(ProximityError::InvalidGeometry)?;
            let vertices = mesh
                .triangles
                .get(triangle as usize)
                .ok_or(ProximityError::InvalidGeometry)?;
            for &vertex in vertices {
                let p = mesh
                    .vertices
                    .get(vertex as usize)
                    .ok_or(ProximityError::InvalidGeometry)?;
                if (0..3).any(|i| !p[i].is_finite() || p[i] < node.lo[i] || p[i] > node.hi[i]) {
                    return Err(ProximityError::InvalidGeometry);
                }
            }
        }
    }
    Ok(())
}

fn coplanar_neighbors(mesh: &TriMesh, a: u32, b: u32) -> bool {
    let a = mesh.triangles[a as usize].map(|i| mesh.vertices[i as usize]);
    let b = mesh.triangles[b as usize].map(|i| mesh.vertices[i as usize]);
    let face_normal = |v: [[f64; 3]; 3]| {
        let ab = iv_sub(iv(v[1]), iv(v[0]));
        let ac = iv_sub(iv(v[2]), iv(v[0]));
        [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ]
    };
    // Opposite traversal of a shared edge is not enough: overlapping faces
    // can lie on the same side of that edge and have opposite orientations.
    if iv_dot(face_normal(a), face_normal(b)).lo() <= 0.0 {
        return false;
    }
    // Axis-aligned planes are exact directly, including extreme exponents.
    if (0..3).any(|axis| a.iter().chain(b.iter()).all(|v| v[axis] == a[0][axis])) {
        return true;
    }
    // The exact predicate owner's certificate excludes intermediate overflow
    // and underflow. Outside this conservative exponent envelope, retain the
    // triangulation edge instead of asserting a coplanar continuation patch.
    if a.iter()
        .chain(b.iter())
        .flatten()
        .any(|v| *v != 0.0 && !(1.0e-60..=1.0e60).contains(&v.abs()))
    {
        return false;
    }
    for p in a.iter().chain(b.iter()) {
        for q in a.iter().chain(b.iter()) {
            if (0..3).any(|i| {
                let d = (p[i] - q[i]).abs();
                d != 0.0 && !(1.0e-60..=1.0e60).contains(&d)
            }) {
                return false;
            }
        }
    }
    b.into_iter().all(|p| {
        fs_ivl::predicates::orient3d(a[0], a[1], a[2], p) == fs_ivl::predicates::Sign::Zero
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn with_cx(operation: impl FnOnce(&Cx<'_>)) {
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 37,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            operation(&cx);
        });
    }

    #[test]
    fn g1_triangle_distance_enclosures_cover_face_edge_vertex_and_coplanar_exterior() {
        let triangle = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        for (point, expected) in [
            (Point3::new(0.2, 0.2, 2.0), 2.0),
            (Point3::new(1.0, 1.0, 1.0), 1.5_f64.sqrt()),
            (Point3::new(-1.0, 0.0, 0.0), 1.0),
            (Point3::new(0.8, 0.8, 0.0), 0.6 / 2.0_f64.sqrt()),
        ] {
            let (lower, witness) = triangle_proximity(point, triangle, 0).unwrap();
            assert!(
                lower <= expected && witness.distance_upper_m >= expected,
                "{lower} <= {expected} <= {}",
                witness.distance_upper_m
            );
            assert!(witness.distance_upper_m - lower < 1.0e-12);
        }
    }

    #[test]
    fn g1_thin_separated_faces_keep_positive_non_tunneling_clearance() {
        with_cx(|cx| {
            let mesh = TriMesh::new(
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0e-9],
                    [1.0, 0.0, 1.0e-9],
                    [0.0, 1.0, 1.0e-9],
                ],
                vec![[0, 1, 2], [3, 5, 4]],
            );
            let q = mesh
                .nearest_surface_with_cx(
                    cx,
                    Point3::new(0.2, 0.2, 0.4e-9),
                    &mut ProximityBudget::new(100),
                )
                .unwrap();
            assert!(q.clearance_m > 0.39e-9 && q.clearance_m <= 0.4e-9);
            assert!(q.witness.unwrap().distance_upper_m >= 0.4e-9);
        });
    }

    #[test]
    fn g1_large_face_clearance_stays_sharp_at_curved_arrival_tolerance() {
        let triangle = [[-10.0, -10.0, 1.0], [10.0, -10.0, 1.0], [0.0, 10.0, 1.0]];
        for gap in [1.0e-6, 1.0e-8, 1.0e-10] {
            let point = Point3::new(0.083531731, 0.0, 1.0 - gap);
            let expected = 1.0 - point.z;
            let (lower, witness) = triangle_proximity(point, triangle, 0).unwrap();
            assert!(
                lower > expected * 0.999,
                "collapsed gap {expected:e}: {lower:e}"
            );
            assert!(lower <= expected && witness.distance_upper_m >= expected);
        }
    }

    #[test]
    fn g3_nearest_bvh_prunes_distant_triangles_and_retains_canonical_ties() {
        with_cx(|cx| {
            let mut vertices = Vec::new();
            let mut triangles = Vec::new();
            for i in 0..512u32 {
                let x = f64::from(i) * 4.0;
                vertices.extend([[x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.0]]);
                triangles.push([3 * i, 3 * i + 1, 3 * i + 2]);
            }
            let mesh = TriMesh::new(vertices, triangles);
            let prepared =
                MeshPlanarPatches::try_new(&mesh, cx, &mut ProximityBudget::new(100_000)).unwrap();
            let mut budget = ProximityBudget::new(10_000);
            let q = prepared
                .nearest_surface_with_cx(cx, Point3::new(0.25, 0.25, 0.5), &mut budget)
                .unwrap();
            assert_eq!(q.witness.unwrap().triangle_index, Some(0));
            assert!(budget.used() < 128, "pruning work {}", budget.used());
            let duplicated = TriMesh::new(mesh.vertices[..3].to_vec(), vec![[0, 1, 2], [0, 1, 2]]);
            let q = duplicated
                .nearest_surface_with_cx(cx, Point3::new(0.25, 0.25, 0.5), &mut budget)
                .unwrap();
            assert_eq!(q.witness.unwrap().triangle_index, Some(0));
        });
    }

    #[test]
    fn g1_coplanar_patch_crosses_internal_diagonal_but_respects_thin_opposite_face() {
        with_cx(|cx| {
            for opposite in [false, true] {
                let mut vertices = vec![
                    [-1.0, -1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [-1.0, 1.0, 0.0],
                ];
                let mut triangles = vec![[0, 1, 2], [0, 2, 3]];
                if opposite {
                    vertices.extend([
                        [-1.0, -1.0, 1.0e-6],
                        [1.0, -1.0, 1.0e-6],
                        [1.0, 1.0, 1.0e-6],
                        [-1.0, 1.0, 1.0e-6],
                    ]);
                    triangles.extend([[4, 6, 5], [4, 7, 6]]);
                }
                let mesh = TriMesh::new(vertices, triangles);
                let mut budget = ProximityBudget::new(1000);
                let patches = MeshPlanarPatches::try_new(&mesh, cx, &mut budget).unwrap();
                assert_eq!(patches.patch_id(0), patches.patch_id(1));
                let query = patches
                    .nearest_surface_with_cx(cx, Point3::new(0.0, 0.0, -1.0e-8), &mut budget)
                    .unwrap();
                let radius = patches
                    .departure_radius(cx, &query.witness.unwrap(), &mut budget)
                    .unwrap()
                    .unwrap();
                if opposite {
                    assert!(radius > 0.99e-6 && radius <= 1.0e-6);
                } else {
                    assert!(radius > 0.99 && radius <= 1.0);
                }
            }
        });
    }

    #[test]
    fn g1_duplicate_vertex_quad_has_one_geometric_patch() {
        with_cx(|cx| {
            let mesh = TriMesh::new(
                vec![
                    [-1.0, -1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [-1.0, -1.0, -0.0],
                    [1.0, 1.0, -0.0],
                    [-1.0, 1.0, -0.0],
                ],
                vec![[0, 1, 2], [3, 4, 5]],
            );
            let mut budget = ProximityBudget::new(1000);
            let patches = MeshPlanarPatches::try_new(&mesh, cx, &mut budget).unwrap();
            assert_eq!(patches.patch_id(0), patches.patch_id(1));
            let q = patches
                .nearest_surface_with_cx(cx, Point3::new(0.0, 0.0, 1.0e-9), &mut budget)
                .unwrap();
            let radius = patches
                .departure_radius(cx, &q.witness.unwrap(), &mut budget)
                .unwrap()
                .unwrap();
            assert!(radius > 0.99 && radius <= 1.0);
        });
    }

    #[test]
    fn g1_oppositely_oriented_overlapping_faces_are_not_one_patch() {
        with_cx(|cx| {
            let mesh = TriMesh::new(
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.5, 0.5, 0.0],
                ],
                vec![[0, 1, 2], [1, 0, 3]],
            );
            let mut budget = ProximityBudget::new(1000);
            let patches = MeshPlanarPatches::try_new(&mesh, cx, &mut budget).unwrap();
            assert_ne!(patches.patch_id(0), patches.patch_id(1));
            let q = patches
                .nearest_surface_with_cx(cx, Point3::new(0.4, 0.1, 1.0e-8), &mut budget)
                .unwrap();
            assert_eq!(
                patches
                    .departure_radius(cx, &q.witness.unwrap(), &mut budget)
                    .unwrap(),
                None
            );
        });
    }

    #[test]
    fn g4_mutated_authoring_geometry_cannot_reuse_stale_bvh_clearance() {
        with_cx(|cx| {
            let mut mesh = TriMesh::new(
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                vec![[0, 1, 2]],
            );
            mesh.vertices[0][2] = 1.0;
            assert!(matches!(
                MeshPlanarPatches::try_new(&mesh, cx, &mut ProximityBudget::new(100)),
                Err(ProximityError::InvalidGeometry)
            ));
            assert_eq!(
                mesh.nearest_surface_with_cx(
                    cx,
                    Point3::new(0.0, 0.0, 1.1),
                    &mut ProximityBudget::new(100)
                ),
                Err(ProximityError::InvalidGeometry)
            );
        });
    }

    #[test]
    fn g0_distance_chart_has_upper_proximity_but_implicit_residual_does_not() {
        struct Implicit;
        impl Chart for Implicit {
            fn eval(&self, p: Point3, _: &Cx<'_>) -> ChartSample {
                ChartSample {
                    signed_distance: 2.0 * p.z,
                    gradient: Some(Vec3::new(0.0, 0.0, 2.0)),
                    lipschitz: Some(2.0),
                    error: NumericalCertificate::no_claim(),
                }
            }
            fn support(&self) -> Aabb {
                Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
            }
            fn trace_step_claim(&self) -> TraceStepClaim {
                TraceStepClaim::LipschitzImplicit
            }
            fn trace_value_enclosure(
                &self,
                _: Point3,
                s: &ChartSample,
                _: &Cx<'_>,
            ) -> NumericalCertificate {
                NumericalCertificate::exact(s.signed_distance)
            }
            fn name(&self) -> &'static str {
                "proximity-implicit-test"
            }
        }
        with_cx(|cx| {
            let mut budget = ProximityBudget::new(20);
            let q =
                chart_proximity(&Implicit, cx, Point3::new(0.0, 0.0, 0.1), &mut budget).unwrap();
            assert!(q.clearance_m > 0.099 && q.clearance_m <= 0.1);
            assert!(q.witness.is_none());
            let sphere = fs_geom::fixtures::SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            };
            let q = chart_proximity(&sphere, cx, Point3::new(0.0, 0.0, 1.1), &mut budget).unwrap();
            assert!(q.clearance_m <= 0.1 + 1.0e-15);
            assert!(q.witness.unwrap().distance_upper_m >= 0.1);
        });
    }

    #[test]
    fn g4_native_query_refuses_exhausted_work_budget() {
        with_cx(|cx| {
            let mesh = TriMesh::new(
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                vec![[0, 1, 2]],
            );
            assert_eq!(
                mesh.nearest_surface_with_cx(
                    cx,
                    Point3::new(0.0, 0.0, 1.0),
                    &mut ProximityBudget::new(0)
                ),
                Err(ProximityError::WorkBudgetExceeded)
            );
        });
    }
}
