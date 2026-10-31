//! Prepared native boundary geometry for curved optical paths.
//!
//! Interval affine maps account for stored quaternion normalization, inverse
//! evaluation and world-point rounding. Local Euclidean bounds are converted
//! using a singular-value enclosure, not an assumed exact binary64 isometry.
//! The prepared top-level BVH and immutable per-mesh BLAS/planar patches are
//! shared across all queries; one cumulative work budget covers the full path.

use super::*;
pub use crate::charts::proximity::ProximityError;
use crate::charts::proximity::{
    MeshPlanarPatches, ProximityBudget, ProximitySide, SurfaceProximity, SurfaceWitness,
    array_point, box_clearance, chart_normal_cone, chart_proximity, chart_signed_distance,
    point_array,
};
use fs_geom::{Aabb, TraceNormalCone};
use fs_ivl::Interval;
use std::collections::BTreeMap;

/// Bit-affecting native clearance and planar departure semantics.
pub const BOUNDARY_GEOMETRY_SEMANTICS_VERSION: u32 = 2;

type IVec = [Interval; 3];
type IMat = [[Interval; 3]; 3];
fn ivec(p: Point3) -> IVec {
    point_array(p).map(Interval::point)
}
fn isub(a: IVec, b: IVec) -> IVec {
    std::array::from_fn(|i| a[i] - b[i])
}
fn iadd(a: IVec, b: IVec) -> IVec {
    std::array::from_fn(|i| a[i] + b[i])
}
fn idot(a: IVec, b: IVec) -> Interval {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn imul(m: IMat, p: IVec) -> IVec {
    m.map(|row| idot(row, p))
}
fn inorm(v: IVec) -> Interval {
    let v = v.map(Interval::abs);
    idot(v, v).sqrt()
}
fn finite(v: IVec) -> bool {
    v.iter().all(|v| v.lo().is_finite() && v.hi().is_finite())
}
fn uncertainty(v: IVec, p: Point3) -> f64 {
    inorm(isub(v, ivec(p))).hi()
}
fn lower_difference(a: f64, b: f64) -> f64 {
    (Interval::point(a) - Interval::point(b)).lo().max(0.0)
}

#[derive(Clone, Copy)]
struct PreparedPose {
    matrix: IMat,
    inverse: IMat,
    translation: IVec,
    min_scale: f64,
    max_scale: f64,
}
impl PreparedPose {
    fn new(transform: RigidTransform) -> Result<Self, ProximityError> {
        let [x, y, z, w] = transform.rotation_xyzw().map(Interval::point);
        let one = Interval::point(1.0);
        let two = Interval::point(2.0);
        let m = [
            [
                one - two * (y * y + z * z),
                two * (x * y - z * w),
                two * (x * z + y * w),
            ],
            [
                two * (x * y + z * w),
                one - two * (x * x + z * z),
                two * (y * z - x * w),
            ],
            [
                two * (x * z - y * w),
                two * (y * z + x * w),
                one - two * (x * x + y * y),
            ],
        ];
        let cof = [
            [
                m[1][1] * m[2][2] - m[1][2] * m[2][1],
                m[1][2] * m[2][0] - m[1][0] * m[2][2],
                m[1][0] * m[2][1] - m[1][1] * m[2][0],
            ],
            [
                m[0][2] * m[2][1] - m[0][1] * m[2][2],
                m[0][0] * m[2][2] - m[0][2] * m[2][0],
                m[0][1] * m[2][0] - m[0][0] * m[2][1],
            ],
            [
                m[0][1] * m[1][2] - m[0][2] * m[1][1],
                m[0][2] * m[1][0] - m[0][0] * m[1][2],
                m[0][0] * m[1][1] - m[0][1] * m[1][0],
            ],
        ];
        let determinant = idot(m[0], cof[0]);
        if determinant.lo() <= 0.0 || !determinant.hi().is_finite() {
            return Err(ProximityError::NumericalRange);
        }
        let inverse: IMat =
            std::array::from_fn(|i| std::array::from_fn(|j| cof[j][i] / determinant));
        let gram: IMat = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                idot(
                    std::array::from_fn(|k| m[k][i]),
                    std::array::from_fn(|k| m[k][j]),
                )
            })
        });
        let mut lower = f64::INFINITY;
        let mut upper = 0.0_f64;
        for i in 0..3 {
            let radius = (0..3)
                .filter(|&j| j != i)
                .fold(Interval::point(0.0), |sum, j| sum + gram[i][j].abs());
            lower = lower.min((gram[i][i] - radius).lo());
            upper = upper.max((gram[i][i] + radius).hi());
        }
        if !lower.is_finite() || lower <= 0.0 || !upper.is_finite() {
            return Err(ProximityError::NumericalRange);
        }
        let min_scale = Interval::point(lower).sqrt().lo();
        let max_scale = Interval::point(upper).sqrt().hi();
        Ok(Self {
            matrix: m,
            inverse,
            translation: transform.translation_m().map(Interval::point),
            min_scale,
            max_scale,
        })
    }
    fn local_point(self, point: Point3) -> Result<(Point3, f64), ProximityError> {
        if point_array(point).iter().any(|v| !v.is_finite()) {
            return Err(ProximityError::InvalidInput);
        }
        let interval = imul(self.inverse, isub(ivec(point), self.translation));
        if !finite(interval) {
            return Err(ProximityError::NumericalRange);
        }
        let center = array_point(interval.map(Interval::midpoint));
        let error = uncertainty(interval, center);
        if !error.is_finite() {
            return Err(ProximityError::NumericalRange);
        }
        Ok((center, error))
    }
    fn world_point(self, point: Point3) -> Result<(Point3, f64), ProximityError> {
        let interval = iadd(imul(self.matrix, ivec(point)), self.translation);
        if !finite(interval) {
            return Err(ProximityError::NumericalRange);
        }
        let center = array_point(interval.map(Interval::midpoint));
        Ok((center, uncertainty(interval, center)))
    }
    fn world_bounds(self, local: Aabb) -> Result<Aabb, ProximityError> {
        let interval = [
            Interval::new(local.min.x, local.max.x),
            Interval::new(local.min.y, local.max.y),
            Interval::new(local.min.z, local.max.z),
        ];
        let world = iadd(imul(self.matrix, interval), self.translation);
        if !finite(world) {
            return Err(ProximityError::NumericalRange);
        }
        Ok(Aabb {
            min: array_point(world.map(Interval::lo)),
            max: array_point(world.map(Interval::hi)),
        })
    }
    fn world_proximity(
        self,
        local: SurfaceProximity,
        query_error: f64,
    ) -> Result<SurfaceProximity, ProximityError> {
        let clearance = if local.clearance_m.is_infinite() {
            f64::INFINITY
        } else {
            (Interval::point(lower_difference(local.clearance_m, query_error))
                * Interval::point(self.min_scale))
            .lo()
            .max(0.0)
        };
        let witness = if let Some(mut w) = local.witness {
            let (point, point_error) = self.world_point(w.surface_point)?;
            w.surface_point = point;
            w.surface_point_error_m = (Interval::point(w.surface_point_error_m)
                * Interval::point(self.max_scale)
                + Interval::point(point_error))
            .hi();
            w.distance_upper_m = ((Interval::point(w.distance_upper_m)
                + Interval::point(query_error))
                * Interval::point(self.max_scale))
            .hi();
            w.normal = w.normal.and_then(|normal| {
                let inverse_transpose =
                    std::array::from_fn(|i| std::array::from_fn(|j| self.inverse[j][i]));
                let n = imul(
                    inverse_transpose,
                    ivec(Point3::new(normal.x, normal.y, normal.z)),
                )
                .map(Interval::midpoint);
                let v = Vec3::new(n[0], n[1], n[2]);
                let length = v.norm();
                (length.is_finite() && length > 0.0).then(|| v.scale(1.0 / length))
            });
            if local.clearance_m <= query_error {
                w.side = ProximitySide::Unresolved;
            }
            if !w.surface_point_error_m.is_finite() || !w.distance_upper_m.is_finite() {
                return Err(ProximityError::NumericalRange);
            }
            Some(w)
        } else {
            None
        };
        Ok(SurfaceProximity {
            clearance_m: clearance,
            witness,
        })
    }

    fn world_normal_cone(self, local: TraceNormalCone) -> Result<TraceNormalCone, ProximityError> {
        let inverse_transpose =
            std::array::from_fn(|i| std::array::from_fn(|j| self.inverse[j][i]));
        let transformed = imul(
            inverse_transpose,
            ivec(Point3::new(local.axis.x, local.axis.y, local.axis.z)),
        );
        let norm = inorm(transformed);
        if norm.lo() <= 0.0 || !norm.hi().is_finite() {
            return Err(ProximityError::NumericalRange);
        }
        let unit = transformed.map(|v| v / norm);
        let n = unit.map(Interval::midpoint);
        let axis = Vec3::new(n[0], n[1], n[2]);
        let length = axis.norm();
        if !length.is_finite() || length <= 0.0 {
            return Err(ProximityError::NumericalRange);
        }
        let axis = axis.scale(1.0 / length);
        let axis_error = uncertainty(unit, Point3::new(axis.x, axis.y, axis.z));
        // Normalization is 2/min(|A n|)-Lipschitz for these nonzero vectors;
        // A is the inverse transpose, with singular values in
        // [1/max_scale,1/min_scale]. Include its interval arithmetic too.
        let deviation = (Interval::point(axis_error)
            + Interval::point(2.0) * Interval::point(self.max_scale)
                / Interval::point(self.min_scale)
                * Interval::point(local.max_deviation))
        .hi();
        if !deviation.is_finite() {
            return Err(ProximityError::NumericalRange);
        }
        Ok(TraceNormalCone {
            axis,
            max_deviation: deviation,
        })
    }

    fn world_signed_distance(self, local: Interval, query_error: f64) -> Interval {
        let local = local + Interval::new(-query_error, query_error);
        let lower_scale = if local.lo() >= 0.0 {
            self.min_scale
        } else {
            self.max_scale
        };
        let upper_scale = if local.hi() <= 0.0 {
            self.min_scale
        } else {
            self.max_scale
        };
        Interval::new(
            (Interval::point(local.lo()) * Interval::point(lower_scale)).lo(),
            (Interval::point(local.hi()) * Interval::point(upper_scale)).hi(),
        )
    }
}

/// Separate numerical point-query witness carrying native frame identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InstanceProximity {
    pub object_id: u64,
    pub geometry_identity: ContentHash,
    pub frame_identity: ContentHash,
    pub proximity: SurfaceProximity,
    pub local_witness: Option<SurfaceWitness>,
}

fn query_instance(
    instance: &GeometryInstance,
    pose: PreparedPose,
    patches: Option<&MeshPlanarPatches<'_>>,
    cx: &Cx<'_>,
    point: Point3,
    budget: &mut ProximityBudget,
) -> Result<InstanceProximity, ProximityError> {
    budget.charge(cx)?;
    let (local, error) = pose.local_point(point)?;
    let proximity = match instance.geometry() {
        SharedGeometry::Mesh(mesh) => {
            if let Some(patches) = patches {
                patches.nearest_surface_with_cx(cx, local, budget)?
            } else {
                mesh.nearest_surface_with_cx(cx, local, budget)?
            }
        }
        SharedGeometry::Chart(chart) => chart_proximity(chart.as_ref(), cx, local, budget)?,
    };
    Ok(InstanceProximity {
        object_id: instance.object_id(),
        geometry_identity: instance.geometry_identity(),
        frame_identity: instance.frame_identity(),
        local_witness: proximity.witness,
        proximity: pose.world_proximity(proximity, error)?,
    })
}

impl GeometryInstance {
    /// Direct convenience query. A [`BoundaryGeometry`] caches the interval
    /// placement for many queries and also prunes unrelated objects.
    pub fn nearest_surface_with_cx(
        &self,
        cx: &Cx<'_>,
        point: Point3,
        budget: &mut ProximityBudget,
    ) -> Result<InstanceProximity, ProximityError> {
        query_instance(
            self,
            PreparedPose::new(self.transform())?,
            None,
            cx,
            point,
            budget,
        )
    }
}

struct PreparedInstance<'scene> {
    pose: PreparedPose,
    bounds: Aabb,
    patches: Option<Arc<MeshPlanarPatches<'scene>>>,
}
struct Node {
    bounds: Aabb,
    left: u32,
    right: u32,
    start: u32,
    count: u32,
}

/// Immutable, borrowed rigid scene with prepared interval transforms, a
/// deterministic top-level BVH and shared native planar patch metadata.
pub struct BoundaryGeometry<'scene> {
    instances: &'scene [GeometryInstance],
    prepared: Vec<PreparedInstance<'scene>>,
    nodes: Vec<Node>,
    order: Vec<u32>,
}

fn checked_bounds(bounds: Aabb) -> Result<Aabb, ProximityError> {
    if point_array(bounds.min)
        .iter()
        .chain(point_array(bounds.max).iter())
        .any(|v| !v.is_finite())
        || (0..3).any(|i| point_array(bounds.min)[i] > point_array(bounds.max)[i])
    {
        return Err(ProximityError::InvalidGeometry);
    }
    Ok(bounds)
}

impl<'scene> BoundaryGeometry<'scene> {
    pub fn try_new(
        cx: &Cx<'_>,
        instances: &'scene [GeometryInstance],
        maximum_build_work: u64,
    ) -> Result<Self, ProximityError> {
        let mut budget = ProximityBudget::new(maximum_build_work);
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(instances.len())
            .map_err(|_| ProximityError::Allocation)?;
        let mut patch_cache: BTreeMap<usize, Arc<MeshPlanarPatches<'scene>>> = BTreeMap::new();
        for instance in instances {
            budget.charge(cx)?;
            let pose = PreparedPose::new(instance.transform())?;
            let patches = if let SharedGeometry::Mesh(mesh) = instance.geometry() {
                let key = Arc::as_ptr(mesh) as usize;
                if let Some(patches) = patch_cache.get(&key) {
                    Some(Arc::clone(patches))
                } else {
                    let patches = Arc::new(MeshPlanarPatches::try_new(mesh, cx, &mut budget)?);
                    patch_cache.insert(key, Arc::clone(&patches));
                    Some(patches)
                }
            } else {
                None
            };
            let local = match (instance.geometry(), patches.as_deref()) {
                (SharedGeometry::Mesh(_), Some(patches)) => {
                    patches.bounds().ok_or(ProximityError::InvalidGeometry)?
                }
                (SharedGeometry::Chart(chart), _) => chart.support(),
                _ => return Err(ProximityError::InvalidGeometry),
            };
            let bounds = pose.world_bounds(checked_bounds(local)?)?;
            prepared.push(PreparedInstance {
                pose,
                bounds,
                patches,
            });
        }
        let count = u32::try_from(instances.len()).map_err(|_| ProximityError::InvalidInput)?;
        let mut order = Vec::new();
        order
            .try_reserve_exact(instances.len())
            .map_err(|_| ProximityError::Allocation)?;
        order.extend(0..count);
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(
                instances
                    .len()
                    .checked_mul(2)
                    .ok_or(ProximityError::Allocation)?,
            )
            .map_err(|_| ProximityError::Allocation)?;
        if count > 0 {
            build_node(
                cx,
                &mut budget,
                instances,
                &prepared,
                &mut order,
                &mut nodes,
                0,
                count,
            )?;
        }
        cx.checkpoint()?;
        Ok(Self {
            instances,
            prepared,
            nodes,
            order,
        })
    }

    #[must_use]
    pub const fn instances(&self) -> &'scene [GeometryInstance] {
        self.instances
    }

    fn nearest(
        &self,
        cx: &Cx<'_>,
        point: Point3,
        budget: &mut ProximityBudget,
        excluded: Option<usize>,
    ) -> Result<(f64, Option<(usize, InstanceProximity)>), ProximityError> {
        cx.checkpoint()?;
        if point_array(point).iter().any(|v| !v.is_finite()) {
            return Err(ProximityError::InvalidInput);
        }
        let mut clearance = f64::INFINITY;
        let mut best: Option<(usize, InstanceProximity)> = None;
        let mut upper = f64::INFINITY;
        if self.nodes.is_empty() {
            return Ok((clearance, best));
        }
        // Balanced median construction bounds depth by 32 for u32 inputs.
        let mut stack = [0u32; 64];
        let mut len = 1;
        while len > 0 {
            len -= 1;
            let index = stack[len] as usize;
            budget.charge(cx)?;
            let node = &self.nodes[index];
            let lower = box_clearance(point, node.bounds)?;
            if lower > upper {
                continue;
            }
            if node.count > 0 {
                for slot in node.start..node.start + node.count {
                    let index = self.order[slot as usize] as usize;
                    if excluded == Some(index) {
                        continue;
                    }
                    let result = query_instance(
                        &self.instances[index],
                        self.prepared[index].pose,
                        self.prepared[index].patches.as_deref(),
                        cx,
                        point,
                        budget,
                    )?;
                    clearance = clearance.min(result.proximity.clearance_m);
                    if let Some(witness) = result.proximity.witness {
                        let wins = witness.distance_upper_m < upper
                            || (witness.distance_upper_m == upper
                                && best.is_none_or(|(i, _)| {
                                    (self.instances[index].object_id(), index)
                                        < (self.instances[i].object_id(), i)
                                }));
                        if wins {
                            upper = witness.distance_upper_m;
                            best = Some((index, result));
                        }
                    }
                }
            } else {
                let dl = box_clearance(point, self.nodes[node.left as usize].bounds)?;
                let dr = box_clearance(point, self.nodes[node.right as usize].bounds)?;
                let (near, far) = if dl <= dr {
                    (node.left, node.right)
                } else {
                    (node.right, node.left)
                };
                if len + 2 > stack.len() {
                    return Err(ProximityError::InvalidGeometry);
                }
                stack[len] = far;
                stack[len + 1] = near;
                len += 2;
            }
        }
        Ok((clearance, best))
    }

    #[cfg(feature = "tracer")]
    #[must_use]
    pub fn oracle<'query, 'arena>(
        &'query self,
        cx: &'query Cx<'arena>,
        maximum_query_work: u64,
        contact_tolerance_m: f64,
    ) -> NativeBoundaryOracle<'query, 'scene, 'arena> {
        NativeBoundaryOracle {
            geometry: self,
            cx,
            budget: ProximityBudget::new(maximum_query_work),
            contact_tolerance_m,
            last_error: None,
            last_contact: None,
            smooth_departure: None,
        }
    }
}

fn union(a: Aabb, b: Aabb) -> Aabb {
    Aabb {
        min: array_point(std::array::from_fn(|i| {
            point_array(a.min)[i].min(point_array(b.min)[i])
        })),
        max: array_point(std::array::from_fn(|i| {
            point_array(a.max)[i].max(point_array(b.max)[i])
        })),
    }
}
#[allow(clippy::too_many_arguments)]
fn build_node(
    cx: &Cx<'_>,
    budget: &mut ProximityBudget,
    instances: &[GeometryInstance],
    prepared: &[PreparedInstance<'_>],
    order: &mut [u32],
    nodes: &mut Vec<Node>,
    start: u32,
    count: u32,
) -> Result<u32, ProximityError> {
    budget.charge(cx)?;
    let mut bounds = prepared[order[start as usize] as usize].bounds;
    for &index in &order[start as usize..(start + count) as usize] {
        bounds = union(bounds, prepared[index as usize].bounds);
    }
    let id = u32::try_from(nodes.len()).map_err(|_| ProximityError::InvalidInput)?;
    nodes.push(Node {
        bounds,
        left: 0,
        right: 0,
        start,
        count,
    });
    if count <= 2 {
        return Ok(id);
    }
    let axis = (0..3)
        .max_by(|&a, &b| {
            (point_array(bounds.max)[a] - point_array(bounds.min)[a])
                .total_cmp(&(point_array(bounds.max)[b] - point_array(bounds.min)[b]))
        })
        .unwrap();
    order[start as usize..(start + count) as usize].sort_by(|&a, &b| {
        let aa = prepared[a as usize].bounds;
        let bb = prepared[b as usize].bounds;
        f64::midpoint(point_array(aa.min)[axis], point_array(aa.max)[axis])
            .total_cmp(&f64::midpoint(
                point_array(bb.min)[axis],
                point_array(bb.max)[axis],
            ))
            .then(
                instances[a as usize]
                    .object_id()
                    .cmp(&instances[b as usize].object_id()),
            )
            .then(a.cmp(&b))
    });
    let left_count = count / 2;
    let left = build_node(
        cx, budget, instances, prepared, order, nodes, start, left_count,
    )?;
    let right = build_node(
        cx,
        budget,
        instances,
        prepared,
        order,
        nodes,
        start + left_count,
        count - left_count,
    )?;
    nodes[id as usize] = Node {
        bounds,
        left,
        right,
        start: 0,
        count: 0,
    };
    Ok(id)
}

/// One path's cumulative native query budget and typed diagnostics.
#[cfg(feature = "tracer")]
pub struct NativeBoundaryOracle<'query, 'scene, 'arena> {
    geometry: &'query BoundaryGeometry<'scene>,
    cx: &'query Cx<'arena>,
    budget: ProximityBudget,
    contact_tolerance_m: f64,
    last_error: Option<ProximityError>,
    last_contact: Option<(crate::gradient_index::event::BoundaryWitness, f64)>,
    smooth_departure: Option<(crate::gradient_index::event::BoundaryWitness, f64)>,
}

#[cfg(feature = "tracer")]
impl NativeBoundaryOracle<'_, '_, '_> {
    #[must_use]
    pub const fn work_used(&self) -> u64 {
        self.budget.used()
    }
    #[must_use]
    pub const fn last_error(&self) -> Option<ProximityError> {
        self.last_error
    }

    fn query_native(
        &mut self,
        point: Point3,
    ) -> Result<crate::gradient_index::event::BoundaryQuery, ProximityError> {
        use crate::gradient_index::event::{BoundaryQuery, BoundarySide, BoundaryWitness};
        if !self.contact_tolerance_m.is_finite() || self.contact_tolerance_m <= 0.0 {
            return Err(ProximityError::InvalidInput);
        }
        let (clearance, best) = self
            .geometry
            .nearest(self.cx, point, &mut self.budget, None)?;
        let Some((index, proximity)) = best else {
            return Ok(BoundaryQuery {
                clearance_m: clearance,
                witness: None,
            });
        };
        let mut world = proximity
            .proximity
            .witness
            .ok_or(ProximityError::InvalidGeometry)?;
        if world.distance_upper_m > self.contact_tolerance_m {
            return Ok(BoundaryQuery {
                clearance_m: clearance,
                witness: None,
            });
        }
        // A small local residual alone does not resolve a world-space
        // representative at a large translation. Include its interval
        // rounding envelope before claiming a contact/departure plane.
        let representative_distance = inorm(isub(ivec(point), ivec(world.surface_point))).hi();
        let representative_collar = (Interval::point(representative_distance)
            + Interval::point(world.surface_point_error_m))
        .hi();
        if representative_collar > self.contact_tolerance_m {
            return Err(ProximityError::AmbiguousContact);
        }
        world.distance_upper_m = world.distance_upper_m.max(representative_collar);
        let Some(normal) = world.normal else {
            return Ok(BoundaryQuery {
                clearance_m: clearance,
                witness: None,
            });
        };
        let prepared = &self.geometry.prepared[index];
        let (feature_id, mut radius) = if let (Some(patches), Some(local)) =
            (prepared.patches.as_ref(), proximity.local_witness)
        {
            let triangle = local
                .triangle_index
                .ok_or(ProximityError::InvalidGeometry)?;
            let patch = patches
                .patch_id(triangle)
                .ok_or(ProximityError::InvalidGeometry)?;
            let radius = patches.departure_radius(self.cx, &local, &mut self.budget)?;
            let radius = radius.map(|r| {
                lower_difference(
                    (Interval::point(r) * Interval::point(prepared.pose.min_scale)).lo(),
                    world.surface_point_error_m,
                )
            });
            // A normal is not a unique surface event at an unresolved edge,
            // vertex or competing noncoplanar patch. This also prevents an
            // outside finite-detector near miss becoming a numerical hit.
            if radius.is_none_or(|r| r <= 4.0 * self.contact_tolerance_m) {
                return Err(ProximityError::AmbiguousContact);
            }
            (patch as usize, radius)
        } else {
            (0, None)
        };
        let (other, _) =
            self.geometry
                .nearest(self.cx, world.surface_point, &mut self.budget, Some(index))?;
        let contact_collar = (Interval::point(world.surface_point_error_m)
            + Interval::point(2.0) * Interval::point(self.contact_tolerance_m))
        .hi();
        if other <= contact_collar {
            return Err(ProximityError::AmbiguousContact);
        }
        if let Some(value) = radius {
            let bounded = value.min(if other.is_infinite() {
                other
            } else {
                lower_difference(other, world.surface_point_error_m)
            });
            radius = (bounded.is_finite() && bounded > 0.0).then_some(bounded);
        }
        let side = match world.side {
            ProximitySide::Positive => BoundarySide::Positive,
            ProximitySide::Negative => BoundarySide::Negative,
            ProximitySide::Unresolved => BoundarySide::OnSurface,
        };
        let witness = BoundaryWitness {
            boundary_id: index,
            feature_id,
            surface_point: world.surface_point,
            normal,
            distance_upper_m: world.distance_upper_m,
            side,
            departure_radius_m: radius,
        };
        self.last_contact = Some((witness, other));
        Ok(BoundaryQuery {
            clearance_m: clearance,
            witness: Some(witness),
        })
    }

    fn smooth_neighborhood(
        &mut self,
        request: crate::gradient_index::event::DepartureRequest,
    ) -> Result<Option<crate::gradient_index::event::DepartureNeighborhood>, ProximityError> {
        use crate::gradient_index::event::{BoundarySide, DepartureNeighborhood, cone_margin};
        self.budget.charge(self.cx)?;
        self.smooth_departure = None;
        let witness = request.contact.witness;
        let Some((admitted, other_clearance)) = self.last_contact else {
            return Err(ProximityError::UnsupportedDeparture);
        };
        if admitted != witness
            || request.surface_tolerance_m != self.contact_tolerance_m
            || point_array(request.origin).iter().any(|v| !v.is_finite())
            || [
                request.direction.x,
                request.direction.y,
                request.direction.z,
            ]
            .iter()
            .any(|v| !v.is_finite())
            || (request.direction.dot(request.direction) - 1.0).abs() > 2e-10
            || !request.min_approach_cosine.is_finite()
            || !(0.0..1.0).contains(&request.min_approach_cosine)
            || request.min_approach_cosine <= 0.0
        {
            return Err(ProximityError::InvalidInput);
        }
        let sign = match request.contact.side {
            BoundarySide::Positive => 1.0,
            BoundarySide::Negative => -1.0,
            BoundarySide::OnSurface => return Err(ProximityError::UnsupportedDeparture),
        };
        let index = witness.boundary_id;
        let object = self
            .geometry
            .instances
            .get(index)
            .ok_or(ProximityError::InvalidInput)?;
        let SharedGeometry::Chart(chart) = object.geometry() else {
            return Ok(None);
        };
        if chart.trace_step_claim() != fs_geom::TraceStepClaim::ExactDistance {
            return Ok(None);
        }
        let pose = self.geometry.prepared[index].pose;
        let (local_center, center_error) = pose.local_point(witness.surface_point)?;
        let proximity = inorm(isub(ivec(request.origin), ivec(witness.surface_point))).hi();
        if proximity > request.surface_tolerance_m {
            return Err(ProximityError::AmbiguousContact);
        }
        let approach = cone_margin(witness.normal, 0.0, sign, request.direction);
        if approach <= 0.0 {
            return Err(ProximityError::AmbiguousContact);
        }
        // Start with enough room for the complete numerical collar. Shrinking
        // a normal cone is bounded; an unresolved grazing contact cannot gain
        // acceptance by shrinking below the distance needed to leave it.
        let mut radius =
            (64.0 * request.surface_tolerance_m / approach).min(if other_clearance.is_finite() {
                (0.875 * other_clearance).next_down()
            } else {
                f64::INFINITY
            });
        if !radius.is_finite() || radius <= 0.0 {
            return Err(ProximityError::NumericalRange);
        }
        let mut supports_neighborhood = false;
        for _ in 0..32 {
            let local_radius = (Interval::point(radius) / Interval::point(pose.min_scale)
                + Interval::point(center_error))
            .hi();
            if let Some(local_cone) = chart_normal_cone(
                chart.as_ref(),
                self.cx,
                local_center,
                local_radius,
                &mut self.budget,
            )? {
                supports_neighborhood = true;
                let cone = pose.world_normal_cone(local_cone)?;
                let margin = cone_margin(cone.axis, cone.max_deviation, sign, request.direction);
                if cone.max_deviation < 1.0 && margin >= request.min_approach_cosine {
                    let required = (Interval::point(proximity)
                        + Interval::point(8.0) * Interval::point(request.surface_tolerance_m)
                            / Interval::point(margin))
                    .hi();
                    if radius > required {
                        self.smooth_departure = Some((witness, radius));
                        return Ok(Some(DepartureNeighborhood::Smooth {
                            radius_m: radius,
                            normal_axis: cone.axis,
                            max_normal_deviation: cone.max_deviation,
                        }));
                    }
                }
            }
            if radius <= 4.0 * request.surface_tolerance_m {
                break;
            }
            radius = (0.5 * radius).next_down();
        }
        if supports_neighborhood {
            Err(ProximityError::AmbiguousContact)
        } else {
            Ok(None)
        }
    }

    fn signed_departure_native(
        &mut self,
        witness: crate::gradient_index::event::BoundaryWitness,
        point: Point3,
    ) -> Result<crate::gradient_index::event::SignedBoundaryDistance, ProximityError> {
        use crate::gradient_index::event::SignedBoundaryDistance;
        self.budget.charge(self.cx)?;
        let (admitted, radius) = self
            .smooth_departure
            .ok_or(ProximityError::UnsupportedDeparture)?;
        if admitted != witness || point_array(point).iter().any(|v| !v.is_finite()) {
            return Err(ProximityError::InvalidInput);
        }
        if inorm(isub(ivec(point), ivec(witness.surface_point))).hi() > radius {
            return Err(ProximityError::AmbiguousContact);
        }
        let instance = &self.geometry.instances[witness.boundary_id];
        let SharedGeometry::Chart(chart) = instance.geometry() else {
            return Err(ProximityError::UnsupportedDeparture);
        };
        let pose = self.geometry.prepared[witness.boundary_id].pose;
        let (local, error) = pose.local_point(point)?;
        let signed = chart_signed_distance(chart.as_ref(), self.cx, local, &mut self.budget)?;
        let world = pose.world_signed_distance(signed, error);
        if !world.lo().is_finite() || !world.hi().is_finite() {
            return Err(ProximityError::NumericalRange);
        }
        Ok(SignedBoundaryDistance {
            lower_m: world.lo(),
            upper_m: world.hi(),
        })
    }

    fn map_error(&mut self, error: ProximityError) -> crate::gradient_index::GradientIndexError {
        use crate::gradient_index::GradientIndexError;
        self.last_error = Some(error);
        match error {
            ProximityError::Cancelled => GradientIndexError::Cancelled,
            ProximityError::WorkBudgetExceeded => GradientIndexError::BoundaryQueryBudgetExceeded,
            ProximityError::AmbiguousContact => GradientIndexError::BoundaryResolutionLimit,
            ProximityError::UnsupportedDeparture => {
                GradientIndexError::UnsupportedBoundaryDeparture
            }
            _ => GradientIndexError::InvalidBoundaryQuery,
        }
    }
}

#[cfg(feature = "tracer")]
impl crate::gradient_index::event::BoundaryOracle for NativeBoundaryOracle<'_, '_, '_> {
    fn query(
        &mut self,
        point: Point3,
    ) -> Result<
        crate::gradient_index::event::BoundaryQuery,
        crate::gradient_index::GradientIndexError,
    > {
        self.query_native(point)
            .map_err(|error| self.map_error(error))
    }

    fn prepare_departure(
        &mut self,
        request: crate::gradient_index::event::DepartureRequest,
    ) -> Result<
        crate::gradient_index::event::DepartureNeighborhood,
        crate::gradient_index::GradientIndexError,
    > {
        if request.contact.witness.departure_radius_m.is_some() {
            self.budget
                .charge(self.cx)
                .map_err(|error| self.map_error(error))?;
            return crate::gradient_index::event::DepartureNeighborhood::planar(
                request.contact.witness,
            );
        }
        self.smooth_neighborhood(request)
            .and_then(|neighborhood| neighborhood.ok_or(ProximityError::UnsupportedDeparture))
            .map_err(|error| self.map_error(error))
    }

    fn admit_arrival(
        &mut self,
        request: crate::gradient_index::event::DepartureRequest,
    ) -> Result<
        Option<crate::gradient_index::event::DepartureNeighborhood>,
        crate::gradient_index::GradientIndexError,
    > {
        if request.contact.witness.departure_radius_m.is_some() {
            return Ok(None);
        }
        self.smooth_neighborhood(request)
            .map_err(|error| self.map_error(error))
    }

    fn signed_departure(
        &mut self,
        contact: crate::gradient_index::event::BoundaryWitness,
        point: Point3,
    ) -> Result<
        crate::gradient_index::event::SignedBoundaryDistance,
        crate::gradient_index::GradientIndexError,
    > {
        self.signed_departure_native(contact, point)
            .map_err(|error| self.map_error(error))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::charts::TriMesh;
    use fs_blake3::hash_domain;

    fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
        let gate = fs_exec::CancelGate::new();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 71,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            f(&cx, &gate);
        });
    }

    fn triangle() -> TriMesh {
        TriMesh::new(
            vec![[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
    }

    fn instance(id: u64, geometry: SharedGeometry, transform: RigidTransform) -> GeometryInstance {
        GeometryInstance::try_new(
            id,
            hash_domain("native-proximity-test", &id.to_le_bytes()),
            geometry,
            transform,
        )
        .unwrap()
    }

    #[test]
    fn g1_rotated_large_translation_bounds_agree_with_independent_world_mesh() {
        with_cx(|cx, _| {
            let angle = 0.73_f64;
            let axis = 14.0_f64.sqrt();
            let q = [
                (angle / 2.0).sin() / axis,
                2.0 * (angle / 2.0).sin() / axis,
                3.0 * (angle / 2.0).sin() / axis,
                (angle / 2.0).cos(),
            ];
            for translation in [[0.0; 3], [1.0e9, -2.0e9, 3.0e9]] {
                let transform = RigidTransform::try_new(q, translation).unwrap();
                let local = triangle();
                let pose = PreparedPose::new(transform).unwrap();
                let mut rounded_world_vertices = Vec::new();
                let mut rounding_envelope = 0.0_f64;
                for &vertex in &local.vertices {
                    let exact_world = iadd(
                        imul(pose.matrix, vertex.map(Interval::point)),
                        pose.translation,
                    );
                    let rounded = transform.transform_point(array_point(vertex));
                    rounding_envelope = rounding_envelope.max(uncertainty(exact_world, rounded));
                    rounded_world_vertices.push(point_array(rounded));
                }
                // This second BLAS contains explicitly transformed binary64
                // vertices. Its difference from the exact affine image is
                // bounded by the maximum endpoint-rounding envelope above.
                let world = TriMesh::new(rounded_world_vertices, vec![[0, 1, 2]]);
                let object = instance(1, SharedGeometry::mesh(local), transform);
                for query in [
                    Point3::new(0.13, 0.2, 0.7),
                    Point3::new(1.4, -0.6, 0.03),
                    Point3::new(-2.0, -1.0, 0.1),
                ] {
                    let point = transform.transform_point(query);
                    let placed = object
                        .nearest_surface_with_cx(cx, point, &mut ProximityBudget::new(1000))
                        .unwrap();
                    let direct = world
                        .nearest_surface_with_cx(cx, point, &mut ProximityBudget::new(1000))
                        .unwrap();
                    let upper = placed.proximity.witness.unwrap().distance_upper_m;
                    let reference_upper = direct.witness.unwrap().distance_upper_m;
                    assert!(placed.proximity.clearance_m <= reference_upper + rounding_envelope);
                    assert!(upper + rounding_envelope >= direct.clearance_m);
                    assert!(upper - placed.proximity.clearance_m < 1.0e-4);
                    assert_eq!(placed.object_id, object.object_id());
                    assert_eq!(placed.frame_identity, object.frame_identity());
                    let normal = placed.proximity.witness.unwrap().normal.unwrap();
                    assert!((normal.norm() - 1.0).abs() < 2e-14);
                }
            }
        });
    }

    #[test]
    fn g3_prepared_scene_prunes_instances_and_shares_mesh_validation() {
        with_cx(|cx, _| {
            let mut mesh = triangle();
            // Unreferenced authoring vertices are validated once, but are not
            // surfaces and must not inflate every placed instance's BLAS box.
            mesh.vertices
                .extend((0..2048).map(|i| [1.0e9 + f64::from(i), 0.0, 0.0]));
            let geometry = SharedGeometry::mesh(mesh);
            let objects = (0..512)
                .map(|index| {
                    instance(
                        index + 1,
                        geometry.clone(),
                        RigidTransform::try_new(
                            [0.0, 0.0, 0.0, 1.0],
                            [index as f64 * 10.0, 0.0, 0.0],
                        )
                        .unwrap(),
                    )
                })
                .collect::<Vec<_>>();
            // Enough work to validate this mesh once, not per placement.
            let prepared = BoundaryGeometry::try_new(cx, &objects, 5000).unwrap();
            assert!(prepared.prepared[0].bounds.max.x < 1.000001);
            assert!(Arc::ptr_eq(
                prepared.prepared[0].patches.as_ref().unwrap(),
                prepared.prepared[511].patches.as_ref().unwrap()
            ));
            let mut budget = ProximityBudget::new(1000);
            let (clearance, closest) = prepared
                .nearest(cx, Point3::new(0.0, 0.0, 0.5), &mut budget, None)
                .unwrap();
            assert_eq!(closest.unwrap().0, 0);
            assert!(clearance > 0.499999 && clearance <= 0.5);
            assert!(
                budget.used() < 128,
                "top-level traversal used {}",
                budget.used()
            );
        });
    }

    #[test]
    #[cfg(feature = "tracer")]
    fn g1_resolved_thin_other_object_caps_departure_ball() {
        use crate::gradient_index::event::BoundaryOracle;
        with_cx(|cx, _| {
            let geometry = SharedGeometry::mesh(triangle());
            let objects = vec![
                instance(8, geometry.clone(), RigidTransform::identity()),
                instance(
                    3,
                    geometry,
                    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0e-5]).unwrap(),
                ),
            ];
            let prepared = BoundaryGeometry::try_new(cx, &objects, 1000).unwrap();
            let mut oracle = prepared.oracle(cx, 1000, 1.0e-8);
            let q = oracle.query(Point3::new(0.0, 0.0, -0.5e-8)).unwrap();
            let witness = q.witness.unwrap();
            assert_eq!(witness.boundary_id, 0);
            let radius = witness.departure_radius_m.unwrap();
            assert!(radius > 0.999e-5 && radius <= 1.0e-5);
        });
    }

    #[test]
    #[cfg(feature = "tracer")]
    fn g4_coincident_objects_and_finite_edge_near_miss_refuse_contact() {
        use crate::gradient_index::{GradientIndexError, event::BoundaryOracle};
        with_cx(|cx, _| {
            let geometry = SharedGeometry::mesh(triangle());
            for gap in [0.0, 0.5e-8] {
                let objects = vec![
                    instance(2, geometry.clone(), RigidTransform::identity()),
                    instance(
                        1,
                        geometry.clone(),
                        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, gap]).unwrap(),
                    ),
                ];
                let prepared = BoundaryGeometry::try_new(cx, &objects, 1000).unwrap();
                let mut oracle = prepared.oracle(cx, 1000, 1.0e-8);
                assert_eq!(
                    oracle.query(Point3::new(0.0, 0.0, -0.5e-8)),
                    Err(GradientIndexError::BoundaryResolutionLimit)
                );
                assert_eq!(oracle.last_error(), Some(ProximityError::AmbiguousContact));
            }
            let objects = vec![instance(
                1,
                SharedGeometry::mesh(TriMesh::new(
                    vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    vec![[0, 1, 2]],
                )),
                RigidTransform::identity(),
            )];
            let prepared = BoundaryGeometry::try_new(cx, &objects, 1000).unwrap();
            let mut oracle = prepared.oracle(cx, 1000, 1.0e-8);
            assert_eq!(
                oracle.query(Point3::new(-0.5e-8, 0.25, -0.2e-8)),
                Err(GradientIndexError::BoundaryResolutionLimit)
            );
            assert_eq!(oracle.last_error(), Some(ProximityError::AmbiguousContact));
        });
    }

    #[test]
    #[cfg(feature = "tracer")]
    fn g4_contact_resolution_cannot_exceed_large_world_coordinate_precision() {
        use crate::gradient_index::{GradientIndexError, event::BoundaryOracle};
        with_cx(|cx, _| {
            let objects = vec![instance(
                1,
                SharedGeometry::mesh(triangle()),
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0e9]).unwrap(),
            )];
            let prepared = BoundaryGeometry::try_new(cx, &objects, 1000).unwrap();
            let mut oracle = prepared.oracle(cx, 1000, 1.0e-9);
            assert_eq!(
                oracle.query(Point3::new(0.0, 0.0, 1.0e9)),
                Err(GradientIndexError::BoundaryResolutionLimit)
            );
            assert_eq!(oracle.last_error(), Some(ProximityError::AmbiguousContact));
            let mut oracle = prepared.oracle(cx, 1000, 1.0e-5);
            assert!(
                oracle
                    .query(Point3::new(0.0, 0.0, 1.0e9))
                    .unwrap()
                    .witness
                    .is_some()
            );
        });
    }

    #[test]
    #[cfg(feature = "tracer")]
    fn g0_exact_distance_chart_arrival_has_no_invented_planar_departure() {
        use crate::gradient_index::event::BoundaryOracle;
        with_cx(|cx, _| {
            let objects = vec![instance(
                1,
                SharedGeometry::chart(fs_geom::fixtures::SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.0,
                }),
                RigidTransform::identity(),
            )];
            let prepared = BoundaryGeometry::try_new(cx, &objects, 1000).unwrap();
            let mut oracle = prepared.oracle(cx, 1000, 1.0e-7);
            let q = oracle.query(Point3::new(0.0, 0.0, 1.0 + 0.5e-7)).unwrap();
            assert!(q.witness.is_some());
            assert_eq!(q.witness.unwrap().departure_radius_m, None);
        });
    }

    #[test]
    #[cfg(feature = "tracer")]
    fn g4_scene_query_budget_is_cumulative_and_cancellation_is_immediate() {
        use crate::gradient_index::{GradientIndexError, event::BoundaryOracle};
        with_cx(|cx, gate| {
            let objects = vec![instance(
                1,
                SharedGeometry::mesh(triangle()),
                RigidTransform::identity(),
            )];
            let prepared = BoundaryGeometry::try_new(cx, &objects, 1000).unwrap();
            let mut oracle = prepared.oracle(cx, 9, 1.0e-8);
            assert!(oracle.query(Point3::new(0.0, 0.0, 0.1)).is_ok());
            assert!(oracle.query(Point3::new(0.0, 0.0, 0.1)).is_ok());
            assert_eq!(
                oracle.query(Point3::new(0.0, 0.0, 0.1)),
                Err(GradientIndexError::BoundaryQueryBudgetExceeded)
            );
            assert_eq!(
                oracle.last_error(),
                Some(ProximityError::WorkBudgetExceeded)
            );
            let mut oracle = prepared.oracle(cx, 1000, 1.0e-8);
            gate.request();
            assert_eq!(
                oracle.query(Point3::new(0.0, 0.0, 0.1)),
                Err(GradientIndexError::Cancelled)
            );
            assert_eq!(oracle.last_error(), Some(ProximityError::Cancelled));
        });
    }
}
