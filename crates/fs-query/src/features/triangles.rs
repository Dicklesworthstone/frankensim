//! Point/triangle distance on an immutable topology with a refittable BVH.
//!
//! Preparation reuses the feature query owner's deterministic median-split
//! hierarchy. A refit updates boxes in linear work without rediscovering edges
//! or sorting triangles. Each bounded query returns both the nearest triangle
//! and EVERY triangle within the requested distance. Distances and barycentric
//! witnesses are floating-point geometry, not interval distance certificates;
//! the outward boxes only make broad-phase pruning conservative.

use super::{Bvh, BvhKind, Feature, FeatureBox, FeatureComplex, MAX_COMPLEX_FEATURES};
use fs_exec::Cx;

/// Conservative admission of declared linear sphere/triangle paths.
pub mod sweep;

/// The exact material feature selected by the barycentric witness.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TriangleFeature {
    Vertex(usize),
    Edge(usize, usize),
    Face,
}

/// Closest point on one nondegenerate triangle, in the caller's coordinate unit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrianglePoint {
    pub point: [f64; 3],
    pub barycentric: [f64; 3],
    pub distance: f64,
    /// Unit normal following the supplied triangle winding.
    pub unit_normal: [f64; 3],
    pub feature: TriangleFeature,
}

/// Input/range/geometry and explicit work/output refusals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TriangleQueryError {
    Invalid(&'static str),
    DegenerateTriangle,
    NumericalRange,
    WorkBudgetExhausted { maximum_visits: usize },
    CandidateBudgetExhausted { maximum_candidates: usize },
    Cancelled,
}
impl core::fmt::Display for TriangleQueryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "triangle distance query refused: {self:?}")
    }
}
impl std::error::Error for TriangleQueryError {}

/// Allocation-free closed-triangle closest point, including edges and vertices.
///
/// Edges are scaled before dot/cross products, so physical lengths are never
/// squared. A dimensionless double-area below `64*EPSILON` refuses rather than
/// returning unstable barycentric material coordinates. Edge projection clamps
/// to the actual closed segment; face coordinates are never tolerance-snapped.
pub fn closest_triangle_point(
    point: [f64; 3],
    triangle: [[f64; 3]; 3],
) -> Result<TrianglePoint, TriangleQueryError> {
    vector(point)?;
    for p in triangle {
        vector(p)?;
    }
    let ab = difference(triangle[1], triangle[0])?;
    let ac = difference(triangle[2], triangle[0])?;
    let bc = difference(triangle[2], triangle[1])?;
    let scale = norm(ab)?.max(norm(ac)?).max(norm(bc)?);
    if scale == 0.0 {
        return Err(TriangleQueryError::DegenerateTriangle);
    }
    let e = ab.map(|x| x / scale);
    let f = ac.map(|x| x / scale);
    let n = vector(cross(e, f))?;
    let area = norm(n)?;
    if area <= 64.0 * f64::EPSILON {
        return Err(TriangleQueryError::DegenerateTriangle);
    }
    let unit_normal = n.map(|x| x / area);
    let offset = difference(point, triangle[0])?;
    let q = vector(offset.map(|x| x / scale))?;
    // Dotting a cross product with the UNIT normal avoids squaring the small
    // scaled area. The component perpendicular to the face cancels directly.
    let u = finite(dot(cross(q, f), unit_normal) / area)?;
    let v = finite(dot(cross(e, q), unit_normal) / area)?;
    let w = finite(1.0 - u - v)?;
    if w >= 0.0 && u >= 0.0 && v >= 0.0 {
        return witness(point, triangle, [w, u, v], unit_normal);
    }
    let mut best: Option<TrianglePoint> = None;
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        let edge = difference(triangle[b], triangle[a])?.map(|x| x / scale);
        let r = vector(difference(point, triangle[a])?.map(|x| x / scale))?;
        let denominator = finite(dot(edge, edge))?;
        if denominator <= 0.0 {
            return Err(TriangleQueryError::DegenerateTriangle);
        }
        let t = finite(dot(r, edge) / denominator)?.clamp(0.0, 1.0);
        let mut weights = [0.0; 3];
        weights[a] = 1.0 - t;
        weights[b] = t;
        let candidate = witness(point, triangle, weights, unit_normal)?;
        if best.is_none_or(|current| candidate.distance < current.distance) {
            best = Some(candidate);
        }
    }
    best.ok_or(TriangleQueryError::NumericalRange)
}

fn witness(
    query: [f64; 3],
    triangle: [[f64; 3]; 3],
    barycentric: [f64; 3],
    unit_normal: [f64; 3],
) -> Result<TrianglePoint, TriangleQueryError> {
    let point = vector(core::array::from_fn(|axis| {
        triangle[0][axis]
            + barycentric[1] * (triangle[1][axis] - triangle[0][axis])
            + barycentric[2] * (triangle[2][axis] - triangle[0][axis])
    }))?;
    let distance = norm(difference(query, point)?)?;
    let feature = if barycentric[0] == 1.0 {
        TriangleFeature::Vertex(0)
    } else if barycentric[1] == 1.0 {
        TriangleFeature::Vertex(1)
    } else if barycentric[2] == 1.0 {
        TriangleFeature::Vertex(2)
    } else if barycentric[2] == 0.0 {
        TriangleFeature::Edge(0, 1)
    } else if barycentric[1] == 0.0 {
        TriangleFeature::Edge(0, 2)
    } else if barycentric[0] == 0.0 {
        TriangleFeature::Edge(1, 2)
    } else {
        TriangleFeature::Face
    };
    Ok(TrianglePoint {
        point,
        barycentric,
        distance,
        unit_normal,
        feature,
    })
}

/// Immutable triangle incidence and spatial partition, prepared once.
pub struct PreparedTriangleSurface {
    vertex_count: usize,
    triangles: Vec<[u32; 3]>,
    bvh: Bvh,
}

impl core::fmt::Debug for PreparedTriangleSurface {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedTriangleSurface")
            .field("vertex_count", &self.vertex_count)
            .field("triangles", &self.triangles.len())
            .finish_non_exhaustive()
    }
}

impl PreparedTriangleSurface {
    /// Build using only face boxes; no duplicated vertex/edge candidates.
    /// Original ordering is retained in every returned triangle identifier.
    pub fn new(
        positions: &[[f64; 3]],
        triangles: &[[u32; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, TriangleQueryError> {
        poll(cx)?;
        if positions.is_empty()
            || triangles.is_empty()
            || positions.len().saturating_add(triangles.len()) > MAX_COMPLEX_FEATURES
        {
            return Err(TriangleQueryError::Invalid(
                "nonempty bounded triangle mesh required",
            ));
        }
        for p in positions {
            poll(cx)?;
            vector(*p)?;
        }
        let boxes = face_boxes(positions, triangles, cx)?;
        let features = (0..triangles.len())
            .map(|i| Feature::Face(u32::try_from(i).expect("admitted feature ceiling")))
            .collect();
        let bvh = Bvh::build(&FeatureComplex { features, boxes }, 0.0);
        poll(cx)?;
        Ok(Self {
            vertex_count: positions.len(),
            triangles: triangles.to_vec(),
            bvh,
        })
    }

    #[must_use]
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// Refit every node from current coordinates without sorting or changing
    /// topology. No partially updated geometry is published on failure.
    pub fn refit<'s, 'p>(
        &'s self,
        positions: &'p [[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<RefitTriangleSurface<'s, 'p>, TriangleQueryError> {
        poll(cx)?;
        if positions.len() != self.vertex_count {
            return Err(TriangleQueryError::Invalid(
                "refit vertex count must match prepared topology",
            ));
        }
        let boxes = face_boxes(positions, &self.triangles, cx)?;
        let bounds = self.refit_bounds(&boxes, cx)?;
        Ok(RefitTriangleSurface {
            prepared: self,
            positions,
            boxes,
            bounds,
        })
    }

    /// Enclose straight vertex interpolation between two coordinate images.
    /// Endpoint-union boxes enclose every intermediate triangle, including
    /// rotation/deformation of its plane. They do not certify orientation or
    /// nonintersection during that interpolation.
    pub fn refit_swept(
        &self,
        start: &[[f64; 3]],
        end: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<SweptTriangleSurface<'_>, TriangleQueryError> {
        poll(cx)?;
        if start.len() != self.vertex_count || end.len() != self.vertex_count {
            return Err(TriangleQueryError::Invalid(
                "swept vertices must match prepared topology",
            ));
        }
        let mut boxes = face_boxes(start, &self.triangles, cx)?;
        let end_boxes = face_boxes(end, &self.triangles, cx)?;
        for (start, end) in boxes.iter_mut().zip(end_boxes) {
            poll(cx)?;
            *start = start.merged(&end);
        }
        let bounds = self.refit_bounds(&boxes, cx)?;
        Ok(SweptTriangleSurface {
            prepared: self,
            boxes,
            bounds,
        })
    }

    fn refit_bounds(
        &self,
        boxes: &[FeatureBox],
        cx: &Cx<'_>,
    ) -> Result<Vec<FeatureBox>, TriangleQueryError> {
        let mut bounds: Vec<_> = self.bvh.nodes.iter().map(|n| n.bounds).collect();
        for (id, node) in self.bvh.nodes.iter().enumerate().rev() {
            poll(cx)?;
            bounds[id] = match node.kind {
                BvhKind::Leaf { start, len } => {
                    let indices = &self.bvh.order[start as usize..(start + len) as usize];
                    let mut bound = boxes[indices[0] as usize];
                    for &index in &indices[1..] {
                        bound = bound.merged(&boxes[index as usize]);
                    }
                    bound
                }
                BvhKind::Internal { left, right } => {
                    bounds[left as usize].merged(&bounds[right as usize])
                }
            };
        }
        Ok(bounds)
    }
}

/// One current coordinate image sharing the original hierarchy and topology.
pub struct RefitTriangleSurface<'s, 'p> {
    prepared: &'s PreparedTriangleSurface,
    positions: &'p [[f64; 3]],
    boxes: Vec<FeatureBox>,
    bounds: Vec<FeatureBox>,
}

/// Cached hierarchy refitted to straight endpoint motion enclosures.
pub struct SweptTriangleSurface<'s> {
    prepared: &'s PreparedTriangleSurface,
    boxes: Vec<FeatureBox>,
    bounds: Vec<FeatureBox>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SweptTriangleCandidates {
    /// Complete potentially overlapping faces, in original triangle order.
    pub triangles: Vec<usize>,
    pub visited_nodes: usize,
    pub tested_triangles: usize,
}

impl SweptTriangleSurface<'_> {
    /// All triangle boxes overlapping a sphere's straight centre segment.
    /// `radius_with_padding` is the complete geometric/enclosure radius;
    /// the caller retains its provenance. Both sides round bounds outward.
    pub fn sphere_candidates(
        &self,
        start: [f64; 3],
        end: [f64; 3],
        radius_with_padding: f64,
        maximum_visits: usize,
        maximum_candidates: usize,
        cx: &Cx<'_>,
    ) -> Result<SweptTriangleCandidates, TriangleQueryError> {
        poll(cx)?;
        vector(start)?;
        vector(end)?;
        if !radius_with_padding.is_finite() || radius_with_padding < 0.0 {
            return Err(TriangleQueryError::Invalid(
                "finite nonnegative swept sphere extent required",
            ));
        }
        let sphere = FeatureBox::of_points(&[start, end]).inflated(radius_with_padding);
        vector(sphere.min)?;
        vector(sphere.max)?;
        let mut triangles = Vec::new();
        let mut stack = vec![0usize];
        let mut visited_nodes = 0usize;
        let mut tested_triangles = 0usize;
        while let Some(id) = stack.pop() {
            poll(cx)?;
            consume(visited_nodes, tested_triangles, maximum_visits)?;
            visited_nodes += 1;
            if !sphere.overlaps(&self.bounds[id]) {
                continue;
            }
            match self.prepared.bvh.nodes[id].kind {
                BvhKind::Leaf { start, len } => {
                    for &index in &self.prepared.bvh.order[start as usize..(start + len) as usize] {
                        consume(visited_nodes, tested_triangles, maximum_visits)?;
                        tested_triangles += 1;
                        let index = index as usize;
                        if sphere.overlaps(&self.boxes[index]) {
                            if triangles.len() == maximum_candidates {
                                return Err(TriangleQueryError::CandidateBudgetExhausted {
                                    maximum_candidates,
                                });
                            }
                            triangles.push(index);
                        }
                    }
                }
                BvhKind::Internal { left, right } => {
                    stack.push(right as usize);
                    stack.push(left as usize);
                }
            }
        }
        triangles.sort_unstable();
        poll(cx)?;
        Ok(SweptTriangleCandidates {
            triangles,
            visited_nodes,
            tested_triangles,
        })
    }
}

/// Complete nearest/nearby result; nearby entries follow original face order.
#[derive(Clone, Debug)]
pub struct TrianglePointQuery {
    pub nearest_triangle: usize,
    pub nearest: TrianglePoint,
    pub nearby: Vec<(usize, TrianglePoint)>,
    pub visited_nodes: usize,
    pub tested_triangles: usize,
}

impl RefitTriangleSurface<'_, '_> {
    /// Outward-rounded box of the complete current surface.
    #[must_use]
    pub fn bounds(&self) -> [[f64; 3]; 2] {
        [self.bounds[0].min, self.bounds[0].max]
    }

    /// Query the nearest point and all triangles at distance <= `within`.
    /// Each visited node AND tested triangle consumes one unit of the shared
    /// work cap. Output exhaustion refuses without publishing a partial set.
    pub fn point_query(
        &self,
        point: [f64; 3],
        within: f64,
        maximum_visits: usize,
        maximum_candidates: usize,
        cx: &Cx<'_>,
    ) -> Result<TrianglePointQuery, TriangleQueryError> {
        poll(cx)?;
        vector(point)?;
        if !within.is_finite() || within < 0.0 {
            return Err(TriangleQueryError::Invalid(
                "finite nonnegative query distance required",
            ));
        }
        let mut nearest: Option<(usize, TrianglePoint)> = None;
        let mut nearby = Vec::new();
        let mut stack = vec![0usize];
        let mut visited_nodes = 0usize;
        let mut tested_triangles = 0usize;
        while let Some(id) = stack.pop() {
            poll(cx)?;
            consume(visited_nodes, tested_triangles, maximum_visits)?;
            visited_nodes += 1;
            let limit = nearest.map_or(f64::INFINITY, |(_, p)| p.distance.max(within));
            if box_distance_lower(self.bounds[id], point)? > limit {
                continue;
            }
            match self.prepared.bvh.nodes[id].kind {
                BvhKind::Leaf { start, len } => {
                    for &index in &self.prepared.bvh.order[start as usize..(start + len) as usize] {
                        let index = index as usize;
                        let limit = nearest.map_or(f64::INFINITY, |(_, p)| p.distance.max(within));
                        if box_distance_lower(self.boxes[index], point)? > limit {
                            continue;
                        }
                        consume(visited_nodes, tested_triangles, maximum_visits)?;
                        tested_triangles += 1;
                        let triangle =
                            self.prepared.triangles[index].map(|i| self.positions[i as usize]);
                        let candidate = closest_triangle_point(point, triangle)?;
                        if nearest.is_none_or(|(i, p)| {
                            candidate.distance < p.distance
                                || (candidate.distance == p.distance && index < i)
                        }) {
                            nearest = Some((index, candidate));
                        }
                        if candidate.distance <= within {
                            if nearby.len() == maximum_candidates {
                                return Err(TriangleQueryError::CandidateBudgetExhausted {
                                    maximum_candidates,
                                });
                            }
                            nearby.push((index, candidate));
                        }
                    }
                }
                BvhKind::Internal { left, right } => {
                    let left = left as usize;
                    let right = right as usize;
                    let dl = box_distance_lower(self.bounds[left], point)?;
                    let dr = box_distance_lower(self.bounds[right], point)?;
                    if dl <= dr {
                        stack.push(right);
                        stack.push(left);
                    } else {
                        stack.push(left);
                        stack.push(right);
                    }
                }
            }
        }
        nearby.sort_unstable_by_key(|(index, _)| *index);
        poll(cx)?;
        let (nearest_triangle, nearest) = nearest.ok_or(TriangleQueryError::NumericalRange)?;
        Ok(TrianglePointQuery {
            nearest_triangle,
            nearest,
            nearby,
            visited_nodes,
            tested_triangles,
        })
    }
}

fn face_boxes(
    positions: &[[f64; 3]],
    triangles: &[[u32; 3]],
    cx: &Cx<'_>,
) -> Result<Vec<FeatureBox>, TriangleQueryError> {
    let mut boxes = Vec::with_capacity(triangles.len());
    for triangle in triangles {
        poll(cx)?;
        if triangle.iter().any(|&i| i as usize >= positions.len()) {
            return Err(TriangleQueryError::Invalid(
                "triangle vertex is out of range",
            ));
        }
        let points = triangle.map(|i| positions[i as usize]);
        // Admission shares precisely the narrow phase's scale/range checks.
        closest_triangle_point(points[0], points)?;
        let bounds = FeatureBox::of_points(&points);
        vector(bounds.min)?;
        vector(bounds.max)?;
        boxes.push(bounds);
    }
    Ok(boxes)
}

// The maximum of downward-rounded axis gaps is a lower bound on Euclidean
// distance, without an unproven directed-rounding claim for sqrt/hypot.
fn box_distance_lower(bounds: FeatureBox, p: [f64; 3]) -> Result<f64, TriangleQueryError> {
    let mut distance = 0.0_f64;
    for (axis, value) in p.into_iter().enumerate() {
        let gap = if value < bounds.min[axis] {
            finite(bounds.min[axis] - value)?.next_down()
        } else if value > bounds.max[axis] {
            finite(value - bounds.max[axis])?.next_down()
        } else {
            0.0
        };
        distance = distance.max(gap);
    }
    Ok(distance)
}
fn consume(
    nodes: usize,
    triangles: usize,
    maximum_visits: usize,
) -> Result<(), TriangleQueryError> {
    if nodes.saturating_add(triangles) >= maximum_visits {
        Err(TriangleQueryError::WorkBudgetExhausted { maximum_visits })
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, TriangleQueryError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TriangleQueryError::NumericalRange)
    }
}
fn vector(value: [f64; 3]) -> Result<[f64; 3], TriangleQueryError> {
    for v in value {
        finite(v)?;
    }
    Ok(value)
}
fn difference(a: [f64; 3], b: [f64; 3]) -> Result<[f64; 3], TriangleQueryError> {
    vector(core::array::from_fn(|i| a[i] - b[i]))
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}
fn norm(v: [f64; 3]) -> Result<f64, TriangleQueryError> {
    finite(v[0].hypot(v[1]).hypot(v[2]))
}
fn poll(cx: &Cx<'_>) -> Result<(), TriangleQueryError> {
    cx.checkpoint().map_err(|_| TriangleQueryError::Cancelled)
}
