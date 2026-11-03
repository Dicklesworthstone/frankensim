//! Six-node topology on the original feature BVH, with linear-time refits.

use super::{
    CellBudget, QuadraticQueryControls, QuadraticQueryError, QuadraticTrianglePoint,
    finite_interval, hull, relative_bernstein_controls, solve, validate_controls,
};
use super::super::{
    Bvh, BvhKind, Feature, FeatureBox, FeatureComplex, MAX_COMPLEX_FEATURES,
    SweptTriangleCandidates, TriangleQueryError, box_distance_lower, consume, poll, vector,
};
use fs_exec::Cx;
use fs_ivl::Interval;

/// Immutable six-node incidence and the feature owner's deterministic BVH.
/// Refit changes only geometric boxes, retaining original patch identifiers.
pub struct PreparedQuadraticTriangleSurface {
    vertex_count: usize,
    triangles: Vec<[u32; 6]>,
    bvh: Bvh,
}

impl core::fmt::Debug for PreparedQuadraticTriangleSurface {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedQuadraticTriangleSurface")
            .field("vertex_count", &self.vertex_count)
            .field("triangles", &self.triangles.len())
            .finish_non_exhaustive()
    }
}

impl PreparedQuadraticTriangleSurface {
    pub fn new(
        positions: &[[f64; 3]], triangles: &[[u32; 6]], cx: &Cx<'_>,
    ) -> Result<Self, QuadraticQueryError> {
        poll(cx)?;
        if positions.is_empty() || triangles.is_empty()
            || positions.len().saturating_add(triangles.len()) > MAX_COMPLEX_FEATURES
        {
            return Err(QuadraticQueryError::Invalid("nonempty bounded six-node triangle mesh required"));
        }
        let boxes = face_boxes(positions, triangles, cx)?;
        let features = (0..triangles.len()).map(|i| Feature::Face(u32::try_from(i).expect("admitted feature ceiling"))).collect();
        let bvh = Bvh::build(&FeatureComplex { features, boxes }, 0.0);
        poll(cx)?;
        Ok(Self { vertex_count: positions.len(), triangles: triangles.to_vec(), bvh })
    }

    #[must_use]
    pub fn triangles(&self) -> &[[u32; 6]] { &self.triangles }

    /// Refit the exact quadratic Bernstein hulls, not the corner triangles.
    /// All geometry is validated before a new coordinate image is published.
    pub fn refit<'s, 'p>(
        &'s self, positions: &'p [[f64; 3]], cx: &Cx<'_>,
    ) -> Result<RefitQuadraticTriangleSurface<'s, 'p>, QuadraticQueryError> {
        if positions.len() != self.vertex_count {
            return Err(QuadraticQueryError::Invalid("refit vertex count must match prepared topology"));
        }
        let boxes = face_boxes(positions, &self.triangles, cx)?;
        let bounds = self.refit_bounds(&boxes, cx)?;
        Ok(RefitQuadraticTriangleSurface { prepared: self, positions, boxes, bounds })
    }

    /// Endpoint hull unions enclose every declared straight nodal interpolation
    /// because the Bernstein conversion is linear in the six node positions.
    pub fn refit_swept<'s>(
        &'s self, start: &[[f64; 3]], end: &[[f64; 3]], cx: &Cx<'_>,
    ) -> Result<SweptQuadraticTriangleSurface<'s>, QuadraticQueryError> {
        if start.len() != self.vertex_count || end.len() != self.vertex_count {
            return Err(QuadraticQueryError::Invalid("swept vertex count must match prepared topology"));
        }
        let mut boxes = face_boxes(start, &self.triangles, cx)?;
        let end_boxes = face_boxes(end, &self.triangles, cx)?;
        for (a, b) in boxes.iter_mut().zip(end_boxes) { *a = a.merged(&b); }
        let bounds = self.refit_bounds(&boxes, cx)?;
        Ok(SweptQuadraticTriangleSurface { prepared: self, boxes, bounds })
    }

    fn refit_bounds(&self, boxes: &[FeatureBox], cx: &Cx<'_>) -> Result<Vec<FeatureBox>, QuadraticQueryError> {
        let mut bounds: Vec<_> = self.bvh.nodes.iter().map(|n| n.bounds).collect();
        for (id, node) in self.bvh.nodes.iter().enumerate().rev() {
            poll(cx)?;
            bounds[id] = match node.kind {
                BvhKind::Leaf { start, len } => {
                    let indices = &self.bvh.order[start as usize..(start + len) as usize];
                    let mut bound = boxes[indices[0] as usize];
                    for &index in &indices[1..] { bound = bound.merged(&boxes[index as usize]); }
                    bound
                }
                BvhKind::Internal { left, right } => bounds[left as usize].merged(&bounds[right as usize]),
            };
        }
        Ok(bounds)
    }
}

pub struct RefitQuadraticTriangleSurface<'s, 'p> {
    prepared: &'s PreparedQuadraticTriangleSurface,
    positions: &'p [[f64; 3]],
    boxes: Vec<FeatureBox>,
    bounds: Vec<FeatureBox>,
}

/// Global minimum enclosure and a conservative complete nearby-patch set.
/// Each nearby patch contributes one position-localized witness, not a claim
/// that every stationary point or every sphere intersection was enumerated.
#[derive(Clone, Debug, PartialEq)]
pub struct QuadraticTrianglePointQuery {
    pub nearest_triangle: usize,
    pub nearest: QuadraticTrianglePoint,
    pub distance_lower: f64,
    pub distance_upper: f64,
    /// Encloses oriented side on EVERY patch still eligible to be globally
    /// nearest, including near ties. Positive lower bound certifies that all
    /// such nearest material points face the query with their supplied winding.
    pub nearest_oriented_side_bounds: [f64; 2],
    /// Every patch whose true minimum distance is <= the requested band is
    /// included. Interval uncertainty may include additional patches.
    pub nearby: Vec<(usize, QuadraticTrianglePoint)>,
    pub visited_nodes: usize,
    pub tested_triangles: usize,
    pub examined_cells: usize,
}

impl RefitQuadraticTriangleSurface<'_, '_> {
    #[must_use]
    pub fn bounds(&self) -> [[f64; 3]; 2] { [self.bounds[0].min, self.bounds[0].max] }

    /// Traverse the original BVH once, sharing one global cell budget across
    /// narrow-phase solves. `maximum_visits` counts nodes plus tested patches.
    /// `maximum_candidates` caps the retained union of nearby and potentially
    /// nearest patches (so at least one slot is required). No partial result
    /// is returned on any work, storage, localization, or output refusal.
    pub fn point_query(
        &self, point: [f64; 3], within: f64, maximum_visits: usize,
        maximum_candidates: usize, controls: QuadraticQueryControls, cx: &Cx<'_>,
    ) -> Result<QuadraticTrianglePointQuery, QuadraticQueryError> {
        poll(cx)?;
        vector(point)?;
        validate_controls(controls)?;
        if !within.is_finite() || within < 0.0 {
            return Err(QuadraticQueryError::Invalid("finite nonnegative query distance required"));
        }
        let mut budget = CellBudget { used: 0, maximum: controls.maximum_cells };
        let mut nearest: Option<(usize, QuadraticTrianglePoint)> = None;
        let mut retained: Vec<(usize, QuadraticTrianglePoint)> = Vec::new();
        let mut stack = vec![0usize];
        let mut visited_nodes = 0;
        let mut tested_triangles = 0;
        while let Some(id) = stack.pop() {
            poll(cx)?;
            consume(visited_nodes, tested_triangles, maximum_visits)?;
            visited_nodes += 1;
            let limit = nearest.map_or(f64::INFINITY, |(_, p)| p.distance_upper.max(within));
            if box_distance_lower(self.bounds[id], point)? > limit { continue; }
            match self.prepared.bvh.nodes[id].kind {
                BvhKind::Leaf { start, len } => {
                    for &index in &self.prepared.bvh.order[start as usize..(start + len) as usize] {
                        poll(cx)?;
                        let index = index as usize;
                        let limit = nearest.map_or(f64::INFINITY, |(_, p)| p.distance_upper.max(within));
                        if box_distance_lower(self.boxes[index], point)? > limit { continue; }
                        consume(visited_nodes, tested_triangles, maximum_visits)?;
                        tested_triangles += 1;
                        let nodes = self.prepared.triangles[index].map(|i| self.positions[i as usize]);
                        let candidate = solve(point, nodes, controls, &mut budget, cx)?;
                        if nearest.is_none_or(|(i, p)| candidate.distance_upper < p.distance_upper
                            || (candidate.distance_upper == p.distance_upper && index < i))
                        {
                            nearest = Some((index, candidate));
                            let limit = candidate.distance_upper.max(within);
                            retained.retain(|(_, p)| p.distance_lower <= limit);
                        }
                        let limit = nearest.map_or(f64::INFINITY, |(_, p)| p.distance_upper.max(within));
                        if candidate.distance_lower <= limit {
                            if retained.len() == maximum_candidates {
                                return Err(TriangleQueryError::CandidateBudgetExhausted { maximum_candidates }.into());
                            }
                            retained.push((index, candidate));
                        }
                    }
                }
                BvhKind::Internal { left, right } => {
                    let left = left as usize; let right = right as usize;
                    let dl = box_distance_lower(self.bounds[left], point)?;
                    let dr = box_distance_lower(self.bounds[right], point)?;
                    if dl <= dr { stack.push(right); stack.push(left); }
                    else { stack.push(left); stack.push(right); }
                }
            }
        }
        let (nearest_triangle, nearest) = nearest.ok_or(QuadraticQueryError::NumericalRange)?;
        let mut distance_lower = nearest.distance_upper;
        let mut nearest_oriented_side_bounds = [f64::INFINITY, f64::NEG_INFINITY];
        for (_, candidate) in &retained {
            distance_lower = distance_lower.min(candidate.distance_lower);
            if candidate.distance_lower <= nearest.distance_upper {
                nearest_oriented_side_bounds[0] = nearest_oriented_side_bounds[0].min(candidate.oriented_side_bounds[0]);
                nearest_oriented_side_bounds[1] = nearest_oriented_side_bounds[1].max(candidate.oriented_side_bounds[1]);
            }
        }
        retained.retain(|(_, p)| p.distance_lower <= within);
        retained.sort_unstable_by_key(|(i, _)| *i);
        poll(cx)?;
        Ok(QuadraticTrianglePointQuery {
            nearest_triangle, nearest, distance_lower, distance_upper: nearest.distance_upper,
            nearest_oriented_side_bounds, nearby: retained, visited_nodes, tested_triangles,
            examined_cells: budget.used,
        })
    }
}

pub struct SweptQuadraticTriangleSurface<'s> {
    prepared: &'s PreparedQuadraticTriangleSurface,
    boxes: Vec<FeatureBox>,
    bounds: Vec<FeatureBox>,
}

impl SweptQuadraticTriangleSurface<'_> {
    /// Complete broad-phase candidate set for a straight sphere-centre path
    /// against straight nodal interpolation. Exact path admission belongs to
    /// `quadratic::sweep`; overlapping hulls alone do not establish collision.
    pub fn sphere_candidates(
        &self, start: [f64; 3], end: [f64; 3], radius: f64,
        maximum_visits: usize, maximum_candidates: usize, cx: &Cx<'_>,
    ) -> Result<SweptTriangleCandidates, QuadraticQueryError> {
        poll(cx)?;
        vector(start)?; vector(end)?;
        if !radius.is_finite() || radius < 0.0 {
            return Err(QuadraticQueryError::Invalid("finite nonnegative swept radius required"));
        }
        let query = FeatureBox::of_points(&[start, end]).inflated(radius);
        vector(query.min)?; vector(query.max)?;
        let mut triangles = Vec::new();
        let mut stack = vec![0usize];
        let mut visited_nodes = 0; let mut tested_triangles = 0;
        while let Some(id) = stack.pop() {
            poll(cx)?;
            consume(visited_nodes, tested_triangles, maximum_visits)?;
            visited_nodes += 1;
            if !self.bounds[id].overlaps(&query) { continue; }
            match self.prepared.bvh.nodes[id].kind {
                BvhKind::Leaf { start, len } => {
                    for &index in &self.prepared.bvh.order[start as usize..(start + len) as usize] {
                        consume(visited_nodes, tested_triangles, maximum_visits)?;
                        tested_triangles += 1;
                        if self.boxes[index as usize].overlaps(&query) {
                            if triangles.len() == maximum_candidates {
                                return Err(TriangleQueryError::CandidateBudgetExhausted { maximum_candidates }.into());
                            }
                            triangles.push(index as usize);
                        }
                    }
                }
                BvhKind::Internal { left, right } => { stack.push(right as usize); stack.push(left as usize); }
            }
        }
        triangles.sort_unstable();
        poll(cx)?;
        Ok(SweptTriangleCandidates { triangles, visited_nodes, tested_triangles })
    }
}

fn face_boxes(
    positions: &[[f64; 3]], triangles: &[[u32; 6]], cx: &Cx<'_>,
) -> Result<Vec<FeatureBox>, QuadraticQueryError> {
    for point in positions { poll(cx)?; vector(*point)?; }
    let mut boxes = Vec::with_capacity(triangles.len());
    for triangle in triangles {
        poll(cx)?;
        if triangle.iter().any(|&i| i as usize >= positions.len()) {
            return Err(QuadraticQueryError::Invalid("quadratic triangle node is out of range"));
        }
        let nodes = triangle.map(|i| positions[i as usize]);
        let relative = hull(relative_bernstein_controls(nodes, nodes[0])?);
        let mut min = [0.0; 3]; let mut max = min;
        for axis in 0..3 {
            let bound = relative[axis] + Interval::point(nodes[0][axis]);
            finite_interval(bound)?;
            min[axis] = bound.lo(); max[axis] = bound.hi();
        }
        boxes.push(FeatureBox { min, max });
    }
    Ok(boxes)
}
