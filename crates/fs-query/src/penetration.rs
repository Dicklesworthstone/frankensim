//! Certified penetration-depth brackets for proven-overlapping convex sets.
//!
//! A strictly positive common ball proves that the origin is interior to the
//! Minkowski difference `A - B`. From that proof, this module seeds an inner
//! octahedron and expands its convex hull with support points (EPA style).
//! The closest inner-hull face supplies a monotone lower bound; a support
//! plane in that face direction supplies a monotone upper bound.

use std::collections::BTreeMap;

use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

use crate::{ConvexSupportMap, QueryError};

/// Default support-expansion budget for penetration brackets.
pub const CONVEX_PENETRATION_DEFAULT_ITERATIONS: u32 = 128;

/// Hard work ceiling for one penetration query.
pub const CONVEX_PENETRATION_MAX_ITERATIONS: u32 = 1_024;

const CHECKPOINT_STRIDE: u32 = 8;
const MAX_EPA_FACES: usize = 16_384;
// Conservative relative deflation for cross products, face planes, rounded
// Minkowski subtraction, and visibility classification. This is deliberately
// far larger than binary64 epsilon; it weakens tightness, never authority.
const GEOMETRY_GUARD: f64 = 9.094_947_017_729_282e-13; // 2^-40

/// A strictly positive common ball that can admit penetration analysis.
///
/// Fields are private so a scalar or raw boolean cannot masquerade as an
/// overlap proof. The penetration routine revalidates the ball against both
/// support maps, which also prevents a witness from an unrelated pair from
/// authorizing a claim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConvexOverlapWitness {
    center: Point3,
    inradius_lower: f64,
}

impl ConvexOverlapWitness {
    pub(crate) fn from_common_ball(center: Point3, inradius_lower: f64) -> Option<Self> {
        let finite_center = center.x.is_finite() && center.y.is_finite() && center.z.is_finite();
        (finite_center && inradius_lower.is_finite() && inradius_lower > 0.0).then_some(Self {
            center,
            inradius_lower,
        })
    }

    /// Center of the certified ball contained in both bodies.
    #[must_use]
    pub const fn center(self) -> Point3 {
        self.center
    }

    /// Certified lower bound on the common ball's radius.
    #[must_use]
    pub const fn inradius_lower(self) -> f64 {
        self.inradius_lower
    }
}

/// A certified bracket on minimum translational penetration depth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConvexPenetration {
    /// Monotone certified lower bound, strictly positive.
    pub lo: f64,
    /// Monotone certified upper bound from an admitted support direction.
    pub hi: f64,
    /// Support-expansion iterations actually executed.
    pub iterations: u32,
    /// Minkowski support evaluations actually executed.
    pub support_evaluations: u32,
    /// Whether the retained bracket closed to the rounding-scale tolerance.
    pub converged: bool,
    /// Deterministic unit-ish direction furnishing the retained upper bound.
    pub normal: [f64; 3],
}

/// Prove a strictly positive common ball at `center` for two support maps.
///
/// Each map must implement [`ConvexSupportMap::contained_ball_radius`]. The
/// smaller admitted radius is sealed into the returned witness.
///
/// # Errors
/// [`QueryError::ConvexOverlapUnproven`] if either map lacks a positive
/// containment certificate at the requested point.
pub fn convex_overlap_witness(
    a: &dyn ConvexSupportMap,
    b: &dyn ConvexSupportMap,
    center: Point3,
) -> Result<ConvexOverlapWitness, QueryError> {
    let finite = center.x.is_finite() && center.y.is_finite() && center.z.is_finite();
    if !finite {
        return Err(overlap_refusal("common-ball center is non-finite"));
    }
    let radius_a = admitted_radius(a, center)
        .ok_or_else(|| overlap_refusal("input a cannot prove a positive contained ball"))?;
    let radius_b = admitted_radius(b, center)
        .ok_or_else(|| overlap_refusal("input b cannot prove a positive contained ball"))?;
    ConvexOverlapWitness::from_common_ball(center, radius_a.min(radius_b))
        .ok_or_else(|| overlap_refusal("common-ball radius is not strictly positive"))
}

/// Certified EPA-style penetration depth for a proven-overlapping convex pair.
///
/// The witness is revalidated against `a` and `b`. Its common ball implies a
/// centered ball inside `A - B`, which supplies the initial positive lower
/// bound and an inner octahedron. Every expansion adds an actual Minkowski
/// support point, so the polytope remains inside `A - B`. Closest-face
/// distances can only raise `lo`; directional support planes can only lower
/// `hi`. Exhausting the bounded work budget returns the honest retained
/// bracket with `converged == false`.
///
/// # Errors
/// [`QueryError::ConvexOverlapUnproven`] for a missing, touching-only, or
/// pair-mismatched witness; [`QueryError::ConvexInvalidShape`] for a zero work
/// budget; [`QueryError::ConvexInvalidSupport`] for malformed support or hull
/// arithmetic; [`QueryError::Cancelled`] on cancellation.
pub fn convex_penetration_depth(
    a: &dyn ConvexSupportMap,
    b: &dyn ConvexSupportMap,
    witness: &ConvexOverlapWitness,
    max_iterations: u32,
    cx: &Cx<'_>,
) -> Result<ConvexPenetration, QueryError> {
    if max_iterations == 0 {
        return Err(QueryError::ConvexInvalidShape {
            reason: "penetration iteration budget must be positive",
        });
    }
    let witness = revalidate_witness(a, b, *witness)?;
    let inner_radius = twice_lower(witness.inradius_lower)?;
    if !(inner_radius.is_finite() && inner_radius > 0.0) {
        return Err(overlap_refusal(
            "common ball is too small to retain a positive Minkowski inball",
        ));
    }
    let slack = combined_slack(a, b)?;
    let mut vertices = initial_octahedron(inner_radius);
    let mut faces = initial_faces(&vertices)?;
    let budget = max_iterations.min(CONVEX_PENETRATION_MAX_ITERATIONS);
    let mut best_lo = inner_radius;
    let mut best_hi = f64::INFINITY;
    let mut best_normal = Vec3::new(1.0, 0.0, 0.0);
    let mut iterations = 0;
    let mut support_evaluations = 0;
    let mut converged = false;

    for iteration in 0..budget {
        if iteration % CHECKPOINT_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let selected = closest_face(&faces).ok_or_else(hull_refusal)?;
        let hull_lo = faces
            .iter()
            .map(|face| face.distance_lo)
            .fold(f64::INFINITY, f64::min);
        if !hull_lo.is_finite() || hull_lo < 0.0 {
            return Err(hull_refusal());
        }
        best_lo = best_lo.max(hull_lo);

        let direction = faces[selected].normal;
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let support = support_difference(a, b, direction)?;
        support_evaluations += 1;
        iterations = iteration + 1;
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let support_hi =
            directional_support_upper(direction, support.point, slack, support.roundoff)?;
        if support_hi < best_hi {
            best_hi = support_hi;
            best_normal = direction;
        }
        if !(best_hi.is_finite() && best_lo.is_finite() && best_lo <= best_hi) {
            return Err(QueryError::ConvexInvalidSupport {
                at: [best_lo, best_hi, slack],
            });
        }
        let tolerance = (best_hi.max(1.0) * 64.0 * f64::EPSILON).next_up();
        if best_hi - best_lo <= tolerance {
            converged = true;
            break;
        }
        if duplicate_vertex(&vertices, support.point) {
            break;
        }
        let Some(expanded) = expand_hull(&vertices, &faces, support.point, support.roundoff)?
        else {
            break;
        };
        if expanded.len() > MAX_EPA_FACES {
            break;
        }
        vertices.push(support.point);
        faces = expanded;
    }

    if !best_hi.is_finite() {
        return Err(hull_refusal());
    }
    Ok(ConvexPenetration {
        lo: best_lo,
        hi: best_hi,
        iterations,
        support_evaluations,
        converged,
        normal: [best_normal.x, best_normal.y, best_normal.z],
    })
}

fn admitted_radius(map: &dyn ConvexSupportMap, center: Point3) -> Option<f64> {
    map.contained_ball_radius(center)
        .filter(|radius| radius.is_finite() && *radius > 0.0)
}

fn revalidate_witness(
    a: &dyn ConvexSupportMap,
    b: &dyn ConvexSupportMap,
    witness: ConvexOverlapWitness,
) -> Result<ConvexOverlapWitness, QueryError> {
    let radius_a = admitted_radius(a, witness.center)
        .ok_or_else(|| overlap_refusal("input a does not contain the witness ball"))?;
    let radius_b = admitted_radius(b, witness.center)
        .ok_or_else(|| overlap_refusal("input b does not contain the witness ball"))?;
    ConvexOverlapWitness::from_common_ball(
        witness.center,
        witness.inradius_lower.min(radius_a).min(radius_b),
    )
    .ok_or_else(|| overlap_refusal("revalidated common-ball radius is not positive"))
}

fn overlap_refusal(reason: &'static str) -> QueryError {
    QueryError::ConvexOverlapUnproven { reason }
}

fn hull_refusal() -> QueryError {
    QueryError::ConvexInvalidSupport {
        at: [f64::NAN, f64::NAN, f64::NAN],
    }
}

fn twice_lower(radius: f64) -> Result<f64, QueryError> {
    let twice = radius * 2.0;
    if twice.is_finite() && twice > 0.0 {
        Ok(twice.next_down())
    } else {
        Err(QueryError::ConvexInvalidSupport {
            at: [radius, 2.0, twice],
        })
    }
}

fn combined_slack(a: &dyn ConvexSupportMap, b: &dyn ConvexSupportMap) -> Result<f64, QueryError> {
    let lhs = a.support_slack();
    let rhs = b.support_slack();
    if !(lhs.is_finite() && rhs.is_finite() && lhs >= 0.0 && rhs >= 0.0) {
        return Err(QueryError::ConvexInvalidSupport {
            at: [lhs, rhs, 0.0],
        });
    }
    let sum = lhs + rhs;
    if !sum.is_finite() {
        return Err(QueryError::ConvexInvalidSupport {
            at: [lhs, rhs, sum],
        });
    }
    Ok(if sum == 0.0 { 0.0 } else { sum.next_up() })
}

fn initial_octahedron(radius: f64) -> Vec<Vec3> {
    vec![
        Vec3::new(radius, 0.0, 0.0),
        Vec3::new(-radius, 0.0, 0.0),
        Vec3::new(0.0, radius, 0.0),
        Vec3::new(0.0, -radius, 0.0),
        Vec3::new(0.0, 0.0, radius),
        Vec3::new(0.0, 0.0, -radius),
    ]
}

fn initial_faces(vertices: &[Vec3]) -> Result<Vec<EpaFace>, QueryError> {
    let indices = [
        [0, 2, 4],
        [0, 4, 3],
        [0, 3, 5],
        [0, 5, 2],
        [1, 4, 2],
        [1, 3, 4],
        [1, 5, 3],
        [1, 2, 5],
    ];
    indices
        .into_iter()
        .map(|indices| EpaFace::new(vertices, indices).ok_or_else(hull_refusal))
        .collect()
}

#[derive(Debug, Clone, Copy)]
struct EpaFace {
    indices: [usize; 3],
    normal: Vec3,
    plane: f64,
    distance_lo: f64,
}

impl EpaFace {
    fn new(vertices: &[Vec3], mut indices: [usize; 3]) -> Option<Self> {
        let p = *vertices.get(indices[0])?;
        let q = *vertices.get(indices[1])?;
        let r = *vertices.get(indices[2])?;
        let mut raw = cross(sub(q, p), sub(r, p));
        let raw_norm = raw.norm();
        if !(raw_norm.is_finite() && raw_norm > 0.0) {
            return None;
        }
        if raw.dot(p) < 0.0 {
            indices.swap(1, 2);
            raw = raw.scale(-1.0);
        }
        let normal = raw.scale(1.0 / raw_norm);
        let normal_norm_hi = norm_upper(normal);
        if !(normal.x.is_finite()
            && normal.y.is_finite()
            && normal.z.is_finite()
            && normal_norm_hi.is_finite()
            && normal_norm_hi > 0.0)
        {
            return None;
        }
        let plane = normal.dot(p);
        let coordinate_scale = vertices
            .iter()
            .flat_map(|v| [v.x.abs(), v.y.abs(), v.z.abs()])
            .fold(1.0_f64, f64::max);
        let guard = (coordinate_scale * GEOMETRY_GUARD).next_up();
        let projected_lo = dot_lower(normal, p);
        let distance_lo = ((projected_lo / normal_norm_hi).next_down() - guard)
            .next_down()
            .max(0.0);
        if !(plane.is_finite() && plane > 0.0 && distance_lo.is_finite()) {
            return None;
        }
        Some(Self {
            indices,
            normal,
            plane,
            distance_lo,
        })
    }

    fn visibility(&self, point: Vec3, scale: f64) -> bool {
        let guard = (scale * GEOMETRY_GUARD).next_up();
        self.normal.dot(point) > self.plane + guard
    }
}

fn closest_face(faces: &[EpaFace]) -> Option<usize> {
    faces
        .iter()
        .enumerate()
        .min_by(|(_, lhs), (_, rhs)| {
            lhs.distance_lo
                .total_cmp(&rhs.distance_lo)
                .then_with(|| lhs.indices.cmp(&rhs.indices))
        })
        .map(|(index, _)| index)
}

#[derive(Debug, Clone, Copy)]
struct SupportVertex {
    point: Vec3,
    roundoff: f64,
}

fn support_difference(
    a: &dyn ConvexSupportMap,
    b: &dyn ConvexSupportMap,
    direction: Vec3,
) -> Result<SupportVertex, QueryError> {
    let pa = a.support_point(direction);
    let pb = b.support_point(direction.scale(-1.0));
    for point in [pa, pb] {
        if !(point.x.is_finite() && point.y.is_finite() && point.z.is_finite()) {
            return Err(QueryError::ConvexInvalidSupport {
                at: [point.x, point.y, point.z],
            });
        }
    }
    let support = pa.delta_from(pb);
    if support.x.is_finite() && support.y.is_finite() && support.z.is_finite() {
        let roundoff = norm_upper(Vec3::new(
            subtraction_roundoff(support.x),
            subtraction_roundoff(support.y),
            subtraction_roundoff(support.z),
        ));
        if roundoff.is_finite() {
            Ok(SupportVertex {
                point: support,
                roundoff,
            })
        } else {
            Err(QueryError::ConvexInvalidSupport {
                at: [support.x, support.y, support.z],
            })
        }
    } else {
        Err(QueryError::ConvexInvalidSupport {
            at: [support.x, support.y, support.z],
        })
    }
}

fn directional_support_upper(
    direction: Vec3,
    support: Vec3,
    slack: f64,
    subtraction_roundoff: f64,
) -> Result<f64, QueryError> {
    let norm_lo = norm_lower(direction);
    let norm_hi = norm_upper(direction);
    if !(norm_lo.is_finite() && norm_hi.is_finite() && norm_lo > 0.0) {
        return Err(QueryError::ConvexInvalidSupport {
            at: [direction.x, direction.y, direction.z],
        });
    }
    let total_error = (slack + subtraction_roundoff).next_up();
    let error_term = (total_error * norm_hi).next_up();
    let numerator = (dot_upper(direction, support) + error_term).next_up();
    if !(total_error.is_finite() && error_term.is_finite() && numerator.is_finite()) {
        return Err(QueryError::ConvexInvalidSupport {
            at: [slack, subtraction_roundoff, norm_hi],
        });
    }
    let upper = (numerator / norm_lo).next_up();
    if upper.is_finite() && upper > 0.0 {
        Ok(upper)
    } else {
        Err(QueryError::ConvexInvalidSupport {
            at: [support.x, support.y, support.z],
        })
    }
}

fn duplicate_vertex(vertices: &[Vec3], candidate: Vec3) -> bool {
    let scale = vertices
        .iter()
        .flat_map(|v| [v.x.abs(), v.y.abs(), v.z.abs()])
        .chain([candidate.x.abs(), candidate.y.abs(), candidate.z.abs()])
        .fold(1.0_f64, f64::max);
    let tolerance = scale * GEOMETRY_GUARD;
    vertices.iter().any(|vertex| {
        let delta = sub(*vertex, candidate);
        delta.dot(delta) <= tolerance * tolerance
    })
}

fn expand_hull(
    vertices: &[Vec3],
    faces: &[EpaFace],
    support: Vec3,
    support_roundoff: f64,
) -> Result<Option<Vec<EpaFace>>, QueryError> {
    let scale = vertices
        .iter()
        .flat_map(|v| [v.x.abs(), v.y.abs(), v.z.abs()])
        .chain([support.x.abs(), support.y.abs(), support.z.abs()])
        .fold(1.0_f64, f64::max);
    let visible: Vec<bool> = faces
        .iter()
        .map(|face| face.visibility(support, scale))
        .collect();
    if !visible.iter().any(|is_visible| *is_visible) {
        return Ok(None);
    }

    let mut edges: BTreeMap<(usize, usize), ((usize, usize), u8)> = BTreeMap::new();
    for (face, is_visible) in faces.iter().zip(&visible) {
        if !is_visible {
            continue;
        }
        for (from, to) in [
            (face.indices[0], face.indices[1]),
            (face.indices[1], face.indices[2]),
            (face.indices[2], face.indices[0]),
        ] {
            let key = (from.min(to), from.max(to));
            let entry = edges.entry(key).or_insert(((from, to), 0));
            entry.1 = entry.1.saturating_add(1);
        }
    }
    if edges.values().any(|(_, count)| *count > 2) {
        return Err(hull_refusal());
    }
    let horizon: Vec<(usize, usize)> = edges
        .into_values()
        .filter_map(|(edge, count)| (count == 1).then_some(edge))
        .collect();
    if horizon.len() < 3 {
        return Ok(None);
    }

    let mut expanded: Vec<EpaFace> = faces
        .iter()
        .zip(&visible)
        .filter_map(|(face, is_visible)| (!is_visible).then_some(*face))
        .collect();
    let new_index = vertices.len();
    let mut all_vertices = vertices.to_vec();
    all_vertices.push(support);
    for (from, to) in horizon {
        let face = EpaFace::new(&all_vertices, [from, to, new_index]).ok_or_else(hull_refusal)?;
        expanded.push(face);
    }
    if expanded.is_empty() {
        return Ok(None);
    }
    let guard = (scale * GEOMETRY_GUARD * 4.0 + support_roundoff).next_up();
    let valid = expanded.iter().all(|face| {
        all_vertices
            .iter()
            .all(|vertex| face.normal.dot(*vertex) <= face.plane + guard)
    });
    Ok(valid.then_some(expanded))
}

fn sub(lhs: Vec3, rhs: Vec3) -> Vec3 {
    Vec3::new(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z)
}

fn cross(lhs: Vec3, rhs: Vec3) -> Vec3 {
    Vec3::new(
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x,
    )
}

fn subtraction_roundoff(value: f64) -> f64 {
    if value == 0.0 {
        f64::from_bits(1)
    } else {
        let up = (value.next_up() - value).abs();
        let down = (value - value.next_down()).abs();
        up.max(down)
    }
}

fn dot_lower(lhs: Vec3, rhs: Vec3) -> f64 {
    let sum = ((lhs.x * rhs.x).next_down() + (lhs.y * rhs.y).next_down()).next_down();
    (sum + (lhs.z * rhs.z).next_down()).next_down()
}

fn dot_upper(lhs: Vec3, rhs: Vec3) -> f64 {
    let sum = ((lhs.x * rhs.x).next_up() + (lhs.y * rhs.y).next_up()).next_up();
    (sum + (lhs.z * rhs.z).next_up()).next_up()
}

fn norm_lower(vector: Vec3) -> f64 {
    let square_lower = |value: f64| {
        let square = value * value;
        if square == 0.0 {
            0.0
        } else {
            square.next_down().max(0.0)
        }
    };
    let sum = (square_lower(vector.x) + square_lower(vector.y))
        .next_down()
        .max(0.0);
    let sum = (sum + square_lower(vector.z)).next_down().max(0.0);
    if sum == 0.0 {
        0.0
    } else {
        sum.sqrt().next_down()
    }
}

fn norm_upper(vector: Vec3) -> f64 {
    let square_upper = |value: f64| {
        let square = value * value;
        if square == 0.0 { 0.0 } else { square.next_up() }
    };
    let sum = (square_upper(vector.x) + square_upper(vector.y)).next_up();
    let sum = (sum + square_upper(vector.z)).next_up();
    if sum == 0.0 {
        0.0
    } else {
        sum.sqrt().next_up()
    }
}
