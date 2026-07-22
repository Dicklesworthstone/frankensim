//! Self-intersection freedom as a PROOF over ADMITTED soups: sweep-and-
//! prune broad phase, EXACT narrow phase — plane-separation early exits,
//! then exact edge-vs-triangle tests (four `orient3d` signs each;
//! complete for non-coplanar pairs because every intersection-segment
//! endpoint lies on some edge), with exact 2D `orient2d` handling for the
//! coplanar case.
//!
//! Faces that share a vertex INDEX are not skipped: they are decided
//! against their SHARED FEATURE (the shared vertex, or the shared edge),
//! so a face that pierces a neighbour it happens to touch is caught. The
//! shared-feature decision is exact too — see
//! [`shared_feature_intersect`] for the argument.
//!
//! What a PASS proves is bounded by ADMISSION. The certificate refuses,
//! before any predicate runs, on out-of-range vertex indices, on
//! non-finite coordinates (`orient3d`'s own finiteness assertion is
//! unreachable behind a NaN-dropping broad phase, so the guard has to be
//! here), and on exactly-degenerate faces (a zero-area face has no plane,
//! and every plane-based argument below silently degenerates on one). A
//! refusal is recorded, never absorbed: [`SelfIntersectReport::proven_free`]
//! is false whenever any refusal is present.
//!
//! Exactness makes a false PASS impossible on the pairs actually decided;
//! configurations in exact contact (shared plane touching, coincident
//! patches) are reported CONSERVATIVELY as intersections of kind
//! `Touching` — the bounded, listed false-FAIL class the acceptance
//! contract allows.

use fs_geom::Point3;
use fs_ivl::{Sign, orient2d, orient3d};
use fs_rep_mesh::Soup;
use std::collections::BTreeSet;

/// How a flagged pair intersects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectKind {
    /// Interiors cross (strict intersection, exact).
    Crossing,
    /// Exact contact / coplanar overlap — conservative flag.
    Touching,
}

/// Why the certificate refused part of its input. A report carrying any
/// refusal proves NOTHING about the soup.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelfIntersectRefusal {
    /// A triangle references a vertex index outside `positions`. Nothing
    /// was tested.
    VertexIndexOutOfRange {
        /// Face index.
        face: usize,
        /// The out-of-range vertex index.
        vertex: u32,
    },
    /// A referenced vertex position carries a non-finite coordinate.
    /// Nothing was tested: NaN silently survives an AABB broad phase
    /// (`min`/`max` drop it) and would cull every pair containing it.
    NonFiniteVertex {
        /// The offending vertex index.
        vertex: u32,
    },
    /// A face has EXACTLY zero area (a repeated index, or three exactly
    /// collinear corners). It has no supporting plane, so no pair
    /// containing it was tested.
    DegenerateFace {
        /// Face index.
        face: usize,
    },
}

/// The certificate.
#[derive(Debug, Clone)]
pub struct SelfIntersectReport {
    /// Flagged face pairs. Empty AND `refusals` empty ⟺ PROVEN free of
    /// self-intersections outside the faces' shared vertices/edges
    /// (exact arithmetic).
    pub intersections: Vec<(usize, usize, IntersectKind)>,
    /// Candidate pairs the exact narrow phase examined.
    pub pairs_tested: u64,
    /// The subset of `pairs_tested` that share at least one vertex index
    /// and were decided against their shared feature rather than skipped.
    pub shared_feature_pairs_tested: u64,
    /// Inputs the certificate REFUSED, in deterministic order. Non-empty
    /// ⟹ no proof is claimed for any part of the soup.
    pub refusals: Vec<SelfIntersectRefusal>,
}

impl SelfIntersectReport {
    /// True when the whole soup was ADMITTED: in-range indices, finite
    /// coordinates, no exactly-degenerate face. A report that is not
    /// admitted carries no claim.
    #[must_use]
    pub fn admitted(&self) -> bool {
        self.refusals.is_empty()
    }

    /// True when non-intersection is PROVEN.
    ///
    /// The proof covers every face pair whose bounding boxes overlap,
    /// including pairs that share a vertex or an edge: those are decided
    /// against their shared feature, so a legitimate adjacency contact
    /// passes while a neighbour-piercing face does not. It is false
    /// whenever the soup was not [admitted](Self::admitted).
    #[must_use]
    pub fn proven_free(&self) -> bool {
        self.refusals.is_empty() && self.intersections.is_empty()
    }
}

fn p3(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

fn finite_point(p: Point3) -> bool {
    p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
}

/// Drop coordinate `drop` (0 → keep `(y, z)`, 1 → `(z, x)`, 2 → `(x, y)`).
fn project(p: [f64; 3], drop: usize) -> [f64; 2] {
    match drop {
        0 => [p[1], p[2]],
        1 => [p[2], p[0]],
        _ => [p[0], p[1]],
    }
}

/// Which coordinate to DROP for an exact, INJECTIVE planar projection of
/// `t`; `None` ⟺ `t` has exactly zero area.
///
/// Each candidate answers one component of the triangle's normal, and
/// every component of `(t₁−t₀) × (t₂−t₀)` is literally the 2D orientation
/// determinant of the corresponding axis-aligned projection — so
/// `orient2d` decides both the degeneracy and the projection choice
/// EXACTLY, with no floating-point magnitude heuristic. A nonzero normal
/// component along the dropped axis is exactly the condition for the
/// projection to be injective on the triangle's plane.
fn projection_axis(t: &[[f64; 3]; 3]) -> Option<usize> {
    (0..3).find(|&drop| {
        orient2d(
            project(t[0], drop),
            project(t[1], drop),
            project(t[2], drop),
        ) != Sign::Zero
    })
}

/// Exact sign classification of `t2`'s plane against `t1`'s vertices.
fn plane_signs(t_plane: &[[f64; 3]; 3], pts: &[[f64; 3]; 3]) -> [Sign; 3] {
    core::array::from_fn(|i| orient3d(t_plane[0], t_plane[1], t_plane[2], pts[i]))
}

fn all(signs: [Sign; 3], s: Sign) -> bool {
    signs.iter().all(|&x| x == s)
}

/// Exact triangle-triangle intersection (closed triangles). `None`
/// means PROVEN disjoint; `Some(kind)` localizes the contact class.
///
/// Every coordinate must be FINITE — `orient3d` asserts it. An exactly
/// degenerate (zero-area) triangle has no supporting plane, so the
/// plane-based argument cannot prove disjointness; such a pair is
/// reported CONSERVATIVELY as `Touching` rather than passed.
#[must_use]
pub fn tri_tri_intersect(t1: [Point3; 3], t2: [Point3; 3]) -> Option<IntersectKind> {
    let a = [p3(t1[0]), p3(t1[1]), p3(t1[2])];
    let b = [p3(t2[0]), p3(t2[1]), p3(t2[2])];
    tri_tri_intersect_raw(&a, &b)
}

fn tri_tri_intersect_raw(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> Option<IntersectKind> {
    if projection_axis(a).is_none() || projection_axis(b).is_none() {
        // No plane: the separation argument below is vacuous on a
        // degenerate face. Fail closed.
        return Some(IntersectKind::Touching);
    }
    let sa = plane_signs(b, a); // T1's vertices vs plane(T2)
    if all(sa, Sign::Positive) || all(sa, Sign::Negative) {
        return None; // strictly separated by plane(T2): PROVEN
    }
    let sb = plane_signs(a, b);
    if all(sb, Sign::Positive) || all(sb, Sign::Negative) {
        return None;
    }
    if sa == [Sign::Zero; 3] {
        // Coplanar: exact 2D overlap test.
        return coplanar_overlap(a, b);
    }
    // General case: for non-coplanar triangles, any intersection
    // segment ends on an edge of one of them — so T1 ∩ T2 ≠ ∅ iff
    // some edge of T1 meets T2 or some edge of T2 meets T1. Each
    // edge-triangle test is four exact orient3d signs.
    let mut touching = false;
    for i in 0..3 {
        match segment_triangle(a[i], a[(i + 1) % 3], b) {
            Some(IntersectKind::Crossing) => return Some(IntersectKind::Crossing),
            Some(IntersectKind::Touching) => touching = true,
            None => {}
        }
        match segment_triangle(b[i], b[(i + 1) % 3], a) {
            Some(IntersectKind::Crossing) => return Some(IntersectKind::Crossing),
            Some(IntersectKind::Touching) => touching = true,
            None => {}
        }
    }
    touching.then_some(IntersectKind::Touching)
}

/// Decide two faces that share at least one vertex INDEX against their
/// SHARED FEATURE. `None` ⟺ PROVEN that the two closed triangles meet in
/// exactly that shared feature — the legitimate adjacency contact.
/// `Some(kind)` ⟹ they meet somewhere else too.
///
/// Both triangles must be exactly non-degenerate with finite coordinates
/// (the certificate admits that before calling). Write `F` for the shared
/// feature and `π₁`, `π₂` for the two supporting planes.
///
/// - **Three shared vertices.** The faces have identical corner sets, so
///   they coincide as point sets: `Touching`.
/// - **Two shared vertices `P`, `Q`, distinct planes.** Both planes
///   contain the line `PQ`, so `T₁ ∩ T₂ ⊆ π₁ ∩ π₂ = line(PQ)`, and a
///   non-degenerate triangle meets the line through one of its own edges
///   in exactly that edge. Hence `T₁ ∩ T₂ = PQ = F` — no test needed.
/// - **Two shared vertices, coplanar.** The overlap has positive area iff
///   the two opposite corners are on the SAME side of line `PQ`; strictly
///   opposite sides give `T₁ ∩ T₂ = PQ`. One `orient2d` each.
/// - **One shared vertex `P`, distinct planes.** `T₁ ∩ T₂ ⊆ L = π₁ ∩ π₂`.
///   `Iₖ = Tₖ ∩ L` is a chord of `Tₖ` with `P` as an ENDPOINT (a chord
///   through a vertex of a convex set ends there), so `J = I₁ ∩ I₂` is a
///   segment from `P` and `J ⊋ {P}` iff the far endpoint of the shorter
///   chord lies in the other triangle. The far endpoint of `I₁` is
///   exactly `edge(T₁ opposite P) ∩ π₂`; therefore `J ⊋ {P}` iff the edge
///   of `T₁` opposite `P` meets `T₂`, or the edge of `T₂` opposite `P`
///   meets `T₁`. Those are two ordinary exact segment-triangle tests —
///   and, crucially, the edges INCIDENT to `P` (whose contact at `P` is
///   the legitimate adjacency) are never tested.
/// - **One shared vertex, coplanar.** Near `P` each triangle is its
///   corner cone, so `T₁ ∩ T₂ ⊋ {P}` iff the two 2D cones share a ray,
///   which for cones of opening angle `< π` happens iff one cone's
///   generator lies in the other cone. Four exact `orient2d` cone tests.
fn shared_feature_intersect(
    ti: [u32; 3],
    tj: [u32; 3],
    a: &[[f64; 3]; 3],
    b: &[[f64; 3]; 3],
) -> Option<IntersectKind> {
    let mut shared_a = [0usize; 3];
    let mut shared_b = [0usize; 3];
    let mut n_a = 0;
    let mut n_b = 0;
    for k in 0..3 {
        if tj.contains(&ti[k]) {
            shared_a[n_a] = k;
            n_a += 1;
        }
        if ti.contains(&tj[k]) {
            shared_b[n_b] = k;
            n_b += 1;
        }
    }
    match n_a.min(n_b) {
        0 => tri_tri_intersect_raw(a, b),
        3 => Some(IntersectKind::Touching), // identical corner sets
        2 => shared_edge_intersect(&shared_a[..2], &shared_b[..2], a, b),
        _ => shared_vertex_intersect(shared_a[0], shared_b[0], a, b),
    }
}

fn shared_edge_intersect(
    shared_a: &[usize],
    shared_b: &[usize],
    a: &[[f64; 3]; 3],
    b: &[[f64; 3]; 3],
) -> Option<IntersectKind> {
    let a_other = (0..3)
        .find(|k| !shared_a.contains(k))
        .expect("a triangle has three corners");
    let b_other = (0..3)
        .find(|k| !shared_b.contains(k))
        .expect("a triangle has three corners");
    if orient3d(b[0], b[1], b[2], a[a_other]) != Sign::Zero {
        // Distinct planes: the intersection is exactly the shared edge.
        return None;
    }
    // Coplanar. The two triangles hinge on the shared edge; they overlap
    // in positive area unless their opposite corners are STRICTLY on
    // opposite sides of it.
    let Some(drop) = projection_axis(a) else {
        return Some(IntersectKind::Touching);
    };
    let p = project(a[shared_a[0]], drop);
    let q = project(a[shared_a[1]], drop);
    let side_a = orient2d(p, q, project(a[a_other], drop));
    let side_b = orient2d(p, q, project(b[b_other], drop));
    let strictly_opposite = matches!(
        (side_a, side_b),
        (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
    );
    if strictly_opposite {
        None
    } else {
        Some(IntersectKind::Touching)
    }
}

fn shared_vertex_intersect(
    pa: usize,
    pb: usize,
    a: &[[f64; 3]; 3],
    b: &[[f64; 3]; 3],
) -> Option<IntersectKind> {
    let (a1, a2) = (a[(pa + 1) % 3], a[(pa + 2) % 3]);
    let (b1, b2) = (b[(pb + 1) % 3], b[(pb + 2) % 3]);
    let coplanar = orient3d(b[0], b[1], b[2], a1) == Sign::Zero
        && orient3d(b[0], b[1], b[2], a2) == Sign::Zero;
    if !coplanar {
        // Only the edges OPPOSITE the shared vertex carry the question.
        let mut touching = false;
        for probe in [segment_triangle(a1, a2, b), segment_triangle(b1, b2, a)] {
            match probe {
                Some(IntersectKind::Crossing) => return Some(IntersectKind::Crossing),
                Some(IntersectKind::Touching) => touching = true,
                None => {}
            }
        }
        return touching.then_some(IntersectKind::Touching);
    }
    let Some(drop) = projection_axis(a) else {
        return Some(IntersectKind::Touching);
    };
    let apex = project(a[pa], drop);
    let (qa1, qa2) = (project(a1, drop), project(a2, drop));
    let (qb1, qb2) = (project(b1, drop), project(b2, drop));
    let cones_meet = in_cone(apex, qa1, qa2, qb1)
        || in_cone(apex, qa1, qa2, qb2)
        || in_cone(apex, qb1, qb2, qa1)
        || in_cone(apex, qb1, qb2, qa2);
    if cones_meet {
        Some(IntersectKind::Touching)
    } else {
        None
    }
}

/// Exact 2D test: is the ray `apex → r` inside the closed cone spanned by
/// `apex → u` and `apex → v`? The cone's opening angle is `< π` (its
/// generators come from a non-degenerate triangle corner), so the cone is
/// the intersection of two half-planes and two `orient2d` signs decide it.
fn in_cone(apex: [f64; 2], u: [f64; 2], v: [f64; 2], r: [f64; 2]) -> bool {
    let (u, v) = match orient2d(apex, u, v) {
        Sign::Positive => (u, v),
        Sign::Negative => (v, u),
        // A flat corner cannot happen behind the degeneracy admission;
        // fail closed if it ever does.
        Sign::Zero => return true,
    };
    orient2d(apex, u, r) != Sign::Negative && orient2d(apex, r, v) != Sign::Negative
}

/// Exact segment-vs-triangle: `(p, q)` against `(a, b, c)`.
/// Strict crossing needs the endpoints strictly on opposite sides of
/// the plane AND the segment's line passing strictly inside the
/// triangle; any on-boundary sign yields the conservative `Touching`.
fn segment_triangle(p: [f64; 3], q: [f64; 3], t: &[[f64; 3]; 3]) -> Option<IntersectKind> {
    let s1 = orient3d(t[0], t[1], t[2], p);
    let s2 = orient3d(t[0], t[1], t[2], q);
    if s1 == s2 && s1 != Sign::Zero {
        return None; // both endpoints strictly on one side
    }
    if s1 == Sign::Zero && s2 == Sign::Zero {
        return None; // collinear-with-plane handled by the coplanar path
    }
    // Side volumes: the segment's line passes through the triangle iff
    // the three tetrahedra (p,q,edge) agree in orientation.
    let v1 = orient3d(p, q, t[0], t[1]);
    let v2 = orient3d(p, q, t[1], t[2]);
    let v3 = orient3d(p, q, t[2], t[0]);
    let signs = [v1, v2, v3];
    let pos = signs.iter().filter(|&&s| s == Sign::Positive).count();
    let neg = signs.iter().filter(|&&s| s == Sign::Negative).count();
    if pos > 0 && neg > 0 {
        return None; // the line misses the triangle: PROVEN
    }
    let boundary = signs.contains(&Sign::Zero) || s1 == Sign::Zero || s2 == Sign::Zero;
    if boundary {
        Some(IntersectKind::Touching)
    } else {
        Some(IntersectKind::Crossing)
    }
}

/// Exact coplanar overlap: any 2D edge pair crosses, or one triangle's
/// vertex lies inside the other (all via exact `orient2d`).
fn coplanar_overlap(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> Option<IntersectKind> {
    // Projection axis: drop a coordinate whose normal component is
    // EXACTLY nonzero, so the projection is injective on the shared
    // plane. `projection_axis` decides that with `orient2d`, not with a
    // floating-point magnitude comparison.
    let Some(drop) = projection_axis(a).or_else(|| projection_axis(b)) else {
        return Some(IntersectKind::Touching);
    };
    let q = |p: [f64; 3]| -> [f64; 2] { project(p, drop) };
    let qa = [q(a[0]), q(a[1]), q(a[2])];
    let qb = [q(b[0]), q(b[1]), q(b[2])];
    // Segment-pair crossings.
    for i in 0..3 {
        for j in 0..3 {
            let (p1, p2) = (qa[i], qa[(i + 1) % 3]);
            let (p3v, p4) = (qb[j], qb[(j + 1) % 3]);
            let d1 = orient2d(p3v, p4, p1);
            let d2 = orient2d(p3v, p4, p2);
            let d3 = orient2d(p1, p2, p3v);
            let d4 = orient2d(p1, p2, p4);
            let opposite = |x: Sign, y: Sign| {
                matches!(
                    (x, y),
                    (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
                )
            };
            if opposite(d1, d2) && opposite(d3, d4) {
                return Some(IntersectKind::Touching); // coplanar contact class
            }
        }
    }
    // Containment either way (strict interior via consistent signs).
    let inside = |p: [f64; 2], t: &[[f64; 2]; 3]| -> bool {
        let s0 = orient2d(t[0], t[1], p);
        let s1 = orient2d(t[1], t[2], p);
        let s2 = orient2d(t[2], t[0], p);
        (s0 == s1 && s1 == s2 && s0 != Sign::Zero)
            || (s0 != Sign::Zero || s1 != Sign::Zero || s2 != Sign::Zero)
                && [s0, s1, s2]
                    .iter()
                    .filter(|&&s| s != Sign::Zero)
                    .collect::<Vec<_>>()
                    .windows(2)
                    .all(|w| w[0] == w[1])
    };
    if qa.iter().any(|&p| inside(p, &qb)) || qb.iter().any(|&p| inside(p, &qa)) {
        return Some(IntersectKind::Touching);
    }
    None
}

/// Admission: everything that must hold before an exact predicate may
/// run. Returns the refusals in deterministic order, plus the per-face
/// degeneracy mask (only meaningful when indices and coordinates passed).
fn admit(soup: &Soup) -> (Vec<SelfIntersectRefusal>, Vec<[[f64; 3]; 3]>, Vec<bool>) {
    let mut refusals = Vec::new();
    let vertices = soup.positions.len();
    for (face, t) in soup.triangles.iter().enumerate() {
        for &vertex in t {
            if vertex as usize >= vertices {
                refusals.push(SelfIntersectRefusal::VertexIndexOutOfRange { face, vertex });
            }
        }
    }
    if !refusals.is_empty() {
        return (refusals, Vec::new(), Vec::new());
    }
    let mut non_finite: BTreeSet<u32> = BTreeSet::new();
    for t in &soup.triangles {
        for &vertex in t {
            if !finite_point(soup.positions[vertex as usize]) {
                non_finite.insert(vertex);
            }
        }
    }
    if !non_finite.is_empty() {
        refusals.extend(
            non_finite
                .into_iter()
                .map(|vertex| SelfIntersectRefusal::NonFiniteVertex { vertex }),
        );
        return (refusals, Vec::new(), Vec::new());
    }
    let corners: Vec<[[f64; 3]; 3]> = soup
        .triangles
        .iter()
        .map(|t| t.map(|v| p3(soup.positions[v as usize])))
        .collect();
    let mut degenerate = vec![false; corners.len()];
    for (face, t) in soup.triangles.iter().enumerate() {
        let repeated = t[0] == t[1] || t[1] == t[2] || t[0] == t[2];
        if repeated || projection_axis(&corners[face]).is_none() {
            degenerate[face] = true;
            refusals.push(SelfIntersectRefusal::DegenerateFace { face });
        }
    }
    (refusals, corners, degenerate)
}

/// Prove a soup free of self-intersections. Sweep-and-prune on x, then
/// an exact narrow phase; faces sharing a vertex or an edge are decided
/// against their SHARED FEATURE rather than skipped.
///
/// The soup is ADMITTED first: out-of-range indices, non-finite
/// coordinates, and exactly-degenerate faces are recorded as
/// [`SelfIntersectRefusal`]s and make [`SelfIntersectReport::proven_free`]
/// false. Non-finite coordinates stop the run outright — `orient3d`
/// asserts finiteness, and a NaN would otherwise be dropped by the
/// broad-phase `min`/`max` into an inverted box that culls every pair.
#[must_use]
pub fn self_intersection_certificate(soup: &Soup) -> SelfIntersectReport {
    let (refusals, corners, degenerate) = admit(soup);
    if corners.is_empty() && !soup.triangles.is_empty() {
        // Admission stopped before any geometry was read.
        return SelfIntersectReport {
            intersections: Vec::new(),
            pairs_tested: 0,
            shared_feature_pairs_tested: 0,
            refusals,
        };
    }
    let nf = corners.len();
    // AABBs + sweep order. Every coordinate is finite here (admitted).
    let mut boxes = Vec::with_capacity(nf);
    for ps in &corners {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in ps {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        boxes.push((lo, hi));
    }
    let mut order: Vec<usize> = (0..nf).collect();
    order.sort_by(|&i, &j| {
        boxes[i].0[0]
            .partial_cmp(&boxes[j].0[0])
            .expect("admission proved every coordinate finite")
            .then(i.cmp(&j))
    });
    let mut intersections = Vec::new();
    let mut pairs_tested = 0u64;
    let mut shared_feature_pairs_tested = 0u64;
    for (oi, &i) in order.iter().enumerate() {
        if degenerate[i] {
            continue;
        }
        for &j in order.iter().skip(oi + 1) {
            if boxes[j].0[0] > boxes[i].1[0] {
                break; // sweep axis separation: no further overlaps
            }
            if degenerate[j] {
                continue; // already refused; it has no plane to argue with
            }
            // Remaining axes.
            if boxes[j].0[1] > boxes[i].1[1]
                || boxes[i].0[1] > boxes[j].1[1]
                || boxes[j].0[2] > boxes[i].1[2]
                || boxes[i].0[2] > boxes[j].1[2]
            {
                continue;
            }
            let ti = soup.triangles[i];
            let tj = soup.triangles[j];
            let shares_vertex = ti.iter().any(|v| tj.contains(v));
            pairs_tested += 1;
            let kind = if shares_vertex {
                shared_feature_pairs_tested += 1;
                shared_feature_intersect(ti, tj, &corners[i], &corners[j])
            } else {
                tri_tri_intersect_raw(&corners[i], &corners[j])
            };
            if let Some(kind) = kind {
                let (lo, hi) = (i.min(j), i.max(j));
                intersections.push((lo, hi, kind));
            }
        }
    }
    intersections.sort_unstable_by_key(|&(a, b, _)| (a, b));
    SelfIntersectReport {
        intersections,
        pairs_tested,
        shared_feature_pairs_tested,
        refusals,
    }
}
