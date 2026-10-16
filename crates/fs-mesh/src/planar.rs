//! Constrained planar outlines, hole carving and bounded Delaunay refinement.
//!
//! Orientations and Delaunay decisions use fs-ivl exact signs. Metric choices
//! (circumcentres, angles, sizing) use f64. Success checks all nonexempt angles;
//! work exhaustion is a refusal, never a partially qualified mesh. Acute input
//! vertices retain their coordinates and are explicitly listed. This is not a
//! theorem that every finite-precision PSLG admits the requested quality.

use fs_exec::Cx;
use fs_ivl::{Sign, orient2d};
use std::collections::{BTreeMap, BTreeSet};

#[cfg(test)]
mod tests;
mod triangulation;
mod uniform;

/// Closed curve on t in [0,1], with an explicit bound on |d²p/dt²| [m].
/// The endpoint positions must be identical; NURBS sampling belongs in a caller
/// adapter that supplies a valid derivative bound, not an fs-mesh dependency.
pub trait ParametricCurve {
    fn position(&self, t: f64) -> [f64; 2];
    fn second_derivative_bound_m(&self) -> f64;
}

/// An unclosed point list (the last-to-first segment is implicit), or a curve.
pub enum Loop<'a> {
    Polyline(&'a [[f64; 2]]),
    Parametric {
        curve: &'a dyn ParametricCurve,
        chord_tolerance_m: f64,
    },
}
/// One connected domain. Orientation of supplied loops is immaterial.
pub struct PlanarOutline<'a> {
    pub outer: Loop<'a>,
    pub holes: Vec<Loop<'a>>,
}
/// Caller-owned ceilings, including the number of refinement insertions.
#[derive(Debug, Clone, Copy)]
pub struct PlanarBudget {
    pub max_vertices: usize,
    pub max_triangles: usize,
    pub max_refinement_steps: usize,
}
/// Global edge ceiling and optional local ceiling; neither changes geometry.
pub struct PlanarOptions<'a> {
    pub maximum_edge_m: f64,
    pub minimum_angle_deg: f64,
    pub size_field: Option<&'a dyn Fn([f64; 2]) -> f64>,
    pub budget: PlanarBudget,
}
/// A refusal leaves no partially admitted result.
#[derive(Debug, Clone, PartialEq)]
pub enum PlanarError {
    InvalidInput(&'static str),
    SelfIntersection {
        loop_index: usize,
    },
    HoleOutside {
        hole: usize,
    },
    OverlappingHoles {
        first: usize,
        second: usize,
    },
    Budget {
        what: &'static str,
        requested: usize,
        limit: usize,
    },
    Precision(&'static str),
    Topology(&'static str),
    Cancelled,
}
impl std::fmt::Display for PlanarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FS-MESH-PLANAR: {self:?}")
    }
}
impl std::error::Error for PlanarError {}

/// Boundary loops and source-segment chains refer to the returned vertex IDs.
#[derive(Debug, Clone)]
pub struct PlanarReceipt {
    pub boundary_loops: Vec<Vec<usize>>,
    pub source_segment_chains: Vec<Vec<usize>>,
    pub acute_input_vertices: Vec<usize>,
    pub exempt_triangles: Vec<usize>,
    pub minimum_angle_deg: f64,
    pub polygon_area_m2: f64,
    pub mesh_area_m2: f64,
    pub chord_tolerances_m: Vec<Option<f64>>,
    pub refinement_insertions: usize,
}
#[derive(Debug, Clone)]
pub struct PlanarMesh {
    pub vertices: Vec<[f64; 2]>,
    pub triangles: Vec<[usize; 3]>,
    pub receipt: PlanarReceipt,
}
#[derive(Debug, Clone)]
struct Segment {
    a: usize,
    b: usize,
    source: usize,
    loop_id: usize,
}

fn cap(what: &'static str, requested: usize, limit: usize) -> Result<(), PlanarError> {
    if requested > limit {
        Err(PlanarError::Budget {
            what,
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), PlanarError> {
    cx.checkpoint().map_err(|_| PlanarError::Cancelled)
}
fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn distance(a: [f64; 2], b: [f64; 2]) -> f64 {
    let d = sub(a, b);
    fs_math::det::sqrt(dot(d, d))
}
fn edge(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}
fn opposite(a: Sign, b: Sign) -> bool {
    matches!(
        (a, b),
        (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
    )
}
fn on_segment(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> bool {
    orient2d(a, b, p) == Sign::Zero
        && (0..2).all(|k| p[k] >= a[k].min(b[k]) && p[k] <= a[k].max(b[k]))
}
fn intersects(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    (opposite(orient2d(a, b, c), orient2d(a, b, d))
        && opposite(orient2d(c, d, a), orient2d(c, d, b)))
        || on_segment(a, b, c)
        || on_segment(a, b, d)
        || on_segment(c, d, a)
        || on_segment(c, d, b)
}
fn inside(p: [f64; 2], ring: &[[f64; 2]]) -> bool {
    let mut winding = 0i32;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if on_segment(a, b, p) {
            return true;
        }
        if a[1] <= p[1] && b[1] > p[1] && orient2d(a, b, p) == Sign::Positive {
            winding += 1;
        }
        if a[1] > p[1] && b[1] <= p[1] && orient2d(a, b, p) == Sign::Negative {
            winding -= 1;
        }
    }
    winding != 0
}
fn area(ring: &[[f64; 2]]) -> f64 {
    let origin = ring[0];
    (1..ring.len() - 1)
        .map(|i| {
            let a = sub(ring[i], origin);
            let b = sub(ring[i + 1], origin);
            (a[0] * b[1] - a[1] * b[0]) / 2.0
        })
        .sum::<f64>()
}
fn sample(input: &Loop<'_>, remaining: usize) -> Result<(Vec<[f64; 2]>, Option<f64>), PlanarError> {
    let (points, eps) = match input {
        Loop::Polyline(p) => {
            cap("vertices", p.len(), remaining)?;
            (p.to_vec(), None)
        }
        Loop::Parametric {
            curve,
            chord_tolerance_m: eps,
        } => {
            let bound = curve.second_derivative_bound_m();
            if !eps.is_finite() || *eps <= 0.0 || !bound.is_finite() || bound <= 0.0 {
                return Err(PlanarError::InvalidInput(
                    "positive finite chord tolerance and curvature bound required",
                ));
            }
            if curve.position(0.0) != curve.position(1.0) {
                return Err(PlanarError::InvalidInput("curve is not closed"));
            }
            // Linear interpolation error <= max|p''| Delta(t)^2 / 8.
            let count = fs_math::det::sqrt(bound / (8.0 * eps)).ceil().max(3.0);
            if !count.is_finite() || count > remaining as f64 {
                return Err(PlanarError::Budget {
                    what: "curve vertices",
                    requested: remaining.saturating_add(1),
                    limit: remaining,
                });
            }
            let n = count as usize;
            (
                (0..n)
                    .map(|i| curve.position(i as f64 / n as f64))
                    .collect(),
                Some(*eps),
            )
        }
    };
    if points.len() < 3
        || points
            .iter()
            .flatten()
            .any(|v| !v.is_finite() || v.abs() > 1e50)
    {
        return Err(PlanarError::InvalidInput(
            "loop needs at least three finite, representable coordinates",
        ));
    }
    Ok((points, eps))
}

/// Mesh without changing supplied vertices or replacing a hole by material.
/// Exact predicates are applied to represented coordinates. If a required
/// segment split has no representable collinear midpoint, precision refuses.
pub fn mesh_outline(
    outline: &PlanarOutline<'_>,
    options: &PlanarOptions<'_>,
    cx: &Cx<'_>,
) -> Result<PlanarMesh, PlanarError> {
    checkpoint(cx)?;
    if !options.maximum_edge_m.is_finite()
        || options.maximum_edge_m <= 0.0
        || !options.minimum_angle_deg.is_finite()
        || !(0.0..=20.7).contains(&options.minimum_angle_deg)
    {
        return Err(PlanarError::InvalidInput(
            "positive size and angle in [0,20.7] degrees required",
        ));
    }
    let mut rings = Vec::new();
    let mut tolerances = Vec::new();
    let mut count = 0;
    for input in std::iter::once(&outline.outer).chain(outline.holes.iter()) {
        let (p, eps) = sample(input, options.budget.max_vertices.saturating_sub(count))?;
        count += p.len();
        rings.push(p);
        tolerances.push(eps);
    }
    for (r, ring) in rings.iter().enumerate() {
        checkpoint(cx)?;
        let n = ring.len();
        for i in 0..n {
            if ring[i] == ring[(i + 1) % n] {
                return Err(PlanarError::InvalidInput("duplicate loop vertex"));
            }
            for j in i + 1..n {
                if j == i + 1 || (i == 0 && j == n - 1) {
                    continue;
                }
                if intersects(ring[i], ring[(i + 1) % n], ring[j], ring[(j + 1) % n]) {
                    return Err(PlanarError::SelfIntersection { loop_index: r });
                }
            }
        }
        if area(ring) == 0.0 {
            return Err(PlanarError::InvalidInput("zero-area loop"));
        }
    }
    for h in 1..rings.len() {
        if !rings[h].iter().all(|&p| inside(p, &rings[0])) {
            return Err(PlanarError::HoleOutside { hole: h - 1 });
        }
        for other in 0..h {
            let crosses = (0..rings[h].len()).any(|i| {
                (0..rings[other].len()).any(|j| {
                    intersects(
                        rings[h][i],
                        rings[h][(i + 1) % rings[h].len()],
                        rings[other][j],
                        rings[other][(j + 1) % rings[other].len()],
                    )
                })
            });
            if other == 0 && crosses {
                return Err(PlanarError::HoleOutside { hole: h - 1 });
            }
            if other > 0
                && (crosses
                    || inside(rings[h][0], &rings[other])
                    || inside(rings[other][0], &rings[h]))
            {
                return Err(PlanarError::OverlappingHoles {
                    first: other - 1,
                    second: h - 1,
                });
            }
        }
    }
    let polygon_area =
        area(&rings[0]).abs() - rings.iter().skip(1).map(|p| area(p).abs()).sum::<f64>();
    // A triangulated polygon with k holes and boundary-only vertices has
    // F=V+2k-2. Refuse this necessary lower bound before the initial CDT.
    let minimum_triangles = count
        .checked_add(
            outline
                .holes
                .len()
                .checked_mul(2)
                .ok_or(PlanarError::InvalidInput("hole count overflow"))?,
        )
        .and_then(|n| n.checked_sub(2))
        .ok_or(PlanarError::InvalidInput("triangle count overflow"))?;
    cap("triangles", minimum_triangles, options.budget.max_triangles)?;
    let mut vertices = Vec::new();
    let mut segments = Vec::new();
    let mut acute = Vec::new();
    let mut originals = Vec::new();
    for (r, ring) in rings.iter().enumerate() {
        let base = vertices.len();
        vertices.extend_from_slice(ring);
        for i in 0..ring.len() {
            let a = base + i;
            let b = base + (i + 1) % ring.len();
            originals.push((a, b));
            segments.push(Segment {
                a,
                b,
                source: originals.len() - 1,
                loop_id: r,
            });
            let u = sub(ring[(i + ring.len() - 1) % ring.len()], ring[i]);
            let v = sub(ring[(i + 1) % ring.len()], ring[i]);
            if dot(u, v) > 0.5 * fs_math::det::sqrt(dot(u, u) * dot(v, v)) {
                acute.push(a);
            }
        }
    }
    let mut triangles = triangulation::initial(&vertices, &segments, &rings, cx)?;
    cap("triangles", triangles.len(), options.budget.max_triangles)?;
    let mut insertions = 0_usize;
    loop {
        checkpoint(cx)?;
        let adjacency = triangulation::edges(&triangles);
        if segments
            .iter()
            .any(|s| !adjacency.contains_key(&edge(s.a, s.b)))
        {
            return Err(PlanarError::Topology(
                "a constrained boundary edge was lost",
            ));
        }
        let protected: BTreeSet<_> = segments.iter().map(|s| edge(s.a, s.b)).collect();
        // Only vertices visible through retained material can encroach. Carve
        // holes BEFORE refinement; points across a removed hole cannot trigger
        // an endless mutual segment-bisection cascade.
        let encroached = segments
            .iter()
            .enumerate()
            .find(|(_, s)| {
                adjacency[&edge(s.a, s.b)].iter().any(|&(t, _)| {
                    let c = *triangles[t]
                        .iter()
                        .find(|&&v| v != s.a && v != s.b)
                        .unwrap();
                    dot(
                        sub(vertices[c], vertices[s.a]),
                        sub(vertices[c], vertices[s.b]),
                    ) < 0.0
                        && !acute.contains(&c)
                })
            })
            .map(|(i, _)| i);
        let mut bad = None;
        for (t, tri) in triangles.iter().enumerate() {
            let p = tri.map(|i| vertices[i]);
            let centre = [
                (p[0][0] + p[1][0] + p[2][0]) / 3.0,
                (p[0][1] + p[1][1] + p[2][1]) / 3.0,
            ];
            let size = options
                .size_field
                .map(|f| f(centre))
                .unwrap_or(options.maximum_edge_m)
                .min(options.maximum_edge_m);
            if !size.is_finite() || size <= 0.0 {
                return Err(PlanarError::InvalidInput(
                    "size field must be positive finite",
                ));
            }
            let longest = (0..3)
                .map(|i| distance(p[i], p[(i + 1) % 3]))
                .fold(0.0_f64, f64::max);
            let exempt = tri.iter().any(|i| acute.contains(i));
            if longest > size * (1.0 + 1e-12)
                || (!exempt && minimum_angle(p) < options.minimum_angle_deg * (1.0 - 1e-12))
            {
                bad = Some(t);
                break;
            }
        }
        if encroached.is_none() && bad.is_none() {
            break;
        }
        cap(
            "refinement steps",
            insertions.saturating_add(1),
            options.budget.max_refinement_steps,
        )?;
        cap(
            "vertices",
            vertices.len().saturating_add(1),
            options.budget.max_vertices,
        )?;
        let mut split = encroached;
        let candidate = if let Some(t) = bad {
            triangulation::circumcentre(triangles[t].map(|i| vertices[i]))?
        } else {
            [0.0; 2]
        };
        if split.is_none() {
            split = segments.iter().position(|s| {
                dot(sub(candidate, vertices[s.a]), sub(candidate, vertices[s.b])) < 0.0
            });
        }
        let point = if let Some(i) = split {
            let s = &segments[i];
            let point = [
                f64::midpoint(vertices[s.a][0], vertices[s.b][0]),
                f64::midpoint(vertices[s.a][1], vertices[s.b][1]),
            ];
            if point == vertices[s.a]
                || point == vertices[s.b]
                || !on_segment(vertices[s.a], vertices[s.b], point)
            {
                return Err(PlanarError::Precision(
                    "boundary midpoint is not exactly collinear and interior",
                ));
            }
            point
        } else {
            candidate
        };
        if point.iter().any(|x| !x.is_finite()) || vertices.contains(&point) {
            return Err(PlanarError::Precision(
                "refinement point duplicates a vertex or is nonfinite",
            ));
        }
        let id = vertices.len();
        let triangle_delta = if split.is_some() { 1 } else { 2 };
        cap(
            "triangles",
            triangles.len().saturating_add(triangle_delta),
            options.budget.max_triangles,
        )?;
        vertices.push(point);
        triangulation::insert(&vertices, &mut triangles, id, &protected, cx)?;
        if let Some(i) = split {
            let s = segments[i].clone();
            segments[i].b = id;
            segments.push(Segment { a: id, ..s });
        }
        let protected = segments.iter().map(|s| edge(s.a, s.b)).collect();
        triangulation::legalize(&vertices, &mut triangles, &protected, cx)?;
        insertions += 1;
    }
    let mut chains = Vec::new();
    for (source, &(a, b)) in originals.iter().enumerate() {
        let mut chain = vec![a];
        let mut current = a;
        while current != b {
            let s = segments
                .iter()
                .find(|s| s.source == source && s.a == current)
                .ok_or(PlanarError::Topology("broken source-segment chain"))?;
            chain.push(s.b);
            current = s.b;
            if chain.len() > segments.len() + 1 {
                return Err(PlanarError::Topology("cyclic source-segment chain"));
            }
        }
        chains.push(chain);
    }
    let mut loops = vec![Vec::new(); rings.len()];
    for (source, chain) in chains.iter().enumerate() {
        let r = segments
            .iter()
            .find(|s| s.source == source)
            .unwrap()
            .loop_id;
        loops[r].extend_from_slice(&chain[..chain.len() - 1]);
    }
    let exempt_triangles: Vec<_> = triangles
        .iter()
        .enumerate()
        .filter(|(_, t)| t.iter().any(|v| acute.contains(v)))
        .map(|(i, _)| i)
        .collect();
    let min_angle = triangles
        .iter()
        .map(|t| minimum_angle(t.map(|i| vertices[i])))
        .fold(f64::INFINITY, f64::min);
    let mesh_area = triangles
        .iter()
        .map(|t| area(&t.map(|i| vertices[i])))
        .sum::<f64>();
    let edges = triangulation::edges(&triangles);
    if vertices.len() as isize - edges.len() as isize + triangles.len() as isize
        != 1 - outline.holes.len() as isize
        || (mesh_area - polygon_area).abs() > 1e-12 * polygon_area
    {
        return Err(PlanarError::Topology(
            "Euler characteristic or area conservation failed",
        ));
    }
    Ok(PlanarMesh {
        vertices,
        triangles,
        receipt: PlanarReceipt {
            boundary_loops: loops,
            source_segment_chains: chains,
            acute_input_vertices: acute,
            exempt_triangles,
            minimum_angle_deg: min_angle,
            polygon_area_m2: polygon_area,
            mesh_area_m2: mesh_area,
            chord_tolerances_m: tolerances,
            refinement_insertions: insertions,
        },
    })
}

/// Smallest represented triangle angle [degrees], evaluated by atan2(cross,dot).
pub fn minimum_angle(p: [[f64; 2]; 3]) -> f64 {
    (0..3)
        .map(|i| {
            let a = sub(p[(i + 1) % 3], p[i]);
            let b = sub(p[(i + 2) % 3], p[i]);
            fs_math::det::atan2((a[0] * b[1] - a[1] * b[0]).abs(), dot(a, b)) * 180.0
                / core::f64::consts::PI
        })
        .fold(f64::INFINITY, f64::min)
}
