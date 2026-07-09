//! Radius-edge quality refinement with the Ruppert policy floor (bead
//! uee3): insert circumcenters of the worst tets; the SMALL-INPUT-ANGLE
//! POLICY is a minimum-edge floor from the input's closest-pair
//! spacing — any insertion that would land closer than the floor to an
//! existing vertex YIELDS (counted as protected) instead of cascading;
//! escapes are skipped and counted. `split_hull_facets` adds the classical
//! Ruppert rule via DIAMETRAL ENCROACHMENT BALLS (`facet_diametral_ball`): a
//! circumcenter that lands inside a hull facet's minimum-enclosing sphere splits
//! that facet instead of being inserted (the sliver-making insertion is
//! replaced by a lower-dimensional split); an escaping circumcenter that
//! encroaches nothing is an unfixable boundary sliver and is skipped, not split.
//! This encroachment protection is exact-audit-clean and deterministic and
//! MEASURABLY shrinks the convex-hull regression (~8×: the ledgered worst-case
//! radius-edge fell from ~2.8e18 to ~3.5e17 on the tmesh-011 fixture), but does
//! NOT eliminate it: the residual slivers come from near-boundary INTERIOR
//! vertices making thin tets with facet-split points, so true full-Ruppert
//! quality stays coupled to boundary-layer / constrained-recovery refinement,
//! exactly as the classical termination theory says. The flag defaults OFF.
//! Deterministic: worst-first with a canonical tie-break, sequential inserts
//! through the exact-predicate kernel; the exact audit is the trip-wire.

use crate::delaunay::{GHOST, Mesh, MeshError, Tetrahedralization};
use fs_exec::Cx;

/// Refinement policy.
#[derive(Debug, Clone, Copy)]
pub struct RefineOptions {
    /// Radius-edge target (2.0 is the classical safe bound).
    pub max_radius_edge: f64,
    /// Steiner-point budget.
    pub max_steiner: u32,
    /// Split encroached hull facets instead of skipping offenders
    /// whose circumcenters escape (the full-Ruppert upgrade).
    pub split_hull_facets: bool,
    /// Small-angle policy: minimum new-edge floor as a fraction of the
    /// input closest-pair spacing (insertions below it yield).
    pub min_edge_factor: f64,
}

impl Default for RefineOptions {
    fn default() -> Self {
        RefineOptions {
            max_radius_edge: 2.0,
            max_steiner: 2000,
            split_hull_facets: false,
            min_edge_factor: 0.25,
        }
    }
}

/// What refinement did (ledger evidence).
#[derive(Debug, Clone, Copy, Default)]
pub struct RefineStats {
    /// Circumcenters inserted.
    pub steiner_inserted: u32,
    /// Offenders skipped because the circumcenter left the hull.
    pub skipped_outside_hull: u32,
    /// Worst radius-edge ratio before.
    pub worst_before: f64,
    /// Worst radius-edge ratio after.
    pub worst_after: f64,
    /// Offenders left whose action was blocked (policy-protected or
    /// facet-splitting disabled) — COUNTED, not hidden.
    pub unrefinable_remaining: u32,
    /// Offenders left that were still refinable when budgets ran out.
    pub refinable_remaining: u32,
    /// Encroached hull facets split (the full-Ruppert path).
    pub hull_facets_split: u32,
    /// Insertions yielded by the minimum-edge (small-angle) policy.
    pub protected_by_policy: u32,
}

impl RefineStats {
    /// Canonical JSON object.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"steiner_inserted\":{},\"skipped_outside_hull\":{},\
             \"worst_before\":{:.4},\"worst_after\":{:.4},\
             \"unrefinable_remaining\":{},\"refinable_remaining\":{},\
             \"hull_facets_split\":{},\"protected_by_policy\":{}}}",
            self.steiner_inserted,
            self.skipped_outside_hull,
            self.worst_before,
            self.worst_after,
            self.unrefinable_remaining,
            self.refinable_remaining,
            self.hull_facets_split,
            self.protected_by_policy
        )
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Circumcenter of a positively oriented tet (f64 — the point is a
/// STEINER point, not a certificate; exactness lives in the predicates
/// that re-triangulate around it).
fn circumcenter(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> Option<[f64; 3]> {
    let (u, v, w) = (sub(b, a), sub(c, a), sub(d, a));
    let m = [u, v, w];
    let rhs = [0.5 * dot(u, u), 0.5 * dot(v, v), 0.5 * dot(w, w)];
    let den = det3(m);
    if den.abs() < 1e-300 {
        return None;
    }
    let col = |k: usize| {
        let mut mm = m;
        for (row, r) in mm.iter_mut().zip(rhs) {
            row[k] = r;
        }
        det3(mm) / den
    };
    let x = [col(0), col(1), col(2)];
    Some([a[0] + x[0], a[1] + x[1], a[2] + x[2]])
}

/// Radius-edge ratio of a tet (`circumradius / shortest edge`).
fn radius_edge(pts: &[[f64; 3]; 4]) -> Option<f64> {
    let cc = circumcenter(pts[0], pts[1], pts[2], pts[3])?;
    let r = dot(sub(cc, pts[0]), sub(cc, pts[0])).sqrt();
    let mut shortest = f64::INFINITY;
    for i in 0..4 {
        for j in (i + 1)..4 {
            let e = sub(pts[i], pts[j]);
            shortest = shortest.min(dot(e, e).sqrt());
        }
    }
    (shortest > 0.0).then(|| r / shortest)
}

/// The hull-facet split point: the in-plane circumcenter when it lies
/// INSIDE the triangle, else the longest edge's midpoint — an obtuse
/// facet's circumcenter sits outside the facet (still on the hull
/// plane) and inserting it manufactures degenerate flat tets (the
/// battery measured radius-edge ratios of 1e16 before this cascade
/// rule; Ruppert's lower-dimensional split is the classical cure).
fn facet_split_point(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Option<[f64; 3]> {
    let (u, v) = (sub(b, a), sub(c, a));
    let uu = dot(u, u);
    let vv = dot(v, v);
    let uv = dot(u, v);
    let den = 2.0 * (uu * vv - uv * uv);
    if den.abs() < 1e-300 {
        return None;
    }
    let s = (vv * (uu - uv)) / den;
    let t = (uu * (vv - uv)) / den;
    let eps = 1e-6;
    let candidate = if s >= eps && t >= eps && s + t <= 1.0 - eps {
        [
            a[0] + s * u[0] + t * v[0],
            a[1] + s * u[1] + t * v[1],
            a[2] + s * u[2] + t * v[2],
        ]
    } else {
        // Longest-edge midpoint (the lower-dimensional cascade).
        let w = sub(c, b);
        let (lu, lv, lw) = (uu, vv, dot(w, w));
        let (p, q) = if lu >= lv && lu >= lw {
            (a, b)
        } else if lv >= lw {
            (a, c)
        } else {
            (b, c)
        };
        [
            f64::midpoint(p[0], q[0]),
            f64::midpoint(p[1], q[1]),
            f64::midpoint(p[2], q[2]),
        ]
    };
    // Pull strictly into the facet interior (a point EXACTLY on a hull
    // edge is collinear-degenerate for the kernel — the audit went red
    // before this blend; 1/32 is dyadic, preserving G3 behavior).
    let centroid = [
        (a[0] + b[0] + c[0]) / 3.0,
        (a[1] + b[1] + c[1]) / 3.0,
        (a[2] + b[2] + c[2]) / 3.0,
    ];
    let blend = 1.0 / 32.0;
    Some([
        candidate[0] + blend * (centroid[0] - candidate[0]),
        candidate[1] + blend * (centroid[1] - candidate[1]),
        candidate[2] + blend * (centroid[2] - candidate[2]),
    ])
}

/// Minimum enclosing sphere of a triangle facet — `(center, radius²)`. Acute
/// facet → its circumsphere; obtuse → the diametral sphere of the longest edge.
/// This is Shewchuk's subfacet encroachment ball: a vertex inside it would make
/// a boundary sliver (f64 — a Steiner heuristic; the exact audit is the trip-wire).
fn facet_diametral_ball(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Option<([f64; 3], f64)> {
    let (u, v) = (sub(b, a), sub(c, a));
    let uu = dot(u, u);
    let vv = dot(v, v);
    let uv = dot(u, v);
    let den = 2.0 * (uu * vv - uv * uv);
    if den.abs() < 1e-300 {
        return None;
    }
    let s = (vv * (uu - uv)) / den;
    let t = (uu * (vv - uv)) / den;
    if s >= 0.0 && t >= 0.0 && s + t <= 1.0 {
        // Acute: the in-plane circumcenter lies inside → its circumsphere.
        let center = [
            a[0] + s * u[0] + t * v[0],
            a[1] + s * u[1] + t * v[1],
            a[2] + s * u[2] + t * v[2],
        ];
        Some((center, dot(sub(center, a), sub(center, a))))
    } else {
        // Obtuse: the longest edge's diametral sphere encloses the facet.
        let w = sub(c, b);
        let (lu, lv, lw) = (uu, vv, dot(w, w));
        let (p, q) = if lu >= lv && lu >= lw {
            (a, b)
        } else if lv >= lw {
            (a, c)
        } else {
            (b, c)
        };
        let center = [
            f64::midpoint(p[0], q[0]),
            f64::midpoint(p[1], q[1]),
            f64::midpoint(p[2], q[2]),
        ];
        Some((center, dot(sub(center, p), sub(center, p))))
    }
}

/// The first hull (ghost) facet whose diametral ball STRICTLY contains `p` — the
/// boundary facet `p` encroaches. Scans all ghost tets (O(hull) — the 1e7 perf
/// lane owns a spatial index); deterministic in tet-slot order.
fn encroached_hull_facet(mesh: &Mesh, p: [f64; 3]) -> Option<[u32; 3]> {
    for t in 0..mesh.tets.len() {
        if !mesh.alive[t] || !mesh.is_ghost(t as u32) {
            continue;
        }
        let tv = mesh.tets[t];
        let f = [tv[0], tv[1], tv[2]];
        let pts = &mesh.points;
        if let Some((c, r2)) =
            facet_diametral_ball(pts[f[0] as usize], pts[f[1] as usize], pts[f[2] as usize])
        {
            let d = sub(p, c);
            if dot(d, d) < r2 {
                return Some(f);
            }
        }
    }
    None
}

/// The input's closest-pair spacing (deterministic O(n²) at fixture
/// scale — the perf lane owns bigger inputs).
fn closest_pair(points: &[[f64; 3]], upto: usize) -> f64 {
    let mut best = f64::INFINITY;
    for i in 0..upto.min(points.len()) {
        for j in (i + 1)..upto.min(points.len()) {
            let d = sub(points[i], points[j]);
            best = best.min(dot(d, d));
        }
    }
    best.sqrt()
}

/// Refine in place until the radius-edge bound holds or budgets run out.
/// Steiner points append after [`Tetrahedralization::steiner_from`].
///
/// # Errors
/// [`MeshError::Cancelled`] between insertions.
#[allow(clippy::too_many_lines)] // one worst-first loop with full accounting
pub fn refine(
    tetra: &mut Tetrahedralization,
    opts: RefineOptions,
    cx: &Cx<'_>,
) -> Result<RefineStats, MeshError> {
    let mut stats = RefineStats::default();
    let worst_ratio = |t: &Tetrahedralization| -> Vec<(f64, [u32; 4])> {
        let pts = &t.mesh.points;
        let mut offenders: Vec<(f64, [u32; 4])> = t
            .tets()
            .into_iter()
            .filter_map(|tet| {
                let q: [[f64; 3]; 4] = core::array::from_fn(|k| pts[tet[k] as usize]);
                radius_edge(&q).map(|r| (r, tet))
            })
            .filter(|&(r, _)| r > opts.max_radius_edge)
            .collect();
        // Worst first; canonical vertex tuple breaks ties (P2).
        offenders.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(b.1.cmp(&a.1)));
        offenders
    };
    let canon = |t: [u32; 4]| {
        let mut s = t;
        s.sort_unstable();
        s
    };
    let live_set = |t: &Tetrahedralization| -> std::collections::BTreeSet<[u32; 4]> {
        t.tets().into_iter().map(canon).collect()
    };
    let initial = worst_ratio(tetra);
    stats.worst_before = initial.first().map_or(0.0, |o| o.0);
    let edge_floor =
        opts.min_edge_factor * closest_pair(&tetra.mesh.points, tetra.steiner_from as usize);
    // A skipped tet's identity never recurs after it dies (points only
    // accumulate), so the skip set is PERMANENT — each unrefinable tet
    // is counted once.
    let mut skipped: std::collections::BTreeSet<[u32; 4]> = std::collections::BTreeSet::new();
    'rounds: while stats.steiner_inserted < opts.max_steiner {
        cx.checkpoint()?;
        let offenders = worst_ratio(tetra);
        let mut live = live_set(tetra);
        let mut progressed = false;
        for &(_, tet) in &offenders {
            if stats.steiner_inserted >= opts.max_steiner {
                break;
            }
            let key = canon(tet);
            if skipped.contains(&key) || !live.contains(&key) {
                continue; // unrefinable, or killed earlier this round
            }
            let pts = &tetra.mesh.points;
            let q: [[f64; 3]; 4] = core::array::from_fn(|k| pts[tet[k] as usize]);
            let insertable =
                circumcenter(q[0], q[1], q[2], q[3]).filter(|cc| cc.iter().all(|x| x.is_finite()));
            let Some(cc) = insertable else {
                skipped.insert(key);
                continue;
            };
            // Where does the circumcenter land? Inside → interior
            // insertion; on a ghost → the FULL-RUPPERT path splits the
            // encroached hull facet instead of skipping.
            let new_idx = tetra.mesh.points.len() as u32;
            let seed = tetra.mesh.locate(cc, new_idx);
            let escaped = tetra.mesh.tets[seed as usize][3] == GHOST;
            // Full-Ruppert rule: split a hull facet IFF the circumcenter
            // encroaches its diametral ball; otherwise insert the circumcenter
            // (interior) or skip it (an escaping circumcenter that encroaches
            // nothing is an unfixable boundary sliver — domain extension is out
            // of scope). Splitting a facet a far circumcenter merely escaped
            // through was what manufactured the boundary slivers.
            let encroached = opts
                .split_hull_facets
                .then(|| encroached_hull_facet(&tetra.mesh, cc))
                .flatten();
            let mut split_facet = false;
            let candidate = if let Some(f) = encroached {
                split_facet = true;
                let pts = &tetra.mesh.points;
                facet_split_point(pts[f[0] as usize], pts[f[1] as usize], pts[f[2] as usize])
            } else if !escaped {
                Some(cc)
            } else {
                stats.skipped_outside_hull += 1;
                skipped.insert(key);
                continue;
            };
            let Some(point) = candidate else {
                skipped.insert(key);
                continue;
            };
            // Small-angle policy: yield if the new point would create
            // an edge below the floor.
            let nearest = {
                let mut best = f64::INFINITY;
                for p in &tetra.mesh.points {
                    let d = sub(*p, point);
                    best = best.min(dot(d, d));
                }
                best.sqrt()
            };
            if nearest < edge_floor {
                stats.protected_by_policy += 1;
                skipped.insert(key);
                continue;
            }
            tetra.mesh.points.push(point);
            if tetra.mesh.insert(new_idx) {
                stats.steiner_inserted += 1;
                if split_facet {
                    stats.hull_facets_split += 1;
                }
                progressed = true;
                live = live_set(tetra);
            } else {
                skipped.insert(key);
            }
        }
        if !progressed {
            break 'rounds;
        }
    }
    let remaining = worst_ratio(tetra);
    stats.worst_after = remaining.first().map_or(0.0, |o| o.0);
    stats.unrefinable_remaining = remaining
        .iter()
        .filter(|(_, t)| skipped.contains(&canon(*t)))
        .count() as u32;
    stats.refinable_remaining = remaining.len() as u32 - stats.unrefinable_remaining;
    tetra.mesh.stats.tets_final = tetra.tets().len() as u64;
    Ok(stats)
}
