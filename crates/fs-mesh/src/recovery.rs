//! Constrained boundary recovery, conforming-Delaunay slice (bead
//! uee3 item 1): every PLC SEGMENT becomes a union of mesh edges by
//! recursive midpoint Steiner insertion — if (a, b) is not an edge of
//! the current tetrahedralization, insert the midpoint and recurse on
//! the halves (the classic stitching argument: sub-segments shorter
//! than the local feature size have empty diametral balls and are
//! Delaunay edges). The BOUNDARY CORRESPONDENCE table maps every
//! recovered sub-edge back to its parent segment BY CONSTRUCTION
//! (the recursion knows which points it created for which segment) —
//! and the battery re-verifies each recorded sub-edge against the
//! finished mesh anyway. Depth/budget caps are counted honestly
//! (`unrecovered`), never silently dropped. Segment Steiner midpoints
//! snap onto the original parent chord so a general-position segment
//! stays on the line under bisection. A second diagonal that crosses
//! the first at its centre adopts that existing vertex when bitwise
//! midpoints no longer match after a general-position rotation.
//!
//! INTERIOR FACET recovery ([`recover_facets`], the uee3/iw3l successor
//! slice): every SIMPLE planar PLC facet becomes a union of mesh FACES
//! by longest-edge midpoint bisection of a valid triangulation — the
//! 2D analogue of the segment argument (sub-triangles below the local
//! feature size have empty min-enclosing balls and are Delaunay faces).
//! The triangulation is EXACT-PREDICATE ROBUST: convex facets keep the
//! cheap fan (unchanged); NON-CONVEX facets are ear-clipped in the
//! facet plane using `orient2d` for the ear and containment tests
//! (bead iw3l item (a): interior/non-convex PLC facets) — a loop that
//! is not a simple non-degenerate polygon is counted `unrecovered`,
//! never faked. Steiner midpoints on INTERIOR triangulation edges are
//! snapped onto the parent facet's supporting plane (Newell normal
//! through the first loop vertex) so a planar facet stays planar under
//! bisection even when it is not axis-aligned. A point already on the
//! plane (`t == 0`) is left bitwise unchanged, which keeps the
//! axis-aligned identity. Constraint edges — the original loop and
//! every bisection descendant — keep the raw midpoint so two
//! non-coplanar parent facets that share a crease adopt the same
//! Steiner vertex.

use crate::delaunay::{GHOST, MeshError, Tetrahedralization};
use fs_exec::Cx;
use fs_ivl::{Sign, orient2d};
use std::collections::BTreeSet;

/// Recovery policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecoveryOptions {
    /// Bisection depth cap per segment (2^depth sub-edges at worst).
    pub max_depth: u32,
    /// Total Steiner budget.
    pub max_steiner: u32,
}

impl Default for RecoveryOptions {
    fn default() -> Self {
        RecoveryOptions {
            max_depth: 12,
            max_steiner: 4000,
        }
    }
}

/// Recovery evidence.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct RecoveryStats {
    /// Segments requested.
    pub segments_in: u64,
    /// Segments fully recovered as edge chains.
    pub recovered: u64,
    /// Segments abandoned at a cap (HONESTY counter — must be zero
    /// for a pass).
    pub unrecovered: u64,
    /// Steiner points inserted on segments.
    pub steiner_inserted: u64,
    /// Deepest bisection level used.
    pub max_depth_used: u32,
    /// Sub-edges in the correspondence table.
    pub sub_edges: u64,
}

impl RecoveryStats {
    /// Canonical JSON ledger row.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"segments_in\":{},\"recovered\":{},\"unrecovered\":{},\
             \"steiner_inserted\":{},\"max_depth_used\":{},\"sub_edges\":{}}}",
            self.segments_in,
            self.recovered,
            self.unrecovered,
            self.steiner_inserted,
            self.max_depth_used,
            self.sub_edges
        )
    }
}

/// The boundary correspondence: every recovered sub-edge (sorted
/// vertex pair) with its parent segment index — the DWR mapping back
/// to source charts.
#[derive(Debug, Clone, Default)]
pub struct Correspondence {
    /// (sub-edge, parent segment) rows, deterministic order.
    pub rows: Vec<([u32; 2], u32)>,
}

/// Live mesh edge set, for the volumetric layer's diagnostics.
pub(crate) fn edge_set_of(tetra: &Tetrahedralization) -> BTreeSet<[u32; 2]> {
    edge_set(tetra)
}

/// The crate's segment membership test (parameter strictly inside `(a, b)`
/// within the chord tolerance), for the volumetric layer's Steiner
/// perturbation: a point on a chord must be moved ALONG it.
pub(crate) fn chord_parameter(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> Option<f64> {
    parameter_on_segment(p, a, b)
}

/// `snap_to_line` for the volumetric layer.
pub(crate) fn snap_point_to_line(point: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    snap_to_line(point, a, b)
}

/// `snap_to_plane` for the volumetric layer.
pub(crate) fn snap_point_to_plane(point: [f64; 3], origin: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    snap_to_plane(point, origin, normal)
}

/// Live mesh edge set (sorted vertex pairs of live real tets).
fn edge_set(tetra: &Tetrahedralization) -> BTreeSet<[u32; 2]> {
    let mut edges = BTreeSet::new();
    for tet in tetra.tets() {
        for i in 0..4 {
            for j in (i + 1)..4 {
                let (a, b) = (tet[i], tet[j]);
                if a == GHOST || b == GHOST {
                    continue;
                }
                edges.insert(if a < b { [a, b] } else { [b, a] });
            }
        }
    }
    edges
}

/// Recover every PLC segment as a chain of mesh edges. Segment
/// endpoints are indices into the ORIGINAL input points (before any
/// Steiner insertion).
///
/// # Errors
/// [`MeshError::Cancelled`] between insertions.
///
/// # Panics
/// Only on kernel programmer contracts.
pub fn recover_segments(
    tetra: &mut Tetrahedralization,
    segments: &[[u32; 2]],
    opts: RecoveryOptions,
    cx: &Cx<'_>,
) -> Result<(RecoveryStats, Correspondence), MeshError> {
    let mut stats = RecoveryStats {
        segments_in: segments.len() as u64,
        ..RecoveryStats::default()
    };
    let mut table = Correspondence::default();
    let mut edges = edge_set(tetra);
    // Coordinate-bits index: a bisection midpoint that ALREADY exists
    // as a vertex (segments crossing at a shared midpoint — the four
    // body diagonals of a box all meet at its center) is ADOPTED, not
    // abandoned: bitwise equality to the exact midpoint of on-segment
    // endpoints puts the twin on the segment by construction.
    let mut by_bits: std::collections::BTreeMap<[u64; 3], u32> = tetra
        .mesh
        .points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            (
                [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()],
                u32::try_from(i).expect("point count fits u32"),
            )
        })
        .collect();
    for (sid, &[a, b]) in segments.iter().enumerate() {
        cx.checkpoint()?;
        // Chain of on-segment vertices, kept in parameter order: the
        // recursion only ever SPLITS an interval, so a sorted list of
        // (dyadic parameter, vertex) is the whole bookkeeping.
        let (oa, ob) = (tetra.mesh.points[a as usize], tetra.mesh.points[b as usize]);
        let mut chain: Vec<(f64, u32)> = vec![(0.0, a), (1.0, b)];
        // Vertices already ON the chord (a re-mesh after Steiner
        // perturbation, or a refinement re-recovery) are part of the chain
        // from the start: a sub-edge that skipped one could never be a mesh
        // edge, and bisection would only mint midpoints beside it until the
        // depth cap. A fresh PLC has none, so its chain is unchanged.
        for v in chain_on_chord(&tetra.mesh.points, a, b) {
            if v != a
                && v != b
                && let Some(t) = parameter_on_segment(tetra.mesh.points[v as usize], oa, ob)
            {
                chain.push((t, v));
            }
        }
        chain.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        // Work stack of open sub-intervals (param lo, vert lo, param
        // hi, vert hi, depth).
        let mut stack: Vec<(f64, u32, f64, u32, u32)> = chain
            .windows(2)
            .map(|w| (w[0].0, w[0].1, w[1].0, w[1].1, 0))
            .collect();
        let mut failed = false;
        while let Some((tlo, vlo, thi, vhi, depth)) = stack.pop() {
            let key = if vlo < vhi { [vlo, vhi] } else { [vhi, vlo] };
            if edges.contains(&key) {
                continue;
            }
            if depth >= opts.max_depth || stats.steiner_inserted >= u64::from(opts.max_steiner) {
                failed = true;
                continue;
            }
            // Midpoint Steiner point (exact halving of the parameter;
            // coordinates via f64::midpoint per axis, then snapped onto
            // the ORIGINAL parent segment so later generations cannot
            // walk off a general-position chord).
            let (pa, pb) = (
                tetra.mesh.points[vlo as usize],
                tetra.mesh.points[vhi as usize],
            );
            let mid = snap_to_line(
                [
                    f64::midpoint(pa[0], pb[0]),
                    f64::midpoint(pa[1], pb[1]),
                    f64::midpoint(pa[2], pb[2]),
                ],
                oa,
                ob,
            );
            let bits = [mid[0].to_bits(), mid[1].to_bits(), mid[2].to_bits()];
            let split = if let Some(&twin) = by_bits.get(&bits) {
                // Adopt the existing on-segment vertex.
                Some(twin)
            } else if let Some(twin) = adopt_on_segment(&tetra.mesh.points, oa, ob, mid) {
                // Crossing PLC diagonals after a general-position
                // rotation do not share bitwise midpoints. The first
                // diagonal's centre still lies on the second chord.
                Some(twin)
            } else {
                let new_idx = u32::try_from(tetra.mesh.points.len()).expect("point count fits u32");
                tetra.mesh.points.push(mid);
                if tetra.mesh.insert(new_idx) {
                    stats.steiner_inserted += 1;
                    stats.max_depth_used = stats.max_depth_used.max(depth + 1);
                    by_bits.insert(bits, new_idx);
                    edges = edge_set(tetra);
                    if std::env::var_os("FS_MESH_TRACE_MIDPOINTS").is_some() && edges.contains(&key)
                    {
                        eprintln!(
                            "TRACE recovery: segment {sid} midpoint {new_idx} of sub-edge {key:?} left the edge alive"
                        );
                    }
                    Some(new_idx)
                } else {
                    // A vertex with different stored bits collided in
                    // the kernel's duplicate guard — cannot happen when
                    // the bits index is complete; count honestly.
                    None
                }
            };
            if let Some(v) = split {
                let tmid = f64::midpoint(tlo, thi);
                let pos = chain
                    .binary_search_by(|(t, _)| t.partial_cmp(&tmid).expect("finite"))
                    .unwrap_err();
                chain.insert(pos, (tmid, v));
                stack.push((tlo, vlo, tmid, v, depth + 1));
                stack.push((tmid, v, thi, vhi, depth + 1));
            } else {
                failed = true;
            }
            if stats.steiner_inserted.is_multiple_of(64) {
                cx.checkpoint()?;
            }
        }
        // Verify the finished chain edge-by-edge against the mesh and
        // record the correspondence.
        let mut all_edges = true;
        let sid32 = u32::try_from(sid).expect("segment count fits u32");
        for w in chain.windows(2) {
            let (u, v) = (w[0].1, w[1].1);
            let key = if u < v { [u, v] } else { [v, u] };
            if edges.contains(&key) {
                table.rows.push((key, sid32));
                stats.sub_edges += 1;
            } else {
                all_edges = false;
            }
        }
        if all_edges && !failed {
            stats.recovered += 1;
        } else {
            stats.unrecovered += 1;
            if std::env::var_os("FS_MESH_TRACE_RECOVERY").is_some() {
                let missing: Vec<String> = chain
                    .windows(2)
                    .filter(|w| {
                        let (u, v) = (w[0].1, w[1].1);
                        !edges.contains(&if u < v { [u, v] } else { [v, u] })
                    })
                    .map(|w| format!("[{:.6}..{:.6}]", w[0].0, w[1].0))
                    .collect();
                eprintln!(
                    "TRACE segments: segment {sid} ({a}-{b}) FAILED: chain {} vertices, depth-capped {failed}, missing sub-edges {}: {} (steiner so far {}, duplicates_skipped {}, exhaustive_locates {})",
                    chain.len(),
                    missing.len(),
                    missing.join(" "),
                    stats.steiner_inserted,
                    tetra.mesh.stats.duplicates_skipped,
                    tetra.mesh.stats.exhaustive_locates
                );
                // Anatomy of the first missing sub-edge: how far apart its
                // endpoints are, how many live tets each touches, whether
                // either sits on the hull, and whether any live tet holds
                // both (which would contradict `edges`).
                if let Some(w) = chain.windows(2).find(|w| {
                    let (u, v) = (w[0].1, w[1].1);
                    !edges.contains(&if u < v { [u, v] } else { [v, u] })
                }) {
                    let (u, v) = (w[0].1, w[1].1);
                    let (pu, pv) = (tetra.mesh.points[u as usize], tetra.mesh.points[v as usize]);
                    let dist = ((pu[0] - pv[0]).powi(2)
                        + (pu[1] - pv[1]).powi(2)
                        + (pu[2] - pv[2]).powi(2))
                    .sqrt();
                    let mut touch_u = 0;
                    let mut touch_v = 0;
                    let mut both = 0;
                    let mut hull_u = false;
                    let mut hull_v = false;
                    for (ti, tet) in tetra.mesh.tets.iter().enumerate() {
                        if !tetra.mesh.alive[ti] {
                            continue;
                        }
                        let ghost = tet[3] == GHOST;
                        let has_u = tet.contains(&u);
                        let has_v = tet.contains(&v);
                        if ghost {
                            hull_u |= has_u;
                            hull_v |= has_v;
                            continue;
                        }
                        touch_u += usize::from(has_u);
                        touch_v += usize::from(has_v);
                        both += usize::from(has_u && has_v);
                    }
                    let along = parameter_on_segment(pv, oa, ob)
                        .or_else(|| parameter_on_segment(pu, oa, ob));
                    eprintln!(
                        "TRACE segments:   first missing ({u},{v}): |uv| {dist:.3e} m, live tets touching u {touch_u} v {touch_v} both {both}, hull u {hull_u} v {hull_v}, on-segment parameter {along:?}"
                    );
                }
            }
        }
    }
    Ok((stats, table))
}

/// Facet-recovery evidence.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FacetRecoveryStats {
    /// Facets requested.
    pub facets_in: u64,
    /// Facets fully recovered as face unions.
    pub recovered: u64,
    /// Facets abandoned at a cap (HONESTY counter).
    pub unrecovered: u64,
    /// Steiner points inserted on facets.
    pub steiner_inserted: u64,
    /// Bisection rounds used (worst facet).
    pub rounds_used: u32,
    /// Sub-faces in the correspondence table.
    pub sub_faces: u64,
}

impl FacetRecoveryStats {
    /// Canonical JSON ledger row.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"facets_in\":{},\"recovered\":{},\"unrecovered\":{},\
             \"steiner_inserted\":{},\"rounds_used\":{},\"sub_faces\":{}}}",
            self.facets_in,
            self.recovered,
            self.unrecovered,
            self.steiner_inserted,
            self.rounds_used,
            self.sub_faces
        )
    }
}

/// The facet correspondence: every recovered sub-face (sorted vertex
/// triple) with its parent facet index.
#[derive(Debug, Clone, Default)]
pub struct FacetCorrespondence {
    /// (sub-face, parent facet) rows, deterministic order.
    pub rows: Vec<([u32; 3], u32)>,
}

/// Live mesh face set (sorted vertex triples of live real tets).
/// Live real faces → the apex of each incident real tet (`GHOST` in the
/// second slot of a hull face). Which side of a facet's plane those apexes
/// lie on is what tells a tiling apart from a double cover
/// (`coplanar_tiling`).
type FaceApexes = std::collections::BTreeMap<[u32; 3], [u32; 2]>;

fn face_set(tetra: &Tetrahedralization) -> FaceApexes {
    let mut faces = FaceApexes::new();
    for tet in tetra.tets() {
        for skip in 0..4 {
            let mut f = [0u32; 3];
            let mut j = 0;
            for i in 0..4 {
                if i != skip {
                    f[j] = tet[i];
                    j += 1;
                }
            }
            f.sort_unstable();
            let slots = faces.entry(f).or_insert([GHOST, GHOST]);
            if slots[0] == GHOST {
                slots[0] = tet[skip];
            } else {
                slots[1] = tet[skip];
            }
        }
    }
    faces
}

/// Project a facet loop to 2D by dropping the DOMINANT-normal axis and keeping
/// the two remaining coordinates VERBATIM (exact sub-coordinates, so `orient2d`
/// on them stays exact). Newell's method gives a robust normal for any
/// near-planar simple loop.
/// Newell's polygon normal. Used both to pick the 2-D drop axis and
/// to snap Steiner midpoints back onto the supporting plane.
fn newell_normal(points: &[[f64; 3]], loop_verts: &[u32]) -> [f64; 3] {
    let m = loop_verts.len();
    let mut nrm = [0.0f64; 3];
    for i in 0..m {
        let a = points[loop_verts[i] as usize];
        let b = points[loop_verts[(i + 1) % m] as usize];
        nrm[0] += (a[1] - b[1]) * (a[2] + b[2]);
        nrm[1] += (a[2] - b[2]) * (a[0] + b[0]);
        nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    nrm
}

/// Existing vertex on `ab` closest to `target`, if it is within
/// `1e-12` of the chord length of both the line and the target.
fn adopt_on_segment(
    points: &[[f64; 3]],
    a: [f64; 3],
    b: [f64; 3],
    target: [f64; 3],
) -> Option<u32> {
    let scale2 = (b[0] - a[0]).mul_add(
        b[0] - a[0],
        (b[1] - a[1]).mul_add(b[1] - a[1], (b[2] - a[2]) * (b[2] - a[2])),
    );
    let tol2 = 1e-24 * scale2.max(1.0);
    let mut best: Option<(u32, f64)> = None;
    for (i, p) in points.iter().enumerate() {
        if parameter_on_segment(*p, a, b).is_none() {
            continue;
        }
        let d2 = (p[0] - target[0]).mul_add(
            p[0] - target[0],
            (p[1] - target[1]).mul_add(p[1] - target[1], (p[2] - target[2]) * (p[2] - target[2])),
        );
        if d2 > tol2 {
            continue;
        }
        let better = match best {
            None => true,
            Some((_, bd)) => d2 < bd,
        };
        if better {
            best = Some((u32::try_from(i).expect("point count fits u32"), d2));
        }
    }
    best.map(|(i, _)| i)
}

/// Existing vertex within the segment tolerance (`1e-12` of the chord
/// `ab`) of `target`, closest first — the facet-recovery twin of
/// [`adopt_on_segment`] for midpoints that need not lie on a segment.
fn adopt_near(points: &[[f64; 3]], a: [f64; 3], b: [f64; 3], target: [f64; 3]) -> Option<u32> {
    let scale2 = (b[0] - a[0]).mul_add(
        b[0] - a[0],
        (b[1] - a[1]).mul_add(b[1] - a[1], (b[2] - a[2]) * (b[2] - a[2])),
    );
    let tol2 = 1e-24 * scale2.max(1.0);
    let mut best: Option<(u32, f64)> = None;
    for (i, p) in points.iter().enumerate() {
        let d2 = (p[0] - target[0]).mul_add(
            p[0] - target[0],
            (p[1] - target[1]).mul_add(p[1] - target[1], (p[2] - target[2]) * (p[2] - target[2])),
        );
        if d2 > tol2 {
            continue;
        }
        if best.is_none_or(|(_, bd)| d2 < bd) {
            best = Some((u32::try_from(i).expect("point count fits u32"), d2));
        }
    }
    best.map(|(i, _)| i)
}

/// Parameter `t ∈ (0, 1)` if `p` lies on the open segment `ab`
/// within a relative residual of `1e-12` of the chord length.
fn parameter_on_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> Option<f64> {
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let dd = d[0].mul_add(d[0], d[1].mul_add(d[1], d[2] * d[2]));
    if dd == 0.0 || !dd.is_finite() {
        return None;
    }
    let w = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let t = w[0].mul_add(d[0], w[1].mul_add(d[1], w[2] * d[2])) / dd;
    if t <= 0.0 || t >= 1.0 {
        return None;
    }
    let q = [
        t.mul_add(d[0], a[0]),
        t.mul_add(d[1], a[1]),
        t.mul_add(d[2], a[2]),
    ];
    let res2 = (p[0] - q[0]).mul_add(
        p[0] - q[0],
        (p[1] - q[1]).mul_add(p[1] - q[1], (p[2] - q[2]) * (p[2] - q[2])),
    );
    if res2 > 1e-24 * dd {
        return None;
    }
    Some(t)
}

fn sorted2(a: u32, b: u32) -> [u32; 2] {
    if a < b { [a, b] } else { [b, a] }
}

/// Vertex `i` lies strictly between the chord endpoints `a` and `b` within
/// the segment tolerance. The chord is evaluated with its endpoints in index
/// order so that the two facets sharing a segment — which walk it in
/// opposite directions — reach the same verdict on every vertex near it.
fn on_chord(points: &[[f64; 3]], i: u32, a: u32, b: u32) -> bool {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    parameter_on_segment(points[i as usize], points[lo as usize], points[hi as usize]).is_some()
}

/// Every mesh vertex on the chord `(a, b)` — its endpoints and each vertex
/// within the segment tolerance of it — in parameter order from the lower
/// index to the higher. Consecutive pairs are the sub-edges any tiling that
/// touches this segment must use as its free edges: a chain that skipped a
/// vertex would pass a rounding hair from it, and the facet across the
/// segment, tiled through that vertex, would then not meet this one.
fn chain_on_chord(points: &[[f64; 3]], a: u32, b: u32) -> Vec<u32> {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    let (pa, pb) = (points[lo as usize], points[hi as usize]);
    let d = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
    let pad = 1e-12 * d[0].mul_add(d[0], d[1].mul_add(d[1], d[2] * d[2])).sqrt();
    let min = [
        pa[0].min(pb[0]) - pad,
        pa[1].min(pb[1]) - pad,
        pa[2].min(pb[2]) - pad,
    ];
    let max = [
        pa[0].max(pb[0]) + pad,
        pa[1].max(pb[1]) + pad,
        pa[2].max(pb[2]) + pad,
    ];
    let mut on: Vec<(f64, u32)> = vec![(0.0, lo), (1.0, hi)];
    for (i, p) in points.iter().enumerate() {
        let i = u32::try_from(i).expect("point count fits u32");
        if i == lo
            || i == hi
            || p[0] < min[0]
            || p[0] > max[0]
            || p[1] < min[1]
            || p[1] > max[1]
            || p[2] < min[2]
            || p[2] > max[2]
        {
            continue;
        }
        if let Some(t) = parameter_on_segment(*p, pa, pb) {
            on.push((t, i));
        }
    }
    on.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
    on.into_iter().map(|(_, v)| v).collect()
}

/// The free edges every tiling of the facet must have: the consecutive
/// pairs of each loop edge's chord chain.
fn required_boundary(points: &[[f64; 3]], loop_verts: &[u32]) -> BTreeSet<[u32; 2]> {
    let n = loop_verts.len();
    let mut required = BTreeSet::new();
    for k in 0..n {
        let chain = chain_on_chord(points, loop_verts[k], loop_verts[(k + 1) % n]);
        for w in chain.windows(2) {
            required.insert(sorted2(w[0], w[1]));
        }
    }
    required
}

/// Edges used by exactly one of `tiles`; `None` when some edge is used more
/// than twice (not a tiling of anything).
fn free_edges(tiles: &[[u32; 3]]) -> Option<BTreeSet<[u32; 2]>> {
    let mut edge_use: std::collections::BTreeMap<[u32; 2], u32> = std::collections::BTreeMap::new();
    for t in tiles {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *edge_use.entry(sorted2(x, y)).or_insert(0) += 1;
        }
    }
    let mut free = BTreeSet::new();
    for (edge, uses) in edge_use {
        match uses {
            1 => {
                free.insert(edge);
            }
            2 => {}
            _ => return None,
        }
    }
    Some(free)
}

/// Orthogonal projection onto the line through `a` and `b`. A point
/// already on the line is returned unchanged so axis-aligned bitwise
/// identity and concurrent-segment twins survive.
fn snap_to_line(point: [f64; 3], a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let dd = d[0].mul_add(d[0], d[1].mul_add(d[1], d[2] * d[2]));
    if dd == 0.0 || !dd.is_finite() {
        return point;
    }
    let w = [point[0] - a[0], point[1] - a[1], point[2] - a[2]];
    let t = w[0].mul_add(d[0], w[1].mul_add(d[1], w[2] * d[2])) / dd;
    let q = [
        t.mul_add(d[0], a[0]),
        t.mul_add(d[1], a[1]),
        t.mul_add(d[2], a[2]),
    ];
    if q[0] == point[0] && q[1] == point[1] && q[2] == point[2] {
        return point;
    }
    q
}

/// Orthogonal projection onto the plane through `origin` with
/// unnormalized `normal`. A point already on the plane (`t == 0`) is
/// returned unchanged so axis-aligned bitwise identity survives.
fn snap_to_plane(point: [f64; 3], origin: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let nn = normal[0].mul_add(
        normal[0],
        normal[1].mul_add(normal[1], normal[2] * normal[2]),
    );
    if nn == 0.0 || !nn.is_finite() {
        return point;
    }
    let w0 = point[0] - origin[0];
    let w1 = point[1] - origin[1];
    let w2 = point[2] - origin[2];
    let t = w0.mul_add(normal[0], w1.mul_add(normal[1], w2 * normal[2])) / nn;
    if t == 0.0 {
        return point;
    }
    [
        point[0] - t * normal[0],
        point[1] - t * normal[1],
        point[2] - t * normal[2],
    ]
}

fn project_facet(points: &[[f64; 3]], loop_verts: &[u32]) -> Vec<[f64; 2]> {
    let (u, v) = plane_axes(newell_normal(points, loop_verts));
    loop_verts
        .iter()
        .map(|&idx| {
            let p = points[idx as usize];
            [p[u], p[v]]
        })
        .collect()
}

/// The two coordinate axes a plane with normal `nrm` projects onto: drop
/// the dominant axis (largest |component|) and keep the other two in
/// ascending axis order (a fixed, deterministic choice).
fn plane_axes(nrm: [f64; 3]) -> (usize, usize) {
    let mut ax = 0usize;
    let mut best = nrm[0].abs();
    for (a, &c) in nrm.iter().enumerate().skip(1) {
        if c.abs() > best {
            best = c.abs();
            ax = a;
        }
    }
    match ax {
        0 => (1usize, 2usize),
        1 => (0usize, 2usize),
        _ => (0usize, 1usize),
    }
}

/// True iff the projected simple polygon is convex (every turn one way;
/// collinear vertices allowed). Convex facets keep the exact fan triangulation.
fn is_convex(proj: &[[f64; 2]]) -> bool {
    let n = proj.len();
    let mut sign = 0i8;
    for i in 0..n {
        let a = proj[(i + n - 1) % n];
        let b = proj[i];
        let c = proj[(i + 1) % n];
        match orient2d(a, b, c) {
            Sign::Positive => {
                if sign < 0 {
                    return false;
                }
                sign = 1;
            }
            Sign::Negative => {
                if sign > 0 {
                    return false;
                }
                sign = -1;
            }
            Sign::Zero => {}
        }
    }
    true
}

/// Closed-region containment via exact `orient2d`: is `p` inside triangle
/// `(a, b, c)` (oriented `ccw`), boundary included?
fn in_triangle(p: [f64; 2], a: [f64; 2], b: [f64; 2], c: [f64; 2], ccw: bool) -> bool {
    for (x, y) in [(a, b), (b, c), (c, a)] {
        match orient2d(x, y, p) {
            Sign::Zero => {}
            Sign::Positive => {
                if !ccw {
                    return false;
                }
            }
            Sign::Negative => {
                if ccw {
                    return false;
                }
            }
        }
    }
    true
}

/// Ear-clipping triangulation of a SIMPLE (convex or non-convex) planar facet
/// loop, exact-predicate robust: `orient2d` decides both the convex-corner test
/// and the ear-emptiness test, and the scan clips the FIRST valid ear in index
/// order (deterministic). Returns the triangle vertex triples (original point
/// indices), or `None` if the loop is not a simple non-degenerate polygon (an
/// HONEST failure — the caller counts it `unrecovered`).
fn ear_clip(proj: &[[f64; 2]], loop_verts: &[u32]) -> Option<Vec<[u32; 3]>> {
    let m = loop_verts.len();
    if m < 3 {
        return None;
    }
    // Signed area (shoelace) → orientation; zero area is degenerate.
    let mut area2 = 0.0f64;
    for i in 0..m {
        let a = proj[i];
        let b = proj[(i + 1) % m];
        area2 += a[0].mul_add(b[1], -(b[0] * a[1]));
    }
    if area2 == 0.0 {
        return None;
    }
    let ccw = area2 > 0.0;
    let mut poly: Vec<usize> = (0..m).collect();
    let mut tris: Vec<[u32; 3]> = Vec::with_capacity(m - 2);
    // Each successful pass removes one vertex; `m` passes is the hard ceiling.
    for _ in 0..m {
        if poly.len() == 3 {
            break;
        }
        let n = poly.len();
        let mut clipped = false;
        for i in 0..n {
            let ip = poly[(i + n - 1) % n];
            let ic = poly[i];
            let inx = poly[(i + 1) % n];
            let convex = matches!(
                (orient2d(proj[ip], proj[ic], proj[inx]), ccw),
                (Sign::Positive, true) | (Sign::Negative, false)
            );
            if !convex {
                continue; // reflex or collinear corner is never an ear
            }
            let empty = poly.iter().all(|&k| {
                k == ip
                    || k == ic
                    || k == inx
                    || !in_triangle(proj[k], proj[ip], proj[ic], proj[inx], ccw)
            });
            if empty {
                tris.push([loop_verts[ip], loop_verts[ic], loop_verts[inx]]);
                poly.remove(i);
                clipped = true;
                break;
            }
        }
        if !clipped {
            return None; // no ear found → not a simple polygon
        }
    }
    if poly.len() != 3 {
        return None;
    }
    tris.push([
        loop_verts[poly[0]],
        loop_verts[poly[1]],
        loop_verts[poly[2]],
    ]);
    Some(tris)
}

/// Coplanar mesh faces that TILE the triangular facet `loop_verts`,
/// whichever diagonals the kernel chose.
///
/// The bisection driver above asks for one particular triangulation of
/// the facet (its own sub-triangles). That is a SUFFICIENT condition,
/// not the definition: a facet is recovered when SOME set of mesh faces
/// in its plane tiles it. The two differ exactly on co-circular ties —
/// an axis-aligned rectangle's four corners share a circle, the
/// kernel's symbolic perturbation picks one diagonal, and midpoint
/// bisection only manufactures smaller squares with the same tie, so
/// the old test could starve at the round cap on geometry the mesh had
/// in fact conformed to. Every facet edge is a recovered segment
/// (a chain of mesh edges), so no mesh face in the plane can straddle
/// the facet boundary; hence the faces with all three vertices inside
/// or on the facet tile it iff every edge used by exactly one of them
/// lies on the facet boundary. Vertex membership is by provenance and
/// the crate's segment tolerance (see the classifier below); a vertex
/// that is not recognised simply does not count, which can only delay
/// recognition, never fake it. Non-triangular loops return `None` and
/// keep the sub-triangle path.
/// A vertex IN the facet's plane and STRICTLY inside its triangle, by the
/// same relative tolerance the segment classifier uses (1e-12 of the facet
/// scale). Constrained refinement inserts facet split points that no facet
/// minted; a valid PLC has no other vertex there, so recognising such a
/// point by geometry lets the next recovery pass adopt it into the tiling
/// instead of calling the facet broken.
fn inside_facet_interior(p: [f64; 3], corners: &[[f64; 3]; 3]) -> bool {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let u = sub(corners[1], corners[0]);
    let v = sub(corners[2], corners[0]);
    let w = sub(p, corners[0]);
    let n = cross(u, v);
    let nn = dot(n, n);
    if nn == 0.0 || !nn.is_finite() {
        return false;
    }
    let scale2 = dot(u, u).max(dot(v, v));
    let h = dot(n, w);
    if h * h > 1e-24 * nn * scale2 {
        return false;
    }
    let inv = 1.0 / nn;
    let s = dot(cross(w, v), n) * inv;
    let t = dot(cross(u, w), n) * inv;
    let eps = 1e-12;
    s > eps && t > eps && s + t < 1.0 - eps
}

/// Diagnostic twin of [`coplanar_tiling`] (trace only): lists the faces whose
/// three vertices are geometrically in the facet's plane and within its
/// triangle (loose test), each vertex's classification, and the free edges
/// that stop the tiling from closing.
fn explain_tiling(
    points: &[[f64; 3]],
    faces: &FaceApexes,
    loop_verts: &[u32],
    interior: &BTreeSet<u32>,
) -> String {
    if loop_verts.len() != 3 {
        return "non-triangular loop".to_string();
    }
    let corners = [loop_verts[0], loop_verts[1], loop_verts[2]];
    let cp = [
        points[corners[0] as usize],
        points[corners[1] as usize],
        points[corners[2] as usize],
    ];
    let classify = |i: u32| -> Option<u8> {
        let mut bits = 0u8;
        for k in 0..3 {
            let a = corners[k];
            let b = corners[(k + 1) % 3];
            if i == a
                || i == b
                || parameter_on_segment(points[i as usize], cp[k], cp[(k + 1) % 3]).is_some()
            {
                bits |= 1 << k;
            }
        }
        if bits != 0 {
            Some(bits)
        } else if interior.contains(&i) || inside_facet_interior(points[i as usize], &cp) {
            Some(0)
        } else {
            None
        }
    };
    // Loose geometric membership: in the plane (1e-9 relative) and inside the
    // triangle (barycentric within 1e-6).
    let loose = |i: u32| -> bool {
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let cross = |a: [f64; 3], b: [f64; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let u = sub(cp[1], cp[0]);
        let v = sub(cp[2], cp[0]);
        let w = sub(points[i as usize], cp[0]);
        let n = cross(u, v);
        let nn = dot(n, n);
        let scale2 = dot(u, u).max(dot(v, v));
        let h = dot(n, w);
        if h * h > 1e-18 * nn * scale2 {
            return false;
        }
        let s = dot(cross(w, v), n) / nn;
        let tt = dot(cross(u, w), n) / nn;
        s > -1e-6 && tt > -1e-6 && s + tt < 1.0 + 1e-6
    };
    let mut out = String::new();
    let mut tiles = 0usize;
    let mut edge_use: std::collections::BTreeMap<[u32; 2], (u32, u8)> =
        std::collections::BTreeMap::new();
    for face in faces.keys() {
        if !(loose(face[0]) && loose(face[1]) && loose(face[2])) {
            continue;
        }
        let c = [classify(face[0]), classify(face[1]), classify(face[2])];
        if c.iter().any(Option::is_none) {
            out.push_str(&format!(" face{:?} unclassified {:?};", face, c));
            continue;
        }
        tiles += 1;
        for (x, bx, y, by) in [
            (face[0], c[0].unwrap(), face[1], c[1].unwrap()),
            (face[1], c[1].unwrap(), face[2], c[2].unwrap()),
            (face[0], c[0].unwrap(), face[2], c[2].unwrap()),
        ] {
            let key = if x < y { [x, y] } else { [y, x] };
            let e = edge_use.entry(key).or_insert((0, bx & by));
            e.0 += 1;
        }
    }
    let free: Vec<String> = edge_use
        .iter()
        .filter(|(_, (uses, bits))| *uses == 1 && *bits == 0)
        .map(|(k, _)| format!("{:?}", k))
        .collect();
    let over: Vec<String> = edge_use
        .iter()
        .filter(|(_, (uses, _))| *uses > 2)
        .map(|(k, (uses, _))| format!("{:?}x{uses}", k))
        .collect();
    format!(
        "{tiles} in-plane tiles; free-off-boundary edges {:?}; overused {:?};{out}",
        free, over
    )
}

/// The tiles of the facet's plane on `side`: the in-plane faces with
/// exactly one incident tet whose apex lies strictly on that side. That is
/// the boundary of the tets on that side, restricted to the plane — a
/// manifold sheet by construction, whatever stack of zero-volume tets the
/// kernel left between the two triangulations of each co-circular quad.
fn sheet_on_side(
    tiles: &[([u32; 3], [u32; 2])],
    side: Sign,
    plane_side: &dyn Fn(u32) -> Option<Sign>,
) -> Vec<[u32; 3]> {
    tiles
        .iter()
        .filter(|(_, apexes)| {
            (plane_side(apexes[0]) == Some(side)) != (plane_side(apexes[1]) == Some(side))
        })
        .map(|(face, _)| *face)
        .collect()
}

fn coplanar_tiling(
    points: &[[f64; 3]],
    faces: &FaceApexes,
    loop_verts: &[u32],
    interior: &BTreeSet<u32>,
    required: &BTreeSet<[u32; 2]>,
    trace: bool,
) -> Option<Vec<[u32; 3]>> {
    if loop_verts.len() != 3 {
        return None;
    }
    let corners = [loop_verts[0], loop_verts[1], loop_verts[2]];
    let corner_points = [
        points[corners[0] as usize],
        points[corners[1] as usize],
        points[corners[2] as usize],
    ];
    // Vertex → edge-incidence bits (bit k = on facet edge k, between
    // corner k and corner k+1) for vertices ON the facet; `None` for
    // every other vertex. Membership is by PROVENANCE and the crate's
    // segment tolerance, not by exact predicates: boundary Steiner
    // points on non-dyadic coordinates are collinear with their chord
    // only to an ulp, and an exact `orient2d` would call them outside.
    // Corners are corners; a vertex within the 1e-12 chord residual of
    // edge k (`parameter_on_segment`) is on edge k; a vertex this facet
    // minted on one of its own interior edges is interior; nothing
    // else is on the facet (a valid PLC has no foreign vertex inside a
    // facet, and other facets' Steiner points lie in their own facets).
    let mut cache: std::collections::BTreeMap<u32, Option<u8>> = std::collections::BTreeMap::new();
    let mut classify = |i: u32| -> Option<u8> {
        if let Some(&known) = cache.get(&i) {
            return known;
        }
        let mut bits = 0u8;
        for k in 0..3 {
            let a = corners[k];
            let b = corners[(k + 1) % 3];
            if i == a || i == b || on_chord(points, i, a, b) {
                bits |= 1 << k;
            }
        }
        let verdict = if bits != 0 {
            Some(bits)
        } else if interior.contains(&i) || inside_facet_interior(points[i as usize], &corner_points)
        {
            Some(0)
        } else {
            None
        };
        cache.insert(i, verdict);
        verdict
    };
    let mut tiles: Vec<([u32; 3], [u32; 2])> = Vec::new();
    let mut edge_use: std::collections::BTreeMap<[u32; 2], (u32, u8)> =
        std::collections::BTreeMap::new();
    for (face, apexes) in faces {
        let (Some(b0), Some(b1), Some(b2)) =
            (classify(face[0]), classify(face[1]), classify(face[2]))
        else {
            continue;
        };
        let edges = [
            (face[0], b0, face[1], b1),
            (face[1], b1, face[2], b2),
            (face[0], b0, face[2], b2),
        ];
        // A face with an edge along a facet edge that is not one of that
        // edge's chain sub-edges skips a vertex sitting on the segment (a
        // midpoint a rounding hair off its chord, kept alive beside the
        // edge by a needle tet). The facet across the segment is tiled
        // through that vertex, so this face would leave a hole between
        // the two tilings: it is not a tile of anything.
        if edges
            .iter()
            .any(|&(x, bx, y, by)| bx & by != 0 && !required.contains(&sorted2(x, y)))
        {
            continue;
        }
        tiles.push((*face, *apexes));
        for (x, bx, y, by) in edges {
            let entry = edge_use.entry(sorted2(x, y)).or_insert((0, bx & by));
            entry.0 += 1;
        }
    }
    if tiles.is_empty() {
        return None;
    }
    // Fast path — the in-plane faces are already one tiling: every free
    // edge (used once) lies on ONE facet edge (both endpoints on that
    // edge, inside the facet ⇒ the edge is on it), no edge is used more
    // than twice, and the free edges together are exactly the boundary
    // chains — fewer is a hole on the boundary, a different set a skipped
    // vertex. This is the historical answer, byte-identical.
    let mut free = BTreeSet::new();
    let mut single_cover = true;
    for (&edge, &(uses, shared_bits)) in &edge_use {
        match uses {
            2 => {}
            1 => {
                if shared_bits == 0 {
                    single_cover = false;
                    break;
                }
                free.insert(edge);
            }
            _ => {
                single_cover = false;
                break;
            }
        }
    }
    if single_cover && free == *required {
        return Some(tiles.iter().map(|(face, _)| *face).collect());
    }
    // Otherwise the plane holds both triangulations of some quads (the
    // kernel's zero-volume tets between them, `repair_flat_tets`): take
    // the sheet on one side of the plane — the boundary of that side's
    // tets. Either side is a tiling of a conformed facet; a hull layer has
    // real tets on one side only, so both are tried.
    let normal = {
        let u = [
            corner_points[1][0] - corner_points[0][0],
            corner_points[1][1] - corner_points[0][1],
            corner_points[1][2] - corner_points[0][2],
        ];
        let v = [
            corner_points[2][0] - corner_points[0][0],
            corner_points[2][1] - corner_points[0][1],
            corner_points[2][2] - corner_points[0][2],
        ];
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
        let scale2 =
            (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).max(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        (n, nn * scale2 * 1e-24)
    };
    let plane_side = |q: u32| -> Option<Sign> {
        if q == GHOST {
            return None;
        }
        let p = points[q as usize];
        let (n, tol2) = normal;
        let h = n[0] * (p[0] - corner_points[0][0])
            + n[1] * (p[1] - corner_points[0][1])
            + n[2] * (p[2] - corner_points[0][2]);
        if h * h <= tol2 {
            None
        } else if h > 0.0 {
            Some(Sign::Positive)
        } else {
            Some(Sign::Negative)
        }
    };
    for side in [Sign::Positive, Sign::Negative] {
        let sheet = sheet_on_side(&tiles, side, &plane_side);
        if sheet.is_empty() {
            continue;
        }
        match free_edges(&sheet) {
            Some(free) if free == *required => return Some(sheet),
            outcome => {
                if trace {
                    let free_len = outcome.as_ref().map_or(0, BTreeSet::len);
                    let off: Vec<[u32; 2]> = outcome
                        .as_ref()
                        .map(|f| f.symmetric_difference(required).copied().take(8).collect())
                        .unwrap_or_default();
                    eprintln!(
                        "TRACE sheet: side {side:?}: {} tiles, {} free edges vs {} required (overused: {}); first differences {off:?}{}",
                        sheet.len(),
                        free_len,
                        required.len(),
                        outcome.is_none(),
                        if required.len() <= 20 {
                            format!(
                                "; free {:?} required {:?} tiles {:?}",
                                outcome
                                    .as_ref()
                                    .map(|f| f.iter().copied().collect::<Vec<_>>()),
                                required.iter().copied().collect::<Vec<_>>(),
                                sheet
                            )
                        } else {
                            String::new()
                        }
                    );
                }
            }
        }
    }
    None
}

/// Recover every SIMPLE planar PLC facet (vertex loop into the
/// ORIGINAL points) as a union of mesh faces. Passes repeat until a
/// whole sweep changes nothing (later insertions can undo earlier
/// facets); only the verification against the finished mesh decides.
///
/// # Errors
/// [`MeshError::Cancelled`] between insertions.
///
/// # Panics
/// Only on kernel programmer contracts (and facets with < 3 vertices).
pub fn recover_facets(
    tetra: &mut Tetrahedralization,
    facets: &[Vec<u32>],
    opts: RecoveryOptions,
    cx: &Cx<'_>,
) -> Result<(FacetRecoveryStats, FacetCorrespondence), MeshError> {
    recover_facets_with_points(tetra, facets, &[], &[], opts, cx)
}

/// Split `tris` (a facet's own triangulation) at mesh vertex `v`, which
/// lies in the facet's plane: strictly inside one triangle → that triangle
/// becomes three; on an edge shared by two → both become two. A vertex the
/// tolerant test places in no triangle leaves `tris` alone (the coplanar
/// tiling classifier still recognises it by geometry). Returns whether the
/// triangulation changed.
fn split_tris_at_vertex(points: &[[f64; 3]], tris: &mut Vec<[u32; 3]>, v: u32) -> bool {
    let p = points[v as usize];
    let bary = |t: [u32; 3]| -> Option<(f64, f64)> {
        let a = points[t[0] as usize];
        let b = points[t[1] as usize];
        let c = points[t[2] as usize];
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let q = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
        if nn == 0.0 || !nn.is_finite() {
            return None;
        }
        let cross = |x: [f64; 3], y: [f64; 3]| {
            [
                x[1] * y[2] - x[2] * y[1],
                x[2] * y[0] - x[0] * y[2],
                x[0] * y[1] - x[1] * y[0],
            ]
        };
        let dot = |x: [f64; 3], y: [f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
        Some((dot(cross(q, w), n) / nn, dot(cross(u, q), n) / nn))
    };
    let eps = 1e-9;
    // Strictly inside one triangle.
    for i in 0..tris.len() {
        let t = tris[i];
        if t.contains(&v) {
            return false;
        }
        if let Some((s, r)) = bary(t)
            && s > eps
            && r > eps
            && s + r < 1.0 - eps
        {
            tris.swap_remove(i);
            tris.push([t[0], t[1], v]);
            tris.push([t[1], t[2], v]);
            tris.push([t[2], t[0], v]);
            return true;
        }
    }
    // On an edge: split every triangle that has it.
    for i in 0..tris.len() {
        let t = tris[i];
        let Some((s, r)) = bary(t) else {
            continue;
        };
        if s < -eps || r < -eps || s + r > 1.0 + eps {
            continue;
        }
        let edge = if r.abs() <= eps {
            Some([t[0], t[1]])
        } else if (1.0 - s - r).abs() <= eps {
            Some([t[1], t[2]])
        } else if s.abs() <= eps {
            Some([t[2], t[0]])
        } else {
            None
        };
        let Some(edge) = edge else {
            continue;
        };
        let key = if edge[0] < edge[1] {
            [edge[0], edge[1]]
        } else {
            [edge[1], edge[0]]
        };
        let mut next = Vec::with_capacity(tris.len() + 2);
        let mut changed = false;
        for tt in tris.iter() {
            let has = |x: u32, y: u32| {
                let k = if x < y { [x, y] } else { [y, x] };
                k == key
            };
            if has(tt[0], tt[1]) {
                next.push([tt[0], v, tt[2]]);
                next.push([v, tt[1], tt[2]]);
                changed = true;
            } else if has(tt[1], tt[2]) {
                next.push([tt[0], tt[1], v]);
                next.push([tt[0], v, tt[2]]);
                changed = true;
            } else if has(tt[2], tt[0]) {
                next.push([tt[0], tt[1], v]);
                next.push([v, tt[1], tt[2]]);
                changed = true;
            } else {
                next.push(*tt);
            }
        }
        if changed {
            *tris = next;
            return true;
        }
    }
    false
}

/// INCREMENTAL [`recover_facets`]: a facet whose sub-faces from an earlier
/// recovery are given in `seed_tiles` (`(sorted sub-face, facet index)`,
/// i.e. the previous correspondence rows) starts from THAT triangulation
/// rather than from its loop triangle, and vertices that already lie inside
/// a facet (constrained refinement's wall split points,
/// `(facet index, vertex)`) are stitched into it before the passes begin and
/// counted as interior. Bisection then repairs only the sub-faces an
/// insertion actually destroyed. MEASURED 2026-09-02: restarting from the
/// loop triangle after 23 wall splits on the two-fin comb burned the whole
/// 4,000-point recovery budget and left 30 of 60 facets open.
///
/// # Errors
/// [`MeshError::Cancelled`] between insertions.
#[allow(clippy::too_many_lines)] // one facet-recovery narrative: rounds loop + verification
pub fn recover_facets_with_points(
    tetra: &mut Tetrahedralization,
    facets: &[Vec<u32>],
    points_on_facets: &[(u32, u32)],
    seed_tiles: &[([u32; 3], u32)],
    opts: RecoveryOptions,
    cx: &Cx<'_>,
) -> Result<(FacetRecoveryStats, FacetCorrespondence), MeshError> {
    let mut stats = FacetRecoveryStats {
        facets_in: facets.len() as u64,
        ..FacetRecoveryStats::default()
    };
    let mut table = FacetCorrespondence::default();
    let mut by_bits: std::collections::BTreeMap<[u64; 3], u32> = tetra
        .mesh
        .points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            (
                [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()],
                u32::try_from(i).expect("point count fits u32"),
            )
        })
        .collect();

    /// Per-facet refinement state, persistent across passes.
    struct FacetWork {
        /// The facet's own edge-conforming triangulation (fan or ear
        /// clip, bisected as recovery proceeds).
        tris: Vec<[u32; 3]>,
        /// The original loop edges and every bisection descendant.
        constraint_edges: BTreeSet<[u32; 2]>,
        /// Steiner vertices this facet minted on its own interior edges.
        interior: BTreeSet<u32>,
        origin: [f64; 3],
        normal: [f64; 3],
        rounds: u32,
        failed: bool,
    }

    /// Rows proving the facet is a union of CURRENT mesh faces: its own
    /// sub-triangles when they are all faces (the historical rows,
    /// byte-identical), else any coplanar tiling.
    fn satisfied(
        points: &[[f64; 3]],
        faces: &FaceApexes,
        loop_verts: &[u32],
        work: &FacetWork,
    ) -> Option<Vec<[u32; 3]>> {
        let required = required_boundary(points, loop_verts);
        let own: Vec<[u32; 3]> = work
            .tris
            .iter()
            .map(|t| {
                let mut k = *t;
                k.sort_unstable();
                k
            })
            .collect();
        // The facet's own sub-triangles count only while their boundary is
        // the current chains: a neighbour may have minted a vertex on a
        // shared edge that the mesh kept beside the edge (see
        // `chain_on_chord`), and a tiling that skips it is not closed
        // against the neighbour's.
        if own.iter().all(|k| faces.contains_key(k)) && free_edges(&own).as_ref() == Some(&required)
        {
            return Some(own);
        }
        coplanar_tiling(points, faces, loop_verts, &work.interior, &required, false)
    }

    let mut work: Vec<Option<FacetWork>> = Vec::with_capacity(facets.len());
    for (fid, loop_verts) in facets.iter().enumerate() {
        cx.checkpoint()?;
        assert!(loop_verts.len() >= 3, "facet needs at least 3 vertices");
        // Triangulate the facet: CONVEX → the exact fan (unchanged — cheapest,
        // and keeps the axis-aligned convex path bit-for-bit); NON-CONVEX →
        // exact-predicate ear-clipping in the facet plane (bead iw3l item (a)).
        // A loop that is not a simple non-degenerate polygon is counted
        // `unrecovered`, never faked.
        let proj = project_facet(&tetra.mesh.points, loop_verts);
        let tris: Vec<[u32; 3]> = if is_convex(&proj) {
            (1..loop_verts.len() - 1)
                .map(|i| [loop_verts[0], loop_verts[i], loop_verts[i + 1]])
                .collect()
        } else if let Some(t) = ear_clip(&proj, loop_verts) {
            t
        } else {
            work.push(None);
            continue;
        };
        let origin = tetra.mesh.points[loop_verts[0] as usize];
        let normal = newell_normal(&tetra.mesh.points, loop_verts);
        let mut constraint_edges: BTreeSet<[u32; 2]> = BTreeSet::new();
        let nloop = loop_verts.len();
        for i in 0..nloop {
            let a = loop_verts[i];
            let b = loop_verts[(i + 1) % nloop];
            constraint_edges.insert(if a < b { [a, b] } else { [b, a] });
        }
        let mut interior = BTreeSet::new();
        let mut tris = tris;
        let mut constraint_edges = constraint_edges;
        // Incremental start: the previous recovery's sub-faces of this facet,
        // their boundary edges (used once) as the constraint edges, and every
        // vertex not on the loop boundary as interior.
        let seeded: Vec<[u32; 3]> = seed_tiles
            .iter()
            .filter(|(_, parent)| *parent as usize == fid)
            .map(|(face, _)| *face)
            .collect();
        if !seeded.is_empty() {
            let mut edge_use: std::collections::BTreeMap<[u32; 2], u32> =
                std::collections::BTreeMap::new();
            for t in &seeded {
                for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    *edge_use
                        .entry(if x < y { [x, y] } else { [y, x] })
                        .or_insert(0) += 1;
                }
            }
            constraint_edges = edge_use
                .into_iter()
                .filter(|(_, uses)| *uses == 1)
                .map(|(edge, _)| edge)
                .collect();
            let on_loop = |v: u32| {
                loop_verts.contains(&v)
                    || (0..nloop).any(|i| {
                        parameter_on_segment(
                            tetra.mesh.points[v as usize],
                            tetra.mesh.points[loop_verts[i] as usize],
                            tetra.mesh.points[loop_verts[(i + 1) % nloop] as usize],
                        )
                        .is_some()
                    })
            };
            for t in &seeded {
                for &v in t {
                    if !on_loop(v) {
                        interior.insert(v);
                    }
                }
            }
            tris = seeded;
        }
        for &(facet_id, v) in points_on_facets {
            if facet_id as usize == fid {
                split_tris_at_vertex(&tetra.mesh.points, &mut tris, v);
                interior.insert(v);
            }
        }
        work.push(Some(FacetWork {
            tris,
            constraint_edges,
            interior,
            origin,
            normal,
            rounds: 0,
            failed: false,
        }));
    }

    // PASSES. One pass gives every unsatisfied facet one batch round of
    // longest-edge bisection against the current mesh. A single sweep
    // is not a proof: a later insertion re-tetrahedralizes its cavity
    // and can flip away faces an earlier facet had already conformed
    // to, so sweeps repeat until every facet is satisfied or no facet
    // can take another round (the per-facet round cap and the Steiner
    // budget are the honest caps; MEASURED: a four-fin comb needed 13+
    // passes but only 8 rounds on its worst facet), and only the
    // verification against the FINISHED mesh below decides.
    // Trace bookkeeping: where the Steiner points went.
    let mut constraint_splits = 0u64;
    let mut interior_splits = 0u64;
    'passes: loop {
        cx.checkpoint()?;
        let faces = face_set(tetra);
        let pending: Vec<usize> = work
            .iter()
            .enumerate()
            .filter_map(|(fid, w)| {
                let w = w.as_ref()?;
                if w.failed || satisfied(&tetra.mesh.points, &faces, &facets[fid], w).is_some() {
                    None
                } else {
                    Some(fid)
                }
            })
            .collect();
        if pending.is_empty() {
            break;
        }
        if std::env::var_os("FS_MESH_TRACE_RECOVERY").is_some() {
            let detail: Vec<String> = pending
                .iter()
                .take(6)
                .map(|&fid| {
                    let w = work[fid].as_ref().expect("pending facet has work");
                    let missing = w
                        .tris
                        .iter()
                        .filter(|t| {
                            let mut k = **t;
                            k.sort_unstable();
                            !faces.contains_key(&k)
                        })
                        .count();
                    format!(
                        "f{fid}(r{} tris{} miss{} int{} pts{} tile:{})",
                        w.rounds,
                        w.tris.len(),
                        missing,
                        w.interior.len(),
                        points_on_facets
                            .iter()
                            .filter(|(f, _)| *f as usize == fid)
                            .count(),
                        coplanar_tiling(
                            &tetra.mesh.points,
                            &faces,
                            &facets[fid],
                            &w.interior,
                            &required_boundary(&tetra.mesh.points, &facets[fid]),
                            false
                        )
                        .map_or("none", |_| "found")
                    )
                })
                .collect();
            eprintln!(
                "TRACE recovery pass: pending {} steiner {} (constraint {} interior {}) :: {}",
                pending.len(),
                stats.steiner_inserted,
                constraint_splits,
                interior_splits,
                detail.join(" ")
            );
            if let Some(&fid) = pending
                .iter()
                .find(|&&fid| points_on_facets.iter().any(|(f, _)| *f as usize == fid))
            {
                let w = work[fid].as_ref().expect("pending facet has work");
                eprintln!(
                    "TRACE recovery explain f{fid}: {}",
                    explain_tiling(&tetra.mesh.points, &faces, &facets[fid], &w.interior)
                );
            }
        }
        // A pass that gives no facet a round has nothing left to try
        // (every pending facet is capped); stop rather than spin. Total
        // passes are bounded by facets × `max_depth` rounds.
        let mut progressed = false;
        for fid in pending {
            let loop_verts = &facets[fid];
            let Some(w) = work[fid].as_mut() else {
                continue;
            };
            let faces = face_set(tetra);
            if satisfied(&tetra.mesh.points, &faces, loop_verts, w).is_some() {
                continue; // a neighbour's round conformed it meanwhile
            }
            if w.rounds >= opts.max_depth {
                w.failed = true;
                if std::env::var_os("FS_MESH_TRACE_RECOVERY").is_some() {
                    let faces = face_set(tetra);
                    let missing = w
                        .tris
                        .iter()
                        .filter(|t| {
                            let mut k = **t;
                            k.sort_unstable();
                            !faces.contains_key(&k)
                        })
                        .count();
                    eprintln!(
                        "TRACE recovery: facet {fid} FAILED at round cap {}: {} own tris, {} missing, {} interior points, {} stitched split points, tiling {}",
                        w.rounds,
                        w.tris.len(),
                        missing,
                        w.interior.len(),
                        points_on_facets
                            .iter()
                            .filter(|(f, _)| *f as usize == fid)
                            .count(),
                        coplanar_tiling(
                            &tetra.mesh.points,
                            &faces,
                            loop_verts,
                            &w.interior,
                            &required_boundary(&tetra.mesh.points, loop_verts),
                            false
                        )
                        .map_or("none", |_| "found")
                    );
                }
                continue;
            }
            if stats.steiner_inserted >= u64::from(opts.max_steiner) {
                w.failed = true;
                break 'passes;
            }
            w.rounds += 1;
            progressed = true;
            stats.rounds_used = stats.rounds_used.max(w.rounds);
            // Own sub-triangles to refine: those that are not mesh faces,
            // and those that ARE faces but skip a vertex a neighbour put on
            // one of this facet's edges (a midpoint the mesh kept beside its
            // edge in a needle tet, see `chain_on_chord`): such a face exists
            // but no tiling through it closes against the neighbour's.
            // Longest-edge bisection then adopts the vertex when its edge
            // comes up — raw midpoints of dyadic chains coincide bitwise and
            // `adopt_near` covers the ulp — while keeping the sub-triangles
            // balanced (MEASURED 2026-09-02: fanning a triangle over a long
            // chain instead made interior splits run away, 3393 of 3703
            // Steiner points on the four-fin comb that needs 136).
            let chains: Vec<std::collections::BTreeMap<u32, usize>> = (0..loop_verts.len())
                .map(|k| {
                    chain_on_chord(
                        &tetra.mesh.points,
                        loop_verts[k],
                        loop_verts[(k + 1) % loop_verts.len()],
                    )
                    .into_iter()
                    .enumerate()
                    .map(|(i, v)| (v, i))
                    .collect()
                })
                .collect();
            let skips_chain = |t: &[u32; 3]| -> bool {
                chains.iter().any(|position| {
                    [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]
                        .iter()
                        .any(|&(x, y)| {
                            matches!(
                                (position.get(&x), position.get(&y)),
                                (Some(&i), Some(&j)) if i.abs_diff(j) > 1
                            )
                        })
                })
            };
            let missing: Vec<usize> = w
                .tris
                .iter()
                .enumerate()
                .filter(|(_, t)| {
                    let mut k = **t;
                    k.sort_unstable();
                    !faces.contains_key(&k) || skips_chain(t)
                })
                .map(|(i, _)| i)
                .collect();
            // Batch: the LONGEST edge of EVERY missing triangle
            // (deterministic: BTreeSet of sorted pairs), split at
            // midpoints this round. One-split-per-round was MEASURED
            // to starve at the rounds cap (12 rounds, facet still
            // open); batching converges in a handful of rounds.
            let mut split_edges: BTreeSet<[u32; 2]> = BTreeSet::new();
            for &mi in &missing {
                let t = w.tris[mi];
                let pts = &tetra.mesh.points;
                let mut best: Option<([u32; 2], f64)> = None;
                for (u, v) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    let (pu, pv) = (pts[u as usize], pts[v as usize]);
                    let d2 =
                        (pu[0] - pv[0]).powi(2) + (pu[1] - pv[1]).powi(2) + (pu[2] - pv[2]).powi(2);
                    let key = if u < v { [u, v] } else { [v, u] };
                    let better = match &best {
                        None => true,
                        Some((bk, bd)) => d2 > *bd || (d2.to_bits() == bd.to_bits() && key < *bk),
                    };
                    if better {
                        best = Some((key, d2));
                    }
                }
                split_edges.insert(best.expect("triangle has edges").0);
            }
            for [u, v] in split_edges {
                if stats.steiner_inserted >= u64::from(opts.max_steiner) {
                    w.failed = true;
                    break 'passes;
                }
                let (pu, pv) = (tetra.mesh.points[u as usize], tetra.mesh.points[v as usize]);
                let raw = [
                    f64::midpoint(pu[0], pv[0]),
                    f64::midpoint(pu[1], pv[1]),
                    f64::midpoint(pu[2], pv[2]),
                ];
                let edge = if u < v { [u, v] } else { [v, u] };
                let on_constraint = w.constraint_edges.contains(&edge);
                let mid = if on_constraint {
                    // Shared crease edges must agree across parent
                    // planes; the raw midpoint of two constraint
                    // vertices is the twin key.
                    raw
                } else {
                    snap_to_plane(raw, w.origin, w.normal)
                };
                let bits = [mid[0].to_bits(), mid[1].to_bits(), mid[2].to_bits()];
                let m = if let Some(&twin) = by_bits.get(&bits) {
                    twin
                } else if let Some(near) = adopt_near(&tetra.mesh.points, pu, pv, mid) {
                    // The same geometric point reached by a different
                    // bisection path (e.g. a segment-recovery midpoint of
                    // non-dyadic endpoints) differs by an ulp; minting a
                    // twin an ulp away breeds sliver tets and breaks the
                    // boundary chain, so adopt it exactly as segment
                    // recovery does.
                    near
                } else {
                    let new_idx =
                        u32::try_from(tetra.mesh.points.len()).expect("point count fits u32");
                    tetra.mesh.points.push(mid);
                    if tetra.mesh.insert(new_idx) {
                        stats.steiner_inserted += 1;
                        by_bits.insert(bits, new_idx);
                        if on_constraint {
                            constraint_splits += 1;
                        } else {
                            interior_splits += 1;
                        }
                        // Rebuilding the edge set per insertion is O(tets): its
                        // own switch, so the pass trace stays usable on big runs.
                        if std::env::var_os("FS_MESH_TRACE_MIDPOINTS").is_some()
                            && edge_set(tetra).contains(&edge)
                        {
                            eprintln!(
                                "TRACE recovery: facet {fid} midpoint {new_idx} of edge {edge:?} (constraint {on_constraint}) left the edge alive"
                            );
                        }
                        new_idx
                    } else {
                        w.failed = true;
                        break;
                    }
                };
                if on_constraint {
                    w.constraint_edges
                        .insert(if u < m { [u, m] } else { [m, u] });
                    w.constraint_edges
                        .insert(if v < m { [v, m] } else { [m, v] });
                } else {
                    w.interior.insert(m);
                }
                // Split EVERY facet triangle sharing edge (u, v) so
                // the facet triangulation stays edge-conforming.
                let mut next: Vec<[u32; 3]> = Vec::with_capacity(w.tris.len() + 2);
                for tt in &w.tris {
                    if tt.contains(&u) && tt.contains(&v) {
                        let third = *tt
                            .iter()
                            .find(|&&x| x != u && x != v)
                            .expect("third vertex");
                        next.push([u, m, third]);
                        next.push([m, v, third]);
                    } else {
                        next.push(*tt);
                    }
                }
                w.tris = next;
            }
            cx.checkpoint()?;
        }
        if !progressed {
            break;
        }
    }

    // Verify against the FINISHED mesh and record correspondence.
    let faces = face_set(tetra);
    if std::env::var_os("FS_MESH_TRACE_RECOVERY").is_some() {
        let report = tetra.audit(false);
        eprintln!(
            "TRACE recovery: finished mesh: {} points, kernel stats {}, exact audit violations {} (first: {:?})",
            tetra.mesh.points.len(),
            tetra.mesh.stats.to_json(),
            report.violations.len(),
            report.violations.first()
        );
        if let Some(path) = std::env::var_os("FS_MESH_DUMP_MESH") {
            // Debug aid: the finished mesh as text (`p x y z` per point in
            // index order, `t a b c d` per live real tet, `f v0 v1 v2 ...`
            // per facet loop) for offline analysis of a failed recovery.
            let mut out = String::new();
            for p in &tetra.mesh.points {
                out.push_str(&format!("p {:?} {:?} {:?}\n", p[0], p[1], p[2]));
            }
            for tet in tetra.tets() {
                if tet[3] != GHOST {
                    out.push_str(&format!("t {} {} {} {}\n", tet[0], tet[1], tet[2], tet[3]));
                }
            }
            for loop_verts in facets {
                let ids: Vec<String> = loop_verts.iter().map(u32::to_string).collect();
                out.push_str(&format!("f {}\n", ids.join(" ")));
            }
            if let Err(err) = std::fs::write(&path, out) {
                eprintln!("TRACE recovery: mesh dump to {path:?} failed: {err}");
            }
        }
    }
    for (fid, loop_verts) in facets.iter().enumerate() {
        let fid32 = u32::try_from(fid).expect("facet count fits u32");
        let rows = work[fid]
            .as_ref()
            .and_then(|w| satisfied(&tetra.mesh.points, &faces, loop_verts, w));
        match rows {
            Some(rows) => {
                for k in rows {
                    table.rows.push((k, fid32));
                    stats.sub_faces += 1;
                }
                stats.recovered += 1;
            }
            None => {
                stats.unrecovered += 1;
                if std::env::var_os("FS_MESH_TRACE_RECOVERY").is_some() {
                    let corners: Vec<[f64; 3]> = loop_verts
                        .iter()
                        .map(|&v| tetra.mesh.points[v as usize])
                        .collect();
                    let (rounds, tris, missing, interior) =
                        work[fid].as_ref().map_or((0, 0, 0, 0), |w| {
                            let missing = w
                                .tris
                                .iter()
                                .filter(|t| {
                                    let mut k = **t;
                                    k.sort_unstable();
                                    !faces.contains_key(&k)
                                })
                                .count();
                            (w.rounds, w.tris.len(), missing, w.interior.len())
                        });
                    let explain = work[fid].as_ref().map_or_else(String::new, |w| {
                        explain_tiling(&tetra.mesh.points, &faces, loop_verts, &w.interior)
                    });
                    let explain: String = explain.chars().take(700).collect();
                    eprintln!(
                        "TRACE recovery: facet {fid} UNRECOVERED {loop_verts:?} corners {corners:?}: rounds {rounds} own tris {tris} missing {missing} interior {interior}; {explain}"
                    );
                    if let Some(w) = work[fid].as_ref() {
                        // The sheet extraction's own account of the failure.
                        let _ = coplanar_tiling(
                            &tetra.mesh.points,
                            &faces,
                            loop_verts,
                            &w.interior,
                            &required_boundary(&tetra.mesh.points, loop_verts),
                            true,
                        );
                    }
                }
            }
        }
    }
    Ok((stats, table))
}

#[cfg(test)]
mod tests {
    use super::{ear_clip, is_convex, project_facet, sheet_on_side};
    use crate::delaunay::GHOST;
    use fs_ivl::Sign;

    /// The sheet on one side of a facet's plane is the boundary of that
    /// side's tets: a square carrying both triangulations (a zero-volume
    /// tet between them, apexes in the plane) contributes exactly the
    /// diagonal facing that side, a face with real tets on both sides is in
    /// both sheets, and a hull face is in the sheet of its one real tet.
    #[test]
    fn sheet_on_side_is_the_boundary_of_that_sides_tets() {
        // Vertices 0..=8 in the plane; 9 above, 10 below.
        let plane_side = |q: u32| -> Option<Sign> {
            match q {
                9 => Some(Sign::Positive),
                10 => Some(Sign::Negative),
                _ => None,
            }
        };
        let tiles: Vec<([u32; 3], [u32; 2])> = vec![
            ([5, 6, 7], [9, 8]),     // top of the flat [5,6,7,8]
            ([5, 7, 8], [9, 6]),     // top of the flat
            ([5, 6, 8], [7, 10]),    // bottom of the flat
            ([6, 7, 8], [5, 10]),    // bottom of the flat
            ([0, 3, 5], [9, 10]),    // ordinary sheet face
            ([1, 2, 7], [9, GHOST]), // hull face, real tet above
        ];
        let above = sheet_on_side(&tiles, Sign::Positive, &plane_side);
        assert_eq!(above, vec![[5, 6, 7], [5, 7, 8], [0, 3, 5], [1, 2, 7]]);
        let below = sheet_on_side(&tiles, Sign::Negative, &plane_side);
        assert_eq!(below, vec![[5, 6, 8], [6, 7, 8], [0, 3, 5]]);
    }

    fn poly_area(p: &[[f64; 2]]) -> f64 {
        let m = p.len();
        let mut a = 0.0;
        for i in 0..m {
            let u = p[i];
            let v = p[(i + 1) % m];
            a += u[0].mul_add(v[1], -(v[0] * u[1]));
        }
        a.abs() * 0.5
    }

    fn tri_area(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
        ((b[0] - a[0]).mul_add(c[1] - a[1], -((c[0] - a[0]) * (b[1] - a[1])))).abs() * 0.5
    }

    /// Ear-clipping tiles simple polygons EXACTLY: `m − 2` triangles built only
    /// from original vertices whose areas sum to the polygon area (no overlap,
    /// no gap, nothing outside the boundary). Covers a non-convex, NON-star
    /// U-shape (the fan triangulation is wrong here — ear-clipping is required),
    /// a non-convex star-shaped L, and a convex pentagon.
    #[test]
    fn ear_clip_tiles_simple_polygons() {
        let u_shape = vec![
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [2.0, 3.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 3.0],
            [0.0, 3.0],
        ];
        let l_shape = vec![
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ];
        let convex = vec![[0.0, 0.0], [2.0, 0.0], [3.0, 1.0], [1.5, 2.5], [0.0, 1.5]];
        assert!(!is_convex(&u_shape), "U-shape is non-convex");
        assert!(!is_convex(&l_shape), "L-shape is non-convex");
        assert!(is_convex(&convex), "pentagon is convex");
        for poly in [&u_shape, &l_shape, &convex] {
            let m = poly.len();
            let loop_verts: Vec<u32> = (0..m as u32).collect();
            let tris = ear_clip(poly, &loop_verts).expect("simple polygon triangulates");
            assert_eq!(tris.len(), m - 2, "a simple m-gon yields m − 2 triangles");
            let want = poly_area(poly);
            let got: f64 = tris
                .iter()
                .map(|t| {
                    tri_area(
                        poly[t[0] as usize],
                        poly[t[1] as usize],
                        poly[t[2] as usize],
                    )
                })
                .sum();
            assert!(
                (got - want).abs() < 1e-12,
                "triangulation area {got} != polygon area {want}"
            );
            // Only original loop vertices — ear-clipping adds no Steiner points.
            assert!(tris.iter().flatten().all(|&v| (v as usize) < m));
        }
        // A degenerate all-collinear loop is refused, not faked.
        let line = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]];
        let lv: Vec<u32> = (0..4).collect();
        assert!(
            ear_clip(&line, &lv).is_none(),
            "collinear loop is not a polygon"
        );
    }

    /// The 3D projection keeps the dominant-normal axis dropped so the loop
    /// stays a non-degenerate 2D polygon, and convexity survives projection.
    #[test]
    fn project_facet_drops_dominant_axis() {
        // A convex quad in the z = 0.5 plane → projects to the (x, y) plane.
        let points = vec![
            [0.0, 0.0, 0.5],
            [1.0, 0.0, 0.5],
            [1.0, 1.0, 0.5],
            [0.0, 1.0, 0.5],
        ];
        let proj = project_facet(&points, &[0, 1, 2, 3]);
        assert_eq!(proj, vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        assert!(is_convex(&proj));
        // A quad in the x = 2 plane → dominant axis is x, projects to (y, z).
        let yz = vec![
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [2.0, 2.0, 2.0],
            [2.0, 0.0, 2.0],
        ];
        let projx = project_facet(&yz, &[0, 1, 2, 3]);
        assert_eq!(projx, vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
    }
}
