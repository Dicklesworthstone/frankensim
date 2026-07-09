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
//! (`unrecovered`), never silently dropped.
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
//! never faked. General-position (non-axis-aligned) planes remain the
//! recorded successor: f64 midpoints stay EXACTLY coplanar only when
//! the plane is axis-aligned (bitwise-equal coordinate), so the
//! correspondence verification measures the residual rather than
//! assuming it.

use crate::delaunay::{GHOST, MeshError, Tetrahedralization};
use fs_exec::Cx;
use fs_ivl::{Sign, orient2d};
use std::collections::BTreeSet;

/// Recovery policy.
#[derive(Debug, Clone, Copy)]
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
#[derive(Debug, Clone, Copy, Default)]
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
        let mut chain: Vec<(f64, u32)> = vec![(0.0, a), (1.0, b)];
        // Work stack of open sub-intervals (param lo, vert lo, param
        // hi, vert hi, depth).
        let mut stack: Vec<(f64, u32, f64, u32, u32)> = vec![(0.0, a, 1.0, b, 0)];
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
            // coordinates via f64::midpoint per axis).
            let (pa, pb) = (
                tetra.mesh.points[vlo as usize],
                tetra.mesh.points[vhi as usize],
            );
            let mid = [
                f64::midpoint(pa[0], pb[0]),
                f64::midpoint(pa[1], pb[1]),
                f64::midpoint(pa[2], pb[2]),
            ];
            let bits = [mid[0].to_bits(), mid[1].to_bits(), mid[2].to_bits()];
            let split = if let Some(&twin) = by_bits.get(&bits) {
                // Adopt the existing on-segment vertex.
                Some(twin)
            } else {
                let new_idx = u32::try_from(tetra.mesh.points.len()).expect("point count fits u32");
                tetra.mesh.points.push(mid);
                if tetra.mesh.insert(new_idx) {
                    stats.steiner_inserted += 1;
                    stats.max_depth_used = stats.max_depth_used.max(depth + 1);
                    by_bits.insert(bits, new_idx);
                    edges = edge_set(tetra);
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
        }
    }
    Ok((stats, table))
}

/// Facet-recovery evidence.
#[derive(Debug, Clone, Copy, Default)]
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
fn face_set(tetra: &Tetrahedralization) -> BTreeSet<[u32; 3]> {
    let mut faces = BTreeSet::new();
    for tet in tetra.tets() {
        for skip in 0..4 {
            let mut f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| tet[i]).collect();
            if f.contains(&GHOST) {
                continue;
            }
            f.sort_unstable();
            faces.insert([f[0], f[1], f[2]]);
        }
    }
    faces
}

/// Project a facet loop to 2D by dropping the DOMINANT-normal axis and keeping
/// the two remaining coordinates VERBATIM (exact sub-coordinates, so `orient2d`
/// on them stays exact). Newell's method gives a robust normal for any
/// near-planar simple loop.
fn project_facet(points: &[[f64; 3]], loop_verts: &[u32]) -> Vec<[f64; 2]> {
    let m = loop_verts.len();
    let mut nrm = [0.0f64; 3];
    for i in 0..m {
        let a = points[loop_verts[i] as usize];
        let b = points[loop_verts[(i + 1) % m] as usize];
        nrm[0] += (a[1] - b[1]) * (a[2] + b[2]);
        nrm[1] += (a[2] - b[2]) * (a[0] + b[0]);
        nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    // Dominant axis = largest |component|; keep the other two in ascending
    // axis order (a fixed, deterministic choice).
    let mut ax = 0usize;
    let mut best = nrm[0].abs();
    for (a, &c) in nrm.iter().enumerate().skip(1) {
        if c.abs() > best {
            best = c.abs();
            ax = a;
        }
    }
    let (u, v) = match ax {
        0 => (1usize, 2usize),
        1 => (0usize, 2usize),
        _ => (0usize, 1usize),
    };
    loop_verts
        .iter()
        .map(|&idx| {
            let p = points[idx as usize];
            [p[u], p[v]]
        })
        .collect()
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
    tris.push([loop_verts[poly[0]], loop_verts[poly[1]], loop_verts[poly[2]]]);
    Some(tris)
}

/// Recover every SIMPLE planar PLC facet (vertex loop into the
/// ORIGINAL points) as a union of mesh faces.
///
/// # Errors
/// [`MeshError::Cancelled`] between insertions.
///
/// # Panics
/// Only on kernel programmer contracts (and facets with < 3 vertices).
#[allow(clippy::too_many_lines)] // one facet-recovery narrative: rounds loop + verification
pub fn recover_facets(
    tetra: &mut Tetrahedralization,
    facets: &[Vec<u32>],
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
    for (fid, loop_verts) in facets.iter().enumerate() {
        cx.checkpoint()?;
        assert!(loop_verts.len() >= 3, "facet needs at least 3 vertices");
        // Triangulate the facet: CONVEX → the exact fan (unchanged — cheapest,
        // and keeps the axis-aligned convex path bit-for-bit); NON-CONVEX →
        // exact-predicate ear-clipping in the facet plane (bead iw3l item (a)).
        // A loop that is not a simple non-degenerate polygon is counted
        // `unrecovered`, never faked.
        let proj = project_facet(&tetra.mesh.points, loop_verts);
        let mut tris: Vec<[u32; 3]> = if is_convex(&proj) {
            (1..loop_verts.len() - 1)
                .map(|i| [loop_verts[0], loop_verts[i], loop_verts[i + 1]])
                .collect()
        } else if let Some(t) = ear_clip(&proj, loop_verts) {
            t
        } else {
            stats.unrecovered += 1;
            continue;
        };
        let mut failed = false;
        let mut round = 0u32;
        loop {
            let faces = face_set(tetra);
            let missing: Vec<usize> = tris
                .iter()
                .enumerate()
                .filter(|(_, t)| {
                    let mut k = **t;
                    k.sort_unstable();
                    !faces.contains(&k)
                })
                .map(|(i, _)| i)
                .collect();
            if missing.is_empty() {
                break;
            }
            if round >= opts.max_depth {
                failed = true;
                break;
            }
            if stats.steiner_inserted >= u64::from(opts.max_steiner) {
                failed = true;
                break;
            }
            round += 1;
            // Batch: the LONGEST edge of EVERY missing triangle
            // (deterministic: BTreeSet of sorted pairs), split at
            // midpoints this round. One-split-per-round was MEASURED
            // to starve at the rounds cap (12 rounds, facet still
            // open); batching converges in a handful of rounds.
            let mut split_edges: BTreeSet<[u32; 2]> = BTreeSet::new();
            for &mi in &missing {
                let t = tris[mi];
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
                    failed = true;
                    break;
                }
                let (pu, pv) = (tetra.mesh.points[u as usize], tetra.mesh.points[v as usize]);
                let mid = [
                    f64::midpoint(pu[0], pv[0]),
                    f64::midpoint(pu[1], pv[1]),
                    f64::midpoint(pu[2], pv[2]),
                ];
                let bits = [mid[0].to_bits(), mid[1].to_bits(), mid[2].to_bits()];
                let m = if let Some(&twin) = by_bits.get(&bits) {
                    twin
                } else {
                    let new_idx =
                        u32::try_from(tetra.mesh.points.len()).expect("point count fits u32");
                    tetra.mesh.points.push(mid);
                    if tetra.mesh.insert(new_idx) {
                        stats.steiner_inserted += 1;
                        by_bits.insert(bits, new_idx);
                        new_idx
                    } else {
                        failed = true;
                        break;
                    }
                };
                // Split EVERY facet triangle sharing edge (u, v) so
                // the facet triangulation stays edge-conforming.
                let mut next: Vec<[u32; 3]> = Vec::with_capacity(tris.len() + 2);
                for tt in &tris {
                    if tt.contains(&u) && tt.contains(&v) {
                        let w = *tt
                            .iter()
                            .find(|&&x| x != u && x != v)
                            .expect("third vertex");
                        next.push([u, m, w]);
                        next.push([m, v, w]);
                    } else {
                        next.push(*tt);
                    }
                }
                tris = next;
            }
            if failed {
                break;
            }
            stats.rounds_used = stats.rounds_used.max(round);
            cx.checkpoint()?;
        }
        // Verify against the finished mesh and record correspondence.
        let faces = face_set(tetra);
        let fid32 = u32::try_from(fid).expect("facet count fits u32");
        let mut all_faces = true;
        for t in &tris {
            let mut k = *t;
            k.sort_unstable();
            if faces.contains(&k) {
                table.rows.push((k, fid32));
                stats.sub_faces += 1;
            } else {
                all_faces = false;
            }
        }
        if all_faces && !failed {
            stats.recovered += 1;
        } else {
            stats.unrecovered += 1;
        }
    }
    Ok((stats, table))
}
