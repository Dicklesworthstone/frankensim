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
//! (`unrecovered`), never silently dropped. CONSTRAINED-Delaunay
//! facet recovery (interior/non-convex facets) is the recorded
//! successor; convex hull-facet conformity is gated test-side.

use crate::delaunay::{GHOST, MeshError, Tetrahedralization};
use fs_exec::Cx;
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
