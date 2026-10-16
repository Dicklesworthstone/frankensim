//! The h-vs-p DECISION signal (mechanism 3 of 4): where the solution
//! is locally smooth, raising the polynomial order beats splitting
//! cells; where it kinks or layers, h wins. The classifier is the
//! scaled second-to-first-derivative ratio `s_K = h·|H|_F / (|∇u| + δ)`
//! — large means the local Taylor series is not converging at this
//! resolution (h-refine), small means it is (p-enrich).
//!
//! EXECUTION of local p-enrichment awaits the high-order FEEC element
//! families (recorded no-claim); this module emits the decisions the
//! DWR loop will route when they land.

use crate::estimate::q1;
use fs_cutfem::Quadtree;
use std::collections::BTreeMap;

/// The per-cell routing decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Decision {
    /// Split the cell (non-smooth: kink/layer at this resolution).
    HRefine,
    /// Raise the local order (smooth: spectral convergence available).
    PEnrich,
}

/// Classify every leaf of a UNIFORM grid.
///
/// # Panics
/// If the grid is not uniform.
#[must_use]
pub fn h_vs_p(
    grid: &Quadtree,
    nodal: &BTreeMap<(u32, u32), f64>,
    threshold: f64,
) -> BTreeMap<(u32, u32, u32), Decision> {
    let level = grid.leaves().next().expect("nonempty grid").0;
    assert!(
        grid.leaves().all(|c| c.0 == level),
        "h-vs-p v1 needs a uniform grid"
    );
    let h = 1.0 / f64::from(1u32 << level);
    let s = 1u32 << (grid.max_level() - level);
    let ext = grid.node_extent();
    let val = |gi: i64, gj: i64| -> f64 {
        let gi = gi.clamp(0, i64::from(ext)) as u32;
        let gj = gj.clamp(0, i64::from(ext)) as u32;
        nodal.get(&(gi, gj)).copied().unwrap_or(0.0)
    };
    let mut out = BTreeMap::new();
    for c in grid.leaves() {
        let (lo, hi) = grid.rect(c);
        let corners = grid.corner_nodes(c);
        // Gradient at the cell center from the bilinear.
        let center = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];
        let (_, gr) = q1(lo, hi, center);
        let mut g = [0.0f64; 2];
        for a in 0..4 {
            let v = nodal.get(&corners[a]).copied().unwrap_or(0.0);
            g[0] += gr[a][0] * v;
            g[1] += gr[a][1] * v;
        }
        let gnorm = g[0].hypot(g[1]);
        // Second differences across the cell (same stencil as aniso).
        let (gi, gj) = (i64::from(corners[0].0), i64::from(corners[0].1));
        let st = i64::from(s);
        let hxx = (val(gi + 2 * st, gj) - 2.0 * val(gi + st, gj) + val(gi, gj)) / (h * h);
        let hyy = (val(gi, gj + 2 * st) - 2.0 * val(gi, gj + st) + val(gi, gj)) / (h * h);
        let hxy =
            (val(gi + st, gj + st) - val(gi + st, gj) - val(gi, gj + st) + val(gi, gj)) / (h * h);
        let hf = (hxx * hxx + hyy * hyy + 2.0 * hxy * hxy).sqrt();
        let s_k = h * hf / (gnorm + 1e-12);
        out.insert(
            c,
            if s_k > threshold {
                Decision::HRefine
            } else {
                Decision::PEnrich
            },
        );
    }
    out
}
