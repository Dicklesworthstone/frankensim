//! Anisotropic metric synthesis (mechanism 2 of 4): the continuous-
//! mesh model. Directional information comes from the recovered
//! HESSIAN of the primal (second differences on the nodal lattice),
//! importance from the adjoint weight per cell; the metric is the
//! absolute Hessian rescaled so the implied cell count Σ√det(M)·|K|
//! meets a target complexity — fs-mesh's remesher consumes it as a
//! `MetricField` (the battery drives that path end-to-end on a planar
//! sheet).
//!
//! v1 surface: UNIFORM grids (second differences need a regular
//! stencil); graded-tree recovery is a recorded follow-up.

use fs_cutfem::Quadtree;
use std::collections::BTreeMap;

/// Per-cell 2×2 metric tensors, complexity-normalized.
///
/// `weight` is the per-cell adjoint importance (|z| or |η| mass);
/// `target_cells` is the complexity budget the metric should imply.
///
/// # Panics
/// If the grid is not uniform (all leaves at one level).
#[must_use]
pub fn synthesize_metric(
    grid: &Quadtree,
    nodal: &BTreeMap<(u32, u32), f64>,
    weight: &BTreeMap<(u32, u32, u32), f64>,
    target_cells: f64,
) -> BTreeMap<(u32, u32, u32), [[f64; 2]; 2]> {
    let level = grid.leaves().next().expect("nonempty grid").0;
    assert!(
        grid.leaves().all(|c| c.0 == level),
        "metric synthesis v1 needs a uniform grid"
    );
    let h = 1.0 / f64::from(1u32 << level);
    let s = 1u32 << (grid.max_level() - level);
    let ext = grid.node_extent();
    let val = |gi: i64, gj: i64| -> f64 {
        let gi = gi.clamp(0, i64::from(ext)) as u32;
        let gj = gj.clamp(0, i64::from(ext)) as u32;
        nodal.get(&(gi, gj)).copied().unwrap_or(0.0)
    };
    let mut raw: BTreeMap<(u32, u32, u32), [[f64; 2]; 2]> = BTreeMap::new();
    let mut mass;
    let mut _pre = 0.0f64;
    for c in grid.leaves() {
        let corners = grid.corner_nodes(c);
        // Cell-centered second differences from the corner stencil.
        let (gi, gj) = (i64::from(corners[0].0), i64::from(corners[0].1));
        let st = i64::from(s);
        // Forward AND backward second differences, larger magnitude
        // wins: a one-sided stencil straddling an odd layer's
        // inflection point vanishes by antisymmetry and would zero the
        // metric exactly where the layer is (measured pathology).
        let maxabs = |a: f64, b: f64| if a.abs() >= b.abs() { a } else { b };
        let hxx_f = (val(gi + 2 * st, gj) - 2.0 * val(gi + st, gj) + val(gi, gj)
            + val(gi + 2 * st, gj + st)
            - 2.0 * val(gi + st, gj + st)
            + val(gi, gj + st))
            / (2.0 * h * h);
        let hxx_b = (val(gi + st, gj) - 2.0 * val(gi, gj) + val(gi - st, gj)
            + val(gi + st, gj + st)
            - 2.0 * val(gi, gj + st)
            + val(gi - st, gj + st))
            / (2.0 * h * h);
        let hxx = maxabs(hxx_f, hxx_b);
        let hyy_f = (val(gi, gj + 2 * st) - 2.0 * val(gi, gj + st) + val(gi, gj)
            + val(gi + st, gj + 2 * st)
            - 2.0 * val(gi + st, gj + st)
            + val(gi + st, gj))
            / (2.0 * h * h);
        let hyy_b = (val(gi, gj + st) - 2.0 * val(gi, gj) + val(gi, gj - st)
            + val(gi + st, gj + st)
            - 2.0 * val(gi + st, gj)
            + val(gi + st, gj - st))
            / (2.0 * h * h);
        let hyy = maxabs(hyy_f, hyy_b);
        let hxy = (val(gi + st, gj + st) - val(gi + st, gj) - val(gi, gj + st) + val(gi, gj))
            / (h * h);
        // Absolute Hessian: |H| via closed-form 2×2 spectral abs.
        let tr = hxx + hyy;
        let det = hxx * hyy - hxy * hxy;
        let disc = (0.25 * tr * tr - det).max(0.0).sqrt();
        let (l1, l2) = (0.5 * tr + disc, 0.5 * tr - disc);
        // Eigenvector of l1.
        let (ex, ey) = if hxy.abs() > 1e-30 {
            let n = (l1 - hyy).hypot(hxy);
            ((l1 - hyy) / n, hxy / n)
        } else if hxx >= hyy {
            (1.0, 0.0)
        } else {
            (0.0, 1.0)
        };
        // The metric's strong axis is the LARGEST-|λ| eigenvector —
        // l1 is the largest SIGNED one, so a dominant negative
        // curvature (the concave side of a layer) lives on l2 and the
        // axes must swap.
        let (ex, ey, a1, a2) = if l2.abs() > l1.abs() {
            (-ey, ex, l2.abs().max(1e-12), l1.abs().max(1e-12))
        } else {
            (ex, ey, l1.abs().max(1e-12), l2.abs().max(1e-12))
        };
        // Anisotropy cap (100:1) BEFORE normalization — flooring after
        // scaling would inflate the implied complexity by orders.
        let a2 = a2.max(a1 * 1e-4);
        let w = weight.get(&c).copied().unwrap_or(0.0).abs().max(1e-12);
        // M = w·(a1 e⊗e + a2 e⊥⊗e⊥).
        let m = [
            [
                w * (a1 * ex * ex + a2 * ey * ey),
                w * (a1 - a2) * ex * ey,
            ],
            [
                w * (a1 - a2) * ex * ey,
                w * (a1 * ey * ey + a2 * ex * ex),
            ],
        ];
        raw.insert(c, m);
    }
    // Complexity normalization: cells implied = Σ√det(sM)·|K| = s·mass
    // (2D scaling) → s = target / mass.
    // Flat-region floor relative to the global scale, BEFORE the
    // complexity normalization (mass recomputed after flooring).
    let amax = raw
        .values()
        .map(|m| m[0][0].max(m[1][1]))
        .fold(0.0f64, f64::max)
        .max(1e-300);
    mass = 0.0;
    for m in raw.values_mut() {
        let floor = amax * 1e-4;
        m[0][0] = m[0][0].max(floor);
        m[1][1] = m[1][1].max(floor);
        let dm = (m[0][0] * m[1][1] - m[0][1] * m[1][0]).max(0.0);
        mass += dm.sqrt() * h * h;
    }
    let scale = if mass > 0.0 { target_cells / mass } else { 1.0 };
    for m in raw.values_mut() {
        for row in m.iter_mut() {
            for v in row.iter_mut() {
                *v *= scale;
            }
        }
    }
    raw
}
