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
    let mut mass = 0.0f64;
    for c in grid.leaves() {
        let corners = grid.corner_nodes(c);
        // Cell-centered second differences from the corner stencil.
        let (gi, gj) = (i64::from(corners[0].0), i64::from(corners[0].1));
        let st = i64::from(s);
        let hxx = (val(gi + 2 * st, gj) - 2.0 * val(gi + st, gj) + val(gi, gj)
            + val(gi + 2 * st, gj + st)
            - 2.0 * val(gi + st, gj + st)
            + val(gi, gj + st))
            / (2.0 * h * h);
        let hyy = (val(gi, gj + 2 * st) - 2.0 * val(gi, gj + st) + val(gi, gj)
            + val(gi + st, gj + 2 * st)
            - 2.0 * val(gi + st, gj + st)
            + val(gi + st, gj))
            / (2.0 * h * h);
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
        let (a1, a2) = (l1.abs().max(1e-12), l2.abs().max(1e-12));
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
        let dm = (m[0][0] * m[1][1] - m[0][1] * m[1][0]).max(0.0);
        mass += dm.sqrt() * h * h;
        raw.insert(c, m);
    }
    // Complexity normalization: cells implied = Σ√det(sM)·|K| = s·mass
    // (2D scaling) → s = target / mass.
    let scale = if mass > 0.0 { target_cells / mass } else { 1.0 };
    for m in raw.values_mut() {
        for row in m.iter_mut() {
            for v in row.iter_mut() {
                *v *= scale;
            }
        }
        // Anisotropy cap 100:1 and a floor keep the remesher sane.
        let tr = m[0][0] + m[1][1];
        let floor = tr.max(1e-12) * 1e-2 * 0.5;
        m[0][0] = m[0][0].max(floor);
        m[1][1] = m[1][1].max(floor);
    }
    raw
}
