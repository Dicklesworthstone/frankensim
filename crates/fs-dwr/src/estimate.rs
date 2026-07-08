//! The DWR core: `J(u) − J(u_h) ≈ r(z − I_h z)` with the enriched
//! adjoint z solved on the once-refined grid. Indicators are the
//! SIGNED per-coarse-cell contributions of the full discrete residual
//! — interior `∫ f·w − ∇u_h·∇w` (Q1 is Laplacian-free, so the interior
//! strong residual is exactly f) plus the Nitsche interface terms of
//! fs-cutfem's form. The coarse ghost-penalty contribution is an
//! O(γh)-scaled correction deliberately absorbed into effectivity
//! (measured by the battery's documented band).

use fs_cutfem::quad::{cut_cell_rules, tensor_gauss};
use fs_cutfem::{CellClass, CutFemError, CutSdf, FemParams, Quadtree, Space};
use std::collections::BTreeMap;

/// A volumetric goal functional `J(u) = ∫ jw·u` (region averages,
/// windowed integrals — the localized-QoI family).
pub struct GoalContext<'a> {
    /// The goal weight field jw.
    pub weight: &'a dyn Fn(f64, f64) -> f64,
}

/// The DWR output for one grid.
#[derive(Debug, Clone)]
pub struct DwrEstimate {
    /// Signed estimate Σ η_K ≈ J(u) − J(u_h).
    pub eta_signed: f64,
    /// Σ |η_K| (the marking mass).
    pub eta_abs: f64,
    /// Signed indicator per coarse leaf.
    pub indicators: BTreeMap<(u32, u32, u32), f64>,
    /// J(u_h).
    pub j_primal: f64,
    /// Primal free-DOF count.
    pub dofs: usize,
    /// Primal nodal solution.
    pub nodal: BTreeMap<(u32, u32), f64>,
}

/// Q1 shapes on an axis-aligned cell (fs-cutfem corner order).
pub(crate) fn q1(lo: [f64; 2], hi: [f64; 2], p: [f64; 2]) -> ([f64; 4], [[f64; 2]; 4]) {
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = (p[0] - lo[0]) / hx;
    let et = (p[1] - lo[1]) / hy;
    (
        [
            (1.0 - xi) * (1.0 - et),
            xi * (1.0 - et),
            xi * et,
            (1.0 - xi) * et,
        ],
        [
            [-(1.0 - et) / hx, -(1.0 - xi) / hy],
            [(1.0 - et) / hx, -xi / hy],
            [et / hx, xi / hy],
            [-et / hx, (1.0 - xi) / hy],
        ],
    )
}

/// Evaluate a nodal field and its gradient on one cell at a point.
fn eval_cell(
    grid: &Quadtree,
    cell: (u32, u32, u32),
    nodal: &BTreeMap<(u32, u32), f64>,
    p: [f64; 2],
) -> (f64, [f64; 2]) {
    let (lo, hi) = grid.rect(cell);
    let corners = grid.corner_nodes(cell);
    let (n, g) = q1(lo, hi, p);
    let mut v = 0.0;
    let mut gr = [0.0f64; 2];
    for a in 0..4 {
        let val = nodal.get(&corners[a]).copied().unwrap_or(0.0);
        v += n[a] * val;
        gr[0] += g[a][0] * val;
        gr[1] += g[a][1] * val;
    }
    (v, gr)
}

/// The bulk quadrature rule for one cell of a built space.
fn bulk_rule(
    space: &Space<'_>,
    grid: &Quadtree,
    sdf: &dyn CutSdf,
    cell: (u32, u32, u32),
    depth: u32,
) -> Vec<([f64; 2], f64)> {
    let (lo, hi) = grid.rect(cell);
    match space.class_of(cell) {
        Some(CellClass::Inside) => {
            let mut v = Vec::with_capacity(9);
            tensor_gauss(lo, hi, &mut v);
            v
        }
        Some(CellClass::Cut) => cut_cell_rules(sdf, lo, hi, depth).bulk,
        _ => Vec::new(),
    }
}

/// `J(u_h) = ∫ jw·u_h` over the active domain.
#[must_use]
pub fn goal_value(
    space: &Space<'_>,
    grid: &Quadtree,
    sdf: &dyn CutSdf,
    nodal: &BTreeMap<(u32, u32), f64>,
    goal: &GoalContext<'_>,
    depth: u32,
) -> f64 {
    let mut j = 0.0;
    for cell in grid.leaves() {
        for (p, w) in bulk_rule(space, grid, sdf, cell, depth) {
            let (u, _) = eval_cell(grid, cell, nodal, p);
            j += w * (goal.weight)(p[0], p[1]) * u;
        }
    }
    j
}

/// Run the DWR estimate on one grid: primal solve, enriched adjoint on
/// the once-refined grid, signed per-cell indicators.
///
/// # Errors
/// fs-cutfem build/solve teaching errors.
#[allow(clippy::too_many_lines)] // primal + adjoint + weighting is one narrative
pub fn estimate(
    grid: &Quadtree,
    sdf: &dyn CutSdf,
    params: FemParams,
    f: &dyn Fn(f64, f64) -> f64,
    g: &dyn Fn(f64, f64) -> f64,
    goal: &GoalContext<'_>,
) -> Result<DwrEstimate, CutFemError> {
    let space = Space::build(grid, sdf, params)?;
    let sol = space.solve(f, g)?;
    let j_primal = goal_value(&space, grid, sdf, &sol.nodal, goal, params.quad_depth);
    // Enriched adjoint: one-level-finer solve, homogeneous data.
    let fine = grid.refined_once();
    let fspace = Space::build(&fine, sdf, params)?;
    let adj = fspace.solve(goal.weight, &|_, _| 0.0)?;
    // Indicators: loop coarse leaves; integrate on their fine children
    // with w = z_fine − I_h z_fine (coarse-node interpolant of z).
    let mut indicators: BTreeMap<(u32, u32, u32), f64> = BTreeMap::new();
    for cell in grid.leaves() {
        if space.class_of(cell) == Some(CellClass::Outside) || space.class_of(cell).is_none() {
            continue;
        }
        let (clo, chi) = grid.rect(cell);
        let ccorners = grid.corner_nodes(cell);
        // Coarse-node values of z (coarse lattice ⊂ fine lattice, ×2).
        let zc: [f64; 4] = core::array::from_fn(|a| {
            let n = ccorners[a];
            adj.nodal.get(&(2 * n.0, 2 * n.1)).copied().unwrap_or(0.0)
        });
        let h = grid.cell_h(cell);
        let pen = params.nitsche_beta / h;
        let mut eta = 0.0f64;
        let (lv, i, j) = cell;
        for di in 0..2u32 {
            for dj in 0..2u32 {
                let child = (lv + 1, 2 * i + di, 2 * j + dj);
                // Bulk: ∫ f·w − ∇u_h·∇w on the child's rule.
                for (p, w) in bulk_rule(&fspace, &fine, sdf, child, params.quad_depth) {
                    let (zf, gzf) = eval_cell(&fine, child, &adj.nodal, p);
                    let (nc, gc) = q1(clo, chi, p);
                    let mut zi = 0.0;
                    let mut gzi = [0.0f64; 2];
                    for a in 0..4 {
                        zi += nc[a] * zc[a];
                        gzi[0] += gc[a][0] * zc[a];
                        gzi[1] += gc[a][1] * zc[a];
                    }
                    let wgt = zf - zi;
                    let gw = [gzf[0] - gzi[0], gzf[1] - gzi[1]];
                    let (_, gu) = eval_cell(grid, cell, &sol.nodal, p);
                    eta += w * (f(p[0], p[1]) * wgt - (gu[0] * gw[0] + gu[1] * gw[1]));
                }
                // Nitsche interface terms of the coarse form:
                // r_Γ(w) = ∫_Γ ∂n u_h·w + (∂n w + pen·w)(u_h − g)…
                // with the sign convention of fs-cutfem's assembly.
                if fspace.class_of(child) == Some(CellClass::Cut) {
                    let (flo, fhi) = fine.rect(child);
                    for (p, w, nrm) in cut_cell_rules(sdf, flo, fhi, params.quad_depth).iface {
                        let (zf, gzf) = eval_cell(&fine, child, &adj.nodal, p);
                        let (nc, gc) = q1(clo, chi, p);
                        let mut zi = 0.0;
                        let mut gzi = [0.0f64; 2];
                        for a in 0..4 {
                            zi += nc[a] * zc[a];
                            gzi[0] += gc[a][0] * zc[a];
                            gzi[1] += gc[a][1] * zc[a];
                        }
                        let wgt = zf - zi;
                        let dnw = (gzf[0] - gzi[0]) * nrm[0] + (gzf[1] - gzi[1]) * nrm[1];
                        let (u, gu) = eval_cell(grid, cell, &sol.nodal, p);
                        let dnu = gu[0] * nrm[0] + gu[1] * nrm[1];
                        let gv = g(p[0], p[1]);
                        eta += w * (dnu * wgt + (dnw - pen * wgt) * (u - gv));
                    }
                }
            }
        }
        indicators.insert(cell, eta);
    }
    let eta_signed: f64 = indicators.values().sum();
    let eta_abs: f64 = indicators.values().map(|v| v.abs()).sum();
    Ok(DwrEstimate {
        eta_signed,
        eta_abs,
        indicators,
        j_primal,
        dofs: space.dof_count(),
        nodal: sol.nodal,
    })
}
