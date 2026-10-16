//! The DWR core: `J(u) − J(u_h) ≈ r(z_{h/2} − z_h)` with coarse and
//! enriched adjoints solved on the two active spaces. Indicators are the
//! SIGNED per-coarse-cell contributions of the full discrete residual
//! — interior `∫ f·w − ∇u_h·∇w` (Q1 is Laplacian-free, so the interior
//! strong residual is exactly f) plus the Nitsche interface terms of
//! fs-cutfem's form. The coarse ghost-penalty contribution is an
//! O(γh)-scaled correction deliberately absorbed into effectivity
//! (measured by the battery's documented band).

use fs_cutfem::quad::tensor_gauss;
use fs_cutfem::{CellClass, CellKey, CutFemError, CutSdf, FemParams, NodeKey, Quadtree, Space};
use std::collections::{BTreeMap, BTreeSet};

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
    space: &Space<'_>,
    cell: CellKey,
    nodal: &BTreeMap<NodeKey, f64>,
    p: [f64; 2],
) -> Result<(f64, [f64; 2]), CutFemError> {
    if !space.active_cells().contains(&cell) {
        return Err(CutFemError::InvalidFemInput {
            what: format!("cannot evaluate inactive scalar CutFEM cell {cell:?}"),
        });
    }
    if p.iter().any(|coordinate| !coordinate.is_finite()) {
        return Err(CutFemError::InvalidFemInput {
            what: format!("scalar CutFEM evaluation point {p:?} must be finite"),
        });
    }
    let grid = space.grid();
    let (lo, hi) = grid.rect(cell);
    if (0..2).any(|axis| p[axis] < lo[axis] || p[axis] > hi[axis]) {
        return Err(CutFemError::InvalidFemInput {
            what: format!(
                "scalar CutFEM evaluation point {p:?} lies outside cell {cell:?} rectangle {lo:?}--{hi:?}"
            ),
        });
    }
    let corners = grid.corner_nodes(cell);
    let (n, g) = q1(lo, hi, p);
    if n.iter().any(|shape| !shape.is_finite())
        || g.iter().flatten().any(|gradient| !gradient.is_finite())
    {
        return Err(CutFemError::InvalidFemInput {
            what: format!("scalar Q1 evaluation is non-finite on cell {cell:?} at {p:?}"),
        });
    }
    let mut v = 0.0;
    let mut gr = [0.0f64; 2];
    for a in 0..4 {
        let val = nodal
            .get(&corners[a])
            .copied()
            .ok_or_else(|| CutFemError::InvalidFemInput {
                what: format!(
                    "scalar CutFEM cell {cell:?} is missing nodal corner {:?}",
                    corners[a]
                ),
            })?;
        if !val.is_finite() {
            return Err(CutFemError::InvalidFemInput {
                what: format!(
                    "scalar CutFEM nodal value at corner {:?} is non-finite",
                    corners[a]
                ),
            });
        }
        v += n[a] * val;
        gr[0] += g[a][0] * val;
        gr[1] += g[a][1] * val;
    }
    if !v.is_finite() || gr.iter().any(|gradient| !gradient.is_finite()) {
        return Err(CutFemError::InvalidFemInput {
            what: format!("scalar CutFEM field evaluation is non-finite on cell {cell:?}"),
        });
    }
    Ok((v, gr))
}

/// The bulk quadrature rule for one cell of a built space.
fn bulk_rule(space: &Space<'_>, cell: CellKey) -> Result<Vec<([f64; 2], f64)>, CutFemError> {
    let grid = space.grid();
    let (lo, hi) = grid.rect(cell);
    match space.class_of(cell) {
        Some(CellClass::Inside) => {
            let mut v = Vec::with_capacity(9);
            tensor_gauss(lo, hi, &mut v);
            Ok(v)
        }
        Some(CellClass::Cut) => space
            .cut_rules()
            .get(&cell)
            .map(|rules| rules.bulk.clone())
            .ok_or_else(|| CutFemError::InvalidFemInput {
                what: format!("active cut cell {cell:?} has no retained quadrature rule"),
            }),
        Some(CellClass::Outside) | None => Err(CutFemError::InvalidFemInput {
            what: format!("cell {cell:?} is not active in the scalar CutFEM space"),
        }),
    }
}

/// `J(u_h) = ∫ jw·u_h` over the active domain.
///
/// # Errors
/// Refuses missing/non-finite nodal data, non-finite goal weights, or
/// topology/rule mismatches instead of inventing a zero contribution.
pub fn goal_value(
    space: &Space<'_>,
    nodal: &BTreeMap<NodeKey, f64>,
    goal: &GoalContext<'_>,
) -> Result<f64, CutFemError> {
    let mut j = 0.0;
    for &cell in space.active_cells() {
        for (p, w) in bulk_rule(space, cell)? {
            let (u, _) = eval_cell(space, cell, nodal, p)?;
            let weight = (goal.weight)(p[0], p[1]);
            if !weight.is_finite() {
                return Err(CutFemError::InvalidFemInput {
                    what: format!("goal weight is non-finite at {p:?}"),
                });
            }
            j += w * weight * u;
            if !j.is_finite() {
                return Err(CutFemError::InvalidFemInput {
                    what: "scalar goal accumulation is non-finite".to_string(),
                });
            }
        }
    }
    Ok(j)
}

fn validate_active_hierarchy(
    coarse: &BTreeSet<CellKey>,
    fine: &BTreeSet<CellKey>,
) -> Result<(), CutFemError> {
    for &(level, i, j) in fine {
        let cell = (level, i, j);
        if level == 0 {
            return Err(CutFemError::InvalidFemInput {
                what: format!("enriched active cell {cell:?} has no parent level"),
            });
        }
        let parent = (level - 1, i / 2, j / 2);
        if !coarse.contains(&parent) {
            return Err(CutFemError::InvalidFemInput {
                what: format!(
                    "enriched active cell {cell:?} has inactive coarse parent {parent:?}"
                ),
            });
        }
    }
    for &(level, i, j) in coarse {
        let cell = (level, i, j);
        let has_active_child = (0..2u32)
            .any(|di| (0..2u32).any(|dj| fine.contains(&(level + 1, 2 * i + di, 2 * j + dj))));
        if !has_active_child {
            return Err(CutFemError::InvalidFemInput {
                what: format!("coarse active cell {cell:?} has no enriched active child"),
            });
        }
    }
    Ok(())
}

/// Run the DWR estimate on one grid: primal and adjoint solves on the
/// coarse space, an enriched adjoint on the once-refined space, and signed
/// per-cell indicators.
///
/// # Errors
/// Propagates fs-cutfem build/solve teaching errors and refuses non-nested
/// active coverage or missing/non-finite field evidence.
#[allow(clippy::too_many_lines)] // primal + adjoint + weighting is one narrative
pub fn estimate(
    grid: &Quadtree,
    sdf: &dyn CutSdf,
    params: FemParams,
    f: &dyn Fn(f64, f64) -> f64,
    g: &dyn Fn(f64, f64) -> f64,
    goal: &GoalContext<'_>,
) -> Result<DwrEstimate, CutFemError> {
    if grid.max_level() >= 16 {
        return Err(CutFemError::InvalidFemInput {
            what: "scalar DWR enrichment requires one lattice level below the level-16 cap"
                .to_string(),
        });
    }
    let space = Space::build(grid, sdf, params)?;
    let sol = space.solve(f, g)?;
    let j_primal = goal_value(&space, &sol.nodal, goal)?;
    // Enriched adjoint: one-level-finer solve, homogeneous data.
    let fine = grid.refined_once();
    let fspace = Space::build(&fine, sdf, params)?;
    validate_active_hierarchy(space.active_cells(), fspace.active_cells())?;
    let coarse_adj = space.solve(goal.weight, &|_, _| 0.0)?;
    let adj = fspace.solve(goal.weight, &|_, _| 0.0)?;
    // Indicators: loop coarse leaves; integrate on their fine children
    // with the standard two-level weight w = z_fine − z_coarse.
    let mut indicators: BTreeMap<(u32, u32, u32), f64> = BTreeMap::new();
    for &cell in space.active_cells() {
        let h = grid.cell_h(cell);
        let pen = params.nitsche_beta / h;
        let mut eta = 0.0f64;
        let (lv, i, j) = cell;
        for di in 0..2u32 {
            for dj in 0..2u32 {
                let child = (lv + 1, 2 * i + di, 2 * j + dj);
                // Bulk: ∫ f·w − ∇u_h·∇w on the child's rule.
                if fspace.active_cells().contains(&child) {
                    for (p, w) in bulk_rule(&fspace, child)? {
                        let (zf, gzf) = eval_cell(&fspace, child, &adj.nodal, p)?;
                        let (zc, gzc) = eval_cell(&space, cell, &coarse_adj.nodal, p)?;
                        let wgt = zf - zc;
                        let gw = [gzf[0] - gzc[0], gzf[1] - gzc[1]];
                        let (_, gu) = eval_cell(&space, cell, &sol.nodal, p)?;
                        let source = f(p[0], p[1]);
                        if !source.is_finite() {
                            return Err(CutFemError::InvalidFemInput {
                                what: format!("scalar DWR source is non-finite at {p:?}"),
                            });
                        }
                        eta += w * (source * wgt - (gu[0] * gw[0] + gu[1] * gw[1]));
                        if !eta.is_finite() {
                            return Err(CutFemError::InvalidFemInput {
                                what: format!(
                                    "scalar DWR bulk indicator is non-finite on coarse cell {cell:?}"
                                ),
                            });
                        }
                    }
                }
                // Nitsche interface terms of the coarse form:
                // r_Γ(w) = ∫_Γ ∂n u_h·w + (∂n w + pen·w)(u_h − g)…
                // with the sign convention of fs-cutfem's assembly.
                if fspace.class_of(child) == Some(CellClass::Cut) {
                    let rules = fspace.cut_rules().get(&child).ok_or_else(|| {
                        CutFemError::InvalidFemInput {
                            what: format!(
                                "enriched cut cell {child:?} has no retained interface rule"
                            ),
                        }
                    })?;
                    for &(p, w, nrm) in &rules.iface {
                        let (zf, gzf) = eval_cell(&fspace, child, &adj.nodal, p)?;
                        let (zc, gzc) = eval_cell(&space, cell, &coarse_adj.nodal, p)?;
                        let wgt = zf - zc;
                        let dnw = (gzf[0] - gzc[0]) * nrm[0] + (gzf[1] - gzc[1]) * nrm[1];
                        let (u, gu) = eval_cell(&space, cell, &sol.nodal, p)?;
                        let dnu = gu[0] * nrm[0] + gu[1] * nrm[1];
                        let gv = g(p[0], p[1]);
                        if !gv.is_finite() {
                            return Err(CutFemError::InvalidFemInput {
                                what: format!("scalar DWR boundary data is non-finite at {p:?}"),
                            });
                        }
                        eta += w * (dnu * wgt + (dnw - pen * wgt) * (u - gv));
                        if !eta.is_finite() {
                            return Err(CutFemError::InvalidFemInput {
                                what: format!(
                                    "scalar DWR interface indicator is non-finite on coarse cell {cell:?}"
                                ),
                            });
                        }
                    }
                }
            }
        }
        indicators.insert(cell, eta);
    }
    let eta_signed: f64 = indicators.values().sum();
    let eta_abs: f64 = indicators.values().map(|v| v.abs()).sum();
    if !eta_signed.is_finite() || !eta_abs.is_finite() {
        return Err(CutFemError::InvalidFemInput {
            what: "scalar DWR indicator totals are non-finite".to_string(),
        });
    }
    Ok(DwrEstimate {
        eta_signed,
        eta_abs,
        indicators,
        j_primal,
        dofs: space.dof_count(),
        nodal: sol.nodal,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_cutfem::{Circle, HalfPlane};

    fn all_inside() -> HalfPlane {
        HalfPlane {
            normal: [1.0, 0.0],
            offset: 2.0,
        }
    }

    fn unit_nodal(space: &Space<'_>) -> BTreeMap<NodeKey, f64> {
        let mut nodal = BTreeMap::new();
        for &cell in space.active_cells() {
            for node in space.grid().corner_nodes(cell) {
                nodal.insert(node, 1.0);
            }
        }
        nodal
    }

    fn invalid_what(error: CutFemError) -> String {
        match error {
            CutFemError::InvalidFemInput { what } => what,
            other => panic!("expected InvalidFemInput, got {other:?}"),
        }
    }

    #[test]
    fn goal_value_refuses_a_missing_active_corner() {
        let grid = Quadtree::uniform(1);
        let domain = all_inside();
        let space = Space::build(&grid, &domain, FemParams::default()).expect("space");
        let mut nodal = unit_nodal(&space);
        let cell = *space.active_cells().first().expect("active cell");
        let missing = grid.corner_nodes(cell)[0];
        nodal.remove(&missing);
        let goal = GoalContext {
            weight: &|_, _| 1.0,
        };

        let error = goal_value(&space, &nodal, &goal).expect_err("missing corner must refuse");
        assert!(invalid_what(error).contains("missing nodal corner"));
    }

    #[test]
    fn goal_value_refuses_non_finite_nodal_and_weight_evidence() {
        let grid = Quadtree::uniform(1);
        let domain = all_inside();
        let space = Space::build(&grid, &domain, FemParams::default()).expect("space");
        let mut nodal = unit_nodal(&space);
        let node = *nodal.first_key_value().expect("nodal value").0;
        nodal.insert(node, f64::NAN);
        let finite_goal = GoalContext {
            weight: &|_, _| 1.0,
        };
        let error = goal_value(&space, &nodal, &finite_goal)
            .expect_err("non-finite nodal value must refuse");
        assert!(invalid_what(error).contains("nodal value"));

        let nodal = unit_nodal(&space);
        let non_finite_goal = GoalContext {
            weight: &|_, _| f64::INFINITY,
        };
        let error = goal_value(&space, &nodal, &non_finite_goal)
            .expect_err("non-finite goal weight must refuse");
        assert!(invalid_what(error).contains("goal weight"));
    }

    #[test]
    fn enriched_evaluation_never_zero_fills_a_missing_corner() {
        let grid = Quadtree::uniform(1);
        let domain = all_inside();
        let space = Space::build(&grid, &domain, FemParams::default()).expect("space");
        let cell = (1, 0, 0);
        let nodal = grid
            .corner_nodes(cell)
            .into_iter()
            .take(3)
            .map(|node| (node, 1.0))
            .collect();
        let (lo, hi) = grid.rect(cell);
        let point = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];

        let error = eval_cell(&space, cell, &nodal, point)
            .expect_err("missing enriched field corner must refuse");
        assert!(invalid_what(error).contains("missing nodal corner"));
    }

    #[test]
    fn active_hierarchy_is_checked_in_both_directions() {
        let coarse = BTreeSet::from([(1, 0, 0)]);
        let one_child = BTreeSet::from([(2, 0, 0)]);
        validate_active_hierarchy(&coarse, &one_child)
            .expect("one active child is sufficient coverage");

        let orphan_child = BTreeSet::from([(2, 2, 2)]);
        let error = validate_active_hierarchy(&coarse, &orphan_child)
            .expect_err("orphan enriched child must refuse");
        assert!(invalid_what(error).contains("inactive coarse parent"));

        let coarse_with_orphan = BTreeSet::from([(1, 0, 0), (1, 1, 1)]);
        let error = validate_active_hierarchy(&coarse_with_orphan, &one_child)
            .expect_err("orphan coarse cell must refuse");
        assert!(invalid_what(error).contains("no enriched active child"));
    }

    #[test]
    fn valid_disk_cut_with_absent_mapped_fine_nodes_still_estimates() {
        let grid = Quadtree::uniform(4);
        let domain = Circle {
            center: [0.5, 0.5],
            radius: 0.35,
        };
        let params = FemParams::default();
        let coarse_space = Space::build(&grid, &domain, params).expect("coarse space");
        let fine_grid = grid.refined_once();
        let fine_space = Space::build(&fine_grid, &domain, params).expect("fine space");
        let cell = (4, 2, 5);
        assert!(coarse_space.active_cells().contains(&cell));
        let active_children = (0..2u32)
            .flat_map(|di| (0..2u32).map(move |dj| (5, 4 + di, 10 + dj)))
            .filter(|child| fine_space.active_cells().contains(child))
            .count();
        assert!((1..4).contains(&active_children));

        let one = |_: f64, _: f64| 1.0;
        let zero = |_: f64, _: f64| 0.0;
        let fine_nodes: BTreeSet<NodeKey> = fine_space
            .active_cells()
            .iter()
            .flat_map(|&fine_cell| fine_grid.corner_nodes(fine_cell))
            .collect();
        let absent_mapped_nodes = grid
            .corner_nodes(cell)
            .map(|node| (2 * node.0, 2 * node.1))
            .into_iter()
            .filter(|node| !fine_nodes.contains(node))
            .count();
        assert!(absent_mapped_nodes > 0);

        let goal = GoalContext { weight: &one };
        let estimate = estimate(&grid, &domain, params, &zero, &zero, &goal)
            .expect("coarse adjoint makes the valid non-nested cut estimable");
        assert!(estimate.eta_signed.is_finite());
        assert!(estimate.eta_abs.is_finite());
    }

    #[test]
    fn level_sixteen_estimate_refuses_before_refinement() {
        let grid = Quadtree::with_room(0, 16);
        let domain = all_inside();
        let goal = GoalContext {
            weight: &|_, _| 1.0,
        };
        let zero = |_: f64, _: f64| 0.0;

        let error = estimate(&grid, &domain, FemParams::default(), &zero, &zero, &goal)
            .expect_err("lattice-cap enrichment must refuse");
        assert!(invalid_what(error).contains("level-16 cap"));
    }

    #[test]
    fn level_fifteen_retains_one_enrichment_level() {
        let grid = Quadtree::with_room(0, 15);
        let domain = Circle {
            center: [0.5, 0.5],
            radius: 0.35,
        };
        let goal = GoalContext {
            weight: &|_, _| 1.0,
        };
        let zero = |_: f64, _: f64| 0.0;

        let estimate = estimate(&grid, &domain, FemParams::default(), &zero, &zero, &goal)
            .expect("level 15 has exactly one safe enriched lattice level");
        assert!(estimate.eta_signed.is_finite());
        assert!(estimate.eta_abs.is_finite());
    }
}
