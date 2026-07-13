//! The octree h-refinement loop (mechanism 1 of 4): solve → estimate →
//! Dörfler-mark → split → rebalance → restore the uniform cut band
//! fs-cutfem's ghost penalty requires. The accuracy-per-DOF trajectory
//! is the ledgered evidence — goal-oriented refinement must beat
//! uniform on localized QoIs or the estimator is decoration.

use crate::estimate::{DwrEstimate, GoalContext, estimate};
use crate::mark::dorfler;
use fs_cutfem::{CutFemError, CutSdf, FemParams, Quadtree};
use std::fmt::Write as _;

/// One adaptive iteration's evidence.
#[derive(Debug, Clone)]
pub struct AdaptStep {
    /// Primal free DOFs at this step.
    pub dofs: usize,
    /// J(u_h).
    pub j: f64,
    /// Signed estimate.
    pub eta_signed: f64,
    /// Marking mass Σ|η_K|.
    pub eta_abs: f64,
    /// Cells marked (0 on the final, estimate-only step).
    pub marked: usize,
}

impl AdaptStep {
    /// Ledger-style JSON row.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"dofs\":{},\"j\":{:.10e},\"eta_signed\":{:.4e},\
             \"eta_abs\":{:.4e},\"marked\":{}}}",
            self.dofs, self.j, self.eta_signed, self.eta_abs, self.marked
        );
        s
    }
}

fn refinement_ceiling(max_level: u32) -> Result<u32, CutFemError> {
    max_level
        .checked_sub(1)
        .ok_or_else(|| CutFemError::InvalidFemInput {
            what: "scalar DWR adaptivity requires refinement headroom when another marked iteration remains"
                .to_string(),
        })
}

/// Run `iters` adaptive cycles (the last records without refining).
/// The grid must carry enough `with_room` headroom for the splits.
///
/// # Errors
/// Propagates fs-cutfem build/solve errors and the DWR estimator's structured
/// refusals for non-nested active coverage or missing/non-finite field
/// evidence. A non-final iteration with marked cells also refuses when the
/// grid has no reserved refinement headroom.
#[allow(clippy::too_many_arguments)] // the PDE problem statement is the argument list
pub fn adapt_loop(
    grid: &mut Quadtree,
    sdf: &dyn CutSdf,
    params: FemParams,
    f: &dyn Fn(f64, f64) -> f64,
    g: &dyn Fn(f64, f64) -> f64,
    goal: &GoalContext<'_>,
    theta: f64,
    iters: usize,
) -> Result<(Vec<AdaptStep>, DwrEstimate), CutFemError> {
    let mut steps = Vec::new();
    loop {
        let est = estimate(grid, sdf, params, f, g, goal)?;
        let last = steps.len() + 1 >= iters;
        if last {
            steps.push(AdaptStep {
                dofs: est.dofs,
                j: est.j_primal,
                eta_signed: est.eta_signed,
                eta_abs: est.eta_abs,
                marked: 0,
            });
            return Ok((steps, est));
        }
        let marked = dorfler(&est.indicators, theta);
        steps.push(AdaptStep {
            dofs: est.dofs,
            j: est.j_primal,
            eta_signed: est.eta_signed,
            eta_abs: est.eta_abs,
            marked: marked.len(),
        });
        if !marked.is_empty() {
            let ceiling = refinement_ceiling(grid.max_level())?;
            for c in &marked {
                if grid.is_leaf(*c) && c.0 < ceiling {
                    grid.split(*c);
                }
            }
        }
        grid.balance();
        // Restore the uniform interface band at the finest level any
        // cut-adjacent cell reached (the ghost-penalty precondition).
        let mut band_level = 0u32;
        for c in grid.leaves().collect::<Vec<_>>() {
            let (lo, hi) = grid.rect(c);
            let h = hi[0] - lo[0];
            let ilo = [(lo[0] - h).max(0.0), (lo[1] - h).max(0.0)];
            let ihi = [(hi[0] + h).min(1.0), (hi[1] + h).min(1.0)];
            if sdf.enclose(ilo, ihi).contains_zero() {
                band_level = band_level.max(c.0);
            }
        }
        grid.refine_toward_interface(sdf, band_level);
    }
}

#[cfg(test)]
mod tests {
    use super::{adapt_loop, refinement_ceiling};
    use crate::estimate::{GoalContext, estimate};
    use crate::mark::dorfler;
    use fs_cutfem::{Circle, CutFemError, FemParams, Quadtree};

    #[test]
    fn adapt_loop_refuses_zero_level_marked_headroom_without_bypassing_safe_paths() {
        let error = refinement_ceiling(0).expect_err("zero-level grid has no refinement headroom");
        assert!(matches!(
            error,
            CutFemError::InvalidFemInput { what }
                if what.contains("refinement headroom")
        ));
        assert_eq!(
            refinement_ceiling(1).expect("one level has a zero ceiling"),
            0
        );
        assert_eq!(refinement_ceiling(4).expect("positive headroom"), 3);

        let sdf = Circle {
            center: [0.5, 0.5],
            radius: 0.35,
        };
        let params = FemParams {
            ghost_gamma: 0.0,
            ..FemParams::default()
        };
        let source = |x: f64, y: f64| 1.0 + x + 2.0 * y;
        let boundary = |_: f64, _: f64| 0.0;
        let weight = |x: f64, y: f64| 1.0 + x * x + y;
        let goal = GoalContext { weight: &weight };
        let probe_grid = Quadtree::uniform(0);
        let probe = estimate(&probe_grid, &sdf, params, &source, &boundary, &goal)
            .expect("level-zero fixture has a valid DWR estimate");
        assert!(
            !dorfler(&probe.indicators, 0.5).is_empty(),
            "nonzero fixture must reach the marked headroom gate"
        );

        let mut final_only = Quadtree::uniform(0);
        let (steps, _) = adapt_loop(
            &mut final_only,
            &sdf,
            params,
            &source,
            &boundary,
            &goal,
            0.5,
            1,
        )
        .expect("the final estimate-only iteration requires no headroom");
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].marked, 0);

        let mut no_headroom = Quadtree::uniform(0);
        let error = adapt_loop(
            &mut no_headroom,
            &sdf,
            params,
            &source,
            &boundary,
            &goal,
            0.5,
            2,
        )
        .expect_err("a marked continuation must refuse zero-level headroom");
        assert!(matches!(
            error,
            CutFemError::InvalidFemInput { what }
                if what.contains("refinement headroom")
        ));

        let zero = |_: f64, _: f64| 0.0;
        let empty_goal = GoalContext { weight: &zero };
        let mut empty_marking = Quadtree::uniform(0);
        let (steps, _) = adapt_loop(
            &mut empty_marking,
            &sdf,
            params,
            &zero,
            &zero,
            &empty_goal,
            0.5,
            2,
        )
        .expect("an empty marked set requires no refinement headroom");
        assert_eq!(steps.len(), 2);
        assert!(steps.iter().all(|step| step.marked == 0));
    }
}
