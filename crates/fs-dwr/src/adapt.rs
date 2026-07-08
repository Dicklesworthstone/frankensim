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

/// Run `iters` adaptive cycles (the last records without refining).
/// The grid must carry enough `with_room` headroom for the splits.
///
/// # Errors
/// fs-cutfem build/solve teaching errors.
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
        for c in &marked {
            if grid.is_leaf(*c) && c.0 < grid.max_level() - 1 {
                grid.split(*c);
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
