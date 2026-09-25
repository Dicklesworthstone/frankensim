//! Goal-controlled primal work. The existing PCG recurrence proposes fields;
//! only an outward bound on the returned physical field admits success.

use fs_solver::{CgState, CsrOp};

use super::{ConductionError, Cx, LinearGoalAnalysis, LinearGoalAnalyzer, invalid, poll};

/// Explicit accuracy and work policy for a fixed linear thermal goal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearGoalSolveConfig {
    /// Absolute error in the declared goal's units, finite and positive.
    /// A caller may allocate, for example, 0.1 times its discretization
    /// allowance. This driver does not invent or certify that allowance.
    pub absolute_tolerance: f64,
    /// Shared primal iteration budget across all correction attempts.
    /// Zero admits an initial check but no primal work.
    pub max_primal_iterations: usize,
    /// PCG iterations between goal checks, in 1..=32.
    pub check_every: usize,
    /// Extra defect solves after the initial correction attempt. Zero
    /// disables retries, not the goal enclosure or returned-field check.
    pub max_defect_corrections: usize,
}

/// Why the bounded goal-controlled driver stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinearGoalStop {
    /// The returned field's FULL algebraic goal enclosure met the tolerance.
    GoalTolerance,
    /// No more primal iterations were admitted.
    IterationBudget,
    /// The initial field had no full inverse/dual-error goal bound.
    BoundUnavailable,
    /// The checked goal bound did not improve, or a rounded zero defect
    /// left a nonzero outward bound. This is not proof of an accuracy floor.
    NoProgress,
    /// A recurrence ended without meeting the goal and no further defect
    /// solve was admitted. A small recursive residual is not success.
    DefectCorrectionBudget,
}

/// A field accepted ONLY for the named linear goal. This intentionally is
/// not a ConductionSolution: small mean/selected-node error cannot stand in
/// for whole-field convergence, energy balance or a moving maximum.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearGoalSolve {
    /// Best checked physical temperature field; prescribed values unchanged.
    pub temperature: Vec<f64>,
    /// The enclosure of this exact returned field, not of an inner correction.
    pub analysis: LinearGoalAnalysis,
    /// Goal success or an explicit unresolved stopping reason.
    pub stop: LinearGoalStop,
    /// All primal iterations spent, including rejected candidates and retries.
    pub primal_iterations: usize,
    /// Extra defect solves actually started.
    pub defect_corrections: usize,
    /// Completed outward checks, including the initial field.
    pub goal_checks: usize,
}

impl LinearGoalAnalyzer<'_> {
    /// Improve an admissible field until its absolute algebraic GOAL error is
    /// small enough, rather than oversolving to an unrelated residual target.
    ///
    /// Reuses this analyzer's assembled operator, dual and checked stability
    /// proposal. The dual/stability work reported in `analysis` is preparation
    /// work, not repeated inside the primal loop. A loose goal may require no
    /// primal solve at all. Missing inverse evidence never triggers a fallback
    /// that silently claims the requested goal accuracy.
    ///
    /// Each defect RHS is normalized before PCG, then the correction is added
    /// to the physical temperature with FMA. The FULL returned field is checked
    /// after that rounding. Only improving goal enclosures replace the best
    /// iterate; retries share the original explicit iteration budget.
    ///
    /// # Errors
    /// Invalid policy, field/assembly-domain failures, nonfinite arithmetic
    /// and cancellation. Budget or numerical nonconvergence returns an explicit
    /// non-success state with the best checked field, never a forged success.
    pub fn solve_to_goal(
        &self,
        cx: &Cx<'_>,
        initial_temperature: &[f64],
        config: LinearGoalSolveConfig,
    ) -> Result<LinearGoalSolve, ConductionError> {
        self.solve_to_goal_observed(cx, initial_temperature, config, |_, _| {})
    }

    /// Same solve with a callback after each completed candidate enclosure.
    /// A callback may report progress or request cancellation through its Cx
    /// gate. Cancellation is checked AFTER the callback and before publication,
    /// even when that candidate would otherwise satisfy the goal. Candidates
    /// need not improve; the returned field is the best checked one.
    ///
    /// # Errors
    /// The same refusals as [`Self::solve_to_goal`].
    #[allow(clippy::too_many_lines)] // Keep the shared budget, recurrence lifetime and publication gate visible together.
    pub fn solve_to_goal_observed(
        &self,
        cx: &Cx<'_>,
        initial_temperature: &[f64],
        config: LinearGoalSolveConfig,
        mut observe: impl FnMut(usize, &LinearGoalAnalysis),
    ) -> Result<LinearGoalSolve, ConductionError> {
        poll(cx, 0)?;
        if !(config.absolute_tolerance.is_finite() && config.absolute_tolerance > 0.0)
            || !(1..=32).contains(&config.check_every)
        {
            return Err(invalid("goal solve requires a finite positive absolute tolerance and check_every in 1..=32"));
        }
        let analysis = self.analyze(cx, initial_temperature)?;
        observe(0, &analysis);
        poll(cx, 0)?;
        let mut result = LinearGoalSolve {
            temperature: initial_temperature.to_vec(), analysis,
            stop: LinearGoalStop::IterationBudget, primal_iterations: 0,
            defect_corrections: 0, goal_checks: 1,
        };
        if result.analysis.meets_absolute_tolerance(config.absolute_tolerance) {
            result.stop = LinearGoalStop::GoalTolerance;
            poll(cx, 0)?;
            return Ok(result);
        }
        let Some(bound) = result.analysis.enclosure.goal_error() else {
            result.stop = LinearGoalStop::BoundUnavailable;
            poll(cx, 0)?;
            return Ok(result);
        };
        if config.max_primal_iterations == 0 {
            poll(cx, 0)?;
            return Ok(result);
        }
        let mut best_bound = bound.magnitude_upper();
        let op = CsrOp::symmetric(self.response.matrix.clone());
        let pre = crate::solve::spd_preconditioner(&self.response.matrix);
        'attempt: loop {
            poll(cx, result.primal_iterations)?;
            let base = self.response.dofs.gather(&result.temperature);
            let mut defect = Vec::with_capacity(base.len());
            let mut scale = 0.0_f64;
            let mut visited = 0_usize;
            for (row, &rhs) in self.rhs.iter().enumerate() {
                if row % 512 == 0 { poll(cx, result.primal_iterations)?; }
                let (columns, values) = self.response.matrix.row(row);
                let mut r = rhs;
                for (&column, &a) in columns.iter().zip(values) {
                    if visited % 512 == 0 { poll(cx, result.primal_iterations)?; }
                    visited = visited.wrapping_add(1);
                    r = finite((-a).mul_add(base[column], r))?;
                }
                scale = scale.max(r.abs());
                defect.push(r);
            }
            if scale == 0.0 {
                result.stop = LinearGoalStop::NoProgress;
                break;
            }
            for (index, value) in defect.iter_mut().enumerate() {
                if index % 512 == 0 { poll(cx, result.primal_iterations)?; }
                *value /= scale;
            }
            let start_bound = best_bound;
            let mut state = CgState::new(&op, &pre, &defect);
            loop {
                poll(cx, result.primal_iterations)?;
                let before = state.iters;
                let batch = config.check_every.min(config.max_primal_iterations - result.primal_iterations);
                // This cutoff only ends an exhausted recurrence. It cannot
                // admit a goal, a field, or a residual-convergence claim.
                state.run(&op, &pre, f64::MIN_POSITIVE, batch);
                result.primal_iterations += state.iters - before;
                poll(cx, result.primal_iterations)?;
                if state.x.iter().all(|v| v.is_finite()) {
                    let mut candidate = self.response.dofs.prescribed().to_vec();
                    for (i, &vertex) in self.response.dofs.free().iter().enumerate() {
                        if i % 512 == 0 { poll(cx, result.primal_iterations)?; }
                        candidate[vertex] = finite(scale.mul_add(state.x[i], base[i]))?;
                    }
                    let checked = self.analyze(cx, &candidate)?;
                    result.goal_checks = result.goal_checks.saturating_add(1);
                    observe(result.primal_iterations, &checked);
                    poll(cx, result.primal_iterations)?;
                    if let Some(candidate_bound) = checked.enclosure.goal_error() {
                        if candidate_bound.magnitude_upper() < best_bound {
                            best_bound = candidate_bound.magnitude_upper();
                            result.temperature = candidate;
                            result.analysis = checked;
                        }
                    }
                    if result.analysis.meets_absolute_tolerance(config.absolute_tolerance) {
                        result.stop = LinearGoalStop::GoalTolerance;
                        break 'attempt;
                    }
                } else {
                    // The previous best field is still finite and enclosed.
                    // Do not publish this broken inner iterate.
                    break;
                }
                if result.primal_iterations == config.max_primal_iterations {
                    result.stop = LinearGoalStop::IterationBudget;
                    break 'attempt;
                }
                if state.iters == before || !state.rel_residual().is_finite()
                    || state.rel_residual() <= f64::MIN_POSITIVE
                {
                    break;
                }
            }
            if result.primal_iterations == config.max_primal_iterations {
                result.stop = LinearGoalStop::IterationBudget;
                break;
            }
            if best_bound >= start_bound {
                result.stop = LinearGoalStop::NoProgress;
                break;
            }
            if result.defect_corrections == config.max_defect_corrections {
                result.stop = LinearGoalStop::DefectCorrectionBudget;
                break;
            }
            result.defect_corrections += 1;
        }
        poll(cx, result.primal_iterations)?;
        Ok(result)
    }
}

fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) }
    else { Err(ConductionError::NonFinite { field: "goal-controlled primal", bits: value.to_bits() }) }
}
