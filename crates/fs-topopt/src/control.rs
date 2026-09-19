//! Bounded work for the real topology filter and elasticity solves.
//!
//! Reuses resumable `fs_solver::CgState`; no competing linear solver or
//! checkpoint format. Polls before allocation, every 32 Krylov iterations and
//! before publishing an evaluation. A matrix application is not preemptible.

use std::ops::ControlFlow;

use fs_solver::{CgState, LinearOp};
use fs_sparse::precond::IdentityPrecond;

/// Limits shared by every filter, load solve and pullback in one computation.
#[derive(Debug, Clone, Copy)]
pub struct SolveBudget {
    /// Maximum iterations in any one linear solve. Zero admits only zero RHS.
    pub per_solve_iterations: usize,
    /// Cumulative Krylov iterations, including rejected optimization trials.
    pub total_iterations: usize,
}

impl Default for SolveBudget {
    fn default() -> Self {
        Self { per_solve_iterations: 50_000, total_iterations: usize::MAX }
    }
}

/// Actual work, not a count of accepted optimization steps.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SolveWork {
    /// Started linear solves, including zero-RHS solves and interrupted trials.
    pub linear_solves: usize,
    /// Completed Krylov iterations across all started solves.
    pub linear_iterations: usize,
}

/// Context provided to cancellation/wall-time callbacks.
#[derive(Debug, Clone, Copy)]
pub struct SolveProgress {
    /// `filter`, `filter-transpose`, `elasticity`, or a named algebraic stage.
    pub stage: &'static str,
    /// Iterations in the active solve; zero outside a linear solve.
    pub solve_iterations: usize,
    /// Cumulative work, including rejected candidates.
    pub work: SolveWork,
}

/// A stopped evaluation never supplies a partial gradient as a usable result.
#[derive(Debug, Clone, PartialEq)]
pub enum EvaluationStop {
    /// Callback requested a stop, for example from an asupersync cancellation
    /// context or an expired wall budget.
    Cancelled,
    /// This solve used its full iteration allowance without convergence.
    LinearBudget { stage: &'static str, iterations: usize, residual_estimate: f64 },
    /// The computation's cumulative iteration allowance was consumed.
    TotalBudget { stage: &'static str },
    /// Non-finite arithmetic or failure of the linear recurrence.
    Breakdown { stage: &'static str },
}

impl std::fmt::Display for EvaluationStop {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "topology evaluation stopped: {self:?}")
    }
}

impl std::error::Error for EvaluationStop {}

/// One budget/accounting context, reused through a whole optimization.
/// The callback may bind `fs_exec::Cx::checkpoint`; no global cancel state.
/// Model admission retains the component APIs' panic-on-invalid-input contract.
pub struct SolveControl<'a> {
    budget: SolveBudget,
    work: SolveWork,
    callback: &'a mut dyn FnMut(SolveProgress) -> ControlFlow<()>,
}

impl<'a> SolveControl<'a> {
    /// Start fresh work accounting under explicit linear-iteration limits.
    pub fn new(
        budget: SolveBudget,
        callback: &'a mut dyn FnMut(SolveProgress) -> ControlFlow<()>,
    ) -> Self {
        Self { budget, work: SolveWork::default(), callback }
    }

    /// Work consumed so far, also available after a stopped evaluation.
    #[must_use]
    pub fn work(&self) -> SolveWork { self.work }

    pub(crate) fn checkpoint(&mut self, stage: &'static str) -> Result<(), EvaluationStop> {
        self.poll(stage, 0)
    }

    fn poll(&mut self, stage: &'static str, solve_iterations: usize) -> Result<(), EvaluationStop> {
        match (self.callback)(SolveProgress { stage, solve_iterations, work: self.work }) {
            ControlFlow::Continue(()) => Ok(()),
            ControlFlow::Break(()) => Err(EvaluationStop::Cancelled),
        }
    }

    pub(crate) fn solve(
        &mut self,
        op: &impl LinearOp,
        rhs: &[f64],
        tolerance: f64,
        component_limit: usize,
        stage: &'static str,
    ) -> Result<Vec<f64>, EvaluationStop> {
        self.poll(stage, 0)?;
        assert_eq!(rhs.len(), op.n(), "linear RHS shape mismatch");
        if !rhs.iter().all(|value| value.is_finite()) {
            return Err(EvaluationStop::Breakdown { stage });
        }
        self.work.linear_solves = self.work.linear_solves.checked_add(1)
            .ok_or(EvaluationStop::TotalBudget { stage })?;
        if rhs.iter().all(|value| *value == 0.0) {
            let result = vec![0.0; op.n()];
            self.poll(stage, 0)?;
            return Ok(result);
        }
        let limit = component_limit.min(self.budget.per_solve_iterations);
        if self.work.linear_iterations == self.budget.total_iterations {
            return Err(EvaluationStop::TotalBudget { stage });
        }
        if limit == 0 {
            return Err(EvaluationStop::LinearBudget {
                stage, iterations: 0, residual_estimate: 1.0,
            });
        }
        let mut state = CgState::new(op, &IdentityPrecond, rhs);
        loop {
            self.poll(stage, state.iters)?;
            let residual = state.rel_residual();
            if !residual.is_finite() || !state.x.iter().all(|value| value.is_finite()) {
                return Err(EvaluationStop::Breakdown { stage });
            }
            if residual < tolerance { return Ok(state.x); }
            if state.iters == limit {
                return Err(EvaluationStop::LinearBudget {
                    stage, iterations: state.iters, residual_estimate: residual,
                });
            }
            let remaining = self.budget.total_iterations - self.work.linear_iterations;
            if remaining == 0 { return Err(EvaluationStop::TotalBudget { stage }); }
            let batch = 32.min(limit - state.iters).min(remaining);
            let before = state.iters;
            let _ = state.run(op, &IdentityPrecond, tolerance, batch);
            self.work.linear_iterations += state.iters - before;
            // CgState's recurrence does not depend on its diagnostic history.
            // Do not repeatedly clone an ever-growing history in each batch.
            state.history.clear();
            if state.iters == before {
                return Err(EvaluationStop::Breakdown { stage });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn operator() -> fs_solver::CsrOp {
        let mut coo = fs_sparse::Coo::new(7, 7);
        for i in 0..7 { coo.push(i, i, 1.0 + i as f64); }
        fs_solver::CsrOp::symmetric(coo.assemble())
    }

    #[test]
    fn g5_batched_cg_matches_the_existing_uninterrupted_solver() {
        let op = operator();
        let rhs = [1.0; 7];
        let mut reference = CgState::new(&op, &IdentityPrecond, &rhs);
        assert!(reference.run(&op, &IdentityPrecond, 1e-12, 100).converged);
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let actual = control.solve(&op, &rhs, 1e-12, 100, "test").unwrap();
        assert_eq!(actual, reference.x);
        assert_eq!(control.work.linear_iterations, reference.iters);
    }

    #[test]
    fn g4_iteration_allowances_count_actual_work_across_solves() {
        let op = operator();
        let rhs = [1.0; 7];
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget {
            per_solve_iterations: 100, total_iterations: 8,
        }, &mut callback);
        control.solve(&op, &rhs, 1e-12, 100, "first").unwrap();
        assert_eq!(control.work.linear_iterations, 7);
        assert_eq!(control.solve(&op, &rhs, 1e-12, 100, "second"),
            Err(EvaluationStop::TotalBudget { stage: "second" }));
        assert_eq!(control.work.linear_iterations, 8);
    }

    #[test]
    fn g4_zero_budget_and_zero_rhs_have_distinct_results() {
        let op = operator();
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget {
            per_solve_iterations: 0, total_iterations: 100,
        }, &mut callback);
        assert_eq!(control.solve(&op, &[0.0; 7], 1e-12, 100, "zero").unwrap(), vec![0.0; 7]);
        assert!(matches!(control.solve(&op, &[1.0; 7], 1e-12, 100, "nonzero"),
            Err(EvaluationStop::LinearBudget { iterations: 0, .. })));
        assert_eq!(control.work.linear_iterations, 0);
    }

    #[test]
    fn g4_cancellation_after_real_krylov_work_does_not_publish_an_iterate() {
        let op = operator();
        let mut callback = |progress: SolveProgress| {
            if progress.solve_iterations > 0 { ControlFlow::Break(()) }
            else { ControlFlow::Continue(()) }
        };
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        assert_eq!(control.solve(&op, &[1.0; 7], 1e-12, 100, "test"), Err(EvaluationStop::Cancelled));
        assert!(control.work.linear_iterations > 0);
        assert!(control.work.linear_iterations <= 32);
    }
}
