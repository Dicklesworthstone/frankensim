//! Feasibility-first DE on the existing serial/parallel search worker.
//!
//! Every inequality is a finite, caller-scaled g(x) <= 0. A fixed explicit
//! tolerance determines feasibility. Feasible candidates outrank ALL infeasible
//! candidates regardless of objective magnitude. Among infeasible candidates,
//! minimize the largest positive scaled violation, then the original objective.
//! Among feasible candidates, minimize only the objective. There is no penalty
//! coefficient, objective rescaling, discarded constraint or invented fitness.
//!
//! This is a declared max-violation policy, not a sum-violation formulation.
//! Row scaling matters and must remain fixed throughout the search. Feasibility
//! concerns only the supplied rows; it is not a continuum/error certificate.

/// Save and restore drained searches, retaining their feasibility policy.
pub mod snapshot;

use super::{DifferentialEvolutionConfig, DifferentialEvolutionError as Error,
    DifferentialEvolutionPoint, DifferentialEvolutionReport, DifferentialEvolutionRunError,
    DifferentialEvolutionScore, RankedDifferentialEvolutionState, finite, invalid};
use super::parallel::DifferentialEvolutionParallelError;
use fs_alloc::OperationMemoryLease;
use fs_exec::{Cx, TilePool};
use std::cmp::Ordering;

/// One successful physical callback. All declared rows must be present, in a
/// fixed order and scaling. A physical failure belongs in Result::Err, not here.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstrainedEvaluation {
    pub objective: f64,
    pub inequalities: Vec<f64>,
}

/// Original objective and constraint summary, never a scalar penalty score.
/// Constructed only by the search's immutable feasibility policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FeasibilityScore {
    objective: f64,
    maximum_violation: f64,
    feasible: bool,
}
impl FeasibilityScore {
    pub const fn objective(self) -> f64 { self.objective }
    /// Actual max(0,g_i), before subtracting or applying the tolerance.
    pub const fn max_scaled_violation(self) -> f64 { self.maximum_violation }
    pub const fn is_feasible(self) -> bool { self.feasible }
}
impl DifferentialEvolutionScore for FeasibilityScore {
    fn checked(self) -> Result<Self, Error> {
        finite(self.objective)?;
        finite(self.maximum_violation)?;
        if self.maximum_violation < 0.0 { return Err(invalid("DE violation must be nonnegative")); }
        Ok(self)
    }
    fn rank(self, other: Self) -> Ordering {
        match (self.feasible, other.feasible) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            (true, true) => self.objective.total_cmp(&other.objective),
            (false, false) => self.maximum_violation.total_cmp(&other.maximum_violation)
                .then_with(|| self.objective.total_cmp(&other.objective)),
        }
    }
    // Equal ranks accept the trial, preserving exploration on flat responses.
    fn accepts(self, parent: Self) -> bool { !self.rank(parent).is_gt() }
    fn reaches(self, target: f64) -> bool { self.feasible && self.objective <= target }
}

/// Report always separates objective, violation and feasible/infeasible status.
/// An exhausted or paused best candidate may still be infeasible.
pub type ConstrainedDifferentialEvolutionReport = DifferentialEvolutionReport<FeasibilityScore>;
pub type ConstrainedDifferentialEvolutionPoint = DifferentialEvolutionPoint<FeasibilityScore>;

/// Penalty-free constrained selection using the original DE candidate stream,
/// frozen generations, callback accounting, cancellation and parallel cache.
///
/// The number of rows and their scaled feasibility tolerance are fixed at
/// construction. Caller-owned model/row semantics must also remain unchanged.
/// An objective target is reached ONLY by a feasible candidate. Budget
/// exhaustion proves neither optimality nor absence of a feasible design.
/// Checkpoints retain the original feasibility policy and completed callbacks.
/// Their distinct format refuses scalar DE snapshots and changed requirements.
#[derive(Debug)]
pub struct ConstrainedDifferentialEvolutionState {
    state: RankedDifferentialEvolutionState<FeasibilityScore>,
    constraint_count: usize,
    feasibility_tolerance: f64,
}
impl ConstrainedDifferentialEvolutionState {
    /// Admit original box/search controls plus a positive row count and finite
    /// nonnegative tolerance. No objective is invoked. Every callback must
    /// return exactly constraint_count rows. The count also bounds the scalar
    /// reduction work; callback-owned vectors/physics memory need their own cap.
    pub fn new(cx: &Cx<'_>, bounds: &[[f64; 2]], initial: Option<&[f64]>,
        config: DifferentialEvolutionConfig, constraint_count: usize,
        feasibility_tolerance: f64) -> Result<Self, Error>
    {
        Self::new_with(bounds, initial, config, constraint_count, feasibility_tolerance,
            &mut || cx.checkpoint().is_ok())
    }

    fn new_with(bounds: &[[f64; 2]], initial: Option<&[f64]>,
        config: DifferentialEvolutionConfig, constraint_count: usize,
        feasibility_tolerance: f64, poll: &mut impl FnMut() -> bool) -> Result<Self, Error>
    {
        super::checkpoint(poll)?;
        if constraint_count == 0 || !feasibility_tolerance.is_finite() || feasibility_tolerance < 0.0 {
            return Err(invalid("constrained DE needs a positive row count and finite nonnegative feasibility tolerance"));
        }
        Ok(Self { state: RankedDifferentialEvolutionState::new_with(bounds, initial, config, poll)?,
            constraint_count, feasibility_tolerance })
    }

    pub const fn constraint_count(&self) -> usize { self.constraint_count }
    pub const fn feasibility_tolerance(&self) -> f64 { self.feasibility_tolerance }
    pub const fn evaluations(&self) -> usize { self.state.evaluations() }
    pub const fn evaluation_limit(&self) -> usize { self.state.evaluation_limit() }
    pub const fn generations(&self) -> usize { self.state.generations() }
    pub fn best(&self) -> Option<&ConstrainedDifferentialEvolutionPoint> { self.state.best() }
    pub fn population(&self) -> &[ConstrainedDifferentialEvolutionPoint] { self.state.population() }
    pub fn extend_evaluation_budget(&mut self, total: usize) -> Result<(), Error> {
        self.state.extend_evaluation_budget(total)
    }

    /// Recheck a physical replay with this SAME count/tolerance. Does not spend
    /// a callback or change accepted search state. Reduction is one bounded
    /// O(constraint_count) unit so a successful callback can be cached even if
    /// it requests cancellation just before returning its finite result.
    pub fn score(&self, objective: f64, inequalities: &[f64]) -> Result<FeasibilityScore, Error> {
        score(objective, inequalities, self.constraint_count, self.feasibility_tolerance)
    }

    /// Serial callbacks share the same work caps and retained-candidate retries
    /// as scalar DE. Wrong row count or a nonfinite value consumes the attempt
    /// but publishes no score. No simulation error is converted into infeasibility.
    pub fn advance<E>(&mut self, cx: &Cx<'_>, additional_evaluations: usize,
        evaluate: &mut impl FnMut(&[f64]) -> Result<ConstrainedEvaluation, E>)
        -> Result<ConstrainedDifferentialEvolutionReport, DifferentialEvolutionRunError<E>>
    {
        self.advance_with(additional_evaluations, evaluate, &mut || cx.checkpoint().is_ok())
    }

    fn advance_with<E>(&mut self, additional: usize,
        evaluate: &mut impl FnMut(&[f64]) -> Result<ConstrainedEvaluation, E>,
        poll: &mut impl FnMut() -> bool)
        -> Result<ConstrainedDifferentialEvolutionReport, DifferentialEvolutionRunError<E>>
    {
        let (count, tolerance) = (self.constraint_count, self.feasibility_tolerance);
        self.state.advance_with(additional, &mut |point| {
            let value = evaluate(point).map_err(ScoreError::Callback)?;
            score(value.objective, &value.inequalities, count, tolerance).map_err(ScoreError::Invalid)
        }, poll).map_err(|error| match error {
            DifferentialEvolutionRunError::Solver(error) => DifferentialEvolutionRunError::Solver(error),
            DifferentialEvolutionRunError::Evaluation(ScoreError::Callback(error)) =>
                DifferentialEvolutionRunError::Evaluation(error),
            DifferentialEvolutionRunError::Evaluation(ScoreError::Invalid(error)) =>
                DifferentialEvolutionRunError::Solver(error),
        })
    }

    /// Production TilePool execution with ordered feasibility-first selection.
    /// Completed feasible AND infeasible results survive a failing sibling;
    /// serial or parallel continuation retries only missing scores. Target
    /// stopping is between waves and every started physical callback is charged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_parallel<E: Send, F>(&mut self, cx: &Cx<'_>, pool: &TilePool,
        lease: &OperationMemoryLease, additional_evaluations: usize, max_wave: usize,
        evaluate: &F)
        -> Result<ConstrainedDifferentialEvolutionReport, DifferentialEvolutionParallelError<E>>
    where F: Fn(&Cx<'_>, &[f64]) -> Result<ConstrainedEvaluation, E> + Sync,
    {
        let (count, tolerance) = (self.constraint_count, self.feasibility_tolerance);
        self.state.advance_parallel(cx, pool, lease, additional_evaluations, max_wave, &|worker, point| {
            let value = evaluate(worker, point).map_err(ScoreError::Callback)?;
            score(value.objective, &value.inequalities, count, tolerance).map_err(ScoreError::Invalid)
        }).map_err(|error| match error {
            DifferentialEvolutionParallelError::Solver(error) => DifferentialEvolutionParallelError::Solver(error),
            DifferentialEvolutionParallelError::Executor(error) => DifferentialEvolutionParallelError::Executor(error),
            DifferentialEvolutionParallelError::Evaluation { member, source: ScoreError::Callback(source) } =>
                DifferentialEvolutionParallelError::Evaluation { member, source },
            DifferentialEvolutionParallelError::Evaluation { source: ScoreError::Invalid(error), .. } =>
                DifferentialEvolutionParallelError::Solver(error),
        })
    }
}

fn score(objective: f64, inequalities: &[f64], count: usize, tolerance: f64) -> Result<FeasibilityScore, Error> {
    finite(objective)?;
    if inequalities.len() != count { return Err(invalid("constrained DE callback changed its inequality row count")); }
    let mut maximum_violation = 0.0_f64;
    for &row in inequalities { maximum_violation = maximum_violation.max(finite(row)?); }
    Ok(FeasibilityScore { objective, maximum_violation, feasible: maximum_violation <= tolerance })
}

enum ScoreError<E> { Callback(E), Invalid(Error) }

#[cfg(test)]
mod tests;
