//! Feasibility-ranked population fitting through the original thermal model.
//!
//! Own the immutable physical requirements and the search together. A fitting
//! target never bypasses an energy/state/hot-spot limit. Failed physics remains
//! a typed error, rather than an infeasible score or a fabricated penalty.
//! One physical-evaluation slot is held back for replay of the selected design.

/// Whole-study restart, including the independent physical-work account.
pub mod restart;

use fs_ascent::scan_design::constraints::value::{ScanConstraintValues, ValueConstrainedScanStudy};
use fs_ascent::scan_design::value::ScanPowerValueEvaluation;
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::ScanSensitivityError;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::{
    ScanObservation, spatial::ScanObservationSpec,
};
use fs_conduction::bc::history::march::moving::MovingHeatSource;
use fs_dfo::differential_evolution::{DifferentialEvolutionConfig, DifferentialEvolutionError,
    DifferentialEvolutionRunError, DifferentialEvolutionStop};
use fs_dfo::differential_evolution::parallel::DifferentialEvolutionParallelError;
use fs_dfo::differential_evolution::constrained::{ConstrainedDifferentialEvolutionPoint,
    ConstrainedDifferentialEvolutionReport, ConstrainedDifferentialEvolutionState,
    ConstrainedEvaluation, FeasibilityScore};
use fs_exec::{Cx, TilePool};
use fs_alloc::OperationMemoryLease;

pub type PopulationScanError = DifferentialEvolutionRunError<ScanSensitivityError>;

/// Preserve admission and drained parallel execution failures separately.
#[derive(Debug)]
pub enum PopulationScanParallelError {
    Admission(ScanSensitivityError),
    Search(DifferentialEvolutionParallelError<ScanSensitivityError>),
}
impl std::fmt::Display for PopulationScanParallelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Admission(error) => write!(f, "parallel scan admission: {error}"),
            Self::Search(error) => std::fmt::Display::fmt(error, f),
        }
    }
}
impl std::error::Error for PopulationScanParallelError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Admission(error) => Some(error), Self::Search(error) => Some(error) }
    }
}

/// Search status and physical work are distinct; callback admission can fail
/// before a physical reservation. Neither a pause nor budget exhaustion is a
/// constrained-stationarity or infeasibility certificate.
#[derive(Debug)]
pub struct PopulationScanReport {
    pub search: ConstrainedDifferentialEvolutionReport,
    pub physical_evaluations: usize,
    pub physical_evaluation_limit: usize,
    /// Excludes the one remaining slot reserved for a final physical replay.
    pub remaining_physical_search_evaluations: usize,
}
impl PopulationScanReport {
    pub fn physical_budget_exhausted(&self) -> bool {
        self.search.stop == DifferentialEvolutionStop::Paused
            && self.remaining_physical_search_evaluations == 0
    }
}

/// An independently repeated physical evaluation, not a cached fitness value.
/// Feasibility and target mismatch remain separate. The original model's mesh,
/// time grid and material declarations still bound every physical claim.
#[derive(Debug)]
pub struct PopulationScanReplay {
    pub coordinates: Vec<f64>,
    pub score: FeasibilityScore,
    pub objective_target: Option<f64>,
    pub physical: ScanPowerValueEvaluation,
    pub constraints: ScanConstraintValues,
    pub sources: Vec<MovingHeatSource>,
    pub search_evaluations: usize,
    pub physical_evaluations: usize,
}
impl PopulationScanReplay {
    pub fn target_recovered(&self) -> bool {
        self.score.is_feasible()
            && self.objective_target.is_some_and(|target| self.score.objective() <= target)
    }
}

/// In-process constrained search with serial and explicitly shared execution.
/// The serial API accepts the original FnMut physical factory. Parallel calls
/// additionally require Fn + Sync and Sync observations: a mutable-only factory
/// is never silently cloned or protected by a global physical-solver mutex.
///
/// Search boxes use the study's optimizer coordinates, not automatically SI.
/// For nonnegative-linear powers with reference_w=1 they are watts. Source
/// reconstruction admits both corners before any physical callback. Footprint
/// support and material-range failures still propagate from the real solver.
/// No public mutable study/search access can replace requirements during a run.
pub struct PopulationScanStudy<F, O = ScanObservation> {
    study: ValueConstrainedScanStudy<F, O>,
    search: ConstrainedDifferentialEvolutionState,
    target: Option<f64>,
}
impl<F, O: ScanObservationSpec> PopulationScanStudy<F, O> {
    /// Constructor performs no physics. Existing study probes remain charged.
    /// The remaining physical allowance must cover every admitted search
    /// attempt plus one final replay, even when a probe later fails early.
    pub fn new(cx: &Cx<'_>, study: ValueConstrainedScanStudy<F, O>, bounds: &[[f64; 2]],
        config: DifferentialEvolutionConfig, feasibility_tolerance: f64)
        -> Result<Self, PopulationScanError>
    {
        check(cx)?;
        admit_work(0, config.max_evaluations, study.evaluations(), study.evaluation_limit())?;
        if bounds.len() != study.initial_coordinates().len() {
            return Err(invalid("population box must have one interval per selected scan coordinate").into());
        }
        let target = config.target;
        let search = ConstrainedDifferentialEvolutionState::new(cx, bounds,
            Some(study.initial_coordinates()), config, study.constraint_count(), feasibility_tolerance)?;
        let lower: Vec<_> = bounds.iter().map(|bound| bound[0]).collect();
        let upper: Vec<_> = bounds.iter().map(|bound| bound[1]).collect();
        study.sources_at(cx, &lower).map_err(DifferentialEvolutionRunError::Evaluation)?;
        study.sources_at(cx, &upper).map_err(DifferentialEvolutionRunError::Evaluation)?;
        check(cx)?;
        Ok(Self { study, search, target })
    }

    pub fn best(&self) -> Option<&ConstrainedDifferentialEvolutionPoint> { self.search.best() }
    pub fn population(&self) -> &[ConstrainedDifferentialEvolutionPoint] { self.search.population() }
    pub fn search_evaluations(&self) -> usize { self.search.evaluations() }
    pub fn physical_evaluations(&self) -> usize { self.study.evaluations() }
    pub fn search_evaluation_limit(&self) -> usize { self.search.evaluation_limit() }
    pub fn physical_evaluation_limit(&self) -> usize { self.study.evaluation_limit() }
    pub fn feasibility_tolerance(&self) -> f64 { self.search.feasibility_tolerance() }

    /// Keep partial initialization/generations and all spent attempts. Validate
    /// both lifetime limits before changing either one. Increasing just one is
    /// allowed if the other already covers the requested work and replay.
    pub fn extend_evaluation_budgets(&mut self, search_total: usize, physical_total: usize)
        -> Result<(), PopulationScanError>
    {
        if search_total < self.search.evaluation_limit() || physical_total < self.study.evaluation_limit()
            || (search_total == self.search.evaluation_limit() && physical_total == self.study.evaluation_limit())
        { return Err(invalid("population allowances may only increase, with at least one changed").into()); }
        admit_work(self.search.evaluations(), search_total, self.study.evaluations(), physical_total)?;
        // The original extension owners only refuse a non-increasing limit;
        // those conditions have already been checked, with no fallible work after.
        if physical_total > self.study.evaluation_limit() {
            self.study.extend_evaluation_budget(physical_total)
                .map_err(DifferentialEvolutionRunError::Evaluation)?;
        }
        if search_total > self.search.evaluation_limit() { self.search.extend_evaluation_budget(search_total)?; }
        Ok(())
    }

    /// Run original physical trajectories and pass ALL scaled inequalities to
    /// the existing feasibility-first comparator. Objective magnitude cannot
    /// compensate for violating a constraint. Callback and physical budgets
    /// remain separate and failures consume their original reservations.
    pub fn advance<'h, 'm, 'c>(&mut self, cx: &Cx<'_>, additional_evaluations: usize)
        -> Result<PopulationScanReport, PopulationScanError>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        check(cx)?;
        let allowance = additional_evaluations.min(self.remaining_search_work());
        let study = &mut self.study;
        let search = self.search.advance(cx, allowance, &mut |point| {
            let (physical, constraints) = study.evaluate_at(cx, point)?;
            Ok::<_, ScanSensitivityError>(ConstrainedEvaluation { objective: physical.physical.value,
                inequalities: constraints.residuals })
        })?;
        Ok(PopulationScanReport { search, physical_evaluations: self.study.evaluations(),
            physical_evaluation_limit: self.study.evaluation_limit(),
            remaining_physical_search_evaluations: self.remaining_search_work() })
    }

    /// Evaluate bounded waves on the original TilePool with all physical
    /// requirements attached. One borrowed batch keeps source controls, model,
    /// observations and constraint scales fixed, while each worker constructs
    /// its own run from its own Cx. Started physical reservations settle back
    /// into the original study before success OR an error returns.
    ///
    /// Successful sibling scores remain in the DE cache on failure/cancellation;
    /// either advance API can continue without repeating them. Target checks
    /// occur between waves. The reserved serial replay is never dispatched.
    /// The executor lease does not cap callback-owned physical heap allocations.
    pub fn advance_parallel<'h, 'm, 'c>(&mut self, cx: &Cx<'_>, pool: &TilePool,
        lease: &OperationMemoryLease, additional_evaluations: usize, max_wave: usize)
        -> Result<PopulationScanReport, PopulationScanParallelError>
    where F: Fn(&Cx<'_>, Vec<MovingHeatSource>)
            -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError> + Sync,
        O: Sync,
    {
        let allowance = additional_evaluations.min(self.remaining_search_work())
            .min(self.search.evaluation_limit().saturating_sub(self.search.evaluations()));
        let attempted = {
            let batch = self.study.value_batch(cx, allowance)
                .map_err(PopulationScanParallelError::Admission)?;
            self.search.advance_parallel(cx, pool, lease, allowance, max_wave, &|worker, point| {
                let (physical, constraints) = batch.evaluate_at(worker, point)?;
                Ok::<_, ScanSensitivityError>(ConstrainedEvaluation {
                    objective: physical.physical.value, inequalities: constraints.residuals,
                })
            })
            // The pool has drained. Drop settles even failed reservations here,
            // before the outer Result propagates and before another advance.
        };
        let search = attempted.map_err(PopulationScanParallelError::Search)?;
        Ok(PopulationScanReport { search, physical_evaluations: self.study.evaluations(),
            physical_evaluation_limit: self.study.evaluation_limit(),
            remaining_physical_search_evaluations: self.remaining_search_work() })
    }

    /// Re-run every attached physical requirement at the best-ranked design.
    /// A best-so-far may still be infeasible; its report says so. The repeated
    /// score must match exactly before any sources or report are returned.
    /// This spends a physical evaluation without refunding search work. Extra
    /// replays can reduce future search capacity; extend budgets explicitly.
    pub fn replay_best<'h, 'm, 'c>(&mut self, cx: &Cx<'_>)
        -> Result<PopulationScanReplay, PopulationScanError>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        check(cx)?;
        let best = self.search.best().ok_or_else(|| invalid("no evaluated population design to replay"))?;
        let coordinates = best.coordinates.clone();
        let expected = best.value;
        let (physical, constraints) = self.study.evaluate_at(cx, &coordinates)
            .map_err(DifferentialEvolutionRunError::Evaluation)?;
        let score = self.search.score(physical.physical.value, &constraints.residuals)?;
        if score.objective().to_bits() != expected.objective().to_bits()
            || score.max_scaled_violation().to_bits() != expected.max_scaled_violation().to_bits()
            || score.is_feasible() != expected.is_feasible()
        { return Err(invalid("physical replay changed the selected objective or feasibility score").into()); }
        let sources = self.study.sources_at(cx, &coordinates)
            .map_err(DifferentialEvolutionRunError::Evaluation)?;
        check(cx)?;
        Ok(PopulationScanReplay { coordinates, score, objective_target: self.target,
            physical, constraints, sources, search_evaluations: self.search.evaluations(),
            physical_evaluations: self.study.evaluations() })
    }

    fn remaining_search_work(&self) -> usize {
        self.study.evaluation_limit().saturating_sub(self.study.evaluations()).saturating_sub(1)
    }
}

fn admit_work(search_spent: usize, search_limit: usize, physical_spent: usize, physical_limit: usize)
    -> Result<(), DifferentialEvolutionError>
{
    let needed = search_limit.checked_sub(search_spent).and_then(|remaining| remaining.checked_add(1));
    let available = physical_limit.checked_sub(physical_spent);
    match (needed, available) {
        (Some(needed), Some(available)) if needed <= available => Ok(()),
        _ => Err(invalid("physical budget must cover remaining population attempts and one replay")),
    }
}
fn invalid(message: &'static str) -> DifferentialEvolutionError { DifferentialEvolutionError::Invalid(message) }
fn check(cx: &Cx<'_>) -> Result<(), PopulationScanError> {
    cx.checkpoint().map_err(|_| DifferentialEvolutionError::Cancelled.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_replay_is_reserved_without_refunding_prior_physical_probes() {
        assert!(admit_work(0, 10, 3, 14).is_ok());
        assert!(admit_work(0, 10, 3, 13).is_err());
        assert!(admit_work(10, 10, 13, 14).is_ok());
        assert!(admit_work(10, 10, 14, 14).is_err());
    }
    #[test]
    fn g0_spent_counters_and_overflow_cannot_manufacture_search_capacity() {
        assert!(admit_work(11, 10, 0, 100).is_err());
        assert!(admit_work(0, 0, 11, 10).is_err());
        assert!(admit_work(0, usize::MAX, 0, usize::MAX).is_err());
        assert!(admit_work(usize::MAX, usize::MAX, usize::MAX - 1, usize::MAX).is_ok());
    }
}
