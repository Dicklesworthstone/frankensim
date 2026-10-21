//! Optimize one design against every declared physical case, without a surrogate.
use fs_ascent::scan_design::constraints::value::scenarios::{ScanScenarioStudy,
    ScanScenarioEvaluation, ScanScenarioError};
use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{ScanSensitivityError,
    objective::{ScanObservation, spatial::ScanObservationSpec}};
use fs_conduction::bc::history::march::moving::MovingHeatSource;
use fs_dfo::differential_evolution::{DifferentialEvolutionConfig, DifferentialEvolutionError,
    DifferentialEvolutionRunError};
use fs_dfo::differential_evolution::constrained::{ConstrainedDifferentialEvolutionState,
    ConstrainedDifferentialEvolutionReport, ConstrainedEvaluation, FeasibilityScore};
use fs_dfo::differential_evolution::parallel::DifferentialEvolutionParallelError;
use fs_alloc::OperationMemoryLease;
use fs_exec::{Cx, TilePool};

pub type ScenarioPopulationError = DifferentialEvolutionRunError<ScanScenarioError>;
#[derive(Debug)]
pub enum ScenarioPopulationParallelError {
    Admission(ScanScenarioError),
    Search(DifferentialEvolutionParallelError<ScanScenarioError>),
}
impl std::fmt::Display for ScenarioPopulationParallelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { Self::Admission(e) => write!(f,"{e}"), Self::Search(e) => write!(f,"{e}") }
    }
}
impl std::error::Error for ScenarioPopulationParallelError {}

#[derive(Debug)]
pub struct ScenarioPopulationReport {
    pub search: ConstrainedDifferentialEvolutionReport,
    /// One physical trajectory reservation per case, not one per ensemble.
    pub physical_evaluations: Vec<usize>,
    pub physical_evaluation_limits: Vec<usize>,
    /// Capacity after reserving one final replay in EVERY case.
    pub remaining_physical_search_candidates: usize,
}

#[derive(Debug)]
pub struct ScenarioPopulationReplay {
    pub coordinates: Vec<f64>,
    pub score: FeasibilityScore,
    pub objective_target: Option<f64>,
    pub evaluation: ScanScenarioEvaluation,
    pub sources: Vec<MovingHeatSource>,
    pub search_evaluations: usize,
    pub physical_evaluations: Vec<usize>,
}
impl ScenarioPopulationReplay {
    /// This is the search's declared aggregate target, not necessarily every
    /// individual case target under a weighted-mean objective.
    pub fn aggregate_target_met(&self) -> bool {
        self.score.is_feasible() && self.objective_target.is_some_and(|t| self.score.objective() <= t)
    }
    pub fn every_case_target_met(&self) -> bool {
        self.score.is_feasible() && self.objective_target.is_some_and(|t|
            self.evaluation.cases.iter().all(|c| c.physical.physical.value <= t))
    }
}

/// The existing constrained DE worker owns all mutation, selection, caching and
/// scheduling. This adapter supplies one original physical multi-case callback.
/// Requirements cannot change during continuation; failed cases stay errors.
pub struct ScenarioPopulationStudy<F, O = ScanObservation> {
    study: ScanScenarioStudy<F, O>,
    search: ConstrainedDifferentialEvolutionState,
    target: Option<f64>,
}
impl<F, O: ScanObservationSpec> ScenarioPopulationStudy<F, O> {
    pub fn new(cx: &Cx<'_>, study: ScanScenarioStudy<F, O>, bounds: &[[f64;2]],
        config: DifferentialEvolutionConfig, feasibility_tolerance: f64)
        -> Result<Self, ScenarioPopulationError>
    {
        check(cx)?;
        require_capacity(config.max_evaluations, study.remaining_evaluations())?;
        let target=config.target;
        let search=ConstrainedDifferentialEvolutionState::new(cx,bounds,Some(study.initial_coordinates()),
            config,study.constraint_count(),feasibility_tolerance)?;
        study.sources_at(cx,&bounds.iter().map(|b|b[0]).collect::<Vec<_>>())
            .map_err(DifferentialEvolutionRunError::Evaluation)?;
        study.sources_at(cx,&bounds.iter().map(|b|b[1]).collect::<Vec<_>>())
            .map_err(DifferentialEvolutionRunError::Evaluation)?;
        Ok(Self {study,search,target})
    }

    pub fn case_names(&self) -> &[String] { self.study.case_names() }
    pub fn search_evaluations(&self) -> usize { self.search.evaluations() }
    pub fn physical_evaluations(&self) -> Vec<usize> { self.study.evaluations() }
    pub fn remaining_search_candidates(&self) -> usize { self.study.remaining_evaluations().saturating_sub(1) }

    /// Preserve every spent counter and partial generation. All case totals
    /// and the search ceiling are checked before any allowance changes.
    pub fn extend_evaluation_budgets(&mut self, search_total: usize, physical_totals: &[usize])
        -> Result<(), ScenarioPopulationError>
    {
        let limits=self.study.evaluation_limits();
        let spent=self.study.evaluations();
        if search_total<self.search.evaluation_limit() || physical_totals.len()!=limits.len()
            || physical_totals.iter().zip(&limits).any(|(n,old)|n<old)
            || (search_total==self.search.evaluation_limit() && physical_totals==limits.as_slice()) {
            return Err(invalid("multi-case work allowances may only increase and must name every case").into());
        }
        let remaining=search_total.checked_sub(self.search.evaluations())
            .ok_or_else(||invalid("multi-case search counter exceeds its proposed limit"))?;
        for (&total,&used) in physical_totals.iter().zip(&spent) {
            require_capacity(remaining,total.checked_sub(used)
                .ok_or_else(||invalid("case counter exceeds its proposed limit"))?)?;
        }
        if physical_totals!=limits.as_slice() {
            self.study.extend_evaluation_budgets(physical_totals).map_err(DifferentialEvolutionRunError::Evaluation)?;
        }
        if search_total>self.search.evaluation_limit() { self.search.extend_evaluation_budget(search_total)?; }
        Ok(())
    }

    pub fn advance<'h,'m,'c>(&mut self,cx:&Cx<'_>,additional:usize)
        -> Result<ScenarioPopulationReport,ScenarioPopulationError>
    where F:FnMut(&Cx<'_>,Vec<MovingHeatSource>)->Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        check(cx)?;
        let allowance=additional.min(self.remaining_search_candidates());
        let study=&mut self.study;
        let search=self.search.advance(cx,allowance,&mut |point| {
            let result=study.evaluate_at(cx,point)?;
            Ok::<_,ScanScenarioError>(ConstrainedEvaluation {objective:result.objective,inequalities:result.residuals})
        })?;
        Ok(self.report(search))
    }

    pub fn advance_parallel<'h,'m,'c>(&mut self,cx:&Cx<'_>,pool:&TilePool,lease:&OperationMemoryLease,
        additional:usize,max_wave:usize)->Result<ScenarioPopulationReport,ScenarioPopulationParallelError>
    where F:Fn(&Cx<'_>,Vec<MovingHeatSource>)->Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>+Sync,
        O:Sync,
    {
        let allowance=additional.min(self.remaining_search_candidates())
            .min(self.search.evaluation_limit().saturating_sub(self.search.evaluations()));
        let attempted={
            let batch=self.study.value_batch(cx,allowance).map_err(ScenarioPopulationParallelError::Admission)?;
            self.search.advance_parallel(cx,pool,lease,allowance,max_wave,&|worker,point| {
                let result=batch.evaluate_at(worker,point)?;
                Ok::<_,ScanScenarioError>(ConstrainedEvaluation {objective:result.objective,inequalities:result.residuals})
            })
            // Workers have drained; unused reservations settle in every case,
            // including on a failed candidate, before propagating the Result.
        };
        let search=attempted.map_err(ScenarioPopulationParallelError::Search)?;
        Ok(self.report(search))
    }

    /// Repeat all original cases at the selected design. A score match is not
    /// a stationarity or uncertainty-coverage certificate. The immutable model
    /// must remain deterministic; cached DE scores do not retain all old fields.
    pub fn replay_best<'h,'m,'c>(&mut self,cx:&Cx<'_>)
        -> Result<ScenarioPopulationReplay,ScenarioPopulationError>
    where F:FnMut(&Cx<'_>,Vec<MovingHeatSource>)->Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        check(cx)?;
        let best=self.search.best().ok_or_else(||invalid("no evaluated multi-case design to replay"))?;
        let coordinates=best.coordinates.clone();
        let expected=best.value;
        let evaluation=self.study.evaluate_at(cx,&coordinates).map_err(DifferentialEvolutionRunError::Evaluation)?;
        let score=self.search.score(evaluation.objective,&evaluation.residuals)?;
        if score.objective().to_bits()!=expected.objective().to_bits()
            || score.max_scaled_violation().to_bits()!=expected.max_scaled_violation().to_bits()
            || score.is_feasible()!=expected.is_feasible() {
            return Err(invalid("multi-case physical replay changed the selected score").into());
        }
        let sources=self.study.sources_at(cx,&coordinates).map_err(DifferentialEvolutionRunError::Evaluation)?;
        check(cx)?;
        Ok(ScenarioPopulationReplay {coordinates,score,objective_target:self.target,evaluation,sources,
            search_evaluations:self.search.evaluations(),physical_evaluations:self.study.evaluations()})
    }
    fn report(&self,search:ConstrainedDifferentialEvolutionReport)->ScenarioPopulationReport {
        ScenarioPopulationReport {search,physical_evaluations:self.study.evaluations(),
            physical_evaluation_limits:self.study.evaluation_limits(),
            remaining_physical_search_candidates:self.remaining_search_candidates()}
    }
}
fn require_capacity(search_remaining:usize,physical_remaining:usize)->Result<(),DifferentialEvolutionError> {
    if search_remaining.checked_add(1).is_none_or(|n|n>physical_remaining) {
        Err(invalid("each case must cover remaining search attempts and one physical replay"))
    } else {Ok(())}
}
fn invalid(s:&'static str)->DifferentialEvolutionError {DifferentialEvolutionError::Invalid(s)}
fn check(cx:&Cx<'_>)->Result<(),ScenarioPopulationError> {
    cx.checkpoint().map_err(|_|DifferentialEvolutionError::Cancelled.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_every_case_needs_an_independent_replay_reservation() {
        assert!(require_capacity(10,11).is_ok());
        assert!(require_capacity(10,10).is_err());
        assert!(require_capacity(0,1).is_ok());
        assert!(require_capacity(0,0).is_err());
        assert!(require_capacity(usize::MAX,usize::MAX).is_err());
    }
}
