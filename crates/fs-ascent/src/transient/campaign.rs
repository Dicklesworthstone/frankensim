//! Joint transient calibration across experiments sharing decision coordinates.
//!
//! Each experiment retains its own initial conditions, forcing, sensors, time
//! interval and numerical policy. Its forward solve and sampled discrete
//! adjoint are delegated to evaluate_transient. The weighted sum is passed to
//! the existing SQP engine, not a second optimizer or a concatenated artificial
//! ODE. Experiments run sequentially in stable ID order with bounded storage.

use super::{TransientConfig, TransientError, TransientEvaluation, TransientFamily,
    admit, box_sample, evaluate_transient};
use crate::sqp::{SqpError, SqpRunReport, SqpState, SqpStop};

pub struct Experiment<'a, F> {
    pub id: u64,
    pub family: &'a F,
    pub config: TransientConfig,
    /// Fixed positive scalar multiplying this experiment's entire objective.
    pub weight: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CampaignError {
    Invalid(&'static str),
    Experiment { id: u64, source: TransientError },
    TrialLimit,
    ExperimentLimit,
    NonFiniteSum,
    Cancelled,
}
impl std::fmt::Display for CampaignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "transient campaign failed: {self:?}")
    }
}
impl std::error::Error for CampaignError {}
pub type CampaignStudyError = SqpError<CampaignError>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CampaignWork { pub trials: usize, pub experiments: usize }

/// Caller-owned, non-cloneable cumulative allowance, retained even when study
/// initialization fails. Attempts are charged BEFORE entering a model factory;
/// failed/cancelled attempts are not refunded. Limits can increase, never reset.
/// These are callback counts, not seconds/FLOPs. Each charged experiment is
/// additionally bounded by its own forward/replay configuration.
pub struct CampaignControl { work: CampaignWork, max_trials: usize, max_experiments: usize }
impl CampaignControl {
    pub fn new(max_trials: usize, max_experiments: usize) -> Self {
        Self { work: CampaignWork { trials: 0, experiments: 0 }, max_trials, max_experiments }
    }
    pub fn work(&self) -> CampaignWork { self.work }
    pub fn extend(&mut self, max_trials: usize, max_experiments: usize) -> Result<(), CampaignError> {
        if max_trials < self.max_trials || max_experiments < self.max_experiments {
            return Err(CampaignError::Invalid("campaign allowances can only increase"));
        }
        self.max_trials = max_trials; self.max_experiments = max_experiments; Ok(())
    }
    fn begin(&mut self, experiments: usize) -> Result<(), CampaignError> {
        if self.work.trials >= self.max_trials { return Err(CampaignError::TrialLimit); }
        // Admit a complete trial before spending on its first experiment. A
        // numerical failure can still stop early; only attempted cases count.
        if experiments > self.max_experiments - self.work.experiments {
            return Err(CampaignError::ExperimentLimit);
        }
        self.work.trials += 1; Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExperimentEvaluation { pub id: u64, pub weight: f64, pub result: TransientEvaluation }
#[derive(Debug, Clone, PartialEq)]
pub struct CampaignEvaluation {
    pub point: Vec<f64>,
    pub value: f64,
    pub gradient: Vec<f64>,
    /// Unweighted results, in ascending ID order, at exactly `point`.
    pub experiments: Vec<ExperimentEvaluation>,
}

fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), CampaignError> {
    if cancelled() { Err(CampaignError::Cancelled) } else { Ok(()) }
}
fn case_error(id: u64, source: TransientError) -> CampaignError {
    if source == TransientError::Cancelled { CampaignError::Cancelled }
    else { CampaignError::Experiment { id, source } }
}
fn normalize(error: CampaignStudyError) -> CampaignStudyError {
    match error { SqpError::Evaluation(CampaignError::Cancelled) => SqpError::Cancelled, other => other }
}
fn same_bounds(a: &[[f64; 2]], b: &[[f64; 2]]) -> bool {
    a.len() == b.len() && a.iter().flatten().zip(b.iter().flatten()).all(|(a,b)| a.to_bits() == b.to_bits())
}
#[derive(Clone, Copy, Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, x: f64) -> Result<(), CampaignError> {
        let next = self.value + x;
        self.correction += if self.value.abs() >= x.abs() { (self.value-next)+x } else { (x-next)+self.value };
        self.value = next;
        if !self.value.is_finite() || !self.correction.is_finite() || !self.total().is_finite() {
            return Err(CampaignError::NonFiniteSum);
        }
        Ok(())
    }
    fn total(self) -> f64 { self.value + self.correction }
}

/// All families must use the SAME physical meaning and ordering for decisions,
/// not merely equal numerical boxes. That semantic mapping is caller-supplied;
/// exact bound equality is checked here. Families must remain pure/unchanged.
pub struct TransientCampaign<'a, F> {
    experiments: Vec<Experiment<'a, F>>,
    bounds: Vec<[f64; 2]>,
    max_kkt_dimension: usize,
}
impl<'a, F: TransientFamily> TransientCampaign<'a, F> {
    /// Cap the experiment count before traversal. `max_retained_components`
    /// bounds scalar entries in ONE retained CampaignEvaluation: its aggregate
    /// point/gradient plus each case's point/gradient and declared terminal
    /// state cap. Accepted and trial results may coexist; RK workspace, inputs,
    /// records, allocator metadata and callback-owned memory are separate.
    pub fn new(mut experiments: Vec<Experiment<'a, F>>, max_experiments: usize,
        max_retained_components: usize) -> Result<Self, CampaignError>
    {
        if experiments.is_empty() || experiments.len() > max_experiments {
            return Err(CampaignError::Invalid("nonempty experiment list must fit its cap"));
        }
        let first = &experiments[0];
        admit(first.family, &first.config).map_err(|e| case_error(first.id,e))?;
        let n = first.family.bounds().len();
        let mut retained = n.checked_mul(2).ok_or(CampaignError::Invalid("retained-size overflow"))?;
        let mut max_kkt_dimension = first.config.max_kkt_dimension;
        for experiment in &experiments {
            if !experiment.weight.is_finite() || experiment.weight <= 0.0 {
                return Err(CampaignError::Invalid("experiment weights must be finite and positive"));
            }
            if !same_bounds(first.family.bounds(),experiment.family.bounds()) {
                return Err(CampaignError::Invalid("experiments must share the exact decision box"));
            }
            admit(experiment.family, &experiment.config).map_err(|e| case_error(experiment.id,e))?;
            retained = retained.checked_add(2*n).and_then(|v| v.checked_add(experiment.config.max_state_components))
                .ok_or(CampaignError::Invalid("retained-size overflow"))?;
            if retained > max_retained_components {
                return Err(CampaignError::Invalid("campaign result exceeds retained scalar cap"));
            }
            max_kkt_dimension = max_kkt_dimension.min(experiment.config.max_kkt_dimension);
        }
        let bounds = first.family.bounds().to_vec();
        experiments.sort_by_key(|case| case.id);
        if experiments.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(CampaignError::Invalid("experiment IDs must be unique"));
        }
        Ok(Self { experiments, bounds, max_kkt_dimension })
    }
    pub fn bounds(&self) -> &[[f64; 2]] { &self.bounds }
    pub fn experiments(&self) -> &[Experiment<'a, F>] { &self.experiments }

    /// No partial case family is returned on failure. Ordinary box refusals
    /// spend no physical work; SQP separately counts those callback attempts.
    /// A failed partial family is recomputed on retry, with prior work retained.
    pub fn evaluate<C: FnMut() -> bool>(&self, point: &[f64], control: &mut CampaignControl,
        cancelled: &mut C) -> Result<Option<CampaignEvaluation>, CampaignError>
    {
        poll(cancelled)?;
        let n = self.bounds.len();
        if point.len() != n || point.iter().any(|v| !v.is_finite()) {
            return Err(CampaignError::Invalid("decision point must be finite and dimension-matched"));
        }
        if point.iter().zip(&self.bounds).any(|(x,b)| *x < b[0] || *x > b[1]) { return Ok(None); }
        for e in &self.experiments {
            poll(cancelled)?;
            if !same_bounds(e.family.bounds(), &self.bounds) {
                return Err(CampaignError::Invalid("experiment decision box changed"));
            }
            admit(e.family, &e.config).map_err(|error| case_error(e.id,error))?;
        }
        control.begin(self.experiments.len())?;
        let mut value = Sum::default(); let mut gradient = vec![Sum::default(); n];
        let mut experiments = Vec::with_capacity(self.experiments.len());
        for e in &self.experiments {
            poll(cancelled)?;
            control.work.experiments += 1;
            let result = evaluate_transient(e.family, &e.config, point, cancelled)
                .map_err(|error| case_error(e.id,error))?
                .ok_or(CampaignError::Invalid("experiment domain changed during evaluation"))?;
            value.add(e.weight * result.value)?;
            for (sum,g) in gradient.iter_mut().zip(&result.gradient) { sum.add(e.weight * g)?; }
            experiments.push(ExperimentEvaluation { id: e.id, weight: e.weight, result });
        }
        poll(cancelled)?;
        Ok(Some(CampaignEvaluation { point: point.to_vec(), value: value.total(),
            gradient: gradient.into_iter().map(Sum::total).collect(), experiments }))
    }
}

/// Shared-parameter SQP using one immutable campaign and an exclusively borrowed
/// cumulative work allowance. Local KKT convergence is not a proof of physical
/// validity, parameter identifiability, noise independence, or global optimality.
pub struct CampaignStudy<'a, 'work, F: TransientFamily> {
    campaign: &'a TransientCampaign<'a, F>,
    control: &'work mut CampaignControl,
    state: SqpState,
    accepted: CampaignEvaluation,
}
impl<'a, 'work, F: TransientFamily> CampaignStudy<'a, 'work, F> {
    pub fn new<C: FnMut() -> bool>(campaign: &'a TransientCampaign<'a,F>, point: &[f64],
        control: &'work mut CampaignControl, cancelled: &mut C) -> Result<Self, CampaignStudyError>
    {
        let mut accepted = None;
        let state = SqpState::try_new(point, campaign.max_kkt_dimension, &mut |x| {
            let Some(result) = campaign.evaluate(x, control, cancelled)? else { return Ok(None); };
            let sample = box_sample(&campaign.bounds, &result.point, result.value, &result.gradient);
            accepted = Some(result); Ok(Some(sample))
        }, None).map_err(normalize)?;
        poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
        let accepted = accepted.ok_or(SqpError::Invalid("missing complete initial campaign"))?;
        Ok(Self { campaign, control, state, accepted })
    }
    pub fn optimizer(&self) -> &SqpState { &self.state }
    pub fn accepted(&self) -> &CampaignEvaluation { &self.accepted }
    pub fn work(&self) -> CampaignWork { self.control.work() }
    pub fn extend_work(&mut self, trials: usize, experiments: usize) -> Result<(), CampaignError> {
        self.control.extend(trials, experiments)
    }
    pub fn run<C: FnMut() -> bool>(&mut self, tolerance: f64, additional_iterations: usize,
        maximum_evaluations: usize, cancelled: &mut C) -> Result<SqpRunReport, CampaignStudyError>
    {
        poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
        let mut report = self.advance(tolerance, 0, maximum_evaluations, cancelled)?;
        for _ in 0..additional_iterations {
            if report.stop != SqpStop::IterationLimit { break; }
            poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
            report = self.advance(tolerance, 1, maximum_evaluations, cancelled)?;
        }
        poll(cancelled).map_err(SqpError::Evaluation).map_err(normalize)?;
        Ok(report)
    }
    fn advance<C: FnMut() -> bool>(&mut self, tolerance: f64, steps: usize,
        maximum_evaluations: usize, cancelled: &mut C) -> Result<SqpRunReport, CampaignStudyError>
    {
        let campaign = self.campaign; let control = &mut *self.control; let mut candidate = None;
        let outcome = self.state.try_run(&mut |point| {
            let Some(result) = campaign.evaluate(point, control, cancelled)? else { return Ok(None); };
            let sample = box_sample(&campaign.bounds, &result.point, result.value, &result.gradient);
            candidate = Some(result); Ok(Some(sample))
        }, tolerance, steps, maximum_evaluations, None);
        if let Some(result) = candidate {
            if result.point.as_slice() == self.state.point() { self.accepted = result; }
        }
        outcome.map_err(normalize)
    }
}

#[cfg(test)]
mod tests;
