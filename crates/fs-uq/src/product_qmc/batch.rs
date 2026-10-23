//! Concurrent QMC without changing the original net or between-scramble estimator.
//! Drivers own scheduling and must drain all jobs before returning their slots.
mod checkpoint;

use super::{QmcConfig, QmcExecution, QmcReport};
use crate::{UqBatchInput, UqPlan, UqStatus};
use core::fmt::Display;

/// Paid work is distinct from inference. Only the complete fixed design has a result.
#[derive(Debug, Clone, PartialEq)]
pub struct QmcBatchReport {
    pub config: QmcConfig,
    pub samples_planned: usize,
    pub evaluations_attempted: usize,
    pub evaluations_accepted: usize,
    pub contiguous_observations: usize,
    pub retained_out_of_order: usize,
    pub status: UqStatus,
    pub rejection_reason: Option<String>,
    /// Original QMC report, never an iid pointwise Monte Carlo standard error.
    pub result: Option<QmcReport>,
}

/// Bounded, recoverable waves over one original replicated QMC experiment.
///
/// Every ordinal is replicate-major and sampled by `QmcExecution::parameters`.
/// Responses are published to that same owner's reducer in ordinal order. A
/// scheduler may finish points in any order but must keep each outcome in its
/// job's slot. No new wave opens while a previous wave has missing values.
/// Complete siblings survive interruption, cloning and durable restoration.
///
/// At most 256 pending slots and the original one-million-call cap are retained.
/// The callback must represent the same deterministic physical model on every
/// invocation. Worker identity and completion time are not model parameters.
/// Solver resource admission and cancellation remain the driver's responsibility.
/// No shortened-net inference, optional stopping, accuracy or speedup is claimed.
#[derive(Debug, Clone)]
pub struct QmcBatchExecution {
    scalar: QmcExecution,
    pending: Vec<Option<f64>>,
}

impl QmcBatchExecution {
    pub fn new(plan: &UqPlan, config: QmcConfig) -> Result<Self, &'static str> {
        let mut scalar = QmcExecution::new(plan, config)?;
        scalar.values.try_reserve_exact(plan.budget_max_samples)
            .map_err(|_| "cannot reserve QMC batch observations")?;
        let mut pending = Vec::new();
        pending.try_reserve_exact(256).map_err(|_| "cannot reserve QMC batch slots")?;
        Ok(Self { scalar, pending })
    }

    pub fn plan(&self) -> &UqPlan { self.scalar.plan() }
    pub const fn config(&self) -> QmcConfig { self.scalar.config() }
    pub const fn evaluations_attempted(&self) -> usize { self.scalar.attempted }
    /// Contiguous original-ordinal prefix, not a selected set of completed jobs.
    pub fn observations(&self) -> &[f64] { &self.scalar.values }
    pub fn observation(&self, ordinal: usize) -> Option<f64> {
        let n = self.scalar.values.len();
        if ordinal < n { Some(self.scalar.values[ordinal]) }
        else { self.pending.get(ordinal - n).copied().flatten() }
    }

    pub fn report(&self) -> QmcBatchReport {
        let retained = self.pending.iter().filter(|x| x.is_some()).count();
        QmcBatchReport {
            config: self.scalar.config,
            samples_planned: self.scalar.plan.budget_max_samples,
            evaluations_attempted: self.scalar.attempted,
            evaluations_accepted: self.scalar.values.len() + retained,
            contiguous_observations: self.scalar.values.len(),
            retained_out_of_order: retained,
            status: self.scalar.status,
            rejection_reason: self.scalar.failure.clone(),
            result: (self.scalar.status == UqStatus::Complete).then(|| self.scalar.report()),
        }
    }

    /// Dispatch at most `allowance` missing physical evaluations in ONE wave.
    ///
    /// `None` means unfinished/unlaunched, never a discarded inconvenient result.
    /// `Some(Err(_))` and nonfinite responses permanently refuse the experiment.
    /// All finite siblings are retained before reporting the lowest failed ordinal
    /// or cancellation. Drivers must join/drain work before returning, converting
    /// caught worker panics to errors. Uncaught panics propagate.
    ///
    /// Width may change after interruption. Zero allowance and terminal runs do
    /// not poll or dispatch. Interrupted costs are owned by the physical solver;
    /// only completed or rejected callbacks spend an evaluation here.
    pub fn advance_batch<E: Display>(
        &mut self, allowance: usize, mut cancelled: impl FnMut() -> bool,
        evaluate: impl FnOnce(&[UqBatchInput], &mut [Option<Result<f64, E>>]),
    ) -> Result<QmcBatchReport, &'static str> {
        if allowance > 256 { return Err("QMC batch allowance must be in 0..=256"); }
        if allowance == 0 || matches!(self.scalar.status, UqStatus::Complete | UqStatus::Refused) {
            return Ok(self.report());
        }
        if cancelled() { return Ok(self.pause()); }
        if self.pending.is_empty() {
            let count = allowance.min(self.scalar.plan.budget_max_samples - self.scalar.values.len());
            self.pending.try_reserve_exact(count).map_err(|_| "cannot reserve QMC pending slots")?;
            self.pending.resize(count, None);
        }
        // Clones need not keep spare capacity. Admit publication before physics.
        self.scalar.values.try_reserve_exact(self.pending.len())
            .map_err(|_| "cannot reserve QMC batch publication")?;
        let mut jobs = Vec::new();
        let mut outcomes = Vec::new();
        jobs.try_reserve_exact(allowance).map_err(|_| "cannot reserve QMC batch inputs")?;
        outcomes.try_reserve_exact(allowance).map_err(|_| "cannot reserve QMC batch outcomes")?;
        for (offset, value) in self.pending.iter().enumerate() {
            if value.is_some() { continue; }
            if cancelled() { return Ok(self.pause()); }
            let ordinal = self.scalar.values.len() + offset;
            let parameters = match self.scalar.parameters(ordinal) {
                Ok(parameters) => parameters,
                Err(reason) => return Ok(self.fail(format!("QMC sample {ordinal}: {reason}"))),
            };
            jobs.push(UqBatchInput { ordinal, parameters });
            outcomes.push(None);
            if jobs.len() == allowance { break; }
        }
        if cancelled() { return Ok(self.pause()); }
        evaluate(&jobs, &mut outcomes);
        let mut failure = None;
        let mut unfinished = false;
        for (job, outcome) in jobs.iter().zip(outcomes) {
            let Some(outcome) = outcome else { unfinished = true; continue; };
            self.scalar.attempted += 1;
            match outcome {
                Ok(value) if value.is_finite() => {
                    self.pending[job.ordinal - self.scalar.values.len()] = Some(value);
                }
                Ok(_) => {
                    if failure.is_none() { failure = Some(format!("QMC model evaluation {} returned a non-finite QoI", job.ordinal)); }
                }
                Err(error) => {
                    if failure.is_none() { failure = Some(format!("QMC model evaluation {} refused: {error}", job.ordinal)); }
                }
            }
        }
        let prefix = self.pending.iter().take_while(|x| x.is_some()).count();
        self.scalar.values.extend(self.pending.drain(..prefix).map(|x| x.expect("finite QMC prefix")));
        if let Some(reason) = failure { return Ok(self.fail(reason)); }
        if self.scalar.values.len() == self.scalar.plan.budget_max_samples {
            self.scalar.status = UqStatus::Complete;
        } else if unfinished || cancelled() {
            return Ok(self.pause());
        } else {
            self.scalar.status = UqStatus::BudgetTruncated;
        }
        Ok(self.report())
    }

    fn pause(&mut self) -> QmcBatchReport {
        self.scalar.status = UqStatus::Cancelled;
        self.report()
    }
    fn fail(&mut self, reason: String) -> QmcBatchReport {
        self.scalar.status = UqStatus::Refused;
        self.scalar.failure = Some(reason);
        self.report()
    }
}
