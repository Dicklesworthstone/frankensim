//! Parallel physical sensitivity without changing the original pick-freeze design.
//! Drivers drain all jobs and return each outcome in their original ordinal slot.
mod checkpoint;

use super::{SobolExecution, SobolReport};
use crate::{UqBatchInput, UqPlan, UqStatus};
use core::fmt::Display;

/// Paid work is not an estimate. Only the entire declared design has a result.
#[derive(Debug, Clone, PartialEq)]
pub struct SobolBatchReport {
    pub evaluations_planned: usize,
    pub evaluations_attempted: usize,
    pub evaluations_accepted: usize,
    pub contiguous_observations: usize,
    pub retained_out_of_order: usize,
    pub base_samples: usize,
    pub evaluations_per_row: usize,
    pub completed_rows: usize,
    pub status: UqStatus,
    pub rejection_reason: Option<String>,
    /// The original Jansen report, including its undefined-variance diagnostic.
    pub result: Option<SobolReport>,
}

/// Bounded waves of A, B and hybrid physical solves under one immutable plan.
///
/// Sampling, cached A/B rows and final Jansen reductions belong to the original
/// SobolExecution. Rows may cross wave boundaries; neither worker identity nor
/// completion order enters the design. A successful later sibling survives an
/// earlier interruption, cloning and durable restoration. Only missing slots
/// are dispatched, and a new wave cannot bypass an unfinished one.
///
/// The original one-million-call limit and at most 256 pending slots bound
/// storage. Drivers own solver resources, cancellation and draining. They must
/// evaluate the same deterministic physical model after every resumption. This
/// adds no correlated-input interpretation, stopping rule or speedup guarantee.
#[derive(Debug, Clone)]
pub struct SobolBatchExecution {
    scalar: SobolExecution,
    pending: Vec<Option<f64>>,
}

impl SobolBatchExecution {
    pub fn new(plan: &UqPlan) -> Result<Self, &'static str> {
        let scalar = SobolExecution::new(plan)?;
        let mut pending = Vec::new();
        pending.try_reserve_exact(256).map_err(|_| "cannot reserve Sobol batch slots")?;
        Ok(Self { scalar, pending })
    }

    pub const fn plan(&self) -> &UqPlan { self.scalar.plan() }
    pub const fn evaluations_attempted(&self) -> usize { self.scalar.attempted }
    /// Contiguous row-major observations, never a selected set of fast workers.
    pub fn observations(&self) -> &[f64] { &self.scalar.values }
    pub fn observation(&self, ordinal: usize) -> Option<f64> {
        let n = self.scalar.values.len();
        if ordinal < n { Some(self.scalar.values[ordinal]) }
        else { self.pending.get(ordinal - n).copied().flatten() }
    }

    pub fn report(&self) -> SobolBatchReport {
        let retained = self.pending.iter().filter(|x| x.is_some()).count();
        SobolBatchReport {
            evaluations_planned: self.scalar.plan.budget_max_samples,
            evaluations_attempted: self.scalar.attempted,
            evaluations_accepted: self.scalar.values.len() + retained,
            contiguous_observations: self.scalar.values.len(),
            retained_out_of_order: retained,
            base_samples: self.scalar.rows,
            evaluations_per_row: self.scalar.width,
            completed_rows: self.scalar.values.len() / self.scalar.width,
            status: self.scalar.status,
            rejection_reason: self.scalar.failure.clone(),
            result: (self.scalar.status == UqStatus::Complete).then(|| self.scalar.report()),
        }
    }

    /// Dispatch at most allowance missing physical evaluations in ONE wave.
    ///
    /// Each outcome slot is initially None. Some(Ok(finite)) retains a response;
    /// None means unfinished/unlaunched work, not a discarded observed result.
    /// Errors and nonfinite values permanently refuse the experiment. Retain
    /// all finite siblings before reporting the lowest failed ordinal or a
    /// cancellation. Drivers must drain before returning and convert caught
    /// worker panics to errors; unhandled driver panics propagate.
    ///
    /// Allowance is in 0..=256 and can change after resumption. Zero allowance
    /// and terminal runs do not dispatch or poll. Interrupted attempts do not
    /// consume ordinals; their costs must be metered by the physical solver.
    pub fn advance_batch<E: Display>(
        &mut self, allowance: usize, mut cancelled: impl FnMut() -> bool,
        evaluate: impl FnOnce(&[UqBatchInput], &mut [Option<Result<f64, E>>]),
    ) -> Result<SobolBatchReport, &'static str> {
        if allowance > 256 { return Err("Sobol batch allowance must be in 0..=256"); }
        if allowance == 0 || matches!(self.scalar.status, UqStatus::Complete | UqStatus::Refused) {
            return Ok(self.report());
        }
        if cancelled() { return Ok(self.pause()); }
        if self.pending.is_empty() {
            let count = allowance.min(self.scalar.plan.budget_max_samples - self.scalar.values.len());
            self.pending.try_reserve_exact(count).map_err(|_| "cannot reserve Sobol pending slots")?;
            self.pending.resize(count, None);
        }
        // A clone need not preserve spare capacity. Admit publication before physics.
        self.scalar.values.try_reserve_exact(self.pending.len())
            .map_err(|_| "cannot reserve Sobol batch publication")?;
        let mut jobs = Vec::new();
        let mut outcomes = Vec::new();
        jobs.try_reserve_exact(allowance).map_err(|_| "cannot reserve Sobol batch inputs")?;
        outcomes.try_reserve_exact(allowance).map_err(|_| "cannot reserve Sobol batch outcomes")?;
        for offset in 0..self.pending.len() {
            if self.pending[offset].is_some() { continue; }
            if cancelled() { return Ok(self.pause()); }
            let ordinal = self.scalar.values.len() + offset;
            let mut parameters = Vec::new();
            parameters.try_reserve_exact(self.scalar.plan.parameters.len())
                .map_err(|_| "cannot reserve Sobol physical input")?;
            match self.scalar.parameters(ordinal) {
                Ok(values) => parameters.extend_from_slice(values),
                Err(reason) => return Ok(self.fail(format!("Sobol sample {ordinal}: {reason}"))),
            }
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
                    if failure.is_none() { failure = Some(format!("Sobol model evaluation {} returned a non-finite QoI", job.ordinal)); }
                }
                Err(error) => {
                    if failure.is_none() { failure = Some(format!("Sobol model evaluation {} refused: {error}", job.ordinal)); }
                }
            }
        }
        let prefix = self.pending.iter().take_while(|x| x.is_some()).count();
        self.scalar.values.extend(self.pending.drain(..prefix).map(|x| x.expect("finite Sobol prefix")));
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

    fn pause(&mut self) -> SobolBatchReport {
        self.scalar.status = UqStatus::Cancelled;
        self.report()
    }
    fn fail(&mut self, reason: String) -> SobolBatchReport {
        self.scalar.status = UqStatus::Refused;
        self.scalar.failure = Some(reason);
        self.report()
    }
}
