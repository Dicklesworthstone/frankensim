//! Concurrent physical Shapley evaluations using the original conditional law.
//! Completion order never changes the row, physical coalition or pair half.
mod checkpoint;

use super::{ShapleyExecution, ShapleyReport};
use crate::{UqBatchInput, UqPlan, UqStatus};
use core::fmt::Display;

/// Paid physical work, separate from a complete-design sensitivity estimate.
#[derive(Debug, Clone, PartialEq)]
pub struct ShapleyBatchReport {
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
    /// Original physical Shapley report, including undefined normalization.
    /// Absent until the ENTIRE fixed design succeeds, not just a few pairs.
    pub result: Option<ShapleyReport>,
}

/// Bounded concurrent conditional pairs for the original Gaussian experiment.
///
/// The scalar owner admits every permuted conditional factor, samples every
/// physical input and reduces the complete design. No independent-input swap,
/// covariance inverse or new germ stream replaces the declared dependence.
/// Waves may split pairs and cross coalition/row boundaries. Successful later
/// calls survive interruption, cloning and disk recovery in their original slots.
/// Only missing slots are dispatched; no later wave bypasses a missing one.
///
/// Storage retains the original one-million-response cap and <=256 pending
/// slots. Drivers own physical solver resources and MUST drain before returning.
/// They must evaluate the same deterministic model and match outcomes to jobs.
/// No partial attribution, adaptive stopping or measured speedup is implied.
#[derive(Debug, Clone)]
pub struct ShapleyBatchExecution {
    scalar: ShapleyExecution,
    pending: Vec<Option<f64>>,
}

impl ShapleyBatchExecution {
    pub fn new(plan: &UqPlan) -> Result<Self, &'static str> {
        let scalar = ShapleyExecution::new(plan)?;
        let mut pending = Vec::new();
        pending.try_reserve_exact(256).map_err(|_| "cannot reserve Shapley batch slots")?;
        Ok(Self { scalar, pending })
    }

    pub const fn plan(&self) -> &UqPlan { self.scalar.plan() }
    pub const fn evaluations_attempted(&self) -> usize { self.scalar.attempted }
    /// Contiguous row/subset/half prefix, never a selected set of fast workers.
    pub fn observations(&self) -> &[f64] { &self.scalar.values }
    pub fn observation(&self, ordinal: usize) -> Option<f64> {
        let n = self.scalar.values.len();
        if ordinal < n { Some(self.scalar.values[ordinal]) }
        else { self.pending.get(ordinal - n).copied().flatten() }
    }

    pub fn report(&self) -> ShapleyBatchReport {
        let retained = self.pending.iter().filter(|value| value.is_some()).count();
        ShapleyBatchReport {
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

    /// Dispatch at most allowance missing evaluations in ONE wave (0..=256).
    ///
    /// Slots start as None. Some(Ok(finite)) retains an original response; None
    /// means unfinished/unlaunched work, never an inconvenient observed result.
    /// Some(Err) and nonfinite responses permanently refuse the experiment.
    /// All finite siblings are retained before the lowest failed ordinal or
    /// cancellation is reported. Caught worker panics must become errors;
    /// unhandled driver panics propagate. Callback side effects must tolerate
    /// retry of unfinished work, whose cost the driver must meter separately.
    ///
    /// Cancellation is polled while preparing inputs and after the drained wave.
    /// Zero allowance and terminal states neither poll nor dispatch. Width may
    /// change after resumption without changing the design or refunding work.
    pub fn advance_batch<E: Display>(
        &mut self, allowance: usize, mut cancelled: impl FnMut() -> bool,
        evaluate: impl FnOnce(&[UqBatchInput], &mut [Option<Result<f64, E>>]),
    ) -> Result<ShapleyBatchReport, &'static str> {
        if allowance > 256 { return Err("Shapley batch allowance must be in 0..=256"); }
        if allowance == 0 || matches!(self.scalar.status, UqStatus::Complete | UqStatus::Refused) {
            return Ok(self.report());
        }
        if cancelled() { return Ok(self.pause()); }
        if self.pending.is_empty() {
            let count = allowance.min(self.scalar.plan.budget_max_samples - self.scalar.values.len());
            self.pending.try_reserve_exact(count).map_err(|_| "cannot reserve Shapley pending slots")?;
            self.pending.resize(count, None);
        }
        // Clone may discard spare capacity. Admit publication before physics.
        self.scalar.values.try_reserve_exact(self.pending.len())
            .map_err(|_| "cannot reserve Shapley batch publication")?;
        let mut jobs = Vec::new();
        let mut outcomes = Vec::new();
        jobs.try_reserve_exact(allowance).map_err(|_| "cannot reserve Shapley batch inputs")?;
        outcomes.try_reserve_exact(allowance).map_err(|_| "cannot reserve Shapley batch outcomes")?;
        for offset in 0..self.pending.len() {
            if self.pending[offset].is_some() { continue; }
            if cancelled() { return Ok(self.pause()); }
            let ordinal = self.scalar.values.len() + offset;
            let mut parameters = Vec::new();
            parameters.try_reserve_exact(self.scalar.plan.parameters.len())
                .map_err(|_| "cannot reserve Shapley physical input")?;
            match self.scalar.parameters_at(ordinal) {
                Ok(values) => parameters.extend_from_slice(values),
                Err(reason) => return Ok(self.fail(format!("Shapley sample {ordinal}: {reason}"))),
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
                    if failure.is_none() { failure = Some(format!("Shapley physical evaluation {} returned a non-finite QoI", job.ordinal)); }
                }
                Err(error) => {
                    if failure.is_none() { failure = Some(format!("Shapley physical evaluation {} refused: {error}", job.ordinal)); }
                }
            }
        }
        let prefix = self.pending.iter().take_while(|value| value.is_some()).count();
        self.scalar.values.extend(self.pending.drain(..prefix).map(|value| value.expect("finite Shapley prefix")));
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

    fn pause(&mut self) -> ShapleyBatchReport {
        self.scalar.status = UqStatus::Cancelled;
        self.report()
    }
    fn fail(&mut self, reason: String) -> ShapleyBatchReport {
        self.scalar.status = UqStatus::Refused;
        self.scalar.failure = Some(reason);
        self.report()
    }
}
