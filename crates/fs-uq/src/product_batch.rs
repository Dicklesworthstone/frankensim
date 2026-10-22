//! Ordered Monte Carlo observations from bounded, externally scheduled waves.
//! The original scalar owner supplies admission, sampling and final statistics.
//! No scheduler, worker identity or completion order enters the probability law.
mod checkpoint;

use core::fmt::Display;
use crate::product_execution::sample_parameters;
use crate::{UqExecution, UqPlan, UqResult, UqStatus};

/// One immutable physical input; outcome slot j belongs to job j, not its finish order.
#[derive(Debug, Clone, PartialEq)]
pub struct UqBatchInput {
    pub ordinal: usize,
    pub parameters: Vec<f64>,
}

/// Progress is separate from inference: only the whole fixed design has statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct UqBatchReport {
    pub samples_planned: usize,
    /// Completed or rejected calls. Unfinished/unlaunched slots do not spend an ordinal.
    pub evaluations_attempted: usize,
    /// Includes successful later ordinals retained behind an unfinished earlier call.
    pub evaluations_accepted: usize,
    pub contiguous_observations: usize,
    pub retained_out_of_order: usize,
    pub status: UqStatus,
    pub rejection_reason: Option<String>,
    pub result: Option<UqResult>,
}

/// Bounded waves for one immutable fixed-count Monte Carlo experiment.
///
/// A driver may use a caller-owned TilePool, scoped workers, or remote solvers.
/// It MUST drain its work before returning and put each outcome in its matching
/// slot. No new wave is opened until every member of the current wave completes.
/// Successful siblings survive interruption, cloning and checkpoint restoration.
/// Errors/nonfinite results permanently refuse the experiment, never get skipped.
///
/// One million observations and at most 256 pending slots/parameter vectors.
/// Callback-owned solver resources need their own admission. The callback must
/// remain deterministic from the original physical inputs, not worker identity.
/// No partial distribution, optional-stopping inference or speedup is claimed.
#[derive(Debug, Clone)]
pub struct UqBatchExecution {
    scalar: UqExecution,
    pending: Vec<Option<f64>>,
}

impl UqBatchExecution {
    pub fn new(plan: &UqPlan) -> Result<Self, &'static str> {
        let mut scalar = UqExecution::new(plan)?;
        // Reserve admitted storage before the first physical evaluation.
        scalar.values.try_reserve_exact(plan.budget_max_samples)
            .map_err(|_| "cannot reserve batch observation storage")?;
        let mut pending = Vec::new();
        pending.try_reserve_exact(256).map_err(|_| "cannot reserve batch pending storage")?;
        Ok(Self { scalar, pending })
    }

    pub const fn plan(&self) -> &UqPlan { &self.scalar.plan }
    pub fn evaluations_attempted(&self) -> usize { self.scalar.attempted }
    /// Contiguous original sample prefix, not a selected completed-worker sample.
    pub fn observations(&self) -> &[f64] { &self.scalar.values }
    /// Accepted value at an original ordinal, including a retained later sibling.
    pub fn observation(&self, ordinal: usize) -> Option<f64> {
        let n = self.scalar.values.len();
        if ordinal < n { Some(self.scalar.values[ordinal]) }
        else { self.pending.get(ordinal - n).copied().flatten() }
    }

    pub fn report(&self) -> UqBatchReport {
        let retained = self.pending.iter().filter(|x| x.is_some()).count();
        let result = (self.scalar.status == UqStatus::Complete).then(|| self.scalar.report());
        // The canonical scalar reduction may refuse unrepresentable full statistics.
        let refused = result.as_ref().filter(|r| r.status == UqStatus::Refused);
        UqBatchReport {
            samples_planned: self.scalar.plan.budget_max_samples,
            evaluations_attempted: self.scalar.attempted,
            evaluations_accepted: self.scalar.values.len() + retained,
            contiguous_observations: self.scalar.values.len(), retained_out_of_order: retained,
            status: refused.map_or(self.scalar.status, |_| UqStatus::Refused),
            rejection_reason: refused.map_or_else(|| self.scalar.failure.clone(), |r| r.rejection_reason.clone()),
            result: result.filter(|r| r.status != UqStatus::Refused),
        }
    }

    /// Dispatch at most allowance missing samples in ONE wave (0..=256).
    ///
    /// Each initially empty outcome slot accepts Some(Ok(finite_value)), a
    /// terminal Some(Err(error)), or remains None for unfinished/unlaunched work.
    /// None must never censor a result or a numerical refusal. A driver that
    /// catches worker panics must report them as Err; unhandled panics propagate.
    /// Drivers own internal cancellation, draining and unfinished-call costs.
    ///
    /// Retain ALL returned finite siblings before considering cancellation or
    /// errors. Report the lowest failed ordinal, independent of finish order.
    /// Cancellation is polled during input preparation and after the drained
    /// wave. Zero allowance and terminal executions never invoke the driver.
    /// Width can change on resume; only missing ordinals are resubmitted.
    pub fn advance_batch<E: Display>(
        &mut self, allowance: usize, mut cancelled: impl FnMut() -> bool,
        evaluate: impl FnOnce(&[UqBatchInput], &mut [Option<Result<f64, E>>]),
    ) -> Result<UqBatchReport, &'static str> {
        if allowance > 256 { return Err("batch allowance must be in 0..=256"); }
        if allowance == 0 || matches!(self.scalar.status, UqStatus::Complete | UqStatus::Refused) {
            return Ok(self.report());
        }
        if cancelled() { return Ok(self.pause()); }
        if self.pending.is_empty() {
            let count = allowance.min(self.scalar.plan.budget_max_samples - self.scalar.values.len());
            self.pending.try_reserve_exact(count).map_err(|_| "cannot reserve batch pending slots")?;
            self.pending.resize(count, None);
        }
        // Clone need not retain spare capacity. Re-admit publication storage
        // before dispatch, including every successful sibling that may drain.
        self.scalar.values.try_reserve_exact(self.pending.len())
            .map_err(|_| "cannot reserve completed batch publication")?;
        let mut jobs = Vec::new();
        let mut outcomes = Vec::new();
        jobs.try_reserve_exact(allowance).map_err(|_| "cannot reserve batch inputs")?;
        outcomes.try_reserve_exact(allowance).map_err(|_| "cannot reserve batch outcomes")?;
        for (offset, value) in self.pending.iter().enumerate() {
            if value.is_some() { continue; }
            if cancelled() { return Ok(self.pause()); }
            let ordinal = self.scalar.values.len() + offset;
            let parameters = match sample_parameters(&self.scalar.plan, self.scalar.factor.as_deref(), ordinal) {
                Ok(parameters) => parameters,
                Err(reason) => return Ok(self.fail(format!("sample {ordinal}: {reason}"))),
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
                    if failure.is_none() { failure = Some(format!("model evaluation {} returned a non-finite QoI", job.ordinal)); }
                }
                Err(error) => {
                    if failure.is_none() { failure = Some(format!("model evaluation {} refused: {error}", job.ordinal)); }
                }
            }
        }
        // Publish only in original ordinal order. Keep successes behind every hole.
        let prefix = self.pending.iter().take_while(|x| x.is_some()).count();
        self.scalar.values.extend(self.pending.drain(..prefix).map(|x| x.expect("finite prefix")));
        if let Some(reason) = failure { return Ok(self.fail(reason)); }
        if self.scalar.values.len() == self.scalar.plan.budget_max_samples {
            self.scalar.status = UqStatus::Complete;
            let report = self.report();
            if report.status == UqStatus::Refused {
                self.scalar.status = UqStatus::Refused;
                self.scalar.failure = report.rejection_reason.clone();
            }
            return Ok(report);
        }
        if unfinished || cancelled() { return Ok(self.pause()); }
        self.scalar.status = UqStatus::BudgetTruncated;
        Ok(self.report())
    }

    fn pause(&mut self) -> UqBatchReport {
        self.scalar.status = UqStatus::Cancelled;
        self.report()
    }
    fn fail(&mut self, reason: String) -> UqBatchReport {
        self.scalar.status = UqStatus::Refused;
        self.scalar.failure = Some(reason);
        self.report()
    }
}
