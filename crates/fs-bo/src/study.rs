//! Resumable noisy design studies for caller-owned simulation jobs.
//!
//! `ask` reads outstanding work, `tell` records one result without numerical
//! work, and `advance` closes a complete batch and selects the next one. Arrival
//! order never becomes GP training order. This is deterministic BATCH BO, not
//! asynchronous fantasy/pending-point acquisition. External job dispatch and
//! cancellation belong to the caller; this module owns no threads or I/O.

use crate::gp::Gp;
use crate::hyper::zero_noise_duplicates;
use crate::noisy::{NoisyBoConfig, NoisyBoReport, NoisyObservation};

/// Stable logical job identity and the exact requested coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct NoisyRequest {
    /// Caller-assigned study key. Use a distinct key for independent studies.
    /// This is an association check, not authentication.
    pub study_id: u64,
    /// Zero-based position in the eventual ordered observation history.
    pub evaluation: usize,
    /// Coordinates in the original input units; return these unchanged.
    pub x: Vec<f64>,
}

/// One selected job and its retained response, if received.
#[derive(Debug, Clone, PartialEq)]
pub struct PendingEvaluation {
    /// Immutable through the study's public API.
    pub request: NoisyRequest,
    /// Original observation, including its declared variance.
    pub observation: Option<NoisyObservation>,
}

/// Result of recording an external simulation response.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TellDisposition {
    /// First accepted response for this job.
    Recorded,
    /// An identical response was already retained. No state or work changed.
    Replayed,
}

/// Failures leave pending observations and the last complete report intact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoisyStudyError {
    /// The request names a different caller-assigned study key.
    ForeignStudy,
    /// No job with this index has been issued.
    UnknownEvaluation,
    /// Coordinates differ bitwise from the issued job.
    RequestMismatch,
    /// Non-finite observation/centered value or negative/non-finite variance.
    InvalidObservation,
    /// A previously accepted job has a different value or variance.
    ConflictingObservation,
    /// Not every selected job has a response yet.
    WaitingForObservations,
    /// The completed data do not admit a finite fixed-noise GP.
    InvalidModel,
    /// The continuation hook requested a pause.
    Cancelled,
    /// A cumulative numerical-attempt counter cannot be represented.
    WorkOverflow,
}

impl std::fmt::Display for NoisyStudyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "noisy study: {self:?}")
    }
}
impl std::error::Error for NoisyStudyError {}

/// Numerical ATTEMPTS, including cancelled or failed advances. Not flop counts
/// or objective-evaluation counts: simulations are owned by the caller.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NoisyStudyWork {
    /// Attempts to fit a complete observation prefix.
    pub model_fit_attempts: usize,
    /// Attempts to select a whole new batch, including discarded partial ones.
    pub acquisition_batch_attempts: usize,
}

/// All state needed to pause an in-memory study, including partial responses.
/// `clone()` is an independent checkpoint; it never repeats received physical
/// evaluations. No persistence format or cross-version replay is claimed.
#[derive(Debug, Clone)]
pub struct NoisyStudy {
    study_id: u64,
    config: NoisyBoConfig,
    dim: usize,
    iters: usize,
    report: NoisyBoReport,
    pending: Vec<PendingEvaluation>,
    complete: bool,
    work: NoisyStudyWork,
}

fn same_observation(a: NoisyObservation, b: NoisyObservation) -> bool {
    a.value.to_bits() == b.value.to_bits()
        && a.noise_variance.to_bits() == b.noise_variance.to_bits()
}

fn check_request(expected: &[f64], actual: &[f64]) -> Result<(), NoisyStudyError> {
    if expected.len() != actual.len()
        || expected.iter().zip(actual).any(|(a, b)| a.to_bits() != b.to_bits())
    {
        return Err(NoisyStudyError::RequestMismatch);
    }
    Ok(())
}

fn requests(study_id: u64, first: usize, points: Vec<Vec<f64>>) -> Vec<PendingEvaluation> {
    points.into_iter().enumerate().map(|(slot, x)| PendingEvaluation {
        request: NoisyRequest { study_id, evaluation: first + slot, x }, observation: None,
    }).collect()
}

impl NoisyStudy {
    /// Prepare the initial Sobol jobs without invoking a simulation or fitting
    /// a model. At most `n_init + iters * config.q` distinct jobs will be issued.
    ///
    /// # Panics
    /// Invalid configuration uses the same admission checks as `minimize_noisy`.
    pub fn new(study_id: u64, dim: usize, n_init: usize, iters: usize,
        config: &NoisyBoConfig) -> Self
    {
        config.validate_study(dim, n_init, iters);
        Self { study_id, config: config.clone(), dim, iters,
            report: NoisyBoReport { x: Vec::new(), observations: Vec::new(), incumbent_trace: Vec::new() },
            pending: requests(study_id, 0, config.initial_design(dim, n_init)),
            complete: false, work: NoisyStudyWork::default() }
    }

    /// Outstanding jobs in stable logical order. Repeated calls do no work and
    /// return the same IDs; the caller must not dispatch an in-flight ID twice.
    /// Empty while a complete received batch awaits `advance`, or when finished.
    #[must_use]
    pub fn ask(&self) -> Vec<NoisyRequest> {
        self.pending.iter().filter(|p| p.observation.is_none()).map(|p| p.request.clone()).collect()
    }

    /// All selected jobs, including already received responses. A failed
    /// simulation is retried by the caller using the same unanswered request;
    /// there is no fake objective value for failure or timeout.
    #[must_use]
    pub fn pending(&self) -> &[PendingEvaluation] { &self.pending }

    /// Last fully modeled prefix. Received results in the current batch remain
    /// in `pending()` until a successful advance, even after cancellation.
    #[must_use]
    pub fn report(&self) -> &NoisyBoReport { &self.report }

    /// True only after the final complete batch has been modeled successfully.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.complete }

    /// Cumulative attempts are not rolled back by cancellation or model failure.
    #[must_use]
    pub fn work(&self) -> NoisyStudyWork { self.work }

    /// Accept results in any order. Bit-identical delivery retries, even for
    /// past batches, are idempotent; conflicting retries never replace data.
    /// Validation finishes before any mutation. This does no model/solver work.
    pub fn tell(&mut self, request: &NoisyRequest, observation: NoisyObservation)
        -> Result<TellDisposition, NoisyStudyError>
    {
        if request.study_id != self.study_id { return Err(NoisyStudyError::ForeignStudy); }
        let committed = self.report.x.len();
        let (point, previous) = if request.evaluation < committed {
            (&self.report.x[request.evaluation], Some(self.report.observations[request.evaluation]))
        } else {
            let pending = self.pending.get(request.evaluation - committed)
                .ok_or(NoisyStudyError::UnknownEvaluation)?;
            (&pending.request.x, pending.observation)
        };
        check_request(point, &request.x)?;
        if !observation.value.is_finite() || !(observation.value - self.config.prior_mean).is_finite()
            || !observation.noise_variance.is_finite() || observation.noise_variance < 0.0
        {
            return Err(NoisyStudyError::InvalidObservation);
        }
        if let Some(previous) = previous {
            return if same_observation(previous, observation) { Ok(TellDisposition::Replayed) }
                else { Err(NoisyStudyError::ConflictingObservation) };
        }
        self.pending[request.evaluation - committed].observation = Some(observation);
        Ok(TellDisposition::Recorded)
    }

    /// Close a complete batch and atomically publish either the next batch or
    /// final report. No objective callbacks are invoked. Repeating after final
    /// completion is a no-op; an incomplete batch refuses without numerical work.
    pub fn advance(&mut self) -> Result<(), NoisyStudyError> {
        self.advance_controlled(&mut || true)
    }

    /// Pause before/after the dense fit and between greedy acquisition slots.
    /// The last report and ALL received responses survive a refused advance;
    /// retrying regenerates the same next jobs. Discarded numerical work remains
    /// counted in `work()`. One dense fit or a slot's bounded CMA-ES searches
    /// cannot be preempted here. Adapt a Cx checkpoint to this continuation hook.
    ///
    /// # Panics
    /// Inherits existing acquisition/posterior arithmetic assertions. Such a
    /// panic, like cancellation, precedes committing history or new requests.
    pub fn advance_controlled(&mut self, keep_going: &mut dyn FnMut() -> bool)
        -> Result<(), NoisyStudyError>
    {
        if self.complete { return Ok(()); }
        if self.pending.iter().any(|p| p.observation.is_none()) {
            return Err(NoisyStudyError::WaitingForObservations);
        }
        if !keep_going() { return Err(NoisyStudyError::Cancelled); }
        let mut next = self.report.clone();
        for pending in &self.pending {
            next.x.push(pending.request.x.clone());
            next.observations.push(pending.observation.expect("complete batch"));
        }
        self.work.model_fit_attempts = self.work.model_fit_attempts.checked_add(1)
            .ok_or(NoisyStudyError::WorkOverflow)?;
        let values: Vec<f64> = next.observations.iter().map(|o| o.value - self.config.prior_mean).collect();
        let variances: Vec<f64> = next.observations.iter().map(|o| o.noise_variance).collect();
        if zero_noise_duplicates(&next.x, &variances) { return Err(NoisyStudyError::InvalidModel); }
        let gp = Gp::try_fit_diag(&next.x, &values, self.config.kernel.clone(), &variances)
            .filter(|gp| gp.lml.is_finite()).ok_or(NoisyStudyError::InvalidModel)?;
        if !keep_going() { return Err(NoisyStudyError::Cancelled); }
        next.incumbent_trace.push(self.config.study_incumbent(&gp, &next.x));
        let stage = self.report.incumbent_trace.len();
        let complete = stage == self.iters;
        let pending = if complete { Vec::new() } else {
            self.work.acquisition_batch_attempts = self.work.acquisition_batch_attempts.checked_add(1)
                .ok_or(NoisyStudyError::WorkOverflow)?;
            let points = self.config.select_batch(&gp, &next.x, self.dim, stage, keep_going)
                .ok_or(NoisyStudyError::Cancelled)?;
            requests(self.study_id, next.x.len(), points)
        };
        if !keep_going() { return Err(NoisyStudyError::Cancelled); }
        self.report = next;
        self.pending = pending;
        self.complete = complete;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
