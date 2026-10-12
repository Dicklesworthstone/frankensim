//! External multi-output simulations with no complete-batch barrier.
//! Reuses the constrained model, feasibility, and acquisition implementations.
use super::bo::{condition_controlled, recommend};
use super::{ConstrainedBoConfig, ConstrainedBoError, ConstrainedBoReport,
    ConstrainedObservation, ConstrainedProposalError, ConstrainedRecommendation};
use crate::study::{NoisyRequest, NoisyStudyError, NoisyStudyWork, TellDisposition};

/// A simulation's complete objective and constraints, retained together.
#[derive(Debug, Clone, PartialEq)]
pub struct PendingConstrainedEvaluation {
    /// Stable job identity and exact input coordinates.
    pub request: NoisyRequest,
    /// None until ALL declared outputs arrive in one valid response.
    pub observation: Option<ConstrainedObservation>,
}

/// The received-data set behind the last successful model update.
#[derive(Debug, Clone, PartialEq)]
pub struct AsyncConstrainedModelReport {
    /// Logical job IDs in model row order; pending jobs are never training rows.
    pub observed_evaluations: Vec<usize>,
    /// The observation_index indexes observed_evaluations, NOT all issued jobs.
    /// None means no received design meets the model probability policy. It is
    /// not proof of physical infeasibility, and does not stop further sampling.
    pub recommendation: Option<ConstrainedRecommendation>,
}

/// Refusals retain all results, reservations, and the last modeled data set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AsyncConstrainedError {
    /// Invalid total-job budget or in-flight allowance.
    InvalidPolicy,
    /// Async learning cadence is not declared here; a learning request refuses
    /// rather than silently freezing or rescheduling any output's kernel.
    LearningUnsupported,
    /// Request identity or duplicate-delivery conflict.
    Delivery(NoisyStudyError),
    /// Invalid config, output count/value, or a refused output model/posterior.
    Model(ConstrainedBoError),
    /// At least one received simulation is needed before adaptive selection.
    WaitingForObservations,
    /// A paused transition publishes no partial models or reservations.
    Cancelled,
    /// Bounded selection returned an already observed/reserved design.
    NoNovelCandidate,
    /// Attempt accounting cannot be represented.
    WorkOverflow,
}
impl std::fmt::Display for AsyncConstrainedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "asynchronous constrained study: {self:?}")
    }
}
impl std::error::Error for AsyncConstrainedError {}
impl From<ConstrainedBoError> for AsyncConstrainedError {
    fn from(error: ConstrainedBoError) -> Self { Self::Model(error) }
}
impl From<ConstrainedProposalError> for AsyncConstrainedError {
    fn from(error: ConstrainedProposalError) -> Self {
        match error {
            ConstrainedProposalError::InvalidRequest => Self::InvalidPolicy,
            ConstrainedProposalError::Cancelled => Self::Cancelled,
            ConstrainedProposalError::NoNovelCandidate => Self::NoNovelCandidate,
        }
    }
}

/// Pending-aware constrained design search over caller-owned simulation jobs.
///
/// Each GP retains a fixed declared kernel/prior and per-observation variance.
/// Candidate evaluations may violate constraints: the recommendation threshold
/// is not a safe-exploration guarantee. Outputs are modeled independently.
///
/// Clone is an in-memory checkpoint. Replay requires identical observations
/// and availability sets at refill time, not just the optimizer seed. This
/// adds no filesystem, worker threads, durable format, or job-lease protocol.
#[derive(Debug, Clone)]
pub struct AsyncConstrainedStudy {
    study_id: u64,
    config: ConstrainedBoConfig,
    initial: Vec<Vec<f64>>,
    max_evaluations: usize,
    max_in_flight: usize,
    jobs: Vec<PendingConstrainedEvaluation>,
    model_report: Option<AsyncConstrainedModelReport>,
    work: NoisyStudyWork,
}

fn jobs(study_id: u64, first: usize, points: &[Vec<f64>]) -> Vec<PendingConstrainedEvaluation> {
    points.iter().enumerate().map(|(slot, x)| PendingConstrainedEvaluation {
        request: NoisyRequest { study_id, evaluation: first + slot, x: x.clone() }, observation: None,
    }).collect()
}
fn identical(a: &ConstrainedObservation, b: &ConstrainedObservation) -> bool {
    a.constraints.len() == b.constraints.len()
        && std::iter::once(&a.objective).chain(&a.constraints)
            .zip(std::iter::once(&b.objective).chain(&b.constraints))
            .all(|(x, y)| x.value.to_bits() == y.value.to_bits()
                && x.noise_variance.to_bits() == y.noise_variance.to_bits())
}

impl AsyncConstrainedStudy {
    /// Admit a fixed-kernel campaign and issue the first bounded initial wave.
    ///
    /// The total budget includes initial jobs and every reserved replacement.
    /// max_in_flight counts unanswered issued jobs, including jobs awaiting an
    /// external retry. search.q separately caps adaptive jobs per refill.
    /// Invalid policies (including any requested kernel learning) refuse before
    /// any job generation. No simulation or model fit runs in this constructor.
    pub fn new(study_id: u64, dim: usize, n_init: usize, max_evaluations: usize,
        max_in_flight: usize, config: &ConstrainedBoConfig) -> Result<Self, AsyncConstrainedError>
    {
        if n_init == 0 || max_evaluations < n_init || max_in_flight == 0
            || u32::try_from(max_evaluations).is_err()
        { return Err(AsyncConstrainedError::InvalidPolicy); }
        config.validate(dim, max_evaluations, 0)?;
        if config.objective_learning.is_some() || config.constraints.iter().any(|c| c.learning.is_some()) {
            return Err(AsyncConstrainedError::LearningUnsupported);
        }
        let initial = config.search.initial_design(dim, n_init);
        let issued = jobs(study_id, 0, &initial[..n_init.min(max_in_flight)]);
        Ok(Self { study_id, config: config.clone(), initial, max_evaluations,
            max_in_flight, jobs: issued, model_report: None, work: NoisyStudyWork::default() })
    }

    /// Read unanswered jobs. This is not a dispatch/lease operation: callers
    /// must deduplicate outstanding IDs before launching physical simulations.
    #[must_use]
    pub fn ask(&self) -> Vec<NoisyRequest> {
        self.jobs.iter().filter(|j| j.observation.is_none()).map(|j| j.request.clone()).collect()
    }
    /// Every issued request and unmodified complete response, in issuance order.
    #[must_use]
    pub fn jobs(&self) -> &[PendingConstrainedEvaluation] { &self.jobs }
    /// Last fully published model update; later received data remain in jobs().
    #[must_use]
    pub fn model_report(&self) -> Option<&AsyncConstrainedModelReport> { self.model_report.as_ref() }
    /// All unanswered reservations, whether dispatched yet or not.
    #[must_use]
    pub fn in_flight(&self) -> usize { self.jobs.iter().filter(|j| j.observation.is_none()).count() }
    /// A model_fit_attempt counts one MULTI-OUTPUT conditioning pass, including
    /// refused/cancelled passes; acquisition_batch_attempts includes discarded
    /// proposals. These are neither covariance-factorization nor physics counts.
    #[must_use]
    pub fn work(&self) -> NoisyStudyWork { self.work }
    /// All reserved results must receive a final model update. A valid final
    /// model can still have no recommended feasible design.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.jobs.len() == self.max_evaluations && self.in_flight() == 0
            && self.model_report.as_ref().is_some_and(|r|
                r.observed_evaluations.len() == self.max_evaluations)
    }

    /// Accept ALL outcomes in declaration order without doing numerical work.
    /// Incorrect constraints, variances, IDs or coordinates never release a
    /// slot. Identical retries are idempotent; changed objective OR constraint
    /// values refuse even after finalization. Never use fake outcomes for a
    /// failed/timeout simulation: it remains an unanswered request for retry.
    pub fn tell(&mut self, request: &NoisyRequest, observation: ConstrainedObservation)
        -> Result<TellDisposition, AsyncConstrainedError>
    {
        let delivery = AsyncConstrainedError::Delivery;
        if request.study_id != self.study_id { return Err(delivery(NoisyStudyError::ForeignStudy)); }
        let job = self.jobs.get(request.evaluation).ok_or(delivery(NoisyStudyError::UnknownEvaluation))?;
        if job.request.x.len() != request.x.len() || job.request.x.iter().zip(&request.x)
            .any(|(a, b)| a.to_bits() != b.to_bits())
        { return Err(delivery(NoisyStudyError::RequestMismatch)); }
        self.config.validate_observation(&observation, request.evaluation)?;
        if let Some(previous) = &job.observation {
            return if identical(previous, &observation) { Ok(TellDisposition::Replayed) }
                else { Err(delivery(NoisyStudyError::ConflictingObservation)) };
        }
        self.jobs[request.evaluation].observation = Some(observation);
        Ok(TellDisposition::Recorded)
    }

    /// Refit received outcomes and issue replacements without waiting for all
    /// outstanding jobs. Returns the number of newly reserved evaluations.
    pub fn refill(&mut self) -> Result<usize, AsyncConstrainedError> {
        self.refill_controlled(&mut || true)
    }

    /// Publish all output models/recommendations and all new jobs together.
    ///
    /// Initial Sobol jobs are released in order as slots free up. Once all have
    /// been issued, adaptation requires one received complete output vector,
    /// not a complete batch. All unanswered designs enter the joint feasible
    /// acquisition. Final-budget refills still incorporate late results.
    ///
    /// Continuation checks occur around and between output fits and between
    /// greedy acquisition slots. A single fit or slot search is non-preemptible.
    /// Cancellation/failure retains all received results and the previous model
    /// report; attempts remain counted. Selection ordinals depend on the number
    /// of already issued adaptive jobs, never failed proposal attempts.
    ///
    /// # Panics
    /// Inherits acquisition arithmetic assertions, before state publication.
    pub fn refill_controlled(&mut self, keep_going: &mut dyn FnMut() -> bool)
        -> Result<usize, AsyncConstrainedError>
    {
        if !keep_going() { return Err(AsyncConstrainedError::Cancelled); }
        let issued = self.jobs.len();
        let free = (self.max_in_flight - self.in_flight()).min(self.max_evaluations - issued);
        if issued < self.initial.len() {
            let count = free.min(self.initial.len() - issued);
            if count == 0 { return Ok(0); }
            let next = jobs(self.study_id, issued, &self.initial[issued..issued + count]);
            if !keep_going() { return Err(AsyncConstrainedError::Cancelled); }
            self.jobs.extend(next);
            return Ok(count);
        }
        let received: Vec<usize> = self.jobs.iter().enumerate()
            .filter_map(|(i, j)| j.observation.as_ref().map(|_| i)).collect();
        let count = free.min(self.config.search.q);
        if count == 0 && self.model_report.as_ref().is_some_and(|r| r.observed_evaluations == received) {
            return Ok(0);
        }
        if received.is_empty() {
            return if count == 0 { Ok(0) } else { Err(AsyncConstrainedError::WaitingForObservations) };
        }
        // A trial data set, never an asynchronous completion-order history.
        let mut trial = ConstrainedBoReport {
            x: received.iter().map(|&i| self.jobs[i].request.x.clone()).collect(),
            observations: received.iter().map(|&i| self.jobs[i].observation.clone().expect("received result")).collect(),
            recommendations: Vec::new(), fits: Vec::new(), likelihood_evaluations: 0, posterior_only_fits: 0,
        };
        self.work.model_fit_attempts = self.work.model_fit_attempts.checked_add(1)
            .ok_or(AsyncConstrainedError::WorkOverflow)?;
        let mut kernels = self.config.initial_kernels();
        let models = condition_controlled(&mut trial, &self.config, &mut kernels, 0, keep_going)?
            .ok_or(AsyncConstrainedError::Cancelled)?;
        let report = AsyncConstrainedModelReport { observed_evaluations: received,
            recommendation: recommend(&models, &trial.x, &self.config)? };
        let next = if count == 0 { Vec::new() } else {
            let pending: Vec<_> = self.jobs.iter().filter(|j| j.observation.is_none())
                .map(|j| j.request.x.clone()).collect();
            self.work.acquisition_batch_attempts = self.work.acquisition_batch_attempts.checked_add(1)
                .ok_or(AsyncConstrainedError::WorkOverflow)?;
            let points = self.config.propose_pending(&models, &trial.x, &pending, count,
                issued - self.initial.len(), keep_going)?;
            jobs(self.study_id, issued, &points)
        };
        if !keep_going() { return Err(AsyncConstrainedError::Cancelled); }
        let count = next.len();
        self.jobs.reserve(count); // allocation, if needed, precedes publication
        self.jobs.extend(next);
        self.model_report = Some(report);
        Ok(count)
    }
}

#[cfg(test)]
mod tests;
