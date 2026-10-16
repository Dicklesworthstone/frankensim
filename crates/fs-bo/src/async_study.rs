//! Asynchronous fixed-kernel noisy BO for caller-owned simulation workers.
//!
//! Unlike `NoisyStudy`, a refill does not wait for a complete batch. It fits
//! received results in logical job order and treats every unanswered job as a
//! pending latent design in the acquisition. The caller owns actual execution,
//! in-flight deduplication and failure retries. No threads or I/O live here.

use crate::gp::Gp;
use crate::gp::prepared::PreparedError;
use crate::cached::CachedProposalError;
use crate::hyper::zero_noise_duplicates;
use crate::noisy::{NoisyBoConfig, NoisyObservation};
use crate::pending::PendingProposalError;
use crate::study::{NoisyRequest, NoisyStudyError, NoisyStudyWork, PendingEvaluation, TellDisposition};

/// Policy/advance failures preserve all received results and issued requests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsyncStudyError {
    /// Zero/overflowing budget or concurrency, or initial count above budget.
    InvalidPolicy,
    /// No received observation is available to condition a replacement on.
    WaitingForObservations,
    /// The received data do not admit a finite fixed-noise posterior.
    InvalidModel,
    /// The caller paused a numerical transition before publication.
    Cancelled,
    /// The bounded search did not return distinct unreserved designs.
    NoNovelCandidate,
    /// Explicit fixed-jitter cached acquisition refused its posterior.
    CachedPosterior(PreparedError),
    /// Attempt counters cannot be represented.
    WorkOverflow,
}

impl std::fmt::Display for AsyncStudyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "asynchronous noisy study: {self:?}")
    }
}
impl std::error::Error for AsyncStudyError {}

impl From<PendingProposalError> for AsyncStudyError {
    fn from(error: PendingProposalError) -> Self {
        match error {
            PendingProposalError::InvalidRequest => Self::InvalidPolicy,
            PendingProposalError::Cancelled => Self::Cancelled,
            PendingProposalError::NoNovelCandidate => Self::NoNovelCandidate,
        }
    }
}

impl From<CachedProposalError> for AsyncStudyError {
    fn from(error: CachedProposalError) -> Self {
        match error {
            CachedProposalError::InvalidRequest => Self::InvalidPolicy,
            CachedProposalError::Cancelled => Self::Cancelled,
            CachedProposalError::Posterior(e) => Self::CachedPosterior(e),
            CachedProposalError::NoNovelCandidate => Self::NoNovelCandidate,
        }
    }
}

/// Last successfully published model recommendation, with its exact data set.
/// Newly received results remain in `jobs()` until another successful refill
/// updates this record. No pending value is represented as an observation.
#[derive(Debug, Clone, PartialEq)]
pub struct AsyncModelReport {
    /// Indices into `jobs()`, in stable issuance order (not arrival order).
    pub observed_evaluations: Vec<usize>,
    /// Logical job index of the lowest posterior-mean RECEIVED design.
    pub recommended_evaluation: usize,
    /// Latent posterior mean in original objective units.
    pub mean: f64,
    /// Latent posterior variance, excluding measurement noise.
    pub variance: f64,
}

/// Sliding-window evaluation campaign with a fixed total job budget.
///
/// Clone retains every received result and unanswered job. Different result
/// sets at refill time may intentionally produce different designs. Replay
/// requires the same observation data AND decision-time availability sets;
/// merely permuting arrivals without a refill between them has no effect.
/// This API has no durable checkpoint wire format; batch-study checkpoints
/// and their semantics are unchanged.
#[derive(Debug, Clone)]
pub struct AsyncNoisyStudy {
    study_id: u64,
    config: NoisyBoConfig,
    initial: Vec<Vec<f64>>,
    max_evaluations: usize,
    max_in_flight: usize,
    jobs: Vec<PendingEvaluation>,
    model_report: Option<AsyncModelReport>,
    work: NoisyStudyWork,
}

fn jobs(study_id: u64, first: usize, points: &[Vec<f64>]) -> Vec<PendingEvaluation> {
    points.iter().enumerate().map(|(offset, x)| PendingEvaluation {
        request: NoisyRequest { study_id, evaluation: first + offset, x: x.clone() },
        observation: None,
    }).collect()
}

impl AsyncNoisyStudy {
    /// Prepare at most `max_in_flight` initial jobs without numerical solves.
    /// More Sobol initial jobs are issued as slots free up. After the entire
    /// initial design has been ISSUED, adaptive jobs may start after just one
    /// received result, even if some initial jobs are still running.
    ///
    /// `max_evaluations` includes the initial design and every reserved job,
    /// not just completed results. `config.q` caps adaptive jobs per refill;
    /// `max_in_flight` independently caps all unanswered issued jobs.
    ///
    /// # Panics
    /// The underlying BO configuration uses its existing admission assertions.
    /// Budget/concurrency errors return before any job generation.
    pub fn new(study_id: u64, dim: usize, n_init: usize, max_evaluations: usize,
        max_in_flight: usize, config: &NoisyBoConfig) -> Result<Self, AsyncStudyError>
    {
        if n_init == 0 || max_evaluations < n_init || max_in_flight == 0
            || u32::try_from(max_evaluations).is_err()
            || max_evaluations.checked_mul(max_evaluations).is_none()
            || max_evaluations.checked_mul(dim).is_none()
            || max_evaluations.checked_mul(config.mc_samples).is_none()
        {
            return Err(AsyncStudyError::InvalidPolicy);
        }
        config.validate_study(dim, max_evaluations, 0);
        let initial = config.initial_design(dim, n_init);
        let issued = jobs(study_id, 0, &initial[..n_init.min(max_in_flight)]);
        Ok(Self { study_id, config: config.clone(), initial, max_evaluations,
            max_in_flight, jobs: issued, model_report: None, work: NoisyStudyWork::default() })
    }

    /// Read unanswered jobs; repeated calls do no computation or dispatch.
    /// Maintain an external in-flight set so polling does not launch duplicates.
    #[must_use]
    pub fn ask(&self) -> Vec<NoisyRequest> {
        self.jobs.iter().filter(|j| j.observation.is_none()).map(|j| j.request.clone()).collect()
    }

    /// All issued requests and retained original observations, in logical order.
    #[must_use]
    pub fn jobs(&self) -> &[PendingEvaluation] { &self.jobs }

    /// Last successfully modeled availability set and its recommendation.
    #[must_use]
    pub fn model_report(&self) -> Option<&AsyncModelReport> { self.model_report.as_ref() }

    /// Count of unanswered issued jobs, including jobs awaiting external retry.
    #[must_use]
    pub fn in_flight(&self) -> usize {
        self.jobs.iter().filter(|j| j.observation.is_none()).count()
    }

    /// Counts include numerical attempts discarded by cancellation/refusal.
    #[must_use]
    pub fn work(&self) -> NoisyStudyWork { self.work }

    /// Completion requires all budgeted results AND their final model update.
    /// Call `refill` after final results even when no new jobs can be issued.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.jobs.len() == self.max_evaluations && self.in_flight() == 0
            && self.model_report.as_ref().is_some_and(|r|
                r.observed_evaluations.len() == self.max_evaluations)
    }

    /// Record a result without fitting. Requests must retain exact coordinates.
    /// Identical retries are idempotent, including after finalization; invalid
    /// or conflicting results never consume a slot or overwrite a response.
    /// No numerical value is substituted for a failed external simulation.
    pub fn tell(&mut self, request: &NoisyRequest, observation: NoisyObservation)
        -> Result<TellDisposition, NoisyStudyError>
    {
        if request.study_id != self.study_id { return Err(NoisyStudyError::ForeignStudy); }
        let job = self.jobs.get_mut(request.evaluation).ok_or(NoisyStudyError::UnknownEvaluation)?;
        if job.request.x.len() != request.x.len() || job.request.x.iter().zip(&request.x)
            .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            return Err(NoisyStudyError::RequestMismatch);
        }
        if !observation.value.is_finite()
            || !(observation.value - self.config.prior_mean).is_finite()
            || !observation.noise_variance.is_finite() || observation.noise_variance < 0.0
        {
            return Err(NoisyStudyError::InvalidObservation);
        }
        if let Some(previous) = job.observation {
            return if previous.value.to_bits() == observation.value.to_bits()
                && previous.noise_variance.to_bits() == observation.noise_variance.to_bits()
            { Ok(TellDisposition::Replayed) } else { Err(NoisyStudyError::ConflictingObservation) };
        }
        job.observation = Some(observation);
        Ok(TellDisposition::Recorded)
    }

    /// Fill available slots without waiting on stragglers. Returns the number
    /// of NEW jobs issued; use `ask` to discover their stable request IDs.
    pub fn refill(&mut self) -> Result<usize, AsyncStudyError> {
        self.refill_controlled(&mut || true)
    }

    /// Opt-in cached acquisition with a fixed posterior-jitter variance.
    /// The baseline factor and latent minima are computed once per refill,
    /// not for every CMA-ES probe. Existing `refill` behavior is unchanged.
    /// Retain this explicit jitter and choice of API when replaying decisions.
    pub fn refill_cached(&mut self, diagonal_jitter: f64) -> Result<usize, AsyncStudyError> {
        self.refill_cached_controlled(diagonal_jitter, &mut || true)
    }

    /// Cancellable cached refill. Jitter is validated before any job issuance
    /// or fit; cached numerical failures retain observations, model and jobs.
    /// Sample-row/extension checks improve cancellation granularity, but the
    /// current bounded CMA-ES call may finish with inert callbacks on refusal.
    pub fn refill_cached_controlled(&mut self, diagonal_jitter: f64,
        keep_going: &mut dyn FnMut() -> bool) -> Result<usize, AsyncStudyError>
    {
        if !diagonal_jitter.is_finite() || diagonal_jitter < 0.0 {
            return Err(AsyncStudyError::InvalidPolicy);
        }
        self.refill_impl(keep_going, Some(diagonal_jitter))
    }

    /// Atomically publish a modeled availability set and its replacement jobs.
    /// Pending designs enter the joint acquisition, never the training values.
    /// Once total issuance reaches the budget, later refills only update the
    /// model as remaining responses arrive. An unchanged full window is a no-op.
    ///
    /// Cancellation or numerical refusal retains all observations, previous
    /// model report and reservations. Retry from the same data uses the same
    /// selection seed, keyed by the number of adaptive jobs already issued.
    /// One dense fit or one greedy slot remains non-preemptible. Successful
    /// issue counts respect both the total budget and the in-flight cap.
    ///
    /// # Panics
    /// Existing GP/acquisition arithmetic assertions remain. All such work
    /// precedes publication of jobs or a new model report.
    pub fn refill_controlled(&mut self, keep_going: &mut dyn FnMut() -> bool)
        -> Result<usize, AsyncStudyError>
    {
        self.refill_impl(keep_going, None)
    }

    fn refill_impl(&mut self, keep_going: &mut dyn FnMut() -> bool, cached_jitter: Option<f64>)
        -> Result<usize, AsyncStudyError>
    {
        if !keep_going() { return Err(AsyncStudyError::Cancelled); }
        let issued = self.jobs.len();
        let free = (self.max_in_flight - self.in_flight()).min(self.max_evaluations - issued);
        // Complete the declared initial design first, but do not require its
        // outstanding results before adaptive selection once it is issued.
        if issued < self.initial.len() {
            let count = free.min(self.initial.len() - issued);
            if count == 0 { return Ok(0); }
            let new = jobs(self.study_id, issued, &self.initial[issued..issued + count]);
            if !keep_going() { return Err(AsyncStudyError::Cancelled); }
            self.jobs.extend(new);
            return Ok(count);
        }
        let observed_evaluations: Vec<usize> = self.jobs.iter().enumerate()
            .filter_map(|(i, j)| j.observation.map(|_| i)).collect();
        let count = free.min(self.config.q);
        if count == 0 && self.model_report.as_ref().is_some_and(|r|
            r.observed_evaluations == observed_evaluations) { return Ok(0); }
        if observed_evaluations.is_empty() {
            return if count == 0 { Ok(0) } else { Err(AsyncStudyError::WaitingForObservations) };
        }
        let mut x = Vec::with_capacity(observed_evaluations.len());
        let mut values = Vec::with_capacity(observed_evaluations.len());
        let mut variances = Vec::with_capacity(observed_evaluations.len());
        for &i in &observed_evaluations {
            let job = &self.jobs[i];
            let o = job.observation.expect("received result");
            x.push(job.request.x.clone());
            values.push(o.value - self.config.prior_mean);
            variances.push(o.noise_variance);
        }
        self.work.model_fit_attempts = self.work.model_fit_attempts.checked_add(1)
            .ok_or(AsyncStudyError::WorkOverflow)?;
        if zero_noise_duplicates(&x, &variances) { return Err(AsyncStudyError::InvalidModel); }
        let gp = Gp::try_fit_diag(&x, &values, self.config.kernel.clone(), &variances)
            .filter(|g| g.lml.is_finite()).ok_or(AsyncStudyError::InvalidModel)?;
        if !keep_going() { return Err(AsyncStudyError::Cancelled); }
        let best = self.config.study_incumbent(&gp, &x);
        let report = AsyncModelReport {
            recommended_evaluation: observed_evaluations[best.observation_index],
            observed_evaluations, mean: best.mean, variance: best.variance,
        };
        let new = if count == 0 { Vec::new() } else {
            let pending: Vec<Vec<f64>> = self.jobs.iter().filter(|j| j.observation.is_none())
                .map(|j| j.request.x.clone()).collect();
            self.work.acquisition_batch_attempts = self.work.acquisition_batch_attempts.checked_add(1)
                .ok_or(AsyncStudyError::WorkOverflow)?;
            let points = if let Some(jitter) = cached_jitter {
                self.config.propose_pending_cached(&gp, &x, &pending, count,
                    issued - self.initial.len(), jitter, keep_going)?
            } else {
                self.config.propose_pending(&gp, &x, &pending, count,
                    issued - self.initial.len(), keep_going)?
            };
            jobs(self.study_id, issued, &points)
        };
        if !keep_going() { return Err(AsyncStudyError::Cancelled); }
        self.jobs.reserve(new.len()); // allocation, if needed, precedes publication
        self.jobs.extend(new);
        self.model_report = Some(report);
        Ok(count)
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod cached_tests;
