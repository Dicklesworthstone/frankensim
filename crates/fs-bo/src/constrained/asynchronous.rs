//! External multi-output simulations with no complete-batch barrier.
//! Reuses the constrained model, feasibility, and acquisition implementations.
use super::bo::{condition_controlled, recommend};
use super::{ConstrainedBoConfig, ConstrainedBoError, ConstrainedBoReport,
    ConstrainedObservation, ConstrainedProposalError, ConstrainedRecommendation};
use crate::study::{NoisyRequest, NoisyStudyError, NoisyStudyWork, TellDisposition};
use crate::gp::Kernel;

mod checkpoint;

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
    /// The fixed-kernel constructor refuses a learning request. Select
    /// new_with_learning to explicitly use model-update-based refit cadence.
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
/// Each GP retains its declared prior and per-observation variance. Kernels
/// stay fixed by default; new_with_learning opts into bounded per-output fits.
/// Candidate evaluations may violate constraints: the recommendation threshold
/// is not a safe-exploration guarantee. Outputs are modeled independently.
///
/// Clone is an in-memory checkpoint. Replay requires identical observations
/// and availability sets at refill time, not just the optimizer seed. This
/// owns no filesystem, worker threads, or job-lease protocol. checkpoint_bytes
/// and from_checkpoint also retain this state across process restarts.
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
    kernels: Vec<Kernel>,
    model_updates: usize,
    likelihood_evaluations: usize,
    posterior_only_fits: usize,
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
        Self::initialize(study_id, dim, n_init, max_evaluations, max_in_flight, config, false)
    }

    /// Opt into per-output kernel learning without a complete-batch barrier.
    ///
    /// Each configured refit_every counts successful model updates on NEW
    /// received-job availability sets: first update is zero, then one, two,
    /// etc. It does not count arrivals, issued batches or refused attempts.
    /// Repeated refills on the same data may reserve more jobs, but only
    /// recondition the accepted kernels; they do not retrain or advance cadence.
    ///
    /// Each output warm-starts its own last accepted kernel and keeps its
    /// declared mean and original noise variances. Kernels and the new report
    /// commit only together with a complete proposal. Clone retains this state.
    /// No pending result is used as a training target, even for refitting.
    pub fn new_with_learning(study_id: u64, dim: usize, n_init: usize,
        max_evaluations: usize, max_in_flight: usize, config: &ConstrainedBoConfig)
        -> Result<Self, AsyncConstrainedError>
    {
        Self::initialize(study_id, dim, n_init, max_evaluations, max_in_flight, config, true)
    }

    fn initialize(study_id: u64, dim: usize, n_init: usize, max_evaluations: usize,
        max_in_flight: usize, config: &ConstrainedBoConfig, allow_learning: bool)
        -> Result<Self, AsyncConstrainedError>
    {
        if n_init == 0 || max_evaluations < n_init || max_in_flight == 0
            || u32::try_from(max_evaluations).is_err()
        { return Err(AsyncConstrainedError::InvalidPolicy); }
        config.validate(dim, max_evaluations, 0)?;
        if !allow_learning && (config.objective_learning.is_some()
            || config.constraints.iter().any(|c| c.learning.is_some())) {
            return Err(AsyncConstrainedError::LearningUnsupported);
        }
        // There can be at most one distinct nonempty received set per result.
        // Admit the worst-case total learning work without pretending q jobs
        // arrive together or multiplying the physical-job budget by q.
        let mut fit_work = 0usize;
        for learning in std::iter::once(config.objective_learning.as_ref())
            .chain(config.constraints.iter().map(|c| c.learning.as_ref())).flatten()
        {
            let stages = 1 + (max_evaluations - 1) / learning.refit_every;
            fit_work = stages.checked_mul(learning.fit.max_evaluations)
                .and_then(|n| fit_work.checked_add(n))
                .ok_or(AsyncConstrainedError::InvalidPolicy)?;
        }
        let initial = config.search.initial_design(dim, n_init);
        let issued = jobs(study_id, 0, &initial[..n_init.min(max_in_flight)]);
        Ok(Self { study_id, config: config.clone(), initial, max_evaluations,
            max_in_flight, jobs: issued, model_report: None, work: NoisyStudyWork::default(),
            kernels: config.initial_kernels(), model_updates: 0,
            likelihood_evaluations: 0, posterior_only_fits: 0 })
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
    /// Last accepted objective kernel followed by constraint kernels. A failed
    /// transition cannot expose partially refitted models through this view.
    #[must_use]
    pub fn kernels(&self) -> &[Kernel] { &self.kernels }
    /// Number of accepted, distinct received-data revisions. The next new set
    /// uses this number as its zero-based learning-stage ordinal.
    #[must_use]
    pub fn model_updates(&self) -> usize { self.model_updates }
    /// Likelihood probes from COMMITTED learning stages only; discarded or
    /// cancelled passes remain visible through work().model_fit_attempts.
    #[must_use]
    pub fn likelihood_evaluations(&self) -> usize { self.likelihood_evaluations }
    /// Ordinary conditioning fits from committed passes, including fixed
    /// outputs and extra reservations made without receiving any new data.
    #[must_use]
    pub fn posterior_only_fits(&self) -> usize { self.posterior_only_fits }
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
    /// Continuation checks occur around and between output fits, inside learned
    /// gradient/search work and between greedy acquisition slots. A dense
    /// factorization or a slot search remains non-preemptible.
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
        let new_data = self.model_report.as_ref().is_none_or(|r| r.observed_evaluations != received);
        if count == 0 && !new_data {
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
        let mut kernels = self.kernels.clone();
        let mut conditioning = self.config.clone();
        if !new_data {
            conditioning.objective_learning = None;
            for output in &mut conditioning.constraints { output.learning = None; }
        }
        let models = condition_controlled(&mut trial, &conditioning, &mut kernels,
            self.model_updates, keep_going)?
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
        let model_updates = self.model_updates.checked_add(usize::from(new_data))
            .ok_or(AsyncConstrainedError::WorkOverflow)?;
        let likelihood_evaluations = self.likelihood_evaluations.checked_add(trial.likelihood_evaluations)
            .ok_or(AsyncConstrainedError::WorkOverflow)?;
        let posterior_only_fits = self.posterior_only_fits.checked_add(trial.posterior_only_fits)
            .ok_or(AsyncConstrainedError::WorkOverflow)?;
        self.jobs.reserve(count); // allocation, if needed, precedes publication
        self.jobs.extend(next);
        self.model_report = Some(report);
        self.kernels = kernels;
        self.model_updates = model_updates;
        self.likelihood_evaluations = likelihood_evaluations;
        self.posterior_only_fits = posterior_only_fits;
        Ok(count)
    }
}

#[cfg(test)]
mod tests;
