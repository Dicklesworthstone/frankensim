//! Exact-bit durable state for constrained asynchronous workers. No numerical
//! work runs while restoring: the last modeled availability set is retained
//! separately from more recently received observations. Trusted storage only.

use super::*;
use crate::constrained::OutcomeConstraint;
use crate::gp::Matern;
use crate::hyper::HeteroFitConfig;
use crate::learning::NoisyLearningConfig;
use crate::noisy::{NoisyBoConfig, NoisyObservation};
use crate::study::{CheckpointError, CheckpointLimits};

const MAGIC: &[u8; 8] = b"FSACBO\x01\0";
type Result<T> = std::result::Result<T, CheckpointError>;
fn require(ok: bool) -> Result<()> {
    if ok { Ok(()) } else { Err(CheckpointError::InvalidState) }
}
fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325u64,
        |h, b| (h ^ u64::from(*b)).wrapping_mul(0x100_0000_01b3))
}

struct Writer(Vec<u8>);
impl Writer {
    fn word(&mut self, v: u64) { self.0.extend_from_slice(&v.to_le_bytes()); }
    fn size(&mut self, v: usize) { self.word(u64::try_from(v).expect("64-bit checkpoint size")); }
    fn float(&mut self, v: f64) { self.word(v.to_bits()); }
    fn floats(&mut self, vs: &[f64]) { for v in vs { self.float(*v); } }
    fn kernel(&mut self, k: &Kernel) {
        self.word(match k.family { Matern::Half => 0, Matern::ThreeHalves => 1, Matern::FiveHalves => 2 });
        self.float(k.signal); self.floats(&k.lengthscales);
    }
    fn learning(&mut self, p: Option<&NoisyLearningConfig>) {
        self.word(u64::from(p.is_some()));
        if let Some(p) = p {
            self.size(p.refit_every);
            for (lo, hi) in &p.fit.lengthscale_bounds { self.float(*lo); self.float(*hi); }
            self.float(p.fit.signal_bounds.0); self.float(p.fit.signal_bounds.1);
            self.size(p.fit.starts); self.size(p.fit.max_iterations); self.size(p.fit.max_evaluations);
            self.float(p.fit.gradient_tolerance); self.word(p.fit.seed);
        }
    }
    fn observation(&mut self, o: &NoisyObservation) {
        self.float(o.value); self.float(o.noise_variance);
    }
}

struct Reader<'a>(&'a [u8]);
impl Reader<'_> {
    fn word(&mut self) -> Result<u64> {
        let b = self.0.get(..8).ok_or(CheckpointError::Truncated)?;
        let v = u64::from_le_bytes(b.try_into().map_err(|_| CheckpointError::Truncated)?);
        self.0 = &self.0[8..]; Ok(v)
    }
    fn size(&mut self) -> Result<usize> {
        usize::try_from(self.word()?).map_err(|_| CheckpointError::LimitExceeded)
    }
    fn float(&mut self) -> Result<f64> {
        let v = f64::from_bits(self.word()?); require(v.is_finite())?; Ok(v)
    }
    fn flag(&mut self) -> Result<bool> {
        match self.word()? { 0 => Ok(false), 1 => Ok(true), _ => Err(CheckpointError::InvalidState) }
    }
    // The byte frame, not only caller-supplied counts, bounds EVERY allocation.
    fn space(&self, count: usize, words_each: usize) -> Result<()> {
        if count > self.0.len() / 8 / words_each { Err(CheckpointError::Truncated) } else { Ok(()) }
    }
    fn count(&mut self, max: usize, words_each: usize) -> Result<usize> {
        let n = self.size()?;
        if n > max { return Err(CheckpointError::LimitExceeded); }
        self.space(n, words_each)?; Ok(n)
    }
    fn floats(&mut self, n: usize) -> Result<Vec<f64>> {
        self.space(n, 1)?; (0..n).map(|_| self.float()).collect()
    }
    fn kernel(&mut self, dim: usize) -> Result<Kernel> {
        let family = match self.word()? {
            0 => Matern::Half, 1 => Matern::ThreeHalves, 2 => Matern::FiveHalves,
            _ => return Err(CheckpointError::InvalidState),
        };
        let signal = self.float()?; let lengthscales = self.floats(dim)?;
        require(signal > 0.0 && lengthscales.iter().all(|v| *v > 0.0))?;
        Ok(Kernel { family, signal, lengthscales })
    }
    fn learning(&mut self, dim: usize) -> Result<Option<NoisyLearningConfig>> {
        if !self.flag()? { return Ok(None); }
        let refit_every = self.size()?;
        self.space(dim, 2)?;
        let lengthscale_bounds = (0..dim).map(|_| Ok((self.float()?, self.float()?)))
            .collect::<Result<Vec<_>>>()?;
        Ok(Some(NoisyLearningConfig { refit_every, fit: HeteroFitConfig {
            lengthscale_bounds, signal_bounds: (self.float()?, self.float()?),
            starts: self.size()?, max_iterations: self.size()?, max_evaluations: self.size()?,
            gradient_tolerance: self.float()?, seed: self.word()?,
        }}))
    }
    fn point(&mut self, dim: usize, bounds: (f64, f64)) -> Result<Vec<f64>> {
        let x = self.floats(dim)?;
        require(x.iter().all(|v| (bounds.0..=bounds.1).contains(v)))?; Ok(x)
    }
    fn observation(&mut self) -> Result<NoisyObservation> {
        let value = self.float()?; let noise_variance = self.float()?;
        require(noise_variance >= 0.0)?; Ok(NoisyObservation { value, noise_variance })
    }
}

impl AsyncConstrainedStudy {
    /// Encode fixed or learned state, including pending jobs, late responses,
    /// accepted kernels and refit clocks. This performs no fit or simulation.
    /// Persist atomically in caller-owned storage. The v1 checksum detects
    /// accidental corruption, NOT forgery. Physics identity and worker leases
    /// must be bound by the application; this is not a dispatch journal.
    #[must_use]
    pub fn checkpoint_bytes(&self) -> Vec<u8> {
        let mut w = Writer(MAGIC.to_vec());
        let c = &self.config; let s = &c.search;
        w.word(self.study_id); w.size(s.kernel.lengthscales.len());
        w.size(self.max_evaluations); w.size(self.max_in_flight);
        w.size(c.constraints.len() + 1);
        w.float(s.bounds.0); w.float(s.bounds.1); w.kernel(&s.kernel); w.float(s.prior_mean);
        for n in [s.q, s.mc_samples, s.acq_starts, s.acq_evals] { w.size(n); }
        w.word(s.seed); w.learning(c.objective_learning.as_ref());
        for output in &c.constraints {
            w.kernel(&output.kernel); w.float(output.prior_mean); w.float(output.upper_bound);
            w.learning(output.learning.as_ref());
        }
        w.float(c.reference); w.float(c.recommendation_probability);
        w.size(self.initial.len()); for x in &self.initial { w.floats(x); }
        w.size(self.jobs.len());
        for job in &self.jobs {
            w.floats(&job.request.x); w.word(u64::from(job.observation.is_some()));
            if let Some(o) = &job.observation {
                w.observation(&o.objective);
                for constraint in &o.constraints { w.observation(constraint); }
            }
        }
        w.word(u64::from(self.model_report.is_some()));
        if let Some(r) = &self.model_report {
            w.size(r.observed_evaluations.len()); for id in &r.observed_evaluations { w.size(*id); }
            w.word(u64::from(r.recommendation.is_some()));
            if let Some(best) = &r.recommendation {
                w.size(best.observation_index); w.float(best.objective_mean); w.float(best.objective_variance);
                w.floats(&best.constraint_means); w.floats(&best.constraint_variances);
                w.float(best.joint_feasibility);
            }
        }
        for k in &self.kernels { w.kernel(k); }
        for n in [self.model_updates, self.likelihood_evaluations, self.posterior_only_fits,
            self.work.model_fit_attempts, self.work.acquisition_batch_attempts] { w.size(n); }
        w.word(checksum(&w.0)); w.0
    }

    /// Restore trusted v1 state without generating jobs, fitting models or
    /// changing the last modeled availability set. `max_evaluations` caps the
    /// ENTIRE reservation budget; `max_bytes` bounds the frame and allocations,
    /// including output-model counts. Malformed input never yields a partial
    /// study. Numerical results are structurally checked, not re-verified.
    ///
    /// Replay requires the same numerical build/profile, external physics and
    /// subsequent availability sets. Already received results stay idempotent;
    /// unanswered IDs are unchanged. Reconcile them with the external dispatcher
    /// before retrying any job: restoring a file does not acquire worker leases.
    pub fn from_checkpoint(bytes: &[u8], limits: CheckpointLimits) -> Result<Self> {
        if bytes.len() > limits.max_bytes { return Err(CheckpointError::LimitExceeded); }
        if bytes.len() < 16 { return Err(CheckpointError::Truncated); }
        if &bytes[..8] != MAGIC { return Err(CheckpointError::IncompatibleVersion); }
        let split = bytes.len() - 8;
        if Reader(&bytes[split..]).word()? != checksum(&bytes[..split]) {
            return Err(CheckpointError::ChecksumMismatch);
        }
        let mut r = Reader(&bytes[8..split]);
        let study_id = r.word()?; let dim = r.size()?;
        require((1..=fs_rand::qmc::MAX_SOBOL_DIM).contains(&dim))?;
        let max_evaluations = r.size()?; let max_in_flight = r.size()?;
        if max_evaluations > limits.max_evaluations { return Err(CheckpointError::LimitExceeded); }
        require(max_evaluations > 0 && u32::try_from(max_evaluations).is_ok() && max_in_flight > 0)?;
        let outputs = r.count(usize::MAX, dim + 4)?; require(outputs > 0)?;
        let bounds = (r.float()?, r.float()?); let kernel = r.kernel(dim)?; let prior_mean = r.float()?;
        let search = NoisyBoConfig { bounds, kernel, prior_mean, q: r.size()?, mc_samples: r.size()?,
            acq_starts: r.size()?, acq_evals: r.size()?, seed: r.word()? };
        let objective_learning = r.learning(dim)?;
        let constraints = (1..outputs).map(|_| Ok(OutcomeConstraint {
            kernel: r.kernel(dim)?, prior_mean: r.float()?, upper_bound: r.float()?, learning: r.learning(dim)?,
        })).collect::<Result<Vec<_>>>()?;
        let config = ConstrainedBoConfig { search, objective_learning, constraints,
            reference: r.float()?, recommendation_probability: r.float()? };
        config.validate(dim, max_evaluations, 0).map_err(|_| CheckpointError::InvalidState)?;
        let n_initial = r.count(max_evaluations, dim)?; require(n_initial > 0)?;
        let initial = (0..n_initial).map(|_| r.point(dim, bounds)).collect::<Result<Vec<_>>>()?;
        let issued = r.count(max_evaluations, dim + 1)?;
        let mut jobs = Vec::with_capacity(issued);
        for evaluation in 0..issued {
            let x = r.point(dim, bounds)?;
            let observation = if r.flag()? {
                r.space(outputs, 2)?;
                let o = ConstrainedObservation { objective: r.observation()?,
                    constraints: (1..outputs).map(|_| r.observation()).collect::<Result<Vec<_>>>()? };
                config.validate_observation(&o, evaluation).map_err(|_| CheckpointError::InvalidState)?;
                Some(o)
            } else { None };
            jobs.push(PendingConstrainedEvaluation { request: NoisyRequest { study_id, evaluation, x }, observation });
        }
        let model_report = if r.flag()? {
            let n = r.count(issued, 1)?;
            let observed_evaluations = (0..n).map(|_| r.size()).collect::<Result<Vec<_>>>()?;
            let recommendation = if r.flag()? {
                Some(ConstrainedRecommendation { observation_index: r.size()?, objective_mean: r.float()?,
                    objective_variance: r.float()?, constraint_means: r.floats(outputs - 1)?,
                    constraint_variances: r.floats(outputs - 1)?, joint_feasibility: r.float()? })
            } else { None };
            Some(AsyncConstrainedModelReport { observed_evaluations, recommendation })
        } else { None };
        r.space(outputs, dim + 2)?;
        let kernels = (0..outputs).map(|_| r.kernel(dim)).collect::<Result<Vec<_>>>()?;
        let model_updates = r.size()?; let likelihood_evaluations = r.size()?; let posterior_only_fits = r.size()?;
        let work = NoisyStudyWork { model_fit_attempts: r.size()?, acquisition_batch_attempts: r.size()? };
        require(r.0.is_empty())?;
        let study = Self { study_id, config, initial, max_evaluations, max_in_flight, jobs,
            model_report, work, kernels, model_updates, likelihood_evaluations, posterior_only_fits };
        study.validate_restored()?;
        Ok(study)
    }

    fn validate_restored(&self) -> Result<()> {
        let issued = self.jobs.len(); let n_initial = self.initial.len();
        require(issued >= n_initial.min(self.max_in_flight) && self.in_flight() <= self.max_in_flight)?;
        for (job, original) in self.jobs.iter().zip(&self.initial) {
            require(job.request.x.iter().zip(original).all(|(a, b)| a.to_bits() == b.to_bits()))?;
        }
        require(self.work.acquisition_batch_attempts <= self.work.model_fit_attempts)?;
        let adaptive = issued.saturating_sub(n_initial);
        let q = self.config.search.q;
        let required_proposals = adaptive / q + usize::from(adaptive % q != 0);
        require(required_proposals <= self.work.acquisition_batch_attempts)?;
        let modeled = if let Some(r) = &self.model_report {
            let ids = &r.observed_evaluations;
            require(issued >= n_initial && !ids.is_empty() && ids.windows(2).all(|w| w[0] < w[1]))?;
            require(ids.iter().all(|&i| self.jobs.get(i).is_some_and(|j| j.observation.is_some())))?;
            require(self.model_updates > 0 && self.model_updates <= ids.len()
                && self.model_updates <= self.work.model_fit_attempts)?;
            if let Some(best) = &r.recommendation {
                require(best.observation_index < ids.len() && best.objective_variance >= 0.0
                    && best.constraint_variances.iter().all(|v| *v >= 0.0)
                    && best.joint_feasibility >= self.config.recommendation_probability
                    && best.joint_feasibility <= 1.0)?;
            }
            true
        } else {
            require(adaptive == 0 && self.model_updates == 0 && self.likelihood_evaluations == 0
                && self.posterior_only_fits == 0)?;
            if issued < n_initial { require(self.work == NoisyStudyWork::default())?; }
            false
        };
        let c = &self.config;
        let policies = std::iter::once((&c.search.kernel, c.objective_learning.as_ref()))
            .chain(c.constraints.iter().map(|o| (&o.kernel, o.learning.as_ref())));
        let mut minimum_probes = 0usize; let mut maximum_probes = 0usize;
        let mut planned_probes = 0usize;
        for ((original, learning), accepted) in policies.zip(&self.kernels) {
            require(original.family == accepted.family)?;
            if !modeled || learning.is_none() {
                require(original.signal.to_bits() == accepted.signal.to_bits()
                    && original.lengthscales.iter().zip(&accepted.lengthscales)
                        .all(|(a, b)| a.to_bits() == b.to_bits()))?;
            }
            if let Some(p) = learning {
                p.fit.validate(accepted).map_err(|_| CheckpointError::InvalidState)?;
                let planned = 1 + (self.max_evaluations - 1) / p.refit_every;
                planned_probes = planned.checked_mul(p.fit.max_evaluations)
                    .and_then(|n| planned_probes.checked_add(n)).ok_or(CheckpointError::LimitExceeded)?;
                let fits = if self.model_updates == 0 { 0 } else { 1 + (self.model_updates - 1) / p.refit_every };
                minimum_probes = minimum_probes.checked_add(fits).ok_or(CheckpointError::LimitExceeded)?;
                maximum_probes = fits.checked_mul(p.fit.max_evaluations)
                    .and_then(|n| maximum_probes.checked_add(n)).ok_or(CheckpointError::LimitExceeded)?;
            }
        }
        require((minimum_probes..=maximum_probes).contains(&self.likelihood_evaluations)
            && self.likelihood_evaluations <= planned_probes)?;
        let outputs = self.kernels.len();
        let required_conditioning = self.model_updates.checked_mul(outputs)
            .and_then(|n| n.checked_sub(minimum_probes)).ok_or(CheckpointError::InvalidState)?;
        let conditioning_passes = self.posterior_only_fits / outputs
            + usize::from(self.posterior_only_fits % outputs != 0);
        require(self.posterior_only_fits >= required_conditioning && conditioning_passes <= self.work.model_fit_attempts)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
