//! Versioned, exact-bit checkpoint bytes for fixed and learned external studies.
//! The v1 frame is eight-byte magic, little-endian u64 words (f64 bit patterns
//! for floating values), then FNV-1a64 of the preceding bytes. That checksum is
//! accidental-corruption detection, NOT authentication. Use trusted storage.
//! Request IDs are reconstructed from the ordered history, never renumbered.

use super::*;
use crate::gp::Matern;
use crate::hyper::HeteroFitConfig;
use crate::noisy::NoisyIncumbent;

const MAGIC: &[u8; 8] = b"FSNBO\0\x01\0";

/// Restore limits apply before allocation or numerical work. They bound the
/// frame and the full planned job count, not just already received results.
#[derive(Debug, Clone, Copy)]
pub struct CheckpointLimits {
    /// Largest accepted encoded frame, including header/checksum.
    pub max_bytes: usize,
    /// Largest permitted `n_init + batches * q` for a restored study.
    pub max_evaluations: usize,
}

/// A rejected checkpoint never mutates an existing study or issues jobs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckpointError {
    /// Caller limits or a nonrepresentable length were exceeded.
    LimitExceeded,
    /// Magic/schema version does not match this reader.
    IncompatibleVersion,
    /// A declared field or list extends past the available bytes.
    Truncated,
    /// Accidental corruption detected; this check is not a signature.
    ChecksumMismatch,
    /// Invalid policy, values, state relationships, flags, or trailing data.
    InvalidState,
}

impl std::fmt::Display for CheckpointError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "noisy study checkpoint: {self:?}")
    }
}
impl std::error::Error for CheckpointError {}
type Result<T> = std::result::Result<T, CheckpointError>;

fn require(ok: bool) -> Result<()> {
    if ok { Ok(()) } else { Err(CheckpointError::InvalidState) }
}
fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, b| (h ^ u64::from(*b)).wrapping_mul(0x100_0000_01b3))
}

struct Writer(Vec<u8>);
impl Writer {
    fn word(&mut self, value: u64) { self.0.extend_from_slice(&value.to_le_bytes()); }
    fn size(&mut self, value: usize) { self.word(u64::try_from(value).expect("64-bit checkpoint index")); }
    fn float(&mut self, value: f64) { self.word(value.to_bits()); }
    fn kernel(&mut self, kernel: &Kernel) {
        self.word(match kernel.family { Matern::Half => 0, Matern::ThreeHalves => 1, Matern::FiveHalves => 2 });
        self.float(kernel.signal);
        for v in &kernel.lengthscales { self.float(*v); }
    }
    fn observation(&mut self, observation: NoisyObservation) {
        self.float(observation.value); self.float(observation.noise_variance);
    }
}

struct Reader<'a>(&'a [u8]);
impl Reader<'_> {
    fn word(&mut self) -> Result<u64> {
        let bytes: [u8; 8] = self.0.get(..8).ok_or(CheckpointError::Truncated)?
            .try_into().map_err(|_| CheckpointError::Truncated)?;
        self.0 = &self.0[8..];
        Ok(u64::from_le_bytes(bytes))
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
    fn count(&mut self, max: usize, min_words: usize) -> Result<usize> {
        let n = self.size()?;
        if n > max { return Err(CheckpointError::LimitExceeded); }
        if n > self.0.len() / 8 / min_words { return Err(CheckpointError::Truncated); }
        Ok(n)
    }
    fn floats(&mut self, n: usize) -> Result<Vec<f64>> {
        if n > self.0.len() / 8 { return Err(CheckpointError::Truncated); }
        (0..n).map(|_| self.float()).collect()
    }
    fn kernel(&mut self, dim: usize) -> Result<Kernel> {
        let family = match self.word()? {
            0 => Matern::Half, 1 => Matern::ThreeHalves, 2 => Matern::FiveHalves,
            _ => return Err(CheckpointError::InvalidState),
        };
        let signal = self.float()?;
        let lengthscales = self.floats(dim)?;
        require(signal > 0.0 && lengthscales.iter().all(|v| *v > 0.0))?;
        Ok(Kernel { family, signal, lengthscales })
    }
    fn point(&mut self, dim: usize, bounds: (f64, f64)) -> Result<Vec<f64>> {
        let point = self.floats(dim)?;
        require(point.iter().all(|v| (bounds.0..=bounds.1).contains(v)))?;
        Ok(point)
    }
    fn observation(&mut self, prior: f64) -> Result<NoisyObservation> {
        let value = self.float()?; let noise_variance = self.float()?;
        require((value - prior).is_finite() && noise_variance >= 0.0)?;
        Ok(NoisyObservation { value, noise_variance })
    }
}

impl NoisyStudy {
    /// Encode all accepted and pending state without fitting, selecting jobs or
    /// rerunning simulations. Values retain their exact IEEE-754 bits. Persist
    /// these bytes atomically using caller-owned file/database operations.
    /// Replay requires the same numerical build/profile and callback outcomes;
    /// wire compatibility alone is not a cross-build determinism guarantee.
    #[must_use]
    pub fn checkpoint_bytes(&self) -> Vec<u8> {
        let mut w = Writer(MAGIC.to_vec());
        let stages = self.report.incumbent_trace.len();
        let initial = if stages == 0 { self.pending.len() }
            else { self.report.x.len() - (stages - 1) * self.config.q };
        w.word(self.study_id); w.size(self.dim); w.size(initial); w.size(self.iters);
        let c = &self.config;
        w.float(c.bounds.0); w.float(c.bounds.1); w.kernel(&c.kernel); w.float(c.prior_mean);
        for n in [c.q, c.mc_samples, c.acq_starts, c.acq_evals] { w.size(n); }
        w.word(c.seed);
        w.word(u64::from(self.learning.is_some()));
        if let Some(p) = &self.learning {
            w.size(p.refit_every);
            for (lo, hi) in &p.fit.lengthscale_bounds { w.float(*lo); w.float(*hi); }
            w.float(p.fit.signal_bounds.0); w.float(p.fit.signal_bounds.1);
            for n in [p.fit.starts, p.fit.max_iterations, p.fit.max_evaluations] { w.size(n); }
            w.float(p.fit.gradient_tolerance); w.word(p.fit.seed);
        }
        w.size(self.work.model_fit_attempts); w.size(self.work.acquisition_batch_attempts);
        w.size(stages); w.size(self.report.x.len());
        for (x, observation) in self.report.x.iter().zip(&self.report.observations) {
            for v in x { w.float(*v); } w.observation(*observation);
        }
        for incumbent in &self.report.incumbent_trace {
            w.size(incumbent.observation_index); w.float(incumbent.mean); w.float(incumbent.variance);
        }
        w.size(self.pending.len());
        for pending in &self.pending {
            for v in &pending.request.x { w.float(*v); }
            w.word(u64::from(pending.observation.is_some()));
            if let Some(o) = pending.observation { w.observation(o); }
        }
        w.size(self.fits.len());
        for fit in &self.fits {
            w.size(fit.after_batches); w.size(fit.observation_count); w.kernel(&fit.kernel);
            w.float(fit.initial_lml); w.float(fit.lml); w.size(fit.evaluations);
            w.float(fit.projected_gradient_inf); w.word(u64::from(fit.evaluation_limit_reached)); w.word(fit.seed);
        }
        w.word(checksum(&w.0));
        w.0
    }

    /// Restore a v1 checkpoint without invoking a GP, optimizer, or callback.
    /// All list lengths are bounded by both caller limits and available bytes
    /// before allocation. Reject malformed progress, policies, observations,
    /// refit cadence and budget records instead of constructing a partial study.
    /// Failed-model pending batches remain restorable for inspection/retry.
    ///
    /// This is structural admission plus accidental-corruption detection, NOT
    /// authentication or a numerical re-verification of previously stored fits.
    /// Only restore trusted checkpoints; the declared numerical work policy is
    /// retained and executes when the caller later requests `advance`.
    pub fn from_checkpoint(bytes: &[u8], limits: CheckpointLimits) -> Result<Self> {
        if bytes.len() > limits.max_bytes { return Err(CheckpointError::LimitExceeded); }
        if bytes.len() < 16 { return Err(CheckpointError::Truncated); }
        if &bytes[..8] != MAGIC { return Err(CheckpointError::IncompatibleVersion); }
        let split = bytes.len() - 8;
        let mut tail = Reader(&bytes[split..]);
        if tail.word()? != checksum(&bytes[..split]) { return Err(CheckpointError::ChecksumMismatch); }
        let mut r = Reader(&bytes[8..split]);
        let study_id = r.word()?; let dim = r.size()?; let initial = r.size()?; let iters = r.size()?;
        require((1..=fs_rand::qmc::MAX_SOBOL_DIM).contains(&dim))?;
        require(initial > 0 && u32::try_from(initial).is_ok())?;
        let bounds = (r.float()?, r.float()?); let kernel = r.kernel(dim)?; let prior_mean = r.float()?;
        let q = r.size()?; let mc_samples = r.size()?; let acq_starts = r.size()?; let acq_evals = r.size()?;
        let seed = r.word()?;
        let span = bounds.1 - bounds.0;
        require(span.is_finite() && span > 0.0 && 0.2 * span > 0.0 && q > 0
            && mc_samples > 0 && u32::try_from(mc_samples).is_ok()
            && acq_starts > 0 && u32::try_from(acq_starts).is_ok() && acq_evals > 0)?;
        let total = iters.checked_mul(q).and_then(|n| n.checked_add(initial))
            .ok_or(CheckpointError::LimitExceeded)?;
        if total > limits.max_evaluations { return Err(CheckpointError::LimitExceeded); }
        for n in [dim, total, mc_samples] {
            total.checked_mul(n).ok_or(CheckpointError::LimitExceeded)?;
        }
        let config = NoisyBoConfig { bounds, kernel, prior_mean, q, mc_samples, acq_starts, acq_evals, seed };
        let learning = if r.flag()? {
            let refit_every = r.size()?;
            if dim > r.0.len() / 16 { return Err(CheckpointError::Truncated); }
            let lengthscale_bounds = (0..dim).map(|_| Ok((r.float()?, r.float()?)))
                .collect::<Result<Vec<_>>>()?;
            let signal_bounds = (r.float()?, r.float()?);
            let starts = r.size()?; let max_iterations = r.size()?; let max_evaluations = r.size()?;
            let gradient_tolerance = r.float()?; let seed = r.word()?;
            let fit = HeteroFitConfig { lengthscale_bounds, signal_bounds, starts, max_iterations,
                max_evaluations, gradient_tolerance, seed };
            require(refit_every > 0 && fit.validate(&config.kernel).is_ok())?;
            (iters / refit_every).checked_add(1).and_then(|n| n.checked_mul(max_evaluations))
                .ok_or(CheckpointError::LimitExceeded)?;
            Some(NoisyLearningConfig { refit_every, fit })
        } else { None };
        let work = NoisyStudyWork { model_fit_attempts: r.size()?, acquisition_batch_attempts: r.size()? };
        let stages = r.size()?;
        require(stages <= iters + 1)?; // initial >= 1 and q >= 1 made iters + 1 representable
        let complete = stages == iters + 1;
        let count = r.count(total, dim + 2)?;
        let expected = if stages == 0 { 0 } else { initial + (stages - 1) * q };
        require(count == expected && work.model_fit_attempts >= stages
            && work.acquisition_batch_attempts >= stages.min(iters)
            && work.acquisition_batch_attempts <= work.model_fit_attempts)?;
        let mut report = NoisyBoReport { x: Vec::new(), observations: Vec::new(), incumbent_trace: Vec::new() };
        for _ in 0..count {
            report.x.push(r.point(dim, bounds)?); report.observations.push(r.observation(prior_mean)?);
        }
        if stages > r.0.len() / 24 { return Err(CheckpointError::Truncated); }
        for stage in 0..stages {
            let observation_index = r.size()?; let mean = r.float()?; let variance = r.float()?;
            require(observation_index < initial + stage * q && variance >= 0.0)?;
            report.incumbent_trace.push(NoisyIncumbent { observation_index, mean, variance });
        }
        let pending_count = r.count(total - count, dim + 1)?;
        require(pending_count == if complete { 0 } else if stages == 0 { initial } else { q })?;
        let mut pending = Vec::new();
        for slot in 0..pending_count {
            let x = r.point(dim, bounds)?;
            let observation = if r.flag()? { Some(r.observation(prior_mean)?) } else { None };
            pending.push(PendingEvaluation { request: NoisyRequest { study_id, evaluation: count + slot, x }, observation });
        }
        let fit_count = r.count(stages, dim + 10)?;
        let expected_fits = learning.as_ref().map_or(0, |p|
            if stages == 0 { 0 } else { 1 + (stages - 1) / p.refit_every });
        require(fit_count == expected_fits)?;
        let mut fits = Vec::new();
        for index in 0..fit_count {
            let policy = learning.as_ref().ok_or(CheckpointError::InvalidState)?;
            let after_batches = r.size()?; let observation_count = r.size()?; let kernel = r.kernel(dim)?;
            let initial_lml = r.float()?; let lml = r.float()?; let evaluations = r.size()?;
            let projected_gradient_inf = r.float()?; let evaluation_limit_reached = r.flag()?; let seed = r.word()?;
            require(after_batches == index * policy.refit_every && observation_count == initial + after_batches * q
                && kernel.family == config.kernel.family && policy.fit.validate(&kernel).is_ok()
                && lml >= initial_lml && projected_gradient_inf >= 0.0
                && evaluations > 0 && evaluations <= policy.fit.max_evaluations
                && evaluation_limit_reached == (evaluations == policy.fit.max_evaluations)
                && seed == (policy.fit.seed ^ config.seed ^ (after_batches as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)))?;
            fits.push(NoisyFitRecord { after_batches, observation_count, kernel, initial_lml, lml,
                evaluations, projected_gradient_inf, evaluation_limit_reached, seed });
        }
        require(r.0.is_empty())?;
        Ok(Self { study_id, config, dim, iters, report, pending, complete, work, learning, fits })
    }
}
