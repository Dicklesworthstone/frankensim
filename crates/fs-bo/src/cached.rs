//! Cached noisy acquisition: factor and sample the baseline once per search.
//!
//! All baseline points remain in the model: this does not prune observations or
//! approximate the covariance by independent marginals. Fixed posterior jitter
//! is explicitly declared; the existing adaptive-jitter entry points are intact.

use crate::gp::Gp;
use crate::gp::prepared::{PreparedError, PreparedPosterior};
use crate::noisy::{NoisyBoConfig, joint_normal_bank};

/// Cache a distinct baseline factor AND the minimum of each fixed latent draw.
/// Candidate scores only extend factor rows and sample the candidate values.
///
/// Bank layout is original baseline columns followed by `candidate_capacity`
/// columns. Duplicate baseline/candidate coordinates, including signed zero,
/// use their first occurrence's column. Baseline duplicates consume no new
/// random variable. Pending designs may be included in the baseline, without
/// adding their unknown outcomes to the fitted GP.
pub struct CachedNoisyImprovement<'gp, 'bank> {
    posterior: PreparedPosterior<'gp>,
    bank: &'bank [f64],
    width: usize,
    baseline_len: usize,
    columns: Vec<usize>,
    minima: Vec<f64>,
    capacity: usize,
}

fn check(keep_going: &mut dyn FnMut() -> bool) -> Result<(), PreparedError> {
    if keep_going() { Ok(()) } else { Err(PreparedError::Cancelled) }
}

impl<'gp, 'bank> CachedNoisyImprovement<'gp, 'bank> {
    /// Prepare one fixed bank for all probes and greedy slots in a search.
    pub fn new(gp: &'gp Gp, baseline: &[Vec<f64>], candidate_capacity: usize,
        bank: &'bank [f64], diagonal_jitter: f64) -> Result<Self, PreparedError>
    {
        Self::new_controlled(gp, baseline, candidate_capacity, bank, diagonal_jitter, &mut || true)
    }

    /// Polls between baseline projections/covariance rows and sampled minima.
    /// A dense Cholesky and individual triangular solves are not preemptible.
    pub fn new_controlled(gp: &'gp Gp, baseline: &[Vec<f64>], candidate_capacity: usize,
        bank: &'bank [f64], diagonal_jitter: f64, keep_going: &mut dyn FnMut() -> bool)
        -> Result<Self, PreparedError>
    {
        let width = baseline.len().checked_add(candidate_capacity).ok_or(PreparedError::ShapeOverflow)?;
        let dim = gp.kernel.lengthscales.len();
        if baseline.is_empty() || candidate_capacity == 0 || bank.is_empty()
            || bank.len() % width != 0 || bank.iter().any(|v| !v.is_finite())
            || baseline.iter().any(|p| p.len() != dim || p.iter().any(|v| !v.is_finite())) {
            return Err(PreparedError::InvalidInput);
        }
        let mut unique = Vec::new();
        let mut columns = Vec::new();
        for (i, point) in baseline.iter().enumerate() {
            check(keep_going)?;
            if !unique.contains(point) { unique.push(point.clone()); columns.push(i); }
        }
        let posterior = PreparedPosterior::new_controlled(gp, &unique, diagonal_jitter, keep_going)?;
        let mut minima = Vec::with_capacity(bank.len() / width);
        for row in bank.chunks_exact(width) {
            check(keep_going)?;
            let mut best = f64::INFINITY;
            for i in 0..unique.len() {
                let mut value = posterior.mean()[i];
                for j in 0..=i { value = posterior.factor().l(i, j).mul_add(row[columns[j]], value); }
                if !value.is_finite() { return Err(PreparedError::NonFinite); }
                best = best.min(value);
            }
            minima.push(best);
        }
        Ok(Self { posterior, bank, width, baseline_len: baseline.len(), columns, minima, capacity: candidate_capacity })
    }

    /// Number of distinct points factored once, including any pending designs.
    #[must_use]
    pub fn baseline_size(&self) -> usize { self.columns.len() }

    /// Score a candidate set using the leading candidate slots of this bank.
    pub fn score(&self, candidates: &[Vec<f64>]) -> Result<f64, PreparedError> {
        self.score_controlled(candidates, &mut || true)
    }

    /// Neither a score nor its refusal mutates the baseline cache. Identity
    /// duplicates have exactly zero incremental gain, not jitter-induced gain.
    pub fn score_controlled(&self, candidates: &[Vec<f64>], keep_going: &mut dyn FnMut() -> bool)
        -> Result<f64, PreparedError>
    {
        let dim = self.posterior.baseline()[0].len();
        if candidates.is_empty() || candidates.len() > self.capacity
            || candidates.iter().any(|p| p.len() != dim || p.iter().any(|v| !v.is_finite())) {
            return Err(PreparedError::InvalidInput);
        }
        check(keep_going)?;
        let mut novel = Vec::new();
        let mut columns = Vec::new();
        for (i, p) in candidates.iter().enumerate() {
            if !self.posterior.baseline().contains(p) && !novel.contains(p) {
                novel.push(p.clone()); columns.push(self.baseline_len + i);
            }
        }
        if novel.is_empty() { return Ok(0.0); }
        let extension = self.posterior.extend_controlled(&novel, keep_going)?;
        let b = self.columns.len(); let q = novel.len(); let stride = b + q;
        let mut average = 0.0;
        for (sample, row) in self.bank.chunks_exact(self.width).enumerate() {
            check(keep_going)?;
            let mut best = f64::INFINITY;
            for i in 0..q {
                let lower = &extension.lower_rows[i * stride..(i + 1) * stride];
                let mut value = extension.mean[i];
                for j in 0..b { value = lower[j].mul_add(row[self.columns[j]], value); }
                for j in 0..=i { value = lower[b + j].mul_add(row[columns[j]], value); }
                if !value.is_finite() { return Err(PreparedError::NonFinite); }
                best = best.min(value);
            }
            let gain = (self.minima[sample] - best).max(0.0);
            if !gain.is_finite() { return Err(PreparedError::NonFinite); }
            average += (gain - average) / (sample + 1) as f64;
        }
        Ok(average)
    }
}

/// An unsuccessful cached proposal publishes no partial candidate set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CachedProposalError {
    /// Invalid shape, count, coordinates or regularization declaration.
    InvalidRequest,
    /// Caller paused the search.
    Cancelled,
    /// Fixed-jitter posterior preparation/extension failed.
    Posterior(PreparedError),
    /// Bounded search selected an already evaluated/reserved design.
    NoNovelCandidate,
}
impl std::fmt::Display for CachedProposalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "cached noisy proposal: {self:?}")
    }
}
impl std::error::Error for CachedProposalError {}
impl From<PreparedError> for CachedProposalError {
    fn from(e: PreparedError) -> Self {
        if e == PreparedError::Cancelled { Self::Cancelled } else { Self::Posterior(e) }
    }
}

impl NoisyBoConfig {
    /// Opt-in replacement selection with one cached observed/pending baseline.
    /// Reuses production CMA-ES, seed mapping and common normal-bank columns.
    ///
    /// `diagonal_jitter` is a fixed nonnegative POSTERIOR variance, in squared
    /// objective units, not measurement noise. Unlike `propose_pending`, this
    /// path does not adapt jitter or fall back to independent draws. It can
    /// change selected designs and refuse ill-conditioned trials; choose this
    /// policy explicitly and retain it when replaying a decision.
    ///
    /// A numerical failure or cancellation poisons the complete proposal. The
    /// current bounded CMA-ES call may finish bookkeeping with inert callbacks,
    /// but no further posterior work or external simulation is performed.
    /// Preparation is per call (the GP and pending set may change on refills).
    /// No baseline pruning or additional objective calls are introduced.
    ///
    /// # Panics
    /// Invalid base optimizer settings retain existing BO admission assertions.
    pub fn propose_pending_cached(&self, gp: &Gp, observed: &[Vec<f64>], pending: &[Vec<f64>],
        count: usize, selection: usize, diagonal_jitter: f64, keep_going: &mut dyn FnMut() -> bool)
        -> Result<Vec<Vec<f64>>, CachedProposalError>
    {
        let dim = self.kernel.lengthscales.len();
        let baseline_len = observed.len().checked_add(pending.len()).ok_or(CachedProposalError::InvalidRequest)?;
        let width = baseline_len.checked_add(count).ok_or(CachedProposalError::InvalidRequest)?;
        if observed.is_empty() || count == 0 || count > self.q
            || gp.kernel.lengthscales.len() != dim || width.checked_mul(width).is_none()
            || width.checked_mul(self.mc_samples).is_none()
            || !diagonal_jitter.is_finite() || diagonal_jitter < 0.0
            || observed.iter().chain(pending).any(|p| p.len() != dim || p.iter().any(|v|
                !v.is_finite() || *v < self.bounds.0 || *v > self.bounds.1)) {
            return Err(CachedProposalError::InvalidRequest);
        }
        let mut policy = self.clone(); policy.q = count;
        policy.validate_study(dim, baseline_len, 1);
        check(keep_going)?;
        let seed = self.seed ^ 0x4E45_4942 ^ (selection as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        let bank = joint_normal_bank(self.mc_samples, width, seed);
        let baseline: Vec<_> = observed.iter().chain(pending).cloned().collect();
        let cached = CachedNoisyImprovement::new_controlled(gp, &baseline, count, &bank, diagonal_jitter, keep_going)?;
        let (lo, hi) = self.bounds;
        let mut batch: Vec<Vec<f64>> = Vec::new();
        for slot in 0..count {
            let slot_seed = seed ^ (slot as u64).wrapping_mul(0xD1B5_4A32_D192_ED03);
            let sobol = fs_rand::qmc::Sobol::scrambled(dim, slot_seed);
            let mut point = vec![0.0; dim];
            let mut best: Option<(f64, Vec<f64>)> = None;
            for start in 0..self.acq_starts {
                check(keep_going)?;
                sobol.point(u32::try_from(start + 1).expect("admitted starts"), &mut point);
                let initial: Vec<_> = point.iter().map(|u| (hi - lo).mul_add(*u, lo)).collect();
                let mut error = None;
                let mut objective = |x: &[f64]| {
                    if error.is_some() { return 0.0; }
                    let mut trial = batch.clone();
                    trial.push(x.iter().map(|v| v.clamp(lo, hi)).collect());
                    match cached.score_controlled(&trial, keep_going) {
                        Ok(value) => -value,
                        Err(e) => { error = Some(e); 0.0 }
                    }
                };
                let params = fs_dfo::CmaParams::standard(dim, 0.2 * (hi - lo), self.acq_evals, f64::NEG_INFINITY);
                let report = fs_dfo::cmaes(&mut objective, &initial, &params,
                    slot_seed ^ (start as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15));
                if let Some(error) = error { return Err(error.into()); }
                if !report.f_best.is_finite() { return Err(PreparedError::NonFinite.into()); }
                if best.as_ref().is_none_or(|(value, _)| -report.f_best > *value) {
                    best = Some((-report.f_best, report.x_best.iter().map(|v| v.clamp(lo, hi)).collect()));
                }
            }
            let candidate = best.expect("admitted starts").1;
            if baseline.contains(&candidate) || batch.contains(&candidate) {
                return Err(CachedProposalError::NoNovelCandidate);
            }
            batch.push(candidate);
        }
        check(keep_going)?;
        Ok(batch)
    }
}

#[cfg(test)]
mod tests;
