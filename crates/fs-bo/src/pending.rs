//! Pending-aware noisy improvement and bounded replacement-batch selection.
//!
//! For minimization, integrate E[(min f(observed union pending) - min f(new))_+]
//! under ONE latent posterior conditioned only on received observations. This
//! is the marginal gain beyond designs already being evaluated, not a fake
//! observation or a posterior-mean fantasy. Equivalently it is the difference
//! between total joint improvement with and without the new designs, using
//! the same draws. See BoTorch's incremental qLogNEI semantics; no logarithmic
//! smoothing or new dependency is introduced here.

use crate::noisy::{NoisyBoConfig, q_noisy_expected_improvement};
use crate::gp::Gp;

/// A refused pending-aware proposal never issues external jobs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PendingProposalError {
    /// Empty observed set, zero/oversized count, incompatible dimensions,
    /// out-of-box/non-finite coordinates, or unrepresentable work shapes.
    InvalidRequest,
    /// Caller requested a pause before publication of a complete proposal.
    Cancelled,
    /// The bounded acquisition search returned an already reserved design.
    /// No replicate value-of-information is inferred from zero improvement.
    NoNovelCandidate,
}

impl std::fmt::Display for PendingProposalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "pending-aware noisy proposal: {self:?}")
    }
}
impl std::error::Error for PendingProposalError {}

/// Incremental q-NEI beyond currently running designs, in objective units.
///
/// Bank columns are observed, pending, then candidate coordinates. An empty
/// pending set exactly recovers ordinary q-NEI. Duplicate coordinates share
/// one latent draw; repeating a pending design therefore has zero gain even
/// if its unknown outcome has high variance. Pending results MUST NOT be put
/// in the GP training data until they arrive. The model and known baseline
/// are caller supplied: this function cannot inspect their provenance.
///
/// Uses the existing joint-posterior numerical jitter/fallback and finite-bank
/// integration. These remain model estimates, not confidence certificates or
/// estimates of information gained by replicating a noisy measurement.
///
/// # Panics
/// Requires at least one observed design and the same finite, rectangular
/// bank/coordinate admission as `q_noisy_expected_improvement`.
#[must_use]
pub fn q_pending_noisy_expected_improvement(
    gp: &Gp,
    observed: &[Vec<f64>],
    pending: &[Vec<f64>],
    candidates: &[Vec<f64>],
    bank: &[f64],
) -> f64 {
    assert!(!observed.is_empty(), "pending q-NEI needs received observations");
    let baseline: Vec<Vec<f64>> = observed.iter().chain(pending).cloned().collect();
    q_noisy_expected_improvement(gp, &baseline, candidates, bank)
}

impl NoisyBoConfig {
    /// Select up to `self.q` replacement jobs while older jobs are unfinished.
    ///
    /// Reuses the existing CMA-ES greedy q-NEI engine, treating pending designs
    /// as uncertain baseline values, never training targets. `selection` is a
    /// caller-owned logical decision ordinal: retry a refused decision with
    /// the same ordinal, observed set and pending set. Replay also requires
    /// the same GP and numerical profile. Point order fixes the normal-bank
    /// column identity. Cancellation is checked at the existing slot boundaries.
    ///
    /// An entire proposal is rejected if any returned coordinate equals an
    /// observed, pending, or earlier selected design (including signed zero).
    /// A flat/degenerate search may thus refuse despite other unsearched novel
    /// designs existing; no unbudgeted search or replication is substituted.
    ///
    /// # Panics
    /// Invalid base BO configuration retains `minimize_noisy`'s admission
    /// assertions. Arithmetic failures inside the GP/acquisition also retain
    /// their existing assertions. External dispatch is never performed here.
    pub fn propose_pending(
        &self,
        gp: &Gp,
        observed: &[Vec<f64>],
        pending: &[Vec<f64>],
        count: usize,
        selection: usize,
        keep_going: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<Vec<f64>>, PendingProposalError> {
        let dim = self.kernel.lengthscales.len();
        let baseline_len = observed.len().checked_add(pending.len())
            .ok_or(PendingProposalError::InvalidRequest)?;
        let width = baseline_len.checked_add(count)
            .ok_or(PendingProposalError::InvalidRequest)?;
        if observed.is_empty() || count == 0 || count > self.q
            || gp.kernel.lengthscales.len() != dim
            || width.checked_mul(width).is_none()
            || width.checked_mul(self.mc_samples).is_none()
            || observed.iter().chain(pending).any(|p| p.len() != dim
                || p.iter().any(|v| !v.is_finite()
                    || *v < self.bounds.0 || *v > self.bounds.1))
        {
            return Err(PendingProposalError::InvalidRequest);
        }
        // Admit all ordinary optimizer knobs before any acquisition work.
        let mut config = self.clone();
        config.q = count;
        config.validate_study(dim, baseline_len, 1);
        if !keep_going() { return Err(PendingProposalError::Cancelled); }
        let baseline: Vec<Vec<f64>> = observed.iter().chain(pending).cloned().collect();
        let candidates = config.select_batch(gp, &baseline, dim, selection, keep_going)
            .ok_or(PendingProposalError::Cancelled)?;
        for (i, candidate) in candidates.iter().enumerate() {
            if baseline.iter().chain(&candidates[..i]).any(|p| p == candidate) {
                return Err(PendingProposalError::NoNovelCandidate);
            }
        }
        Ok(candidates)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::{Kernel, Matern};
    use crate::noisy::joint_normal_bank;

    fn config() -> NoisyBoConfig {
        NoisyBoConfig { bounds: (-2.0, 2.0),
            kernel: Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.4] },
            prior_mean: 0.0, q: 2, mc_samples: 32, acq_starts: 1, acq_evals: 16, seed: 73 }
    }

    fn model() -> Gp {
        Gp::fit(&[vec![0.0]], &[0.2], config().kernel, 0.3)
    }

    #[test]
    fn no_pending_recovers_existing_qnei_and_search_bits() {
        let gp = model();
        let observed = vec![vec![0.0]];
        let candidates = vec![vec![0.7], vec![-0.8]];
        let bank = joint_normal_bank(256, 3, 13);
        assert_eq!(q_pending_noisy_expected_improvement(&gp, &observed, &[], &candidates, &bank).to_bits(),
            q_noisy_expected_improvement(&gp, &observed, &candidates, &bank).to_bits());
        let c = config();
        let expected = c.select_batch(&gp, &observed, 1, 4, &mut || true).unwrap();
        let actual = c.propose_pending(&gp, &observed, &[], 2, 4, &mut || true).unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn pending_nei_matches_joint_draw_marginal_not_plugin_fantasy() {
        let gp = model();
        let points = vec![vec![0.0], vec![0.6], vec![-0.8]];
        let bank = joint_normal_bank(2048, 3, 79);
        let (mean, lower) = gp.predict_joint(&points);
        let mut expected = 0.0;
        let mut total = 0.0;
        let mut already_reserved = 0.0;
        for row in bank.chunks_exact(3) {
            let mut f = mean.clone();
            for i in 0..3 {
                for j in 0..=i { f[i] = lower[i * 3 + j].mul_add(row[j], f[i]); }
            }
            expected += (f[0].min(f[1]) - f[2]).max(0.0);
            total += (f[0] - f[1].min(f[2])).max(0.0);
            already_reserved += (f[0] - f[1]).max(0.0);
        }
        let actual = q_pending_noisy_expected_improvement(&gp, &points[..1], &points[1..2], &points[2..], &bank);
        assert!(actual > 0.01);
        assert!((actual - expected / 2048.0).abs() < 1e-12);
        assert!((actual - (total - already_reserved) / 2048.0).abs() < 1e-12);
    }

    #[test]
    fn pending_duplicates_have_no_gain_and_nearby_points_stay_correlated() {
        let gp = model();
        let bank = joint_normal_bank(2048, 3, 83);
        assert_eq!(q_pending_noisy_expected_improvement(&gp, &[vec![0.0]],
            &[vec![1.0]], &[vec![1.0]], &bank), 0.0);
        let marginal = q_pending_noisy_expected_improvement(&gp, &[vec![0.0]],
            &[vec![1.0]], &[vec![1.00001]], &bank);
        let without_pending = q_noisy_expected_improvement(&gp, &[vec![0.0]],
            &[vec![1.00001]], &joint_normal_bank(2048, 2, 83));
        assert!(without_pending > 0.1);
        assert!(marginal < 0.001 * without_pending);
    }

    #[test]
    fn selection_is_bounded_replayable_and_never_reserves_a_duplicate() {
        let gp = model();
        let c = config();
        let observed = vec![vec![0.0]];
        let pending = vec![vec![0.4]];
        let a = c.propose_pending(&gp, &observed, &pending, 2, 9, &mut || true).unwrap();
        let b = c.propose_pending(&gp, &observed, &pending, 2, 9, &mut || true).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.len(), 2);
        assert!(a.iter().flatten().all(|v| (-2.0..=2.0).contains(v)));
        assert!(a.iter().all(|p| !observed.contains(p) && !pending.contains(p)));
        assert_ne!(a[0], a[1]);
        assert_eq!(c.propose_pending(&gp, &observed, &pending, 2, 9, &mut || false),
            Err(PendingProposalError::Cancelled));
        for count in [0, 3, usize::MAX] {
            assert_eq!(c.propose_pending(&gp, &observed, &pending, count, 9, &mut || true),
                Err(PendingProposalError::InvalidRequest));
        }
        assert_eq!(c.propose_pending(&gp, &observed, &[vec![3.0]], 1, 9, &mut || true),
            Err(PendingProposalError::InvalidRequest));
    }
}
