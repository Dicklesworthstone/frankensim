//! Incremental feasible utility beyond received AND unfinished designs.
use super::bo::select_batch;
use super::ConstrainedBoConfig;
use crate::gp::Gp;

/// A refused proposal never reserves an external simulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstrainedProposalError {
    /// Invalid config, model/output count, coordinates or allocation shape.
    InvalidRequest,
    /// Paused at a greedy-slot boundary; no partial proposal is returned.
    Cancelled,
    /// The bounded search returned an already received/reserved design.
    /// This does not prove that no unsearched novel design exists.
    NoNovelCandidate,
}
impl std::fmt::Display for ConstrainedProposalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "constrained pending proposal: {self:?}")
    }
}
impl std::error::Error for ConstrainedProposalError {}

impl ConstrainedBoConfig {
    /// Propose replacements while other objective/constraint evaluations run.
    ///
    /// `models[0]` is the centered objective GP and remaining models correspond
    /// to the declared constraints. Condition them ONLY on received results.
    /// Pending designs enter each output's joint latent draw, never its targets.
    /// The score is E[U(received union pending union proposed) -
    /// U(received union pending)], using the existing reference-capped feasible
    /// utility. An infeasible pending design therefore cannot block a feasible
    /// replacement merely by having an attractive objective.
    ///
    /// Uses the synchronous driver's bounded searches, seed policy and bank
    /// ordering. `selection` is a caller-owned logical decision ordinal; retries
    /// need the same ordinal, models and ordered availability sets. At most
    /// `search.q` designs are returned. Duplicates (including signed zero) refuse
    /// the whole proposal instead of manufacturing improvement through jitter.
    ///
    /// No model fitting, external dispatch, extra fallback search or safety
    /// guarantee is added. Candidate evaluations may violate constraints. The
    /// recommendation probability is NOT an admission rule for physical jobs.
    /// One slot's bounded search is non-preemptible; the continuation hook is
    /// checked between slots and before returning the complete proposal.
    ///
    /// # Panics
    /// Inherits GP/acquisition assertions for non-finite numerical results.
    pub fn propose_pending(&self, models: &[Gp], received: &[Vec<f64>],
        pending: &[Vec<f64>], count: usize, selection: usize,
        keep_going: &mut dyn FnMut() -> bool)
        -> Result<Vec<Vec<f64>>, ConstrainedProposalError>
    {
        let bad = ConstrainedProposalError::InvalidRequest;
        let dim = self.search.kernel.lengthscales.len();
        let baseline_len = received.len().checked_add(pending.len()).ok_or(bad)?;
        let points = baseline_len.checked_add(count).ok_or(bad)?;
        let outputs = self.constraints.len().checked_add(1).ok_or(bad)?;
        if received.is_empty() || count == 0 || count > self.search.q
            || models.len() != outputs
            || models.iter().any(|m| m.kernel.lengthscales.len() != dim)
            || self.validate(dim, points, 0).is_err()
            || points.checked_mul(points).is_none()
            || points.checked_mul(outputs).and_then(|v| v.checked_mul(self.search.mc_samples)).is_none()
            || received.iter().chain(pending).any(|p| p.len() != dim
                || p.iter().any(|v| !v.is_finite()
                    || *v < self.search.bounds.0 || *v > self.search.bounds.1))
        { return Err(bad); }
        if !keep_going() { return Err(ConstrainedProposalError::Cancelled); }
        let baseline: Vec<_> = received.iter().chain(pending).cloned().collect();
        let selected = select_batch(models, &baseline, count, selection, self, keep_going)
            .ok_or(ConstrainedProposalError::Cancelled)?;
        for (i, point) in selected.iter().enumerate() {
            if baseline.iter().chain(&selected[..i]).any(|p| p == point) {
                return Err(ConstrainedProposalError::NoNovelCandidate);
            }
        }
        Ok(selected)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Kernel, Matern};
    use crate::noisy::{NoisyBoConfig, joint_normal_bank};
    use crate::constrained::{ConstraintModel, OutcomeConstraint, q_feasible_noisy_improvement};

    fn config() -> ConstrainedBoConfig {
        let k = Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] };
        ConstrainedBoConfig {
            search: NoisyBoConfig { bounds: (0.0, 1.0), kernel: k.clone(), prior_mean: 0.0,
                q: 2, mc_samples: 64, acq_starts: 2, acq_evals: 32, seed: 91 },
            objective_learning: None,
            constraints: vec![OutcomeConstraint { kernel: k, prior_mean: 0.0,
                upper_bound: 0.0, learning: None }],
            reference: 2.0, recommendation_probability: 0.9,
        }
    }
    fn models(c: &ConstrainedBoConfig) -> Vec<Gp> {
        vec![Gp::fit(&[vec![0.2]], &[0.3], c.search.kernel.clone(), 0.05),
            Gp::fit(&[vec![0.2]], &[-1.0], c.constraints[0].kernel.clone(), 0.05)]
    }

    #[test]
    fn g3_no_pending_preserves_the_shared_batch_search() {
        let c = config(); let m = models(&c); let x = vec![vec![0.2]];
        let old = select_batch(&m, &x, 2, 7, &c, &mut || true).unwrap();
        let new = c.propose_pending(&m, &x, &[], 2, 7, &mut || true).unwrap();
        assert_eq!(new, old);
        let mut large_cap = c.clone(); large_cap.search.q = usize::MAX;
        // Allocation admission covers the ACTUAL requested count, not a
        // hypothetical full wave of the configured maximum batch size.
        assert_eq!(large_cap.propose_pending(&m, &x, &[], 2, 7, &mut || true).unwrap(), old);
    }

    #[test]
    fn g0_incremental_feasible_gain_matches_one_joint_draw_not_pending_means() {
        let c = config(); let m = models(&c);
        let x = vec![vec![0.2], vec![0.6], vec![0.9]];
        let bank = joint_normal_bank(2048, 6, 731);
        let objective = m[0].predict_joint(&x);
        let constraint = m[1].predict_joint(&x);
        let mut expected = 0.0;
        let mut total = 0.0;
        let mut pending_gain = 0.0;
        for row in bank.chunks_exact(6) {
            let mut utility = [0.0f64; 3];
            for i in 0..3 {
                let mut f = objective.0[i]; let mut g = constraint.0[i];
                for j in 0..=i {
                    f = objective.1[i * 3 + j].mul_add(row[2 * j], f);
                    g = constraint.1[i * 3 + j].mul_add(row[2 * j + 1], g);
                }
                if g <= 0.0 { utility[i] = (c.reference - f).max(0.0); }
            }
            let baseline = utility[0].max(utility[1]);
            expected += (utility[2] - baseline).max(0.0);
            total += utility[0].max(utility[1]).max(utility[2]) - utility[0];
            pending_gain += utility[0].max(utility[1]) - utility[0];
        }
        let constraints = [ConstraintModel { gp: &m[1], upper_bound: 0.0 }];
        let actual = q_feasible_noisy_improvement(&m[0], &constraints, &x[..2], &x[2..], c.reference, &bank);
        assert!(actual > 1e-3, "active replacement surface");
        assert!((actual - expected / 2048.0).abs() < 1e-12);
        assert!((actual - (total - pending_gain) / 2048.0).abs() < 1e-12);
    }

    #[test]
    fn g0_infeasible_pending_objective_does_not_suppress_feasible_gain() {
        let mut c = config(); c.search.kernel.lengthscales = vec![1e-5];
        c.constraints[0].kernel.lengthscales = vec![1e-5];
        let x = vec![vec![0.0], vec![0.5], vec![1.0]];
        let f = Gp::fit(&x, &[1.0, -100.0, 0.0], c.search.kernel.clone(), 1e-10);
        let g = Gp::fit(&x, &[-1.0, 10.0, -1.0], c.constraints[0].kernel.clone(), 1e-10);
        let models = [ConstraintModel { gp: &g, upper_bound: 0.0 }];
        // Zero bank isolates the exact stated mean draws. The running design
        // is infeasible, despite an objective lower than every received one.
        let gain = q_feasible_noisy_improvement(&f, &models, &x[..2], &x[2..], c.reference, &[0.0; 6]);
        assert!((gain - 1.0).abs() < 1e-6);
        assert_eq!(q_feasible_noisy_improvement(&f, &models, &x[..2], &x[1..2], c.reference, &[0.0; 6]), 0.0);
    }

    #[test]
    fn g4_proposal_admission_and_cancellation_publish_no_partial_batch() {
        let c = config(); let m = models(&c); let x = vec![vec![0.2]];
        for count in [0, 3, usize::MAX] {
            assert_eq!(c.propose_pending(&m, &x, &[], count, 0, &mut || true),
                Err(ConstrainedProposalError::InvalidRequest));
        }
        assert_eq!(c.propose_pending(&m[..1], &x, &[], 1, 0, &mut || true),
            Err(ConstrainedProposalError::InvalidRequest));
        assert_eq!(c.propose_pending(&m, &x, &[vec![f64::NAN]], 1, 0, &mut || true),
            Err(ConstrainedProposalError::InvalidRequest));
        let mut checks = 0;
        assert_eq!(c.propose_pending(&m, &x, &[vec![0.5]], 2, 0,
            &mut || { checks += 1; checks < 4 }), Err(ConstrainedProposalError::Cancelled));
        assert_eq!(checks, 4); // after the first slot, before the second
    }
}
