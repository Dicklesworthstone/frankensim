//! Opt-in kernel learning for the existing noisy BO acquisition/history loop.
//! Declared observation noise and prior mean remain unchanged throughout.

use crate::gp::{Gp, Kernel};
use crate::hyper::{HeteroFitConfig, HeteroFitError, fit_heteroscedastic, zero_noise_duplicates};
use crate::noisy::{NoisyBoConfig, NoisyBoReport, NoisyObservation};

/// Explicit refit policy. The original fixed-kernel entry point is unchanged.
#[derive(Debug, Clone)]
pub struct NoisyLearningConfig {
    /// Fit after initialization and every this many COMPLETED batches.
    /// The posterior still receives new observations after every batch.
    pub refit_every: usize,
    /// Per-refit kernel-learning bounds, restarts and likelihood-work budget.
    /// Each refit starts from the previous selected kernel. Seeds combine this
    /// seed, the BO seed and the completed-batch ordinal.
    pub fit: HeteroFitConfig,
}

/// Actual model learned from one prefix of the observation history.
#[derive(Debug, Clone)]
pub struct NoisyFitRecord {
    /// Zero means initialization; otherwise the completed-batch ordinal.
    pub after_batches: usize,
    /// Number of observations in this fit, including repeated observations.
    pub observation_count: usize,
    /// Selected covariance in the original coordinate/objective units.
    pub kernel: Kernel,
    /// LML of the warm-start kernel on this SAME data prefix.
    pub initial_lml: f64,
    /// LML of the selected kernel on this data prefix; at least initial_lml.
    pub lml: f64,
    /// Actual likelihood probes including the warm start and rejected probes.
    pub evaluations: usize,
    /// Unit projected log-gradient mapping at the selected kernel.
    pub projected_gradient_inf: f64,
    /// True when this refit's exact likelihood allowance was consumed.
    pub evaluation_limit_reached: bool,
    /// Resolved seed for replaying this refit's restart points.
    pub seed: u64,
}

/// Objective history, recommendations and measured model-learning work.
#[derive(Debug, Clone)]
pub struct NoisyLearnedReport {
    /// Unmodified callback observations and posterior-mean recommendations.
    pub report: NoisyBoReport,
    /// One record per scheduled refit, in completed-batch order.
    pub fits: Vec<NoisyFitRecord>,
    /// Sum of `fits[*].evaluations`. Does not include acquisition-posterior
    /// factorizations or the ordinary conditioning fits counted below.
    pub likelihood_evaluations: usize,
    /// Fits using the last learned kernel between refit stages. Each uses all
    /// currently available observations and exactly one training factorization.
    pub posterior_only_fits: usize,
}

/// Run noisy BO with bounded, warm-started heteroscedastic kernel learning.
///
/// On success there are exactly `n_init + iters * config.q` objective calls,
/// `1 + iters / learning.refit_every` learning stages and at most that many
/// times `learning.fit.max_evaluations` likelihood probes. A learning-stage
/// budget exhaustion retains its best evaluated kernel and is reported; it
/// does not terminate physical sampling or claim model convergence.
///
/// Every stage conditions on the full observed history with the ORIGINAL
/// variances. Only kernel signal variance and lengthscales are learned; the
/// supplied prior mean, noise variances and acquisition seed policy stay fixed.
/// Each greedy batch is selected entirely before its callbacks, using the
/// exact same acquisition engine as the fixed-kernel entry point.
///
/// Learning policy errors return before any objective call. Invalid initial
/// or later training covariances return an error before the next batch. Exact
/// duplicate zero-noise observations refuse, including signed-zero coordinates.
/// Noise-bearing replicates remain distinct observations in the likelihood.
///
/// # Panics
/// Invalid base BO configuration, observation values/variances and acquisition
/// arithmetic retain `noisy::minimize_noisy`'s existing panic semantics.
///
/// This driver remains synchronous/dense and has no within-acquisition or
/// objective-callback cancellation. It is point-estimate kernel learning, not
/// hyperparameter marginalization, a statistical certificate or a guaranteed
/// improvement over fixed-kernel optimization.
pub fn minimize_noisy_with_learning(
    f: &mut dyn FnMut(&[f64]) -> NoisyObservation,
    dim: usize,
    n_init: usize,
    iters: usize,
    config: &NoisyBoConfig,
    learning: &NoisyLearningConfig,
) -> Result<NoisyLearnedReport, HeteroFitError> {
    learning.fit.validate(&config.kernel)?;
    if learning.refit_every == 0 { return Err(HeteroFitError::InvalidConfig); }
    let stages = (iters / learning.refit_every).checked_add(1)
        .ok_or(HeteroFitError::InvalidConfig)?;
    stages.checked_mul(learning.fit.max_evaluations).ok_or(HeteroFitError::InvalidConfig)?;
    let mut kernel = config.kernel.clone();
    let mut fits = Vec::new();
    let mut likelihood_evaluations = 0usize;
    let mut posterior_only_fits = 0usize;
    let report = config.run_with_model(f, dim, n_init, iters, &mut |x, observations, stage| {
        let values: Vec<f64> = observations.iter().map(|o| o.value - config.prior_mean).collect();
        let variances: Vec<f64> = observations.iter().map(|o| o.noise_variance).collect();
        if zero_noise_duplicates(x, &variances) { return Err(HeteroFitError::InvalidInitialModel); }
        if stage % learning.refit_every == 0 {
            let mut policy = learning.fit.clone();
            policy.seed ^= config.seed ^ (stage as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
            let trained = fit_heteroscedastic(x, &values, &variances, &kernel, &policy)?;
            likelihood_evaluations += trained.evaluations; // admitted total bound
            kernel = trained.model.kernel.clone();
            fits.push(NoisyFitRecord { after_batches: stage, observation_count: x.len(),
                kernel: kernel.clone(), initial_lml: trained.initial_lml, lml: trained.model.lml,
                evaluations: trained.evaluations, projected_gradient_inf: trained.projected_gradient_inf,
                evaluation_limit_reached: trained.evaluation_limit_reached, seed: policy.seed });
            // Reuse the already evaluated winning GP. No unbudgeted refit.
            Ok(trained.model)
        } else {
            let model = Gp::try_fit_diag(x, &values, kernel.clone(), &variances)
                .filter(|gp| gp.lml.is_finite()).ok_or(HeteroFitError::InvalidInitialModel)?;
            posterior_only_fits += 1;
            Ok(model)
        }
    })?;
    Ok(NoisyLearnedReport { report, fits, likelihood_evaluations, posterior_only_fits })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::Matern;
    use crate::noisy::minimize_noisy;

    fn configs() -> (NoisyBoConfig, NoisyLearningConfig) {
        let bo = NoisyBoConfig { bounds: (0.0, 1.0),
            kernel: Kernel { family: Matern::FiveHalves, signal: 1.0, lengthscales: vec![0.3] },
            prior_mean: 0.5, q: 2, mc_samples: 16, acq_starts: 1, acq_evals: 16, seed: 97 };
        let learning = NoisyLearningConfig { refit_every: 2,
            fit: HeteroFitConfig { lengthscale_bounds: vec![(0.05, 2.0)], signal_bounds: (0.01, 5.0),
                starts: 2, max_iterations: 12, max_evaluations: 20, gradient_tolerance: 1e-7, seed: 13 } };
        (bo, learning)
    }

    fn observation(x: &[f64]) -> NoisyObservation {
        NoisyObservation { value: (x[0] - 0.37).powi(2), noise_variance: 0.01 + 0.03 * x[0] }
    }

    #[test]
    fn g3_one_evaluation_learning_matches_fixed_kernel_entire_trajectory() {
        let (c, mut learning) = configs();
        learning.fit.max_evaluations = 1;
        let learned = minimize_noisy_with_learning(&mut observation, 1, 3, 2, &c, &learning).unwrap();
        let fixed = minimize_noisy(&mut observation, 1, 3, 2, &c);
        assert_eq!(learned.report, fixed);
        assert_eq!(learned.likelihood_evaluations, 2);
        assert_eq!(learned.posterior_only_fits, 1);
        assert!(learned.fits.iter().all(|r| r.evaluation_limit_reached));
    }

    #[test]
    fn g3_refit_cadence_keeps_all_observations_and_warm_start_identity() {
        let (c, learning) = configs();
        let mut calls = Vec::new();
        let result = minimize_noisy_with_learning(&mut |x| { calls.push(x.to_vec()); observation(x) },
            1, 3, 3, &c, &learning).unwrap();
        assert_eq!(calls, result.report.x);
        assert_eq!(calls.len(), 9);
        assert_eq!(result.fits.iter().map(|f| (f.after_batches, f.observation_count)).collect::<Vec<_>>(),
            vec![(0, 3), (2, 7)]);
        assert_eq!(result.posterior_only_fits, 2);
        assert_eq!(result.likelihood_evaluations, result.fits.iter().map(|f| f.evaluations).sum::<usize>());
        assert!(result.likelihood_evaluations <= 2 * learning.fit.max_evaluations);
        let mut previous = c.kernel.clone();
        for fit in &result.fits {
            let n = fit.observation_count;
            let values: Vec<f64> = result.report.observations[..n].iter().map(|o| o.value - c.prior_mean).collect();
            let noise: Vec<f64> = result.report.observations[..n].iter().map(|o| o.noise_variance).collect();
            let old = Gp::try_fit_diag(&result.report.x[..n], &values, previous, &noise).unwrap();
            assert_eq!(old.lml.to_bits(), fit.initial_lml.to_bits());
            assert!(fit.lml >= fit.initial_lml);
            previous = fit.kernel.clone();
        }
        // Reconstruct every recommendation, including non-refit stages, from
        // the full prefix and its most recently learned covariance.
        for (stage, incumbent) in result.report.incumbent_trace.iter().enumerate() {
            let n = 3 + stage * c.q;
            let chosen = result.fits.iter().rev().find(|f| f.after_batches <= stage).unwrap();
            let y: Vec<f64> = result.report.observations[..n].iter().map(|o| o.value - c.prior_mean).collect();
            let noise: Vec<f64> = result.report.observations[..n].iter().map(|o| o.noise_variance).collect();
            let gp = Gp::try_fit_diag(&result.report.x[..n], &y, chosen.kernel.clone(), &noise).unwrap();
            let mut best = 0;
            for i in 1..n {
                if gp.predict(&result.report.x[i]).0 < gp.predict(&result.report.x[best]).0 { best = i; }
            }
            let (mean, variance) = gp.predict(&result.report.x[best]);
            assert_eq!(incumbent.observation_index, best);
            assert_eq!(incumbent.mean, mean + c.prior_mean);
            assert_eq!(incumbent.variance, variance);
        }
        for (x, o) in result.report.x.iter().zip(&result.report.observations) { assert_eq!(*o, observation(x)); }
    }

    #[test]
    fn g5_learned_study_and_resolved_models_replay() {
        let (c, learning) = configs();
        let a = minimize_noisy_with_learning(&mut observation, 1, 3, 2, &c, &learning).unwrap();
        let b = minimize_noisy_with_learning(&mut observation, 1, 3, 2, &c, &learning).unwrap();
        assert_eq!(a.report, b.report);
        assert_eq!(a.likelihood_evaluations, b.likelihood_evaluations);
        for (fa, fb) in a.fits.iter().zip(&b.fits) {
            assert_eq!(fa.kernel.lengthscales, fb.kernel.lengthscales);
            assert_eq!(fa.kernel.signal.to_bits(), fb.kernel.signal.to_bits());
            assert_eq!(fa.lml.to_bits(), fb.lml.to_bits());
            assert_eq!(fa.seed, fb.seed);
            assert_eq!(fa.projected_gradient_inf.to_bits(), fb.projected_gradient_inf.to_bits());
        }
    }

    #[test]
    fn g0_invalid_training_policy_refuses_before_callbacks() {
        let (c, policy) = configs();
        for variant in 0..4 {
            let mut bad = policy.clone();
            match variant {
                0 => bad.refit_every = 0,
                1 => bad.fit.max_evaluations = 0,
                2 => bad.fit.signal_bounds = (2.0, 5.0),
                _ => bad.fit.max_evaluations = usize::MAX,
            }
            let mut calls = 0;
            let result = minimize_noisy_with_learning(&mut |x| { calls += 1; observation(x) }, 1, 3, 2, &c, &bad);
            assert!(matches!(result, Err(HeteroFitError::InvalidConfig)));
            assert_eq!(calls, 0);
        }
    }

    #[test]
    fn g0_initial_only_study_learns_without_extra_objective_calls() {
        let (c, learning) = configs();
        let result = minimize_noisy_with_learning(&mut observation, 1, 3, 0, &c, &learning).unwrap();
        assert_eq!(result.report.x.len(), 3);
        assert_eq!(result.fits.len(), 1);
        assert_eq!(result.posterior_only_fits, 0);
        assert!(result.fits[0].lml > result.fits[0].initial_lml, "learning must be active");
    }

    #[test]
    fn g0_exact_constraint_aliases_are_rejected_but_noisy_replicates_are_not() {
        assert!(zero_noise_duplicates(&[vec![0.0], vec![-0.0]], &[0.0, 0.0]));
        assert!(!zero_noise_duplicates(&[vec![0.0], vec![-0.0]], &[0.0, 0.1]));
        assert!(!zero_noise_duplicates(&[vec![0.0], vec![1.0]], &[0.0, 0.0]));
    }
}
