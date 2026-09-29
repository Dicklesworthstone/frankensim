//! Closed-loop constrained q-NEI with independent, optionally learned outputs.
use super::{ConstraintModel, q_feasible_noisy_improvement};
use crate::gp::{Gp, Kernel};
use crate::hyper::{HeteroFitError, fit_heteroscedastic, zero_noise_duplicates};
use crate::learning::{NoisyFitRecord, NoisyLearningConfig};
use crate::noisy::{NoisyBoConfig, NoisyObservation, joint_normal_bank};

/// One physical evaluation; constraint entries follow configuration order.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstrainedObservation {
    /// Observed objective and its observation-noise variance.
    pub objective: NoisyObservation,
    /// Observed constraint outcomes and THEIR individual noise variances.
    pub constraints: Vec<NoisyObservation>,
}

/// Independent GP policy for an upper-bounded physical outcome.
#[derive(Debug, Clone)]
pub struct OutcomeConstraint {
    /// Initial covariance, with signal variance in squared outcome units.
    pub kernel: Kernel,
    /// Fixed prior mean in this outcome's original units.
    pub prior_mean: f64,
    /// Feasible means `outcome <= upper_bound`, in original units.
    pub upper_bound: f64,
    /// Optional bounded kernel learning; observation noise is never learned.
    pub learning: Option<NoisyLearningConfig>,
}

/// Explicit objective, constraints, acquisition and recommendation policy.
#[derive(Debug, Clone)]
pub struct ConstrainedBoConfig {
    /// Existing noisy-BO search limits, seed and INITIAL objective GP/prior.
    pub search: NoisyBoConfig,
    /// Optional objective kernel learning, independently scheduled.
    pub objective_learning: Option<NoisyLearningConfig>,
    /// Independent constraint models; an empty list is allowed.
    pub constraints: Vec<OutcomeConstraint>,
    /// Finite zero-utility objective reference, in original objective units.
    /// This affects acquisition, not the recommendation's feasibility test.
    pub reference: f64,
    /// Recommend only evaluated points whose estimated joint latent feasibility
    /// probability reaches this value in (0,1]. Not a frequentist bound.
    pub recommendation_probability: f64,
}

/// Lowest posterior-mean evaluated objective meeting the probability policy.
/// This is a model recommendation, not a certified feasible physical design.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstrainedRecommendation {
    /// Index into the aligned input/observation histories.
    pub observation_index: usize,
    /// Objective posterior mean and latent variance in original units.
    pub objective_mean: f64,
    /// Measurement noise is not re-added to this variance.
    pub objective_variance: f64,
    /// Constraint posterior means in declaration order, in original units.
    pub constraint_means: Vec<f64>,
    /// Corresponding latent variances.
    pub constraint_variances: Vec<f64>,
    /// Product of marginal Gaussian feasibility probabilities: outputs are
    /// modeled independently. Uses the existing approximate normal CDF.
    pub joint_feasibility: f64,
}

/// One scheduled learning stage for one output.
#[derive(Debug, Clone)]
pub struct OutputFitRecord {
    /// Zero is objective; j+1 is constraint j.
    pub output: usize,
    /// Refit data prefix, resolved seed, kernel, likelihood and actual work.
    pub fit: NoisyFitRecord,
}

/// Complete successful constrained study, including unresolved recommendations.
#[derive(Debug, Clone)]
pub struct ConstrainedBoReport {
    /// All inputs in actual callback order.
    pub x: Vec<Vec<f64>>,
    /// Original, unmodified objective and constraint observations.
    pub observations: Vec<ConstrainedObservation>,
    /// After initialization and every batch. `None` means no evaluated point
    /// meets the MODEL probability policy; it does not prove infeasibility.
    pub recommendations: Vec<Option<ConstrainedRecommendation>>,
    /// Scheduled kernel fits for all outputs, in stage/output order.
    pub fits: Vec<OutputFitRecord>,
    /// Actual likelihood probes, including rejected learning trials.
    pub likelihood_evaluations: usize,
    /// Ordinary conditioning fits between refits, or for fixed outputs.
    /// Acquisition joint factorizations are not counted in either fit counter.
    pub posterior_only_fits: usize,
}

/// Fail-closed admission, observation and model failures.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstrainedBoError {
    /// Invalid configuration or an unrepresentable declared shape/work bound.
    InvalidConfig,
    /// Wrong number of constraint observations at a callback.
    OutputCount { evaluation: usize, expected: usize, actual: usize },
    /// Non-finite value/centered value, or invalid observation variance.
    Observation { evaluation: usize, output: usize },
    /// A fixed or learned output model refused; no subsequent batch is run.
    Model { output: usize, source: HeteroFitError },
    /// Non-finite posterior mean/variance while forming recommendations.
    Posterior { output: usize },
}
impl std::fmt::Display for ConstrainedBoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "constrained noisy BO: {self:?}")
    }
}
impl std::error::Error for ConstrainedBoError {}

struct OutputPolicy<'a> {
    initial: &'a Kernel,
    prior: f64,
    learning: Option<&'a NoisyLearningConfig>,
}

impl ConstrainedBoConfig {
    fn outputs(&self) -> Vec<OutputPolicy<'_>> {
        std::iter::once(OutputPolicy { initial: &self.search.kernel, prior: self.search.prior_mean,
            learning: self.objective_learning.as_ref() }).chain(self.constraints.iter().map(|c|
                OutputPolicy { initial: &c.kernel, prior: c.prior_mean, learning: c.learning.as_ref() })).collect()
    }

    fn validate(&self, dim: usize, n_init: usize, iters: usize) -> Result<(), ConstrainedBoError> {
        let bad = ConstrainedBoError::InvalidConfig;
        let s = &self.search;
        let (lo, hi) = s.bounds;
        if !(1..=fs_rand::qmc::MAX_SOBOL_DIM).contains(&dim) || n_init == 0
            || u32::try_from(n_init).is_err() || s.q == 0
            || !lo.is_finite() || !hi.is_finite() || !(hi - lo).is_finite() || 0.2 * (hi - lo) <= 0.0
            || s.mc_samples == 0 || u32::try_from(s.mc_samples).is_err()
            || s.acq_starts == 0 || u32::try_from(s.acq_starts).is_err() || s.acq_evals == 0
            || !self.reference.is_finite() || !(self.reference - s.prior_mean).is_finite()
            || !self.recommendation_probability.is_finite()
            || self.recommendation_probability <= 0.0 || self.recommendation_probability > 1.0
        { return Err(bad); }
        let total = iters.checked_mul(s.q).and_then(|v| v.checked_add(n_init)).ok_or(bad.clone())?;
        let stages = iters.checked_add(1).ok_or(bad.clone())?;
        let outputs = self.constraints.len().checked_add(1).ok_or(bad.clone())?;
        total.checked_mul(total).ok_or(bad.clone())?;
        total.checked_mul(dim).ok_or(bad.clone())?;
        total.checked_mul(outputs).and_then(|v| v.checked_mul(s.mc_samples)).ok_or(bad.clone())?;
        stages.checked_mul(outputs).ok_or(bad.clone())?;
        let mut work = 0usize;
        for output in self.outputs() {
            let k = output.initial;
            if !output.prior.is_finite() || !k.signal.is_finite() || k.signal <= 0.0
                || k.lengthscales.len() != dim || k.lengthscales.iter().any(|v| !v.is_finite() || *v <= 0.0)
            { return Err(bad); }
            if let Some(policy) = output.learning {
                if policy.refit_every == 0 || policy.fit.validate(k).is_err() { return Err(bad); }
                let fits = 1 + iters / policy.refit_every;
                work = fits.checked_mul(policy.fit.max_evaluations).and_then(|v| work.checked_add(v)).ok_or(bad.clone())?;
            }
        }
        if self.constraints.iter().any(|c| !c.upper_bound.is_finite() || !(c.upper_bound - c.prior_mean).is_finite()) {
            return Err(bad);
        }
        Ok(())
    }
}

fn observe(f: &mut dyn FnMut(&[f64]) -> ConstrainedObservation, x: &[f64], c: &ConstrainedBoConfig,
    evaluation: usize) -> Result<ConstrainedObservation, ConstrainedBoError>
{
    let observed = f(x);
    if observed.constraints.len() != c.constraints.len() {
        return Err(ConstrainedBoError::OutputCount { evaluation, expected: c.constraints.len(), actual: observed.constraints.len() });
    }
    for (output, (value, policy)) in std::iter::once(&observed.objective).chain(&observed.constraints).zip(c.outputs()).enumerate() {
        if !value.value.is_finite() || !(value.value - policy.prior).is_finite()
            || !value.noise_variance.is_finite() || value.noise_variance < 0.0 {
            return Err(ConstrainedBoError::Observation { evaluation, output });
        }
    }
    Ok(observed)
}

fn condition(report: &mut ConstrainedBoReport, c: &ConstrainedBoConfig, kernels: &mut [Kernel], stage: usize)
    -> Result<Vec<Gp>, ConstrainedBoError>
{
    let mut models = Vec::new();
    for (output, policy) in c.outputs().into_iter().enumerate() {
        let observed: Vec<_> = report.observations.iter().map(|o|
            if output == 0 { o.objective } else { o.constraints[output - 1] }).collect();
        let y: Vec<_> = observed.iter().map(|o| o.value - policy.prior).collect();
        let noise: Vec<_> = observed.iter().map(|o| o.noise_variance).collect();
        let failure = |source| ConstrainedBoError::Model { output, source };
        if zero_noise_duplicates(&report.x, &noise) { return Err(failure(HeteroFitError::InvalidInitialModel)); }
        let model = if let Some(learning) = policy.learning.filter(|p| stage % p.refit_every == 0) {
            let mut fit = learning.fit.clone();
            fit.seed ^= c.search.seed ^ (stage as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)
                ^ (output as u64).wrapping_mul(0xD1B5_4A32_D192_ED03);
            let trained = fit_heteroscedastic(&report.x, &y, &noise, &kernels[output], &fit).map_err(failure)?;
            kernels[output] = trained.model.kernel.clone();
            report.likelihood_evaluations += trained.evaluations;
            report.fits.push(OutputFitRecord { output, fit: NoisyFitRecord {
                after_batches: stage, observation_count: report.x.len(), kernel: kernels[output].clone(),
                initial_lml: trained.initial_lml, lml: trained.model.lml, evaluations: trained.evaluations,
                projected_gradient_inf: trained.projected_gradient_inf,
                evaluation_limit_reached: trained.evaluation_limit_reached, seed: fit.seed,
            }});
            trained.model
        } else {
            let model = Gp::try_fit_diag(&report.x, &y, kernels[output].clone(), &noise)
                .filter(|m| m.lml.is_finite()).ok_or_else(|| failure(HeteroFitError::InvalidInitialModel))?;
            report.posterior_only_fits += 1;
            model
        };
        models.push(model);
    }
    Ok(models)
}

fn recommend(models: &[Gp], x: &[Vec<f64>], c: &ConstrainedBoConfig)
    -> Result<Option<ConstrainedRecommendation>, ConstrainedBoError>
{
    let mut best: Option<ConstrainedRecommendation> = None;
    for (index, point) in x.iter().enumerate() {
        let (centered, objective_variance) = models[0].predict(point);
        let objective_mean = centered + c.search.prior_mean;
        if !objective_mean.is_finite() || !objective_variance.is_finite() || objective_variance < 0.0 {
            return Err(ConstrainedBoError::Posterior { output: 0 });
        }
        let mut joint_feasibility = 1.0;
        let mut constraint_means = Vec::new();
        let mut constraint_variances = Vec::new();
        for (j, constraint) in c.constraints.iter().enumerate() {
            let (mean, variance) = models[j + 1].predict(point);
            let physical = mean + constraint.prior_mean;
            if !physical.is_finite() || !variance.is_finite() || variance < 0.0 {
                return Err(ConstrainedBoError::Posterior { output: j + 1 });
            }
            let bound = constraint.upper_bound - constraint.prior_mean;
            let p = if variance == 0.0 {
                if mean <= bound { 1.0 } else { 0.0 }
            } else {
                crate::acq::phi_cdf((bound - mean) / fs_math::det::sqrt(variance)).clamp(0.0, 1.0)
            };
            if !p.is_finite() { return Err(ConstrainedBoError::Posterior { output: j + 1 }); }
            joint_feasibility *= p;
            constraint_means.push(physical);
            constraint_variances.push(variance);
        }
        if joint_feasibility >= c.recommendation_probability
            && best.as_ref().is_none_or(|old| objective_mean < old.objective_mean) {
            best = Some(ConstrainedRecommendation { observation_index: index, objective_mean, objective_variance,
                constraint_means, constraint_variances, joint_feasibility });
        }
    }
    Ok(best)
}

// Same production CMA-ES engine as noisy BO, now scoring a multi-output utility.
fn argmax(acq: &dyn Fn(&[f64]) -> f64, dim: usize, s: &NoisyBoConfig, seed: u64) -> Vec<f64> {
    let (lo, hi) = s.bounds;
    let sobol = fs_rand::qmc::Sobol::scrambled(dim, seed);
    let mut point = vec![0.0; dim];
    let mut best: Option<(f64, Vec<f64>)> = None;
    for start in 0..s.acq_starts {
        sobol.point(u32::try_from(start + 1).expect("admitted starts"), &mut point);
        let initial: Vec<_> = point.iter().map(|u| (hi - lo).mul_add(*u, lo)).collect();
        let mut objective = |x: &[f64]| {
            let bounded: Vec<_> = x.iter().map(|v| v.clamp(lo, hi)).collect();
            let value = acq(&bounded);
            assert!(value.is_finite() && value >= 0.0, "invalid constrained acquisition");
            -value
        };
        let params = fs_dfo::CmaParams::standard(dim, 0.2 * (hi - lo), s.acq_evals, f64::NEG_INFINITY);
        let r = fs_dfo::cmaes(&mut objective, &initial, &params, seed ^ (start as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15));
        assert!(r.f_best.is_finite(), "constrained acquisition search must be finite");
        if best.as_ref().is_none_or(|(value, _)| -r.f_best > *value) {
            best = Some((-r.f_best, r.x_best.iter().map(|v| v.clamp(lo, hi)).collect()));
        }
    }
    best.expect("admitted nonzero starts").1
}

/// Optimize an objective subject to independently modeled noisy upper outcomes.
///
/// On success exactly `n_init + iters * search.q` callbacks are made. Every
/// callback supplies ALL outputs and their noise variances in original units.
/// All candidates of a greedy batch are selected before that batch's callbacks;
/// one point-major/output-minor normal bank supplies consistent row prefixes.
/// Every output conditions on the full observed history at each batch boundary.
/// Optional learning uses explicit per-output cadence, bounds and fit budgets.
///
/// No recommendation is fabricated when feasibility is unresolved. The threshold
/// controls only reported recommendations: candidate evaluations can violate
/// constraints. This is NOT safe BO. The reference defines capped improvement
/// and must cover objective values worth exploring; an overly optimistic
/// reference can suppress search. Recommendations may still exceed that reference.
///
/// Invalid policies refuse before callbacks. Invalid observations/model fits
/// stop before any subsequent evaluation. Existing objective callback panics
/// and non-finite acquisition arithmetic propagate; no successful report is
/// returned for them. The driver is synchronous and dense, without within-fit
/// or callback cancellation. It assumes independent outputs and observation
/// noise; correlated outputs/measurement noise require a different model.
/// Gauntlet tests cover bounded numerical fixtures, not physical validation.
pub fn minimize_constrained(
    f: &mut dyn FnMut(&[f64]) -> ConstrainedObservation,
    dim: usize,
    n_init: usize,
    iters: usize,
    config: &ConstrainedBoConfig,
) -> Result<ConstrainedBoReport, ConstrainedBoError> {
    config.validate(dim, n_init, iters)?;
    let s = &config.search;
    let (lo, hi) = s.bounds;
    let mut kernels: Vec<Kernel> = config.outputs().iter().map(|p| (*p.initial).clone()).collect();
    let mut report = ConstrainedBoReport { x: Vec::new(), observations: Vec::new(), recommendations: Vec::new(),
        fits: Vec::new(), likelihood_evaluations: 0, posterior_only_fits: 0 };
    let sobol = fs_rand::qmc::Sobol::scrambled(dim, s.seed);
    let mut point = vec![0.0; dim];
    for index in 0..n_init {
        sobol.point(u32::try_from(index + 1).expect("admitted initial count"), &mut point);
        let input: Vec<_> = point.iter().map(|u| (hi - lo).mul_add(*u, lo)).collect();
        report.observations.push(observe(f, &input, config, index)?);
        report.x.push(input);
    }
    let mut models = condition(&mut report, config, &mut kernels, 0)?;
    report.recommendations.push(recommend(&models, &report.x, config)?);
    let outputs = config.constraints.len() + 1;
    for iteration in 0..iters {
        let seed = s.seed ^ 0x434E_4549 ^ (iteration as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        let width = (report.x.len() + s.q) * outputs;
        let bank = joint_normal_bank(s.mc_samples, width, seed);
        let constraints: Vec<_> = config.constraints.iter().enumerate().map(|(j, c)|
            ConstraintModel { gp: &models[j + 1], upper_bound: c.upper_bound - c.prior_mean }).collect();
        let mut batch: Vec<Vec<f64>> = Vec::new();
        for slot in 0..s.q {
            let active_width = (report.x.len() + slot + 1) * outputs;
            let mut prefix = Vec::with_capacity(s.mc_samples * active_width);
            for row in bank.chunks_exact(width) { prefix.extend_from_slice(&row[..active_width]); }
            let candidate = argmax(&|input| {
                let mut trial = batch.clone(); trial.push(input.to_vec());
                q_feasible_noisy_improvement(&models[0], &constraints, &report.x, &trial,
                    config.reference - s.prior_mean, &prefix)
            }, dim, s, seed ^ (slot as u64).wrapping_mul(0xD1B5_4A32_D192_ED03));
            batch.push(candidate);
        }
        for input in batch {
            report.observations.push(observe(f, &input, config, report.x.len())?);
            report.x.push(input);
        }
        models = condition(&mut report, config, &mut kernels, iteration + 1)?;
        report.recommendations.push(recommend(&models, &report.x, config)?);
    }
    Ok(report)
}

#[cfg(test)]
mod tests;
