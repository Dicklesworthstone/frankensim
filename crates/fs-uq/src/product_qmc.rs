//! Replicated randomized QMC for the same physical plans admitted by Monte Carlo.
//!
//! Reuses fs-rand's random-access Owen-scrambled Sobol nets and fs-bo's normal
//! quantile. Replicate means, NOT dependent points within a net, are the units
//! used for a descriptive standard error. No optional-stopping bound is implied.

#[path = "product_qmc_checkpoint.rs"]
mod checkpoint;

use core::fmt::Display;
use fs_evidence::Color;
use fs_rand::qmc::{MAX_SOBOL_DIM, Sobol};
use crate::{PropagationMethod, UncertaintyKind, UqExecution, UqPlan, UqStatus};

/// Explicit quadrature layout. The product must equal the plan's sample budget:
/// no silently discarded remainder, thinning, skipped origin, or partial net.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QmcConfig {
    /// Independently keyed scrambles, in 2..=256. Need not be a power of two.
    pub replicates: usize,
    /// Points per replicate: a power of two, at least two.
    pub samples_per_replicate: usize,
}

/// Mean over equally sized COMPLETE scrambles and its between-replicate error.
#[derive(Debug, Clone, PartialEq)]
pub struct QmcEstimate {
    /// Average of the retained replicate means.
    pub mean: f64,
    /// sqrt(sum((replicate_mean - mean)^2) / (R*(R-1))).
    /// Absent until two complete replicates exist; never a confidence interval.
    pub standard_error: Option<f64>,
}

/// Report from one retained execution. Incomplete net points count as paid work
/// but do not enter an estimate. Refused executions publish no estimates at all.
#[derive(Debug, Clone, PartialEq)]
pub struct QmcReport {
    /// The plan's physical quantity of interest.
    pub qoi_name: String,
    /// Exact, immutable replicate and within-replicate layout.
    pub config: QmcConfig,
    /// Original total evaluator allowance.
    pub samples_planned: usize,
    /// Completed evaluator calls, including a terminal failed/nonfinite call.
    /// An interrupted call is retried and not counted, like UqExecution.
    pub samples_evaluated: usize,
    /// All accepted points, including an unfinished replicate.
    pub samples_accepted: usize,
    /// Full replicates retained (even on failure, for work accounting only).
    pub completed_replicates: usize,
    /// Complete replicate means in their deterministic input order. Empty on failure.
    pub replicate_means: Vec<f64>,
    /// QoI mean and its descriptive randomized-quadrature standard error.
    pub estimate: Option<QmcEstimate>,
    /// Estimate for QoI <= threshold, only when a threshold was declared.
    /// Uses replicate proportions, not a Bernoulli iid formula on net points.
    pub compliance: Option<QmcEstimate>,
    /// Estimated only, never verified physical or interval evidence.
    pub evidence_color: Color,
    /// Complete means the entire predeclared sample budget was consumed.
    pub status: UqStatus,
    /// Original terminal sampling/model failure, without a replacement observation.
    pub rejection_reason: Option<String>,
}

/// Resumable bounded propagation for 1..=10 parameters using the existing Sobol
/// table. No pseudo-random tail is silently substituted beyond that limit.
///
/// Clone is an in-memory checkpoint; [`Self::checkpoint`] persists it.
/// The owned plan and layout cannot change on
/// resume, and the callback must keep the same model meaning. Sample ordinals
/// are replicate-major and random-access; interruption retries the exact point.
/// A rejected sample is terminal: continuing after discarding it would censor
/// the distribution. All accepted values are retained within the original cap.
///
/// Uniforms use the midpoints of fs-rand's 32-bit cells, so Gaussian quantiles
/// are finite WITHOUT skipping the first Sobol point or clamping its tails.
/// The resulting finite-grid quadrature and approximate inverse-normal map can
/// introduce bias. Between-replicate error does NOT bound that bias, solver or
/// physical-model error. Nor is zero replicate variation a proof of exactness.
#[derive(Debug, Clone)]
pub struct QmcExecution {
    plan: UqPlan,
    config: QmcConfig,
    factor: Option<Vec<Vec<f64>>>,
    sobol: Sobol,
    values: Vec<f64>,
    attempted: usize,
    status: UqStatus,
    failure: Option<String>,
}

impl QmcExecution {
    /// Admit the entire probability model and net layout before any evaluation.
    /// Only QuasiMonteCarlo plans are accepted. The existing Monte Carlo
    /// admission owns marginal, joint-Gaussian PSD, name, unit and total-budget
    /// validation; this adapter does not implement an alternative copula.
    pub fn new(plan: &UqPlan, config: QmcConfig) -> Result<Self, &'static str> {
        if plan.method != PropagationMethod::QuasiMonteCarlo {
            return Err("replicated QMC requires a QuasiMonteCarlo plan");
        }
        if !(1..=MAX_SOBOL_DIM).contains(&plan.parameters.len()) {
            return Err("replicated QMC supports 1..=10 declared parameters; no MC tail fallback");
        }
        if !(2..=256).contains(&config.replicates)
            || config.samples_per_replicate < 2
            || !config.samples_per_replicate.is_power_of_two()
            || config.replicates.checked_mul(config.samples_per_replicate)
                != Some(plan.budget_max_samples)
        {
            return Err("QMC needs 2..=256 replicates times a power-of-two point count >=2, exactly matching the sample budget");
        }
        // This temporary plan is for SHARED ADMISSION ONLY. No MC samples are
        // drawn, and the stored/reported plan remains QuasiMonteCarlo.
        let mut admission = plan.clone();
        admission.method = PropagationMethod::MonteCarlo;
        let admitted = UqExecution::new(&admission)?;
        Ok(Self {
            plan: plan.clone(), config, factor: admitted.factor,
            sobol: Sobol::new(plan.parameters.len()), values: Vec::new(),
            attempted: 0, status: UqStatus::BudgetTruncated, failure: None,
        })
    }

    /// Immutable model declaration and lifetime sample ceiling.
    #[must_use]
    pub fn plan(&self) -> &UqPlan { &self.plan }
    /// Immutable replicate layout.
    #[must_use]
    pub const fn config(&self) -> QmcConfig { self.config }
    /// Accepted observations, including the unfinished net, never sorted.
    #[must_use]
    pub fn observations(&self) -> &[f64] { &self.values }
    /// Completed/terminally rejected callback count; interruptions do not count.
    #[must_use]
    pub const fn evaluations_attempted(&self) -> usize { self.attempted }

    /// Report without invoking the model. Partial reports include only FULL
    /// replicates: they cannot mistake an interrupted net for a balanced rule.
    /// Calling this repeatedly costs O(accepted samples) and leaves state intact.
    #[must_use]
    pub fn report(&self) -> QmcReport {
        let complete = self.values.len() / self.config.samples_per_replicate;
        let mut means = Vec::new();
        let mut probabilities = Vec::new();
        if self.failure.is_none() {
            for block in self.values.chunks_exact(self.config.samples_per_replicate) {
                means.push(scaled_mean(block));
                if let Some(threshold) = self.plan.compliance_threshold {
                    probabilities.push(block.iter().filter(|&&x| x <= threshold).count() as f64
                        / self.config.samples_per_replicate as f64);
                }
            }
        }
        let estimate = summarize_replicates(&means);
        let compliance = summarize_replicates(&probabilities);
        let dispersion = estimate.as_ref().and_then(|s| s.standard_error).unwrap_or(0.0);
        let label = if self.failure.is_some() {
            "refused; no randomized-QMC estimate"
        } else if complete < 2 {
            "randomized-QMC; fewer than two complete scrambles; standard error unavailable"
        } else {
            "randomized-QMC between-scramble standard error; no stopping, discretization or model-error bound"
        };
        QmcReport {
            qoi_name: self.plan.target_qoi.clone(), config: self.config,
            samples_planned: self.plan.budget_max_samples, samples_evaluated: self.attempted,
            samples_accepted: self.values.len(), completed_replicates: complete,
            replicate_means: means, estimate, compliance,
            evidence_color: Color::Estimated { estimator: label.into(), dispersion },
            status: self.status, rejection_reason: self.failure.clone(),
        }
    }

    /// Execute at most this many additional evaluations. A zero allowance is a
    /// no-op. Poll cancellation between samples; completed values remain paid
    /// work. Complete and refused executions are terminal. Panics are not caught.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F) -> QmcReport
    where F: FnMut(&[f64]) -> Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }

    /// `Ok(None)` means an UNFINISHED model evaluation, never a discarded result.
    /// Stop and retry the same point on resume. `Err` and nonfinite values refuse
    /// permanently. Callback side effects must tolerate retry. Interrupted work
    /// is not counted as a completed evaluation; track its cost in the solver.
    pub fn advance_interruptible<F, E, C>(
        &mut self, allowance: usize, mut cancelled: C, mut evaluate: F,
    ) -> QmcReport
    where F: FnMut(&[f64]) -> Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        if allowance == 0 || matches!(self.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let count = allowance.min(self.plan.budget_max_samples - self.values.len());
        for _ in 0..count {
            if cancelled() { self.status = UqStatus::Cancelled; return self.report(); }
            let parameters = match self.parameters(self.values.len()) {
                Ok(x) => x,
                Err(reason) => return self.fail(reason.into()),
            };
            let evaluated = evaluate(&parameters);
            if matches!(&evaluated, Ok(None)) {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            self.attempted += 1;
            match evaluated {
                Ok(Some(x)) if x.is_finite() => self.values.push(x),
                Ok(Some(_)) => return self.fail("model returned a non-finite QoI".into()),
                Err(error) => return self.fail(format!("model evaluation {} refused: {error}", self.values.len())),
                Ok(None) => unreachable!("interruption returned without consuming the ordinal"),
            }
        }
        self.status = if self.values.len() == self.plan.budget_max_samples {
            UqStatus::Complete
        } else { UqStatus::BudgetTruncated };
        self.report()
    }

    fn fail(&mut self, reason: String) -> QmcReport {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }

    fn parameters(&self, ordinal: usize) -> Result<Vec<f64>, &'static str> {
        let replicate = ordinal / self.config.samples_per_replicate;
        let index = ordinal % self.config.samples_per_replicate;
        // Odd multiply plus addition is injective modulo 2^64. No replicate
        // seed collisions in the admitted range; no dependence on chunk sizes.
        let seed = self.plan.seed.wrapping_add((replicate as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15));
        let mut uniforms = vec![0.0; self.plan.parameters.len()];
        self.sobol.point_with_owen_scramble(index as u32, &mut uniforms, seed);
        for u in &mut uniforms { *u += 1.0 / 8_589_934_592.0; }
        let normals: Vec<f64> = if self.factor.is_some() {
            uniforms.iter().map(|&u| fs_bo::phi_inv(u)).collect()
        } else { Vec::new() };
        let mut values = Vec::with_capacity(uniforms.len());
        for (i, p) in self.plan.parameters.iter().enumerate() {
            let value = match p.kind {
                UncertaintyKind::AleatoryUniform { lo, hi } => {
                    if lo == hi { lo } else { (1.0 - uniforms[i]) * lo + uniforms[i] * hi }
                }
                UncertaintyKind::AleatoryGaussian { mean, std_dev } => {
                    if std_dev == 0.0 { mean } else {
                        let z = self.factor.as_ref().map_or_else(
                            || fs_bo::phi_inv(uniforms[i]),
                            |l| (0..=i).map(|j| l[i][j] * normals[j]).sum(),
                        );
                        mean + std_dev * z
                    }
                }
                _ => unreachable!("shared admission requires implemented probability measures"),
            };
            if !value.is_finite() { return Err("sampled QMC parameter overflowed"); }
            values.push(value);
        }
        Ok(values)
    }
}

// Scale BEFORE summing to avoid overflow for finite QoIs. A mean cannot lie
// outside [-scale,scale]; clamping only that rounding overshoot keeps constants
// at f64::MAX representable without changing the probability domain.
fn scaled_mean(values: &[f64]) -> f64 {
    let scale = values.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    if scale == 0.0 { return 0.0; }
    (values.iter().map(|x| x / scale).sum::<f64>() / values.len() as f64).clamp(-1.0, 1.0) * scale
}

fn summarize_replicates(values: &[f64]) -> Option<QmcEstimate> {
    if values.is_empty() { return None; }
    let mean = scaled_mean(values);
    let scale = values.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let standard_error = if values.len() < 2 { None } else if scale == 0.0 { Some(0.0) } else {
        let center = mean / scale;
        let squares: f64 = values.iter().map(|x| (x / scale - center).powi(2)).sum();
        let n = values.len() as f64;
        // Divide before rescaling: the error of the mean may be representable
        // even when the sample standard deviation would overflow.
        let error = fs_math::det::sqrt(squares / (n * (n - 1.0))) * scale;
        error.is_finite().then_some(error)
    };
    Some(QmcEstimate { mean, standard_error })
}
