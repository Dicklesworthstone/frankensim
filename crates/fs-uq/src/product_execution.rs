//! Resumable, bounded execution of admitted Monte Carlo plans.
//!
//! A sample is addressed by its original plan seed and ordinal, not the chunk
//! size or process lifetime. Accepted observations are retained in ordinal order
//! so resumed statistics use exactly the same reductions as an uninterrupted run.
//! Cancellation is polled BETWEEN evaluations; a long-running solver must also
//! implement its own cancellation. Statistical stopping is not claimed here.

use core::fmt::Display;

use fs_evidence::Color;

use crate::product_plan::{
    UncertaintyKind, UqPlan, UqResult, UqStatus, admit_plan, refused,
};

/// An admitted Monte Carlo execution with an immutable lifetime sample budget.
///
/// The evaluator must represent the same deterministic model on every call.
/// A rejected or non-finite evaluation permanently refuses this execution:
/// dropping failed samples and continuing would bias the reported distribution.
/// Storage is bounded by the plan's admitted maximum of one million observations.
#[derive(Debug, Clone)]
pub struct UqExecution {
    pub(crate) plan: UqPlan,
    pub(crate) factor: Option<Vec<Vec<f64>>>,
    pub(crate) values: Vec<f64>,
    pub(crate) attempted: usize,
    pub(crate) status: UqStatus,
    pub(crate) failure: Option<String>,
}

impl UqExecution {
    /// Admit the complete plan before any model evaluation or sample allocation.
    ///
    /// # Errors
    /// Returns the same admission diagnostic as [`crate::UqPropagator::run`].
    pub fn new(plan: &UqPlan) -> Result<Self, &'static str> {
        let factor = admit_plan(plan)?;
        Ok(Self {
            plan: plan.clone(),
            factor,
            values: Vec::new(),
            attempted: 0,
            status: UqStatus::BudgetTruncated,
            failure: None,
        })
    }

    /// The immutable admitted plan, including its original lifetime budget.
    #[must_use]
    pub const fn plan(&self) -> &UqPlan {
        &self.plan
    }

    /// Accepted observations, in deterministic sample order (never sorted in place).
    #[must_use]
    pub fn observations(&self) -> &[f64] {
        &self.values
    }

    /// Number of evaluator calls, including a terminal rejected evaluation.
    #[must_use]
    pub const fn evaluations_attempted(&self) -> usize {
        self.attempted
    }

    /// Current report without generating a sample or invoking the evaluator.
    ///
    /// For zero observations, statistics are absent. For one observation (or an unrepresentable partial dispersion),
    /// variance and standard error are unavailable: inspect
    /// `std_dev`, not the required numeric `sampling_error` placeholder, first.
    /// Partial-run standard errors are descriptive, not optional-stopping bounds.
    #[must_use]
    pub fn report(&self) -> UqResult {
        if let Some(reason) = &self.failure {
            return refused(&self.plan, reason.clone(), self.attempted);
        }
        summarize(&self.plan, &self.values, self.status)
    }

    /// Execute at most `max_evaluations` additional samples, within the ORIGINAL
    /// plan budget. A zero allowance is a no-op. Cancelled and budget-truncated
    /// executions can continue; completed and refused executions are terminal.
    ///
    /// `cancelled` is checked before each sample. A completed evaluation is
    /// retained before the next check, so cancellation never discards paid work.
    /// Evaluator errors retain their diagnostic and never become successful,
    /// zero-valued, or silently skipped observations. Panics are not caught.
    pub fn advance<F, E, C>(
        &mut self,
        max_evaluations: usize,
        mut cancelled: C,
        mut evaluator: F,
    ) -> UqResult
    where
        F: FnMut(&[f64]) -> Result<f64, E>,
        E: Display,
        C: FnMut() -> bool,
    {
        if matches!(self.status, UqStatus::Complete | UqStatus::Refused)
            || max_evaluations == 0
        {
            return self.report();
        }
        let count = max_evaluations.min(self.plan.budget_max_samples - self.values.len());
        for _ in 0..count {
            if cancelled() {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            let parameters = match sample_parameters(
                &self.plan,
                self.factor.as_deref(),
                self.values.len(),
            ) {
                Ok(parameters) => parameters,
                Err(reason) => return self.fail(reason.into()),
            };
            self.attempted += 1;
            let value = match evaluator(&parameters) {
                Ok(value) if value.is_finite() => value,
                Ok(_) => return self.fail("model returned a non-finite QoI".into()),
                Err(error) => {
                    return self.fail(format!(
                        "model evaluation {} refused: {error}",
                        self.values.len()
                    ));
                }
            };
            self.values.push(value);
        }
        self.status = if self.values.len() == self.plan.budget_max_samples {
            UqStatus::Complete
        } else {
            UqStatus::BudgetTruncated
        };
        let report = self.report();
        if report.status == UqStatus::Refused {
            self.status = UqStatus::Refused;
            self.failure = report.rejection_reason.clone();
        }
        report
    }

    fn fail(&mut self, reason: String) -> UqResult {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }
}

fn sample_parameters(
    plan: &UqPlan,
    factor: Option<&[Vec<f64>]>,
    ordinal: usize,
) -> Result<Vec<f64>, &'static str> {
    // Admission caps the lifetime budget at 1_000_000, below u32::MAX.
    let mut stream = fs_rand::StreamKey {
        seed: plan.seed,
        kernel: 0x0517,
        tile: ordinal as u32,
    }
    .stream();
    let normals: Vec<f64> = if factor.is_some() {
        (0..plan.parameters.len())
            .map(|_| stream.next_normal())
            .collect()
    } else {
        Vec::new()
    };
    let mut parameters = Vec::with_capacity(plan.parameters.len());
    for (i, parameter) in plan.parameters.iter().enumerate() {
        let value = match parameter.kind {
            UncertaintyKind::AleatoryGaussian { mean, std_dev } => {
                let z = factor.map_or_else(
                    || stream.next_normal(),
                    |l| (0..=i).map(|j| l[i][j] * normals[j]).sum(),
                );
                mean + std_dev * z
            }
            UncertaintyKind::AleatoryUniform { lo, hi } => {
                let u = stream.next_f64();
                (1.0 - u) * lo + u * hi
            }
            _ => unreachable!("admission requires an implemented probability measure"),
        };
        if !value.is_finite() {
            return Err("sampled parameter overflowed");
        }
        parameters.push(value);
    }
    Ok(parameters)
}

fn summarize(plan: &UqPlan, values: &[f64], status: UqStatus) -> UqResult {
    let n = values.len();
    if n == 0 {
        let mut result = refused(plan, "no observations", 0);
        result.status = status;
        result.rejection_reason = None;
        result.evidence_color = Color::Estimated {
            estimator: "no observations; statistics unavailable".into(),
            dispersion: 0.0,
        };
        return result;
    }
    // Preserve the existing full-run reduction order and overflow protection.
    let scale = values.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let divisor = if scale == 0.0 { 1.0 } else { scale };
    let normalized_mean = values.iter().map(|x| x / divisor).sum::<f64>() / n as f64;
    let mean = normalized_mean * scale;
    let raw_std_dev = (n > 1).then(|| {
        let variance = values
            .iter()
            .map(|x| (x / divisor - normalized_mean).powi(2))
            .sum::<f64>()
            / (n - 1) as f64;
        variance.sqrt() * scale
    });
    if !mean.is_finite()
        || (status == UqStatus::Complete && raw_std_dev.is_some_and(|value| !value.is_finite()))
    {
        return refused(plan, "QoI statistics exceed finite f64 range", n);
    }
    // A partial sample standard deviation may exceed f64::MAX even when the
    // final one does not. Do not make work partitioning change execution.
    let std_dev = raw_std_dev.filter(|value| value.is_finite());
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let percentiles = [0.05, 0.50, 0.95].map(|p| sorted[((n as f64 * p) as usize).min(n - 1)]);
    let probability_of_compliance = plan.compliance_threshold.map(|threshold| {
        values.iter().filter(|&&value| value <= threshold).count() as f64 / n as f64
    });
    let sampling_error = std_dev.map_or(0.0, |sd| sd / (n as f64).sqrt());
    let estimator = if n == 1 {
        "one observation; variance and standard error unavailable"
    } else if std_dev.is_none() {
        "partial empirical-monte-carlo; dispersion exceeds finite f64 range"
    } else if status == UqStatus::Complete {
        "empirical-monte-carlo; no confidence or model-error bound"
    } else {
        "partial empirical-monte-carlo; no optional-stopping confidence or model-error bound"
    };
    UqResult {
        qoi_name: plan.target_qoi.clone(),
        method_used: plan.method,
        samples_evaluated: n,
        mean: Some(mean),
        std_dev,
        percentiles: Some(percentiles),
        interval_bounds: [sorted[0], sorted[n - 1]],
        probability_of_compliance,
        sampling_error,
        evidence_color: Color::Estimated {
            estimator: estimator.into(),
            dispersion: sampling_error,
        },
        status,
        rejection_reason: None,
    }
}
