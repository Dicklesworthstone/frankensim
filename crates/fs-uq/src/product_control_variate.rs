//! Frozen linear control variates for the mean of the actual model responses.
//!
//! For declared parameter means m and a gradient g fixed before sampling,
//! average Y_i - g^T(X_i - m). The control has known expectation zero; an
//! approximate adjoint can change variance, but does not replace the model or
//! change the target mean. Means come from the admitted probability law, never
//! from caller-supplied nominal values or the observed sample mean.
//!
//! This is NOT a distribution of the physical QoI. Quantiles, maxima, CVaR and
//! compliance indicators must continue to use the unadjusted observations.
//! Standard errors here are descriptive fixed-sample estimates, not confidence
//! sequences, physical error bounds, or an assertion of variance reduction.

mod checkpoint;

use core::fmt;

use crate::product_execution::sample_parameters;
use crate::{UncertaintyKind, UqExecution, UqPlan, UqStatus};

/// Admission, interruption, or numerical failure of a mean-only assessment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UqControlError {
    /// One finite coefficient per parameter, in exact declaration order, required.
    InvalidGradient,
    /// Coefficients cannot be installed after this execution has observed data.
    AlreadySampled,
    /// The frozen control belongs to a different complete plan.
    PlanMismatch,
    /// A failed sample cannot be dropped to obtain a more convenient estimate.
    RefusedExecution,
    /// Sampling, centering or a statistic exceeded finite f64 arithmetic.
    NumericalRange,
    /// The caller interrupted replay/statistics; no assessment was published.
    Cancelled,
}

impl fmt::Display for UqControlError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidGradient => "control requires one finite derivative per declared parameter",
            Self::AlreadySampled => "freeze the control before any model observations",
            Self::PlanMismatch => "control and execution have different uncertainty plans",
            Self::RefusedExecution => "a refused execution cannot publish a controlled mean",
            Self::NumericalRange => "control-variate arithmetic exceeds finite f64 range",
            Self::Cancelled => "control-variate assessment cancelled before publication",
        })
    }
}
impl std::error::Error for UqControlError {}

/// A fixed linear predictor bound to the entire admitted Monte Carlo plan.
///
/// Each derivative is in QoI units per corresponding declared parameter unit.
/// The caller owns the adjoint's parameter binding, numerical validity and model
/// identity. A wrong but fixed coefficient may INCREASE variance; it must not be
/// fitted to these observations or selected afterwards for a favorable mean.
/// Even a late-created fresh execution cannot legitimize such data-dependent
/// selection. The constructor's early-freeze check is not a proof of independence.
///
/// To recover, use [`Self::checkpoint`] and [`Self::restore`] to retain the
/// frozen coefficients and RAW responses together. Assessment deterministically
/// regenerates their parameters without rerunning either nominal or sample physics.
#[derive(Debug, Clone)]
pub struct LinearControlVariate {
    plan: UqPlan,
    means: Vec<f64>,
    gradient: Vec<f64>,
}

impl LinearControlVariate {
    /// Distribution means in parameter declaration order, not a sample estimate.
    #[must_use]
    pub fn parameter_means(&self) -> &[f64] { &self.means }

    /// The coefficients fixed before sampling, in declaration order.
    #[must_use]
    pub fn gradient(&self) -> &[f64] { &self.gradient }
}

/// Mean-only comparison on the SAME observations, always Estimated evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearControlEstimate {
    /// Number of actual model observations, including every accepted ordinal.
    pub n: usize,
    /// Unadjusted sample mean of the physical/model observable.
    pub raw_mean: f64,
    /// Mean of Y - g^T(X - E[X]); estimates the same declared model expectation.
    pub mean: f64,
    /// Unadjusted sample standard deviation; absent with fewer than two samples.
    pub raw_std_dev: Option<f64>,
    /// Sample standard deviation of the ADJUSTED estimator, not of physical Y.
    pub std_dev: Option<f64>,
    /// Unadjusted descriptive standard error, not an optional-stopping bound.
    pub raw_standard_error: Option<f64>,
    /// Adjusted descriptive standard error, not an optional-stopping bound.
    pub standard_error: Option<f64>,
    /// Adjusted/raw sample variance. Values above one are retained, not clipped.
    /// None for fewer than two samples, zero raw dispersion, or overflow.
    pub variance_ratio: Option<f64>,
}

impl UqExecution {
    /// Means of the admitted marginals, including a joint Gaussian's marginals.
    ///
    /// These are suitable for a deterministic nominal linearization. Correlation
    /// changes the joint sampler, not these marginal expectations. Interval-only
    /// and unspecified inputs were already refused by UqExecution construction.
    #[must_use]
    pub fn parameter_means(&self) -> Vec<f64> {
        self.plan.parameters.iter().map(|parameter| match parameter.kind {
            UncertaintyKind::AleatoryGaussian { mean, .. } => mean,
            UncertaintyKind::AleatoryUniform { lo, hi } => lo.midpoint(hi),
            _ => unreachable!("UqExecution admits probability measures only"),
        }).collect()
    }

    /// Freeze a linear mean control BEFORE this execution has observed a sample.
    ///
    /// No externally supplied expectation is accepted. The mean of the centered
    /// control is known from the exact declared marginals, not from a nominal
    /// design point which may differ from them. Does not alter the executor.
    ///
    /// # Errors
    /// Refuses failed/nonempty executions, wrong coefficient count or nonfinite
    /// coefficients. Choosing coefficients using OTHER observations remains the
    /// caller's responsibility and is not made valid by constructing a fresh run.
    pub fn freeze_linear_control_variate(
        &self, gradient: &[f64],
    ) -> Result<LinearControlVariate, UqControlError> {
        if self.status == UqStatus::Refused || self.failure.is_some() {
            return Err(UqControlError::RefusedExecution);
        }
        if self.attempted != 0 || !self.values.is_empty() {
            return Err(UqControlError::AlreadySampled);
        }
        if gradient.len() != self.plan.parameters.len() || !gradient.iter().all(|g| g.is_finite()) {
            return Err(UqControlError::InvalidGradient);
        }
        Ok(LinearControlVariate {
            plan: self.plan.clone(), means: self.parameter_means(), gradient: gradient.to_vec(),
        })
    }

    /// Assess a frozen control without invoking any model or changing observations.
    ///
    /// # Errors
    /// As [`Self::assess_linear_control_variate_interruptible`], without interruption.
    pub fn assess_linear_control_variate(
        &self, control: &LinearControlVariate,
    ) -> Result<Option<LinearControlEstimate>, UqControlError> {
        self.assess_linear_control_variate_interruptible(control, || false)
    }

    /// Replay the existing sampler and compute raw/controlled mean statistics.
    ///
    /// O(n*d) work and O(n+d) additional storage under the existing admitted sample
    /// and dimension caps. Polls before each regenerated sample and every 512
    /// entries in reductions. Raw observation bits, compliance results, sample
    /// counts and checkpoint bytes are unchanged, including after interruption.
    /// A partial execution's result is descriptive only. Do not stop a mean study
    /// on these standard errors or transfer them to a probability/quantile claim.
    ///
    /// # Errors
    /// Refuses plan substitution, failed executions, nonrepresentable arithmetic
    /// and caller interruption. No model evaluation or additional random draw is
    /// consumed, and no partial assessment is returned.
    pub fn assess_linear_control_variate_interruptible(
        &self, control: &LinearControlVariate, mut cancelled: impl FnMut() -> bool,
    ) -> Result<Option<LinearControlEstimate>, UqControlError> {
        if self.plan != control.plan { return Err(UqControlError::PlanMismatch); }
        if self.status == UqStatus::Refused || self.failure.is_some() {
            return Err(UqControlError::RefusedExecution);
        }
        poll(&mut cancelled)?;
        if self.values.is_empty() { return Ok(None); }
        let mut adjusted = Vec::with_capacity(self.values.len());
        for (ordinal, &value) in self.values.iter().enumerate() {
            poll(&mut cancelled)?;
            let parameters = sample_parameters(&self.plan, self.factor.as_deref(), ordinal)
                .map_err(|_| UqControlError::NumericalRange)?;
            let mut correction = 0.0;
            for ((&x, &mean), &gradient) in parameters.iter().zip(&control.means).zip(&control.gradient) {
                // A zero coefficient must not create 0*infinity in an unused term.
                if gradient != 0.0 {
                    correction = finite(gradient.mul_add(finite(x - mean)?, correction))?;
                }
            }
            adjusted.push(finite(value - correction)?);
        }
        let (raw_mean, raw_std_dev) = moments(&self.values, &mut cancelled)?;
        let (mean, std_dev) = moments(&adjusted, &mut cancelled)?;
        let root_n = (self.values.len() as f64).sqrt();
        let variance_ratio = raw_std_dev.zip(std_dev).and_then(|(raw, controlled)| {
            if raw == 0.0 { return None; }
            let ratio = (controlled / raw).powi(2);
            ratio.is_finite().then_some(ratio)
        });
        poll(&mut cancelled)?;
        Ok(Some(LinearControlEstimate {
            n: self.values.len(), raw_mean, mean, raw_std_dev, std_dev,
            raw_standard_error: raw_std_dev.map(|s| s / root_n),
            standard_error: std_dev.map(|s| s / root_n), variance_ratio,
        }))
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), UqControlError> {
    if cancelled() { Err(UqControlError::Cancelled) } else { Ok(()) }
}
fn finite(value: f64) -> Result<f64, UqControlError> {
    if value.is_finite() { Ok(value) } else { Err(UqControlError::NumericalRange) }
}

// Same scaled, ordered sample moments as product_execution's raw summary, with
// no sorting/quantile work: the adjusted samples are NOT a physical distribution.
fn moments(values: &[f64], cancelled: &mut impl FnMut() -> bool)
    -> Result<(f64, Option<f64>), UqControlError>
{
    let mut scale = 0.0_f64;
    for (i, value) in values.iter().enumerate() {
        if i % 512 == 0 { poll(cancelled)?; }
        scale = scale.max(value.abs());
    }
    let divisor = if scale == 0.0 { 1.0 } else { scale };
    let mut total = 0.0;
    for (i, value) in values.iter().enumerate() {
        if i % 512 == 0 { poll(cancelled)?; }
        total += value / divisor;
    }
    let center = total / values.len() as f64;
    let mean = finite(center * scale)?;
    let std_dev = if values.len() > 1 {
        let mut square_sum = 0.0;
        for (i, value) in values.iter().enumerate() {
            if i % 512 == 0 { poll(cancelled)?; }
            square_sum += (value / divisor - center).powi(2);
        }
        Some(finite((square_sum / (values.len() - 1) as f64).sqrt() * scale)?)
    } else { None };
    Ok((mean, std_dev))
}

#[cfg(test)]
mod tests;
