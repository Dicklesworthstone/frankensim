//! Defensive Gaussian importance sampling of physical failure probabilities.
//!
//! In independent normal germs, q(z)=epsilon*phi(z)+(1-epsilon)*phi(z-mu).
//! The target/proposal likelihood ratio is bounded by 1/epsilon. Multiplying
//! a failure indicator by epsilon times that ratio gives a [0,1] observation
//! for the existing bounded-mean confidence-sequence owner.
//! Bead: `frankensim-ext-uq-reliability-upgrades-pwui`.
pub mod proposal;

use core::fmt::Display;
use fs_eproc::GaussianMixtureCs;
use fs_evidence::Color;
use fs_math::det;

use crate::product_execution::sample_parameters;
use crate::product_plan::admit_plan;
use crate::product_sensitivity::Sum;
use crate::{
    AnytimeEstimate, CorrelationModel, UncertaintyKind, UqComplianceError, UqPlan, UqStatus,
};

/// A fixed proposal, chosen before observing this execution's outcomes.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportanceConfig {
    /// Shift in independent STANDARD-NORMAL germs, before the admitted target
    /// Cholesky transform; these are not physical parameter displacements.
    pub germ_shift: Vec<f64>,
    /// Original-target mixture mass, in (0,1], with a finite reciprocal.
    /// Larger values protect unexplored failure modes but reduce targeting.
    pub defensive_fraction: f64,
}

/// One paid physical evaluation and its original, unnormalized weight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImportanceObservation {
    pub qoi: f64,
    /// Failure means QoI > the immutable plan's compliance ceiling.
    pub failed: bool,
    /// Target/proposal ratio; never normalized by observed total weight.
    pub likelihood_weight: f64,
    /// epsilon * likelihood_weight, evaluated directly in log space in (0,1].
    pub scaled_weight: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImportanceEstimate {
    /// Mean weighted failure indicator. Finite-sample estimates can exceed one
    /// and are not clipped; this is not the proposal's observed failure rate.
    pub failure_probability: f64,
    /// Descriptive standard error, absent before two observations or outside
    /// finite numerical range. Zero error after no failures is not zero risk.
    pub standard_error: Option<f64>,
    /// (sum failure weights)^2 / sum squared failure weights; absent without
    /// observed failures. This diagnostic is not a confidence guarantee.
    pub failure_effective_samples: Option<f64>,
    /// Diagnostic E_q[target/proposal] estimate, whose population value is one.
    pub mean_likelihood_weight: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ImportanceReport {
    pub qoi_name: String,
    pub threshold: f64,
    pub evaluations_planned: usize,
    pub evaluations_attempted: usize,
    pub evaluations_accepted: usize,
    pub proposal_failures: usize,
    /// Descriptive statistics for the retained prefix; never supplied after a
    /// terminal failure, since that prefix could be selected by model refusal.
    pub estimate: Option<ImportanceEstimate>,
    pub unavailable_reason: Option<&'static str>,
    pub evidence_color: Color,
    pub status: UqStatus,
    pub rejection_reason: Option<String>,
}

/// Budgeted rare-event execution under an explicit Gaussian target law.
///
/// Physical samples retain the original units, order and admitted correlation,
/// including singular covariance. Weights live in independent germ space, so
/// no physical covariance inverse or determinant is needed. The evaluator and
/// threshold must remain fixed. Choose the proposal independently of this
/// run's observations. `proposal` can learn a shift from a completed pilot for
/// a SEPARATE execution; it never retroactively reweights this run.
///
/// Clone retains paid calls. Cancellation is checked between calls; expensive
/// solvers own their internal cancellation. The original one-million-call cap
/// bounds retained observations. No durable checkpoint or population-model
/// validation is implied by the estimate or its sampling-confidence sequence.
#[derive(Debug, Clone)]
pub struct ImportanceExecution {
    plan: UqPlan,
    germ_plan: UqPlan,
    config: ImportanceConfig,
    factor: Option<Vec<Vec<f64>>>,
    norm_squared: f64,
    unchanged: bool,
    values: Vec<ImportanceObservation>,
    attempted: usize,
    status: UqStatus,
    failure: Option<String>,
}

impl ImportanceExecution {
    pub fn new(plan: &UqPlan, config: ImportanceConfig) -> Result<Self, &'static str> {
        let factor = admit_plan(plan)?;
        if plan.compliance_threshold.is_none() {
            return Err("importance sampling requires a fixed compliance ceiling");
        }
        if !matches!(
            plan.correlation,
            CorrelationModel::Independent | CorrelationModel::JointGaussian { .. }
        ) || plan
            .parameters
            .iter()
            .any(|p| !matches!(p.kind, UncertaintyKind::AleatoryGaussian { .. }))
        {
            return Err(
                "importance sampling requires an explicit Independent or JointGaussian Gaussian law",
            );
        }
        let epsilon = config.defensive_fraction;
        if !epsilon.is_finite()
            || epsilon <= 0.0
            || epsilon > 1.0
            || !epsilon.recip().is_finite()
            || config.germ_shift.len() != plan.parameters.len()
        {
            return Err(
                "provide one finite germ shift per input and defensive fraction in (0,1] with finite reciprocal",
            );
        }
        let mut norm = Sum::default();
        for &shift in &config.germ_shift {
            if !shift.is_finite() || !(shift * shift).is_finite() {
                return Err("importance proposal shift norm is not finite");
            }
            norm.add(shift * shift);
            if !norm.value.is_finite() {
                return Err("importance proposal shift norm is not finite");
            }
        }
        let unchanged = epsilon == 1.0 || config.germ_shift.iter().all(|x| *x == 0.0);
        let mut germ_plan = plan.clone();
        germ_plan.correlation = CorrelationModel::Independent;
        for p in &mut germ_plan.parameters {
            p.kind = UncertaintyKind::AleatoryGaussian {
                mean: 0.0,
                std_dev: 1.0,
            };
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(plan.budget_max_samples)
            .map_err(|_| "cannot reserve admitted importance observation storage")?;
        Ok(Self {
            plan: plan.clone(),
            germ_plan,
            config,
            factor,
            norm_squared: norm.value,
            unchanged,
            values,
            attempted: 0,
            status: UqStatus::BudgetTruncated,
            failure: None,
        })
    }

    pub const fn plan(&self) -> &UqPlan {
        &self.plan
    }
    pub const fn config(&self) -> &ImportanceConfig {
        &self.config
    }
    pub fn observations(&self) -> &[ImportanceObservation] {
        &self.values
    }

    pub fn report(&self) -> ImportanceReport {
        let (estimate, unavailable_reason) =
            if self.status == UqStatus::Refused || self.values.is_empty() {
                (None, None)
            } else {
                match summarize(&self.values, self.config.defensive_fraction) {
                    Some(estimate) => (Some(estimate), None),
                    None => (
                        None,
                        Some("importance statistics exceed the finite numerical range"),
                    ),
                }
            };
        ImportanceReport {
            qoi_name: self.plan.target_qoi.clone(),
            threshold: self.plan.compliance_threshold.expect("admitted threshold"),
            evaluations_planned: self.plan.budget_max_samples,
            evaluations_attempted: self.attempted,
            evaluations_accepted: self.values.len(),
            proposal_failures: self.values.iter().filter(|v| v.failed).count(),
            evidence_color: Color::Estimated {
                estimator:
                    "defensive Gaussian importance sampling; descriptive error, declared model only"
                        .into(),
                dispersion: estimate
                    .and_then(|e| e.standard_error)
                    .unwrap_or(f64::INFINITY),
            },
            estimate,
            unavailable_reason,
            status: self.status,
            rejection_reason: self.failure.clone(),
        }
    }

    pub fn advance<F, E, C>(
        &mut self,
        allowance: usize,
        cancelled: C,
        mut evaluate: F,
    ) -> ImportanceReport
    where
        F: FnMut(&[f64]) -> Result<f64, E>,
        E: Display,
        C: FnMut() -> bool,
    {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }

    /// Ok(None) pauses an unfinished call; resume retries the SAME mixture
    /// component, germs, physical inputs and weight. Errors/nonfinite outputs
    /// permanently refuse, with no skipping or replacement. Panics propagate;
    /// callback side effects must tolerate retry of unfinished work.
    pub fn advance_interruptible<F, E, C>(
        &mut self,
        allowance: usize,
        mut cancelled: C,
        mut evaluate: F,
    ) -> ImportanceReport
    where
        F: FnMut(&[f64]) -> Result<Option<f64>, E>,
        E: Display,
        C: FnMut() -> bool,
    {
        if allowance == 0 || matches!(self.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let count = allowance.min(self.plan.budget_max_samples - self.values.len());
        for _ in 0..count {
            if cancelled() {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            let (parameters, scaled_weight, likelihood_weight) = match self.draw(self.values.len())
            {
                Ok(sample) => sample,
                Err(reason) => return self.fail(reason.into()),
            };
            let result = evaluate(&parameters);
            if matches!(&result, Ok(None)) {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            self.attempted += 1;
            let qoi = match result {
                Ok(Some(value)) if value.is_finite() => value,
                Ok(Some(_)) => {
                    return self.fail("importance physical model returned a non-finite QoI".into());
                }
                Err(error) => {
                    return self.fail(format!(
                        "importance physical evaluation {} refused: {error}",
                        self.values.len()
                    ));
                }
                Ok(None) => unreachable!("unfinished call did not consume an ordinal"),
            };
            self.values.push(ImportanceObservation {
                qoi,
                failed: qoi > self.plan.compliance_threshold.expect("admitted threshold"),
                likelihood_weight,
                scaled_weight,
            });
        }
        self.status = if self.values.len() == self.plan.budget_max_samples {
            UqStatus::Complete
        } else {
            UqStatus::BudgetTruncated
        };
        self.report()
    }

    /// Anytime sampling-confidence sequence for P(QoI > threshold).
    ///
    /// epsilon*w*failure lies in [0,1] under the declared mixture law, with
    /// conditional mean epsilon*p. Replays every paid value through the
    /// existing sigma=1/2, rho=1 Gaussian-mixture CS, then rescales by epsilon.
    /// When the proposal equals the target, known unit weights instead retain
    /// the ordinary-MC indicator sequence and its unscaled confidence radius.
    /// This conservative range-based bound can be MUCH wider than descriptive
    /// importance standard error, especially for rare events. Neither zero
    /// observed failures nor a small empirical error proves low target risk.
    ///
    /// Fix alpha before inspecting results. Multiple models/proposals/seeds or
    /// confidence choices need their own selection control. The continuous-law
    /// mathematical guarantee does not bound finite-grid transform/rounding,
    /// solver, model-discrepancy or physical-population error; no outward-rounded
    /// certificate is minted. Nonfinite or probability-incompatible results
    /// refuse. Clipping to [0,1] never makes `converged` artificially true.
    pub fn assess_failure_probability(
        &self,
        alpha: f64,
        target_half_width: f64,
    ) -> Result<Option<AnytimeEstimate>, UqComplianceError> {
        if !alpha.is_finite() || alpha <= 0.0 || alpha >= 1.0 || !alpha.recip().is_finite() {
            return Err(UqComplianceError::InvalidAlpha);
        }
        if !target_half_width.is_finite() || target_half_width < 0.0 {
            return Err(UqComplianceError::InvalidHalfWidth);
        }
        if self.status == UqStatus::Refused {
            return Err(UqComplianceError::RefusedExecution);
        }
        let mut confidence = GaussianMixtureCs::new(0.5, 1.0, alpha);
        for value in &self.values {
            confidence.observe(if !value.failed {
                0.0
            } else if self.unchanged {
                1.0
            } else {
                value.scaled_weight
            });
        }
        let Some((center, radius)) = confidence.interval() else {
            return Ok(None);
        };
        let epsilon = if self.unchanged {
            1.0
        } else {
            self.config.defensive_fraction
        };
        let center = center / epsilon;
        let radius = radius / epsilon;
        if !center.is_finite()
            || !radius.is_finite()
            || radius < 0.0
            || (center == 0.0 && self.values.iter().any(|value| value.failed))
        {
            return Err(UqComplianceError::NumericalRange);
        }
        let lo = (center - radius).max(0.0);
        let hi = (center + radius).min(1.0);
        if lo > hi {
            return Err(UqComplianceError::NumericalRange);
        }
        Ok(Some(AnytimeEstimate {
            mean: center,
            lo,
            hi,
            n: confidence.len(),
            converged: radius <= target_half_width,
        }))
    }

    fn fail(&mut self, reason: String) -> ImportanceReport {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }

    fn draw(&self, ordinal: usize) -> Result<(Vec<f64>, f64, f64), &'static str> {
        let epsilon = self.config.defensive_fraction;
        if self.unchanged {
            // Preserve the original target sampler BIT FOR BIT, including
            // singleton marginals, factor arithmetic and normal addressing.
            return Ok((
                sample_parameters(&self.plan, self.factor.as_deref(), ordinal)?,
                epsilon,
                1.0,
            ));
        }
        let germs = self.draw_germs(ordinal)?;
        let mut product = Sum::default();
        for (&z, &shift) in germs.iter().zip(&self.config.germ_shift) {
            let term = z * shift;
            if !z.is_finite() || !term.is_finite() {
                return Err("importance germ or log-density product overflowed");
            }
            product.add(term);
            if !product.value.is_finite() {
                return Err("importance log-density reduction overflowed");
            }
        }
        let ratio_log = product.value - 0.5 * self.norm_squared;
        let a = det::ln(epsilon);
        let b = det::ln(1.0 - epsilon) + ratio_log;
        if !b.is_finite() {
            return Err("importance log-density ratio is not finite");
        }
        let pivot = a.max(b);
        let denominator_log = pivot + det::ln(det::exp(a - pivot) + det::exp(b - pivot));
        let scaled = det::exp(a - denominator_log);
        let weight = scaled / epsilon;
        if !scaled.is_finite()
            || scaled <= 0.0
            || scaled > 1.0
            || !weight.is_finite()
            || weight <= 0.0
        {
            return Err("importance weight underflowed, overflowed or left its admitted range");
        }
        let mut physical = Vec::with_capacity(germs.len());
        for (i, parameter) in self.plan.parameters.iter().enumerate() {
            let z = self
                .factor
                .as_ref()
                .map_or_else(|| germs[i], |l| (0..=i).map(|j| l[i][j] * germs[j]).sum());
            let UncertaintyKind::AleatoryGaussian { mean, std_dev } = parameter.kind else {
                unreachable!("Gaussian admission");
            };
            let value = mean + std_dev * z;
            if !value.is_finite() {
                return Err("importance physical parameter overflowed");
            }
            physical.push(value);
        }
        Ok((physical, scaled, weight))
    }

    // Single owner for the actual mixture germs, including pilot replay.
    // No inversion of a possibly singular physical Cholesky map is involved.
    fn draw_germs(&self, ordinal: usize) -> Result<Vec<f64>, &'static str> {
        let mut selector = fs_rand::StreamKey {
            seed: self.plan.seed,
            kernel: 0x4953_4d58,
            tile: ordinal as u32,
        }.stream();
        let shifted = !self.unchanged
            && selector.next_f64() >= self.config.defensive_fraction;
        let mut germs = sample_parameters(&self.germ_plan, None, ordinal)?;
        if shifted {
            for (z, &shift) in germs.iter_mut().zip(&self.config.germ_shift) {
                *z += shift;
                if !z.is_finite() {
                    return Err("importance germ or log-density product overflowed");
                }
            }
        }
        Ok(germs)
    }
}

fn summarize(values: &[ImportanceObservation], epsilon: f64) -> Option<ImportanceEstimate> {
    let n = values.len() as f64;
    let event_scale = values
        .iter()
        .filter(|v| v.failed)
        .fold(0.0_f64, |s, v| s.max(v.scaled_weight));
    let weight_scale = values.iter().fold(0.0_f64, |s, v| s.max(v.scaled_weight));
    let mut weights = Sum::default();
    let mut events = Sum::default();
    let mut squares = Sum::default();
    for value in values {
        weights.add(value.scaled_weight / weight_scale);
        if value.failed {
            let y = value.scaled_weight / event_scale;
            events.add(y);
            squares.add(y * y);
        }
    }
    let average = events.value / n;
    let scale = event_scale / epsilon;
    let probability = average * scale;
    let mean_weight = (weights.value / n) * (weight_scale / epsilon);
    let standard_error = if values.len() < 2 {
        None
    } else {
        let mut variance = Sum::default();
        for value in values {
            let y = if value.failed {
                value.scaled_weight / event_scale
            } else {
                0.0
            };
            variance.add((y - average) * (y - average));
        }
        let error = det::sqrt(variance.value / (n - 1.0) / n) * scale;
        (error.is_finite() && (variance.value == 0.0 || error > 0.0)).then_some(error)
    };
    if !probability.is_finite()
        || (event_scale > 0.0 && probability == 0.0)
        || !mean_weight.is_finite()
        || mean_weight == 0.0
    {
        return None;
    }
    Some(ImportanceEstimate {
        failure_probability: probability,
        standard_error,
        failure_effective_samples: (event_scale > 0.0)
            .then(|| events.value * events.value / squares.value),
        mean_likelihood_weight: mean_weight,
    })
}
