//! Product polynomial chaos with physical callbacks and held-out data.
//!
//! A fixed Philox design fits a Hermite/Legendre product surrogate; disjoint later
//! ordinals measure its error against the physical model. Coefficient moments
//! describe the surrogate, not a confidence interval or a certified model bound.
//! Bead: `frankensim-extreal-program-f85xj.6.7`.

use core::fmt::Display;
use fs_evidence::Color;
use fs_la::factor::{Qr, lu, qr};

use crate::pce::{total_degree_basis_size, total_degree_indices};
use crate::product_execution::sample_parameters;
use crate::product_pce_model::{PceBasis, ProductPceModel};
use crate::product_plan::admit_plan;
use crate::product_qmc::summarize_replicates;
use crate::{CorrelationModel, PropagationMethod, UncertaintyKind, UqPlan, UqStatus};

/// Fixed total-degree regression and a declared numerical conditioning limit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PceConfig {
    /// Total polynomial degree, in 0..=8. All declared parameters supply germs.
    pub degree: usize,
    /// Prefix used for fitting, at least twice the total-degree basis size.
    /// Remaining calls in the original plan are held out; at least two remain.
    pub training_samples: usize,
    /// Maximum estimated 1-norm condition of the scaled design's R factor.
    /// A numerical screen, not a rigorous condition-number bound.
    pub max_condition_estimate: f64,
}

impl PceConfig {
    /// Use the fixed degree/count with a conservative numerical screen of 1e8.
    #[must_use]
    pub const fn new(degree: usize, training_samples: usize) -> Self {
        Self {
            degree,
            training_samples,
            max_condition_estimate: 1e8,
        }
    }
}

/// Statistics available only after all fixed training and validation calls.
#[derive(Debug, Clone, PartialEq)]
pub struct PceSummary {
    /// Constant coefficient: expectation under the retained independent germ laws.
    pub surrogate_mean: f64,
    /// Sum of squared nonconstant coefficients, absent if not representable.
    /// This is response variance of the surrogate, not estimator variance.
    pub surrogate_variance: Option<f64>,
    /// Surrogate standard deviation, retained when its square is unrepresentable.
    pub surrogate_std_dev: Option<f64>,
    /// Empirical RMS physical-minus-surrogate residual on the training prefix.
    pub training_rmse: f64,
    /// Empirical RMS physical-minus-surrogate residual on disjoint held-out data.
    /// Neither this diagnostic nor its square bounds error at other inputs.
    pub validation_rmse: f64,
    /// Mean physical-minus-surrogate residual on the held-out observations.
    pub validation_mean_error: f64,
    /// Physical QoI <= threshold frequency using held-out observations only.
    /// No surrogate indicators or training observations enter this estimate.
    pub probability_of_compliance: Option<f64>,
    /// Descriptive iid standard error of the held-out compliance frequency.
    /// No confidence sequence, optional-stopping or rare-event guarantee.
    pub compliance_standard_error: Option<f64>,
}

/// Paid work, declared design, and completed fit diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct PceReport {
    /// Physical quantity of interest from the immutable input plan.
    pub qoi_name: String,
    /// Immutable regression and conditioning choices.
    pub config: PceConfig,
    /// Independent germ dimension, in physical declaration order.
    pub germ_dimension: usize,
    /// Number of graded-lex total-degree basis terms.
    pub basis_terms: usize,
    /// Numerical R-factor condition estimate computed before physical work.
    pub r_condition_estimate: f64,
    /// Lifetime physical-call budget, including training and held-out data.
    pub samples_planned: usize,
    /// Completed or terminally rejected calls. Unfinished calls consume no ordinal.
    pub samples_evaluated: usize,
    /// Paid finite training observations retained in their original order.
    pub training_samples_accepted: usize,
    /// Paid finite held-out observations retained in their original order.
    pub validation_samples_accepted: usize,
    /// Absent until the full original allocation completes, and on refusal.
    pub summary: Option<PceSummary>,
    /// Always Estimated with infinite dispersion: no total error bound is claimed.
    pub evidence_color: Color,
    /// Complete, paused by allowance/cancellation, or terminally refused.
    pub status: UqStatus,
    /// Original terminal diagnostic; observations are never skipped or replaced.
    pub rejection_reason: Option<String>,
}

/// Resumable physical-model execution for a law-aware product surrogate.
///
/// Admission supports 1..=16 Gaussian/uniform parameters, degree <= 8, <= 128 basis
/// terms and <= 65,536 total physical calls. Training work is capped at
/// 16,777,216 row/term/term products. Independent mixtures and explicitly joint
/// Gaussian laws reuse product admission, including singular laws and constant
/// marginals. Unknown dependence and non-Gaussian correlated laws still refuse.
///
/// The independent germs use the existing Philox seed/kernel/ordinal sampler.
/// Gaussian germs use Hermite polynomials; uniform [0,1] germs use shifted
/// Legendre polynomials. Physical Gaussian inputs are `mean_i + sigma_i *
/// (L * germ)_i`; uniform inputs use the original convex endpoint interpolation.
/// The admitted correlation factor L is identity for independence. Germ sensitivity
/// of the returned model therefore need not describe physical-coordinate effects.
/// The evaluator and its physical meaning must stay fixed through resumption.
///
/// Clone retains all paid observations. Cancellation is polled between physical
/// calls and before/after the final fit; solvers must poll their own scopes.
/// Bounded design construction, QR and final diagnostics are synchronous. There
/// is no durable checkpoint, degree adaptation or statistical stopping rule.
#[derive(Debug, Clone)]
pub struct PceExecution {
    plan: UqPlan,
    config: PceConfig,
    factor: Option<Vec<Vec<f64>>>,
    basis: Vec<PceBasis>,
    germs: Vec<Vec<f64>>,
    indices: Vec<Vec<usize>>,
    design_scale: f64,
    regression: Qr,
    condition: f64,
    values: Vec<f64>,
    attempted: usize,
    model: Option<ProductPceModel>,
    summary: Option<PceSummary>,
    status: UqStatus,
    failure: Option<String>,
}

impl PceExecution {
    /// Admit probability, storage, oversampling and numerical conditioning
    /// before any physical-model callback. The final two or more ordinals are
    /// reserved for validation and never influence regression coefficients.
    pub fn new(plan: &UqPlan, config: PceConfig) -> Result<Self, &'static str> {
        if plan.method != PropagationMethod::PolynomialChaos {
            return Err("PCE execution requires a PolynomialChaos plan");
        }
        let mut admitted = plan.clone();
        admitted.method = PropagationMethod::MonteCarlo;
        let factor = admit_plan(&admitted)?;
        let dim = plan.parameters.len();
        if dim > 16 || config.degree > 8 || plan.budget_max_samples > 65_536 {
            return Err("PCE supports 1..=16 Gaussian/uniform parameters, degree <=8 and <=65536 calls");
        }
        let basis: Vec<_> = plan.parameters.iter().map(|p| match p.kind {
            UncertaintyKind::AleatoryGaussian { .. } => PceBasis::StandardNormal,
            UncertaintyKind::AleatoryUniform { .. } => PceBasis::Uniform01,
            _ => unreachable!("product admission requires Gaussian or uniform marginals"),
        }).collect();
        if !config.max_condition_estimate.is_finite()
            || !(1.0..=1e8).contains(&config.max_condition_estimate)
        {
            return Err("PCE condition limit must be finite and in [1, 1e8]");
        }
        let terms =
            total_degree_basis_size(dim, config.degree).ok_or("PCE basis size overflows")?;
        if terms > 128
            || config.training_samples < 2 * terms
            || config.training_samples > plan.budget_max_samples.saturating_sub(2)
        {
            return Err("PCE needs <=128 terms, >=2*basis training samples and >=2 held-out calls");
        }
        if config
            .training_samples
            .checked_mul(terms)
            .and_then(|v| v.checked_mul(terms))
            .is_none_or(|work| work > 16_777_216)
        {
            return Err("PCE dense regression exceeds the bounded work envelope");
        }
        // Preserve each distribution's draw order using the canonical sampler
        // on standardized independent laws. No inverse transform or division
        // by a zero standard deviation/interval width is needed.
        let mut germ_plan = admitted;
        germ_plan.correlation = CorrelationModel::Independent;
        for (parameter, family) in germ_plan.parameters.iter_mut().zip(&basis) {
            parameter.kind = match family {
                PceBasis::StandardNormal => UncertaintyKind::AleatoryGaussian {
                    mean: 0.0, std_dev: 1.0,
                },
                PceBasis::Uniform01 => UncertaintyKind::AleatoryUniform { lo: 0.0, hi: 1.0 },
            };
        }
        let mut germs = Vec::new();
        germs
            .try_reserve_exact(plan.budget_max_samples)
            .map_err(|_| "PCE germ storage allocation refused")?;
        for ordinal in 0..plan.budget_max_samples {
            germs.push(sample_parameters(&germ_plan, None, ordinal)?);
        }
        let indices = total_degree_indices(dim, config.degree);
        let mut design = Vec::new();
        design
            .try_reserve_exact(config.training_samples * terms)
            .map_err(|_| "PCE design storage allocation refused")?;
        for germ in &germs[..config.training_samples] {
            for index in &indices {
                let value = index
                    .iter()
                    .zip(germ)
                    .zip(&basis)
                    .map(|((&degree, &x), &family)| family.value(degree, x))
                    .product::<f64>();
                if !value.is_finite() {
                    return Err("PCE design basis overflowed");
                }
                design.push(value);
            }
        }
        // A common scale preserves the least-squares problem and keeps the
        // existing QR kernel's squared column norms away from overflow.
        let design_scale = design.iter().fold(1.0_f64, |scale, x| scale.max(x.abs()));
        for value in &mut design {
            *value /= design_scale;
        }
        let regression = qr(&design, config.training_samples, terms);
        let mut r = vec![0.0; terms * terms];
        for i in 0..terms {
            for j in i..terms {
                r[i * terms + j] = regression.r(i, j);
            }
        }
        if r.iter().any(|x| !x.is_finite()) {
            return Err("PCE QR factor is non-finite");
        }
        let scale = (0..terms)
            .map(|i| r[i * terms + i].abs())
            .fold(0.0_f64, f64::max);
        if (0..terms).any(|i| r[i * terms + i].abs() <= 1e-12 * scale) {
            return Err("PCE regression design has a numerically deficient pivot");
        }
        let condition = lu(&r, terms)
            .map_err(|_| "PCE regression design is singular")?
            .condition_1(&r);
        if !condition.is_finite() || condition > config.max_condition_estimate {
            return Err("PCE regression exceeds the declared R-factor condition estimate");
        }
        let mut values = Vec::new();
        values
            .try_reserve_exact(plan.budget_max_samples)
            .map_err(|_| "PCE observation storage allocation refused")?;
        Ok(Self {
            plan: plan.clone(),
            config,
            factor,
            basis,
            germs,
            indices,
            design_scale,
            regression,
            condition,
            values,
            attempted: 0,
            model: None,
            summary: None,
            status: UqStatus::BudgetTruncated,
            failure: None,
        })
    }

    /// Immutable physical plan, retaining the PolynomialChaos method.
    #[must_use]
    pub const fn plan(&self) -> &UqPlan {
        &self.plan
    }

    /// Fixed independent germs; training is the prefix in `config.training_samples`.
    #[must_use]
    pub fn germs(&self) -> &[Vec<f64>] {
        &self.germs
    }

    /// Reference laws and polynomial families, available before any model call.
    #[must_use]
    pub fn basis(&self) -> &[PceBasis] {
        &self.basis
    }

    /// Accepted physical QoIs, in original ordinal order, including held-out data.
    #[must_use]
    pub fn observations(&self) -> &[f64] {
        &self.values
    }

    /// Completed or terminally rejected physical calls; interruption is unfinished.
    #[must_use]
    pub const fn evaluations_attempted(&self) -> usize {
        self.attempted
    }

    /// Completed surrogate bound to its independent germ laws. The returned
    /// `ProductPceModel` supports evaluation, gradients and coefficient Sobol
    /// indices; none is a certified derivative or sensitivity of physics.
    #[must_use]
    pub fn model(&self) -> Option<&ProductPceModel> {
        (self.status == UqStatus::Complete)
            .then_some(self.model.as_ref())
            .flatten()
    }

    /// Map a finite germ vector through the admitted physical probability law.
    /// This is the exact transformation used for every callback, also for
    /// singular joint laws; no inverse covariance or whitening is required.
    /// Uniform coordinates must lie in [0,1], including for constant marginals.
    pub fn physical_parameters(&self, germ: &[f64]) -> Result<Vec<f64>, &'static str> {
        if germ.len() != self.plan.parameters.len() || germ.iter().any(|x| !x.is_finite()) {
            return Err("PCE germ must have the declared dimension and finite coordinates");
        }
        self.plan
            .parameters
            .iter()
            .enumerate()
            .map(|(i, parameter)| {
                let value = match parameter.kind {
                    UncertaintyKind::AleatoryGaussian { mean, std_dev } => {
                        let z = self.factor.as_ref()
                            .map_or(germ[i], |l| (0..=i).map(|j| l[i][j] * germ[j]).sum());
                        if std_dev == 0.0 { mean } else { mean + std_dev * z }
                    }
                    UncertaintyKind::AleatoryUniform { lo, hi } => {
                        let u = germ[i];
                        if !(0.0..=1.0).contains(&u) {
                            return Err("uniform PCE germs must lie in [0,1]");
                        }
                        if lo == hi { lo } else { (1.0 - u) * lo + u * hi }
                    }
                    _ => unreachable!("admission requires Gaussian or uniform marginals"),
                };
                if value.is_finite() {
                    Ok(value)
                } else {
                    Err("PCE sampled physical parameter overflowed")
                }
            })
            .collect()
    }

    /// Summarize paid work without callbacks or refitting. Partial and refused
    /// states expose counts only, never a provisional final surrogate estimate.
    #[must_use]
    pub fn report(&self) -> PceReport {
        let complete = self.status == UqStatus::Complete;
        PceReport {
            qoi_name: self.plan.target_qoi.clone(), config: self.config,
            germ_dimension: self.plan.parameters.len(), basis_terms: self.indices.len(),
            r_condition_estimate: self.condition, samples_planned: self.plan.budget_max_samples,
            samples_evaluated: self.attempted,
            training_samples_accepted: self.values.len().min(self.config.training_samples),
            validation_samples_accepted: self.values.len().saturating_sub(self.config.training_samples),
            summary: if complete { self.summary.clone() } else { None },
            evidence_color: Color::Estimated {
                estimator: if complete && self.basis.iter().all(|&b| b == PceBasis::StandardNormal) {
                    "Gaussian Hermite surrogate with held-out physical diagnostics; no error or confidence bound"
                } else if complete {
                    "Hermite-Legendre product surrogate with held-out physical diagnostics; no error or confidence bound"
                } else { "incomplete or refused PCE; no surrogate estimate" }.into(),
                dispersion: f64::INFINITY,
            }, status: self.status, rejection_reason: self.failure.clone(),
        }
    }

    /// Complete at most this many additional physical calls. Zero is a no-op;
    /// complete and refused states are terminal. No observations are replaced.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F) -> PceReport
    where
        F: FnMut(&[f64]) -> Result<f64, E>,
        E: Display,
        C: FnMut() -> bool,
    {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }

    /// `Ok(None)` pauses unfinished work and retries the identical physical input
    /// on resume without consuming a sample ordinal. The solver owns interrupted
    /// costs. Errors and nonfinite values are terminal; callback panics are not
    /// caught and side effects must tolerate retry of unfinished calls.
    pub fn advance_interruptible<F, E, C>(
        &mut self,
        allowance: usize,
        mut cancelled: C,
        mut evaluate: F,
    ) -> PceReport
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
            let parameters = match self.physical_parameters(&self.germs[self.values.len()]) {
                Ok(parameters) => parameters,
                Err(reason) => return self.fail(reason.into()),
            };
            let result = evaluate(&parameters);
            if matches!(&result, Ok(None)) {
                self.status = UqStatus::Cancelled;
                return self.report();
            }
            self.attempted += 1;
            match result {
                Ok(Some(value)) if value.is_finite() => self.values.push(value),
                Ok(Some(_)) => {
                    return self.fail("PCE physical model returned a non-finite QoI".into());
                }
                Ok(None) => unreachable!("interruption did not consume the ordinal"),
                Err(error) => {
                    return self.fail(format!(
                        "PCE physical evaluation {} refused: {error}",
                        self.values.len()
                    ));
                }
            }
        }
        if self.values.len() != self.plan.budget_max_samples {
            self.status = UqStatus::BudgetTruncated;
            return self.report();
        }
        if cancelled() {
            self.status = UqStatus::Cancelled;
            return self.report();
        }
        if self.model.is_none() {
            if let Err(reason) = self.finish_fit() {
                return self.fail(reason.into());
            }
        }
        self.status = if cancelled() {
            UqStatus::Cancelled
        } else {
            UqStatus::Complete
        };
        self.report()
    }

    fn fail(&mut self, reason: String) -> PceReport {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }

    fn finish_fit(&mut self) -> Result<(), &'static str> {
        let n = self.config.training_samples;
        let response_scale = self.values[..n].iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        let rhs: Vec<_> = self.values[..n]
            .iter()
            .map(|&y| {
                if response_scale == 0.0 {
                    0.0
                } else {
                    y / response_scale
                }
            })
            .collect();
        let coefficients: Vec<_> = self
            .regression
            .solve_ls(&rhs)
            .into_iter()
            .map(|c| c / self.design_scale * response_scale)
            .collect();
        if coefficients.iter().any(|c| !c.is_finite()) {
            return Err("PCE fitted coefficients are not representable");
        }
        let model = ProductPceModel::new(
            self.basis.clone(),
            self.indices.clone(),
            coefficients,
        )?;
        let (training_rmse, _) = residuals(&model, &self.germs[..n], &self.values[..n])?;
        let (validation_rmse, validation_mean_error) =
            residuals(&model, &self.germs[n..], &self.values[n..])?;
        let surrogate_std_dev = scaled_norm(&model.coefficients()[1..], 1.0);
        let surrogate_variance = surrogate_std_dev.and_then(|sd| {
            let variance = sd * sd;
            (variance.is_finite() && (variance > 0.0 || sd == 0.0)).then_some(variance)
        });
        let compliance = self.plan.compliance_threshold.and_then(|threshold| {
            let indicators: Vec<_> = self.values[n..]
                .iter()
                .map(|&value| f64::from(u8::from(value <= threshold)))
                .collect();
            summarize_replicates(&indicators)
        });
        self.summary = Some(PceSummary {
            surrogate_mean: model.mean(),
            surrogate_variance,
            surrogate_std_dev,
            training_rmse,
            validation_rmse,
            validation_mean_error,
            probability_of_compliance: compliance.as_ref().map(|c| c.mean),
            compliance_standard_error: compliance.and_then(|c| c.standard_error),
        });
        self.model = Some(model);
        Ok(())
    }
}

fn residuals(
    model: &ProductPceModel,
    germs: &[Vec<f64>],
    values: &[f64],
) -> Result<(f64, f64), &'static str> {
    let mut errors = Vec::with_capacity(values.len());
    for (germ, &value) in germs.iter().zip(values) {
        let error = value - model.try_eval(germ)?;
        if !error.is_finite() {
            return Err("PCE prediction or physical residual is not representable");
        }
        errors.push(error);
    }
    let rmse =
        scaled_norm(&errors, errors.len() as f64).ok_or("PCE residual RMS is not representable")?;
    let mean = summarize_replicates(&errors)
        .ok_or("PCE residual samples are missing")?
        .mean;
    Ok((rmse, mean))
}

fn scaled_norm(values: &[f64], divisor: f64) -> Option<f64> {
    let scale = values.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let result = if scale == 0.0 {
        0.0
    } else {
        fs_math::det::sqrt(values.iter().map(|x| (x / scale).powi(2)).sum::<f64>() / divisor)
            * scale
    };
    result.is_finite().then_some(result)
}
