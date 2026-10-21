//! One bounded physical-mean study: fit, freeze, independently correct.
use core::fmt::Display;
use fs_evidence::Color;
use super::{PceControlError, PceMeanControl, PceMeanEstimate, admit_replay};
use crate::{PceConfig, PceExecution, PceReport, PropagationMethod, UqExecution, UqPlan, UqResult, UqStatus};

const PRODUCTION_KEY_XOR: u64 = 0x5043_454d_4541_4e31;

/// Choices fixed before either phase evaluates physics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PceMeanConfig {
    /// Fitting plus held-out calls, taken from the ORIGINAL lifetime budget.
    pub fit_evaluations: usize,
    /// Training prefix, degree and numerical conditioning policy.
    pub regression: PceConfig,
    /// Fixed control coefficient; one subtracts the centered fitted polynomial.
    pub control_weight: f64,
}
/// Current phase, independent of cancellation or terminal numerical refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PceMeanStage {
    /// Fitting/validation calls or completing their retained QR regression.
    Fitting,
    /// Independent physical production, with its frozen nonlinear control.
    Production,
    /// Physics is paid; only bounded control replay/statistics remain.
    Assessment,
    /// The full fixed study, including its controlled-mean assessment, completed.
    Complete,
}
/// Physical-call accounting and separate surrogate/physics results.
#[derive(Debug, Clone, PartialEq)]
pub struct PceMeanReport {
    /// Original physical quantity of interest.
    pub qoi_name: String,
    /// Total allowance including all fitting, validation and production calls.
    pub model_calls_planned: usize,
    /// Completed or terminally rejected calls across both original executors.
    pub model_calls_evaluated: usize,
    /// Fitting key (also the original study key).
    pub fit_seed: u64,
    /// Independently addressed production key, original key XOR the fixed tag.
    pub production_seed: u64,
    /// Surrogate diagnostics only; never the physical-mean estimate.
    pub fit: PceReport,
    /// Unadjusted production diagnostics, including raw compliance when declared.
    /// Absent before the control is frozen. May remain valid after control refusal.
    pub production: Option<UqResult>,
    /// Available only when the entire fixed experiment and assessment complete.
    pub estimate: Option<PceMeanEstimate>,
    /// Mean-control stage, retained on pause/refusal for diagnosis.
    pub stage: PceMeanStage,
    /// Complete, paused or terminally refused. No statistical early stopping.
    pub status: UqStatus,
    /// Original fit, physical-model or control-assessment diagnostic.
    pub rejection_reason: Option<String>,
    /// Always Estimated; dispersion describes production sampling error only.
    pub evidence_color: Color,
}

/// Owns the existing PCE and Monte Carlo executors under one fixed call budget.
///
/// No fitting or held-out observation enters the production mean. Both phases
/// use the caller's fixed model; changing its meaning on resumption is invalid.
/// Stream keys differ within this run, under the existing pseudorandom model;
/// cross-study callers must avoid key reuse. Cloning retains paid work, not an
/// independent experiment. No adaptive counts/weights or durable codec is added.
#[derive(Debug, Clone)]
pub struct PceMeanExecution {
    plan: UqPlan,
    config: PceMeanConfig,
    fit: PceExecution,
    production: UqExecution,
    control: Option<PceMeanControl>,
    estimate: Option<PceMeanEstimate>,
    stage: PceMeanStage,
    status: UqStatus,
    failure: Option<String>,
}
impl PceMeanExecution {
    /// Admit a MonteCarlo mean-study plan whose sample budget covers BOTH phases.
    /// PCE's training/holdout restrictions apply; at least two production calls
    /// must remain. The production replay work limit is checked before physics.
    pub fn new(plan: &UqPlan, config: PceMeanConfig) -> Result<Self, &'static str> {
        if plan.method != PropagationMethod::MonteCarlo || plan.budget_max_samples > 1_000_000 {
            return Err("PCE-controlled mean needs a MonteCarlo plan with at most one million total calls");
        }
        if !config.control_weight.is_finite() { return Err("PCE control weight must be finite"); }
        let production_count = plan.budget_max_samples.checked_sub(config.fit_evaluations)
            .filter(|&n| n >= 2).ok_or("reserve at least two production calls after the PCE fit")?;
        let mut fit_plan = plan.clone();
        fit_plan.method = PropagationMethod::PolynomialChaos;
        fit_plan.budget_max_samples = config.fit_evaluations;
        let fit = PceExecution::new(&fit_plan, config.regression)?;
        admit_replay(production_count, fit.report().basis_terms, plan.parameters.len())
            .map_err(|_| "PCE production control replay exceeds the work envelope")?;
        let mut production_plan = plan.clone();
        production_plan.budget_max_samples = production_count;
        production_plan.seed ^= PRODUCTION_KEY_XOR;
        let production = UqExecution::new(&production_plan)?;
        Ok(Self {
            plan: plan.clone(), config, fit, production, control: None, estimate: None,
            stage: PceMeanStage::Fitting, status: UqStatus::BudgetTruncated, failure: None,
        })
    }
    /// Immutable original plan, whose budget includes fitting and production.
    #[must_use]
    pub const fn plan(&self) -> &UqPlan { &self.plan }
    /// Immutable regression/allocation/control choices.
    #[must_use]
    pub const fn config(&self) -> PceMeanConfig { self.config }
    /// Read-only fitting owner, retaining all training and held-out observations.
    #[must_use]
    pub const fn fit(&self) -> &PceExecution { &self.fit }
    /// Read-only raw production owner; never includes any fitting observations.
    #[must_use]
    pub const fn production(&self) -> &UqExecution { &self.production }
    /// Frozen control, available after a successful fit/phase transition.
    #[must_use]
    pub fn control(&self) -> Option<&PceMeanControl> { self.control.as_ref() }
    /// Paid/terminally rejected calls, excluding unfinished callback attempts.
    #[must_use]
    pub fn evaluations_attempted(&self) -> usize {
        self.fit.evaluations_attempted() + self.production.evaluations_attempted()
    }
    /// No callbacks or control replay. Raw summaries use their existing owner.
    #[must_use]
    pub fn report(&self) -> PceMeanReport {
        let estimate = if self.status == UqStatus::Complete { self.estimate.clone() } else { None };
        let dispersion = estimate.as_ref().and_then(|e| e.standard_error).unwrap_or(f64::INFINITY);
        PceMeanReport {
            qoi_name: self.plan.target_qoi.clone(), model_calls_planned: self.plan.budget_max_samples,
            model_calls_evaluated: self.evaluations_attempted(), fit_seed: self.fit.plan().seed,
            production_seed: self.production.plan().seed, fit: self.fit.report(),
            production: self.control.as_ref().map(|_| self.production.report()), estimate,
            stage: self.stage, status: self.status, rejection_reason: self.failure.clone(),
            evidence_color: Color::Estimated {
                estimator: "independent PCE-controlled physical mean; no confidence or model-error bound".into(),
                dispersion,
            },
        }
    }
    /// Spend at most this many additional physical calls across both phases.
    /// Zero allowance and completed/refused executions are no-ops.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F) -> PceMeanReport
    where F: FnMut(&[f64]) -> Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }
    /// Unfinished `Ok(None)` calls retry identical inputs without consuming a
    /// paid ordinal. Fit-boundary or final-assessment cancellation never repeats
    /// completed physics. A positive allowance can finish pure computation even
    /// when every physical call is already paid. Long callbacks own cancellation.
    /// Errors/nonfinite values permanently refuse; no failed sample is replaced.
    pub fn advance_interruptible<F, E, C>(&mut self, allowance: usize, mut cancelled: C, mut evaluate: F)
        -> PceMeanReport
    where F: FnMut(&[f64]) -> Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        if allowance == 0 || matches!(self.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let mut remaining = allowance.min(self.plan.budget_max_samples - self.evaluations_attempted());
        if self.control.is_none() {
            if cancelled() { return self.pause(); }
            let before = self.fit.evaluations_attempted();
            let fit = self.fit.advance_interruptible(remaining, &mut cancelled, &mut evaluate);
            remaining -= self.fit.evaluations_attempted() - before;
            if fit.status == UqStatus::Refused {
                return self.fail(fit.rejection_reason.unwrap_or_else(|| "PCE fit refused".into()));
            }
            if fit.status != UqStatus::Complete { self.status = fit.status; return self.report(); }
            if cancelled() { return self.pause(); }
            match PceMeanControl::new(&self.production, &self.fit, self.config.control_weight) {
                Ok(control) => self.control = Some(control),
                Err(error) => return self.fail(error.to_string()),
            }
            self.stage = PceMeanStage::Production;
        }
        if remaining > 0 {
            let production = self.production.advance_interruptible(remaining, &mut cancelled, &mut evaluate);
            if production.status == UqStatus::Refused {
                return self.fail(production.rejection_reason.unwrap_or_else(|| "production physics refused".into()));
            }
            if production.status == UqStatus::Cancelled { return self.pause(); }
        }
        if self.production.evaluations_attempted() != self.production.plan().budget_max_samples {
            self.status = UqStatus::BudgetTruncated;
            return self.report();
        }
        self.stage = PceMeanStage::Assessment;
        let Some(control) = self.control.as_ref() else {
            return self.fail("missing frozen PCE control".into());
        };
        match control.assess_interruptible(&self.production, &mut cancelled) {
            Ok(Some(estimate)) => self.estimate = Some(estimate),
            Ok(None) => return self.fail("completed production has no observations".into()),
            Err(PceControlError::Cancelled) => return self.pause(),
            Err(error) => return self.fail(error.to_string()),
        }
        self.status = UqStatus::Complete;
        self.stage = PceMeanStage::Complete;
        self.report()
    }
    fn pause(&mut self) -> PceMeanReport {
        self.status = UqStatus::Cancelled;
        self.report()
    }
    fn fail(&mut self, reason: String) -> PceMeanReport {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.estimate = None;
        self.report()
    }
}
