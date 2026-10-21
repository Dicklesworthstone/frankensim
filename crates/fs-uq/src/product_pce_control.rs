//! Frozen nonlinear mean controls learned by PCE, assessed on separate physics.
//!
//! For a polynomial P fitted independently of production, E[Y - b(P-E[P])]
//! equals E[Y] under the declared ideal germ law. A fixed approximation changes
//! variance, not that target. Finite PRNG grids, arithmetic and physical-model
//! errors are not bounded here. No confidence sequence or accuracy certificate.
//! Bead: `frankensim-extreal-program-f85xj.6.7`; plan sections 8.8 and 9.7.

mod execution;
pub use execution::{PceMeanConfig, PceMeanExecution, PceMeanReport, PceMeanStage};

use crate::product_execution::sample_parameters;
use crate::product_qmc::summarize_replicates;
use crate::{CorrelationModel, PceBasis, PceExecution, ProductPceModel, UqExecution, UqPlan, UqStatus, UncertaintyKind};

/// Invalid binding, bounded-work refusal or interrupted assessment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PceControlError {
    /// The original fit must have completed, including its held-out calls.
    IncompleteFit,
    /// Physical names, units, marginals, dependence and QoI must match the fit.
    LawMismatch,
    /// Production cannot reuse the fit's Philox key.
    ReusedSeed,
    /// The coefficient must be finite and fixed before production.
    InvalidWeight,
    /// Controls cannot be installed after any production call.
    AlreadySampled,
    /// Assessment must use the complete original production plan.
    PlanMismatch,
    /// Failed physics cannot yield a selected-prefix controlled estimate.
    RefusedExecution,
    /// Declared replay work or storage cannot be admitted.
    Budget,
    /// Polynomial, adjusted response or summary arithmetic is unrepresentable.
    NumericalRange,
    /// Cancellation publishes no partial assessment and changes no observations.
    Cancelled,
}
impl std::fmt::Display for PceControlError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::IncompleteFit => "PCE control requires a completed fit",
            Self::LawMismatch => "PCE control and production need identical physical laws and QoI",
            Self::ReusedSeed => "PCE production must use a different Philox key from fitting",
            Self::InvalidWeight => "PCE control weight must be finite and frozen",
            Self::AlreadySampled => "freeze PCE control before production observations",
            Self::PlanMismatch => "PCE control belongs to a different production plan",
            Self::RefusedExecution => "refused physics cannot publish a controlled mean",
            Self::Budget => "PCE control replay exceeds its admitted work or storage",
            Self::NumericalRange => "PCE controlled-mean arithmetic is not representable",
            Self::Cancelled => "PCE control assessment cancelled",
        })
    }
}
impl std::error::Error for PceControlError {}

/// Mean-only comparison on the same physical observations. The adjusted
/// observations are NOT a physical distribution: do not use them for quantiles,
/// CVaR, extrema or threshold indicators. Sampling errors are descriptive only.
#[derive(Debug, Clone, PartialEq)]
pub struct PceMeanEstimate {
    /// All accepted production observations, in their original order.
    pub samples: usize,
    /// Ordinary mean of the physical QoI.
    pub raw_mean: f64,
    /// Mean of Y - weight * (P(germ) - E[P]); not the surrogate expectation.
    pub mean: f64,
    /// Ordinary descriptive standard error; unavailable with one observation.
    pub raw_standard_error: Option<f64>,
    /// Adjusted descriptive standard error, never an optional-stopping bound.
    pub standard_error: Option<f64>,
    /// Adjusted/raw sample variance; values above one are deliberately retained.
    /// Unavailable for a zero denominator or unrepresentable ratio.
    pub variance_ratio: Option<f64>,
}

/// Immutable fitted control bound to one complete Monte Carlo production plan.
///
/// Construction requires a completed fit, an empty production execution, and
/// different Philox keys. These checks prevent accidental within-run reuse;
/// they cannot prove the user did not select a model/weight using other views
/// of production data. Fix the control without those outcomes. Numerical
/// conditioning, model quality and surrogate diagnostics are not bias bounds.
#[derive(Debug, Clone)]
pub struct PceMeanControl {
    plan: UqPlan,
    germ_plan: UqPlan,
    centered: ProductPceModel,
    weight: f64,
    fit_seed: u64,
}
impl PceMeanControl {
    /// Freeze a completed Gaussian/uniform product surrogate before production.
    /// Joint Gaussian laws, including singular factors, need no inverse map:
    /// assessment regenerates the ORIGINAL independent germs by sample ordinal.
    /// Replay is capped at 67,108,864 term-coordinate products per assessment,
    /// using the full production allowance, before any production work is done.
    pub fn new(production: &UqExecution, fit: &PceExecution, weight: f64) -> Result<Self, PceControlError> {
        if production.status == UqStatus::Refused || production.failure.is_some() {
            return Err(PceControlError::RefusedExecution);
        }
        if production.evaluations_attempted() != 0 || !production.observations().is_empty() {
            return Err(PceControlError::AlreadySampled);
        }
        let model = fit.model().ok_or(PceControlError::IncompleteFit)?;
        let plan = production.plan();
        if plan.parameters != fit.plan().parameters || plan.correlation != fit.plan().correlation
            || plan.target_qoi != fit.plan().target_qoi
        {
            return Err(PceControlError::LawMismatch);
        }
        if plan.seed == fit.plan().seed { return Err(PceControlError::ReusedSeed); }
        if !weight.is_finite() { return Err(PceControlError::InvalidWeight); }
        admit_replay(plan.budget_max_samples, model.indices().len(), model.dimension())?;
        // Strip the constant COEFFICIENT, not two large evaluated values. This
        // keeps small nonconstant controls next to a large physical offset.
        let mut coefficients = model.coefficients().to_vec();
        coefficients[0] = 0.0;
        let centered = ProductPceModel::new(model.basis().to_vec(), model.indices().to_vec(), coefficients)
            .map_err(|_| PceControlError::NumericalRange)?;
        let mut germ_plan = plan.clone();
        germ_plan.correlation = CorrelationModel::Independent;
        for (parameter, family) in germ_plan.parameters.iter_mut().zip(model.basis()) {
            parameter.kind = match family {
                PceBasis::StandardNormal => UncertaintyKind::AleatoryGaussian { mean: 0.0, std_dev: 1.0 },
                PceBasis::Uniform01 => UncertaintyKind::AleatoryUniform { lo: 0.0, hi: 1.0 },
            };
        }
        Ok(Self { plan: plan.clone(), germ_plan, centered, weight, fit_seed: fit.plan().seed })
    }

    /// The entire immutable production plan, including its lifetime allowance.
    #[must_use]
    pub const fn plan(&self) -> &UqPlan { &self.plan }
    /// Key used by the source fit; never reused by this production execution.
    #[must_use]
    pub const fn fit_seed(&self) -> u64 { self.fit_seed }
    /// Frozen coefficient. One uses the fitted polynomial without reweighting.
    #[must_use]
    pub const fn weight(&self) -> f64 { self.weight }

    /// Assess retained production responses without any new physics evaluations.
    pub fn assess(&self, production: &UqExecution) -> Result<Option<PceMeanEstimate>, PceControlError> {
        self.assess_interruptible(production, || false)
    }

    /// Replay only independent germs and the cheap centered polynomial. Raw
    /// response bits, compliance statistics and checkpoint bytes are unchanged.
    /// Polls before each replayed sample and around the shared synchronous moment
    /// reductions. Partial-run errors are descriptive, not valid stopping rules.
    /// Refusal or cancellation returns no partial assessment; retry is read-only.
    pub fn assess_interruptible(&self, production: &UqExecution, mut cancelled: impl FnMut() -> bool)
        -> Result<Option<PceMeanEstimate>, PceControlError>
    {
        if production.plan() != &self.plan { return Err(PceControlError::PlanMismatch); }
        if production.status == UqStatus::Refused || production.failure.is_some() {
            return Err(PceControlError::RefusedExecution);
        }
        poll(&mut cancelled)?;
        let values = production.observations();
        if values.is_empty() { return Ok(None); }
        let mut adjusted = Vec::new();
        adjusted.try_reserve_exact(values.len()).map_err(|_| PceControlError::Budget)?;
        for (ordinal, &value) in values.iter().enumerate() {
            poll(&mut cancelled)?;
            let corrected = if self.weight == 0.0 { value } else {
                let germ = sample_parameters(&self.germ_plan, None, ordinal)
                    .map_err(|_| PceControlError::NumericalRange)?;
                let control = self.centered.try_eval(&germ).map_err(|_| PceControlError::NumericalRange)?;
                (-self.weight).mul_add(control, value)
            };
            if !corrected.is_finite() { return Err(PceControlError::NumericalRange); }
            adjusted.push(corrected);
        }
        poll(&mut cancelled)?;
        let raw = summarize_replicates(values).ok_or(PceControlError::NumericalRange)?;
        poll(&mut cancelled)?;
        let estimate = summarize_replicates(&adjusted).ok_or(PceControlError::NumericalRange)?;
        let variance_ratio = raw.standard_error.zip(estimate.standard_error).and_then(|(raw, controlled)| {
            if raw == 0.0 { return None; }
            let ratio = (controlled / raw).powi(2);
            (ratio.is_finite() && (ratio > 0.0 || controlled == 0.0)).then_some(ratio)
        });
        poll(&mut cancelled)?;
        Ok(Some(PceMeanEstimate {
            samples: values.len(), raw_mean: raw.mean, mean: estimate.mean,
            raw_standard_error: raw.standard_error, standard_error: estimate.standard_error, variance_ratio,
        }))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PceControlError> {
    if cancelled() { Err(PceControlError::Cancelled) } else { Ok(()) }
}
fn admit_replay(samples: usize, terms: usize, dimension: usize) -> Result<(), PceControlError> {
    if samples.checked_mul(terms).and_then(|n| n.checked_mul(dimension))
        .is_none_or(|n| n > 67_108_864)
    { Err(PceControlError::Budget) } else { Ok(()) }
}
