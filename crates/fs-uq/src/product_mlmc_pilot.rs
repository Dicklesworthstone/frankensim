//! Pilot-allocated MLMC, reusing the existing coupled physical-call executor.
//!
//! A fixed pilot learns correction dispersion. A disjoint keyed production
//! stream then runs ONE immutable allocation. Pilot observations never enter
//! the production estimator and production observations never change its size.
//! This avoids reusing the data that chose N as though N had been fixed ahead
//! of those observations. The finest declared level is always the target.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.7`.
use core::fmt::Display;
use fs_evidence::Color;
use crate::product_mlmc_allocation::{physical_calls, validate_policy};
use crate::product_plan::admit_plan;
use crate::{allocate_mlmc, MlmcAllocation, MlmcAllocationError, MlmcAllocationPolicy,
    MlmcEstimate, MlmcExecution, MlmcExecutionReport, PropagationMethod, UqPlan, UqStatus};

// A distinct Philox key, not an unrelated sampler or reuse of pilot ordinals
// under the same key. The two keys are recorded in each report for exact replay.
const PRODUCTION_KEY_DOMAIN: u64 = 0x4d4c_4d43_5052_4f44;

/// Fixed pilot, physical cost model, and production accuracy/allocation policy.
#[derive(Debug, Clone, PartialEq)]
pub struct PilotMlmcConfig {
    /// At least two correction samples on EACH declared level before allocation.
    pub pilot_samples_per_level: usize,
    /// Cost of one physical evaluation per resolution, in common arbitrary units.
    pub physical_model_costs: Vec<f64>,
    pub target_standard_error: f64,
    pub minimum_production_samples_per_level: usize,
    /// Positive QoI-unit safeguard used only in allocation, not final statistics.
    pub correction_std_dev_floor: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PilotMlmcStage { Pilot, Production }

/// Complete lifetime accounting and separately identified pilot/production data.
#[derive(Debug, Clone, PartialEq)]
pub struct PilotMlmcReport {
    pub stage: PilotMlmcStage,
    pub status: UqStatus,
    /// Original cap includes ALL pilot and production physical calls.
    pub maximum_model_calls: usize,
    pub model_calls_evaluated: usize,
    /// Exact pilot-plus-production allocation, known only after pilot admission.
    pub model_calls_planned: Option<usize>,
    pub pilot_seed: u64,
    pub production_seed: u64,
    pub target_standard_error: f64,
    /// Diagnostic pilot statistics are NOT the final finest-level estimator.
    pub pilot: MlmcExecutionReport,
    pub allocation: Option<MlmcAllocation>,
    pub production: Option<MlmcExecutionReport>,
    /// Production only, and only when the entire selected allocation completes.
    pub estimate: Option<MlmcEstimate>,
    /// Production physical indicators only. Allocation targets QoI variance,
    /// not compliance variance; no compliance accuracy target is implied.
    pub compliance: Option<MlmcEstimate>,
    /// Descriptive comparison after completion; None when SE is unavailable.
    /// False does not trigger more sampling or get relabeled Complete accuracy.
    pub observed_sampling_target_met: Option<bool>,
    pub evidence_color: Color,
    pub allocation_refusal: Option<MlmcAllocationError>,
    pub rejection_reason: Option<String>,
}

/// One owning, resumable pilot/production experiment on a caller-fixed ladder.
///
/// `plan.budget_max_samples` is a lifetime cap, not a required spend. Admission
/// reserves enough calls for the full pilot and the minimum production design.
/// An infeasible pilot-derived target refuses BEFORE any production callback.
/// Cancellation, chunk boundaries and Clone preserve all accepted corrections,
/// paid fine half-pairs, and the chosen allocation. No failed sample is replaced.
///
/// Both phases use `MlmcExecution`, including its law admission, exact paired
/// parameters, stable statistics and terminal failure semantics. Pilot uses the
/// original seed; production uses seed XOR 0x4d4c4d4350524f44. Statistical stream
/// independence is the existing Philox sampling assumption, not a new proof.
/// Costs must describe the same deterministic model throughout both phases.
/// There is no level/bias adaptation, confidence sequence, durable checkpoint,
/// continuum error certificate, or guarantee that pilot variance predicts future
/// variance. Long model calls must poll their own cancellation scope.
#[derive(Debug, Clone)]
pub struct PilotMlmcExecution {
    plan: UqPlan,
    config: PilotMlmcConfig,
    policy: MlmcAllocationPolicy,
    pilot: MlmcExecution,
    production: Option<MlmcExecution>,
    allocation: Option<MlmcAllocation>,
    status: UqStatus,
    allocation_refusal: Option<MlmcAllocationError>,
    failure: Option<String>,
}

impl PilotMlmcExecution {
    /// Validate the lifetime cap, fixed pilot, minimum production design, and
    /// complete probability law before any physical evaluation. Construction
    /// and <=64-level allocation reductions are bounded synchronous operations.
    pub fn new(plan: &UqPlan, config: PilotMlmcConfig) -> Result<Self, MlmcAllocationError> {
        if plan.method != PropagationMethod::MultilevelMonteCarlo {
            return Err(MlmcAllocationError::Invalid("pilot MLMC requires a MultilevelMonteCarlo plan"));
        }
        let dim = config.physical_model_costs.len();
        if dim == 0 || dim > 64 || config.pilot_samples_per_level < 2 {
            return Err(MlmcAllocationError::Invalid("pilot MLMC needs 1..=64 levels and at least two pilot corrections per level"));
        }
        let mut admission = plan.clone();
        admission.method = PropagationMethod::MonteCarlo;
        admit_plan(&admission).map_err(MlmcAllocationError::Invalid)?;
        let pilot_counts = vec![config.pilot_samples_per_level; dim];
        let pilot_calls = physical_calls(&pilot_counts).ok_or(MlmcAllocationError::Invalid("pilot physical-call count overflows"))?;
        let remaining = plan.budget_max_samples.checked_sub(pilot_calls)
            .ok_or(MlmcAllocationError::BudgetExceeded { maximum_model_calls: plan.budget_max_samples })?;
        let policy = MlmcAllocationPolicy {
            target_standard_error: config.target_standard_error,
            minimum_samples_per_level: config.minimum_production_samples_per_level,
            correction_std_dev_floor: config.correction_std_dev_floor,
            maximum_model_calls: remaining,
        };
        validate_policy(&config.physical_model_costs, policy)?;
        let mut pilot_plan = plan.clone();
        pilot_plan.budget_max_samples = pilot_calls;
        let pilot = MlmcExecution::new(&pilot_plan, &pilot_counts).map_err(MlmcAllocationError::Invalid)?;
        Ok(Self { plan: plan.clone(), config, policy, pilot, production: None,
            allocation: None, status: UqStatus::BudgetTruncated, allocation_refusal: None, failure: None })
    }

    #[must_use]
    pub const fn plan(&self) -> &UqPlan { &self.plan }
    #[must_use]
    pub const fn config(&self) -> &PilotMlmcConfig { &self.config }
    /// Read-only access prevents bypassing phase/lifetime accounting.
    #[must_use]
    pub const fn pilot(&self) -> &MlmcExecution { &self.pilot }
    #[must_use]
    pub fn production(&self) -> Option<&MlmcExecution> { self.production.as_ref() }
    #[must_use]
    pub fn evaluations_attempted(&self) -> usize {
        self.pilot.evaluations_attempted() + self.production.as_ref().map_or(0, MlmcExecution::evaluations_attempted)
    }

    #[must_use]
    pub fn report(&self) -> PilotMlmcReport {
        let pilot = self.pilot.report();
        let production = self.production.as_ref().map(MlmcExecution::report);
        let complete = self.status == UqStatus::Complete;
        let estimate = if complete { production.as_ref().and_then(|r| r.estimate.clone()) } else { None };
        let compliance = if complete { production.as_ref().and_then(|r| r.compliance.clone()) } else { None };
        let observed_sampling_target_met = estimate.as_ref().and_then(|e| e.standard_error)
            .map(|se| se <= self.config.target_standard_error);
        let dispersion = estimate.as_ref().and_then(|e| e.standard_error).unwrap_or(f64::INFINITY);
        PilotMlmcReport {
            stage: if self.production.is_some() { PilotMlmcStage::Production } else { PilotMlmcStage::Pilot },
            status: self.status, maximum_model_calls: self.plan.budget_max_samples,
            model_calls_evaluated: self.evaluations_attempted(),
            model_calls_planned: self.allocation.as_ref().map(|a| pilot.model_calls_planned + a.model_calls),
            pilot_seed: self.plan.seed, production_seed: self.plan.seed ^ PRODUCTION_KEY_DOMAIN,
            target_standard_error: self.config.target_standard_error,
            pilot, allocation: self.allocation.clone(), production, estimate, compliance,
            observed_sampling_target_met,
            evidence_color: Color::Estimated {
                estimator: "pilot-allocated fixed-level MLMC; disjoint production keys; descriptive sampling error, no bias or confidence bound".into(),
                dispersion,
            },
            allocation_refusal: self.allocation_refusal.clone(), rejection_reason: self.failure.clone(),
        }
    }

    /// Spend at most `allowance` additional physical calls across BOTH phases.
    /// Zero allowance and terminal runs are no-ops. Completing the pilot can
    /// admit the production allocation without spending a production call.
    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F) -> PilotMlmcReport
    where F: FnMut(usize, &[f64]) -> Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |level, x| evaluate(level, x).map(Some))
    }

    /// `Ok(None)` retains paid work and retries the same unfinished call.
    /// Errors/nonfinite responses are terminal in either phase. Cancellation
    /// brackets the phase transition as well as each physical evaluation.
    pub fn advance_interruptible<F, E, C>(
        &mut self, allowance: usize, mut cancelled: C, mut evaluate: F,
    ) -> PilotMlmcReport
    where F: FnMut(usize, &[f64]) -> Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        if allowance == 0 || matches!(self.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let mut remaining = allowance.min(self.plan.budget_max_samples - self.evaluations_attempted());
        if self.production.is_none() {
            let before = self.pilot.evaluations_attempted();
            let report = self.pilot.advance_interruptible(remaining, &mut cancelled, &mut evaluate);
            remaining -= self.pilot.evaluations_attempted() - before;
            if report.status != UqStatus::Complete {
                self.status = report.status;
                self.failure = report.rejection_reason;
                return self.report();
            }
            // The pilot may already be complete after a cancelled transition;
            // its original executor then performs no repeated physical calls.
            if cancelled() { self.status = UqStatus::Cancelled; return self.report(); }
            if let Err(error) = self.prepare_production(&report) {
                self.failure = Some(error.to_string());
                self.allocation_refusal = Some(error);
                self.status = UqStatus::Refused;
                return self.report();
            }
            if cancelled() { self.status = UqStatus::Cancelled; return self.report(); }
        }
        if remaining == 0 { self.status = UqStatus::BudgetTruncated; return self.report(); }
        let report = self.production.as_mut().expect("successful admission installed production")
            .advance_interruptible(remaining, cancelled, evaluate);
        self.status = report.status;
        self.failure = report.rejection_reason;
        self.report()
    }

    fn prepare_production(&mut self, pilot: &MlmcExecutionReport) -> Result<(), MlmcAllocationError> {
        let mut deviations = Vec::with_capacity(pilot.levels.len());
        for level in &pilot.levels {
            let se = level.correction.as_ref().and_then(|c| c.standard_error)
                .ok_or(MlmcAllocationError::Numerical("pilot correction dispersion is unavailable"))?;
            let sd = se * fs_math::det::sqrt(level.samples_accepted as f64);
            if !sd.is_finite() {
                return Err(MlmcAllocationError::Numerical("pilot correction standard deviation is not representable"));
            }
            deviations.push(sd);
        }
        let allocation = allocate_mlmc(&deviations, &self.config.physical_model_costs, self.policy)?;
        let mut plan = self.plan.clone();
        plan.budget_max_samples = allocation.model_calls;
        plan.seed ^= PRODUCTION_KEY_DOMAIN;
        let production = MlmcExecution::new(&plan, &allocation.samples_per_level)
            .map_err(MlmcAllocationError::Invalid)?;
        // Publish both together; no partially installed allocation can execute.
        self.production = Some(production);
        self.allocation = Some(allocation);
        Ok(())
    }
}
