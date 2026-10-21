//! Fixed-budget adaptive pilots followed by independent, frozen production.
use core::fmt::Display;
use fs_evidence::Color;
use super::{ImportanceFitConfig, ImportanceFitError, ImportanceProposalFit};
use crate::{ImportanceConfig, ImportanceEstimate, ImportanceExecution, ImportanceReport,
    UqPlan, UqStatus};

/// All counts, fitting choices and defensive mass are fixed before physics.
#[derive(Debug, Clone, PartialEq)]
pub struct PilotImportanceConfig {
    /// Complete adaptation rounds in 1..=16. No outcome-dependent early exit.
    pub pilot_rounds: usize,
    /// At least two model calls per round; these never enter production results.
    pub samples_per_pilot: usize,
    pub initial_proposal: ImportanceConfig,
    pub fit: ImportanceFitConfig,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PilotImportanceStage {
    Pilot { round: usize },
    Production,
    Complete,
}

/// Retained fitting diagnostics, not a pooled failure-probability estimate.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportancePilotRound {
    pub seed: u64,
    pub sampled_proposal: ImportanceConfig,
    pub fit: ImportanceProposalFit,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PilotImportanceReport {
    pub qoi_name: String,
    pub model_calls_planned: usize,
    pub model_calls_evaluated: usize,
    pub stage: PilotImportanceStage,
    pub completed_pilots: Vec<ImportancePilotRound>,
    pub current_seed: u64,
    pub current_proposal: ImportanceConfig,
    pub current_evaluations: usize,
    /// Raw report from ONLY the independently keyed, fixed final proposal.
    /// Partial production statistics are descriptive diagnostics, not stopping rules.
    pub production: Option<ImportanceReport>,
    /// Available only after the complete fixed production design; absent on refusal.
    pub estimate: Option<ImportanceEstimate>,
    pub status: UqStatus,
    pub rejection_reason: Option<String>,
    pub evidence_color: Color,
}

/// Uses the existing importance executor for every pilot and production call.
///
/// The original plan allowance covers all rounds plus at least two production
/// samples, with the canonical one-million-call cap. Each phase has a distinct
/// Philox key; production observations never choose their own proposal, sample
/// count, or threshold. Completed pilots are discarded only after their fit and
/// next owner are ready. History retains O(rounds*dimension) diagnostics rather
/// than all old observations. Clone resumes the SAME experiment.
///
/// This estimates a fixed declared model's failure probability, not model
/// validity. A single learned mean may stagnate or miss disconnected modes.
/// No bias/variance improvement, achieved accuracy, or optimal proposal is
/// promised. Defensive mass, unit covariance and all physical laws are retained.
#[derive(Debug, Clone)]
pub struct PilotImportanceExecution {
    plan: UqPlan,
    config: PilotImportanceConfig,
    production_samples: usize,
    completed: Vec<ImportancePilotRound>,
    current: ImportanceExecution,
    status: UqStatus,
    failure: Option<String>,
}

impl PilotImportanceExecution {
    pub fn new(plan: &UqPlan, config: PilotImportanceConfig) -> Result<Self, &'static str> {
        if !(1..=16).contains(&config.pilot_rounds) || config.samples_per_pilot < 2
            || plan.budget_max_samples > 1_000_000
        {
            return Err("importance learning needs 1..=16 pilots, >=2 calls each and <=1000000 total calls");
        }
        let pilot_calls = config.pilot_rounds.checked_mul(config.samples_per_pilot)
            .ok_or("importance pilot call count overflows")?;
        let production_samples = plan.budget_max_samples.checked_sub(pilot_calls)
            .filter(|&n| n >= 2).ok_or("reserve at least two independent production calls")?;
        config.fit.admit(config.samples_per_pilot, plan.parameters.len())
            .map_err(|_| "importance pilot fitting policy or replay work is invalid")?;
        if pilot_calls.checked_mul(plan.parameters.len()).is_none_or(|n| n > 67_108_864) {
            return Err("aggregate importance pilot replay exceeds the work envelope");
        }
        let mut first = plan.clone();
        first.budget_max_samples = config.samples_per_pilot;
        // ImportanceExecution admits the original method, Gaussian law, threshold,
        // proposal and numerical factor before any physical callback.
        let current = ImportanceExecution::new(&first, config.initial_proposal.clone())?;
        let mut completed = Vec::new();
        completed.try_reserve_exact(config.pilot_rounds)
            .map_err(|_| "cannot reserve importance pilot history")?;
        Ok(Self { plan: plan.clone(), config, production_samples, completed, current,
            status: UqStatus::BudgetTruncated, failure: None })
    }

    pub const fn plan(&self) -> &UqPlan { &self.plan }
    pub const fn config(&self) -> &PilotImportanceConfig { &self.config }
    pub fn pilot_history(&self) -> &[ImportancePilotRound] { &self.completed }
    /// Original current phase, for read-only inspection of paid observations.
    pub const fn current(&self) -> &ImportanceExecution { &self.current }
    /// Fixed final producer, absent until all pilot fits have completed. Its
    /// existing sampling-confidence assessment can be applied to production
    /// alone, with alpha chosen in advance and the original numerical/model caveats.
    pub fn production(&self) -> Option<&ImportanceExecution> {
        (self.completed.len() == self.config.pilot_rounds).then_some(&self.current)
    }
    pub fn evaluations_attempted(&self) -> usize {
        self.completed.len() * self.config.samples_per_pilot + self.current.attempted
    }
    pub fn report(&self) -> PilotImportanceReport {
        let production = self.production().map(ImportanceExecution::report);
        let estimate = if self.status == UqStatus::Complete {
            production.as_ref().and_then(|r| r.estimate)
        } else { None };
        let stage = if self.completed.len() < self.config.pilot_rounds {
            PilotImportanceStage::Pilot { round: self.completed.len() }
        } else if self.status == UqStatus::Complete { PilotImportanceStage::Complete }
        else { PilotImportanceStage::Production };
        PilotImportanceReport {
            qoi_name: self.plan.target_qoi.clone(), model_calls_planned: self.plan.budget_max_samples,
            model_calls_evaluated: self.evaluations_attempted(), stage,
            completed_pilots: self.completed.clone(), current_seed: self.current.plan.seed,
            current_proposal: self.current.config.clone(), current_evaluations: self.current.attempted,
            production, estimate, status: self.status, rejection_reason: self.failure.clone(),
            evidence_color: Color::Estimated {
                estimator: "pilot-learned fixed production importance sampling; no model-error or accuracy guarantee".into(),
                dispersion: estimate.and_then(|e| e.standard_error).unwrap_or(f64::INFINITY),
            },
        }
    }

    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut evaluate: F)
        -> PilotImportanceReport
    where F: FnMut(&[f64]) -> Result<f64, E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |x| evaluate(x).map(Some))
    }

    /// Spend a total additional physical-call allowance across phase boundaries.
    /// Zero allowance is a no-op. Positive allowance may finish a paid pilot's
    /// pure fitting/transition or publish completed production with no new call.
    /// Interrupted calls retry identical inputs and weights. Cancellation also
    /// brackets each fit, owner transition and final publication. Completed model
    /// calls are retained, and errors/nonfinite responses permanently refuse.
    pub fn advance_interruptible<F, E, C>(
        &mut self, allowance: usize, mut cancelled: C, mut evaluate: F,
    ) -> PilotImportanceReport
    where F: FnMut(&[f64]) -> Result<Option<f64>, E>, E: Display, C: FnMut() -> bool {
        if allowance == 0 || matches!(self.status, UqStatus::Complete | UqStatus::Refused) {
            return self.report();
        }
        let mut remaining = allowance.min(self.plan.budget_max_samples - self.evaluations_attempted());
        loop {
            if cancelled() { return self.pause(); }
            if self.current.status != UqStatus::Complete {
                if remaining == 0 {
                    self.status = UqStatus::BudgetTruncated;
                    return self.report();
                }
                let before = self.current.attempted;
                let report = self.current.advance_interruptible(remaining, &mut cancelled, &mut evaluate);
                remaining -= self.current.attempted - before;
                match report.status {
                    UqStatus::Refused => return self.fail(report.rejection_reason
                        .unwrap_or_else(|| "importance physical execution refused".into())),
                    UqStatus::Cancelled => return self.pause(),
                    UqStatus::BudgetTruncated => {
                        self.status = UqStatus::BudgetTruncated;
                        return self.report();
                    }
                    UqStatus::Complete => {}
                }
            }
            if cancelled() { return self.pause(); }
            if self.completed.len() == self.config.pilot_rounds {
                self.status = UqStatus::Complete;
                return self.report();
            }
            let fit = match self.current.fit_proposal(self.config.fit, &mut cancelled) {
                Ok(fit) => fit,
                Err(ImportanceFitError::Cancelled) => return self.pause(),
                Err(error) => return self.fail(error.to_string()),
            };
            if cancelled() { return self.pause(); }
            let next_index = self.completed.len() + 1;
            let mut next_plan = self.plan.clone();
            next_plan.seed = phase_seed(self.plan.seed, next_index);
            next_plan.budget_max_samples = if next_index == self.config.pilot_rounds {
                self.production_samples
            } else { self.config.samples_per_pilot };
            let next = match ImportanceExecution::new(&next_plan, fit.proposal.clone()) {
                Ok(next) => next,
                Err(error) => return self.fail(error.into()),
            };
            if cancelled() { return self.pause(); }
            // Commit together: all paid pilot data remain live until both the
            // proposal fit and the next execution owner have succeeded.
            self.completed.push(ImportancePilotRound {
                seed: self.current.plan.seed,
                sampled_proposal: self.current.config.clone(), fit,
            });
            self.current = next;
        }
    }

    fn pause(&mut self) -> PilotImportanceReport {
        self.status = UqStatus::Cancelled;
        self.report()
    }
    fn fail(&mut self, reason: String) -> PilotImportanceReport {
        self.status = UqStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }
}

// Odd multiplier makes the bounded phase indices injective modulo 2^64,
// including wraparound. Distinct keys model independent pseudorandom streams;
// callers remain responsible for avoiding key reuse between separate studies.
fn phase_seed(seed: u64, phase: usize) -> u64 {
    seed.wrapping_add((phase as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15))
}
