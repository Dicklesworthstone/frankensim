//! ANYTIME + REFUSAL SEMANTICS (addendum Proposal 8, bead lmp4.17;
//! ships behind `ladder-planner` because it drives the planner, but the
//! CONTRACT here survives even if the planner is frozen by its kill
//! criterion — this is the product win): every query whose first budget can
//! fund a solve returns IMMEDIATELY with a wide certified interval and tightens
//! monotonically as budget is spent; every interval carries its color
//! and a "WHAT WOULD TIGHTEN THIS" hint that prices the next move; and
//! when the budget cannot discharge the query, the system SAYS SO with
//! the interval it DID achieve and the price of the gap — never a
//! silent best-effort number dressed as an answer.
//!
//! Determinism (G5): the planner underneath is deterministic, so a
//! replayed query reproduces the same interval trajectory.

use crate::planner::{
    AnswerCache, CostTable, PlanControl, PlanError, PlanObserver, PlanOp, PlanOutcome,
    ProblemFamily, VerifierCertificate, plan_observed, validate_bound, validate_budget,
    validate_finite, validate_positive_finite, validate_rung_cells,
};
use fs_evidence::Color;

/// Maximum entries in one operational budget ladder.
pub const MAX_BUDGET_RUNGS: usize = 4_096;

/// Whether the caller permits the planner to execute beyond an emitted step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnytimeControl {
    /// Continue toward later budget rungs.
    Continue,
    /// Return the emitted deterministic prefix without later side effects.
    Stop,
}

/// One point on the anytime trajectory.
#[derive(Debug, Clone, PartialEq)]
pub struct IntervalStep {
    /// The budget this step was allowed (cells).
    pub budget: f64,
    /// Cumulative solved-cell work completed when this step was emitted.
    pub spent: f64,
    /// The certified half-width achieved.
    pub bound: f64,
    /// The Proposal-3 color of the interval (equilibrated enclosures
    /// are VERIFIED; the operator always knows what they hold).
    pub color: Color,
    /// Stable family of the verifier that minted `color`.
    pub verifier_family: &'static str,
    /// Reconstructed-flux identity bound to the verifier certificate.
    pub flux_hash: u64,
    /// The "what would tighten this" hint, priced.
    pub hint: String,
    /// True when the query discharged at this budget.
    pub discharged: bool,
}

/// The anytime result: the trajectory plus the final verdict.
#[derive(Debug, Clone, PartialEq)]
pub struct AnytimeReport {
    /// The interval trajectory, one entry per budget rung that produced a
    /// certificate. Unfunded rungs cannot honestly contribute a colored step.
    pub trajectory: Vec<IntervalStep>,
    /// The refusal note when the final budget could not discharge —
    /// the achieved interval AND the price of the gap, teaching.
    pub refusal: Option<String>,
    /// True when the operational observer stopped the cumulative execution.
    pub stopped_by_observer: bool,
}

impl AnytimeReport {
    /// The final certified bound.
    #[must_use]
    pub fn final_bound(&self) -> f64 {
        self.trajectory.last().map_or(f64::INFINITY, |s| s.bound)
    }

    /// Did the query discharge within the final budget?
    #[must_use]
    pub fn discharged(&self) -> bool {
        self.trajectory.last().is_some_and(|s| s.discharged)
    }

    /// Did the caller stop the operational execution at an emitted rung?
    #[must_use]
    pub fn stopped_by_observer(&self) -> bool {
        self.stopped_by_observer
    }
}

/// The "what would tighten this" hint: extrapolate the price of closing
/// the gap from the achieved bound and spend (O(h) energy convergence:
/// cells scale like bound/tol), and NAME where the money goes (the
/// residual's hot region when the planner exposes it, else the
/// operator menu's next move). Cold telemetry degrades to a generic
/// but still-priced hint.
///
/// # Errors
/// Returns [`PlanError`] for malformed inputs or non-finite extrapolation.
pub fn tighten_hint(
    bound: f64,
    tol: f64,
    cells_spent: f64,
    costs: &CostTable,
    hot_region: Option<(f64, f64)>,
) -> Result<String, PlanError> {
    validate_bound(bound, "achieved_bound")?;
    validate_positive_finite(tol, "tolerance")?;
    validate_finite(cells_spent, "cells_spent")?;
    if cells_spent < 0.0 {
        return Err(PlanError::InvalidScalar {
            field: "cells_spent",
            reason: "must be non-negative",
        });
    }
    if let Some((lo, hi)) = hot_region
        && (!lo.is_finite() || !hi.is_finite() || lo < 0.0 || lo >= hi || hi > 1.0)
    {
        return Err(PlanError::InvalidScalar {
            field: "hot_region",
            reason: "must be a finite ordered subinterval of [0,1]",
        });
    }
    if bound <= tol {
        return Ok("already at tolerance — spend nothing".to_string());
    }
    let factor = (bound / tol).max(1.0);
    let projected = cells_spent * factor;
    if !factor.is_finite() || !projected.is_finite() {
        return Err(PlanError::NumericalFailure {
            stage: "anytime gap-price extrapolation",
        });
    }
    let extra = (projected - cells_spent).max(1.0);
    let next_op = if costs.predict(PlanOp::DwrRefine) <= costs.predict(PlanOp::Climb) {
        PlanOp::DwrRefine
    } else {
        PlanOp::Climb
    };
    Ok(match hot_region {
        Some((lo, hi)) => format!(
            "closing ±{bound:.3e} to ±{tol:.3e} needs ~{extra:.0} more cells via \
             {}, mostly on the region x ∈ [{lo:.2}, {hi:.2}]",
            next_op.name()
        ),
        None => format!(
            "closing ±{bound:.3e} to ±{tol:.3e} needs ~{extra:.0} more cells via {}",
            next_op.name()
        ),
    })
}

/// Collect an operational anytime execution over an increasing budget ladder.
///
/// This is the convenience wrapper around [`run_anytime_observed`]. Steps are
/// produced while one cumulative planner execution is running, not reconstructed
/// from a completed final-budget run. The wrapper simply continues after every
/// emitted step and returns the collected trajectory.
///
/// # Errors
/// Returns [`PlanError`] for an invalid query, budget/rung ladder, cache
/// replay, cost table, or planner computation.
pub fn run_anytime(
    family: &ProblemFamily,
    theta: f64,
    tol: f64,
    budget_ladder: &[f64],
    rung_cells: &[usize],
    cache: &mut dyn AnswerCache,
    costs: &mut CostTable,
) -> Result<AnytimeReport, PlanError> {
    run_anytime_observed(
        family,
        theta,
        tol,
        budget_ladder,
        rung_cells,
        cache,
        costs,
        |_| AnytimeControl::Continue,
    )
}

/// Run one cumulative planner execution and expose each certified budget rung
/// before work for any later rung begins.
///
/// Returning [`AnytimeControl::Stop`] prevents every later planner operation,
/// allocation, cache insertion, and telemetry update. A replay that returns the
/// same controls observes the same deterministic prefix.
///
/// # Errors
/// Returns [`PlanError`] for an invalid query, resource-driving ladder, cache
/// replay, cost table, observer-time hint, or planner computation.
#[allow(clippy::too_many_arguments)]
pub fn run_anytime_observed<F>(
    family: &ProblemFamily,
    theta: f64,
    tol: f64,
    budget_ladder: &[f64],
    rung_cells: &[usize],
    cache: &mut dyn AnswerCache,
    costs: &mut CostTable,
    observer: F,
) -> Result<AnytimeReport, PlanError>
where
    F: FnMut(&IntervalStep) -> AnytimeControl,
{
    validate_finite(theta, "theta")?;
    validate_positive_finite(tol, "tolerance")?;
    validate_rung_cells(rung_cells)?;
    validate_budget_ladder(budget_ladder)?;
    let final_budget = budget_ladder
        .last()
        .copied()
        .ok_or(PlanError::EmptySequence {
            field: "budget_ladder",
        })?;
    let mut driver = AnytimeDriver::new(budget_ladder, tol, observer)?;
    let observed = plan_observed(
        family,
        theta,
        tol,
        final_budget,
        rung_cells,
        cache,
        costs,
        &mut driver,
    )?;
    let planned_cost = outcome_cost(&observed.outcome);
    validate_bound(planned_cost, "planned_cost")?;
    if planned_cost > final_budget {
        return Err(PlanError::NumericalFailure {
            stage: "anytime planner exceeded its admitted budget",
        });
    }
    if !observed.stopped
        && let Some((certificate, mesh)) = outcome_certificate(&observed.outcome)
    {
        driver.finish(planned_cost, certificate, mesh, costs)?;
    }

    let stopped_by_observer = observed.stopped || driver.stopped;
    let refusal = if driver.trajectory.last().is_some_and(|step| step.discharged) {
        None
    } else if let Some(step) = driver.trajectory.last() {
        let disposition = if stopped_by_observer {
            "STOPPED by the operational observer"
        } else {
            "REFUSED at the requested tolerance"
        };
        let spend = if stopped_by_observer {
            format!(
                "after {:.1} cells spent at budget rung {:.1}",
                step.spent, step.budget
            )
        } else {
            format!("within the final {final_budget:.1}-cell budget")
        };
        Some(format!(
            "{disposition}: achieved a certified ±{:.3e} (verified) {spend}; {}. \
             No best-effort point estimate is returned.",
            step.bound, step.hint
        ))
    } else {
        let reason = outcome_reason(&observed.outcome)
            .unwrap_or("the admitted execution produced no independently verified candidate");
        Some(format!(
            "REFUSED without a certified interval: {reason}. No best-effort point \
             estimate or evidence color is returned."
        ))
    };
    Ok(AnytimeReport {
        trajectory: driver.trajectory,
        refusal,
        stopped_by_observer,
    })
}

fn outcome_cost(outcome: &PlanOutcome) -> f64 {
    match outcome {
        PlanOutcome::Discharged { cost, .. }
        | PlanOutcome::RefusedWithBest { cost, .. }
        | PlanOutcome::RefusedWithoutAnswer { cost, .. } => *cost,
    }
}

fn outcome_certificate(outcome: &PlanOutcome) -> Option<(&VerifierCertificate, &[f64])> {
    match outcome {
        PlanOutcome::Discharged {
            certificate, mesh, ..
        } => Some((certificate, mesh)),
        PlanOutcome::RefusedWithBest {
            best_certificate,
            best_mesh,
            ..
        } => Some((best_certificate, best_mesh)),
        PlanOutcome::RefusedWithoutAnswer { .. } => None,
    }
}

fn outcome_reason(outcome: &PlanOutcome) -> Option<&str> {
    match outcome {
        PlanOutcome::Discharged { .. } => None,
        PlanOutcome::RefusedWithBest { reason, .. }
        | PlanOutcome::RefusedWithoutAnswer { reason, .. } => Some(reason),
    }
}

struct AnytimeDriver<'a, F> {
    budgets: &'a [f64],
    tolerance: f64,
    next_budget: usize,
    trajectory: Vec<IntervalStep>,
    callback: F,
    stopped: bool,
    discharged: bool,
}

impl<'a, F> AnytimeDriver<'a, F>
where
    F: FnMut(&IntervalStep) -> AnytimeControl,
{
    fn new(budgets: &'a [f64], tolerance: f64, callback: F) -> Result<Self, PlanError> {
        let mut trajectory = Vec::new();
        trajectory
            .try_reserve_exact(budgets.len())
            .map_err(|_| PlanError::AllocationFailed {
                stage: "anytime trajectory",
                requested: budgets.len(),
            })?;
        Ok(Self {
            budgets,
            tolerance,
            next_budget: 0,
            trajectory,
            callback,
            stopped: false,
            discharged: false,
        })
    }

    fn emit(
        &mut self,
        spent: f64,
        certificate: &VerifierCertificate,
        mesh: &[f64],
        costs: &CostTable,
    ) -> Result<PlanControl, PlanError> {
        if self.stopped || self.discharged || self.next_budget >= self.budgets.len() {
            return Ok(if self.stopped {
                PlanControl::Stop
            } else {
                PlanControl::Continue
            });
        }
        let budget = self.budgets[self.next_budget];
        if spent > budget {
            return Err(PlanError::NumericalFailure {
                stage: "anytime emitted certificate exceeds rung budget",
            });
        }
        let bound = certificate.bound();
        validate_bound(bound, "anytime_certificate_bound")?;
        let discharged = bound <= self.tolerance;
        self.trajectory.push(IntervalStep {
            budget,
            spent,
            bound,
            color: certificate.color().clone(),
            verifier_family: certificate.verifier_family(),
            flux_hash: certificate.flux_hash(),
            hint: tighten_hint(bound, self.tolerance, spent, costs, hot_region_of(mesh))?,
            discharged,
        });
        self.next_budget += 1;
        self.discharged = discharged;
        let emitted = self.trajectory.last().ok_or(PlanError::NumericalFailure {
            stage: "anytime emitted-step retention",
        })?;
        if (self.callback)(emitted) == AnytimeControl::Stop {
            self.stopped = true;
            Ok(PlanControl::Stop)
        } else {
            Ok(PlanControl::Continue)
        }
    }

    fn finish(
        &mut self,
        spent: f64,
        certificate: &VerifierCertificate,
        mesh: &[f64],
        costs: &CostTable,
    ) -> Result<(), PlanError> {
        while !self.stopped && !self.discharged && self.next_budget < self.budgets.len() {
            self.emit(spent, certificate, mesh, costs)?;
        }
        Ok(())
    }
}

impl<F> PlanObserver for AnytimeDriver<'_, F>
where
    F: FnMut(&IntervalStep) -> AnytimeControl,
{
    fn before_work(
        &mut self,
        spent: f64,
        next_cost: f64,
        best: Option<(&VerifierCertificate, &[f64])>,
        costs: &CostTable,
    ) -> Result<PlanControl, PlanError> {
        validate_bound(spent, "anytime_spent")?;
        validate_positive_finite(next_cost, "anytime_next_cost")?;
        let next_total = spent + next_cost;
        if !next_total.is_finite() {
            return Err(PlanError::NumericalFailure {
                stage: "anytime next-operation cost accumulation",
            });
        }
        while !self.stopped
            && !self.discharged
            && self
                .budgets
                .get(self.next_budget)
                .is_some_and(|budget| *budget < next_total)
        {
            if let Some((certificate, mesh)) = best {
                if self.emit(spent, certificate, mesh, costs)? == PlanControl::Stop {
                    return Ok(PlanControl::Stop);
                }
            } else {
                self.next_budget += 1;
            }
        }
        Ok(if self.stopped {
            PlanControl::Stop
        } else {
            PlanControl::Continue
        })
    }

    fn certified(
        &mut self,
        spent: f64,
        best: (&VerifierCertificate, &[f64]),
        costs: &CostTable,
    ) -> Result<PlanControl, PlanError> {
        while self
            .budgets
            .get(self.next_budget)
            .is_some_and(|budget| *budget < spent)
        {
            self.next_budget += 1;
        }
        if self
            .budgets
            .get(self.next_budget)
            .is_some_and(|budget| spent <= *budget)
        {
            self.emit(spent, best.0, best.1, costs)
        } else {
            Ok(PlanControl::Continue)
        }
    }
}

fn validate_budget_ladder(budget_ladder: &[f64]) -> Result<(), PlanError> {
    if budget_ladder.is_empty() {
        return Err(PlanError::EmptySequence {
            field: "budget_ladder",
        });
    }
    if budget_ladder.len() > MAX_BUDGET_RUNGS {
        return Err(PlanError::ResourceLimit {
            field: "budget_ladder",
            requested: budget_ladder.len(),
            limit: MAX_BUDGET_RUNGS,
        });
    }
    for (index, &budget) in budget_ladder.iter().enumerate() {
        validate_budget(budget, "budget_ladder").map_err(|_| PlanError::InvalidSequenceEntry {
            field: "budget_ladder",
            index,
            reason: "budgets must be finite, positive, and within exact cell-accounting range",
        })?;
        if index > 0 && budget <= budget_ladder[index - 1] {
            return Err(PlanError::InvalidSequenceEntry {
                field: "budget_ladder",
                index,
                reason: "budgets must be strictly increasing",
            });
        }
    }
    Ok(())
}

/// The densest-mesh window (where refinement concentrated): the hint's
/// "where the money goes". `None` on uniform meshes.
fn hot_region_of(mesh: &[f64]) -> Option<(f64, f64)> {
    if mesh.len() < 3 {
        return None;
    }
    let mut min_h = f64::INFINITY;
    let mut max_h = 0.0f64;
    for e in 0..mesh.len() - 1 {
        let h = mesh[e + 1] - mesh[e];
        min_h = min_h.min(h);
        max_h = max_h.max(h);
    }
    if max_h < 2.0 * min_h {
        return None; // effectively uniform: no hot region to name
    }
    // The window spanned by the finest quartile of elements.
    let cutoff = 2.0 * min_h;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for e in 0..mesh.len() - 1 {
        if mesh[e + 1] - mesh[e] <= cutoff {
            lo = lo.min(mesh[e]);
            hi = hi.max(mesh[e + 1]);
        }
    }
    (lo < hi).then_some((lo, hi))
}
