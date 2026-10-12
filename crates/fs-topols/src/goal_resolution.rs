//! Signed-DWR resolution of proposed same-area compliance improvements.
//!
//! This is an optional admission hook on `ProjectedOptimizer`, not another
//! optimizer. A coarse win must survive an enriched solve, and both the change
//! in the signed DWR correction and the measured change in the improvement
//! must fit the caller's declared fraction of the coarse gain. These are
//! Estimated numerical screens, NOT continuum-error bounds or optimality tests.
//! In particular, absolute cell-indicator mass is never used as an error radius.

use std::cell::RefCell;
use std::ops::ControlFlow;

use fs_cutfem::CutFemError;
use crate::projected::{ProjectedOptimizer, ProjectedProgress, ProjectedStage};
use crate::{ComplianceDwrAssessment, ComplianceDwrStage, EvaluatedFinalState,
    assess_compliance_dwr_controlled};

/// Additional DWR work, separate from the existing proposal/projection solves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DwrGatePolicy {
    /// Largest permitted numerical change as a fraction of the coarse gain.
    /// Must be strictly between zero and one; passing is not a guarantee.
    pub max_error_fraction: f64,
    /// Two baseline solves plus two per possible candidate must fit BEFORE work.
    /// This allowance applies to one method invocation. An interrupted attempt
    /// must be charged by the caller and a retry explicitly admitted anew.
    pub max_solves_per_attempt: usize,
}

fn invalid(what: impl Into<String>) -> CutFemError {
    CutFemError::InvalidElasticityInput { what: what.into() }
}

impl DwrGatePolicy {
    /// Validate the entire bounded candidate family before starting any solve.
    pub fn validate(self, level: u32, max_candidates: usize) -> Result<(), CutFemError> {
        if !(1..=5).contains(&level) || !(1..=64).contains(&max_candidates)
            || !self.max_error_fraction.is_finite()
            || self.max_error_fraction <= 0.0 || self.max_error_fraction >= 1.0
        {
            return Err(invalid("DWR gate requires level 1..=5, 1..=64 candidates and an error fraction strictly between zero and one"));
        }
        let required = max_candidates.checked_add(1).and_then(|n| n.checked_mul(2))
            .ok_or_else(|| invalid("DWR solve admission overflow"))?;
        if self.max_solves_per_attempt < required || self.max_solves_per_attempt > 130 {
            return Err(invalid(format!("DWR gate requires {required}..=130 additional solves per attempt")));
        }
        Ok(())
    }
}

/// Bounded scalar projection of real DWR output; not an independently verified
/// certificate. The producer below binds it to the accepted-state mechanics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DwrMeasurement {
    pub snapshot: u64,
    pub level: u32,
    pub area: f64,
    pub coarse_compliance: f64,
    pub enriched_compliance: f64,
    pub eta_signed: f64,
    pub coarse_relative_residual: f64,
    pub enriched_relative_residual: f64,
}

impl DwrMeasurement {
    /// Check finite numerical evidence and the canonical true-residual regime.
    /// This validates representation, not the provenance of caller-built data.
    pub fn validate(self) -> Result<(), CutFemError> {
        if !(1..=5).contains(&self.level) || !self.area.is_finite() || self.area <= 0.0
            || !self.coarse_compliance.is_finite() || self.coarse_compliance < 0.0
            || !self.enriched_compliance.is_finite() || self.enriched_compliance < 0.0
            || !self.eta_signed.is_finite()
            || [self.coarse_relative_residual, self.enriched_relative_residual].iter()
                .any(|r| !r.is_finite() || *r < 0.0 || *r >= 1e-12)
        {
            return Err(invalid("DWR gate requires finite same-geometry measurements and recomputed residuals below 1e-12"));
        }
        Ok(())
    }

    fn from_assessment(a: ComplianceDwrAssessment, expected: EvaluatedFinalState,
        level: u32) -> Result<Self, CutFemError>
    {
        let measured = Self { snapshot: a.snapshot, level: a.level, area: a.volume,
            coarse_compliance: a.estimate.j_primal, enriched_compliance: a.estimate.j_enriched,
            eta_signed: a.estimate.eta_signed,
            coarse_relative_residual: a.estimate.primal_relative_residual,
            enriched_relative_residual: a.estimate.enriched_relative_residual };
        measured.validate()?;
        if measured.level != level || measured.snapshot != expected.snapshot
            || measured.area.to_bits() != expected.volume.to_bits()
            || measured.coarse_compliance.to_bits() != expected.compliance.to_bits()
        {
            return Err(invalid("DWR coarse solve does not reproduce the exact accepted/projected mechanics"));
        }
        Ok(measured)
    }
}

/// Recomputable arithmetic of the Estimated improvement screen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DwrGainDecision {
    pub coarse_gain: f64,
    pub enriched_gain: f64,
    pub signed_correction_change: f64,
    pub observed_gain_change: f64,
    pub allowance: f64,
    /// None means this numerical screen passed; it does not mint a certificate.
    pub refusal: Option<&'static str>,
}

/// Compare numerical evidence at identical levels. Public inputs are data, not
/// authority. The optimizer only calls this on its own bound, actual solves.
/// A shared bias can survive this screen: no continuum guarantee is inferred.
pub fn compare_dwr_gains(baseline: DwrMeasurement, candidate: DwrMeasurement,
    policy: DwrGatePolicy, min_relative_improvement: f64) -> Result<DwrGainDecision, CutFemError>
{
    baseline.validate()?;
    candidate.validate()?;
    policy.validate(baseline.level, 1)?;
    if candidate.level != baseline.level || !min_relative_improvement.is_finite()
        || !(0.0..1.0).contains(&min_relative_improvement)
    {
        return Err(invalid("DWR comparison requires identical levels and relative decrease in [0,1)"));
    }
    let coarse_gain = baseline.coarse_compliance - candidate.coarse_compliance;
    let enriched_gain = baseline.enriched_compliance - candidate.enriched_compliance;
    let signed_correction_change = (baseline.eta_signed - candidate.eta_signed).abs();
    let observed_gain_change = (enriched_gain - coarse_gain).abs();
    let allowance = policy.max_error_fraction * coarse_gain.max(0.0);
    if [coarse_gain, enriched_gain, signed_correction_change, observed_gain_change, allowance]
        .iter().any(|v| !v.is_finite())
    {
        return Err(invalid("DWR improvement comparison overflow"));
    }
    let refusal = if !(candidate.coarse_compliance
        < baseline.coarse_compliance * (1.0 - min_relative_improvement)) {
        Some("coarse-decrease-not-met")
    } else if !(candidate.enriched_compliance
        < baseline.enriched_compliance * (1.0 - min_relative_improvement)) {
        Some("enriched-decrease-not-met")
    } else if signed_correction_change.max(observed_gain_change) > allowance {
        Some("discretization-dominates-gain")
    } else { None };
    Ok(DwrGainDecision { coarse_gain, enriched_gain, signed_correction_change,
        observed_gain_change, allowance, refusal })
}

/// Phase boundaries of DWR and the original optimizer. Each two-solve DWR call
/// remains indivisible between its BeforeEstimate/BeforePublish callbacks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DwrCheckStage {
    Baseline(ComplianceDwrStage),
    Candidate { index: usize, stage: ComplianceDwrStage },
    Optimizer(ProjectedStage),
}

#[derive(Debug, Clone)]
pub struct DwrCandidateCheck {
    pub index: usize,
    pub measurement: Option<DwrMeasurement>,
    pub decision: Option<DwrGainDecision>,
    pub refusal: Option<String>,
}

#[derive(Debug, Clone)]
pub enum DwrCheckedProgress {
    /// No further update is permitted; no DWR solves were started.
    IterationLimit,
    Searched {
        baseline: DwrMeasurement,
        candidates: Vec<DwrCandidateCheck>,
        /// Upper work debit: two for each assessment started, even if it failed.
        /// This is not the measured number of completed linear solves.
        reserved_solves: usize,
        progress: ProjectedProgress,
    },
}

impl ProjectedOptimizer {
    /// Gate a same-material step using the existing candidate admission hook.
    /// Refusal contracts the next proposal; cancellation or any error leaves
    /// geometry, ordinal, multiplier and accepted mechanics unchanged. Only
    /// otherwise admissible coarse candidates start extra DWR work.
    ///
    /// At most `2 * (1 + max_candidates)` extra solves are reserved. Assembly,
    /// both DWR solves and residual integration are indivisible; the existing
    /// proposal/evaluation CG polls are preserved. The caller owns memory and
    /// wall-time admission, including interrupted-attempt charges on resume.
    pub fn advance_one_dwr_controlled<B>(&mut self, policy: DwrGatePolicy,
        mut control: impl FnMut(DwrCheckStage) -> ControlFlow<B>)
        -> Result<ControlFlow<B, DwrCheckedProgress>, CutFemError>
    {
        let settings = self.checkpoint().settings();
        let fixture = self.checkpoint().fixture();
        let search = self.controls();
        policy.validate(settings.level, search.max_candidates)?;
        if self.checkpoint().is_complete() {
            return Ok(ControlFlow::Continue(DwrCheckedProgress::IterationLimit));
        }
        let baseline = match assess_compliance_dwr_controlled(
            self.checkpoint().geometry(), fixture, settings,
            |stage| control(DwrCheckStage::Baseline(stage)),
        )? {
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(a) => DwrMeasurement::from_assessment(a, self.current(), settings.level)?,
        };
        let mut reserved_solves = 2usize;
        let mut candidates = Vec::new();
        candidates.try_reserve_exact(search.max_candidates)
            .map_err(|_| invalid("DWR candidate evidence allocation refused"))?;
        let controller = RefCell::new(&mut control);
        let progress = self.advance_one_admitted_controlled(|index, geometry, state, _| {
            // Full-family admission above makes this debit bounded before work.
            reserved_solves += 2;
            let assessed = assess_compliance_dwr_controlled(geometry, fixture, settings,
                |stage| (controller.borrow_mut())(DwrCheckStage::Candidate { index, stage }));
            let measured = match assessed {
                Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                Ok(ControlFlow::Continue(a)) => DwrMeasurement::from_assessment(a, state, settings.level),
                Err(error) => Err(error),
            };
            let (measurement, decision, refusal) = match measured {
                Ok(m) => match compare_dwr_gains(baseline, m, policy, search.min_relative_improvement) {
                    Ok(d) => (Some(m), Some(d), d.refusal.map(str::to_string)),
                    Err(error) => (Some(m), None, Some(error.to_string())),
                },
                Err(error) => (None, None, Some(error.to_string())),
            };
            candidates.push(DwrCandidateCheck { index, measurement, decision, refusal: refusal.clone() });
            Ok(ControlFlow::Continue(refusal))
        }, |stage| (controller.borrow_mut())(DwrCheckStage::Optimizer(stage)))?;
        match progress {
            ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(progress) => Ok(ControlFlow::Continue(DwrCheckedProgress::Searched {
                baseline, candidates, reserved_solves, progress,
            })),
        }
    }
}

#[cfg(test)]
#[path = "goal_resolution/tests.rs"]
mod tests;
