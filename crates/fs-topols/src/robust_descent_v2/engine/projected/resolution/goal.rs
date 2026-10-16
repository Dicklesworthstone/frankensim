//! Weighted DWR uses the original candidate hook and lifetime solve allowance.
use super::super::*;
use crate::robust_resolution::goal::{MultiLoadDwrCandidate, MultiLoadDwrMeasurement,
    MultiLoadDwrPolicy, MultiLoadDwrProgress, MultiLoadDwrStage, compare_weighted_dwr_gains};
use crate::{WeightedComplianceDwrStage, assess_weighted_compliance_dwr_controlled};

// Each entered DWR case reserves two capped solves before its indivisible kernel.
// Failed cases can finish fewer; reservations remain in the original checkpoint.
fn assess<B>(geometry: &GridSdf, state: &MultiLoadProjectedState,
    loads: &[RobustLoadCase], settings: OptimizeSettings, spent: &mut usize, limit: usize,
    mut control: impl FnMut(WeightedComplianceDwrStage) -> ControlFlow<B>)
    -> Result<ControlFlow<B, MultiLoadDwrMeasurement>, CutFemError>
{
    let count = loads.len().checked_mul(2).ok_or_else(|| invalid("DWR work overflow"))?;
    if spent.checked_add(count).is_none_or(|total| total > limit) {
        return Err(invalid("complete DWR load family exceeds the remaining solve allowance"));
    }
    let result = assess_weighted_compliance_dwr_controlled(geometry, loads, settings, |stage| {
        if let ControlFlow::Break(reason) = control(stage) { return ControlFlow::Break(reason); }
        if matches!(stage, WeightedComplianceDwrStage::BeforeCase { .. }) { *spent += 2; }
        ControlFlow::Continue(())
    })?;
    match result {
        ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
        ControlFlow::Continue(a) => MultiLoadDwrMeasurement::from_assessment(a, loads, state, settings.level)
            .map(ControlFlow::Continue),
    }
}

impl MultiLoadProjectedOptimizer {
    /// Check weighted-sum updates using independently solved coarse/enriched
    /// responses and non-cancelling per-case error changes relative to the gain.
    /// The existing proposal, area, stress and publication rules remain active.
    /// Zero-weight cases are still solved and checked against any stress limit.
    ///
    /// Each entered DWR case charges two slots in the existing lifetime budget.
    /// Thus `solves_started` includes conservative reservations, not just actual
    /// completed solves, when this method is used. A baseline requires 2*N slots;
    /// a candidate requires 3*N including its original N equilibrium solves.
    ///
    /// This method supports weighted sums without stress restoration. Policy is
    /// retained by the caller alongside the checkpoint. Cancellation preserves
    /// accepted geometry, stress, ordinal and multiplier but not spent work.
    /// Original CG polls remain; each two-solve DWR case is indivisible.
    pub fn advance_one_dwr_polling<B>(&mut self, policy: MultiLoadDwrPolicy, poll_iters: usize,
        mut control: impl FnMut(MultiLoadDwrStage) -> ControlFlow<B>)
        -> Result<ControlFlow<B, MultiLoadDwrProgress>, CutFemError>
    {
        let settings = self.kernel.settings;
        policy.validate(settings.level, self.kernel.load_cases.len())?;
        if self.kernel.aggregate != RobustAggregate::WeightedSum || self.restoration_reduction.is_some()
            || poll_iters == 0 {
            return Err(invalid("weighted DWR descent requires weighted-sum, positive polling and no stress restoration"));
        }
        if self.next_iteration == settings.iterations {
            return Ok(ControlFlow::Continue(MultiLoadDwrProgress::IterationLimit));
        }
        let baseline_cost = 2 * self.kernel.load_cases.len();
        if self.solves_started.checked_add(baseline_cost).is_none_or(|n| n > self.controls.max_solves) {
            return Ok(ControlFlow::Continue(MultiLoadDwrProgress::SolveBudget));
        }
        if let ControlFlow::Break(reason) = control(MultiLoadDwrStage::Prepare) {
            return Ok(ControlFlow::Break(reason));
        }
        let loads = self.kernel.load_cases.clone();
        let current = self.current();
        let limit = self.controls.max_solves;
        let baseline = match assess(&self.current.phi, &current, &loads, settings,
            &mut self.solves_started, limit, |stage| control(MultiLoadDwrStage::Baseline(stage)))? {
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(measured) => measured,
        };
        let decrease = self.controls.min_relative_improvement;
        let controller = std::cell::RefCell::new(&mut control);
        let mut candidates = Vec::new();
        candidates.try_reserve_exact(self.controls.max_candidates)
            .map_err(|_| invalid("weighted DWR evidence allocation refused"))?;
        let progress = self.advance_one_admitted_scheduled(Some(poll_iters), 3,
            |index, candidate, _, spent, _| {
                let measured = assess(&candidate.phi, &MultiLoadProjectedState::of(candidate),
                    &loads, settings, spent, limit,
                    |stage| (controller.borrow_mut())(MultiLoadDwrStage::Candidate { index, stage }));
                let (measurement, decision, refusal) = match measured {
                    Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                    Ok(ControlFlow::Continue(m)) => match compare_weighted_dwr_gains(
                        &baseline, &m, &loads, policy, decrease) {
                        Ok(d) => (Some(m), Some(d), d.refusal.map(str::to_string)),
                        Err(error) => (Some(m), None, Some(error.to_string())),
                    },
                    Err(error) => (None, None, Some(format!("weighted DWR solve: {error}"))),
                };
                candidates.push(MultiLoadDwrCandidate { index, measurement, decision, refusal: refusal.clone() });
                Ok(ControlFlow::Continue(refusal))
            }, |stage| (controller.borrow_mut())(MultiLoadDwrStage::Optimizer(stage)))?;
        match progress {
            ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(progress) => Ok(ControlFlow::Continue(MultiLoadDwrProgress::Searched {
                baseline, candidates, progress,
            })),
        }
    }
}

#[cfg(test)]
#[path = "goal_tests.rs"]
mod tests;
