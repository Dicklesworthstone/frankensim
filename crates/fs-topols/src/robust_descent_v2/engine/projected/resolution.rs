//! All-grid/all-load admission in the existing projected search, not a new solver.
use super::*;
use crate::refinement::prolongate_level_set;
use crate::resolution::ResolutionPolicy;
use crate::robust_resolution::{MultiLoadMeshCandidate, MultiLoadMeshProgress, MultiLoadMeshStage,
    MultiLoadResolutionReport, MultiLoadResolutionRung, MultiLoadResolutionStage};
use crate::SampledStressEvaluation;

fn rung(level: u32, stress: &RobustSampledStressEvaluation) -> MultiLoadResolutionRung {
    MultiLoadResolutionRung { level, cases: (0..stress.case_compliances.len()).map(|case|
        SampledStressEvaluation {
            compliance: stress.case_compliances[case], volume: stress.volume,
            sampled_max_von_mises: stress.case_sampled_max_von_mises[case],
            max_location: stress.case_max_locations[case], sample_count: stress.case_sample_counts[case],
            snapshot: stress.snapshot,
        }).collect() }
}

// Coarse fields have ALREADY been solved and sampled by the original owner.
// Only finer solves are new work. Callers reserve all of them before starting.
#[allow(clippy::too_many_arguments)]
fn assess<B>(geometry: &GridSdf, coarse: &RobustSampledStressEvaluation,
    loads: &[RobustLoadCase], settings: OptimizeSettings, aggregate: RobustAggregate,
    policy: ResolutionPolicy, poll_iters: usize, spent: &mut usize,
    mut control: impl FnMut(MultiLoadResolutionStage) -> ControlFlow<B>)
    -> Result<ControlFlow<B, MultiLoadResolutionReport>, CutFemError>
{
    let mut rungs = vec![rung(settings.level, coarse)];
    let mut field = geometry.clone();
    for level in settings.level + 1..=settings.level + policy.extra_levels {
        if let ControlFlow::Break(reason) = control(MultiLoadResolutionStage::Refine(level)) {
            return Ok(ControlFlow::Break(reason));
        }
        field = prolongate_level_set(&field, &[])?.geometry;
        let kernel = Kernel::new(&field, loads, OptimizeSettings { level, ..settings }, aggregate)?;
        let state = match kernel.evaluate_scheduled(field.clone(), Some(poll_iters), |case, progress| {
            let stage = match progress {
                CaseProgress::Start => MultiLoadResolutionStage::CaseSolve { level, case, complete: false },
                CaseProgress::Complete => MultiLoadResolutionStage::CaseSolve { level, case, complete: true },
                CaseProgress::Iterations(iterations) => MultiLoadResolutionStage::CaseIterations { level, case, iterations },
            };
            if let ControlFlow::Break(reason) = control(stage) { return ControlFlow::Break(reason); }
            // Before the actual solve; numerical refusal and interruption spend it.
            // The complete family fits the checked remaining allowance, so no overflow.
            if matches!(progress, CaseProgress::Start) { *spent += 1; }
            ControlFlow::Continue(())
        })? {
            ControlFlow::Continue(state) => state,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        let stress = match stress::sample_state_controlled(&kernel, &state,
            |case, cell| control(MultiLoadResolutionStage::StressCell { level, case, cell }))?
        {
            ControlFlow::Continue(stress) => stress,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        rungs.push(rung(level, &stress));
    }
    if let ControlFlow::Break(reason) = control(MultiLoadResolutionStage::Publish) {
        return Ok(ControlFlow::Break(reason));
    }
    Ok(ControlFlow::Continue(MultiLoadResolutionReport { rungs }))
}

impl MultiLoadProjectedOptimizer {
    /// Require every independent load to pass the ORIGINAL stress limit and
    /// observed resolution policy on all declared grids, and the selected
    /// aggregate to improve on every matching grid, before accepting a design.
    ///
    /// Uses existing CutFEM equilibrium, stress samples, projection and proposal
    /// code. Fine displacement fields are independently solved, never interpolated;
    /// the level set is prolonged without changing loads, material or area target.
    /// Zero-weight loads remain mandatory. A bad baseline returns measurements
    /// without geometry search. This lane requires an installed stress limit and
    /// does not combine with stress restoration.
    ///
    /// ALL actual case starts, including failed/interrupted fine solves, spend
    /// `max_solves`. A complete grid/load family must fit before each candidate.
    /// Cancellation never replaces accepted geometry, stress, ordinal or search
    /// multiplier, but spent work is retained in the existing checkpoint format.
    /// `policy` is caller-owned and must be retained beside that checkpoint.
    /// As elsewhere, assembly, prolongation and one cell's probes are indivisible.
    pub fn advance_one_resolution_polling<B>(&mut self, policy: ResolutionPolicy, poll_iters: usize,
        mut control: impl FnMut(MultiLoadMeshStage) -> ControlFlow<B>)
        -> Result<ControlFlow<B, MultiLoadMeshProgress>, CutFemError>
    {
        let settings = self.kernel.settings;
        policy.validate(settings.level)?;
        let limit = self.stress_limit.ok_or_else(|| invalid("multi-load mesh checks require an installed sampled stress limit"))?;
        if poll_iters == 0 || self.restoration_reduction.is_some() {
            return Err(invalid("mesh-checked multi-load descent requires positive polling and no stress-restoration policy"));
        }
        if self.next_iteration == settings.iterations {
            return Ok(ControlFlow::Continue(MultiLoadMeshProgress::IterationLimit));
        }
        let loads = self.kernel.load_cases.clone();
        let aggregate = self.kernel.aggregate;
        let extra = loads.len().checked_mul(policy.extra_levels as usize)
            .ok_or_else(|| invalid("mesh/load solve count overflowed"))?;
        if self.controls.max_solves - self.solves_started < extra {
            return Ok(ControlFlow::Continue(MultiLoadMeshProgress::SolveBudget));
        }
        let coarse = self.current_stress.as_ref().ok_or_else(|| invalid("missing complete baseline stress family"))?;
        let baseline = match assess(&self.current.phi, coarse, &loads, settings, aggregate,
            policy, poll_iters, &mut self.solves_started,
            |stage| control(MultiLoadMeshStage::Baseline(stage)))?
        {
            ControlFlow::Continue(report) => report,
            ControlFlow::Break(reason) => return Ok(ControlFlow::Break(reason)),
        };
        if let Some(reason) = baseline.refusal(policy, settings.volfrac, limit, &loads)? {
            return Ok(ControlFlow::Continue(MultiLoadMeshProgress::UnresolvedBaseline { baseline, reason }));
        }
        let decrease = self.controls.min_relative_improvement;
        let controller = std::cell::RefCell::new(&mut control);
        let mut candidates = Vec::new();
        let progress = self.advance_one_admitted_scheduled(Some(poll_iters), policy.extra_levels as usize + 1,
            |index, candidate, coarse, spent, _| {
                let coarse = coarse.ok_or_else(|| invalid("missing complete candidate stress family"))?;
                let result = assess(&candidate.phi, coarse, &loads, settings, aggregate, policy, poll_iters, spent,
                    |stage| (controller.borrow_mut())(MultiLoadMeshStage::Candidate { index, stage }));
                let (report, refusal) = match result {
                    Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                    Ok(ControlFlow::Continue(report)) => {
                        let refusal = report.comparison_refusal(&baseline, policy, settings.volfrac,
                            limit, &loads, aggregate, decrease)?;
                        (Some(report), refusal)
                    }
                    Err(error) => (None, Some(format!("independent-load mesh solve: {error}"))),
                };
                candidates.push(MultiLoadMeshCandidate { index, report, refusal: refusal.clone() });
                Ok(ControlFlow::Continue(refusal))
            }, |stage| (controller.borrow_mut())(MultiLoadMeshStage::Optimizer(stage)))?;
        match progress {
            ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(progress) => Ok(ControlFlow::Continue(MultiLoadMeshProgress::Searched {
                baseline, candidates, progress,
            })),
        }
    }
}

#[cfg(test)]
mod tests;
