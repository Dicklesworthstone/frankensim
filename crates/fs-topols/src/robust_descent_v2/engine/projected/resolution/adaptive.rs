//! Adaptive comparison inside the original projected candidate transaction.
use super::super::*;
use crate::adaptive_dwr::gain::{AdaptiveDescentCandidate, AdaptiveDescentStage,
    AdaptiveDescentUpdate, AdaptiveGainPolicy, AdaptiveGainReport, AdaptiveGainStage,
    AdaptiveGainStop, compare_adaptive_designs_controlled};

fn bind(report: &AdaptiveGainReport, baseline: &MultiLoadProjectedState,
    candidate: &MultiLoadProjectedState) -> Result<(), CutFemError> {
    let first = report.rows.first().ok_or_else(|| invalid("adaptive pair has no original-grid row"))?;
    for (side, state) in [baseline, candidate].into_iter().enumerate() {
        let measured = &report.initial_case_compliances[side];
        let (row, area) = if side == 0 {
            (&first.baseline, first.baseline_areas[0])
        } else { (&first.candidate, first.candidate_areas[0]) };
        if report.snapshots[side] != state.snapshot || measured.len() != state.case_compliances.len()
            || measured.iter().zip(&state.case_compliances).any(|(a, b)| a.to_bits() != b.to_bits())
            || row.coarse_compliance.to_bits() != state.objective.to_bits()
            || area.to_bits() != state.volume.to_bits() {
            return Err(invalid("adaptive comparison does not reproduce every original projected response"));
        }
    }
    Ok(())
}

impl MultiLoadProjectedOptimizer {
    /// Resolve each coarse-admissible gain using a jointly refined analysis grid
    /// BEFORE replacing geometry, stress, cached solutions, ordinal or multiplier.
    /// Refused comparisons contract through the ORIGINAL bounded candidate search.
    /// Existing area and unweighted sampled-stress checks remain mandatory.
    ///
    /// Each candidate reserves N*(1 + 4*max_rounds) slots in the existing lifetime
    /// solve allowance before proposal work: N original equilibria and both
    /// designs' two-level assessments. Early success spends only entered cases;
    /// each entered estimator conservatively charges two slots even on failure.
    /// Spent work survives cancellation and the existing checkpoint/restore path.
    ///
    /// The extra policy is caller-owned and must accompany that checkpoint.
    /// This method supports weighted sums, including one case, with or without
    /// an original sampled-stress limit; it does not support stress restoration.
    /// Stress is NOT certified or rechecked on adaptive grids by this method.
    /// Original CG polls remain; each two-solve DWR case remains indivisible.
    pub fn advance_one_adaptive_dwr_polling<B>(
        &mut self, policy: AdaptiveGainPolicy, poll_iters: usize,
        mut control: impl FnMut(AdaptiveDescentStage) -> ControlFlow<B>,
    ) -> Result<ControlFlow<B, AdaptiveDescentUpdate>, CutFemError> {
        let settings = self.kernel.settings;
        let count = self.kernel.load_cases.len();
        policy.validate(settings.level, count)?;
        if self.kernel.aggregate != RobustAggregate::WeightedSum
            || self.restoration_reduction.is_some() || poll_iters == 0 {
            return Err(invalid("adaptive gain descent needs weighted-sum, positive polling and no stress restoration"));
        }
        if self.next_iteration == settings.iterations {
            return Ok(ControlFlow::Continue(AdaptiveDescentUpdate {
                progress: MultiLoadProjectedProgress::IterationLimit, candidates: Vec::new(),
            }));
        }
        let levels = policy.refinement.max_rounds.checked_mul(4).and_then(|n| n.checked_add(1))
            .ok_or_else(|| invalid("adaptive candidate work overflow"))?;
        let required = count.checked_mul(levels).ok_or_else(|| invalid("adaptive case family overflow"))?;
        if self.solves_started.checked_add(required).is_none_or(|total| total > self.controls.max_solves) {
            return Ok(ControlFlow::Continue(AdaptiveDescentUpdate {
                progress: MultiLoadProjectedProgress::SolveBudget(Vec::new()), candidates: Vec::new(),
            }));
        }
        let baseline = self.current();
        let geometry = self.current.phi.clone();
        let loads = self.kernel.load_cases.clone();
        let decrease = self.controls.min_relative_improvement;
        let mut candidates = Vec::new();
        candidates.try_reserve_exact(self.controls.max_candidates)
            .map_err(|_| invalid("adaptive candidate history allocation refused"))?;
        let controller = std::cell::RefCell::new(&mut control);
        let progress = self.advance_one_admitted_scheduled(Some(poll_iters), levels,
            |index, candidate, _, spent, _| {
                let result = compare_adaptive_designs_controlled(&geometry, &candidate.phi,
                    &loads, settings, policy, decrease, |stage| {
                        if let ControlFlow::Break(reason) = (controller.borrow_mut())(
                            AdaptiveDescentStage::Candidate { index, stage },
                        ) { return ControlFlow::Break(reason); }
                        // The owner admitted the complete worst-case family before
                        // evolving this candidate. Failed/late-cancelled work stays spent.
                        if matches!(stage, AdaptiveGainStage::BeforeCase { .. }) { *spent += 2; }
                        ControlFlow::Continue(())
                    });
                let (report, refusal) = match result {
                    Ok(ControlFlow::Break(reason)) => return Ok(ControlFlow::Break(reason)),
                    Ok(ControlFlow::Continue(report)) => {
                        bind(&report, &baseline, &MultiLoadProjectedState::of(candidate))?;
                        let passed = report.stop == AdaptiveGainStop::ResolvedGain
                            && report.rows.last().is_some_and(|row| row.admissible);
                        let refusal = (!passed).then(|| format!("adaptive gain: {}", report.stop.as_str()));
                        (Some(report), refusal)
                    }
                    Err(error) => (None, Some(format!("adaptive gain solve: {error}"))),
                };
                candidates.push(AdaptiveDescentCandidate { index, report, refusal: refusal.clone() });
                Ok(ControlFlow::Continue(refusal))
            }, |stage| (controller.borrow_mut())(AdaptiveDescentStage::Optimizer(stage)))?;
        match progress {
            ControlFlow::Break(reason) => Ok(ControlFlow::Break(reason)),
            ControlFlow::Continue(progress) => Ok(ControlFlow::Continue(AdaptiveDescentUpdate {
                progress, candidates,
            })),
        }
    }
}

#[cfg(test)]
#[path = "adaptive_tests.rs"]
mod tests;
