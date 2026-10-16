//! Measured stress restoration on every declared grid, using the same immutable
//! limit and independent load family as ordinary mesh-checked descent.
use super::*;

fn worst(rung: &MultiLoadResolutionRung) -> f64 {
    rung.cases.iter().map(|case| case.sampled_max_von_mises).fold(0.0, f64::max)
}

impl MultiLoadResolutionReport {
    /// Admit a restoration update only when the worst UNWEIGHTED sampled stress
    /// excess strictly decreases on every still-infeasible grid, by more than
    /// `reduction`, or reaches the unchanged bound. Already-feasible grids must
    /// remain feasible. Compliance may increase, but each load still passes its
    /// declared compliance-resolution and material-area checks.
    ///
    /// A candidate feasible on the optimization grid must be feasible on ALL
    /// declared grids before publication, because the existing optimizer then
    /// leaves restoration. No stress gradient or continuum bound is implied.
    /// This is also the retained-report reader's acceptance rule.
    #[allow(clippy::too_many_arguments)]
    pub fn restoration_comparison_refusal(&self, baseline: &Self,
        policy: ResolutionPolicy, target: f64, limit: SampledStressLimit,
        loads: &[RobustLoadCase], reduction: f64) -> Result<Option<String>, CutFemError>
    {
        let limit = SampledStressLimit::new(limit.max_von_mises, limit.absolute_tolerance)?;
        if !reduction.is_finite() || !(0.0..1.0).contains(&reduction) {
            return Err(invalid("invalid mesh-checked stress restoration reduction"));
        }
        if let Some(reason) = baseline.resolution_refusal(policy, target, loads)? {
            return Ok(Some(reason));
        }
        if let Some(reason) = self.resolution_refusal(policy, target, loads)? {
            return Ok(Some(reason));
        }
        if baseline.rungs[0].level != self.rungs[0].level {
            return Err(invalid("stress restoration comparison changed grid levels"));
        }
        let bound = limit.admitted_max();
        if worst(&baseline.rungs[0]) <= bound {
            return Err(invalid("mesh-checked restoration requires an infeasible optimization-grid baseline"));
        }
        let crossing = worst(&self.rungs[0]) <= bound;
        for (old, new) in baseline.rungs.iter().zip(&self.rungs) {
            if (old.cases[0].volume - new.cases[0].volume).abs() > policy.area_tolerance {
                return Ok(Some(format!("level {} changes material area beyond the allowance", new.level)));
            }
            let previous = worst(old);
            let candidate = worst(new);
            if candidate <= bound { continue; }
            if crossing {
                return Ok(Some(format!("level {} remains overstressed at the proposed restoration exit", new.level)));
            }
            if previous <= bound {
                return Ok(Some(format!("level {} loses sampled stress feasibility during restoration", new.level)));
            }
            if !(candidate - bound < (previous - bound) * (1.0 - reduction)) {
                return Ok(Some(format!("level {} does not reproduce the required worst sampled stress excess reduction", new.level)));
            }
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DesignBoxEdge;

    fn loads() -> Vec<RobustLoadCase> {
        vec![
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 4.0], 0.0).unwrap(),
        ]
    }
    fn policy() -> ResolutionPolicy {
        ResolutionPolicy { extra_levels: 1, absolute_compliance_tolerance: 1.0,
            relative_compliance_tolerance: 0.0, area_tolerance: 0.01 }
    }
    fn report(stress: [[f64; 2]; 2], compliance: f64) -> MultiLoadResolutionReport {
        MultiLoadResolutionReport { rungs: stress.iter().enumerate().map(|(level, stresses)|
            MultiLoadResolutionRung { level: 3 + level as u32,
                cases: stresses.iter().map(|&sampled_max_von_mises| SampledStressEvaluation {
                    compliance, volume: 0.6, sampled_max_von_mises,
                    max_location: [0.5, 0.5], sample_count: 10, snapshot: 7,
                }).collect() }).collect() }
    }
    fn compare(a: &MultiLoadResolutionReport, b: &MultiLoadResolutionReport)
        -> Result<Option<String>, CutFemError>
    {
        b.restoration_comparison_refusal(a, policy(), 0.6,
            SampledStressLimit::new(10.0, 0.0).unwrap(), &loads(), 0.1)
    }

    #[test]
    fn projected_restoration_checks_all_grids_and_zero_weight_cases() {
        let a = report([[2.0, 20.0], [4.0, 30.0]], 1.0);
        let b = report([[3.0, 17.0], [3.0, 24.0]], 2.0);
        // Compliance increases: restoration is not mislabeled compliance descent.
        assert!(compare(&a, &b).unwrap().is_none());
        let mut bad = b.clone();
        bad.rungs[1].cases[1].sampled_max_von_mises = 29.0;
        assert!(compare(&a, &bad).unwrap().unwrap().contains("level 4"));
        bad.rungs[1].cases[1].sampled_max_von_mises = 28.0; // equality is not strict progress
        assert!(compare(&a, &bad).unwrap().is_some());
        assert!(compare(&a, &a).unwrap().is_some());
        bad.rungs[1].cases.pop();
        assert!(compare(&a, &bad).is_err());
    }

    #[test]
    fn projected_restoration_exit_requires_all_grids_feasible() {
        let a = report([[2.0, 20.0], [4.0, 30.0]], 1.0);
        let mut b = report([[2.0, 9.0], [4.0, 15.0]], 2.0);
        // Even a large fine-grid improvement cannot end restoration prematurely.
        assert!(compare(&a, &b).unwrap().unwrap().contains("restoration exit"));
        b.rungs[1].cases[1].sampled_max_von_mises = 10.0;
        assert!(compare(&a, &b).unwrap().is_none());
        assert!(compare(&b, &b).is_err()); // feasible designs use compliance descent
    }

    #[test]
    fn projected_restoration_preserves_already_feasible_grids_and_active_case_switches() {
        let a = report([[20.0, 15.0], [4.0, 9.0]], 1.0);
        let mut b = report([[17.0, 18.0], [4.0, 10.0]], 2.0);
        assert!(compare(&a, &b).unwrap().is_none());
        b.rungs[1].cases[1].sampled_max_von_mises = 11.0;
        assert!(compare(&a, &b).unwrap().unwrap().contains("loses sampled stress feasibility"));
    }

    #[test]
    fn projected_restoration_keeps_numerical_area_and_measurement_gates() {
        let a = report([[2.0, 20.0], [4.0, 30.0]], 1.0);
        let b = report([[3.0, 17.0], [3.0, 24.0]], 2.0);
        let mut bad = b.clone();
        bad.rungs[1].cases[1].compliance = 4.0;
        assert!(compare(&a, &bad).unwrap().unwrap().contains("case 1"));
        bad = b.clone();
        for case in &mut bad.rungs[1].cases { case.volume = 0.7; }
        assert!(compare(&a, &bad).unwrap().is_some());
        bad = b.clone();
        bad.rungs[1].cases[1].sampled_max_von_mises = f64::NAN;
        assert!(compare(&a, &bad).is_err());
        bad = b.clone();
        bad.rungs[1].cases[1].sample_count = 0;
        assert!(compare(&a, &bad).is_err());
        let limit = SampledStressLimit::new(10.0, 0.0).unwrap();
        assert!(a.refusal(policy(), 0.6, limit, &loads()).unwrap().is_some());
        assert!(a.resolution_refusal(policy(), 0.6, &loads()).unwrap().is_none());
        for reduction in [-1.0, 1.0, f64::NAN] {
            assert!(b.restoration_comparison_refusal(&a, policy(), 0.6, limit,
                &loads(), reduction).is_err());
        }
    }
}
