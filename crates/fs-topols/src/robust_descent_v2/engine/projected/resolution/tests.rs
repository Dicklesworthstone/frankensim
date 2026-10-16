use super::*;
use crate::{DesignBoxEdge, evaluate_robust_sampled_stress};

fn loads() -> Vec<RobustLoadCase> {
    vec![
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
        // A physically larger, opposite load with NO objective weight.
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 4.0], 0.0).unwrap(),
    ]
}
fn policy() -> ResolutionPolicy {
    ResolutionPolicy { extra_levels: 1, absolute_compliance_tolerance: 1e8,
        relative_compliance_tolerance: 0.0, area_tolerance: 0.1 }
}
fn owner(max_solves: usize, aggregate: RobustAggregate) -> MultiLoadProjectedOptimizer {
    let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = geometry.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % 9 == 0 || index % 9 == 8).collect();
    MultiLoadProjectedOptimizer::new(geometry, &loads(),
        OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default() }, aggregate, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        MultiLoadProjectedSettings { max_candidates: 6, max_solves, ..MultiLoadProjectedSettings::default() },
    ).unwrap().with_sampled_stress_limit(SampledStressLimit::new(1e8, 0.0).unwrap()).unwrap()
}
fn measured(owner: &MultiLoadProjectedOptimizer) -> MultiLoadResolutionReport {
    let mut spent = owner.solves_started();
    let result = assess(owner.geometry(), owner.current_stress().unwrap(), &loads(),
        owner.kernel.settings, owner.kernel.aggregate, policy(), 1, &mut spent,
        |_| ControlFlow::<()>::Continue(())).unwrap();
    assert_eq!(spent, owner.solves_started() + 2);
    let ControlFlow::Continue(report) = result else { panic!("not cancelled") };
    report
}
fn advance(owner: &mut MultiLoadProjectedOptimizer) -> MultiLoadMeshProgress {
    let result = owner.advance_one_resolution_polling(policy(), 1,
        |_| ControlFlow::<()>::Continue(())).unwrap();
    let ControlFlow::Continue(result) = result else { panic!("not cancelled") };
    result
}

#[test]
fn projected_mesh_fine_fields_are_real_independent_solves_including_zero_weight_loads() {
    let owner = owner(1000, RobustAggregate::WeightedSum);
    let report = measured(&owner);
    let fine = prolongate_level_set(owner.geometry(), &[]).unwrap().geometry;
    let replay = evaluate_robust_sampled_stress(&fine, &loads(),
        OptimizeSettings { level: 4, ..owner.kernel.settings }, RobustAggregate::WeightedSum).unwrap();
    assert_eq!(report.rungs[0], rung(3, owner.current_stress().unwrap()));
    assert_eq!(report.rungs[1], rung(4, &replay));
    assert_ne!(report.rungs[0].cases[0].compliance.to_bits(), report.rungs[1].cases[0].compliance.to_bits());
    for grid in &report.rungs {
        let a = &grid.cases[0]; let b = &grid.cases[1];
        assert!(a.compliance > 0.0 && a.sampled_max_von_mises > 0.0);
        assert!((b.compliance / a.compliance - 16.0).abs() < 1e-7);
        assert!((b.sampled_max_von_mises / a.sampled_max_von_mises - 4.0).abs() < 1e-7);
        assert!(a.sample_count > 0 && b.sample_count > 0);
    }
}

#[test]
fn projected_mesh_finer_zero_weight_overstress_and_partial_families_cannot_pass() {
    let owner = owner(1000, RobustAggregate::WeightedSum);
    let mut report = measured(&owner);
    let limit = owner.stress_limit().unwrap();
    assert!(report.refusal(policy(), 0.6, limit, &loads()).unwrap().is_none());
    // Mutate measured public data to isolate the all-grid/all-load gate, not to
    // impersonate a new physics oracle. The zero-weight fine load must be read.
    report.rungs[1].cases[1].sampled_max_von_mises = 2.0 * limit.admitted_max();
    let reason = report.refusal(policy(), 0.6, limit, &loads()).unwrap().unwrap();
    assert!(reason.contains("case 1") && reason.contains("level 4"));
    report.rungs[1].cases.pop();
    assert!(report.refusal(policy(), 0.6, limit, &loads()).is_err());
}

#[test]
fn projected_mesh_uses_matching_aggregates_and_individual_load_resolution_not_weighted_cancellation() {
    let state = |c| SampledStressEvaluation { compliance: c, volume: 0.6,
        sampled_max_von_mises: 1.0, max_location: [0.5, 0.5], sample_count: 10, snapshot: 7 };
    let baseline = MultiLoadResolutionReport { rungs: vec![
        MultiLoadResolutionRung { level: 3, cases: vec![state(2.0), state(32.0)] },
        MultiLoadResolutionRung { level: 4, cases: vec![state(2.0), state(32.0)] },
    ] };
    let mut candidate = baseline.clone();
    for grid in &mut candidate.rungs { grid.cases[0].compliance = 1.8; }
    let limit = SampledStressLimit::new(10.0, 0.0).unwrap();
    for aggregate in [RobustAggregate::WeightedSum, RobustAggregate::WorstWeightedCase] {
        assert!(candidate.comparison_refusal(&baseline, policy(), 0.6, limit, &loads(), aggregate, 0.0).unwrap().is_none());
        candidate.rungs[1].cases[0].compliance = 2.1;
        assert!(candidate.comparison_refusal(&baseline, policy(), 0.6, limit, &loads(), aggregate, 0.0)
            .unwrap().unwrap().contains("level 4"));
        candidate.rungs[1].cases[0].compliance = 1.8;
    }
    candidate.rungs[1].cases[1].compliance = 40.0;
    let tight = ResolutionPolicy { absolute_compliance_tolerance: 1.0, ..policy() };
    assert!(candidate.refusal(tight, 0.6, limit, &loads()).unwrap().unwrap().contains("case 1"));
    candidate.rungs[1].cases[1].compliance = f64::NAN;
    assert!(candidate.refusal(tight, 0.6, limit, &loads()).is_err());
}

#[test]
fn projected_mesh_unresolved_baseline_never_proposes_or_changes_accepted_state() {
    let mut owner = owner(1000, RobustAggregate::WeightedSum);
    let before = owner.current();
    let stress = owner.current_stress().unwrap().clone();
    let multiplier = owner.search_multiplier().to_bits();
    let strict = ResolutionPolicy { absolute_compliance_tolerance: 1e-30, ..policy() };
    let result = owner.advance_one_resolution_polling(strict, 1, |stage| {
        assert!(!matches!(stage, MultiLoadMeshStage::Optimizer(_)));
        ControlFlow::<()>::Continue(())
    }).unwrap();
    assert!(matches!(result, ControlFlow::Continue(MultiLoadMeshProgress::UnresolvedBaseline { .. })));
    assert_eq!(owner.current(), before);
    assert_eq!(owner.current_stress(), Some(&stress));
    assert_eq!(owner.search_multiplier().to_bits(), multiplier);
    assert_eq!(owner.next_iteration(), 0);
    assert_eq!(owner.solves_started(), 4);
}

#[test]
fn projected_mesh_reserves_complete_families_in_the_original_solve_allowance() {
    let mut short = owner(3, RobustAggregate::WeightedSum);
    let result = short.advance_one_resolution_polling(policy(), 1, |_| -> ControlFlow<()> {
        panic!("an unfunded baseline must do no work");
    }).unwrap();
    assert!(matches!(result, ControlFlow::Continue(MultiLoadMeshProgress::SolveBudget)));
    assert_eq!(short.solves_started(), 2);
    let mut enough_baseline = owner(5, RobustAggregate::WeightedSum);
    let result = enough_baseline.advance_one_resolution_polling(policy(), 1, |stage| {
        assert!(!matches!(stage, MultiLoadMeshStage::Optimizer(_)));
        ControlFlow::<()>::Continue(())
    }).unwrap();
    assert!(matches!(result, ControlFlow::Continue(MultiLoadMeshProgress::Searched {
        progress: MultiLoadProjectedProgress::SolveBudget(_), .. })));
    assert_eq!(enough_baseline.solves_started(), 4);
    assert_eq!(enough_baseline.next_iteration(), 0);
    // A completed owner never starts another baseline grid assessment.
    enough_baseline.next_iteration = enough_baseline.kernel.settings.iterations;
    assert!(matches!(enough_baseline.advance_one_resolution_polling(policy(), 1,
        |_| -> ControlFlow<()> { panic!("completed owner is inert") }).unwrap(),
        ControlFlow::Continue(MultiLoadMeshProgress::IterationLimit)));
}

#[test]
fn projected_mesh_cancelled_fine_case_and_candidate_keep_state_spend_work_and_retry_exactly() {
    for stop in 0..3 {
        let mut interrupted = owner(1000, RobustAggregate::WeightedSum);
        let before = interrupted.current();
        let stress = interrupted.current_stress().unwrap().clone();
        let multiplier = interrupted.search_multiplier().to_bits();
        let result = interrupted.advance_one_resolution_polling(policy(), 1, |stage| {
            let hit = match (stop, stage) {
                (0, MultiLoadMeshStage::Baseline(MultiLoadResolutionStage::CaseIterations { case: 1, .. })) => true,
                (1, MultiLoadMeshStage::Candidate { stage: MultiLoadResolutionStage::StressCell { case: 1, .. }, .. }) => true,
                (2, MultiLoadMeshStage::Candidate { stage: MultiLoadResolutionStage::Publish, .. }) => true,
                _ => false,
            };
            if hit { ControlFlow::Break("stop") } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break("stop")), "must reach the intended numerical boundary");
        assert_eq!(interrupted.current(), before);
        assert_eq!(interrupted.current_stress(), Some(&stress));
        assert_eq!(interrupted.search_multiplier().to_bits(), multiplier);
        assert_eq!(interrupted.next_iteration(), 0);
        let spent = interrupted.solves_started() - 2;
        assert!(spent > 0);
        let retry = advance(&mut interrupted);
        let mut fresh = owner(1000, RobustAggregate::WeightedSum);
        let clean = advance(&mut fresh);
        assert_eq!(std::mem::discriminant(&retry), std::mem::discriminant(&clean));
        assert_eq!(interrupted.geometry().nodes(), fresh.geometry().nodes());
        assert_eq!(interrupted.current_stress(), fresh.current_stress());
        assert_eq!(interrupted.next_iteration(), fresh.next_iteration());
        assert_eq!(interrupted.solves_started(), fresh.solves_started() + spent);
        if let MultiLoadMeshProgress::Searched { baseline, candidates, progress: MultiLoadProjectedProgress::Accepted(step) } = retry {
            let accepted = candidates.iter().find(|c| c.index == step.attempts.last().unwrap().index).unwrap();
            assert!(accepted.refusal.is_none());
            assert!(accepted.report.as_ref().unwrap().comparison_refusal(&baseline, policy(), 0.6,
                interrupted.stress_limit().unwrap(), &loads(), RobustAggregate::WeightedSum,
                interrupted.controls.min_relative_improvement).unwrap().is_none());
        }
    }
}
