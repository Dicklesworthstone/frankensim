use super::*;
use crate::adaptive_dwr::{AdaptiveDwrPolicy, gain::ComparedDesign};
use crate::{DesignBoxEdge, evaluate_robust_sampled_stress};

fn policy(strict: bool) -> AdaptiveGainPolicy {
    AdaptiveGainPolicy { refinement: AdaptiveDwrPolicy { max_level: 5, max_rounds: 2,
        max_cells: 512, absolute_tolerance: if strict { 1e-30 } else { 1e6 },
        relative_tolerance: 0.0, mark_fraction: 0.25 },
        max_error_fraction: 0.99, area_tolerance: 0.05 }
}
fn owner(max_solves: usize) -> MultiLoadProjectedOptimizer {
    let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = phi.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    let loads = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 4.0], 0.0).unwrap(),
    ];
    MultiLoadProjectedOptimizer::new(phi, &loads, OptimizeSettings { level: 3,
        iterations: 2, volfrac: 0.6, move_cells: 0.1, nucleation_period: 0,
        ..OptimizeSettings::default() }, RobustAggregate::WeightedSum, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        MultiLoadProjectedSettings { max_candidates: 3, max_solves,
            ..MultiLoadProjectedSettings::default() }).unwrap()
        .with_sampled_stress_limit(SampledStressLimit::new(1e8, 0.0).unwrap()).unwrap()
}
fn run(owner: &mut MultiLoadProjectedOptimizer, strict: bool) -> AdaptiveDescentUpdate {
    let ControlFlow::Continue(update) = owner.advance_one_adaptive_dwr_polling(policy(strict), 1,
        |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("not cancelled") };
    update
}

#[test]
fn projected_adaptive_descent_reserves_the_complete_pair_before_proposal_work() {
    // Two baseline slots, then 2 original + 16 adaptive slots per full candidate.
    let mut short = owner(19);
    let before = short.checkpoint_bytes();
    let result = short.advance_one_adaptive_dwr_polling(policy(true), 1,
        |_| -> ControlFlow<()> { panic!("unfunded candidate must not start") }).unwrap();
    assert!(matches!(result, ControlFlow::Continue(AdaptiveDescentUpdate {
        progress: MultiLoadProjectedProgress::SolveBudget(_), .. })));
    assert_eq!(short.checkpoint_bytes(), before);
    for budget in [20, 21] {
        let mut admitted = owner(budget);
        let result = admitted.advance_one_adaptive_dwr_polling(policy(true), 1, |stage| {
            assert!(matches!(stage, AdaptiveDescentStage::Optimizer(MultiLoadProjectedStage::Direction)));
            ControlFlow::Break(17)
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break(17)));
        assert_eq!(admitted.solves_started(), 2);
    }
    for bad in [0, 1] {
        let mut state = owner(1000);
        if bad == 0 { state.kernel.aggregate = RobustAggregate::WorstWeightedCase; }
        else { state.restoration_reduction = Some(0.0); }
        assert!(state.advance_one_adaptive_dwr_polling(policy(false), 1,
            |_| -> ControlFlow<()> { panic!("unsupported mode started work") }).is_err());
    }
}

#[test]
fn projected_adaptive_descent_really_refines_before_refusing_an_unresolved_update() {
    let mut state = owner(1000);
    let current = state.current();
    let stress = state.current_stress().cloned();
    let bits = state.geometry().nodes().to_vec();
    let ell = state.search_multiplier().to_bits();
    let update = run(&mut state, true);
    assert!(matches!(update.progress, MultiLoadProjectedProgress::NoDescent(_)));
    assert!(!update.candidates.is_empty(), "an actual coarse-admissible comparison is required");
    let comparisons: Vec<_> = update.candidates.iter().filter_map(|c| c.report.as_ref()).collect();
    assert!(!comparisons.is_empty(), "at least one real pair must be solved");
    assert!(comparisons.iter().any(|r| r.grid.leaf_count() > 64 && r.rows.len() == 2));
    for report in comparisons { assert_ne!(report.stop, AdaptiveGainStop::ResolvedGain); }
    assert_eq!(state.current(), current);
    assert_eq!(state.current_stress(), stress.as_ref());
    assert_eq!(state.geometry().nodes(), bits);
    assert_eq!(state.search_multiplier().to_bits(), ell);
    assert_eq!(state.next_iteration(), 0);
    assert!(state.solves_started() > 2);
}

#[test]
fn projected_adaptive_descent_cancelled_pair_retains_work_and_restores_the_exact_search() {
    let mut state = owner(1000);
    let current = state.current();
    let stress = state.current_stress().cloned();
    let ell = state.search_multiplier().to_bits();
    let result = state.advance_one_adaptive_dwr_polling(policy(true), 1, |stage| {
        if matches!(stage, AdaptiveDescentStage::Candidate {
            stage: AdaptiveGainStage::BeforeCase { round: 0, design: ComparedDesign::Candidate, case: 1 }, ..
        }) { ControlFlow::Break(17) } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(matches!(result, ControlFlow::Break(17)));
    assert_eq!(state.current(), current);
    assert_eq!(state.current_stress(), stress.as_ref());
    assert_eq!(state.search_multiplier().to_bits(), ell);
    assert_eq!(state.next_iteration(), 0);
    let spent = state.solves_started() - 2;
    assert!(spent >= 8, "original family and three entered DWR cases stay charged");
    let mut recovery = 4;
    let mut restored = MultiLoadProjectedOptimizer::restore_checkpoint(&state.checkpoint_bytes(), &mut recovery).unwrap();
    assert_eq!(recovery, 0);
    assert_eq!(restored.solves_started(), state.solves_started());
    let resumed = run(&mut restored, true);
    let mut fresh = owner(1000);
    let clean = run(&mut fresh, true);
    assert_eq!(format!("{resumed:?}"), format!("{clean:?}"));
    assert_eq!(restored.current(), fresh.current());
    assert_eq!(restored.current_stress(), fresh.current_stress());
    assert_eq!(restored.solves_started(), fresh.solves_started() + spent);
}

#[test]
fn projected_adaptive_descent_acceptance_keeps_original_stress_and_independent_case_binding() {
    let mut state = owner(1000);
    let previous = state.current();
    let update = run(&mut state, false);
    assert!(!update.candidates.is_empty());
    if let MultiLoadProjectedProgress::Accepted(step) = &update.progress {
        let checked = update.candidates.last().unwrap().report.as_ref().unwrap();
        assert_eq!(checked.stop, AdaptiveGainStop::ResolvedGain);
        assert!(checked.rows.last().unwrap().admissible);
        bind(checked, &step.previous, &step.state).unwrap();
        let actual = evaluate_robust_sampled_stress(state.geometry(), state.load_cases(),
            state.settings(), state.aggregate()).unwrap();
        assert_eq!(state.current_stress(), Some(&actual));
        assert!(state.current().objective < previous.objective);
    } else {
        assert_eq!(state.current(), previous);
        assert!(update.candidates.iter().all(|check| check.refusal.is_some()));
    }
    let check = update.candidates.iter().find(|c| c.report.is_some()).unwrap();
    let attempts = match &update.progress {
        MultiLoadProjectedProgress::Accepted(step) => &step.attempts,
        MultiLoadProjectedProgress::NoDescent(attempts)
            | MultiLoadProjectedProgress::SolveBudget(attempts) => attempts,
        MultiLoadProjectedProgress::IterationLimit => panic!("search must run"),
    };
    let original = attempts.iter().find(|a| a.index == check.index)
        .and_then(|a| a.state.as_ref()).unwrap();
    let first = check.report.as_ref().unwrap();
    bind(first, &previous, original).unwrap();
    let mut changed = first.clone();
    changed.initial_case_compliances[0][1] *= 2.0;
    // The zero-weight column cannot be hidden by an unchanged aggregate.
    assert!(bind(&changed, &previous, original).is_err());
}
