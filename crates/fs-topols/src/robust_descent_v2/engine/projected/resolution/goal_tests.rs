use super::*;
use crate::DesignBoxEdge;

fn loads() -> Vec<RobustLoadCase> {
    vec![RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 4.0], 0.0).unwrap()]
}
fn policy() -> MultiLoadDwrPolicy { MultiLoadDwrPolicy { max_error_fraction: 0.9 } }
fn owner(limit: usize) -> MultiLoadProjectedOptimizer {
    let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = geometry.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    MultiLoadProjectedOptimizer::new(geometry, &loads(), OptimizeSettings {
        level: 3, iterations: 2, volfrac: 0.6, move_cells: 0.1, nucleation_period: 0,
        ..OptimizeSettings::default()
    }, RobustAggregate::WeightedSum, fixed, VolumeProjectionSettings {
        target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
    }, MultiLoadProjectedSettings { max_candidates: 6, max_solves: limit,
        ..MultiLoadProjectedSettings::default() }).unwrap()
        .with_sampled_stress_limit(SampledStressLimit::new(1e8, 0.0).unwrap()).unwrap()
}
fn run(owner: &mut MultiLoadProjectedOptimizer) -> MultiLoadDwrProgress {
    let ControlFlow::Continue(result) = owner.advance_one_dwr_polling(policy(), 1,
        |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("not cancelled") };
    result
}

#[test]
fn weighted_dwr_real_baseline_preserves_independent_zero_weight_equilibria() {
    let mut owner = owner(6); // Two original + four DWR slots; no candidate fits.
    let current = owner.current();
    let MultiLoadDwrProgress::Searched { baseline, candidates, progress } = run(&mut owner)
        else { panic!("baseline fits") };
    assert!(matches!(progress, MultiLoadProjectedProgress::SolveBudget(_)));
    assert!(candidates.is_empty());
    baseline.matches_state(&loads(), &current, 3).unwrap();
    assert!((baseline.cases[1].1.coarse_compliance / baseline.cases[0].1.coarse_compliance - 16.0).abs() < 1e-7);
    assert!((baseline.cases[1].1.enriched_compliance / baseline.cases[0].1.enriched_compliance - 16.0).abs() < 1e-7);
    assert_eq!(owner.solves_started(), 6);
    assert_eq!(owner.current(), current);
    assert_eq!(owner.next_iteration(), 0);
}

#[test]
fn weighted_dwr_budget_and_unsupported_modes_refuse_before_work() {
    let mut short = owner(5);
    assert!(matches!(short.advance_one_dwr_polling(policy(), 1,
        |_| -> ControlFlow<()> { panic!("unfunded baseline") }).unwrap(),
        ControlFlow::Continue(MultiLoadDwrProgress::SolveBudget)));
    assert_eq!(short.solves_started(), 2);
    short.kernel.aggregate = RobustAggregate::WorstWeightedCase;
    assert!(short.advance_one_dwr_polling(policy(), 1,
        |_| -> ControlFlow<()> { panic!("unsupported objective") }).is_err());
    short.kernel.aggregate = RobustAggregate::WeightedSum;
    short.restoration_reduction = Some(0.0);
    assert!(short.advance_one_dwr_polling(policy(), 1,
        |_| -> ControlFlow<()> { panic!("unsupported restoration") }).is_err());
}

#[test]
fn weighted_dwr_cancelled_case_preserves_state_and_charges_lifetime_work() {
    for stop in 0..3 {
        let mut interrupted = owner(1000);
        let current = interrupted.current();
        let stress = interrupted.current_stress().unwrap().clone();
        let bits: Vec<_> = interrupted.geometry().nodes().iter().map(|x| x.to_bits()).collect();
        let multiplier = interrupted.search_multiplier().to_bits();
        let result = interrupted.advance_one_dwr_polling(policy(), 1, |stage| {
            let hit = match (stop, stage) {
                (0, MultiLoadDwrStage::Prepare) => true,
                (1, MultiLoadDwrStage::Baseline(WeightedComplianceDwrStage::BeforeCase { case: 1 })) => true,
                (2, MultiLoadDwrStage::Candidate { stage: WeightedComplianceDwrStage::BeforePublish, .. }) => true,
                _ => false,
            };
            if hit { ControlFlow::Break(17) } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break(17)), "real checkpoint must be reached");
        assert_eq!(interrupted.current(), current);
        assert_eq!(interrupted.current_stress(), Some(&stress));
        assert_eq!(interrupted.geometry().nodes().iter().map(|x| x.to_bits()).collect::<Vec<_>>(), bits);
        assert_eq!(interrupted.search_multiplier().to_bits(), multiplier);
        assert_eq!(interrupted.next_iteration(), 0);
        let spent = interrupted.solves_started() - 2;
        if stop == 0 { assert_eq!(spent, 0); } else { assert!(spent > 0); }
        let retry = run(&mut interrupted);
        let mut fresh = owner(1000);
        let clean = run(&mut fresh);
        assert_eq!(format!("{retry:?}"), format!("{clean:?}"));
        assert_eq!(interrupted.current(), fresh.current());
        assert_eq!(interrupted.current_stress(), fresh.current_stress());
        assert_eq!(interrupted.solves_started(), fresh.solves_started() + spent);
    }
}
