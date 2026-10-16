use super::*;
use crate::robust::evaluate_robust_design;

fn loads(scale: f64) -> Vec<RobustLoadCase> {
    vec![
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -scale], 0.7).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.5 * scale, 0.0], 0.3).unwrap(),
    ]
}

fn setup(
    cases: &[RobustLoadCase], aggregate: RobustAggregate, controls: MultiLoadProjectedSettings,
) -> Result<MultiLoadProjectedOptimizer, CutFemError> {
    let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = geometry.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % 9 == 0 || index % 9 == 8).collect();
    MultiLoadProjectedOptimizer::new(
        geometry, cases,
        OptimizeSettings {
            level: 3, iterations: 3, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default()
        },
        aggregate, fixed,
        VolumeProjectionSettings {
            target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
        }, controls,
    )
}

fn optimizer(aggregate: RobustAggregate) -> MultiLoadProjectedOptimizer {
    setup(&loads(1.0), aggregate, MultiLoadProjectedSettings::default()).unwrap()
}

fn require_accepted(progress: MultiLoadProjectedProgress) -> Box<MultiLoadProjectedStep> {
    match progress {
        MultiLoadProjectedProgress::Accepted(step) => step,
        other => panic!("fixture must exercise a real accepted geometry update: {other:?}"),
    }
}

#[test]
fn complete_baseline_matches_the_independent_public_evaluator() {
    for aggregate in [RobustAggregate::WeightedSum, RobustAggregate::WorstWeightedCase] {
        let optimizer = optimizer(aggregate);
        let replay = evaluate_robust_design(
            optimizer.geometry(), &loads(1.0), optimizer.kernel.settings, aggregate,
        ).unwrap();
        assert!((optimizer.baseline().volume - 0.6).abs() <= 1e-4);
        assert_eq!(optimizer.current().snapshot, replay.snapshot);
        assert_eq!(optimizer.current().case_compliances, replay.case_compliances);
        assert_eq!(optimizer.current().objective.to_bits(), replay.objective.to_bits());
        assert_eq!(optimizer.solves_started(), 2);
        assert_eq!(optimizer.next_iteration(), 0);
    }
}

#[test]
fn opposite_cases_do_not_cancel_and_source_scaling_reaches_physics() {
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.5).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 1.0], 0.5).unwrap(),
    ];
    let optimizer = setup(&cases, RobustAggregate::WeightedSum, MultiLoadProjectedSettings::default()).unwrap();
    let state = optimizer.current();
    assert!(state.objective > 0.0);
    assert_eq!(state.case_compliances[0].to_bits(), state.case_compliances[1].to_bits());
    assert_eq!(state.objective.to_bits(), state.case_compliances[0].to_bits());
    let original = setup(&loads(1.0), RobustAggregate::WeightedSum, MultiLoadProjectedSettings::default()).unwrap();
    let scaled = setup(&loads(2.0), RobustAggregate::WeightedSum, MultiLoadProjectedSettings::default()).unwrap();
    assert_eq!(original.geometry().nodes(), scaled.geometry().nodes());
    for (&a, &b) in original.current().case_compliances.iter().zip(&scaled.current().case_compliances) {
        assert!((b - 4.0 * a).abs() <= 1e-8 * a.abs());
    }
}

#[test]
fn accepted_state_is_same_area_and_replays_under_every_case() {
    let mut optimizer = optimizer(RobustAggregate::WeightedSum);
    let baseline = optimizer.baseline().clone();
    let fixed = optimizer.fixed.clone();
    let step = require_accepted(optimizer.advance_one().unwrap());
    assert_eq!(step.previous, baseline);
    assert_eq!(step.iteration, 0);
    assert_eq!(optimizer.next_iteration(), 1);
    assert!(step.state.objective < step.previous.objective);
    assert_ne!(step.state.snapshot, step.previous.snapshot);
    assert!((step.state.volume - 0.6).abs() <= 1e-4);
    for (index, value) in fixed {
        assert_eq!(optimizer.geometry().nodes()[index].to_bits(), value.to_bits());
    }
    let replay = evaluate_robust_design(
        optimizer.geometry(), &loads(1.0), optimizer.kernel.settings, RobustAggregate::WeightedSum,
    ).unwrap();
    assert_eq!(step.state.snapshot, replay.snapshot);
    assert_eq!(step.state.case_compliances, replay.case_compliances);
    assert_eq!(step.state.objective.to_bits(), replay.objective.to_bits());
    assert_eq!(optimizer.baseline(), &baseline);
}

#[test]
fn late_case_cancellation_retains_state_but_spends_work_and_retry_replays() {
    let mut interrupted = optimizer(RobustAggregate::WeightedSum);
    let before = interrupted.current();
    let field = interrupted.geometry().nodes().to_vec();
    let ell = interrupted.ell.to_bits();
    let result = interrupted.advance_one_controlled(|stage| {
        if matches!(stage, MultiLoadProjectedStage::CaseSolve { case: 1, complete: true, .. }) {
            ControlFlow::Break("after last scenario")
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(matches!(result, ControlFlow::Break("after last scenario")));
    assert_eq!(interrupted.solves_started(), 4);
    assert_eq!(interrupted.current(), before);
    assert_eq!(interrupted.geometry().nodes(), field.as_slice());
    assert_eq!(interrupted.ell.to_bits(), ell);
    assert_eq!(interrupted.next_iteration(), 0);
    let resumed = require_accepted(interrupted.advance_one().unwrap());
    let mut clean = optimizer(RobustAggregate::WeightedSum);
    let uninterrupted = require_accepted(clean.advance_one().unwrap());
    assert_eq!(resumed.state, uninterrupted.state);
    assert_eq!(interrupted.geometry().nodes(), clean.geometry().nodes());
    assert_eq!(interrupted.ell.to_bits(), clean.ell.to_bits());
    assert_eq!(interrupted.solves_started(), clean.solves_started() + 2);
}

#[test]
fn publication_cancellation_does_not_commit_an_already_solved_candidate() {
    let mut optimizer = optimizer(RobustAggregate::WeightedSum);
    let before = optimizer.current();
    let geometry = optimizer.geometry().nodes().to_vec();
    let result = optimizer.advance_one_controlled(|stage| {
        if matches!(stage, MultiLoadProjectedStage::Publish(_)) {
            ControlFlow::Break(17)
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(matches!(result, ControlFlow::Break(17)), "must reach the publication boundary");
    assert_eq!(optimizer.current(), before);
    assert_eq!(optimizer.geometry().nodes(), geometry.as_slice());
    assert_eq!(optimizer.next_iteration(), 0);
    assert!(optimizer.solves_started() >= 4);
}

#[test]
fn a_partial_family_allowance_cannot_start_work_or_change_the_baseline() {
    let controls = MultiLoadProjectedSettings { max_solves: 3, ..MultiLoadProjectedSettings::default() };
    let mut optimizer = setup(&loads(1.0), RobustAggregate::WeightedSum, controls).unwrap();
    let before = optimizer.current();
    let result = optimizer.advance_one_controlled(|_| -> ControlFlow<()> {
        panic!("no numerical stage may start without a full-family allowance")
    }).unwrap();
    assert!(matches!(result, ControlFlow::Continue(MultiLoadProjectedProgress::SolveBudget(_))));
    assert_eq!(optimizer.solves_started(), 2);
    assert_eq!(optimizer.current(), before);
    assert_eq!(optimizer.next_iteration(), 0);
}

#[test]
fn rejecting_candidates_spends_the_exact_started_solve_count() {
    let controls = MultiLoadProjectedSettings {
        max_candidates: 2, min_relative_improvement: 0.999_999,
        max_solves: 6, ..MultiLoadProjectedSettings::default()
    };
    let mut optimizer = setup(&loads(1.0), RobustAggregate::WeightedSum, controls).unwrap();
    let before = optimizer.current();
    let mut starts = 0;
    let progress = optimizer.advance_one_controlled(|stage| {
        if matches!(stage, MultiLoadProjectedStage::CaseSolve { complete: false, .. }) { starts += 1; }
        ControlFlow::<()>::Continue(())
    }).unwrap();
    assert!(matches!(progress, ControlFlow::Continue(MultiLoadProjectedProgress::NoDescent(_))));
    assert!(starts >= 2, "must execute at least one complete candidate family");
    assert_eq!(optimizer.solves_started(), 2 + starts);
    assert_eq!(optimizer.current(), before);
}

#[test]
fn zero_weight_late_case_refusal_does_not_produce_a_partial_objective() {
    let geometry = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let cases = [
        loads(1.0)[0],
        RobustLoadCase::new(DesignBoxEdge::Top, 0.375, 0.625, [1.0, 0.0], 0.0).unwrap(),
    ];
    let settings = OptimizeSettings { level: 3, ..OptimizeSettings::default() };
    let kernel = Kernel::new(&geometry, &cases, settings, RobustAggregate::WeightedSum).unwrap();
    let mut boundaries = Vec::new();
    let result = kernel.evaluate_controlled(geometry, |case, complete| {
        boundaries.push((case, complete));
        ControlFlow::<()>::Continue(())
    });
    assert!(result.is_err(), "top-edge traction lies outside material even at zero objective weight");
    assert_eq!(boundaries, vec![(0, false), (0, true), (1, false)]);
}

#[test]
fn controls_that_cannot_admit_a_baseline_refuse() {
    assert!(setup(&loads(1.0), RobustAggregate::WeightedSum,
        MultiLoadProjectedSettings { max_solves: 1, ..MultiLoadProjectedSettings::default() }).is_err());
    assert!(setup(&loads(1.0), RobustAggregate::WeightedSum,
        MultiLoadProjectedSettings { max_candidates: 0, ..MultiLoadProjectedSettings::default() }).is_err());
    assert!(setup(&[], RobustAggregate::WeightedSum, MultiLoadProjectedSettings::default()).is_err());
}
