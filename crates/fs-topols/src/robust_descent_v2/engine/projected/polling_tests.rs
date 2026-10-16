//! G3/G4/G5 tests of real candidate solves and accepted-state continuation.
use super::*;

fn optimizer(aggregate: RobustAggregate, budget: usize) -> MultiLoadProjectedOptimizer {
    let field = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = field.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % 9 == 0 || index % 9 == 8).collect();
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.7).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.5, 0.0], 0.3).unwrap(),
    ];
    MultiLoadProjectedOptimizer::new(
        field, &cases,
        OptimizeSettings {
            level: 3, iterations: 2, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default()
        },
        aggregate, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        MultiLoadProjectedSettings { max_candidates: 8, max_solves: budget, ..Default::default() },
    ).expect("real area-feasible baseline")
}

fn continuing(stage: MultiLoadProjectedStage) -> ControlFlow<()> {
    let _ = stage;
    ControlFlow::Continue(())
}

#[test]
fn polling_batches_preserve_both_aggregates_stress_and_restart_bytes() {
    for aggregate in [RobustAggregate::WeightedSum, RobustAggregate::WorstWeightedCase] {
        for stress_enabled in [false, true] {
            let mut reference = optimizer(aggregate, 256);
            if stress_enabled {
                reference = reference.with_sampled_stress_limit(
                    SampledStressLimit::new(1e8, 0.0).unwrap(),
                ).unwrap();
            }
            let initial = reference.checkpoint_bytes();
            let expected = reference.advance_one().unwrap();
            assert!(matches!(&expected, MultiLoadProjectedProgress::Accepted(_)),
                "fixture must accept a real update: {expected:?}");
            for batch in [1, 7, 64] {
                let mut allowance = 4;
                let mut replay = MultiLoadProjectedOptimizer::restore_checkpoint(&initial, &mut allowance).unwrap();
                assert_eq!(allowance, 0);
                let mut previous = [0; 2];
                let mut polls = 0;
                let result = replay.advance_one_polling(batch, |stage| {
                    if let MultiLoadProjectedStage::CaseIterations { case, iterations, .. } = stage {
                        // Reset at each candidate family; cumulative CG iterations
                        // within a case include all true-residual correction passes.
                        if iterations == 0 { previous[case] = 0; }
                        assert!(iterations >= previous[case]);
                        assert!(iterations - previous[case] <= batch);
                        previous[case] = iterations;
                        polls += 1;
                    }
                    ControlFlow::<()>::Continue(())
                }).unwrap();
                let ControlFlow::Continue(actual) = result else { panic!("unexpected stop") };
                assert!(polls > 2, "must enter the actual CG iteration loop");
                assert_eq!(format!("{actual:?}"), format!("{expected:?}"));
                assert_eq!(replay.current.solutions, reference.current.solutions);
                assert_eq!(replay.current_stress(), reference.current_stress());
                assert_eq!(replay.checkpoint_bytes(), reference.checkpoint_bytes());
            }
        }
    }
}

#[test]
fn cancellation_inside_each_case_discards_partial_family_and_survives_recovery() {
    for stop_case in 0..2 {
        let mut state = optimizer(RobustAggregate::WeightedSum, 256)
            .with_sampled_stress_limit(SampledStressLimit::new(1e8, 0.0).unwrap()).unwrap();
        let current = state.current();
        let nodes = state.geometry().nodes().to_vec();
        let stress = state.current_stress().cloned();
        let ell = state.ell.to_bits();
        let mut saw_iterations = false;
        let stop = state.advance_one_polling(3, |stage| {
            if let MultiLoadProjectedStage::CaseIterations { case, iterations, .. } = stage {
                if case == stop_case && iterations > 0 {
                    saw_iterations = true;
                    return ControlFlow::Break("stopped inside CG");
                }
            }
            ControlFlow::Continue(())
        }).unwrap();
        assert!(saw_iterations);
        assert!(matches!(stop, ControlFlow::Break("stopped inside CG")));
        assert_eq!(state.current(), current);
        assert_eq!(state.geometry().nodes(), nodes);
        assert_eq!(state.current_stress(), stress.as_ref());
        assert_eq!(state.ell.to_bits(), ell);
        assert_eq!(state.next_iteration(), 0);
        assert_eq!(state.solves_started(), 2 + stop_case + 1);
        let checkpoint = state.checkpoint_bytes();
        let mut recovery = 4;
        let mut restored = MultiLoadProjectedOptimizer::restore_checkpoint(&checkpoint, &mut recovery).unwrap();
        assert_eq!(recovery, 0);
        assert_eq!(restored.checkpoint_bytes(), checkpoint);
        let expected = state.advance_one().unwrap();
        let ControlFlow::Continue(actual) = restored.advance_one_polling(7, continuing).unwrap()
            else { panic!("retry refused") };
        assert!(matches!(&actual, MultiLoadProjectedProgress::Accepted(_)));
        assert_eq!(format!("{actual:?}"), format!("{expected:?}"));
        assert_eq!(restored.checkpoint_bytes(), state.checkpoint_bytes());
    }
}

#[test]
fn stress_and_final_publication_stops_keep_the_previous_complete_solution() {
    for stop_at_publication in [false, true] {
        let mut state = optimizer(RobustAggregate::WeightedSum, 256)
            .with_sampled_stress_limit(SampledStressLimit::new(1e8, 0.0).unwrap()).unwrap();
        let before = state.current();
        let before_stress = state.current_stress().cloned();
        let before_fields = state.current.solutions.clone();
        let mut stopped = false;
        let result = state.advance_one_polling(7, |stage| {
            let at_stop = if stop_at_publication {
                matches!(stage, MultiLoadProjectedStage::Publish(_))
            } else {
                matches!(stage, MultiLoadProjectedStage::StressCell { case: 1, cell: 1, .. })
            };
            if at_stop { stopped = true; ControlFlow::Break(()) }
            else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(stopped && matches!(result, ControlFlow::Break(())));
        assert_eq!(state.current(), before);
        assert_eq!(state.current_stress(), before_stress.as_ref());
        assert_eq!(state.current.solutions, before_fields);
        assert_eq!(state.next_iteration(), 0);
        assert!(state.solves_started() >= 4);
    }
}

#[test]
fn restoration_can_cancel_inside_cg_without_losing_its_phase_or_work() {
    let mut state = optimizer(RobustAggregate::WeightedSum, 256)
        .with_stress_restoration(SampledStressLimit::new(1e-30, 0.0).unwrap(), 0.01).unwrap();
    let before = state.current();
    let before_stress = state.current_stress().cloned();
    assert!(state.is_restoring_stress());
    let result = state.advance_one_polling(1, |stage| {
        if matches!(stage, MultiLoadProjectedStage::CaseIterations { iterations: 1.., .. }) {
            ControlFlow::Break(123)
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(matches!(result, ControlFlow::Break(123)));
    assert_eq!(state.current(), before);
    assert_eq!(state.current_stress(), before_stress.as_ref());
    assert_eq!(state.restoration_updates(), 0);
    assert!(state.is_restoring_stress());
    assert_eq!(state.solves_started(), 3);
    let bytes = state.checkpoint_bytes();
    assert!(bytes.starts_with(b"fs-topols/projected-multiload/checkpoint/2\n"));
    let mut allowance = 4;
    let restored = MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut allowance).unwrap();
    assert_eq!(restored.checkpoint_bytes(), bytes);
    assert!(restored.is_restoring_stress());
}

#[test]
fn zero_poll_interval_and_exhausted_solve_budget_do_not_start_work() {
    let mut state = optimizer(RobustAggregate::WeightedSum, 2);
    let before = state.checkpoint_bytes();
    assert!(state.advance_one_polling(0, continuing).is_err());
    let mut polls = 0;
    let result = state.advance_one_polling(1, |_| {
        polls += 1;
        ControlFlow::<()>::Continue(())
    }).unwrap();
    assert!(matches!(result, ControlFlow::Continue(MultiLoadProjectedProgress::SolveBudget(_))));
    assert_eq!(polls, 0);
    assert_eq!(state.checkpoint_bytes(), before);
}

#[test]
fn numerical_refusal_remains_a_refusal_and_is_charged_not_a_partial_solution() {
    let mut state = optimizer(RobustAggregate::WeightedSum, 4);
    // Private unit-test fault: the public material card would reject this
    // before constructing a study. Exercise a failed scheduled assembly.
    state.kernel.material.youngs = f64::NAN;
    let before = state.current();
    let result = state.advance_one_polling(1, continuing).unwrap();
    let ControlFlow::Continue(MultiLoadProjectedProgress::SolveBudget(attempts)) = result
        else { panic!("refused starts must exhaust the fixed allowance") };
    assert_eq!(state.current(), before);
    assert_eq!(state.next_iteration(), 0);
    // One failed first case starts, leaving too little for another full family.
    assert_eq!(state.solves_started(), 3);
    assert_eq!(attempts.len(), 1);
    assert!(attempts[0].state.is_none());
    assert!(attempts[0].refusal.is_some());
}
