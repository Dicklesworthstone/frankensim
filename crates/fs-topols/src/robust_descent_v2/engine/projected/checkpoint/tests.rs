use super::*;

fn optimizer(aggregate: RobustAggregate, stress: bool) -> MultiLoadProjectedOptimizer {
    let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = phi.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % 9 == 0 || index % 9 == 8).collect();
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.7).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.5, 0.0], 0.3).unwrap(),
    ];
    let optimizer = MultiLoadProjectedOptimizer::new(phi, &cases,
        OptimizeSettings { level: 3, iterations: 3, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 2, hole_radius_cells: 1.0, ..OptimizeSettings::default() },
        aggregate, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        MultiLoadProjectedSettings::default(),
    ).unwrap();
    if stress { optimizer.with_sampled_stress_limit(SampledStressLimit::new(1e8, 0.0).unwrap()).unwrap() }
    else { optimizer }
}

fn accepted(optimizer: &mut MultiLoadProjectedOptimizer) {
    assert!(matches!(optimizer.advance_one().unwrap(), MultiLoadProjectedProgress::Accepted(_)),
        "restart fixture must make a real accepted geometry update");
}

#[test]
fn exact_round_trip_preserves_all_policies_and_separately_charges_recovery() {
    for aggregate in [RobustAggregate::WeightedSum, RobustAggregate::WorstWeightedCase] {
        for stress in [false, true] {
            let original = optimizer(aggregate, stress);
            let bytes = original.checkpoint_bytes();
            let mut recovery = 5;
            let restored = MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut recovery).unwrap();
            assert_eq!(recovery, 1);
            assert_eq!(restored.solves_started(), original.solves_started());
            assert_eq!(restored.checkpoint_bytes(), bytes);
            assert_eq!(restored.stress_limit(), original.stress_limit());
            assert_eq!(restored.current_stress(), original.current_stress());
        }
    }
}

#[test]
fn restart_after_accepted_update_reproduces_the_uninterrupted_tail() {
    for stress in [false, true] {
        let mut original = optimizer(RobustAggregate::WeightedSum, stress);
        accepted(&mut original);
        let initial = original.baseline_geometry().nodes().to_vec();
        assert_ne!(original.geometry().nodes(), initial);
        let bytes = original.checkpoint_bytes();
        let mut recovery = 4;
        let mut resumed = MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut recovery).unwrap();
        assert_eq!(resumed.next_iteration(), 1);
        assert_eq!(resumed.baseline_geometry().nodes(), initial);
        // Continue through the global nucleation ordinal rather than restarting it.
        for _ in 0..3 {
            let a = original.advance_one().unwrap();
            let b = resumed.advance_one().unwrap();
            assert_eq!(format!("{a:?}"), format!("{b:?}"));
            assert_eq!(resumed.checkpoint_bytes(), original.checkpoint_bytes());
            if !matches!(a, MultiLoadProjectedProgress::Accepted(_)) { break; }
        }
    }
}

#[test]
fn malformed_or_underfunded_recovery_starts_no_pde_work() {
    let original = optimizer(RobustAggregate::WeightedSum, false);
    let bytes = original.checkpoint_bytes();
    for length in [0, MAGIC.len() - 1, MAGIC.len(), 100, bytes.len() - 1] {
        let mut budget = 4;
        assert!(MultiLoadProjectedOptimizer::restore_checkpoint(&bytes[..length], &mut budget).is_err());
        assert_eq!(budget, 4);
    }
    let mut extra = bytes.clone(); extra.push(0);
    let mut budget = 4;
    assert!(MultiLoadProjectedOptimizer::restore_checkpoint(&extra, &mut budget).is_err());
    assert_eq!(budget, 4);
    let mut budget = 3;
    assert!(MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut budget).is_err());
    assert_eq!(budget, 3);
}

#[test]
fn corrupt_numerical_evidence_refuses_after_replay_without_refunding_work() {
    let original = optimizer(RobustAggregate::WeightedSum, true);
    let mut bytes = original.checkpoint_bytes();
    // Corrupt the retained stress snapshot, leaving declarations and fields valid.
    *bytes.last_mut().unwrap() ^= 1;
    let mut budget = 6;
    let result = MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut budget);
    assert!(result.is_err());
    assert_eq!(budget, 2);
    assert_eq!(original.solves_started(), 2);
}

#[test]
fn interrupted_candidate_work_survives_checkpoint_recovery() {
    let mut original = optimizer(RobustAggregate::WeightedSum, false);
    let before = original.current();
    let mut seen = false;
    let stopped = original.advance_one_controlled(|stage| {
        if matches!(stage, MultiLoadProjectedStage::CaseSolve { case: 1, complete: true, .. }) {
            seen = true; ControlFlow::Break("pause")
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(seen && matches!(stopped, ControlFlow::Break("pause")));
    assert_eq!(original.current(), before);
    assert!(original.solves_started() > 2);
    let bytes = original.checkpoint_bytes();
    let mut budget = 4;
    let resumed = MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut budget).unwrap();
    assert_eq!(resumed.checkpoint_bytes(), bytes);
    assert_eq!(resumed.solves_started(), original.solves_started());
    assert_eq!(resumed.next_iteration(), 0);
}
