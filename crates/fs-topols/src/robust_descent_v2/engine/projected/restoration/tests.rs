use super::*;

fn optimizer(max_solves: usize) -> MultiLoadProjectedOptimizer {
    // A narrow central neck creates a genuine geometry-dependent demand.
    let phi = GridSdf::from_fn(8, &|x, y| (y - 0.5).abs() - (0.15 + 0.3 * (2.0*x - 1.0).powi(2)));
    let fixed = phi.nodes().iter().copied().enumerate().filter(|(i, _)| {
        i % 9 == 0 || i % 9 == 8 || i / 9 == 0 || i / 9 == 8
    }).collect();
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.4375, 0.5625, [0.0, -1.0], 1.0).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.4375, 0.5625, [0.5, 0.0], 0.0).unwrap(),
    ];
    MultiLoadProjectedOptimizer::new(phi, &cases, OptimizeSettings {
        level: 3, iterations: 3, volfrac: 0.5, move_cells: 0.1,
        nucleation_period: 0, ..OptimizeSettings::default()
    }, RobustAggregate::WeightedSum, fixed, VolumeProjectionSettings {
        target: 0.5, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
    }, MultiLoadProjectedSettings { max_solves, max_candidates: 16,
        ..MultiLoadProjectedSettings::default() }).expect("area-feasible neck fixture")
}

fn stress_of(optimizer: &MultiLoadProjectedOptimizer) -> RobustSampledStressEvaluation {
    match optimizer.sample_stress_controlled(&optimizer.current, |_, _| ControlFlow::<()>::Continue(())).unwrap() {
        ControlFlow::Continue(value) => value,
        ControlFlow::Break(()) => panic!("unrequested cancellation"),
    }
}

#[test]
fn g0_restoration_gate_is_strict_and_switches_at_the_exact_bound() {
    assert!(reduces_violation(12.0, 10.0, 10.0, 0.99));
    assert!(reduces_violation(12.0, 11.7, 10.0, 0.1));
    for next in [12.0, 13.0, f64::NAN, f64::INFINITY, -1.0] {
        assert!(!reduces_violation(12.0, next, 10.0, 0.0));
    }
    assert!(!reduces_violation(12.0, 11.0, 10.0, 0.5)); // strict equality boundary
    assert!(!reduces_violation(10.0, 9.0, 10.0, 0.0)); // already feasible
    for bad in [-0.1, 1.0, f64::NAN, f64::INFINITY] { assert!(validate_reduction(bad).is_err()); }
}

#[test]
fn g0_partial_repair_is_not_mislabeled_as_stress_feasibility() {
    let plain = optimizer(2);
    let before = plain.checkpoint_bytes();
    let bound = SampledStressLimit::new(stress_of(&plain).worst_sampled_von_mises * 0.5, 0.0).unwrap();
    let mut restoring = plain.with_stress_restoration(bound, 0.01).unwrap();
    assert!(restoring.is_restoring_stress());
    assert_eq!(restoring.restoration_updates(), 0);
    assert!(matches!(restoring.advance_one().unwrap(), MultiLoadProjectedProgress::SolveBudget(_)));
    assert!(restoring.is_restoring_stress());
    assert_eq!(restoring.solves_started(), 2);
    let state = restoring.current();
    let mut candidate = stress_of(&restoring);
    candidate.worst_sampled_von_mises *= 0.9;
    assert!(restoring.require_candidate(&state, Some(&candidate)).is_ok());
    assert!(candidate.worst_sampled_von_mises > bound.admitted_max());
    assert!(restoring.require_candidate(&state, None).is_err());
    assert_ne!(before, restoring.checkpoint_bytes()); // explicit policy identity
}

#[test]
fn g3_zero_weight_governing_case_has_a_real_unweighted_proposal() {
    let mut run = optimizer(64);
    let settings = run.settings();
    let force = RobustLoadCase::new(DesignBoxEdge::Right, 0.4375, 0.5625, [0.0, -10.0], 0.0).unwrap();
    // Rebuild the declared problem rather than changing a live evaluated load.
    let cases = [run.load_cases()[0], force];
    run = MultiLoadProjectedOptimizer::new(run.geometry().clone(), &cases, settings,
        run.aggregate(), run.fixed_nodes().to_vec(), run.projection_settings(), run.controls()).unwrap();
    let measured = stress_of(&run);
    assert_eq!(measured.worst_stress_case, 1);
    let direction = run.kernel.direction_for_case(&run.current, 0, Some(1)).unwrap();
    assert!(direction.smooth.iter().any(|value| *value > 0.0));
    assert_eq!(run.load_cases()[1].weight(), 0.0);
    assert!(run.kernel.direction_for_case(&run.current, 0, Some(2)).is_err());
}

#[test]
fn g3_actual_restoration_update_reduces_stress_and_replays_from_checkpoint() {
    let plain = optimizer(128);
    let measured = stress_of(&plain).worst_sampled_von_mises;
    let mut run = plain.with_stress_restoration(
        SampledStressLimit::new(measured * (1.0 - 1e-6), 0.0).unwrap(), 0.0,
    ).unwrap();
    let initial = run.checkpoint_bytes();
    let MultiLoadProjectedProgress::Accepted(step) = run.advance_one().unwrap() else {
        panic!("neck fixture must produce an actual stress-reducing update; no zero-step pass");
    };
    assert!(step.restoration);
    assert!(step.stress.as_ref().unwrap().worst_sampled_von_mises < measured);
    assert!((step.state.volume - 0.5).abs() <= 1e-4);
    assert_eq!(run.restoration_updates(), 1);
    let independent = crate::evaluate_robust_sampled_stress(run.geometry(), run.load_cases(),
        run.settings(), run.aggregate()).unwrap();
    assert_eq!(independent, *run.current_stress().unwrap());
    let mut recovery = 4;
    let mut replay = MultiLoadProjectedOptimizer::restore_checkpoint(&initial, &mut recovery).unwrap();
    assert_eq!(recovery, 0);
    assert!(matches!(replay.advance_one().unwrap(), MultiLoadProjectedProgress::Accepted(_)));
    assert_eq!(replay.checkpoint_bytes(), run.checkpoint_bytes());
    let mut recovery = 4;
    let restored = MultiLoadProjectedOptimizer::restore_checkpoint(&run.checkpoint_bytes(), &mut recovery).unwrap();
    assert_eq!(restored.checkpoint_bytes(), run.checkpoint_bytes());
    assert_eq!(restored.is_restoring_stress(), run.is_restoring_stress());
}

#[test]
fn g4_late_stress_sampling_cancellation_preserves_phase_geometry_and_progress() {
    let plain = optimizer(128);
    let bound = SampledStressLimit::new(stress_of(&plain).worst_sampled_von_mises * 0.5, 0.0).unwrap();
    let mut run = plain.with_stress_restoration(bound, 0.0).unwrap();
    let before = run.current();
    let stress = run.current_stress().unwrap().clone();
    let outcome = run.advance_one_controlled(|stage| {
        if matches!(stage, MultiLoadProjectedStage::StressCell { case: 1, .. }) {
            ControlFlow::Break("last-case-sampling")
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(matches!(outcome, ControlFlow::Break("last-case-sampling")));
    assert_eq!(run.current(), before);
    assert_eq!(run.current_stress(), Some(&stress));
    assert_eq!(run.next_iteration(), 0);
    assert_eq!(run.restoration_updates(), 0);
    assert_eq!(run.solves_started(), 4); // cancelled candidate work is not refunded
    let mut recovery = 4;
    let restored = MultiLoadProjectedOptimizer::restore_checkpoint(&run.checkpoint_bytes(), &mut recovery).unwrap();
    assert_eq!(restored.checkpoint_bytes(), run.checkpoint_bytes());
}

#[test]
fn g5_strict_mode_bytes_and_feasible_acceptance_policy_remain_unchanged() {
    let plain = optimizer(64);
    assert_eq!(plain.checkpoint_bytes(), plain.checkpoint_v1_bytes());
    let bound = SampledStressLimit::new(stress_of(&plain).worst_sampled_von_mises * 2.0, 0.0).unwrap();
    let ordinary = optimizer(64).with_sampled_stress_limit(bound).unwrap();
    let restorative = plain.with_stress_restoration(bound, 0.01).unwrap();
    assert!(!restorative.is_restoring_stress());
    let state = restorative.current();
    let stress = restorative.current_stress().unwrap();
    assert!(restorative.require_candidate(&state, Some(stress)).is_err()); // no improvement
    let mut improved = state.clone();
    improved.objective *= 0.9;
    assert!(restorative.require_candidate(&improved, Some(stress)).is_ok());
    let mut overloaded = stress.clone();
    overloaded.worst_sampled_von_mises = bound.admitted_max() * 1.01;
    assert!(restorative.require_candidate(&improved, Some(&overloaded)).is_err());
    assert_eq!(ordinary.checkpoint_bytes(), ordinary.checkpoint_v1_bytes());
    assert!(restorative.with_sampled_stress_limit(bound).is_err());
}

#[test]
fn g0_versioned_restoration_header_refuses_invalid_policy_before_solves() {
    let plain = optimizer(2);
    let mut bytes = MAGIC.to_vec();
    bytes.extend_from_slice(&f64::NAN.to_bits().to_le_bytes());
    bytes.extend_from_slice(&0u64.to_le_bytes());
    bytes.extend_from_slice(&plain.checkpoint_bytes());
    let mut recovery = 4;
    assert!(MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut recovery).is_err());
    assert_eq!(recovery, 4);
    for length in MAGIC.len()..MAGIC.len() + 16 {
        assert!(checkpoint_payload(&bytes[..length]).is_err());
    }
    let mut missing_limit = MAGIC.to_vec();
    missing_limit.extend_from_slice(&0.01f64.to_bits().to_le_bytes());
    missing_limit.extend_from_slice(&0u64.to_le_bytes());
    missing_limit.extend_from_slice(&plain.checkpoint_bytes());
    assert!(MultiLoadProjectedOptimizer::restore_checkpoint(&missing_limit, &mut recovery).is_err());
    assert_eq!(recovery, 4);
}

#[test]
fn g3_refinement_preserves_restoration_policy_without_relaxing_the_limit() {
    let coarse = optimizer(2).with_stress_restoration(SampledStressLimit::new(1e-30, 0.0).unwrap(), 0.01).unwrap();
    let source = coarse.checkpoint_bytes();
    let (fine, _) = crate::refinement::refine_projected_study(&coarse, 1, 2).unwrap();
    assert_eq!(fine.stress_restoration_reduction(), Some(0.01));
    assert_eq!(fine.stress_limit(), coarse.stress_limit());
    assert!(fine.is_restoring_stress());
    assert_eq!(fine.restoration_updates(), 0);
    assert_eq!(coarse.checkpoint_bytes(), source);
}
