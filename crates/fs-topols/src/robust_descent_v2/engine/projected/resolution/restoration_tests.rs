use super::*;
use crate::DesignBoxEdge;

fn owner(max_solves: usize) -> MultiLoadProjectedOptimizer {
    let phi = GridSdf::from_fn(8, &|x, y|
        (y - 0.5).abs() - (0.15 + 0.3 * (2.0 * x - 1.0).powi(2)));
    let fixed = phi.nodes().iter().copied().enumerate().filter(|(i, _)|
        i % 9 == 0 || i % 9 == 8 || i / 9 == 0 || i / 9 == 8).collect();
    let loads = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.4375, 0.5625, [0.0, -1.0], 1.0).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.4375, 0.5625, [0.0, -10.0], 0.0).unwrap(),
    ];
    MultiLoadProjectedOptimizer::new(phi, &loads, OptimizeSettings {
        level: 3, iterations: 2, volfrac: 0.5, move_cells: 0.1,
        nucleation_period: 0, ..OptimizeSettings::default()
    }, RobustAggregate::WeightedSum, fixed, VolumeProjectionSettings {
        target: 0.5, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
    }, MultiLoadProjectedSettings { max_solves, max_candidates: 16,
        ..MultiLoadProjectedSettings::default() }).unwrap()
        .with_stress_restoration(SampledStressLimit::new(1e-8, 0.0).unwrap(), 0.0).unwrap()
}
fn policy() -> ResolutionPolicy {
    ResolutionPolicy { extra_levels: 1, absolute_compliance_tolerance: 1e8,
        relative_compliance_tolerance: 0.0, area_tolerance: 0.1 }
}
fn advance(run: &mut MultiLoadProjectedOptimizer) -> MultiLoadMeshProgress {
    match run.advance_one_resolution_polling(policy(), 1, |_| ControlFlow::<()>::Continue(())).unwrap() {
        ControlFlow::Continue(progress) => progress,
        ControlFlow::Break(()) => panic!("unrequested cancellation"),
    }
}

#[test]
fn projected_restoration_mesh_baseline_can_be_overstressed_but_not_underresolved() {
    let mut run = owner(5);
    assert!(run.is_restoring_stress());
    let before = run.current();
    assert!(matches!(advance(&mut run), MultiLoadMeshProgress::Searched {
        progress: MultiLoadProjectedProgress::SolveBudget(_), .. }));
    assert_eq!(run.current(), before);
    assert_eq!(run.solves_started(), 4);
    let mut strict = owner(128);
    let result = strict.advance_one_resolution_polling(ResolutionPolicy {
        absolute_compliance_tolerance: 1e-30, ..policy()
    }, 1, |stage| {
        assert!(!matches!(stage, MultiLoadMeshStage::Optimizer(_)));
        ControlFlow::<()>::Continue(())
    }).unwrap();
    assert!(matches!(result, ControlFlow::Continue(MultiLoadMeshProgress::UnresolvedBaseline { .. })));
    assert_eq!(strict.next_iteration(), 0);
}

#[test]
fn projected_restoration_mesh_real_update_and_checkpoint_replay() {
    let mut run = owner(128);
    assert_eq!(run.current_stress().unwrap().worst_stress_case, 1);
    let initial = run.checkpoint_bytes();
    let MultiLoadMeshProgress::Searched { baseline, candidates,
        progress: MultiLoadProjectedProgress::Accepted(step) } = advance(&mut run)
    else { panic!("neck fixture must actually restore sampled stress on both grids") };
    assert!(step.restoration);
    assert_eq!(run.restoration_updates(), 1);
    assert!(run.is_restoring_stress()); // partial repair is not feasibility
    let candidate = candidates.iter().find(|c| c.refusal.is_none()).unwrap().report.as_ref().unwrap();
    assert!(candidate.restoration_comparison_refusal(&baseline, policy(), 0.5,
        run.stress_limit().unwrap(), run.load_cases(), 0.0).unwrap().is_none());
    let fine = prolongate_level_set(run.geometry(), &[]).unwrap().geometry;
    let independent = crate::evaluate_robust_sampled_stress(&fine, run.load_cases(),
        OptimizeSettings { level: 4, ..run.settings() }, run.aggregate()).unwrap();
    assert_eq!(candidate.rungs[1], rung(4, &independent));
    let mut recovery = 4;
    let mut replay = MultiLoadProjectedOptimizer::restore_checkpoint(&initial, &mut recovery).unwrap();
    assert!(matches!(advance(&mut replay), MultiLoadMeshProgress::Searched {
        progress: MultiLoadProjectedProgress::Accepted(_), .. }));
    assert_eq!(run.checkpoint_bytes(), replay.checkpoint_bytes());
    let mut recovery = 4;
    let restored = MultiLoadProjectedOptimizer::restore_checkpoint(&run.checkpoint_bytes(), &mut recovery).unwrap();
    assert_eq!(restored.checkpoint_bytes(), run.checkpoint_bytes());
}

#[test]
fn projected_restoration_mesh_cancellation_keeps_phase_and_charges_actual_case_starts() {
    for late in [false, true] {
        let mut run = owner(128);
        let before = run.current();
        let stress = run.current_stress().unwrap().clone();
        let result = run.advance_one_resolution_polling(policy(), 1, |stage| {
            let stop = if late {
                matches!(stage, MultiLoadMeshStage::Candidate { stage: MultiLoadResolutionStage::Publish, .. })
            } else {
                matches!(stage, MultiLoadMeshStage::Baseline(MultiLoadResolutionStage::CaseIterations { case: 1, .. }))
            };
            if stop { ControlFlow::Break("stop") } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break("stop")));
        assert_eq!(run.current(), before);
        assert_eq!(run.current_stress(), Some(&stress));
        assert_eq!(run.next_iteration(), 0);
        assert_eq!(run.restoration_updates(), 0);
        assert!(run.is_restoring_stress());
        assert!(run.solves_started() >= if late { 8 } else { 4 });
        let mut recovery = 4;
        let mut restored = MultiLoadProjectedOptimizer::restore_checkpoint(&run.checkpoint_bytes(), &mut recovery).unwrap();
        assert_eq!(restored.checkpoint_bytes(), run.checkpoint_bytes());
        advance(&mut run);
        advance(&mut restored);
        assert_eq!(restored.checkpoint_bytes(), run.checkpoint_bytes());
    }
}
