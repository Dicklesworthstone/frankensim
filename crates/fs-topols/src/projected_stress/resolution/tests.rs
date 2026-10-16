use super::*;
use crate::evaluated::DesignEvaluationStage;

fn policy() -> ResolutionPolicy {
    ResolutionPolicy { extra_levels: 1, absolute_compliance_tolerance: 1e6,
        relative_compliance_tolerance: 0.5, area_tolerance: 1.0 }
}
fn area_owner() -> ProjectedOptimizer {
    let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = phi.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    ProjectedOptimizer::new(phi, Cantilever { load: 1.0, band: 0.125 },
        OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default() }, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        ProjectedSettings { max_candidates: 8, poll_iters: 1, ..ProjectedSettings::default() }).unwrap()
}
fn owner() -> ProjectedStressOptimizer {
    ProjectedStressOptimizer::new(area_owner(), SampledStressLimit::new(1e12, 0.0).unwrap()).unwrap()
}
fn assess(geometry: &GridSdf, fixture: Cantilever, settings: OptimizeSettings) -> StressResolutionReport {
    let ControlFlow::Continue(report) = assess_stress_resolution_controlled(geometry,
        fixture, settings, policy(), 1, |_| ControlFlow::<()>::Continue(())).unwrap()
        else { panic!("uninterrupted assessment") };
    report
}
fn unchanged(a: &ProjectedStressOptimizer, b: &ProjectedStressOptimizer) {
    assert_eq!(a.current(), b.current());
    assert_eq!(a.baseline(), b.baseline());
    assert_eq!(a.checkpoint().geometry().nodes(), b.checkpoint().geometry().nodes());
    assert_eq!(a.checkpoint().ell().to_bits(), b.checkpoint().ell().to_bits());
    assert_eq!(a.checkpoint().next_iteration(), b.checkpoint().next_iteration());
}

#[test]
fn projected_stress_resolution_matches_independent_fine_mechanics_and_physical_load_scaling() {
    let optimizer = owner();
    let cp = optimizer.checkpoint();
    let report = assess(cp.geometry(), cp.fixture(), cp.settings());
    assert_eq!(&report.rungs[0].evaluation, optimizer.current());
    let fine = prolongate_level_set(cp.geometry(), &[]).unwrap().geometry;
    let exact = evaluate_sampled_stress(&fine, cp.fixture(), OptimizeSettings { level: 4, ..cp.settings() }).unwrap();
    assert_eq!(report.rungs[1].evaluation, exact);
    assert!(exact.sample_count > report.rungs[0].evaluation.sample_count);
    let doubled = assess(cp.geometry(), Cantilever { load: 2.0, ..cp.fixture() }, cp.settings());
    for (a, b) in report.rungs.iter().zip(&doubled.rungs) {
        assert_eq!(a.evaluation.snapshot, b.evaluation.snapshot);
        assert_eq!(a.evaluation.volume, b.evaluation.volume);
        assert!((b.evaluation.compliance / a.evaluation.compliance - 4.0).abs() < 1e-7);
        assert!((b.evaluation.sampled_max_von_mises / a.evaluation.sampled_max_von_mises - 2.0).abs() < 1e-7);
    }
}

fn measurements(coarse: f64, fine: f64, stress: f64) -> StressResolutionReport {
    StressResolutionReport { rungs: [coarse, fine].iter().enumerate().map(|(i, &compliance)|
        StressResolutionRung { level: 3+i as u32, evaluation: SampledStressEvaluation {
            compliance, volume: 0.6, sampled_max_von_mises: stress, max_location: [0.5, 0.5],
            sample_count: 16*(i+1), snapshot: 42+i as u64,
        }}).collect() }
}

#[test]
fn projected_stress_resolution_rejects_fine_overstress_without_enlarging_the_limit() {
    let baseline = measurements(10.0, 11.0, 9.0);
    let limit = SampledStressLimit::new(10.0, 0.25).unwrap();
    let mut candidate = measurements(9.0, 10.0, 8.0);
    candidate.rungs[1].evaluation.sampled_max_von_mises = 10.25;
    assert!(candidate.comparison_refusal(&baseline, policy(), 0.6, limit, 0.0).unwrap().is_none());
    candidate.rungs[1].evaluation.sampled_max_von_mises = 10.26;
    let why = candidate.comparison_refusal(&baseline, policy(), 0.6, limit, 0.0).unwrap().unwrap();
    assert!(why.contains("level 4 sampled von Mises stress"));
    assert_eq!(limit.admitted_max(), 10.25);
    candidate.rungs[1].evaluation.sampled_max_von_mises = 8.0;
    candidate.rungs[1].evaluation.compliance = 11.1;
    assert!(candidate.comparison_refusal(&baseline, policy(), 0.6, limit, 0.0).unwrap().unwrap().contains("level 4"));
    candidate.rungs.pop();
    assert!(candidate.refusal(policy(), 0.6, limit).is_err());
}

#[test]
fn projected_stress_resolution_unresolved_baseline_never_starts_a_geometry_proposal() {
    let mut optimizer = owner();
    let before = optimizer.clone();
    let strict = ResolutionPolicy { absolute_compliance_tolerance: 1e-14,
        relative_compliance_tolerance: 0.0, ..policy() };
    let result = optimizer.advance_one_resolution_controlled(strict, |stage| {
        assert!(matches!(stage, MeshCheckStage::Baseline(_)));
        ControlFlow::<()>::Continue(())
    }).unwrap();
    let ControlFlow::Continue(MeshCheckedStressProgress::UnresolvedBaseline { baseline, .. }) = result
        else { panic!("the beam must expose its measured discretization change") };
    assert_eq!(baseline.rungs.len(), 2);
    unchanged(&optimizer, &before);
}

#[test]
fn projected_stress_resolution_cancellation_in_fine_cg_or_sampling_preserves_the_owner() {
    for in_sampling in [false, true] {
        let mut optimizer = owner();
        let before = optimizer.clone();
        let mut reached = false;
        let result = optimizer.advance_one_resolution_controlled(policy(), |stage| {
            let stop = match stage {
                MeshCheckStage::Baseline(ResolutionStage::Evaluate { level: 4, stage: DesignEvaluationStage::Solve(n) }) => !in_sampling && n > 0,
                MeshCheckStage::Baseline(ResolutionStage::Evaluate { level: 4, stage: DesignEvaluationStage::StressCell(_) }) => in_sampling,
                _ => false,
            };
            if stop { reached = true; ControlFlow::Break("stop") } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(reached);
        assert!(matches!(result, ControlFlow::Break("stop")));
        unchanged(&optimizer, &before);
    }
}

#[test]
fn projected_stress_resolution_candidate_checks_are_transactional_and_retry_is_exact() {
    let mut optimizer = owner();
    let before = optimizer.clone();
    let mut reached = false;
    let result = optimizer.advance_one_resolution_controlled(policy(), |stage| {
        if matches!(stage, MeshCheckStage::Candidate { stage: ResolutionStage::Publish, .. }) {
            reached = true; ControlFlow::Break("stop")
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(reached, "exercise an actual candidate, not only baseline checks");
    assert!(matches!(result, ControlFlow::Break("stop")));
    unchanged(&optimizer, &before);
    let mut whole = before;
    let a = optimizer.advance_one_resolution_controlled(policy(), |_| ControlFlow::<()>::Continue(())).unwrap();
    let b = whole.advance_one_resolution_controlled(policy(), |_| ControlFlow::<()>::Continue(())).unwrap();
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    unchanged(&optimizer, &whole);
    let ControlFlow::Continue(MeshCheckedStressProgress::Searched { baseline, candidates, update }) = a
        else { panic!("bounded baseline-admitted search must execute") };
    assert!(!candidates.is_empty());
    match update.progress {
        ProjectedProgress::Accepted(_) => {
            let accepted = candidates.last().unwrap().report.as_ref().unwrap();
            assert!(accepted.comparison_refusal(&baseline, policy(), 0.6, optimizer.limit(),
                optimizer.optimizer.controls().min_relative_improvement).unwrap().is_none());
            assert_eq!(&accepted.rungs[0].evaluation, optimizer.current());
        }
        ProjectedProgress::NoDescent(_) => assert_eq!(optimizer.checkpoint().next_iteration(), 0),
        ProjectedProgress::IterationLimit => panic!("fresh owner cannot be complete"),
    }
}

#[test]
fn projected_stress_resolution_validates_policy_and_does_not_solve_a_completed_owner() {
    let mut optimizer = owner();
    let before = optimizer.clone();
    assert!(optimizer.advance_one_resolution_controlled(ResolutionPolicy { extra_levels: 3, ..policy() },
        |_| -> ControlFlow<()> { panic!("invalid policy must refuse before solving") }).is_err());
    unchanged(&optimizer, &before);
    let cp = optimizer.checkpoint();
    let completed = OptimizeCheckpoint::restore(cp.geometry().clone(), cp.fixture(), cp.settings(),
        cp.settings().iterations, cp.ell()).unwrap();
    let mut terminal = ProjectedStressOptimizer::from_checkpoint(completed,
        optimizer.optimizer.fixed_nodes().to_vec(), optimizer.optimizer.projection_settings(),
        optimizer.optimizer.controls(), optimizer.limit()).unwrap();
    assert!(matches!(terminal.advance_one_resolution_controlled(policy(), |_| -> ControlFlow<()> {
        panic!("completed owner must do no assessment work")
    }).unwrap(), ControlFlow::Continue(MeshCheckedStressProgress::IterationLimit)));
}
