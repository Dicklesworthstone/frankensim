use super::*;
use crate::projected::ProjectedSettings;
use crate::volume::VolumeProjectionSettings;

fn policy() -> ResolutionPolicy {
    ResolutionPolicy { extra_levels: 1, absolute_compliance_tolerance: 1e6,
        relative_compliance_tolerance: 0.0, area_tolerance: 0.5 }
}
fn owner() -> ProjectedOptimizer {
    let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = phi.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    let settings = OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6,
        move_cells: 0.1, nucleation_period: 0, ..OptimizeSettings::default() };
    ProjectedOptimizer::new(phi, Cantilever { load: 1.0, band: 0.125 }, settings, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        ProjectedSettings { max_candidates: 8, poll_iters: 1, ..ProjectedSettings::default() }).unwrap()
}
fn assess(owner: &ProjectedOptimizer, load: f64) -> ResolutionReport {
    let ControlFlow::Continue(report) = assess_compliance_resolution_controlled(
        owner.checkpoint().geometry(), Cantilever { load, ..owner.checkpoint().fixture() },
        owner.checkpoint().settings(), policy(), 1, |_| ControlFlow::<()>::Continue(())).unwrap()
    else { panic!("uninterrupted assessment") };
    report
}
fn unchanged(before: &ProjectedOptimizer, after: &ProjectedOptimizer) {
    assert_eq!(before.checkpoint().geometry().nodes(), after.checkpoint().geometry().nodes());
    assert_eq!(before.checkpoint().next_iteration(), after.checkpoint().next_iteration());
    assert_eq!(before.checkpoint().ell().to_bits(), after.checkpoint().ell().to_bits());
    assert_eq!(before.current().compliance.to_bits(), after.current().compliance.to_bits());
}

#[test]
fn projected_resolution_rungs_are_actual_same_geometry_solves_and_loads_scale() {
    let owner = owner();
    let before = owner.checkpoint().geometry().nodes().to_vec();
    let report = assess(&owner, 1.0);
    assert_eq!(report.rungs.len(), 2);
    let fine = prolongate_level_set(owner.checkpoint().geometry(), &[]).unwrap().geometry;
    let oracle = crate::evaluate_compliance_design(&fine, owner.checkpoint().fixture(),
        OptimizeSettings { level: 4, ..owner.checkpoint().settings() }).unwrap();
    assert_eq!(report.rungs[1].snapshot, oracle.snapshot);
    assert_eq!(report.rungs[1].compliance.to_bits(), oracle.compliance.to_bits());
    assert_eq!(report.rungs[1].area.to_bits(), oracle.volume.to_bits());
    let doubled = assess(&owner, 2.0);
    for (a, b) in report.rungs.iter().zip(&doubled.rungs) {
        assert!((b.compliance / a.compliance - 4.0).abs() < 1e-8);
        assert_eq!(a.area.to_bits(), b.area.to_bits());
        assert_eq!(a.snapshot, b.snapshot);
    }
    assert_eq!(owner.checkpoint().geometry().nodes(), before);
}

#[test]
fn projected_resolution_tight_tolerance_stops_before_any_geometry_proposal() {
    let mut owner = owner();
    let before = owner.clone();
    let strict = ResolutionPolicy { absolute_compliance_tolerance: f64::MIN_POSITIVE, ..policy() };
    let ControlFlow::Continue(progress) = owner.advance_one_resolution_controlled(strict, |stage| {
        assert!(matches!(stage, MeshCheckStage::Baseline(_)));
        ControlFlow::<()>::Continue(())
    }).unwrap() else { panic!("uninterrupted check") };
    let MeshCheckedProgress::UnresolvedBaseline { baseline, .. } = progress
        else { panic!("the cut-beam grid difference must be nonzero") };
    assert!(baseline.max_compliance_change > strict.absolute_compliance_tolerance);
    assert!(!baseline.agrees);
    unchanged(&before, &owner);
}

#[test]
fn projected_resolution_cancellation_during_fine_cg_leaves_state_unchanged() {
    let mut owner = owner();
    let before = owner.clone();
    let result = owner.advance_one_resolution_controlled(policy(), |stage| {
        if matches!(stage, MeshCheckStage::Baseline(ResolutionStage::Evaluate {
            level: 4, stage: DesignEvaluationStage::Solve(n) }) if n > 0) {
            ControlFlow::Break("fine CG cancelled")
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(matches!(result, ControlFlow::Break("fine CG cancelled")));
    unchanged(&before, &owner);
}

#[test]
fn projected_resolution_cancel_after_candidate_assessment_does_not_publish() {
    let mut owner = owner();
    let before = owner.clone();
    let result = owner.advance_one_resolution_controlled(policy(), |stage| {
        if matches!(stage, MeshCheckStage::Candidate { stage: ResolutionStage::Publish, .. }) {
            ControlFlow::Break(71)
        } else { ControlFlow::Continue(()) }
    }).unwrap();
    assert!(matches!(result, ControlFlow::Break(71)));
    unchanged(&before, &owner);
    let mut fresh = before;
    let a = owner.advance_one_resolution_controlled(policy(), |_| ControlFlow::<()>::Continue(())).unwrap();
    let b = fresh.advance_one_resolution_controlled(policy(), |_| ControlFlow::<()>::Continue(())).unwrap();
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    unchanged(&fresh, &owner);
}

fn measured(values: [f64; 2], areas: [f64; 2]) -> ResolutionReport {
    summarize((0..2).map(|i| ResolutionRung { level: 3+i as u32,
        compliance: values[i], area: areas[i], snapshot: i as u64 }).collect(), policy()).unwrap()
}

#[test]
fn projected_resolution_rejects_coarse_only_wins_and_added_material() {
    let baseline = measured([100.0, 99.0], [0.6; 2]);
    let reversed = measured([90.0, 101.0], [0.6; 2]);
    assert!(candidate_refusal(&baseline, &reversed, policy(), 0.6, 0.0).unwrap().contains("level 4"));
    let better = measured([90.0, 89.0], [0.6; 2]);
    assert!(candidate_refusal(&baseline, &better, policy(), 0.6, 0.01).is_none());
    let added = measured([90.0, 89.0], [0.6, 0.65]);
    assert!(candidate_refusal(&baseline, &added,
        ResolutionPolicy { area_tolerance: 0.01, ..policy() }, 0.6, 0.0).is_some());
    let mut wrong_grid = better;
    wrong_grid.rungs[1].level = 5;
    assert!(candidate_refusal(&baseline, &wrong_grid, policy(), 0.6, 0.0).is_some());
}

#[test]
fn projected_resolution_policy_refuses_before_work_and_checks_all_declared_levels() {
    for p in [ResolutionPolicy { extra_levels: 0, ..policy() },
        ResolutionPolicy { extra_levels: 3, ..policy() },
        ResolutionPolicy { area_tolerance: f64::NAN, ..policy() },
        ResolutionPolicy { absolute_compliance_tolerance: 0.0, ..policy() },
        ResolutionPolicy { relative_compliance_tolerance: 1.0, ..policy() }] {
        assert!(p.validate(3).is_err());
    }
    assert!(policy().validate(7).is_err());
    let p = ResolutionPolicy { extra_levels: 2, absolute_compliance_tolerance: 1.0, ..policy() };
    let report = summarize(vec![ResolutionRung { level: 2, compliance: 100.0, area: 0.6, snapshot: 1 },
        ResolutionRung { level: 3, compliance: 90.0, area: 0.6, snapshot: 2 },
        ResolutionRung { level: 4, compliance: 90.5, area: 0.6, snapshot: 3 }], p).unwrap();
    assert!(!report.agrees);
    assert_eq!(report.max_compliance_change, 10.0);
}
