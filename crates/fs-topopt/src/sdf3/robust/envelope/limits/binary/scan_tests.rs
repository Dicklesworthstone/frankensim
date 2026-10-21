use super::*;
use crate::sdf3::robust::envelope::tests::fixture;
use crate::{SolveBudget, SolveProgress};
use std::ops::ControlFlow;
fn caps() -> ComplianceCaps3 { ComplianceCaps3::new(vec![1e9; 2], 0.02, 2).unwrap() }

#[test]
fn scan_matches_independent_endpoints_and_filters_only_once() {
    let (mut study, y, z) = fixture(2.0);
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.35 + 0.045*i as f64).collect();
    let thresholds = [0.2, 0.45, 0.55, 0.8];
    let options = BinaryComplianceScanOptions3::default();
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), &thresholds, options, &mut control).unwrap();
    assert_eq!(report.stop, BinaryComplianceScanStop3::ScheduleComplete);
    assert_eq!(report.candidates.len(), thresholds.len());
    assert_eq!(report.work.linear_solves, 1 + 2*thresholds.len());
    let mut selected = None;
    let mut best_volume = f64::INFINITY;
    for (i, &threshold) in thresholds.iter().enumerate() {
        let fresh = study.evaluate_binary_compliance_caps(&raw, &[&y, &z], &caps(), options.endpoint(threshold), &mut control).unwrap();
        assert_eq!(report.candidates[i].case_compliances, fresh.case_compliances);
        assert_eq!(report.candidates[i].dilated_volume, fresh.dilated.volume_fraction);
        if fresh.feasible(options.smooth_tolerance) && fresh.dilated.volume_fraction < best_volume {
            selected = Some(i); best_volume = fresh.dilated.volume_fraction;
        }
    }
    assert_eq!(report.selected_index, selected);
    assert_eq!(report.best.unwrap().dilated.volume_fraction, best_volume);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn equal_volume_uses_first_entry_and_volume_cap_can_reject_all() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.8; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let options = BinaryComplianceScanOptions3::default();
    let report = study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), &[0.1, 0.1], options, &mut control).unwrap();
    assert_eq!(report.selected_index, Some(0));
    let report = study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), &[0.1, 0.1],
        BinaryComplianceScanOptions3 { dilated_volume_limit: 0.0, ..options }, &mut control).unwrap();
    assert_eq!(report.stop, BinaryComplianceScanStop3::ScheduleComplete);
    assert!(report.best.is_none() && report.selected_index.is_none());
    assert!(report.candidates.iter().all(|c| c.physical_caps_satisfied && !c.feasible));
}

#[test]
fn cancelled_scan_retains_the_last_completed_incumbent() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.8; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-binary-cap-scan-accepted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), &[0.1, 0.5],
        BinaryComplianceScanOptions3::default(), &mut control).unwrap();
    assert_eq!(report.stop, BinaryComplianceScanStop3::EvaluationStopped);
    assert_eq!(report.candidates.len(), 1); assert_eq!(report.selected_index, Some(0));
    assert!(report.best.as_ref().unwrap().feasible(1e-6));
    assert_eq!(report.work.linear_solves, 3);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn whole_schedule_is_admitted_before_work_and_unmet_caps_do_not_fall_back_to_gray() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.8; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let options = BinaryComplianceScanOptions3::default();
    for schedule in [&[][..], &[0.2, f64::NAN][..], &[0.2, 0.0][..]] {
        assert!(study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), schedule, options, &mut control).is_err());
    }
    assert!(study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), &[0.2, 0.3],
        BinaryComplianceScanOptions3 { max_candidates: 1, ..options }, &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
    let impossible = ComplianceCaps3::new(vec![1e-10; 2], 0.02, 2).unwrap();
    let report = study.scan_binary_compliance_caps(&raw, &[&y, &z], &impossible, &[0.1, 0.5], options, &mut control).unwrap();
    assert_eq!(report.stop, BinaryComplianceScanStop3::ScheduleComplete);
    assert!(report.best.is_none()); assert!(report.candidates.iter().all(|c| !c.physical_caps_satisfied));
}

#[test]
fn shared_budget_exhaustion_in_next_candidate_keeps_completed_fields() {
    let (mut study, y, z) = fixture(2.0);
    let raw = vec![0.8; study.cells()];
    let scales = study.operator().scales().to_vec();
    let options = BinaryComplianceScanOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let first = study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), &[0.1], options, &mut control).unwrap();
    assert_eq!(first.stop, BinaryComplianceScanStop3::ScheduleComplete);
    let budget = first.work.linear_iterations + 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut poll);
    let report = study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), &[0.1, 0.1], options, &mut control).unwrap();
    assert_eq!(report.stop, BinaryComplianceScanStop3::EvaluationStopped);
    assert!(matches!(report.error, Some(BinaryComplianceError3::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(report.work.linear_iterations, budget);
    assert_eq!(report.candidates, first.candidates);
    assert_eq!(report.selected_index, Some(0));
    assert_eq!(report.best.unwrap().displacements, first.best.unwrap().displacements);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn cancelling_an_unpublished_load_family_keeps_only_complete_candidates() {
    let (mut study, y, z) = fixture(2.0);
    let raw = vec![0.8; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut candidates = 0;
    let mut poll = |p: SolveProgress| {
        if p.stage == "sdf3-binary-cap-scan-candidate" { candidates += 1; }
        if candidates == 2 && p.stage == "sdf3-binary-cap-elasticity" && p.solve_iterations > 0 {
            ControlFlow::Break(())
        } else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = study.scan_binary_compliance_caps(&raw, &[&y, &z], &caps(), &[0.1, 0.5],
        BinaryComplianceScanOptions3::default(), &mut control).unwrap();
    assert_eq!(report.stop, BinaryComplianceScanStop3::EvaluationStopped);
    assert_eq!(report.error, Some(BinaryComplianceError3::Evaluation(EvaluationStop::Cancelled)));
    assert_eq!(report.candidates.len(), 1);
    assert_eq!(report.selected_index, Some(0));
    let retained = report.best.unwrap();
    assert_eq!(retained.displacements.len(), 2);
    assert_eq!(report.candidates[0].case_compliances, retained.case_compliances);
    assert!(report.work.linear_iterations > retained.work.linear_iterations);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn exhaustive_scan_enumerates_every_distinct_three_field_material_state() {
    let (mut study, y, z) = fixture(2.0);
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.31 + 0.052*i as f64).collect();
    let options = BinaryComplianceScanOptions3 { max_candidates: 128, ..Default::default() };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let source = study.binary_cap_source(&raw, options.robustness, &mut control).unwrap();
    let mut expected = Vec::new();
    for design in &source {
        expected.extend(design.projected.iter().copied().filter(|v| *v > 0.0));
    }
    expected.push(1.0);
    expected.sort_by(f64::total_cmp);
    expected.dedup_by(|a,b| *a == *b);
    let report = study.scan_all_binary_compliance_caps(&raw, &[&y,&z], &caps(), options, &mut control).unwrap();
    let actual: Vec<_> = report.candidates.iter().map(|c| c.threshold).collect();
    assert_eq!(actual, expected);
    assert!(actual.windows(2).all(|w| w[0] < w[1]));
    assert_eq!(actual.last(), Some(&1.0));

    let mut limited = options;
    limited.max_candidates = expected.len()-1;
    let before = control.work();
    assert!(study.scan_all_binary_compliance_caps(&raw, &[&y,&z], &caps(), limited, &mut control).is_err());
    assert!(control.work().linear_solves > before.linear_solves);
}
