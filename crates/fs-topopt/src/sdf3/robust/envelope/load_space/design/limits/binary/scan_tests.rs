use super::*;
use crate::sdf3::robust::envelope::tests::fixture;
use crate::{EvaluationStop, SolveBudget, SolveProgress};
use std::ops::ControlFlow;

fn loose() -> LoadSpaceCap3 { LoadSpaceCap3::new(2, 1e9, 0.02).unwrap() }

#[test]
fn scan_replays_each_full_load_set_and_uses_one_source_filter() {
    let (mut study, y, z) = fixture(2.0);
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.35+0.045*i as f64).collect();
    let basis = [&y[..], &z[..]];
    let thresholds = [0.2, 0.45, 0.55, 0.8];
    let options = BinaryLoadSpaceScanOptions3::default();
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = study.scan_binary_load_space_cap(&raw, &basis, &loose(), &thresholds, options, &mut control).unwrap();
    assert_eq!(report.stop, BinaryLoadSpaceScanStop3::ScheduleComplete);
    assert_eq!(report.candidates.len(), thresholds.len());
    assert_eq!(report.work.linear_solves, 1+2*thresholds.len());
    let mut selected = None;
    let mut least = f64::INFINITY;
    for (i, &threshold) in thresholds.iter().enumerate() {
        let fresh = study.evaluate_binary_load_space_cap(&raw, &basis, &loose(), options.endpoint(threshold), &mut control).unwrap();
        let row = &report.candidates[i];
        assert_eq!(row.worst_compliance, fresh.worst_compliance);
        assert_eq!(row.maximum, fresh.spectrum.maximum);
        assert_eq!(row.smooth_maximum, fresh.spectrum.smooth_maximum);
        assert_eq!(row.basis_compliances, fresh.physical.case_compliances);
        assert_eq!(row.dilated_volume, fresh.physical.dilated.volume_fraction);
        let feasible = fresh.feasible(options.smooth_tolerance)
            && fresh.physical.dilated.volume_fraction <= options.dilated_volume_limit;
        assert_eq!(row.feasible, feasible);
        if feasible && row.dilated_volume < least { selected = Some(i); least = row.dilated_volume; }
    }
    assert_eq!(report.selected_index, selected);
    let best = report.best.unwrap();
    assert_eq!(best.physical.dilated.volume_fraction, least);
    let fresh = study.evaluate_binary_load_space_cap(&raw, &basis, &loose(),
        options.endpoint(thresholds[selected.unwrap()]), &mut control).unwrap();
    assert_eq!(best.worst_displacement, fresh.worst_displacement);
    assert_eq!(best.spectrum.gram, fresh.spectrum.gram);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn individually_feasible_columns_cannot_be_selected_as_a_combined_load_pass() {
    let (mut study, y, _) = fixture(2.0);
    let raw = vec![0.8; study.cells()]; let basis = [&y[..], &y[..]];
    let options = BinaryLoadSpaceScanOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let e = study.evaluate_binary_load_space_cap(&raw, &basis, &loose(), options.endpoint(0.1), &mut control).unwrap();
    let cap = LoadSpaceCap3::new(2, 1.25*e.physical.case_compliances[0], 0.02).unwrap();
    let report = study.scan_binary_load_space_cap(&raw, &basis, &cap, &[0.1, 0.2], options, &mut control).unwrap();
    assert_eq!(report.stop, BinaryLoadSpaceScanStop3::ScheduleComplete);
    assert!(report.best.is_none() && report.selected_index.is_none());
    for row in report.candidates {
        assert!(row.basis_compliances.iter().all(|c| *c <= cap.cap()));
        assert!(row.worst_compliance > cap.cap());
        assert!(!row.load_set_feasible && !row.feasible);
    }
}

#[test]
fn tied_material_keeps_first_entry_and_an_independent_volume_cap_can_reject_all() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.8; study.cells()];
    let basis = [&y[..], &z[..]]; let options = BinaryLoadSpaceScanOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = study.scan_binary_load_space_cap(&raw, &basis, &loose(), &[0.1, 0.1], options, &mut control).unwrap();
    assert_eq!(report.selected_index, Some(0));
    let report = study.scan_binary_load_space_cap(&raw, &basis, &loose(), &[0.1, 0.1],
        BinaryLoadSpaceScanOptions3 { dilated_volume_limit: 0.0, ..options }, &mut control).unwrap();
    assert_eq!(report.stop, BinaryLoadSpaceScanStop3::ScheduleComplete);
    assert!(report.candidates.iter().all(|c| c.load_set_feasible && !c.feasible));
    assert!(report.best.is_none() && report.selected_index.is_none());
}

#[test]
fn entire_schedule_basis_and_policy_refuse_before_physical_work() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.6; study.cells()];
    let basis = [&y[..], &z[..]]; let options = BinaryLoadSpaceScanOptions3::default();
    let mut poll = |_| -> ControlFlow<()> { panic!("invalid input reached physical work") };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    for schedule in [&[][..], &[0.2, f64::NAN][..], &[0.2, 0.0][..], &[0.2, 1.1][..]] {
        assert!(study.scan_binary_load_space_cap(&raw, &basis, &loose(), schedule, options, &mut control).is_err());
    }
    for options in [BinaryLoadSpaceScanOptions3 { max_candidates: 1, ..options },
        BinaryLoadSpaceScanOptions3 { smooth_tolerance: f64::NAN, ..options },
        BinaryLoadSpaceScanOptions3 { dilated_volume_limit: -0.1, ..options },
        BinaryLoadSpaceScanOptions3 { max_load_cases: 1, ..options }] {
        assert!(study.scan_binary_load_space_cap(&raw, &basis, &loose(), &[0.2, 0.4], options, &mut control).is_err());
    }
    assert!(study.scan_binary_load_space_cap(&raw, &[&y], &loose(), &[0.2], options, &mut control).is_err());
    assert_eq!(control.work().linear_solves, 0);
}

#[test]
fn a_late_stop_preserves_completed_spectra_not_a_partly_assessed_candidate() {
    for stop_after_completed in [true, false] {
        let (mut study, y, z) = fixture(2.0); let raw = vec![0.8; study.cells()];
        let scales = study.operator().scales().to_vec();
        let mut candidates = 0;
        let mut poll = |p: SolveProgress| {
            if p.stage == "sdf3-binary-load-space-scan-candidate" { candidates += 1; }
            if (stop_after_completed && p.stage == "sdf3-binary-load-space-scan-accepted")
                || (!stop_after_completed && candidates == 2 && p.stage == "sdf3-binary-load-space-spectrum") {
                ControlFlow::Break(())
            } else { ControlFlow::Continue(()) }
        };
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let report = study.scan_binary_load_space_cap(&raw, &[&y, &z], &loose(), &[0.1, 0.5],
            BinaryLoadSpaceScanOptions3::default(), &mut control).unwrap();
        assert_eq!(report.stop, BinaryLoadSpaceScanStop3::EvaluationStopped);
        assert_eq!(report.error, Some(BinaryComplianceError3::Evaluation(EvaluationStop::Cancelled)));
        assert_eq!(report.candidates.len(), 1);
        assert_eq!(report.selected_index, Some(0));
        let best = report.best.unwrap();
        assert!(best.feasible(1e-6));
        assert_eq!(best.worst_compliance, report.candidates[0].worst_compliance);
        assert_eq!(report.work.linear_solves, if stop_after_completed { 3 } else { 5 });
        assert_eq!(study.operator().scales(), scales);
    }
}

#[test]
fn cumulative_krylov_limit_retains_the_complete_first_candidate() {
    let (mut study, y, z) = fixture(2.0); let raw = vec![0.8; study.cells()];
    let basis = [&y[..], &z[..]]; let options = BinaryLoadSpaceScanOptions3::default();
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let first = study.scan_binary_load_space_cap(&raw, &basis, &loose(), &[0.1], options, &mut control).unwrap();
    let budget = first.work.linear_iterations+1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut poll);
    let limited = study.scan_binary_load_space_cap(&raw, &basis, &loose(), &[0.1, 0.1], options, &mut control).unwrap();
    assert_eq!(limited.stop, BinaryLoadSpaceScanStop3::EvaluationStopped);
    assert!(matches!(limited.error, Some(BinaryComplianceError3::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(limited.work.linear_iterations, budget);
    assert_eq!(limited.candidates, first.candidates);
    assert_eq!(limited.selected_index, Some(0));
    let (a,b) = (limited.best.unwrap(), first.best.unwrap());
    assert_eq!(a.spectrum.gram, b.spectrum.gram);
    assert_eq!(a.worst_force, b.worst_force);
    assert_eq!(a.worst_displacement, b.worst_displacement);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn exhaustive_load_space_scan_uses_every_distinct_material_state() {
    let (mut study, y, z) = fixture(2.0);
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.31 + 0.052*i as f64).collect();
    let basis = [&y[..], &z[..]];
    let options = BinaryLoadSpaceScanOptions3 { max_candidates: 128, ..Default::default() };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let source = study.binary_cap_source(&raw, options.robustness, &mut control).unwrap();
    let expected = crate::sdf3::robust::envelope::limits::binary::critical_thresholds(
        &source, options.max_candidates).unwrap();
    let report = study.scan_all_binary_load_space_cap(&raw, &basis, &loose(), options, &mut control).unwrap();
    let actual: Vec<_> = report.candidates.iter().map(|c| c.threshold).collect();
    assert_eq!(actual, expected);
    assert_eq!(actual.last(), Some(&1.0));
    assert!(actual.windows(2).all(|w| w[0] < w[1]));

    for row in &report.candidates {
        let fresh = study.evaluate_binary_load_space_cap(&raw, &basis, &loose(),
            options.endpoint(row.threshold), &mut control).unwrap();
        assert_eq!(row.worst_compliance, fresh.worst_compliance);
        assert_eq!(row.dilated_volume, fresh.physical.dilated.volume_fraction);
        assert_eq!(row.feasible, fresh.feasible(options.smooth_tolerance)
            && fresh.physical.dilated.volume_fraction <= options.dilated_volume_limit);
    }
}
