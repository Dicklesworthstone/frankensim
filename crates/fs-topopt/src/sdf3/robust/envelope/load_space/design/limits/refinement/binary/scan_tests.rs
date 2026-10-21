use super::*;
use super::super::super::tests::{build, trees, params, laws, options};
use crate::{EvaluationStop, SolveBudget, SolveProgress};

fn policy() -> BinaryRefinedScanOptions3 {
    BinaryRefinedScanOptions3 { enrichment:options(), change_tolerance:1.0, ..Default::default() }
}

#[test]
fn g1_selection_matches_independent_two_grid_checks_and_retains_matching_fields() {
    let (tree,fine) = trees();
    let mut study = CutDensityStudy3::new(build(&tree),0.15,SimpParams::default());
    let raw: Vec<_> = (0..study.cells()).map(|i|0.35+0.05*i as f64).collect();
    let incoming = study.params(); let scales = study.operator().scales().to_vec();
    let cap = LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let thresholds = [0.2,0.5,0.65,0.65];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let report = study.scan_binary_load_space_cap_enrichment(&raw,&laws(),params(),&cap,
        &thresholds,policy(),&mut control,|_|Ok(build(&fine))).unwrap();
    assert_eq!(report.stop,BinaryLoadSpaceScanStop3::ScheduleComplete,"{report:?}");
    assert_eq!(report.attempted_candidates,thresholds.len());
    assert_eq!(report.candidates.len(),thresholds.len());
    let mut expected: Option<(usize,BinaryLoadSpaceCapRefinement3)> = None;
    for row in &report.candidates {
        let direct = study.estimate_binary_load_space_cap_enrichment(build(&fine),&laws(),&raw,
            params(),&cap,policy().selection.endpoint(row.threshold),options(),&mut control).unwrap();
        assert_eq!(row.coarse_worst_compliance,direct.coarse.worst_compliance);
        assert_eq!(row.fine_worst_compliance,direct.fine.worst_compliance);
        assert_eq!(row.normalized_gram_change,direct.normalized_gram_change);
        assert_eq!(row.two_grid_passed,direct.stable_and_feasible(policy().selection.smooth_tolerance,1.0));
        if row.feasible && expected.as_ref().is_none_or(|(_,best)|
            direct.physical.dilated.volume_fraction < best.physical.dilated.volume_fraction) {
            expected = Some((row.index,direct));
        }
    }
    let (index,direct) = expected.expect("loose authored cap admits an endpoint");
    assert_eq!(report.selected_index,Some(index));
    assert_ne!(report.selected_index,Some(3)); // equal-material duplicate keeps first
    let best = report.best.unwrap();
    assert_eq!(best.physical.dilated.solid,direct.physical.dilated.solid);
    assert_eq!(best.fine.spectrum.gram,direct.fine.spectrum.gram);
    assert_eq!(best.fine.worst_displacement,direct.fine.worst_displacement);
    assert_eq!(study.operator().scales(),scales);
    assert_eq!(study.params().beta,incoming.beta);
}

#[test]
fn g0_coarse_feasibility_cannot_waive_grid_change_or_material_limits() {
    let (tree,fine) = trees();
    let mut study = CutDensityStudy3::new(build(&tree),0.15,params());
    let raw = vec![0.8;study.cells()];
    let cap = LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let mut p = policy(); p.change_tolerance = 0.0;
    let failed = study.scan_binary_load_space_cap_enrichment(&raw,&laws(),params(),&cap,
        &[0.2],p,&mut control,|_|Ok(build(&fine))).unwrap();
    assert!(failed.candidates[0].coarse_worst_compliance < cap.cap());
    assert!(failed.candidates[0].normalized_gram_change > 0.0);
    assert!(!failed.candidates[0].two_grid_passed);
    assert!(failed.best.is_none() && failed.selected_index.is_none());
    p = policy(); p.selection.dilated_volume_limit = 0.0;
    let material = study.scan_binary_load_space_cap_enrichment(&raw,&laws(),params(),&cap,
        &[0.2],p,&mut control,|_|Ok(build(&fine))).unwrap();
    assert!(material.candidates[0].two_grid_passed);
    assert!(!material.candidates[0].feasible && material.best.is_none());
}

#[test]
fn g4_late_stops_and_failed_builds_keep_only_completed_two_grid_endpoints() {
    let (tree,fine) = trees();
    let mut study = CutDensityStudy3::new(build(&tree),0.15,params());
    let raw = vec![0.8;study.cells()]; let scales = study.operator().scales().to_vec();
    let cap = LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut poll = |p:SolveProgress| if p.stage == "sdf3-binary-refined-scan-accepted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let first = study.scan_binary_load_space_cap_enrichment(&raw,&laws(),params(),&cap,
        &[0.2,0.3],policy(),&mut control,|_|Ok(build(&fine))).unwrap();
    assert_eq!(first.stop,BinaryLoadSpaceScanStop3::EvaluationStopped);
    assert_eq!(first.candidates.len(),1); assert_eq!(first.selected_index,Some(0));
    assert!(matches!(first.error,Some(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let mut builds = 0;
    let failed = study.scan_binary_load_space_cap_enrichment(&raw,&laws(),params(),&cap,
        &[0.2,0.3],policy(),&mut control,|_| {
            builds += 1;
            if builds == 2 { Err(GoalRefinementError3::Invalid("test build refusal")) }
            else { Ok(build(&fine)) }
        }).unwrap();
    assert_eq!(failed.attempted_candidates,2);
    assert_eq!(failed.stop,BinaryLoadSpaceScanStop3::EvaluationStopped);
    assert_eq!(failed.candidates,first.candidates);
    assert_eq!(failed.best.unwrap().fine.worst_displacement,first.best.unwrap().fine.worst_displacement);
    assert_eq!(study.operator().scales(),scales);
}

#[test]
fn g4_invalid_schedules_refuse_before_builds_and_shared_work_is_never_reset() {
    let (tree,fine) = trees();
    let mut study = CutDensityStudy3::new(build(&tree),0.15,params());
    let raw = vec![0.8;study.cells()];
    let cap = LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    for thresholds in [vec![],vec![0.2,f64::NAN],vec![0.0]] {
        assert!(study.scan_binary_load_space_cap_enrichment(&raw,&laws(),params(),&cap,
            &thresholds,policy(),&mut control,|_|panic!("invalid schedule invoked builder")).is_err());
    }
    assert_eq!(control.work().linear_solves,0);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations:1,..Default::default() },&mut poll);
    let stopped = study.scan_binary_load_space_cap_enrichment(&raw,&laws(),params(),&cap,
        &[0.2,0.3],policy(),&mut control,|_|Ok(build(&fine))).unwrap();
    assert_eq!(stopped.stop,BinaryLoadSpaceScanStop3::EvaluationStopped);
    assert_eq!(stopped.work.linear_iterations,1);
    assert!(stopped.candidates.is_empty() && stopped.best.is_none());
}


#[test]
fn g1_exhaustive_refined_scan_matches_every_distinct_source_state() {
    let (tree,fine) = trees();
    let mut study = CutDensityStudy3::new(build(&tree),0.15,SimpParams::default());
    let raw: Vec<_> = (0..study.cells()).map(|i|0.31+0.052*i as f64).collect();
    let cap = LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut p = policy(); p.selection.max_candidates = 128;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);

    let incoming = study.params();
    study.params = params();
    let source = study.binary_cap_source(&raw,p.selection.robustness,&mut control).unwrap();
    study.params = incoming;
    let expected = crate::sdf3::robust::envelope::limits::binary::critical_thresholds(
        &source,p.selection.max_candidates).unwrap();

    let mut builds = 0;
    let report = study.scan_all_binary_load_space_cap_enrichment(
        &raw,&laws(),params(),&cap,p,&mut control,|_| {
            builds += 1; Ok(build(&fine))
        }).unwrap();
    assert_eq!(report.stop,BinaryLoadSpaceScanStop3::ScheduleComplete,"{report:?}");
    assert_eq!(builds,expected.len());
    assert_eq!(report.attempted_candidates,expected.len());
    assert_eq!(report.candidates.iter().map(|c|c.threshold).collect::<Vec<_>>(),expected);
    assert_eq!(expected.last(),Some(&1.0));
    assert!(expected.windows(2).all(|w|w[0]<w[1]));
    for row in &report.candidates {
        let direct = study.estimate_binary_load_space_cap_enrichment(
            build(&fine),&laws(),&raw,params(),&cap,p.selection.endpoint(row.threshold),
            p.enrichment,&mut control).unwrap();
        assert_eq!(row.coarse_worst_compliance,direct.coarse.worst_compliance);
        assert_eq!(row.fine_worst_compliance,direct.fine.worst_compliance);
        assert_eq!(row.normalized_gram_change,direct.normalized_gram_change);
    }
}

#[test]
fn g4_exhaustive_state_limit_refuses_before_any_enriched_build() {
    let (tree,fine) = trees();
    let mut study = CutDensityStudy3::new(build(&tree),0.15,params());
    let raw: Vec<_> = (0..study.cells()).map(|i|0.31+0.052*i as f64).collect();
    let cap = LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut p = policy(); p.selection.max_candidates = 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let before = control.work();
    let result = study.scan_all_binary_load_space_cap_enrichment(
        &raw,&laws(),params(),&cap,p,&mut control,|_| -> Result<AdaptiveElasticity3,GoalRefinementError3> {
            panic!("over-budget exhaustive scan invoked enriched builder")
        });
    assert!(result.is_err());
    assert_eq!(control.work().linear_solves,before.linear_solves+1,
        "only the bounded source-projection filter may run before refusal");
    let _ = fine; // keep fixture symmetry explicit
}
