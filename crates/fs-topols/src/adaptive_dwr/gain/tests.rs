use super::*;
use crate::DesignBoxEdge;

fn policy() -> AdaptiveGainPolicy {
    AdaptiveGainPolicy { refinement: AdaptiveDwrPolicy { max_level: 5, max_rounds: 2,
        max_cells: 512, absolute_tolerance: 1e-30, relative_tolerance: 0.0, mark_fraction: 0.25 },
        max_error_fraction: 0.5, area_tolerance: 1e-3 }
}
fn fixture() -> (GridSdf, GridSdf, Vec<RobustLoadCase>, OptimizeSettings) {
    (GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35),
        GridSdf::from_fn(8, &|x, y| (y - 0.5).abs() - (0.35 + 0.05 * (x - 0.5))),
        vec![RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 2.0], 0.0).unwrap()],
        OptimizeSettings { level: 3, volfrac: 0.7, ..OptimizeSettings::default() })
}
fn run(p: AdaptiveGainPolicy) -> AdaptiveGainReport {
    let (a, b, cases, settings) = fixture();
    let ControlFlow::Continue(report) = compare_adaptive_designs_controlled(
        &a, &b, &cases, settings, p, 0.0, |_| ControlFlow::<()>::Continue(())).unwrap()
        else { panic!("not cancelled") };
    report
}
fn measurement(coarse: f64, enriched: f64) -> AdaptiveDwrRow {
    AdaptiveDwrRow { cells: 64, max_cell_level: 3, coarse_compliance: coarse,
        enriched_compliance: enriched, eta_signed: 0.0, correction_magnitude: 0.0,
        observed_change: 0.0, tolerance: 10.0, coarse_dofs: 100, enriched_dofs: 300,
        relative_residual: 1e-14 }
}

#[test]
fn projected_adaptive_gain_requires_both_area_and_resolution_even_with_a_large_gain() {
    let check = |old, new, areas, c, e| row(old, new, [0.6; 2], areas,
        c, e, 0.6, 0.0, policy()).unwrap();
    let old = measurement(10.0, 11.0);
    let new = measurement(9.0, 10.0);
    assert!(check(old.clone(), new.clone(), [0.6; 2], 0.5, 0.5).admissible);
    assert!(!check(old.clone(), new.clone(), [0.6; 2], 0.5f64.next_up(), 0.0).admissible);
    let bad_area = check(old.clone(), new.clone(), [0.6, 0.602], 0.0, 0.0);
    assert!(bad_area.decreasing && !bad_area.areas_match && !bad_area.admissible);
    let unresolved = check(AdaptiveDwrRow { correction_magnitude: 11.0, ..old.clone() },
        new, [0.6; 2], 0.0, 0.0);
    assert!(!unresolved.endpoints_resolved && !unresolved.admissible);
    assert!(!check(old.clone(), old, [0.6; 2], 0.0, 0.0).admissible);
}

#[test]
fn projected_adaptive_gain_independent_error_changes_cannot_cancel() {
    let (correction, observed) = changes([
        (1.0, [10.0, 11.0, 1.0], [9.0, 12.0, 3.0]),
        (1.0, [10.0, 11.0, 1.0], [9.0, 8.0, -1.0]),
    ].into_iter()).unwrap();
    assert_eq!(correction, 4.0);
    assert_eq!(observed, 4.0);
    let measured = row(measurement(20.0, 22.0), measurement(18.0, 20.0),
        [0.6; 2], [0.6; 2], correction, observed, 0.6, 0.0, policy()).unwrap();
    assert!(measured.decreasing);
    assert!(!measured.admissible);
    assert!(changes([(1.0, [1.0, 1.0, f64::MAX], [1.0, 1.0, -f64::MAX])].into_iter()).is_err());
    assert!(changes([(0.0, [1.0, f64::NAN, 0.0], [1.0; 3])].into_iter()).is_err());
}

#[test]
fn projected_adaptive_gain_resolves_both_designs_on_the_actual_shared_mixed_grid() {
    let report = run(policy());
    assert_eq!(report.stop, AdaptiveGainStop::RoundBudget);
    assert_eq!(report.rows.len(), 2);
    assert_eq!(report.solves, 16);
    assert!(report.grid.leaf_count() > 64 && report.grid.leaf_count() < 256);
    assert!(report.grid.leaves().any(|c| c.0 == 3));
    assert!(report.grid.leaves().any(|c| c.0 == 4));
    let (a, b, cases, settings) = fixture();
    let material = prepare(&a, &cases, settings, policy().refinement).unwrap();
    let initial = Quadtree::with_room(3, 5);
    for (side, (geometry, estimates)) in [(&a, &report.baseline_estimates),
        (&b, &report.candidate_estimates)].into_iter().enumerate() {
        assert_eq!(report.snapshots[side], snapshot(geometry));
        for (index, case) in cases.iter().enumerate() {
            let first = estimate_case(&initial, geometry, &material, *case).unwrap();
            assert_eq!(report.initial_case_compliances[side][index].to_bits(), first.j_primal.to_bits());
            let replay = estimate_case(&report.grid, geometry, &material, *case).unwrap();
            assert_eq!(replay.indicators, estimates[index].indicators);
            assert_eq!(replay.j_primal.to_bits(), estimates[index].j_primal.to_bits());
        }
        assert!((estimates[1].j_enriched / estimates[0].j_enriched - 4.0).abs() < 1e-7,
            "opposite zero-weight load retains its independently solved response");
    }
    for cell in report.grid.leaves() {
        for direction in fs_cutfem::FaceDirection::ALL {
            report.grid.face_neighbors(cell, direction).unwrap();
        }
    }
    assert_eq!(report.rows, run(policy()).rows);
}

#[test]
fn projected_adaptive_gain_budget_stops_retain_a_fully_solved_pair() {
    for (refinement, expected) in [
        (AdaptiveDwrPolicy { max_cells: 64, ..policy().refinement }, AdaptiveGainStop::CellBudget),
        (AdaptiveDwrPolicy { max_rounds: 1, ..policy().refinement }, AdaptiveGainStop::RoundBudget),
        (AdaptiveDwrPolicy { max_level: 3, ..policy().refinement }, AdaptiveGainStop::LevelLimit),
    ] {
        let report = run(AdaptiveGainPolicy { refinement, ..policy() });
        assert_eq!(report.stop, expected);
        assert_eq!(report.rows.len(), 1);
        assert_eq!(report.solves, 8);
        assert_eq!(report.grid.leaf_count(), 64);
        assert!(report.baseline_estimates.iter().chain(&report.candidate_estimates)
            .all(|e| e.indicators.keys().all(|c| report.grid.is_leaf(*c))));
    }
}

#[test]
fn projected_adaptive_gain_cancellation_returns_no_partial_pair_and_retries_exactly() {
    let (a, b, cases, settings) = fixture();
    let bits = |g: &GridSdf| g.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    let before = [bits(&a), bits(&b)];
    for stop in [AdaptiveGainStage::BeforeCase { round: 0, design: ComparedDesign::Candidate, case: 1 },
        AdaptiveGainStage::BeforeRefine(0), AdaptiveGainStage::BeforePublish] {
        let mut debit = 0;
        let result = compare_adaptive_designs_controlled(&a, &b, &cases, settings, policy(), 0.0, |stage| {
            if stage == stop { return ControlFlow::Break(17); }
            if matches!(stage, AdaptiveGainStage::BeforeCase { .. }) { debit += 2; }
            ControlFlow::Continue(())
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break(17)));
        assert!(debit > 0 && debit <= policy().max_solves(cases.len()).unwrap());
        assert_eq!([bits(&a), bits(&b)], before);
    }
    assert_eq!(run(policy()).rows, run(policy()).rows);
}

#[test]
fn projected_adaptive_gain_invalid_controls_refuse_before_any_callback() {
    let (a, b, cases, settings) = fixture();
    for bad in [AdaptiveGainPolicy { max_error_fraction: f64::NAN, ..policy() },
        AdaptiveGainPolicy { max_error_fraction: 1.0, ..policy() },
        AdaptiveGainPolicy { area_tolerance: -1.0, ..policy() },
        AdaptiveGainPolicy { refinement: AdaptiveDwrPolicy { max_cells: 63,
            ..policy().refinement }, ..policy() }] {
        assert!(compare_adaptive_designs_controlled(&a, &b, &cases, settings, bad, 0.0,
            |_| -> ControlFlow<()> { panic!("invalid controls started work") }).is_err());
    }
    assert_eq!(policy().max_solves(2), Some(16));
}
