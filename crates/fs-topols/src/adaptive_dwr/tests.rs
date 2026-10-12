use super::*;
use crate::{DesignBoxEdge, assess_weighted_compliance_dwr};

fn policy() -> AdaptiveDwrPolicy {
    AdaptiveDwrPolicy { max_level: 5, max_rounds: 2, max_cells: 512,
        absolute_tolerance: 1e-30, relative_tolerance: 0.0, mark_fraction: 0.25 }
}
fn fixture() -> (GridSdf, Vec<RobustLoadCase>, OptimizeSettings) {
    (GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35),
        vec![RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 2.0], 0.0).unwrap()],
        OptimizeSettings { level: 3, ..OptimizeSettings::default() })
}
fn run(p: AdaptiveDwrPolicy) -> AdaptiveDwrReport {
    let (phi, cases, settings) = fixture();
    let ControlFlow::Continue(report) = assess_adaptive_dwr_controlled(&phi, &cases, settings, p,
        |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("not cancelled") };
    report
}

#[test]
fn projected_adaptive_dwr_starts_from_canonical_mechanics_and_retains_zero_weight_cases() {
    let (phi, cases, settings) = fixture();
    let before: Vec<_> = phi.nodes().iter().map(|x| x.to_bits()).collect();
    let direct = assess_weighted_compliance_dwr(&phi, &cases, settings).unwrap();
    let p = AdaptiveDwrPolicy { absolute_tolerance: 1e6, ..policy() };
    let ControlFlow::Continue(report) = assess_adaptive_dwr_controlled(&phi, &cases, settings, p,
        |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("not cancelled") };
    assert_eq!(report.stop, AdaptiveDwrStop::EstimatedTolerance);
    assert_eq!(report.rows.len(), 1);
    assert_eq!(report.rows[0].cells, 64);
    assert_eq!(report.rows[0].coarse_compliance.to_bits(), direct.coarse_compliance.to_bits());
    assert_eq!(report.original_area.to_bits(), direct.volume.to_bits());
    assert_eq!(report.snapshot, direct.snapshot);
    assert_eq!(report.solves, 4);
    for (measured, canonical) in report.estimates.iter().zip(&direct.cases) {
        assert_eq!(measured.j_primal.to_bits(), canonical.estimate.j_primal.to_bits());
        assert_eq!(measured.j_enriched.to_bits(), canonical.estimate.j_enriched.to_bits());
    }
    assert!((report.estimates[1].j_primal / report.estimates[0].j_primal - 4.0).abs() < 1e-8);
    assert_eq!(phi.nodes().iter().map(|x| x.to_bits()).collect::<Vec<_>>(), before);
}

#[test]
fn projected_adaptive_dwr_really_refines_locally_and_replays_the_complete_family() {
    let report = run(policy());
    assert_eq!(report.stop, AdaptiveDwrStop::RoundBudget);
    assert_eq!(report.rows.len(), 2);
    assert_eq!(report.solves, 8);
    let first = &report.rows[0]; let last = &report.rows[1];
    assert!(last.cells > first.cells);
    assert!(last.cells < 4 * first.cells, "bulk marking must not masquerade as uniform refinement");
    assert_eq!(last.cells, report.grid.leaf_count());
    assert_eq!(last.max_cell_level, 4);
    assert!(report.grid.leaves().any(|c| c.0 == 3));
    assert!(last.enriched_dofs > last.coarse_dofs);
    assert!(last.relative_residual < 1e-12);
    for cell in report.grid.leaves() {
        for direction in fs_cutfem::FaceDirection::ALL {
            report.grid.face_neighbors(cell, direction).unwrap();
        }
    }
    let replay = run(policy());
    assert_eq!(report.rows, replay.rows);
    assert_eq!(report.grid.leaves().collect::<Vec<_>>(), replay.grid.leaves().collect::<Vec<_>>());
    for (a, b) in report.estimates.iter().zip(&replay.estimates) {
        assert_eq!(a.indicators, b.indicators);
        assert_eq!(a.j_primal.to_bits(), b.j_primal.to_bits());
    }
}

#[test]
fn projected_adaptive_dwr_limits_retain_the_last_actually_solved_grid() {
    for (p, stop) in [
        (AdaptiveDwrPolicy { max_rounds: 1, ..policy() }, AdaptiveDwrStop::RoundBudget),
        (AdaptiveDwrPolicy { max_cells: 64, ..policy() }, AdaptiveDwrStop::CellBudget),
        (AdaptiveDwrPolicy { max_level: 3, ..policy() }, AdaptiveDwrStop::LevelLimit),
    ] {
        let report = run(p);
        assert_eq!(report.stop, stop);
        assert_eq!(report.rows.len(), 1);
        assert_eq!(report.grid.leaf_count(), 64);
        assert_eq!(report.solves, 4);
        assert!(report.estimates.iter().all(|e| e.indicators.keys().all(|c| report.grid.is_leaf(*c))));
    }
}

#[test]
fn projected_adaptive_dwr_cancelled_cases_or_refinement_never_publish_a_partial_report() {
    let (phi, cases, settings) = fixture();
    let bits: Vec<_> = phi.nodes().iter().map(|x| x.to_bits()).collect();
    for stop in [AdaptiveDwrStage::BeforeCase { round: 0, case: 1 },
        AdaptiveDwrStage::BeforeRefine(0), AdaptiveDwrStage::BeforePublish] {
        let result = assess_adaptive_dwr_controlled(&phi, &cases, settings, policy(), |stage| {
            if stage == stop { ControlFlow::Break(17) } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break(17)), "stop {stop:?} must be reached");
        assert_eq!(phi.nodes().iter().map(|x| x.to_bits()).collect::<Vec<_>>(), bits);
    }
}

#[test]
fn projected_adaptive_dwr_marks_signed_magnitudes_and_caps_balance_before_allocation() {
    let grid = Quadtree::with_room(2, 5);
    let indicators = BTreeMap::from([((2, 0, 0), 5.0), ((2, 3, 3), 5.0)]);
    assert_eq!(mark(&grid, &indicators, 0.5).unwrap(), vec![(2, 0, 0)]);
    assert!(refine(&grid, &[(2, 0, 0)], 18).is_none());
    let first = refine(&grid, &[(2, 0, 0)], 19).unwrap();
    assert_eq!(first.leaf_count(), 19);
    // This split touches an unsplit level-2 neighbor. The direct +3 cells fit
    // in 22, but the required balancing split does not; nothing is published.
    assert!(refine(&first, &[(3, 1, 0)], 22).is_none());
    let balanced = refine(&first, &[(3, 1, 0)], 100).unwrap();
    assert!(balanced.leaf_count() > 22);
    assert_eq!(first.leaf_count(), 19);
    assert_eq!(grid.leaf_count(), 16);
}

#[test]
fn projected_adaptive_dwr_invalid_admission_starts_no_solve() {
    let (phi, cases, settings) = fixture();
    for p in [AdaptiveDwrPolicy { max_cells: 63, ..policy() },
        AdaptiveDwrPolicy { max_level: 7, ..policy() },
        AdaptiveDwrPolicy { max_rounds: 0, ..policy() },
        AdaptiveDwrPolicy { mark_fraction: 0.0, ..policy() },
        AdaptiveDwrPolicy { absolute_tolerance: f64::NAN, ..policy() }] {
        assert!(assess_adaptive_dwr_controlled(&phi, &cases, settings, p,
            |_| -> ControlFlow<()> { panic!("invalid admission started work") }).is_err());
    }
    assert_eq!(policy().max_solves(cases.len()), Some(8));
}
