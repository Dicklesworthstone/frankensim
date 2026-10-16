use super::*;
use crate::robust_descent::{MultiLoadProjectedProgress, MultiLoadProjectedSettings};
use crate::volume::VolumeProjectionSettings;
use crate::{OptimizeSettings, RobustAggregate, RobustLoadCase, SampledStressLimit};
use fs_cutfem::DesignBoxEdge;

#[test]
fn dyadic_transfer_preserves_old_bits_and_the_piecewise_bilinear_field() {
    let mut coarse = GridSdf::from_fn(8, &|x, y| (7.0 * x).sin() + x * y - y * y);
    *coarse.node_mut(3, 2) = -0.0;
    let fine = prolongate_level_set(&coarse, &[]).unwrap().geometry;
    for j in 0..=8 {
        for i in 0..=8 {
            assert_eq!(coarse.node(i, j).to_bits(), fine.node(2 * i, 2 * j).to_bits());
        }
    }
    for j in 0..=37 {
        for i in 0..=41 {
            let point = [f64::from(i) / 41.0, f64::from(j) / 37.0];
            let original = coarse.value_at(point);
            let transferred = fine.value_at(point);
            assert!((original - transferred).abs() <= 16.0 * f64::EPSILON * (1.0 + original.abs()));
        }
    }
}

#[test]
fn fixed_edges_regions_and_isolated_nodes_have_exact_interpolation_support() {
    let coarse = GridSdf::from_fn(4, &|x, y| x + y - 0.5);
    let fixed: Vec<_> = coarse.nodes().iter().copied().enumerate()
        .filter(|(index, _)| index % 5 == 0 || *index == 12).collect();
    let fine = prolongate_level_set(&coarse, &fixed).unwrap();
    assert_eq!(fine.fixed_nodes.len(), 10); // Nine left-edge nodes plus the isolated node.
    assert!(fine.fixed_nodes.iter().all(|&(i, v)|
        (i % 9 == 0 || i == 40) && v.to_bits() == fine.geometry.nodes()[i].to_bits()));
    let square = [0, 1, 5, 6].map(|i| (i, coarse.nodes()[i]));
    let fine = prolongate_level_set(&coarse, &square).unwrap();
    assert_eq!(fine.fixed_nodes.len(), 9);
    assert!(fine.fixed_nodes.iter().all(|&(i, _)| i % 9 <= 2 && i / 9 <= 2));
}

#[test]
fn transfer_refuses_bad_declarations_without_changing_the_coarse_field() {
    let coarse = GridSdf::from_fn(4, &|x, y| x - y);
    let original: Vec<_> = coarse.nodes().iter().map(|v| v.to_bits()).collect();
    for fixed in [vec![(0, 1.0)], vec![(25, 0.0)], vec![(0, 0.0), (0, 0.0)],
                  vec![(1, 0.25), (0, 0.0)], vec![(0, f64::NAN)]] {
        assert!(prolongate_level_set(&coarse, &fixed).is_err());
    }
    assert_eq!(original, coarse.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>());
    assert!(prolongate_level_set(&GridSdf::from_fn(3, &|_, _| -1.0), &[]).is_err());
    assert!(prolongate_level_set(&GridSdf::from_fn(256, &|_, _| -1.0), &[]).is_err());
    assert!(prolongate_level_set(&GridSdf::from_fn(2, &|_, _| f64::NAN), &[]).is_err());
}

#[test]
fn extreme_finite_midpoints_do_not_overflow() {
    let coarse = GridSdf::from_fn(2, &|x, y| if x + y < 1.0 { f64::MAX } else { -f64::MAX });
    let fine = prolongate_level_set(&coarse, &[]).unwrap();
    assert!(fine.geometry.nodes().iter().all(|v| v.is_finite()));
    let maximum = GridSdf::from_fn(128, &|x, y| x - y);
    assert_eq!(prolongate_level_set(&maximum, &[]).unwrap().geometry.n(), 256);
}

fn fixture(aggregate: RobustAggregate) -> MultiLoadProjectedOptimizer {
    let field = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = field.nodes().iter().copied().enumerate().filter(|(index, _)| {
        index % 9 == 0 || index % 9 == 8 || index / 9 == 0 || index / 9 == 8
    }).collect();
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.75).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 1.0], 0.25).unwrap(),
    ];
    MultiLoadProjectedOptimizer::new(field, &cases,
        OptimizeSettings { level: 3, volfrac: 0.6, iterations: 2, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default() },
        aggregate, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        MultiLoadProjectedSettings { max_candidates: 8, max_solves: 64,
            ..MultiLoadProjectedSettings::default() },
    ).unwrap()
}

#[test]
fn finer_study_rebaselines_at_fixed_area_and_replays_through_the_existing_checkpoint() {
    for aggregate in [RobustAggregate::WeightedSum, RobustAggregate::WorstWeightedCase] {
        let coarse = fixture(aggregate);
        let original = coarse.checkpoint_bytes();
        let (fine, report) = refine_projected_study(&coarse, 3, 96).unwrap();
        assert_eq!(coarse.checkpoint_bytes(), original);
        assert_eq!((report.coarse_level, report.fine_level), (3, 4));
        assert_eq!(report.coarse_endpoint, coarse.current());
        assert_eq!(fine.settings().iterations, 3);
        assert_eq!(fine.next_iteration(), 0);
        assert_eq!(fine.max_solves(), 96);
        assert_eq!(fine.solves_started(), 2);
        assert_eq!(fine.load_cases(), coarse.load_cases());
        assert_eq!(fine.aggregate(), aggregate);
        assert_eq!(fine.settings().youngs.to_bits(), coarse.settings().youngs.to_bits());
        assert_eq!(fine.settings().poisson.to_bits(), coarse.settings().poisson.to_bits());
        assert!((fine.current().volume - 0.6).abs() <= 1e-4);
        assert_eq!(fine.current(), *fine.baseline());
        assert_eq!(report.fine_baseline, fine.current());
        let independent = crate::evaluate_robust_design(fine.geometry(), fine.load_cases(),
            fine.settings(), fine.aggregate()).unwrap();
        assert_eq!(independent.objective.to_bits(), fine.current().objective.to_bits());
        assert_eq!(independent.case_compliances, fine.current().case_compliances);
        assert_eq!(independent.snapshot, fine.current().snapshot);
        let bytes = fine.checkpoint_bytes();
        let restored = MultiLoadProjectedOptimizer::restore_checkpoint(&bytes, &mut 4).unwrap();
        assert_eq!(restored.checkpoint_bytes(), bytes);
        assert!(fine.current_stress().is_none());
    }
}

#[test]
fn stress_policy_survives_refinement_and_spent_source_work_is_not_refunded() {
    let mut coarse = fixture(RobustAggregate::WeightedSum)
        .with_sampled_stress_limit(SampledStressLimit::new(1e6, 0.0).unwrap()).unwrap();
    // Either a genuine step or a bounded terminal attempt must be retained in
    // the source; the test does not mislabel an exhausted search as acceptance.
    let progress = coarse.advance_one().unwrap();
    assert!(matches!(progress, MultiLoadProjectedProgress::Accepted(_)
        | MultiLoadProjectedProgress::NoDescent(_) | MultiLoadProjectedProgress::SolveBudget(_)));
    let original = coarse.checkpoint_bytes();
    let (mut fine, _) = refine_projected_study(&coarse, 2, 2).unwrap();
    assert_eq!(coarse.checkpoint_bytes(), original);
    assert_eq!(fine.stress_limit(), coarse.stress_limit());
    assert!(fine.current_stress().is_some());
    let independent = crate::evaluate_robust_sampled_stress(fine.geometry(), fine.load_cases(),
        fine.settings(), fine.aggregate()).unwrap();
    assert_eq!(fine.current_stress(), Some(&independent));
    assert!(matches!(fine.advance_one().unwrap(), MultiLoadProjectedProgress::SolveBudget(_)));
    assert!(refine_projected_study(&coarse, 0, 2).is_err());
    assert!(refine_projected_study(&coarse, 2, 1).is_err());
    assert_eq!(coarse.checkpoint_bytes(), original);
}
