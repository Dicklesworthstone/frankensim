use super::*;
use crate::refinement::prolongate_level_set;
use crate::volume::{VolumeProjectionSettings, project_material_volume};

fn rectangle(phase: DesignPhase, lo: [f64; 2], hi: [f64; 2], margin: f64) -> DesignRegion {
    DesignRegion::new(phase, lo, hi, margin).unwrap()
}
fn bits(field: &GridSdf) -> Vec<u64> { field.nodes().iter().map(|v| v.to_bits()).collect() }

#[test]
fn whole_cell_coverage_retains_a_subcell_clearance_and_original_input() {
    let input = GridSdf::from_fn(8, &|_, _| -0.1);
    let before = bits(&input);
    let region = rectangle(DesignPhase::Void, [0.30, 0.30], [0.31, 0.31], 0.02);
    let result = prepare_design_regions(&input, &[], &[region]).unwrap();
    assert_eq!(bits(&input), before);
    assert_eq!(result.covered_cells, [[2, 2, 3, 3]]);
    assert_eq!(result.fixed_nodes.len(), 4);
    assert_eq!(result.changed_nodes, 4);
    assert_eq!(result.void_nodes, 4);
    for i in 0..=20 {
        for j in 0..=20 {
            let x = 0.25 + f64::from(i) / 160.0;
            let y = 0.25 + f64::from(j) / 160.0;
            assert!(result.geometry.value_at([x, y]) >= 0.02 - 1e-16);
        }
    }
    for (index, value) in input.nodes().iter().enumerate() {
        if !result.fixed_nodes.iter().any(|&(node, _)| node == index) {
            assert_eq!(result.geometry.nodes()[index].to_bits(), value.to_bits());
        }
    }
}

#[test]
fn strongest_same_phase_margin_is_order_independent_and_fixed_bits_survive() {
    let input = GridSdf::from_fn(8, &|_, _| -0.1);
    let a = rectangle(DesignPhase::Material, [0.125, 0.25], [0.5, 0.5], 0.2);
    let b = rectangle(DesignPhase::Material, [0.25, 0.25], [0.625, 0.5], 0.3);
    let fixed = [(0, -0.1), (80, -0.1)];
    let first = prepare_design_regions(&input, &fixed, &[a, b]).unwrap();
    let second = prepare_design_regions(&input, &fixed, &[b, a]).unwrap();
    assert_eq!(bits(&first.geometry), bits(&second.geometry));
    assert_eq!(first.fixed_nodes, second.fixed_nodes);
    assert_eq!(first.geometry.node(3, 3), -0.3);
    assert_eq!(first.fixed_nodes.first(), Some(&(0, -0.1)));
    assert_eq!(first.fixed_nodes.last(), Some(&(80, -0.1)));
}

#[test]
fn opposite_phases_and_incompatible_existing_fixed_values_refuse_transactionally() {
    let input = GridSdf::from_fn(8, &|_, _| -0.1);
    let before = bits(&input);
    let a = rectangle(DesignPhase::Material, [0.125, 0.25], [0.25, 0.5], 0.02);
    let b = rectangle(DesignPhase::Void, [0.25, 0.25], [0.375, 0.5], 0.02);
    // No positive-area rectangle overlap, but a shared Q1 node cannot have both signs.
    assert!(prepare_design_regions(&input, &[], &[a, b]).is_err());
    assert!(prepare_design_regions(&input, &[(2 + 2 * 9, -0.1)], &[b]).is_err());
    let stronger = rectangle(DesignPhase::Material, [0.0, 0.0], [0.125, 0.125], 0.2);
    assert!(prepare_design_regions(&input, &[(0, -0.1)], &[stronger]).is_err());
    assert_eq!(bits(&input), before);
}

#[test]
fn malformed_requests_and_empty_identity_have_no_authoring_side_effect() {
    for (lo, hi, margin) in [
        ([0.0, 0.0], [0.0, 0.2], 0.1), ([f64::NAN, 0.0], [0.2, 0.2], 0.1),
        ([-0.1, 0.0], [0.2, 0.2], 0.1), ([0.0, 0.0], [1.1, 0.2], 0.1),
        ([0.0, 0.0], [0.2, 0.2], 0.0), ([0.0, 0.0], [0.2, 0.2], f64::INFINITY),
    ] { assert!(DesignRegion::new(DesignPhase::Void, lo, hi, margin).is_err()); }
    let input = GridSdf::from_fn(8, &|x, y| x - y);
    let empty = prepare_design_regions(&input, &[(0, 0.0)], &[]).unwrap();
    assert_eq!(bits(&empty.geometry), bits(&input));
    assert_eq!(empty.fixed_nodes, [(0, 0.0)]);
    assert_eq!(empty.changed_nodes, 0);
    assert!(prepare_design_regions(&input, &[(0, -0.0)], &[]).is_err());
    assert!(prepare_design_regions(&input, &[(0, 0.0), (0, 0.0)], &[]).is_err());
    assert!(prepare_design_regions(&input, &[(81, 0.0)], &[]).is_err());
    let all = rectangle(DesignPhase::Material, [0.0, 0.0], [1.0, 1.0], 0.01);
    assert!(prepare_design_regions(&input, &[], &[all; 64]).is_ok());
    assert!(prepare_design_regions(&input, &[], &[all; 65]).is_err());
    assert!(prepare_design_regions(&GridSdf::from_fn(3, &|_, _| 0.1), &[], &[]).is_err());
}

#[test]
fn every_rasterization_and_publication_boundary_cancels_without_mutation() {
    let input = GridSdf::from_fn(8, &|_, _| -0.1);
    let before = bits(&input);
    let region = rectangle(DesignPhase::Void, [0.25, 0.25], [0.5, 0.5], 0.02);
    let mut stages = Vec::new();
    prepare_design_regions_controlled(&input, &[], &[region], |stage| {
        stages.push(stage); ControlFlow::<()>::Continue(())
    }).unwrap();
    for stop in stages {
        let result = prepare_design_regions_controlled(&input, &[], &[region], |stage| {
            if stage == stop { ControlFlow::Break(stop) } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(matches!(result, ControlFlow::Break(reason) if reason == stop));
        assert_eq!(bits(&input), before);
    }
}

#[test]
fn projection_and_refinement_preserve_whole_solid_and_void_cells() {
    let input = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.42);
    let regions = [
        rectangle(DesignPhase::Material, [0.125, 0.375], [0.25, 0.625], 0.04),
        rectangle(DesignPhase::Void, [0.5, 0.5], [0.625, 0.625], 0.04),
    ];
    let mut prepared = prepare_design_regions(&input, &[], &regions).unwrap();
    project_material_volume(&mut prepared.geometry, 3, &prepared.fixed_nodes,
        VolumeProjectionSettings { target: 0.7, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 }).unwrap();
    for &(index, value) in &prepared.fixed_nodes {
        assert_eq!(prepared.geometry.nodes()[index].to_bits(), value.to_bits());
    }
    let fine = prolongate_level_set(&prepared.geometry, &prepared.fixed_nodes).unwrap();
    for (phase, lo, hi) in [
        (DesignPhase::Material, [0.125, 0.375], [0.25, 0.625]),
        (DesignPhase::Void, [0.5, 0.5], [0.625, 0.625]),
    ] {
        for i in 0..=10 {
            for j in 0..=10 {
                let p = [lo[0] + (hi[0] - lo[0]) * f64::from(i) / 10.0,
                    lo[1] + (hi[1] - lo[1]) * f64::from(j) / 10.0];
                let v = fine.geometry.value_at(p);
                assert!(match phase { DesignPhase::Material => v < 0.0, DesignPhase::Void => v > 0.0 });
            }
        }
    }
}

#[test]
fn actual_projected_updates_and_recovered_studies_keep_the_authored_clearance() {
    use crate::robust_descent::{MultiLoadProjectedOptimizer, MultiLoadProjectedProgress, MultiLoadProjectedSettings};
    use crate::{OptimizeSettings, RobustAggregate, RobustLoadCase};
    use fs_cutfem::DesignBoxEdge;
    let input = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.42);
    let fixed: Vec<_> = input.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8 || i / 9 == 0 || i / 9 == 8).collect();
    let prepared = prepare_design_regions(&input, &fixed, &[
        rectangle(DesignPhase::Void, [0.5, 0.5], [0.625, 0.625], 0.02),
    ]).unwrap();
    let expected = prepared.fixed_nodes.clone();
    let cases = [RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap()];
    let mut optimizer = MultiLoadProjectedOptimizer::new(prepared.geometry, &cases,
        OptimizeSettings { level: 3, iterations: 1, volfrac: 0.7, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default() }, RobustAggregate::WeightedSum,
        prepared.fixed_nodes, VolumeProjectionSettings {
            target: 0.7, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64,
        }, MultiLoadProjectedSettings { max_candidates: 16, max_solves: 17, ..MultiLoadProjectedSettings::default() }).unwrap();
    assert!(matches!(optimizer.advance_one().unwrap(), MultiLoadProjectedProgress::Accepted(_)),
        "regression requires a real optimized endpoint, not only baseline preservation");
    let mut recovery = 2;
    let restored = MultiLoadProjectedOptimizer::restore_checkpoint(&optimizer.checkpoint_bytes(), &mut recovery).unwrap();
    assert_eq!(restored.fixed_nodes(), expected.as_slice());
    assert_eq!(bits(restored.geometry()), bits(optimizer.geometry()));
    for &(index, value) in &expected {
        assert_eq!(restored.geometry().nodes()[index].to_bits(), value.to_bits());
    }
}
