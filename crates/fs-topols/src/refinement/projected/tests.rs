use super::*;
use crate::{Cantilever, GridSdf, evaluate_compliance_design};
use crate::evaluated::DesignEvaluationStage;
use crate::projected::{ProjectedProgress, ProjectedSettings};
use crate::volume::VolumeProjectionSettings;

fn coarse(load: f64) -> ProjectedOptimizer {
    let phi = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
    let fixed = phi.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    ProjectedOptimizer::new(phi, Cantilever { load, band: 0.125 },
        OptimizeSettings { level: 3, iterations: 2, volfrac: 0.6, move_cells: 0.1,
            nucleation_period: 0, ..OptimizeSettings::default() }, fixed,
        VolumeProjectionSettings { target: 0.6, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        ProjectedSettings { max_candidates: 8, poll_iters: 1, ..ProjectedSettings::default() }).unwrap()
}
fn refine(coarse: &ProjectedOptimizer) -> (ProjectedOptimizer, VolumeRefinementReport) {
    let ControlFlow::Continue(result) = refine_controlled(coarse, 3,
        |_| ControlFlow::<()>::Continue(())).unwrap() else { panic!("unexpected stop") };
    result
}
fn same_source(a: &ProjectedOptimizer, b: &ProjectedOptimizer) {
    assert_eq!(a.checkpoint().geometry().nodes(), b.checkpoint().geometry().nodes());
    assert_eq!(a.checkpoint().next_iteration(), b.checkpoint().next_iteration());
    assert_eq!(a.checkpoint().ell().to_bits(), b.checkpoint().ell().to_bits());
    assert_eq!(a.current().compliance.to_bits(), b.current().compliance.to_bits());
}

#[test]
fn projected_refinement_uses_the_accepted_endpoint_and_solves_a_new_feasible_baseline() {
    let mut source = coarse(1.0);
    assert!(matches!(source.advance_one().unwrap(), ProjectedProgress::Accepted(_)));
    let before = source.clone();
    let transfer = prolongate_level_set(source.checkpoint().geometry(), source.fixed_nodes()).unwrap();
    let (fine, report) = refine(&source);
    same_source(&source, &before);
    assert_eq!((report.coarse_level, report.fine_level, report.source_updates), (3, 4, 1));
    assert_eq!(fine.checkpoint().next_iteration(), 0);
    assert_eq!(fine.checkpoint().settings().iterations, 3);
    assert_eq!(fine.checkpoint().ell().to_bits(), source.checkpoint().settings().ell0.to_bits());
    assert_eq!(fine.fixed_nodes(), transfer.fixed_nodes);
    for &(node, value) in fine.fixed_nodes() {
        assert_eq!(fine.checkpoint().geometry().nodes()[node].to_bits(), value.to_bits());
    }
    assert!((fine.current().volume - 0.6).abs() <= 1e-4);
    let oracle = evaluate_compliance_design(fine.checkpoint().geometry(),
        fine.checkpoint().fixture(), fine.checkpoint().settings()).unwrap();
    assert_eq!(oracle.snapshot, report.fine_baseline.snapshot);
    assert_eq!(oracle.compliance.to_bits(), report.fine_baseline.compliance.to_bits());
    assert_eq!(fine.baseline().compliance.to_bits(), fine.current().compliance.to_bits());
    assert_eq!(report.transferred_area.to_bits(),
        material_volume(&Quadtree::uniform(4), &transfer.geometry).to_bits());
}

#[test]
fn projected_refinement_cancellation_inside_fine_cg_and_before_return_preserves_source() {
    let source = coarse(1.0);
    let before = source.clone();
    for inside in [true, false] {
        let mut reached = false;
        let result = refine_controlled(&source, 3, |stage| {
            let stop = if inside {
                matches!(stage, VolumeRefinementStage::Baseline(ProjectedSetupStage::Evaluation(
                    DesignEvaluationStage::Solve(n))) if n > 0)
            } else { stage == VolumeRefinementStage::Publish };
            if stop { reached = true; ControlFlow::Break("stop") } else { ControlFlow::Continue(()) }
        }).unwrap();
        assert!(reached);
        assert!(matches!(result, ControlFlow::Break("stop")));
        same_source(&source, &before);
    }
    let (retry, _) = refine(&source);
    let (reference, _) = refine(&before);
    same_source(&retry, &reference);
}

#[test]
fn projected_refinement_keeps_canonical_checkpoint_resume_on_the_fine_grid() {
    let (mut fine, _) = refine(&coarse(1.0));
    let mut resumed = ProjectedOptimizer::from_checkpoint(fine.checkpoint().clone(),
        fine.fixed_nodes().to_vec(), fine.projection_settings(), fine.controls()).unwrap();
    let a = fine.advance_one().unwrap();
    let b = resumed.advance_one().unwrap();
    assert_eq!(std::mem::discriminant(&a), std::mem::discriminant(&b));
    same_source(&fine, &resumed);
}

#[test]
fn projected_refinement_refuses_invalid_new_work_before_callbacks() {
    let source = coarse(1.0);
    for updates in [0, 10_001, usize::MAX] {
        let result = refine_controlled(&source, updates, |_| -> ControlFlow<()> {
            panic!("invalid work must be refused before any numerical stage")
        });
        assert!(result.is_err());
    }
}

#[test]
fn projected_refinement_preserves_load_scaling_in_actual_fine_mechanics() {
    let (one, _) = refine(&coarse(1.0));
    let (two, _) = refine(&coarse(2.0));
    assert_eq!(one.checkpoint().geometry().nodes(), two.checkpoint().geometry().nodes());
    assert!((two.current().compliance / one.current().compliance - 4.0).abs() < 1e-8);
}
