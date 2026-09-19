//! G0/G3/G5 checks of the real filtered, independently loaded elasticity path.
use std::ops::ControlFlow;

use fs_topopt::multi_load::{
    MultiLoadOcOptions, MultiLoadOcTermination, multi_load_optimality_criteria,
};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::{DensityElasticity, DensityFilter, DesignPipeline, SimpParams};

fn fixture() -> (DesignPipeline, DensityElasticity, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (complex, positions) = fs_feec::kuhn_cube(2);
    let elasticity = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p| p[0] < 1e-12);
    let mut force = vec![0.0; elasticity.n()];
    for (v, p) in positions.iter().enumerate() {
        if p[0] > 1.0 - 1e-12 { force[3 * v + 2] = -1.0; }
    }
    let rho = vec![0.5; elasticity.cells()];
    let volumes = fs_feec::element_geometry(&complex, &positions)
        .vol_signed.iter().map(|v| v.abs()).collect();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params: SimpParams::default(),
    };
    (pipeline, elasticity, rho, force, volumes)
}

fn options() -> MultiLoadOcOptions {
    MultiLoadOcOptions { max_iterations: 5, change_tolerance: 0.0, ..Default::default() }
}

#[test]
fn g3_independent_opposite_loads_do_not_cancel_during_optimization() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let opposite: Vec<f64> = force.iter().map(|f| -f).collect();
    let run = multi_load_optimality_criteria(&pipeline, &mut elasticity, &[
        LoadCase { force: &force, weight: 0.25 },
        LoadCase { force: &opposite, weight: 0.75 },
    ], &rho, &volumes, options(), || ControlFlow::Continue(()));
    assert!(run.history.len() > 1);
    assert!(run.history[0].compliance > 0.0);
    assert!(run.history.last().unwrap().compliance < run.history[0].compliance);
    for row in &run.history {
        assert!(row.volume_fraction <= 0.5 + 1e-8);
        assert_eq!(row.case_compliances.len(), 2);
        assert!((row.case_compliances[0] - row.case_compliances[1]).abs()
            <= 1e-8 * row.case_compliances[0].max(1.0));
    }
    for rows in run.history.windows(2) {
        assert!(rows[1].compliance <= rows[0].compliance);
        assert!(rows[1].max_change <= options().move_limit + 1e-14);
    }
    assert!(run.rho.iter().all(|r| (1e-3..=1.0).contains(r)));
}

#[test]
fn g0_final_history_describes_returned_design_and_operator() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let loads = [LoadCase { force: &force, weight: 1.0 }];
    let run = multi_load_optimality_criteria(&pipeline, &mut elasticity, &loads,
        &rho, &volumes, options(), || ControlFlow::Continue(()));
    let (_, projected, expected_moduli) = pipeline.forward(&run.rho);
    assert_eq!(run.projected_rho, projected);
    assert_eq!(elasticity.moduli, expected_moduli);
    let actual = pipeline.multi_load_compliance_and_gradient(&mut elasticity, &run.rho, &loads);
    assert_eq!(run.displacements, actual.displacements);
    let (volume, _) = pipeline.volume_and_gradient(&run.rho, &volumes);
    let last = run.history.last().unwrap();
    assert_eq!(last.compliance.to_bits(), actual.compliance.to_bits());
    assert!((last.volume_fraction - volume).abs() < 1e-12);
    assert_eq!(last.iteration + 1, run.history.len());
}

#[test]
fn g5_repeated_studies_replay_bitwise() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let loads = [LoadCase { force: &force, weight: 1.0 }];
    let a = multi_load_optimality_criteria(&pipeline, &mut elasticity, &loads,
        &rho, &volumes, options(), || ControlFlow::Continue(()));
    let b = multi_load_optimality_criteria(&pipeline, &mut elasticity, &loads,
        &rho, &volumes, options(), || ControlFlow::Continue(()));
    assert_eq!(a.termination, b.termination);
    assert_eq!(a.rho, b.rho);
    assert_eq!(a.projected_rho, b.projected_rho);
    assert_eq!(a.displacements, b.displacements);
    assert_eq!(a.history.len(), b.history.len());
    for (a, b) in a.history.iter().zip(&b.history) {
        assert_eq!(a.compliance.to_bits(), b.compliance.to_bits());
        assert_eq!(a.volume_fraction.to_bits(), b.volume_fraction.to_bits());
        assert_eq!(a.max_change.to_bits(), b.max_change.to_bits());
        assert_eq!(a.case_compliances, b.case_compliances);
    }
}

#[test]
fn g4_cancel_before_analysis_does_not_fabricate_a_solved_design() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let original_moduli = elasticity.moduli.clone();
    let run = multi_load_optimality_criteria(&pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes, options(),
        || ControlFlow::Break(()));
    assert_eq!(run.termination, MultiLoadOcTermination::Cancelled);
    assert!(run.history.is_empty());
    assert!(run.projected_rho.is_empty());
    assert!(run.displacements.is_empty());
    assert_eq!(run.rho, rho);
    assert_eq!(elasticity.moduli, original_moduli);
}

#[test]
fn g4_cancel_during_multiplier_search_retains_only_solved_state() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let mut calls = 0;
    let run = multi_load_optimality_criteria(&pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes, options(), || {
            calls += 1;
            if calls == 6 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
        });
    assert_eq!(run.termination, MultiLoadOcTermination::Cancelled);
    assert_eq!(run.history.len(), 1);
    assert_eq!(run.rho, rho);
    assert_eq!(elasticity.moduli, pipeline.forward(&rho).2);
    let solved = pipeline.multi_load_compliance_and_gradient(&mut elasticity, &rho,
        &[LoadCase { force: &force, weight: 1.0 }]);
    assert_eq!(run.displacements, solved.displacements);
    assert_eq!(run.projected_rho, pipeline.forward(&rho).1);
}

#[test]
fn g0_zero_update_budget_still_reports_real_baseline() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let run = multi_load_optimality_criteria(&pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
        MultiLoadOcOptions { max_iterations: 0, ..options() }, || ControlFlow::Continue(()));
    assert_eq!(run.termination, MultiLoadOcTermination::IterationBudget);
    assert_eq!(run.history.len(), 1);
    assert_eq!(run.rho, rho);
    assert!(run.history[0].compliance > 0.0);
}

#[test]
#[should_panic(expected = "starting projected design exceeds")]
fn g0_full_material_is_not_admitted_as_a_half_volume_baseline() {
    let (pipeline, mut elasticity, mut rho, force, volumes) = fixture();
    rho.fill(1.0);
    multi_load_optimality_criteria(&pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
        options(), || ControlFlow::Continue(()));
}
