//! Real P1/filter checks of interruption and atomic evaluation, not mock solves.
use std::ops::ControlFlow;
use fs_topopt::{DensityElasticity, DensityFilter, DesignPipeline, EvaluationStop,
    SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::pipeline::LoadCase;

fn fixture() -> (DesignPipeline, DensityElasticity, Vec<f64>, Vec<f64>) {
    let (complex, positions) = fs_feec::kuhn_cube(2);
    let elasticity = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p| p[0] < 1e-12);
    let mut force = vec![0.0; elasticity.n()];
    for (i, p) in positions.iter().enumerate() {
        if p[0] > 1.0 - 1e-12 { force[3 * i + 2] = -1.0; }
    }
    let rho = (0..elasticity.cells()).map(|i| 0.3 + 0.02 * (i % 17) as f64).collect();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15), params: SimpParams::default(),
    };
    (pipeline, elasticity, rho, force)
}

#[test]
fn g3_controlled_multi_load_matches_independent_single_load_solves() {
    let (pipeline, mut elasticity, rho, force) = fixture();
    let opposite: Vec<f64> = force.iter().map(|f| -2.0 * f).collect();
    let (c, u, g) = pipeline.compliance_and_gradient(&mut elasticity, &rho, &force);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let result = pipeline.try_multi_load_compliance_and_gradient(&mut elasticity, &rho, &[
        LoadCase { force: &force, weight: 0.25 }, LoadCase { force: &opposite, weight: 0.75 },
    ], &mut control).unwrap();
    assert!((result.compliance - 3.25 * c).abs() < 1e-9 * c);
    assert_eq!(result.displacements[0], u);
    for (actual, expected) in result.gradient.iter().zip(&g) {
        assert!((actual - 3.25 * expected).abs() < 1e-8 * expected.abs().max(1.0));
    }
    assert_eq!(control.work().linear_solves, 4); // filter, two loads, transpose
    assert!(control.work().linear_iterations > 0);
}

#[test]
fn g4_stops_inside_elasticity_and_adjoint_restore_the_previous_operator() {
    for stage in ["elasticity", "filter-transpose"] {
        let (pipeline, mut elasticity, rho, force) = fixture();
        let previous = elasticity.moduli.clone();
        let mut callback = |progress: SolveProgress| {
            if progress.stage == stage && progress.solve_iterations > 0 { ControlFlow::Break(()) }
            else { ControlFlow::Continue(()) }
        };
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let result = pipeline.try_multi_load_compliance_and_gradient(&mut elasticity, &rho,
            &[LoadCase { force: &force, weight: 1.0 }], &mut control);
        assert!(matches!(result, Err(EvaluationStop::Cancelled)), "stage={stage}");
        assert_eq!(elasticity.moduli, previous);
        assert!(control.work().linear_iterations > 0);
    }
}

#[test]
fn g4_total_iteration_budget_is_a_recoverable_stop_not_a_partial_gradient() {
    let (pipeline, mut elasticity, rho, force) = fixture();
    let previous = elasticity.moduli.clone();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget {
        per_solve_iterations: 100, total_iterations: 1,
    }, &mut callback);
    let result = pipeline.try_multi_load_compliance_and_gradient(&mut elasticity, &rho,
        &[LoadCase { force: &force, weight: 1.0 }], &mut control);
    assert!(matches!(result, Err(EvaluationStop::TotalBudget { .. })));
    assert_eq!(elasticity.moduli, previous);
    assert_eq!(control.work().linear_iterations, 1);
}
