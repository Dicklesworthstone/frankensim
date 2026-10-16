//! Physical pencil and design-gradient regressions against the dense reference.
use std::ops::ControlFlow;
use fs_solver::LinearOp;
use fs_topopt::{DensityElasticity, DensityFilter, DesignPipeline, SimpParams,
    SolveBudget, SolveControl, SolveProgress, EvaluationStop};
use fs_topopt::modal::{MatrixFreeEigenOptions, EigenfrequencyObjectiveOptions, MatrixFreeEigenError,
    controlled_matrix_free_eigenpairs, controlled_matrix_free_eigenfrequency_objective};

fn settings() -> MatrixFreeEigenOptions {
    MatrixFreeEigenOptions { relative_tolerance: 1e-9, linear_tolerance: 1e-12, ..Default::default() }
}
fn norm(v: &[f64]) -> f64 { v.iter().map(|x| x * x).sum::<f64>().sqrt() }

#[test]
fn physical_consistent_mass_modes_match_dense_and_original_residuals() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let mut elasticity = DensityElasticity::new(&mesh, &positions, 17.0, 0.3, &|p| p[0] == 0.0);
    let rho: Vec<f64> = (0..elasticity.cells()).map(|i| 0.4 + 0.03 * (i % 7) as f64).collect();
    elasticity.moduli = rho.iter().map(|r| 1e-6 + (1.0 - 1e-6) * r * r * r).collect();
    let mass: Vec<f64> = rho.iter().map(|r| 2.7 * r).collect();
    let (reference, _, _) = fs_topopt::eigenfreq::lowest_eigenpairs(&elasticity, &mass, 3);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let result = controlled_matrix_free_eigenpairs(&elasticity, &mass, settings(), &mut control).unwrap();
    assert_eq!(result.values.len(), 3);
    for i in 0..3 {
        assert!((result.values[i] - reference[i]).abs() / reference[i] < 1e-7);
        let phi = &result.modes[i];
        let mut kphi = vec![0.0; elasticity.n()];
        let mut mphi = vec![0.0; elasticity.n()];
        elasticity.apply(phi, &mut kphi);
        elasticity.apply_mass(&mass, phi, &mut mphi);
        let residual: Vec<f64> = kphi.iter().zip(&mphi).map(|(k, m)| k - result.values[i] * m).collect();
        let relative = norm(&residual) / (norm(&kphi) + result.values[i] * norm(&mphi));
        assert!(relative <= settings().relative_tolerance * 1.01);
        assert!((relative - result.relative_residuals[i]).abs() < 1e-12);
        assert!((phi.iter().zip(&mphi).map(|(p, m)| p * m).sum::<f64>() - 1.0).abs() < 1e-9);
        for (d, free) in elasticity.free().iter().enumerate() {
            if !free { assert_eq!(phi[d], 0.0); }
        }
    }
    assert!(result.work.linear_solves > 0);
}

#[test]
fn frequency_gradient_crosses_filter_projection_simp_and_physical_mass() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let mut elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.3, &|p| p[0] == 0.0);
    let direction: Vec<f64> = (0..elasticity.cells()).map(|i| 0.1 + (i % 5) as f64 / 10.0).collect();
    let options = EigenfrequencyObjectiveOptions { solver: settings(), aggregation_beta: 2.0, reference_density: 2.7 };
    for (base, variation, beta) in [(0.4, 0.02, 2.0), (0.4, 0.02, 7.0), (0.06, 0.002, 2.0)] {
        let rho: Vec<f64> = (0..elasticity.cells())
            .map(|i| base + variation * (i % 9) as f64).collect();
        let pipeline = DesignPipeline { filter: DensityFilter::new(&mesh, &positions, 0.1),
            params: SimpParams { beta, ..Default::default() } };
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let center = controlled_matrix_free_eigenfrequency_objective(&pipeline, &mut elasticity, &rho,
            options, &mut control).unwrap();
        let step = 1e-5;
        let plus: Vec<f64> = rho.iter().zip(&direction).map(|(r, d)| r + step * d).collect();
        let minus: Vec<f64> = rho.iter().zip(&direction).map(|(r, d)| r - step * d).collect();
        let high = controlled_matrix_free_eigenfrequency_objective(&pipeline, &mut elasticity, &plus,
            options, &mut control).unwrap().aggregate;
        let low = controlled_matrix_free_eigenfrequency_objective(&pipeline, &mut elasticity, &minus,
            options, &mut control).unwrap().aggregate;
        let actual: f64 = center.gradient.iter().zip(&direction).map(|(g, d)| g * d).sum();
        let fd = (high - low) / (2.0 * step);
        assert!((actual - fd).abs() / actual.abs().max(fd.abs()).max(1e-8) < 2e-4,
            "beta={beta}: analytic {actual:e}, FD {fd:e}");
    }
}

#[test]
fn cancelled_modal_solve_restores_moduli_without_refunding_work() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let mut elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.3, &|p| p[0] == 0.0);
    let previous = elasticity.moduli.clone();
    let rho = vec![0.55; elasticity.cells()];
    let pipeline = DesignPipeline { filter: DensityFilter::new(&mesh, &positions, 0.1),
        params: SimpParams::default() };
    let mut callback = |p: SolveProgress| {
        if p.stage == "eigen-elasticity" && p.solve_iterations > 0 { ControlFlow::Break(()) }
        else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let result = controlled_matrix_free_eigenfrequency_objective(&pipeline, &mut elasticity, &rho,
        EigenfrequencyObjectiveOptions::default(), &mut control);
    assert!(matches!(result, Err(MatrixFreeEigenError::Operator(EvaluationStop::Cancelled))));
    assert_eq!(elasticity.moduli, previous);
    assert!(control.work().linear_iterations > 0);
}

#[test]
fn empty_mass_and_linear_budget_exhaustion_are_not_valid_mode_reports() {
    let (mesh, positions) = fs_feec::kuhn_cube(2);
    let elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.3, &|p| p[0] == 0.0);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { per_solve_iterations: 0, total_iterations: 100 }, &mut callback);
    let empty = vec![0.0; elasticity.cells()];
    assert!(matches!(controlled_matrix_free_eigenpairs(&elasticity, &empty, settings(), &mut control),
        Err(MatrixFreeEigenError::InvalidInput(_))));
    let mass = vec![1.0; elasticity.cells()];
    assert!(matches!(controlled_matrix_free_eigenpairs(&elasticity, &mass, settings(), &mut control),
        Err(MatrixFreeEigenError::Operator(EvaluationStop::LinearBudget { iterations: 0, .. }))));
    assert_eq!(control.work().linear_iterations, 0);
}
