//! Source material density reaches both modal front doors and their derivatives.

use std::ops::ControlFlow;

use fs_rep_mesh::TetComplex;
use fs_solid::linear3::{TetElasticMaterial, TetMaterialField};
use fs_topopt::modal::{
    EigenfrequencyObjectiveOptions, MatrixFreeEigenError, MatrixFreeEigenObjective,
    MatrixFreeEigenOptions, controlled_matrix_free_eigenfrequency_objective,
};
use fs_topopt::{
    DensityElasticity, DensityFilter, DesignPipeline, SimpParams, SolveBudget, SolveControl,
    eigenfrequency_objective,
};

fn mesh() -> (TetComplex, Vec<[f64; 3]>) {
    let (mesh, mut positions) = fs_feec::kuhn_cube(1);
    // Unequal side lengths separate the bending modes on this six-cell fixture.
    for point in &mut positions {
        point[0] *= 1.6;
        point[1] *= 0.9;
        point[2] *= 0.7;
    }
    (mesh, positions)
}

fn pipeline(mesh: &TetComplex, positions: &[[f64; 3]], beta: f64) -> DesignPipeline {
    DesignPipeline {
        filter: DensityFilter::new(mesh, positions, 0.12),
        params: SimpParams {
            e_min: 1e-6,
            penal: 3.0,
            beta,
            eta: 0.5,
        },
    }
}

fn material(density: f64, youngs: f64, poisson: f64, identity: u8) -> TetElasticMaterial {
    TetElasticMaterial::try_new(density, youngs, poisson, [identity; 32]).unwrap()
}

fn matrix_free(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    rho: &[f64],
    count: usize,
) -> MatrixFreeEigenObjective {
    let options = EigenfrequencyObjectiveOptions {
        solver: MatrixFreeEigenOptions {
            count,
            // The complete twelve-dof space keeps this material regression tiny.
            oversampling: 12,
            max_iterations: 32,
            relative_tolerance: 1e-9,
            linear_tolerance: 1e-12,
            ..Default::default()
        },
        aggregation_beta: 2.0,
        ..Default::default()
    };
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(
        SolveBudget {
            per_solve_iterations: 512,
            total_iterations: 20_000,
        },
        &mut callback,
    );
    controlled_matrix_free_eigenfrequency_objective(
        pipeline,
        elasticity,
        rho,
        options,
        &mut control,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(1e-12);
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, relative tolerance={tolerance:e}"
    );
}

fn vectors_close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let error = actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    let scale = expected
        .iter()
        .map(|v| v * v)
        .sum::<f64>()
        .sqrt()
        .max(1e-12);
    assert!(
        error.is_finite() && error <= tolerance * scale,
        "vector relative error={:e}, tolerance={tolerance:e}",
        error / scale
    );
}

#[test]
fn doubling_material_density_preserves_compliance_and_halves_squared_frequency() {
    let (mesh, positions) = mesh();
    let light = material(2.4, 17.0, 0.28, 1);
    let heavy = material(4.8, 17.0, 0.28, 2);
    let fixed = |p: [f64; 3]| p[0] == 0.0;
    let mut first = DensityElasticity::with_material_field(
        &mesh,
        &positions,
        TetMaterialField::Uniform(&light),
        &fixed,
    )
    .unwrap();
    let mut second = DensityElasticity::with_material_field(
        &mesh,
        &positions,
        TetMaterialField::Uniform(&heavy),
        &fixed,
    )
    .unwrap();
    let pipeline = pipeline(&mesh, &positions, 3.0);
    let rho: Vec<f64> = (0..first.cells())
        .map(|i| 0.42 + 0.035 * i as f64)
        .collect();
    let mut force = vec![0.0; first.n()];
    for (node, point) in positions.iter().enumerate() {
        if point[0] > 1.0 {
            force[3 * node + 2] = -0.25;
        }
    }
    let (compliance, _, gradient) = pipeline.compliance_and_gradient(&mut first, &rho, &force);
    let (heavy_compliance, _, heavy_gradient) =
        pipeline.compliance_and_gradient(&mut second, &rho, &force);
    close(heavy_compliance, compliance, 1e-12);
    vectors_close(&heavy_gradient, &gradient, 1e-12);
    let first_mass: f64 = first.reference_cell_masses_kg().unwrap().iter().sum();
    let second_mass: f64 = second.reference_cell_masses_kg().unwrap().iter().sum();
    close(second_mass, 2.0 * first_mass, 1e-14);

    let (lambda, gradient) = eigenfrequency_objective(&pipeline, &mut first, &rho, 1, 2.0);
    let (heavy_lambda, heavy_gradient) =
        eigenfrequency_objective(&pipeline, &mut second, &rho, 1, 2.0);
    close(heavy_lambda, 0.5 * lambda, 1e-9);
    let half_gradient: Vec<f64> = gradient.iter().map(|g| 0.5 * g).collect();
    vectors_close(&heavy_gradient, &half_gradient, 1e-8);
    let first_modes = matrix_free(&pipeline, &mut first, &rho, 1);
    let second_modes = matrix_free(&pipeline, &mut second, &rho, 1);
    close(first_modes.modal.values[0], lambda, 1e-7);
    close(second_modes.modal.values[0], 0.5 * lambda, 1e-7);

    // Low-level mass weights remain absolute kg/m3, so the source density must
    // not also be folded into their geometric unit-density element blocks.
    let explicit_mass = vec![7.0; first.cells()];
    let (_, mass_first, _) = first.assemble_dense(&explicit_mass);
    let (_, mass_second, _) = second.assemble_dense(&explicit_mass);
    assert_eq!(mass_first, mass_second);
}

#[test]
fn heterogeneous_material_modes_and_gradients_match_dense_and_finite_differences() {
    let (mesh, positions) = mesh();
    let materials: Vec<_> = (0..mesh.tets.len())
        .map(|cell| {
            if cell % 2 == 0 {
                material(2.4, 17.0, 0.28, 3)
            } else {
                material(7.1, 31.0, 0.22, 4)
            }
        })
        .collect();
    let mut elasticity = DensityElasticity::with_material_field(
        &mesh,
        &positions,
        TetMaterialField::PerElement(&materials),
        &|p| p[0] == 0.0,
    )
    .unwrap();
    drop(materials);
    let Some(TetMaterialField::PerElement(retained)) = elasticity.material_field() else {
        panic!("per-cell materials must remain attached to the operator");
    };
    for (cell, retained) in retained.iter().enumerate() {
        assert_eq!(
            retained.material_state_identity,
            [if cell % 2 == 0 { 3 } else { 4 }; 32]
        );
        close(
            elasticity.material_density_kg_m3(cell).unwrap(),
            if cell % 2 == 0 { 2.4 } else { 7.1 },
            1e-14,
        );
    }

    let direction = [0.8, -0.3, 0.5, -0.6, 0.2, 0.7];
    assert_eq!(elasticity.cells(), direction.len());
    // Exercise both ordinary and rho^6 void mass interpolation after a
    // nontrivial Helmholtz filter and Heaviside projection.
    for (base, variation, beta) in [(0.42, 0.035, 3.0), (0.055, 0.003, 2.0)] {
        let pipeline = pipeline(&mesh, &positions, beta);
        let rho: Vec<f64> = (0..elasticity.cells())
            .map(|i| base + variation * i as f64)
            .collect();
        let (dense, dense_gradient) =
            eigenfrequency_objective(&pipeline, &mut elasticity, &rho, 3, 2.0);
        let native = matrix_free(&pipeline, &mut elasticity, &rho, 3);
        close(native.aggregate, dense, 1e-7);
        vectors_close(&native.gradient, &dense_gradient, 2e-6);

        let step = 1e-5;
        let plus: Vec<f64> = rho
            .iter()
            .zip(&direction)
            .map(|(r, d)| r + step * d)
            .collect();
        let minus: Vec<f64> = rho
            .iter()
            .zip(&direction)
            .map(|(r, d)| r - step * d)
            .collect();
        let high = eigenfrequency_objective(&pipeline, &mut elasticity, &plus, 3, 2.0).0;
        let low = eigenfrequency_objective(&pipeline, &mut elasticity, &minus, 3, 2.0).0;
        let dense_fd = (high - low) / (2.0 * step);
        let high = matrix_free(&pipeline, &mut elasticity, &plus, 3).aggregate;
        let low = matrix_free(&pipeline, &mut elasticity, &minus, 3).aggregate;
        let native_fd = (high - low) / (2.0 * step);
        let dense_directional: f64 = dense_gradient
            .iter()
            .zip(&direction)
            .map(|(g, d)| g * d)
            .sum();
        let native_directional: f64 = native
            .gradient
            .iter()
            .zip(&direction)
            .map(|(g, d)| g * d)
            .sum();
        close(dense_directional, dense_fd, 2e-4);
        close(native_directional, native_fd, 2e-4);
    }
}

#[test]
fn bound_material_density_cannot_be_rescaled_by_modal_options() {
    let (mesh, positions) = mesh();
    let source = material(2.4, 17.0, 0.28, 5);
    let mut elasticity = DensityElasticity::with_material_field(
        &mesh,
        &positions,
        TetMaterialField::Uniform(&source),
        &|p| p[0] == 0.0,
    )
    .unwrap();
    let previous = elasticity.moduli.clone();
    let pipeline = pipeline(&mesh, &positions, 3.0);
    let rho = vec![0.5; elasticity.cells()];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let result = controlled_matrix_free_eigenfrequency_objective(
        &pipeline,
        &mut elasticity,
        &rho,
        EigenfrequencyObjectiveOptions {
            reference_density: 2.4,
            ..Default::default()
        },
        &mut control,
    );
    assert!(matches!(result, Err(MatrixFreeEigenError::InvalidInput(_))));
    assert_eq!(elasticity.moduli, previous);
    assert_eq!(control.work().linear_solves, 0);
}
