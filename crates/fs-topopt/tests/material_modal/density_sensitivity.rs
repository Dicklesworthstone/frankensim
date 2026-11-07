//! Physical density is an independent material coordinate, not SIMP occupancy.

use super::*;
use fs_topopt::eigenfreq::{
    eigenvalue_reference_density_gradient, lowest_eigenpairs, mass_interp, scatter_full, smooth_min,
};
use fs_topopt::elasticity::ReferenceDensitySensitivityError;

fn dense_pullback(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    occupancy: &[f64],
    count: usize,
) -> (f64, Vec<f64>, f64) {
    let (_, projected, moduli) = pipeline.forward(occupancy);
    elasticity.moduli = moduli;
    let mass: Vec<_> = projected
        .iter()
        .enumerate()
        .map(|(cell, &r)| elasticity.reference_density_kg_m3(cell).unwrap() * mass_interp(r))
        .collect();
    let (values, modes, free) = lowest_eigenpairs(elasticity, &mass, count);
    let (aggregate, weights) = smooth_min(&values, 2.0);
    let mut gradient = vec![0.0; elasticity.n() / 3];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(
        SolveBudget {
            per_solve_iterations: 0,
            total_iterations: 0,
        },
        &mut callback,
    );
    for ((&value, mode), &weight) in values.iter().zip(&modes).zip(&weights) {
        let mode = scatter_full(mode, &free, elasticity.n());
        let derivative = eigenvalue_reference_density_gradient(
            elasticity,
            &projected,
            value,
            &mode,
            &mut control,
        )
        .unwrap();
        for (g, d) in gradient.iter_mut().zip(derivative) {
            *g += weight * d;
        }
    }
    assert_eq!(control.work().linear_solves, 0);
    assert_eq!(control.work().linear_iterations, 0);
    let homogeneity = -values.iter().zip(weights).map(|(v, w)| v * w).sum::<f64>();
    (aggregate, gradient, homogeneity)
}

#[test]
fn physical_density_modal_gradients_match_both_solvers_and_original_field_perturbations() {
    let (mesh, positions) = mesh();
    let source = material(99.0, 17.0, 0.28, 12);
    for (uniform, base, beta) in [(true, 0.47, 3.0), (false, 0.055, 2.0)] {
        let density: Vec<_> = positions
            .iter()
            .map(|p| {
                if uniform {
                    4.0
                } else {
                    2.0 + 3.0 * p[0] + p[1]
                }
            })
            .collect();
        let build = |values: &[f64]| {
            DensityElasticity::with_material_field_and_nodal_reference_density(
                &mesh,
                &positions,
                TetMaterialField::Uniform(&source),
                values,
                &|p| p[0] == 0.0,
            )
            .unwrap()
        };
        let pipeline = pipeline(&mesh, &positions, beta);
        let occupancy: Vec<_> = (0..mesh.tets.len())
            .map(|i| base + 0.002 * i as f64)
            .collect();
        let mut elasticity = build(&density);
        let (value, gradient, homogeneous) =
            dense_pullback(&pipeline, &mut elasticity, &occupancy, 3);
        let native = matrix_free(&pipeline, &mut elasticity, &occupancy, 3);
        let actual = native.reference_density_gradient.as_ref().unwrap();
        close(native.aggregate, value, 1e-7);
        vectors_close(actual, &gradient, 3e-6);
        assert_eq!(actual.len(), positions.len());
        assert!(actual.iter().all(|v| *v <= 0.0));
        close(
            gradient.iter().zip(&density).map(|(g, r)| g * r).sum(),
            homogeneous,
            2e-9,
        );

        // Independent original-field perturbations rebuild the true mass tiles,
        // so this detects differentiating the mean while freezing those tiles.
        let direction: Vec<_> = (0..density.len())
            .map(|i| {
                if i % 2 == 0 {
                    0.4 + 0.03 * i as f64
                } else {
                    -0.3
                }
            })
            .collect();
        let step = 1e-4;
        let plus: Vec<_> = density
            .iter()
            .zip(&direction)
            .map(|(r, d)| r + step * d)
            .collect();
        let minus: Vec<_> = density
            .iter()
            .zip(&direction)
            .map(|(r, d)| r - step * d)
            .collect();
        let high = eigenfrequency_objective(&pipeline, &mut build(&plus), &occupancy, 3, 2.0).0;
        let low = eigenfrequency_objective(&pipeline, &mut build(&minus), &occupancy, 3, 2.0).0;
        let analytic: f64 = actual.iter().zip(direction).map(|(g, d)| g * d).sum();
        close(analytic, (high - low) / (2.0 * step), 3e-4);
        assert_eq!(elasticity.nodal_reference_density_kg_m3().unwrap(), density);
        assert_eq!(elasticity.material_density_kg_m3(0), Some(99.0));
    }
}

#[test]
fn physical_inventory_and_quadratic_mass_derivatives_include_clamped_density_nodes() {
    let (mesh, positions) = mesh();
    let source = material(99.0, 17.0, 0.28, 13);
    let density: Vec<_> = positions.iter().map(|p| 2.0 + p[0] + 0.3 * p[2]).collect();
    let build = |values: &[f64]| {
        DensityElasticity::with_material_field_and_nodal_reference_density(
            &mesh,
            &positions,
            TetMaterialField::Uniform(&source),
            values,
            &|p| p[0] == 0.0,
        )
        .unwrap()
    };
    let elasticity = build(&density);
    let occupancy: Vec<_> = (0..mesh.tets.len())
        .map(|i| 0.02 + 0.008 * i as f64)
        .collect();
    let translation: Vec<_> = positions.iter().flat_map(|_| [1.0, 0.0, 0.0]).collect();
    let mut motion: Vec<_> = positions
        .iter()
        .flat_map(|p| [1.0 + p[0], 2.0 * p[1] - p[2], p[0] - p[1]])
        .collect();
    // cell_kinetic is the reduced physical form and masks fixed coordinates.
    // Use that same field for its independent mass-tile comparison below.
    for (dof, &free) in elasticity.free().iter().enumerate() {
        if !free {
            motion[dof] = 0.0;
        }
    }
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let inventory = elasticity
        .reference_density_inventory_gradient(&occupancy, &mut control)
        .unwrap();
    let translated = elasticity
        .reference_density_mass_gradient(&occupancy, &translation, &mut control)
        .unwrap();
    vectors_close(&inventory, &translated, 2e-14);
    for (node, p) in positions.iter().enumerate() {
        if p[0] == 0.0 {
            assert!(inventory[node] > 0.0);
        }
    }
    let total: f64 = elasticity
        .reference_cell_masses_kg()
        .unwrap()
        .iter()
        .zip(&occupancy)
        .map(|(m, r)| m * r)
        .sum();
    close(
        inventory.iter().zip(&density).map(|(g, r)| g * r).sum(),
        total,
        2e-14,
    );
    let gradient = elasticity
        .reference_density_mass_gradient(&occupancy, &motion, &mut control)
        .unwrap();
    let form = |el: &DensityElasticity| {
        el.cell_kinetic(&motion)
            .iter()
            .enumerate()
            .map(|(cell, q)| q * el.reference_density_kg_m3(cell).unwrap() * occupancy[cell])
            .sum::<f64>()
    };
    close(
        gradient.iter().zip(&density).map(|(g, r)| g * r).sum(),
        form(&elasticity),
        3e-14,
    );
    for node in 0..density.len() {
        let step = 1e-4;
        let mut plus = density.clone();
        let mut minus = density.clone();
        plus[node] += step;
        minus[node] -= step;
        close(
            gradient[node],
            (form(&build(&plus)) - form(&build(&minus))) / (2.0 * step),
            2e-9,
        );
    }
    assert_eq!(control.work().linear_solves, 0);
}

#[test]
fn absent_density_and_invalid_fields_refuse_and_cancellation_restores_modal_state() {
    let (mesh, positions) = mesh();
    let source = material(2.4, 17.0, 0.28, 14);
    let legacy = DensityElasticity::with_material_field(
        &mesh,
        &positions,
        TetMaterialField::Uniform(&source),
        &|p| p[0] == 0.0,
    )
    .unwrap();
    let density = vec![2.4; positions.len()];
    let mut elasticity = DensityElasticity::with_material_field_and_nodal_reference_density(
        &mesh,
        &positions,
        TetMaterialField::Uniform(&source),
        &density,
        &|p| p[0] == 0.0,
    )
    .unwrap();
    let occupancy = vec![0.5; mesh.tets.len()];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    assert_eq!(
        legacy.reference_density_inventory_gradient(&occupancy, &mut control),
        Err(ReferenceDensitySensitivityError::MissingNodalDensity)
    );
    assert!(matches!(
        elasticity.reference_density_inventory_gradient(&[1.1; 6], &mut control),
        Err(ReferenceDensitySensitivityError::InvalidInput(_))
    ));
    assert!(matches!(
        elasticity.reference_density_mass_gradient(&occupancy, &[0.0], &mut control),
        Err(ReferenceDensitySensitivityError::InvalidInput(_))
    ));
    assert!(matches!(
        elasticity.reference_density_mass_gradient(
            &occupancy,
            &vec![f64::NAN; elasticity.n()],
            &mut control
        ),
        Err(ReferenceDensitySensitivityError::Nonfinite { .. })
    ));
    assert_eq!(control.work().linear_solves, 0);

    elasticity.moduli.fill(0.91);
    let previous = elasticity.moduli.clone();
    let mass = elasticity.reference_cell_masses_kg().unwrap().to_vec();
    let mut stop = |progress: fs_topopt::SolveProgress| {
        if progress.stage == "reference-density-mass" {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut stop);
    let result = controlled_matrix_free_eigenfrequency_objective(
        &pipeline(&mesh, &positions, 3.0),
        &mut elasticity,
        &occupancy,
        EigenfrequencyObjectiveOptions {
            solver: MatrixFreeEigenOptions {
                count: 1,
                oversampling: 12,
                max_iterations: 32,
                ..Default::default()
            },
            ..Default::default()
        },
        &mut control,
    );
    assert!(matches!(
        result,
        Err(MatrixFreeEigenError::Operator(
            fs_topopt::EvaluationStop::Cancelled
        ))
    ));
    assert_eq!(elasticity.moduli, previous);
    assert_eq!(elasticity.reference_cell_masses_kg().unwrap(), mass);
    assert_eq!(elasticity.nodal_reference_density_kg_m3().unwrap(), density);
}
