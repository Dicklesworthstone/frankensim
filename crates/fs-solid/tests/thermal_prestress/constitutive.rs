//! Constitutive elimination changes the derivative, not the accepted stress.
use super::*;

#[test]
fn replacement_tangent_keeps_original_stress_and_inertial_inventory() {
    with_cx(|cx| {
        let mat = material();
        let thermal = strain(0.01);
        let p = problem(&mat, &TETS, &FIXED);
        let field = TetThermalStrainField::Uniform(&thermal);
        let u = p
            .solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
            .unwrap();
        let original = p
            .assemble_thermal_tangent(&u, field, ContentHash([1; 32]), 0.1, cx)
            .unwrap();
        let tangent = TetElasticMaterial::try_new(6.0, 2000.0, 0.25, [42; 32]).unwrap();
        let assemble = || {
            p.assemble_thermal_tangent_with_constitutive_materials(
                &u,
                field,
                &[],
                TetMaterialField::Uniform(&tangent),
                ContentHash([1; 32]),
                0.1,
                cx,
            )
            .unwrap()
        };
        let result = assemble();
        assert_eq!(result.stress(), original.stress());
        assert_eq!(
            result.geometric_stiffness().to_dense(),
            original.geometric_stiffness().to_dense()
        );
        assert_eq!(
            result.pencil().mass.to_dense(),
            original.pencil().mass.to_dense()
        );
        assert_eq!(result.pencil().free_dofs, original.pencil().free_dofs);
        close(result.pencil().total_mass_kg, 1.0);
        let n = result.pencil().free_dofs.len();
        for row in 0..n {
            for column in 0..n {
                let material = 2.0 * original.elastic_stiffness().get(row, column);
                close(result.elastic_stiffness().get(row, column), material);
                close(
                    result.pencil().stiffness.get(row, column),
                    material + original.geometric_stiffness().get(row, column),
                );
            }
        }
        assert_ne!(result.identity(), original.identity());
        assert_eq!(result.identity(), assemble().identity());
        // Reusing a source label cannot hide different numerical derivatives.
        let different = TetElasticMaterial::try_new(6.0, 1700.0, 0.25, [42; 32]).unwrap();
        let different = p
            .assemble_thermal_tangent_with_constitutive_materials(
                &u,
                field,
                &[],
                TetMaterialField::Uniform(&different),
                ContentHash([1; 32]),
                0.1,
                cx,
            )
            .unwrap();
        assert_ne!(result.identity(), different.identity());
    });
}

// Independent initial-stress potential in Green strain, evaluated under an
// affine virtual displacement on the unit tet. Its Hessian at the origin has
// the declared small-displacement K_material + K_geometric terms. This oracle
// never calls either assembly kernel or its material matrix.
fn virtual_energy(scale: f64, gradient: [[f64; 3]; 3], young: f64, nu: f64, sigma_xx: f64) -> f64 {
    let green: [[f64; 3]; 3] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            0.5 * scale * (gradient[i][j] + gradient[j][i])
                + 0.5 * scale * scale * (0..3).map(|k| gradient[k][i] * gradient[k][j]).sum::<f64>()
        })
    });
    let trace = (0..3).map(|i| green[i][i]).sum::<f64>();
    let lambda = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = young / (2.0 * (1.0 + nu));
    (sigma_xx * green[0][0]
        + 0.5 * lambda * trace * trace
        + mu * green
            .iter()
            .flatten()
            .map(|value| value * value)
            .sum::<f64>())
        / 6.0
}

#[test]
fn full_pencil_matches_an_independent_initial_stress_potential_hessian() {
    with_cx(|cx| {
        let mat = material();
        let thermal = strain(0.01);
        let p = problem(&mat, &TETS, &FIXED);
        let field = TetThermalStrainField::Uniform(&thermal);
        let u = p
            .solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
            .unwrap();
        let tangent = TetElasticMaterial::try_new(6.0, 1500.0, 0.1, [43; 32]).unwrap();
        let result = p
            .assemble_thermal_tangent_with_constitutive_materials(
                &u,
                field,
                &[],
                TetMaterialField::Uniform(&tangent),
                ContentHash([1; 32]),
                0.1,
                cx,
            )
            .unwrap();
        let full = [
            [0.3, -0.1, 0.7],
            [0.2, 0.9, -0.4],
            [-0.3, 0.2, 0.6],
            [0.1, -0.6, 0.8],
        ];
        let mut admitted = [[0.0; 3]; 4];
        let q: Vec<_> = result
            .pencil()
            .free_dofs
            .iter()
            .map(|&dof| {
                admitted[dof / 3][dof % 3] = full[dof / 3][dof % 3];
                full[dof / 3][dof % 3]
            })
            .collect();
        let gradient =
            core::array::from_fn(|i| core::array::from_fn(|j| admitted[j + 1][i] - admitted[0][i]));
        let energy = |s| virtual_energy(s, gradient, 1500.0, 0.1, -10.0);
        let h = 0.001;
        let expected = (-energy(-2.0 * h) + 16.0 * energy(-h) - 30.0 * energy(0.0)
            + 16.0 * energy(h)
            - energy(2.0 * h))
            / (12.0 * h * h);
        let mut applied = vec![0.0; q.len()];
        result.pencil().stiffness.spmv(&q, &mut applied);
        let actual: f64 = q.iter().zip(applied).map(|(a, b)| a * b).sum();
        close(actual, expected);
    });
}

#[test]
fn compressive_instability_stays_signed_and_bad_material_or_history_refuses() {
    with_cx(|cx| {
        let mat = material();
        let thermal = strain(0.01);
        let p = problem(&mat, &TETS, &FIXED);
        let field = TetThermalStrainField::Uniform(&thermal);
        let u = p
            .solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
            .unwrap();
        // Explicit authored softened tangent at the same prestressed state.
        let soft = TetElasticMaterial::try_new(6.0, 0.1, 0.25, [44; 32]).unwrap();
        let result = p
            .assemble_thermal_tangent_with_constitutive_materials(
                &u,
                field,
                &[],
                TetMaterialField::Uniform(&soft),
                ContentHash([1; 32]),
                0.1,
                cx,
            )
            .unwrap();
        assert!(eigen(&result).iter().any(|value| *value < -1.0));
        let wrong_density = TetElasticMaterial::try_new(7.0, 1000.0, 0.25, [45; 32]).unwrap();
        assert!(
            p.assemble_thermal_tangent_with_constitutive_materials(
                &u,
                field,
                &[],
                TetMaterialField::Uniform(&wrong_density),
                ContentHash([1; 32]),
                0.1,
                cx,
            )
            .is_err()
        );
        assert!(matches!(
            p.assemble_thermal_tangent_with_constitutive_materials(
                &u,
                field,
                &[],
                TetMaterialField::PerElement(&[]),
                ContentHash([1; 32]),
                0.1,
                cx,
            ),
            Err(TetElasticError::MaterialCountMismatch { .. })
        ));
        assert!(matches!(
            p.assemble_thermal_tangent_with_constitutive_materials(
                &u,
                field,
                &[(4, 0.25)],
                TetMaterialField::Uniform(&soft),
                ContentHash([1; 32]),
                0.1,
                cx,
            ),
            Err(TetElasticError::StressRecoveryInputMismatch)
        ));
        let changed = strain(0.02);
        assert!(
            p.assemble_thermal_tangent_with_constitutive_materials(
                &u,
                TetThermalStrainField::Uniform(&changed),
                &[],
                TetMaterialField::Uniform(&soft),
                ContentHash([1; 32]),
                0.1,
                cx,
            )
            .is_err()
        );
        let retry = p
            .assemble_thermal_tangent_with_constitutive_materials(
                &u,
                field,
                &[],
                TetMaterialField::Uniform(&soft),
                ContentHash([1; 32]),
                0.1,
                cx,
            )
            .unwrap();
        assert_eq!(retry.identity(), result.identity());
    });
}
