//! Native material tensors/densities share the existing topology operator.
use fs_feec::whitney::ElementGeometryError;
use fs_rep_mesh::TetComplex;
use fs_solid::linear3::{TetElasticError, TetElasticMaterial, TetMaterialField};
use fs_solver::LinearOp;
use fs_topopt::elasticity::{DensityElasticity, MaterialElasticityError};

fn single_tet() -> (TetComplex, Vec<[f64; 3]>) {
    (
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
}

fn close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1e-12);
    assert!(
        (actual - expected).abs() <= 2e-12 * scale,
        "actual {actual:e}, expected {expected:e}"
    );
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn affine_displacement(positions: &[[f64; 3]]) -> Vec<f64> {
    let gradient = [[0.01, 0.02, -0.03], [0.04, -0.01, 0.05], [0.06, 0.07, 0.02]];
    positions
        .iter()
        .flat_map(|position| gradient.map(|row| dot(&row, position)))
        .collect()
}

#[test]
fn full_coupled_mandel_tensor_preserves_physical_stress_and_energy() {
    let (complex, positions) = single_tet();
    // Strictly diagonally dominant and symmetric: normal/shear as well as
    // shear/shear coupling exposes scaling only the rows or only the columns.
    let diagonal = [12.0, 14.0, 16.0, 8.0, 9.0, 10.0];
    let mut stiffness = [[0.0; 6]; 6];
    for (i, row) in stiffness.iter_mut().enumerate() {
        row[i] = diagonal[i];
    }
    for (i, j, value) in [(0, 3, 1.4), (1, 5, -0.8), (2, 4, 0.6), (3, 4, 0.7)] {
        stiffness[i][j] = value;
        stiffness[j][i] = value;
    }
    let material = TetElasticMaterial::try_new_mandel(6.0, stiffness, [9; 32]).unwrap();
    let operator = DensityElasticity::with_material_field(
        &complex,
        &positions,
        TetMaterialField::Uniform(&material),
        &|_| false,
    )
    .unwrap();
    let displacement = affine_displacement(&positions);
    let inv_sqrt_two = core::f64::consts::FRAC_1_SQRT_2;
    let strain_mandel = [
        0.01,
        -0.01,
        0.02,
        0.06 * inv_sqrt_two,
        0.12 * inv_sqrt_two,
        0.03 * inv_sqrt_two,
    ];
    let stress_mandel: [f64; 6] = core::array::from_fn(|row| dot(&stiffness[row], &strain_mandel));
    let actual_stress = operator.cell_stress(&displacement)[0];
    for component in 0..6 {
        let expected = stress_mandel[component] * if component < 3 { 1.0 } else { inv_sqrt_two };
        close(actual_stress[component], expected);
    }
    let expected_quadratic = dot(&strain_mandel, &stress_mandel) / 6.0;
    close(operator.cell_energies(&displacement)[0], expected_quadratic);
    let mut force = vec![0.0; displacement.len()];
    operator.apply(&displacement, &mut force);
    close(dot(&displacement, &force), expected_quadratic);
    match operator.material_field().unwrap() {
        TetMaterialField::Uniform(retained) => assert_eq!(retained, &material),
        TetMaterialField::PerElement(_) => panic!("uniform field must remain uniform"),
    }
}

#[test]
fn isotropic_path_matches_legacy_and_mass_blocks_remain_unit_density() {
    let (complex, positions) = single_tet();
    let material = TetElasticMaterial::try_new(6.0, 120.0, 0.27, [7; 32]).unwrap();
    let legacy = DensityElasticity::new(&complex, &positions, 120.0, 0.27, &|_| false);
    let operator = DensityElasticity::with_material_field(
        &complex,
        &positions,
        TetMaterialField::Uniform(&material),
        &|_| false,
    )
    .unwrap();
    let displacement = affine_displacement(&positions);
    let mut actual = vec![0.0; operator.n()];
    let mut expected = vec![0.0; operator.n()];
    operator.apply(&displacement, &mut actual);
    legacy.apply(&displacement, &mut expected);
    for (&a, &b) in actual.iter().zip(&expected) {
        close(a, b);
    }
    for (&a, &b) in operator.cell_stress(&displacement)[0]
        .iter()
        .zip(&legacy.cell_stress(&displacement)[0])
    {
        close(a, b);
    }
    let translation: Vec<f64> = (0..operator.n())
        .map(|dof| if dof % 3 == 0 { 1.0 } else { 0.0 })
        .collect();
    close(operator.cell_kinetic(&translation)[0], 1.0 / 6.0);
    operator.apply_mass(&[1.0], &translation, &mut actual);
    close(dot(&translation, &actual), 1.0 / 6.0);
    operator.apply_mass(&[material.density_kg_m3], &translation, &mut actual);
    close(dot(&translation, &actual), 1.0);
    close(operator.reference_cell_masses_kg().unwrap()[0], 1.0);
    assert_eq!(operator.material_density_kg_m3(0), Some(6.0));
    assert!(legacy.material_field().is_none());
    assert!(legacy.material_density_kg_m3(0).is_none());
    assert!(legacy.reference_cell_masses_kg().is_none());
}

#[test]
fn ordered_field_is_owned_and_full_reference_mass_includes_fixed_nodes() {
    let (complex, positions) = fs_feec::kuhn_cube(1);
    let mut materials: Vec<_> = (0..complex.tets.len())
        .map(|cell| {
            TetElasticMaterial::try_new(
                2.0 + cell as f64,
                20.0 + cell as f64,
                0.2,
                [cell as u8 + 1; 32],
            )
            .unwrap()
        })
        .collect();
    let originals = materials.clone();
    let operator = DensityElasticity::with_material_field(
        &complex,
        &positions,
        TetMaterialField::PerElement(&materials),
        &|point| point[0] == 0.0,
    )
    .unwrap();
    materials[0].density_kg_m3 = 999.0;
    materials[0].material_state_identity = [99; 32];
    assert!(operator.free().iter().any(|free| !free));
    match operator.material_field().unwrap() {
        TetMaterialField::PerElement(retained) => assert_eq!(retained, originals.as_slice()),
        TetMaterialField::Uniform(_) => panic!("ordered heterogeneous field must remain ordered"),
    }
    let masses = operator.reference_cell_masses_kg().unwrap();
    for (cell, material) in originals.iter().enumerate() {
        assert_eq!(
            operator.material_density_kg_m3(cell),
            Some(material.density_kg_m3)
        );
        close(masses[cell], material.density_kg_m3 / 6.0);
    }
    close(
        masses.iter().sum(),
        originals.iter().map(|m| m.density_kg_m3 / 6.0).sum(),
    );
    assert_eq!(operator.material_density_kg_m3(operator.cells()), None);
}

#[test]
fn incomplete_and_mutated_material_fields_are_refused() {
    let (complex, positions) = single_tet();
    let material = TetElasticMaterial::try_new(6.0, 120.0, 0.27, [7; 32]).unwrap();
    assert!(matches!(
        DensityElasticity::with_material_field(
            &complex,
            &positions,
            TetMaterialField::PerElement(&[]),
            &|_| false,
        ),
        Err(MaterialElasticityError::Solid(
            TetElasticError::MaterialCountMismatch {
                expected: 1,
                got: 0,
            }
        ))
    ));
    for density in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
        let mut invalid = material.clone();
        invalid.density_kg_m3 = density;
        assert!(matches!(
            DensityElasticity::with_material_field(
                &complex,
                &positions,
                TetMaterialField::Uniform(&invalid),
                &|_| false,
            ),
            Err(MaterialElasticityError::Solid(
                TetElasticError::InvalidMaterial { .. }
            ))
        ));
    }
    let mut unidentified = material;
    unidentified.material_state_identity = [0; 32];
    assert!(matches!(
        DensityElasticity::with_material_field(
            &complex,
            &positions,
            TetMaterialField::Uniform(&unidentified),
            &|_| false,
        ),
        Err(MaterialElasticityError::Solid(
            TetElasticError::InvalidMaterial { .. }
        ))
    ));
}

#[test]
fn malformed_mesh_and_unrepresentable_physical_mass_are_refused() {
    let (complex, positions) = single_tet();
    let material = TetElasticMaterial::try_new(6.0, 120.0, 0.27, [7; 32]).unwrap();
    let build = |mesh: &TetComplex, points: &[[f64; 3]]| {
        DensityElasticity::with_material_field(
            mesh,
            points,
            TetMaterialField::Uniform(&material),
            &|_| false,
        )
    };
    assert!(matches!(
        build(&complex, &positions[..3]),
        Err(MaterialElasticityError::Geometry(
            ElementGeometryError::Positions
        ))
    ));
    let mut nonfinite = positions.clone();
    nonfinite[0][0] = f64::NAN;
    assert!(matches!(
        build(&complex, &nonfinite),
        Err(MaterialElasticityError::Geometry(
            ElementGeometryError::Positions
        ))
    ));
    for connectivity in [[0, 1, 1, 3], [0, 1, 2, 4]] {
        let mut malformed = complex.clone();
        malformed.tets[0] = connectivity;
        assert!(matches!(
            build(&malformed, &positions),
            Err(MaterialElasticityError::Geometry(
                ElementGeometryError::Cell(0)
            ))
        ));
    }
    let mut singular = positions.clone();
    singular[3] = [0.25, 0.25, 0.0];
    assert!(matches!(
        build(&complex, &singular),
        Err(MaterialElasticityError::Geometry(
            ElementGeometryError::Cell(0)
        ))
    ));
    let mut unreferenced = complex.clone();
    unreferenced.vertex_count += 1;
    let mut extra_node = positions.clone();
    extra_node.push([2.0, 2.0, 2.0]);
    assert!(matches!(
        build(&unreferenced, &extra_node),
        Err(MaterialElasticityError::Solid(
            TetElasticError::UnreferencedNode { node: 4 }
        ))
    ));
    assert!(matches!(
        build(&TetComplex::from_tets(0, vec![]), &[]),
        Err(MaterialElasticityError::Solid(TetElasticError::EmptyMesh))
    ));
    let huge_density = TetElasticMaterial::try_new(f64::MAX, 120.0, 0.27, [8; 32]).unwrap();
    let large_positions: Vec<_> = positions.iter().map(|p| p.map(|v| v * 10.0)).collect();
    assert!(matches!(
        DensityElasticity::with_material_field(
            &complex,
            &large_positions,
            TetMaterialField::Uniform(&huge_density),
            &|_| false,
        ),
        Err(MaterialElasticityError::UnrepresentableElement { cell: 0 })
    ));
}
