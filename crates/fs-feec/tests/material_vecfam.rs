//! Physical graded inertia in original H(curl)/H(div) moment coordinates.
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_feec::highorder::vecfam::{Family, VecSpace, VectorMaterialMass};
use fs_feec::{
    CellWeight, WeightedAssemblyLimits, WeightedError, element_geometry, single_tet, two_tets,
};
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;

const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const LIMITS: WeightedAssemblyLimits = WeightedAssemblyLimits {
    max_cells: 100,
    max_dofs: 10_000,
    max_triplets: 1_000_000,
};
const TENSOR: [[f64; 3]; 3] = [[2.0, 0.3, -0.2], [0.3, 3.0, 0.4], [-0.2, 0.4, 4.0]];

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-8 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn apply(mass: &VectorMaterialMass<'_, '_>, input: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; mass.ndof()];
    mass.apply_into(input, &mut out, &mut vec![0.0; mass.ndof()], &mut || false)
        .unwrap();
    out
}
fn profile(complex: &TetComplex, positions: &[[f64; 3]]) -> Vec<CellWeight> {
    complex
        .tets
        .iter()
        .map(|tet| CellWeight::AffineScalar(tet.map(|v| 3.0 + 12.0 * positions[v as usize][0])))
        .collect()
}
fn affine_field([x, y, z]: [f64; 3]) -> [f64; 3] {
    [1.0 + x, 2.0 - y, -0.5 + 2.0 * z]
}

#[test]
fn graded_vector_mass_preserves_exact_inventory_centroid_and_energy() {
    let (complex, positions) = single_tet();
    let weights = profile(&complex, &positions);
    for family in [Family::Nedelec, Family::Rt] {
        for r in 1..=4 {
            let space = VecSpace::new(&complex, &positions, r, family);
            let mass = space
                .prepare_material_mass(
                    &positions,
                    &weights,
                    DENSITY,
                    LIMITS,
                    1_000_000,
                    &mut || false,
                )
                .unwrap();
            let constant = space.interpolate(&positions, &|_| [1.0, 2.0, -0.5]);
            // integral rho = 1 on the unit tetrahedron; |constant|² = 21/4.
            near(dot(&constant, &apply(&mass, &constant)), 21.0 / 4.0);
            assert_eq!(
                mass.energy_dims(Dims([1, 0, -1, 0, 0, 0])).unwrap(),
                Dims([2, 1, -2, 0, 0, 0])
            );
            if r >= 2 {
                let one = space.interpolate(&positions, &|_| [1.0, 0.0, 0.0]);
                let load = apply(&mass, &one);
                for axis in 0..3 {
                    let coordinate = space.interpolate(&positions, &|p| [p[axis], 0.0, 0.0]);
                    near(
                        dot(&coordinate, &load),
                        [13.0 / 40.0, 9.0 / 40.0, 9.0 / 40.0][axis],
                    );
                }
                let field = space.interpolate(&positions, &affine_field);
                // Independent simplex monomial moments of rho*|field|².
                near(dot(&field, &apply(&mass, &field)), 307.0 / 60.0);
            }
        }
    }
}

#[test]
fn polynomial_action_and_diagonal_match_consistent_scalar_affine_and_tensor_assembly() {
    let (complex, positions) = two_tets();
    for family in [Family::Nedelec, Family::Rt] {
        for r in [1, 3, 4] {
            let space = VecSpace::new(&complex, &positions, r, family);
            let sets = [
                [CellWeight::Scalar(1.0); 2],
                [
                    CellWeight::AffineScalar([2.0, 3.0, 7.0, 5.0]),
                    CellWeight::Scalar(11.0),
                ],
                [CellWeight::Tensor(TENSOR), CellWeight::Scalar(0.7)],
            ];
            for weights in sets {
                let mass = space
                    .prepare_material_mass(
                        &positions,
                        &weights,
                        DENSITY,
                        LIMITS,
                        1_000_000,
                        &mut || false,
                    )
                    .unwrap();
                let matrix = mass.assemble(&mut || false).unwrap();
                let dense = matrix.to_dense();
                for (i, &diagonal) in mass.diagonal().iter().enumerate() {
                    assert_eq!(diagonal.to_bits(), dense[i * space.ndof + i].to_bits());
                }
                let input: Vec<_> = (0..space.ndof)
                    .map(|i| (i % 13) as f64 / 7.0 - 0.6)
                    .collect();
                let actual = apply(&mass, &input);
                let mut expected = vec![0.0; space.ndof];
                matrix.spmv(&input, &mut expected);
                let scale = expected.iter().fold(1.0_f64, |s, v| s.max(v.abs()));
                for (a, b) in actual.iter().zip(&expected) {
                    assert!(
                        (a - b).abs() < 3e-9 * scale,
                        "{a:e} vs {b:e}, scale {scale:e}"
                    );
                }
                assert!(dot(&input, &actual) > 0.0);
                assert_eq!(actual, apply(&mass, &input));
                if weights == [CellWeight::Scalar(1.0); 2] {
                    let original = space.mass(&positions).to_dense();
                    for (a, b) in dense.iter().zip(original) {
                        near(*a, b);
                    }
                }
            }
        }
    }
}

#[test]
fn spatial_force_density_recovers_physical_acceleration_in_native_moment_solve() {
    let (complex, positions) = two_tets();
    let geometry = element_geometry(&complex, &positions);
    // A material jump at the shared face must weight both forcing and inertia.
    let weights: Vec<_> = complex
        .tets
        .iter()
        .enumerate()
        .map(|(cell, tet)| {
            CellWeight::AffineScalar(
                tet.map(|v| (cell + 1) as f64 * (3.0 + 12.0 * positions[v as usize][0])),
            )
        })
        .collect();
    for family in [Family::Nedelec, Family::Rt] {
        let space = VecSpace::new(&complex, &positions, 2, family);
        let mass = space
            .prepare_material_mass(
                &positions,
                &weights,
                DENSITY,
                LIMITS,
                1_000_000,
                &mut || false,
            )
            .unwrap();
        let matrix = mass.assemble(&mut || false).unwrap();
        let mut force = vec![0.0; space.ndof];
        // Independent direct load quadrature, using the physical force rho*a,
        // not the matrix applied to the expected discrete acceleration.
        for (cell, tet) in complex.tets.iter().enumerate() {
            let dofs = space.element_dofs(cell);
            let element = &space.elements[cell];
            for (lambda, weight) in duffy_quadrature(7) {
                let point = core::array::from_fn(|c| {
                    (0..4)
                        .map(|a| lambda[a] * positions[tet[a] as usize][c])
                        .sum()
                });
                let acceleration = affine_field(point);
                let rho = (cell + 1) as f64 * (3.0 + 12.0 * point[0]);
                let scale = weight * 6.0 * geometry.vol_signed[cell].abs() * rho;
                for (i, &global) in dofs.iter().enumerate() {
                    let shape =
                        element.funcs[i].eval_local(&element.monos, element.chart.local(point));
                    force[global] += scale * dot(&shape, &acceleration);
                }
            }
        }
        let mut solved = force.clone();
        fs_la::factor::cholesky(&matrix.to_dense(), space.ndof)
            .unwrap()
            .solve(&mut solved);
        let mut residual = vec![0.0; space.ndof];
        matrix.spmv(&solved, &mut residual);
        for (a, b) in residual.iter().zip(&force) {
            near(*a, *b);
        }
        for cell in 0..complex.tets.len() {
            for lambda in [[0.25; 4], [0.1, 0.2, 0.3, 0.4]] {
                let point = core::array::from_fn(|c| {
                    (0..4)
                        .map(|a| lambda[a] * positions[complex.tets[cell][a] as usize][c])
                        .sum()
                });
                let observed = space.eval_in(cell, &solved, point);
                for (a, b) in observed.into_iter().zip(affine_field(point)) {
                    near(a, b);
                }
            }
        }
    }
}

#[test]
fn native_moment_units_match_geometric_scaling_of_every_matrix_block() {
    let (complex, positions) = single_tet();
    let doubled: Vec<_> = positions.iter().map(|p| p.map(|x| 2.0 * x)).collect();
    let weights = [CellWeight::Scalar(3.0)];
    let velocity = Dims([1, 0, -1, 0, 0, 0]);
    for family in [Family::Nedelec, Family::Rt] {
        let space = VecSpace::new(&complex, &positions, 3, family);
        let scaled_space = VecSpace::new(&complex, &doubled, 3, family);
        let mass = space
            .prepare_material_mass(
                &positions,
                &weights,
                DENSITY,
                LIMITS,
                1_000_000,
                &mut || false,
            )
            .unwrap();
        let scaled_mass = scaled_space
            .prepare_material_mass(&doubled, &weights, DENSITY, LIMITS, 1_000_000, &mut || {
                false
            })
            .unwrap();
        let original = mass.assemble(&mut || false).unwrap().to_dense();
        let scaled = scaled_mass.assemble(&mut || false).unwrap().to_dense();
        for i in 0..space.ndof {
            let power_i = i8::from(family == Family::Nedelec && i < space.cell_off);
            assert_eq!(
                mass.coordinate_dims(i, velocity).unwrap(),
                Dims([1 + power_i, 0, -1, 0, 0, 0])
            );
            assert_eq!(
                mass.output_dims(i, velocity).unwrap(),
                Dims([1 - power_i, 1, -1, 0, 0, 0])
            );
            for j in 0..space.ndof {
                let power_j = i8::from(family == Family::Nedelec && j < space.cell_off);
                assert_eq!(
                    mass.entry_dims(i, j).unwrap(),
                    Dims([-power_i - power_j, 1, 0, 0, 0, 0])
                );
                near(
                    scaled[i * space.ndof + j],
                    original[i * space.ndof + j] * 2.0_f64.powi(i32::from(3 - power_i - power_j)),
                );
            }
        }
        assert!(matches!(
            mass.output_dims(space.ndof, velocity),
            Err(WeightedError::InvalidGeometry { .. })
        ));
    }
}

#[test]
fn graded_physical_energy_survives_vertex_relabeling_and_anisotropic_materials() {
    let (complex, positions) = two_tets();
    let permutation = [3, 0, 4, 1, 2];
    let mut relabeled_positions = vec![[0.0; 3]; positions.len()];
    for (old, &new) in permutation.iter().enumerate() {
        relabeled_positions[new] = positions[old];
    }
    let relabeled = TetComplex::from_tets(
        positions.len(),
        complex
            .tets
            .iter()
            .map(|t| t.map(|v| permutation[v as usize] as u32))
            .collect(),
    );
    let weights = profile(&complex, &positions);
    let relabeled_weights = profile(&relabeled, &relabeled_positions);
    for family in [Family::Nedelec, Family::Rt] {
        let a = VecSpace::new(&complex, &positions, 3, family);
        let b = VecSpace::new(&relabeled, &relabeled_positions, 3, family);
        let ma = a
            .prepare_material_mass(
                &positions,
                &weights,
                DENSITY,
                LIMITS,
                1_000_000,
                &mut || false,
            )
            .unwrap();
        let mb = b
            .prepare_material_mass(
                &relabeled_positions,
                &relabeled_weights,
                DENSITY,
                LIMITS,
                1_000_000,
                &mut || false,
            )
            .unwrap();
        let ua = a.interpolate(&positions, &affine_field);
        let ub = b.interpolate(&relabeled_positions, &affine_field);
        near(dot(&ua, &apply(&ma, &ua)), dot(&ub, &apply(&mb, &ub)));
        let tensor_weights = [CellWeight::Tensor(TENSOR); 2];
        let tensor_mass = a
            .prepare_material_mass(
                &positions,
                &tensor_weights,
                DENSITY,
                LIMITS,
                1_000_000,
                &mut || false,
            )
            .unwrap();
        let value = [1.0, 2.0, -0.5];
        let constant = a.interpolate(&positions, &|_| value);
        let volume: f64 = element_geometry(&complex, &positions)
            .vol_signed
            .iter()
            .map(|v| v.abs())
            .sum();
        let quadratic: f64 = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| value[i] * TENSOR[i][j] * value[j])
                    .sum::<f64>()
            })
            .sum();
        near(
            dot(&constant, &apply(&tensor_mass, &constant)),
            volume * quadratic,
        );
    }
}

#[test]
fn bounded_preparation_retains_matrix_free_action_and_atomic_refusal() {
    let (complex, positions) = single_tet();
    let mut space = VecSpace::new(&complex, &positions, 3, Family::Nedelec);
    let weights = profile(&complex, &positions);
    assert!(matches!(
        space.prepare_material_mass(&positions, &weights, DENSITY, LIMITS, 0, &mut || false),
        Err(WeightedError::LimitExceeded { .. })
    ));
    {
        let mass = space
            .prepare_material_mass(
                &positions,
                &weights,
                DENSITY,
                WeightedAssemblyLimits {
                    max_triplets: 0,
                    ..LIMITS
                },
                1_000_000,
                &mut || false,
            )
            .unwrap();
        assert!(matches!(
            mass.assemble(&mut || false),
            Err(WeightedError::LimitExceeded {
                resource: "triplets",
                ..
            })
        ));
        let mut input = space.interpolate(&positions, &|_| [1.0, 0.0, 0.0]);
        near(dot(&input, &apply(&mass, &input)), 1.0);
        let mut output = vec![42.0; space.ndof];
        let original = output.clone();
        let mut scratch = vec![0.0; space.ndof];
        let mut polls = 0;
        assert_eq!(
            mass.apply_into(&input, &mut output, &mut scratch, &mut || {
                polls += 1;
                polls == 8
            }),
            Err(WeightedError::Cancelled)
        );
        assert_eq!(output, original);
        input[0] = f64::NAN;
        assert_eq!(
            mass.apply_into(&input, &mut output, &mut scratch, &mut || false),
            Err(WeightedError::NonFiniteEntry)
        );
        assert_eq!(output, original);
    }
    for weight in [
        CellWeight::Scalar(0.0),
        CellWeight::AffineScalar([1.0, 2.0, -1.0, 4.0]),
        CellWeight::Tensor([[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    ] {
        assert!(matches!(
            space.prepare_material_mass(
                &positions,
                &[weight],
                DENSITY,
                LIMITS,
                1_000_000,
                &mut || false
            ),
            Err(WeightedError::InvalidCoefficient { .. })
        ));
    }
    let _ = space.elements[0].funcs[0].comp[0].pop();
    assert!(matches!(
        space.prepare_material_mass(
            &positions,
            &weights,
            DENSITY,
            LIMITS,
            1_000_000,
            &mut || false
        ),
        Err(WeightedError::InvalidGeometry { .. })
    ));
}
