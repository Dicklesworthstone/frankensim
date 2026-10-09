//! G0/G1/G3/G4: real material-dependent Whitney pairings and elliptic solves.

use fs_feec::{
    CellWeight, WeightedAssemblyLimits, WeightedError, deram0, deram1, deram2, element_geometry,
    incidence_to_csr, kuhn_cube, mass_matrix, on_unit_cube_boundary, single_tet, stiffness,
    two_tets, weighted_mass_matrix,
};
use fs_qty::Dims;
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;

const LIMITS: WeightedAssemblyLimits = WeightedAssemblyLimits {
    max_cells: 30_000,
    max_dofs: 100_000,
    max_triplets: 1_100_000,
};

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "actual {actual:e}, expected {expected:e}"
    );
}

fn quadratic(matrix: &Csr, values: &[f64]) -> f64 {
    let mut product = vec![0.0; values.len()];
    matrix.spmv(values, &mut product);
    values.iter().zip(product).map(|(a, b)| a * b).sum()
}

fn tensor_energy(tensor: [[f64; 3]; 3], vector: [f64; 3]) -> f64 {
    (0..3)
        .map(|i| {
            (0..3)
                .map(|j| vector[i] * tensor[i][j] * vector[j])
                .sum::<f64>()
        })
        .sum()
}

#[test]
fn unit_weights_reproduce_all_degrees_and_exact_symmetry() {
    let (complex, positions) = two_tets();
    let geo = element_geometry(&complex, &positions);
    for degree in 0..=3 {
        let weighted = weighted_mass_matrix(
            &complex,
            &geo,
            degree,
            &[CellWeight::Scalar(1.0); 2],
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap();
        let expected = mass_matrix(&complex, &geo, degree).to_dense();
        let actual = weighted.matrix().to_dense();
        for (&a, &b) in actual.iter().zip(&expected) {
            near(a, b);
        }
        let n = weighted.matrix().nrows();
        for i in 0..n {
            for j in 0..n {
                assert_eq!(actual[i * n + j].to_bits(), actual[j * n + i].to_bits());
            }
        }
        assert_eq!(weighted.degree(), degree);
        assert_eq!(
            weighted.entry_dims(),
            Dims([3 - 2 * degree as i8, 0, 0, 0, 0, 0])
        );
    }
}

#[test]
fn scalar_jumps_preserve_integrated_constant_fields_and_units() {
    let (complex, positions) = two_tets();
    let geo = element_geometry(&complex, &positions);
    let weights = [CellWeight::Scalar(2.0), CellWeight::Scalar(11.0)];
    let weighted_volume = 2.0 * geo.vol_signed[0].abs() + 11.0 * geo.vol_signed[1].abs();
    let coefficient_dims = Dims([1, 1, -3, -1, 0, 0]); // W/(m K)
    let vector = [0.7, -1.3, 2.1];
    let norm_squared: f64 = vector.iter().map(|v| v * v).sum();
    for degree in 0..=3 {
        let matrix = weighted_mass_matrix(
            &complex,
            &geo,
            degree,
            &weights,
            coefficient_dims,
            LIMITS,
            &mut || false,
        )
        .unwrap();
        let values = match degree {
            0 => vec![1.0; positions.len()],
            1 => deram1(&complex, &positions, &|_| vector),
            2 => deram2(&complex, &positions, &|_| vector),
            _ => geo.vol_signed.iter().map(|v| v.abs()).collect(),
        };
        near(
            quadratic(matrix.matrix(), &values),
            weighted_volume
                * if matches!(degree, 1 | 2) {
                    norm_squared
                } else {
                    1.0
                },
        );
        assert_eq!(matrix.coefficient_dims(), coefficient_dims);
        assert_eq!(
            matrix.output_dims(Dims([0, 0, 0, 1, 0, 0])).unwrap(),
            Dims([4 - 2 * degree as i8, 1, -3, 0, 0, 0])
        );
    }
}

#[test]
fn anisotropic_one_and_two_forms_match_independent_volume_energy() {
    let (complex, positions) = two_tets();
    let geo = element_geometry(&complex, &positions);
    let tensors = [
        [[2.0, 0.3, -0.2], [0.3, 3.0, 0.4], [-0.2, 0.4, 4.0]],
        [[7.0, -1.0, 0.7], [-1.0, 5.0, -0.3], [0.7, -0.3, 2.0]],
    ];
    let weights = tensors.map(CellWeight::Tensor);
    let vector = [0.7, -1.3, 2.1];
    let expected: f64 = tensors
        .iter()
        .zip(&geo.vol_signed)
        .map(|(&tensor, volume)| volume.abs() * tensor_energy(tensor, vector))
        .sum();
    for degree in [1, 2] {
        let mass = weighted_mass_matrix(
            &complex,
            &geo,
            degree,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap();
        let values = if degree == 1 {
            deram1(&complex, &positions, &|_| vector)
        } else {
            deram2(&complex, &positions, &|_| vector)
        };
        near(quadratic(mass.matrix(), &values), expected);
    }
}

#[test]
fn stored_cell_orientation_and_coordinate_scaling_preserve_pairings() {
    let (complex, positions) = single_tet();
    let reversed = TetComplex::from_tets(4, vec![[0, 2, 1, 3]]);
    let geo = element_geometry(&complex, &positions);
    let reversed_geo = element_geometry(&reversed, &positions);
    let scaled_positions: Vec<_> = positions.iter().map(|p| p.map(|x| 2.0 * x)).collect();
    let scaled_geo = element_geometry(&complex, &scaled_positions);
    for degree in 0..=3 {
        let weights = [CellWeight::Scalar(3.0)];
        let original = weighted_mass_matrix(
            &complex,
            &geo,
            degree,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap()
        .matrix()
        .to_dense();
        let reverse = weighted_mass_matrix(
            &reversed,
            &reversed_geo,
            degree,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap()
        .matrix()
        .to_dense();
        let scaled = weighted_mass_matrix(
            &complex,
            &scaled_geo,
            degree,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap()
        .matrix()
        .to_dense();
        for ((&a, &b), &c) in original.iter().zip(&reverse).zip(&scaled) {
            near(a, b);
            near(c, a * 2.0_f64.powi(3 - 2 * i32::from(degree)));
        }
    }
    // Weights never mutate the exact integer incidence sequence.
    assert!(
        complex
            .d1()
            .apply(&complex.d0().apply(&[1, 2, 4, 8]))
            .iter()
            .all(|&x| x == 0)
    );
}

#[test]
fn coefficient_geometry_and_extent_refusals_are_typed() {
    let (complex, positions) = single_tet();
    let geo = element_geometry(&complex, &positions);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            weighted_mass_matrix(
                &complex,
                &geo,
                1,
                &[CellWeight::Scalar(bad)],
                Dims::NONE,
                LIMITS,
                &mut || false
            ),
            Err(WeightedError::InvalidCoefficient { cell: 0, .. })
        ));
    }
    for bad in [
        [[1.0, 0.5, 0.0], [0.4, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, f64::NAN]],
    ] {
        assert!(matches!(
            weighted_mass_matrix(
                &complex,
                &geo,
                1,
                &[CellWeight::Tensor(bad)],
                Dims::NONE,
                LIMITS,
                &mut || false
            ),
            Err(WeightedError::InvalidCoefficient { .. })
        ));
    }
    let identity = CellWeight::Tensor([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    for degree in [0, 3] {
        assert!(matches!(
            weighted_mass_matrix(
                &complex,
                &geo,
                degree,
                &[identity],
                Dims::NONE,
                LIMITS,
                &mut || false
            ),
            Err(WeightedError::InvalidCoefficient { .. })
        ));
    }
    let weights = [CellWeight::Scalar(1.0)];
    assert!(matches!(
        weighted_mass_matrix(&complex, &geo, 4, &weights, Dims::NONE, LIMITS, &mut || {
            false
        }),
        Err(WeightedError::InvalidDegree(4))
    ));
    assert!(matches!(
        weighted_mass_matrix(&complex, &geo, 0, &[], Dims::NONE, LIMITS, &mut || false),
        Err(WeightedError::CoefficientCount { .. })
    ));
    for limits in [
        WeightedAssemblyLimits {
            max_cells: 0,
            ..LIMITS
        },
        WeightedAssemblyLimits {
            max_dofs: 3,
            ..LIMITS
        },
        WeightedAssemblyLimits {
            max_triplets: 15,
            ..LIMITS
        },
    ] {
        assert!(matches!(
            weighted_mass_matrix(&complex, &geo, 0, &weights, Dims::NONE, limits, &mut || {
                false
            }),
            Err(WeightedError::LimitExceeded { .. })
        ));
    }
    assert_eq!(
        weighted_mass_matrix(
            &complex,
            &geo,
            0,
            &weights,
            Dims([127, 0, 0, 0, 0, 0]),
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::DimensionOverflow)
    );
    let mut malformed_geo = element_geometry(&complex, &positions);
    malformed_geo.grads.clear();
    assert!(matches!(
        weighted_mass_matrix(
            &complex,
            &malformed_geo,
            1,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::InvalidGeometry { .. })
    ));
    let mut malformed_complex = complex.clone();
    malformed_complex.edges.clear();
    assert!(matches!(
        weighted_mass_matrix(
            &malformed_complex,
            &geo,
            1,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::InvalidGeometry { .. })
    ));
}

#[test]
fn scaled_tensor_admission_does_not_overflow_representable_mass() {
    let (complex, positions) = single_tet();
    let geo = element_geometry(&complex, &positions);
    for scale in [1e-250, 1e250, f64::MAX] {
        let weight = CellWeight::Tensor([[scale, 0.0, 0.0], [0.0, scale, 0.0], [0.0, 0.0, scale]]);
        let mass = weighted_mass_matrix(
            &complex,
            &geo,
            1,
            &[weight],
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap();
        for (&actual, expected) in mass
            .matrix()
            .to_dense()
            .iter()
            .zip(mass_matrix(&complex, &geo, 1).to_dense())
        {
            assert!(actual.is_finite());
            near(actual / scale, expected);
        }
    }
}

#[test]
fn cancellation_during_assembly_sort_and_publication_is_retryable() {
    let (complex, positions) = two_tets();
    let geo = element_geometry(&complex, &positions);
    let weights = [CellWeight::Scalar(2.0), CellWeight::Scalar(7.0)];
    let mut polls = 0;
    let expected =
        weighted_mass_matrix(&complex, &geo, 1, &weights, Dims::NONE, LIMITS, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, 2, 4, 6, 10, polls / 3, polls / 2, polls - 1, polls] {
        let mut seen = 0;
        assert_eq!(
            weighted_mass_matrix(&complex, &geo, 1, &weights, Dims::NONE, LIMITS, &mut || {
                seen += 1;
                seen == stop
            }),
            Err(WeightedError::Cancelled)
        );
        assert_eq!(
            weighted_mass_matrix(&complex, &geo, 1, &weights, Dims::NONE, LIMITS, &mut || {
                false
            })
            .unwrap(),
            expected
        );
    }
    let empty = TetComplex::from_tets(0, vec![]);
    let empty_geo = fs_feec::ElementGeometry {
        vol_signed: vec![],
        grads: vec![],
        gram: vec![],
    };
    assert_eq!(
        weighted_mass_matrix(
            &empty,
            &empty_geo,
            0,
            &[],
            Dims::NONE,
            WeightedAssemblyLimits {
                max_cells: 0,
                max_dofs: 0,
                max_triplets: 0
            },
            &mut || false
        )
        .unwrap()
        .matrix()
        .nrows(),
        0
    );
}

#[test]
fn anisotropic_diffusion_manufactured_solution_converges_at_second_order() {
    // G1: -div(K grad u) = 6*pi^2*u, K=diag(1,2,3), zero Dirichlet,
    // u=sin(pi*x)sin(pi*y)sin(pi*z). This solves the real FEEC composition.
    let pi = std::f64::consts::PI;
    let exact = |p: [f64; 3]| p.into_iter().map(|x| (pi * x).sin()).product::<f64>();
    let mut errors = Vec::new();
    for n in [3, 6, 12] {
        let (complex, positions) = kuhn_cube(n);
        let geo = element_geometry(&complex, &positions);
        let weights = vec![
            CellWeight::Tensor([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]);
            complex.tets.len()
        ];
        let weighted = weighted_mass_matrix(
            &complex,
            &geo,
            1,
            &weights,
            Dims([1, 1, -3, -1, 0, 0]),
            LIMITS,
            &mut || false,
        )
        .unwrap();
        let stiffness = stiffness(&incidence_to_csr(&complex.d0()), weighted.matrix());
        let mass = mass_matrix(&complex, &geo, 0);
        let forcing = deram0(&positions, &|p| 6.0 * pi * pi * exact(p));
        let mut rhs = vec![0.0; positions.len()];
        mass.spmv(&forcing, &mut rhs);
        let interior: Vec<_> = (0..positions.len())
            .filter(|&v| !on_unit_cube_boundary(positions[v]))
            .collect();
        let mut slot = vec![usize::MAX; positions.len()];
        for (i, &v) in interior.iter().enumerate() {
            slot[v] = i;
        }
        let mut reduced = fs_sparse::Coo::new(interior.len(), interior.len());
        for (i, &v) in interior.iter().enumerate() {
            let (cols, values) = stiffness.row(v);
            for (&col, &value) in cols.iter().zip(values) {
                if slot[col] != usize::MAX {
                    reduced.push(i, slot[col], value);
                }
            }
        }
        let rhs: Vec<_> = interior.iter().map(|&v| rhs[v]).collect();
        let mut solution = vec![0.0; interior.len()];
        let report = fs_sparse::precond::pcg(
            &reduced.assemble(),
            &rhs,
            &mut solution,
            &fs_sparse::precond::IdentityPrecond,
            1e-12,
            10_000,
        );
        assert!(report.converged, "n={n}: {report:?}");
        let mut error = vec![0.0; positions.len()];
        for (i, &v) in interior.iter().enumerate() {
            error[v] = solution[i] - exact(positions[v]);
        }
        errors.push(quadratic(&mass, &error).sqrt());
    }
    let orders = [
        (errors[0] / errors[1]).log2(),
        (errors[1] / errors[2]).log2(),
    ];
    assert!(
        orders[0] > 1.6 && (orders[1] - 2.0).abs() < 0.25,
        "G1 anisotropic diffusion errors {errors:?}, orders {orders:?}"
    );
}
