//! G0/G1/G3/G4: material-weighted high-order H1 operators and actual solves.
use fs_feec::highorder::material::{H1Pairing, MAX_MATERIAL_ORDER, MaterialOperator};
use fs_feec::highorder::simplex::SimplexSpace;
use fs_feec::{
    CellWeight, WeightedAssemblyLimits, WeightedError, element_geometry, kuhn_cube, single_tet,
    two_tets,
};
use fs_qty::Dims;
use fs_sparse::{
    Csr,
    precond::{IdentityPrecond, pcg},
};

const LIMITS: WeightedAssemblyLimits = WeightedAssemblyLimits {
    max_cells: 10_000,
    max_dofs: 100_000,
    max_triplets: 2_000_000,
};
const TENSOR: [[f64; 3]; 3] = [[2.0, 0.3, -0.2], [0.3, 3.0, 0.4], [-0.2, 0.4, 4.0]];
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-11 * b.abs().max(1.0), "{a:e} != {b:e}");
}
fn energy(matrix: &Csr, values: &[f64]) -> f64 {
    let mut product = vec![0.0; values.len()];
    matrix.spmv(values, &mut product);
    values.iter().zip(product).map(|(u, v)| u * v).sum()
}
fn tensor_energy(tensor: [[f64; 3]; 3], v: [f64; 3]) -> f64 {
    (0..3)
        .map(|i| (0..3).map(|j| v[i] * tensor[i][j] * v[j]).sum::<f64>())
        .sum()
}

#[test]
fn prepared_material_action_matches_consistent_graded_and_tensor_operators() {
    let (complex, positions) = two_tets();
    let geo = element_geometry(&complex, &positions);
    for r in [1, 3, 6] {
        let space = SimplexSpace::new(&complex, r);
        for (pairing, weights) in [
            (
                H1Pairing::Mass,
                [
                    CellWeight::AffineScalar([2.0, 5.0, 7.0, 11.0]),
                    CellWeight::Scalar(13.0),
                ],
            ),
            (
                H1Pairing::Diffusion,
                [
                    CellWeight::AffineScalar([2.0, 5.0, 7.0, 11.0]),
                    CellWeight::Tensor(TENSOR),
                ],
            ),
        ] {
            let op = MaterialOperator::new(
                &space,
                &geo,
                pairing,
                &weights,
                Dims::NONE,
                LIMITS,
                &mut || false,
            )
            .unwrap();
            let assembled = op.assemble(&mut || false).unwrap();
            let action = op.prepare_action(1_000_000, &mut || false).unwrap();
            assert_eq!(action.ndof(), space.ndof);
            let dense = assembled.matrix().to_dense();
            for (i, &d) in action.diagonal().iter().enumerate() {
                assert_eq!(d.to_bits(), dense[i * space.ndof + i].to_bits());
            }
            let input: Vec<_> = (0..space.ndof)
                .map(|i| (i % 13) as f64 / 7.0 - 0.8)
                .collect();
            let mut expected = vec![0.0; space.ndof];
            assembled.matrix().spmv(&input, &mut expected);
            let mut actual = vec![0.0; space.ndof];
            let mut scratch = vec![0.0; space.ndof];
            action
                .apply_into(&input, &mut actual, &mut scratch, &mut || false)
                .unwrap();
            for (&a, &b) in actual.iter().zip(&expected) {
                near(a, b);
            }
            // Reuse the same prepared plan and scratch; prior scratch must not
            // contribute to the next physical load or residual.
            let first = actual.clone();
            action
                .apply_into(&input, &mut actual, &mut scratch, &mut || false)
                .unwrap();
            assert_eq!(first, actual);
        }
    }
}

#[test]
fn matrix_free_density_load_retains_mass_centroid_without_triplet_storage() {
    let (complex, positions) = single_tet();
    let geo = element_geometry(&complex, &positions);
    let space = SimplexSpace::new(&complex, 6);
    // rho = 3 + 12*x on the unit tetrahedron. Exact inventory is 1 kg;
    // its centroid is (13/40, 9/40, 9/40), independently of basis order.
    let weights = [CellWeight::AffineScalar(
        complex.tets[0].map(|v| 3.0 + 12.0 * positions[v as usize][0]),
    )];
    let op = MaterialOperator::new(
        &space,
        &geo,
        H1Pairing::Mass,
        &weights,
        Dims([-3, 1, 0, 0, 0, 0]),
        WeightedAssemblyLimits {
            max_triplets: 0,
            ..LIMITS
        },
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        op.assemble(&mut || false),
        Err(WeightedError::LimitExceeded {
            resource: "triplets",
            ..
        })
    ));
    let action = op.prepare_action(1_000_000, &mut || false).unwrap();
    let mut constant = vec![0.0; space.ndof];
    assert_eq!(action.entry_dims(), Dims([0, 1, 0, 0, 0, 0]));
    assert_eq!(
        action.output_dims(Dims([1, 0, -2, 0, 0, 0])).unwrap(),
        Dims([1, 1, -2, 0, 0, 0])
    );
    constant[..positions.len()].fill(1.0);
    let mut load = vec![0.0; space.ndof];
    let mut scratch = vec![0.0; space.ndof];
    action
        .apply_into(&constant, &mut load, &mut scratch, &mut || false)
        .unwrap();
    near(load[..positions.len()].iter().sum(), 1.0);
    for c in 0..3 {
        let moment: f64 = positions
            .iter()
            .enumerate()
            .map(|(i, p)| p[c] * load[i])
            .sum();
        near(moment, [13.0 / 40.0, 9.0 / 40.0, 9.0 / 40.0][c]);
    }
}

#[test]
fn material_action_refusal_is_atomic_and_storage_is_bounded() {
    let (complex, positions) = single_tet();
    let geo = element_geometry(&complex, &positions);
    let space = SimplexSpace::new(&complex, 3);
    let weights = [CellWeight::AffineScalar([2.0, 5.0, 7.0, 11.0])];
    let op = MaterialOperator::new(
        &space,
        &geo,
        H1Pairing::Mass,
        &weights,
        Dims::NONE,
        LIMITS,
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        op.prepare_action(0, &mut || false),
        Err(WeightedError::LimitExceeded {
            resource: "prepared material action bytes",
            ..
        })
    ));
    let action = op.prepare_action(1_000_000, &mut || false).unwrap();
    let mut input = vec![1.0; space.ndof];
    let mut output = vec![23.0; space.ndof];
    let mut scratch = vec![0.0; space.ndof];
    let mut polls = 0;
    assert!(matches!(
        action.apply_into(&input, &mut output, &mut scratch, &mut || {
            polls += 1;
            polls > 30
        }),
        Err(WeightedError::Cancelled)
    ));
    assert_eq!(output, vec![23.0; space.ndof]);
    input[space.ndof / 2] = f64::NAN;
    assert!(matches!(
        action.apply_into(&input, &mut output, &mut scratch, &mut || false),
        Err(WeightedError::NonFiniteEntry)
    ));
    assert_eq!(output, vec![23.0; space.ndof]);
    assert!(
        action
            .apply_into(
                &input[..space.ndof - 1],
                &mut output,
                &mut scratch,
                &mut || false
            )
            .is_err()
    );
    assert_eq!(output, vec![23.0; space.ndof]);
}

#[test]
fn affine_density_highorder_mass_and_diffusion_match_exact_whitney_restriction() {
    use fs_feec::{incidence_to_csr, stiffness, weighted_mass_matrix};
    let (complex, positions) = two_tets();
    let geo = element_geometry(&complex, &positions);
    let weights = [
        CellWeight::AffineScalar([2.0, 5.0, 7.0, 11.0]),
        CellWeight::AffineScalar([13.0, 3.0, 17.0, 2.0]),
    ];
    for pairing in [H1Pairing::Mass, H1Pairing::Diffusion] {
        let degree = if pairing == H1Pairing::Mass { 0 } else { 1 };
        let whitney = weighted_mass_matrix(
            &complex,
            &geo,
            degree,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap();
        let reference = if pairing == H1Pairing::Mass {
            whitney.matrix().to_dense()
        } else {
            stiffness(&incidence_to_csr(&complex.d0()), whitney.matrix()).to_dense()
        };
        for r in [1, 2, 4, 6] {
            let space = SimplexSpace::new(&complex, r);
            let matrix = MaterialOperator::new(
                &space,
                &geo,
                pairing,
                &weights,
                Dims::NONE,
                LIMITS,
                &mut || false,
            )
            .unwrap()
            .assemble(&mut || false)
            .unwrap();
            let actual = matrix.matrix().to_dense();
            for i in 0..positions.len() {
                for j in 0..positions.len() {
                    near(
                        actual[i * space.ndof + j],
                        reference[i * positions.len() + j],
                    );
                }
            }
            let trial: Vec<_> = (0..space.ndof)
                .map(|i| 0.1 + (i % 7) as f64 / 7.0)
                .collect();
            assert!(energy(matrix.matrix(), &trial) > 0.0);
        }
    }
}

#[test]
fn affine_highorder_weight_samples_the_original_higher_modes() {
    use fs_feec::highorder::simplex::duffy_quadrature;
    let (complex, positions) = single_tet();
    let geo = element_geometry(&complex, &positions);
    let density = [2.0, 5.0, 7.0, 11.0];
    for r in [2, 4, 6] {
        let space = SimplexSpace::new(&complex, r);
        let weights = [CellWeight::AffineScalar(density)];
        let matrix = MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Mass,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap()
        .assemble(&mut || false)
        .unwrap();
        let trial: Vec<_> = (0..space.ndof)
            .map(|i| ((3 * i + 1) % 11) as f64 / 11.0)
            .collect();
        let dofs = space.element_dofs(0);
        let expected: f64 = duffy_quadrature(r + 4)
            .into_iter()
            .map(|(lambda, w)| {
                let rho: f64 = (0..4).map(|a| density[a] * lambda[a]).sum();
                let value: f64 = dofs
                    .iter()
                    .map(|(i, basis)| trial[*i] * basis.eval(lambda, r))
                    .sum();
                w * 6.0 * geo.vol_signed[0].abs() * rho * value * value
            })
            .sum();
        near(energy(matrix.matrix(), &trial), expected);
        let constant = [CellWeight::AffineScalar([5.0; 4])];
        let old = [CellWeight::Scalar(5.0)];
        let assemble = |weights| {
            MaterialOperator::new(
                &space,
                &geo,
                H1Pairing::Mass,
                weights,
                Dims::NONE,
                LIMITS,
                &mut || false,
            )
            .unwrap()
            .assemble(&mut || false)
            .unwrap()
        };
        assert_eq!(assemble(&constant), assemble(&old));
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                MaterialOperator::new(
                    &space,
                    &geo,
                    H1Pairing::Mass,
                    &[CellWeight::AffineScalar([1.0, bad, 2.0, 3.0])],
                    Dims::NONE,
                    LIMITS,
                    &mut || false,
                ),
                Err(WeightedError::InvalidCoefficient { cell: 0, .. })
            ));
        }
    }
}

#[test]
fn unit_materials_reproduce_existing_highorder_operators() {
    let (complex, positions) = two_tets();
    let geo = element_geometry(&complex, &positions);
    for r in 1..=MAX_MATERIAL_ORDER {
        let space = SimplexSpace::new(&complex, r);
        for pairing in [H1Pairing::Mass, H1Pairing::Diffusion] {
            let op = MaterialOperator::new(
                &space,
                &geo,
                pairing,
                &[CellWeight::Scalar(1.0); 2],
                Dims::NONE,
                LIMITS,
                &mut || false,
            )
            .unwrap();
            let matrix = op.assemble(&mut || false).unwrap();
            let reference = match pairing {
                H1Pairing::Mass => space.mass(&positions),
                H1Pairing::Diffusion => space.stiffness(&positions),
            }
            .to_dense();
            let actual = matrix.matrix().to_dense();
            for (&a, &b) in actual.iter().zip(&reference) {
                near(a, b);
            }
            for i in 0..space.ndof {
                for j in 0..space.ndof {
                    assert_eq!(
                        actual[i * space.ndof + j].to_bits(),
                        actual[j * space.ndof + i].to_bits()
                    );
                }
            }
            assert_eq!(matrix.pairing(), pairing);
        }
    }
}

#[test]
fn discontinuous_anisotropy_has_correct_linear_energy_and_units() {
    let (complex, positions) = two_tets();
    let geo = element_geometry(&complex, &positions);
    let coefficient_dims = Dims([1, 1, -3, -1, 0, 0]);
    let weights = [CellWeight::Tensor(TENSOR), CellWeight::Scalar(7.0)];
    let gradient = [0.7, -1.3, 2.1];
    let expected = geo.vol_signed[0].abs() * tensor_energy(TENSOR, gradient)
        + geo.vol_signed[1].abs() * 7.0 * gradient.iter().map(|v| v * v).sum::<f64>();
    for r in [1, 2, 4, 6] {
        let space = SimplexSpace::new(&complex, r);
        let op = MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Diffusion,
            &weights,
            coefficient_dims,
            LIMITS,
            &mut || false,
        )
        .unwrap();
        let matrix = op.assemble(&mut || false).unwrap();
        let mut values = vec![0.0; space.ndof];
        for (i, p) in positions.iter().enumerate() {
            values[i] = 2.0 + (0..3).map(|j| p[j] * gradient[j]).sum::<f64>();
        }
        near(energy(matrix.matrix(), &values), expected);
        assert_eq!(matrix.entry_dims(), Dims([2, 1, -3, -1, 0, 0]));
        assert_eq!(matrix.coefficient_dims(), coefficient_dims);
        assert_eq!(
            matrix.output_dims(Dims([0, 0, 0, 1, 0, 0])).unwrap(),
            Dims([2, 1, -3, 0, 0, 0])
        );
        values.fill(0.0);
        values[..positions.len()].fill(1.0);
        near(energy(matrix.matrix(), &values), 0.0);
    }
}

#[test]
fn storage_mass_preserves_jumps_constant_fields_and_scaling() {
    let (complex, positions) = two_tets();
    let weights = [CellWeight::Scalar(2.0), CellWeight::Scalar(11.0)];
    let space = SimplexSpace::new(&complex, 3);
    let mut constant = vec![0.0; space.ndof];
    constant[..positions.len()].fill(1.0);
    let mut energies = Vec::new();
    for scale in [1.0, 3.0] {
        let scaled: Vec<_> = positions.iter().map(|p| p.map(|v| scale * v)).collect();
        let geo = element_geometry(&complex, &scaled);
        let op = MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Mass,
            &weights,
            Dims([-1, 1, -2, -1, 0, 0]),
            LIMITS,
            &mut || false,
        )
        .unwrap();
        let mass = op.assemble(&mut || false).unwrap();
        let actual = energy(mass.matrix(), &constant);
        near(
            actual,
            2.0 * geo.vol_signed[0].abs() + 11.0 * geo.vol_signed[1].abs(),
        );
        assert_eq!(mass.entry_dims(), Dims([2, 1, -2, -1, 0, 0]));
        energies.push(actual);
    }
    near(energies[1], 27.0 * energies[0]);
}

#[test]
fn malformed_materials_layouts_and_exact_budget_edges_refuse() {
    let (complex, positions) = single_tet();
    let geo = element_geometry(&complex, &positions);
    let mut space = SimplexSpace::new(&complex, 2);
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            MaterialOperator::new(
                &space,
                &geo,
                H1Pairing::Mass,
                &[CellWeight::Scalar(value)],
                Dims::NONE,
                LIMITS,
                &mut || false
            ),
            Err(WeightedError::InvalidCoefficient { .. })
        ));
    }
    assert!(matches!(
        MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Mass,
            &[CellWeight::Tensor(TENSOR)],
            Dims::NONE,
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::InvalidCoefficient { .. })
    ));
    let bad = [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(matches!(
        MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Diffusion,
            &[CellWeight::Tensor(bad)],
            Dims::NONE,
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::InvalidCoefficient { .. })
    ));
    assert!(matches!(
        MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Mass,
            &[],
            Dims::NONE,
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::CoefficientCount { .. })
    ));
    assert!(matches!(
        MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Mass,
            &[CellWeight::Scalar(1.0)],
            Dims([127, 0, 0, 0, 0, 0]),
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::DimensionOverflow)
    ));
    let weights = [CellWeight::Scalar(1.0)];
    for cap in [99, 100] {
        let limits = WeightedAssemblyLimits {
            max_cells: 1,
            max_dofs: space.ndof,
            max_triplets: cap,
        };
        let op = MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Mass,
            &weights,
            Dims::NONE,
            limits,
            &mut || false,
        )
        .unwrap();
        assert_eq!(op.assemble(&mut || false).is_ok(), cap == 100);
    }
    space.ndof += 1;
    assert!(matches!(
        MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Mass,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::InvalidGeometry { .. })
    ));
    space.r = usize::MAX;
    assert!(matches!(
        MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Mass,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::InvalidGeometry { .. })
    ));
}

#[test]
fn cancellation_during_quadrature_and_publication_has_deterministic_retry() {
    let (complex, positions) = single_tet();
    let geo = element_geometry(&complex, &positions);
    let space = SimplexSpace::new(&complex, 2);
    let weights = [CellWeight::Tensor(TENSOR)];
    let op = MaterialOperator::new(
        &space,
        &geo,
        H1Pairing::Diffusion,
        &weights,
        Dims::NONE,
        LIMITS,
        &mut || false,
    )
    .unwrap();
    let mut count = 0;
    let expected = op
        .assemble(&mut || {
            count += 1;
            false
        })
        .unwrap();
    for stop in [1, 10, count / 2, count - 1, count] {
        let mut seen = 0;
        assert_eq!(
            op.assemble(&mut || {
                seen += 1;
                seen == stop
            }),
            Err(WeightedError::Cancelled)
        );
        assert_eq!(op.assemble(&mut || false).unwrap(), expected);
    }
}

#[test]
fn anisotropic_manufactured_diffusion_has_cubic_l2_convergence() {
    let pi = std::f64::consts::PI;
    let solution = |p: [f64; 3]| (pi * p[0]).sin() * (pi * p[1]).sin() * (pi * p[2]).sin();
    let forcing = |p: [f64; 3]| {
        let s = p.map(|x| (pi * x).sin());
        let c = p.map(|x| (pi * x).cos());
        pi * pi
            * ((TENSOR[0][0] + TENSOR[1][1] + TENSOR[2][2]) * s[0] * s[1] * s[2]
                - 2.0
                    * (TENSOR[0][1] * c[0] * c[1] * s[2]
                        + TENSOR[0][2] * c[0] * s[1] * c[2]
                        + TENSOR[1][2] * s[0] * c[1] * c[2]))
    };
    let mut errors = Vec::new();
    for m in [2, 4] {
        let (complex, positions) = kuhn_cube(m);
        let geo = element_geometry(&complex, &positions);
        let space = SimplexSpace::new(&complex, 2);
        let weights = vec![CellWeight::Tensor(TENSOR); complex.tets.len()];
        let op = MaterialOperator::new(
            &space,
            &geo,
            H1Pairing::Diffusion,
            &weights,
            Dims::NONE,
            LIMITS,
            &mut || false,
        )
        .unwrap();
        let stiffness = op.assemble(&mut || false).unwrap();
        let load = space.load(&positions, &forcing);
        let boundary = space.boundary_mask();
        let interior: Vec<_> = (0..space.ndof).filter(|&i| !boundary[i]).collect();
        let mut slot = vec![usize::MAX; space.ndof];
        for (i, &d) in interior.iter().enumerate() {
            slot[d] = i;
        }
        let mut reduced = fs_sparse::Coo::new(interior.len(), interior.len());
        for (i, &d) in interior.iter().enumerate() {
            let (cols, vals) = stiffness.matrix().row(d);
            for (&j, &v) in cols.iter().zip(vals) {
                if slot[j] != usize::MAX {
                    reduced.push(i, slot[j], v);
                }
            }
        }
        let rhs: Vec<_> = interior.iter().map(|&d| load[d]).collect();
        let mut values = vec![0.0; interior.len()];
        let report = pcg(
            &reduced.assemble(),
            &rhs,
            &mut values,
            &IdentityPrecond,
            1e-12,
            20_000,
        );
        assert!(report.converged, "{report:?}");
        let mut full = vec![0.0; space.ndof];
        for (i, &d) in interior.iter().enumerate() {
            full[d] = values[i];
        }
        errors.push(space.l2_error(&positions, &full, &solution));
    }
    let order = (errors[0] / errors[1]).log2();
    println!("anisotropic P2 diffusion: L2={errors:?}, order={order}");
    assert!(
        order > 2.7,
        "expected cubic convergence, got {order}: {errors:?}"
    );
}
