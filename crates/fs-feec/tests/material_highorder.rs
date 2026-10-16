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
