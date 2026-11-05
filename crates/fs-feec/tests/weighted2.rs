//! G0/G1/G3/G4 regressions for weighted two-dimensional FEEC assembly.

use fs_feec::{
    CellWeight2, WeightedAssemblyLimits, WeightedError, WeightedMass, incidence_to_csr, stiffness,
    weighted_mass_matrix_2d,
};
use fs_qty::Dims;
use fs_rep_mesh::{Metric2, TriComplex2};

const LIMITS: WeightedAssemblyLimits = WeightedAssemblyLimits {
    max_cells: 1000,
    max_dofs: 3000,
    max_triplets: 9000,
};

fn triangle(metric: Metric2) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "weighted-triangle",
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        vec![[0, 1, 2]],
        metric,
    )
    .unwrap()
}

fn assemble(c: &TriComplex2, k: u8, w: &[CellWeight2]) -> WeightedMass {
    weighted_mass_matrix_2d(c, k, w, Dims::NONE, LIMITS, &mut || false).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn energy(matrix: &fs_sparse::Csr, values: &[f64]) -> f64 {
    values
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            values
                .iter()
                .enumerate()
                .map(|(j, &y)| x * matrix.get(i, j) * y)
                .sum::<f64>()
        })
        .sum()
}

#[test]
fn g1_graded_density_integrates_planar_and_axisymmetric_fields() {
    for axisymmetric in [false, true] {
        let offset = if axisymmetric { 1.0 } else { 0.0 };
        let metric = if axisymmetric {
            Metric2::axisymmetric(2.0).unwrap()
        } else {
            Metric2::planar(2.0).unwrap()
        };
        let c = TriComplex2::from_indexed_triangles(
            "graded-density",
            vec![[offset, 0.0], [offset + 1.0, 0.0], [offset, 1.0]],
            vec![[0, 1, 2]],
            metric,
        )
        .unwrap();
        let density = [CellWeight2::AffineScalar([2.0, 5.0, 7.0])];
        // Independent Cartesian polynomial integrals on x>=0,y>=0,x+y<=1:
        // rho=2+3x+5y; scalar=1+2x+3y; vector=(1-y,2+r).
        // The physical volume measure is 2 dA or 2(1+x) dA.
        let mass = assemble(&c, 0, &density);
        let expected_mass = if axisymmetric { 25.0 / 4.0 } else { 14.0 / 3.0 };
        close(energy(mass.matrix(), &[1.0; 3]), expected_mass);
        let expected_scalar = if axisymmetric {
            154.0 / 3.0
        } else {
            383.0 / 10.0
        };
        close(energy(mass.matrix(), &[1.0, 3.0, 4.0]), expected_scalar);
        let velocity: Vec<_> = c
            .edges()
            .iter()
            .map(|&[a, b]| {
                let a = c.vertices()[a as usize];
                let b = c.vertices()[b as usize];
                (1.0 - f64::midpoint(a[1], b[1])) * (b[0] - a[0])
                    + (2.0 + f64::midpoint(a[0], b[0])) * (b[1] - a[1])
            })
            .collect();
        let expected_vector = if axisymmetric {
            449.0 / 6.0
        } else {
            279.0 / 10.0
        };
        close(
            energy(assemble(&c, 1, &density).matrix(), &velocity),
            expected_vector,
        );
        // A unit top-form field has the original embedded area as its dof.
        close(
            energy(assemble(&c, 2, &density).matrix(), &[0.5]),
            expected_mass,
        );
        for k in 0..=2 {
            assert_eq!(
                assemble(&c, k, &[CellWeight2::AffineScalar([5.0; 3])]).matrix(),
                assemble(&c, k, &[CellWeight2::Scalar(5.0)]).matrix(),
            );
        }
    }
}

#[test]
fn g3_graded_density_face_order_and_admission() {
    let vertices = vec![[1.0, 0.0], [2.0, 0.0], [1.0, 1.0]];
    let global_density = [2.0, 5.0, 7.0];
    let mut reference: [Option<Vec<f64>>; 3] = [None, None, None];
    for face in [[0, 1, 2], [2, 1, 0], [1, 2, 0]] {
        let c = TriComplex2::from_indexed_triangles(
            "graded-orientation",
            vertices.clone(),
            vec![face],
            Metric2::axisymmetric(2.0).unwrap(),
        )
        .unwrap();
        let density = c.faces()[0].map(|v| global_density[v as usize]);
        for degree in 0..=2 {
            let matrix = assemble(&c, degree, &[CellWeight2::AffineScalar(density)])
                .matrix()
                .to_dense();
            if let Some(expected) = &reference[degree as usize] {
                for (&a, &b) in matrix.iter().zip(expected) {
                    close(a, b);
                }
            } else {
                reference[degree as usize] = Some(matrix);
            }
            for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                assert!(matches!(
                    weighted_mass_matrix_2d(
                        &c,
                        degree,
                        &[CellWeight2::AffineScalar([1.0, bad, 2.0])],
                        Dims::NONE,
                        LIMITS,
                        &mut || false,
                    ),
                    Err(WeightedError::InvalidCoefficient { cell: 0, .. })
                ));
            }
        }
    }
}

#[test]
fn g0_planar_normalization_and_physical_units() {
    let c = triangle(Metric2::planar(2.0).unwrap());
    let dims = Dims([-1, 1, -2, -1, 0, 0]); // J/(m^3 K).
    let m = weighted_mass_matrix_2d(
        &c,
        0,
        &[CellWeight2::Scalar(3.0)],
        dims,
        LIMITS,
        &mut || false,
    )
    .unwrap();
    assert_eq!(m.coefficient_dims(), dims);
    assert_eq!(m.entry_dims(), Dims([2, 1, -2, -1, 0, 0])); // J/K.
    assert_eq!(m.degree(), 0);
    for i in 0..3 {
        for j in 0..3 {
            close(m.matrix().get(i, j), if i == j { 0.5 } else { 0.25 });
        }
    }
    close(energy(m.matrix(), &[1.0; 3]), 3.0);
    let top = assemble(&c, 2, &[CellWeight2::Scalar(3.0)]);
    close(top.matrix().get(0, 0), 12.0);
    assert_eq!(top.entry_dims(), Dims([-1, 0, 0, 0, 0, 0]));
}

#[test]
fn g1_anisotropic_stiffness_matches_affine_energy() {
    let c = triangle(Metric2::planar(2.0).unwrap());
    let m = assemble(&c, 1, &[CellWeight2::Tensor([[3.0, 0.5], [0.5, 2.0]])]);
    let k = stiffness(&incidence_to_csr(&c.d0()), m.matrix());
    // T=7+2x+3y: grad(T)^T W grad(T)=36, volume=1.
    close(energy(&k, &[7.0, 9.0, 10.0]), 36.0);
    close(energy(&k, &[1.0; 3]), 0.0);
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(
                m.matrix().get(i, j).to_bits(),
                m.matrix().get(j, i).to_bits()
            );
        }
    }
    assert_eq!(c.d1().apply(&c.d0().apply(&[3, -2, 8])), vec![0]);
}

#[test]
fn g1_axisymmetric_mass_uses_exact_radial_moments() {
    let c = TriComplex2::from_indexed_triangles(
        "weighted-meridian",
        vec![[1.0, 0.0], [2.0, 0.0], [1.0, 1.0]],
        vec![[0, 1, 2]],
        Metric2::axisymmetric(2.0).unwrap(),
    )
    .unwrap();
    let m = assemble(&c, 0, &[CellWeight2::Scalar(1.0)]);
    // r=1+lambda_1; exact cubic simplex moments, angle=2, area=1/2.
    let exact = [
        [1.0 / 5.0, 7.0 / 60.0, 1.0 / 10.0],
        [7.0 / 60.0, 4.0 / 15.0, 7.0 / 60.0],
        [1.0 / 10.0, 7.0 / 60.0, 1.0 / 5.0],
    ];
    for (i, row) in exact.iter().enumerate() {
        for (j, &x) in row.iter().enumerate() {
            close(m.matrix().get(i, j), x);
        }
    }
    assert!(m.matrix().get(1, 1) > m.matrix().get(0, 0)); // Fails for centroid-scaled mass.
    close(energy(m.matrix(), &[1.0; 3]), 4.0 / 3.0);
    let edges = assemble(&c, 1, &[CellWeight2::Tensor([[3.0, 0.5], [0.5, 2.0]])]);
    let k = stiffness(&incidence_to_csr(&c.d0()), edges.matrix());
    close(energy(&k, &[2.0, 4.0, 5.0]), 48.0);
    close(
        assemble(&c, 2, &[CellWeight2::Scalar(1.0)])
            .matrix()
            .get(0, 0),
        16.0 / 3.0,
    );
}

#[test]
fn g3_coefficient_jumps_share_dofs_and_preserve_topology() {
    let c = TriComplex2::from_indexed_triangles(
        "weighted-interface",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        vec![[0, 1, 2], [0, 2, 3]],
        Metric2::planar(1.0).unwrap(),
    )
    .unwrap();
    let w = [CellWeight2::Scalar(2.0), CellWeight2::Scalar(10.0)];
    let m = assemble(&c, 1, &w);
    assert_eq!(m.matrix().nrows(), 5);
    let k = stiffness(&incidence_to_csr(&c.d0()), m.matrix());
    close(energy(&k, &[0.0, 1.0, 2.0, 1.0]), 12.0);
    close(energy(assemble(&c, 0, &w).matrix(), &[1.0; 4]), 6.0);
    assert_eq!(c.d1().apply(&c.d0().apply(&[4, 9, -1, 3])), vec![0, 0]);
}

#[test]
fn g3_orientation_and_reindexing_preserve_energy() {
    for face in [[0, 1, 2], [2, 1, 0], [1, 2, 0]] {
        let c = TriComplex2::from_indexed_triangles(
            "weighted-orientation",
            vec![[0.0, 1.0], [0.0, 0.0], [1.0, 0.0]],
            vec![face],
            Metric2::planar(2.0).unwrap(),
        )
        .unwrap();
        let m = assemble(&c, 1, &[CellWeight2::Tensor([[3.0, 0.5], [0.5, 2.0]])]);
        let k = stiffness(&incidence_to_csr(&c.d0()), m.matrix());
        close(energy(&k, &[10.0, 7.0, 9.0]), 36.0);
    }
}

#[test]
fn g0_scaled_tensor_avoids_overflow_in_representable_entries() {
    let c = triangle(Metric2::planar(2.0).unwrap());
    let unit = assemble(&c, 1, &[CellWeight2::Scalar(1.0)]);
    let scaled = assemble(
        &c,
        1,
        &[CellWeight2::Tensor([[f64::MAX, 0.0], [0.0, f64::MAX]])],
    );
    for i in 0..3 {
        for j in 0..3 {
            let value = scaled.matrix().get(i, j);
            assert!(value.is_finite());
            close(value / f64::MAX, unit.matrix().get(i, j));
        }
    }
}

#[test]
fn g0_invalid_inputs_are_typed_refusals() {
    let c = triangle(Metric2::planar(1.0).unwrap());
    for w in [
        CellWeight2::Scalar(0.0),
        CellWeight2::Scalar(-1.0),
        CellWeight2::Scalar(f64::NAN),
        CellWeight2::Scalar(f64::INFINITY),
        CellWeight2::Tensor([[1.0, 0.0], [0.1, 1.0]]),
        CellWeight2::Tensor([[1.0, 2.0], [2.0, 1.0]]),
        CellWeight2::Tensor([[1.0, 1.0], [1.0, 1.0]]),
        CellWeight2::Tensor([[1.0, 0.0], [0.0, f64::INFINITY]]),
    ] {
        assert!(matches!(
            weighted_mass_matrix_2d(&c, 1, &[w], Dims::NONE, LIMITS, &mut || false),
            Err(WeightedError::InvalidCoefficient { cell: 0, .. })
        ));
    }
    for k in [0, 2] {
        assert!(matches!(
            weighted_mass_matrix_2d(
                &c,
                k,
                &[CellWeight2::Tensor([[1.0, 0.0], [0.0, 1.0]])],
                Dims::NONE,
                LIMITS,
                &mut || false
            ),
            Err(WeightedError::InvalidCoefficient { .. })
        ));
    }
    assert!(matches!(
        weighted_mass_matrix_2d(
            &c,
            3,
            &[CellWeight2::Scalar(1.0)],
            Dims::NONE,
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::InvalidDegree(3))
    ));
    assert!(matches!(
        weighted_mass_matrix_2d(&c, 0, &[], Dims::NONE, LIMITS, &mut || false),
        Err(WeightedError::CoefficientCount {
            expected: 1,
            actual: 0
        })
    ));
    assert!(matches!(
        weighted_mass_matrix_2d(
            &c,
            0,
            &[CellWeight2::Scalar(1.0)],
            Dims([127, 0, 0, 0, 0, 0]),
            LIMITS,
            &mut || false
        ),
        Err(WeightedError::DimensionOverflow)
    ));
}

#[test]
fn g0_exact_limits_and_one_beyond() {
    let c = triangle(Metric2::planar(1.0).unwrap());
    let limits = WeightedAssemblyLimits {
        max_cells: 1,
        max_dofs: 3,
        max_triplets: 9,
    };
    weighted_mass_matrix_2d(
        &c,
        1,
        &[CellWeight2::Scalar(1.0)],
        Dims::NONE,
        limits,
        &mut || false,
    )
    .unwrap();
    for refused in [
        WeightedAssemblyLimits {
            max_cells: 0,
            ..limits
        },
        WeightedAssemblyLimits {
            max_dofs: 2,
            ..limits
        },
        WeightedAssemblyLimits {
            max_triplets: 8,
            ..limits
        },
    ] {
        assert!(matches!(
            weighted_mass_matrix_2d(
                &c,
                1,
                &[CellWeight2::Scalar(1.0)],
                Dims::NONE,
                refused,
                &mut || false
            ),
            Err(WeightedError::LimitExceeded { .. })
        ));
    }
}

#[test]
fn g4_every_checkpoint_cancels_atomically_and_retry_matches() {
    let c = triangle(Metric2::planar(1.0).unwrap());
    let w = [CellWeight2::Scalar(2.0)];
    let mut checkpoints = 0;
    let expected = weighted_mass_matrix_2d(&c, 1, &w, Dims::NONE, LIMITS, &mut || {
        checkpoints += 1;
        false
    })
    .unwrap();
    for stop in 1..=checkpoints {
        let mut polls = 0;
        let result = weighted_mass_matrix_2d(&c, 1, &w, Dims::NONE, LIMITS, &mut || {
            polls += 1;
            polls == stop
        });
        assert!(
            matches!(result, Err(WeightedError::Cancelled)),
            "checkpoint {stop}"
        );
    }
    assert_eq!(assemble(&c, 1, &w).matrix(), expected.matrix());
}
