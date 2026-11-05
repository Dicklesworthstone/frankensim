//! Spatial density in the original vector moment coordinates (G0/G1/G3).
use std::cell::Cell;

use fs_feec::highorder::vecfam::{Family, VecSpace, VectorMaterialMass};
use fs_feec::{CellWeight, WeightedAssemblyLimits, WeightedError, single_tet, two_tets};
use fs_qty::Dims;

const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const LIMITS: WeightedAssemblyLimits = WeightedAssemblyLimits {
    max_cells: 10,
    max_dofs: 1000,
    max_triplets: 100_000,
};

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8 * b.abs().max(1.0), "{a:e} != {b:e}");
}

fn action(mass: &VectorMaterialMass<'_, '_>, input: &[f64]) -> Vec<f64> {
    let mut output = vec![0.0; mass.ndof()];
    mass.apply_into(input, &mut output, &mut vec![0.0; mass.ndof()], &mut || {
        false
    })
    .unwrap();
    output
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn quadratic_density_retains_inventory_centroid_energy_and_source_snapshot() {
    let (complex, positions) = single_tet();
    let weights = [CellWeight::Scalar(1.0)];
    for family in [Family::Nedelec, Family::Rt] {
        for r in [2, 4] {
            let space = VecSpace::new(&complex, &positions, r, family);
            let calls = Cell::new(0);
            let source_scale = Cell::new(1.0);
            let mass = space
                .prepare_material_mass_with_scalar_field(
                    &positions,
                    &weights,
                    Dims::NONE,
                    DENSITY,
                    r + 3,
                    LIMITS,
                    1_000_000,
                    &mut |cell, l| -> Result<f64, WeightedError> {
                        calls.set(calls.get() + 1);
                        let p: [f64; 3] = core::array::from_fn(|axis| {
                            (0..4)
                                .map(|a| l[a] * positions[complex.tets[cell][a] as usize][axis])
                                .sum()
                        });
                        Ok(source_scale.get() * (1.0 + 4.0 * p[0] * p[0] + 3.0 * p[1]))
                    },
                    &mut || false,
                )
                .unwrap();
            assert_eq!(calls.get(), (r + 3).pow(3));
            source_scale.set(99.0);
            let one = space.interpolate(&positions, &|_| [1.0, 0.0, 0.0]);
            let load = action(&mass, &one);
            // Independent factorial moments: integral x^a y^b z^c =
            // a! b! c! / (a+b+c+3)! on the unit tetrahedron.
            near(dot(&one, &load), 43.0 / 120.0);
            for axis in 0..3 {
                let x = space.interpolate(&positions, &|p| [p[axis], 0.0, 0.0]);
                near(dot(&x, &load), [1.0 / 10.0, 37.0 / 360.0, 7.0 / 90.0][axis]);
            }
            let velocity =
                space.interpolate(&positions, &|[x, y, z]| [1.0 + x, 2.0 - y, -0.5 + 2.0 * z]);
            near(dot(&velocity, &action(&mass, &velocity)), 17299.0 / 10080.0);
            assert_eq!(mass.coefficient_dims(), DENSITY);
            assert_eq!(
                mass.energy_dims(Dims([1, 0, -1, 0, 0, 0])).unwrap(),
                Dims([2, 1, -2, 0, 0, 0])
            );
            let assembled = mass.assemble(&mut || false).unwrap();
            let mut assembled_load = vec![0.0; space.ndof];
            assembled.spmv(&one, &mut assembled_load);
            for (&a, &b) in load.iter().zip(&assembled_load) {
                near(a, b);
            }
            assert_eq!(calls.get(), (r + 3).pow(3));
        }
    }
}

#[test]
fn spatial_multiplier_and_anisotropy_share_csr_action_and_diagonal() {
    let (complex, positions) = two_tets();
    let tensor = [[2.0, 0.25, 0.1], [0.25, 3.0, -0.2], [0.1, -0.2, 4.0]];
    let weights = [
        CellWeight::Tensor(tensor),
        CellWeight::AffineScalar([2.0, 3.0, 7.0, 5.0]),
    ];
    for family in [Family::Nedelec, Family::Rt] {
        let space = VecSpace::new(&complex, &positions, 3, family);
        let mass = space
            .prepare_material_mass_with_scalar_field(
                &positions,
                &weights,
                Dims::NONE,
                DENSITY,
                7,
                LIMITS,
                1_000_000,
                &mut |cell, l| -> Result<f64, WeightedError> {
                    Ok((cell + 1) as f64 * (1.0 + l[0] * l[0] + 2.0 * l[2] * l[3]))
                },
                &mut || false,
            )
            .unwrap();
        let matrix = mass.assemble(&mut || false).unwrap();
        let dense = matrix.to_dense();
        for (i, &d) in mass.diagonal().iter().enumerate() {
            assert_eq!(d.to_bits(), dense[i * space.ndof + i].to_bits());
        }
        let input: Vec<_> = (0..space.ndof)
            .map(|i| (i % 17) as f64 / 11.0 - 0.5)
            .collect();
        let applied = action(&mass, &input);
        let mut expected = vec![0.0; space.ndof];
        matrix.spmv(&input, &mut expected);
        let scale = expected.iter().fold(1.0_f64, |s, x| s.max(x.abs()));
        for (&a, &b) in applied.iter().zip(&expected) {
            assert!((a - b).abs() < 3e-9 * scale);
        }
        assert!(dot(&input, &applied) > 0.0);
        assert_eq!(applied, action(&mass, &input));
    }
}

#[test]
fn nonpolynomial_inverse_density_converges_with_the_declared_quadrature() {
    let (complex, positions) = single_tet();
    let weights = [CellWeight::Scalar(1.0)];
    let a = 20.0_f64;
    // Integrate (1-x)^2 / (2*(1+a*x)) over x in [0,1].
    let exact =
        (a + 1.0).powi(2) * (1.0 + a).ln() / (2.0 * a.powi(3)) - (3.0 * a + 2.0) / (4.0 * a * a);
    for family in [Family::Nedelec, Family::Rt] {
        let space = VecSpace::new(&complex, &positions, 1, family);
        let one = space.interpolate(&positions, &|_| [1.0, 0.0, 0.0]);
        let mut previous = f64::INFINITY;
        for q in [4, 8, 16] {
            let mass = space
                .prepare_material_mass_with_scalar_field(
                    &positions,
                    &weights,
                    Dims::NONE,
                    Dims([3, -1, 0, 0, 0, 0]),
                    q,
                    LIMITS,
                    1_000_000,
                    &mut |cell, l| -> Result<f64, WeightedError> {
                        let x: f64 = (0..4)
                            .map(|i| l[i] * positions[complex.tets[cell][i] as usize][0])
                            .sum();
                        Ok(1.0 / (1.0 + a * x))
                    },
                    &mut || false,
                )
                .unwrap();
            let error = (dot(&one, &action(&mass, &one)) - exact).abs();
            assert!(
                error < previous / 10.0,
                "quadrature {q}: {error:e} after {previous:e}"
            );
            previous = error;
        }
        assert!(previous < 1e-7);
    }
}

#[test]
fn independent_density_scales_do_not_overflow_a_representable_operator() {
    let (complex, positions) = single_tet();
    let positions: Vec<_> = positions.iter().map(|p| p.map(|x| x * 1e-40)).collect();
    let weights = [CellWeight::Scalar(1e160)];
    for family in [Family::Nedelec, Family::Rt] {
        let space = VecSpace::new(&complex, &positions, 1, family);
        let mass = space
            .prepare_material_mass_with_scalar_field(
                &positions,
                &weights,
                Dims::NONE,
                DENSITY,
                4,
                LIMITS,
                1_000_000,
                &mut |_, _| -> Result<f64, WeightedError> { Ok(1e160) },
                &mut || false,
            )
            .unwrap();
        let one = space.interpolate(&positions, &|_| [1.0, 0.0, 0.0]);
        near(dot(&one, &action(&mass, &one)) / 1e200, 1.0 / 6.0);
    }
}

#[derive(Debug, PartialEq)]
enum SourceError {
    Weight(WeightedError),
    Unavailable,
}
impl From<WeightedError> for SourceError {
    fn from(error: WeightedError) -> Self {
        Self::Weight(error)
    }
}

#[test]
fn source_errors_admission_and_cancelled_action_do_not_publish_partial_output() {
    let (complex, positions) = single_tet();
    let weights = [CellWeight::Scalar(1.0)];
    let space = VecSpace::new(&complex, &positions, 2, Family::Rt);
    let calls = Cell::new(0);
    for (q, cap) in [(4, 1_000_000), (17, 1_000_000), (5, 1)] {
        let result = space.prepare_material_mass_with_scalar_field(
            &positions,
            &weights,
            Dims::NONE,
            DENSITY,
            q,
            LIMITS,
            cap,
            &mut |_, _| -> Result<f64, WeightedError> {
                calls.set(calls.get() + 1);
                Ok(1.0)
            },
            &mut || false,
        );
        assert!(result.is_err());
        assert_eq!(calls.get(), 0);
    }
    let result = space.prepare_material_mass_with_scalar_field(
        &positions,
        &weights,
        Dims::NONE,
        DENSITY,
        5,
        LIMITS,
        1_000_000,
        &mut |_, _| -> Result<f64, SourceError> {
            calls.set(calls.get() + 1);
            if calls.get() == 3 {
                Err(SourceError::Unavailable)
            } else {
                Ok(1.0)
            }
        },
        &mut || false,
    );
    assert!(matches!(result, Err(SourceError::Unavailable)));
    assert_eq!(calls.get(), 3);
    for invalid in [0.0, -1.0, f64::NAN] {
        let result = space.prepare_material_mass_with_scalar_field(
            &positions,
            &weights,
            Dims::NONE,
            DENSITY,
            5,
            LIMITS,
            1_000_000,
            &mut |_, _| -> Result<f64, WeightedError> { Ok(invalid) },
            &mut || false,
        );
        assert!(matches!(
            result,
            Err(WeightedError::InvalidCoefficient { .. })
        ));
    }
    let no_csr = WeightedAssemblyLimits {
        max_triplets: 0,
        ..LIMITS
    };
    let mass = space
        .prepare_material_mass_with_scalar_field(
            &positions,
            &weights,
            Dims::NONE,
            DENSITY,
            5,
            no_csr,
            1_000_000,
            &mut |_, l| -> Result<f64, WeightedError> { Ok(1.0 + l[0] * l[0]) },
            &mut || false,
        )
        .unwrap();
    assert!(mass.assemble(&mut || false).is_err());
    let input = vec![1.0; space.ndof];
    let mut output = vec![1234.0; space.ndof];
    let mut scratch = vec![0.0; space.ndof];
    let mut polls = 0;
    assert!(
        mass.apply_into(&input, &mut output, &mut scratch, &mut || {
            polls += 1;
            polls == 10
        })
        .is_err()
    );
    assert_eq!(output, vec![1234.0; space.ndof]);
    mass.apply_into(&input, &mut output, &mut scratch, &mut || false)
        .unwrap();
    assert!(output.iter().all(|x| x.is_finite()));
}
