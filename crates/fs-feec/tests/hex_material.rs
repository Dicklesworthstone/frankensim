//! Consistent density on the native sum-factorized Q_r layout.
use fs_feec::WeightedError;
use fs_feec::highorder::hex::TensorSpace;
use fs_feec::highorder::quad1d::{gauss_legendre, lobatto_shapes};
use fs_qty::Dims;

const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-12 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn linear_coefficients(space: &TensorSpace, field: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
    let mut values = vec![0.0; space.ndof()];
    for x in 0..=space.m() {
        for y in 0..=space.m() {
            for z in 0..=space.m() {
                let point = [x, y, z].map(|i| i as f64 / space.m() as f64);
                values[space.gid(x * space.r(), y * space.r(), z * space.r())] = field(point);
            }
        }
    }
    values
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn density([x, y, z]: [f64; 3]) -> f64 {
    1.0 + x * y + y * z + z * x
}

#[test]
fn nonseparable_density_retains_inventory_centroid_and_kinetic_energy() {
    for (m, r) in [(1, 1), (2, 3), (2, 6)] {
        let space = TensorSpace::new(m, r);
        let mut calls = 0;
        let mass = space
            .prepare_weighted_mass(
                DENSITY,
                r + 2,
                1_000_000,
                &mut |_, point| {
                    calls += 1;
                    Ok::<_, WeightedError>(density(point))
                },
                &mut || false,
            )
            .unwrap();
        assert_eq!(calls, m.pow(3) * (r + 2).pow(3));
        assert_eq!(mass.entry_dims(), Dims([0, 1, 0, 0, 0, 0]));
        assert_eq!(
            mass.output_dims(Dims([1, 0, -2, 0, 0, 0])).unwrap(),
            Dims([1, 1, -2, 0, 0, 0])
        );
        assert_eq!(mass.workspace_len(), space.ndof() + 2 * (r + 2).pow(3));
        let one = linear_coefficients(&space, |_| 1.0);
        let x = linear_coefficients(&space, |p| p[0]);
        let y = linear_coefficients(&space, |p| p[1]);
        let z = linear_coefficients(&space, |p| p[2]);
        let mut load = vec![0.0; space.ndof()];
        let mut workspace = vec![0.0; mass.workspace_len()];
        mass.apply_into(&one, &mut load, &mut workspace, &mut || false)
            .unwrap();
        // Exact monomial integrals on [0,1]^3, independent of basis and rule.
        near(dot(&one, &load), 7.0 / 4.0);
        for coordinate in [&x, &y, &z] {
            near(dot(coordinate, &load), 23.0 / 24.0);
        }
        mass.apply_into(&x, &mut load, &mut workspace, &mut || false)
            .unwrap();
        near(dot(&x, &load), 2.0 / 3.0);
        near(dot(&y, &load), 19.0 / 36.0);
        let first = load.clone();
        workspace.fill(-1234.0);
        mass.apply_into(&x, &mut load, &mut workspace, &mut || false)
            .unwrap();
        assert_eq!(first, load);
        assert_eq!(calls, m.pow(3) * (r + 2).pow(3));
    }
}

#[test]
fn constant_density_matches_the_existing_assembled_1d_mass_product() {
    for (m, r) in [(2, 2), (1, 6)] {
        let space = TensorSpace::new(m, r);
        let mass = space
            .prepare_weighted_mass(
                DENSITY,
                r + 2,
                1_000_000,
                &mut |_, _| Ok::<_, WeightedError>(4.5),
                &mut || false,
            )
            .unwrap();
        let input: Vec<_> = (0..space.ndof())
            .map(|i| (i % 17) as f64 / 9.0 - 0.6)
            .collect();
        let mut actual = vec![0.0; space.ndof()];
        let mut workspace = vec![0.0; mass.workspace_len()];
        mass.apply_into(&input, &mut actual, &mut workspace, &mut || false)
            .unwrap();
        let (one_d, _) = space.assembled_1d();
        let n = space.n1();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let mut expected = 0.0;
                    for x in 0..n {
                        for y in 0..n {
                            for z in 0..n {
                                expected += 4.5
                                    * one_d[i * n + x]
                                    * one_d[j * n + y]
                                    * one_d[k * n + z]
                                    * input[space.gid(x, y, z)];
                            }
                        }
                    }
                    near(actual[space.gid(i, j, k)], expected);
                    near(
                        mass.diagonal()[space.gid(i, j, k)],
                        4.5 * one_d[i * n + i] * one_d[j * n + j] * one_d[k * n + k],
                    );
                }
            }
        }
    }
}

#[test]
fn rectangular_contractions_match_direct_consistent_quadrature_including_bubbles() {
    let space = TensorSpace::new(2, 2);
    let q = 4;
    let mass = space
        .prepare_weighted_mass(
            DENSITY,
            q,
            1_000_000,
            &mut |_, point| Ok::<_, WeightedError>(density(point)),
            &mut || false,
        )
        .unwrap();
    let input: Vec<_> = (0..space.ndof())
        .map(|i| (i % 11) as f64 / 5.0 - 0.9)
        .collect();
    let mut actual = vec![0.0; space.ndof()];
    mass.apply_into(
        &input,
        &mut actual,
        &mut vec![0.0; mass.workspace_len()],
        &mut || false,
    )
    .unwrap();
    let mut expected = vec![0.0; space.ndof()];
    let mut diagonal = vec![0.0; space.ndof()];
    let (nodes, weights) = gauss_legendre(q);
    let basis: Vec<_> = nodes
        .iter()
        .map(|&x| lobatto_shapes(space.r(), x).0)
        .collect();
    let h = 1.0 / space.m() as f64;
    for cx in 0..space.m() {
        for cy in 0..space.m() {
            for cz in 0..space.m() {
                for i in 0..q {
                    for j in 0..q {
                        for k in 0..q {
                            let point = [
                                (cx as f64 + (nodes[i] + 1.0) / 2.0) * h,
                                (cy as f64 + (nodes[j] + 1.0) / 2.0) * h,
                                (cz as f64 + (nodes[k] + 1.0) / 2.0) * h,
                            ];
                            let weight = density(point)
                                * weights[i]
                                * weights[j]
                                * weights[k]
                                * (h / 2.0).powi(3);
                            let mut local = Vec::new();
                            for x in 0..=space.r() {
                                for y in 0..=space.r() {
                                    for z in 0..=space.r() {
                                        let global = space.gid(
                                            space.lat1(cx, x),
                                            space.lat1(cy, y),
                                            space.lat1(cz, z),
                                        );
                                        local.push((
                                            global,
                                            basis[i][x] * basis[j][y] * basis[k][z],
                                        ));
                                    }
                                }
                            }
                            let trial: f64 = local.iter().map(|&(g, b)| input[g] * b).sum();
                            for &(global, b) in &local {
                                expected[global] += weight * b * trial;
                                diagonal[global] += weight * b * b;
                            }
                        }
                    }
                }
            }
        }
    }
    for ((&a, &b), (&d, &e)) in actual
        .iter()
        .zip(&expected)
        .zip(mass.diagonal().iter().zip(&diagonal))
    {
        near(a, b);
        near(d, e);
    }
    assert!(dot(&input, &actual) > 0.0);
}

#[test]
fn discontinuous_cell_density_is_not_smoothed_across_interfaces() {
    let space = TensorSpace::new(2, 3);
    let mass = space
        .prepare_weighted_mass(
            DENSITY,
            5,
            1_000_000,
            &mut |cell, _| Ok::<_, WeightedError>(if cell[0] == 0 { 2.0 } else { 8.0 }),
            &mut || false,
        )
        .unwrap();
    let one = linear_coefficients(&space, |_| 1.0);
    let x = linear_coefficients(&space, |p| p[0]);
    let mut load = vec![0.0; mass.ndof()];
    mass.apply_into(
        &one,
        &mut load,
        &mut vec![0.0; mass.workspace_len()],
        &mut || false,
    )
    .unwrap();
    near(dot(&one, &load), 5.0);
    near(dot(&x, &load), 13.0 / 4.0);
}

#[test]
fn finite_extreme_density_is_scaled_after_integration() {
    let space = TensorSpace::new(1, 1);
    let mass = space
        .prepare_weighted_mass(
            DENSITY,
            3,
            1_000_000,
            &mut |_, _| Ok::<_, WeightedError>(1e308),
            &mut || false,
        )
        .unwrap();
    for &value in mass.diagonal() {
        near(value / 1e308, 1.0 / 27.0);
    }
    let one = linear_coefficients(&space, |_| 1.0);
    let mut load = vec![0.0; mass.ndof()];
    mass.apply_into(
        &one,
        &mut load,
        &mut vec![0.0; mass.workspace_len()],
        &mut || false,
    )
    .unwrap();
    near(dot(&one, &load) / 1e308, 1.0);
}

#[derive(Debug, PartialEq)]
enum SourceError {
    Native(WeightedError),
    MissingMaterial,
}
impl From<WeightedError> for SourceError {
    fn from(error: WeightedError) -> Self {
        Self::Native(error)
    }
}

#[test]
fn admission_preserves_source_errors_and_application_refusal_is_atomic() {
    let space = TensorSpace::new(2, 2);
    let mut calls = 0;
    let refused = space.prepare_weighted_mass(
        DENSITY,
        4,
        0,
        &mut |_, _| {
            calls += 1;
            Ok::<_, WeightedError>(1.0)
        },
        &mut || false,
    );
    assert!(matches!(refused, Err(WeightedError::LimitExceeded { .. })));
    assert_eq!(calls, 0);
    assert!(matches!(
        space.prepare_weighted_mass(
            DENSITY,
            2,
            1_000_000,
            &mut |_, _| Ok::<_, WeightedError>(1.0),
            &mut || false
        ),
        Err(WeightedError::InvalidGeometry { .. })
    ));
    assert!(matches!(
        space.prepare_weighted_mass(
            DENSITY,
            4,
            1_000_000,
            &mut |_, _| Err::<f64, _>(SourceError::MissingMaterial),
            &mut || false
        ),
        Err(SourceError::MissingMaterial)
    ));
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            space.prepare_weighted_mass(
                DENSITY,
                4,
                1_000_000,
                &mut |_, _| Ok::<_, WeightedError>(value),
                &mut || false
            ),
            Err(WeightedError::InvalidCoefficient { .. })
        ));
    }
    let mass = space
        .prepare_weighted_mass(
            DENSITY,
            4,
            1_000_000,
            &mut |_, _| Ok::<_, WeightedError>(1.0),
            &mut || false,
        )
        .unwrap();
    let mut input = linear_coefficients(&space, |_| 1.0);
    let mut output = vec![42.0; mass.ndof()];
    let unchanged = output.clone();
    let mut workspace = vec![0.0; mass.workspace_len()];
    let mut polls = 0;
    assert_eq!(
        mass.apply_into(&input, &mut output, &mut workspace, &mut || {
            polls += 1;
            polls == 5
        }),
        Err(WeightedError::Cancelled)
    );
    assert_eq!(output, unchanged);
    input[4] = f64::NAN;
    assert_eq!(
        mass.apply_into(&input, &mut output, &mut workspace, &mut || false),
        Err(WeightedError::NonFiniteEntry)
    );
    assert_eq!(output, unchanged);
}
