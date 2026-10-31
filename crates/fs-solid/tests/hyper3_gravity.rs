//! G0/G1/G3: original reference mass, current-position gravity derivatives,
//! conservative equilibrium and quadrature refinement on native hyper3.
#[path = "tet10/support.rs"]
mod support;
use support::{block, with_cx};

use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::gravity::GravitySample3;
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_solid::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};

fn material() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 3.0,
        },
        5.0,
    )
    .unwrap()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    law: &'a Hyperelastic,
    prescribed: &'a [(usize, f64)],
) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        material: HyperelasticLawField::Uniform(law),
        prescribed_m: prescribed,
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: Default::default(),
        settings: HyperTetSettings::default(),
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance={tolerance:e}"
    );
}

fn affine(x: [f64; 3]) -> Result<GravitySample3, HyperTetError> {
    let gradient = [[-0.4, 0.1, 0.0], [0.1, -0.7, 0.05], [0.0, 0.05, -1.2]];
    let constant = [0.3, -0.5, -2.0];
    let acceleration =
        std::array::from_fn(|i| constant[i] + (0..3).map(|j| gradient[i][j] * x[j]).sum::<f64>());
    let potential = -(0..3)
        .map(|i| 0.5 * x[i] * (acceleration[i] + constant[i]))
        .sum::<f64>();
    Ok(GravitySample3 {
        potential_j_kg: potential,
        acceleration_m_s2: acceleration,
        acceleration_gradient_s2: gradient,
    })
}

#[test]
fn g0_affine_gravity_uses_full_heterogeneous_reference_mass_at_current_positions() {
    let mesh = block([1; 3], [1.0; 3]);
    let law = material();
    let p = problem(&mesh.nodes, &mesh.tets, &law, &[(0, 0.0)]);
    let density: Vec<_> = (0..mesh.tets.len()).map(|i| 2.0 + i as f64).collect();
    let u: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|x| {
            [
                0.2 * x[0] + 0.1 * x[2] + 0.7,
                -0.1 * x[1] - 0.2,
                0.15 * x[2] + 0.05 * x[0] + 0.3,
            ]
        })
        .collect();
    let positions: Vec<_> = mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(node, x)| std::array::from_fn(|i| x[i] + u[3 * node + i]))
        .collect();
    let acceleration: Vec<_> = positions
        .iter()
        .flat_map(|&x| affine(x).unwrap().acceleration_m_s2)
        .collect();
    with_cx(|cx| {
        let mass = p.reference_mass_matrix(&density, cx).unwrap();
        let prepared = p.prepare(cx).unwrap();
        let load = prepared.gravity_load(&density, &u, &affine, cx).unwrap();
        let mut expected = vec![0.0; u.len()];
        mass.spmv(&acceleration, &mut expected);
        for (&actual, &expected) in load.force_n.iter().zip(&expected) {
            near(actual, expected, 2e-14);
        }
        assert!(
            load.force_n[0].abs() > 0.01,
            "supported loads must remain in the full vector"
        );
        let gradient = affine([0.0; 3]).unwrap().acceleration_gradient_s2;
        for row in 0..u.len() {
            for column in 0..u.len() {
                near(
                    load.force_tangent_n_m.get(row, column),
                    mass.get(3 * (row / 3), 3 * (column / 3)) * gradient[row % 3][column % 3],
                    2e-14,
                );
            }
        }
        let constant = affine([0.0; 3]).unwrap().acceleration_m_s2;
        let expected_potential: f64 = (0..u.len())
            .map(|row| {
                let row_mass: f64 = mass.row(row).1.iter().sum();
                -0.5 * positions[row / 3][row % 3] * (expected[row] + row_mass * constant[row % 3])
            })
            .sum();
        near(load.potential_energy_j, expected_potential, 2e-13);
        let mut doubled = u.clone();
        for (node, x) in mesh.nodes.iter().enumerate() {
            for i in 0..3 {
                doubled[3 * node + i] += x[i];
            }
        }
        let uniform = |x: [f64; 3]| {
            Ok(GravitySample3 {
                potential_j_kg: -0.3 * x[0] + 0.4 * x[1] + 9.8 * x[2],
                acceleration_m_s2: [0.3, -0.4, -9.8],
                acceleration_gradient_s2: [[0.0; 3]; 3],
            })
        };
        assert_eq!(
            prepared
                .gravity_load(&density, &doubled, &uniform, cx)
                .unwrap()
                .force_n,
            prepared
                .gravity_load(&density, &u, &uniform, cx)
                .unwrap()
                .force_n,
            "stretching changes current volume, never the reference weight"
        );
    });
}

fn quartic(x: [f64; 3]) -> Result<GravitySample3, HyperTetError> {
    let radius2: f64 = x.iter().map(|v| v * v).sum();
    Ok(GravitySample3 {
        potential_j_kg: 0.25 * radius2 * radius2,
        acceleration_m_s2: x.map(|v| -radius2 * v),
        acceleration_gradient_s2: std::array::from_fn(|i| {
            std::array::from_fn(|j| -2.0 * x[i] * x[j] - if i == j { radius2 } else { 0.0 })
        }),
    })
}

#[test]
fn g0_g3_nonlinear_force_tangent_and_potential_share_one_current_configuration() {
    let mesh = block([1; 3], [1.0; 3]);
    let law = material();
    let p = problem(&mesh.nodes, &mesh.tets, &law, &[]);
    let density = vec![2.0; mesh.tets.len()];
    let u: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|x| [0.1 + 0.2 * x[1], -0.05 * x[0], 0.1 * x[2]])
        .collect();
    with_cx(|cx| {
        let prepared = p.prepare(cx).unwrap();
        let load = prepared.gravity_load(&density, &u, &quartic, cx).unwrap();
        let mut moment = [0.0; 3];
        for (node, reference) in mesh.nodes.iter().enumerate() {
            let x: [f64; 3] = std::array::from_fn(|i| reference[i] + u[3 * node + i]);
            let f = &load.force_n[3 * node..3 * node + 3];
            for i in 0..3 {
                moment[i] += x[(i + 1) % 3] * f[(i + 2) % 3] - x[(i + 2) % 3] * f[(i + 1) % 3];
            }
        }
        for value in moment {
            near(value, 0.0, 2e-13);
        }
        for column in 0..u.len() {
            let mut plus = u.clone();
            let mut minus = u.clone();
            plus[column] += 1e-6;
            minus[column] -= 1e-6;
            let a = prepared
                .gravity_load(&density, &plus, &quartic, cx)
                .unwrap();
            let b = prepared
                .gravity_load(&density, &minus, &quartic, cx)
                .unwrap();
            near(
                -(a.potential_energy_j - b.potential_energy_j) / 2e-6,
                load.force_n[column],
                2e-9,
            );
            for row in 0..u.len() {
                near(
                    (a.force_n[row] - b.force_n[row]) / 2e-6,
                    load.force_tangent_n_m.get(row, column),
                    2e-9,
                );
                near(
                    load.force_tangent_n_m.get(row, column),
                    load.force_tangent_n_m.get(column, row),
                    2e-14,
                );
            }
        }
        let rotate = |x: [f64; 3]| [-x[1], x[0], x[2]];
        let nodes: Vec<_> = mesh.nodes.iter().copied().map(rotate).collect();
        let turned_u: Vec<_> = u
            .chunks_exact(3)
            .flat_map(|x| rotate([x[0], x[1], x[2]]))
            .collect();
        let turned = problem(&nodes, &mesh.tets, &law, &[])
            .gravity_load(&density, &turned_u, &quartic, cx)
            .unwrap();
        near(turned.potential_energy_j, load.potential_energy_j, 2e-14);
        for (actual, original) in turned
            .force_n
            .chunks_exact(3)
            .zip(load.force_n.chunks_exact(3))
        {
            for (&actual, expected) in
                actual
                    .iter()
                    .zip(rotate([original[0], original[1], original[2]]))
            {
                near(actual, expected, 2e-14);
            }
        }
    });
}

#[test]
fn g1_native_newton_retains_gravity_gradient_pressure_and_supported_weight() {
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let tets = [[0, 1, 2, 3]];
    let law = material();
    let prescribed: Vec<_> = (0..12).filter(|&d| d != 3).map(|d| (d, 0.0)).collect();
    let p = HyperTetProblem {
        settings: HyperTetSettings {
            load_steps: 3,
            ..HyperTetSettings::default()
        },
        ..problem(&nodes, &tets, &law, &prescribed)
    };
    let stretch = 1.25_f64;
    let gradient = 3.0;
    let piola = 2.0 * (stretch - 1.0 / stretch) + 3.0 * stretch.ln() / stretch;
    with_cx(|cx| {
        let surface =
            ClosedTetPressureSurface::new(&nodes, &tets, Default::default(), 4, cx).unwrap();
        for pressure in [0.0, 0.2] {
            for datum in [0.0, 1e12] {
                // Exact consistent mass projection on this single stretching tet:
                // f_1x/V = a/4 - k*stretch/10, with rho=1.
                let a = 4.0 * (piola + pressure) + 0.4 * gradient * stretch;
                let field = |x: [f64; 3]| {
                    Ok(GravitySample3 {
                        potential_j_kg: 0.5 * gradient * x[0] * x[0] - a * x[0] + datum,
                        acceleration_m_s2: [a - gradient * x[0], 0.0, 0.0],
                        acceleration_gradient_s2: [[-gradient, 0.0, 0.0], [0.0; 3], [0.0; 3]],
                    })
                };
                let solution = if pressure == 0.0 {
                    p.solve_with_gravity(&[1.0], &field, cx).unwrap()
                } else {
                    p.solve_with_gravity_and_pressure(
                        &[1.0],
                        &field,
                        &surface,
                        TetPressure::Uniform(pressure),
                        cx,
                    )
                    .unwrap()
                };
                near(solution.displacement_m[1][0], stretch - 1.0, 2e-10);
                let u: Vec<_> = solution.displacement_m.iter().flatten().copied().collect();
                let gravity = p.gravity_load(&[1.0], &u, &field, cx).unwrap();
                let pressure_force = surface.force_n(&u, pressure, cx).unwrap();
                for i in 0..3 {
                    let reaction: f64 = solution
                        .reactions_n
                        .iter()
                        .filter(|(d, _)| d % 3 == i)
                        .map(|(_, r)| r)
                        .sum();
                    near(
                        reaction
                            + (0..4)
                                .map(|node| {
                                    gravity.force_n[3 * node + i] + pressure_force[3 * node + i]
                                })
                                .sum::<f64>(),
                        0.0,
                        2e-10,
                    );
                }
                assert!(
                    solution
                        .steps
                        .last()
                        .unwrap()
                        .residual_history
                        .last()
                        .unwrap()
                        < &1e-10
                );
            }
        }
    });
}

#[test]
fn g1_nonlinear_gravity_quadrature_refines_with_original_material_mesh() {
    let law = material();
    let mut errors = Vec::new();
    for n in [1, 2, 4] {
        let mesh = block([n; 3], [1.0; 3]);
        let p = problem(&mesh.nodes, &mesh.tets, &law, &[]);
        let field = |x: [f64; 3]| {
            let x2 = x[0] * x[0];
            let x4 = x2 * x2;
            Ok(GravitySample3 {
                potential_j_kg: x4 * x2 / 6.0,
                acceleration_m_s2: [-x4 * x[0], 0.0, 0.0],
                acceleration_gradient_s2: [[-5.0 * x4, 0.0, 0.0], [0.0; 3], [0.0; 3]],
            })
        };
        with_cx(|cx| {
            let load = p
                .gravity_load(
                    &vec![1.0; mesh.tets.len()],
                    &vec![0.0; 3 * mesh.nodes.len()],
                    &field,
                    cx,
                )
                .unwrap();
            let force: f64 = load.force_n.chunks_exact(3).map(|f| f[0]).sum();
            errors.push([
                (force + 1.0 / 6.0).abs(),
                (load.potential_energy_j - 1.0 / 42.0).abs(),
            ]);
        });
    }
    for pair in errors.windows(2) {
        for i in 0..2 {
            assert!(pair[1][i] < pair[0][i] / 3.0, "{errors:?}");
        }
    }
}

#[test]
fn g1_potential_datum_cannot_reject_a_resolved_newton_decrement() {
    // Native Neo-Hookean SIMPLE SHEAR has exactly quadratic energy: with
    // mu=6 and V=1/6, E=u²/2. Reference mass=4 and g_x=1 give the sole free
    // node force 1 and U=4*datum-1-u. Every ramp step has true decrement
    // -0.5/123², despite larger cancelling elastic/gravity increments. At
    // step 57 the rounded quadrature gives dE≈0.00373455 and dUg≈-0.00373357:
    // a false positive total change, with both force and tangent unchanged.
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let tets = [[0, 1, 2, 3]];
    let law = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 6.0,
            lambda: 3.0,
        },
        5.0,
    )
    .unwrap();
    let prescribed: Vec<_> = (0..12).filter(|&d| d != 6).map(|d| (d, 0.0)).collect();
    let p = HyperTetProblem {
        settings: HyperTetSettings {
            load_steps: 123,
            max_backtracks: 0,
            ..HyperTetSettings::default()
        },
        ..problem(&nodes, &tets, &law, &prescribed)
    };
    with_cx(|cx| {
        for datum in [0.0, 2.5e11] {
            let field = |x: [f64; 3]| {
                Ok(GravitySample3 {
                    potential_j_kg: datum - x[0],
                    acceleration_m_s2: [1.0, 0.0, 0.0],
                    acceleration_gradient_s2: [[0.0; 3]; 3],
                })
            };
            let solution = p.solve_with_gravity(&[24.0], &field, cx).unwrap();
            near(solution.displacement_m[2][0], 1.0, 2e-12);
            near(solution.strain_energy_j, 0.5, 2e-12);
            assert!(solution.steps.iter().all(|step| step.backtracks == 0));
        }
    });
}

#[test]
fn g0_invalid_fields_refuse_without_rebinding_the_original_problem() {
    let mesh = block([1; 3], [1.0; 3]);
    let law = material();
    let p = problem(&mesh.nodes, &mesh.tets, &law, &[]);
    let density = vec![2.0; mesh.tets.len()];
    let zero = vec![0.0; 3 * mesh.nodes.len()];
    with_cx(|cx| {
        let expected = p.gravity_load(&density, &zero, &affine, cx).unwrap();
        assert!(p.gravity_load(&[], &zero, &affine, cx).is_err());
        assert!(
            p.gravity_load(&vec![-1.0; density.len()], &zero, &affine, cx)
                .is_err()
        );
        assert!(p.gravity_load(&density, &[], &affine, cx).is_err());
        for kind in 0..3 {
            let field = |x| {
                let mut s = affine(x)?;
                match kind {
                    0 => s.potential_j_kg = f64::NAN,
                    1 => s.acceleration_m_s2[0] = f64::INFINITY,
                    _ => s.acceleration_gradient_s2[0][1] += 1.0,
                }
                Ok(s)
            };
            assert!(p.gravity_load(&density, &zero, &field, cx).is_err());
        }
        let outside = |_| {
            Err(HyperTetError::InvalidInput {
                what: "outside source domain".into(),
            })
        };
        assert!(p.gravity_load(&density, &zero, &outside, cx).is_err());
        let cancelled = |_| Err(HyperTetError::Cancelled);
        assert!(matches!(
            p.gravity_load(&density, &zero, &cancelled, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert_eq!(
            p.gravity_load(&density, &zero, &affine, cx)
                .unwrap()
                .force_n,
            expected.force_n
        );
        let unloaded = |_| {
            Ok(GravitySample3 {
                potential_j_kg: 0.0,
                acceleration_m_s2: [0.0; 3],
                acceleration_gradient_s2: [[0.0; 3]; 3],
            })
        };
        let baseline = p.solve(cx).unwrap();
        let with_zero = p.solve_with_gravity(&density, &unloaded, cx).unwrap();
        assert_eq!(with_zero.displacement_m, baseline.displacement_m);
        assert_eq!(with_zero.reactions_n, baseline.reactions_n);
    });
}
