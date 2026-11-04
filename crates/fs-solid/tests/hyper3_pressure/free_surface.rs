//! Independent clipped-volume references for prescribed planar immersion.

use super::*;

fn field(pressure_pa: f64, gradient_pa_m: [f64; 3]) -> TetPressure<'static> {
    TetPressure::FreeSurface {
        origin_m: [0.0; 3],
        pressure_pa,
        gradient_pa_m,
    }
}

fn displacement(points: [[f64; 3]; 4]) -> Vec<f64> {
    points
        .iter()
        .zip(NODES)
        .flat_map(|(point, reference)| (0..3).map(move |axis| point[axis] - reference[axis]))
        .collect()
}

/// Truncated-power divided differences on a simplex, independently of the
/// production wet-triangle/tetrahedron decompositions. Distinct pressure nodes
/// are chosen for this oracle; repeated levels use separate analytic fixtures.
/// E[p_+] = sum_i (p_i)_+^4 / (4 product_{j!=i}(p_i-p_j)).
fn positive_moment_oracle(
    points: [[f64; 3]; 4],
    pressure: f64,
    gradient: [f64; 3],
) -> (f64, f64, [[f64; 3]; 4]) {
    let values = points.map(|point| pressure + dot(gradient, point));
    let denominators: [f64; 4] = std::array::from_fn(|i| {
        (0..4)
            .filter(|&j| j != i)
            .map(|j| {
                let difference = values[i] - values[j];
                assert!(difference.abs() > 1.0e-4);
                difference
            })
            .product()
    });
    let terms: [f64; 4] = std::array::from_fn(|i| values[i].max(0.0).powi(4) / denominators[i]);
    let mean = terms.into_iter().sum::<f64>() / 4.0;
    let fraction = (0..4)
        .map(|i| values[i].max(0.0).powi(3) / denominators[i])
        .sum::<f64>();
    let (volume, mut force) = tet_oracle(points, mean);
    for node in 0..4 {
        let mut derivative = values[node].max(0.0).powi(3) / denominators[node];
        for other in 0..4 {
            if other != node {
                derivative -= terms[node] / (4.0 * (values[node] - values[other]));
                derivative += terms[other] / (4.0 * (values[other] - values[node]));
            }
        }
        for axis in 0..3 {
            force[node][axis] -= volume * derivative * gradient[axis];
        }
    }
    (volume * mean, volume * fraction, force)
}

#[test]
fn clipped_force_tangent_and_work_match_independent_volume_divided_differences() {
    with_cx(|_, cx| {
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let points = [
            [0.1, -0.2, 0.3],
            [1.3, -0.1, 0.35],
            [0.25, 0.7, 0.25],
            [0.15, -0.05, 1.4],
        ];
        let u = displacement(points);
        let gradient = [-0.7, 1.1, -1.8];
        for (wet_count, pressure) in [(1, 0.1), (2, 1.2), (3, 2.0)] {
            assert_eq!(
                points
                    .iter()
                    .filter(|&&point| pressure + dot(gradient, point) > 0.0)
                    .count(),
                wet_count
            );
            let loading = field(pressure, gradient);
            let (potential, wet_volume, expected) =
                positive_moment_oracle(points, pressure, gradient);
            let force = surface.force_n(&u, loading, cx).unwrap();
            for (actual, expected) in force.iter().zip(expected.iter().flatten()) {
                near(*actual, *expected, 2.0e-13);
            }
            for axis in 0..3 {
                near(
                    (0..4).map(|node| force[3 * node + axis]).sum(),
                    -wet_volume * gradient[axis],
                    2.0e-13,
                );
            }
            let tangent = surface.force_tangent(&u, loading, cx).unwrap();
            let epsilon = 1.0e-6;
            for column in 0..12 {
                let mut plus = points;
                let mut minus = points;
                plus[column / 3][column % 3] += epsilon;
                minus[column / 3][column % 3] -= epsilon;
                let (forward_potential, _, forward) =
                    positive_moment_oracle(plus, pressure, gradient);
                let (backward_potential, _, backward) =
                    positive_moment_oracle(minus, pressure, gradient);
                near(
                    force[column],
                    -(forward_potential - backward_potential) / (2.0 * epsilon),
                    5.0e-9,
                );
                for row in 0..12 {
                    near(
                        tangent.get(row, column),
                        (forward[row / 3][row % 3] - backward[row / 3][row % 3]) / (2.0 * epsilon),
                        5.0e-9,
                    );
                    near(tangent.get(row, column), tangent.get(column, row), 3.0e-14);
                }
            }
            let work = surface.work_j(&[0.0; 12], &u, loading, cx).unwrap();
            near(
                work,
                positive_moment_oracle(NODES, pressure, gradient).0 - potential,
                2.0e-13,
            );
            near(
                surface.work_j(&u, &[0.0; 12], loading, cx).unwrap(),
                -work,
                2.0e-14,
            );
            let middle: Vec<f64> = u
                .iter()
                .enumerate()
                .map(|(i, value)| 0.3 * value + if i % 3 == 1 { 0.25 } else { 0.0 })
                .collect();
            near(
                surface.work_j(&[0.0; 12], &middle, loading, cx).unwrap()
                    + surface.work_j(&middle, &u, loading, cx).unwrap(),
                work,
                2.0e-13,
            );
        }
    });
}

/// Unit-tet integrals below z=depth: wet volume, first x moment, integral depth.
fn horizontal_integrals(depth: f64) -> (f64, f64, f64) {
    if depth <= 0.0 {
        return (0.0, 0.0, 0.0);
    }
    if depth >= 1.0 {
        return (1.0 / 6.0, 1.0 / 24.0, depth / 6.0 - 1.0 / 24.0);
    }
    let volume = depth * (3.0 - 3.0 * depth + depth * depth) / 6.0;
    let moment_x = depth * (4.0 - 6.0 * depth + 4.0 * depth * depth - depth.powi(3)) / 24.0;
    let potential = depth * depth * (6.0 - 4.0 * depth + depth * depth) / 24.0;
    (volume, moment_x, potential)
}

#[test]
fn immersion_sweep_has_exact_archimedes_force_torque_and_crossing_work() {
    with_cx(|_, cx| {
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let gradient = 7.3;
        let level = 0.6;
        let loading = field(gradient * level, [0.0, 0.0, -gradient]);
        let translated = |depth| -> Vec<f64> {
            (0..12)
                .map(|i| if i % 3 == 2 { level - depth } else { 0.0 })
                .collect()
        };
        let dry = translated(-0.3);
        for depth in [-0.3, 0.0, 0.125, 0.5, 0.999, 1.0, 1.4] {
            let u = translated(depth);
            let force = surface.force_n(&u, loading, cx).unwrap();
            let (volume, moment, potential) = horizontal_integrals(depth);
            let resultant: [f64; 3] =
                std::array::from_fn(|axis| (0..4).map(|node| force[3 * node + axis]).sum());
            near(resultant[0], 0.0, 2.0e-14);
            near(resultant[1], 0.0, 2.0e-14);
            near(resultant[2], gradient * volume, 2.0e-14);
            let mut torque = [0.0; 3];
            for node in 0..4 {
                let current = std::array::from_fn(|axis| NODES[node][axis] + u[3 * node + axis]);
                let moment = cross(current, std::array::from_fn(|axis| force[3 * node + axis]));
                for axis in 0..3 {
                    torque[axis] += moment[axis];
                }
            }
            near(torque[0], gradient * moment, 2.0e-14);
            near(torque[1], -gradient * moment, 2.0e-14);
            near(torque[2], 0.0, 2.0e-14);
            near(
                surface.work_j(&dry, &u, loading, cx).unwrap(),
                -gradient * potential,
                2.0e-14,
            );
        }
        // A step crossing the plane is not a cubic-force Simpson interval.
        let after = translated(0.8);
        let midpoint: Vec<f64> = dry
            .iter()
            .zip(&after)
            .map(|(a, b)| 0.5 * a + 0.5 * b)
            .collect();
        let f0 = surface.force_n(&dry, loading, cx).unwrap();
        let fm = surface.force_n(&midpoint, loading, cx).unwrap();
        let f1 = surface.force_n(&after, loading, cx).unwrap();
        let simpson = (0..12)
            .map(|i| (f0[i] + 4.0 * fm[i] + f1[i]) * (after[i] - dry[i]) / 6.0)
            .sum::<f64>();
        let work = surface.work_j(&dry, &after, loading, cx).unwrap();
        assert!((simpson - work).abs() > 1.0e-3);
        near(work, -gradient * horizontal_integrals(0.8).2, 2.0e-14);
    });
}

#[test]
fn one_and_two_wet_vertex_caps_remain_accurate_at_small_immersion() {
    with_cx(|_, cx| {
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let gradient = 2.7;
        for height in [0.1, 0.5, 0.9, 0.999_999] {
            let cap = 1.0 - height;
            for two_wet in [false, true] {
                let direction = if two_wet {
                    [gradient, gradient, 0.0]
                } else {
                    [0.0, 0.0, gradient]
                };
                let loading = field(-gradient * height, direction);
                let volume = if two_wet {
                    cap * cap * (1.0 + 2.0 * height) / 6.0
                } else {
                    cap.powi(3) / 6.0
                };
                let potential = if two_wet {
                    gradient * cap.powi(3) * (1.0 + height) / 12.0
                } else {
                    gradient * cap.powi(4) / 24.0
                };
                let force = surface.force_n(&[0.0; 12], loading, cx).unwrap();
                for axis in 0..3 {
                    near(
                        (0..4).map(|node| force[3 * node + axis]).sum(),
                        -volume * direction[axis],
                        1.0e-9 * gradient * volume.max(1.0e-30),
                    );
                }
                let dry: Vec<f64> = (0..12)
                    .map(|i| {
                        if (two_wet && i % 3 == 0) || (!two_wet && i % 3 == 2) {
                            -1.0
                        } else {
                            0.0
                        }
                    })
                    .collect();
                near(
                    surface.work_j(&dry, &[0.0; 12], loading, cx).unwrap(),
                    -potential,
                    1.0e-9 * potential.max(1.0e-30),
                );
            }
        }
    });
}

#[test]
fn free_surface_limits_admission_and_coplanar_tangent_are_explicit() {
    with_cx(|gate, cx| {
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let u: Vec<f64> = (0..12).map(|i| 0.003 * i as f64).collect();
        let gradient = [0.3, -0.2, 0.4];
        let wet = field(3.0, gradient);
        let affine = TetPressure::Affine {
            origin_m: [0.0; 3],
            pressure_pa: 3.0,
            gradient_pa_m: gradient,
        };
        assert_eq!(
            surface.force_n(&u, wet, cx).unwrap(),
            surface.force_n(&u, affine, cx).unwrap()
        );
        let wet_tangent = surface.force_tangent(&u, wet, cx).unwrap();
        let affine_tangent = surface.force_tangent(&u, affine, cx).unwrap();
        for row in 0..12 {
            for column in 0..12 {
                assert_eq!(
                    wet_tangent.get(row, column),
                    affine_tangent.get(row, column)
                );
            }
        }
        near(
            surface.work_j(&[0.0; 12], &u, wet, cx).unwrap(),
            surface.work_j(&[0.0; 12], &u, affine, cx).unwrap(),
            2.0e-14,
        );
        for pressure in [0.0, 2.3, -2.3] {
            let constant = field(pressure, [0.0; 3]);
            let expected = pressure.max(0.0);
            assert_eq!(
                surface.validate_pressure(constant, cx).unwrap(),
                expected == 0.0
            );
            assert_eq!(
                surface.force_n(&u, constant, cx).unwrap(),
                surface.force_n(&u, expected, cx).unwrap()
            );
            assert_eq!(
                surface.work_j(&[0.0; 12], &u, constant, cx).unwrap(),
                surface.work_j(&[0.0; 12], &u, expected, cx).unwrap()
            );
        }
        let touching = field(0.0, [0.0, 0.0, -2.0]);
        assert!(
            !surface.validate_pressure(touching, cx).unwrap(),
            "a currently dry body may enter the fluid"
        );
        assert_eq!(
            surface.force_n(&[0.0; 12], touching, cx).unwrap(),
            vec![0.0; 12]
        );
        let tangent = surface.force_tangent(&[0.0; 12], touching, cx).unwrap();
        for row in 0..12 {
            for column in 0..12 {
                assert_eq!(tangent.get(row, column), 0.0);
            }
        }
        let entered: Vec<f64> = (0..12)
            .map(|i| if i % 3 == 2 { -0.01 } else { 0.0 })
            .collect();
        let force = surface.force_n(&entered, touching, cx).unwrap();
        near(
            (0..4).map(|node| force[3 * node + 2]).sum(),
            2.0 * horizontal_integrals(0.01).0,
            2.0e-14,
        );
        near(
            surface.work_j(&[0.0; 12], &entered, touching, cx).unwrap(),
            -2.0 * horizontal_integrals(0.01).2,
            2.0e-14,
        );
        for bad in [
            TetPressure::FreeSurface {
                origin_m: [f64::NAN, 0.0, 0.0],
                pressure_pa: 0.0,
                gradient_pa_m: gradient,
            },
            field(f64::INFINITY, gradient),
            field(0.0, [0.0, f64::NEG_INFINITY, 0.0]),
        ] {
            assert!(matches!(
                surface.force_n(&u, bad, cx),
                Err(HyperTetError::InvalidInput { .. })
            ));
            assert!(matches!(
                surface.force_tangent(&u, bad, cx),
                Err(HyperTetError::InvalidInput { .. })
            ));
            assert!(matches!(
                surface.work_j(&[0.0; 12], &u, bad, cx),
                Err(HyperTetError::InvalidInput { .. })
            ));
        }
        gate.request();
        assert!(matches!(
            surface.force_n(&u, wet, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            surface.force_tangent(&u, wet, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            surface.work_j(&[0.0; 12], &u, wet, cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}

#[test]
fn clipping_opposite_large_pressures_does_not_overflow_the_edge_difference() {
    with_cx(|_, cx| {
        let mut nodes = NODES;
        nodes[1][0] = 2.0;
        let surface =
            ClosedTetPressureSurface::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let gradient = 1.7e308;
        let loading = TetPressure::FreeSurface {
            origin_m: [1.0, 0.0, 0.0],
            pressure_pa: 0.0,
            gradient_pa_m: [gradient, 0.0, 0.0],
        };
        let force = surface.force_n(&[0.0; 12], loading, cx).unwrap();
        assert!(force.iter().all(|value| value.is_finite()));
        near(
            (0..4).map(|node| force[3 * node]).sum(),
            -gradient / 24.0,
            gradient * 2.0e-15,
        );
        let tangent = surface.force_tangent(&[0.0; 12], loading, cx).unwrap();
        for row in 0..12 {
            for column in 0..12 {
                assert!(tangent.get(row, column).is_finite());
            }
        }
        let mut dry = [0.0; 12];
        dry[3] = -1.0;
        near(
            surface.work_j(&dry, &[0.0; 12], loading, cx).unwrap(),
            -gradient / 96.0,
            gradient * 2.0e-15,
        );
    });
}

#[test]
fn native_newton_equilibrium_crosses_a_wet_vertex_and_retains_support_reactions() {
    with_cx(|_, cx| {
        let (mu, lambda) = (2.0, 3.0);
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let points = [
            [0.0; 3],
            [0.92, 0.0, 0.0],
            [0.05, 0.95, 0.0],
            [0.02, -0.03, 0.94],
        ];
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let free = [3, 6, 7, 9, 10, 11];
        let jacobian = dot(points[1], cross(points[2], points[3]));
        let coefficient = (lambda * jacobian.ln() - mu) / jacobian;
        let cofactors = [
            cross(points[2], points[3]),
            cross(points[3], points[1]),
            cross(points[1], points[2]),
        ];
        let mut internal = [0.0; 12];
        for column in 0..3 {
            for axis in 0..3 {
                let force =
                    (mu * points[column + 1][axis] + coefficient * cofactors[column][axis]) / 6.0;
                internal[3 * (column + 1) + axis] = force;
                internal[axis] -= force;
            }
        }
        let (pressure, gradient) = (0.19, [0.15, -0.2, -0.25]);
        assert!(pressure + dot(gradient, NODES[2]) < 0.0);
        assert!(pressure + dot(gradient, points[2]) > 0.0);
        let applied = positive_moment_oracle(points, pressure, gradient).2;
        let dead: Vec<(usize, f64)> = (0..12)
            .map(|dof| {
                (
                    dof,
                    if free.contains(&dof) {
                        internal[dof] - applied[dof / 3][dof % 3]
                    } else {
                        0.13
                    },
                )
            })
            .collect();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let problem = HyperTetProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            material: HyperelasticLawField::Uniform(&material),
            prescribed_m: &fixed,
            nodal_forces_n: &dead,
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings {
                load_steps: 4,
                ..HyperTetSettings::default()
            },
        };
        let solution = problem
            .solve_with_pressure_field(&surface, field(pressure, gradient), cx)
            .unwrap();
        for (node, u) in solution.displacement_m.iter().enumerate() {
            for axis in 0..3 {
                near(u[axis], points[node][axis] - NODES[node][axis], 3.0e-9);
            }
        }
        for (dof, reaction) in solution.reactions_n {
            near(
                reaction,
                internal[dof] - applied[dof / 3][dof % 3] - dead[dof].1,
                3.0e-9,
            );
        }
        near(solution.min_det_f, jacobian, 3.0e-9);
    });
}
