//! G0/G3: stratified fluid loading through native solid pressure owners.
use super::*;
use fs_solid::hyper3::pressure::{HydrostaticDensityLayer, StratifiedHydrostaticPressure};

const LAYERS: [HydrostaticDensityLayer; 3] = [
    HydrostaticDensityLayer {
        top_depth_m: 0.0,
        density_kg_m3: 1.3,
    },
    HydrostaticDensityLayer {
        top_depth_m: 0.2,
        density_kg_m3: 2.1,
    },
    // An admitted prescribed density inversion catches reversed clipping.
    HydrostaticDensityLayer {
        top_depth_m: 0.55,
        density_kg_m3: 0.9,
    },
];
const ORIGIN: [f64; 3] = [0.13, -0.17, 0.28];
const GRAVITY: [f64; 3] = [0.31, -0.57, -1.09];

fn profile(cx: &Cx<'_>) -> StratifiedHydrostaticPressure {
    StratifiedHydrostaticPressure::new(ORIGIN, GRAVITY, &LAYERS, 3, cx).unwrap()
}
fn surface(cx: &Cx<'_>) -> ClosedTetPressureSurface {
    ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap()
}
fn state(f: [f64; 9], translation: [f64; 3]) -> [f64; 12] {
    std::array::from_fn(|d| {
        let (node, axis) = (d / 3, d % 3);
        translation[axis]
            + (0..3)
                .map(|j| f[3 * axis + j] * NODES[node][j])
                .sum::<f64>()
            - NODES[node][axis]
    })
}
fn endpoints() -> ([f64; 12], [f64; 12]) {
    (
        state(
            [1.08, 0.06, -0.04, 0.02, 0.94, 0.05, -0.03, 0.04, 1.11],
            [0.1, -0.2, 0.55],
        ),
        state(
            [0.91, 0.1, 0.02, -0.04, 1.06, -0.08, 0.02, -0.06, 0.98],
            [-0.05, 0.04, -0.65],
        ),
    )
}

// Independent truncated-power divided difference for integral_tet p(x) dV.
// Production uses positive clipped simplices and current surface tractions.
// The test's oblique field has distinct vertex depths at both configurations.
fn potential(u: &[f64; 12]) -> f64 {
    let x: [[f64; 3]; 4] = std::array::from_fn(|node| {
        std::array::from_fn(|axis| NODES[node][axis] + u[3 * node + axis])
    });
    let volume = dot(sub(x[1], x[0]), cross(sub(x[2], x[0]), sub(x[3], x[0]))) / 6.0;
    let magnitude = dot(GRAVITY, GRAVITY).sqrt();
    let q = x.map(|point| dot(GRAVITY, sub(point, ORIGIN)));
    let mut mean = 0.0;
    let mut previous_density = 0.0;
    for layer in LAYERS {
        let values = q.map(|value| value - magnitude * layer.top_depth_m);
        let clipped_mean = (0..4)
            .map(|i| {
                let denominator = (0..4)
                    .filter(|&j| i != j)
                    .map(|j| values[i] - values[j])
                    .product::<f64>();
                values[i].max(0.0).powi(4) / (4.0 * denominator)
            })
            .sum::<f64>();
        mean += (layer.density_kg_m3 - previous_density) * clipped_mean;
        previous_density = layer.density_kg_m3;
    }
    volume * mean
}

#[test]
fn stratified_archimedes_resultant_torque_and_world_frame_are_consistent() {
    with_cx(|_, cx| {
        let nodes = [
            [0.2, -0.3, -0.8],
            [1.2, -0.3, -0.8],
            [1.2, 0.7, -0.8],
            [0.2, 0.7, -0.8],
            [0.2, -0.3, 0.2],
            [1.2, -0.3, 0.2],
            [1.2, 0.7, 0.2],
            [0.2, 0.7, 0.2],
        ];
        let tets = [
            [0, 1, 2, 6],
            [0, 2, 3, 6],
            [0, 3, 7, 6],
            [0, 7, 4, 6],
            [0, 4, 5, 6],
            [0, 5, 1, 6],
        ];
        let layers = [
            HydrostaticDensityLayer {
                top_depth_m: 0.0,
                density_kg_m3: 800.0,
            },
            HydrostaticDensityLayer {
                top_depth_m: 0.35,
                density_kg_m3: 1200.0,
            },
            HydrostaticDensityLayer {
                top_depth_m: 0.7,
                density_kg_m3: 950.0,
            },
        ];
        let origin = [0.0, 0.0, 0.1];
        let gravity = [0.0, 0.0, -9.0];
        let field = StratifiedHydrostaticPressure::new(origin, gravity, &layers, 3, cx).unwrap();
        let surface =
            ClosedTetPressureSurface::new(&nodes, &tets, TetAssemblyBudget::standard(), 24, cx)
                .unwrap();
        let force = surface.force_n(&[0.0; 24], field.loading(), cx).unwrap();
        let expected = 9.0 * (800.0 * 0.35 + 1200.0 * 0.35 + 950.0 * 0.2);
        let mut resultant = [0.0; 3];
        let mut torque = [0.0; 3];
        for node in 0..8 {
            let f = std::array::from_fn(|axis| force[3 * node + axis]);
            let moment = cross(nodes[node], f);
            for axis in 0..3 {
                resultant[axis] += f[axis];
                torque[axis] += moment[axis];
            }
        }
        near(resultant[0], 0.0, 2e-11);
        near(resultant[1], 0.0, 2e-11);
        near(resultant[2], expected, 2e-11);
        near(torque[0], 0.2 * expected, 2e-11);
        near(torque[1], -0.7 * expected, 2e-11);
        near(torque[2], 0.0, 2e-11);
        // Independent pressure integral over depth, including the density drop.
        near(
            field.pressure_pa([0.0, 0.0, -0.8], cx).unwrap(),
            expected,
            2e-11,
        );
        assert_eq!(field.pressure_pa([0.0, 0.0, 0.3], cx).unwrap(), 0.0);
        let rotate = |x: [f64; 3]| [x[2], x[0], x[1]];
        let rotated = ClosedTetPressureSurface::new(
            &nodes.map(rotate),
            &tets,
            TetAssemblyBudget::standard(),
            24,
            cx,
        )
        .unwrap();
        let field_rotated =
            StratifiedHydrostaticPressure::new(rotate(origin), rotate(gravity), &layers, 3, cx)
                .unwrap();
        let moved = rotated
            .force_n(&[0.0; 24], field_rotated.loading(), cx)
            .unwrap();
        for node in 0..8 {
            let expected = rotate(std::array::from_fn(|axis| force[3 * node + axis]));
            for axis in 0..3 {
                near(moved[3 * node + axis], expected[axis], 2e-11);
            }
        }
    });
}

#[test]
fn stratified_force_tangent_and_interface_crossing_work_match_volume_potential() {
    with_cx(|_, cx| {
        let surface = surface(cx);
        let profile = profile(cx);
        let field = profile.loading();
        let (before, after) = endpoints();
        let force = surface.force_n(&after, field, cx).unwrap();
        let tangent = surface.force_tangent(&after, field, cx).unwrap();
        let epsilon = 2e-6;
        for column in 0..12 {
            let mut plus = after;
            let mut minus = after;
            plus[column] += epsilon;
            minus[column] -= epsilon;
            near(
                force[column],
                -(potential(&plus) - potential(&minus)) / (2.0 * epsilon),
                2e-8,
            );
            let forward = surface.force_n(&plus, field, cx).unwrap();
            let backward = surface.force_n(&minus, field, cx).unwrap();
            for row in 0..12 {
                near(
                    tangent.get(row, column),
                    (forward[row] - backward[row]) / (2.0 * epsilon),
                    2e-8,
                );
                near(tangent.get(row, column), tangent.get(column, row), 2e-13);
            }
        }
        let work = surface.work_j(&before, &after, field, cx).unwrap();
        near(work, potential(&before) - potential(&after), 2e-13);
        near(
            surface.work_j(&after, &before, field, cx).unwrap(),
            -work,
            2e-13,
        );
        let middle: [f64; 12] = std::array::from_fn(|i| 0.4 * before[i] + 0.6 * after[i]);
        near(
            surface.work_j(&before, &middle, field, cx).unwrap()
                + surface.work_j(&middle, &after, field, cx).unwrap(),
            work,
            2e-13,
        );
        for factor in [0.0, 0.25, 2.0] {
            let scaled = TetPressure::Stratified {
                profile: &profile,
                load_factor: factor,
            };
            let load = surface.force_n(&after, scaled, cx).unwrap();
            for i in 0..12 {
                near(load[i], factor * force[i], 2e-13);
            }
            near(
                surface.work_j(&before, &after, scaled, cx).unwrap(),
                factor * work,
                2e-13,
            );
        }
    });
}

#[test]
fn stratified_native_discrete_gradient_closes_work_and_endpoint_tangent() {
    with_cx(|_, cx| {
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 3.0,
            },
            5.0,
        )
        .unwrap();
        let problem = HyperTetProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            material: HyperelasticLawField::Uniform(&material),
            prescribed_m: &[],
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings::default(),
        };
        let prepared = problem.prepare(cx).unwrap();
        let surface = surface(cx);
        let profile = profile(cx);
        let (before, after) = endpoints();
        let force = prepared
            .energy_momentum_force_with_pressure(&before, &after, &surface, profile.loading(), cx)
            .unwrap();
        let reverse = prepared
            .energy_momentum_force_with_pressure(&after, &before, &surface, profile.loading(), cx)
            .unwrap();
        let work = (0..12)
            .map(|i| force[i] * (after[i] - before[i]))
            .sum::<f64>();
        let energy = problem.stored_energy_j(&after, cx).unwrap()
            - problem.stored_energy_j(&before, cx).unwrap();
        near(work, energy + potential(&after) - potential(&before), 3e-12);
        for i in 0..12 {
            near(force[i], reverse[i], 3e-12);
        }
        let tangent = prepared
            .energy_momentum_tangent_with_pressure(&before, &after, &surface, profile.loading(), cx)
            .unwrap();
        let epsilon = 2e-6;
        for column in 0..12 {
            let mut plus = after;
            let mut minus = after;
            plus[column] += epsilon;
            minus[column] -= epsilon;
            let forward = prepared
                .energy_momentum_force_with_pressure(
                    &before,
                    &plus,
                    &surface,
                    profile.loading(),
                    cx,
                )
                .unwrap();
            let backward = prepared
                .energy_momentum_force_with_pressure(
                    &before,
                    &minus,
                    &surface,
                    profile.loading(),
                    cx,
                )
                .unwrap();
            for row in 0..12 {
                near(
                    tangent.get(row, column),
                    (forward[row] - backward[row]) / (2.0 * epsilon),
                    3e-7,
                );
            }
        }
    });
}

#[test]
fn stratified_single_layer_equal_density_and_zero_gravity_limits() {
    with_cx(|_, cx| {
        let surface = surface(cx);
        let density = 1.3;
        let layer = [HydrostaticDensityLayer {
            top_depth_m: 0.0,
            density_kg_m3: density,
        }];
        let field = StratifiedHydrostaticPressure::new(ORIGIN, GRAVITY, &layer, 1, cx).unwrap();
        let equal = StratifiedHydrostaticPressure::new(
            ORIGIN,
            GRAVITY,
            &[
                layer[0],
                HydrostaticDensityLayer {
                    top_depth_m: 0.4,
                    density_kg_m3: density,
                },
            ],
            2,
            cx,
        )
        .unwrap();
        let uniform = TetPressure::FreeSurface {
            origin_m: ORIGIN,
            pressure_pa: 0.0,
            gradient_pa_m: GRAVITY.map(|g| density * g),
        };
        let (before, after) = endpoints();
        let a = surface.force_n(&after, field.loading(), cx).unwrap();
        let b = surface.force_n(&after, uniform, cx).unwrap();
        assert_eq!(a, surface.force_n(&after, equal.loading(), cx).unwrap());
        for i in 0..12 {
            near(a[i], b[i], 2e-14);
        }
        let a = surface.force_tangent(&after, field.loading(), cx).unwrap();
        let b = surface.force_tangent(&after, uniform, cx).unwrap();
        for i in 0..12 {
            for j in 0..12 {
                near(a.get(i, j), b.get(i, j), 2e-14);
            }
        }
        near(
            surface
                .work_j(&before, &after, field.loading(), cx)
                .unwrap(),
            surface.work_j(&before, &after, uniform, cx).unwrap(),
            2e-14,
        );
        let zero = StratifiedHydrostaticPressure::new(ORIGIN, [0.0; 3], &LAYERS, 3, cx).unwrap();
        assert_eq!(
            surface.force_n(&after, zero.loading(), cx).unwrap(),
            vec![0.0; 12]
        );
        assert_eq!(
            surface.work_j(&before, &after, zero.loading(), cx).unwrap(),
            0.0
        );
    });
}

#[test]
fn stratified_native_equilibrium_continuation_reaches_manufactured_state() {
    with_cx(|_, cx| {
        let (mu, lambda) = (2.0, 3.0);
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let x = [
            [0.0; 3],
            [0.92, 0.0, 0.0],
            [0.05, 0.95, 0.0],
            [0.02, -0.03, 0.94],
        ];
        let u: [f64; 12] = std::array::from_fn(|d| x[d / 3][d % 3] - NODES[d / 3][d % 3]);
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let free = [3, 6, 7, 9, 10, 11];
        let jacobian = dot(x[1], cross(x[2], x[3]));
        let coefficient = (lambda * jacobian.ln() - mu) / jacobian;
        let cofactors = [cross(x[2], x[3]), cross(x[3], x[1]), cross(x[1], x[2])];
        let mut internal = [0.0; 12];
        for column in 0..3 {
            for axis in 0..3 {
                let value =
                    (mu * x[column + 1][axis] + coefficient * cofactors[column][axis]) / 6.0;
                internal[3 * (column + 1) + axis] = value;
                internal[axis] -= value;
            }
        }
        let surface = surface(cx);
        let profile = StratifiedHydrostaticPressure::new(
            [0.1, -0.2, 0.3],
            [0.013, -0.019, -0.031],
            &LAYERS,
            3,
            cx,
        )
        .unwrap();
        // Surface-force accuracy is independently checked above. Here the
        // imposed equilibrium exercises native load scaling and globalization.
        let applied = surface.force_n(&u, profile.loading(), cx).unwrap();
        let dead: Vec<_> = (0..12)
            .map(|dof| {
                (
                    dof,
                    if free.contains(&dof) {
                        internal[dof] - applied[dof]
                    } else {
                        0.13
                    },
                )
            })
            .collect();
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
            .solve_with_pressure_field(&surface, profile.loading(), cx)
            .unwrap();
        for node in 0..4 {
            for axis in 0..3 {
                near(
                    solution.displacement_m[node][axis],
                    u[3 * node + axis],
                    3e-9,
                );
            }
        }
        for (dof, reaction) in solution.reactions_n {
            near(reaction, internal[dof] - applied[dof] - dead[dof].1, 3e-9);
        }
    });
}

#[test]
fn stratified_admission_budgets_and_cancellation_are_explicit() {
    with_cx(|gate, cx| {
        assert!(StratifiedHydrostaticPressure::new(ORIGIN, GRAVITY, &[], 3, cx).is_err());
        assert!(StratifiedHydrostaticPressure::new(ORIGIN, GRAVITY, &LAYERS, 2, cx).is_err());
        for bad in [
            [
                HydrostaticDensityLayer {
                    top_depth_m: 0.1,
                    ..LAYERS[0]
                },
                LAYERS[1],
                LAYERS[2],
            ],
            [
                LAYERS[0],
                HydrostaticDensityLayer {
                    top_depth_m: 0.0,
                    ..LAYERS[1]
                },
                LAYERS[2],
            ],
            [
                LAYERS[0],
                HydrostaticDensityLayer {
                    density_kg_m3: 0.0,
                    ..LAYERS[1]
                },
                LAYERS[2],
            ],
            [
                LAYERS[0],
                LAYERS[1],
                HydrostaticDensityLayer {
                    top_depth_m: f64::INFINITY,
                    ..LAYERS[2]
                },
            ],
        ] {
            assert!(matches!(
                StratifiedHydrostaticPressure::new(ORIGIN, GRAVITY, &bad, 3, cx),
                Err(HyperTetError::InvalidInput { .. })
            ));
        }
        assert!(
            StratifiedHydrostaticPressure::new(ORIGIN, [f64::NAN, 0.0, 0.0], &LAYERS, 3, cx)
                .is_err()
        );
        let surface = surface(cx);
        let profile = profile(cx);
        let (before, after) = endpoints();
        for factor in [-1.0, f64::INFINITY, f64::NAN] {
            let bad = TetPressure::Stratified {
                profile: &profile,
                load_factor: factor,
            };
            assert!(surface.force_n(&after, bad, cx).is_err());
            assert!(surface.force_tangent(&after, bad, cx).is_err());
            assert!(surface.work_j(&before, &after, bad, cx).is_err());
        }
        gate.request();
        assert!(matches!(
            surface.force_n(&after, profile.loading(), cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            surface.work_j(&before, &after, profile.loading(), cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            profile.pressure_pa(ORIGIN, cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            StratifiedHydrostaticPressure::new(ORIGIN, GRAVITY, &LAYERS, 3, cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}
