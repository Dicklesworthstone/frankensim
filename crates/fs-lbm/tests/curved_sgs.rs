//! Analytic closure/geometry and native transport checks, not turbulent-flow validation.

use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, CollisionModel3, CurvedWall3, CurvedWallLink3, E3, FaceBoundary3,
    TILE, equilibrium3,
};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse, eddy_viscosity};

const ORIGIN: [f64; 3] = [0.7, -1.3, 2.1];
const SHEAR: f64 = 0.004;
const LOWER: f64 = 0.25;
const UPPER: f64 = 0.8;
const MODELS: [CollisionModel3; 5] = [
    CollisionModel3::Bgk { tau: 0.83 },
    CollisionModel3::Trt {
        tau: 0.83,
        magic: 3.0 / 16.0,
    },
    CollisionModel3::CentralMoment {
        second_order_rate: 1.2,
        higher_order_rate: 1.1,
    },
    CollisionModel3::CentralMomentBulk {
        shear_rate: 1.2,
        bulk_rate: 1.05,
        higher_order_rate: 1.1,
    },
    CollisionModel3::ReducedCumulant {
        second_order_rate: 1.2,
        third_order_rate: 1.1,
        fourth_order_rate: 1.05,
    },
];

fn settings(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.4,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, difference {:e}",
        actual - expected,
    );
}

fn native_index([nx, ny, _]: [usize; 3], [x, y, z]: [usize; 3]) -> usize {
    let tile = (z / TILE * (ny / TILE) + y / TILE) * (nx / TILE) + x / TILE;
    let lane = (z % TILE * TILE + y % TILE) * TILE + x % TILE;
    tile * TILE.pow(3) + lane
}

fn bits(grid: &BoundaryGrid3) -> Vec<u64> {
    let [nx, ny, nz] = grid.dimensions();
    let mut result = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if !grid.is_solid(x, y, z) {
                    result.extend(grid.populations(x, y, z).map(f64::to_bits));
                }
            }
        }
    }
    result
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn inventory(grid: &BoundaryGrid3) -> (f64, [f64; 3], [f64; 3]) {
    let [nx, ny, nz] = grid.dimensions();
    let (mut mass, mut momentum, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                let arm = [
                    x as f64 + 0.5 - ORIGIN[0],
                    y as f64 + 0.5 - ORIGIN[1],
                    z as f64 + 0.5 - ORIGIN[2],
                ];
                for (q, f) in grid.populations(x, y, z).into_iter().enumerate() {
                    let e = E3[q];
                    let p = [f * f64::from(e.0), f * f64::from(e.1), f * f64::from(e.2)];
                    let torque = cross(arm, p);
                    mass += f;
                    for a in 0..3 {
                        momentum[a] += p[a];
                        angular[a] += torque[a];
                    }
                }
            }
        }
    }
    (mass, momentum, angular)
}

fn channel(
    axis: usize,
    quadratic: f64,
    forced: bool,
    model: CollisionModel3,
) -> (BoundaryGrid3, CurvedWall3, Vec<[f64; 3]>) {
    let component = (axis + 1) % 3;
    let profile = |s: f64| 0.01 + SHEAR * s + quadratic * s * s;
    let length = 7.0 + LOWER + UPPER;
    let mut faces = [FaceBoundary3::Periodic; 6];
    let mut bottom = [0.0; 3];
    let mut top = [0.0; 3];
    bottom[component] = profile(0.0);
    top[component] = profile(length);
    faces[2 * axis] = FaceBoundary3::Wall { velocity: bottom };
    faces[2 * axis + 1] = FaceBoundary3::Wall { velocity: top };
    let mut grid =
        BoundaryGrid3::with_collision_model(8, 8, 8, model, [0.0; 3], BoundarySpec3::new(faces));
    if forced {
        grid.set_force_field(|p| [1e-5 * (p[1] + 1) as f64, -2e-6, 3e-6]);
    }
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let p = [x, y, z];
                let mut physical = [0.0; 3];
                physical[component] = profile(p[axis] as f64 + LOWER);
                let force = grid.force_density(x, y, z);
                let raw = core::array::from_fn(|a| physical[a] - 0.5 * force[a]);
                let mut f = equilibrium3(1.0, raw);
                // Non-equilibrium shear with zero mass and momentum change:
                // the native collision must actually use the prepared rates.
                f[7] += 2e-5;
                f[8] += 2e-5;
                f[9] -= 2e-5;
                f[10] -= 2e-5;
                grid.set_populations(x, y, z, &f);
            }
        }
    }
    let links: Vec<_> = grid
        .boundary_links()
        .into_iter()
        .map(|link| CurvedWallLink3 {
            fraction: if link.cell[axis] == 0 { LOWER } else { UPPER },
            link,
        })
        .collect();
    let wall = CurvedWall3::new(&grid, &links, &CancelGate::new_clock_free()).unwrap();
    let velocities = wall
        .links()
        .map(|link| {
            if link.link.cell[axis] == 0 {
                bottom
            } else {
                top
            }
        })
        .collect();
    (grid, wall, velocities)
}

#[test]
fn displaced_wall_profiles_use_actual_distances_and_the_guo_velocity_on_every_axis() {
    let gate = CancelGate::new_clock_free();
    let config = settings(SgsModel::Smagorinsky { cs: 0.2 });
    for axis in 0..3 {
        for forced in [false, true] {
            for quadratic in [0.0, 0.00025] {
                let (grid, wall, velocities) = channel(axis, quadratic, forced, MODELS[0]);
                let before = bits(&grid);
                let mut workspace = grid.sgs_workspace().unwrap();
                let allocation = workspace.collisions().as_ptr();
                let report = grid
                    .prepare_curved_sgs(&wall, &velocities, &gate, config, &mut workspace)
                    .unwrap();
                let mut expected_power = 0.0;
                for z in 0..8 {
                    for y in 0..8 {
                        for x in 0..8 {
                            let p = [x, y, z];
                            let shear = SHEAR + 2.0 * quadratic * (p[axis] as f64 + LOWER);
                            let eddy = (0.2_f64 * 1.4).powi(2) * shear;
                            let collision = workspace.collisions()[native_index([8; 3], p)];
                            near(
                                collision.kinematic_viscosity() - grid.viscosity(),
                                eddy,
                                5e-16,
                            );
                            expected_power += eddy * shear * shear;
                        }
                    }
                }
                near(report.sgs_dissipation_lu, expected_power, 3e-17);
                assert_eq!(report.fluid_cells, 512);
                assert_eq!(allocation, workspace.collisions().as_ptr());
                for closure in [SgsModel::Wale { cw: 0.5 }, SgsModel::Vreman { cv: 0.07 }] {
                    let laminar = grid
                        .prepare_curved_sgs(
                            &wall,
                            &velocities,
                            &gate,
                            settings(closure),
                            &mut workspace,
                        )
                        .unwrap();
                    assert!(laminar.maximum_eddy_viscosity_lu < 1e-15);
                }
                assert_eq!(bits(&grid), before);
            }
        }
    }
}

#[test]
fn unequal_wall_distances_use_the_center_value_in_a_two_wall_stencil() {
    let gate = CancelGate::new_clock_free();
    let mut grid = BoundaryGrid3::new(
        4,
        4,
        4,
        0.83,
        [0.0; 3],
        BoundarySpec3::new([FaceBoundary3::Periodic; 6]),
    );
    grid.voxelize_sdf(|p| if p[1] == 2.5 { 1.0 } else { -1.0 });
    for z in 0..4 {
        for x in 0..4 {
            grid.set_populations(x, 2, z, &equilibrium3(1.0, [0.01, 0.0, 0.0]));
        }
    }
    let links: Vec<_> = grid
        .boundary_links()
        .into_iter()
        .map(|link| CurvedWallLink3 {
            fraction: if E3[link.direction].1 > 0 { 0.7 } else { 0.9 },
            link,
        })
        .collect();
    let wall = CurvedWall3::new(&grid, &links, &gate).unwrap();
    let velocities: Vec<_> = wall
        .links()
        .map(|link| {
            let y = link.intersection()[1] - 2.5;
            [0.01 + SHEAR * y + 0.001 * y * y, 0.0, 0.0]
        })
        .collect();
    let mut workspace = grid.sgs_workspace().unwrap();
    let report = grid
        .prepare_curved_sgs(
            &wall,
            &velocities,
            &gate,
            settings(SgsModel::Smagorinsky { cs: 0.2 }),
            &mut workspace,
        )
        .unwrap();
    near(
        report.maximum_eddy_viscosity_lu,
        (0.2_f64 * 1.4).powi(2) * SHEAR,
        3e-16,
    );
    assert_eq!(report.fluid_cells, 16);
}

#[test]
fn rotating_wall_samples_retain_the_full_gradient_for_wale_and_vreman() {
    let gate = CancelGate::new_clock_free();
    let omega = [0.001, -0.002, 0.003];
    let velocity = |p: [f64; 3]| {
        let rotation = cross(omega, p.map(|v| v - 4.0));
        [rotation[0] + 0.01, rotation[1] - 0.005, rotation[2] + 0.003]
    };
    let mut grid = BoundaryGrid3::new(
        8,
        8,
        8,
        0.83,
        [0.0; 3],
        BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]),
    );
    grid.set_force_field(|p| [1e-5 * (p[1] + 1) as f64, -2e-6, 3e-6]);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let u = velocity([x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5]);
                let force = grid.force_density(x, y, z);
                grid.set_populations(
                    x,
                    y,
                    z,
                    &equilibrium3(1.0, core::array::from_fn(|a| u[a] - 0.5 * force[a])),
                );
            }
        }
    }
    let links: Vec<_> = grid
        .boundary_links()
        .into_iter()
        .map(|link| CurvedWallLink3 {
            link,
            fraction: 0.8,
        })
        .collect();
    let wall = CurvedWall3::new(&grid, &links, &gate).unwrap();
    let velocities: Vec<_> = wall
        .links()
        .map(|link| velocity(link.intersection()))
        .collect();
    let gradient = [
        [0.0, -omega[2], omega[1]],
        [omega[2], 0.0, -omega[0]],
        [-omega[1], omega[0], 0.0],
    ];
    let mut workspace = grid.sgs_workspace().unwrap();
    for model in [SgsModel::Wale { cw: 0.5 }, SgsModel::Vreman { cv: 0.07 }] {
        let expected = eddy_viscosity(model, 1.4, gradient).unwrap();
        assert!(expected > 0.0);
        let report = grid
            .prepare_curved_sgs(&wall, &velocities, &gate, settings(model), &mut workspace)
            .unwrap();
        for collision in workspace.collisions() {
            near(
                collision.kinematic_viscosity() - grid.viscosity(),
                expected,
                5e-16,
            );
        }
        assert!(report.sgs_dissipation_lu < 1e-25);
    }
}

#[test]
fn modeled_collision_uses_the_native_curved_transport_and_complete_load_account() {
    let gate = CancelGate::new_clock_free();
    for model in MODELS {
        let (mut actual, wall, velocities) = channel(1, 0.00025, true, model);
        let (mut reference, _, _) = channel(1, 0.00025, true, model);
        let (mut molecular, _, _) = channel(1, 0.00025, true, model);
        let config = settings(SgsModel::Smagorinsky { cs: 0.2 });
        let before = inventory(&actual);
        let mut workspace = actual.sgs_workspace().unwrap();
        let expected_report = actual
            .prepare_curved_sgs(&wall, &velocities, &gate, config, &mut workspace)
            .unwrap();
        assert!(expected_report.sgs_dissipation_lu > 0.0);
        let mut tau = Vec::new();
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..8 {
                    let collision = workspace.collisions()[native_index([8; 3], [x, y, z])];
                    tau.push(match collision {
                        CollisionModel3::Bgk { tau } | CollisionModel3::Trt { tau, .. } => tau,
                        CollisionModel3::CentralMoment {
                            second_order_rate, ..
                        }
                        | CollisionModel3::ReducedCumulant {
                            second_order_rate, ..
                        } => 1.0 / second_order_rate,
                        CollisionModel3::CentralMomentBulk { shear_rate, .. } => 1.0 / shear_rate,
                    });
                }
            }
        }
        let mut sampled = velocities.iter();
        let expected_loads = reference
            .step_with_curved_walls_and_relaxation(
                &wall,
                |_| *sampled.next().unwrap(),
                ORIGIN,
                &tau,
                &gate,
            )
            .unwrap();
        let (loads, report) = actual
            .step_with_curved_walls_and_sgs(
                &wall,
                &velocities,
                ORIGIN,
                &gate,
                config,
                &mut workspace,
            )
            .unwrap();
        assert_eq!(report, expected_report);
        assert_eq!(actual.collision_model(), model);
        for (a, b) in bits(&actual).into_iter().zip(bits(&reference)) {
            near(f64::from_bits(a), f64::from_bits(b), 4e-16);
        }
        for i in 0..3 {
            near(
                loads.curved.wall_impulse[i],
                expected_loads.curved.wall_impulse[i],
                2e-12,
            );
            near(
                loads.curved.wall_angular_impulse[i],
                expected_loads.curved.wall_angular_impulse[i],
                2e-11,
            );
        }
        near(
            loads.curved.wall_work,
            expected_loads.curved.wall_work,
            2e-13,
        );
        let mut sampled = velocities.iter();
        let _ = molecular
            .step_with_curved_walls(&wall, |_| *sampled.next().unwrap(), ORIGIN, &gate)
            .unwrap();
        assert_ne!(bits(&actual), bits(&molecular));
        let after = inventory(&actual);
        near(
            after.0 - before.0,
            loads.curved.fluid_mass_change
                + loads.other_walls.fluid_mass_change
                + loads.open.mass_change,
            2e-10,
        );
        for i in 0..3 {
            near(
                after.1[i] - before.1[i],
                loads.curved.fluid_population_impulse[i]
                    + loads.other_walls.fluid_population_impulse[i]
                    + loads.open.impulse[i]
                    + loads.body_force_impulse[i],
                3e-12,
            );
            near(
                after.2[i] - before.2[i],
                loads.curved.fluid_population_angular_impulse[i]
                    + loads.other_walls.fluid_population_angular_impulse[i]
                    + loads.open.angular_impulse[i]
                    + loads.periodic.angular_momentum_change[i]
                    + loads.body_force_angular_impulse[i],
                3e-11,
            );
        }
    }
}

#[test]
fn zero_coefficient_preserves_curved_step_bits_and_refusals_leave_replay_intact() {
    let gate = CancelGate::new_clock_free();
    for model in MODELS {
        let (mut actual, wall, velocities) = channel(1, 0.00025, true, model);
        let (mut reference, _, _) = channel(1, 0.00025, true, model);
        let mut workspace = actual.sgs_workspace().unwrap();
        let initial = bits(&actual);
        let nonzero = settings(SgsModel::Smagorinsky { cs: 0.2 });
        let mut limited = nonzero;
        limited.maximum_eddy_to_molecular = 0.0;
        assert!(
            actual
                .step_with_curved_walls_and_sgs(
                    &wall,
                    &velocities,
                    ORIGIN,
                    &gate,
                    limited,
                    &mut workspace,
                )
                .is_err()
        );
        let mut invalid = velocities.clone();
        *invalid.last_mut().unwrap() = [f64::NAN; 3];
        assert!(
            actual
                .step_with_curved_walls_and_sgs(
                    &wall,
                    &invalid,
                    ORIGIN,
                    &gate,
                    nonzero,
                    &mut workspace,
                )
                .is_err()
        );
        assert!(
            actual
                .step_with_curved_walls_and_sgs(
                    &wall,
                    &velocities[..velocities.len() - 1],
                    ORIGIN,
                    &gate,
                    nonzero,
                    &mut workspace,
                )
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            actual
                .step_with_curved_walls_and_sgs(
                    &wall,
                    &velocities,
                    ORIGIN,
                    &cancelled,
                    nonzero,
                    &mut workspace,
                )
                .is_err()
        );
        assert_eq!(bits(&actual), initial);
        for closure in [
            SgsModel::Smagorinsky { cs: 0.0 },
            SgsModel::Wale { cw: 0.0 },
            SgsModel::Vreman { cv: 0.0 },
        ] {
            let mut sampled = velocities.iter();
            let expected = reference
                .step_with_curved_walls(&wall, |_| *sampled.next().unwrap(), ORIGIN, &gate)
                .unwrap();
            let (loads, report) = actual
                .step_with_curved_walls_and_sgs(
                    &wall,
                    &velocities,
                    ORIGIN,
                    &gate,
                    settings(closure),
                    &mut workspace,
                )
                .unwrap();
            assert_eq!(loads, expected);
            assert_eq!(report.sgs_dissipation_lu, 0.0);
            assert_eq!(bits(&actual), bits(&reference));
        }
    }
}
