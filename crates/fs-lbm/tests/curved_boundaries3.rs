use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, CollisionModel3, CurvedWall3, CurvedWallLink3, E3, FaceBoundary3,
    OPP3, Q3, W3, collide_cell3, equilibrium3,
};

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, error {:e}",
        (a - b).abs()
    );
}
fn near3(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for i in 0..3 {
        near(a[i], b[i], tolerance);
    }
}
fn closed() -> BoundarySpec3 {
    BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6])
}
fn periodic() -> BoundarySpec3 {
    BoundarySpec3::new([FaceBoundary3::Periodic; 6])
}
fn channel(u: [f64; 3]) -> BoundarySpec3 {
    let p = FaceBoundary3::Periodic;
    BoundarySpec3::new([
        p,
        p,
        FaceBoundary3::Wall { velocity: u },
        FaceBoundary3::Wall { velocity: u },
        p,
        p,
    ])
}
fn snapshots(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let [nx, ny, nz] = grid.dimensions();
    let mut values = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if !grid.is_solid(x, y, z) {
                    values.push(grid.populations(x, y, z));
                }
            }
        }
    }
    values
}
fn inventory(grid: &BoundaryGrid3, origin: [f64; 3]) -> (f64, [f64; 3], [f64; 3]) {
    let [nx, ny, nz] = grid.dimensions();
    let (mut mass, mut linear, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                let r = [
                    x as f64 + 0.5 - origin[0],
                    y as f64 + 0.5 - origin[1],
                    z as f64 + 0.5 - origin[2],
                ];
                for (q, f) in grid.populations(x, y, z).into_iter().enumerate() {
                    let c = E3[q];
                    let p = [f * f64::from(c.0), f * f64::from(c.1), f * f64::from(c.2)];
                    let l = cross(r, p);
                    mass += f;
                    for i in 0..3 {
                        linear[i] += p[i];
                        angular[i] += l[i];
                    }
                }
            }
        }
    }
    (mass, linear, angular)
}
fn seed(grid: &mut BoundaryGrid3) {
    let [nx, ny, nz] = grid.dimensions();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if grid.is_solid(x, y, z) {
                    continue;
                }
                let rho = 0.98 + 0.002 * (x + 2 * y + z) as f64;
                let u = [
                    0.001 * (y as f64 - 2.0),
                    -0.002 * (x as f64 - 3.0),
                    0.0007 * (z as f64 - 4.0),
                ];
                let f = equilibrium3(rho, u);
                grid.set_populations(x, y, z, &f);
            }
        }
    }
}

const MODELS: [CollisionModel3; 4] = [
    CollisionModel3::Bgk { tau: 0.83 },
    CollisionModel3::Trt {
        tau: 0.83,
        magic: 3.0 / 16.0,
    },
    CollisionModel3::CentralMoment {
        second_order_rate: 1.2,
        higher_order_rate: 1.1,
    },
    CollisionModel3::ReducedCumulant {
        second_order_rate: 1.2,
        third_order_rate: 1.1,
        fourth_order_rate: 1.05,
    },
];

#[test]
fn halfway_recovers_original_populations_and_exchange_for_every_collision_family() {
    let u = [0.017, 0.0, -0.009];
    let gate = CancelGate::new();
    for model in MODELS {
        let mut original =
            BoundaryGrid3::with_collision_model(4, 8, 4, model, [1e-5, 0.0, -2e-5], channel(u));
        let mut curved =
            BoundaryGrid3::with_collision_model(4, 8, 4, model, [1e-5, 0.0, -2e-5], channel(u));
        seed(&mut original);
        seed(&mut curved);
        let links = original.boundary_links();
        let binding = CurvedWall3::new(
            &curved,
            &links
                .iter()
                .map(|&link| CurvedWallLink3 {
                    link,
                    fraction: 0.5,
                })
                .collect::<Vec<_>>(),
            &gate,
        )
        .unwrap();
        for _ in 0..3 {
            let a = original.step_with_wall_momentum(&links, [0.7, -1.3, 2.1]);
            let b = curved
                .step_with_curved_walls(&binding, |_| u, [0.7, -1.3, 2.1], &gate)
                .unwrap();
            assert_eq!(a, b.curved);
            assert_eq!(b.other_walls.measured_links, 0);
            for (a, b) in snapshots(&original)
                .into_iter()
                .flatten()
                .zip(snapshots(&curved).into_iter().flatten())
            {
                assert_eq!(a.to_bits(), b.to_bits(), "{model:?}");
            }
        }
    }
}

#[test]
fn both_interpolation_branches_use_the_original_postcollision_population() {
    let gate = CancelGate::new();
    for fraction in [0.25, 0.5, 0.8, 1.0] {
        for normal in 0..3 {
            let force = [1e-5, -3e-5, 2e-5];
            let mut grid = BoundaryGrid3::new(8, 8, 8, 0.9, force, closed());
            seed(&mut grid);
            let link = grid
                .boundary_links()
                .into_iter()
                .find(|link| {
                    link.cell[normal] == 0
                        && (0..3).filter(|&i| i != normal).all(|i| link.cell[i] == 3)
                        && {
                            let c = E3[link.direction];
                            let e = [c.0, c.1, c.2];
                            e[normal] == 1 && e[(normal + 1) % 3] == 1
                        }
                })
                .unwrap();
            let binding =
                CurvedWall3::new(&grid, &[CurvedWallLink3 { link, fraction }], &gate).unwrap();
            let [x, y, z] = link.cell;
            let q = link.direction;
            let c = E3[q];
            let e = [f64::from(c.0), f64::from(c.1), f64::from(c.2)];
            let post =
                collide_cell3(grid.populations(x, y, z), grid.collision_model(), force).unwrap();
            let far_cell = core::array::from_fn::<_, 3, _>(|i| {
                (link.cell[i] as i32 + [c.0, c.1, c.2][i]) as usize
            });
            let far = collide_cell3(
                grid.populations(far_cell[0], far_cell[1], far_cell[2]),
                grid.collision_model(),
                force,
            )
            .unwrap();
            let u = [0.013, -0.009, 0.007];
            let rho = post.iter().sum::<f64>();
            let correction = 6.0 * W3[q] * rho * (e[0] * u[0] + e[1] * u[1] + e[2] * u[2]);
            let expected = if fraction < 0.5 {
                2.0 * fraction * post[OPP3[q]] + (1.0 - 2.0 * fraction) * far[OPP3[q]] + correction
            } else {
                (post[OPP3[q]] + correction) / (2.0 * fraction)
                    + (1.0 - 1.0 / (2.0 * fraction)) * post[q]
            };
            let before = inventory(&grid, [0.0; 3]);
            let report = grid
                .step_with_curved_walls(&binding, |_| u, [0.0; 3], &gate)
                .unwrap();
            near(grid.populations(x, y, z)[q], expected, 3e-17);
            let incoming = grid.populations(x, y, z)[q];
            let outgoing = post[OPP3[q]];
            let j = core::array::from_fn(|i| (-e[i] - u[i]) * outgoing - (e[i] - u[i]) * incoming);
            near3(report.curved.wall_impulse, j, 3e-17);
            let point = CurvedWallLink3 { link, fraction }.intersection();
            near3(report.curved.wall_angular_impulse, cross(point, j), 3e-16);
            near(
                report.curved.wall_work,
                (0..3).map(|i| j[i] * u[i]).sum(),
                3e-18,
            );
            let after = inventory(&grid, [0.0; 3]);
            near(
                after.0 - before.0,
                report.curved.fluid_mass_change + report.other_walls.fluid_mass_change,
                3e-10,
            );
            for i in 0..3 {
                near(
                    after.1[i] - before.1[i],
                    report.curved.fluid_population_impulse[i]
                        + report.other_walls.fluid_population_impulse[i]
                        + report.body_force_impulse[i],
                    2e-12,
                );
                near(
                    after.2[i] - before.2[i],
                    report.curved.fluid_population_angular_impulse[i]
                        + report.other_walls.fluid_population_angular_impulse[i]
                        + report.body_force_angular_impulse[i],
                    2e-11,
                );
            }
        }
    }
}

#[test]
fn generated_sphere_retains_actual_torque_work_mass_and_periodic_transport() {
    let gate = CancelGate::new();
    let center = [4.0; 3];
    let radius = 1.67_f64;
    let sdf =
        |p: [f64; 3]| (0..3).map(|i| (p[i] - center[i]).powi(2)).sum::<f64>() - radius * radius;
    let origin = [-1.7, 2.3, 0.4];
    let omega = [0.003, -0.002, 0.005];
    for model in MODELS {
        let mut grid =
            BoundaryGrid3::with_collision_model(8, 8, 8, model, [2e-6, -1e-6, 3e-6], periodic());
        grid.voxelize_sdf(sdf);
        let wall = CurvedWall3::from_level_set(&grid, sdf, 48, 1e-12, &gate).unwrap();
        assert!(wall.links().any(|x| x.fraction < 0.5));
        assert!(wall.links().any(|x| x.fraction > 0.5));
        for link in wall.links() {
            near(sdf(link.intersection()), 0.0, 8e-12);
        }
        seed(&mut grid);
        for _ in 0..4 {
            let before = inventory(&grid, origin);
            let result = grid
                .step_with_curved_walls(
                    &wall,
                    |link| {
                        cross(
                            omega,
                            core::array::from_fn(|i| link.intersection()[i] - center[i]),
                        )
                    },
                    origin,
                    &gate,
                )
                .unwrap();
            let after = inventory(&grid, origin);
            assert_eq!(result.other_walls.measured_links, 0);
            assert!(result.periodic.face_crossings > 0);
            near(after.0 - before.0, result.curved.fluid_mass_change, 3e-10);
            let center_arm = core::array::from_fn(|i| center[i] - origin[i]);
            let orbital = cross(center_arm, result.curved.wall_impulse);
            let torque_center = core::array::from_fn::<_, 3, _>(|i| {
                result.curved.wall_angular_impulse[i] - orbital[i]
            });
            near(
                result.curved.wall_work,
                (0..3).map(|i| omega[i] * torque_center[i]).sum(),
                3e-15,
            );
            for i in 0..3 {
                near(
                    after.1[i] - before.1[i],
                    result.curved.fluid_population_impulse[i] + result.body_force_impulse[i],
                    3e-12,
                );
                near(
                    after.2[i] - before.2[i],
                    result.curved.fluid_population_angular_impulse[i]
                        + result.body_force_angular_impulse[i]
                        + result.periodic.angular_momentum_change[i],
                    3e-11,
                );
                near(
                    result.curved.wall_impulse[i] + result.curved.fluid_population_impulse[i],
                    result.curved.wall_velocity_mass_impulse[i],
                    3e-14,
                );
                near(
                    result.curved.wall_angular_impulse[i]
                        + result.curved.fluid_population_angular_impulse[i],
                    result.curved.wall_velocity_mass_angular_impulse[i],
                    2e-13,
                );
            }
        }
    }
}

#[test]
fn fractional_wall_couette_reaches_the_physical_surface_velocities() {
    let gate = CancelGate::new();
    let (lower, upper) = (1.2, 11.3);
    let sdf = |p: [f64; 3]| (p[1] - lower).min(upper - p[1]);
    let mut grid = BoundaryGrid3::new(4, 12, 4, 0.8, [0.0; 3], channel([0.0; 3]));
    grid.voxelize_sdf(sdf);
    let wall = CurvedWall3::from_level_set(&grid, sdf, 48, 1e-12, &gate).unwrap();
    for _ in 0..2000 {
        let _ = grid
            .step_with_curved_walls(
                &wall,
                |link| {
                    if link.intersection()[1] < 2.0 {
                        [-0.01, 0.0, 0.0]
                    } else {
                        [0.03, 0.0, 0.0]
                    }
                },
                [0.0; 3],
                &gate,
            )
            .unwrap();
    }
    for y in 1..11 {
        let expected = -0.01 + 0.04 * (y as f64 + 0.5 - lower) / (upper - lower);
        let u = grid.velocity(2, y, 2);
        near(u[0], expected, 2e-10);
        near(u[1], 0.0, 2e-13);
        near(u[2], 0.0, 2e-13);
    }
}

fn circular_couette_error(scale: usize) -> f64 {
    let n = 12 * scale;
    let center = 6.0 * scale as f64;
    let inner = 2.3 * scale as f64;
    let outer = 5.2 * scale as f64;
    let omega = 0.006 / scale as f64;
    let sdf = |p: [f64; 3]| {
        let r = ((p[0] - center).powi(2) + (p[1] - center).powi(2)).sqrt();
        (r - inner).min(outer - r)
    };
    let p = FaceBoundary3::Periodic;
    let w = FaceBoundary3::stationary_wall();
    let mut grid = BoundaryGrid3::new(
        n,
        n,
        4,
        1.0,
        [0.0; 3],
        BoundarySpec3::new([w, w, w, w, p, p]),
    );
    grid.voxelize_sdf(sdf);
    let gate = CancelGate::new();
    let wall = CurvedWall3::from_level_set(&grid, sdf, 48, 1e-12, &gate).unwrap();
    for _ in 0..(600 * scale * scale) {
        let _ = grid
            .step_with_curved_walls(
                &wall,
                |link| {
                    let p = link.intersection();
                    let r = ((p[0] - center).powi(2) + (p[1] - center).powi(2)).sqrt();
                    if r < (inner + outer) * 0.5 {
                        [-omega * (p[1] - center), omega * (p[0] - center), 0.0]
                    } else {
                        [0.0; 3]
                    }
                },
                [center, center, 0.0],
                &gate,
            )
            .unwrap();
    }
    let (mut error, mut norm) = (0.0, 0.0);
    for y in 0..n {
        for x in 0..n {
            if grid.is_solid(x, y, 1) {
                continue;
            }
            let dx = x as f64 + 0.5 - center;
            let dy = y as f64 + 0.5 - center;
            let r2 = dx * dx + dy * dy;
            let angular = omega * inner * inner / (outer * outer - inner * inner)
                * (outer * outer / r2 - 1.0);
            let expected = [-angular * dy, angular * dx, 0.0];
            let u = grid.velocity(x, y, 1);
            for i in 0..3 {
                error += (u[i] - expected[i]).powi(2);
                norm += expected[i] * expected[i];
            }
        }
    }
    (error / norm).sqrt()
}

#[test]
fn curved_circular_couette_refines_toward_the_cylindrical_solution() {
    // Resolve the inner radius by 4.6 and 9.2 cells. The 2.3-cell radius is
    // deliberately outside this accuracy envelope (about 12% relative error).
    let coarse = circular_couette_error(2);
    let fine = circular_couette_error(4);
    eprintln!("curved Couette relative errors: coarse={coarse:e}, fine={fine:e}");
    assert!(coarse < 0.02, "coarse curved-flow error {coarse:e}");
    assert!(fine < 0.01, "fine curved-flow error {fine:e}");
    assert!(fine < 0.6 * coarse, "refinement {coarse:e} -> {fine:e}");
}

#[test]
fn geometry_donor_and_late_failure_refuse_without_fluid_publication() {
    let gate = CancelGate::new();
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.8, [0.0; 3], closed());
    let link = grid.boundary_links()[0];
    assert!(
        CurvedWall3::new(
            &grid,
            &[CurvedWallLink3 {
                link,
                fraction: 0.0
            }],
            &gate
        )
        .is_err()
    );
    let same = CurvedWallLink3 {
        link,
        fraction: 0.5,
    };
    assert!(CurvedWall3::new(&grid, &[same, same], &gate).is_err());
    let geometry = grid
        .boundary_links()
        .into_iter()
        .map(|link| CurvedWallLink3 {
            link,
            fraction: 0.5,
        })
        .collect::<Vec<_>>();
    let wall = CurvedWall3::new(&grid, &geometry, &gate).unwrap();
    let initial = snapshots(&grid);
    let mut calls = 0;
    assert!(
        grid.step_with_curved_walls(
            &wall,
            |_| {
                calls += 1;
                if calls == geometry.len() {
                    [f64::NAN, 0.0, 0.0]
                } else {
                    [0.0; 3]
                }
            },
            [0.0; 3],
            &gate
        )
        .is_err()
    );
    assert_eq!(snapshots(&grid), initial);
    let late_gate = CancelGate::new();
    let mut calls = 0;
    assert!(
        grid.step_with_curved_walls(
            &wall,
            |_| {
                calls += 1;
                if calls == geometry.len() {
                    late_gate.request();
                }
                [0.0; 3]
            },
            [0.0; 3],
            &late_gate
        )
        .is_err()
    );
    assert_eq!(snapshots(&grid), initial);
    let _ = grid
        .step_with_curved_walls(&wall, |_| [0.0; 3], [0.0; 3], &gate)
        .unwrap();
    let mut trapped = BoundaryGrid3::new(4, 4, 4, 0.8, [0.0; 3], closed());
    trapped.voxelize_sdf(|p| if p == [1.5; 3] { 1.0 } else { -1.0 });
    let selected = trapped.boundary_links()[0];
    assert!(
        CurvedWall3::new(
            &trapped,
            &[CurvedWallLink3 {
                link: selected,
                fraction: 0.25
            }],
            &gate
        )
        .is_err()
    );
    assert!(CurvedWall3::from_level_set(&grid, |_| -1.0, 48, 1e-12, &gate).is_err());
    assert!(CurvedWall3::from_level_set(&grid, |_| f64::NAN, 48, 1e-12, &gate).is_err());
}

#[test]
fn open_halfway_retains_native_populations_for_all_axes_and_collision_models() {
    let gate = CancelGate::new();
    for model in MODELS {
        for axis in 0..3 {
            for pressure_inlet in [false, true] {
                let mut faces = [FaceBoundary3::stationary_wall(); 6];
                let mut inlet = [0.0; 3];
                inlet[axis] = 0.014;
                faces[2 * axis] = if pressure_inlet {
                    FaceBoundary3::Pressure { density: 1.003 }
                } else {
                    FaceBoundary3::Velocity { velocity: inlet }
                };
                faces[2 * axis + 1] = FaceBoundary3::Pressure { density: 0.999 };
                let make = || {
                    let mut grid = BoundaryGrid3::with_collision_model(
                        8,
                        8,
                        8,
                        model,
                        [1e-6, -2e-6, 3e-6],
                        BoundarySpec3::new(faces),
                    );
                    seed(&mut grid);
                    grid
                };
                let (mut original, mut curved) = (make(), make());
                let links = original.boundary_links();
                let wall = CurvedWall3::new(
                    &curved,
                    &links
                        .iter()
                        .map(|&link| CurvedWallLink3 {
                            link,
                            fraction: 0.5,
                        })
                        .collect::<Vec<_>>(),
                    &gate,
                )
                .unwrap();
                for _ in 0..2 {
                    let native = original.step_with_wall_momentum(&links, [0.7, -1.3, 2.1]);
                    let result = curved
                        .step_with_curved_walls(&wall, |_| [0.0; 3], [0.7, -1.3, 2.1], &gate)
                        .unwrap();
                    assert_eq!(native, result.curved);
                    assert_eq!(result.other_walls.measured_links, 0);
                    assert!(result.open.open_links > 0);
                    assert_eq!(result.open.reconstructed_cells, 2 * 8 * 8);
                    for (a, b) in snapshots(&original)
                        .into_iter()
                        .flatten()
                        .zip(snapshots(&curved).into_iter().flatten())
                    {
                        assert_eq!(a.to_bits(), b.to_bits(), "{model:?}, axis {axis}");
                    }
                }
            }
        }
    }
}

#[test]
fn rotating_curved_body_in_open_flow_closes_separate_population_accounts() {
    let gate = CancelGate::new();
    let center = [6.0; 3];
    let origin = [0.7, -1.3, 2.1];
    let omega = [0.002, -0.003, 0.004];
    let sdf = |p: [f64; 3]| {
        ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2) + (p[2] - center[2]).powi(2))
            .sqrt()
            - 2.3
    };
    let p = FaceBoundary3::Periodic;
    let w = FaceBoundary3::stationary_wall();
    for model in MODELS {
        let mut grid = BoundaryGrid3::with_collision_model(
            12,
            12,
            12,
            model,
            [1e-6, -2e-6, 3e-6],
            BoundarySpec3::new([
                FaceBoundary3::Velocity {
                    velocity: [0.014, 0.0, 0.0],
                },
                FaceBoundary3::Pressure { density: 0.999 },
                p,
                p,
                w,
                w,
            ]),
        );
        grid.voxelize_sdf(sdf);
        seed(&mut grid);
        grid.rebind_force_field(|cell| {
            [
                2e-7 * cell[2] as f64,
                -1e-7 * cell[0] as f64,
                3e-7 * cell[1] as f64,
            ]
        });
        let wall = CurvedWall3::from_level_set(&grid, sdf, 48, 1e-12, &gate).unwrap();
        assert!(wall.links().any(|l| l.fraction < 0.5));
        assert!(wall.links().any(|l| l.fraction > 0.5));
        for _ in 0..3 {
            let before = inventory(&grid, origin);
            let result = grid
                .step_with_curved_walls(
                    &wall,
                    |link| {
                        cross(
                            omega,
                            core::array::from_fn(|i| link.intersection()[i] - center[i]),
                        )
                    },
                    origin,
                    &gate,
                )
                .unwrap();
            let after = inventory(&grid, origin);
            assert!(result.open.open_links > 0);
            assert!(result.periodic.face_crossings > 0);
            assert!(result.other_walls.measured_links > 0);
            near(
                after.0 - before.0,
                result.curved.fluid_mass_change
                    + result.other_walls.fluid_mass_change
                    + result.open.mass_change,
                3e-9,
            );
            for i in 0..3 {
                near(
                    after.1[i] - before.1[i],
                    result.curved.fluid_population_impulse[i]
                        + result.other_walls.fluid_population_impulse[i]
                        + result.open.impulse[i]
                        + result.body_force_impulse[i],
                    3e-11,
                );
                near(
                    after.2[i] - before.2[i],
                    result.curved.fluid_population_angular_impulse[i]
                        + result.other_walls.fluid_population_angular_impulse[i]
                        + result.open.angular_impulse[i]
                        + result.periodic.angular_momentum_change[i]
                        + result.body_force_angular_impulse[i],
                    3e-10,
                );
                near(
                    result.curved.wall_impulse[i] + result.curved.fluid_population_impulse[i],
                    result.curved.wall_velocity_mass_impulse[i],
                    3e-13,
                );
            }
            // The reservoir owns the measured input, not a force assigned to
            // the rotating obstacle. The original prescribed moment survives.
            near3(grid.velocity(0, 6, 6), [0.014, 0.0, 0.0], 3e-14);
        }
    }
}

#[test]
fn fractional_open_duct_preserves_wall_corners_and_develops_pressure_driven_flow() {
    let gate = CancelGate::new();
    let (nx, ny) = (96, 16);
    let (lower, upper) = (0.7, ny as f64 - 0.7);
    let sdf = |p: [f64; 3]| (p[1] - lower).min(upper - p[1]);
    let p = FaceBoundary3::Periodic;
    let w = FaceBoundary3::stationary_wall();
    let (inlet, outlet, tau) = (1.001, 0.999, 0.8);
    let mut grid = BoundaryGrid3::new(
        nx,
        ny,
        4,
        tau,
        [0.0; 3],
        BoundarySpec3::new([
            FaceBoundary3::Pressure { density: inlet },
            FaceBoundary3::Pressure { density: outlet },
            w,
            w,
            p,
            p,
        ]),
    );
    grid.voxelize_sdf(sdf);
    let wall = CurvedWall3::from_level_set(&grid, sdf, 48, 1e-12, &gate).unwrap();
    let corner = wall
        .links()
        .find(|link| link.link.cell == [0, 1, 1] && E3[link.link.direction] == (-1, 1, 0))
        .unwrap();
    let [x, y, z] = corner.link.cell;
    let q = corner.link.direction;
    let post = collide_cell3(grid.populations(x, y, z), grid.collision_model(), [0.0; 3]).unwrap();
    let width = 2.0 * corner.fraction;
    let expected = post[OPP3[q]] / width + (width - 1.0) / width * post[q];
    let result = grid
        .step_with_curved_walls(&wall, |_| [0.0; 3], [0.0; 3], &gate)
        .unwrap();
    assert_eq!(grid.populations(x, y, z)[q].to_bits(), expected.to_bits());
    assert!(grid.open_link_mask(x, y, z) & (1 << q) == 0);
    assert_eq!(result.open.reconstructed_cells, 2 * (ny - 2) * 4);
    for _ in 0..9600 {
        let _ = grid
            .step_with_curved_walls(&wall, |_| [0.0; 3], [0.0; 3], &gate)
            .unwrap();
    }
    // On-site pressure locations are x=0.5 and x=nx-0.5; physical walls lie
    // at the level-set intersections, not the voxel halfway planes.
    let viscosity = (tau - 0.5) / 3.0;
    let pressure_gradient = (inlet - outlet) / (3.0 * (nx - 1) as f64);
    let (mut error, mut norm) = (0.0, 0.0);
    for y in 1..ny - 1 {
        let position = y as f64 + 0.5;
        let expected =
            pressure_gradient * (position - lower) * (upper - position) / (2.0 * viscosity);
        let actual = grid.velocity(nx / 2, y, 2)[0];
        error += (actual - expected).powi(2);
        norm += expected * expected;
    }
    let relative = (error / norm).sqrt();
    eprintln!("fractional open-duct velocity error: {relative:e}");
    assert!(
        relative < 0.03,
        "resolved open-duct velocity error {relative:e}"
    );
}

#[test]
fn late_open_reconstruction_overflow_does_not_publish_populations() {
    let gate = CancelGate::new();
    let p = FaceBoundary3::Periodic;
    let mut grid = BoundaryGrid3::new(
        4,
        4,
        4,
        0.8,
        [0.0; 3],
        BoundarySpec3::new([
            FaceBoundary3::Pressure { density: 1e200 },
            FaceBoundary3::Pressure { density: 1.0 },
            p,
            p,
            p,
            p,
        ]),
    );
    let wall = CurvedWall3::new(&grid, &[], &gate).unwrap();
    let before = snapshots(&grid);
    let error = grid
        .step_with_curved_walls(&wall, |_| unreachable!(), [0.0, 1e200, 0.0], &gate)
        .unwrap_err();
    assert!(matches!(
        error,
        fs_lbm::d3q19::BoundaryStepError3::WallMomentum {
            reason: "open-boundary momentum account overflowed"
        }
    ));
    assert_eq!(snapshots(&grid), before);
}
