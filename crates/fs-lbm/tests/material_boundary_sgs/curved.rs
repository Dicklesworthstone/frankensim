use super::*;
use fs_lbm::d3q19::{CurvedWall3, CurvedWallLink3, CurvedWallStep3};

const LOWER: f64 = 0.25;
const UPPER: f64 = 0.8;
const SHEAR: f64 = 0.004;

fn collisions(tau: f64) -> [CollisionModel3; 5] {
    [
        model(tau, 0),
        model(tau, 1),
        model(tau, 2),
        CollisionModel3::CentralMomentBulk {
            shear_rate: 1.0 / tau,
            bulk_rate: 0.93,
            higher_order_rate: 1.27,
        },
        model(tau, 3),
    ]
}

fn channel(
    axis: usize,
    quadratic: f64,
    collision: CollisionModel3,
) -> (BoundaryGrid3, CurvedWall3, Vec<[f64; 3]>) {
    let component = (axis + 1) % 3;
    let profile = |s: f64| 0.01 + SHEAR * s + quadratic * s * s;
    let mut bottom = [0.0; 3];
    let mut top = [0.0; 3];
    bottom[component] = profile(0.0);
    top[component] = profile(3.0 + LOWER + UPPER);
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[2 * axis] = FaceBoundary3::Wall { velocity: bottom };
    faces[2 * axis + 1] = FaceBoundary3::Wall { velocity: top };
    let mut grid = BoundaryGrid3::with_collision_model(
        4,
        4,
        4,
        collision,
        [1e-6, -2e-6, 0.0],
        BoundarySpec3::new(faces),
    );
    grid.set_force_field(|p| [2e-6 * (p[1] + 1) as f64, -3e-6, 1e-6]);
    seed(&mut grid, |p| {
        let mut velocity = [0.0; 3];
        velocity[component] = profile(p[axis] as f64 + LOWER);
        velocity
    });
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

fn inventory(grid: &BoundaryGrid3) -> (f64, [f64; 3], [f64; 3]) {
    let [nx, ny, nz] = grid.dimensions();
    let (mut mass, mut momentum, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let arm = [
                    x as f64 + 0.5 - ORIGIN[0],
                    y as f64 + 0.5 - ORIGIN[1],
                    z as f64 + 0.5 - ORIGIN[2],
                ];
                for (q, f) in grid.populations(x, y, z).into_iter().enumerate() {
                    let e = E3[q];
                    let p = [f * f64::from(e.0), f * f64::from(e.1), f * f64::from(e.2)];
                    mass += f;
                    for axis in 0..3 {
                        let b = (axis + 1) % 3;
                        let c = (axis + 2) % 3;
                        momentum[axis] += p[axis];
                        angular[axis] += arm[b] * p[c] - arm[c] * p[b];
                    }
                }
            }
        }
    }
    (mass, momentum, angular)
}

fn check_exchange(
    before: (f64, [f64; 3], [f64; 3]),
    after: (f64, [f64; 3], [f64; 3]),
    receipt: CurvedWallStep3,
) {
    close(
        after.0 - before.0,
        receipt.curved.fluid_mass_change
            + receipt.other_walls.fluid_mass_change
            + receipt.open.mass_change,
        3e-11,
    );
    for axis in 0..3 {
        close(
            after.1[axis] - before.1[axis],
            receipt.curved.fluid_population_impulse[axis]
                + receipt.other_walls.fluid_population_impulse[axis]
                + receipt.open.impulse[axis]
                + receipt.body_force_impulse[axis],
            3e-12,
        );
        close(
            after.2[axis] - before.2[axis],
            receipt.curved.fluid_population_angular_impulse[axis]
                + receipt.other_walls.fluid_population_angular_impulse[axis]
                + receipt.open.angular_impulse[axis]
                + receipt.periodic.angular_momentum_change[axis]
                + receipt.body_force_angular_impulse[axis],
            3e-11,
        );
        close(
            receipt.curved.wall_impulse[axis] + receipt.curved.fluid_population_impulse[axis],
            receipt.curved.wall_velocity_mass_impulse[axis],
            3e-12,
        );
        close(
            receipt.curved.wall_angular_impulse[axis]
                + receipt.curved.fluid_population_angular_impulse[axis],
            receipt.curved.wall_velocity_mass_angular_impulse[axis],
            3e-11,
        );
    }
}

#[test]
fn displaced_quadratic_walls_supply_the_same_analytic_shear_to_material_and_sgs() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut config = settings(SgsModel::Smagorinsky { cs: 0.2 });
    config.filter_width_lu = 1.4;
    for axis in 0..3 {
        for quadratic in [0.0, 0.00025] {
            for collision in collisions(1.37) {
                let (mut grid, wall, velocities) = channel(axis, quadratic, collision);
                let mut expected = grid.clone();
                let before = inventory(&grid);
                let mut tau = Vec::new();
                let mut molecular_power = 0.0;
                let mut modeled_power = 0.0;
                let mut maximum_ratio = 0.0_f64;
                for z in 0..4 {
                    for y in 0..4 {
                        for x in 0..4 {
                            let cell = [x, y, z];
                            let shear = SHEAR + 2.0 * quadratic * (cell[axis] as f64 + LOWER);
                            let mu = 0.4 + 0.02 * x as f64 + 0.1 * shear / DT;
                            let molecular = mu / 1000.0 * DT / (DX * DX);
                            let eddy = (0.2_f64 * 1.4).powi(2) * shear;
                            tau.push(0.5 + 3.0 * (molecular + eddy));
                            molecular_power += molecular * shear * shear;
                            modeled_power += eddy * shear * shear;
                            maximum_ratio = maximum_ratio.max(eddy / molecular);
                        }
                    }
                }
                let mut velocity = velocities.iter();
                let expected_loads = expected
                    .step_with_curved_walls_and_relaxation(
                        &wall,
                        |_| *velocity.next().unwrap(),
                        ORIGIN,
                        &tau,
                        &gate,
                    )
                    .unwrap();
                let mut work = scaling
                    .boundary_rheology_sgs_workspace(&grid, 64, &gate)
                    .unwrap();
                let bytes = work.workspace_bytes();
                let mut calls = 0;
                let (loads, report) = scaling
                    .step_curved_walls_with_material_rheology_and_sgs(
                        &mut grid,
                        &wall,
                        &velocities,
                        ORIGIN,
                        &gate,
                        config,
                        &mut work,
                        |point| {
                            calls += 1;
                            let shear = SHEAR + 2.0 * quadratic * (point.cell[axis] as f64 + LOWER);
                            close(point.shear_rate_per_s, shear / DT, 2e-13);
                            Ok(fluid(
                                0.4 + 0.02 * point.cell[0] as f64 + 0.1 * point.shear_rate_per_s,
                            ))
                        },
                    )
                    .unwrap();
                assert_eq!(calls, 64);
                assert_eq!(report.fluid_cells, 64);
                assert_eq!(grid.collision_model(), collision);
                assert_eq!(work.workspace_bytes(), bytes);
                let accepted = work.last_report().unwrap();
                assert_eq!(accepted.sgs, report);
                assert_eq!(accepted.wall_momentum, Some(loads.curved));
                for sample in work.applied() {
                    let mu = 0.4
                        + 0.02 * sample.point.cell[0] as f64
                        + 0.1 * sample.point.shear_rate_per_s;
                    close(sample.transport.dynamic_viscosity_pa_s(), mu, 2e-14);
                    close(
                        sample.shear_relaxation_time,
                        0.5 + 3.0 * mu / 1000.0 * DT / (DX * DX),
                        2e-14,
                    );
                }
                close(report.resolved_dissipation_lu, molecular_power, 3e-15);
                close(report.sgs_dissipation_lu, modeled_power, 3e-15);
                close(report.maximum_eddy_to_molecular, maximum_ratio, 3e-14);
                for (actual, expected) in bits(&grid).into_iter().zip(bits(&expected)) {
                    for (a, b) in actual.into_iter().zip(expected) {
                        close(f64::from_bits(a), f64::from_bits(b), 3e-15);
                    }
                }
                for axis in 0..3 {
                    close(
                        loads.curved.wall_impulse[axis],
                        expected_loads.curved.wall_impulse[axis],
                        2e-12,
                    );
                    close(
                        loads.curved.wall_angular_impulse[axis],
                        expected_loads.curved.wall_angular_impulse[axis],
                        2e-11,
                    );
                }
                close(
                    loads.curved.wall_work,
                    expected_loads.curved.wall_work,
                    2e-13,
                );
                check_exchange(before, inventory(&grid), loads);
            }
        }
    }
}

#[test]
fn constant_source_recovers_curved_sgs_bits_and_keeps_one_read_per_applied_epoch() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    for (original, molecular) in collisions(1.37)
        .into_iter()
        .zip(collisions(scaling.flow().tau))
    {
        let (mut actual, wall, velocities) = channel(1, 0.00025, original);
        let (mut expected, _, _) = channel(1, 0.00025, molecular);
        let mut work = scaling
            .boundary_rheology_sgs_workspace(&actual, 64, &gate)
            .unwrap();
        let mut native = expected.sgs_workspace().unwrap();
        let bytes = work.workspace_bytes();
        for closure in [
            SgsModel::Smagorinsky { cs: 0.0 },
            SgsModel::Wale { cw: 0.5 },
            SgsModel::Vreman { cv: 0.07 },
        ] {
            let native_result = expected
                .step_with_curved_walls_and_sgs(
                    &wall,
                    &velocities,
                    ORIGIN,
                    &gate,
                    settings(closure),
                    &mut native,
                )
                .unwrap();
            let actual_result = scaling
                .step_curved_walls_with_material_rheology_and_sgs(
                    &mut actual,
                    &wall,
                    &velocities,
                    ORIGIN,
                    &gate,
                    settings(closure),
                    &mut work,
                    |point| {
                        Ok(SingleRead {
                            point,
                            reads: Cell::new(0),
                        })
                    },
                )
                .unwrap();
            assert_eq!(actual_result, native_result);
            assert_eq!(bits(&actual), bits(&expected));
            assert_eq!(actual.collision_model(), original);
            assert_eq!(work.workspace_bytes(), bytes);
            assert_eq!(work.applied().len(), 64);
            for sample in work.applied() {
                assert_eq!(sample.point, sample.material.point);
                assert_eq!(sample.material.reads.get(), 1);
                assert_eq!(sample.shear_relaxation_time, scaling.flow().tau);
            }
        }
    }
}

#[test]
fn curved_source_wall_ratio_and_late_cancellation_refusals_preserve_replay() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    let (mut grid, wall, velocities) = channel(1, 0.00025, model(1.37, 0));
    let config = settings(SgsModel::Smagorinsky { cs: 0.2 });
    let mut work = scaling
        .boundary_rheology_sgs_workspace(&grid, 64, &gate)
        .unwrap();
    scaling
        .step_curved_walls_with_material_rheology_and_sgs(
            &mut grid,
            &wall,
            &velocities,
            ORIGIN,
            &gate,
            config,
            &mut work,
            |_| Ok(fluid(0.8)),
        )
        .unwrap();
    let before = bits(&grid);
    let applied = work.applied().to_vec();
    let accepted = work.last_report();
    let mut replay = grid.clone();
    let mut replay_work = scaling
        .boundary_rheology_sgs_workspace(&replay, 64, &gate)
        .unwrap();
    let source_error = MaterialError::State {
        what: "curved source refused its state".into(),
    };
    let mut calls = 0;
    let failed = scaling.step_curved_walls_with_material_rheology_and_sgs(
        &mut grid,
        &wall,
        &velocities,
        ORIGIN,
        &gate,
        config,
        &mut work,
        |_| {
            calls += 1;
            if calls == 64 {
                Err(source_error.clone())
            } else {
                Ok(fluid(1.2))
            }
        },
    );
    assert!(
        matches!(failed, Err(MaterialBoundarySgsError3::Source { source: MaterialLbmError::Transport(e), .. }) if e == source_error)
    );
    assert_eq!(calls, 64);
    assert_eq!(bits(&grid), before);
    assert_eq!(work.applied(), applied);
    let mut invalid = velocities.clone();
    *invalid.last_mut().unwrap() = [f64::NAN; 3];
    assert!(
        scaling
            .step_curved_walls_with_material_rheology_and_sgs(
                &mut grid,
                &wall,
                &invalid,
                ORIGIN,
                &gate,
                config,
                &mut work,
                |_| -> Result<FluidTransportState, MaterialError> {
                    panic!("invalid wall sampled a material")
                },
            )
            .is_err()
    );
    let mut capped = config;
    capped.maximum_eddy_to_molecular = 1.0;
    // The constructor's nu is about 0.29; this source's nu is 1e-6.
    // Only the actual material viscosity can trigger this ratio refusal.
    assert!(matches!(
        scaling.step_curved_walls_with_material_rheology_and_sgs(
            &mut grid,
            &wall,
            &velocities,
            ORIGIN,
            &gate,
            capped,
            &mut work,
            |_| Ok(fluid(1e-5)),
        ),
        Err(MaterialBoundarySgsError3::Sgs(BoundarySgsError3::Closure(
            SgsError::ViscosityRatioExceeded { .. }
        )))
    ));
    assert_eq!(bits(&grid), before);
    assert_eq!(work.applied(), applied);
    let cancelled = CancelGate::new_clock_free();
    calls = 0;
    assert!(matches!(
        scaling.step_curved_walls_with_material_rheology_and_sgs(
            &mut grid,
            &wall,
            &velocities,
            ORIGIN,
            &cancelled,
            config,
            &mut work,
            |_| {
                calls += 1;
                if calls == 64 {
                    cancelled.request();
                }
                Ok(fluid(1.2))
            },
        ),
        Err(MaterialBoundarySgsError3::Sgs(BoundarySgsError3::Boundary(
            BoundaryStepError3::Cancelled
        )))
    ));
    assert_eq!(calls, 64);
    assert_eq!(bits(&grid), before);
    assert_eq!(work.applied(), applied);
    assert_eq!(work.last_report(), accepted);
    let actual = scaling
        .step_curved_walls_with_material_rheology_and_sgs(
            &mut grid,
            &wall,
            &velocities,
            ORIGIN,
            &gate,
            config,
            &mut work,
            |_| Ok(fluid(1.2)),
        )
        .unwrap();
    let expected = scaling
        .step_curved_walls_with_material_rheology_and_sgs(
            &mut replay,
            &wall,
            &velocities,
            ORIGIN,
            &gate,
            config,
            &mut replay_work,
            |_| Ok(fluid(1.2)),
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(bits(&grid), bits(&replay));
    assert_eq!(work.applied(), replay_work.applied());
    assert_eq!(work.last_report(), replay_work.last_report());
}
