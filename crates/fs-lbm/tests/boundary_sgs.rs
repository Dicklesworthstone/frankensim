//! Analytic and replay controls for explicit SGS on the native boundary grid.
//! These small fixtures check closure integration, not turbulent DNS accuracy.

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySgsError3, BoundarySpec3, BoundaryStepError3, CollisionModel3,
    FaceBoundary3, TILE, equilibrium3,
};
use fs_lbm::sgs3::{SgsConfig, SgsError, SgsModel, SgsUse, SgsWorkspace3};

const SEED: u64 = 0xB035_0535;
const ORIGIN: [f64; 3] = [2.5, -1.0, 3.0];

fn config(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn fluid_bits(grid: &BoundaryGrid3) -> Vec<u64> {
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

fn flat_index([nx, ny, _]: [usize; 3], [x, y, z]: [usize; 3]) -> usize {
    let tile = (z / TILE * (ny / TILE) + y / TILE) * (nx / TILE) + x / TILE;
    let lane = (z % TILE * TILE + y % TILE) * TILE + x % TILE;
    tile * TILE * TILE * TILE + lane
}

fn fixture(
    dimensions: [usize; 3],
    boundaries: BoundarySpec3,
    collision: CollisionModel3,
    obstacle: bool,
) -> BoundaryGrid3 {
    let [nx, ny, nz] = dimensions;
    let mut grid = BoundaryGrid3::with_collision_model(nx, ny, nz, collision, [0.0; 3], boundaries);
    if obstacle {
        grid.voxelize_sdf(|p| {
            let r = core::array::from_fn::<_, 3, _>(|a| p[a] - dimensions[a] as f64 * 0.5);
            r.iter().map(|v| v * v).sum::<f64>() - 4.0
        });
    }
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if !grid.is_solid(x, y, z) {
                    let u = [
                        0.015 * (y as f64 + 0.5) / ny as f64,
                        -0.01 * (z as f64 + 0.5) / nz as f64,
                        0.007 * (x as f64 + 0.5) / nx as f64,
                    ];
                    grid.set_populations(x, y, z, &equilibrium3(1.0, u));
                }
            }
        }
    }
    grid
}

#[test]
fn zero_sgs_preserves_populations_and_actual_wall_loads() {
    let gate = CancelGate::new_clock_free();
    let laws = [
        CollisionModel3::Bgk { tau: 0.8 },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.25,
            higher_order_rate: 1.1,
        },
        CollisionModel3::ReducedCumulant {
            second_order_rate: 1.25,
            third_order_rate: 1.1,
            fourth_order_rate: 1.0,
        },
    ];
    for boundaries in [
        BoundarySpec3::lid_cavity([0.025, 0.0, 0.01]),
        BoundarySpec3::velocity_pressure_duct_z([0.0, 0.0, 0.02], 1.0),
        BoundarySpec3::periodic(),
    ] {
        for collision in laws {
            let mut ordinary = fixture([8; 3], boundaries, collision, false);
            let mut controlled = fixture([8; 3], boundaries, collision, false);
            let links = ordinary.boundary_links();
            let mut workspace = controlled.sgs_workspace().unwrap();
            let allocation = workspace.collisions().as_ptr();
            for model in [
                SgsModel::Smagorinsky { cs: 0.0 },
                SgsModel::Wale { cw: 0.0 },
            ] {
                for _ in 0..2 {
                    let expected = ordinary.step_with_wall_momentum(&links, ORIGIN);
                    let actual = controlled
                        .step_with_sgs(&gate, config(model), &mut workspace, Some((&links, ORIGIN)))
                        .unwrap();
                    assert_eq!(actual.sgs.sgs_dissipation_lu, 0.0);
                    assert_eq!(actual.wall_momentum, Some(expected));
                    assert_eq!(fluid_bits(&controlled), fluid_bits(&ordinary));
                    assert_eq!(controlled.collision_model(), collision);
                    assert_eq!(workspace.collisions().as_ptr(), allocation);
                }
            }
        }
    }
}

const COUETTE_HEIGHT: usize = 8;
const COUETTE_SPEED: f64 = 0.04;

fn couette(with_force: bool) -> BoundaryGrid3 {
    let periodic = FaceBoundary3::Periodic;
    let bottom = [0.012, 0.0, -0.006];
    let top = [bottom[0] + COUETTE_SPEED, 0.0, bottom[2]];
    let boundaries = BoundarySpec3::new([
        periodic,
        periodic,
        FaceBoundary3::Wall { velocity: bottom },
        FaceBoundary3::Wall { velocity: top },
        periodic,
        periodic,
    ]);
    let uniform_force = if with_force {
        [1e-4, 0.0, 0.0]
    } else {
        [0.0; 3]
    };
    let mut grid = BoundaryGrid3::new(4, COUETTE_HEIGHT, 4, 0.8, uniform_force, boundaries);
    let field = |p: [usize; 3]| [2e-4 * (p[1] as f64 + 0.5), 0.0, 0.0];
    if with_force {
        grid.set_force_field(field);
    }
    for z in 0..4 {
        for y in 0..COUETTE_HEIGHT {
            for x in 0..4 {
                let physical = [
                    bottom[0] + COUETTE_SPEED * (y as f64 + 0.5) / COUETTE_HEIGHT as f64,
                    0.0,
                    bottom[2],
                ];
                // Guo forcing exposes u = (momentum + F/2) / rho. Seed the
                // same physical shear with and without a spatial force field.
                let raw = core::array::from_fn(|a| {
                    physical[a]
                        - if with_force {
                            0.5 * (uniform_force[a] + field([x, y, z])[a])
                        } else {
                            0.0
                        }
                });
                grid.set_populations(x, y, z, &equilibrium3(1.0, raw));
                for (actual, expected) in grid.velocity(x, y, z).into_iter().zip(physical) {
                    assert!((actual - expected).abs() < 6e-17);
                }
            }
        }
    }
    grid
}

#[test]
fn halfway_couette_uses_physical_shear_and_wale_vanishes() {
    let gate = CancelGate::new_clock_free();
    let shear = COUETTE_SPEED / COUETTE_HEIGHT as f64;
    let mut settings = config(SgsModel::Smagorinsky { cs: 0.2 });
    settings.filter_width_lu = 1.5;
    let expected = (0.2_f64 * 1.5).powi(2) * shear;
    for with_force in [false, true] {
        let grid = couette(with_force);
        let before = fluid_bits(&grid);
        let mut workspace = grid.sgs_workspace().unwrap();
        let report = grid.prepare_sgs(&gate, settings, &mut workspace).unwrap();
        assert_eq!(report.fluid_cells, 4 * COUETTE_HEIGHT * 4);
        assert!((report.maximum_eddy_viscosity_lu - expected).abs() < 2e-16);
        // This includes the first/last cell centers a half-cell from a wall.
        for (cell, collision) in workspace.collisions().iter().enumerate() {
            let observed = collision.kinematic_viscosity() - grid.viscosity();
            assert!(
                (observed - expected).abs() < 2e-16,
                "cell={cell} forced={with_force}: {observed} != {expected}"
            );
        }
        let expected_dissipation = report.fluid_cells as f64 * expected * shear * shear;
        assert!((report.sgs_dissipation_lu - expected_dissipation).abs() < 1e-16);
        let wale = grid
            .prepare_sgs(&gate, config(SgsModel::Wale { cw: 0.5 }), &mut workspace)
            .unwrap();
        assert!(wale.maximum_eddy_viscosity_lu < 1e-20);
        assert_eq!(
            fluid_bits(&grid),
            before,
            "preparation cannot advance populations"
        );
    }
}

#[test]
fn nonzero_sgs_replays_populations_and_wall_loads_for_every_worker_count() {
    let make_grid = || {
        // Four pooled kernel groups, moving-wall corners, stationary voxel
        // links, and nonuniform Guo forcing all use the native stepping path.
        let mut grid = fixture(
            [16, 8, 16],
            BoundarySpec3::lid_cavity([0.025, 0.0, 0.01]),
            CollisionModel3::Bgk { tau: 0.8 },
            true,
        );
        grid.set_force_field(|p| [1e-6 * (p[1] as f64 + 1.0), -2e-6, 0.0]);
        grid
    };
    for model in [
        SgsModel::Smagorinsky { cs: 0.17 },
        SgsModel::Wale { cw: 0.5 },
    ] {
        let gate = CancelGate::new_clock_free();
        let mut serial = make_grid();
        let links = serial.boundary_links();
        let mut workspace = serial.sgs_workspace().unwrap();
        let mut expected = Vec::new();
        for _ in 0..3 {
            let report = serial
                .step_with_sgs(&gate, config(model), &mut workspace, Some((&links, ORIGIN)))
                .unwrap();
            assert!(report.sgs.maximum_eddy_viscosity_lu > 0.0);
            assert_eq!(report.wall_momentum.unwrap().measured_links, links.len());
            expected.push((report, fluid_bits(&serial)));
        }
        for workers in [1, 2, 4] {
            let pool = TilePool::for_host(workers, SEED);
            let mut pooled = make_grid();
            let mut workspace = pooled.sgs_workspace().unwrap();
            for (reference, populations) in &expected {
                let actual = pooled
                    .step_pooled_with_sgs(
                        &pool,
                        &gate,
                        config(model),
                        &mut workspace,
                        Some((&links, ORIGIN)),
                    )
                    .unwrap();
                assert_eq!(actual.sgs, reference.sgs, "workers={workers}");
                assert_eq!(
                    actual.wall_momentum, reference.wall_momentum,
                    "workers={workers}"
                );
                assert_eq!(&fluid_bits(&pooled), populations, "workers={workers}");
            }
        }
    }
}

#[test]
fn invalid_closure_cap_and_cancellation_leave_a_retryable_state() {
    let pool = TilePool::for_host(2, SEED);
    for pooled in [false, true] {
        let step = |grid: &mut BoundaryGrid3,
                    gate: &CancelGate,
                    settings: SgsConfig,
                    workspace: &mut SgsWorkspace3<CollisionModel3>| {
            if pooled {
                grid.step_pooled_with_sgs(&pool, gate, settings, workspace, None)
            } else {
                grid.step_with_sgs(gate, settings, workspace, None)
            }
        };
        let gate = CancelGate::new_clock_free();
        let mut grid = couette(false);
        let before = fluid_bits(&grid);
        let mut workspace = grid.sgs_workspace().unwrap();
        assert!(matches!(
            step(
                &mut grid,
                &gate,
                config(SgsModel::Smagorinsky { cs: -0.1 }),
                &mut workspace
            ),
            Err(BoundarySgsError3::Closure(SgsError::InvalidCoefficient))
        ));
        assert_eq!(fluid_bits(&grid), before);
        let settings = config(SgsModel::Smagorinsky { cs: 0.17 });
        let mut capped = settings;
        capped.maximum_eddy_to_molecular = 0.0;
        assert!(matches!(
            step(&mut grid, &gate, capped, &mut workspace),
            Err(BoundarySgsError3::Closure(
                SgsError::ViscosityRatioExceeded { .. }
            ))
        ));
        assert_eq!(fluid_bits(&grid), before);
        gate.request();
        assert!(matches!(
            step(&mut grid, &gate, settings, &mut workspace),
            Err(BoundarySgsError3::Boundary(BoundaryStepError3::Cancelled))
        ));
        assert_eq!(fluid_bits(&grid), before);
        let fresh = CancelGate::new_clock_free();
        let actual = step(&mut grid, &fresh, settings, &mut workspace).unwrap();
        let mut reference = couette(false);
        let mut reference_workspace = reference.sgs_workspace().unwrap();
        let expected = reference
            .step_with_sgs(&fresh, settings, &mut reference_workspace, None)
            .unwrap();
        assert_eq!(actual.sgs, expected.sgs);
        assert_eq!(actual.wall_momentum, None);
        assert_eq!(fluid_bits(&grid), fluid_bits(&reference));
    }
}

#[test]
fn open_faces_use_onsite_velocity_and_extrapolated_pressure_gradients() {
    let periodic = FaceBoundary3::Periodic;
    let boundaries = BoundarySpec3::new([
        FaceBoundary3::Velocity {
            velocity: [0.01, 0.0, 0.0],
        },
        FaceBoundary3::Pressure { density: 1.0 },
        periodic,
        periodic,
        periodic,
        periodic,
    ]);
    let dimensions = [8, 4, 4];
    let mut grid = BoundaryGrid3::new(8, 4, 4, 0.8, [0.0; 3], boundaries);
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..8 {
                let s = x as f64;
                grid.set_populations(
                    x,
                    y,
                    z,
                    &equilibrium3(1.0, [0.01 + 0.001 * s + 0.0002 * s * s, 0.0, 0.0]),
                );
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    let mut workspace = grid.sgs_workspace().unwrap();
    grid.prepare_sgs(
        &gate,
        config(SgsModel::Smagorinsky { cs: 0.2 }),
        &mut workspace,
    )
    .unwrap();
    for x in 0..8 {
        // The prescribed velocity lives at x=0 (not x=-1/2). A second-order
        // one-sided derivative is exact for this quadratic; the pressure
        // outlet instead extrapolates velocity with zero normal derivative.
        let derivative = if x == 7 {
            0.0
        } else {
            0.001 + 0.0004 * x as f64
        };
        let expected = 0.2 * 0.2 * fs_math::det::sqrt(2.0) * derivative;
        let collision = workspace.collisions()[flat_index(dimensions, [x, 1, 1])];
        let actual = collision.kinematic_viscosity() - grid.viscosity();
        assert!(
            (actual - expected).abs() < 2e-16,
            "x={x}: {actual} != {expected}"
        );
    }
}
