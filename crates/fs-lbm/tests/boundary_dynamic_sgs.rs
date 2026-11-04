//! Native dynamic-LES integration controls, not a developed-turbulence claim.

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::{
    BoundaryDynamicSgsConfig3, BoundaryDynamicSgsError3, BoundaryDynamicSgsWorkspace3,
    BoundaryGrid3, BoundarySgsError3, BoundarySpec3, BoundaryStepError3, CollisionModel3, E3,
    FaceBoundary3, TILE, equilibrium3,
};
use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicSgs, ChannelDynamicSgsConfig};
use fs_lbm::sgs3::dynamic::{
    DynamicAveraging, DynamicSgsConfig, DynamicSgsError, PeriodicDynamicSgs,
};
use fs_lbm::sgs3::{SgsError, SgsUse};

const SHAPE: [usize; 3] = [8, 12, 8];
const ORIGIN: [f64; 3] = [2.5, -1.0, 3.0];
const LAW: CollisionModel3 = CollisionModel3::ReducedCumulant {
    second_order_rate: 1.25,
    third_order_rate: 1.1,
    fourth_order_rate: 0.95,
};

fn periodic_config(averaging: DynamicAveraging) -> DynamicSgsConfig {
    DynamicSgsConfig {
        averaging,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn channel_config() -> ChannelDynamicSgsConfig {
    ChannelDynamicSgsConfig {
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn noise(cell: usize, salt: u64) -> f64 {
    let mut h = ((cell as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ salt;
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / 9_007_199_254_740_992.0 * 2.0 - 1.0
}

fn xyz(shape: [usize; 3], cell: usize) -> [usize; 3] {
    [
        cell % shape[0],
        cell / shape[0] % shape[1],
        cell / (shape[0] * shape[1]),
    ]
}

fn native(shape: [usize; 3], [x, y, z]: [usize; 3]) -> usize {
    let tile = (z / TILE * (shape[1] / TILE) + y / TILE) * (shape[0] / TILE) + x / TILE;
    let lane = (z % TILE * TILE + y % TILE) * TILE + x % TILE;
    tile * TILE * TILE * TILE + lane
}

fn walls(normal: usize, boost: f64) -> [[f64; 3]; 2] {
    let mut low = [0.0; 3];
    low[(normal + 1) % 3] = 0.012 + boost;
    low[(normal + 2) % 3] = -0.006;
    let mut high = low;
    high[(normal + 1) % 3] += 0.04;
    [low, high]
}

fn fixture(
    shape: [usize; 3],
    normal: Option<usize>,
    boost: f64,
    varying: bool,
    forced: bool,
    collision: CollisionModel3,
) -> BoundaryGrid3 {
    let mut faces = [FaceBoundary3::Periodic; 6];
    let wall_velocity = normal.map(|axis| {
        let wall = walls(axis, boost);
        faces[2 * axis] = FaceBoundary3::Wall { velocity: wall[0] };
        faces[2 * axis + 1] = FaceBoundary3::Wall { velocity: wall[1] };
        wall
    });
    let [nx, ny, nz] = shape;
    let force = if forced {
        [2e-4, -1e-4, 1.5e-4]
    } else {
        [0.0; 3]
    };
    let mut grid = BoundaryGrid3::with_collision_model(
        nx,
        ny,
        nz,
        collision,
        force,
        BoundarySpec3::new(faces),
    );
    if forced {
        grid.set_force_field(|[x, y, z]| {
            let cell = (z * ny + y) * nx + x;
            core::array::from_fn(|a| 4e-4 * noise(cell, 71 + a as u64))
        });
    }
    for cell in 0..nx * ny * nz {
        let p = xyz(shape, cell);
        let eta = normal.map_or(0.5, |a| (p[a] as f64 + 0.5) / shape[a] as f64);
        let envelope = 4.0 * eta * (1.0 - eta);
        let density = if varying {
            1.0 + 0.08 * noise(cell, 911)
        } else {
            1.0
        };
        let force = grid.force_density(p[0], p[1], p[2]);
        let raw = core::array::from_fn(|a| {
            let mean =
                wall_velocity.map_or(0.0, |wall| (1.0 - eta) * wall[0][a] + eta * wall[1][a]);
            let perturbation = if varying {
                0.02 * envelope * noise(cell, 7 + 31 * a as u64)
            } else {
                0.0
            };
            mean + perturbation - 0.5 * force[a] / density
        });
        let mut f = equilibrium3(density, raw);
        if varying {
            // A prescribed anisotropic stress with zero mass/first moment
            // makes relaxation and wall traction sensitive on the first step.
            let stress = 2e-4 * noise(cell, 101);
            f[1] += stress;
            f[2] += stress;
            f[3] -= stress;
            f[4] -= stress;
        }
        grid.set_populations(p[0], p[1], p[2], &f);
    }
    grid
}

fn fields(grid: &BoundaryGrid3) -> (Vec<f64>, Vec<[f64; 3]>) {
    (0..grid.dimensions().iter().product())
        .map(|cell| {
            let [x, y, z] = xyz(grid.dimensions(), cell);
            (grid.density(x, y, z), grid.velocity(x, y, z))
        })
        .unzip()
}

fn bits(grid: &BoundaryGrid3) -> Vec<u64> {
    (0..grid.dimensions().iter().product())
        .flat_map(|cell| {
            let [x, y, z] = xyz(grid.dimensions(), cell);
            grid.populations(x, y, z).map(f64::to_bits)
        })
        .collect()
}

fn physical_momentum_and_force(grid: &BoundaryGrid3) -> ([f64; 3], [f64; 3]) {
    let mut momentum = [0.0; 3];
    let mut total_force = [0.0; 3];
    for cell in 0..grid.dimensions().iter().product() {
        let [x, y, z] = xyz(grid.dimensions(), cell);
        let rho = grid.density(x, y, z);
        let u = grid.velocity(x, y, z);
        let force = grid.force_density(x, y, z);
        for a in 0..3 {
            momentum[a] += rho * u[a];
            total_force[a] += force[a];
        }
    }
    (momentum, total_force)
}

#[test]
fn periodic_fit_uses_physical_force_frame_and_maps_native_collisions() {
    let gate = CancelGate::new_clock_free();
    let grid = fixture(SHAPE, None, 0.0, true, true, LAW);
    let (rho, u) = fields(&grid);
    let before = bits(&grid);
    let mut work = grid.dynamic_sgs_workspace(&gate, usize::MAX).unwrap();
    let mut oracle = PeriodicDynamicSgs::new(SHAPE, 1.0, usize::MAX).unwrap();
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        let config = periodic_config(averaging);
        let expected = oracle.prepare(config, &rho, &u, || false).unwrap();
        let actual = grid
            .prepare_dynamic_sgs(
                &gate,
                BoundaryDynamicSgsConfig3::Periodic(config),
                &mut work,
            )
            .unwrap();
        assert_eq!(actual.fit, expected);
        assert_eq!(work.eddy_viscosities(), oracle.viscosities());
        if averaging == DynamicAveraging::LocalBox {
            assert!(actual.fit.maximum_eddy_viscosity > 0.0);
            assert!(actual.dissipation.sgs_dissipation_lu > 0.0);
        }
        for (cell, eddy) in oracle.viscosities().iter().copied().enumerate() {
            let collision = work.collisions()[native(SHAPE, xyz(SHAPE, cell))];
            assert!((collision.kinematic_viscosity() - grid.viscosity() - eddy).abs() < 1e-16);
            assert!(matches!(
                collision,
                CollisionModel3::ReducedCumulant {
                    third_order_rate: 1.1,
                    fourth_order_rate: 0.95,
                    ..
                }
            ));
        }
        assert_eq!(bits(&grid), before);
        // A raw-momentum fit is a different model: guard against dropping F/2.
        let raw: Vec<_> = (0..rho.len())
            .map(|cell| {
                let [x, y, z] = xyz(SHAPE, cell);
                let f = grid.populations(x, y, z);
                let mut m = [0.0; 3];
                for (population, e) in f.into_iter().zip(E3) {
                    for (a, component) in [e.0, e.1, e.2].into_iter().enumerate() {
                        m[a] += population * f64::from(component);
                    }
                }
                m.map(|v| v / rho[cell])
            })
            .collect();
        oracle.prepare(config, &rho, &raw, || false).unwrap();
        if averaging == DynamicAveraging::LocalBox {
            assert!(
                work.eddy_viscosities()
                    .iter()
                    .zip(oracle.viscosities())
                    .any(|(a, b)| (a - b).abs() > 1e-8)
            );
        }
    }
}

#[test]
fn channel_fit_reuses_current_moving_walls_in_every_orientation() {
    let gate = CancelGate::new_clock_free();
    for normal in 0..3 {
        let grid = fixture(SHAPE, Some(normal), 0.0, true, true, LAW);
        let mut work = grid.dynamic_sgs_workspace(&gate, usize::MAX).unwrap();
        let collisions_ptr = work.collisions().as_ptr();
        let eddy_ptr = work.eddy_viscosities().as_ptr();
        let bytes = work.workspace_bytes();
        let mut oracle = ChannelDynamicSgs::new(SHAPE, 1.0, normal, usize::MAX).unwrap();
        // Reuse the workspace on the same topology with changed physical walls.
        for boost in [0.0, 0.017] {
            let grid = fixture(SHAPE, Some(normal), boost, true, true, LAW);
            let (rho, u) = fields(&grid);
            let expected = oracle
                .prepare(channel_config(), &rho, &u, walls(normal, boost), || false)
                .unwrap();
            let actual = grid
                .prepare_dynamic_sgs(
                    &gate,
                    BoundaryDynamicSgsConfig3::Channel(channel_config()),
                    &mut work,
                )
                .unwrap();
            assert_eq!(actual.fit, expected);
            assert_eq!(work.eddy_viscosities(), oracle.viscosities());
            assert!(actual.fit.maximum_eddy_viscosity > 0.0, "normal={normal}");
            assert_eq!(work.collisions().as_ptr(), collisions_ptr);
            assert_eq!(work.eddy_viscosities().as_ptr(), eddy_ptr);
            assert_eq!(work.workspace_bytes(), bytes);
            oracle
                .prepare(channel_config(), &rho, &u, [[0.0; 3]; 2], || false)
                .unwrap();
            assert!(
                work.eddy_viscosities()
                    .iter()
                    .zip(oracle.viscosities())
                    .any(|(a, b)| (a - b).abs() > 1e-8),
                "wall motion must enter the fit on axis {normal}"
            );
        }
    }
}

#[test]
fn laminar_couette_has_no_modeled_transport_and_preserves_original_step() {
    let gate = CancelGate::new_clock_free();
    for collision in [
        CollisionModel3::Bgk { tau: 0.8 },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.25,
            higher_order_rate: 1.1,
        },
        LAW,
    ] {
        for normal in 0..3 {
            let mut flow = fixture([8; 3], Some(normal), 0.0, false, false, collision);
            let mut original = flow.clone();
            let links = flow.boundary_links();
            let mut work = flow.dynamic_sgs_workspace(&gate, usize::MAX).unwrap();
            let expected = original.step_with_wall_momentum(&links, ORIGIN);
            let actual = flow
                .step_with_dynamic_sgs(
                    &gate,
                    BoundaryDynamicSgsConfig3::Channel(channel_config()),
                    &mut work,
                    Some((&links, ORIGIN)),
                )
                .unwrap();
            assert!(actual.dynamic.fit.maximum_eddy_viscosity < 1e-20);
            assert!(actual.dynamic.dissipation.sgs_dissipation_lu < 1e-20);
            assert_eq!(work.collisions(), vec![collision; 8 * 8 * 8]);
            assert_eq!(bits(&flow), bits(&original));
            assert_eq!(actual.wall_momentum, Some(expected));
        }
    }
}

#[test]
fn native_serial_and_pooled_dynamic_steps_close_actual_wall_momentum() {
    let gate = CancelGate::new_clock_free();
    let config = BoundaryDynamicSgsConfig3::Channel(channel_config());
    let initial = fixture(SHAPE, Some(1), 0.0, true, true, LAW);
    let links = initial.boundary_links();
    let mut ordinary = initial.clone();
    let ordinary_loads = ordinary.step_with_wall_momentum(&links, ORIGIN);
    for workers in [1, 2, 3] {
        let pool = TilePool::for_host(workers, 0xBD51_59D3);
        let mut serial = initial.clone();
        let mut pooled = initial.clone();
        let mut work = serial.dynamic_sgs_workspace(&gate, usize::MAX).unwrap();
        let mut parallel = pooled.dynamic_sgs_workspace(&gate, usize::MAX).unwrap();
        for iteration in 0..2 {
            let (before, force) = physical_momentum_and_force(&serial);
            let actual = serial
                .step_with_dynamic_sgs(&gate, config, &mut work, Some((&links, ORIGIN)))
                .unwrap();
            let replay = pooled
                .step_pooled_with_dynamic_sgs(
                    &pool,
                    &gate,
                    config,
                    &mut parallel,
                    Some((&links, ORIGIN)),
                )
                .unwrap();
            assert_eq!(actual, replay);
            assert_eq!(bits(&serial), bits(&pooled));
            assert!(actual.dynamic.fit.maximum_eddy_viscosity > 0.0);
            let loads = actual.wall_momentum.unwrap();
            let (after, _) = physical_momentum_and_force(&serial);
            for a in 0..3 {
                assert!(
                    (after[a] - before[a] - force[a] - loads.fluid_population_impulse[a]).abs()
                        < 2e-12
                );
                assert!(
                    (loads.wall_impulse[a] + loads.fluid_population_impulse[a]
                        - loads.wall_velocity_mass_impulse[a])
                        .abs()
                        < 2e-12
                );
            }
            assert_eq!(loads.measured_links, links.len());
            if iteration == 0 {
                assert_ne!(bits(&serial), bits(&ordinary));
                assert_ne!(loads.wall_impulse, ordinary_loads.wall_impulse);
            }
        }
    }
}

#[test]
fn caps_acceptance_and_cancellation_refuse_without_advancing_and_allow_retry() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0xBD51_59D3);
    for pooled in [false, true] {
        let mut flow = fixture(SHAPE, None, 0.0, true, true, LAW);
        let initial = flow.clone();
        let before = bits(&flow);
        let config = periodic_config(DynamicAveraging::LocalBox);
        let mut work = flow.dynamic_sgs_workspace(&gate, usize::MAX).unwrap();
        let step = |flow: &mut BoundaryGrid3,
                    gate: &CancelGate,
                    config,
                    work: &mut BoundaryDynamicSgsWorkspace3| {
            if pooled {
                flow.step_pooled_with_dynamic_sgs(
                    &pool,
                    gate,
                    BoundaryDynamicSgsConfig3::Periodic(config),
                    work,
                    None,
                )
            } else {
                flow.step_with_dynamic_sgs(
                    gate,
                    BoundaryDynamicSgsConfig3::Periodic(config),
                    work,
                    None,
                )
            }
        };
        assert!(matches!(
            step(
                &mut flow,
                &gate,
                DynamicSgsConfig {
                    maximum_coefficient_squared: 0.0,
                    ..config
                },
                &mut work
            ),
            Err(BoundaryDynamicSgsError3::Dynamic(
                DynamicSgsError::Coefficient { .. }
            ))
        ));
        assert_eq!(bits(&flow), before);
        assert!(matches!(
            step(
                &mut flow,
                &gate,
                DynamicSgsConfig {
                    maximum_eddy_to_molecular: 0.0,
                    ..config
                },
                &mut work
            ),
            Err(BoundaryDynamicSgsError3::Boundary(
                BoundarySgsError3::Closure(SgsError::ViscosityRatioExceeded { .. })
            ))
        ));
        assert_eq!(bits(&flow), before);
        assert!(matches!(
            step(
                &mut flow,
                &gate,
                DynamicSgsConfig {
                    usage: SgsUse::AcceptanceRender,
                    ..config
                },
                &mut work
            ),
            Err(BoundaryDynamicSgsError3::Dynamic(DynamicSgsError::Invalid(
                _
            )))
        ));
        assert_eq!(bits(&flow), before);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            step(&mut flow, &cancelled, config, &mut work),
            Err(BoundaryDynamicSgsError3::Boundary(
                BoundarySgsError3::Boundary(BoundaryStepError3::Cancelled)
            ))
        ));
        assert_eq!(bits(&flow), before);
        let mut reference = initial;
        let mut fresh = reference.dynamic_sgs_workspace(&gate, usize::MAX).unwrap();
        let expected = reference
            .step_with_dynamic_sgs(
                &gate,
                BoundaryDynamicSgsConfig3::Periodic(config),
                &mut fresh,
                None,
            )
            .unwrap();
        assert_eq!(step(&mut flow, &gate, config, &mut work).unwrap(), expected);
        assert_eq!(bits(&flow), bits(&reference));
    }
}

#[test]
fn workspace_binds_shape_and_topology_and_enforces_complete_byte_budget() {
    let gate = CancelGate::new_clock_free();
    let mut flow = fixture(SHAPE, Some(1), 0.0, true, false, LAW);
    let count = SHAPE.iter().product();
    let required = BoundaryDynamicSgsWorkspace3::required_bytes(count).unwrap();
    assert!(
        matches!(flow.dynamic_sgs_workspace(&gate, required - 1), Err(BoundaryDynamicSgsError3::Dynamic(DynamicSgsError::WorkspaceBytes { required: r, maximum: m })) if r == required && m == required - 1)
    );
    let mut work = flow.dynamic_sgs_workspace(&gate, required).unwrap();
    assert_eq!(work.workspace_bytes(), required);
    assert!(BoundaryDynamicSgsWorkspace3::required_bytes(usize::MAX).is_err());
    let config = BoundaryDynamicSgsConfig3::Channel(channel_config());
    for other in [
        fixture([12, 8, 8], Some(1), 0.0, true, false, LAW),
        fixture(SHAPE, Some(0), 0.0, true, false, LAW),
        fixture(SHAPE, None, 0.0, true, false, LAW),
    ] {
        assert!(matches!(
            other.prepare_dynamic_sgs(&gate, config, &mut work),
            Err(BoundaryDynamicSgsError3::Boundary(
                BoundarySgsError3::Workspace
            ))
        ));
    }
    assert!(matches!(
        flow.prepare_dynamic_sgs(
            &gate,
            BoundaryDynamicSgsConfig3::Periodic(periodic_config(DynamicAveraging::Domain)),
            &mut work
        ),
        Err(BoundaryDynamicSgsError3::Dynamic(DynamicSgsError::Invalid(
            _
        )))
    ));
    flow.voxelize_sdf(|p| {
        if p[0] < 1.0 && p[1] < 1.0 && p[2] < 1.0 {
            -1.0
        } else {
            1.0
        }
    });
    assert!(matches!(
        flow.prepare_dynamic_sgs(&gate, config, &mut work),
        Err(BoundaryDynamicSgsError3::Dynamic(DynamicSgsError::Invalid(
            _
        )))
    ));
    for boundaries in [
        BoundarySpec3::lid_cavity([0.01, 0.0, 0.0]),
        BoundarySpec3::velocity_pressure_duct_z([0.0, 0.0, 0.01], 1.0),
    ] {
        let unsupported = BoundaryGrid3::new(8, 8, 8, 0.8, [0.0; 3], boundaries);
        assert!(matches!(
            unsupported.dynamic_sgs_workspace(&gate, usize::MAX),
            Err(BoundaryDynamicSgsError3::Dynamic(DynamicSgsError::Invalid(
                _
            )))
        ));
    }
}
