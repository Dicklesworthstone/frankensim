//! Native channel/history transactions and physical wall-frame regressions.

use super::*;
use crate::d3q19::{CollisionModel3, equilibrium3};
use crate::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use crate::lattice3::{D3q19, GridBudget3};
use crate::sgs3::SgsUse;
use fs_exec::{RunError, RunReport, TileKernel, TilePool};
use std::sync::atomic::{AtomicUsize, Ordering};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 503,
                kernel_id: 136,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn config() -> LagrangianSgsConfig {
    LagrangianSgsConfig {
        memory_time: 4.0,
        maximum_departure_cells: 1.0,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn cumulant() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Cumulant {
            rates: Rates27 {
                omega: [1.8, 1.13, 1.21, 1.34, 1.47, 1.58, 1.66, 1.72, 1.39, 1.51],
            },
        },
        maximum_mach: 0.2,
    }
}

fn eq27(rho: f64, velocity: [f64; 3]) -> [f64; 27] {
    equilibrium27(rho, velocity, Equilibrium27::Factorized)
}

fn random(cell: usize, component: usize) -> f64 {
    let mut hash =
        (cell as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ (31 * component as u64 + 7);
    hash ^= hash >> 30;
    hash = hash.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    hash ^= hash >> 27;
    hash = hash.wrapping_mul(0x94d0_49bb_1331_11eb);
    hash ^= hash >> 31;
    (hash >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}

fn flow<L: SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
    normal: Option<usize>,
    force: [f64; 3],
    walls: [[f64; 3]; 2],
    varying: bool,
) -> Grid3<L, Q> {
    let axis = normal.unwrap_or(2);
    Grid3::try_new(
        [5; 3],
        collision,
        force,
        core::array::from_fn(|i| normal != Some(i)),
        GridBudget3 {
            maximum_cells: 125,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let axes = [(axis + 1) % 3, (axis + 2) % 3, axis];
            let [x, y, z] = axes.map(|i| p[i]);
            let cell = (z * 5 + y) * 5 + x;
            let h = (p[axis] as f64 + 0.5) / 5.0;
            let rho = 1.0
                + if varying {
                    0.002 * random(cell, 3)
                } else {
                    0.0
                };
            let mut velocity = core::array::from_fn(|i| walls[0][i] * (1.0 - h) + walls[1][i] * h);
            if varying {
                for i in 0..3 {
                    velocity[axes[i]] += 0.007 * 4.0 * h * (1.0 - h) * random(cell, i);
                }
            }
            for i in 0..3 {
                velocity[i] -= 0.5 * force[i] / rho;
            }
            equilibrium(rho, velocity)
        },
    )
    .unwrap()
}

fn drive<L: Lattice3<Q>, const Q: usize>(
    flow: &Grid3<L, Q>,
    cx: &Cx<'_>,
    normal: usize,
    walls: [[f64; 3]; 2],
) -> (WallMomentumProbe3<L, Q>, Vec<[f64; 3]>) {
    let probe = flow.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
    let velocity = probe
        .wall_links()
        .map(|link| walls[usize::from(link.exterior_faces[2 * normal + 1])])
        .collect();
    (probe, velocity)
}

fn reference<L: SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
    moving: bool,
) where
    L::Collision: core::fmt::Debug + PartialEq,
{
    for normal in 0..3 {
        let tangent = (normal + 1) % 3;
        let mut force = [0.0; 3];
        force[tangent] = 2e-5;
        let mut walls = [[0.0; 3]; 2];
        if moving {
            walls[0][tangent] = 0.004;
            walls[1][tangent] = 0.012;
        }
        let new = || flow::<L, Q>(collision, equilibrium, Some(normal), force, walls, true);
        let mut serial = new();
        let mut native = new();
        let mut pooled = new();
        let mut work = serial.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut parallel = pooled.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut model =
            LagrangianDynamicSgs::new_channel([5; 3], [1.0; 3], normal, usize::MAX).unwrap();
        let molecular: Vec<_> = (0..125)
            .map(|cell| L::with_shear_relaxation(collision, 0.7 + 0.0005 * cell as f64).unwrap())
            .collect();
        let original = molecular.clone();
        let bytes = work.workspace_bytes();
        let pointer = work.collisions().as_ptr();
        assert_eq!(work.normal_axis(), Some(normal));
        let pool = TilePool::for_host(2, 503);
        let gate = CancelGate::new_clock_free();
        let mut saw_transport = false;
        with_cx(|cx| {
            let (probe, velocities) = drive(&serial, cx, normal, walls);
            for step in 0..3 {
                let macros: Vec<_> = (0..125)
                    .map(|cell| native.macroscopic(cell).unwrap())
                    .collect();
                let density: Vec<_> = macros.iter().map(|m| m.0).collect();
                let velocity: Vec<_> = macros.iter().map(|m| m.1).collect();
                let trial = model
                    .prepare_with_walls(config(), step as f64, &density, &velocity, walls, || false)
                    .unwrap();
                let rates: Vec<_> = trial
                    .viscosities()
                    .iter()
                    .enumerate()
                    .map(|(cell, &eddy)| {
                        let base = molecular[cell];
                        if eddy == 0.0 {
                            base
                        } else {
                            let nu = L::molecular_viscosity(base).unwrap();
                            L::with_shear_relaxation(base, 0.5 + (nu + eddy) / L::CS2).unwrap()
                        }
                    })
                    .collect();
                let expected = trial.report();
                let (actual, pooled_report) = if moving {
                    let native_load = native
                        .step_serial_with_moving_wall_momentum(
                            cx,
                            &probe,
                            &velocities,
                            Some(&rates),
                        )
                        .unwrap();
                    let actual = serial
                        .step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                            cx,
                            config(),
                            Some(&molecular),
                            &probe,
                            &velocities,
                            &mut work,
                        )
                        .unwrap();
                    let parallel = pooled
                        .step_pooled_with_channel_lagrangian_sgs_and_moving_walls(
                            &pool,
                            &gate,
                            config(),
                            Some(&molecular),
                            &probe,
                            &velocities,
                            &mut parallel,
                        )
                        .unwrap();
                    assert_eq!(actual.1, native_load);
                    assert_eq!(parallel.1, native_load);
                    (actual.0, parallel.0)
                } else {
                    native.step_serial_with_collisions(cx, &rates).unwrap();
                    (
                        serial
                            .step_serial_with_lagrangian_sgs(
                                cx,
                                config(),
                                Some(&molecular),
                                &mut work,
                            )
                            .unwrap(),
                        pooled
                            .step_pooled_with_lagrangian_sgs(
                                &pool,
                                &gate,
                                config(),
                                Some(&molecular),
                                &mut parallel,
                            )
                            .unwrap(),
                    )
                };
                trial.commit();
                assert_eq!(actual.fit, expected);
                assert_eq!(actual, pooled_report);
                assert_eq!(serial.populations(), native.populations());
                assert_eq!(serial.populations(), pooled.populations());
                assert_eq!(work.collisions(), rates);
                assert_eq!(work.accepted_fit_time(), Some(step as f64));
                assert_eq!(parallel.accepted_fit_time(), Some(step as f64));
                assert_eq!(serial.completed_steps(), step + 1);
                assert_eq!(work.workspace_bytes(), bytes);
                assert_eq!(work.collisions().as_ptr(), pointer);
                saw_transport |= actual.fit.maximum_eddy_viscosity > 0.0;
            }
        });
        assert!(saw_transport);
        assert_eq!(molecular, original);
    }
}

#[test]
fn stationary_channels_match_independent_history_and_native_transport() {
    reference::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3, false);
    reference::<D3q27, 27>(cumulant(), eq27, false);
}

#[test]
fn moving_channels_match_independent_history_wall_receipts_and_worker_counts() {
    reference::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3, true);
    reference::<D3q27, 27>(cumulant(), eq27, true);
}

#[test]
fn workspace_binds_periodicity_normal_and_actual_storage() {
    let law = CollisionModel3::Bgk { tau: 0.8 };
    let new = |normal| flow::<D3q19, 19>(law, equilibrium3, normal, [0.0; 3], [[0.0; 3]; 2], true);
    let mut channel = new(Some(2));
    let mut other = new(Some(1));
    let mut periodic = new(None);
    let mut work = channel
        .channel_lagrangian_sgs_workspace(usize::MAX)
        .unwrap();
    let mut periodic_work = periodic.lagrangian_sgs_workspace(usize::MAX).unwrap();
    let bytes = work.workspace_bytes();
    assert!(channel.channel_lagrangian_sgs_workspace(bytes - 1).is_err());
    assert_eq!(
        channel
            .channel_lagrangian_sgs_workspace(bytes)
            .unwrap()
            .workspace_bytes(),
        bytes
    );
    assert!(
        periodic
            .channel_lagrangian_sgs_workspace(usize::MAX)
            .is_err()
    );
    assert!(channel.lagrangian_sgs_workspace(usize::MAX).is_err());
    let before = channel.populations().to_vec();
    with_cx(|cx| {
        assert!(
            channel
                .step_serial_with_lagrangian_sgs(cx, config(), None, &mut periodic_work)
                .is_err()
        );
        assert!(
            periodic
                .step_serial_with_lagrangian_sgs(cx, config(), None, &mut work)
                .is_err()
        );
        assert!(
            other
                .step_serial_with_lagrangian_sgs(cx, config(), None, &mut work)
                .is_err()
        );
    });
    assert_eq!(work.accepted_fit_time(), None);
    assert_eq!(periodic_work.accepted_fit_time(), None);
    assert_eq!(channel.populations(), before);
    assert_eq!(channel.completed_steps(), 0);
    assert_eq!(periodic.completed_steps(), 0);
    assert_eq!(other.completed_steps(), 0);
}

fn laminar<L: SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) {
    for normal in 0..3 {
        let tangent = (normal + 1) % 3;
        for difference in [0.0, 0.02] {
            let mut walls = [[0.0; 3]; 2];
            walls[0][tangent] = 0.01;
            walls[1][tangent] = 0.01 + difference;
            let mut grid =
                flow::<L, Q>(collision, equilibrium, Some(normal), [0.0; 3], walls, false);
            let mut work = grid.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
            with_cx(|cx| {
                let (probe, velocity) = drive(&grid, cx, normal, walls);
                let (report, load) = grid
                    .step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                        cx,
                        config(),
                        None,
                        &probe,
                        &velocity,
                        &mut work,
                    )
                    .unwrap();
                let expected =
                    125.0 * L::molecular_viscosity(collision).unwrap() * (difference / 5.0).powi(2);
                assert!(report.fit.maximum_eddy_viscosity < 1e-14);
                assert!(report.dissipation.sgs_dissipation_lu < 1e-16);
                if difference == 0.0 {
                    assert!(report.dissipation.resolved_dissipation_lu < 1e-26);
                    assert!(load.wall_impulse[tangent].abs() < 1e-13);
                    assert!(load.wall_work.abs() < 1e-14);
                } else {
                    assert!(
                        (report.dissipation.resolved_dissipation_lu - expected).abs()
                            < 1e-11 * expected
                    );
                }
                assert_eq!(work.accepted_fit_time(), Some(0.0));
            });
        }
    }
}

#[test]
fn couette_and_comotion_use_physical_wall_strain_on_both_lattices() {
    laminar::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
    laminar::<D3q27, 27>(cumulant(), eq27);
}

struct RefuseAfterStreaming {
    pool: TilePool,
    calls: AtomicUsize,
}
impl KernelRunner for RefuseAfterStreaming {
    fn workers(&self) -> usize {
        self.pool.workers()
    }
    fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        let result = self.pool.run_with_gate(kernel, gate);
        if self.calls.fetch_add(1, Ordering::AcqRel) == 1 {
            gate.request();
        }
        result
    }
}

#[test]
fn invalid_drive_and_late_cancellation_leave_accepted_history_retryable() {
    let law = cumulant();
    let walls = [[0.004, 0.0, 0.0], [0.012, 0.0, 0.0]];
    let new = || flow::<D3q27, 27>(law, eq27, Some(2), [2e-5, 0.0, 0.0], walls, true);
    let mut grid = new();
    let mut control = new();
    let mut work = grid.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
    let mut clean = control
        .channel_lagrangian_sgs_workspace(usize::MAX)
        .unwrap();
    with_cx(|cx| {
        let (probe, velocity) = drive(&grid, cx, 2, walls);
        grid.step_serial_with_channel_lagrangian_sgs_and_moving_walls(
            cx,
            config(),
            None,
            &probe,
            &velocity,
            &mut work,
        )
        .unwrap();
        control
            .step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &velocity,
                &mut clean,
            )
            .unwrap();
        let before = grid.populations().to_vec();
        let partial = grid
            .wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[5])
            .unwrap();
        assert!(
            grid.step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &partial,
                &vec![walls[1]; partial.measured_links()],
                &mut work,
            )
            .is_err()
        );
        let nonuniform: Vec<_> = probe
            .wall_links()
            .map(|link| [0.001 * (1.0 + link.fluid_cell[0] as f64), 0.0, 0.0])
            .collect();
        assert!(
            grid.step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &nonuniform,
                &mut work,
            )
            .is_err()
        );
        let gate = CancelGate::new_clock_free();
        let refusing = RefuseAfterStreaming {
            pool: TilePool::for_host(2, 503),
            calls: AtomicUsize::new(0),
        };
        assert!(
            grid.step_pooled_with_channel_lagrangian_sgs_and_moving_walls(
                &refusing,
                &gate,
                config(),
                None,
                &probe,
                &velocity,
                &mut work,
            )
            .is_err()
        );
        assert_eq!(refusing.calls.load(Ordering::Acquire), 2);
        assert_eq!(grid.populations(), before);
        assert_eq!(grid.completed_steps(), 1);
        assert_eq!(work.accepted_fit_time(), Some(0.0));
        let retry = grid
            .step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &velocity,
                &mut work,
            )
            .unwrap();
        let expected = control
            .step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &velocity,
                &mut clean,
            )
            .unwrap();
        assert_eq!(retry, expected);
        assert_eq!(grid.populations(), control.populations());
        assert_eq!(work.accepted_fit_time(), Some(1.0));
    });
}

#[test]
fn late_native_mach_refusal_does_not_publish_a_prepared_wall_history() {
    let mut grid = flow::<D3q19, 19>(
        CollisionModel3::Bgk { tau: 0.8 },
        equilibrium3,
        Some(2),
        [0.002, 0.0, 0.0],
        [[0.0; 3]; 2],
        false,
    );
    grid.budget.maximum_mach = 0.001;
    let mut work = grid.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
    let before = grid.populations().to_vec();
    with_cx(|cx| {
        let (probe, velocity) = drive(&grid, cx, 2, [[0.0003, 0.0, 0.0]; 2]);
        let (trial, _) = grid
            .prepare_lagrangian_sgs_with_motion(
                config(),
                None,
                &mut work.local,
                &mut work.filter,
                None,
                Some(WallMotion3 {
                    probe: &probe,
                    velocities: &velocity,
                }),
                || Ok(()),
            )
            .unwrap();
        drop(trial);
        assert!(matches!(
            grid.step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &velocity,
                &mut work,
            ),
            Err(SgsStepError::Grid(GridError3::Mach { .. }))
        ));
    });
    assert_eq!(grid.populations(), before);
    assert_eq!(grid.completed_steps(), 0);
    assert_eq!(work.accepted_fit_time(), None);
}
