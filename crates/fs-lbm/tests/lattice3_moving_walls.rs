//! G0/G1/G3/G4/G5: native moving-wall transport, momentum, work and rollback.

use fs_exec::{CancelGate, Cx, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{CollisionModel3, equilibrium3};
use fs_lbm::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, GridError3, Lattice3};
use std::sync::atomic::{AtomicUsize, Ordering};

const SEED: u64 = 0xD327_B0A7;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: SEED,
                kernel_id: 131,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn bits<L: Lattice3<Q>, const Q: usize>(grid: &Grid3<L, Q>) -> Vec<u64> {
    grid.populations()
        .iter()
        .flatten()
        .map(|v| v.to_bits())
        .collect()
}

fn config(model: CollisionModel27) -> CollisionConfig27 {
    CollisionConfig27 {
        model,
        maximum_mach: 0.2,
    }
}

fn bgk(tau: f64) -> CollisionConfig27 {
    config(CollisionModel27::Bgk {
        tau,
        equilibrium: Equilibrium27::Factorized,
    })
}

fn couette<L: Lattice3<Q>, const Q: usize>(
    collision: L::Collision,
    tau: f64,
    tangent: usize,
    speed: f64,
    equilibrium: impl Fn(f64, [f64; 3]) -> [f64; Q],
) {
    let height = 6;
    let normal = (tangent + 1) % 3;
    let rho = 1.4;
    let mut dims = [1; 3];
    dims[normal] = height;
    let mut periodic = [true; 3];
    periodic[normal] = false;
    let mut grid = Grid3::<L, Q>::try_new(
        dims,
        collision,
        [0.0; 3],
        periodic,
        GridBudget3 {
            maximum_cells: height,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| equilibrium(rho, [0.0; 3]),
    )
    .unwrap();
    with_cx(|cx| {
        // Only the driven wall is selected; the opposite wall remains native
        // stationary bounce-back and is excluded from this wall's load.
        let probe = grid
            .wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[2 * normal + 1])
            .unwrap();
        let mut u_wall = [0.0; 3];
        u_wall[tangent] = speed;
        let velocities = vec![u_wall; probe.measured_links()];
        let first = grid
            .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
            .unwrap();
        // Independent half-space lattice moment: sum(w e_t^2)=1/18 on either
        // D3Q19 or D3Q27. A wall first acting on resting fluid receives -rho*u/3.
        assert!((first.wall_impulse[tangent] + rho * speed / 3.0).abs() < 1e-13);
        assert!((first.wall_work + rho * speed * speed / 3.0).abs() < 1e-14);
        let nu = (tau - 0.5) / 3.0;
        let steps = (20.0 * (height * height) as f64 / (std::f64::consts::PI.powi(2) * nu)).ceil()
            as usize
            + 80;
        let mut load = first;
        for _ in 1..steps {
            load = grid
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
                .unwrap();
        }
        for cell in 0..height {
            let (density, velocity) = grid.macroscopic(cell).unwrap();
            let expected = speed * (cell as f64 + 0.5) / height as f64;
            assert!((velocity[tangent] - expected).abs() < 2e-8 * speed.abs());
            assert!((density - rho).abs() < 1e-10);
            for a in 0..3 {
                if a != tangent {
                    assert!(velocity[a].abs() < 1e-12);
                }
            }
        }
        let shear = rho * nu * speed / height as f64; // area = 1
        assert!((load.wall_impulse[tangent] + shear).abs() < 2e-8 * shear.abs());
        assert!((load.wall_work + shear * speed).abs() < 2e-8 * (shear * speed).abs());
        assert!(load.wall_work < 0.0);
        assert!(load.fluid_mass_change.abs() < 1e-13);
        for a in 0..3 {
            assert!(
                (load.wall_impulse[a] + load.fluid_population_impulse[a]
                    - load.wall_velocity_mass_impulse[a])
                    .abs()
                    < 1e-13
            );
        }
    });
}

#[test]
fn couette_profile_viscous_drag_and_boundary_work_match_the_analytic_solution() {
    couette::<D3q19, 19>(
        CollisionModel3::Bgk { tau: 0.8 },
        0.8,
        0,
        0.015,
        equilibrium3,
    );
    couette::<D3q19, 19>(
        CollisionModel3::Trt {
            tau: 1.2,
            magic: 3.0 / 16.0,
        },
        1.2,
        2,
        -0.015,
        equilibrium3,
    );
    for (axis, model, tau) in [
        (
            0,
            CollisionModel27::Bgk {
                tau: 0.8,
                equilibrium: Equilibrium27::Quadratic,
            },
            0.8,
        ),
        (
            1,
            CollisionModel27::Bgk {
                tau: 1.2,
                equilibrium: Equilibrium27::Factorized,
            },
            1.2,
        ),
        (
            2,
            CollisionModel27::Trt {
                tau: 0.8,
                magic: 3.0 / 16.0,
                equilibrium: Equilibrium27::Factorized,
            },
            0.8,
        ),
        (
            0,
            CollisionModel27::CentralMoment {
                rates: Rates27 {
                    omega: [1.25, 1.0, 0.8, 0.9, 1.0, 1.1, 1.2, 0.7, 0.8, 0.9],
                },
            },
            0.8,
        ),
        (
            2,
            CollisionModel27::Cumulant {
                rates: Rates27 {
                    omega: [1.25, 1.0, 0.8, 0.9, 1.0, 1.1, 1.2, 0.7, 0.8, 0.9],
                },
            },
            0.8,
        ),
    ] {
        couette::<D3q27, 27>(
            config(model),
            tau,
            axis,
            if axis == 1 { -0.015 } else { 0.015 },
            |r, u| equilibrium27(r, u, Equilibrium27::Factorized),
        );
    }
}

#[test]
fn co_moving_fluid_and_surfaces_do_not_create_drag_or_wall_power() {
    with_cx(|cx| {
        let velocity = [0.04, 0.0, 0.0];
        let mut grid = Grid3::<D3q27, 27>::try_new(
            [2, 6, 2],
            bgk(0.8),
            [0.0; 3],
            [true, false, true],
            GridBudget3 {
                maximum_cells: 24,
                maximum_mach: 0.2,
            },
            |_| false,
            |_| equilibrium27(1.3, velocity, Equilibrium27::Factorized),
        )
        .unwrap();
        let probe = grid.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
        let velocities = vec![velocity; probe.measured_links()];
        for _ in 0..30 {
            let load = grid
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
                .unwrap();
            assert!(load.wall_impulse[0].abs() < 1e-13);
            assert!(load.wall_work.abs() < 1e-14);
        }
        for cell in 0..24 {
            let (rho, u) = grid.macroscopic(cell).unwrap();
            assert!((rho - 1.3).abs() < 1e-12);
            for a in 0..3 {
                assert!((u[a] - velocity[a]).abs() < 1e-13);
            }
        }
    });
}

fn fixture() -> Grid3<D3q27, 27> {
    Grid3::try_new(
        [7, 6, 3],
        bgk(0.8),
        [2e-6, 0.0, -1e-6],
        [true, false, true],
        GridBudget3 {
            maximum_cells: 126,
            maximum_mach: 0.2,
        },
        |p| p[1] == 0,
        |p| {
            equilibrium27(
                1.1 + 1e-4 * p[0] as f64,
                [0.003, 0.0, -0.002],
                Equilibrium27::Factorized,
            )
        },
    )
    .unwrap()
}

fn population_moments(grid: &Grid3<D3q27, 27>) -> [f64; 4] {
    let mut value = [0.0; 4];
    for (cell, f) in grid.populations().iter().enumerate() {
        if (cell / 7) % 6 == 0 {
            continue;
        }
        for (q, &population) in f.iter().enumerate() {
            value[0] += population;
            for a in 0..3 {
                value[a + 1] += f64::from(D3q27::E[q][a]) * population;
            }
        }
    }
    value
}

#[test]
fn changing_surface_drive_and_material_collision_field_preserve_population_balances() {
    with_cx(|cx| {
        let mut grid = fixture();
        let mut ordinary = fixture();
        let probe = grid
            .wall_momentum_probe(cx, [2.1, 0.2, -0.4], |_| true)
            .unwrap();
        let collisions: Vec<_> = (0..126).map(|i| bgk(0.7 + 0.04 * (i % 7) as f64)).collect();
        for step in 0..8 {
            let velocities: Vec<_> = probe
                .wall_links()
                .map(|link| {
                    if link.solid_cell.is_some() {
                        [-0.003, 0.0, 0.001]
                    } else {
                        [0.004 + 0.001 * step as f64, 0.0, -0.002]
                    }
                })
                .collect();
            let before = population_moments(&grid);
            let load = grid
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, Some(&collisions))
                .unwrap();
            ordinary
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, Some(&collisions))
                .unwrap();
            let after = population_moments(&grid);
            assert!((after[0] - before[0] - load.fluid_mass_change).abs() < 2e-11);
            assert!(load.fluid_mass_change.abs() < 1e-13);
            for (a, force) in [2e-6, 0.0, -1e-6].into_iter().enumerate() {
                assert!(
                    (after[a + 1]
                        - before[a + 1]
                        - 105.0 * force
                        - load.fluid_population_impulse[a])
                        .abs()
                        < 2e-12
                );
                assert!(
                    (load.wall_impulse[a] + load.fluid_population_impulse[a]
                        - load.wall_velocity_mass_impulse[a])
                        .abs()
                        < 2e-13
                );
                assert!(
                    (load.wall_angular_impulse[a] + load.fluid_population_angular_impulse[a]
                        - load.wall_velocity_mass_angular_impulse[a])
                        .abs()
                        < 2e-12
                );
            }
        }
        // Neither drive nor local viscosity is latched into the next step.
        let zero = vec![[0.0; 3]; probe.measured_links()];
        grid.step_serial_with_moving_wall_momentum(cx, &probe, &zero, None)
            .unwrap();
        ordinary.step_serial(cx).unwrap();
        assert_eq!(bits(&grid), bits(&ordinary));
    });
}

#[test]
fn zero_wall_motion_and_stationary_measurements_retain_identical_bits() {
    with_cx(|cx| {
        let mut stationary = fixture();
        let mut moving = fixture();
        let probe = moving
            .wall_momentum_probe(cx, [1.1, -0.7, 0.3], |_| true)
            .unwrap();
        let zero = vec![[0.0; 3]; probe.measured_links()];
        for _ in 0..5 {
            let expected = stationary
                .step_serial_with_wall_momentum(cx, &probe, None)
                .unwrap();
            let actual = moving
                .step_serial_with_moving_wall_momentum(cx, &probe, &zero, None)
                .unwrap();
            assert_eq!(bits(&moving), bits(&stationary));
            assert_eq!(
                actual.wall_impulse.map(f64::to_bits),
                expected.wall_impulse.map(f64::to_bits)
            );
            assert_eq!(
                actual.wall_angular_impulse.map(f64::to_bits),
                expected.wall_angular_impulse.map(f64::to_bits)
            );
            assert_eq!(actual.wall_work, 0.0);
            assert_eq!(actual.fluid_mass_change, 0.0);
        }
    });
}

#[test]
fn moving_wall_tiles_and_loads_are_identical_across_worker_counts() {
    with_cx(|cx| {
        let mut serial = fixture();
        let probe = serial
            .wall_momentum_probe(cx, [0.5, 0.1, 0.7], |_| true)
            .unwrap();
        let velocities: Vec<_> = probe
            .wall_links()
            .map(|link| {
                if link.solid_cell.is_some() {
                    [0.0; 3]
                } else {
                    [0.01, 0.0, -0.004]
                }
            })
            .collect();
        let laws: Vec<_> = (0..126).map(|i| bgk(0.7 + 0.04 * (i % 7) as f64)).collect();
        let mut receipts = Vec::new();
        for _ in 0..4 {
            receipts.push(
                serial
                    .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, Some(&laws))
                    .unwrap(),
            );
        }
        for workers in [1, 2, 7] {
            let pool = TilePool::for_host(workers, SEED);
            let gate = CancelGate::new_clock_free();
            let mut grid = fixture();
            for expected in &receipts {
                let load = grid
                    .step_pooled_with_moving_wall_momentum(
                        &pool,
                        &gate,
                        &probe,
                        &velocities,
                        Some(&laws),
                    )
                    .unwrap();
                assert_eq!(load, *expected);
            }
            assert_eq!(bits(&grid), bits(&serial));
        }
    });
}

#[test]
fn fixed_mask_refuses_normal_motion_and_incomplete_stencils_without_publication() {
    with_cx(|cx| {
        let mut grid = fixture();
        let before = bits(&grid);
        let probe = grid
            .wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[3])
            .unwrap();
        let count = probe.measured_links();
        for velocities in [
            vec![[0.0; 3]; count - 1],
            vec![[f64::NAN, 0.0, 0.0]; count],
            vec![[0.2, 0.0, 0.0]; count],
            vec![[0.0, 0.01, 0.0]; count],
        ] {
            assert!(
                grid.step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
                    .is_err()
            );
            assert_eq!(bits(&grid), before);
        }
        let incomplete = grid
            .wall_momentum_probe(cx, [0.0; 3], |link| {
                link.exterior_faces[3] && link.incoming_direction == [1, -1, 0]
            })
            .unwrap();
        let velocities = vec![[0.01, 0.0, 0.0]; incomplete.measured_links()];
        assert!(matches!(
            grid.step_serial_with_moving_wall_momentum(cx, &incomplete, &velocities, None),
            Err(GridError3::Input(
                "fixed moving-wall links must conserve local mass"
            ))
        ));
        assert_eq!(bits(&grid), before);
        assert_eq!(grid.completed_steps(), 0);
    });
}

#[test]
fn cancellation_after_moving_stream_keeps_state_and_retry_exact() {
    struct CancelAfterStream {
        pool: TilePool,
        calls: AtomicUsize,
    }
    impl KernelRunner for CancelAfterStream {
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
    with_cx(|cx| {
        let mut grid = fixture();
        let before = bits(&grid);
        let probe = grid.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
        let velocities: Vec<_> = probe
            .wall_links()
            .map(|link| {
                if link.solid_cell.is_some() {
                    [0.0; 3]
                } else {
                    [0.01, 0.0, 0.0]
                }
            })
            .collect();
        let runner = CancelAfterStream {
            pool: TilePool::for_host(2, SEED),
            calls: AtomicUsize::new(0),
        };
        let gate = CancelGate::new_clock_free();
        assert!(matches!(
            grid.step_pooled_with_moving_wall_momentum(&runner, &gate, &probe, &velocities, None),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(bits(&grid), before);
        assert_eq!(grid.completed_steps(), 0);
        let actual = grid
            .step_pooled_with_moving_wall_momentum(
                &runner.pool,
                &CancelGate::new_clock_free(),
                &probe,
                &velocities,
                None,
            )
            .unwrap();
        let mut reference = fixture();
        let expected = reference
            .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, None)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(bits(&grid), bits(&reference));
    });
}
