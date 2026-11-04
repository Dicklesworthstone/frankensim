//! G0/G1/G3/G4/G5: actual finite-wall drag, inertia, forcing and native transport.

use fs_exec::{CancelGate, Cx, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{CollisionModel3, Face3, equilibrium3};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, equilibrium27};
use fs_lbm::lattice3::finite_wall::{FinitePlanarWall3, FiniteWallBudget3, PlanarWallState3};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, GridError3, Lattice3};
use std::sync::atomic::{AtomicUsize, Ordering};

const SEED: u64 = 0xf117_ea11;
const DIMS: [usize; 3] = [4, 5, 6];
const RHO: f64 = 1.3;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: SEED,
                kernel_id: 132,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
fn budget() -> FiniteWallBudget3 {
    FiniteWallBudget3 {
        maximum_relative_impulse_defect: 2e-10,
        maximum_relative_work_defect: 2e-10,
    }
}
fn state() -> PlanarWallState3 {
    PlanarWallState3 {
        mass: 3.0,
        axial_inertia: 9.0,
        velocity: [0.009, -0.012, 0.0013],
    }
}
fn center(face: Face3) -> [f64; 3] {
    let n = face as usize / 2;
    let mut c = DIMS.map(|d| d as f64 / 2.0);
    c[n] = if (face as usize).is_multiple_of(2) {
        0.0
    } else {
        DIMS[n] as f64
    };
    c
}
fn grid<L: Lattice3<Q>, const Q: usize>(
    face: Face3,
    collision: L::Collision,
    force: [f64; 3],
    equilibrium: impl Fn(f64, [f64; 3]) -> [f64; Q],
) -> Grid3<L, Q> {
    let mut periodic = [true; 3];
    periodic[face as usize / 2] = false;
    Grid3::try_new(
        DIMS,
        collision,
        force,
        periodic,
        GridBudget3 {
            maximum_cells: DIMS.iter().product(),
            maximum_mach: 0.2,
        },
        |_| false,
        |_| equilibrium(RHO, [0.0; 3]),
    )
    .unwrap()
}
fn d27() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    }
}
fn bits<L: Lattice3<Q>, const Q: usize>(grid: &Grid3<L, Q>) -> Vec<u64> {
    grid.populations()
        .iter()
        .flatten()
        .map(|x| x.to_bits())
        .collect()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} at {tolerance:e}"
    );
}
fn populations<L: Lattice3<Q>, const Q: usize>(grid: &Grid3<L, Q>) -> [f64; 4] {
    let mut out = [0.0; 4];
    for cell in grid.populations() {
        for (q, &f) in cell.iter().enumerate() {
            out[0] += f;
            for a in 0..3 {
                out[a + 1] += f * f64::from(L::E[q][a]);
            }
        }
    }
    out
}

// Independent half-space lattice moments: each complete face cell has
// tangential resistance rho/3; a centered rectangular face has
// mean(r²) = (Na² + Nb² - 2)/12. No production response assembly is used.
fn first_step<L: Lattice3<Q>, const Q: usize>(
    face: Face3,
    collision: L::Collision,
    equilibrium: impl Fn(f64, [f64; 3]) -> [f64; Q],
) {
    with_cx(|cx| {
        let n = face as usize / 2;
        let a = (n + 1) % 3;
        let b = (n + 2) % 3;
        let mut fluid = grid::<L, Q>(face, collision, [0.0; 3], equilibrium);
        let before = populations(&fluid);
        let initial = state();
        let mut wall =
            FinitePlanarWall3::new(&fluid, cx, face, center(face), initial, |_| true).unwrap();
        let drive = [0.012, 0.007, -0.0008];
        let area = (DIMS[a] * DIMS[b]) as f64;
        let k = [
            RHO * area / 3.0,
            RHO * area / 3.0,
            RHO * area * ((DIMS[a] * DIMS[a] + DIMS[b] * DIMS[b] - 2) as f64) / 36.0,
        ];
        let masses = [initial.mass, initial.mass, initial.axial_inertia];
        let report = fluid
            .step_serial_with_finite_wall(cx, &mut wall, drive, budget(), None)
            .unwrap();
        for i in 0..3 {
            let middle =
                (2.0 * masses[i] * initial.velocity[i] + drive[i]) / (2.0 * masses[i] + k[i]);
            near(report.velocity_midpoint[i], middle, 2e-15);
            near(
                report.after.velocity[i],
                2.0 * middle - initial.velocity[i],
                4e-15,
            );
            let measured = if i < 2 {
                report.hydrodynamic.wall_impulse[[a, b][i]]
            } else {
                report.hydrodynamic.wall_angular_impulse[n]
            };
            near(measured, -k[i] * middle, 2e-13);
            near(report.momentum_residual[i], 0.0, 2e-13);
        }
        let expected_work = -(0..3)
            .map(|i| k[i] * report.velocity_midpoint[i].powi(2))
            .sum::<f64>();
        near(report.hydrodynamic.wall_work, expected_work, 2e-15);
        near(report.work_residual, 0.0, 2e-15);
        assert!(report.kinetic_energy_change < 0.0);
        let after = populations(&fluid);
        near(after[0], before[0], 3e-11);
        for axis in 0..3 {
            near(
                after[axis + 1] - before[axis + 1]
                    + report.hydrodynamic.wall_impulse[axis]
                    + report.stationary_walls.wall_impulse[axis],
                0.0,
                4e-12,
            );
            near(
                report.hydrodynamic.wall_angular_impulse[axis]
                    + report.hydrodynamic.fluid_population_angular_impulse[axis]
                    - report.hydrodynamic.wall_velocity_mass_angular_impulse[axis],
                0.0,
                2e-12,
            );
        }
        assert_eq!(wall.state(), report.after);
        assert_eq!(fluid.completed_steps(), 1);
    });
}

#[test]
fn finite_inertia_and_drive_match_independent_half_space_resistance_on_all_faces() {
    for face in Face3::ALL {
        first_step::<D3q19, 19>(face, CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
        first_step::<D3q27, 27>(face, d27(), |rho, u| {
            equilibrium27(rho, u, Equilibrium27::Factorized)
        });
    }
    first_step::<D3q19, 19>(
        Face3::YMax,
        CollisionModel3::Trt {
            tau: 1.1,
            magic: 3.0 / 16.0,
        },
        equilibrium3,
    );
    first_step::<D3q27, 27>(
        Face3::XMin,
        CollisionConfig27 {
            model: CollisionModel27::Trt {
                tau: 1.2,
                magic: 3.0 / 16.0,
                equilibrium: Equilibrium27::Quadratic,
            },
            maximum_mach: 0.2,
        },
        |rho, u| equilibrium27(rho, u, Equilibrium27::Quadratic),
    );
}

#[test]
fn changing_drive_and_viscosity_use_the_same_native_stream_and_close_body_work() {
    with_cx(|cx| {
        let face = Face3::YMax;
        let c = center(face);
        let force = [2e-6, 0.0, -3e-6];
        let make = || {
            grid::<D3q27, 27>(face, d27(), force, |rho, u| {
                equilibrium27(rho, u, Equilibrium27::Factorized)
            })
        };
        let mut fluid = make();
        let mut reference = make();
        let mut wall = FinitePlanarWall3::new(&fluid, cx, face, c, state(), |_| true).unwrap();
        let probe = reference
            .wall_momentum_probe(cx, c, |link| link.exterior_faces[face as usize])
            .unwrap();
        let collisions: Vec<_> = (0..DIMS.iter().product())
            .map(|i| CollisionConfig27 {
                model: CollisionModel27::Bgk {
                    tau: 0.75 + 0.02 * (i % 11) as f64,
                    equilibrium: Equilibrium27::Factorized,
                },
                maximum_mach: 0.2,
            })
            .collect();
        let mut works = [0.0; 3];
        for tick in 0..24 {
            let before = populations(&fluid);
            let drive = if tick < 12 {
                [0.002, -0.003, 0.0007]
            } else {
                [0.0; 3]
            };
            let receipt = fluid
                .step_serial_with_finite_wall(cx, &mut wall, drive, budget(), Some(&collisions))
                .unwrap();
            let v = receipt.velocity_midpoint;
            let velocities: Vec<_> = probe
                .wall_links()
                .map(|link| {
                    let p = link.fluid_cell;
                    [
                        v[1] + v[2] * (p[2] as f64 + 0.5 - c[2]),
                        0.0,
                        v[0] - v[2] * (p[0] as f64 + 0.5 - c[0]),
                    ]
                })
                .collect();
            let expected = reference
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, Some(&collisions))
                .unwrap();
            assert_eq!(receipt.hydrodynamic, expected);
            assert_eq!(bits(&fluid), bits(&reference));
            let after = populations(&fluid);
            near(before[0], after[0], 3e-11);
            let n = face as usize / 2;
            let mut wall_impulse = [0.0; 3];
            for i in 0..2 {
                wall_impulse[(n + 1 + i) % 3] =
                    receipt.after.mass * (receipt.after.velocity[i] - receipt.before.velocity[i]);
            }
            let mut external = [0.0; 3];
            external[(n + 1) % 3] = drive[0];
            external[(n + 2) % 3] = drive[1];
            for a in 0..3 {
                near(
                    after[a + 1] - before[a + 1]
                        + wall_impulse[a]
                        + receipt.stationary_walls.wall_impulse[a]
                        - receipt.support_impulse[a]
                        - external[a]
                        - receipt.fluid_force_impulse[a],
                    0.0,
                    3e-11,
                );
            }
            works[0] += receipt.kinetic_energy_change;
            works[1] += receipt.hydrodynamic.wall_work;
            works[2] += receipt.external_work;
        }
        near(works[0], works[1] + works[2], 2e-13);
    });
}

#[test]
fn off_center_window_couples_translation_and_spin_through_actual_traction() {
    with_cx(|cx| {
        let face = Face3::ZMin;
        let mut fluid = grid::<D3q27, 27>(face, d27(), [0.0; 3], |rho, u| {
            equilibrium27(rho, u, Equilibrium27::Factorized)
        });
        let c = [0.0, 0.0, 0.0];
        let initial = PlanarWallState3 {
            velocity: [0.01, 0.0, 0.0],
            ..state()
        };
        let mut wall =
            FinitePlanarWall3::new(&fluid, cx, face, c, initial, |p| p[0] < 2 && p[1] < 3).unwrap();
        let report = fluid
            .step_serial_with_finite_wall(cx, &mut wall, [0.0; 3], budget(), None)
            .unwrap();
        assert!(report.after.velocity[2] > 0.0);
        let v = report.velocity_midpoint;
        let k = RHO * 6.0 / 3.0;
        let (mx, my, r2) = (1.0, 1.5, 1.25 + 35.0 / 12.0);
        let expected = [
            -k * (v[0] - my * v[2]),
            -k * (v[1] + mx * v[2]),
            -k * (-my * v[0] + mx * v[1] + r2 * v[2]),
        ];
        for i in 0..3 {
            let actual = if i < 2 {
                report.hydrodynamic.wall_impulse[i]
            } else {
                report.hydrodynamic.wall_angular_impulse[2]
            };
            near(actual, expected[i], 3e-13);
        }
        near(report.work_residual, 0.0, 2e-15);
    });
}

#[test]
fn serial_and_pooled_finite_wall_response_are_bit_identical() {
    with_cx(|cx| {
        let face = Face3::ZMax;
        let make = || {
            grid::<D3q27, 27>(face, d27(), [0.0; 3], |rho, u| {
                equilibrium27(rho, u, Equilibrium27::Factorized)
            })
        };
        let mut serial = make();
        let mut source =
            FinitePlanarWall3::new(&serial, cx, face, center(face), state(), |_| true).unwrap();
        let mut expected = Vec::new();
        for _ in 0..12 {
            expected.push(
                serial
                    .step_serial_with_finite_wall(
                        cx,
                        &mut source,
                        [0.002, 0.001, -0.0003],
                        budget(),
                        None,
                    )
                    .unwrap(),
            );
        }
        for workers in [1, 2, 5] {
            let mut fluid = make();
            let mut wall =
                FinitePlanarWall3::new(&fluid, cx, face, center(face), state(), |_| true).unwrap();
            let pool = TilePool::for_host(workers, SEED);
            let gate = CancelGate::new_clock_free();
            for receipt in &expected {
                let actual = fluid
                    .step_pooled_with_finite_wall(
                        &pool,
                        &gate,
                        &mut wall,
                        [0.002, 0.001, -0.0003],
                        budget(),
                        None,
                    )
                    .unwrap();
                assert_eq!(actual, *receipt);
            }
            assert_eq!(bits(&fluid), bits(&serial));
            assert_eq!(wall.state(), source.state());
        }
    });
}

#[test]
fn late_cancellation_and_inadmissible_drives_publish_neither_wall_nor_fluid() {
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
            let out = self.pool.run_with_gate(kernel, gate);
            if self.calls.fetch_add(1, Ordering::AcqRel) == 1 {
                gate.request();
            }
            out
        }
    }
    with_cx(|cx| {
        let face = Face3::YMin;
        let make = || {
            grid::<D3q27, 27>(face, d27(), [0.0; 3], |rho, u| {
                equilibrium27(rho, u, Equilibrium27::Factorized)
            })
        };
        let mut fluid = make();
        let mut wall =
            FinitePlanarWall3::new(&fluid, cx, face, center(face), state(), |_| true).unwrap();
        let saved = bits(&fluid);
        let saved_wall = wall.state();
        for drive in [[f64::NAN, 0.0, 0.0], [1e8, 0.0, 0.0]] {
            assert!(
                fluid
                    .step_serial_with_finite_wall(cx, &mut wall, drive, budget(), None)
                    .is_err()
            );
            assert_eq!(bits(&fluid), saved);
            assert_eq!(wall.state(), saved_wall);
        }
        let pool = CancelAfterStream {
            pool: TilePool::for_host(2, SEED),
            calls: AtomicUsize::new(0),
        };
        let gate = CancelGate::new_clock_free();
        assert!(matches!(
            fluid.step_pooled_with_finite_wall(&pool, &gate, &mut wall, [0.001; 3], budget(), None),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(bits(&fluid), saved);
        assert_eq!(wall.state(), saved_wall);
        assert_eq!(fluid.completed_steps(), 0);
        let actual = fluid
            .step_serial_with_finite_wall(cx, &mut wall, [0.001; 3], budget(), None)
            .unwrap();
        let mut reference = make();
        let mut expected_wall =
            FinitePlanarWall3::new(&reference, cx, face, center(face), state(), |_| true).unwrap();
        let expected = reference
            .step_serial_with_finite_wall(cx, &mut expected_wall, [0.001; 3], budget(), None)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(bits(&fluid), bits(&reference));
    });
}
