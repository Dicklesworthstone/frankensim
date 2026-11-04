//! Native D3Q19/D3Q27 porous-source momentum, work, transport and rollback.

use fs_exec::{CancelGate, Cx, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{CollisionModel3, equilibrium3};
use fs_lbm::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use fs_lbm::lattice3::sources::resistance::{
    DarcyForchheimer3, LinearResistance3, ResistanceGrid3,
};
use fs_lbm::lattice3::sources::{ForcedGrid3, ForcedLattice3};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, GridError3};
use std::sync::atomic::{AtomicUsize, Ordering};

const SEED: u64 = 0xD327_DA2C;
const AXES: [[f64; 3]; 3] = [
    [2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    [-2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0],
];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
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

fn models() -> [CollisionConfig27; 6] {
    let rates = Rates27 {
        omega: [1.25, 1.1, 0.7, 1.3, 0.85, 1.0, 1.4, 0.9, 1.1, 1.2],
    };
    [
        CollisionModel27::Bgk {
            tau: 0.8,
            equilibrium: Equilibrium27::Quadratic,
        },
        CollisionModel27::Bgk {
            tau: 0.8,
            equilibrium: Equilibrium27::Factorized,
        },
        CollisionModel27::Trt {
            tau: 0.8,
            magic: 3.0 / 16.0,
            equilibrium: Equilibrium27::Quadratic,
        },
        CollisionModel27::Trt {
            tau: 0.8,
            magic: 3.0 / 16.0,
            equilibrium: Equilibrium27::Factorized,
        },
        CollisionModel27::CentralMoment { rates },
        CollisionModel27::Cumulant { rates },
    ]
    .map(|model| CollisionConfig27 {
        model,
        maximum_mach: 0.2,
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:e}"
    );
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn cell(index: usize, dims: [usize; 3]) -> [usize; 3] {
    [
        index % dims[0],
        (index / dims[0]) % dims[1],
        index / (dims[0] * dims[1]),
    ]
}
fn force(law: DarcyForchheimer3, velocity: [f64; 3]) -> [f64; 3] {
    let relative = core::array::from_fn(|a| law.linear.skeleton_velocity[a] - velocity[a]);
    let speed = dot(relative, relative).sqrt();
    core::array::from_fn(|a| {
        law.linear.body_force[a]
            + law.quadratic * speed * relative[a]
            + (0..3)
                .map(|k| {
                    law.linear.axes[k][a]
                        * law.linear.principal[k]
                        * dot(law.linear.axes[k], relative)
                })
                .sum::<f64>()
    })
}
fn snapshot<L: ForcedLattice3<Q>, const Q: usize>(
    grid: &ResistanceGrid3<L, Q>,
) -> (Vec<u64>, Vec<u64>, usize) {
    (
        grid.flow()
            .populations()
            .iter()
            .flatten()
            .map(|v| v.to_bits())
            .collect(),
        grid.flow()
            .forces()
            .iter()
            .flatten()
            .map(|v| v.to_bits())
            .collect(),
        grid.flow().completed_steps(),
    )
}
fn uniform<L: ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    pop: [f64; Q],
    law: DarcyForchheimer3,
) -> ResistanceGrid3<L, Q> {
    let grid = Grid3::try_new(
        [2, 2, 2],
        collision,
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: 8,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| pop,
    )
    .unwrap();
    with_cx(|cx| ResistanceGrid3::try_new(grid, cx, vec![law; 8], [0.0; 3]).unwrap())
}

fn linear_oracle<L: ForcedLattice3<Q>, const Q: usize>(collision: L::Collision, pop: [f64; Q]) {
    let law = LinearResistance3 {
        principal: [0.02, 0.08, 0.2],
        axes: AXES,
        skeleton_velocity: [0.003, -0.001, 0.002],
        body_force: [4e-5, -3e-5, 2e-5],
    };
    let mut grid = uniform::<L, Q>(collision, pop, law.into());
    assert_eq!(
        grid.flow().populations()[0].map(f64::to_bits),
        pop.map(f64::to_bits)
    );
    let (rho, mut expected) = grid.flow().macroscopic(0).unwrap();
    for _ in 0..5 {
        let previous = expected;
        let relative = core::array::from_fn(|a| previous[a] - law.skeleton_velocity[a]);
        let principal: [f64; 3] = core::array::from_fn(|k| {
            ((rho - 0.5 * law.principal[k]) * dot(law.axes[k], relative)
                + dot(law.axes[k], law.body_force))
                / (rho + 0.5 * law.principal[k])
        });
        expected = core::array::from_fn(|a| {
            law.skeleton_velocity[a] + (0..3).map(|k| law.axes[k][a] * principal[k]).sum::<f64>()
        });
        let receipt = with_cx(|cx| grid.step_serial(cx, None).unwrap());
        let (density, actual) = grid.flow().macroscopic(0).unwrap();
        close(density, rho, 4e-14);
        for a in 0..3 {
            close(actual[a], expected[a], 2e-15);
            close(
                grid.flow().forces()[0][a],
                force(law.into(), actual)[a],
                3e-15,
            );
            close(
                8.0 * rho * (actual[a] - previous[a]),
                receipt.external_impulse[a] - receipt.skeleton_impulse[a],
                4e-14,
            );
        }
        close(
            receipt.fluid_drag_work + receipt.skeleton_work,
            -receipt.relative_dissipation,
            1e-16,
        );
        assert!(receipt.relative_dissipation >= 0.0);
    }
}

#[test]
fn rotated_linear_relaxation_reuses_accepted_force_for_every_native_collision_family() {
    for model in models() {
        linear_oracle::<D3q27, 27>(
            model,
            equilibrium27(1.3, [0.011, -0.007, 0.005], Equilibrium27::Factorized),
        );
    }
    for model in [
        CollisionModel3::Bgk { tau: 0.8 },
        CollisionModel3::CentralMoment {
            second_order_rate: 1.25,
            higher_order_rate: 0.9,
        },
    ] {
        linear_oracle::<D3q19, 19>(model, equilibrium3(1.3, [0.011, -0.007, 0.005]));
    }
}

#[test]
fn nonlinear_decay_is_second_order_and_resolved_stiff_reversal_dissipates() {
    let initial = [0.02, -0.01, 0.005];
    let speed = dot(initial, initial).sqrt();
    let coefficient = 30.0;
    let errors = [8, 16, 32].map(|steps| {
        let law = DarcyForchheimer3 {
            linear: LinearResistance3::default(),
            quadratic: coefficient / steps as f64,
        };
        let first_force = force(law, initial);
        let raw = core::array::from_fn(|a| initial[a] - 0.5 * first_force[a]);
        let mut grid = uniform::<D3q27, 27>(
            models()[5],
            equilibrium27(1.0, raw, Equilibrium27::Factorized),
            law,
        );
        for a in 0..3 {
            close(grid.flow().macroscopic(0).unwrap().1[a], initial[a], 1e-16);
        }
        for _ in 0..steps {
            let _ = with_cx(|cx| grid.step_serial(cx, None).unwrap());
        }
        let actual = grid.flow().macroscopic(0).unwrap().1;
        (actual[0] - initial[0] / (1.0 + coefficient * speed)).abs()
    });
    for pair in errors.windows(2) {
        assert!((pair[0] / pair[1] - 4.0).abs() < 0.03, "{errors:?}");
    }

    for (principal, quadratic) in [([500.0; 3], 0.0), ([0.01, 0.4, 20.0], 4000.0)] {
        let law = DarcyForchheimer3 {
            linear: LinearResistance3 {
                principal,
                axes: AXES,
                ..LinearResistance3::default()
            },
            quadratic,
        };
        let mut grid = uniform::<D3q27, 27>(
            models()[1],
            equilibrium27(1.0, initial, Equilibrium27::Factorized),
            law,
        );
        let mut old = grid.flow().macroscopic(0).unwrap().1;
        for _ in 0..5 {
            let receipt = with_cx(|cx| grid.step_serial(cx, None).unwrap());
            let new = grid.flow().macroscopic(0).unwrap().1;
            assert!(dot(new, new) <= dot(old, old) * (1.0 + 1e-10));
            if quadratic == 0.0 {
                assert!(dot(new, old) < 0.0);
            }
            close(
                4.0 * (dot(new, new) - dot(old, old)),
                receipt.fluid_drag_work,
                2e-16,
            );
            close(
                receipt.fluid_drag_work,
                -receipt.relative_dissipation,
                2e-16,
            );
            old = new;
        }
    }
}

const DIMS: [usize; 3] = [7, 5, 4];
const ORIGIN: [f64; 3] = [-0.4, 1.3, 0.2];
fn solid(p: [usize; 3]) -> bool {
    p == [0, 2, 0] || p == [3, 2, 2]
}
fn field(p: [usize; 3]) -> DarcyForchheimer3 {
    DarcyForchheimer3 {
        linear: LinearResistance3 {
            principal: [
                0.02 + 0.001 * p[0] as f64,
                0.08 + 0.003 * p[1] as f64,
                0.2 + 0.001 * p[2] as f64,
            ],
            axes: AXES,
            skeleton_velocity: [0.001, -0.002, 0.003],
            body_force: [
                1e-5 * (p[1] + 1) as f64,
                -2e-5 * (p[2] + 1) as f64,
                3e-5 * (p[0] + 1) as f64,
            ],
        },
        quadratic: 1.7,
    }
}
fn fixture<L: ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    periodic: [bool; 3],
    eq: impl Fn(f64, [f64; 3]) -> [f64; Q],
) -> ResistanceGrid3<L, Q> {
    let native = Grid3::try_new(
        DIMS,
        collision,
        [0.0; 3],
        periodic,
        GridBudget3 {
            maximum_cells: 140,
            maximum_mach: 0.2,
        },
        solid,
        |p| {
            eq(
                1.0 + 0.003 * (p[0] + p[1] + 2 * p[2]) as f64,
                [
                    0.01 + 0.001 * p[1] as f64,
                    -0.006 + 0.001 * p[2] as f64,
                    0.005 - 0.001 * p[0] as f64,
                ],
            )
        },
    )
    .unwrap();
    with_cx(|cx| {
        ResistanceGrid3::try_new(
            native,
            cx,
            (0..140).map(|i| field(cell(i, DIMS))).collect(),
            ORIGIN,
        )
        .unwrap()
    })
}
fn totals<L: ForcedLattice3<Q>, const Q: usize>(
    grid: &ResistanceGrid3<L, Q>,
) -> (f64, [f64; 3], [f64; 3]) {
    let mut mass = 0.0;
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (i, f) in grid.flow().populations().iter().enumerate() {
        let p = cell(i, DIMS);
        if solid(p) {
            continue;
        }
        mass += f.iter().sum::<f64>();
        let momentum: [f64; 3] = core::array::from_fn(|a| {
            f.iter()
                .zip(L::E)
                .map(|(f, e)| f * f64::from(e[a]))
                .sum::<f64>()
                + 0.5 * grid.flow().forces()[i][a]
        });
        let moment = cross(
            core::array::from_fn(|a| p[a] as f64 + 0.5 - ORIGIN[a]),
            momentum,
        );
        for a in 0..3 {
            linear[a] += momentum[a];
            angular[a] += moment[a];
        }
    }
    (mass, linear, angular)
}

fn transport_oracle<L: ForcedLattice3<Q>, const Q: usize>(
    grid: &ResistanceGrid3<L, Q>,
    periodic: [bool; 3],
    collisions: &[L::Collision],
) -> (Vec<[f64; Q]>, [f64; 3], usize) {
    let mut next = grid.flow().populations().to_vec();
    let mut angular = [0.0; 3];
    let mut crossings = 0;
    for (i, f) in grid.flow().populations().iter().enumerate() {
        let p = cell(i, DIMS);
        if solid(p) {
            continue;
        }
        let post = L::collide(*f, collisions[i], grid.flow().forces()[i]).unwrap();
        for q in 0..Q {
            let mut target = [0; 3];
            let mut wall = false;
            for a in 0..3 {
                let value = p[a] as i32 + L::E[q][a];
                wall |= !periodic[a] && (value < 0 || value >= DIMS[a] as i32);
                target[a] = value.rem_euclid(DIMS[a] as i32) as usize;
            }
            if wall || solid(target) {
                next[i][L::OPP[q]] = post[q];
            } else {
                let target_index = (target[2] * DIMS[1] + target[1]) * DIMS[0] + target[0];
                next[target_index][q] = post[q];
                let jump = core::array::from_fn(|a| {
                    target[a] as f64 - p[a] as f64 - f64::from(L::E[q][a])
                });
                let moment = cross(jump, L::E[q].map(|e| f64::from(e) * post[q]));
                for a in 0..3 {
                    angular[a] += moment[a];
                    if jump[a] != 0.0 {
                        crossings += 1;
                    }
                }
            }
        }
    }
    (next, angular, crossings)
}

fn spatial<L: ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    eq: impl Fn(f64, [f64; 3]) -> [f64; Q] + Copy,
) {
    for periodic in [[true, false, true], [true; 3]] {
        let mut serial = fixture::<L, Q>(collision, periodic, eq);
        let collisions = vec![collision; 140];
        let mut expected = Vec::new();
        let mut saw_periodic_torque = false;
        for _ in 0..3 {
            let before = totals(&serial);
            let (native, angular, crossings) = transport_oracle(&serial, periodic, &collisions);
            let receipt = with_cx(|cx| serial.step_serial(cx, Some(&collisions)).unwrap());
            assert_eq!(serial.flow().populations(), native);
            let after = totals(&serial);
            close(after.0, before.0, 1e-11);
            for a in 0..3 {
                close(
                    after.1[a] - before.1[a],
                    receipt.external_impulse[a] - receipt.skeleton_impulse[a]
                        + receipt.walls.fluid_population_impulse[a],
                    2e-12,
                );
                close(
                    after.2[a] - before.2[a],
                    receipt.external_angular_impulse[a] - receipt.skeleton_angular_impulse[a]
                        + receipt.walls.fluid_population_angular_impulse[a]
                        + receipt.periodic.angular_momentum_change[a],
                    5e-12,
                );
                close(
                    receipt.periodic.angular_momentum_change[a],
                    angular[a],
                    2e-12,
                );
                saw_periodic_torque |= angular[a].abs() > 1e-6;
            }
            assert_eq!(receipt.periodic.face_crossings, crossings);
            assert!(receipt.walls.measured_links > 0 && crossings > 0);
            close(
                receipt.fluid_drag_work + receipt.skeleton_work,
                -receipt.relative_dissipation,
                1e-14,
            );
            for i in 0..140 {
                if solid(cell(i, DIMS)) {
                    assert_eq!(serial.flow().forces()[i], [0.0; 3]);
                    continue;
                }
                let expected = force(
                    field(cell(i, DIMS)),
                    serial.flow().macroscopic(i).unwrap().1,
                );
                for a in 0..3 {
                    close(serial.flow().forces()[i][a], expected[a], 5e-15);
                }
            }
            expected.push((snapshot(&serial), receipt));
        }
        assert!(saw_periodic_torque);
        for workers in [1, 3] {
            let mut parallel = fixture::<L, Q>(collision, periodic, eq);
            let pool = TilePool::for_host(workers, SEED);
            for (state, receipt) in &expected {
                let actual = parallel
                    .step_pooled(&pool, &CancelGate::new_clock_free(), Some(&collisions))
                    .unwrap();
                assert_eq!(snapshot(&parallel), *state);
                assert_eq!(actual, *receipt);
            }
        }
    }
}

#[test]
fn original_transport_and_wall_periodic_momentum_are_complete_and_worker_invariant() {
    spatial::<D3q27, 27>(models()[1], |r, u| {
        equilibrium27(r, u, Equilibrium27::Factorized)
    });
    spatial::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
}

fn zero_drag<L: ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    eq: impl Fn(f64, [f64; 3]) -> [f64; Q] + Copy,
) {
    let make = || {
        Grid3::<L, Q>::try_new(
            DIMS,
            collision,
            [0.0; 3],
            [true, false, true],
            GridBudget3 {
                maximum_cells: 140,
                maximum_mach: 0.2,
            },
            solid,
            |_| eq(1.2, [0.008, -0.003, 0.002]),
        )
        .unwrap()
    };
    let laws: Vec<DarcyForchheimer3> = (0..140)
        .map(|i| {
            LinearResistance3 {
                body_force: field(cell(i, DIMS)).linear.body_force,
                ..LinearResistance3::default()
            }
            .into()
        })
        .collect();
    with_cx(|cx| {
        let forces = laws
            .iter()
            .enumerate()
            .map(|(i, law)| {
                if solid(cell(i, DIMS)) {
                    [0.0; 3]
                } else {
                    law.linear.body_force
                }
            })
            .collect();
        let mut reference = ForcedGrid3::try_new(make(), cx, forces).unwrap();
        let mut implicit = ResistanceGrid3::try_new(make(), cx, laws, ORIGIN).unwrap();
        let collisions = vec![collision; 140];
        for _ in 0..3 {
            reference.step_serial(cx, &collisions).unwrap();
            let receipt = implicit.step_serial(cx, Some(&collisions)).unwrap();
            assert_eq!(reference.populations(), implicit.flow().populations());
            assert_eq!(reference.forces(), implicit.flow().forces());
            assert_eq!(receipt.skeleton_impulse, [0.0; 3]);
            assert_eq!(receipt.relative_dissipation, 0.0);
        }
    });
}

#[test]
fn zero_resistance_is_exact_native_spatial_forcing_for_both_lattices() {
    for collision in models() {
        zero_drag::<D3q27, 27>(collision, |r, u| {
            equilibrium27(r, u, Equilibrium27::Factorized)
        });
    }
    zero_drag::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
}

#[test]
fn unresolved_stiff_force_refuses_and_resolved_contrast_retains_actual_source_history() {
    let law = LinearResistance3 {
        principal: [1e16; 3],
        ..LinearResistance3::default()
    };
    let grid = Grid3::<D3q27, 27>::try_new(
        [1; 3],
        models()[1],
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: 1,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| equilibrium27(1.0, [0.01, 0.0, 0.0], Equilibrium27::Factorized),
    )
    .unwrap();
    assert!(matches!(
        with_cx(|cx| ResistanceGrid3::try_new(grid, cx, vec![law.into()], [0.0; 3])),
        Err(GridError3::MaterialField { .. })
    ));
    linear_oracle::<D3q27, 27>(
        models()[5],
        equilibrium27(1.3, [0.011, -0.007, 0.005], Equilibrium27::Factorized),
    );
    let law = LinearResistance3 {
        principal: [1e-10, 0.03, 800.0],
        axes: AXES,
        ..LinearResistance3::default()
    };
    let mut grid = uniform::<D3q27, 27>(
        models()[1],
        equilibrium27(1.0, [0.012, -0.008, 0.004], Equilibrium27::Factorized),
        law.into(),
    );
    for _ in 0..4 {
        let (rho0, u0) = grid.flow().macroscopic(0).unwrap();
        let f0 = grid.flow().forces()[0];
        let receipt = with_cx(|cx| grid.step_serial(cx, None).unwrap());
        let (rho1, u1) = grid.flow().macroscopic(0).unwrap();
        for a in 0..3 {
            let applied = 4.0 * (f0[a] + grid.flow().forces()[0][a]);
            close(8.0 * (rho1 * u1[a] - rho0 * u0[a]), applied, 5e-14);
            close(
                receipt.external_impulse[a] - receipt.skeleton_impulse[a],
                applied,
                2e-16,
            );
        }
    }
}

#[test]
fn cancellation_after_either_native_pass_rolls_back_history_and_retries_exactly() {
    struct CancelPass {
        pool: TilePool,
        calls: AtomicUsize,
        after: usize,
    }
    impl KernelRunner for CancelPass {
        fn workers(&self) -> usize {
            self.pool.workers()
        }
        fn run_with_gate<K: TileKernel>(
            &self,
            kernel: &K,
            gate: &CancelGate,
        ) -> (Result<K::Out, RunError>, RunReport) {
            let result = self.pool.run_with_gate(kernel, gate);
            if self.calls.fetch_add(1, Ordering::AcqRel) == self.after {
                gate.request();
            }
            result
        }
    }
    let make = || {
        fixture::<D3q27, 27>(models()[1], [true, false, true], |r, u| {
            equilibrium27(r, u, Equilibrium27::Factorized)
        })
    };
    let mut control = make();
    let expected = with_cx(|cx| control.step_serial(cx, None).unwrap());
    for after in [0, 1] {
        let mut grid = make();
        let before = snapshot(&grid);
        let runner = CancelPass {
            pool: TilePool::for_host(3, SEED),
            calls: AtomicUsize::new(0),
            after,
        };
        assert!(matches!(
            grid.step_pooled(&runner, &CancelGate::new_clock_free(), None),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(snapshot(&grid), before);
        let receipt = grid
            .step_pooled(&runner.pool, &CancelGate::new_clock_free(), None)
            .unwrap();
        assert_eq!(snapshot(&grid), snapshot(&control));
        assert_eq!(receipt, expected);
    }
    let mut grid = make();
    let before = snapshot(&grid);
    let mut collisions = vec![models()[1]; 140];
    collisions[139].model = CollisionModel27::Bgk {
        tau: 0.4,
        equilibrium: Equilibrium27::Factorized,
    };
    assert!(with_cx(|cx| grid.step_serial(cx, Some(&collisions))).is_err());
    assert_eq!(snapshot(&grid), before);
    let receipt = with_cx(|cx| grid.step_serial(cx, None).unwrap());
    assert_eq!(snapshot(&grid), snapshot(&control));
    assert_eq!(receipt, expected);
}

#[test]
fn inadmissible_streamed_endpoint_never_publishes_partial_forces_or_step_count() {
    for pooled in [false, true] {
        let native = Grid3::<D3q27, 27>::try_new(
            [7, 5, 4],
            models()[1],
            [0.0; 3],
            [true; 3],
            GridBudget3 {
                maximum_cells: 140,
                maximum_mach: 0.2,
            },
            |_| false,
            |p| {
                equilibrium27(
                    if p == [5, 3, 2] { 1.0 } else { 1e-6 },
                    [0.0; 3],
                    Equilibrium27::Factorized,
                )
            },
        )
        .unwrap();
        let law = LinearResistance3 {
            principal: [0.01; 3],
            ..LinearResistance3::default()
        };
        let mut grid = with_cx(|cx| {
            ResistanceGrid3::try_new(native, cx, vec![law.into(); 140], [0.0; 3]).unwrap()
        });
        let before = snapshot(&grid);
        let result = if pooled {
            grid.step_pooled(
                &TilePool::for_host(3, SEED),
                &CancelGate::new_clock_free(),
                None,
            )
        } else {
            with_cx(|cx| grid.step_serial(cx, None))
        };
        assert!(matches!(result, Err(GridError3::Mach { .. })));
        assert_eq!(snapshot(&grid), before);
    }
}
