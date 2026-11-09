use super::*;
use fs_exec::{
    Budget, ExecMode, KernelRunner, RunError, RunReport, StreamKey, TileKernel, TilePool,
};
use std::sync::atomic::{AtomicUsize, Ordering};

const DIMS: [usize; 3] = [7, 6, 5];
const CELLS: usize = 210;
const ORIGIN: [f64; 3] = [0.25, -0.7, 1.1];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn field(epoch: usize) -> Vec<[f64; 3]> {
    (0..CELLS)
        .map(|cell| variable_force(cell + epoch, epoch).map(|v| 4.0 * v))
        .collect()
}

fn make<L: EquilibriumLattice3<Q>, const Q: usize>(
    law: L::Collision,
    axis: usize,
    initial_acceleration: Vec<[f64; 3]>,
) -> ForcedOpenBoundaryGrid3<L, Q> {
    let mut wall_velocity = [0.0; 3];
    wall_velocity[axis] = 0.006;
    let mut faces = [FaceBoundary3::Wall {
        velocity: wall_velocity,
    }; 6];
    faces[2 * axis] = FaceBoundary3::Velocity { velocity: VELOCITY };
    faces[2 * axis + 1] = FaceBoundary3::Pressure { density: RHO };
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        DIMS,
        law,
        [0.0; 3],
        [false; 3],
        GridBudget3 {
            maximum_cells: CELLS,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let cell = p[0] + DIMS[0] * (p[1] + DIMS[1] * p[2]);
            (RHO * (1.0 + 0.03 * (0.2 * cell as f64).cos()), VELOCITY)
        },
    )
    .unwrap();
    let mut inner =
        OpenBoundaryGrid3::try_new_with_poll(grid, BoundarySpec3::new(faces), &mut || Ok(()))
            .unwrap();
    inner
        .set_pressure_profile(
            Face3::ALL[2 * axis + 1],
            &CancelGate::new_clock_free(),
            |p| {
                Ok(RHO
                    * (1.0 + 0.004 * p[(axis + 1) % 3] as f64 - 0.003 * p[(axis + 2) % 3] as f64))
            },
        )
        .unwrap();
    inner
        .set_velocity_profile(
            Face3::ALL[2 * axis],
            InletFlux3::Unconstrained,
            &CancelGate::new_clock_free(),
            |p| {
                Ok(core::array::from_fn(|a| {
                    VELOCITY[a] + 0.0003 * (p[(a + 1) % 3] as f64 - 2.0)
                }))
            },
        )
        .unwrap();
    with_cx(|cx| {
        ForcedOpenBoundaryGrid3::try_new_with_body_acceleration_preserving_velocity(
            inner,
            cx,
            initial_acceleration,
        )
        .unwrap()
    })
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn arm(cell: usize) -> [f64; 3] {
    let p = crate::lattice3::coordinates(cell, DIMS);
    core::array::from_fn(|a| p[a] as f64 + 0.5 - ORIGIN[a])
}

fn moments<L: EquilibriumLattice3<Q>, const Q: usize>(
    flow: &ForcedOpenBoundaryGrid3<L, Q>,
) -> (f64, [f64; 3], [f64; 3]) {
    let (mut mass, mut momentum, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for (cell, populations) in flow.populations().iter().enumerate() {
        mass += populations.iter().sum::<f64>();
        let p = core::array::from_fn(|a| {
            populations
                .iter()
                .zip(L::E)
                .map(|(f, e)| f * f64::from(e[a]))
                .sum::<f64>()
                + 0.5 * flow.forces()[cell][a]
        });
        let torque = cross(arm(cell), p);
        for a in 0..3 {
            momentum[a] += p[a];
            angular[a] += torque[a];
        }
    }
    (mass, momentum, angular)
}

fn exercise<L: EquilibriumLattice3<Q>, const Q: usize>(law: L::Collision) {
    for axis in 0..3 {
        let mut actual = make::<L, Q>(law, axis, field(0));
        let mut reference = make::<L, Q>(law, axis, field(0));
        let before = moments(&actual);
        let start = actual.forces().to_vec();
        let initial_acceleration = field(0);
        for cell in 0..CELLS {
            let (rho, velocity) = actual.macroscopic(cell).unwrap();
            for a in 0..3 {
                near(start[cell][a], rho * initial_acceleration[cell][a], 2e-17);
                near(velocity[a], VELOCITY[a], 3e-14);
            }
        }
        let acceleration = field(1);
        let laws = vec![law; CELLS];
        let receipt = with_cx(|cx| {
            actual
                .advance_serial_with_body_acceleration(cx, Some(&laws), &acceleration, ORIGIN)
                .unwrap()
        });

        // Independent original force-density steps converge the nonlinear
        // endpoint consistency, without using the new direct seam closure.
        let mut force: Vec<_> = reference
            .populations()
            .iter()
            .zip(&acceleration)
            .map(|(f, a)| a.map(|v| v * f.iter().sum::<f64>()))
            .collect();
        let mut frozen = Vec::new();
        for iteration in 0..20 {
            reference
                .stage_serial_with_poll(&force, Some(&laws), ORIGIN, &mut || Ok(()))
                .unwrap();
            if iteration == 0 {
                frozen = reference.inner.grid.next.clone();
            }
            force = reference
                .inner
                .grid
                .next
                .iter()
                .zip(&acceleration)
                .map(|(f, a)| a.map(|v| v * f.iter().sum::<f64>()))
                .collect();
        }
        let (expected, step) = reference
            .stage_serial_with_poll(&force, Some(&laws), ORIGIN, &mut || Ok(()))
            .unwrap();
        reference.publish_with_force_swap(step, &mut force);
        let mut frozen_error = 0.0_f64;
        for cell in 0..CELLS {
            let rho = actual.populations()[cell].iter().sum::<f64>();
            assert_eq!(actual.forces()[cell], acceleration[cell].map(|a| rho * a));
            for q in 0..Q {
                near(
                    actual.populations()[cell][q],
                    reference.populations()[cell][q],
                    7e-13,
                );
                frozen_error =
                    frozen_error.max((actual.populations()[cell][q] - frozen[cell][q]).abs());
            }
        }
        assert!(
            frozen_error > 1e-8,
            "freezing endpoint density must fail this fixture"
        );
        let wall = actual.last_wall_exchange().unwrap();
        let expected_wall = reference.last_wall_exchange().unwrap();
        near(wall.wall_work, expected_wall.wall_work, 2e-12);
        near(wall.wall_work, 0.006 * wall.wall_impulse[axis], 2e-13);
        let after = moments(&actual);
        near(
            after.0 - before.0,
            receipt.mass_change + wall.fluid_mass_change,
            2e-11,
        );
        let impulse = body_impulse(&start, actual.forces());
        let mut body_angular = [0.0; 3];
        for cell in 0..CELLS {
            let local = core::array::from_fn(|a| 0.5 * (start[cell][a] + actual.forces()[cell][a]));
            let torque = cross(arm(cell), local);
            for a in 0..3 {
                body_angular[a] += torque[a];
            }
        }
        for a in 0..3 {
            near(receipt.impulse[a], expected.impulse[a], 2e-11);
            near(
                after.1[a] - before.1[a],
                impulse[a] + receipt.impulse[a] + wall.fluid_population_impulse[a],
                2e-11,
            );
            near(
                after.2[a] - before.2[a],
                body_angular[a]
                    + receipt.angular_impulse[a]
                    + wall.fluid_population_angular_impulse[a],
                6e-11,
            );
        }

        // Two source epochs verify that published rho*a drives the next native
        // collision, with bit-identical tiled populations and wall receipts.
        with_cx(|cx| {
            actual
                .advance_serial_with_body_acceleration(cx, Some(&laws), &field(2), ORIGIN)
                .unwrap()
        });
        for workers in [1, 3] {
            let mut pooled = make::<L, Q>(law, axis, field(0));
            let pool = TilePool::for_host(workers, 71);
            for epoch in 1..=2 {
                pooled
                    .advance_pooled_with_body_acceleration(
                        &pool,
                        &CancelGate::new_clock_free(),
                        Some(&laws),
                        &field(epoch),
                        ORIGIN,
                    )
                    .unwrap();
            }
            assert_eq!(pooled.populations(), actual.populations());
            assert_eq!(pooled.forces(), actual.forces());
            assert_eq!(pooled.last_wall_exchange(), actual.last_wall_exchange());
        }
    }
}

#[test]
fn open_acceleration_closes_native_density_momentum_torque_and_wall_work() {
    exercise::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    exercise::<D3q19, 19>(CollisionModel3::ReducedCumulant {
        second_order_rate: 1.2,
        third_order_rate: 1.1,
        fourth_order_rate: 0.7,
    });
    for model in [
        CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.2; 10] },
        },
    ] {
        exercise::<D3q27, 27>(CollisionConfig27 {
            model,
            maximum_mach: 0.2,
        });
    }
}

#[test]
fn zero_acceleration_and_pooled_force_density_retain_original_bits() {
    let law = CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Quadratic,
        },
        maximum_mach: 0.2,
    };
    let mut actual = make::<D3q27, 27>(law, 0, vec![[0.0; 3]; CELLS]);
    let mut reference = make::<D3q27, 27>(law, 0, vec![[0.0; 3]; CELLS]);
    let mut pooled = make::<D3q27, 27>(law, 0, vec![[0.0; 3]; CELLS]);
    let a = with_cx(|cx| {
        actual
            .advance_serial_with_body_acceleration(cx, None, &vec![[0.0; 3]; CELLS], ORIGIN)
            .unwrap()
    });
    let b = with_cx(|cx| reference.step_serial(cx, None, ORIGIN).unwrap());
    let c = pooled
        .advance_pooled(
            &TilePool::for_host(3, 71),
            &CancelGate::new_clock_free(),
            None,
            vec![[0.0; 3]; CELLS],
            ORIGIN,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(b, c);
    assert_eq!(actual.populations(), reference.populations());
    assert_eq!(pooled.populations(), reference.populations());
    let c = pooled
        .step_pooled(
            &TilePool::for_host(3, 71),
            &CancelGate::new_clock_free(),
            None,
            ORIGIN,
        )
        .unwrap();
    let b = with_cx(|cx| reference.step_serial(cx, None, ORIGIN).unwrap());
    assert_eq!(b, c);
    assert_eq!(pooled.populations(), reference.populations());
}

#[test]
fn open_acceleration_refusal_and_late_cancellation_preserve_all_accepted_state() {
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
    let law = CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    };
    let make = || {
        let mut flow = make::<D3q27, 27>(law, 0, field(0));
        with_cx(|cx| {
            flow.advance_serial_with_body_acceleration(cx, None, &field(1), ORIGIN)
                .unwrap()
        });
        flow
    };
    let mut expected = make();
    with_cx(|cx| {
        expected
            .advance_serial_with_body_acceleration(cx, None, &field(2), ORIGIN)
            .unwrap()
    });
    for after in [0, 1] {
        let mut actual = make();
        let before = (
            actual.populations().to_vec(),
            actual.forces().to_vec(),
            actual.last_wall_exchange(),
        );
        let mut bad = field(2);
        bad[150][0] = 1.0;
        with_cx(|cx| {
            assert!(matches!(
                actual.advance_serial_with_body_acceleration(cx, None, &bad, ORIGIN),
                Err(GridError3::Mach { .. })
            ));
            bad[150][0] = f64::NAN;
            assert!(
                actual
                    .advance_serial_with_body_acceleration(cx, None, &bad, ORIGIN)
                    .is_err()
            );
            assert!(
                actual
                    .advance_serial_with_body_acceleration(cx, None, &field(2)[..CELLS - 1], ORIGIN)
                    .is_err()
            );
        });
        let runner = CancelPass {
            pool: TilePool::for_host(3, 71),
            calls: AtomicUsize::new(0),
            after,
        };
        assert!(matches!(
            actual.advance_pooled_with_body_acceleration(
                &runner,
                &CancelGate::new_clock_free(),
                None,
                &field(2),
                ORIGIN
            ),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(actual.populations(), before.0);
        assert_eq!(actual.forces(), before.1);
        assert_eq!(actual.last_wall_exchange(), before.2);
        assert_eq!(actual.completed_steps(), 1);
        actual
            .advance_pooled_with_body_acceleration(
                &runner.pool,
                &CancelGate::new_clock_free(),
                None,
                &field(2),
                ORIGIN,
            )
            .unwrap();
        assert_eq!(actual.populations(), expected.populations());
        assert_eq!(actual.forces(), expected.forces());
        assert_eq!(actual.last_wall_exchange(), expected.last_wall_exchange());
    }
}
