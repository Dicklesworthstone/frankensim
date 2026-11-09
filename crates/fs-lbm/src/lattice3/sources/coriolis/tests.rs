use super::*;
use crate::d3q19::CollisionModel3;
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
use crate::lattice3::{D3q19, GridBudget3};
use fs_exec::{RunError, RunReport, TileKernel, TilePool};
use std::sync::atomic::{AtomicUsize, Ordering};

const ROTATION: Coriolis3 = Coriolis3 {
    angular_velocity: [0.031, -0.023, 0.041],
};
const INITIAL: [f64; 3] = [0.018, -0.012, 0.009];
const SEED: u64 = 0xC081_0115;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: SEED,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn grid<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    dims: [usize; 3],
    periodic: [bool; 3],
    initial: impl Fn([usize; 3]) -> (f64, [f64; 3]),
) -> Grid3<L, Q> {
    Grid3::try_new_equilibrium(
        dims,
        collision,
        [0.0; 3],
        periodic,
        GridBudget3 {
            maximum_cells: dims.iter().product(),
            maximum_mach: 0.2,
        },
        |_| false,
        initial,
    )
    .unwrap()
}

fn bgk27() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.8,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|d| a[d] * b[d]).sum()
}

fn inertial_and_geostrophic<L: EquilibriumLattice3<Q>, const Q: usize>(law: L::Collision) {
    let omega = ROTATION.angular_velocity;
    let norm = dot(omega, omega).sqrt();
    let axis = omega.map(|value| value / norm);
    for geostrophic in [false, true] {
        let base = if geostrophic {
            cross(omega, INITIAL).map(|a| 2.0 * a)
        } else {
            [0.0; 3]
        };
        let mut flow = with_cx(|cx| {
            ForcedGrid3::try_new_with_coriolis_preserving_velocity(
                grid::<L, Q>(law, [4, 3, 2], [true; 3], |_| (1.2, INITIAL)),
                cx,
                ROTATION,
                vec![base; 24],
            )
            .unwrap()
        });
        let (rho0, u0) = flow.macroscopic(0).unwrap();
        for d in 0..3 {
            assert!((u0[d] - INITIAL[d]).abs() < 2e-15, "no initial half kick");
        }
        for step in 1..=24 {
            let (rho_before, before) = flow.macroscopic(0).unwrap();
            let f0 = flow.forces()[0];
            with_cx(|cx| {
                flow.advance_serial_with_coriolis(cx, &vec![law; 24], ROTATION, vec![base; 24])
                    .unwrap()
            });
            let (rho, velocity) = flow.macroscopic(0).unwrap();
            let force = flow.forces()[0];
            // Independent Rodrigues rotation at the analytic Cayley angle.
            let angle = -2.0 * norm.atan() * step as f64;
            let perpendicular = cross(axis, INITIAL);
            let expected: [f64; 3] = core::array::from_fn(|d| {
                if geostrophic {
                    INITIAL[d]
                } else {
                    INITIAL[d] * angle.cos()
                        + perpendicular[d] * angle.sin()
                        + axis[d] * dot(axis, INITIAL) * (1.0 - angle.cos())
                }
            });
            assert!((rho - rho0).abs() < 2e-13);
            for d in 0..3 {
                assert!((velocity[d] - expected[d]).abs() < 2e-12);
                assert!(
                    (rho * velocity[d] - rho_before * before[d] - 0.5 * (f0[d] + force[d])).abs()
                        < 2e-14
                );
            }
            let coriolis = core::array::from_fn(|d| force[d] - rho * base[d]);
            assert!(
                dot(coriolis, velocity).abs() < 2e-16,
                "Coriolis endpoint power"
            );
            if !geostrophic {
                assert!((dot(velocity, velocity) - dot(INITIAL, INITIAL)).abs() < 2e-13);
                assert!((dot(axis, velocity) - dot(axis, INITIAL)).abs() < 2e-13);
                let impulse = core::array::from_fn(|d| 0.5 * (f0[d] + force[d]));
                let midpoint = core::array::from_fn(|d| 0.5 * (before[d] + velocity[d]));
                assert!(
                    dot(impulse, midpoint).abs() < 2e-16,
                    "isolated source does no work"
                );
            }
        }
    }
}

#[test]
fn native_families_preserve_inertial_speed_and_geostrophic_balance() {
    inertial_and_geostrophic::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    inertial_and_geostrophic::<D3q19, 19>(CollisionModel3::ReducedCumulant {
        second_order_rate: 1.25,
        third_order_rate: 1.35,
        fourth_order_rate: 1.72,
    });
    inertial_and_geostrophic::<D3q27, 27>(bgk27());
    inertial_and_geostrophic::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Cumulant {
            rates: Rates27 { omega: [1.25; 10] },
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn transported_density_rotating_force_and_native_moving_wall_work_share_one_step() {
    let make = || {
        with_cx(|cx| {
            ForcedGrid3::try_new_with_coriolis_preserving_velocity(
                grid::<D3q27, 27>(bgk27(), [4, 4, 2], [true, false, true], |p| {
                    (1.0 + 0.03 * p[0] as f64, [0.01, 0.0, 0.0])
                }),
                cx,
                ROTATION,
                vec![[1e-5, 0.0, 0.0]; 32],
            )
            .unwrap()
        })
    };
    let mut actual = make();
    let mut reference = make();
    let probe = with_cx(|cx| actual.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap());
    let velocities: Vec<_> = probe
        .wall_links()
        .map(|link| {
            [
                if link.exterior_faces[3] {
                    0.015
                } else {
                    -0.005
                },
                0.0,
                0.0,
            ]
        })
        .collect();
    let momentum = |flow: &ForcedGrid3<D3q27, 27>| -> [f64; 3] {
        core::array::from_fn(|d| {
            (0..32)
                .map(|i| {
                    let (rho, u) = flow.macroscopic(i).unwrap();
                    rho * u[d]
                })
                .sum()
        })
    };
    let p0 = momentum(&actual);
    let old_forces = actual.forces().to_vec();
    let next: Vec<_> = (0..32)
        .map(|i| [2e-5 * (i as f64).cos(), -1e-5, 1e-5])
        .collect();
    let receipt = with_cx(|cx| {
        actual
            .advance_serial_with_coriolis_and_moving_wall_momentum(
                cx,
                &probe,
                &velocities,
                &vec![bgk27(); 32],
                ROTATION,
                next.clone(),
            )
            .unwrap()
    });
    let mut density_change = 0.0_f64;
    for (i, base) in next.iter().enumerate() {
        let (rho, u) = actual.macroscopic(i).unwrap();
        density_change = density_change.max((rho - reference.macroscopic(i).unwrap().0).abs());
        let skew = cross(ROTATION.angular_velocity, u);
        for d in 0..3 {
            assert!((actual.forces()[i][d] - rho * (base[d] - 2.0 * skew[d])).abs() < 3e-15);
        }
    }
    assert!(density_change > 1e-5);
    let original = with_cx(|cx| {
        reference
            .advance_serial_with_moving_wall_momentum(
                cx,
                &probe,
                &velocities,
                &vec![bgk27(); 32],
                actual.forces().to_vec(),
            )
            .unwrap()
    });
    assert_eq!(actual.populations(), reference.populations());
    assert_eq!(receipt, original);
    assert_ne!(receipt.wall_work, 0.0);
    for d in 0..3 {
        let impulse: f64 = old_forces
            .iter()
            .zip(actual.forces())
            .map(|(a, b)| 0.5 * (a[d] + b[d]))
            .sum();
        assert!(
            (momentum(&actual)[d] - p0[d] - impulse - receipt.fluid_population_impulse[d]).abs()
                < 3e-13
        );
    }
    for workers in [1, 3] {
        let mut pooled = make();
        let reported = pooled
            .advance_pooled_with_coriolis_and_moving_wall_momentum(
                &TilePool::for_host(workers, SEED),
                &CancelGate::new_clock_free(),
                &probe,
                &velocities,
                &vec![bgk27(); 32],
                ROTATION,
                next.clone(),
            )
            .unwrap();
        assert_eq!(reported, receipt);
        assert_eq!(pooled.populations(), actual.populations());
        assert_eq!(pooled.forces(), actual.forces());
    }
}

#[test]
fn zero_rotation_reproduces_acceleration_bits_with_ragged_executor_tiles() {
    let make = || {
        grid::<D3q27, 27>(bgk27(), [13, 5, 2], [true; 3], |p| {
            (1.0 + 0.01 * (p[0] as f64).cos(), INITIAL)
        })
    };
    let start = vec![[1e-5, -2e-5, 3e-5]; 130];
    let mut original = with_cx(|cx| {
        ForcedGrid3::try_new_with_body_acceleration_preserving_velocity(make(), cx, start.clone())
            .unwrap()
    });
    let mut rotating = with_cx(|cx| {
        ForcedGrid3::try_new_with_coriolis_preserving_velocity(
            make(),
            cx,
            Coriolis3::default(),
            start,
        )
        .unwrap()
    });
    assert_eq!(original.populations(), rotating.populations());
    assert_eq!(original.forces(), rotating.forces());
    for epoch in 0..3 {
        let next: Vec<_> = (0..130)
            .map(|i| [1e-5 * ((i + epoch) as f64).sin(), -2e-5, 3e-5])
            .collect();
        with_cx(|cx| {
            original
                .advance_serial_with_body_acceleration(cx, &vec![bgk27(); 130], next.clone())
                .unwrap()
        });
        rotating
            .advance_pooled_with_coriolis(
                &TilePool::for_host(3, SEED),
                &CancelGate::new_clock_free(),
                &vec![bgk27(); 130],
                Coriolis3::default(),
                next,
            )
            .unwrap();
        assert_eq!(original.populations(), rotating.populations());
        assert_eq!(original.forces(), rotating.forces());
    }
}

#[test]
fn unresolved_or_cancelled_rotation_keeps_source_history_and_retries_exactly() {
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
        with_cx(|cx| {
            ForcedGrid3::try_new_with_coriolis_preserving_velocity(
                grid::<D3q27, 27>(bgk27(), [13, 5, 2], [true; 3], |_| (1.1, INITIAL)),
                cx,
                ROTATION,
                vec![[0.0; 3]; 130],
            )
            .unwrap()
        })
    };
    let laws = vec![bgk27(); 130];
    let next = vec![[1e-5, -2e-5, 3e-5]; 130];
    let mut expected = make();
    with_cx(|cx| {
        expected
            .advance_serial_with_coriolis(cx, &laws, ROTATION, next.clone())
            .unwrap()
    });
    for after in [0, 1] {
        let mut flow = make();
        let populations = flow.populations().to_vec();
        let forces = flow.forces().to_vec();
        with_cx(|cx| {
            assert!(
                flow.advance_serial_with_coriolis(
                    cx,
                    &laws,
                    Coriolis3 {
                        angular_velocity: [f64::NAN; 3]
                    },
                    next.clone()
                )
                .is_err()
            );
            let mut excessive = next.clone();
            excessive[129] = [1.0, 0.0, 0.0];
            assert!(matches!(
                flow.advance_serial_with_coriolis(cx, &laws, ROTATION, excessive),
                Err(GridError3::Mach { cell: 129, .. })
            ));
        });
        let runner = CancelPass {
            pool: TilePool::for_host(3, SEED),
            calls: AtomicUsize::new(0),
            after,
        };
        assert!(matches!(
            flow.advance_pooled_with_coriolis(
                &runner,
                &CancelGate::new_clock_free(),
                &laws,
                ROTATION,
                next.clone(),
            ),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.forces(), forces);
        assert_eq!(flow.completed_steps(), 0);
        flow.advance_pooled_with_coriolis(
            &runner.pool,
            &CancelGate::new_clock_free(),
            &laws,
            ROTATION,
            next.clone(),
        )
        .unwrap();
        assert_eq!(flow.populations(), expected.populations());
        assert_eq!(flow.forces(), expected.forces());
        assert_eq!(flow.completed_steps(), 1);
    }
}
