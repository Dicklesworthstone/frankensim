use super::*;
use crate::d3q19::CollisionModel3;
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27};
use crate::lattice3::{D3q19, GridBudget3};
use fs_exec::{RunError, RunReport, TileKernel, TilePool};
use std::sync::atomic::{AtomicUsize, Ordering};

const SEED: u64 = 0xF8A6_E003;
const INITIAL: [f64; 3] = [0.018, -0.012, 0.009];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey { seed: SEED, kernel_id: 3, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn grid<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    dims: [usize; 3],
    periodic: [bool; 3],
    solid: impl Fn([usize; 3]) -> bool,
    initial: impl Fn([usize; 3]) -> (f64, [f64; 3]),
) -> Grid3<L, Q> {
    Grid3::try_new_equilibrium(
        dims,
        collision,
        [0.0; 3],
        periodic,
        GridBudget3 { maximum_cells: dims.iter().product(), maximum_mach: 0.2 },
        solid,
        initial,
    )
    .unwrap()
}

fn bgk27() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Bgk { tau: 0.8, equilibrium: Equilibrium27::Factorized },
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

// Independent component arithmetic for the prescribed physical law.
fn base(frame: RigidFrame3, cell: usize, dims: [usize; 3], gravity: [f64; 3]) -> [f64; 3] {
    let p = [cell % dims[0], (cell / dims[0]) % dims[1], cell / (dims[0] * dims[1])];
    let r = core::array::from_fn(|d| p[d] as f64 + 0.5 - frame.origin[d]);
    let euler = cross(frame.angular_acceleration, r);
    let centrifugal = cross(frame.angular_velocity, cross(frame.angular_velocity, r));
    core::array::from_fn(|d| gravity[d] - frame.origin_acceleration[d] - euler[d] - centrifugal[d])
}

fn frame(tick: usize) -> RigidFrame3 {
    let angular_acceleration = [2e-4, -1.5e-4, 1e-4];
    RigidFrame3 {
        origin: [2.4, 1.2, 2.7],
        angular_velocity: core::array::from_fn(|d| {
            [0.012, -0.009, 0.015][d] + tick as f64 * angular_acceleration[d]
        }),
        angular_acceleration,
        origin_acceleration: [2e-5 + tick as f64 * 1e-5, -4e-5, 3e-5],
    }
}

fn all_terms<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let dims = [5, 3, 5];
    let count = 75;
    let solid_cell = 37; // [2,1,2], excluded from the physical source field.
    let make = || {
        grid::<L, Q>(collision, dims, [false; 3], |p| p == [2, 1, 2], |p| {
            (1.0 + 0.004 * p[0] as f64 + 0.002 * p[1] as f64, INITIAL)
        })
    };
    let body = |tick: usize| -> Vec<[f64; 3]> {
        (0..count)
            .map(|i| {
                if i == solid_cell {
                    [0.0; 3]
                } else {
                    [1e-5 * (i as f64).cos(), -2e-5, 3e-5 + tick as f64 * 1e-5]
                }
            })
            .collect()
    };
    let make_flow = || {
        with_cx(|cx| {
            ForcedGrid3::try_new_with_rigid_frame_preserving_velocity(make(), cx, frame(0), body(0))
                .unwrap()
        })
    };
    let mut actual = make_flow();
    let mut pooled = make_flow();
    let manual_base = |tick: usize| {
        body(tick).into_iter().enumerate()
            .map(|(i, g)| if i == solid_cell { [0.0; 3] } else { base(frame(tick), i, dims, g) })
            .collect()
    };
    let mut original = with_cx(|cx| {
        ForcedGrid3::try_new_with_coriolis_preserving_velocity(
            make(), cx, Coriolis3 { angular_velocity: frame(0).angular_velocity }, manual_base(0),
        )
        .unwrap()
    });
    for (i, gravity) in body(0).iter().enumerate() {
        if i == solid_cell {
            assert_eq!(actual.forces()[i], [0.0; 3]);
            continue;
        }
        let (rho, velocity) = actual.macroscopic(i).unwrap();
        let skew = cross(frame(0).angular_velocity, velocity);
        let non_coriolis = base(frame(0), i, dims, *gravity);
        for d in 0..3 {
            assert!((velocity[d] - INITIAL[d]).abs() < 3e-15, "no startup half kick");
            assert!((actual.forces()[i][d] - rho * (non_coriolis[d] - 2.0 * skew[d])).abs() < 5e-15);
        }
    }
    let initial_density: Vec<_> = (0..count)
        .map(|i| if i == solid_cell { 0.0 } else { actual.macroscopic(i).unwrap().0 })
        .collect();
    let collisions = vec![collision; count];
    let pool = TilePool::for_host(3, SEED);
    for tick in 1..=3 {
        let end = frame(tick);
        with_cx(|cx| {
            actual.advance_serial_with_rigid_frame(cx, &collisions, end, body(tick)).unwrap();
            original.advance_serial_with_coriolis(
                cx, &collisions, Coriolis3 { angular_velocity: end.angular_velocity }, manual_base(tick),
            )
            .unwrap();
        });
        pooled.advance_pooled_with_rigid_frame(
            &pool, &CancelGate::new_clock_free(), &collisions, end, body(tick),
        )
        .unwrap();
        assert_eq!(actual.populations(), pooled.populations());
        assert_eq!(actual.forces(), pooled.forces());
        assert_eq!(actual.completed_steps(), tick);
        let mut density_change = 0.0_f64;
        for (i, gravity) in body(tick).iter().enumerate() {
            if i == solid_cell {
                assert_eq!(actual.forces()[i], [0.0; 3]);
                continue;
            }
            let (rho, velocity) = actual.macroscopic(i).unwrap();
            let (reference_rho, reference_velocity) = original.macroscopic(i).unwrap();
            let non_coriolis = base(end, i, dims, *gravity);
            let skew = cross(end.angular_velocity, velocity);
            density_change = density_change.max((rho - initial_density[i]).abs());
            assert!((rho - reference_rho).abs() < 5e-13);
            for d in 0..3 {
                assert!((velocity[d] - reference_velocity[d]).abs() < 5e-13);
                assert!((actual.forces()[i][d] - rho * (non_coriolis[d] - 2.0 * skew[d])).abs() < 5e-15);
            }
            let coriolis = core::array::from_fn(|d| actual.forces()[i][d] - rho * non_coriolis[d]);
            assert!(dot(coriolis, velocity).abs() < 3e-16);
        }
        assert!(density_change > 1e-5, "endpoint must use transported density");
    }
}

#[test]
fn native_d3q19_and_d3q27_evolve_changing_full_frames_in_actual_physical_velocity() {
    all_terms::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    all_terms::<D3q19, 19>(CollisionModel3::ReducedCumulant {
        second_order_rate: 1.25, third_order_rate: 1.35, fourth_order_rate: 1.72,
    });
    all_terms::<D3q27, 27>(bgk27());
    all_terms::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Cumulant { rates: Rates27 { omega: [1.25; 10] } },
        maximum_mach: 0.2,
    });
}

fn axial_frame(tick: usize) -> RigidFrame3 {
    RigidFrame3 {
        origin: [4.2, 2.3, 0.5],
        angular_velocity: [0.0, 0.0, 0.018 + tick as f64 * 2e-4],
        angular_acceleration: [0.0, 0.0, 2e-4],
        origin_acceleration: [1e-4 + tick as f64 * 1e-5, -2e-4, 3e-5],
    }
}

#[test]
fn full_frame_stream_keeps_relative_wall_work_and_trapezoidal_body_impulse_separate() {
    let dims = [9, 5, 2];
    let count = 90;
    let make = || {
        with_cx(|cx| {
            ForcedGrid3::try_new_with_rigid_frame_preserving_velocity(
                grid::<D3q27, 27>(
                    bgk27(), dims, [false, false, true], |_| false,
                    |p| (1.0 + 0.004 * p[0] as f64 - 0.003 * p[1] as f64, INITIAL),
                ),
                cx, axial_frame(0), vec![[5e-5, 2e-5, -1e-5]; count],
            )
            .unwrap()
        })
    };
    let mut actual = make();
    let mut original = make();
    let probe = with_cx(|cx| actual.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap());
    let velocities: Vec<_> = probe.wall_links()
        .map(|link| [0.0, 0.0, if link.exterior_faces[1] { 0.009 } else { -0.003 }])
        .collect();
    let momentum = |flow: &ForcedGrid3<D3q27, 27>| -> [f64; 3] {
        core::array::from_fn(|d| {
            (0..count).map(|i| {
                let (rho, velocity) = flow.macroscopic(i).unwrap();
                rho * velocity[d]
            }).sum()
        })
    };
    let before = momentum(&actual);
    let old_forces = actual.forces().to_vec();
    let end = axial_frame(1);
    let gravity = [3e-5, -4e-5, 5e-6];
    let collisions = vec![bgk27(); count];
    let receipt = with_cx(|cx| {
        actual.advance_serial_with_rigid_frame_and_moving_wall_momentum(
            cx, &probe, &velocities, &collisions, end, vec![gravity; count],
        )
        .unwrap()
    });
    for i in 0..count {
        let (rho, u) = actual.macroscopic(i).unwrap();
        let rx = (i % dims[0]) as f64 + 0.5 - end.origin[0];
        let ry = ((i / dims[0]) % dims[1]) as f64 + 0.5 - end.origin[1];
        let w = end.angular_velocity[2];
        let alpha = end.angular_acceleration[2];
        // Explicit signs for right-handed axial spin: Euler, centrifugal,
        // translation and Coriolis all contribute to the exported force.
        let expected = [
            gravity[0] - end.origin_acceleration[0] + alpha * ry + w * w * rx + 2.0 * w * u[1],
            gravity[1] - end.origin_acceleration[1] - alpha * rx + w * w * ry - 2.0 * w * u[0],
            gravity[2] - end.origin_acceleration[2],
        ];
        for d in 0..3 {
            assert!((actual.forces()[i][d] - rho * expected[d]).abs() < 5e-15);
        }
    }
    let reference_receipt = with_cx(|cx| {
        original.advance_serial_with_moving_wall_momentum(
            cx, &probe, &velocities, &collisions, actual.forces().to_vec(),
        )
        .unwrap()
    });
    assert_eq!(actual.populations(), original.populations());
    assert_eq!(receipt, reference_receipt);
    assert_ne!(receipt.wall_work, 0.0);
    for d in 0..3 {
        let impulse: f64 = old_forces.iter().zip(actual.forces())
            .map(|(a, b)| 0.5 * (a[d] + b[d])).sum();
        assert!((momentum(&actual)[d] - before[d] - impulse - receipt.fluid_population_impulse[d]).abs() < 8e-13);
    }
    for workers in [1, 3] {
        let mut pooled = make();
        let reported = pooled.advance_pooled_with_rigid_frame_and_moving_wall_momentum(
            &TilePool::for_host(workers, SEED),
            &CancelGate::new_clock_free(),
            &probe, &velocities, &collisions, end, vec![gravity; count],
        )
        .unwrap();
        assert_eq!(pooled.populations(), actual.populations());
        assert_eq!(pooled.forces(), actual.forces());
        assert_eq!(reported, receipt);
    }
}

fn zero_and_free_fall<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let make = || {
        grid::<L, Q>(collision, [13, 5, 2], [true; 3], |_| false, |p| {
            (1.0 + 0.01 * (p[0] as f64).cos(), INITIAL)
        })
    };
    let neutral = RigidFrame3 { origin: [-300.0, 70.0, -120.0], ..RigidFrame3::default() };
    let start = vec![[1e-5, -2e-5, 3e-5]; 130];
    let mut original = with_cx(|cx| {
        ForcedGrid3::try_new_with_body_acceleration_preserving_velocity(make(), cx, start.clone()).unwrap()
    });
    let mut actual = with_cx(|cx| {
        ForcedGrid3::try_new_with_rigid_frame_preserving_velocity(make(), cx, neutral, start).unwrap()
    });
    assert_eq!(actual.populations(), original.populations());
    assert_eq!(actual.forces(), original.forces());
    let pool = TilePool::for_host(3, SEED);
    for tick in 0..3 {
        let next: Vec<_> = (0..130)
            .map(|i| [1e-5 * ((i + tick) as f64).sin(), -2e-5, 3e-5]).collect();
        with_cx(|cx| {
            original.advance_serial_with_body_acceleration(cx, &vec![collision; 130], next.clone()).unwrap()
        });
        actual.advance_pooled_with_rigid_frame(
            &pool, &CancelGate::new_clock_free(), &vec![collision; 130], neutral, next,
        )
        .unwrap();
        assert_eq!(actual.populations(), original.populations());
        assert_eq!(actual.forces(), original.forces());
    }

    // A freely falling, nonrotating periodic frame feels no gravity.
    let uniform = || grid::<L, Q>(
        collision, [4, 3, 2], [true; 3], |_| false, |_| (1.1, INITIAL),
    );
    let translation = |tick: usize| RigidFrame3 {
        origin_acceleration: [1e-4 + tick as f64 * 2e-5, -3e-4, 2e-4],
        ..RigidFrame3::default()
    };
    let mut free_fall = with_cx(|cx| {
        ForcedGrid3::try_new_with_rigid_frame_preserving_velocity(
            uniform(), cx, translation(0), vec![translation(0).origin_acceleration; 24],
        )
        .unwrap()
    });
    let mut unforced = with_cx(|cx| {
        ForcedGrid3::try_new_preserving_velocity(uniform(), cx, vec![[0.0; 3]; 24]).unwrap()
    });
    for tick in 1..=3 {
        let end = translation(tick);
        with_cx(|cx| {
            free_fall.advance_serial_with_rigid_frame(
                cx, &vec![collision; 24], end, vec![end.origin_acceleration; 24],
            )
            .unwrap();
            unforced.step_serial(cx, &vec![collision; 24]).unwrap();
        });
        assert_eq!(free_fall.populations(), unforced.populations());
        assert!(free_fall.forces().iter().all(|force| *force == [0.0; 3]));
        let (_, velocity) = free_fall.macroscopic(0).unwrap();
        for d in 0..3 {
            assert!((velocity[d] - INITIAL[d]).abs() < 5e-14);
        }
    }
}

#[test]
fn zero_frame_reduces_exactly_and_translational_free_fall_is_periodic_in_both_lattices() {
    zero_and_free_fall::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    zero_and_free_fall::<D3q27, 27>(bgk27());
}

#[test]
fn incompatible_periodicity_late_mach_and_cancel_preserve_frame_history_and_wall_replay() {
    struct CancelPass { pool: TilePool, calls: AtomicUsize, after: usize }
    impl KernelRunner for CancelPass {
        fn workers(&self) -> usize { self.pool.workers() }
        fn run_with_gate<K: TileKernel>(
            &self, kernel: &K, gate: &CancelGate,
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
            ForcedGrid3::try_new_with_rigid_frame_preserving_velocity(
                grid::<D3q27, 27>(
                    bgk27(), [13, 5, 2], [false, false, true], |_| false,
                    |p| (1.0 + 0.001 * p[0] as f64, INITIAL),
                ),
                cx, axial_frame(0), vec![[1e-5, -2e-5, 3e-5]; 130],
            )
            .unwrap()
        })
    };
    let mut expected = make();
    let probe = with_cx(|cx| expected.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap());
    let velocities: Vec<_> = probe.wall_links().map(|_| [0.0, 0.0, 0.004]).collect();
    let collisions = vec![bgk27(); 130];
    let next = vec![[2e-5, -3e-5, 1e-5]; 130];
    let end = axial_frame(1);
    let receipt = with_cx(|cx| {
        expected.advance_serial_with_rigid_frame_and_moving_wall_momentum(
            cx, &probe, &velocities, &collisions, end, next.clone(),
        )
        .unwrap()
    });
    for after in [0, 1] {
        let mut actual = make();
        let populations = actual.populations().to_vec();
        let forces = actual.forces().to_vec();
        with_cx(|cx| {
            for invalid in [
                RigidFrame3 { angular_velocity: [1e-3, 0.0, end.angular_velocity[2]], ..end },
                RigidFrame3 { angular_acceleration: [1e-3, 0.0, 2e-4], ..end },
                RigidFrame3 { origin: [f64::NAN; 3], ..end },
            ] {
                assert!(matches!(
                    actual.advance_serial_with_rigid_frame(cx, &collisions, invalid, next.clone()),
                    Err(GridError3::Input(_))
                ));
            }
            let mut excessive = next.clone();
            excessive[129] = [1.0, 0.0, 0.0];
            assert!(matches!(
                actual.advance_serial_with_rigid_frame_and_moving_wall_momentum(
                    cx, &probe, &velocities, &collisions, end, excessive,
                ),
                Err(GridError3::Mach { cell: 129, .. })
            ));
        });
        let runner = CancelPass { pool: TilePool::for_host(3, SEED), calls: AtomicUsize::new(0), after };
        assert!(matches!(
            actual.advance_pooled_with_rigid_frame_and_moving_wall_momentum(
                &runner, &CancelGate::new_clock_free(), &probe, &velocities, &collisions, end, next.clone(),
            ),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(actual.populations(), populations);
        assert_eq!(actual.forces(), forces);
        assert_eq!(actual.completed_steps(), 0);
        let replay = actual.advance_pooled_with_rigid_frame_and_moving_wall_momentum(
            &runner.pool, &CancelGate::new_clock_free(), &probe, &velocities, &collisions, end, next.clone(),
        )
        .unwrap();
        assert_eq!(actual.populations(), expected.populations());
        assert_eq!(actual.forces(), expected.forces());
        assert_eq!(actual.completed_steps(), 1);
        assert_eq!(replay, receipt);
    }
}
