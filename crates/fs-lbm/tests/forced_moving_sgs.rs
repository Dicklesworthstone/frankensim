//! Force-epoch, moving-wall momentum and SGS composition regressions.

use fs_exec::{CancelGate, Cx, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{CollisionModel3, equilibrium3};
use fs_lbm::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use fs_lbm::lattice3::sources::{ForcedGrid3, ForcedLattice3};
use fs_lbm::lattice3::wall_loads::WallMomentumProbe3;
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, GridError3, Lattice3};
use fs_lbm::sgs3::dynamic::channel::ChannelDynamicSgsConfig;
use fs_lbm::sgs3::{SgsConfig, SgsLattice3, SgsModel, SgsStepError, SgsUse};
use std::sync::atomic::{AtomicUsize, Ordering};

const DIMS: [usize; 3] = [5, 4, 5];
const N: usize = 100;
const WALLS: [[f64; 3]; 2] = [[0.01, -0.004, 0.0], [0.018, 0.003, 0.0]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 103,
                kernel_id: 151,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn law() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Cumulant {
            rates: Rates27 {
                omega: [1.97, 1.13, 1.21, 1.34, 1.47, 1.58, 1.66, 1.72, 1.39, 1.51],
            },
        },
        maximum_mach: 0.2,
    }
}

fn noise(i: usize, salt: u64) -> f64 {
    let mut h = (i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
    h ^= h >> 30;
    h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}

fn forces(epoch: f64) -> Vec<[f64; 3]> {
    (0..N)
        .map(|i| {
            [
                1e-5 * epoch * (2.0 + noise(i, 71)),
                1e-5 * epoch * (noise(i, 93) - 0.2),
                0.0,
            ]
        })
        .collect()
}

// All force frames below represent the same physical density and velocity.
fn grid<L: Lattice3<Q>, const Q: usize>(
    collision: L::Collision,
    frame: &[[f64; 3]],
    uniform_force: [f64; 3],
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) -> Grid3<L, Q> {
    Grid3::try_new(
        DIMS,
        collision,
        uniform_force,
        [true, true, false],
        GridBudget3 {
            maximum_cells: N,
            maximum_mach: 0.2,
        },
        |_| false,
        |[x, y, z]| {
            let i = (z * DIMS[1] + y) * DIMS[0] + x;
            let h = (z as f64 + 0.5) / DIMS[2] as f64;
            let envelope = 4.0 * h * (1.0 - h);
            let u = core::array::from_fn(|a| {
                WALLS[0][a]
                    + h * (WALLS[1][a] - WALLS[0][a])
                    + 0.02 * envelope * noise(i, 31 * a as u64 + 7)
                    - 0.5 * frame[i][a]
            });
            equilibrium(1.0, u)
        },
    )
    .unwrap()
}

fn forced(cx: &Cx<'_>, field: Vec<[f64; 3]>) -> ForcedGrid3<D3q27, 27> {
    let populations = grid::<D3q27, 27>(law(), &field, [0.0; 3], |rho, u| {
        equilibrium27(rho, u, Equilibrium27::Factorized)
    });
    ForcedGrid3::try_new(populations, cx, field).unwrap()
}

fn wall_velocities<L: Lattice3<Q>, const Q: usize>(
    probe: &WallMomentumProbe3<L, Q>,
) -> Vec<[f64; 3]> {
    probe
        .wall_links()
        .map(|link| WALLS[usize::from(link.exterior_faces[5])])
        .collect()
}

fn molecular_field() -> Vec<CollisionConfig27> {
    (0..N)
        .map(|i| {
            let mut collision = law();
            if let CollisionModel27::Cumulant { rates } = &mut collision.model {
                rates.omega[0] -= 0.02 * (i % 3) as f64;
                rates.omega[6] += 0.01 * (i % 2) as f64;
            }
            collision
        })
        .collect()
}

fn local(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn dynamic() -> ChannelDynamicSgsConfig {
    ChannelDynamicSgsConfig {
        maximum_coefficient_squared: 10.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn bits<const Q: usize>(field: &[[f64; Q]]) -> Vec<u64> {
    field.iter().flatten().map(|v| v.to_bits()).collect()
}

fn momentum(flow: &ForcedGrid3<D3q27, 27>) -> [f64; 3] {
    let mut total = [0.0; 3];
    for i in 0..N {
        let (rho, u) = flow.macroscopic(i).unwrap();
        for a in 0..3 {
            total[a] += rho * u[a];
        }
    }
    total
}

fn near_collisions(a: &[CollisionConfig27], b: &[CollisionConfig27]) {
    for (a, b) in a.iter().zip(b) {
        let a = D3q27::molecular_viscosity(*a).unwrap();
        let b = D3q27::molecular_viscosity(*b).unwrap();
        assert!((a - b).abs() < 3e-13, "{a} {b}");
    }
}

#[test]
fn native_wall_steps_keep_uniform_force_parity_and_trapezoidal_momentum() {
    fn uniform<L: ForcedLattice3<Q>, const Q: usize>(
        collision: L::Collision,
        equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
    ) {
        with_cx(|cx| {
            let force = [1e-5, -2e-5, 0.0];
            let field = vec![force; N];
            let mut native = grid::<L, Q>(collision, &field, force, equilibrium);
            let mut spatial = ForcedGrid3::try_new(
                grid::<L, Q>(collision, &field, [0.0; 3], equilibrium),
                cx,
                field,
            )
            .unwrap();
            let probe = native.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
            let velocities = wall_velocities(&probe);
            let collisions = vec![collision; N];
            for _ in 0..2 {
                let expected = native
                    .step_serial_with_moving_wall_momentum(
                        cx,
                        &probe,
                        &velocities,
                        Some(&collisions),
                    )
                    .unwrap();
                let observed = spatial
                    .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, &collisions)
                    .unwrap();
                assert_eq!(observed, expected);
                assert_eq!(bits(spatial.populations()), bits(native.populations()));
            }
        });
    }
    uniform::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
    uniform::<D3q19, 19>(
        CollisionModel3::CentralMoment {
            second_order_rate: 1.25,
            higher_order_rate: 1.1,
        },
        equilibrium3,
    );
    uniform::<D3q27, 27>(law(), |rho, u| {
        equilibrium27(rho, u, Equilibrium27::Factorized)
    });

    with_cx(|cx| {
        let mut driven = forced(cx, forces(1.0));
        let mut stationary = forced(cx, forces(1.0));
        let mut old = forced(cx, forces(1.0));
        let probe = driven.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
        let velocities = wall_velocities(&probe);
        let zeros = vec![[0.0; 3]; probe.measured_links()];
        let collisions = molecular_field();
        for epoch in [-3.0, 0.2, 5.0] {
            let next = forces(epoch);
            let before = momentum(&driven);
            let body: [f64; 3] = core::array::from_fn(|a| {
                driven
                    .forces()
                    .iter()
                    .zip(&next)
                    .map(|(a0, a1)| 0.5 * (a0[a] + a1[a]))
                    .sum()
            });
            let receipt = driven
                .advance_serial_with_moving_wall_momentum(
                    cx,
                    &probe,
                    &velocities,
                    &collisions,
                    next.clone(),
                )
                .unwrap();
            let after = momentum(&driven);
            for a in 0..3 {
                assert!(
                    (after[a] - before[a] - body[a] - receipt.fluid_population_impulse[a]).abs()
                        < 3e-12
                );
            }
            assert!(receipt.wall_work.abs() > 1e-10);
            old.advance_serial(cx, &collisions, next.clone()).unwrap();
            stationary
                .advance_serial_with_moving_wall_momentum(
                    cx,
                    &probe,
                    &zeros,
                    &collisions,
                    next.clone(),
                )
                .unwrap();
            assert_eq!(bits(stationary.populations()), bits(old.populations()));
            assert_eq!(driven.forces(), next);
            assert_eq!(stationary.forces(), old.forces());
            assert_eq!(stationary.completed_steps(), old.completed_steps());
        }
    });
}

#[test]
fn moving_wall_sgs_uses_the_incoming_force_frame_and_zero_drives_keep_old_bits() {
    with_cx(|cx| {
        let reference = forced(cx, forces(0.0));
        let shifted = forced(cx, forces(40.0));
        let probe = shifted.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
        let velocities = wall_velocities(&probe);
        let zeros = vec![[0.0; 3]; probe.measured_links()];
        let mut a = reference.sgs_workspace().unwrap();
        let mut b = shifted.sgs_workspace().unwrap();
        for model in [
            SgsModel::Smagorinsky { cs: 0.12 },
            SgsModel::Wale { cw: 0.5 },
            SgsModel::Vreman { cv: 0.07 },
        ] {
            let expected = reference
                .prepare_channel_sgs_with_moving_walls(
                    cx,
                    local(model),
                    None,
                    &probe,
                    &velocities,
                    &mut a,
                )
                .unwrap();
            let observed = shifted
                .prepare_channel_sgs_with_moving_walls(
                    cx,
                    local(model),
                    None,
                    &probe,
                    &velocities,
                    &mut b,
                )
                .unwrap();
            assert!(expected.maximum_eddy_viscosity_lu > 0.0);
            near_collisions(a.collisions(), b.collisions());
            assert!((expected.sgs_dissipation_lu - observed.sgs_dissipation_lu).abs() < 3e-13);

            let stationary = shifted.prepare_sgs(cx, local(model), None, &mut a).unwrap();
            let moving_zero = shifted
                .prepare_channel_sgs_with_moving_walls(
                    cx,
                    local(model),
                    None,
                    &probe,
                    &zeros,
                    &mut b,
                )
                .unwrap();
            assert_eq!(stationary, moving_zero);
            assert_eq!(a.collisions(), b.collisions());
        }
        let mut a = reference.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
        let mut b = shifted.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
        let expected = reference
            .prepare_channel_dynamic_sgs_with_moving_walls(
                cx,
                dynamic(),
                None,
                &probe,
                &velocities,
                &mut a,
            )
            .unwrap();
        let observed = shifted
            .prepare_channel_dynamic_sgs_with_moving_walls(
                cx,
                dynamic(),
                None,
                &probe,
                &velocities,
                &mut b,
            )
            .unwrap();
        assert!(expected.fit.maximum_eddy_viscosity > 0.0);
        near_collisions(a.collisions(), b.collisions());
        assert!(
            (expected.dissipation.sgs_dissipation_lu - observed.dissipation.sgs_dissipation_lu)
                .abs()
                < 3e-13
        );
        let stationary = shifted
            .prepare_channel_dynamic_sgs(cx, dynamic(), None, &mut a)
            .unwrap();
        let moving_zero = shifted
            .prepare_channel_dynamic_sgs_with_moving_walls(
                cx,
                dynamic(),
                None,
                &probe,
                &zeros,
                &mut b,
            )
            .unwrap();
        assert_eq!(stationary, moving_zero);
        assert_eq!(a.collisions(), b.collisions());

        // Changing only F_(n+1) changes the published physical frame, never the
        // incoming fit, kinetic collision, wall reflection or its impulse.
        let mut first = forced(cx, forces(1.0));
        let mut second = forced(cx, forces(1.0));
        let first_report = first
            .advance_serial_with_channel_dynamic_sgs_and_moving_walls(
                cx,
                dynamic(),
                None,
                &probe,
                &velocities,
                &mut a,
                forces(-4.0),
            )
            .unwrap();
        let second_report = second
            .advance_serial_with_channel_dynamic_sgs_and_moving_walls(
                cx,
                dynamic(),
                None,
                &probe,
                &velocities,
                &mut b,
                forces(7.0),
            )
            .unwrap();
        assert_eq!(first_report, second_report);
        assert_eq!(bits(first.populations()), bits(second.populations()));
        assert_ne!(first.forces(), second.forces());
    });
}

#[test]
fn local_and_dynamic_sgs_match_native_moving_steps_and_all_worker_counts() {
    with_cx(|cx| {
        // The same check exercises each public composition, while its oracle
        // explicitly prepares the field and calls the original native owner.
        macro_rules! compare {
            ($config:expr, $active:expr, $workspace:ident($($arg:expr),*),
             $prepare:ident, $held_serial:ident, $held_pooled:ident,
             $serial:ident, $pooled:ident) => {
                for workers in [1, 3] {
                    let mut native = forced(cx, forces(1.0));
                    let mut serial = forced(cx, forces(1.0));
                    let mut pooled = forced(cx, forces(1.0));
                    let probe = native.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
                    let velocities = wall_velocities(&probe);
                    let molecular = molecular_field();
                    let mut a = native.$workspace($($arg),*).unwrap();
                    let mut b = serial.$workspace($($arg),*).unwrap();
                    let mut c = pooled.$workspace($($arg),*).unwrap();
                    let pool = TilePool::for_host(workers, 103);
                    for (step, epoch) in [1.0, -2.0].into_iter().enumerate() {
                        let expected_fit = native.$prepare(
                            cx, $config, Some(&molecular), &probe, &velocities, &mut a,
                        ).unwrap();
                        if $active {
                            if step == 0 {
                                assert!(a.collisions().iter().zip(&molecular).any(|(a, b)| {
                                    D3q27::molecular_viscosity(*a).unwrap()
                                        > D3q27::molecular_viscosity(*b).unwrap()
                                }));
                            }
                        } else {
                            assert_eq!(a.collisions(), molecular.as_slice());
                        }
                        for (actual, original) in a.collisions().iter().zip(&molecular) {
                            let (
                                CollisionModel27::Cumulant { rates: actual },
                                CollisionModel27::Cumulant { rates: original },
                            ) = (actual.model, original.model) else {
                                panic!("SGS changed the molecular collision family");
                            };
                            assert_eq!(actual.omega[1..], original.omega[1..]);
                        }
                        let next = forces(epoch);
                        let expected_wall = if step == 0 {
                            native.step_serial_with_moving_wall_momentum(
                                cx, &probe, &velocities, a.collisions(),
                            ).unwrap()
                        } else {
                            native.advance_serial_with_moving_wall_momentum(
                                cx, &probe, &velocities, a.collisions(), next.clone(),
                            ).unwrap()
                        };
                        let expected = (expected_fit, expected_wall);
                        let actual = if step == 0 {
                            serial.$held_serial(
                                cx, $config, Some(&molecular), &probe, &velocities, &mut b,
                            ).unwrap()
                        } else {
                            serial.$serial(
                                cx, $config, Some(&molecular), &probe, &velocities, &mut b, next.clone(),
                            ).unwrap()
                        };
                        let parallel = if step == 0 {
                            pooled.$held_pooled(
                                &pool, &CancelGate::new_clock_free(), $config, Some(&molecular),
                                &probe, &velocities, &mut c,
                            ).unwrap()
                        } else {
                            pooled.$pooled(
                                &pool, &CancelGate::new_clock_free(), $config, Some(&molecular),
                                &probe, &velocities, &mut c, next.clone(),
                            ).unwrap()
                        };
                        assert_eq!(actual, expected);
                        assert_eq!(parallel, expected);
                        assert_eq!(bits(serial.populations()), bits(native.populations()));
                        assert_eq!(bits(pooled.populations()), bits(native.populations()));
                        assert_eq!(serial.forces(), next);
                        assert_eq!(pooled.forces(), next);
                        assert_eq!(serial.completed_steps(), native.completed_steps());
                        assert_eq!(pooled.completed_steps(), native.completed_steps());
                    }
                }
            };
        }
        for (model, active) in [
            (SgsModel::Smagorinsky { cs: 0.0 }, false),
            (SgsModel::Wale { cw: 0.0 }, false),
            (SgsModel::Vreman { cv: 0.0 }, false),
            (SgsModel::Smagorinsky { cs: 0.12 }, true),
            (SgsModel::Wale { cw: 0.5 }, true),
            (SgsModel::Vreman { cv: 0.07 }, true),
        ] {
            compare!(
                local(model),
                active,
                sgs_workspace(),
                prepare_channel_sgs_with_moving_walls,
                step_serial_with_channel_sgs_and_moving_walls,
                step_pooled_with_channel_sgs_and_moving_walls,
                advance_serial_with_channel_sgs_and_moving_walls,
                advance_pooled_with_channel_sgs_and_moving_walls
            );
        }
        compare!(
            dynamic(),
            true,
            channel_dynamic_sgs_workspace(usize::MAX),
            prepare_channel_dynamic_sgs_with_moving_walls,
            step_serial_with_channel_dynamic_sgs_and_moving_walls,
            step_pooled_with_channel_dynamic_sgs_and_moving_walls,
            advance_serial_with_channel_dynamic_sgs_and_moving_walls,
            advance_pooled_with_channel_dynamic_sgs_and_moving_walls
        );
    });
}

#[test]
fn refused_endpoints_and_post_stream_cancellation_keep_the_whole_force_epoch() {
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
        let current = forces(1.0);
        let mut flow = forced(cx, current.clone());
        let before = bits(flow.populations());
        let probe = flow.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
        let velocities = wall_velocities(&probe);
        let molecular = molecular_field();
        let mut workspace = flow.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
        let mut nonfinite = forces(-2.0);
        nonfinite[N - 1][0] = f64::NAN;
        for next in [vec![[0.0; 3]; N - 1], nonfinite, vec![[1.0, 0.0, 0.0]; N]] {
            assert!(
                flow.advance_serial_with_channel_dynamic_sgs_and_moving_walls(
                    cx,
                    dynamic(),
                    Some(&molecular),
                    &probe,
                    &velocities,
                    &mut workspace,
                    next,
                )
                .is_err()
            );
            assert_eq!(bits(flow.populations()), before);
            assert_eq!(flow.forces(), current);
            assert_eq!(flow.completed_steps(), 0);
        }
        let runner = CancelAfterStream {
            pool: TilePool::for_host(3, 103),
            calls: AtomicUsize::new(0),
        };
        let gate = CancelGate::new_clock_free();
        let next = forces(-2.0);
        assert!(matches!(
            flow.advance_pooled_with_channel_dynamic_sgs_and_moving_walls(
                &runner,
                &gate,
                dynamic(),
                Some(&molecular),
                &probe,
                &velocities,
                &mut workspace,
                next.clone(),
            ),
            Err(SgsStepError::Grid(GridError3::Cancelled))
        ));
        assert_eq!(runner.calls.load(Ordering::Acquire), 2);
        assert_eq!(bits(flow.populations()), before);
        assert_eq!(flow.forces(), current);
        assert_eq!(flow.completed_steps(), 0);
        let observed = flow
            .advance_pooled_with_channel_dynamic_sgs_and_moving_walls(
                &runner.pool,
                &CancelGate::new_clock_free(),
                dynamic(),
                Some(&molecular),
                &probe,
                &velocities,
                &mut workspace,
                next.clone(),
            )
            .unwrap();
        let mut clean = forced(cx, current);
        let mut clean_workspace = clean.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
        let expected = clean
            .advance_serial_with_channel_dynamic_sgs_and_moving_walls(
                cx,
                dynamic(),
                Some(&molecular),
                &probe,
                &velocities,
                &mut clean_workspace,
                next.clone(),
            )
            .unwrap();
        assert_eq!(observed, expected);
        assert_eq!(bits(flow.populations()), bits(clean.populations()));
        assert_eq!(flow.forces(), next);
        assert_eq!(flow.completed_steps(), 1);
    });
}
