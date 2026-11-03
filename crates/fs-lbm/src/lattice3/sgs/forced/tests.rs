//! Native force-frame, impulse and publication regressions.

use super::*;
use crate::d3q19::CollisionModel3;
use crate::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use crate::lattice3::{D3q19, GridBudget3};
use crate::sgs3::dynamic::{DynamicAveraging, DynamicSgsError};
use crate::sgs3::{SgsModel, SgsUse};
use fs_exec::{RunError, RunReport, TileKernel, TilePool};
use std::sync::atomic::{AtomicUsize, Ordering};

const DIMS: [usize; 3] = [5, 4, 5];
const N: usize = 100;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 103,
                kernel_id: 147,
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

fn force_field(epoch: f64) -> Vec<[f64; 3]> {
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

// Different stored force frames represent the SAME physical rho and u.
fn grid(channel: bool, frame: &[[f64; 3]], uniform: [f64; 3]) -> Grid3<D3q27, 27> {
    Grid3::try_new(
        DIMS,
        law(),
        uniform,
        [true, true, !channel],
        GridBudget3 {
            maximum_cells: N,
            maximum_mach: 0.2,
        },
        |_| false,
        |[x, y, z]| {
            let i = (z * DIMS[1] + y) * DIMS[0] + x;
            let h = z as f64 + 0.5;
            let envelope = if channel {
                4.0 * h * (5.0 - h) / 25.0
            } else {
                1.0
            };
            let u = core::array::from_fn(|a| {
                0.02 * envelope * noise(i, 31 * a as u64 + 7) - 0.5 * frame[i][a]
            });
            equilibrium27(1.0, u, Equilibrium27::Factorized)
        },
    )
    .unwrap()
}

fn forced(channel: bool, force: Vec<[f64; 3]>) -> ForcedGrid3<D3q27, 27> {
    with_cx(|cx| ForcedGrid3::try_new(grid(channel, &force, [0.0; 3]), cx, force).unwrap())
}
fn local(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}
fn dynamic() -> DynamicSgsConfig {
    DynamicSgsConfig {
        averaging: DynamicAveraging::LocalBox,
        maximum_coefficient_squared: 10.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}
fn channel() -> ChannelDynamicSgsConfig {
    ChannelDynamicSgsConfig {
        maximum_coefficient_squared: 10.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}
fn bits<const Q: usize>(field: &[[f64; Q]]) -> Vec<u64> {
    field.iter().flatten().map(|v| v.to_bits()).collect()
}
fn near_collisions(a: &[CollisionConfig27], b: &[CollisionConfig27]) {
    for (a, b) in a.iter().zip(b) {
        let x = D3q27::molecular_viscosity(*a).unwrap();
        let y = D3q27::molecular_viscosity(*b).unwrap();
        assert!((x - y).abs() < 3e-13, "{x} {y}");
    }
}

#[test]
fn g1_all_closures_use_the_current_spatial_force_frame() {
    for walls in [false, true] {
        let reference = forced(walls, force_field(0.0));
        let shifted = forced(walls, force_field(40.0));
        for i in 0..N {
            let (_, a) = reference.macroscopic(i).unwrap();
            let (_, b) = shifted.macroscopic(i).unwrap();
            for axis in 0..3 {
                assert!((a[axis] - b[axis]).abs() < 2e-15);
            }
        }
        with_cx(|cx| {
            let mut a = reference.sgs_workspace().unwrap();
            let mut b = shifted.sgs_workspace().unwrap();
            for model in [
                SgsModel::Smagorinsky { cs: 0.12 },
                SgsModel::Wale { cw: 0.5 },
                SgsModel::Vreman { cv: 0.07 },
            ] {
                let expected = reference
                    .prepare_sgs(cx, local(model), None, &mut a)
                    .unwrap();
                let observed = shifted.prepare_sgs(cx, local(model), None, &mut b).unwrap();
                assert!(expected.maximum_eddy_viscosity_lu > 0.0);
                near_collisions(a.collisions(), b.collisions());
                assert!((expected.sgs_dissipation_lu - observed.sgs_dissipation_lu).abs() < 3e-13);
            }
            if walls {
                let mut a = reference.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
                let mut b = shifted.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
                let expected = reference
                    .prepare_channel_dynamic_sgs(cx, channel(), None, &mut a)
                    .unwrap();
                shifted
                    .prepare_channel_dynamic_sgs(cx, channel(), None, &mut b)
                    .unwrap();
                assert!(expected.fit.maximum_eddy_viscosity > 0.0);
                near_collisions(a.collisions(), b.collisions());
            } else {
                let mut a = reference.dynamic_sgs_workspace(usize::MAX).unwrap();
                let mut b = shifted.dynamic_sgs_workspace(usize::MAX).unwrap();
                let expected = reference
                    .prepare_dynamic_sgs(cx, dynamic(), None, &mut a)
                    .unwrap();
                shifted
                    .prepare_dynamic_sgs(cx, dynamic(), None, &mut b)
                    .unwrap();
                assert!(expected.fit.maximum_eddy_viscosity > 0.0);
                near_collisions(a.collisions(), b.collisions());
            }
        });
    }
}

#[test]
fn g0_g5_uniform_forces_keep_the_original_periodic_dynamic_step_bits() {
    let force = [1e-5, -2e-5, 0.0];
    let frame = vec![force; N];
    for workers in [1, 3] {
        let mut reference = grid(false, &frame, force);
        let mut serial = forced(false, frame.clone());
        let mut pooled = forced(false, frame.clone());
        let mut a = reference.dynamic_sgs_workspace(usize::MAX).unwrap();
        let mut b = serial.dynamic_sgs_workspace(usize::MAX).unwrap();
        let mut c = pooled.dynamic_sgs_workspace(usize::MAX).unwrap();
        let pool = TilePool::for_host(workers, 103);
        with_cx(|cx| {
            for _ in 0..2 {
                let expected = reference
                    .step_serial_with_dynamic_sgs(cx, dynamic(), None, &mut a)
                    .unwrap();
                assert_eq!(
                    expected,
                    serial
                        .step_serial_with_dynamic_sgs(cx, dynamic(), None, &mut b)
                        .unwrap()
                );
                assert_eq!(
                    expected,
                    pooled
                        .step_pooled_with_dynamic_sgs(
                            &pool,
                            &CancelGate::new_clock_free(),
                            dynamic(),
                            None,
                            &mut c
                        )
                        .unwrap()
                );
                assert_eq!(bits(reference.populations()), bits(serial.populations()));
                assert_eq!(bits(reference.populations()), bits(pooled.populations()));
            }
        });
    }
}

fn momentum(flow: &ForcedGrid3<D3q27, 27>) -> [f64; 3] {
    let mut result = [0.0; 3];
    for i in 0..N {
        let (rho, u) = flow.macroscopic(i).unwrap();
        for a in 0..3 {
            result[a] += rho * u[a];
        }
    }
    result
}

#[test]
fn g1_g5_active_dynamic_transport_keeps_trapezoidal_impulse_and_worker_identity() {
    for workers in [1, 2, 4] {
        let mut serial = forced(false, force_field(1.0));
        let mut pooled = forced(false, force_field(1.0));
        let mut a = serial.dynamic_sgs_workspace(usize::MAX).unwrap();
        let mut b = pooled.dynamic_sgs_workspace(usize::MAX).unwrap();
        let pool = TilePool::for_host(workers, 103);
        for (step, epoch) in [2.0, -1.0, 0.0].into_iter().enumerate() {
            let before = momentum(&serial);
            let end = force_field(epoch);
            let mut expected = [0.0; 3];
            for (start, end) in serial.forces().iter().zip(&end) {
                for axis in 0..3 {
                    expected[axis] += 0.5 * (start[axis] + end[axis]);
                }
            }
            let report = with_cx(|cx| {
                serial
                    .advance_serial_with_dynamic_sgs(cx, dynamic(), None, &mut a, end.clone())
                    .unwrap()
            });
            let other = pooled
                .advance_pooled_with_dynamic_sgs(
                    &pool,
                    &CancelGate::new_clock_free(),
                    dynamic(),
                    None,
                    &mut b,
                    end,
                )
                .unwrap();
            if step == 0 {
                assert!(report.fit.maximum_eddy_viscosity > 0.0);
            }
            assert_eq!(report, other);
            assert_eq!(bits(serial.populations()), bits(pooled.populations()));
            assert_eq!(bits(serial.forces()), bits(pooled.forces()));
            let after = momentum(&serial);
            for axis in 0..3 {
                assert!((after[axis] - before[axis] - expected[axis]).abs() < 2e-11);
            }
        }
    }
}

#[test]
fn g0_local_zero_control_preserves_both_original_force_owners() {
    fn run<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(collision: L::Collision) {
        let make = || {
            with_cx(|cx| {
                let grid = Grid3::<L, Q>::try_new(
                    DIMS,
                    collision,
                    [0.0; 3],
                    [true; 3],
                    GridBudget3 {
                        maximum_cells: N,
                        maximum_mach: 0.2,
                    },
                    |_| false,
                    |_| L::W,
                )
                .unwrap();
                ForcedGrid3::try_new(grid, cx, force_field(1.0)).unwrap()
            })
        };
        let mut expected = make();
        let mut actual = make();
        let mut workspace = actual.sgs_workspace().unwrap();
        with_cx(|cx| {
            expected
                .advance_serial(cx, &vec![collision; N], force_field(-1.0))
                .unwrap();
            let report = actual
                .advance_serial_with_sgs(
                    cx,
                    local(SgsModel::Smagorinsky { cs: 0.0 }),
                    None,
                    &mut workspace,
                    force_field(-1.0),
                )
                .unwrap();
            assert_eq!(report.maximum_eddy_viscosity_lu, 0.0);
        });
        assert_eq!(bits(expected.populations()), bits(actual.populations()));
        assert_eq!(bits(expected.forces()), bits(actual.forces()));
    }
    run::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    run::<D3q27, 27>(law());
}

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

#[test]
fn g0_g4_channel_refusals_keep_source_history_and_cancelled_passes_retry_exactly() {
    let mut expected = forced(true, force_field(1.0));
    let mut fresh = expected.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
    let report = with_cx(|cx| {
        expected
            .advance_serial_with_channel_dynamic_sgs(
                cx,
                channel(),
                None,
                &mut fresh,
                force_field(-1.0),
            )
            .unwrap()
    });
    assert!(report.fit.maximum_eddy_viscosity > 0.0);
    for after in [0, 1] {
        let mut actual = forced(true, force_field(1.0));
        let initial = bits(actual.populations());
        let history = bits(actual.forces());
        assert!(actual.channel_dynamic_sgs_workspace(0).is_err());
        let mut work = actual.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
        with_cx(|cx| {
            let mut cap = channel();
            cap.maximum_coefficient_squared = 0.0;
            assert!(matches!(
                actual.advance_serial_with_channel_dynamic_sgs(
                    cx,
                    cap,
                    None,
                    &mut work,
                    force_field(-1.0)
                ),
                Err(SgsStepError::Dynamic(DynamicSgsError::Coefficient { .. }))
            ));
            cap = channel();
            cap.maximum_eddy_to_molecular = 0.0;
            assert!(matches!(
                actual.advance_serial_with_channel_dynamic_sgs(
                    cx,
                    cap,
                    None,
                    &mut work,
                    force_field(-1.0)
                ),
                Err(SgsStepError::Closure(
                    SgsError::ViscosityRatioExceeded { .. }
                ))
            ));
            let mut end = force_field(-1.0);
            end[N - 1][0] = 1.0;
            assert!(matches!(
                actual.advance_serial_with_channel_dynamic_sgs(cx, channel(), None, &mut work, end),
                Err(SgsStepError::Grid(GridError3::Mach { .. }))
            ));
        });
        assert_eq!(bits(actual.populations()), initial);
        assert_eq!(bits(actual.forces()), history);
        assert_eq!(actual.completed_steps(), 0);
        let pool = CancelPass {
            pool: TilePool::for_host(2, 103),
            calls: AtomicUsize::new(0),
            after,
        };
        let gate = CancelGate::new_clock_free();
        assert!(matches!(
            actual.advance_pooled_with_channel_dynamic_sgs(
                &pool,
                &gate,
                channel(),
                None,
                &mut work,
                force_field(-1.0)
            ),
            Err(SgsStepError::Grid(GridError3::Cancelled))
        ));
        assert_eq!(bits(actual.populations()), initial);
        assert_eq!(bits(actual.forces()), history);
        assert_eq!(actual.completed_steps(), 0);
        let retry = actual
            .advance_pooled_with_channel_dynamic_sgs(
                &pool.pool,
                &CancelGate::new_clock_free(),
                channel(),
                None,
                &mut work,
                force_field(-1.0),
            )
            .unwrap();
        assert_eq!(retry, report);
        assert_eq!(bits(actual.populations()), bits(expected.populations()));
        assert_eq!(bits(actual.forces()), bits(expected.forces()));
        assert_eq!(actual.completed_steps(), 1);
    }
}
