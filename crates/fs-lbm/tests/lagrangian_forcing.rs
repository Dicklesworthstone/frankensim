//! Current-force Lagrangian fits and joint native source/history publication.

use fs_exec::{CancelGate, Cx, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{CollisionModel3, equilibrium3};
use fs_lbm::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use fs_lbm::lattice3::sources::{ForcedGrid3, ForcedLattice3};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3};
use fs_lbm::sgs3::dynamic::lagrangian::{LagrangianDynamicSgs, LagrangianSgsConfig};
use fs_lbm::sgs3::{SgsLattice3, SgsUse};
use std::sync::atomic::{AtomicUsize, Ordering};

fn config() -> LagrangianSgsConfig {
    LagrangianSgsConfig {
        memory_time: 4.0,
        maximum_departure_cells: 1.0,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey { seed: 501, kernel_id: 136, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn force_field(epoch: f64) -> Vec<[f64; 3]> {
    (0..125)
        .map(|cell| {
            let [x, y, z] = [cell % 5, cell / 5 % 5, cell / 25]
                .map(|i| core::f64::consts::TAU * (i as f64 + 0.5) / 5.0);
            [
                0.001 + 0.0008 * (x + y + epoch).sin(),
                -0.0006 + 0.0006 * (y + z - epoch).cos(),
                0.0003 + 0.0004 * (z + x + epoch).sin(),
            ]
        })
        .collect()
}

fn forced<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) -> ForcedGrid3<L, Q> {
    let grid = Grid3::<L, Q>::try_new(
        [5; 3],
        collision,
        [0.0; 3],
        [true; 3],
        GridBudget3 { maximum_cells: 125, maximum_mach: 0.2 },
        |_| false,
        |p| {
            let [x, y, z] = p.map(|i| core::f64::consts::TAU * (i as f64 + 0.5) / 5.0);
            equilibrium(
                1.0 + 0.002 * x.sin(),
                [
                    0.003 * x.sin() + 0.002 * (y + z).sin(),
                    0.002 * y.cos() + 0.003 * (z + x).sin(),
                    0.002 * z.sin() + 0.003 * (x + y).cos(),
                ],
            )
        },
    )
    .unwrap();
    ForcedGrid3::try_new(grid, cx, force_field(0.0)).unwrap()
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

fn eq27(rho: f64, u: [f64; 3]) -> [f64; 27] {
    equilibrium27(rho, u, Equilibrium27::Factorized)
}

fn reference<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) where
    L::Collision: core::fmt::Debug + PartialEq,
{
    with_cx(|cx| {
        let mut native = forced::<L, Q>(cx, collision, equilibrium);
        let mut independent = forced::<L, Q>(cx, collision, equilibrium);
        let mut pooled = forced::<L, Q>(cx, collision, equilibrium);
        let mut work = native.lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut pooled_work = pooled.lagrangian_sgs_workspace(usize::MAX).unwrap();
        assert!(native.lagrangian_sgs_workspace(work.workspace_bytes() - 1).is_err());
        let mut model = LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX).unwrap();
        let molecular: Vec<_> = (0..125)
            .map(|cell| L::with_shear_relaxation(collision, 0.7 + 0.0005 * cell as f64).unwrap())
            .collect();
        let original = molecular.clone();
        let pool = TilePool::for_host(2, 501);
        let gate = CancelGate::new_clock_free();
        let mut saw_modeled_transport = false;
        for step in 0..4 {
            let macros: Vec<_> = (0..125)
                .map(|cell| independent.macroscopic(cell).unwrap())
                .collect();
            let density: Vec<_> = macros.iter().map(|m| m.0).collect();
            let velocity: Vec<_> = macros.iter().map(|m| m.1).collect();
            let advance = step % 2 == 0;
            let next = if advance {
                force_field(step as f64 + 1.0)
            } else {
                independent.forces().to_vec()
            };
            let trial = model
                .prepare(config(), step as f64, &density, &velocity, || false)
                .unwrap();
            if step == 0 {
                // Explicitly distinguish the current and endpoint force charts.
                // Reinterpreting incoming populations with F_(n+1) must change
                // this fixture's fitted field, so the equality below is useful.
                let wrong_velocity: Vec<_> = (0..125)
                    .map(|cell| core::array::from_fn(|axis| {
                        velocity[cell][axis]
                            + 0.5 * (next[cell][axis] - independent.forces()[cell][axis])
                                / density[cell]
                    }))
                    .collect();
                let mut wrong = LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX).unwrap();
                let wrong_trial = wrong
                    .prepare(config(), 0.0, &density, &wrong_velocity, || false)
                    .unwrap();
                assert!(trial.viscosities().iter().zip(wrong_trial.viscosities())
                    .any(|(a, b)| (a - b).abs() > 1e-12));
            }
            let rates: Vec<_> = trial.viscosities().iter().enumerate()
                .map(|(cell, &eddy)| {
                    let source = molecular[cell];
                    if eddy == 0.0 {
                        source
                    } else {
                        let nu = L::molecular_viscosity(source).unwrap();
                        L::with_shear_relaxation(source, 0.5 + (nu + eddy) / L::CS2).unwrap()
                    }
                })
                .collect();
            let expected_fit = trial.report();
            let (serial, parallel) = if advance {
                independent.advance_serial(cx, &rates, next.clone()).unwrap();
                (
                    native.advance_serial_with_lagrangian_sgs(
                        cx, config(), Some(&molecular), &mut work, next.clone(),
                    ).unwrap(),
                    pooled.advance_pooled_with_lagrangian_sgs(
                        &pool, &gate, config(), Some(&molecular), &mut pooled_work, next.clone(),
                    ).unwrap(),
                )
            } else {
                independent.step_serial(cx, &rates).unwrap();
                (
                    native.step_serial_with_lagrangian_sgs(
                        cx, config(), Some(&molecular), &mut work,
                    ).unwrap(),
                    pooled.step_pooled_with_lagrangian_sgs(
                        &pool, &gate, config(), Some(&molecular), &mut pooled_work,
                    ).unwrap(),
                )
            };
            trial.commit();
            assert_eq!(serial.fit, expected_fit);
            assert_eq!(serial, parallel);
            saw_modeled_transport |= serial.fit.maximum_eddy_viscosity > 0.0;
            assert_eq!(work.collisions(), rates);
            assert_eq!(work.accepted_fit_time(), Some(step as f64));
            assert_eq!(pooled_work.accepted_fit_time(), Some(step as f64));
            assert_eq!(native.populations(), independent.populations());
            assert_eq!(native.populations(), pooled.populations());
            assert_eq!(native.forces(), next);
            assert_eq!(native.forces(), pooled.forces());
            assert_eq!(native.completed_steps(), step + 1);
        }
        assert!(saw_modeled_transport);
        assert_eq!(molecular, original);
    });
}

#[test]
fn d3q19_spatial_and_timed_forces_match_independent_lagrangian_history() {
    reference::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
}

#[test]
fn d3q27_spatial_and_timed_forces_preserve_cumulant_and_worker_identity() {
    reference::<D3q27, 27>(cumulant(), eq27);
}

struct RefuseAfterStreaming {
    pool: TilePool,
    calls: AtomicUsize,
}
impl KernelRunner for RefuseAfterStreaming {
    fn workers(&self) -> usize { self.pool.workers() }
    fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        let result = self.pool.run_with_gate(kernel, gate);
        if self.calls.fetch_add(1, Ordering::AcqRel) == 1 { gate.request(); }
        result
    }
}

fn rollback<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) {
    with_cx(|cx| {
        let mut flow = forced::<L, Q>(cx, collision, equilibrium);
        let mut control = forced::<L, Q>(cx, collision, equilibrium);
        let mut work = flow.lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut control_work = control.lagrangian_sgs_workspace(usize::MAX).unwrap();
        flow.advance_serial_with_lagrangian_sgs(
            cx, config(), None, &mut work, force_field(1.0),
        ).unwrap();
        control.advance_serial_with_lagrangian_sgs(
            cx, config(), None, &mut control_work, force_field(1.0),
        ).unwrap();
        let populations = flow.populations().to_vec();
        let forces = flow.forces().to_vec();
        assert!(flow.advance_serial_with_lagrangian_sgs(
            cx, config(), None, &mut work, vec![[1.0; 3]; 125],
        ).is_err());
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.forces(), forces);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(work.accepted_fit_time(), Some(0.0));
        let gate = CancelGate::new_clock_free();
        let refusing = RefuseAfterStreaming {
            pool: TilePool::for_host(2, 501), calls: AtomicUsize::new(0),
        };
        assert!(flow.advance_pooled_with_lagrangian_sgs(
            &refusing, &gate, config(), None, &mut work, force_field(2.0),
        ).is_err());
        assert_eq!(refusing.calls.load(Ordering::Acquire), 2);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.forces(), forces);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(work.accepted_fit_time(), Some(0.0));
        let retry = flow.advance_serial_with_lagrangian_sgs(
            cx, config(), None, &mut work, force_field(2.0),
        ).unwrap();
        let expected = control.advance_serial_with_lagrangian_sgs(
            cx, config(), None, &mut control_work, force_field(2.0),
        ).unwrap();
        assert_eq!(retry, expected);
        assert_eq!(flow.populations(), control.populations());
        assert_eq!(flow.forces(), control.forces());
        assert_eq!(flow.completed_steps(), 2);
        assert_eq!(work.accepted_fit_time(), Some(1.0));
    });
}

#[test]
fn rejected_endpoint_or_late_cancellation_keeps_force_and_memory_retryable() {
    rollback::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
    rollback::<D3q27, 27>(cumulant(), eq27);
}
