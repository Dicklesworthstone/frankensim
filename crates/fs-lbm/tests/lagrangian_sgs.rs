//! Native history/flow publication and independent collision regressions.

use fs_exec::{CancelGate, Cx, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{CollisionModel3, equilibrium3};
use fs_lbm::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, Lattice3, ViscousLattice3};
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
            fs_exec::StreamKey {
                seed: 501,
                kernel_id: 136,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn grid<L: Lattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
    varying: bool,
) -> Grid3<L, Q> {
    grid_with_force::<L, Q>(collision, equilibrium, varying, [0.0; 3])
}

fn grid_with_force<L: Lattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
    varying: bool,
    force: [f64; 3],
) -> Grid3<L, Q> {
    Grid3::try_new(
        [5; 3],
        collision,
        force,
        [true; 3],
        GridBudget3 {
            maximum_cells: 125,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let [x, y, z] = p.map(|i| core::f64::consts::TAU * (i as f64 + 0.5) / 5.0);
            let u = if varying {
                [
                    0.003 * x.sin() + 0.002 * (y + z).sin(),
                    0.002 * y.cos() + 0.003 * (z + x).sin(),
                    0.002 * z.sin() + 0.003 * (x + y).cos(),
                ]
            } else {
                [0.01, -0.005, 0.002]
            };
            // These are stored populations. For nonzero force the physical
            // velocity is u + F/(2 rho), which the reference checks below.
            equilibrium(1.0 + if varying { 0.002 * x.sin() } else { 0.0 }, u)
        },
    )
    .unwrap()
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

fn native_reference<L: SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
    force: [f64; 3],
) {
    let mut native = grid_with_force::<L, Q>(collision, equilibrium, true, force);
    let mut independent = grid_with_force::<L, Q>(collision, equilibrium, true, force);
    let mut pooled = grid_with_force::<L, Q>(collision, equilibrium, true, force);
    let mut work = native.lagrangian_sgs_workspace(usize::MAX).unwrap();
    let mut pool_work = pooled.lagrangian_sgs_workspace(usize::MAX).unwrap();
    let mut model = LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX).unwrap();
    let pool = TilePool::for_host(2, 501);
    let gate = CancelGate::new_clock_free();
    let nu = L::molecular_viscosity(collision).unwrap();
    with_cx(|cx| {
        for step in 0..4 {
            let macros: Vec<_> = (0..125)
                .map(|i| independent.macroscopic(i).unwrap())
                .collect();
            if force != [0.0; 3] {
                for (populations, &(rho, physical)) in independent.populations().iter().zip(&macros)
                {
                    let raw_momentum: [f64; 3] = core::array::from_fn(|axis| {
                        populations
                            .iter()
                            .zip(L::E)
                            .map(|(&f, e)| f64::from(e[axis]) * f)
                            .sum()
                    });
                    for axis in 0..3 {
                        let raw = raw_momentum[axis] / rho;
                        let half_force = 0.5 * force[axis] / rho;
                        assert!(
                            (physical[axis] - raw - half_force).abs() < 1e-16,
                            "step {step}, axis {axis}: missing physical half-force velocity"
                        );
                        assert!((physical[axis] - raw).abs() > 1e-5);
                    }
                }
            }
            let density: Vec<_> = macros.iter().map(|m| m.0).collect();
            let velocity: Vec<_> = macros.iter().map(|m| m.1).collect();
            let trial = model
                .prepare(config(), step as f64, &density, &velocity, || false)
                .unwrap();
            let rates: Vec<_> = trial
                .viscosities()
                .iter()
                .map(|&eddy| {
                    if eddy == 0.0 {
                        collision
                    } else {
                        L::with_shear_relaxation(collision, 0.5 + (nu + eddy) / L::CS2).unwrap()
                    }
                })
                .collect();
            independent.step_serial_with_collisions(cx, &rates).unwrap();
            let expected_fit = trial.report();
            trial.commit();
            let report = native
                .step_serial_with_lagrangian_sgs(cx, config(), None, &mut work)
                .unwrap();
            let parallel = pooled
                .step_pooled_with_lagrangian_sgs(&pool, &gate, config(), None, &mut pool_work)
                .unwrap();
            assert_eq!(report.fit, expected_fit);
            assert_eq!(report, parallel);
            assert_eq!(work.accepted_fit_time(), Some(step as f64));
            assert_eq!(native.populations(), independent.populations());
            assert_eq!(native.populations(), pooled.populations());
            assert!(report.fit.maximum_eddy_viscosity > 0.0);
        }
    });
}

#[test]
fn d3q19_and_d3q27_match_independently_prepared_history_and_worker_counts() {
    native_reference::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3, [0.0; 3]);
    native_reference::<D3q27, 27>(cumulant(), eq27, [0.0; 3]);
}

#[test]
fn d3q19_and_d3q27_forced_history_uses_guo_centered_velocity() {
    // Variable density makes F/(2 rho) spatially varying even for this uniform
    // force. All four accepted fits and departures must use that velocity,
    // including the history transported after the first local-box fit.
    let force = [0.001, -0.0006, 0.0003];
    native_reference::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3, force);
    native_reference::<D3q27, 27>(cumulant(), eq27, force);
}

#[test]
fn native_cumulant_retains_molecular_field_and_all_nonshear_rates() {
    let law = cumulant();
    let mut flow = grid::<D3q27, 27>(law, eq27, true);
    let mut work = flow.lagrangian_sgs_workspace(usize::MAX).unwrap();
    let molecular: Vec<_> = (0..125)
        .map(|i| D3q27::with_shear_relaxation(law, 0.6 + 0.001 * i as f64).unwrap())
        .collect();
    let original = molecular.clone();
    with_cx(|cx| {
        flow.step_serial_with_lagrangian_sgs(cx, config(), Some(&molecular), &mut work)
            .unwrap();
    });
    for (actual, source) in work.collisions().iter().zip(&molecular) {
        let CollisionModel27::Cumulant { rates: a } = actual.model else {
            panic!("collision changed");
        };
        let CollisionModel27::Cumulant { rates: b } = source.model else {
            unreachable!();
        };
        assert_eq!(a.omega[1..], b.omega[1..]);
        assert!(
            D3q27::molecular_viscosity(*actual).unwrap()
                >= D3q27::molecular_viscosity(*source).unwrap() - 1e-16
        );
    }
    assert_eq!(molecular, original);
}

struct RefuseNativeStep {
    pool: TilePool,
    calls: AtomicUsize,
}
impl KernelRunner for RefuseNativeStep {
    fn workers(&self) -> usize {
        self.pool.workers()
    }
    fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        let outcome = self.pool.run_with_gate(kernel, gate);
        if self.calls.fetch_add(1, Ordering::AcqRel) == 1 {
            gate.request();
        }
        outcome
    }
}

#[test]
fn late_native_refusal_keeps_history_populations_and_clock_retryable() {
    let law = CollisionModel3::Bgk { tau: 0.8 };
    let mut flow = grid::<D3q19, 19>(law, equilibrium3, true);
    let mut control = grid::<D3q19, 19>(law, equilibrium3, true);
    let mut work = flow.lagrangian_sgs_workspace(usize::MAX).unwrap();
    let mut clean_work = control.lagrangian_sgs_workspace(usize::MAX).unwrap();
    with_cx(|cx| {
        flow.step_serial_with_lagrangian_sgs(cx, config(), None, &mut work)
            .unwrap();
        control
            .step_serial_with_lagrangian_sgs(cx, config(), None, &mut clean_work)
            .unwrap();
        let before = flow.populations().to_vec();
        let gate = CancelGate::new_clock_free();
        let refusing = RefuseNativeStep {
            pool: TilePool::for_host(2, 501),
            calls: AtomicUsize::new(0),
        };
        assert!(
            flow.step_pooled_with_lagrangian_sgs(&refusing, &gate, config(), None, &mut work)
                .is_err()
        );
        assert_eq!(flow.populations(), before);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(work.accepted_fit_time(), Some(0.0));
        assert!(
            flow.step_serial_with_lagrangian_sgs(
                cx,
                LagrangianSgsConfig {
                    maximum_eddy_to_molecular: 0.0,
                    ..config()
                },
                None,
                &mut work
            )
            .is_err()
        );
        assert_eq!(work.accepted_fit_time(), Some(0.0));
        let retry = flow
            .step_serial_with_lagrangian_sgs(cx, config(), None, &mut work)
            .unwrap();
        let expected = control
            .step_serial_with_lagrangian_sgs(cx, config(), None, &mut clean_work)
            .unwrap();
        assert_eq!(retry, expected);
        assert_eq!(flow.populations(), control.populations());
        assert_eq!(work.accepted_fit_time(), Some(1.0));
    });
}

#[test]
fn uniform_flow_is_original_step_exact_and_storage_is_admitted() {
    let law = CollisionModel3::Bgk { tau: 0.8 };
    let mut flow = grid::<D3q19, 19>(law, equilibrium3, false);
    let mut original = grid::<D3q19, 19>(law, equilibrium3, false);
    let mut work = flow.lagrangian_sgs_workspace(usize::MAX).unwrap();
    let bytes = work.workspace_bytes();
    assert!(flow.lagrangian_sgs_workspace(bytes - 1).is_err());
    assert_eq!(
        flow.lagrangian_sgs_workspace(bytes)
            .unwrap()
            .workspace_bytes(),
        bytes
    );
    with_cx(|cx| {
        for _ in 0..3 {
            let report = flow
                .step_serial_with_lagrangian_sgs(cx, config(), None, &mut work)
                .unwrap();
            original.step_serial(cx).unwrap();
            assert_eq!(report.fit.maximum_eddy_viscosity, 0.0);
            assert_eq!(flow.populations(), original.populations());
        }
    });
}
