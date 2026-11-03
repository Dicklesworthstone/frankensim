//! G0/G1 native wall-channel checks, not a Reynolds-number validation campaign.

use super::*;
use crate::d3q19::CollisionModel3;
use crate::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use crate::lattice3::{D3q19, GridBudget3};
use crate::sgs3::SgsUse;
use fs_exec::TilePool;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 93,
                kernel_id: 146,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn config() -> ChannelDynamicSgsConfig {
    ChannelDynamicSgsConfig {
        maximum_coefficient_squared: 10.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
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

fn random(i: usize, salt: u64) -> f64 {
    let mut h = (i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
    h ^= h >> 30;
    h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}

fn channel(normal: usize, varying: bool, force_size: f64) -> Grid3<D3q27, 27> {
    let axes = [(normal + 1) % 3, (normal + 2) % 3, normal];
    let mut dims = [5; 3];
    dims[axes[1]] = 4;
    let mut force = [0.0; 3];
    force[axes[0]] = force_size;
    Grid3::try_new(
        dims,
        law(),
        force,
        core::array::from_fn(|axis| axis != normal),
        GridBudget3 {
            maximum_cells: 100,
            maximum_mach: 0.2,
        },
        |_| false,
        |p| {
            let [x, y, z] = axes.map(|axis| p[axis]);
            let active = (z * 4 + y) * 5 + x;
            let h = z as f64 + 0.5;
            let envelope = 4.0 * h * (5.0 - h) / 25.0;
            let mut velocity = [0.0; 3];
            if varying {
                for component in 0..3 {
                    velocity[axes[component]] =
                        0.02 * envelope * random(active, 31 * component as u64 + 7);
                }
            }
            // Initial population frame satisfies rho*u = sum(f*e) + F/2.
            for axis in 0..3 {
                velocity[axis] -= 0.5 * force[axis];
            }
            equilibrium27(1.0, velocity, Equilibrium27::Factorized)
        },
    )
    .unwrap()
}

fn zero<L: SgsLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let new = || {
        Grid3::<L, Q>::try_new(
            [3, 3, 2],
            collision,
            [0.0; 3],
            [true, true, false],
            GridBudget3 {
                maximum_cells: 18,
                maximum_mach: 0.2,
            },
            |_| false,
            |_| L::W,
        )
        .unwrap()
    };
    let mut original = new();
    let mut actual = new();
    let mut work = actual.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
    let molecular: Vec<_> = (0..18)
        .map(|cell| L::with_shear_relaxation(collision, 0.75 + 0.01 * (cell % 3) as f64).unwrap())
        .collect();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        original
            .step_serial_with_collisions(cx, &molecular)
            .unwrap();
        let result = actual
            .step_serial_with_channel_dynamic_sgs(cx, config(), Some(&molecular), &mut work)
            .unwrap();
        assert_eq!(result.fit.maximum_eddy_viscosity, 0.0);
        assert_eq!(result.dissipation.sgs_dissipation_lu, 0.0);
        assert_eq!(actual.populations(), original.populations());
        assert_eq!(actual.completed_steps(), 1);
    });
}

#[test]
fn g0_channel_zero_fit_retains_d3q19_and_cumulant_molecular_steps() {
    zero::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    zero::<D3q27, 27>(law());
}

#[test]
fn g1_forced_channel_cumulants_keep_other_rates_and_original_serial_pooled_transport() {
    for normal in 0..3 {
        let mut serial = channel(normal, true, 2e-5);
        let mut original = channel(normal, true, 2e-5);
        let mut pooled = channel(normal, true, 2e-5);
        let mut work = serial.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
        let mut parallel = pooled.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
        let molecular: Vec<_> = (0..100)
            .map(|cell| {
                let mut collision = law();
                if let CollisionModel27::Cumulant { ref mut rates } = collision.model {
                    rates.omega[0] = 1.90 + 0.005 * (cell % 5) as f64;
                    rates.omega[1] += 0.001 * (cell % 3) as f64;
                }
                collision
            })
            .collect();
        let mass: f64 = serial.populations().iter().flatten().sum();
        let pointer = work.collisions().as_ptr();
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let before = serial
                .prepare_channel_dynamic_sgs(cx, config(), Some(&molecular), &mut work)
                .unwrap();
            assert!(before.fit.maximum_eddy_viscosity > 0.0);
            for (cell, (&augmented, &base)) in work.collisions().iter().zip(&molecular).enumerate()
            {
                let (
                    CollisionModel27::Cumulant { rates: actual },
                    CollisionModel27::Cumulant { rates: expected },
                ) = (augmented.model, base.model)
                else {
                    panic!("collision family changed");
                };
                assert_eq!(actual.omega[1..], expected.omega[1..]);
                assert_eq!(augmented.maximum_mach, base.maximum_mach);
                let nu = actual.viscosity().unwrap();
                assert!(
                    (nu - expected.viscosity().unwrap() - work.eddy_viscosities()[cell]).abs()
                        < 2e-16
                );
                assert_eq!(
                    work.local.velocity[cell],
                    serial.macroscopic(cell).unwrap().1
                );
            }
            original
                .step_serial_with_collisions(cx, work.collisions())
                .unwrap();
            let accepted = serial
                .step_serial_with_channel_dynamic_sgs(cx, config(), Some(&molecular), &mut work)
                .unwrap();
            let pool = TilePool::for_host(2, 93);
            let actual = pooled
                .step_pooled_with_channel_dynamic_sgs(
                    &pool,
                    &CancelGate::new_clock_free(),
                    config(),
                    Some(&molecular),
                    &mut parallel,
                )
                .unwrap();
            assert_eq!(before, accepted);
            assert_eq!(accepted, actual);
            assert_eq!(serial.populations(), original.populations());
            assert_eq!(serial.populations(), pooled.populations());
            assert_eq!(serial.completed_steps(), 1);
            assert_eq!(pooled.completed_steps(), 1);
            assert_eq!(serial.collision, law());
            assert_eq!(pointer, work.collisions().as_ptr());
            let final_mass: f64 = serial.populations().iter().flatten().sum();
            assert!((final_mass - mass).abs() < 2e-11);
        });
    }
}

#[test]
fn g0_channel_caps_topology_workspace_and_cancellation_refuse_then_retry() {
    let mut flow = channel(2, true, 2e-5);
    let before = flow.populations().to_vec();
    let bytes = 100 * (168 + core::mem::size_of::<CollisionConfig27>());
    assert!(matches!(
        flow.channel_dynamic_sgs_workspace(bytes - 1),
        Err(SgsStepError::Dynamic(
            DynamicSgsError::WorkspaceBytes { .. }
        ))
    ));
    let mut work = flow.channel_dynamic_sgs_workspace(bytes).unwrap();
    assert_eq!(work.workspace_bytes(), bytes);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for configuration in [
            ChannelDynamicSgsConfig {
                maximum_coefficient_squared: 0.0,
                ..config()
            },
            ChannelDynamicSgsConfig {
                maximum_eddy_to_molecular: 0.0,
                ..config()
            },
        ] {
            assert!(
                flow.step_serial_with_channel_dynamic_sgs(cx, configuration, None, &mut work)
                    .is_err()
            );
            assert_eq!(flow.populations(), before);
            assert_eq!(flow.completed_steps(), 0);
        }
        // Reusing a buffer with equal shape but a different wall normal refuses.
        flow.periodic = [false, true, true];
        assert!(matches!(
            flow.prepare_channel_dynamic_sgs(cx, config(), None, &mut work),
            Err(SgsStepError::Grid(GridError3::Input(_)))
        ));
        for periodic in [[true; 3], [false, false, true]] {
            flow.periodic = periodic;
            assert!(flow.channel_dynamic_sgs_workspace(bytes).is_err());
        }
        flow.periodic = [true, true, false];
        flow.solid[0] = true;
        assert!(
            flow.prepare_channel_dynamic_sgs(cx, config(), None, &mut work)
                .is_err()
        );
        flow.solid[0] = false;
    });
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        assert!(matches!(
            flow.step_serial_with_channel_dynamic_sgs(cx, config(), None, &mut work),
            Err(SgsStepError::Grid(GridError3::Cancelled))
        ));
    });
    assert_eq!(flow.populations(), before);
    assert_eq!(flow.completed_steps(), 0);
    let mut reference = channel(2, true, 2e-5);
    let mut fresh = reference.channel_dynamic_sgs_workspace(bytes).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let a = flow
            .step_serial_with_channel_dynamic_sgs(cx, config(), None, &mut work)
            .unwrap();
        let b = reference
            .step_serial_with_channel_dynamic_sgs(cx, config(), None, &mut fresh)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(flow.populations(), reference.populations());
    });
}

#[test]
fn g1_channel_late_forced_mach_failure_keeps_populations_and_clock() {
    let mut flow = channel(2, false, 0.002);
    flow.budget.maximum_mach = 0.001;
    let before = flow.populations().to_vec();
    let mut work = flow.channel_dynamic_sgs_workspace(usize::MAX).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        assert!(matches!(
            flow.step_serial_with_channel_dynamic_sgs(cx, config(), None, &mut work),
            Err(SgsStepError::Grid(GridError3::Mach { .. }))
        ));
    });
    assert_eq!(flow.populations(), before);
    assert_eq!(flow.completed_steps(), 0);
}
