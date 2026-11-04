//! G0/G1/G4/G5 local closures driven by the native moving-wall boundary data.

use super::*;
use crate::d3q19::{CollisionModel3, equilibrium3};
use crate::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use crate::lattice3::{D3q19, GridBudget3};
use crate::sgs3::{SgsModel, SgsUse};
use fs_exec::TilePool;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 93,
                kernel_id: 147,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn models() -> [SgsModel; 3] {
    [
        SgsModel::Smagorinsky { cs: 0.16 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ]
}

fn config(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.2,
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

fn channel<L: SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: impl Fn(f64, [f64; 3]) -> [f64; Q],
    normal: usize,
    walls: [[f64; 3]; 2],
    varied: bool,
    force_size: f64,
) -> Grid3<L, Q> {
    let tangent = (normal + 1) % 3;
    let mut force = [0.0; 3];
    force[tangent] = force_size;
    Grid3::try_new(
        [4; 3],
        collision,
        force,
        core::array::from_fn(|axis| axis != normal),
        GridBudget3 {
            maximum_cells: 64,
            maximum_mach: 0.2,
        },
        |_| false,
        |point| {
            let fraction = (point[normal] as f64 + 0.5) / 4.0;
            let cell = (point[2] * 4 + point[1]) * 4 + point[0];
            let mut velocity = core::array::from_fn(|axis| {
                walls[0][axis] + fraction * (walls[1][axis] - walls[0][axis])
            });
            for axis in 0..3 {
                if varied {
                    velocity[axis] += 0.02
                        * fraction
                        * (1.0 - fraction)
                        * (((cell * 17 + axis * 23 + 11) % 31) as f64 / 31.0 - 0.5);
                }
                velocity[axis] -= 0.5 * force[axis];
            }
            equilibrium(1.0, velocity)
        },
    )
    .unwrap()
}

fn equilibrium(rho: f64, velocity: [f64; 3]) -> [f64; 27] {
    equilibrium27(rho, velocity, Equilibrium27::Factorized)
}

fn drive<L: Lattice3<Q>, const Q: usize>(
    flow: &Grid3<L, Q>,
    cx: &Cx<'_>,
    normal: usize,
    walls: [[f64; 3]; 2],
) -> (WallMomentumProbe3<L, Q>, Vec<[f64; 3]>) {
    let probe = flow.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
    let velocities = probe
        .wall_links()
        .map(|link| walls[usize::from(link.exterior_faces[2 * normal + 1])])
        .collect();
    (probe, velocities)
}

fn step<L: SgsLattice3<Q>, const Q: usize>(
    flow: &mut Grid3<L, Q>,
    cx: &Cx<'_>,
    cfg: SgsConfig,
    probe: &WallMomentumProbe3<L, Q>,
    velocities: &[[f64; 3]],
    work: &mut SgsWorkspace3<L::Collision>,
) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
    flow.step_serial_with_channel_sgs_and_moving_walls(cx, cfg, None, probe, velocities, work)
}

fn couette<L: SgsLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: impl Fn(f64, [f64; 3]) -> [f64; Q] + Copy,
) {
    for normal in 0..3 {
        for difference in [0.0, -0.02, 0.02] {
            let tangent = (normal + 1) % 3;
            let mut walls = [[0.0; 3]; 2];
            walls[0][tangent] = 0.01;
            walls[1][tangent] = 0.01 + difference;
            let flow = channel::<L, Q>(collision, equilibrium, normal, walls, false, 2e-5);
            let molecular: Vec<_> = (0..64)
                .map(|i| L::with_shear_relaxation(collision, 0.75 + 0.01 * (i % 3) as f64).unwrap())
                .collect();
            let shear: f64 = difference / 4.0;
            let expected_loss = molecular
                .iter()
                .map(|&law| L::molecular_viscosity(law).unwrap() * shear * shear)
                .sum::<f64>();
            let mut work = flow.sgs_workspace().unwrap();
            with_cx(&CancelGate::new_clock_free(), |cx| {
                let (probe, velocities) = drive(&flow, cx, normal, walls);
                for model in models() {
                    // For rank-one Couette shear, WALE and Vreman vanish;
                    // Smagorinsky gives (Cs Delta)^2 |du_t/dn| independently.
                    let expected_eddy = match model {
                        SgsModel::Smagorinsky { cs } => (cs * 1.2).powi(2) * shear.abs(),
                        _ => 0.0,
                    };
                    let report = flow
                        .prepare_channel_sgs_with_moving_walls(
                            cx,
                            config(model),
                            Some(&molecular),
                            &probe,
                            &velocities,
                            &mut work,
                        )
                        .unwrap();
                    for (&total, &base) in work.collisions().iter().zip(&molecular) {
                        let eddy = L::molecular_viscosity(total).unwrap()
                            - L::molecular_viscosity(base).unwrap();
                        assert!((eddy - expected_eddy).abs() < 3e-15);
                    }
                    assert!((report.maximum_eddy_viscosity_lu - expected_eddy).abs() < 3e-15);
                    assert!((report.resolved_dissipation_lu - expected_loss).abs() < 2e-16);
                    assert!(
                        (report.sgs_dissipation_lu - 64.0 * expected_eddy * shear * shear).abs()
                            < 2e-16
                    );
                    if difference == 0.0 {
                        assert!(report.resolved_dissipation_lu < 1e-26);
                    }
                }
            });
        }
    }
}

#[test]
fn g1_local_moving_channel_couette_and_comotion_for_both_lattices_and_all_normals() {
    couette::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
    couette::<D3q27, 27>(law(), equilibrium);
}

#[test]
fn g1_local_moving_channel_frame_and_native_serial_pooled_equivalence() {
    let pool = TilePool::for_host(2, 93);
    for normal in 0..3 {
        let mut boost = [0.0; 3];
        boost[(normal + 1) % 3] = 0.012;
        boost[(normal + 2) % 3] = -0.006;
        for model in models() {
            let new = || channel::<D3q27, 27>(law(), equilibrium, normal, [boost; 2], true, 2e-5);
            let mut serial = new();
            let mut native = new();
            let mut pooled = new();
            let fixed = channel::<D3q27, 27>(law(), equilibrium, normal, [[0.0; 3]; 2], true, 2e-5);
            let mut work = serial.sgs_workspace().unwrap();
            let mut parallel = pooled.sgs_workspace().unwrap();
            let mut baseline = fixed.sgs_workspace().unwrap();
            let molecular: Vec<_> = (0..64)
                .map(|i| {
                    let mut collision = law();
                    if let CollisionModel27::Cumulant { ref mut rates } = collision.model {
                        rates.omega[0] = 1.90 + 0.005 * (i % 5) as f64;
                        rates.omega[1] += 0.001 * (i % 3) as f64;
                    }
                    collision
                })
                .collect();
            let pointer = work.collisions().as_ptr();
            with_cx(&CancelGate::new_clock_free(), |cx| {
                let reference = fixed
                    .prepare_sgs(cx, config(model), Some(&molecular), &mut baseline)
                    .unwrap();
                let (probe, velocities) = drive(&serial, cx, normal, [boost; 2]);
                let prepared = serial
                    .prepare_channel_sgs_with_moving_walls(
                        cx,
                        config(model),
                        Some(&molecular),
                        &probe,
                        &velocities,
                        &mut work,
                    )
                    .unwrap();
                assert!(prepared.maximum_eddy_viscosity_lu > 0.0);
                assert!((prepared.sgs_dissipation_lu - reference.sgs_dissipation_lu).abs() < 1e-13);
                for ((&total, &base), &fixed) in work
                    .collisions()
                    .iter()
                    .zip(&molecular)
                    .zip(baseline.collisions())
                {
                    let (
                        CollisionModel27::Cumulant { rates },
                        CollisionModel27::Cumulant { rates: original },
                    ) = (total.model, base.model)
                    else {
                        panic!("collision family changed")
                    };
                    assert_eq!(rates.omega[1..], original.omega[1..]);
                    assert_eq!(total.maximum_mach, base.maximum_mach);
                    assert!(
                        (D3q27::molecular_viscosity(total).unwrap()
                            - D3q27::molecular_viscosity(fixed).unwrap())
                        .abs()
                            < 1e-13
                    );
                }
                let expected = native
                    .step_serial_with_moving_wall_momentum(
                        cx,
                        &probe,
                        &velocities,
                        Some(work.collisions()),
                    )
                    .unwrap();
                let accepted = serial
                    .step_serial_with_channel_sgs_and_moving_walls(
                        cx,
                        config(model),
                        Some(&molecular),
                        &probe,
                        &velocities,
                        &mut work,
                    )
                    .unwrap();
                let actual = pooled
                    .step_pooled_with_channel_sgs_and_moving_walls(
                        &pool,
                        &CancelGate::new_clock_free(),
                        config(model),
                        Some(&molecular),
                        &probe,
                        &velocities,
                        &mut parallel,
                    )
                    .unwrap();
                assert_eq!(accepted, (prepared, expected));
                assert_eq!(accepted, actual);
                assert_eq!(serial.populations(), native.populations());
                assert_eq!(serial.populations(), pooled.populations());
                assert_eq!(serial.completed_steps(), 1);
                assert_eq!(serial.collision, law());
                assert_eq!(pointer, work.collisions().as_ptr());
            });
        }
    }
}

#[test]
fn g0_zero_wall_drive_preserves_stationary_local_sgs_bits() {
    for model in models() {
        let new = || channel::<D3q27, 27>(law(), equilibrium, 2, [[0.0; 3]; 2], true, 2e-5);
        let mut ordinary = new();
        let mut moving = new();
        let mut base = ordinary.sgs_workspace().unwrap();
        let mut work = moving.sgs_workspace().unwrap();
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let (probe, zero) = drive(&moving, cx, 2, [[0.0; 3]; 2]);
            let expected = ordinary
                .step_serial_with_sgs(cx, config(model), None, &mut base)
                .unwrap();
            let actual = step(&mut moving, cx, config(model), &probe, &zero, &mut work).unwrap();
            assert_eq!(actual.0, expected);
            assert_eq!(moving.populations(), ordinary.populations());
            assert_eq!(actual.1.wall_work, 0.0);
        });
    }
}

#[test]
fn g4_local_moving_channel_invalid_drive_and_cancel_are_atomic_and_retryable() {
    let walls = [[0.012, 0.0, 0.0]; 2];
    let new = || channel::<D3q27, 27>(law(), equilibrium, 2, walls, true, 2e-5);
    let mut flow = new();
    let before = flow.populations().to_vec();
    let mut work = flow.sgs_workspace().unwrap();
    let cfg = config(models()[0]);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (probe, velocities) = drive(&flow, cx, 2, walls);
        let partial = flow
            .wall_momentum_probe(cx, [0.0; 3], |link| link.exterior_faces[5])
            .unwrap();
        let partial_drive = vec![[0.0; 3]; partial.measured_links()];
        assert!(step(&mut flow, cx, cfg, &partial, &partial_drive, &mut work).is_err());
        let nonuniform: Vec<_> = probe
            .wall_links()
            .map(|link| [0.001 * (1.0 + link.fluid_cell[0] as f64), 0.0, 0.0])
            .collect();
        WallMotion3 {
            probe: &probe,
            velocities: &nonuniform,
        }
        .admit(&flow, || Ok(()))
        .unwrap();
        for invalid in [
            nonuniform,
            vec![[0.0, 0.0, 0.001]; velocities.len()],
            velocities[..velocities.len() - 1].to_vec(),
        ] {
            assert!(step(&mut flow, cx, cfg, &probe, &invalid, &mut work).is_err());
            assert_eq!(flow.populations(), before);
            assert_eq!(flow.completed_steps(), 0);
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_cx(&cancelled, |cx| {
            assert!(matches!(
                step(&mut flow, cx, cfg, &probe, &velocities, &mut work),
                Err(SgsStepError::Grid(GridError3::Cancelled))
            ));
        });
        assert_eq!(flow.populations(), before);
        assert_eq!(flow.completed_steps(), 0);
        let mut reference = new();
        let mut fresh = reference.sgs_workspace().unwrap();
        let expected = step(&mut reference, cx, cfg, &probe, &velocities, &mut fresh).unwrap();
        let retried = step(&mut flow, cx, cfg, &probe, &velocities, &mut work).unwrap();
        assert_eq!(retried, expected);
        assert_eq!(flow.populations(), reference.populations());
    });
}

#[test]
fn g4_local_moving_channel_late_native_mach_refusal_preserves_state() {
    let mut flow = channel::<D3q27, 27>(law(), equilibrium, 2, [[0.0; 3]; 2], false, 0.002);
    flow.budget.maximum_mach = 0.001;
    let before = flow.populations().to_vec();
    let mut work = flow.sgs_workspace().unwrap();
    let cfg = config(models()[0]);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (probe, velocities) = drive(&flow, cx, 2, [[0.0003, 0.0, 0.0]; 2]);
        flow.prepare_channel_sgs_with_moving_walls(cx, cfg, None, &probe, &velocities, &mut work)
            .unwrap();
        assert!(matches!(
            step(&mut flow, cx, cfg, &probe, &velocities, &mut work),
            Err(SgsStepError::Grid(GridError3::Mach { .. }))
        ));
    });
    assert_eq!(flow.populations(), before);
    assert_eq!(flow.completed_steps(), 0);
    assert_eq!(flow.collision, law());
}
