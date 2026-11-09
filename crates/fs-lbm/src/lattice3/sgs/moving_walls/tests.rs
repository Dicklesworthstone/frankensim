//! G0/G1/G4/G5: actual nonuniform wall drives, interior solids and native steps.

use super::*;
use crate::d3q19::{CollisionModel3, equilibrium3};
use crate::d3q27::{
    CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, Rates27, equilibrium27,
};
use crate::lattice3::{D3q19, GridBudget3, ViscousLattice3};
use crate::sgs3::{SgsModel, SgsUse};
use fs_exec::TilePool;

const DIMS: [usize; 3] = [12, 7, 4];
const CELLS: usize = 336;
const FORCE: [f64; 3] = [2e-5, 0.0, -3e-5];

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 107,
                kernel_id: 153,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn config(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.2,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn models() -> [SgsModel; 3] {
    [
        SgsModel::Smagorinsky { cs: 0.16 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ]
}

fn law27() -> CollisionConfig27 {
    CollisionConfig27 {
        model: CollisionModel27::Cumulant {
            rates: Rates27 {
                omega: [1.9, 1.13, 1.21, 1.34, 1.47, 1.58, 1.66, 1.72, 1.39, 1.51],
            },
        },
        maximum_mach: 0.2,
    }
}

fn equilibrium(rho: f64, u: [f64; 3]) -> [f64; 27] {
    equilibrium27(rho, u, Equilibrium27::Factorized)
}

// A stationary solid sheet divides two fluid regions. Both sheet faces slide
// tangentially with different x/z wavelengths; the outer domain walls rest.
// Cross derivatives make the full gradient non-rank-one for WALE and Vreman.
fn physical_velocity(p: [f64; 3]) -> [f64; 3] {
    let h = if p[1] < 3.5 { p[1] } else { 7.0 - p[1] };
    [
        0.001 * (core::f64::consts::TAU * p[2] / 4.0).sin() * h,
        0.0,
        0.0015 * (core::f64::consts::TAU * p[0] / 12.0).sin() * h,
    ]
}

fn grid<L: SgsLattice3<Q>, const Q: usize>(
    law: L::Collision,
    equilibrium: impl Fn(f64, [f64; 3]) -> [f64; Q],
) -> Grid3<L, Q> {
    Grid3::try_new(
        DIMS,
        law,
        FORCE,
        [true, false, true],
        GridBudget3 {
            maximum_cells: CELLS,
            maximum_mach: 0.2,
        },
        |p| p[1] == 3,
        |p| {
            let u = physical_velocity(p.map(|v| v as f64 + 0.5));
            equilibrium(1.0, core::array::from_fn(|a| u[a] - 0.5 * FORCE[a]))
        },
    )
    .unwrap()
}

fn drive<L: Lattice3<Q>, const Q: usize>(
    grid: &Grid3<L, Q>,
    cx: &Cx<'_>,
) -> (WallMomentumProbe3<L, Q>, Vec<[f64; 3]>) {
    let probe = grid
        .wall_momentum_probe(cx, [0.0; 3], |link| link.solid_cell.is_some())
        .unwrap();
    let velocities = probe
        .wall_links()
        .map(|link| physical_velocity(link.position))
        .collect();
    (probe, velocities)
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-14, "{actual} != {expected}");
}

fn wall_sgs<L: SgsLattice3<Q>, const Q: usize>(
    law: L::Collision,
    equilibrium: impl Fn(f64, [f64; 3]) -> [f64; Q] + Copy,
) {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 107);
    let molecular: Vec<_> = (0..CELLS)
        .map(|cell| L::with_shear_relaxation(law, 0.7 + 0.03 * (cell % 4) as f64).unwrap())
        .collect();
    with_cx(&gate, |cx| {
        for model in models() {
            let mut serial = grid::<L, Q>(law, equilibrium);
            let mut pooled = grid::<L, Q>(law, equilibrium);
            let mut native = grid::<L, Q>(law, equilibrium);
            let (probe, velocities) = drive(&serial, cx);
            let mut work = serial.sgs_workspace().unwrap();
            let mut parallel = pooled.sgs_workspace().unwrap();
            let pointer = work.collisions().as_ptr();
            // The old homogeneous-plane entry correctly refuses this geometry.
            assert!(
                serial
                    .prepare_channel_sgs_with_moving_walls(
                        cx,
                        config(model),
                        Some(&molecular),
                        &probe,
                        &velocities,
                        &mut work,
                    )
                    .is_err()
            );
            let prepared = serial
                .prepare_sgs_with_moving_walls(
                    cx,
                    config(model),
                    Some(&molecular),
                    &probe,
                    &velocities,
                    &mut work,
                )
                .unwrap();
            let mut resolved = 0.0;
            let mut modeled = 0.0;
            let mut maximum = 0.0_f64;
            let kx = core::f64::consts::TAU / 12.0;
            let kz = core::f64::consts::TAU / 4.0;
            for cell in 0..CELLS {
                let p = coordinates(cell, DIMS);
                if p[1] == 3 {
                    continue;
                }
                let [x, y, z] = p.map(|v| v as f64 + 0.5);
                let (h, sign) = if p[1] < 3 { (y, 1.0) } else { (7.0 - y, -1.0) };
                // Exact discrete centered derivatives in periodic coordinates;
                // the normal affine profile is exact through the real wall.
                let gradient = [
                    [
                        0.0,
                        0.001 * (kz * z).sin() * sign,
                        0.001 * kz.sin() * (kz * z).cos() * h,
                    ],
                    [0.0; 3],
                    [
                        0.0015 * kx.sin() * (kx * x).cos() * h,
                        0.0015 * (kx * x).sin() * sign,
                        0.0,
                    ],
                ];
                let eddy = eddy_viscosity(model, 1.2, gradient).unwrap();
                let base = L::molecular_viscosity(molecular[cell]).unwrap();
                close(
                    L::molecular_viscosity(work.collisions()[cell]).unwrap(),
                    base + eddy,
                );
                let rate2 = gradient[0][1].powi(2)
                    + gradient[2][1].powi(2)
                    + (gradient[0][2] + gradient[2][0]).powi(2);
                resolved += base * rate2;
                modeled += eddy * rate2;
                maximum = maximum.max(eddy);
            }
            assert_eq!(prepared.fluid_cells, 288);
            close(prepared.resolved_dissipation_lu, resolved);
            close(prepared.sgs_dissipation_lu, modeled);
            close(prepared.maximum_eddy_viscosity_lu, maximum);
            assert!(maximum > 1e-8);

            let expected_walls = native
                .step_serial_with_moving_wall_momentum(
                    cx,
                    &probe,
                    &velocities,
                    Some(work.collisions()),
                )
                .unwrap();
            let accepted = serial
                .step_serial_with_sgs_and_moving_walls(
                    cx,
                    config(model),
                    Some(&molecular),
                    &probe,
                    &velocities,
                    &mut work,
                )
                .unwrap();
            let parallel_step = pooled
                .step_pooled_with_sgs_and_moving_walls(
                    &pool,
                    &gate,
                    config(model),
                    Some(&molecular),
                    &probe,
                    &velocities,
                    &mut parallel,
                )
                .unwrap();
            assert_eq!(accepted, (prepared, expected_walls));
            assert_eq!(accepted, parallel_step);
            assert_eq!(serial.populations(), native.populations());
            assert_eq!(serial.populations(), pooled.populations());
            assert_eq!(serial.completed_steps(), 1);
            assert_eq!(pointer, work.collisions().as_ptr());
        }
    });
}

#[test]
fn g1_g5_spatial_wall_sgs_uses_exact_surface_gradient_and_original_step_on_both_lattices() {
    wall_sgs::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
    wall_sgs::<D3q27, 27>(law27(), equilibrium);
}

#[test]
fn g0_zero_drives_and_disabled_model_preserve_original_native_paths() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for model in models() {
            let fixed = grid::<D3q27, 27>(law27(), equilibrium);
            let mut moving = grid::<D3q27, 27>(law27(), equilibrium);
            let (probe, velocities) = drive(&fixed, cx);
            let zero = vec![[0.0; 3]; velocities.len()];
            let mut base = fixed.sgs_workspace().unwrap();
            let mut work = moving.sgs_workspace().unwrap();
            let expected = fixed
                .prepare_sgs(cx, config(model), None, &mut base)
                .unwrap();
            let actual = moving
                .prepare_sgs_with_moving_walls(cx, config(model), None, &probe, &zero, &mut work)
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(work.collisions(), base.collisions());
            let corrected = moving
                .prepare_sgs_with_moving_walls(
                    cx,
                    config(model),
                    None,
                    &probe,
                    &velocities,
                    &mut work,
                )
                .unwrap();
            assert!(
                (corrected.resolved_dissipation_lu - expected.resolved_dissipation_lu).abs() > 1e-6
            );

            let mut native = grid::<D3q27, 27>(law27(), equilibrium);
            let molecular: Vec<_> = (0..CELLS)
                .map(|cell| {
                    D3q27::with_shear_relaxation(law27(), 0.71 + 0.01 * (cell % 3) as f64).unwrap()
                })
                .collect();
            let disabled = match model {
                SgsModel::Smagorinsky { .. } => SgsModel::Smagorinsky { cs: 0.0 },
                SgsModel::Wale { .. } => SgsModel::Wale { cw: 0.0 },
                SgsModel::Vreman { .. } => SgsModel::Vreman { cv: 0.0 },
            };
            let expected = native
                .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, Some(&molecular))
                .unwrap();
            let actual = moving
                .step_serial_with_sgs_and_moving_walls(
                    cx,
                    config(disabled),
                    Some(&molecular),
                    &probe,
                    &velocities,
                    &mut work,
                )
                .unwrap();
            assert_eq!(actual.0.maximum_eddy_viscosity_lu, 0.0);
            assert_eq!(actual.1, expected);
            assert_eq!(moving.populations(), native.populations());
            assert_eq!(work.collisions(), molecular);
        }
    });
}

#[test]
fn g4_invalid_motion_closure_and_cancellation_preserve_published_state() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut grid = grid::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
        let (probe, velocities) = drive(&grid, cx);
        let before = grid.populations().to_vec();
        let mut work = grid.sgs_workspace().unwrap();
        for invalid in [
            vec![[0.0; 3]; velocities.len() - 1],
            vec![[0.0, 0.01, 0.0]; velocities.len()],
            vec![[f64::NAN, 0.0, 0.0]; velocities.len()],
        ] {
            assert!(
                grid.step_serial_with_sgs_and_moving_walls(
                    cx,
                    config(models()[0]),
                    None,
                    &probe,
                    &invalid,
                    &mut work,
                )
                .is_err()
            );
        }
        let mut limited = config(models()[0]);
        limited.maximum_eddy_to_molecular = 0.0;
        assert!(
            grid.step_serial_with_sgs_and_moving_walls(
                cx,
                limited,
                None,
                &probe,
                &velocities,
                &mut work,
            )
            .is_err()
        );
        gate.request();
        assert!(
            grid.step_serial_with_sgs_and_moving_walls(
                cx,
                config(models()[0]),
                None,
                &probe,
                &velocities,
                &mut work,
            )
            .is_err()
        );
        assert_eq!(grid.populations(), before);
        assert_eq!(grid.completed_steps(), 0);
    });
}

fn forced_grid(cx: &Cx<'_>) -> ForcedGrid3<D3q27, 27> {
    let forces: Vec<_> = (0..CELLS)
        .map(|cell| {
            let p = coordinates(cell, DIMS);
            if p[1] == 3 {
                [0.0; 3]
            } else {
                [1e-5 * (p[0] + 1) as f64, 0.0, -2e-5 * (p[2] + 1) as f64]
            }
        })
        .collect();
    let grid = Grid3::try_new(
        DIMS,
        law27(),
        [0.0; 3],
        [true, false, true],
        GridBudget3 {
            maximum_cells: CELLS,
            maximum_mach: 0.2,
        },
        |p| p[1] == 3,
        |p| {
            let u = physical_velocity(p.map(|v| v as f64 + 0.5));
            let cell = (p[2] * DIMS[1] + p[1]) * DIMS[0] + p[0];
            equilibrium(1.0, core::array::from_fn(|a| u[a] - 0.5 * forces[cell][a]))
        },
    )
    .unwrap();
    ForcedGrid3::try_new(grid, cx, forces).unwrap()
}

#[test]
fn g3_g4_g5_spatial_force_epoch_and_wall_loads_share_one_native_sgs_transaction() {
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(3, 107);
    with_cx(&gate, |cx| {
        for model in models() {
            for advance in [false, true] {
                let mut serial = forced_grid(cx);
                let mut pooled = forced_grid(cx);
                let mut native = forced_grid(cx);
                let reference = grid::<D3q27, 27>(law27(), equilibrium);
                let (probe, velocities) = drive(&reference, cx);
                let mut work = serial.sgs_workspace().unwrap();
                let mut parallel = pooled.sgs_workspace().unwrap();
                let mut baseline = reference.sgs_workspace().unwrap();
                let mut wrong_frame = reference.sgs_workspace().unwrap();
                let expected = reference
                    .prepare_sgs_with_moving_walls(
                        cx,
                        config(model),
                        None,
                        &probe,
                        &velocities,
                        &mut baseline,
                    )
                    .unwrap();
                let prepared = serial
                    .prepare_sgs_with_moving_walls(
                        cx,
                        config(model),
                        None,
                        &probe,
                        &velocities,
                        &mut work,
                    )
                    .unwrap();
                close(
                    prepared.maximum_eddy_viscosity_lu,
                    expected.maximum_eddy_viscosity_lu,
                );
                close(
                    prepared.resolved_dissipation_lu,
                    expected.resolved_dissipation_lu,
                );
                for (&actual, &wanted) in work.collisions().iter().zip(baseline.collisions()) {
                    close(
                        D3q27::molecular_viscosity(actual).unwrap(),
                        D3q27::molecular_viscosity(wanted).unwrap(),
                    );
                }
                // The raw population velocity has a spatially varying half-
                // force offset. Failing to include the accepted F_n creates
                // a false strain and a measurably different turbulent field.
                serial
                    .grid
                    .prepare_sgs_with_moving_walls(
                        cx,
                        config(model),
                        None,
                        &probe,
                        &velocities,
                        &mut wrong_frame,
                    )
                    .unwrap();
                let error: f64 = work
                    .collisions()
                    .iter()
                    .zip(wrong_frame.collisions())
                    .map(|(&a, &b)| {
                        (D3q27::molecular_viscosity(a).unwrap()
                            - D3q27::molecular_viscosity(b).unwrap())
                        .abs()
                    })
                    .sum();
                assert!(error > 1e-9);

                let before = serial.populations().to_vec();
                let forces = serial.forces().to_vec();
                let mut invalid_endpoint = forces.clone();
                invalid_endpoint[0][0] = f64::NAN;
                assert!(
                    serial
                        .advance_serial_with_sgs_and_moving_walls(
                            cx,
                            config(model),
                            None,
                            &probe,
                            &velocities,
                            &mut work,
                            invalid_endpoint,
                        )
                        .is_err()
                );

                assert_eq!(serial.populations(), before);
                assert_eq!(serial.forces(), forces);
                assert_eq!(serial.completed_steps(), 0);
                let endpoint: Vec<_> = forces
                    .iter()
                    .map(|f| [1.5 * f[0], f[1], -0.5 * f[2]])
                    .collect();
                let expected_walls = if advance {
                    native
                        .advance_serial_with_moving_wall_momentum(
                            cx,
                            &probe,
                            &velocities,
                            work.collisions(),
                            endpoint.clone(),
                        )
                        .unwrap()
                } else {
                    native
                        .step_serial_with_moving_wall_momentum(
                            cx,
                            &probe,
                            &velocities,
                            work.collisions(),
                        )
                        .unwrap()
                };
                let accepted = if advance {
                    serial
                        .advance_serial_with_sgs_and_moving_walls(
                            cx,
                            config(model),
                            None,
                            &probe,
                            &velocities,
                            &mut work,
                            endpoint.clone(),
                        )
                        .unwrap()
                } else {
                    serial
                        .step_serial_with_sgs_and_moving_walls(
                            cx,
                            config(model),
                            None,
                            &probe,
                            &velocities,
                            &mut work,
                        )
                        .unwrap()
                };
                let parallel_step = if advance {
                    pooled
                        .advance_pooled_with_sgs_and_moving_walls(
                            &pool,
                            &gate,
                            config(model),
                            None,
                            &probe,
                            &velocities,
                            &mut parallel,
                            endpoint.clone(),
                        )
                        .unwrap()
                } else {
                    pooled
                        .step_pooled_with_sgs_and_moving_walls(
                            &pool,
                            &gate,
                            config(model),
                            None,
                            &probe,
                            &velocities,
                            &mut parallel,
                        )
                        .unwrap()
                };
                assert_eq!(accepted, (prepared, expected_walls));
                assert_eq!(parallel_step, accepted);
                assert_eq!(serial.populations(), native.populations());
                assert_eq!(serial.populations(), pooled.populations());
                let accepted_forces = if advance { endpoint } else { forces };
                assert_eq!(serial.forces(), accepted_forces);
                assert_eq!(pooled.forces(), accepted_forces);
                assert_eq!(serial.completed_steps(), 1);
                assert_eq!(pooled.completed_steps(), 1);
            }
        }
    });
}
