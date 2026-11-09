//! Continuous isotropic absorption with native full-Stokes scattering histories.

use super::*;

fn x_frame() -> PolarizationFrame {
    PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

fn heat<const N: usize>(
    cx: &Cx<'_>,
    medium: &ScatteringGrid<'_>,
    wavelengths: &[f64; N],
    source: &[Stokes; N],
    config: PolarizedDepositionConfig,
) -> ScatteringPowerDeposition<N> {
    medium
        .deposit_polarized_continuous_absorption_power(
            cx,
            Point3::new(-1.0, 0.5, 0.5),
            x_frame(),
            wavelengths,
            source,
            config,
        )
        .unwrap()
}

fn force<const N: usize>(
    cx: &Cx<'_>,
    medium: &ScatteringGrid<'_>,
    wavelengths: &[f64; N],
    source: &[Stokes; N],
    config: PolarizedDepositionConfig,
) -> ScatteringMomentumDeposition<N> {
    medium
        .deposit_polarized_continuous_absorption_momentum(
            cx,
            Point3::new(-1.0, 0.5, 0.5),
            x_frame(),
            wavelengths,
            source,
            config,
        )
        .unwrap()
}

fn force_balance<const N: usize>(report: &ScatteringMomentumDeposition<N>) {
    balance(&report.power);
    for lane in 0..N {
        let scale = report.power.incident_power_w[lane] / 299_792_458.0;
        for axis in 0..3 {
            assert!(report.momentum.balance_residual_n[lane][axis].abs() <= 3e-12 * scale);
            let cells: f64 = report
                .momentum
                .cells
                .iter()
                .map(|cell| cell.force.mean_n[lane][axis])
                .sum();
            assert!((cells - report.momentum.force.mean_n[lane][axis]).abs() <= 3e-12 * scale);
        }
    }
}

#[test]
fn g1_continuous_stokes_sources_keep_exact_beer_and_weak_no_collision_heat() {
    with_cx(|cx, _| {
        let source = [
            Stokes {
                i: 3.0,
                q: 3.0,
                u: 0.0,
                v: 0.0,
            },
            Stokes {
                i: 7.0,
                q: 0.0,
                u: 4.2,
                v: 5.6,
            },
        ];
        let wavelengths = [530.0, 680.0];
        let config = options(32, 8);
        for s in [0.0, 1e-12] {
            let m = medium(
                cx,
                &[1.0, 2.0],
                1e-7,
                s,
                if s == 0.0 {
                    PhaseFunction::Isotropic
                } else {
                    PhaseFunction::Rayleigh
                },
            );
            let report = force(cx, &m, &wavelengths, &source, config);
            let scalar = m
                .deposit_continuous_absorption_momentum(
                    cx,
                    Point3::new(-1.0, 0.5, 0.5),
                    x_frame().direction(),
                    &wavelengths,
                    &[3.0, 7.0],
                    config.transport,
                )
                .unwrap();
            assert_eq!(report, scalar);
            assert_eq!(report.power, heat(cx, &m, &wavelengths, &source, config));
            assert_eq!(
                report.power.absorption_estimator,
                AbsorptionEstimator::Continuous
            );
            assert_eq!(report.power.absorbed.standard_error_w, [0.0; 2]);
            assert_eq!(report.momentum.force.standard_error_n, [[0.0; 3]; 2]);
            for (lane, source) in source.iter().enumerate() {
                let expected = source.i * -(-3e-7_f64).exp_m1();
                assert!((report.power.absorbed.mean_w[lane] / expected - 1.0).abs() < 2e-13);
                assert!(
                    (report.momentum.force.mean_n[lane][0] * 299_792_458.0 / expected - 1.0).abs()
                        < 2e-13
                );
            }
            force_balance(&report);
        }
    });
}

#[test]
fn g2_polarized_continuous_and_collision_estimators_agree_in_multiple_scattering() {
    with_cx(|cx, _| {
        for (phase, state) in [
            (
                PhaseFunction::Rayleigh,
                Stokes {
                    i: 1.0,
                    q: 0.8,
                    u: 0.6,
                    v: 0.0,
                },
            ),
            (
                PhaseFunction::rayleigh_depolarized(0.08).unwrap(),
                Stokes {
                    i: 1.0,
                    q: 0.0,
                    u: 0.0,
                    v: 1.0,
                },
            ),
        ] {
            let m = medium(cx, &[1.0], 0.4, 3.0, phase);
            let continuous = force(cx, &m, &[550.0], &[state], options(8_192, 96));
            let mut independent = options(8_192, 96);
            independent.transport.seed ^= 0x5dca_71ab;
            let collision = m
                .deposit_polarized_multiple_scattering_momentum(
                    cx,
                    Point3::new(-1.0, 0.5, 0.5),
                    x_frame(),
                    &[550.0],
                    &[state],
                    independent,
                )
                .unwrap();
            assert_eq!(continuous.power.unresolved.mean_w, [0.0]);
            assert_eq!(collision.power.unresolved.mean_w, [0.0]);
            for (a, b) in [
                (continuous.power.absorbed, collision.power.absorbed),
                (
                    continuous.power.ballistic_escape,
                    collision.power.ballistic_escape,
                ),
                (
                    continuous.power.diffuse_escape,
                    collision.power.diffuse_escape,
                ),
            ] {
                let combined =
                    (a.standard_error_w[0].powi(2) + b.standard_error_w[0].powi(2)).sqrt();
                assert!((a.mean_w[0] - b.mean_w[0]).abs() <= 6.0 * combined + 2e-13);
            }
            for axis in 0..3 {
                let a = &continuous.momentum.force;
                let b = &collision.momentum.force;
                let combined = (a.standard_error_n[0][axis].powi(2)
                    + b.standard_error_n[0][axis].powi(2))
                .sqrt();
                assert!(
                    (a.mean_n[0][axis] - b.mean_n[0][axis]).abs()
                        <= 6.0 * combined + 2e-13 / 299_792_458.0
                );
            }
            force_balance(&continuous);
            force_balance(&collision);
        }
    });
}

#[test]
fn g1_continuous_orthogonal_dipoles_heat_different_original_receivers() {
    with_cx(|cx, _| {
        let density = [1.0; 27];
        let cell_map: Vec<_> = (0..27).map(|i| u32::from(i == 13)).collect();
        let absorption = AbsorbingGrid::try_new_indexed(
            cx,
            VolumeGrid::new([3; 3], &density, [0.0; 3], [1.0; 3]),
            palette(
                cx,
                vec![
                    MassAbsorptionSpectrum::constant(20.0).unwrap(),
                    MassAbsorptionSpectrum::constant(0.0).unwrap(),
                ],
                &cell_map,
            ),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 27,
                max_steps_per_segment: 32,
            },
        )
        .unwrap();
        let m = ScatteringGrid::try_new_indexed(
            cx,
            absorption,
            palette(
                cx,
                vec![
                    ScatteringOptics {
                        scattering: MassScatteringSpectrum::constant(0.0).unwrap(),
                        phase: PhaseFunction::Isotropic,
                    },
                    ScatteringOptics {
                        scattering: MassScatteringSpectrum::constant(1.0).unwrap(),
                        phase: PhaseFunction::Rayleigh,
                    },
                ],
                &cell_map,
            ),
        )
        .unwrap();
        let mut config = options(8_192, 6);
        config.transport.max_deposition_cells = 27;
        let run = |q| {
            m.deposit_polarized_continuous_absorption_power(
                cx,
                Point3::new(1.5, 1.5, 1.5),
                frame(),
                &[550.0],
                &[Stokes {
                    i: 1.0,
                    q,
                    u: 0.0,
                    v: 0.0,
                }],
                config,
            )
            .unwrap()
        };
        let x_dipole = run(1.0);
        let y_dipole = run(-1.0);
        let sides = |report: &ScatteringPowerDeposition<1>| {
            let mut xy = [0.0; 2];
            for cell in &report.cells {
                let x = cell.cell_index % 3;
                let y = (cell.cell_index / 3) % 3;
                if x != 1 {
                    xy[0] += cell.absorbed.mean_w[0];
                }
                if y != 1 {
                    xy[1] += cell.absorbed.mean_w[0];
                }
            }
            xy
        };
        let x = sides(&x_dipole);
        let y = sides(&y_dipole);
        assert!(x[1] > x[0] + 0.02, "x-directed dipole heat {x:?}");
        assert!(y[0] > y[1] + 0.02, "y-directed dipole heat {y:?}");
        assert!((x[0] - y[1]).abs() < 0.015 && (x[1] - y[0]).abs() < 0.015);
        balance(&x_dipole);
        balance(&y_dipole);
    });
}

#[test]
fn g3_spectral_polarized_continuous_histories_replay_and_separate_elastic_force_from_heat() {
    with_cx(|cx, _| {
        let phase = PhaseFunction::spectral_table(&[
            (380.0, PhaseFunction::Rayleigh),
            (780.0, PhaseFunction::rayleigh_depolarized(0.08).unwrap()),
        ])
        .unwrap();
        let source = [
            Stokes {
                i: 3.0,
                q: 3.0,
                u: 0.0,
                v: 0.0,
            },
            Stokes {
                i: 7.0,
                q: 1.4,
                u: -2.1,
                v: 4.9,
            },
        ];
        let config = options(512, 3);
        for absorption in [0.0, 0.15] {
            let m = medium(cx, &[1.0, 2.0, 0.5], absorption, 0.85, phase.clone());
            let report = force(cx, &m, &[450.0, 680.0], &source, config);
            assert_eq!(report, force(cx, &m, &[450.0, 680.0], &source, config));
            assert_eq!(report.power, heat(cx, &m, &[450.0, 680.0], &source, config));
            assert!(!report.momentum.cells.is_empty());
            for lane in 0..2 {
                assert!(report.power.unresolved.mean_w[lane] > 0.0);
                assert!(report.power.diffuse_escape.mean_w[lane] > 0.0);
            }
            if absorption == 0.0 {
                assert_eq!(report.power.absorbed.mean_w, [0.0; 2]);
                assert!(report.power.cells.is_empty());
            } else {
                assert!(report.power.absorbed.mean_w.iter().all(|&p| p > 0.0));
            }
            force_balance(&report);
        }
    });
}

#[test]
fn g4_continuous_polarized_admission_phase_limits_and_cancellation_are_atomic() {
    with_cx(|cx, gate| {
        let source = [Stokes {
            i: 1.0,
            q: 1.0,
            u: 0.0,
            v: 0.0,
        }];
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let invalid = medium(cx, &[1.0], 0.4, 0.6, PhaseFunction::Isotropic);
        assert_eq!(
            invalid.deposit_polarized_continuous_absorption_power(
                cx,
                origin,
                x_frame(),
                &[550.0],
                &source,
                options(2, 2),
            ),
            Err(PolarizedScatteringError::UnsupportedPhase)
        );
        let m = medium(cx, &[1e6], 0.1, 1.0, PhaseFunction::Rayleigh);
        let identity = m.content_identity();
        assert_eq!(
            m.deposit_polarized_continuous_absorption_power(
                cx,
                origin,
                x_frame(),
                &[550.0],
                &[Stokes {
                    i: 1.0,
                    q: 1.1,
                    u: 0.0,
                    v: 0.0
                }],
                options(2, 2),
            ),
            Err(BeamError::InvalidStokes.into())
        );
        let config = options(32, 1);
        // Momentum needs the final capped event's physical outgoing direction;
        // refusing a phase proposal must never publish its already scored heat.
        assert_eq!(
            m.deposit_polarized_continuous_absorption_momentum(
                cx,
                origin,
                x_frame(),
                &[550.0],
                &source,
                PolarizedDepositionConfig {
                    max_phase_proposals: 1,
                    ..config
                },
            ),
            Err(VolumeTransportError::Budget.into())
        );
        assert_eq!(m.content_identity(), identity);
        gate.request();
        assert_eq!(
            m.deposit_polarized_continuous_absorption_momentum(
                cx,
                origin,
                x_frame(),
                &[550.0],
                &source,
                config,
            ),
            Err(VolumeTransportError::Cancelled.into())
        );
    });
}
