//! G1/G3 full-Stokes angular statistics and native heat/momentum closures.

use super::super::tests::{balance, config, medium, palette, with_cx};
use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, ScatteringOptics};
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};

fn frame() -> PolarizationFrame {
    PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap()
}
fn options(histories: u32, collisions: u32) -> PolarizedDepositionConfig {
    PolarizedDepositionConfig {
        transport: config(histories, collisions),
        max_phase_proposals: 64,
    }
}
fn beam(power: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        frame(),
        550.0,
        DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        ),
        power,
    )
    .unwrap()
}

#[test]
fn g1_matched_dipole_sampling_retains_power_pdf_and_analytic_angular_moments() {
    for (state, moments) in [
        (
            Stokes {
                i: 3.0,
                q: 3.0,
                u: 0.0,
                v: 0.0,
            },
            [0.2, 0.4, 0.4],
        ),
        (
            Stokes {
                i: 3.0,
                q: -3.0,
                u: 0.0,
                v: 0.0,
            },
            [0.4, 0.2, 0.4],
        ),
        (
            Stokes {
                i: 3.0,
                q: 0.0,
                u: 0.0,
                v: 3.0,
            },
            [0.3, 0.3, 0.4],
        ),
    ] {
        let mut rng = StreamKey {
            seed: 1933,
            kernel: 91,
            tile: 7,
        }
        .stream();
        let mut sum = [0.0; 3];
        let samples = 16_384;
        for _ in 0..samples {
            let sample = PhaseFunction::Rayleigh
                .sample_polarized_power_with_checkpoint(
                    &mut || Ok(()),
                    beam(state),
                    64,
                    &mut || rng.next_f64(),
                )
                .unwrap();
            let direction = sample.beam.frame().direction();
            for (sum, v) in sum.iter_mut().zip([direction.x, direction.y, direction.z]) {
                *sum += v * v;
            }
            // Independent electric-dipole radiation, normalized on the sphere.
            let q = state.q / state.i;
            let expected = 3.0 / (8.0 * core::f64::consts::PI)
                * (1.0
                    - 0.5 * (1.0 + q) * direction.x.powi(2)
                    - 0.5 * (1.0 - q) * direction.y.powi(2));
            assert!((sample.pdf - expected).abs() < 3e-14);
            assert_eq!(sample.beam.power().i, state.i);
            assert!((sample.beam.power().degree() - 1.0).abs() < 3e-13);
        }
        for (sum, expected) in sum.into_iter().zip(moments) {
            assert!(
                (sum / f64::from(samples) - expected).abs() < 0.009,
                "dipole moment {} != {expected}",
                sum / f64::from(samples)
            );
        }
    }
}

#[test]
fn g1_dark_proposals_are_rejected_and_budget_or_cancellation_never_returns_a_packet() {
    let input = beam(Stokes {
        i: 7.0,
        q: 7.0,
        u: 0.0,
        v: 0.0,
    });
    // +x is the exact radiation null of an x-directed dipole. +y is bright.
    let mut draws = [0.5, 0.0, 0.0, 0.5, 0.25, 0.0].into_iter();
    let sample = PhaseFunction::Rayleigh
        .sample_polarized_power_with_checkpoint(&mut || Ok(()), input, 2, &mut || {
            draws.next().unwrap()
        })
        .unwrap();
    assert_eq!(draws.next(), None);
    assert_eq!(sample.beam.power().i, 7.0);
    assert!(sample.beam.frame().direction().y > 1.0 - 1e-14);
    let mut draws = [0.5, 0.0, 0.0].into_iter();
    assert_eq!(
        PhaseFunction::Rayleigh.sample_polarized_power_with_checkpoint(
            &mut || Ok(()),
            input,
            1,
            &mut || draws.next().unwrap(),
        ),
        Err(VolumeTransportError::Budget.into())
    );
    assert_eq!(
        PhaseFunction::Rayleigh.sample_polarized_power_with_checkpoint(
            &mut || Err(VolumeTransportError::Cancelled),
            input,
            64,
            &mut || panic!("cancelled"),
        ),
        Err(VolumeTransportError::Cancelled.into())
    );
    let mut draws = [0.5, 0.25, f64::NAN].into_iter();
    assert_eq!(
        PhaseFunction::Rayleigh.sample_polarized_power_with_checkpoint(
            &mut || Ok(()),
            input,
            1,
            &mut || draws.next().unwrap(),
        ),
        Err(VolumeTransportError::Sample.into())
    );
}

#[test]
fn g1_polarized_absorption_only_reuses_exact_beer_and_scalar_phase_is_not_required() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0, 2.0], 0.4, 0.0, PhaseFunction::Isotropic);
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let frame =
            PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let options = options(2, 1);
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
        let polarized = m
            .deposit_polarized_multiple_scattering_power(
                cx,
                origin,
                frame,
                &[530.0, 680.0],
                &source,
                options,
            )
            .unwrap();
        let scalar = m
            .deposit_multiple_scattering_power(
                cx,
                origin,
                frame.direction(),
                &[530.0, 680.0],
                &[3.0, 7.0],
                options.transport,
            )
            .unwrap();
        assert_eq!(polarized, scalar);
        assert_eq!(polarized.sampled_histories, [0; 2]);
        balance(&polarized);
    });
}

#[test]
fn g1_spectral_dipole_endpoint_retains_weak_power_after_frame_rotation() {
    let phase = PhaseFunction::spectral_table(&[
        (380.0, PhaseFunction::Rayleigh),
        (780.0, PhaseFunction::rayleigh_depolarized(0.08).unwrap()),
    ])
    .unwrap();
    let incident = PolarizedBeam::try_new(
        Point3::new(-1.0, 0.5, 0.5),
        PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap(),
        450.0,
        beam(Stokes::unpolarized(1.0)).medium(),
        Stokes {
            i: 1.0,
            q: 1.0,
            u: 0.0,
            v: 0.0,
        },
    )
    .unwrap();
    // Fixed original proposal that exposed roundoff amplified by the weak
    // pure-Rayleigh endpoint, before the depolarized response was mixed in.
    let draws = [0.44655692445231976, 0.7538174234823749];
    let weighted = phase.sample_polarized(incident, draws).unwrap();
    let direction = weighted.beam.frame().direction();
    let dipole = 1.0 - direction.y * direction.y;
    let pi = core::f64::consts::PI;
    let rayleigh = 3.0 * dipole / (8.0 * pi);
    let depolarized = 3.0 * (0.92 * dipole + 0.08) / (4.0 * pi * 2.08);
    let expected = 0.825 * rayleigh + 0.175 * depolarized;
    assert!(expected > 0.0);
    assert!((weighted.beam.power().i * weighted.pdf - expected).abs() < 3e-15);
    let mut draws = [draws[0], draws[1], 0.0].into_iter();
    let matched = phase
        .sample_polarized_power_with_checkpoint(&mut || Ok(()), incident, 1, &mut || {
            draws.next().unwrap()
        })
        .unwrap();
    assert_eq!(matched.beam.power().i, 1.0);
    assert!((matched.pdf - expected).abs() < 3e-15);
    assert!(matched.beam.power().degree() <= 1.0 + 16.0 * f64::EPSILON);

    // The correction preserves intensity: neither exact extinction nor the
    // positive representable parallel dipole channel becomes an error/floor.
    let parallel = Stokes {
        i: 1.0,
        q: -1.0,
        u: 0.0,
        v: 0.0,
    };
    let dark = PhaseFunction::Rayleigh
        .apply_polarized_at(450.0, 0.0, parallel)
        .unwrap();
    assert_eq!(dark, Stokes::unpolarized(0.0));
    let weak = PhaseFunction::Rayleigh
        .apply_polarized_at(450.0, 1e-10, parallel)
        .unwrap();
    let expected_weak = 3.0 / (8.0 * pi) * 1e-20;
    assert!((weak.i / expected_weak - 1.0).abs() < 4.0 * f64::EPSILON);
}

#[test]
fn g1_spectral_polarized_histories_close_heat_force_and_unresolved_radiation() {
    with_cx(|cx, _| {
        let phase = PhaseFunction::spectral_table(&[
            (380.0, PhaseFunction::Rayleigh),
            (780.0, PhaseFunction::rayleigh_depolarized(0.08).unwrap()),
        ])
        .unwrap();
        let m = medium(cx, &[1.0, 2.0, 0.5], 0.15, 0.85, phase);
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let frame =
            PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
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
        let options = options(512, 3);
        let report = m
            .deposit_polarized_multiple_scattering_momentum(
                cx,
                origin,
                frame,
                &[450.0, 680.0],
                &source,
                options,
            )
            .unwrap();
        let power = m
            .deposit_polarized_multiple_scattering_power(
                cx,
                origin,
                frame,
                &[450.0, 680.0],
                &source,
                options,
            )
            .unwrap();
        assert_eq!(power, report.power);
        balance(&power);
        for (lane, source) in source.into_iter().enumerate() {
            assert!(power.unresolved.mean_w[lane] > 0.0);
            for axis in 0..3 {
                assert!(
                    report.momentum.balance_residual_n[lane][axis].abs()
                        < 2e-12 * source.i / 299_792_458.0
                );
            }
        }
        assert_eq!(
            report,
            m.deposit_polarized_multiple_scattering_momentum(
                cx,
                origin,
                frame,
                &[450.0, 680.0],
                &source,
                options,
            )
            .unwrap()
        );
    });
}

#[test]
fn g1_orthogonal_dipoles_heat_different_original_cells() {
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
        let mut options = options(8192, 6);
        options.transport.max_deposition_cells = 27;
        let run = |q| {
            m.deposit_polarized_multiple_scattering_power(
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
                options,
            )
            .unwrap()
        };
        let along_x = run(1.0);
        let along_y = run(-1.0);
        balance(&along_x);
        balance(&along_y);
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
        let x = sides(&along_x);
        let y = sides(&along_y);
        assert!(x[1] > x[0] + 0.02, "x dipole heat {x:?}");
        assert!(y[0] > y[1] + 0.02, "y dipole heat {y:?}");
        assert!((x[0] - y[1]).abs() < 0.015 && (x[1] - y[0]).abs() < 0.015);
    });
}

#[test]
fn g0_invalid_stokes_and_scalar_scattering_refuse_before_publishing_heat() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0], 0.4, 0.6, PhaseFunction::Isotropic);
        assert_eq!(
            m.deposit_polarized_multiple_scattering_power(
                cx,
                Point3::new(0.5, 0.5, 0.5),
                frame(),
                &[550.0],
                &[Stokes::unpolarized(1.0)],
                options(2, 2),
            ),
            Err(PolarizedScatteringError::UnsupportedPhase)
        );
        let m = medium(cx, &[1.0], 0.4, 0.6, PhaseFunction::Rayleigh);
        assert_eq!(
            m.deposit_polarized_multiple_scattering_power(
                cx,
                Point3::new(0.5, 0.5, 0.5),
                frame(),
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
    });
}
