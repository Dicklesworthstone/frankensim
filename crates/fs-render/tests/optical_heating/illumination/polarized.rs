//! Native full-Stokes source sets, finite pulse fluence and conserved heat.

use super::*;
use fs_render::polarization::Stokes;
use fs_render::polarization::dielectric::transport::PolarizationFrame;
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_render::volumes::heating::illumination::{
    IncidentPolarizedIlluminationSource, ScheduledIlluminationSource,
    ScheduledPolarizedIlluminationSource,
};
use fs_render::volumes::heating::scattering_feedback::IncidentPolarizedScatteringBeam;
use fs_render::volumes::transport::VolumeTransportError;
use fs_render::volumes::transport::scattering::deposition::{
    AbsorptionEstimator, PolarizedDepositionConfig,
};
use fs_render::volumes::transport::scattering::polarized::PolarizedScatteringError;

fn a(rows: &[(f64, f64)]) -> TemperatureAbsorptionTable {
    TemperatureAbsorptionTable::try_new(
        ContentHash([91; 32]),
        &rows
            .iter()
            .map(|&(t, value)| (t, MassAbsorptionSpectrum::constant(value).unwrap()))
            .collect::<Vec<_>>(),
        16,
    )
    .unwrap()
}
fn ac(value: f64) -> TemperatureAbsorptionTable {
    a(&[(290.0, value), (800.0, value)])
}
fn s(value: f64) -> TemperatureScatteringTable {
    TemperatureScatteringTable::try_new(
        ContentHash([92; 32]),
        &[
            (290.0, MassScatteringSpectrum::constant(value).unwrap()),
            (800.0, MassScatteringSpectrum::constant(value).unwrap()),
        ],
        2,
    )
    .unwrap()
}

#[allow(clippy::too_many_arguments)]
fn compiled<'g, const N: usize>(
    cx: &Cx<'_>,
    density: &'g [f64],
    dims: [usize; 3],
    thermal: &DepositionSimulation<'_>,
    wavelengths: [f64; N],
    absorption: &[TemperatureAbsorptionTable],
    scattering: &[TemperatureScatteringTable],
    phases: &[PhaseFunction],
    bindings: &[ThermalScatteringCell],
) -> Arc<TemperatureOpticalField<'g, N>> {
    let geometry = AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(dims, density, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(0.0).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: 27,
            max_steps_per_segment: 64,
        },
    )
    .unwrap();
    Arc::new(
        TemperatureOpticalField::try_new(
            geometry,
            wavelengths,
            absorption,
            scattering,
            phases,
            bindings,
            thermal,
            TemperatureOpticalFieldBudget {
                max_grid_cells: 27,
                max_tables: 16,
                max_total_temperature_knots: 128,
            },
            || Ok(()),
        )
        .unwrap(),
    )
}
fn bindings(cells: usize) -> Vec<ThermalScatteringCell> {
    (0..cells)
        .map(|grid_cell| ThermalScatteringCell {
            grid_cell,
            thermal_cell: CellId(11 + grid_cell as u64),
            absorption_table: 0,
            scattering_table: 0,
            phase: 0,
        })
        .collect()
}
fn frame(right: bool) -> PolarizationFrame {
    PolarizationFrame::try_new(
        Vec3::new(if right { -1.0 } else { 1.0 }, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    )
    .unwrap()
}
fn pol_config(histories: u32, max_phase_proposals: u32) -> PolarizedDepositionConfig {
    PolarizedDepositionConfig {
        transport: ScatteringDepositionConfig {
            histories_per_lane: histories,
            max_deposition_cells: 27,
            ..photons()
        },
        max_phase_proposals,
    }
}
fn profile<const N: usize>(
    segments: &[(f64, f64, [f64; N], [f64; N])],
) -> SpectralPowerSchedule<N> {
    SpectralPowerSchedule::try_new(
        &segments
            .iter()
            .map(
                |&(start_s, end_s, start_power_w, end_power_w)| SpectralPowerSegment {
                    start_s,
                    end_s,
                    start_power_w,
                    end_power_w,
                },
            )
            .collect::<Vec<_>>(),
        16,
    )
    .unwrap()
}
fn scheduled<'g, const N: usize>(
    field: Arc<TemperatureOpticalField<'g, N>>,
    sources: &[ScheduledPolarizedIlluminationSource<N>],
    thermal: &DepositionSimulation<'_>,
    histories: u32,
) -> TemperatureIlluminationHeating<'g, N> {
    TemperatureIlluminationHeating::from_scheduled_polarized_field(
        field,
        sources,
        thermal,
        limits(),
        pol_config(histories, 64),
        AbsorptionEstimator::Continuous,
        tolerance(),
        32,
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn short_spectral_polarized_pulses_and_ramps_deliver_exact_native_beer_fluence() {
    with_cx(|cx| {
        let density = [1.0, 2.0];
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(12, 2.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let law = MassAbsorptionSpectrum::table(&[
            (380.0, 0.2),
            (450.0, 0.2),
            (650.0, 0.7),
            (780.0, 0.7),
        ])
        .unwrap();
        let field = compiled(
            cx,
            &density,
            [2, 1, 1],
            &original,
            [450.0, 650.0],
            &[TemperatureAbsorptionTable::try_new(
                ContentHash([93; 32]),
                &[(290.0, law.clone()), (800.0, law)],
                2,
            )
            .unwrap()],
            &[s(0.0)],
            &[PhaseFunction::Isotropic],
            &bindings(2),
        );
        let left = ScheduledPolarizedIlluminationSource {
            source_id: 20,
            origin: Point3::new(-1.0, 0.5, 0.5),
            frame: frame(false),
            polarization: [
                Stokes {
                    i: 1.0,
                    q: 1.0,
                    u: 0.0,
                    v: 0.0,
                },
                Stokes {
                    i: 1.0,
                    q: 0.0,
                    u: 0.0,
                    v: 1.0,
                },
            ],
            power: profile(&[
                (0.203, 0.207, [8.0, 0.0], [8.0, 0.0]),
                (0.4, 0.45, [0.0, 0.0], [12.0, 6.0]),
                (0.45, 0.5, [12.0, 6.0], [0.0, 0.0]),
            ]),
        };
        let right = ScheduledPolarizedIlluminationSource {
            source_id: 10,
            origin: Point3::new(3.0, 0.5, 0.5),
            frame: frame(true),
            polarization: [
                Stokes {
                    i: 1.0,
                    q: 0.0,
                    u: 1.0,
                    v: 0.0,
                },
                Stokes {
                    i: 1.0,
                    q: 0.0,
                    u: 0.0,
                    v: -1.0,
                },
            ],
            power: profile(&[(0.1, 0.3, [2.0, 10.0], [6.0, 2.0])]),
        };
        let sources = [left, right];
        let optical = scheduled(Arc::clone(&field), &sources, &original, 16);
        let scalar_sources: Vec<_> = sources
            .iter()
            .map(|source| ScheduledIlluminationSource {
                source_id: source.source_id,
                origin: source.origin,
                unit_direction: source.frame.direction(),
                power: source.power.clone(),
            })
            .collect();
        let scalar = TemperatureIlluminationHeating::from_scheduled_field_with_estimator(
            field,
            &scalar_sources,
            &original,
            limits(),
            pol_config(16, 64).transport,
            AbsorptionEstimator::Continuous,
            tolerance(),
            32,
            || Ok(()),
        )
        .unwrap();
        let (hot, accepted, audit) = original
            .advance_observed(&optical, 0.8, config(1.0), || Ok(()))
            .unwrap();
        let (reference, scalar, reference_audit) = original
            .advance_observed(&scalar, 0.8, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(hot.cells(), reference.cells());
        assert_eq!(accepted.energy(), scalar.energy());
        assert_eq!(audit, reference_audit);
        let left_j = [
            (0.207 - 0.203) * 8.0 + 0.5 * (0.5 - 0.4) * 12.0,
            0.5 * (0.5 - 0.4) * 6.0,
        ];
        let right_j = [(0.3 - 0.1) * 4.0, (0.3 - 0.1) * 6.0];
        let mut expected = [0.0; 2];
        for (lane, coefficient) in [0.2_f64, 0.7].into_iter().enumerate() {
            let t0 = (-coefficient).exp();
            let t1 = (-2.0 * coefficient).exp();
            let q0 = (left_j[lane] + right_j[lane] * t1) * (1.0 - t0);
            let q1 = (right_j[lane] + left_j[lane] * t0) * (1.0 - t1);
            expected[0] += q0;
            expected[1] += q1;
            close(
                accepted.energy().incident_j[lane],
                left_j[lane] + right_j[lane],
            );
            close(accepted.energy().absorbed_j[lane], q0 + q1);
            close(
                accepted.energy().ballistic_escape_j[lane],
                (left_j[lane] + right_j[lane]) * t0 * t1,
            );
        }
        close(hot.cells()[&CellId(11)].specific_h_j_kg, expected[0]);
        close(hot.cells()[&CellId(12)].specific_h_j_kg, expected[1] / 2.0);
        close(audit.external_heat_j, expected.iter().sum());
        assert!(
            audit.substeps > 1,
            "the native step must reach the short pulse edges"
        );
        assert_eq!(accepted.transport_builds(), 1);
        for id in [10, 20] {
            assert_eq!(accepted.deposition(id).unwrap().sampled_histories, [0; 2]);
            assert_eq!(
                accepted.deposition(id).unwrap().absorption_estimator,
                AbsorptionEstimator::Continuous
            );
        }
    });
}

#[test]
fn alternating_orthogonal_pulses_heat_distinct_native_receivers_and_reuse_photons() {
    with_cx(|cx| {
        let density = [1.0; 27];
        let masses = [1.0, 6.0, 6.0, 12.0, 2.0];
        let original = thermal(
            (11..16)
                .zip(masses)
                .map(|(id, mass)| cell(id, mass, mass, 0.0))
                .collect(),
            &SENSIBLE,
            vec![],
        );
        let bindings: Vec<_> = (0..27)
            .map(|grid_cell| {
                let x = grid_cell % 3;
                let y = (grid_cell / 3) % 3;
                let center = grid_cell == 13;
                let recipient = if center {
                    11
                } else if x != 1 && y == 1 {
                    12
                } else if y != 1 && x == 1 {
                    13
                } else if x != 1 && y != 1 {
                    14
                } else {
                    15
                };
                ThermalScatteringCell {
                    grid_cell,
                    thermal_cell: CellId(recipient),
                    absorption_table: usize::from(center),
                    scattering_table: usize::from(center),
                    phase: usize::from(center),
                }
            })
            .collect();
        let field = compiled(
            cx,
            &density,
            [3; 3],
            &original,
            [550.0],
            &[ac(20.0), ac(0.0)],
            &[s(0.0), s(1.0)],
            &[PhaseFunction::Isotropic, PhaseFunction::Rayleigh],
            &bindings,
        );
        let make = |id, q, start, end| ScheduledPolarizedIlluminationSource {
            source_id: id,
            origin: Point3::new(1.5, 1.5, 1.5),
            frame: PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0))
                .unwrap(),
            polarization: [Stokes {
                i: 1.0,
                q,
                u: 0.0,
                v: 0.0,
            }],
            power: profile(&[(start, end, [4.0], [4.0])]),
        };
        let x = make(10, 1.0, 0.1, 0.2);
        let y = make(20, -1.0, 0.4, 0.5);
        let optical = scheduled(Arc::clone(&field), &[y.clone(), x.clone()], &original, 4096);
        let reordered = scheduled(field, &[x, y], &original, 4096);
        assert_eq!(optical.model_identity(), reordered.model_identity());
        assert_eq!(
            optical
                .polarized_sources()
                .unwrap()
                .iter()
                .map(|s| s.source_id)
                .collect::<Vec<_>>(),
            vec![10, 20]
        );
        let energy = |state: &DepositionSimulation<'_>, id| {
            let cell = &state.cells()[&CellId(id)];
            cell.mass_kg * cell.specific_h_j_kg
        };
        let (after_x, first, _) = original
            .advance_observed(&optical, 0.2, config(1.0), || Ok(()))
            .unwrap();
        assert!(energy(&after_x, 13) > energy(&after_x, 12) + 0.003);
        assert_eq!(first.transport_builds(), 1);
        let (dark, dark_source, _) = after_x
            .advance_observed(
                &first.with_work_limits(1000, 0),
                0.4,
                config(1.0),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(dark.cells(), after_x.cells());
        assert_eq!(dark_source.traced_segments(), 0);
        let (after_y, second, _) = dark
            .advance_observed(&dark_source, 0.5, config(1.0), || Ok(()))
            .unwrap();
        let second_x = energy(&after_y, 12) - energy(&after_x, 12);
        let second_y = energy(&after_y, 13) - energy(&after_x, 13);
        assert!(
            second_x > second_y + 0.003,
            "second dipole heat x={second_x}, y={second_y}"
        );
        assert_eq!(second.transport_builds(), 0);
        assert_eq!(second.traced_segments(), 0);
        let (resumed, last, _) = after_y
            .advance_observed(&second, 0.6, config(1.0), || Ok(()))
            .unwrap();
        let (whole, replay, audit) = original
            .advance_observed(&reordered, 0.6, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(resumed.cells(), whole.cells());
        assert_eq!(last.energy(), replay.energy());
        for id in [10, 20] {
            assert_eq!(last.deposition(id), replay.deposition(id));
            assert_eq!(replay.deposition(id).unwrap().sampled_histories, [4096]);
        }
        close(replay.energy().incident_j[0], 0.8);
        let e = replay.energy();
        close(
            e.incident_j[0],
            e.absorbed_j[0]
                + e.ballistic_escape_j[0]
                + e.diffuse_escape_j[0]
                + e.unresolved_j[0]
                + e.residual_j[0],
        );
        close(audit.external_heat_j, e.thermal_input_j);
        assert!(
            first
                .with_histories_per_lane(128)
                .unwrap()
                .deposition(10)
                .is_none()
        );
    });
}

#[test]
fn fixed_polarized_sources_preserve_native_histories_and_explicit_canonical_identity() {
    with_cx(|cx| {
        let density = [1.0; 2];
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(12, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let field = compiled(
            cx,
            &density,
            [2, 1, 1],
            &original,
            [550.0],
            &[ac(0.3)],
            &[s(1.0)],
            &[PhaseFunction::Rayleigh],
            &bindings(2),
        );
        let sources = [
            IncidentPolarizedIlluminationSource {
                source_id: 20,
                beam: IncidentPolarizedScatteringBeam {
                    origin: Point3::new(-1.0, 0.5, 0.5),
                    frame: frame(false),
                    wavelengths_nm: [550.0],
                    incident_stokes_w: [Stokes {
                        i: 2.0,
                        q: 2.0,
                        u: 0.0,
                        v: 0.0,
                    }],
                },
            },
            IncidentPolarizedIlluminationSource {
                source_id: 10,
                beam: IncidentPolarizedScatteringBeam {
                    origin: Point3::new(3.0, 0.5, 0.5),
                    frame: frame(true),
                    wavelengths_nm: [550.0],
                    incident_stokes_w: [Stokes {
                        i: 3.0,
                        q: 0.0,
                        u: 1.8,
                        v: 2.4,
                    }],
                },
            },
        ];
        let build = |sources: &[IncidentPolarizedIlluminationSource<1>], estimator| {
            TemperatureIlluminationHeating::from_polarized_field(
                Arc::clone(&field),
                sources,
                &original,
                limits(),
                pol_config(128, 64),
                estimator,
                tolerance(),
                || Ok(()),
            )
            .unwrap()
        };
        let mut identities = Vec::new();
        let coefficients = field.initial_state(&original, &mut || Ok(())).unwrap();
        let prepared = field.field(&coefficients, &mut || Ok(())).unwrap();
        for estimator in [
            AbsorptionEstimator::Collision,
            AbsorptionEstimator::Continuous,
        ] {
            let optical = build(&sources, estimator);
            let reordered = build(&[sources[1], sources[0]], estimator);
            assert_eq!(optical.model_identity(), reordered.model_identity());
            identities.push(optical.model_identity());
            let (_, accepted, _) = original
                .advance_observed(&optical, 0.1, config(0.1), || Ok(()))
                .unwrap();
            let (_, replay, _) = original
                .advance_observed(&reordered, 0.1, config(0.1), || Ok(()))
                .unwrap();
            assert_eq!(accepted.energy(), replay.energy());
            for source in &sources {
                let report = accepted.deposition(source.source_id).unwrap();
                let expected = prepared
                    .deposit_polarized_power(
                        cx,
                        source.beam.origin,
                        source.beam.frame,
                        &source.beam.incident_stokes_w,
                        PolarizedDepositionConfig {
                            transport: report.config,
                            max_phase_proposals: 64,
                        },
                        estimator,
                    )
                    .unwrap();
                assert_eq!(report, &expected);
                assert_eq!(replay.deposition(source.source_id), Some(report));
            }
            let trial = optical
                .with_replica_seed(5)
                .with_histories_per_lane(256)
                .unwrap();
            assert!(trial.is_polarized());
            assert_eq!(trial.absorption_estimator(), estimator);
            assert_eq!(
                trial
                    .polarized_sampling_config()
                    .unwrap()
                    .transport
                    .histories_per_lane,
                256
            );
        }
        assert_ne!(identities[0], identities[1]);
        let mut changed = sources;
        changed[0].beam.incident_stokes_w[0].q = -2.0;
        assert_ne!(
            build(&changed, AbsorptionEstimator::Continuous).model_identity(),
            identities[1]
        );
        changed = sources;
        changed[0].beam.frame =
            PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap();
        assert_ne!(
            build(&changed, AbsorptionEstimator::Continuous).model_identity(),
            identities[1]
        );
        let dark = |q| ScheduledPolarizedIlluminationSource {
            source_id: 37,
            origin: Point3::new(-1.0, 0.5, 0.5),
            frame: frame(false),
            polarization: [Stokes {
                i: 1.0,
                q,
                u: 0.0,
                v: 0.0,
            }],
            power: profile(&[]),
        };
        assert_ne!(
            scheduled(Arc::clone(&field), &[dark(0.0)], &original, 16).model_identity(),
            scheduled(field, &[dark(1.0)], &original, 16).model_identity(),
        );
    });
}

#[test]
fn pulsed_continuous_absorption_rebuilds_changing_temperatures_and_converges_in_fluence() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let field = compiled(
            cx,
            &density,
            [1; 3],
            &original,
            [550.0],
            &[a(&[(300.0, 0.2), (320.0, 2.2)])],
            &[s(1e-12)],
            &[PhaseFunction::Rayleigh],
            &bindings(1),
        );
        let sources = [
            ScheduledPolarizedIlluminationSource {
                source_id: 1,
                origin: Point3::new(-1.0, 0.5, 0.5),
                frame: frame(false),
                polarization: [Stokes {
                    i: 1.0,
                    q: 1.0,
                    u: 0.0,
                    v: 0.0,
                }],
                power: profile(&[(0.0, 1.0, [0.0], [16.0])]),
            },
            ScheduledPolarizedIlluminationSource {
                source_id: 2,
                origin: Point3::new(-1.0, 0.5, 0.5),
                frame: frame(false),
                polarization: [Stokes {
                    i: 1.0,
                    q: -1.0,
                    u: 0.0,
                    v: 0.0,
                }],
                power: profile(&[(0.25, 0.5, [8.0], [8.0])]),
            },
        ];
        let optical = scheduled(field, &sources, &original, 32);
        // Independent separable Beer solution in total incident fluence F=10 J:
        // exp[a(T)] - 1 = (exp[a0] - 1) exp[b F / (m cp)].
        let exact = 300.0 + ((1.0 + 0.2_f64.exp_m1() * 0.25_f64.exp()).ln() - 0.2) / 0.1;
        let mut errors = Vec::new();
        for dt in [0.25, 0.125, 0.0625] {
            let (hot, accepted, audit) = original
                .advance_observed(&optical, 1.0, config(dt), || Ok(()))
                .unwrap();
            errors.push(exact - hot.temperature_k(CellId(11)).unwrap());
            close(accepted.energy().incident_j[0], 10.0);
            close(accepted.energy().absorbed_j[0], audit.external_heat_j);
            close(
                hot.cells()[&CellId(11)].specific_h_j_kg,
                accepted.energy().thermal_input_j,
            );
            assert_eq!(accepted.transport_builds(), audit.substeps);
            for id in [1, 2] {
                let report = accepted.deposition(id).unwrap();
                assert_eq!(report.sampled_histories, [32]);
                assert_eq!(report.absorbed.standard_error_w, [0.0]);
                assert_eq!(report.diffuse_escape.mean_w, [0.0]);
            }
        }
        assert!(errors.iter().all(|error| *error > 0.0));
        assert!(
            errors[1] < 0.55 * errors[0] && errors[2] < 0.55 * errors[1],
            "{errors:?}"
        );
        assert_eq!(original.time_s(), 0.0);
        assert!(optical.deposition(1).is_none());
    });
}

#[test]
fn polarized_program_admission_caps_and_late_cancellation_are_transactional() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let field = compiled(
            cx,
            &density,
            [1; 3],
            &original,
            [550.0],
            &[ac(0.2)],
            &[s(2.0)],
            &[PhaseFunction::Rayleigh],
            &bindings(1),
        );
        let make = |id, start, end| ScheduledPolarizedIlluminationSource {
            source_id: id,
            origin: Point3::new(-1.0, 0.5, 0.5),
            frame: frame(false),
            polarization: [Stokes {
                i: 1.0,
                q: 1.0,
                u: 0.0,
                v: 0.0,
            }],
            power: profile(&[(start, end, [4.0], [4.0])]),
        };
        let sources = [make(10, 0.2, 0.3), make(20, 0.4, 0.5)];
        let create = |sources: &[ScheduledPolarizedIlluminationSource<1>], proposals, segments| {
            TemperatureIlluminationHeating::from_scheduled_polarized_field(
                Arc::clone(&field),
                sources,
                &original,
                limits(),
                pol_config(64, proposals),
                AbsorptionEstimator::Continuous,
                tolerance(),
                segments,
                || Ok(()),
            )
        };
        assert!(matches!(
            create(&sources, 64, 1),
            Err(OpticalFeedbackError::Budget)
        ));
        let mut invalid = make(30, 0.0, 0.1);
        invalid.power = profile(&[]);
        invalid.polarization[0].q = 1.1;
        assert!(matches!(
            create(&[invalid], 64, 16),
            Err(OpticalFeedbackError::Invalid(_))
        ));
        let unsupported = compiled(
            cx,
            &density,
            [1; 3],
            &original,
            [550.0],
            &[ac(0.2)],
            &[s(2.0)],
            &[PhaseFunction::Isotropic],
            &bindings(1),
        );
        assert!(matches!(
            TemperatureIlluminationHeating::from_scheduled_polarized_field(
                unsupported,
                &sources,
                &original,
                limits(),
                pol_config(64, 64),
                AbsorptionEstimator::Continuous,
                tolerance(),
                16,
                || Ok(()),
            ),
            Err(OpticalFeedbackError::Polarized(
                PolarizedScatteringError::UnsupportedPhase
            )),
        ));
        let capped = create(&sources, 1, 16).unwrap();
        assert!(matches!(
            original.advance_observed(&capped, 0.6, config(0.1), || Ok(())),
            Err(OpticalFeedbackError::Transport(
                VolumeTransportError::Budget
            )),
        ));
        let optical = create(&sources, 64, 16).unwrap();
        let mut through_first = 0;
        let (_, first, _) = original
            .advance_observed(&optical, 0.3, config(0.1), || {
                through_first += 1;
                Ok(())
            })
            .unwrap();
        assert!(first.energy().thermal_input_j > 0.0);
        let mut checks = 0;
        let refused = original.advance_observed(&optical, 0.6, config(0.1), || {
            checks += 1;
            if checks == through_first + 20 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            refused,
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(checks, through_first + 20);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
        assert_eq!(optical.energy().thermal_input_j, 0.0);
        assert!(optical.deposition(10).is_none());
        assert_eq!(capped.energy().thermal_input_j, 0.0);
    });
}
