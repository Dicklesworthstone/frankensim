//! Original full-Stokes histories on the changing native temperature clock.

use super::*;
use fs_render::polarization::Stokes;
use fs_render::polarization::dielectric::transport::{BeamError, PolarizationFrame};
use fs_render::volumes::heating::scattering_feedback::IncidentPolarizedScatteringBeam;
use fs_render::volumes::heating::temperature_field::{
    TemperatureOpticalField, TemperatureOpticalFieldBudget,
};
use fs_render::volumes::transport::VolumeTransportError;
use fs_render::volumes::transport::scattering::deposition::PolarizedDepositionConfig;
use fs_render::volumes::transport::scattering::polarized::PolarizedScatteringError;
use std::sync::Arc;

#[allow(clippy::too_many_arguments)]
fn compiled<'g, const N: usize>(
    cx: &Cx<'_>,
    density: &'g [f64],
    dims: [usize; 3],
    wavelengths: [f64; N],
    absorption: &[TemperatureAbsorptionTable],
    scattering: &[TemperatureScatteringTable],
    phases: &[PhaseFunction],
    bindings: &[ThermalScatteringCell],
    thermal: &DepositionSimulation<'_>,
) -> Arc<TemperatureOpticalField<'g, N>> {
    let geometry = AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(dims, density, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(99.0).unwrap(),
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

fn polarized_config(histories: u32, max_phase_proposals: u32) -> PolarizedDepositionConfig {
    PolarizedDepositionConfig {
        transport: ScatteringDepositionConfig {
            max_deposition_cells: 27,
            ..photon_config(histories)
        },
        max_phase_proposals,
    }
}

fn source<const N: usize>(
    wavelengths_nm: [f64; N],
    incident_stokes_w: [Stokes; N],
) -> IncidentPolarizedScatteringBeam<N> {
    IncidentPolarizedScatteringBeam {
        origin: Point3::new(-1.0, 0.5, 0.5),
        frame: PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0))
            .unwrap(),
        wavelengths_nm,
        incident_stokes_w,
    }
}

fn observer<'g, const N: usize>(
    field: Arc<TemperatureOpticalField<'g, N>>,
    beam: IncidentPolarizedScatteringBeam<N>,
    thermal: &DepositionSimulation<'_>,
    config: PolarizedDepositionConfig,
    estimator: AbsorptionEstimator,
) -> Result<TemperatureScatteringHeating<'g, N>, OpticalFeedbackError> {
    TemperatureScatteringHeating::from_polarized_field(
        field,
        beam,
        thermal,
        ScatteringFeedbackBudget {
            max_grid_cells: 27,
            ..budget()
        },
        config,
        estimator,
        tolerance(),
        || Ok(()),
    )
}

#[test]
fn polarized_spectral_absorption_only_preserves_exact_temperature_beer_and_split_replay() {
    with_cx(|cx| {
        let density = [1.0, 0.5];
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(12, 2.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let law = TemperatureAbsorptionTable::try_new(
            ContentHash([97; 32]),
            &[
                (
                    300.0,
                    MassAbsorptionSpectrum::table(&[(380.0, 0.2), (780.0, 0.8)]).unwrap(),
                ),
                (
                    800.0,
                    MassAbsorptionSpectrum::table(&[(380.0, 2.0), (780.0, 8.0)]).unwrap(),
                ),
            ],
            2,
        )
        .unwrap();
        let field = compiled(
            cx,
            &density,
            [2, 1, 1],
            [450.0, 650.0],
            &[law],
            &[s(0.0)],
            &[PhaseFunction::Isotropic],
            &[binding(0, 11, 0), binding(1, 12, 0)],
            &original,
        );
        let beam = source(
            [450.0, 650.0],
            [
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
            ],
        );
        for estimator in [
            AbsorptionEstimator::Collision,
            AbsorptionEstimator::Continuous,
        ] {
            let polarized = observer(
                Arc::clone(&field),
                beam,
                &original,
                polarized_config(32, 64),
                estimator,
            )
            .unwrap();
            let scalar = TemperatureScatteringHeating::from_field_with_estimator(
                Arc::clone(&field),
                IncidentScatteringBeam {
                    origin: beam.origin,
                    unit_direction: beam.frame.direction(),
                    wavelengths_nm: beam.wavelengths_nm,
                    incident_power_w: [3.0, 7.0],
                },
                &original,
                budget(),
                polarized_config(32, 64).transport,
                estimator,
                tolerance(),
                || Ok(()),
            )
            .unwrap();
            assert!(polarized.is_polarized());
            assert!(!scalar.is_polarized());
            assert!(Arc::ptr_eq(
                polarized.optical_field(),
                scalar.optical_field()
            ));
            let (hot, accepted, audit) = original
                .advance_observed(&polarized, 0.5, config(0.125), || Ok(()))
                .unwrap();
            let (reference, scalar, reference_audit) = original
                .advance_observed(&scalar, 0.5, config(0.125), || Ok(()))
                .unwrap();
            assert_eq!(hot.cells(), reference.cells());
            assert_eq!(audit, reference_audit);
            assert_eq!(accepted.energy(), scalar.energy());
            assert_eq!(accepted.deposition(), scalar.deposition());
            assert_eq!(accepted.deposition().unwrap().sampled_histories, [0; 2]);
            assert_eq!(
                accepted.deposition().unwrap().absorbed.standard_error_w,
                [0.0; 2]
            );
            let (middle, first, _) = original
                .advance_observed(&polarized, 0.25, config(0.125), || Ok(()))
                .unwrap();
            let (split, last, _) = middle
                .advance_observed(&first, 0.5, config(0.125), || Ok(()))
                .unwrap();
            assert_eq!(split.cells(), hot.cells());
            assert_eq!(last.energy(), accepted.energy());
            assert_eq!(last.deposition(), accepted.deposition());
            assert!(last.is_polarized());
        }
    });
}

#[test]
fn orthogonal_stokes_sources_heat_different_native_receivers_with_temperature_feedback() {
    with_cx(|cx| {
        let density = [1.0; 27];
        let original = thermal(
            (11..16).map(|id| cell(id, 1.0, 1.0, 0.0)).collect(),
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
            [550.0],
            &[
                absorption(&[(290.0, 20.0), (300.0, 20.0), (330.0, 80.0)]),
                a(0.0),
            ],
            &[s(0.0), s(1.0)],
            &[PhaseFunction::Isotropic, PhaseFunction::Rayleigh],
            &bindings,
            &original,
        );
        let make = |q| {
            observer(
                Arc::clone(&field),
                IncidentPolarizedScatteringBeam {
                    origin: Point3::new(1.5, 1.5, 1.5),
                    frame: PolarizationFrame::try_new(
                        Vec3::new(0.0, 0.0, 1.0),
                        Vec3::new(1.0, 0.0, 0.0),
                    )
                    .unwrap(),
                    wavelengths_nm: [550.0],
                    incident_stokes_w: [Stokes {
                        i: 8.0,
                        q,
                        u: 0.0,
                        v: 0.0,
                    }],
                },
                &original,
                polarized_config(4096, 64),
                AbsorptionEstimator::Continuous,
            )
            .unwrap()
        };
        let x_source = make(8.0);
        let y_source = make(-8.0);
        let (x_hot, x_optical, x_audit) = original
            .advance_observed(&x_source, 0.5, config(0.25), || Ok(()))
            .unwrap();
        let (y_hot, y_optical, y_audit) = original
            .advance_observed(&y_source, 0.5, config(0.25), || Ok(()))
            .unwrap();
        let x_faces = [
            x_hot.cells()[&CellId(12)].specific_h_j_kg,
            x_hot.cells()[&CellId(13)].specific_h_j_kg,
        ];
        let y_faces = [
            y_hot.cells()[&CellId(12)].specific_h_j_kg,
            y_hot.cells()[&CellId(13)].specific_h_j_kg,
        ];
        // Electric-dipole radiation is weaker along the incident electric axis.
        assert!(
            x_faces[1] > x_faces[0] + 0.03,
            "x dipole native heat {x_faces:?}"
        );
        assert!(
            y_faces[0] > y_faces[1] + 0.03,
            "y dipole native heat {y_faces:?}"
        );
        assert!((x_faces[0] - y_faces[1]).abs() < 0.1);
        assert!((x_faces[1] - y_faces[0]).abs() < 0.1);
        for (optical, audit) in [(&x_optical, &x_audit), (&y_optical, &y_audit)] {
            let energy = optical.energy();
            close(energy.thermal_input_j, audit.external_heat_j);
            close(energy.absorbed_j[0], energy.thermal_input_j);
            close(
                energy.incident_j[0],
                energy.absorbed_j[0]
                    + energy.ballistic_escape_j[0]
                    + energy.diffuse_escape_j[0]
                    + energy.unresolved_j[0]
                    + energy.residual_j[0],
            );
            assert_eq!(optical.source_evaluations(), 2);
            assert_eq!(optical.deposition().unwrap().sampled_histories, [4096]);
        }
        let (middle, first, _) = original
            .advance_observed(&x_source, 0.25, config(0.25), || Ok(()))
            .unwrap();
        let (split, last, _) = middle
            .advance_observed(&first, 0.5, config(0.25), || Ok(()))
            .unwrap();
        assert_eq!(split.cells(), x_hot.cells());
        assert_eq!(last.energy(), x_optical.energy());
        assert_eq!(last.deposition(), x_optical.deposition());
        assert_ne!(
            first.deposition().unwrap().medium_identity,
            last.deposition().unwrap().medium_identity,
        );
        assert_eq!(original.time_s(), 0.0);
        assert!(x_source.deposition().is_none());
    });
}

#[test]
fn scalar_phase_that_turns_on_after_heating_refuses_the_entire_native_advance() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(vec![cell(11, 1.0, 1.0, 64.0)], &SENSIBLE, vec![]);
        let field = compiled(
            cx,
            &density,
            [1; 3],
            [550.0],
            &[a(0.5)],
            &[scattering(&[
                (290.0, 0.0),
                (300.0, 0.0),
                (301.0, 1.0),
                (800.0, 1.0),
            ])],
            &[PhaseFunction::Isotropic],
            &[binding(0, 11, 0)],
            &original,
        );
        let optical = observer(
            field,
            source(
                [550.0],
                [Stokes {
                    i: 4.0,
                    q: 4.0,
                    u: 0.0,
                    v: 0.0,
                }],
            ),
            &original,
            polarized_config(32, 64),
            AbsorptionEstimator::Continuous,
        )
        .unwrap();
        assert!(optical.has_scattering());
        let refused = original.advance_observed(&optical, 0.125, config(0.0625), || Ok(()));
        assert!(matches!(
            refused,
            Err(OpticalFeedbackError::Polarized(
                PolarizedScatteringError::UnsupportedPhase
            ))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
        assert_eq!(optical.time_s(), 0.0);
        assert!(optical.deposition().is_none());
        assert_eq!(optical.energy().thermal_input_j, 0.0);
        let (warm, first, _) = original
            .advance_observed(&optical, 0.0625, config(0.0625), || Ok(()))
            .unwrap();
        assert!(warm.temperature_k(CellId(11)).unwrap() > 301.0);
        assert_eq!(first.deposition().unwrap().sampled_histories, [0]);
        let retained_energy = *first.energy();
        assert!(matches!(
            warm.advance_observed(&first, 0.125, config(0.0625), || Ok(())),
            Err(OpticalFeedbackError::Polarized(
                PolarizedScatteringError::UnsupportedPhase
            )),
        ));
        assert_eq!(first.energy(), &retained_energy);
        assert_eq!(first.time_s(), 0.0625);
    });
}

#[test]
fn polarized_source_admission_proposal_caps_and_native_cancellation_never_publish_heat() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let field = compiled(
            cx,
            &density,
            [1; 3],
            [550.0],
            &[a(0.2)],
            &[s(2.0)],
            &[PhaseFunction::Rayleigh],
            &[binding(0, 11, 0)],
            &original,
        );
        let beam = source(
            [550.0],
            [Stokes {
                i: 4.0,
                q: 4.0,
                u: 0.0,
                v: 0.0,
            }],
        );
        let invalid = source(
            [550.0],
            [Stokes {
                i: 4.0,
                q: 4.1,
                u: 0.0,
                v: 0.0,
            }],
        );
        assert!(matches!(
            observer(
                Arc::clone(&field),
                invalid,
                &original,
                polarized_config(64, 64),
                AbsorptionEstimator::Continuous
            ),
            Err(OpticalFeedbackError::Polarized(
                PolarizedScatteringError::Beam(BeamError::InvalidStokes)
            )),
        ));
        assert!(matches!(
            observer(
                Arc::clone(&field),
                beam,
                &original,
                polarized_config(64, 0),
                AbsorptionEstimator::Continuous
            ),
            Err(OpticalFeedbackError::Transport(
                VolumeTransportError::Budget
            )),
        ));
        let capped = observer(
            Arc::clone(&field),
            beam,
            &original,
            polarized_config(64, 1),
            AbsorptionEstimator::Continuous,
        )
        .unwrap();
        assert!(matches!(
            original.advance_observed(&capped, 0.25, config(0.25), || Ok(())),
            Err(OpticalFeedbackError::Transport(
                VolumeTransportError::Budget
            )),
        ));
        let optical = observer(
            field,
            beam,
            &original,
            polarized_config(64, 64),
            AbsorptionEstimator::Continuous,
        )
        .unwrap();
        let mut checks = 0;
        let refused = original.advance_observed(&optical, 0.25, config(0.25), || {
            checks += 1;
            if checks == 50 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            refused,
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(checks, 50);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
        assert_eq!(optical.energy().thermal_input_j, 0.0);
        assert!(optical.deposition().is_none());
        assert!(capped.deposition().is_none());
    });
}

#[test]
fn real_polarized_scattering_is_retained_by_complete_thermal_trajectory_refinement() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let field = compiled(
            cx,
            &density,
            [1; 3],
            [550.0],
            &[absorption(&[(300.0, 0.1), (400.0, 1.1)])],
            &[s(1.0)],
            &[PhaseFunction::Rayleigh],
            &[binding(0, 11, 0)],
            &original,
        );
        let optical = observer(
            field,
            source(
                [550.0],
                [Stokes {
                    i: 4.0,
                    q: 0.0,
                    u: 3.2,
                    v: 2.4,
                }],
            ),
            &original,
            polarized_config(64, 64),
            AbsorptionEstimator::Continuous,
        )
        .unwrap();
        let accuracy = ScatteringEnsembleTolerance {
            enthalpy_absolute_j_kg: 1.0,
            temperature_absolute_k: 1.0,
            optical_energy_absolute_j: 1.0,
            relative: 0.0,
        };
        let result = optical
            .advance_ensemble_refined(
                &original,
                1.0,
                ScatteringFeedbackRefinement {
                    initial_max_dt_s: 0.5,
                    minimum_max_dt_s: 0.0625,
                    maximum_halvings: 3,
                    replicas: 2,
                    ensemble_seed: 173,
                    maximum_photon_doublings: 1,
                    standard_error_multiplier: 3.0,
                    temporal: accuracy,
                    sampling: accuracy,
                    photon_sensitivity: accuracy,
                    conservation_tolerance_j: 1e-9,
                    maximum_total_substeps: 1000,
                    maximum_total_source_evaluations: 1000,
                    maximum_total_traced_segments: 1_000_000,
                    maximum_retained_trajectories: 4,
                    maximum_observable_cells: 1,
                },
                || Ok(()),
            )
            .unwrap();
        assert!(result.diagnostics.temporal_refinement_observed);
        assert!(result.diagnostics.photon_sensitivity.is_some());
        assert_eq!(result.diagnostics.histories_per_lane, 128);
        for trajectory in &result.trajectories {
            assert!(trajectory.optical.is_polarized());
            assert_eq!(
                trajectory.optical.absorption_estimator(),
                AbsorptionEstimator::Continuous
            );
            let report = trajectory.optical.deposition().unwrap();
            assert_eq!(report.sampled_histories, [128]);
            assert!(report.diffuse_escape.mean_w[0] > 0.0);
            close(
                trajectory.audit.external_heat_j,
                trajectory.optical.energy().absorbed_j[0],
            );
            close(
                trajectory.thermal.cells()[&CellId(11)].specific_h_j_kg,
                trajectory.optical.energy().thermal_input_j,
            );
        }
        assert_ne!(
            result.trajectories[0]
                .optical
                .deposition()
                .unwrap()
                .config
                .seed,
            result.trajectories[1]
                .optical
                .deposition()
                .unwrap()
                .config
                .seed,
        );
        assert_eq!(optical.time_s(), 0.0);
        assert!(optical.deposition().is_none());
    });
}
