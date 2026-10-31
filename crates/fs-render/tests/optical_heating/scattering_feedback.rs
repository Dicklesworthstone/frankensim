use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::feedback::{
    IncidentBeam, OpticalFeedbackBudget, OpticalFeedbackError, TemperatureAbsorptionTable,
    TemperatureOpticalHeating, ThermalOpticalCell,
};
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::{
    IncidentScatteringBeam, ScatteringFeedbackBudget, TemperatureScatteringHeating,
    TemperatureScatteringTable, ThermalScatteringCell,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;

fn absorption(rows: &[(f64, f64)]) -> TemperatureAbsorptionTable {
    let rows: Vec<_> = rows
        .iter()
        .map(|&(t, a)| (t, MassAbsorptionSpectrum::constant(a).unwrap()))
        .collect();
    TemperatureAbsorptionTable::try_new(ContentHash([71; 32]), &rows, 16).unwrap()
}
fn scattering(rows: &[(f64, f64)]) -> TemperatureScatteringTable {
    let rows: Vec<_> = rows
        .iter()
        .map(|&(t, s)| (t, MassScatteringSpectrum::constant(s).unwrap()))
        .collect();
    TemperatureScatteringTable::try_new(ContentHash([72; 32]), &rows, 16).unwrap()
}
fn a(value: f64) -> TemperatureAbsorptionTable {
    absorption(&[(290.0, value), (800.0, value)])
}
fn s(value: f64) -> TemperatureScatteringTable {
    scattering(&[(290.0, value), (800.0, value)])
}
fn budget() -> ScatteringFeedbackBudget {
    ScatteringFeedbackBudget {
        max_grid_cells: 8,
        max_tables: 16,
        max_total_temperature_knots: 128,
        max_source_evaluations: 1000,
        max_traced_segments: 20_000_000,
    }
}
fn photon_config(histories: u32) -> ScatteringDepositionConfig {
    ScatteringDepositionConfig {
        histories_per_lane: histories,
        max_collisions_per_history: 64,
        max_total_segments: 1_000_000,
        max_deposition_cells: 8,
        seed: 91,
        source_id: 723,
    }
}
fn tolerance() -> ScatteringHeatingTolerance {
    ScatteringHeatingTolerance {
        power_roundoff_w: 1e-10,
        max_unresolved_w: 1.0,
        max_total_standard_error_w: 10.0,
        max_cell_standard_error_w: 10.0,
    }
}
fn binding(grid_cell: usize, thermal_cell: u64, absorption_table: usize) -> ThermalScatteringCell {
    ThermalScatteringCell {
        grid_cell,
        thermal_cell: CellId(thermal_cell),
        absorption_table,
        scattering_table: 0,
        phase: 0,
    }
}
fn geometry<'g>(cx: &Cx<'_>, density: &'g [f64], dims: [usize; 3]) -> AbsorbingGrid<'g> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(dims, density, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(99.0).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: 8,
            max_steps_per_segment: 64,
        },
    )
    .unwrap()
}
#[allow(clippy::too_many_arguments)]
fn observer<'g>(
    cx: &Cx<'_>,
    density: &'g [f64],
    dims: [usize; 3],
    thermal: &DepositionSimulation<'_>,
    power: f64,
    absorption: &[TemperatureAbsorptionTable],
    scattering: &[TemperatureScatteringTable],
    bindings: &[ThermalScatteringCell],
    config: ScatteringDepositionConfig,
    tolerance: ScatteringHeatingTolerance,
) -> Result<TemperatureScatteringHeating<'g, 1>, OpticalFeedbackError> {
    TemperatureScatteringHeating::try_new(
        geometry(cx, density, dims),
        IncidentScatteringBeam {
            origin: Point3::new(-1.0, 0.5, 0.5),
            unit_direction: Vec3::new(1.0, 0.0, 0.0),
            wavelengths_nm: [550.0],
            incident_power_w: [power],
        },
        absorption,
        scattering,
        &[PhaseFunction::Isotropic],
        bindings,
        thermal,
        budget(),
        config,
        tolerance,
        || Ok(()),
    )
}

#[test]
fn joint_absorption_only_is_exact_native_temperature_beer_feedback() {
    with_cx(|cx| {
        let density = [1.0, 0.5];
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(12, 2.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let tables = [absorption(&[(290.0, 0.2), (800.0, 8.0)]), a(1.5)];
        let joint = observer(
            cx,
            &density,
            [2, 1, 1],
            &original,
            12.0,
            &tables,
            &[s(0.0)],
            &[binding(0, 11, 0), binding(1, 12, 1)],
            photon_config(64),
            tolerance(),
        )
        .unwrap();
        assert!(joint.deposition().is_none());
        assert!(!joint.has_scattering());
        let directed_geometry = geometry(cx, &density, [2, 1, 1]);
        let directed = TemperatureOpticalHeating::try_new(
            cx,
            &directed_geometry,
            IncidentBeam {
                origin: Point3::new(-1.0, 0.5, 0.5),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                max_distance_m: 4.0,
                wavelengths_nm: [550.0],
                incident_power_w: [12.0],
            },
            &tables,
            &[
                ThermalOpticalCell {
                    grid_cell: 0,
                    thermal_cell: CellId(11),
                    absorption_table: 0,
                },
                ThermalOpticalCell {
                    grid_cell: 1,
                    thermal_cell: CellId(12),
                    absorption_table: 1,
                },
            ],
            &original,
            OpticalFeedbackBudget {
                max_path_cells: 8,
                max_tables: 8,
                max_total_temperature_knots: 16,
                max_source_evaluations: 1000,
            },
            1e-10,
            || Ok(()),
        )
        .unwrap();
        let (hot_joint, joint, joint_audit) = original
            .advance_observed(&joint, 1.0, config(0.125), || Ok(()))
            .unwrap();
        let (hot_directed, directed, directed_audit) = original
            .advance_observed(&directed, 1.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(hot_joint.cells(), hot_directed.cells());
        assert_eq!(joint_audit, directed_audit);
        assert_eq!(joint.energy().absorbed_j, directed.energy().absorbed_j);
        assert_eq!(
            joint.energy().ballistic_escape_j,
            directed.energy().transmitted_j
        );
        assert_eq!(
            joint.energy().thermal_input_j,
            directed.energy().thermal_input_j
        );
        assert_eq!(joint.energy().diffuse_escape_j, [0.0]);
        assert_eq!(joint.energy().unresolved_j, [0.0]);
        assert_eq!(joint.deposition().unwrap().sampled_histories, [0]);
        assert_eq!(
            joint
                .conditional_source_precision()
                .maximum_absorbed_standard_error_w,
            [0.0]
        );
        assert_eq!(joint.source_evaluations(), joint_audit.substeps);
    });
}

#[test]
fn scattered_photons_heat_real_off_beam_thermal_cells() {
    with_cx(|cx| {
        let density = [1.0; 4];
        let original = thermal(
            (11..15).map(|id| cell(id, 1.0, 1.0, 0.0)).collect(),
            &SENSIBLE,
            vec![],
        );
        let bindings: Vec<_> = (0..4).map(|i| binding(i, 11 + i as u64, 0)).collect();
        let heating = observer(
            cx,
            &density,
            [2, 2, 1],
            &original,
            20.0,
            &[a(1.0)],
            &[s(3.0)],
            &bindings,
            photon_config(2048),
            tolerance(),
        )
        .unwrap();
        let (heated, accepted, audit) = original
            .advance_observed(&heating, 0.25, config(0.25), || Ok(()))
            .unwrap();
        for id in [CellId(13), CellId(14)] {
            assert!(
                heated.temperature_k(id).unwrap() > 300.0,
                "off-beam receiver {id:?}"
            );
            assert!(accepted.sources_w().unwrap()[&id] > 0.0);
        }
        close(audit.external_heat_j, accepted.energy().thermal_input_j);
        close(
            accepted.energy().thermal_input_j,
            accepted.energy().absorbed_j[0],
        );
        let e = accepted.energy();
        close(
            e.incident_j[0],
            e.ballistic_escape_j[0]
                + e.diffuse_escape_j[0]
                + e.absorbed_j[0]
                + e.unresolved_j[0]
                + e.residual_j[0],
        );
        assert!(
            accepted
                .conditional_source_precision()
                .maximum_absorbed_standard_error_w[0]
                > 0.0
        );
        // Every nonvacuum cell can be reached after scattering, even before a
        // particular finite photon sample has visited it.
        assert!(matches!(
            observer(
                cx,
                &density,
                [2, 2, 1],
                &original,
                20.0,
                &[a(1.0)],
                &[s(3.0)],
                &bindings[..2],
                photon_config(64),
                tolerance()
            ),
            Err(OpticalFeedbackError::Heating(
                OpticalHeatingError::UnmappedCell(2)
            ))
        ));
    });
}

#[test]
fn native_heating_can_turn_on_scattering_and_reach_previously_unlit_cells() {
    with_cx(|cx| {
        let density = [1.0; 4];
        let original = thermal(
            vec![
                cell(11, 1.0, 1.0, 64.0),
                cell(12, 1.0, 1.0, 0.0),
                cell(13, 1.0, 1.0, 0.0),
                cell(14, 1.0, 1.0, 0.0),
            ],
            &SENSIBLE,
            vec![],
        );
        let bindings: Vec<_> = (0..4).map(|i| binding(i, 11 + i as u64, 0)).collect();
        let changing = scattering(&[(290.0, 0.0), (300.0, 0.0), (301.0, 4.0), (800.0, 4.0)]);
        let heating = observer(
            cx,
            &density,
            [2, 2, 1],
            &original,
            16.0,
            &[a(0.75)],
            &[changing],
            &bindings,
            photon_config(2048),
            tolerance(),
        )
        .unwrap();
        assert!(heating.has_scattering());
        let (middle, first, _) = original
            .advance_observed(&heating, 0.125, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(first.deposition().unwrap().sampled_histories, [0]);
        assert_eq!(middle.temperature_k(CellId(13)).unwrap(), 300.0);
        let (heated, second, _) = middle
            .advance_observed(&first, 0.25, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(second.deposition().unwrap().sampled_histories, [2048]);
        assert!(heated.temperature_k(CellId(13)).unwrap() > 300.0);
        assert_ne!(
            first.deposition().unwrap().medium_identity,
            second.deposition().unwrap().medium_identity
        );
        assert_ne!(
            first.deposition().unwrap().config.seed,
            second.deposition().unwrap().config.seed
        );
        assert!(second.energy().diffuse_escape_j[0] > 0.0);
    });
}

#[test]
fn temperature_dependent_absorption_changes_downstream_scattered_heat() {
    with_cx(|cx| {
        let density = [1.0; 2];
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 96.0), cell(12, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let bindings = [binding(0, 11, 0), binding(1, 12, 1)];
        let dynamic = [
            absorption(&[(290.0, 0.1), (300.0, 0.1), (303.0, 4.0), (800.0, 4.0)]),
            a(0.5),
        ];
        let changing = observer(
            cx,
            &density,
            [2, 1, 1],
            &original,
            10.0,
            &dynamic,
            &[s(1.0)],
            &bindings,
            photon_config(2048),
            tolerance(),
        )
        .unwrap();
        let frozen = observer(
            cx,
            &density,
            [2, 1, 1],
            &original,
            10.0,
            &[a(0.1), a(0.5)],
            &[s(1.0)],
            &bindings,
            photon_config(2048),
            tolerance(),
        )
        .unwrap();
        let (hot_changing, changing, _) = original
            .advance_observed(&changing, 1.0, config(0.125), || Ok(()))
            .unwrap();
        let (hot_frozen, frozen, _) = original
            .advance_observed(&frozen, 1.0, config(0.125), || Ok(()))
            .unwrap();
        assert!(
            hot_changing.cells()[&CellId(12)].specific_h_j_kg
                < hot_frozen.cells()[&CellId(12)].specific_h_j_kg * 0.5
        );
        assert!(
            hot_changing.cells()[&CellId(11)].specific_h_j_kg
                > hot_frozen.cells()[&CellId(11)].specific_h_j_kg
        );
        assert!(changing.energy().absorbed_j[0] > frozen.energy().absorbed_j[0]);
        // Paired studies share physical-time draws even though their actual
        // accepted local optical fields now differ.
        assert_eq!(
            changing.deposition().unwrap().config.seed,
            frozen.deposition().unwrap().config.seed
        );
        assert_ne!(
            changing.deposition().unwrap().medium_identity,
            frozen.deposition().unwrap().medium_identity
        );
    });
}

#[test]
fn scattering_feedback_respects_latent_storage_events_and_split_advances() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0)],
            &LATENT,
            vec![TimedEvent {
                time_s: 0.5,
                event: Event::Boundary(CellId(11), boundary(4.0)),
            }],
        );
        let heating = observer(
            cx,
            &density,
            [1, 1, 1],
            &original,
            100.0,
            &[a(1.0)],
            &[s(1.0)],
            &[binding(0, 11, 0)],
            photon_config(1024),
            tolerance(),
        )
        .unwrap();
        let (whole, whole_optics, audit) = original
            .advance_observed(&heating, 1.0, config(0.125), || Ok(()))
            .unwrap();
        let (middle, middle_optics, _) = original
            .advance_observed(&heating, 0.5, config(0.125), || Ok(()))
            .unwrap();
        let (split, split_optics, _) = middle
            .advance_observed(&middle_optics, 1.0, config(0.125), || Ok(()))
            .unwrap();
        assert_eq!(whole.temperature_k(CellId(11)).unwrap(), 310.0);
        assert!(whole.cells()[&CellId(11)].specific_h_j_kg > 10.0);
        assert_eq!(split.cells(), whole.cells());
        assert_eq!(split_optics.energy(), whole_optics.energy());
        assert_eq!(split_optics.deposition(), whole_optics.deposition());
        assert_eq!(audit.events, 1);
        close(
            audit.external_heat_j,
            whole_optics.energy().thermal_input_j + 2.0,
        );
        assert_eq!(whole_optics.source_evaluations(), 8);
        assert_eq!(split_optics.source_evaluations(), 4);
    });
}

#[test]
fn whole_advance_optical_work_and_cancellation_refuse_transactionally() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let heating = observer(
            cx,
            &density,
            [1, 1, 1],
            &original,
            8.0,
            &[a(1.0)],
            &[s(1.0)],
            &[binding(0, 11, 0)],
            photon_config(64),
            tolerance(),
        )
        .unwrap();
        let (_, first, _) = original
            .advance_observed(&heating, 0.25, config(0.25), || Ok(()))
            .unwrap();
        let limited = heating.with_work_limits(1, 20_000_000);
        assert!(matches!(
            original.advance_observed(&limited, 0.5, config(0.25), || Ok(())),
            Err(OpticalFeedbackError::Budget)
        ));
        let limited = heating.with_work_limits(100, first.traced_segments());
        assert!(matches!(
            original.advance_observed(&limited, 0.5, config(0.25), || Ok(())),
            Err(OpticalFeedbackError::Budget)
        ));
        let mut checkpoints = 0;
        let refused = original.advance_observed(&heating, 0.5, config(0.25), || {
            checkpoints += 1;
            if checkpoints == 50 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            refused,
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(checkpoints, 50);
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
        assert_eq!(heating.time_s(), 0.0);
        assert_eq!(heating.energy().thermal_input_j, 0.0);
        assert!(heating.deposition().is_none());
    });
}

#[test]
fn source_domains_and_uncertainty_allowances_gate_actual_heat() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(vec![cell(11, 1.0, 1.0, 100.0)], &SENSIBLE, vec![]);
        let heating = observer(
            cx,
            &density,
            [1, 1, 1],
            &original,
            8.0,
            &[a(1.0)],
            &[scattering(&[(300.0, 1.0), (301.0, 1.0)])],
            &[binding(0, 11, 0)],
            photon_config(64),
            tolerance(),
        )
        .unwrap();
        assert!(matches!(
            original.advance_observed(&heating, 0.1, config(0.1), || Ok(())),
            Err(OpticalFeedbackError::TemperatureDomain {
                thermal_cell: CellId(11),
                ..
            })
        ));
        for strict in [
            ScatteringHeatingTolerance {
                max_total_standard_error_w: 0.0,
                ..tolerance()
            },
            ScatteringHeatingTolerance {
                max_cell_standard_error_w: 0.0,
                ..tolerance()
            },
        ] {
            let heating = observer(
                cx,
                &density,
                [1, 1, 1],
                &original,
                8.0,
                &[a(1.0)],
                &[s(1.0)],
                &[binding(0, 11, 0)],
                photon_config(64),
                strict,
            )
            .unwrap();
            assert!(matches!(
                original.advance_observed(&heating, 0.01, config(0.01), || Ok(())),
                Err(OpticalFeedbackError::Heating(OpticalHeatingError::Invalid(
                    _
                )))
            ));
        }
        let truncated = observer(
            cx,
            &density,
            [1, 1, 1],
            &original,
            8.0,
            &[a(0.001)],
            &[s(20.0)],
            &[binding(0, 11, 0)],
            ScatteringDepositionConfig {
                max_collisions_per_history: 1,
                ..photon_config(64)
            },
            ScatteringHeatingTolerance {
                max_unresolved_w: 0.0,
                ..tolerance()
            },
        )
        .unwrap();
        assert!(matches!(
            original.advance_observed(&truncated, 0.01, config(0.01), || Ok(())),
            Err(OpticalFeedbackError::Heating(OpticalHeatingError::Invalid(
                _
            )))
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    });
}

#[test]
fn changing_history_count_preserves_source_stream_prefix_and_replicas_are_distinct() {
    with_cx(|cx| {
        let density = [1.0];
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let heating = observer(
            cx,
            &density,
            [1, 1, 1],
            &original,
            8.0,
            &[a(1.0)],
            &[s(1.0)],
            &[binding(0, 11, 0)],
            photon_config(64),
            tolerance(),
        )
        .unwrap();
        let (_, first, _) = original
            .advance_observed(&heating, 0.25, config(0.25), || Ok(()))
            .unwrap();
        let (_, more, _) = original
            .advance_observed(
                &heating.with_histories_per_lane(128).unwrap(),
                0.25,
                config(0.25),
                || Ok(()),
            )
            .unwrap();
        let (_, replica, _) = original
            .advance_observed(&heating.with_replica_seed(1), 0.25, config(0.25), || Ok(()))
            .unwrap();
        assert_eq!(
            first.deposition().unwrap().config.seed,
            more.deposition().unwrap().config.seed
        );
        assert_eq!(
            first.deposition().unwrap().config.source_id,
            more.deposition().unwrap().config.source_id
        );
        assert_eq!(more.deposition().unwrap().sampled_histories, [128]);
        assert_ne!(
            first.deposition().unwrap().config.seed,
            replica.deposition().unwrap().config.seed
        );
        assert_ne!(
            first.deposition().unwrap().absorbed,
            replica.deposition().unwrap().absorbed
        );
        assert_eq!(
            first.deposition().unwrap().medium_identity,
            replica.deposition().unwrap().medium_identity
        );
    });
}
