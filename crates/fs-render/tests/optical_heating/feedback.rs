use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::feedback::{
    FeedbackRefinement, IncidentBeam, OpticalFeedbackBudget, OpticalFeedbackError,
    TemperatureAbsorptionTable, TemperatureOpticalHeating, ThermalOpticalCell,
};

fn table(
    lo_k: f64,
    hi_k: f64,
    lo_absorption: f64,
    hi_absorption: f64,
) -> TemperatureAbsorptionTable {
    TemperatureAbsorptionTable::try_new(
        ContentHash([31; 32]),
        &[
            (
                lo_k,
                MassAbsorptionSpectrum::constant(lo_absorption).unwrap(),
            ),
            (
                hi_k,
                MassAbsorptionSpectrum::constant(hi_absorption).unwrap(),
            ),
        ],
        2,
    )
    .unwrap()
}

fn budget() -> OpticalFeedbackBudget {
    OpticalFeedbackBudget {
        max_path_cells: 8,
        max_tables: 8,
        max_total_temperature_knots: 16,
        max_source_evaluations: 1000,
    }
}

fn beam(power: f64, count: usize) -> IncidentBeam<1> {
    IncidentBeam {
        origin: Point3::new(-1.0, 0.5, 0.5),
        unit_direction: Vec3::new(1.0, 0.0, 0.0),
        max_distance_m: count as f64 + 2.0,
        wavelengths_nm: [550.0],
        incident_power_w: [power],
    }
}

fn feedback(
    cx: &Cx<'_>,
    original: &DepositionSimulation<'_>,
    power: f64,
    tables: &[TemperatureAbsorptionTable],
    bindings: &[ThermalOpticalCell],
) -> TemperatureOpticalHeating<1> {
    let density = vec![1.0; bindings.len()];
    let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(1.0).unwrap());
    TemperatureOpticalHeating::try_new(
        cx,
        &medium,
        beam(power, density.len()),
        tables,
        bindings,
        original,
        budget(),
        1.0e-11,
        || Ok(()),
    )
    .unwrap()
}

fn binding(grid_cell: usize, thermal_cell: u64, absorption_table: usize) -> ThermalOpticalCell {
    ThermalOpticalCell {
        grid_cell,
        thermal_cell: CellId(thermal_cell),
        absorption_table,
    }
}

fn refinement() -> FeedbackRefinement {
    FeedbackRefinement {
        initial_max_dt_s: 0.5,
        minimum_max_dt_s: 0.001,
        maximum_halvings: 8,
        maximum_total_substeps: 4000,
        enthalpy_absolute_tolerance_j_kg: 0.004,
        temperature_absolute_tolerance_k: 0.001,
        optical_energy_absolute_tolerance_j: 0.004,
        relative_tolerance: 0.0,
        conservation_tolerance_j: 1.0e-9,
    }
}

#[test]
fn temperature_dependent_beer_heating_matches_an_independent_nonlinear_solution() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let heating = feedback(
            cx,
            &original,
            4.0,
            &[table(300.0, 400.0, 0.2, 10.2)],
            &[binding(0, 11, 0)],
        );
        // dT/dt = P/(m cp) * (1 - exp[-a - b(T-T0)]). With
        // y=exp[a+b(T-T0)], y'=b P/(m cp) (y-1), solved independently.
        let exact_t = 300.0 + ((1.0 + 0.2_f64.exp_m1() * 0.5_f64.exp()).ln() - 0.2) / 0.1;
        let mut errors = Vec::new();
        for dt in [0.5, 0.25, 0.125] {
            let (next, optics, audit) = original
                .advance_observed(&heating, 5.0, config(dt), || Ok(()))
                .unwrap();
            errors.push((next.temperature_k(CellId(11)).unwrap() - exact_t).abs());
            close(optics.energy().thermal_input_j, audit.external_heat_j);
            close(
                optics.energy().incident_j[0],
                optics.energy().absorbed_j[0] + optics.energy().transmitted_j[0],
            );
            assert!(optics.power().absorbed_w[0] > heating.power().absorbed_w[0]);
        }
        assert!((1.9..2.1).contains(&(errors[0] / errors[1])));
        assert!((1.9..2.1).contains(&(errors[1] / errors[2])));
        let (next, optics, report) = heating
            .advance_refined(&original, 5.0, refinement(), || Ok(()))
            .unwrap();
        assert!((next.temperature_k(CellId(11)).unwrap() - exact_t).abs() < 0.0011);
        assert!(report.normalized_difference <= 1.0);
        assert!((0.9..1.1).contains(&report.observed_order.unwrap()));
        assert!(report.total_substeps > report.thermal.substeps);
        assert_eq!(report.total_substeps, report.total_source_evaluations);
        close(
            optics.energy().thermal_input_j,
            report.thermal.external_heat_j,
        );
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
        assert_eq!(heating.energy().thermal_input_j, 0.0);
    });
}

#[test]
fn constant_law_refinement_admits_roundoff_without_inventing_convergence_order() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let heating = feedback(
            cx,
            &original,
            3.7,
            &[table(300.0, 400.0, 0.43, 0.43)],
            &[binding(0, 11, 0)],
        );
        let mut controls = refinement();
        controls.initial_max_dt_s = 0.23;
        controls.maximum_halvings = 2;
        controls.enthalpy_absolute_tolerance_j_kg = 1.0e-10;
        controls.temperature_absolute_tolerance_k = 1.0e-10;
        controls.optical_energy_absolute_tolerance_j = 1.0e-10;
        let (next, optics, report) = heating
            .advance_refined(&original, 1.37, controls, || Ok(()))
            .unwrap();
        let expected_heat = 3.7 * -(-0.43_f64).exp_m1() * 1.37;
        close(optics.energy().thermal_input_j, expected_heat);
        close(
            next.temperature_k(CellId(11)).unwrap(),
            300.0 + expected_heat / 4.0,
        );
        assert_eq!(report.halvings, 2);
        assert!(report.roundoff_limited);
        assert_eq!(report.observed_order, None);
        assert!(report.normalized_difference <= 1.0);
    });
}

#[test]
fn local_temperature_changes_upstream_absorption_and_downstream_shielding() {
    with_cx(|cx| {
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(22, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let heating = feedback(
            cx,
            &original,
            8.0,
            &[
                table(300.0, 400.0, 0.1, 20.1),
                table(300.0, 400.0, 1.0, 1.0),
            ],
            &[binding(0, 11, 0), binding(1, 22, 1)],
        );
        let initial_front = heating.sources_w()[&CellId(11)];
        let initial_back = heating.sources_w()[&CellId(22)];
        let (next, optics, audit) = original
            .advance_observed(&heating, 3.0, config(0.01), || Ok(()))
            .unwrap();
        assert!(optics.sources_w()[&CellId(11)] > initial_front * 1.5);
        assert!(optics.sources_w()[&CellId(22)] < initial_back);
        assert!(next.temperature_k(CellId(22)).unwrap() > next.temperature_k(CellId(11)).unwrap());
        close(next.energy_j(), optics.energy().thermal_input_j);
        close(audit.external_heat_j, optics.energy().absorbed_j[0]);
    });
}

#[test]
fn source_refresh_obeys_native_events_latent_storage_and_split_advances() {
    with_cx(|cx| {
        let scheduled = vec![TimedEvent {
            time_s: 0.5,
            event: Event::Boundary(CellId(11), boundary(4.0)),
        }];
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &LATENT, scheduled);
        let heating = feedback(
            cx,
            &original,
            40.0,
            &[table(300.0, 400.0, 0.5, 1.5)],
            &[binding(0, 11, 0)],
        );
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
        close(
            audit.external_heat_j,
            whole_optics.energy().thermal_input_j + 2.0,
        );
        assert_eq!(audit.events, 1);
    });
}

#[test]
fn temperature_law_domain_cancellation_and_whole_refinement_budgets_roll_back() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let narrow = feedback(
            cx,
            &original,
            40.0,
            &[table(300.0, 300.01, 1.0, 1.1)],
            &[binding(0, 11, 0)],
        );
        assert!(matches!(
            original.advance_observed(&narrow, 1.0, config(1.0), || Ok(())),
            Err(OpticalFeedbackError::TemperatureDomain { .. })
        ));
        assert_eq!(narrow.energy().thermal_input_j, 0.0);
        assert_eq!(original.time_s(), 0.0);
        let heating = feedback(
            cx,
            &original,
            4.0,
            &[table(300.0, 400.0, 0.2, 10.2)],
            &[binding(0, 11, 0)],
        );
        let mut calls = 0;
        let cancelled = heating.advance_refined(&original, 5.0, refinement(), || {
            calls += 1;
            if calls > 30 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            cancelled,
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        let mut limited = refinement();
        limited.maximum_total_substeps = 31; // coarse 10 + fine 20 leaves no third trajectory.
        assert!(matches!(
            heating.advance_refined(&original, 5.0, limited, || Ok(())),
            Err(OpticalFeedbackError::Thermal(DepositionError::Budget(_)))
        ));
        let mut strict = refinement();
        strict.maximum_halvings = 2;
        strict.enthalpy_absolute_tolerance_j_kg = 1.0e-12;
        strict.temperature_absolute_tolerance_k = 1.0e-12;
        strict.optical_energy_absolute_tolerance_j = 1.0e-12;
        assert!(matches!(
            heating.advance_refined(&original, 5.0, strict, || Ok(())),
            Err(OpticalFeedbackError::Convergence { .. })
        ));
        assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
        assert_eq!(heating.energy().thermal_input_j, 0.0);
    });
}

#[test]
fn transparent_cells_require_bindings_and_sourced_domains_are_not_extrapolated() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 4.0)], &SENSIBLE, vec![]);
        let source = ContentHash([49; 32]);
        let low =
            MassAbsorptionSpectrum::source_table(&[(900.0, 0.0), (1200.0, 0.0)], source).unwrap();
        let high =
            MassAbsorptionSpectrum::source_table(&[(900.0, 0.5), (1200.0, 1.0)], source).unwrap();
        let law =
            TemperatureAbsorptionTable::try_new(source, &[(300.0, low), (400.0, high)], 2).unwrap();
        let density = [1.0];
        let medium = grid(cx, &density, low);
        let mut infrared = beam(4.0, 1);
        infrared.wavelengths_nm = [1064.0];
        assert!(matches!(
            TemperatureOpticalHeating::try_new(
                cx,
                &medium,
                infrared,
                &[law.clone()],
                &[],
                &original,
                budget(),
                1.0e-10,
                || Ok(())
            ),
            Err(OpticalFeedbackError::Heating(
                OpticalHeatingError::UnmappedCell(0)
            ))
        ));
        let heating = TemperatureOpticalHeating::try_new(
            cx,
            &medium,
            infrared,
            &[law.clone()],
            &[binding(0, 11, 0)],
            &original,
            budget(),
            1.0e-10,
            || Ok(()),
        )
        .unwrap();
        assert_eq!(heating.power().thermal_input_w, 0.0);
        let (_, warmed, _) = original
            .advance_observed(&heating, 1.0, config(0.1), || Ok(()))
            .unwrap();
        assert!(warmed.energy().thermal_input_j > 0.0);
        assert_eq!(warmed.power().wavelength_bounds_nm, (900.0, 1200.0));
        infrared.wavelengths_nm = [550.0];
        assert!(matches!(
            TemperatureOpticalHeating::try_new(
                cx,
                &medium,
                infrared,
                &[law],
                &[binding(0, 11, 0)],
                &original,
                budget(),
                1.0e-10,
                || Ok(())
            ),
            Err(OpticalFeedbackError::Transport(_))
        ));
    });
}
