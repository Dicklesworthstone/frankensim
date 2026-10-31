use super::*;
use fs_render::volumes::heating::scattering::ScatteringHeatingTolerance;
use fs_render::volumes::heating::scattering_feedback::{
    IncidentScatteringBeam, ScatteringFeedbackBudget, TemperatureScatteringHeating,
};
use fs_render::volumes::transport::scattering::deposition::ScatteringDepositionConfig;
use fs_thermal::enthalpy::ThermalLink;
use fs_thermal::observers::CoupledObservers;

#[test]
fn one_material_field_combines_incident_absorption_reabsorption_and_conduction() {
    with_cx(|cx| {
        let original = thermal(
            vec![warm_cell(11, 500.0), warm_cell(12, 350.0)],
            &RADIATIVE_SENSIBLE,
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(11),
                    b: CellId(12),
                    area_m2: 0.2,
                    conductance_w_k: 0.3,
                }),
            }],
        );
        let field = field(
            cx,
            &[1.0; 2],
            &original,
            &[absorption(&[(100.0, 0.2), (1000.0, 1.2)])],
            &[scattering(&[(100.0, 0.8), (1000.0, 0.8)])],
            &[11, 12],
        );
        let incident = TemperatureScatteringHeating::from_field(
            Arc::clone(&field),
            IncidentScatteringBeam {
                origin: Point3::new(-1.0, 0.5, 0.5),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                wavelengths_nm: [10_000.0],
                incident_power_w: [80.0],
            },
            &original,
            ScatteringFeedbackBudget {
                max_grid_cells: 2,
                max_tables: 2,
                max_total_temperature_knots: 4,
                max_source_evaluations: 1000,
                max_traced_segments: 2_000_000,
            },
            ScatteringDepositionConfig {
                histories_per_lane: 256,
                max_collisions_per_history: 96,
                max_total_segments: 100_000,
                max_deposition_cells: 2,
                seed: 48,
                source_id: 719,
            },
            ScatteringHeatingTolerance {
                power_roundoff_w: 1e-9,
                max_unresolved_w: 1e-8,
                max_total_standard_error_w: 20.0,
                max_cell_standard_error_w: 20.0,
            },
            || Ok(()),
        )
        .unwrap();
        let radiation = observer(field, &original, 20.0, blackbody(300.0));
        assert!(Arc::ptr_eq(incident.optical_field(), radiation.field()));
        let pair = CoupledObservers::new(incident, radiation);

        // Reconstruct the first joint enthalpy increment independently from
        // both original source owners plus the authored conductive link.
        // Inspection is a private clone; the actual advance samples its own
        // identical source state on the shared native timeline.
        let mut inspected = pair.clone();
        let incident_w = inspected.first.sources(&original, &mut || Ok(())).unwrap();
        let bath_w = inspected.second.sources(&original, &mut || Ok(())).unwrap();
        let transfers = inspected
            .second
            .heat_transfers(&original, &mut || Ok(()))
            .unwrap();
        let first_dt = 1.0 / 128.0;
        let (first, accepted, audit) = original
            .advance_observed(&pair, first_dt, config(1.0), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 1);
        for (id, other) in [(CellId(11), CellId(12)), (CellId(12), CellId(11))] {
            let mut watts = incident_w.get(&id).copied().unwrap_or(0.0)
                + bath_w.get(&id).copied().unwrap_or(0.0)
                + 0.3
                    * (original.temperature_k(other).unwrap()
                        - original.temperature_k(id).unwrap());
            for transfer in &transfers {
                if transfer.a == id {
                    watts -= transfer.power_w;
                }
                if transfer.b == id {
                    watts += transfer.power_w;
                }
            }
            let before = &original.cells()[&id];
            let expected = before.specific_h_j_kg + first_dt * watts / before.mass_kg;
            assert!((first.cells()[&id].specific_h_j_kg - expected).abs() < 2e-11);
        }
        assert!(
            (audit.external_heat_j
                - accepted.first.energy().thermal_input_j
                - accepted.second.energy().thermal_external_j)
                .abs()
                < 2e-10
        );

        let (final_state, final_pair, _) = first
            .advance_observed(&accepted, 0.125, config(0.01), || Ok(()))
            .unwrap();
        assert!(Arc::ptr_eq(
            final_pair.first.optical_field(),
            final_pair.second.field()
        ));
        assert_eq!(final_pair.first.time_s(), final_state.time_s());
        assert_eq!(final_pair.second.time_s(), final_state.time_s());
        assert!(final_pair.second.operator_builds() > 1);
        assert!(final_pair.first.energy().thermal_input_j > 0.0);
        assert!(final_pair.second.energy().thermal_external_j < 0.0);
        assert!(final_pair.second.energy().internal_transfer_magnitude_j[0] > 0.0);
        let stored = final_state.energy_j() - original.energy_j();
        let external = final_pair.first.energy().thermal_input_j
            + final_pair.second.energy().thermal_external_j;
        assert!((stored - external).abs() < 2e-9);
    });
}
