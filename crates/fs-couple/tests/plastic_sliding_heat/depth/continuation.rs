//! One accepted history through adaptive heating, sourced resistance and export.

use super::*;
use fs_couple::friction_heat::plastic_wear::adaptive::{
    PlasticSlidingMarchControls, PlasticSlidingMarchTolerance,
};
use fs_couple::friction_heat::plastic_wear::depth::DepthSlidingHeatReport;
use fs_couple::friction_heat::plastic_wear::depth::adaptive::DepthSlidingMarchControls;
use fs_thermal::enthalpy::removal::MaterialRemovalTolerances;

#[test]
fn g3_adaptive_depth_resistance_survives_export_and_actual_reversal() {
    with_cx(|cx, _| {
        let data = depth_fixture_with_extra_claim(
            320.0,
            cx,
            Some(resistance::resistance_claim(resistance::R0, true)),
        );
        let thermal = network_with_resistance(&data.caloric, T0, resistance::R0);
        let contact = resistance::bound(&data, &thermal, resistance::R0, 0.3, cx);
        let original = contact
            .contact_resistance()
            .unwrap()
            .original_source()
            .clone();
        let mut velocity = vec![[0.0; 2]; 25];
        velocity[12][0] = 0.02;
        let march = DepthSlidingMarchControls {
            march: PlasticSlidingMarchControls {
                initial_step_s: 0.02,
                minimum_step_s: 0.000625,
                maximum_step_s: 0.02,
                maximum_attempts: 128,
                maximum_accepted_substeps: 256,
                tolerance: PlasticSlidingMarchTolerance {
                    cell_energy_j: 1e-8,
                    interface_temperature_k: 0.02,
                    pressure_pa: 500.0,
                    length_m: 5e-11,
                    work_j: 1e-8,
                    relative: 0.01,
                },
                energy_tolerance_j: 3e-10,
            },
            body_temperature_k: 0.02,
            hardness_pa: 500.0,
        };
        let check_source = |report: &DepthSlidingHeatReport<[f64; 3]>, distances: [f64; 2]| {
            let sample = report.contact_resistance.as_ref().unwrap();
            let pressure = report
                .mechanics
                .candidate()
                .wear_state()
                .contact_state()
                .pressure_pa()[12];
            let temperature = report.heat.interface_temperature_k;
            assert_eq!(sample.mean_loaded_pressure_pa, pressure);
            assert_eq!(
                sample.source.state_coordinate(resistance::PRESSURE),
                Some(pressure)
            );
            assert_eq!(sample.source.state_coordinate("T"), Some(temperature));
            assert_eq!(sample.source.card_identity(), original.card_identity());
            assert_eq!(
                sample.source.properties()[0].answer().receipt.selected,
                original.properties()[0].answer().receipt.selected
            );
            let expected = resistance::R0 * (2.0 - pressure / 1e6 + 0.03 * (temperature - T0));
            close(
                sample.area_resistance_m2_k_w.iter().sum(),
                expected,
                1e-14,
                1e-20,
            );
            for side in 0..2 {
                close(
                    report.heat.conductances_w_k[side],
                    AREA / (distances[side] / CONDUCTIVITY + sample.area_resistance_m2_k_w[side]),
                    1e-14,
                    1e-20,
                );
            }
            assert_eq!(report.mechanics.wear_normals().len(), 2);
        };
        let first = contact
            .march_prescribed_velocity_candidate(
                &thermal,
                0.04,
                &velocity,
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                march,
                T0,
                transport,
                cx,
            )
            .unwrap();
        assert!(first.report().steps.len() >= 4);
        for step in &first.report().steps {
            check_source(&step.report, [0.5 * LENGTH; 2]);
        }
        close(
            first.contact().source().tangential_states()[12].slip_m()[0],
            velocity[12][0] * 0.04,
            0.0,
            1e-15,
        );
        close(
            first.report().thermal_energy_change_j,
            first.report().frictional_work_j,
            0.0,
            march.march.energy_tolerance_j,
        );
        assert!(first.report().combined_energy_residual_j.abs() <= march.march.energy_tolerance_j);

        let pending = first.contact().unexported_wear_volume_m3();
        assert!(pending.into_iter().all(|volume| volume > 0.0));
        let before_export = first.contact().source().clone();
        let split = first.contact().route().area_resistance_m2_k_w();
        let (geometry, mut patches) = removal::export_geometry(pending, first.thermal());
        let distances = geometry.map(|g| 0.5 * (g.end_m[0] - g.start_m[0]));
        patches[0]
            .replacement
            .as_mut()
            .unwrap()
            .link
            .conductance_w_k =
            series_contact_conductance_w_k(AREA, distances, [CONDUCTIVITY; 2], split.iter().sum())
                .unwrap();
        let exported = first
            .contact()
            .export_wear_candidate(
                first.thermal(),
                geometry,
                &patches,
                [CONDUCTIVITY; 2],
                MaterialRemovalTolerances {
                    mass_kg: 1e-20,
                    volume_m3: 1e-23,
                },
                cx,
            )
            .unwrap();
        assert_eq!(exported.contact().source(), &before_export);
        assert_eq!(exported.contact().unexported_wear_volume_m3(), [0.0; 2]);
        assert_eq!(
            exported
                .contact()
                .contact_resistance()
                .unwrap()
                .original_source(),
            &original
        );
        assert_eq!(exported.contact().route().area_resistance_m2_k_w(), split);
        for id in CELLS {
            assert!(exported.thermal().cells()[&id].mass_kg < first.thermal().cells()[&id].mass_kg);
            assert_eq!(
                exported.thermal().temperature_k(id),
                first.thermal().temperature_k(id)
            );
        }

        let saved_cells = exported.thermal().cells().clone();
        velocity[12][0] = -velocity[12][0];
        let last_temperature = first
            .report()
            .steps
            .last()
            .unwrap()
            .report
            .heat
            .interface_temperature_k;
        let reversed = exported
            .contact()
            .march_prescribed_velocity_candidate(
                exported.thermal(),
                0.08,
                &velocity,
                PlasticSlidingHeatBoundary::HeldApproach,
                controls(),
                march,
                last_temperature,
                transport,
                cx,
            )
            .unwrap();
        for step in &reversed.report().steps {
            check_source(&step.report, distances);
        }
        close(
            reversed.contact().source().tangential_states()[12].slip_m()[0],
            0.0,
            0.0,
            2e-15,
        );
        for side in 0..2 {
            assert!(
                reversed.contact().source().wear_state().volumes_m3()[side]
                    > before_export.wear_state().volumes_m3()[side]
            );
            assert!(reversed.contact().unexported_wear_volume_m3()[side] > 0.0);
        }
        let stored: f64 = CELLS
            .into_iter()
            .map(|id| {
                let old = &saved_cells[&id];
                let new = &reversed.thermal().cells()[&id];
                assert_eq!(new.mass_kg, old.mass_kg);
                old.mass_kg * (new.specific_h_j_kg - old.specific_h_j_kg)
            })
            .sum();
        close(
            stored,
            reversed.report().thermal_energy_change_j,
            1e-12,
            1e-16,
        );
        close(
            stored,
            reversed.report().frictional_work_j,
            0.0,
            march.march.energy_tolerance_j,
        );
        assert!(
            reversed.report().combined_energy_residual_j.abs() <= march.march.energy_tolerance_j
        );
        assert_eq!(
            reversed.contact().route().area_resistance_m2_k_w(),
            reversed
                .report()
                .steps
                .last()
                .unwrap()
                .report
                .contact_resistance
                .as_ref()
                .unwrap()
                .area_resistance_m2_k_w
        );
        assert_eq!(exported.thermal().cells(), &saved_cells);
        assert_eq!(exported.contact().source(), &before_export);
        assert_eq!(reversed.contact().time_s(), 0.08);
    });
}
