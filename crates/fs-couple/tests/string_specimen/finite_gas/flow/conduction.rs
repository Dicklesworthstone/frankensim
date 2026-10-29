//! The existing heat, source-work and nozzle owners share one acoustic sample.

use super::super::electrical::{CP, conducting_string, controls};
use super::*;
use fs_couple::acoustic_realize::finite_gas::{
    FiniteGasStringError, exchange::FiniteGasStringExchangeError,
};
use fs_couple::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use fs_material::conductor::OhmicDrive;

fn stored(
    live: &ThermalMaterialStringRuntime,
    gas: &IsochoricGasParcel,
    reservoir: &ClosedGasParcel,
) -> f64 {
    live.mechanical_internal_energy_j().unwrap()
        + live.mechanical().specimen().mass_kg() * live.thermal().specific_enthalpy_j_kg()
        + gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
        + reservoir.mass_kg() * reservoir.state().specific_internal_energy_j_kg()
}

#[test]
fn g3_conduction_then_native_nozzle_preserves_history_in_both_flow_directions() {
    with_string_cx(|cx, _| {
        for (local_pressure, remote_pressure, direction) in [
            (
                100000.0,
                200000.0,
                FiniteGasStringFlowDirection::IntoSurroundings,
            ),
            (
                200000.0,
                100000.0,
                FiniteGasStringFlowDirection::OutOfSurroundings,
            ),
        ] {
            let mut gas = parcel(0.028, local_pressure);
            let (mut live, card, point, key) =
                conducting_string(cx, gas.state().state(), 310.0, true, 9e-6);
            for epoch in 0..8 {
                live.step_with_isochoric_gas(
                    cx,
                    epoch,
                    &[0.0],
                    &mut gas,
                    FiniteGasStringHeat::default(),
                    solve(),
                )
                .unwrap();
            }
            let mut reservoir = ClosedGasParcel::from_isochoric(parcel_named(
                "introduced-gas",
                0.004,
                remote_pressure,
            ));
            let specimen = live.mechanical().specimen();
            let radius = specimen.radius_m();
            let layers = MaterialConductionLayers {
                material: &card,
                context: &point,
                temperature_axis: "T",
                liquid_fraction_axis: None,
                conductivity: &key,
                selection: MaterialPropertySelection::SingleClaimOnly,
                area_m2: 2.0
                    * core::f64::consts::PI
                    * radius
                    * (specimen.string().length_m + radius),
                gas_distance_m: 0.01,
                body_distance_m: radius / 2.0,
            };
            let before = stored(&live, &gas, &reservoir);
            let heat = FiniteGasStringHeat {
                string_heat_j: 1e-9,
                gas_heat_j: 1e-7,
            };
            let mut thermal_only = live.clone();
            let mut thermal_gas = gas.clone();
            let thermal = thermal_only
                .step_with_isochoric_gas_conduction(
                    cx,
                    8,
                    &[0.0],
                    &mut thermal_gas,
                    heat,
                    solve(),
                    controls(),
                    &layers,
                )
                .unwrap();
            let surrounding = ClosedGasParcel::from_isochoric(thermal_gas);
            let (native, expected_flow) = match direction {
                FiniteGasStringFlowDirection::IntoSurroundings => reservoir.through_nozzle_to(
                    &surrounding,
                    nozzle(),
                    1.0 / 48000.0,
                    flow_solve(),
                    || false,
                ),
                FiniteGasStringFlowDirection::OutOfSurroundings => surrounding.through_nozzle_to(
                    &reservoir,
                    nozzle(),
                    1.0 / 48000.0,
                    flow_solve(),
                    || false,
                ),
            }
            .unwrap();
            let (donor, receiver, expected_transfer) = native.into_parts();
            let (expected_gas, expected_reservoir) = match direction {
                FiniteGasStringFlowDirection::IntoSurroundings => (receiver, donor),
                FiniteGasStringFlowDirection::OutOfSurroundings => (donor, receiver),
            };
            let frame = live
                .step_with_isochoric_gas_nozzle_conduction(
                    cx,
                    8,
                    &[0.0],
                    &mut gas,
                    &mut reservoir,
                    heat,
                    solve(),
                    controls(),
                    &layers,
                    nozzle(),
                    flow_solve(),
                )
                .unwrap();
            assert_eq!(frame.direction, direction);
            assert_eq!(frame.flow, expected_flow);
            assert_eq!(frame.transfer, expected_transfer);
            assert_eq!(gas, expected_gas.into_isochoric());
            assert_eq!(reservoir, expected_reservoir);
            assert_eq!(
                frame.thermal.solid_to_gas_heat_j,
                thermal.solid_to_gas_heat_j
            );
            assert_eq!(
                frame.thermal.source.query_point(),
                thermal.source.query_point()
            );
            assert_eq!(frame.thermal.coupled.gas, thermal.coupled.gas);
            assert_ne!(frame.thermal.solid_to_gas_heat_j, 0.0);
            assert_eq!(live.thermal(), thermal_only.thermal());
            assert_eq!(
                live.mechanical().states(),
                thermal_only.mechanical().states()
            );
            assert_eq!(live.mechanical().ambient(), gas.state().state());
            assert_eq!(live.mechanical().accepted_samples(), 9);
            assert_eq!(live.mechanical().epoch(), 9);
            let sample = &frame.thermal.coupled.string.vibration;
            assert_ne!(
                sample.acoustic.observer_pressure_pa,
                thermal
                    .coupled
                    .string
                    .vibration
                    .acoustic
                    .observer_pressure_pa
            );
            let residual = sample.acoustic.input_work_j + heat.string_heat_j + heat.gas_heat_j
                - (stored(&live, &gas, &reservoir) - before)
                - sample.dissipation.authored_loss_j
                - sample.dissipation.roundoff_residual_j;
            assert!(residual.abs() <= frame.energy_tolerance_j);
            assert!(frame.energy_balance_residual_j.abs() <= frame.energy_tolerance_j);
            live.step(cx, 9, &[0.0], 0.0).unwrap();
            thermal_only.step(cx, 9, &[0.0], 0.0).unwrap();
            assert_ne!(
                live.mechanical().states(),
                thermal_only.mechanical().states()
            );
        }
    });
}

#[test]
fn g1_powered_nozzle_retains_closed_limit_and_counts_joule_work_once() {
    with_string_cx(|cx, _| {
        let gas = parcel(0.028, 100000.0);
        let reservoir =
            ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.004, 200000.0));
        let (initial, card, point, key) =
            conducting_string(cx, gas.state().state(), 310.0, true, 9e-6);
        let specimen = initial.mechanical().specimen();
        let radius = specimen.radius_m();
        let layers = MaterialConductionLayers {
            material: &card,
            context: &point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
            conductivity: &key,
            selection: MaterialPropertySelection::SingleClaimOnly,
            area_m2: 2.0 * core::f64::consts::PI * radius * (specimen.string().length_m + radius),
            gas_distance_m: 0.01,
            body_distance_m: radius / 2.0,
        };
        let drive = OhmicDrive::Current(0.05);
        let heat = FiniteGasStringHeat {
            string_heat_j: 1e-9,
            gas_heat_j: 1e-7,
        };
        let mut closed = initial.clone();
        let mut closed_gas = gas.clone();
        let reference = closed
            .step_with_drive_and_isochoric_gas_conduction(
                cx,
                0,
                &[0.0],
                drive,
                &mut closed_gas,
                heat,
                solve(),
                controls(),
                &layers,
            )
            .unwrap();
        for area_m2 in [0.0, nozzle().area_m2] {
            let mut live = initial.clone();
            let mut local = gas.clone();
            let mut source = reservoir.clone();
            let frame = live
                .step_with_drive_and_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    drive,
                    &mut local,
                    &mut source,
                    heat,
                    solve(),
                    controls(),
                    &layers,
                    GasNozzle {
                        area_m2,
                        ..nozzle()
                    },
                    flow_solve(),
                )
                .unwrap();
            let flow = &frame.coupled;
            let thermal = &flow.thermal;
            let sample = &thermal.coupled.string.vibration;
            assert_eq!(frame.joule_heat_j, reference.joule_heat_j);
            assert_eq!(frame.current_a, reference.current_a);
            assert_eq!(frame.final_current_a, reference.final_current_a);
            assert_eq!(frame.external_string_heat_j, heat.string_heat_j);
            assert_eq!(
                thermal.external_string_heat_j,
                heat.string_heat_j + frame.joule_heat_j
            );
            assert_eq!(
                thermal.source.query_point(),
                frame.conductor_after.material().query_point()
            );
            let r0 = 3e-6 * specimen.string().length_m / specimen.area_m2();
            close(frame.joule_heat_j, 0.05_f64.powi(2) * r0 / 48000.0);
            let residual = frame.joule_heat_j
                + heat.string_heat_j
                + heat.gas_heat_j
                + sample.acoustic.input_work_j
                - (stored(&live, &local, &source) - stored(&initial, &gas, &reservoir))
                - sample.dissipation.authored_loss_j
                - sample.dissipation.roundoff_residual_j;
            assert!(residual.abs() <= flow.energy_tolerance_j);
            assert!(flow.energy_balance_residual_j.abs() <= flow.energy_tolerance_j);
            assert_eq!(live.mechanical().ambient(), local.state().state());
            assert_eq!(live.mechanical().epoch(), 1);
            assert_eq!(live.mechanical().accepted_samples(), 1);
            if area_m2 == 0.0 {
                assert_eq!(flow.transfer.transferred_mass_kg, 0.0);
                assert_eq!(thermal.coupled.string, reference.coupled.coupled.string);
                assert_eq!(local, closed_gas);
                assert_eq!(source, reservoir);
            } else {
                assert!(flow.transfer.transferred_mass_kg > 0.0);
                assert_ne!(
                    sample.acoustic.observer_pressure_pa,
                    reference
                        .coupled
                        .coupled
                        .string
                        .vibration
                        .acoustic
                        .observer_pressure_pa
                );
            }
        }
    });
}

#[test]
fn g4_late_conductor_and_final_flow_temperature_refusals_preserve_three_owners() {
    with_string_cx(|cx, _| {
        for electrical in [true, false] {
            let gas = parcel(0.028, 100000.0);
            let reservoir =
                ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.004, 200000.0));
            let (original, card, point, key) = conducting_string(
                cx,
                gas.state().state(),
                if electrical { 300.0 } else { 310.0 },
                !electrical,
                if electrical { 0.0 } else { 9e-6 },
            );
            let specimen = original.mechanical().specimen();
            let radius = specimen.radius_m();
            let layers = MaterialConductionLayers {
                material: &card,
                context: &point,
                temperature_axis: "T",
                liquid_fraction_axis: None,
                conductivity: &key,
                selection: MaterialPropertySelection::SingleClaimOnly,
                // An isolated layer reaches exactly the authored zero-resistivity
                // endpoint; open flow still completes before conductor refusal.
                area_m2: if electrical {
                    0.0
                } else {
                    2.0 * core::f64::consts::PI * radius * (specimen.string().length_m + radius)
                },
                gas_distance_m: 0.01,
                body_distance_m: radius / 2.0,
            };
            let heat = FiniteGasStringHeat {
                string_heat_j: if electrical {
                    specimen.mass_kg() * CP * 40.0
                } else {
                    0.0
                },
                gas_heat_j: 0.0,
            };
            let limits = FiniteGasStringSolve {
                maximum_gas_temperature_increment_k: if electrical { 100.0 } else { 1e-4 },
                ..solve()
            };
            // Both cases admit the complete conduction proposal. Refusal must
            // occur in subsequent flow or final electrical admission.
            let mut thermal_only = original.clone();
            let mut thermal_gas = gas.clone();
            thermal_only
                .step_with_isochoric_gas_conduction(
                    cx,
                    0,
                    &[0.0],
                    &mut thermal_gas,
                    heat,
                    limits,
                    controls(),
                    &layers,
                )
                .unwrap();
            let mut live = original.clone();
            let mut local = gas.clone();
            let mut source = reservoir.clone();
            let refused = live.step_with_drive_and_isochoric_gas_nozzle_conduction(
                cx,
                0,
                &[0.0],
                OhmicDrive::Current(0.0),
                &mut local,
                &mut source,
                heat,
                limits,
                controls(),
                &layers,
                nozzle(),
                flow_solve(),
            );
            if electrical {
                assert!(matches!(
                    refused,
                    Err(FiniteGasStringExchangeError::String(
                        FiniteGasStringError::String(ThermalStringError::Electrical(_))
                    ))
                ));
            } else {
                assert!(matches!(
                    refused,
                    Err(FiniteGasStringExchangeError::String(
                        FiniteGasStringError::Admission(_)
                    ))
                ));
            }
            assert_eq!(local, gas);
            assert_eq!(source, reservoir);
            assert_eq!(live.thermal(), original.thermal());
            assert_eq!(
                live.mechanical().specimen(),
                original.mechanical().specimen()
            );
            assert_eq!(live.mechanical().states(), original.mechanical().states());
            assert_eq!(live.mechanical().ambient(), original.mechanical().ambient());
            assert_eq!(live.mechanical().accepted_samples(), 0);
            assert_eq!(live.mechanical().epoch(), 0);
            let valid = FiniteGasStringHeat {
                string_heat_j: heat.string_heat_j / 40.0,
                ..heat
            };
            let mut replay = original.clone();
            let mut replay_gas = gas.clone();
            let mut replay_source = reservoir.clone();
            let retry = live
                .step_with_drive_and_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    OhmicDrive::Current(0.0),
                    &mut local,
                    &mut source,
                    valid,
                    solve(),
                    controls(),
                    &layers,
                    nozzle(),
                    flow_solve(),
                )
                .unwrap();
            let expected = replay
                .step_with_drive_and_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    OhmicDrive::Current(0.0),
                    &mut replay_gas,
                    &mut replay_source,
                    valid,
                    solve(),
                    controls(),
                    &layers,
                    nozzle(),
                    flow_solve(),
                )
                .unwrap();
            assert_eq!(
                retry.coupled.thermal.coupled.string,
                expected.coupled.thermal.coupled.string
            );
            assert_eq!(retry.coupled.transfer, expected.coupled.transfer);
            assert_eq!(local, replay_gas);
            assert_eq!(source, replay_source);
        }
    });
}
