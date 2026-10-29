//! Conserved species/enthalpy changes must reach the continuing physical voice.
//! Synthetic gas sources test coupling; they are not measured real-gas data.

use super::*;
use fs_couple::acoustic_realize::finite_gas::flow::FiniteGasStringFlowDirection;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleStepSolve,
};

#[path = "flow/conduction.rs"]
mod conduction;

fn nozzle() -> GasNozzle {
    GasNozzle {
        area_m2: 1e-4,
        discharge_coefficient: 0.8,
    }
}

fn flow_solve() -> GasNozzleStepSolve {
    GasNozzleStepSolve {
        nozzle: solve().gas,
        transfer: GasMassTransferSolve {
            donor: solve().gas,
            receiver: solve().gas,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        maximum_donor_fraction: 0.1,
    }
}

#[test]
fn conserved_inlet_mixture_drives_loading_and_observer_without_resetting_history() {
    with_string_cx(|cx, _| {
        let mut gas = parcel(0.028, 100000.0);
        let mut live = runtime(cx, gas.state().state());
        for epoch in 0..16 {
            live.step_with_isochoric_gas(
                cx,
                epoch,
                &[0.0; 2],
                &mut gas,
                FiniteGasStringHeat::default(),
                solve(),
            )
            .unwrap();
        }
        let mut reservoir =
            ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.004, 200000.0));
        let original = live.clone();
        let original_gas = gas.clone();
        let original_reservoir = reservoir.clone();
        let previous = accelerations(&original);
        let heat = FiniteGasStringHeat {
            string_heat_j: 5e-8,
            gas_heat_j: 0.25,
        };

        // Reference composition: each existing physics owner advances once.
        let mut frozen = original.clone();
        let mut heated_gas = original_gas.clone();
        let heated = frozen
            .step_with_isochoric_gas(cx, 16, &[0.0; 2], &mut heated_gas, heat, solve())
            .unwrap();
        let (native, native_flow) = reservoir
            .through_nozzle_to(
                &ClosedGasParcel::from_isochoric(heated_gas),
                nozzle(),
                1.0 / 48000.0,
                flow_solve(),
                || false,
            )
            .unwrap();
        let (expected_reservoir, expected_gas, expected_transfer) = native.into_parts();
        let frame = live
            .step_with_isochoric_gas_nozzle(
                cx,
                16,
                &[0.0; 2],
                &mut gas,
                &mut reservoir,
                heat,
                solve(),
                nozzle(),
                flow_solve(),
            )
            .unwrap();
        assert_eq!(
            frame.direction,
            FiniteGasStringFlowDirection::IntoSurroundings
        );
        assert_eq!(frame.transfer, expected_transfer);
        assert_eq!(frame.flow, native_flow);
        assert_eq!(gas, expected_gas.into_isochoric());
        assert_eq!(reservoir, expected_reservoir);
        assert_eq!(frame.thermal.gas, heated.gas);
        assert_eq!(frame.thermal.drag_heat_j, heated.drag_heat_j);
        assert_eq!(
            frame.thermal.gas.supplied_heat_j,
            heat.gas_heat_j + frame.thermal.drag_heat_j
        );
        assert!(frame.transfer.transferred_mass_kg > 0.0);
        assert_eq!(gas.state().components().len(), 2);
        assert_eq!(gas.point().axes().get("material-lot"), Some(&7.0));
        assert_eq!(live.mechanical().ambient(), gas.state().state());
        assert_eq!(live.mechanical().epoch(), 17);
        assert_eq!(live.mechanical().accepted_samples(), 17);
        assert_eq!(live.mechanical().states(), frozen.mechanical().states());
        assert_eq!(live.thermal(), frozen.thermal());
        assert_ne!(live.mechanical().modes(), frozen.mechanical().modes());

        // Independent endpoint energy includes variable mass in BOTH volumes;
        // transported upstream enthalpy is not a third external heat source.
        let local_storage = gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
            - original_gas.mass_kg() * original_gas.state().specific_internal_energy_j_kg();
        let reservoir_storage = reservoir.mass_kg()
            * reservoir.state().specific_internal_energy_j_kg()
            - original_reservoir.mass_kg()
                * original_reservoir.state().specific_internal_energy_j_kg();
        let solid_storage = live.mechanical().specimen().mass_kg()
            * (live.thermal().specific_enthalpy_j_kg()
                - original.thermal().specific_enthalpy_j_kg());
        let mechanical_storage = live.mechanical_internal_energy_j().unwrap()
            - original.mechanical_internal_energy_j().unwrap();
        let defect = heat.string_heat_j
            + heat.gas_heat_j
            + frame.thermal.string.vibration.acoustic.input_work_j
            - (local_storage
                + reservoir_storage
                + solid_storage
                + mechanical_storage
                + frame.thermal.string.vibration.dissipation.authored_loss_j
                + frame
                    .thermal
                    .string
                    .vibration
                    .dissipation
                    .roundoff_residual_j);
        assert!(defect.abs() <= frame.energy_tolerance_j);
        assert!(frame.energy_balance_residual_j.abs() <= frame.energy_tolerance_j);
        assert!((local_storage - frame.surrounding_energy_change_j).abs() < 1e-10);
        assert!((reservoir_storage - frame.reservoir_energy_change_j).abs() < 1e-10);
        assert!(
            (gas.mass_kg() + reservoir.mass_kg()
                - original_gas.mass_kg()
                - original_reservoir.mass_kg())
            .abs()
                <= flow_solve().transfer.mass_tolerance_kg
        );

        // The dipole uses the previous ACCEPTED acceleration, not the discarded
        // heat-stage observation; both density and sound speed are final ones.
        let string = live.mechanical().specimen().string();
        let medium = gas.state().state();
        let pi = core::f64::consts::PI;
        let acceleration = accelerations(&live);
        let factor = medium.density * string.width_m
            / (4.0 * pi * (string.lin_density_kg_m * string.length_m / 2.0).sqrt());
        let pressure = factor * 2.0 * string.length_m / pi * acceleration[0]
            - factor * string.length_m.powi(2) / (2.0 * pi * medium.sound_speed)
                * (acceleration[1] - previous[1])
                * 48000.0;
        close(
            frame.thermal.string.vibration.acoustic.observer_pressure_pa,
            pressure,
        );
        assert_ne!(
            frame.thermal.string.vibration.acoustic.observer_pressure_pa,
            heated.string.vibration.acoustic.observer_pressure_pa
        );
        live.step(cx, 17, &[0.0; 2], 0.0).unwrap();
        frozen.step(cx, 17, &[0.0; 2], 0.0).unwrap();
        assert_ne!(live.mechanical().states(), frozen.mechanical().states());
    });
}

#[test]
fn source_data_only_species_substitution_changes_the_continuing_motion_and_sound() {
    with_string_cx(|cx, _| {
        let mut light_gas = parcel(0.028, 100000.0);
        let mut heavy_gas = light_gas.clone();
        let mut light = runtime(cx, light_gas.state().state());
        let mut heavy = light.clone();
        // Only the introduced species' authored molar mass differs. Geometry,
        // initial string, heat, nozzle and solver are otherwise identical.
        let mut light_source =
            ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.004, 200000.0));
        let mut heavy_source =
            ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.044, 200000.0));
        let mut final_pressure = (0.0, 0.0);
        for epoch in 0..32 {
            let left = light
                .step_with_isochoric_gas_nozzle(
                    cx,
                    epoch,
                    &[0.0; 2],
                    &mut light_gas,
                    &mut light_source,
                    FiniteGasStringHeat::default(),
                    solve(),
                    nozzle(),
                    flow_solve(),
                )
                .unwrap();
            let right = heavy
                .step_with_isochoric_gas_nozzle(
                    cx,
                    epoch,
                    &[0.0; 2],
                    &mut heavy_gas,
                    &mut heavy_source,
                    FiniteGasStringHeat::default(),
                    solve(),
                    nozzle(),
                    flow_solve(),
                )
                .unwrap();
            if epoch == 0 {
                // Atmosphere changes do not reset or kick the initial motion.
                assert_eq!(light.mechanical().states(), heavy.mechanical().states());
            }
            assert!(left.energy_balance_residual_j.abs() <= left.energy_tolerance_j);
            assert!(right.energy_balance_residual_j.abs() <= right.energy_tolerance_j);
            final_pressure = (
                left.thermal.string.vibration.acoustic.observer_pressure_pa,
                right.thermal.string.vibration.acoustic.observer_pressure_pa,
            );
        }
        assert_eq!(light.mechanical().accepted_samples(), 32);
        assert_eq!(heavy.mechanical().accepted_samples(), 32);
        assert_eq!(light.mechanical().ambient(), light_gas.state().state());
        assert_eq!(heavy.mechanical().ambient(), heavy_gas.state().state());
        assert_ne!(light.mechanical().states(), heavy.mechanical().states());
        assert_ne!(light.mechanical().modes(), heavy.mechanical().modes());
        assert!((final_pressure.0 - final_pressure.1).abs() > 1e-12);
        assert!(light_gas.state().state().sound_speed > heavy_gas.state().state().sound_speed);
    });
}

#[test]
fn reverse_pressure_uses_native_upstream_owner_and_closed_nozzle_matches_heat_only_step() {
    with_string_cx(|cx, _| {
        let gas = parcel(0.028, 200000.0);
        let live = runtime(cx, gas.state().state());
        let reservoir =
            ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.004, 100000.0));
        let mut heated_gas = gas.clone();
        let mut heated_string = live.clone();
        let heated = heated_string
            .step_with_isochoric_gas(
                cx,
                0,
                &[0.0; 2],
                &mut heated_gas,
                FiniteGasStringHeat::default(),
                solve(),
            )
            .unwrap();
        for area in [0.0, nozzle().area_m2] {
            let mut actual = live.clone();
            let mut actual_gas = gas.clone();
            let mut actual_reservoir = reservoir.clone();
            let geometry = GasNozzle {
                area_m2: area,
                ..nozzle()
            };
            let (native, expected_flow) = ClosedGasParcel::from_isochoric(heated_gas.clone())
                .through_nozzle_to(&reservoir, geometry, 1.0 / 48000.0, flow_solve(), || false)
                .unwrap();
            let (expected_gas, expected_reservoir, expected_transfer) = native.into_parts();
            let frame = actual
                .step_with_isochoric_gas_nozzle(
                    cx,
                    0,
                    &[0.0; 2],
                    &mut actual_gas,
                    &mut actual_reservoir,
                    FiniteGasStringHeat::default(),
                    solve(),
                    geometry,
                    flow_solve(),
                )
                .unwrap();
            assert_eq!(
                frame.direction,
                FiniteGasStringFlowDirection::OutOfSurroundings
            );
            assert_eq!(frame.flow, expected_flow);
            assert_eq!(frame.transfer, expected_transfer);
            assert_eq!(actual_gas, expected_gas.into_isochoric());
            assert_eq!(actual_reservoir, expected_reservoir);
            assert_eq!(
                actual.mechanical().states(),
                heated_string.mechanical().states()
            );
            assert_eq!(actual.mechanical().ambient(), actual_gas.state().state());
            if area == 0.0 {
                assert_eq!(frame.transfer.transferred_mass_kg, 0.0);
                assert_eq!(frame.thermal.string, heated.string);
                assert_eq!(actual_gas, heated_gas);
                assert_eq!(actual_reservoir, reservoir);
            } else {
                assert!(actual_gas.mass_kg() < gas.mass_kg());
                assert!(actual_reservoir.mass_kg() > reservoir.mass_kg());
                assert_eq!(actual_gas.state().components().len(), 1);
                assert_eq!(actual_reservoir.state().components().len(), 2);
            }
        }
    });
}

#[test]
fn late_species_flow_and_temperature_refusals_preserve_three_owners_and_observer_history() {
    with_string_cx(|cx, gate| {
        let mut gas = parcel(0.028, 100000.0);
        let mut pristine = runtime(cx, gas.state().state());
        for epoch in 0..16 {
            pristine
                .step_with_isochoric_gas(
                    cx,
                    epoch,
                    &[0.0; 2],
                    &mut gas,
                    FiniteGasStringHeat::default(),
                    solve(),
                )
                .unwrap();
        }
        for case in 0..4 {
            let mut live = pristine.clone();
            let mut local = gas.clone();
            let mut reservoir = ClosedGasParcel::from_isochoric(if case == 0 {
                // Same canonical species, conflicting original caloric law.
                parcel(0.004, 200000.0)
            } else {
                parcel_named("introduced-gas", 0.004, 200000.0)
            });
            let original_reservoir = reservoir.clone();
            let thermal = FiniteGasStringSolve {
                maximum_gas_temperature_increment_k: if case == 2 { 1e-4 } else { 100.0 },
                ..solve()
            };
            let flow = GasNozzleStepSolve {
                maximum_donor_fraction: if case == 1 { 1e-12 } else { 0.1 },
                ..flow_solve()
            };
            assert!(
                live.step_with_isochoric_gas_nozzle(
                    cx,
                    if case == 3 { 17 } else { 16 },
                    &[0.0; 2],
                    &mut local,
                    &mut reservoir,
                    FiniteGasStringHeat::default(),
                    thermal,
                    nozzle(),
                    flow,
                )
                .is_err()
            );
            assert_eq!(local, gas);
            assert_eq!(reservoir, original_reservoir);
            assert_eq!(live.mechanical().epoch(), 16);
            assert_eq!(live.mechanical().accepted_samples(), 16);
            assert_eq!(live.mechanical().states(), pristine.mechanical().states());
            assert_eq!(live.mechanical().ambient(), pristine.mechanical().ambient());
            assert_eq!(live.thermal(), pristine.thermal());
            let mut control = pristine.clone();
            assert_eq!(
                live.step(cx, 16, &[0.0; 2], 0.0).unwrap(),
                control.step(cx, 16, &[0.0; 2], 0.0).unwrap()
            );
        }
        let mut live = pristine.clone();
        let mut local = gas.clone();
        let mut reservoir =
            ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.004, 200000.0));
        let original_reservoir = reservoir.clone();
        gate.request();
        assert!(
            live.step_with_isochoric_gas_nozzle(
                cx,
                16,
                &[0.0; 2],
                &mut local,
                &mut reservoir,
                FiniteGasStringHeat::default(),
                solve(),
                nozzle(),
                flow_solve(),
            )
            .is_err()
        );
        assert_eq!(local, gas);
        assert_eq!(reservoir, original_reservoir);
        assert_eq!(live.mechanical().states(), pristine.mechanical().states());
        assert_eq!(live.thermal(), pristine.thermal());
        assert_eq!(live.mechanical().epoch(), 16);
    });
}
