//! One endpoint closes electrical work, chemical extent and conductive heat.

use super::*;
use super::super::super::reactive::{
    RELEASE, controls as reacting_controls, gas as reacting_gas, law, policy,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::QuantitySpec;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn reservoir(gas: &IsochoricGasParcel) -> ClosedGasParcel {
    ClosedGasParcel::new(
        gas.model().clone(),
        gas.point()
            .clone()
            .with_quantity(
                "pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                200000.0,
            )
            .unwrap(),
        gas.volume_m3(),
    )
    .unwrap()
}

#[test]
fn g1_reacting_current_and_voltage_close_all_three_endpoint_equations() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let gas = reacting_gas(&reaction);
        let (initial, card, point, key) =
            conducting_string(cx, gas.state().state(), 310.0, false, 9e-6);
        let mut layers = layer(&initial, &card, &point, &key, true);
        layers.gas_distance_m = 0.03;
        let cs = initial.mechanical().specimen().mass_kg() * CP;
        let n = gas.mass_kg() / 0.028;
        let cg = n * 2.5 * R;
        let dt = 1.0 / f64::from(RATE);
        let r0 = 1e-6 * initial.mechanical().specimen().string().length_m
            / initial.mechanical().specimen().area_m2();
        let limits = feedback(5.0 * cs);
        for drive in [
            OhmicDrive::Current(0.0),
            OhmicDrive::Current(0.5),
            OhmicDrive::Voltage(0.5 * r0),
        ] {
            let mut live = initial.clone();
            let mut local = gas.clone();
            let result = live
                .step_with_implicit_drive_and_reacting_isochoric_gas_conduction(
                    cx,
                    0,
                    &[0.0],
                    drive,
                    &mut local,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    reacting_controls(),
                    controls(),
                    &layers,
                    limits,
                )
                .unwrap();
            let coupled = &result.electrical.coupled;
            let thermal = &coupled.coupled;
            let ts = live.thermal().temperature_k();
            let tg = local.state().state().temperature;
            let xi = coupled.chemistry.reaction.extent_mol;

            // Closed-form backward-Euler isomer extent at this temperature.
            // Both species have the same mass and cp, differing in formation energy.
            let k = (RELEASE / (R * tg)).exp();
            let a = dt * 200000.0 * R * tg / 100000.0;
            let expected_xi =
                a * n * (0.9 - 0.1 / k) / (1.0 + a * (1.0 + 1.0 / k));
            assert!((xi - expected_xi).abs() < 3e-12);
            let x = 0.9 - xi / n;
            let y = 0.1 + xi / n;

            // Independent equal-molar-mass Wilke/WMS reduction of the two
            // authored transport curves; no production constitutive query.
            let phi_ab = (1.0 + (1.7_f64 / 4.0).sqrt()).powi(2) / 4.0;
            let phi_ba = (1.0 + (4.0_f64 / 1.7).sqrt()).powi(2) / 4.0;
            let kg = (tg + 500.0) / 700.0
                * (x * 0.03 / (x + y * phi_ab) + y * 0.09 / (y + x * phi_ba));
            let conductance = layers.area_m2
                / (layers.gas_distance_m / kg + layers.body_distance_m / 20.0);
            let q = dt * conductance * (ts - tg);
            assert!((thermal.solid_to_gas_heat_j - q).abs() < 2e-12);
            assert!((thermal.conductance_w_k - conductance).abs() < 2e-12);
            assert!(
                (ts - 310.0 - (result.electrical.joule_heat_j - q) / cs).abs()
                    < 2e-6
            );
            assert!((tg - 300.0 - (q + RELEASE * xi) / cg).abs() < 2e-8);
            let resistance = r0 * (1.0 + 0.2 * (ts - 300.0));
            let expected_work = match drive {
                OhmicDrive::Current(current) => dt * current * current * resistance,
                OhmicDrive::Voltage(voltage) => dt * voltage * voltage / resistance,
            };
            assert!((result.source_work_j - expected_work).abs() < 2e-12);
            assert!(result.source_residual_j.abs() <= limits.source_tolerance_j);
            assert!(
                result.energy_balance_residual_j.abs() <= result.energy_tolerance_j
            );
            assert_eq!(thermal.coupled.drag_heat_j, 0.0);
            assert_eq!(live.mechanical().epoch(), 1);
            assert_eq!(live.mechanical().accepted_samples(), 1);

            if drive == OhmicDrive::Current(0.0) {
                let mut direct = initial.clone();
                let mut direct_gas = gas.clone();
                let original = direct
                    .step_with_reacting_isochoric_gas_conduction(
                        cx,
                        0,
                        &[0.0],
                        &mut direct_gas,
                        FiniteGasStringHeat::default(),
                        &reaction,
                        reacting_controls(),
                        controls(),
                        &layers,
                    )
                    .unwrap();
                assert_eq!(result.electrical.joule_heat_j, 0.0);
                assert_eq!(result.source_work_j, 0.0);
                assert_eq!(coupled.chemistry, original.chemistry);
                assert_eq!(thermal.coupled.string, original.coupled.coupled.string);
                assert_eq!(live.thermal(), direct.thermal());
                assert_eq!(local, direct_gas);
            }
        }
    });
}

#[test]
fn g3_reacting_feedback_transports_once_after_convergence_and_audits_both_volumes() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let gas = reacting_gas(&reaction);
        let source = reservoir(&gas);
        let (initial, card, point, key) =
            conducting_string(cx, gas.state().state(), 310.0, true, 9e-6);
        let mut layers = layer(&initial, &card, &point, &key, true);
        layers.gas_distance_m = 0.03;
        let limits = feedback(5.0 * initial.mechanical().specimen().mass_kg() * CP);
        let mut closed = initial.clone();
        let mut heated = gas.clone();
        let reference = closed
            .step_with_implicit_drive_and_reacting_isochoric_gas_conduction(
                cx,
                0,
                &[0.0],
                OhmicDrive::Current(0.5),
                &mut heated,
                FiniteGasStringHeat::default(),
                &reaction,
                reacting_controls(),
                controls(),
                &layers,
                limits,
            )
            .unwrap();
        assert!(reference.trials > 1);
        for area_m2 in [0.0, nozzle().area_m2] {
            let nozzle = GasNozzle {
                area_m2,
                ..nozzle()
            };
            let (native, expected_flow) = source
                .through_nozzle_to(
                    &ClosedGasParcel::from_isochoric(heated.clone()),
                    nozzle,
                    1.0 / f64::from(RATE),
                    flow_solve(),
                    || false,
                )
                .unwrap();
            let (expected_source, expected_gas, expected_transfer) = native.into_parts();
            let mut live = initial.clone();
            let mut local = gas.clone();
            let mut reservoir = source.clone();
            let result = live
                .step_with_implicit_drive_and_reacting_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    OhmicDrive::Current(0.5),
                    &mut local,
                    &mut reservoir,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    reacting_controls(),
                    controls(),
                    &layers,
                    nozzle,
                    flow_solve(),
                    limits,
                )
                .unwrap();
            let flow = &result.electrical.coupled;
            assert_eq!(flow.flow, expected_flow);
            assert_eq!(flow.transfer, expected_transfer);
            assert_eq!(local, expected_gas.into_isochoric());
            assert_eq!(reservoir, expected_source);
            assert_eq!(
                flow.thermal.chemistry,
                reference.electrical.coupled.chemistry
            );
            assert_eq!(result.source_work_j, reference.source_work_j);
            assert_eq!(result.trials, reference.trials);
            assert_eq!(live.mechanical().states(), closed.mechanical().states());
            assert_eq!(live.mechanical().ambient(), local.state().state());
            assert_eq!(live.mechanical().epoch(), 1);
            assert_eq!(live.mechanical().accepted_samples(), 1);

            let gas_storage = local.mass_kg()
                * local.state().specific_internal_energy_j_kg()
                - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
            let reservoir_storage = reservoir.mass_kg()
                * reservoir.state().specific_internal_energy_j_kg()
                - source.mass_kg() * source.state().specific_internal_energy_j_kg();
            let solid_storage = initial.mechanical().specimen().mass_kg()
                * (live.thermal().specific_enthalpy_j_kg()
                    - initial.thermal().specific_enthalpy_j_kg());
            let mechanical_storage = live.mechanical_internal_energy_j().unwrap()
                - initial.mechanical_internal_energy_j().unwrap();
            let loss = &flow.thermal.coupled.coupled.string.vibration.dissipation;
            let defect = result.source_work_j
                - gas_storage
                - reservoir_storage
                - solid_storage
                - mechanical_storage
                - loss.authored_loss_j
                - loss.roundoff_residual_j;
            assert!(defect.abs() <= result.energy_tolerance_j);
            assert!(
                result.energy_balance_residual_j.abs() <= result.energy_tolerance_j
            );
            if area_m2 == 0.0 {
                assert_eq!(local, heated);
                assert_eq!(reservoir, source);
                assert_eq!(
                    flow.thermal.coupled.coupled.string,
                    reference.electrical.coupled.coupled.coupled.string
                );
            } else {
                assert!(flow.transfer.transferred_mass_kg > 0.0);
            }
        }
    });
}

#[test]
fn g4_reacting_source_refusals_preserve_inventory_and_continuing_sound_on_retry() {
    with_string_cx(|cx, _| {
        let reaction = law(RELEASE).with_implicit_mechanism(policy()).unwrap();
        let gas = reacting_gas(&reaction);
        let source = reservoir(&gas);
        let (initial, card, point, key) =
            conducting_string(cx, gas.state().state(), 310.0, true, 9e-6);
        let mut layers = layer(&initial, &card, &point, &key, true);
        layers.gas_distance_m = 0.03;
        let limits = feedback(5.0 * initial.mechanical().specimen().mass_kg() * CP);
        let mut clean = initial.clone();
        let mut clean_gas = gas.clone();
        let mut clean_source = source.clone();
        let expected = clean
            .step_with_implicit_drive_and_reacting_isochoric_gas_nozzle_conduction(
                cx,
                0,
                &[0.0],
                OhmicDrive::Current(0.5),
                &mut clean_gas,
                &mut clean_source,
                FiniteGasStringHeat::default(),
                &reaction,
                reacting_controls(),
                controls(),
                &layers,
                nozzle(),
                flow_solve(),
                limits,
            )
            .unwrap();
        for case in 0..3 {
            let mut live = initial.clone();
            let mut local = gas.clone();
            let mut reservoir = source.clone();
            let mut chemical = reacting_controls();
            let electrical = if case == 0 {
                feedback(1e-10)
            } else {
                limits
            };
            if case == 1 {
                chemical.reaction.maximum_consumed_fraction = 1e-8;
            }
            if case == 2 {
                reservoir = ClosedGasParcel::from_isochoric(parcel_named(
                    "unlisted-inert",
                    0.028,
                    200000.0,
                ));
            }
            let before_source = reservoir.clone();
            assert!(
                live.step_with_implicit_drive_and_reacting_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    OhmicDrive::Current(0.5),
                    &mut local,
                    &mut reservoir,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    chemical,
                    controls(),
                    &layers,
                    nozzle(),
                    flow_solve(),
                    electrical,
                )
                .is_err()
            );
            assert_eq!(local, gas);
            assert_eq!(reservoir, before_source);
            assert_eq!(live.thermal(), initial.thermal());
            assert_eq!(live.mechanical().states(), initial.mechanical().states());
            assert_eq!(live.mechanical().ambient(), initial.mechanical().ambient());
            assert_eq!(live.mechanical().epoch(), 0);
            assert_eq!(live.mechanical().accepted_samples(), 0);

            reservoir = source.clone();
            let retry = live
                .step_with_implicit_drive_and_reacting_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    OhmicDrive::Current(0.5),
                    &mut local,
                    &mut reservoir,
                    FiniteGasStringHeat::default(),
                    &reaction,
                    reacting_controls(),
                    controls(),
                    &layers,
                    nozzle(),
                    flow_solve(),
                    limits,
                )
                .unwrap();
            assert_eq!(local, clean_gas);
            assert_eq!(reservoir, clean_source);
            assert_eq!(live.thermal(), clean.thermal());
            assert_eq!(live.mechanical().states(), clean.mechanical().states());
            assert_eq!(retry.source_work_j, expected.source_work_j);
            assert_eq!(
                retry.electrical.coupled.thermal.chemistry,
                expected.electrical.coupled.thermal.chemistry
            );
            assert_eq!(
                retry.electrical.coupled.thermal.coupled.coupled.string,
                expected.electrical.coupled.thermal.coupled.coupled.string
            );
        }
    });
}
