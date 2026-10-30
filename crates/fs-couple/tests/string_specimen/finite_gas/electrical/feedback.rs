//! Independent endpoint electrical/thermal balances and continuing gas flow.

use super::super::flow::{flow_solve, nozzle};
use super::*;
use fs_couple::acoustic_realize::finite_gas::electrical::feedback::{
    FiniteGasElectricalFeedbackError, FiniteGasElectricalFeedbackSolve,
};
use fs_couple::acoustic_realize::finite_gas::{
    FiniteGasStringError, exchange::FiniteGasStringExchangeError,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;

fn feedback(maximum_joule_work_j: f64) -> FiniteGasElectricalFeedbackSolve {
    FiniteGasElectricalFeedbackSolve {
        maximum_joule_work_j,
        source_tolerance_j: 1e-13,
        maximum_iterations: 64,
    }
}

fn layer<'a>(
    live: &ThermalMaterialStringRuntime,
    card: &'a MaterialCard,
    point: &'a QueryPoint,
    key: &'a PropertyKey,
    conducting: bool,
) -> MaterialConductionLayers<'a> {
    let specimen = live.mechanical().specimen();
    let radius = specimen.radius_m();
    MaterialConductionLayers {
        material: card,
        context: point,
        temperature_axis: "T",
        liquid_fraction_axis: None,
        conductivity: key,
        selection: MaterialPropertySelection::SingleClaimOnly,
        area_m2: if conducting {
            2.0 * core::f64::consts::PI * radius * (specimen.string().length_m + radius)
        } else {
            0.0
        },
        gas_distance_m: 0.01,
        body_distance_m: radius / 2.0,
    }
}

#[test]
fn g1_endpoint_current_and_voltage_match_independent_affine_resistance_roots() {
    with_string_cx(|cx, _| {
        let gas = parcel(0.028, 100000.0);
        let (initial, card, point, key) =
            conducting_string(cx, gas.state().state(), 300.0, false, 9e-6);
        let layers = layer(&initial, &card, &point, &key, false);
        let specimen = initial.mechanical().specimen();
        let cs = specimen.mass_kg() * CP;
        let r0 = 1e-6 * specimen.string().length_m / specimen.area_m2();
        let limits = feedback(5.0 * cs);
        let heat = FiniteGasStringHeat {
            string_heat_j: cs,
            gas_heat_j: 0.0,
        };
        let dt = 1.0 / f64::from(RATE);
        let w0 = dt * r0; // Unit current and matching initial voltage power.
        let beta = 0.2 / cs;
        let a = 1.0 + beta * heat.string_heat_j;
        for drive in [
            OhmicDrive::Current(1.0),
            OhmicDrive::Current(-1.0),
            OhmicDrive::Voltage(r0),
            OhmicDrive::Voltage(-r0),
        ] {
            // R(T)=R0[1+0.2(T-300)], Cs(T-300)=H+W.
            // Current: W=W0(a+beta W). Voltage: W(a+beta W)=W0.
            let expected = match drive {
                OhmicDrive::Current(_) => w0 * a / (1.0 - w0 * beta),
                OhmicDrive::Voltage(_) => 2.0 * w0 / (a + (a * a + 4.0 * beta * w0).sqrt()),
            };
            let mut live = initial.clone();
            let mut next_gas = gas.clone();
            let result = live
                .step_with_implicit_drive_and_isochoric_gas_conduction(
                    cx,
                    0,
                    &[0.0],
                    drive,
                    &mut next_gas,
                    heat,
                    solve(),
                    controls(),
                    &layers,
                    limits,
                )
                .unwrap();
            let frame = &result.electrical;
            assert!((frame.joule_heat_j - expected).abs() <= 4.0 * limits.source_tolerance_j);
            assert!((result.source_work_j - expected).abs() <= 4.0 * limits.source_tolerance_j);
            assert!(
                (live.thermal().temperature_k() - 300.0 - (expected + heat.string_heat_j) / cs)
                    .abs()
                    < 2e-7
            );
            assert!(result.source_residual_j.abs() <= limits.source_tolerance_j);
            assert_eq!(
                result.source_residual_j,
                frame.joule_heat_j - result.source_work_j
            );
            close(
                result.source_work_j,
                dt * frame.final_current_a.powi(2) * frame.conductor_after.resistance_ohm(),
            );
            let stored = cs * (live.thermal().temperature_k() - 300.0);
            assert!(
                (result.source_work_j + heat.string_heat_j - stored).abs()
                    <= result.energy_tolerance_j
            );
            assert!(result.energy_balance_residual_j.abs() <= result.energy_tolerance_j);
            assert_eq!(live.mechanical().accepted_samples(), 1);
            assert_eq!(live.mechanical().epoch(), 1);
            assert!(result.trials > 2);
            let mut explicit = initial.clone();
            let mut explicit_gas = gas.clone();
            let frozen = explicit
                .step_with_drive_and_isochoric_gas_conduction(
                    cx,
                    0,
                    &[0.0],
                    drive,
                    &mut explicit_gas,
                    heat,
                    solve(),
                    controls(),
                    &layers,
                )
                .unwrap();
            close(frozen.joule_heat_j, w0);
            assert!((frame.joule_heat_j - frozen.joule_heat_j).abs() > 1e-8);
        }
    });
}

#[test]
fn g1_electrical_feedback_and_conduction_close_the_independent_two_capacity_endpoint() {
    with_string_cx(|cx, _| {
        let mut gas = parcel(0.028, 100000.0);
        let (mut live, card, point, key) =
            conducting_string(cx, gas.state().state(), 310.0, false, 9e-6);
        let layers = layer(&live, &card, &point, &key, true);
        let cs = live.mechanical().specimen().mass_kg() * CP;
        let cg = gas.mass_kg() * 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
        let r0 = 1e-6 * live.mechanical().specimen().string().length_m
            / live.mechanical().specimen().area_m2();
        let dt = 1.0 / f64::from(RATE);
        // Independent scalar elimination using the manufactured source curves,
        // rather than either production feedback solver or its reported heat.
        let g = |tg: f64| {
            layers.area_m2
                / (layers.gas_distance_m / (0.025 * (tg + 500.0) / 700.0)
                    + layers.body_distance_m / 20.0)
        };
        let gas_endpoint = |ts: f64| {
            let (mut lo, mut hi) = (300.0, ts);
            for _ in 0..80 {
                let tg = 0.5 * (lo + hi);
                if cg * (tg - 300.0) - dt * g(tg) * (ts - tg) > 0.0 {
                    hi = tg;
                } else {
                    lo = tg;
                }
            }
            0.5 * (lo + hi)
        };
        let (mut lo, mut hi) = (300.0, 330.0);
        for _ in 0..80 {
            let ts = 0.5 * (lo + hi);
            let tg = gas_endpoint(ts);
            let work = dt * r0 * (1.0 + 0.2 * (ts - 300.0));
            let residual = cs * (ts - 310.0) + cg * (tg - 300.0) - work;
            if residual > 0.0 {
                hi = ts;
            } else {
                lo = ts;
            }
        }
        let ts = 0.5 * (lo + hi);
        let tg = gas_endpoint(ts);
        let expected_work = dt * r0 * (1.0 + 0.2 * (ts - 300.0));
        let result = live
            .step_with_implicit_drive_and_isochoric_gas_conduction(
                cx,
                0,
                &[0.0],
                OhmicDrive::Current(1.0),
                &mut gas,
                FiniteGasStringHeat::default(),
                solve(),
                controls(),
                &layers,
                feedback(10.0 * cs),
            )
            .unwrap();
        let exchange = &result.electrical.coupled;
        assert!((live.thermal().temperature_k() - ts).abs() < 2e-6);
        assert!((gas.state().state().temperature - tg).abs() < 2e-8);
        assert!((result.source_work_j - expected_work).abs() < 2e-12);
        assert!((exchange.solid_to_gas_heat_j - dt * g(tg) * (ts - tg)).abs() < 2e-12);
        assert!(exchange.solid_to_gas_heat_j > 0.0);
        assert_eq!(
            exchange.source.query_point(),
            result.electrical.conductor_after.material().query_point()
        );
        assert!(exchange.source_residual_j.abs() <= controls().source_tolerance_j);
        assert!(result.source_residual_j.abs() <= feedback(10.0 * cs).source_tolerance_j);
        assert!(result.energy_balance_residual_j.abs() <= result.energy_tolerance_j);
    });
}

#[test]
fn g3_implicit_source_transports_one_native_inventory_after_convergence() {
    with_string_cx(|cx, _| {
        let gas = parcel(0.028, 100000.0);
        let reservoir =
            ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.004, 200000.0));
        let (initial, card, point, key) =
            conducting_string(cx, gas.state().state(), 310.0, true, 9e-6);
        let layers = layer(&initial, &card, &point, &key, true);
        let limits = feedback(10.0 * initial.mechanical().specimen().mass_kg() * CP);
        let mut closed = initial.clone();
        let mut closed_gas = gas.clone();
        let reference = closed
            .step_with_implicit_drive_and_isochoric_gas_conduction(
                cx,
                0,
                &[0.0],
                OhmicDrive::Current(1.0),
                &mut closed_gas,
                FiniteGasStringHeat::default(),
                solve(),
                controls(),
                &layers,
                limits,
            )
            .unwrap();
        assert!(reference.trials > 2);
        for area_m2 in [0.0, nozzle().area_m2] {
            let geometry = GasNozzle {
                area_m2,
                ..nozzle()
            };
            let (native, expected_flow) = reservoir
                .through_nozzle_to(
                    &ClosedGasParcel::from_isochoric(closed_gas.clone()),
                    geometry,
                    1.0 / f64::from(RATE),
                    flow_solve(),
                    || false,
                )
                .unwrap();
            let (expected_source, expected_gas, expected_transfer) = native.into_parts();
            let mut live = initial.clone();
            let mut local = gas.clone();
            let mut source = reservoir.clone();
            let result = live
                .step_with_implicit_drive_and_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    OhmicDrive::Current(1.0),
                    &mut local,
                    &mut source,
                    FiniteGasStringHeat::default(),
                    solve(),
                    controls(),
                    &layers,
                    geometry,
                    flow_solve(),
                    limits,
                )
                .unwrap();
            let flow = &result.electrical.coupled;
            assert_eq!(flow.flow, expected_flow);
            assert_eq!(flow.transfer, expected_transfer);
            assert_eq!(local, expected_gas.into_isochoric());
            assert_eq!(source, expected_source);
            assert_eq!(result.source_work_j, reference.source_work_j);
            assert_eq!(
                result.electrical.joule_heat_j,
                reference.electrical.joule_heat_j
            );
            assert_eq!(live.mechanical().states(), closed.mechanical().states());
            assert_eq!(live.mechanical().ambient(), local.state().state());
            assert_eq!(live.mechanical().accepted_samples(), 1);
            assert_eq!(live.mechanical().epoch(), 1);
            assert!(result.energy_balance_residual_j.abs() <= result.energy_tolerance_j);
            if area_m2 == 0.0 {
                assert_eq!(
                    flow.thermal.coupled.string,
                    reference.electrical.coupled.coupled.string
                );
                assert_eq!(local, closed_gas);
                assert_eq!(source, reservoir);
            } else {
                assert_eq!(local.state().components().len(), 2);
                assert_ne!(
                    flow.thermal
                        .coupled
                        .string
                        .vibration
                        .acoustic
                        .observer_pressure_pa,
                    reference
                        .electrical
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
fn g4_feedback_work_source_and_cancellation_failures_preserve_three_owners_and_retry() {
    with_string_cx(|cx, gate| {
        let gas = parcel(0.028, 100000.0);
        let reservoir =
            ClosedGasParcel::from_isochoric(parcel_named("introduced-gas", 0.004, 200000.0));
        let (initial, card, point, key) =
            conducting_string(cx, gas.state().state(), 310.0, true, 9e-6);
        let layers = layer(&initial, &card, &point, &key, false);
        let cs = initial.mechanical().specimen().mass_kg() * CP;
        let r0 = 3e-6 * initial.mechanical().specimen().string().length_m
            / initial.mechanical().specimen().area_m2();
        let drive = OhmicDrive::Voltage(r0);
        let limits = feedback(5.0 * cs);
        let mut replay = initial.clone();
        let mut replay_gas = gas.clone();
        let mut replay_source = reservoir.clone();
        let expected = replay
            .step_with_implicit_drive_and_isochoric_gas_nozzle_conduction(
                cx,
                0,
                &[0.0],
                drive,
                &mut replay_gas,
                &mut replay_source,
                FiniteGasStringHeat::default(),
                solve(),
                controls(),
                &layers,
                nozzle(),
                flow_solve(),
                limits,
            )
            .unwrap();
        for (case, bad) in [
            feedback(cs * 1e-6),
            FiniteGasElectricalFeedbackSolve {
                maximum_iterations: 1,
                ..limits
            },
            feedback(50.0 * cs),
        ]
        .into_iter()
        .enumerate()
        {
            let mut live = initial.clone();
            let mut local = gas.clone();
            let mut source = reservoir.clone();
            let refused = live.step_with_implicit_drive_and_isochoric_gas_nozzle_conduction(
                cx,
                0,
                &[0.0],
                drive,
                &mut local,
                &mut source,
                FiniteGasStringHeat::default(),
                solve(),
                controls(),
                &layers,
                nozzle(),
                flow_solve(),
                bad,
            );
            match case {
                0 => assert!(matches!(
                    refused,
                    Err(FiniteGasStringExchangeError::String(
                        FiniteGasStringError::ElectricalFeedback(
                            FiniteGasElectricalFeedbackError::OutsideBracket { .. }
                        )
                    ))
                )),
                1 => assert!(matches!(
                    refused,
                    Err(FiniteGasStringExchangeError::String(
                        FiniteGasStringError::ElectricalFeedback(
                            FiniteGasElectricalFeedbackError::NotConverged { .. }
                        )
                    ))
                )),
                _ => assert!(matches!(
                    refused,
                    Err(FiniteGasStringExchangeError::String(
                        FiniteGasStringError::String(ThermalStringError::Admission(_))
                    ))
                )),
            }
            assert_eq!(local, gas);
            assert_eq!(source, reservoir);
            assert_eq!(live.thermal(), initial.thermal());
            assert_eq!(live.mechanical().states(), initial.mechanical().states());
            assert_eq!(
                live.mechanical().specimen(),
                initial.mechanical().specimen()
            );
            assert_eq!(live.mechanical().ambient(), initial.mechanical().ambient());
            assert_eq!(live.mechanical().accepted_samples(), 0);
            assert_eq!(live.mechanical().epoch(), 0);
            let retry = live
                .step_with_implicit_drive_and_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    drive,
                    &mut local,
                    &mut source,
                    FiniteGasStringHeat::default(),
                    solve(),
                    controls(),
                    &layers,
                    nozzle(),
                    flow_solve(),
                    limits,
                )
                .unwrap();
            assert_eq!(
                retry.electrical.coupled.thermal.coupled.string,
                expected.electrical.coupled.thermal.coupled.string
            );
            assert_eq!(retry.source_work_j, expected.source_work_j);
            assert_eq!(local, replay_gas);
            assert_eq!(source, replay_source);
        }
        let mut cancelled = initial.clone();
        let mut cancelled_gas = gas.clone();
        let mut cancelled_source = reservoir.clone();
        gate.request();
        assert!(
            cancelled
                .step_with_implicit_drive_and_isochoric_gas_nozzle_conduction(
                    cx,
                    0,
                    &[0.0],
                    drive,
                    &mut cancelled_gas,
                    &mut cancelled_source,
                    FiniteGasStringHeat::default(),
                    solve(),
                    controls(),
                    &layers,
                    nozzle(),
                    flow_solve(),
                    limits
                )
                .is_err()
        );
        assert_eq!(cancelled_gas, gas);
        assert_eq!(cancelled_source, reservoir);
        assert_eq!(
            cancelled.mechanical().states(),
            initial.mechanical().states()
        );
        assert_eq!(cancelled.thermal(), initial.thermal());
        assert_eq!(cancelled.mechanical().accepted_samples(), 0);
        assert_eq!(cancelled.mechanical().epoch(), 0);
    });
}
