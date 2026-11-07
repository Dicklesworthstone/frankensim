//! Held speed commands change native flow, equipment consumption and inventory.
use super::*;
use fs_airflow::SourceProvenance;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{
    PumpOperatingError, speed::{PumpSpeedCurve, PumpSpeedMap},
};

fn speed_map() -> PumpSpeedMap {
    let curves = [
        (1000.0, 16.0, 0.2, 180000.0, 0.55, 180000.0, (1.9e-6, 3e-6)),
        (3000.0, 24.0, 0.4, 280000.0, 0.90, 320000.0, (1e-6, 2.5e-6)),
    ]
    .into_iter()
    .map(
        |(
            speed_rpm,
            head,
            shaft_idle,
            shaft_slope,
            electrical_idle,
            electrical_slope,
            viscosity,
        )| {
            let head = PumpHeadCurve::try_new(
                [0.0, 4e-5]
                    .map(|q| PumpHeadPoint {
                        flow: VolumetricFlowRate::new(q),
                        head: Length::new(head),
                    })
                    .to_vec(),
                10.0,
                viscosity,
                SourceProvenance::new(
                    "synthetic held-speed coolant head",
                    format!("coolant-head-{speed_rpm}"),
                ),
            )
            .unwrap();
            let power = PumpDriveCurve::try_new(
                [0.0, 4e-5]
                    .map(|q| PumpPowerPoint {
                        flow: VolumetricFlowRate::new(q),
                        shaft_input: Power::new(shaft_idle + shaft_slope * q),
                        electrical_input: Power::new(electrical_idle + electrical_slope * q),
                    })
                    .to_vec(),
                (990.0, 1010.0),
                viscosity,
                PumpLossRouting {
                    pump_to_reservoir: 1.0,
                    drive_to_reservoir: 1.0,
                },
                SourceProvenance::new(
                    "synthetic independent input powers, not equipment measurements",
                    format!("coolant-input-{speed_rpm}"),
                ),
            )
            .unwrap();
            PumpSpeedCurve {
                speed_rpm,
                curve: head.with_drive_curve(power),
            }
        },
    )
    .collect();
    PumpSpeedMap::try_new(
        curves,
        SourceProvenance::new(
            "synthetic same-flow interpolation between supplied speed data",
            "coolant-speed-map",
        ),
    )
    .unwrap()
}

fn with_heated_network<R>(
    text: &str,
    run: impl FnOnce(&Cx<'_>, &ParallelLiquidNetwork<'_, '_, '_>) -> R,
) -> R {
    with_network(text, |cx, original| {
        // Four thousand times the fixture's resistance gives a 160--240 kPa
        // circulating loop with hydraulic and equipment powers both in watts.
        let resistance = MaterialFlowResistance::try_new(
            2e13,
            0.0,
            Area::new(core::f64::consts::PI * 0.01 * 0.01 / 4.0),
            Length::new(0.01),
            (0.0, 2300.0),
        )
        .unwrap();
        let losses =
            ["first", "last"].map(|region| [LiquidResistanceSection { region, resistance }]);
        let branches: Vec<_> = original
            .branches()
            .iter()
            .enumerate()
            .map(|(i, b)| ParallelLiquidBranch {
                resistances: &losses[i],
                pressure_budget: PressureDrivenBudget {
                    maximum_neglected_hydraulic_power_w: 0.0,
                    ..b.pressure_budget
                },
                ..*b
            })
            .collect();
        let network = ParallelLiquidNetwork::try_new(cx, &branches)
            .unwrap()
            .with_hydraulic_dissipation();
        run(cx, &network)
    })
}

fn controls() -> PumpCouplingControls {
    let mut result = coupling();
    result.network.work.maximum_neglected_hydraulic_power_w = 0.0;
    result.operating_point.pressure_bracket_pa = (120000.0, 280000.0);
    result
}

#[test]
fn changing_calibrated_speed_advances_original_flow_and_coolant_energy_history() {
    with_heated_network(FIXED, |cx, network| {
        let map = speed_map();
        let mut tank = network
            .reservoir(cx, 0.1, Temperature::new(300.0), 0.0)
            .unwrap();
        let initial_enthalpy = tank.state().specific_enthalpy_j_kg();
        let step = ReservoirStepControls {
            maximum_neglected_hydraulic_energy_j: 0.0,
            ..step()
        };
        let mut supplied_j = 0.0;
        // Independently: Q=dp/2e10, head=[16,20,24] m at rho=1000,
        // g=10. Each input curve supplies its own watts at the resulting Q.
        for (end, rpm, pressure, flow, electrical) in [
            (10.0, 1000.0, 160000.0, 8e-6, 1.99),
            (20.0, 2000.0, 200000.0, 1e-5, 3.225),
            (30.0, 3000.0, 240000.0, 1.2e-5, 4.74),
        ] {
            let curve = map.at_speed(rpm).unwrap();
            let result = network
                .advance_reservoir_with_pump_coupled_iqn(
                    cx,
                    &mut tank,
                    end,
                    0.0,
                    &[300.0; 2],
                    &curve,
                    &controls(),
                    step,
                    |_, rows| Ok(body(rows, [0.0; 2])),
                )
                .unwrap();
            let supply = result.circuit.supply;
            let drive = supply.drive_power.unwrap();
            assert_eq!(supply.speed_rpm, Some(rpm));
            assert!((result.circuit.network_pressure_pa - pressure).abs() < 1e-4);
            assert!((supply.inlet_volume_flow_m3_s - flow).abs() < 1e-12);
            assert!((drive.electrical_input_w - electrical).abs() < 1e-7);
            assert!((supply.hydraulic_power_w - pressure * flow).abs() < 1e-7);
            assert_eq!(drive.surroundings_heat_w, 0.0);
            assert_eq!(result.neglected_hydraulic_energy_j, 0.0);
            assert!(
                (result.drive_heat_j + result.deposited_hydraulic_energy_j - 10.0 * electrical)
                    .abs()
                    < 2e-6
            );
            assert!((result.stored_enthalpy_change_j - 10.0 * electrical).abs() < 2e-6);
            assert!(result.energy_residual_j.abs() <= step.energy_absolute_tolerance_j);
            supplied_j += 10.0 * electrical;
            assert!(
                (0.1 * (tank.state().specific_enthalpy_j_kg() - initial_enthalpy) - supplied_j)
                    .abs()
                    < 5e-6
            );
            assert_eq!(tank.time_s(), end);
            assert_eq!(tank.state(), result.final_state);
        }
        // 99.55 J over three actual intervals, not three resets of the inlet.
        assert!((tank.state().temperature_k() - 300.49775).abs() < 1e-7);
    });
}

#[test]
fn selected_speed_keeps_transport_feedback_and_unsupported_trials_roll_back() {
    // Native hydraulic deposition requires constant density. Retain that
    // physical assumption while varying both viscosity and heat capacity.
    let transport = VARIABLE
        .replace("250 1050", "250 1000")
        .replace("450 850", "450 1000");
    with_heated_network(&transport, |cx, network| {
        let map = speed_map();
        let curve = map.at_speed(2000.0).unwrap();
        let mut tank = network
            .reservoir(cx, 0.1, Temperature::new(300.0), 0.0)
            .unwrap();
        let step = ReservoirStepControls {
            temperature_bracket_k: (299.0, 301.0),
            maximum_neglected_hydraulic_energy_j: 0.0,
            ..step()
        };
        let result = network
            .advance_reservoir_with_pump_coupled_iqn(
                cx,
                &mut tank,
                10.0,
                0.0,
                &[300.0; 2],
                &curve,
                &controls(),
                step,
                |_, rows| Ok(body(rows, [0.0; 2])),
            )
            .unwrap();
        let supply = result.circuit.supply;
        let drive = supply.drive_power.unwrap();
        assert!(tank.state().temperature_k() > 300.0);
        assert!(supply.inlet_volume_flow_m3_s > 1e-5);
        assert!(
            (drive.shaft_input_w - (0.3 + 230000.0 * supply.inlet_volume_flow_m3_s)).abs() < 1e-9
        );
        assert!(
            (drive.electrical_input_w - (0.725 + 250000.0 * supply.inlet_volume_flow_m3_s)).abs()
                < 1e-9
        );
        assert!((result.stored_enthalpy_change_j - 10.0 * drive.electrical_input_w).abs() < 3e-6);

        // The source liquid is valid at 315 K, but the 2000 rpm command is
        // restricted to the intersection of the two measured viscosity spans.
        // Do not silently use the wider high-speed calibration or keep heating.
        let saved = tank.clone();
        let refused = network.advance_reservoir_with_pump_coupled_iqn(
            cx,
            &mut tank,
            20.0,
            0.0,
            &[315.0; 2],
            &curve,
            &controls(),
            ReservoirStepControls {
                temperature_bracket_k: (315.0, 316.0),
                ..step
            },
            |_, rows| Ok(body(rows, [0.0; 2])),
        );
        assert!(matches!(
            refused,
            Err(RecirculationError::Physics(
                PumpOperatingError::OutsideCurve { .. }
            ))
        ));
        assert_eq!(tank, saved);
    });
}
