//! Source draw limits on original calibrated pump/material/network solves.
use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::ParallelPressureError;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::source::{
    PumpPowerLimitControls, PumpPowerLimitError,
};
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::speed::{
    PumpSpeedCurve, PumpSpeedMap,
};

fn slice(speed_rpm: f64, head: f64, shaft: f64, idle: f64, slope: f64) -> PumpSpeedCurve {
    let head = PumpHeadCurve::try_new(
        [0.0, 1e-5]
            .map(|q| PumpHeadPoint {
                flow: VolumetricFlowRate::new(q),
                head: Length::new(head),
            })
            .to_vec(),
        10.0,
        (1e-6, 4e-6),
        SourceProvenance::new(
            "synthetic constant-head calibration",
            format!("source-head-{speed_rpm}"),
        ),
    )
    .unwrap();
    let input = PumpDriveCurve::try_new(
        [0.0, 1e-5]
            .map(|q| PumpPowerPoint {
                flow: VolumetricFlowRate::new(q),
                shaft_input: Power::new(shaft),
                electrical_input: Power::new(idle + slope * q),
            })
            .to_vec(),
        (990.0, 1010.0),
        (1e-6, 4e-6),
        PumpLossRouting {
            pump_to_reservoir: 0.25,
            drive_to_reservoir: 0.75,
        },
        SourceProvenance::new(
            "independently authored flow-dependent input",
            format!("source-drive-{speed_rpm}"),
        ),
    )
    .unwrap();
    PumpSpeedCurve {
        speed_rpm,
        curve: head.with_drive_curve(input),
    }
}

fn family() -> PumpSpeedMap {
    PumpSpeedMap::try_new(
        vec![
            slice(1000.0, 0.003, 0.002, 0.004, 100.0),
            slice(3000.0, 0.006, 0.006, 0.010, 200.0),
        ],
        SourceProvenance::new(
            "synthetic same-flow linear speed interpolation",
            "source-limit-family",
        ),
    )
    .unwrap()
}

fn source_controls() -> PumpPowerLimitControls {
    PumpPowerLimitControls {
        lower_speed_rpm: 1000.0,
        // At 2000 rpm: dp=45 Pa, Q=4.5e-6 m3/s, Pe=.007+150*Q.
        // Small explicit shortfall prevents floating-point overshoot acceptance.
        maximum_electrical_input_w: 0.007675 + 1e-9,
        absolute_power_tolerance_w: 2e-9,
        relative_power_tolerance: 0.0,
        max_speed_evaluations: 48,
        max_total_network_evaluations: 2048,
    }
}

fn source_work() -> ParallelPressureBudget {
    ParallelPressureBudget {
        max_flow_evaluations: 500000,
        ..work()
    }
}

fn source_coupling() -> PumpCouplingControls {
    let mut controls = coupled_controls();
    controls.network.work = source_work();
    controls.max_total_solid_solves = 20000;
    controls
}

#[test]
fn feasible_command_is_preserved_and_cutback_closes_actual_flow_power() {
    let map = family();
    let unchanged = map.clone();
    with_network(FIXED, |cx, network| {
        let commanded = map.at_speed(3000.0).unwrap();
        let reference = network
            .solve_with_pump(cx, &[300.0; 2], &commanded, local(), source_work(), root())
            .unwrap();
        let accepted = network
            .solve_with_power_limited_pump(
                cx,
                &[300.0; 2],
                &map,
                3000.0,
                PumpPowerLimitControls {
                    maximum_electrical_input_w: 0.02,
                    ..source_controls()
                },
                local(),
                source_work(),
                root(),
            )
            .unwrap();
        assert!(!accepted.curtailed);
        assert_eq!(accepted.speed_evaluations, 1);
        assert_eq!(accepted.curve, commanded);
        assert_eq!(accepted.supply, reference.supply);
        assert_eq!(accepted.state.mixed_outlet, reference.state.mixed_outlet);
        assert_eq!(
            accepted.total_network_evaluations,
            reference.network_evaluations
        );
        assert_eq!(
            accepted.total_flow_evaluations,
            reference.total_flow_evaluations
        );

        let limited = network
            .solve_with_power_limited_pump(
                cx,
                &[300.0; 2],
                &map,
                3000.0,
                source_controls(),
                local(),
                source_work(),
                root(),
            )
            .unwrap();
        assert!(limited.curtailed);
        assert_eq!(limited.speed_evaluations, 3);
        assert_eq!(limited.supply.speed_rpm, Some(2000.0));
        assert_eq!(limited.curve.speed_rpm(), limited.supply.speed_rpm);
        close(limited.network_pressure_pa, 45.0, 1e-6);
        close(limited.supply.inlet_volume_flow_m3_s, 4.5e-6, 1e-12);
        assert!(limited.pressure_residual_pa.abs() < 1e-6);
        let power = limited.supply.drive_power.unwrap();
        close(power.electrical_input_w, 0.007675, 2e-10);
        assert!(power.electrical_input_w <= source_controls().maximum_electrical_input_w);
        assert!(limited.electrical_power_margin_w <= source_controls().absolute_power_tolerance_w);
        close(
            power.electrical_input_w,
            limited.supply.hydraulic_power_w + power.reservoir_heat_w + power.surroundings_heat_w,
            1e-15,
        );
        close(
            power.reservoir_heat_w,
            0.25 * power.pump_loss_w + 0.75 * power.drive_loss_w,
            1e-15,
        );
        assert!(limited.total_network_evaluations > reference.network_evaluations);
        assert!(limited.total_flow_evaluations > limited.state.flow_evaluations);
        assert_eq!(limited.total_solid_solves, 0);
        // Replaying the owned selected curve gives the same original state.
        let replay = network
            .solve_with_pump(
                cx,
                &[300.0; 2],
                &limited.curve,
                local(),
                source_work(),
                root(),
            )
            .unwrap();
        assert_eq!(replay.supply, limited.supply);
        assert_eq!(replay.state.mixed_outlet, limited.state.mixed_outlet);
    });
    assert_eq!(map, unchanged);
}

#[test]
fn material_feedback_changes_curtailed_speed_and_uses_actual_volume_flow() {
    let map = family();
    with_network(VARIABLE, |cx, network| {
        let cold = network
            .solve_with_power_limited_pump(
                cx,
                &[300.0; 2],
                &map,
                3000.0,
                source_controls(),
                local(),
                source_work(),
                root(),
            )
            .unwrap();
        let hot = network
            .solve_with_power_limited_pump(
                cx,
                &[400.0, 300.0],
                &map,
                3000.0,
                source_controls(),
                local(),
                source_work(),
                root(),
            )
            .unwrap();
        assert!(hot.supply.speed_rpm.unwrap() < cold.supply.speed_rpm.unwrap());
        assert!(
            hot.state.branches[0].solution.mass_flow_kg_per_s
                > hot.state.branches[1].solution.mass_flow_kg_per_s
        );
        for point in [&cold, &hot] {
            let fraction = (point.supply.speed_rpm.unwrap() - 1000.0) / 2000.0;
            let actual_q =
                point.state.total_mass_flow_kg_per_s / point.state.inlet_material.density_kg_m3();
            let independent_input =
                0.004 + 0.006 * fraction + (100.0 + 100.0 * fraction) * actual_q;
            close(
                point.supply.drive_power.unwrap().electrical_input_w,
                independent_input,
                1e-15,
            );
            assert!(independent_input <= source_controls().maximum_electrical_input_w);
            assert!(
                point.electrical_power_margin_w <= source_controls().absolute_power_tolerance_w
            );
            assert!(point.pressure_residual_pa.abs() < 1e-6);
            assert!(point.total_flow_evaluations <= source_work().max_flow_evaluations);
        }
    });
}

#[test]
fn nonmonotone_calibration_needs_only_the_explicit_local_bracket() {
    with_network(FIXED, |cx, network| {
        let provenance = family().source().clone();
        // A high intermediate draw is irrelevant when the requested command fits.
        let peak = PumpSpeedMap::try_new(
            vec![
                slice(1000.0, 0.003, 0.002, 0.004, 0.0),
                slice(2000.0, 0.0045, 0.002, 0.012, 0.0),
                slice(3000.0, 0.006, 0.002, 0.006, 0.0),
            ],
            provenance.clone(),
        )
        .unwrap();
        let point = network
            .solve_with_power_limited_pump(
                cx,
                &[300.0; 2],
                &peak,
                3000.0,
                source_controls(),
                local(),
                source_work(),
                root(),
            )
            .unwrap();
        assert_eq!(point.speed_evaluations, 1);
        assert_eq!(point.supply.speed_rpm, Some(3000.0));
        assert!(!point.curtailed);

        // Feasible interior calibration does not make two infeasible endpoints
        // an admitted bracket, nor justify a globally optimal-speed claim.
        let valley = PumpSpeedMap::try_new(
            vec![
                slice(1000.0, 0.003, 0.002, 0.010, 0.0),
                slice(2000.0, 0.0045, 0.002, 0.004, 0.0),
                slice(3000.0, 0.006, 0.002, 0.012, 0.0),
            ],
            provenance,
        )
        .unwrap();
        assert!(matches!(
            network.solve_with_power_limited_pump(
                cx,
                &[300.0; 2],
                &valley,
                3000.0,
                source_controls(),
                local(),
                source_work(),
                root()
            ),
            Err(PumpPowerLimitError::NotBracketed { .. })
        ));
    });
}

#[test]
fn speed_network_and_flow_caps_span_all_complete_trials() {
    let map = family();
    with_network(FIXED, |cx, network| {
        let commanded = map.at_speed(3000.0).unwrap();
        let reference = network
            .solve_with_pump(cx, &[300.0; 2], &commanded, local(), source_work(), root())
            .unwrap();
        assert!(matches!(
            network.solve_with_power_limited_pump(
                cx,
                &[300.0; 2],
                &map,
                3000.0,
                PumpPowerLimitControls {
                    max_speed_evaluations: 2,
                    ..source_controls()
                },
                local(),
                source_work(),
                root()
            ),
            Err(PumpPowerLimitError::NotConverged {
                speed_evaluations: 2,
                ..
            })
        ));
        let network_cap = reference.network_evaluations;
        assert!(
            matches!(network.solve_with_power_limited_pump(cx, &[300.0; 2], &map,
            3000.0, PumpPowerLimitControls { max_total_network_evaluations: network_cap,
                ..source_controls() }, local(), source_work(), root()),
            Err(PumpPowerLimitError::NetworkBudget { performed, maximum })
                if performed == network_cap && maximum == network_cap)
        );
        let flow_cap = reference.total_flow_evaluations;
        assert!(
            matches!(network.solve_with_power_limited_pump(cx, &[300.0; 2], &map,
            3000.0, source_controls(), local(),
            ParallelPressureBudget { max_flow_evaluations: flow_cap, ..source_work() }, root()),
            Err(PumpPowerLimitError::Operating(PumpOperatingError::Network(
                ParallelPressureError::FlowBudget { performed, maximum })))
                if performed == flow_cap && maximum == flow_cap)
        );
    });
}

#[test]
fn common_solid_is_resolved_at_each_speed_and_total_callback_limit_is_not_reset() {
    let map = family();
    with_network(FIXED, |cx, network| {
        let mut callbacks = 0;
        let point = network
            .solve_with_power_limited_pump_coupled_iqn(
                cx,
                &[300.0; 2],
                &map,
                3000.0,
                source_controls(),
                &source_coupling(),
                |_, rows| {
                    callbacks += 1;
                    Ok(common_solid(rows, [12.0, 8.0]))
                },
            )
            .unwrap();
        assert_eq!(point.supply.speed_rpm, Some(2000.0));
        assert_eq!(point.speed_evaluations, 3);
        assert_eq!(callbacks, point.total_solid_solves);
        assert!(callbacks > point.state.iterations);
        assert!(point.total_flow_evaluations <= source_work().max_flow_evaluations);
        close(point.state.fluid.total_heat_rate_w, 20.0, 2e-5);
        close(
            point.supply.drive_power.unwrap().electrical_input_w,
            0.007675,
            2e-10,
        );
        let g = 4.5 * (-(-0.549_f64 / 4.5).exp_m1());
        let mean = 300.0 + 10.0 / g;
        let difference = 4.0 / (g + 2.0);
        close(point.state.solid.state[0], mean + 0.5 * difference, 1e-5);
        close(point.state.solid.state[1], mean - 0.5 * difference, 1e-5);

        let commanded = map.at_speed(3000.0).unwrap();
        let reference = network
            .solve_with_pump_coupled_iqn(
                cx,
                &[300.0; 2],
                &commanded,
                &source_coupling(),
                |_, rows| Ok(common_solid(rows, [12.0, 8.0])),
            )
            .unwrap();
        let cap = reference.total_solid_solves;
        let mut callbacks = 0;
        let error = network
            .solve_with_power_limited_pump_coupled_iqn(
                cx,
                &[300.0; 2],
                &map,
                3000.0,
                source_controls(),
                &PumpCouplingControls {
                    max_total_solid_solves: cap,
                    ..source_coupling()
                },
                |_, rows| {
                    callbacks += 1;
                    Ok(common_solid(rows, [12.0, 8.0]))
                },
            )
            .unwrap_err();
        assert!(
            matches!(error, PumpPowerLimitError::Operating(PumpOperatingError::SolidBudget {
            performed, maximum }) if performed == cap && maximum == cap)
        );
        assert_eq!(callbacks, cap);
    });
}

#[test]
fn source_and_power_refusals_preserve_native_gates_and_original_map() {
    let map = family();
    let unchanged = map.clone();
    with_network(FIXED, |cx, network| {
        for controls in [
            PumpPowerLimitControls {
                maximum_electrical_input_w: f64::NAN,
                ..source_controls()
            },
            PumpPowerLimitControls {
                lower_speed_rpm: 999.0,
                ..source_controls()
            },
            PumpPowerLimitControls {
                max_total_network_evaluations: 1,
                ..source_controls()
            },
        ] {
            assert!(matches!(
                network.solve_with_power_limited_pump_coupled_iqn::<(), _>(
                    cx,
                    &[300.0; 2],
                    &map,
                    3000.0,
                    controls,
                    &source_coupling(),
                    |_, _| { panic!("invalid source controls must precede common-solid work") }
                ),
                Err(PumpPowerLimitError::Invalid(_))
            ));
        }
        let mut curves = map.curves().to_vec();
        for entry in &mut curves {
            entry.curve = entry.curve.clone().with_drive_curve(drive(1e-5, 0.02));
        }
        let impossible = PumpSpeedMap::try_new(curves, map.source().clone()).unwrap();
        assert!(matches!(
            network.solve_with_power_limited_pump(
                cx,
                &[300.0; 2],
                &impossible,
                3000.0,
                source_controls(),
                local(),
                source_work(),
                root()
            ),
            Err(PumpPowerLimitError::Operating(PumpOperatingError::Invalid(
                _
            )))
        ));

        let mut curves = map.curves().to_vec();
        for entry in &mut curves {
            let unsupported = PumpHeadCurve::try_new(
                entry.curve.points().to_vec(),
                entry.curve.gravity_m_per_s2(),
                (3e-6, 4e-6),
                entry.curve.source().clone(),
            )
            .unwrap()
            .with_drive_curve(entry.curve.drive_curve().unwrap().clone());
            entry.curve = unsupported;
        }
        let unsupported = PumpSpeedMap::try_new(curves, map.source().clone()).unwrap();
        assert!(matches!(
            network.solve_with_power_limited_pump_coupled_iqn::<(), _>(
                cx,
                &[300.0; 2],
                &unsupported,
                3000.0,
                source_controls(),
                &source_coupling(),
                |_, _| { panic!("original source refusal must precede common-solid work") }
            ),
            Err(PumpPowerLimitError::Operating(
                PumpOperatingError::OutsideCurve { .. }
            ))
        ));
    });
    assert_eq!(map, unchanged);
}
