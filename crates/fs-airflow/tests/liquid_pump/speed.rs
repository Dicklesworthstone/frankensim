//! Calibrated speed selection through the original source and network owners.
use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::speed::{
    PumpSpeedCurve, PumpSpeedMap,
};

fn calibrated(
    speed_rpm: f64,
    heads: &[(f64, f64)],
    inputs: &[(f64, f64, f64)],
    density: (f64, f64),
    viscosity: (f64, f64),
) -> PumpSpeedCurve {
    let curve = PumpHeadCurve::try_new(
        heads
            .iter()
            .map(|&(q, h)| PumpHeadPoint {
                flow: VolumetricFlowRate::new(q),
                head: Length::new(h),
            })
            .collect(),
        10.0,
        viscosity,
        SourceProvenance::new(
            "synthetic speed calibration, not equipment measurements",
            format!("head-{speed_rpm}"),
        ),
    )
    .unwrap();
    let input = PumpDriveCurve::try_new(
        inputs
            .iter()
            .map(|&(q, shaft, electrical)| PumpPowerPoint {
                flow: VolumetricFlowRate::new(q),
                shaft_input: Power::new(shaft),
                electrical_input: Power::new(electrical),
            })
            .collect(),
        density,
        viscosity,
        PumpLossRouting {
            pump_to_reservoir: 0.25,
            drive_to_reservoir: 0.75,
        },
        SourceProvenance::new(
            "synthetic independently supplied drive data",
            format!("drive-{speed_rpm}"),
        ),
    )
    .unwrap();
    PumpSpeedCurve {
        speed_rpm,
        curve: curve.with_drive_curve(input),
    }
}

fn family() -> PumpSpeedMap {
    PumpSpeedMap::try_new(
        vec![
            calibrated(
                1000.0,
                &[(0.0, 0.008), (4e-6, 0.006), (1e-5, 0.002)],
                &[
                    (0.0, 0.002, 0.004),
                    (3e-6, 0.003, 0.005),
                    (1e-5, 0.006, 0.009),
                ],
                (990.0, 1050.0),
                (1e-6, 4e-6),
            ),
            calibrated(
                3000.0,
                &[(0.0, 0.016), (6e-6, 0.010), (1e-5, 0.004)],
                &[
                    (0.0, 0.006, 0.010),
                    (7e-6, 0.009, 0.015),
                    (1e-5, 0.012, 0.020),
                ],
                (1000.0, 1020.0),
                (1.5e-6, 3e-6),
            ),
        ],
        SourceProvenance::new(
            "synthetic linear interpolation at fixed inlet flow",
            "speed-family",
        ),
    )
    .unwrap()
}

#[test]
fn same_flow_speed_interpolation_preserves_head_and_independent_power_knots() {
    let map = family();
    let selected = map.at_speed(2000.0).unwrap();
    assert_eq!(selected.speed_rpm(), Some(2000.0));
    assert_eq!(
        selected
            .points()
            .iter()
            .map(|p| p.flow.value())
            .collect::<Vec<_>>(),
        [0.0, 4e-6, 6e-6, 1e-5]
    );
    let input = selected.drive_curve().unwrap();
    assert_eq!(
        input
            .points()
            .iter()
            .map(|p| p.flow.value())
            .collect::<Vec<_>>(),
        [0.0, 3e-6, 7e-6, 1e-5]
    );
    assert_eq!(input.density_domain_kg_m3(), (1000.0, 1020.0));
    assert_eq!(input.viscosity_domain_m2_s(), (1.5e-6, 3e-6));
    assert_eq!(selected.kinematic_viscosity_domain_m2_s(), (1.5e-6, 3e-6));
    assert!(selected.source().citation.contains("2000 rpm"));
    assert!(selected.source().citation.contains("head-1000"));
    assert!(input.source().citation.contains("drive-3000"));
    with_network(FIXED, |cx, network| {
        let state = network
            .solve_pressure_driven(cx, &[300.0; 2], Pressure::new(50.0), local(), work())
            .unwrap();
        for q in [0.0, 1.5e-6, 4e-6, 5e-6, 6.5e-6, 9.5e-6, 1e-5] {
            let low = map.curves()[0]
                .curve
                .evaluate(&state.inlet_material, VolumetricFlowRate::new(q))
                .unwrap();
            let high = map.curves()[1]
                .curve
                .evaluate(&state.inlet_material, VolumetricFlowRate::new(q))
                .unwrap();
            let actual = selected
                .evaluate(&state.inlet_material, VolumetricFlowRate::new(q))
                .unwrap();
            close(actual.head_m, 0.5 * (low.head_m + high.head_m), 1e-15);
            let a = low.drive_power.unwrap();
            let b = high.drive_power.unwrap();
            let power = actual.drive_power.unwrap();
            close(
                power.shaft_input_w,
                0.5 * (a.shaft_input_w + b.shaft_input_w),
                1e-15,
            );
            close(
                power.electrical_input_w,
                0.5 * (a.electrical_input_w + b.electrical_input_w),
                1e-15,
            );
            close(
                power.electrical_input_w,
                actual.hydraulic_power_w + power.reservoir_heat_w + power.surroundings_heat_w,
                1e-15,
            );
            assert_eq!(actual.speed_rpm, Some(2000.0));
        }
    });
    // Calibration points retain original samples/domains rather than borrowing
    // a neighbouring slice's narrower support or modifying measured values.
    for original in map.curves() {
        let exact = map.at_speed(original.speed_rpm).unwrap();
        assert_eq!(exact.points(), original.curve.points());
        assert_eq!(exact.drive_curve(), original.curve.drive_curve());
        assert_eq!(exact.source(), original.curve.source());
    }
    assert_eq!(pump(false).speed_rpm(), None);
}

#[test]
fn calibrated_speed_changes_the_actual_network_operating_point_and_input_power() {
    let map = PumpSpeedMap::try_new(
        vec![
            calibrated(
                1000.0,
                &[(0.0, 0.006), (1e-5, 0.0)],
                &[(0.0, 0.002, 0.004), (1e-5, 0.004, 0.008)],
                (990.0, 1010.0),
                (1e-6, 4e-6),
            ),
            calibrated(
                3000.0,
                &[(0.0, 0.012), (1e-5, 0.002)],
                &[(0.0, 0.006, 0.012), (1e-5, 0.010, 0.022)],
                (990.0, 1010.0),
                (1e-6, 4e-6),
            ),
        ],
        SourceProvenance::new(
            "synthetic speed-dependent pump intersections",
            "speed-network",
        ),
    )
    .unwrap();
    with_network(FIXED, |cx, network| {
        // Original network: dp=1e7*Q. Independently intersect each calibrated
        // linear head curve, including the interpolated 2000 rpm command.
        for (rpm, pressure, electrical) in [
            (1000.0, 37.5, 0.0055),
            (2000.0, 50.0, 0.0115),
            (3000.0, 60.0, 0.018),
        ] {
            let curve = map.at_speed(rpm).unwrap();
            let result = network
                .solve_with_pump(cx, &[300.0; 2], &curve, local(), work(), root())
                .unwrap();
            close(result.network_pressure_pa, pressure, 1e-6);
            close(
                result.state.inlet_volume_flow_m3_per_s,
                pressure / 1e7,
                1e-12,
            );
            close(
                result.supply.drive_power.unwrap().electrical_input_w,
                electrical,
                1e-10,
            );
            close(
                result.supply.hydraulic_power_w,
                pressure * pressure / 1e7,
                1e-10,
            );
            assert_eq!(result.supply.speed_rpm, Some(rpm));
            assert_eq!(result.curve.speed_rpm(), Some(rpm));
            assert!(result.pressure_residual_pa.abs() < 1e-6);
        }
    });
}

#[test]
fn interpolation_retains_only_common_flow_support_and_refuses_extrapolation() {
    let mut curves = family().curves().to_vec();
    curves[1] = calibrated(
        3000.0,
        &[(1e-6, 0.014), (6e-6, 0.010), (9e-6, 0.004)],
        &[
            (2e-6, 0.006, 0.010),
            (7e-6, 0.009, 0.015),
            (8e-6, 0.012, 0.020),
        ],
        (1000.0, 1020.0),
        (1.5e-6, 3e-6),
    );
    let map = PumpSpeedMap::try_new(curves, family().source().clone()).unwrap();
    let selected = map.at_speed(2000.0).unwrap();
    assert_eq!(selected.points().first().unwrap().flow.value(), 2e-6);
    assert_eq!(selected.points().last().unwrap().flow.value(), 8e-6);
    assert_eq!(
        selected
            .drive_curve()
            .unwrap()
            .points()
            .first()
            .unwrap()
            .flow
            .value(),
        2e-6
    );
    assert_eq!(
        selected
            .drive_curve()
            .unwrap()
            .points()
            .last()
            .unwrap()
            .flow
            .value(),
        8e-6
    );
    for speed in [-1.0, 999.0, 3001.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            map.at_speed(speed),
            Err(PumpOperatingError::OutsideCurve { .. })
        ));
    }
    with_network(FIXED, |cx, network| {
        let state = network
            .solve_pressure_driven(cx, &[300.0; 2], Pressure::new(50.0), local(), work())
            .unwrap();
        for q in [0.0, 1e-5] {
            assert!(
                selected
                    .evaluate(&state.inlet_material, VolumetricFlowRate::new(q))
                    .is_err()
            );
            assert!(
                map.at_speed(1000.0)
                    .unwrap()
                    .evaluate(&state.inlet_material, VolumetricFlowRate::new(q))
                    .is_ok()
            );
        }
    });
}

#[test]
fn speed_families_refuse_incompatible_sources_and_keep_actual_power_feasibility() {
    let map = family();
    let mut reversed = map.curves().to_vec();
    reversed.reverse();
    assert!(PumpSpeedMap::try_new(reversed, map.source().clone()).is_err());
    let mut missing = map.curves().to_vec();
    missing[0].curve = pump(false);
    assert!(PumpSpeedMap::try_new(missing, map.source().clone()).is_err());
    let mut incompatible = map.curves().to_vec();
    let input = incompatible[1].curve.drive_curve().unwrap();
    let rerouted = PumpDriveCurve::try_new(
        input.points().to_vec(),
        input.density_domain_kg_m3(),
        input.viscosity_domain_m2_s(),
        PumpLossRouting {
            pump_to_reservoir: 1.0,
            drive_to_reservoir: 1.0,
        },
        input.source().clone(),
    )
    .unwrap();
    incompatible[1].curve = incompatible[1].curve.clone().with_drive_curve(rerouted);
    assert!(PumpSpeedMap::try_new(incompatible, map.source().clone()).is_err());
    let mut disjoint = map.curves().to_vec();
    disjoint[1] = calibrated(
        3000.0,
        &[(2e-5, 0.01), (3e-5, 0.0)],
        &[(2e-5, 0.01, 0.02), (3e-5, 0.02, 0.03)],
        (990.0, 1010.0),
        (1e-6, 4e-6),
    );
    assert!(PumpSpeedMap::try_new(disjoint, map.source().clone()).is_err());

    // Admitting sampled inputs does not certify the nonlinear H(Q)*Q output
    // between samples. The selected speed still enforces the native gate.
    let underpowered = PumpSpeedMap::try_new(
        vec![
            calibrated(
                1000.0,
                &[(0.0, 0.01), (1e-5, 0.0)],
                &[(0.0, 0.0001, 0.001), (1e-5, 0.0001, 0.001)],
                (990.0, 1010.0),
                (1e-6, 4e-6),
            ),
            calibrated(
                3000.0,
                &[(0.0, 0.01), (1e-5, 0.0)],
                &[(0.0, 0.0002, 0.002), (1e-5, 0.0002, 0.002)],
                (990.0, 1010.0),
                (1e-6, 4e-6),
            ),
        ],
        map.source().clone(),
    )
    .unwrap()
    .at_speed(2000.0)
    .unwrap();
    with_network(FIXED, |cx, network| {
        assert!(
            network
                .solve_with_pump(cx, &[300.0; 2], &underpowered, local(), work(), root())
                .is_err()
        );
    });
}

#[test]
fn explicitly_calibrated_stopped_drive_retains_standby_heat_without_shaft_work() {
    let stopped = calibrated(
        0.0,
        &[(0.0, 0.0), (1e-5, 0.0)],
        &[(0.0, 0.0, 0.002), (1e-5, 0.0, 0.002)],
        (990.0, 1010.0),
        (1e-6, 4e-6),
    );
    let moving = family().curves()[0].clone();
    let map =
        PumpSpeedMap::try_new(vec![stopped, moving.clone()], family().source().clone()).unwrap();
    with_network(FIXED, |cx, network| {
        let state = network
            .solve_pressure_driven(cx, &[300.0; 2], Pressure::new(50.0), local(), work())
            .unwrap();
        let stopped = map
            .at_speed(0.0)
            .unwrap()
            .evaluate(&state.inlet_material, VolumetricFlowRate::new(0.0))
            .unwrap();
        assert_eq!(stopped.speed_rpm, Some(0.0));
        assert_eq!(stopped.hydraulic_power_w, 0.0);
        let input = stopped.drive_power.unwrap();
        assert_eq!(input.shaft_input_w, 0.0);
        close(input.electrical_input_w, 0.002, 1e-15);
        close(
            input.reservoir_heat_w + input.surroundings_heat_w,
            0.002,
            1e-15,
        );
    });
    let mut impossible = moving.clone();
    impossible.speed_rpm = 0.0;
    assert!(PumpSpeedMap::try_new(vec![impossible, moving], family().source().clone()).is_err());
}
