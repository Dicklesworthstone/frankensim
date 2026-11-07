//! G1/G3/G4: retained magnetic storage, material feedback and atomic publication.
use super::*;
use fs_mbd::rotordynamics::transient::drive::voltage::{
    RotorVoltageControl, RotorWindingMotor, RotorWindingState,
};

const INDUCTANCE: f64 = 0.2;
fn voltage_power(input: f64) -> RotorVoltageControl {
    RotorVoltageControl {
        source: power(input),
        absolute_current_tolerance_a: 1e-10,
    }
}

#[test]
fn inductive_startup_closes_temperature_current_and_both_stored_energies() {
    let cp = [1.0, 1.0];
    let (cards, curves) = sources(cp, ALPHA);
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: winding(),
        bridge_cell: BRIDGE,
    };
    let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = rotor
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(DT, 1e-12, || false)
        .unwrap();
    let native = stepper
        .voltage_motor(
            RotorWindingMotor {
                drive: drive(),
                inductance_h: INDUCTANCE,
            },
            || false,
        )
        .unwrap();
    let thermal = body(&curves, cp, false, 0.0);
    let initial = RotorWindingState {
        rotor: state(0.0),
        current_a: 0.0,
    };
    let motor = network
        .bind_voltage(&thermal, &native, initial.clone(), || false)
        .unwrap();
    near(motor.parameters().drive.winding_resistance_ohm, R0);
    near(native.motor().drive.winding_resistance_ohm, 999.0);
    let result = network
        .step_voltage_candidate(
            &thermal,
            DT,
            &motor,
            &[0.0],
            12.0,
            voltage_power(1000.0),
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let report = result.report();
    let step = &report.electrical.drive;
    // Independent scalar midpoint/caloric solution for this authored fixture:
    // I_mid = h V / (2L+hR(T1)+h²Kt²/(2J)), Cp*(T1-T0)=hR(T1)*I_mid².
    let winding_rise = 0.35201355809956203;
    near(
        result.thermal().temperature_k(WINDING).unwrap(),
        300.0 + winding_rise,
    );
    near(
        result.thermal().temperature_k(BRIDGE).unwrap(),
        300.7954594421233891,
    );
    near(result.motor().state().current_a, 5.303062947489261);
    near(
        result.motor().state().rotor.velocity[0],
        0.13257657368723152,
    );
    near(
        result.motor().parameters().drive.winding_resistance_ohm,
        0.5006864264382942,
    );
    near(step.magnetic_energy_change_j, 2.812247662503349);
    near(step.motion.kinetic_energy_change_j, 0.01757654789064593);
    near(step.winding_heat_j, winding_rise);
    near(step.bridge_heat_j, 0.7954594421233891);
    near(step.electrical_draw_j, 3.977297210616946);
    near(step.electrical_return_j, 0.0);
    let r = result.motor().parameters().drive.winding_resistance_ohm;
    let current_mid = DT * 12.0 / (2.0 * INDUCTANCE + DT * r + DT * DT / 4.0);
    near(step.midpoint_current_a, current_mid);
    near(step.winding_heat_j, DT * r * current_mid * current_mid);
    near(
        step.magnetic_energy_after_j,
        0.5 * INDUCTANCE * result.motor().state().current_a.powi(2),
    );
    near(
        report.electrical.supplied_to_stored_and_heat_j,
        step.motion.kinetic_energy_change_j
            + step.magnetic_energy_change_j
            + report.deposited_joule_j,
    );
    near(report.total_energy_residual_j, 0.0);
    assert_eq!(motor.state(), &initial);
    assert_eq!(thermal.time_s(), 0.0);

    // Finite source limiting is re-evaluated with the accepted hot resistance.
    let limited = network
        .step_voltage_candidate(
            &thermal,
            DT,
            &motor,
            &[0.0],
            12.0,
            voltage_power(10.0),
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let actual = &limited.report().electrical.drive;
    let r = limited.motor().parameters().drive.winding_resistance_ohm;
    assert!(actual.applied_voltage_v < 12.0);
    near(actual.final_source_power_w, 10.0);
    near(
        actual.midpoint_current_a,
        DT * actual.applied_voltage_v / (2.0 * INDUCTANCE + DT * r + DT * DT / 4.0),
    );
    near(
        limited.thermal().temperature_k(WINDING).unwrap() - 300.0,
        DT * r * actual.midpoint_current_a.powi(2),
    );
    near(limited.report().total_energy_residual_j, 0.0);
}

#[test]
fn continued_inductive_regeneration_keeps_all_histories_on_failed_trials() {
    let cp = [1.0, 1.0];
    let (cards, curves) = sources(cp, ALPHA);
    let context = common::point(300.0);
    let cells = bindings(&cards, &curves, &context, None);
    let network = MotorThermalNetwork {
        cells: &cells,
        links: &[],
        winding: winding(),
        bridge_cell: BRIDGE,
    };
    let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = rotor
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(DT, 1e-12, || false)
        .unwrap();
    let native = stepper
        .voltage_motor(
            RotorWindingMotor {
                drive: drive(),
                inductance_h: INDUCTANCE,
            },
            || false,
        )
        .unwrap();
    let mut thermal = body(&curves, cp, false, 0.0);
    let mut motor = network
        .bind_voltage(
            &thermal,
            &native,
            RotorWindingState {
                rotor: state(10.0),
                current_a: -1.0,
            },
            || false,
        )
        .unwrap();
    for index in 0..2 {
        let prior_cells = thermal.cells().clone();
        let prior_state = motor.state().clone();
        let prior_r = motor.parameters().drive.winding_resistance_ohm;
        let target = (index + 1) as f64 * DT;
        let calls = Cell::new(0usize);
        let result = network
            .step_voltage_candidate(
                &thermal,
                target,
                &motor,
                &[0.0],
                5.0,
                voltage_power(0.0),
                controls(),
                &[0.0; 2],
                || {
                    calls.set(calls.get() + 1);
                    false
                },
            )
            .unwrap();
        let step = &result.report().electrical.drive;
        let r = result.motor().parameters().drive.winding_resistance_ohm;
        let i0 = prior_state.current_a;
        let i1 = result.motor().state().current_a;
        let v0 = prior_state.rotor.velocity[0];
        let v1 = result.motor().state().rotor.velocity[0];
        let mean_i = 0.5 * (i0 + i1);
        let mean_speed = 0.5 * (v0 + v1);
        let magnetic_change = 0.5 * INDUCTANCE * (i1 * i1 - i0 * i0);
        let kinetic_change = v1 * v1 - v0 * v0; // J=2 kg m².
        assert!(i1 < i0 && i0 < 0.0);
        assert!(step.electrical_return_j > 0.0);
        near(step.electrical_draw_j, 0.0);
        near(step.applied_voltage_v, 5.0);
        near(INDUCTANCE * (i1 - i0) / DT + r * mean_i + mean_speed, 5.0);
        near(v1 - v0, DT * mean_i / 2.0);
        near(step.magnetic_energy_before_j, 0.5 * INDUCTANCE * i0 * i0);
        near(step.magnetic_energy_change_j, magnetic_change);
        near(step.motion.kinetic_energy_change_j, kinetic_change);
        near(
            step.electrical_return_j,
            -DT * 5.0 * mean_i * drive().regenerative_efficiency,
        );
        for (cell, heat) in [(WINDING, step.winding_heat_j), (BRIDGE, step.bridge_heat_j)] {
            near(
                result.thermal().temperature_k(cell).unwrap()
                    - thermal.temperature_k(cell).unwrap(),
                heat,
            );
        }
        near(
            kinetic_change
                + magnetic_change
                + result.report().deposited_joule_j
                + step.electrical_return_j,
            0.0,
        );
        near(result.report().total_energy_residual_j, 0.0);
        let last = calls.get();
        calls.set(0);
        assert!(matches!(
            network.step_voltage_candidate(
                &thermal,
                target,
                &motor,
                &[0.0],
                5.0,
                voltage_power(0.0),
                controls(),
                &[0.0; 2],
                || {
                    calls.set(calls.get() + 1);
                    calls.get() >= last
                }
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert!(
            network
                .step_voltage_candidate(
                    &thermal,
                    target,
                    &motor,
                    &[0.0],
                    5.0,
                    voltage_power(0.0),
                    controls(),
                    &[1e6; 2],
                    || false
                )
                .is_err()
        );
        assert!(
            network
                .step_voltage_candidate(
                    &thermal,
                    target + DT,
                    &motor,
                    &[0.0],
                    5.0,
                    voltage_power(0.0),
                    controls(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        assert!(matches!(
            network.step_voltage_candidate(
                &thermal,
                target,
                &motor,
                &[0.0],
                f64::NAN,
                voltage_power(0.0),
                controls(),
                &[0.0; 2],
                || false
            ),
            Err(CircuitEnthalpyError::Source(
                MotorThermalError::VoltageDrive(_)
            ))
        ));
        assert_eq!(thermal.cells(), &prior_cells);
        assert_eq!(motor.state(), &prior_state);
        near(motor.parameters().drive.winding_resistance_ohm, prior_r);
        (thermal, motor, _) = result.into_parts();
    }
}
