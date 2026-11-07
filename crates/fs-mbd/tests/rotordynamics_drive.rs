//! Electrical source limits and independently conserved native rotor drive work.
use fs_mbd::rotordynamics::transient::drive::{
    RotorCurrentDrive, RotorDriveError, RotorDrivePowerControl,
};
use fs_mbd::rotordynamics::transient::{RotorPowerError, RotorTransientState};
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::cell::Cell;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}

fn state(speed: f64) -> RotorTransientState {
    RotorTransientState {
        displacement: vec![0.0],
        velocity: vec![speed],
    }
}

fn drive() -> RotorCurrentDrive {
    RotorCurrentDrive {
        shaft_coordinate: 0,
        torque_constant_nm_per_a: 1.0,
        winding_resistance_ohm: 0.0,
        maximum_current_a: 10.0,
        maximum_terminal_voltage_v: 100.0,
        minimum_speed_rad_s: -100.0,
        maximum_speed_rad_s: 100.0,
        motoring_efficiency: 0.8,
        regenerative_efficiency: 0.5,
    }
}

fn limits(input: f64, returned: f64) -> RotorDrivePowerControl {
    RotorDrivePowerControl {
        maximum_source_input_w: input,
        maximum_source_return_w: returned,
        absolute_power_tolerance_w: 1e-10,
        absolute_voltage_tolerance_v: 1e-12,
    }
}

#[test]
fn holding_torque_draws_current_and_heats_the_winding_at_stall() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.5, 1e-12, || false)
        .unwrap();
    let motor = RotorCurrentDrive {
        torque_constant_nm_per_a: 3.0,
        winding_resistance_ohm: 2.0,
        ..drive()
    };
    let input = state(0.0);
    let output = stepper
        .step_current_drive(&input, &[-6.0], 2.0, motor, limits(11.0, 0.0), || false)
        .unwrap();
    assert_eq!(output.motion.state_after, input);
    assert_eq!(output.current_fraction, 1.0);
    assert_eq!(output.current_a, 2.0);
    assert_eq!(output.torque_nm, 6.0);
    assert_eq!(output.initial_terminal_voltage_v, 4.0);
    assert_eq!(output.final_terminal_voltage_v, 4.0);
    assert_eq!(output.initial_source_power_w, 10.0);
    assert_eq!(output.final_source_power_w, 10.0);
    assert_eq!(output.shaft_work_j, 0.0);
    assert_eq!(output.external_work_j, 0.0);
    near(output.electrical_draw_j, 5.0, 1e-14);
    assert_eq!(output.electrical_return_j, 0.0);
    near(output.winding_heat_j, 2.0 * 2.0_f64.powi(2) * 0.5, 1e-14);
    near(output.bridge_heat_j, 1.0, 1e-14);
    near(output.electrical_balance_residual_j, 0.0, 1e-14);
}

#[test]
fn source_voltage_and_current_capabilities_act_on_actual_stall_acceleration() {
    let model = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.25, 1e-12, || false)
        .unwrap();
    let motor = RotorCurrentDrive {
        winding_resistance_ohm: 0.5,
        ..drive()
    };
    for sign in [-1.0, 1.0] {
        // omega1 = alpha, V1 = 5 alpha and I = 8 alpha, so the
        // electrical source peak is 50 alpha² W including R I².
        let power = stepper
            .step_current_drive(
                &state(0.0),
                &[0.0],
                sign * 8.0,
                motor,
                limits(12.5, 0.0),
                || false,
            )
            .unwrap();
        near(power.current_fraction, 0.5, 2e-14);
        near(power.current_a, sign * 4.0, 2e-13);
        near(power.torque_nm, sign * 4.0, 2e-13);
        near(power.motion.state_after.velocity[0], sign * 0.5, 2e-14);
        near(power.initial_source_power_w, 10.0, 2e-12);
        near(power.final_source_power_w, 12.5, 2e-12);
        near(power.shaft_work_j, 0.25, 2e-13);
        near(power.winding_heat_j, 2.0, 2e-13);
        near(power.electrical_draw_j, 2.8125, 2e-13);
        near(power.bridge_heat_j, 0.5625, 2e-13);

        let voltage = stepper
            .step_current_drive(
                &state(0.0),
                &[0.0],
                sign * 8.0,
                RotorCurrentDrive {
                    maximum_terminal_voltage_v: 2.0,
                    ..motor
                },
                limits(12.5, 0.0),
                || false,
            )
            .unwrap();
        near(voltage.current_fraction, 0.4, 2e-14);
        near(voltage.final_terminal_voltage_v, sign * 2.0, 2e-13);
        near(voltage.motion.state_after.velocity[0], sign * 0.4, 2e-14);
        near(voltage.final_source_power_w, 8.0, 2e-12);

        let current = stepper
            .step_current_drive(
                &state(0.0),
                &[0.0],
                sign * 8.0,
                RotorCurrentDrive {
                    maximum_current_a: 2.0,
                    ..motor
                },
                limits(12.5, 0.0),
                || false,
            )
            .unwrap();
        near(current.current_fraction, 0.25, 2e-14);
        assert_eq!(current.current_a, sign * 2.0);
        near(current.motion.state_after.velocity[0], sign * 0.25, 2e-14);
        near(current.final_source_power_w, 3.125, 2e-13);
    }
}

#[test]
fn reversal_integrates_terminal_work_in_both_directions_before_bridge_losses() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let motor = RotorCurrentDrive {
        winding_resistance_ohm: 0.25,
        ..drive()
    };
    let output = stepper
        .step_current_drive(&state(-1.0), &[0.0], 2.0, motor, limits(4.0, 0.5), || false)
        .unwrap();
    assert_eq!(output.current_fraction, 1.0);
    near(output.motion.state_after.velocity[0], 1.0, 1e-14);
    near(output.shaft_work_j, 0.0, 1e-14);
    // omega=-1+2t, I=2, R=.25: V I=-1+4t. Terminal flow reverses
    // at t=.25, while the mechanical port reverses later at t=.5.
    let recovered_terminal_j = 0.5 * 0.25 * 1.0;
    let consumed_terminal_j = 0.5 * 0.75 * 3.0;
    near(output.initial_source_power_w, -0.5, 1e-14);
    near(output.final_source_power_w, 3.75, 1e-14);
    near(output.electrical_draw_j, consumed_terminal_j / 0.8, 1e-14);
    near(
        output.electrical_return_j,
        recovered_terminal_j * 0.5,
        1e-14,
    );
    near(output.winding_heat_j, 1.0, 1e-14);
    near(
        output.bridge_heat_j,
        consumed_terminal_j * 0.25 + recovered_terminal_j * 0.5,
        1e-14,
    );
    near(
        output.electrical_draw_j - output.electrical_return_j,
        1.34375,
        1e-14,
    );
    near(output.electrical_balance_residual_j, 0.0, 1e-14);

    // Entirely braking shaft work can still require electrical input when
    // copper loss exceeds recovered shaft power. Bridge mode follows V I.
    let short = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.1, 1e-12, || false)
        .unwrap();
    let braking = short
        .step_current_drive(
            &state(-0.1),
            &[0.0],
            1.0,
            RotorCurrentDrive {
                winding_resistance_ohm: 2.0,
                ..drive()
            },
            limits(10.0, 0.0),
            || false,
        )
        .unwrap();
    near(braking.shaft_work_j, -0.005, 1e-14);
    near(braking.winding_heat_j, 0.2, 1e-14);
    near(braking.electrical_draw_j, 0.195 / 0.8, 1e-14);
    assert_eq!(braking.electrical_return_j, 0.0);
    assert!(braking.initial_source_power_w > 0.0);
}

#[test]
fn regeneration_respects_source_acceptance_and_zero_efficiency_dissipates_braking() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let limited = stepper
        .step_current_drive(
            &state(-1.0),
            &[0.0],
            4.0,
            drive(),
            limits(100.0, 0.25),
            || false,
        )
        .unwrap();
    near(limited.current_fraction, 0.125, 1e-14);
    near(limited.initial_source_power_w, -0.25, 1e-14);
    near(limited.final_source_power_w, -0.125, 1e-14);
    near(limited.shaft_work_j, -0.375, 1e-14);
    near(limited.electrical_return_j, 0.1875, 1e-14);
    near(limited.bridge_heat_j, 0.1875, 1e-14);

    let dissipative = stepper
        .step_current_drive(
            &state(-2.0),
            &[0.0],
            1.0,
            RotorCurrentDrive {
                regenerative_efficiency: 0.0,
                ..drive()
            },
            limits(0.0, 0.0),
            || false,
        )
        .unwrap();
    assert_eq!(dissipative.current_fraction, 1.0);
    assert_eq!(dissipative.electrical_draw_j, 0.0);
    assert_eq!(dissipative.electrical_return_j, 0.0);
    near(dissipative.shaft_work_j, -1.5, 1e-14);
    near(dissipative.bridge_heat_j, 1.5, 1e-14);
    near(dissipative.electrical_balance_residual_j, 0.0, 1e-14);
}

#[test]
fn both_terminal_power_constraints_survive_disconnected_active_bearing_intervals() {
    let model = RotorBearingModel::new(1, &[1.0], &[-2.5], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let motor = RotorCurrentDrive {
        winding_resistance_ohm: 4.0,
        motoring_efficiency: 1.0,
        regenerative_efficiency: 1.0,
        ..drive()
    };
    // The ORIGINAL active model gives omega1=-.2-4 alpha. Endpoint
    // terminal powers are 4 alpha²-4 alpha and -.2 alpha. Capping only
    // the first then the second would pick alpha=.5 and violate the first.
    let output = stepper
        .step_current_drive(&state(-4.0), &[9.05], 1.0, motor, limits(0.1, 0.1), || {
            false
        })
        .unwrap();
    let expected = 0.5 * (1.0 - 0.9_f64.sqrt());
    near(output.current_fraction, expected, 1e-13);
    near(output.initial_source_power_w, -0.1, 1e-13);
    near(output.final_source_power_w, -0.2 * expected, 1e-13);
    near(
        output.motion.state_after.velocity[0],
        -0.2 - 4.0 * expected,
        1e-13,
    );
    near(output.electrical_balance_residual_j, 0.0, 1e-13);
    assert_eq!(output.bridge_heat_j, 0.0);
}

#[test]
fn voltage_compliance_can_require_nonzero_current_and_refuses_incompatible_source_limits() {
    let model = RotorBearingModel::new(1, &[100.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let motor = RotorCurrentDrive {
        winding_resistance_ohm: 1.0,
        maximum_terminal_voltage_v: 2.0,
        ..drive()
    };
    let input = state(-3.0);
    // V0=-3+4 alpha, so alpha>=.25 is required; zero torque cannot
    // maintain zero winding current against this back-EMF with only 2 V.
    let output = stepper
        .step_current_drive(&input, &[0.0], 4.0, motor, limits(10.0, 10.0), || false)
        .unwrap();
    assert_eq!(output.current_fraction, 1.0);
    near(output.initial_terminal_voltage_v, 1.0, 1e-14);
    near(output.final_terminal_voltage_v, 1.04, 1e-14);
    assert!(matches!(
        stepper.step_current_drive(&input, &[0.0], 4.0, motor, limits(0.0, 0.0), || false),
        Err(RotorDriveError::VoltageCompliance { .. })
    ));
    assert_eq!(input, state(-3.0));
    assert!(matches!(
        stepper.step_current_drive(
            &state(5.0),
            &[0.0],
            1.0,
            RotorCurrentDrive {
                winding_resistance_ohm: 0.0,
                ..motor
            },
            limits(100.0, 100.0),
            || false
        ),
        Err(RotorDriveError::NoVoltageFeasibility)
    ));
}

#[test]
fn electrical_accounting_retains_full_native_tensor_gyro_and_external_load_work() {
    let mass = [2.0, 0.5, 0.5, 1.0];
    let damping = [0.2, 0.1, -0.4, 0.5];
    let stiffness = [3.0, 0.7, -0.1, 2.0];
    let gyro = [0.0, 0.3, -0.3, 0.0];
    let model = RotorBearingModel::new(2, &mass, &damping, &stiffness, &gyro, || false).unwrap();
    let stepper = model
        .at_speed(4.0, || false)
        .unwrap()
        .midpoint_stepper(0.05, 1e-12, || false)
        .unwrap();
    let input = RotorTransientState {
        displacement: vec![0.4, -0.2],
        velocity: vec![0.3, -0.7],
    };
    let external = [1.0, -0.2];
    let output = stepper
        .step_current_drive(
            &input,
            &external,
            3.0,
            RotorCurrentDrive {
                torque_constant_nm_per_a: 2.0,
                winding_resistance_ohm: 0.3,
                ..drive()
            },
            limits(100.0, 100.0),
            || false,
        )
        .unwrap();
    let after = &output.motion.state_after;
    let qmid = [
        0.5 * (input.displacement[0] + after.displacement[0]),
        0.5 * (input.displacement[1] + after.displacement[1]),
    ];
    let vmid = [
        0.5 * (input.velocity[0] + after.velocity[0]),
        0.5 * (input.velocity[1] + after.velocity[1]),
    ];
    let mut kinetic_change = 0.0;
    let mut bearing_work = 0.0;
    let mut stiffness_work = 0.0;
    let mut gyro_work = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            kinetic_change += 0.5
                * mass[2 * i + j]
                * (after.velocity[i] * after.velocity[j] - input.velocity[i] * input.velocity[j]);
            bearing_work += 0.05 * vmid[i] * damping[2 * i + j] * vmid[j];
            stiffness_work += 0.05 * vmid[i] * stiffness[2 * i + j] * qmid[j];
            gyro_work += 0.05 * 4.0 * vmid[i] * gyro[2 * i + j] * vmid[j];
        }
    }
    let external_work = 0.05 * (external[0] * vmid[0] + external[1] * vmid[1]);
    let shaft_work = 0.05 * 6.0 * vmid[0];
    near(output.external_work_j, external_work, 2e-14);
    near(output.shaft_work_j, shaft_work, 2e-14);
    near(
        kinetic_change + bearing_work + stiffness_work + gyro_work - external_work,
        shaft_work,
        2e-13,
    );
    near(
        output.electrical_draw_j - output.electrical_return_j,
        kinetic_change + bearing_work + stiffness_work - external_work
            + output.winding_heat_j
            + output.bridge_heat_j,
        2e-13,
    );
    near(gyro_work, 0.0, 2e-14);
    near(output.load_partition_residual_j, 0.0, 2e-14);
}

#[test]
fn invalid_models_and_unsupported_speeds_refuse_without_silent_extrapolation() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let input = state(0.0);
    for motor in [
        RotorCurrentDrive {
            shaft_coordinate: 1,
            ..drive()
        },
        RotorCurrentDrive {
            torque_constant_nm_per_a: 0.0,
            ..drive()
        },
        RotorCurrentDrive {
            winding_resistance_ohm: -1.0,
            ..drive()
        },
        RotorCurrentDrive {
            motoring_efficiency: 0.0,
            ..drive()
        },
        RotorCurrentDrive {
            motoring_efficiency: 1.1,
            ..drive()
        },
        RotorCurrentDrive {
            regenerative_efficiency: -0.1,
            ..drive()
        },
        RotorCurrentDrive {
            maximum_terminal_voltage_v: f64::NAN,
            ..drive()
        },
        RotorCurrentDrive {
            minimum_speed_rad_s: 101.0,
            ..drive()
        },
    ] {
        assert!(matches!(
            stepper.step_current_drive(&input, &[0.0], 1.0, motor, limits(100.0, 100.0), || false),
            Err(RotorDriveError::Invalid(_))
        ));
    }
    assert!(matches!(
        stepper.step_current_drive(
            &input,
            &[0.0],
            1.0,
            RotorCurrentDrive {
                maximum_speed_rad_s: 0.1,
                ..drive()
            },
            limits(100.0, 100.0),
            || false
        ),
        Err(RotorDriveError::SpeedDomain { .. })
    ));
    assert!(matches!(
        stepper.step_current_drive(
            &state(-0.1),
            &[0.0],
            1.0,
            RotorCurrentDrive {
                minimum_speed_rad_s: 0.0,
                ..drive()
            },
            limits(100.0, 100.0),
            || false
        ),
        Err(RotorDriveError::SpeedDomain { .. })
    ));
    assert_eq!(input, state(0.0));
}

#[test]
fn every_observed_cancellation_boundary_preserves_borrowed_inputs() {
    let model = RotorBearingModel::new(1, &[1.3], &[0.07], &[2.1], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.3, 1e-12, || false)
        .unwrap();
    let input = state(0.37);
    let saved = input.clone();
    let external = [-0.2];
    let motor = RotorCurrentDrive {
        winding_resistance_ohm: 0.2,
        ..drive()
    };
    let control = limits(0.43, 0.2);
    let calls = Cell::new(0);
    let accepted = stepper
        .step_current_drive(&input, &external, 9.0, motor, control, || {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    for stop in 1..=calls.get() {
        let count = Cell::new(0);
        assert!(matches!(
            stepper.step_current_drive(&input, &external, 9.0, motor, control, || {
                count.set(count.get() + 1);
                count.get() == stop
            }),
            Err(RotorDriveError::Power(RotorPowerError::Dynamics(
                RotorError::Cancelled
            )))
        ));
        assert_eq!(input, saved);
        assert_eq!(external, [-0.2]);
    }
    assert_eq!(
        stepper
            .step_current_drive(&input, &external, 9.0, motor, control, || false)
            .unwrap(),
        accepted
    );
}
