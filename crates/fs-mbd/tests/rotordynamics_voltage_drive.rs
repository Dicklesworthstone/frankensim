//! Finite winding storage, original rotor coupling and independent power balance.
use fs_mbd::rotordynamics::transient::drive::voltage::{
    RotorVoltageControl, RotorVoltageError, RotorWindingMotor, RotorWindingState,
};
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

fn state(speed: f64, current: f64) -> RotorWindingState {
    RotorWindingState {
        rotor: RotorTransientState {
            displacement: vec![0.0],
            velocity: vec![speed],
        },
        current_a: current,
    }
}

fn motor() -> RotorWindingMotor {
    RotorWindingMotor {
        drive: RotorCurrentDrive {
            shaft_coordinate: 0,
            torque_constant_nm_per_a: 1.0,
            winding_resistance_ohm: 0.0,
            maximum_current_a: 10.0,
            maximum_terminal_voltage_v: 100.0,
            minimum_speed_rad_s: -100.0,
            maximum_speed_rad_s: 100.0,
            motoring_efficiency: 0.8,
            regenerative_efficiency: 0.5,
        },
        inductance_h: 1.0,
    }
}

fn limits(input: f64, returned: f64) -> RotorVoltageControl {
    RotorVoltageControl {
        source: RotorDrivePowerControl {
            maximum_source_input_w: input,
            maximum_source_return_w: returned,
            absolute_power_tolerance_w: 1e-10,
            absolute_voltage_tolerance_v: 1e-10,
        },
        absolute_current_tolerance_a: 1e-12,
    }
}

#[test]
fn locked_shaft_rl_charge_and_decay_preserve_current_storage_and_midpoint_heat() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.1, 1e-12, || false)
        .unwrap();
    let parameters = RotorWindingMotor {
        drive: RotorCurrentDrive {
            torque_constant_nm_per_a: 2.0,
            winding_resistance_ohm: 2.0,
            motoring_efficiency: 1.0,
            ..motor().drive
        },
        inductance_h: 0.5,
    };
    let prepared = stepper.voltage_motor(parameters, || false).unwrap();
    let mut input = state(0.0, 0.0);
    for voltage in [10.0, 0.0] {
        // Independent closed-form RL midpoint update, with an external shaft
        // reaction equal to the expected motor torque. This gives zero motion
        // without disabling reciprocal back-EMF or bypassing native mechanics.
        let expected_current = ((2.0 * 0.5 - 0.1 * 2.0) * input.current_a + 2.0 * 0.1 * voltage)
            / (2.0 * 0.5 + 0.1 * 2.0);
        let expected_midpoint = 0.5 * (input.current_a + expected_current);
        let output = prepared
            .step(
                &input,
                &[-2.0 * expected_midpoint],
                voltage,
                limits(1000.0, 1000.0),
                || false,
            )
            .unwrap();
        near(output.winding_current_after_a, expected_current, 1e-13);
        near(output.midpoint_current_a, expected_midpoint, 1e-13);
        near(output.motion.state_after.velocity[0], 0.0, 1e-13);
        near(
            output.magnetic_energy_before_j,
            0.25 * input.current_a.powi(2),
            1e-13,
        );
        near(
            output.magnetic_energy_after_j,
            0.25 * expected_current.powi(2),
            1e-13,
        );
        near(
            output.winding_heat_j,
            0.1 * 2.0 * expected_midpoint.powi(2),
            1e-13,
        );
        near(output.shaft_work_j, 0.0, 1e-13);
        near(
            output.electrical_draw_j,
            0.1 * voltage * expected_midpoint,
            1e-13,
        );
        near(output.electrical_balance_residual_j, 0.0, 1e-13);
        if voltage == 0.0 {
            assert!(output.magnetic_energy_change_j < 0.0);
            assert!(output.winding_heat_j > 0.0);
            assert_eq!(output.electrical_draw_j, 0.0);
        } else {
            // The declared midpoint Joule quadrature must not be relabeled as
            // the exact integral of R I(t)^2 along a linear current ramp.
            let linear_ramp_heat = (0.1 * 2.0 / 3.0)
                * (input.current_a.powi(2)
                    + input.current_a * expected_current
                    + expected_current.powi(2));
            assert!(linear_ramp_heat > output.winding_heat_j);
        }
        input = output.into_state();
    }
    near(input.current_a, 10.0 / 9.0, 1e-13);
}

#[test]
fn rl_current_converges_at_second_order_to_the_continuous_solution() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let parameters = RotorWindingMotor {
        drive: RotorCurrentDrive {
            winding_resistance_ohm: 2.0,
            motoring_efficiency: 1.0,
            ..motor().drive
        },
        inductance_h: 0.5,
    };
    let exact = 5.0 * (1.0 - (-2.0_f64).exp());
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let h = 0.5 / f64::from(steps);
        let stepper = model
            .at_speed(0.0, || false)
            .unwrap()
            .midpoint_stepper(h, 1e-12, || false)
            .unwrap();
        let prepared = stepper.voltage_motor(parameters, || false).unwrap();
        let mut input = state(0.0, 0.0);
        for _ in 0..steps {
            let next_current = ((1.0 - 2.0 * h) * input.current_a + 20.0 * h) / (1.0 + 2.0 * h);
            let reaction = -0.5 * (input.current_a + next_current);
            let output = prepared
                .step(&input, &[reaction], 10.0, limits(1000.0, 1000.0), || false)
                .unwrap();
            near(output.motion.state_after.velocity[0], 0.0, 1e-12);
            input = output.into_state();
        }
        errors.push((input.current_a - exact).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            (3.8..4.3).contains(&(pair[0] / pair[1])),
            "RL errors {errors:?}"
        );
    }
}

#[test]
fn zero_terminal_voltage_exchanges_magnetic_and_mechanical_storage_without_energy_creation() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let prepared = stepper.voltage_motor(motor(), || false).unwrap();
    let first = prepared
        .step(&state(0.0, 1.0), &[0.0], 0.0, limits(0.0, 0.0), || false)
        .unwrap();
    // Direct 2x2 midpoint solution of omega'=I and I'=-omega.
    near(first.winding_current_after_a, 0.6, 1e-14);
    near(first.motion.state_after.velocity[0], 0.8, 1e-14);
    near(first.magnetic_energy_change_j, -0.32, 1e-14);
    near(first.shaft_work_j, 0.32, 1e-14);
    assert_eq!(first.electrical_draw_j, 0.0);
    assert_eq!(first.electrical_return_j, 0.0);
    assert_eq!(first.winding_heat_j, 0.0);
    assert_eq!(first.bridge_heat_j, 0.0);
    let mut input = first.into_state();
    for _ in 0..24 {
        let output = prepared
            .step(&input, &[0.0], 0.0, limits(0.0, 0.0), || false)
            .unwrap();
        near(
            0.5 * output.motion.state_after.velocity[0].powi(2) + output.magnetic_energy_after_j,
            0.5,
            1e-12,
        );
        near(output.electrical_balance_residual_j, 0.0, 1e-13);
        input = output.into_state();
    }
    let generating = prepared
        .step(&state(1.0, 0.0), &[0.0], 0.0, limits(0.0, 0.0), || false)
        .unwrap();
    near(generating.winding_current_after_a, -0.8, 1e-14);
    near(generating.motion.state_after.velocity[0], 0.6, 1e-14);
    near(generating.magnetic_energy_change_j, 0.32, 1e-14);
    near(generating.shaft_work_j, -0.32, 1e-14);
}

#[test]
fn finite_voltage_current_and_source_limits_act_on_the_coupled_current_ramp() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    for sign in [-1.0, 1.0] {
        let prepared = stepper.voltage_motor(motor(), || false).unwrap();
        // I1=.8 V and omega1=.4 V. With command 4 V, the linear-current
        // endpoint source peak is 16 alpha² W after bridge conversion.
        let source = prepared
            .step(
                &state(0.0, 0.0),
                &[0.0],
                sign * 4.0,
                limits(4.0, 0.0),
                || false,
            )
            .unwrap();
        near(source.voltage_fraction, 0.5, 1e-14);
        near(source.applied_voltage_v, sign * 2.0, 1e-14);
        near(source.winding_current_after_a, sign * 1.6, 1e-13);
        near(source.motion.state_after.velocity[0], sign * 0.8, 1e-13);
        near(source.final_source_power_w, 4.0, 1e-13);
        near(source.electrical_draw_j, 2.0, 1e-13);
        near(source.magnetic_energy_change_j, 1.28, 1e-13);
        near(source.shaft_work_j, 0.32, 1e-13);
        near(source.bridge_heat_j, 0.4, 1e-13);

        let current_model = RotorWindingMotor {
            drive: RotorCurrentDrive {
                maximum_current_a: 1.0,
                ..motor().drive
            },
            ..motor()
        };
        let current_drive = stepper.voltage_motor(current_model, || false).unwrap();
        let current = current_drive
            .step(
                &state(0.0, 0.0),
                &[0.0],
                sign * 4.0,
                limits(4.0, 0.0),
                || false,
            )
            .unwrap();
        near(current.voltage_fraction, 0.3125, 1e-14);
        near(current.winding_current_after_a, sign * 1.0, 1e-14);
        near(current.motion.state_after.velocity[0], sign * 0.5, 1e-14);
        near(current.final_source_power_w, 1.5625, 1e-13);

        let voltage_model = RotorWindingMotor {
            drive: RotorCurrentDrive {
                maximum_terminal_voltage_v: 1.0,
                ..motor().drive
            },
            ..motor()
        };
        let voltage_drive = stepper.voltage_motor(voltage_model, || false).unwrap();
        let voltage = voltage_drive
            .step(
                &state(0.0, 0.0),
                &[0.0],
                sign * 4.0,
                limits(4.0, 0.0),
                || false,
            )
            .unwrap();
        near(voltage.voltage_fraction, 0.25, 1e-14);
        near(voltage.applied_voltage_v, sign * 1.0, 1e-14);
        near(voltage.winding_current_after_a, sign * 0.8, 1e-14);
        near(voltage.motion.state_after.velocity[0], sign * 0.4, 1e-14);
    }
}

#[test]
fn current_reversal_and_regeneration_keep_magnetic_and_bridge_work_separate() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let prepared = stepper.voltage_motor(motor(), || false).unwrap();
    let reversal = prepared
        .step(&state(0.0, -1.0), &[0.0], 2.0, limits(100.0, 100.0), || {
            false
        })
        .unwrap();
    near(reversal.winding_current_after_a, 1.0, 1e-14);
    near(reversal.midpoint_current_a, 0.0, 1e-14);
    near(reversal.shaft_work_j, 0.0, 1e-14);
    near(reversal.magnetic_energy_change_j, 0.0, 1e-14);
    // V I(t)=-2+4t gives one .5 J triangle in EACH direction.
    near(reversal.electrical_draw_j, 0.5 / 0.8, 1e-14);
    near(reversal.electrical_return_j, 0.5 * 0.5, 1e-14);
    near(reversal.bridge_heat_j, 0.375, 1e-14);
    near(reversal.electrical_balance_residual_j, 0.0, 1e-14);

    let limited = prepared
        .step(&state(0.0, -1.0), &[0.0], 2.0, limits(100.0, 0.25), || {
            false
        })
        .unwrap();
    near(limited.voltage_fraction, 0.25, 1e-14);
    near(limited.initial_source_power_w, -0.25, 1e-14);
    near(limited.final_source_power_w, -0.05, 1e-14);
    near(limited.winding_current_after_a, -0.2, 1e-14);
    near(limited.magnetic_energy_change_j, -0.48, 1e-14);
    // Stored winding energy simultaneously accelerates the shaft and returns
    // source energy. Choosing bridge mode from SHAFT work would be incorrect.
    near(limited.shaft_work_j, 0.18, 1e-14);
    near(limited.electrical_return_j, 0.15, 1e-14);
    near(limited.bridge_heat_j, 0.15, 1e-14);
}

fn determinant(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

#[test]
fn scalar_schur_solution_matches_independent_full_tensor_electromechanical_equations() {
    let mass = [2.0, 0.5, 0.5, 1.0];
    let damping = [0.2, 0.1, -0.4, 0.5];
    let stiffness = [3.0, 0.7, -0.1, 2.0];
    let gyro = [0.0, 0.3, -0.3, 0.0];
    let model = RotorBearingModel::new(2, &mass, &damping, &stiffness, &gyro, || false).unwrap();
    let h = 0.05;
    let stepper = model
        .at_speed(4.0, || false)
        .unwrap()
        .midpoint_stepper(h, 1e-12, || false)
        .unwrap();
    let parameters = RotorWindingMotor {
        drive: RotorCurrentDrive {
            shaft_coordinate: 1,
            torque_constant_nm_per_a: 2.0,
            winding_resistance_ohm: 0.3,
            ..motor().drive
        },
        inductance_h: 0.4,
    };
    let prepared = stepper.voltage_motor(parameters, || false).unwrap();
    let input = RotorWindingState {
        rotor: RotorTransientState {
            displacement: vec![0.4, -0.2],
            velocity: vec![0.3, -0.7],
        },
        current_a: 0.2,
    };
    let external = [1.0, -0.2];
    let voltage = 3.0;
    // Independently assemble the ORIGINAL 3x3 coupled midpoint system in
    // (delta v0, delta v1, delta I), with no scalar-elimination implementation.
    let mut a = [[0.0; 3]; 3];
    let mut b = [0.0; 3];
    for i in 0..2 {
        b[i] = h * external[i];
        for j in 0..2 {
            let velocity_operator = damping[2 * i + j] + 4.0 * gyro[2 * i + j];
            a[i][j] =
                mass[2 * i + j] + 0.5 * h * velocity_operator + 0.25 * h * h * stiffness[2 * i + j];
            b[i] -= h
                * (velocity_operator * input.rotor.velocity[j]
                    + stiffness[2 * i + j]
                        * (input.rotor.displacement[j] + 0.5 * h * input.rotor.velocity[j]));
        }
    }
    a[1][2] = -0.5 * h * 2.0;
    b[1] += h * 2.0 * input.current_a;
    a[2][1] = 0.5 * h * 2.0;
    a[2][2] = 0.4 + 0.5 * h * 0.3;
    b[2] = h * (voltage - 0.3 * input.current_a - 2.0 * input.rotor.velocity[1]);
    let denominator = determinant(a);
    let mut expected_increment = [0.0; 3];
    for column in 0..3 {
        let mut replaced = a;
        for row in 0..3 {
            replaced[row][column] = b[row];
        }
        expected_increment[column] = determinant(replaced) / denominator;
    }
    let output = prepared
        .step(&input, &external, voltage, limits(1000.0, 1000.0), || false)
        .unwrap();
    assert_eq!(output.voltage_fraction, 1.0);
    for i in 0..2 {
        near(
            output.motion.state_after.velocity[i],
            input.rotor.velocity[i] + expected_increment[i],
            2e-13,
        );
        near(
            output.motion.state_after.displacement[i],
            input.rotor.displacement[i]
                + h * (input.rotor.velocity[i] + 0.5 * expected_increment[i]),
            2e-13,
        );
    }
    near(
        output.winding_current_after_a,
        input.current_a + expected_increment[2],
        2e-13,
    );
    let imid = input.current_a + 0.5 * expected_increment[2];
    let shaft_work = h * 2.0 * imid * (input.rotor.velocity[1] + 0.5 * expected_increment[1]);
    let external_work = h
        * (external[0] * (input.rotor.velocity[0] + 0.5 * expected_increment[0])
            + external[1] * (input.rotor.velocity[1] + 0.5 * expected_increment[1]));
    near(output.shaft_work_j, shaft_work, 2e-13);
    near(output.external_work_j, external_work, 2e-13);
    near(output.winding_heat_j, h * 0.3 * imid * imid, 2e-13);
    near(
        output.magnetic_energy_change_j,
        0.4 * imid * expected_increment[2],
        2e-13,
    );
    near(
        output.electrical_draw_j - output.electrical_return_j,
        output.magnetic_energy_change_j + shaft_work + output.winding_heat_j + output.bridge_heat_j,
        2e-13,
    );
    near(output.motion.gyroscopic_work_j, 0.0, 2e-13);
    near(output.load_partition_residual_j, 0.0, 2e-13);
}

#[test]
fn infeasible_source_cutback_cannot_clip_current_or_discard_winding_energy() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let parameters = RotorWindingMotor {
        drive: RotorCurrentDrive {
            maximum_current_a: 1.0,
            ..motor().drive
        },
        ..motor()
    };
    let prepared = stepper.voltage_motor(parameters, || false).unwrap();
    let input = state(-1.0, 1.0);
    let saved = input.clone();
    let accepted = prepared
        .step(&input, &[0.0], -2.0, limits(100.0, 100.0), || false)
        .unwrap();
    near(accepted.winding_current_after_a, -0.2, 1e-14);
    // At V=0 the actual back-EMF would raise I to 1.4 A. At least .5 V
    // of the negative command is required to remain under the current limit.
    // A source that refuses every return forces zero voltage, so REFUSE the
    // incompatible step rather than clipping I to 1 A or deleting magnetic work.
    assert!(matches!(
        prepared.step(&input, &[0.0], -2.0, limits(100.0, 0.0), || false),
        Err(RotorVoltageError::CurrentCompliance { .. })
    ));
    assert!(matches!(
        prepared.step(&input, &[0.0], 0.0, limits(100.0, 100.0), || false),
        Err(RotorVoltageError::NoCurrentFeasibility)
    ));
    assert_eq!(input, saved);
    assert_eq!(
        prepared
            .step(&input, &[0.0], -2.0, limits(100.0, 100.0), || false)
            .unwrap(),
        accepted
    );
}

#[test]
fn invalid_storage_singular_coupling_and_unrepresentable_current_changes_refuse() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    for inductance in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            stepper.voltage_motor(
                RotorWindingMotor {
                    inductance_h: inductance,
                    ..motor()
                },
                || false
            ),
            Err(RotorVoltageError::Invalid(_))
        ));
    }
    let prepared = stepper.voltage_motor(motor(), || false).unwrap();
    assert!(matches!(
        prepared.step(&state(0.0, 11.0), &[0.0], 0.0, limits(100.0, 100.0), || {
            false
        }),
        Err(RotorVoltageError::CurrentCompliance { .. })
    ));
    assert!(matches!(
        prepared.step(
            &state(0.0, f64::NAN),
            &[0.0],
            0.0,
            limits(100.0, 100.0),
            || false
        ),
        Err(RotorVoltageError::Invalid(_))
    ));
    assert!(matches!(
        prepared.step(
            &state(101.0, 0.0),
            &[0.0],
            0.0,
            limits(100.0, 100.0),
            || false
        ),
        Err(RotorVoltageError::Drive(
            RotorDriveError::SpeedDomain { .. }
        ))
    ));

    let active = RotorBearingModel::new(1, &[1.0], &[-4.0], &[0.0], &[0.0], || false).unwrap();
    let active_step = active
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    assert!(matches!(
        active_step.voltage_motor(
            RotorWindingMotor {
                inductance_h: 0.25,
                ..motor()
            },
            || false
        ),
        Err(RotorVoltageError::SingularCoupling)
    ));

    let tiny_step = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1e-16, 1e-12, || false)
        .unwrap();
    let huge_current = RotorWindingMotor {
        drive: RotorCurrentDrive {
            winding_resistance_ohm: 1.0,
            maximum_current_a: 1e20,
            ..motor().drive
        },
        ..motor()
    };
    let prepared = tiny_step.voltage_motor(huge_current, || false).unwrap();
    let input = state(0.0, 1e16);
    // The finite predicted change cannot be stored beside this huge unchanged
    // current. Its actual endpoint must fail the original winding equation.
    assert!(matches!(
        prepared.step(&input, &[0.0], 0.0, limits(0.0, 0.0), || false),
        Err(RotorVoltageError::WindingBalance { .. })
    ));
    assert_eq!(input, state(0.0, 1e16));
}

#[test]
fn preparation_and_every_observed_step_cancellation_preserve_both_physical_states() {
    let model = RotorBearingModel::new(1, &[1.3], &[0.07], &[2.1], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.3, 1e-12, || false)
        .unwrap();
    let parameters = RotorWindingMotor {
        drive: RotorCurrentDrive {
            winding_resistance_ohm: 0.2,
            ..motor().drive
        },
        inductance_h: 0.4,
    };
    let preparations = Cell::new(0);
    let prepared = stepper
        .voltage_motor(parameters, || {
            preparations.set(preparations.get() + 1);
            false
        })
        .unwrap();
    for stop in 1..=preparations.get() {
        let count = Cell::new(0);
        assert!(matches!(
            stepper.voltage_motor(parameters, || {
                count.set(count.get() + 1);
                count.get() == stop
            }),
            Err(RotorVoltageError::Drive(RotorDriveError::Power(
                RotorPowerError::Dynamics(RotorError::Cancelled)
            )))
        ));
    }
    let input = state(0.37, 0.27);
    let saved = input.clone();
    let calls = Cell::new(0);
    let accepted = prepared
        .step(&input, &[-0.2], 9.0, limits(0.43, 0.2), || {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    for stop in 1..=calls.get() {
        let count = Cell::new(0);
        assert!(matches!(
            prepared.step(&input, &[-0.2], 9.0, limits(0.43, 0.2), || {
                count.set(count.get() + 1);
                count.get() == stop
            }),
            Err(RotorVoltageError::Drive(RotorDriveError::Power(
                RotorPowerError::Dynamics(RotorError::Cancelled)
            )))
        ));
        assert_eq!(input, saved);
    }
    assert_eq!(
        prepared
            .step(&input, &[-0.2], 9.0, limits(0.43, 0.2), || false)
            .unwrap(),
        accepted
    );
}

#[test]
fn resistance_rebind_preserves_weak_coupling_matches_fresh_preparation_and_rolls_back() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.3, 1e-12, || false)
        .unwrap();
    let huge_resistance = RotorWindingMotor {
        drive: RotorCurrentDrive {
            winding_resistance_ohm: 1e200,
            ..motor().drive
        },
        ..motor()
    };
    let prepared = stepper.voltage_motor(huge_resistance, || false).unwrap();
    assert_eq!(prepared.motor(), huge_resistance);
    assert_eq!(prepared.duration_s(), 0.3);
    let input = state(0.4, 0.2);
    let saved = input.clone();
    // The nonbinary duration makes h*R_old carry a large multiplication
    // residual. Subtracting it from the total feedback would erase the weak
    // electromagnetic term even with double-double arithmetic.
    for resistance in [0.0, 0.3] {
        let rebound = prepared
            .with_winding_resistance(resistance, || false)
            .unwrap();
        let fresh_parameters = RotorWindingMotor {
            drive: RotorCurrentDrive {
                winding_resistance_ohm: resistance,
                ..huge_resistance.drive
            },
            ..huge_resistance
        };
        let fresh = stepper.voltage_motor(fresh_parameters, || false).unwrap();
        assert_eq!(rebound.motor(), fresh_parameters);
        assert_eq!(rebound.duration_s(), fresh.duration_s());
        assert_eq!(
            rebound
                .step(&input, &[-0.1], 2.0, limits(1000.0, 1000.0), || false)
                .unwrap(),
            fresh
                .step(&input, &[-0.1], 2.0, limits(1000.0, 1000.0), || false)
                .unwrap()
        );
    }
    let reference = prepared.with_winding_resistance(0.3, || false).unwrap();
    let baseline = reference
        .step(&input, &[-0.1], 2.0, limits(1000.0, 1000.0), || false)
        .unwrap();
    for resistance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            reference.with_winding_resistance(resistance, || false),
            Err(RotorVoltageError::Invalid(_))
        ));
    }
    let calls = Cell::new(0);
    reference
        .with_winding_resistance(0.7, || {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    for stop in 1..=calls.get() {
        let count = Cell::new(0);
        assert!(matches!(
            reference.with_winding_resistance(0.7, || {
                count.set(count.get() + 1);
                count.get() == stop
            }),
            Err(RotorVoltageError::Drive(RotorDriveError::Power(
                RotorPowerError::Dynamics(RotorError::Cancelled)
            )))
        ));
    }
    assert_eq!(
        reference
            .step(&input, &[-0.1], 2.0, limits(1000.0, 1000.0), || false)
            .unwrap(),
        baseline
    );
    assert_eq!(prepared.motor(), huge_resistance);
    assert_eq!(input, saved);

    let active = RotorBearingModel::new(1, &[1.0], &[-4.0], &[0.0], &[0.0], || false).unwrap();
    let active_step = active
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let active_motor = active_step
        .voltage_motor(
            RotorWindingMotor {
                drive: RotorCurrentDrive {
                    winding_resistance_ohm: 1.0,
                    ..motor().drive
                },
                inductance_h: 0.125,
            },
            || false,
        )
        .unwrap();
    assert!(matches!(
        active_motor.with_winding_resistance(0.25, || false),
        Err(RotorVoltageError::SingularCoupling)
    ));
    assert_eq!(active_motor.motor().drive.winding_resistance_ohm, 1.0);
}
