//! Finite electrical spin acceleration, native gyro response and complete work.
use fs_ga::{Mat3, So3, Vec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::rigid::{
    RadialBearing, RigidRotor, RotorSpinCurrentCommand, RotorSpinDriveAxis, RotorSpinDriveControl,
    RotorSpinState,
};
use fs_mbd::rotordynamics::transient::RotorTransientState;
use fs_mbd::rotordynamics::transient::drive::{RotorCurrentDrive, RotorDrivePowerControl};
use std::cell::Cell;

fn rotor(bearings: &[RadialBearing]) -> RigidRotor {
    let inertia = SpatialInertia::new(
        4.0,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0],
        },
    )
    .unwrap();
    RigidRotor::new(inertia, So3::identity(), bearings, 1e-12, || false).unwrap()
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
fn control(h: f64) -> RotorSpinDriveControl {
    RotorSpinDriveControl {
        duration_seconds: h,
        relative_tolerance: 1e-11,
        electrical: RotorDrivePowerControl {
            maximum_source_input_w: 100.0,
            maximum_source_return_w: 100.0,
            absolute_power_tolerance_w: 1e-10,
            absolute_voltage_tolerance_v: 1e-12,
        },
    }
}
fn state() -> RotorTransientState {
    RotorTransientState {
        displacement: vec![0.0; 4],
        velocity: vec![0.0; 4],
    }
}
fn spin(speed: f64) -> RotorSpinState {
    RotorSpinState {
        time_seconds: 7.0,
        phase_rad: 0.3,
        spin_rad_s: speed,
    }
}
fn command(current: f64) -> RotorSpinCurrentCommand {
    RotorSpinCurrentCommand {
        current_a: current,
        shaft_load_torque_nm: 0.0,
        lateral_force: [0.0; 4],
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_actual_polar_inertia_and_electrical_peak_limit_set_startup_acceleration() {
    let rotor = rotor(&[]);
    let motor = rotor
        .current_spin_motor(
            RotorCurrentDrive {
                winding_resistance_ohm: 0.5,
                ..drive()
            },
            RotorSpinDriveAxis::Shaft,
            0.0,
            || false,
        )
        .unwrap();
    let mut limits = control(0.3);
    limits.electrical.maximum_source_input_w = 12.0;
    let step = motor
        .step(&state(), spin(0.0), command(8.0), limits, || false)
        .unwrap();
    // Jp=3, h=.3: Omega1=.1 I. Terminal power=(.1+.5)I²,
    // so source peak=.75 I² and the 12 W limit selects I=4 A.
    near(step.motor.current_a, 4.0, 2e-13);
    near(step.motor.current_fraction, 0.5, 2e-14);
    near(step.spin_after.spin_rad_s, 0.4, 2e-14);
    near(step.angular_acceleration_rad_s2, 4.0 / 3.0, 2e-13);
    near(step.spin_after.phase_rad, 0.36, 2e-15);
    near(step.motor.final_source_power_w, 12.0, 2e-12);
    near(step.motor.motion.kinetic_energy_change_j, 0.24, 2e-13);
    near(step.motor.winding_heat_j, 2.4, 2e-12);
    near(step.motor.bridge_heat_j, 0.66, 2e-12);
    near(step.motor.electrical_draw_j, 3.3, 2e-12);
    near(step.balance_residual_j, 0.0, 3e-13);
    assert_eq!(step.lateral.state_after, state());
}

fn accelerating_whirl(steps: usize) -> [f64; 2] {
    let rotor = rotor(&[]);
    let motor = rotor
        .current_spin_motor(
            RotorCurrentDrive {
                torque_constant_nm_per_a: 3.0,
                motoring_efficiency: 1.0,
                regenerative_efficiency: 1.0,
                ..drive()
            },
            RotorSpinDriveAxis::Shaft,
            0.0,
            || false,
        )
        .unwrap();
    let mut lateral = state();
    lateral.velocity[2] = 0.1;
    lateral.velocity[3] = 0.05;
    let mut spin = spin(0.0);
    for _ in 0..steps {
        let step = motor
            .step(
                &lateral,
                spin,
                command(4.0),
                control(0.8 / steps as f64),
                || false,
            )
            .unwrap();
        near(step.lateral.kinetic_energy_change_j, 0.0, 1e-16);
        near(step.lateral.gyroscopic_work_j, 0.0, 1e-16);
        near(step.balance_residual_j, 0.0, 2e-14);
        lateral = step.lateral.state_after;
        spin = step.spin_after;
    }
    near(spin.spin_rad_s, 3.2, 1e-13);
    near(spin.phase_rad, 0.3 + 1.28, 1e-13);
    [lateral.velocity[2], lateral.velocity[3]]
}

#[test]
fn g1_finite_motor_spin_acceleration_changes_whirl_with_second_order_refinement() {
    // For a free balanced rotor, theta_dot rotates by (Jp/Jd)*integral(Omega dt).
    // Constant finite torque Kt*I=12 gives alpha=4, so this oracle is analytic.
    let angle: f64 = 1.5 * 0.5 * 4.0 * 0.8 * 0.8;
    let (s, c) = angle.sin_cos();
    let exact = [0.1 * c - 0.05 * s, 0.1 * s + 0.05 * c];
    let error = |n| {
        let actual = accelerating_whirl(n);
        (actual[0] - exact[0]).hypot(actual[1] - exact[1])
    };
    let coarse = error(20);
    let medium = error(40);
    let fine = error(80);
    assert!(coarse / medium > 3.97, "{coarse:e}/{medium:e}");
    assert!(medium / fine > 3.99, "{medium:e}/{fine:e}");
    assert!(fine < 2.2e-5, "{fine:e}");
}

#[test]
fn g3_stall_reversal_regeneration_and_shaft_drag_keep_distinct_physical_accounts() {
    let rotor = rotor(&[]);
    let stalled = rotor
        .current_spin_motor(
            RotorCurrentDrive {
                winding_resistance_ohm: 1.0,
                ..drive()
            },
            RotorSpinDriveAxis::Shaft,
            0.0,
            || false,
        )
        .unwrap()
        .step(
            &state(),
            spin(0.0),
            RotorSpinCurrentCommand {
                shaft_load_torque_nm: -3.0,
                ..command(3.0)
            },
            control(0.25),
            || false,
        )
        .unwrap();
    assert_eq!(stalled.spin_after.spin_rad_s, 0.0);
    assert_eq!(stalled.motor.shaft_work_j, 0.0);
    near(stalled.motor.winding_heat_j, 2.25, 1e-14);
    near(stalled.motor.electrical_draw_j, 2.8125, 1e-14);
    let reversing = rotor
        .current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, 0.0, || false)
        .unwrap()
        .step(&state(), spin(-1.0), command(6.0), control(1.0), || false)
        .unwrap();
    near(reversing.spin_after.spin_rad_s, 1.0, 1e-14);
    near(reversing.motor.motion.kinetic_energy_change_j, 0.0, 1e-14);
    near(reversing.motor.electrical_draw_j, 1.875, 1e-14);
    near(reversing.motor.electrical_return_j, 0.75, 1e-14);
    near(reversing.motor.bridge_heat_j, 1.125, 1e-14);
    near(reversing.balance_residual_j, 0.0, 1e-14);
    let drag = rotor
        .current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, 3.0, || false)
        .unwrap()
        .step(&state(), spin(2.0), command(0.0), control(0.2), || false)
        .unwrap();
    near(drag.spin_after.spin_rad_s, 18.0 / 11.0, 1e-14);
    near(
        drag.motor.motion.velocity_force_work_j,
        0.5 * 3.0 * (4.0 - (18.0_f64 / 11.0).powi(2)),
        2e-14,
    );
    assert_eq!(drag.motor.electrical_draw_j, 0.0);
    near(drag.balance_residual_j, 0.0, 2e-14);
}

#[test]
fn g3_original_bearing_operator_is_reused_and_refusals_publish_no_partial_step() {
    let rotor = rotor(&[RadialBearing {
        station_m: 0.4,
        stiffness_n_per_m: [[17.0, -3.0], [2.0, 23.0]],
        damping_ns_per_m: [[0.4, -0.3], [0.2, -0.1]],
    }]);
    assert!(
        rotor
            .current_spin_motor(drive(), RotorSpinDriveAxis::Stationary, 0.0, || false)
            .is_err()
    );
    assert!(
        rotor
            .current_spin_motor(
                RotorCurrentDrive {
                    shaft_coordinate: 2,
                    ..drive()
                },
                RotorSpinDriveAxis::Shaft,
                0.0,
                || false
            )
            .is_err()
    );
    assert!(
        rotor
            .current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, -1.0, || false)
            .is_err()
    );
    let motor = rotor
        .current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, 0.2, || false)
        .unwrap();
    let before = RotorTransientState {
        displacement: vec![0.02, -0.01, 0.03, 0.04],
        velocity: vec![-0.03, 0.05, 0.02, -0.01],
    };
    let initial_spin = spin(2.0);
    let load = RotorSpinCurrentCommand {
        lateral_force: [0.3, -0.1, 0.2, 0.05],
        ..command(1.0)
    };
    let calls = Cell::new(0);
    let expected = motor
        .step(&before, initial_spin, load, control(0.02), || {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    let native = rotor
        .model()
        .at_speed(expected.midpoint_spin_rad_s, || false)
        .unwrap()
        .midpoint_stepper(0.02, 1e-11, || false)
        .unwrap()
        .step(&before, &load.lateral_force, || false)
        .unwrap();
    assert_eq!(expected.lateral, native);
    for stop in [1, calls.get()] {
        let polls = Cell::new(0);
        assert!(
            motor
                .step(&before, initial_spin, load, control(0.02), || {
                    polls.set(polls.get() + 1);
                    polls.get() == stop
                })
                .is_err()
        );
        assert_eq!(
            motor
                .step(&before, initial_spin, load, control(0.02), || false)
                .unwrap(),
            expected
        );
    }
    let mut invalid = before.clone();
    invalid.velocity.pop();
    assert!(
        motor
            .step(&invalid, initial_spin, load, control(0.02), || false)
            .is_err()
    );
    for bad in [0.0, -0.1, f64::INFINITY] {
        assert!(
            motor
                .step(&before, initial_spin, load, control(bad), || false)
                .is_err()
        );
    }
    for speed in [f64::NAN, 101.0] {
        assert!(
            motor
                .step(&before, spin(speed), load, control(0.02), || false)
                .is_err()
        );
    }
    assert!(
        motor
            .step(
                &before,
                RotorSpinState {
                    time_seconds: 1e100,
                    ..initial_spin
                },
                load,
                control(0.02),
                || false
            )
            .is_err()
    );
    assert_eq!(
        motor
            .step(&before, initial_spin, load, control(0.02), || false)
            .unwrap(),
        expected
    );
}
