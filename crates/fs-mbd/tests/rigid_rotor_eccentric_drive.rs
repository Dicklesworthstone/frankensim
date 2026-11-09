//! G1/G3/G4: finite electrical sources on the actual eccentric full body.
use fs_ga::{Mat3, Se3, Vec3 as GaVec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::rigid::{
    EccentricRotorCurrentCommand, EccentricRotorDriveControl, EccentricRotorDriveError,
    ElasticEccentricRotor, RadialBearing, RotorSpinDriveAxis,
};
use fs_mbd::rotordynamics::transient::drive::{RotorCurrentDrive, RotorDrivePowerControl};
use fs_mbd::{Gravity, RigidBodyIntegrator, RigidBodyState, Vec3};
use std::cell::Cell;

fn rotor(tensor: [f64; 9], bearings: &[RadialBearing]) -> ElasticEccentricRotor {
    ElasticEccentricRotor::new(
        SpatialInertia::new(3.0, GaVec3::new(0.25, 0.0, 0.0), Mat3 { m: tensor }).unwrap(),
        Se3::identity(),
        bearings,
        || false,
    )
    .unwrap()
}
fn planar(bearings: &[RadialBearing]) -> ElasticEccentricRotor {
    rotor([1.2, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, 0.8], bearings)
}
fn bearing() -> RadialBearing {
    RadialBearing {
        station_m: 0.0,
        stiffness_n_per_m: [[30.0, 5.0], [5.0, 18.0]],
        damping_ns_per_m: [[0.0; 2]; 2],
    }
}
fn integrator() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn state(rotor: &ElasticEccentricRotor, speed: f64) -> RigidBodyState {
    let l = rotor
        .mass_properties()
        .angular_momentum_body_checked(Vec3::new(0.0, 0.0, speed))
        .unwrap();
    RigidBodyState::new(rotor.reference_pose(), Vec3::ZERO, l).unwrap()
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
fn control(h: f64) -> EccentricRotorDriveControl {
    EccentricRotorDriveControl {
        duration_seconds: h,
        relative_tolerance: 1e-11,
        electrical: RotorDrivePowerControl {
            maximum_source_input_w: 100.0,
            maximum_source_return_w: 100.0,
            absolute_power_tolerance_w: 1e-9,
            absolute_voltage_tolerance_v: 1e-12,
        },
        angular_impulse_tolerance_nms: 1e-12,
        shaft_rate_tolerance_rad_s: 1e-11,
        total_energy_tolerance_j: 0.01,
    }
}
fn command(current: f64) -> EccentricRotorCurrentCommand {
    EccentricRotorCurrentCommand {
        current_a: current,
        external_shaft_torque_nm: 0.0,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_nonprincipal_source_uses_original_mobility_and_complete_body_work() {
    let rotor = rotor([2.0, 0.0, 0.5, 0.0, 2.0, 0.0, 0.5, 0.0, 1.0], &[]);
    let motor = rotor
        .current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, || false)
        .unwrap();
    // Inverse of the original x/z block: I^-1_zz=8/7. Using I_zz=1
    // would give the wrong source acceleration and erase cross-axis work.
    near(motor.source_inertia_kg_m2(), 7.0 / 8.0, 2e-14);
    let before = RigidBodyState::new(
        rotor.reference_pose(),
        Vec3::new(30.0, -60.0, 9.0),
        Vec3::new(0.4, -0.2, 0.7),
    )
    .unwrap();
    let step = motor
        .step(integrator(), before, command(0.6), control(0.1), || false)
        .unwrap();
    let source = &step.sources[0];
    let after = source.impulse.state_after;
    assert_eq!(after.pose(), before.pose());
    assert_eq!(
        after.linear_momentum_world(),
        before.linear_momentum_world()
    );
    near(after.angular_momentum_body().x, 0.4, 0.0);
    near(after.angular_momentum_body().y, -0.2, 0.0);
    near(after.angular_momentum_body().z, 0.73, 2e-15);
    let omega0 = rotor
        .mass_properties()
        .angular_velocity_body_checked(before.angular_momentum_body())
        .unwrap();
    let omega1 = rotor
        .mass_properties()
        .angular_velocity_body_checked(after.angular_momentum_body())
        .unwrap();
    near(omega0.z, 1.2 / 1.75, 2e-14);
    near(omega1.x - omega0.x, -0.5 * 0.03 / 1.75, 2e-14);
    near(omega1.z - omega0.z, 0.03 * 8.0 / 7.0, 2e-14);
    let expected_work = (1.2 / 1.75) * 0.03 + 0.5 * (8.0 / 7.0) * 0.03 * 0.03;
    near(
        source.impulse.work.kinetic_energy_change_j,
        expected_work,
        3e-14,
    );
    near(source.electrical.shaft_work_j, expected_work, 3e-14);
    near(source.terminal_voltage_after_v, omega1.z, 3e-14);
    near(source.source_power_after_w, omega1.z * 0.6 / 0.8, 3e-14);
    near(source.source_work_residual_j, 0.0, 3e-14);
    near(
        step.kinetic_energy_change_j,
        before
            .kinetic_energy_change_to(rotor.mass_properties(), step.state_after)
            .unwrap(),
        3e-14,
    );
}

#[test]
fn g1_power_limited_startup_uses_actual_back_emf_between_sources() {
    let rotor = planar(&[]);
    let motor = rotor
        .current_spin_motor(
            RotorCurrentDrive {
                winding_resistance_ohm: 0.5,
                ..drive()
            },
            RotorSpinDriveAxis::Shaft,
            || false,
        )
        .unwrap();
    let mut limits = control(0.3);
    limits.electrical.maximum_source_input_w = 12.0;
    let step = motor
        .step(
            integrator(),
            state(&rotor, 0.0),
            command(8.0),
            limits,
            || false,
        )
        .unwrap();
    // The original terminal law at each source endpoint is
    // (.5 + .15/.8) I² + omega_before I = .8 * 12.
    let a: f64 = 0.5 + 0.15 / 0.8;
    let first_current = (0.8 * 12.0 / a).sqrt();
    let middle_speed = 0.15 * first_current / 0.8;
    let last_current = 2.0 * 0.8 * 12.0
        / (middle_speed + (middle_speed * middle_speed + 4.0 * a * 0.8 * 12.0).sqrt());
    near(step.sources[0].electrical.current_a, first_current, 2e-13);
    near(step.sources[1].electrical.current_a, last_current, 2e-13);
    near(
        step.sources[1].electrical.initial_terminal_voltage_v,
        middle_speed + 0.5 * last_current,
        2e-13,
    );
    near(
        step.sources[1].shaft_rate_after_rad_s,
        0.15 * (first_current + last_current) / 0.8,
        2e-13,
    );
    for source in &step.sources {
        assert!(source.electrical.initial_source_power_w <= 12.0 + 1e-10);
        near(source.electrical.final_source_power_w, 12.0, 2e-12);
        assert!(source.source_power_after_w <= 12.0 + 1e-10);
        near(
            source.terminal_voltage_after_v,
            source.shaft_rate_after_rad_s + 0.5 * source.electrical.current_a,
            2e-13,
        );
    }
    let axis = step
        .state_after
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    near(axis.x, (0.3 * middle_speed).cos(), 2e-14);
    near(axis.y, (0.3 * middle_speed).sin(), 2e-14);
    near(
        step.winding_heat_j,
        0.5 * 0.15 * (first_current.powi(2) + last_current.powi(2)),
        2e-12,
    );
    assert!(step.bridge_heat_j > 0.0 && step.electrical_draw_j > 0.0);
    near(step.total_energy_residual_j, 0.0, 3e-13);
}

#[test]
fn g1_stall_heating_and_regeneration_remain_distinct_from_external_shaft_work() {
    let rotor = planar(&[]);
    let motor = rotor
        .current_spin_motor(
            RotorCurrentDrive {
                winding_resistance_ohm: 1.0,
                ..drive()
            },
            RotorSpinDriveAxis::Shaft,
            || false,
        )
        .unwrap();
    let mut load = command(2.0);
    load.external_shaft_torque_nm = -2.0;
    let held = motor
        .step(integrator(), state(&rotor, 0.0), load, control(0.1), || {
            false
        })
        .unwrap();
    near(held.sources[1].shaft_rate_after_rad_s, 0.0, 2e-14);
    near(held.electrical_shaft_work_j, 0.0, 2e-14);
    near(held.mechanical.shaft_work_j, 0.0, 2e-14);
    near(held.winding_heat_j, 0.4, 2e-14);
    near(held.bridge_heat_j, 0.1, 2e-14);
    near(held.electrical_draw_j, 0.5, 2e-14);
    near(held.total_energy_residual_j, 0.0, 2e-14);
    let generator = rotor
        .current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, || false)
        .unwrap();
    let braking = generator
        .step(
            integrator(),
            state(&rotor, 3.0),
            command(-1.0),
            control(0.1),
            || false,
        )
        .unwrap();
    near(braking.electrical_shaft_work_j, -0.29375, 2e-13);
    near(braking.electrical_return_j, 0.146875, 2e-13);
    near(braking.bridge_heat_j, 0.146875, 2e-13);
    assert_eq!(braking.electrical_draw_j, 0.0);
    assert_eq!(braking.winding_heat_j, 0.0);
    near(braking.total_energy_residual_j, 0.0, 2e-13);
}

// Independent planar COM ODE; the electrical composition must retain the
// bearing's reciprocal shaft torque rather than merely add a radial force.
fn derivative(y: [f64; 6]) -> [f64; 6] {
    let (s, c) = y[2].sin_cos();
    let q = [y[0] - 0.25 * c, y[1] - 0.25 * s];
    let f = [-30.0 * q[0] - 5.0 * q[1], -5.0 * q[0] - 18.0 * q[1]];
    [
        y[3],
        y[4],
        y[5],
        f[0] / 3.0,
        f[1] / 3.0,
        (0.12 - 0.25 * (-s * f[0] + c * f[1])) / 0.8,
    ]
}
fn reference() -> [f64; 6] {
    let mut y = [0.25, 0.0, 0.0, 0.0, 0.0, 3.0];
    let h = 2.0 / 8000.0;
    for _ in 0..8000 {
        let k1 = derivative(y);
        let k2 = derivative(core::array::from_fn(|i| y[i] + h * k1[i] / 2.0));
        let k3 = derivative(core::array::from_fn(|i| y[i] + h * k2[i] / 2.0));
        let k4 = derivative(core::array::from_fn(|i| y[i] + h * k3[i]));
        y = core::array::from_fn(|i| y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}
fn trajectory(steps: usize, reference: [f64; 6]) -> f64 {
    let rotor = planar(&[bearing()]);
    let motor = rotor
        .current_spin_motor(
            RotorCurrentDrive {
                winding_resistance_ohm: 0.25,
                ..drive()
            },
            RotorSpinDriveAxis::Shaft,
            || false,
        )
        .unwrap();
    let mut current = state(&rotor, 3.0);
    let mut heat = 0.0;
    let h = 2.0 / steps as f64;
    for _ in 0..steps {
        let step = motor
            .step(integrator(), current, command(0.12), control(h), || false)
            .unwrap();
        let actual_middle = rotor
            .mass_properties()
            .angular_velocity_body_checked(
                step.mechanical.native.state_after.angular_momentum_body(),
            )
            .unwrap()
            .z;
        near(
            step.sources[1].electrical.initial_terminal_voltage_v,
            actual_middle + 0.25 * 0.12,
            2e-13,
        );
        near(step.sources[0].electrical.current_a, 0.12, 2e-14);
        near(step.sources[1].electrical.current_a, 0.12, 2e-14);
        heat += step.winding_heat_j;
        current = step.state_after;
    }
    near(heat, 0.25 * 0.12 * 0.12 * 2.0, 2e-13);
    let p = current.pose().position_world();
    let v = current
        .center_of_mass_velocity_world(rotor.mass_properties())
        .unwrap();
    let w = rotor
        .mass_properties()
        .angular_velocity_body_checked(current.angular_momentum_body())
        .unwrap()
        .z;
    assert!(
        w < 0.5,
        "bearing recoil must change the unloaded speed of 3.3 rad/s"
    );
    let axis = current
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    (p.x - reference[0]).abs()
        + (p.y - reference[1]).abs()
        + (v.x - reference[3]).abs()
        + (v.y - reference[4]).abs()
        + (w - reference[5]).abs()
        + (axis.x - reference[2].cos()).abs()
        + (axis.y - reference[2].sin()).abs()
}

#[test]
fn g3_eccentric_orbit_reacts_on_motor_speed_with_second_order_refinement() {
    let oracle = reference();
    let coarse = trajectory(50, oracle);
    let medium = trajectory(100, oracle);
    let fine = trajectory(200, oracle);
    assert!(
        coarse / medium > 3.8 && medium / fine > 3.8,
        "{coarse:e} {medium:e} {fine:e}"
    );
    assert!(fine < 1.5e-4, "{fine:e}");
}

#[test]
fn g3_zero_current_reuses_native_motion_and_final_cancellation_is_atomic() {
    let rotor = planar(&[bearing()]);
    let motor = rotor
        .current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, || false)
        .unwrap();
    let initial = state(&rotor, 3.0);
    let mut load = command(0.0);
    load.external_shaft_torque_nm = 0.1;
    let bare = rotor
        .step(integrator(), initial, 0.1, 0.01, || false)
        .unwrap();
    let zero = motor
        .step(integrator(), initial, load, control(0.01), || false)
        .unwrap();
    assert_eq!(zero.state_after, bare.native.state_after);
    assert_eq!(zero.mechanical, bare);
    assert_eq!(zero.electrical_shaft_work_j, 0.0);
    assert_eq!(zero.winding_heat_j, 0.0);
    let polls = Cell::new(0);
    let accepted = motor
        .step(integrator(), initial, command(0.12), control(0.01), || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let last = polls.get();
    polls.set(0);
    let refused = motor.step(integrator(), initial, command(0.12), control(0.01), || {
        polls.set(polls.get() + 1);
        polls.get() == last
    });
    assert_eq!(refused, Err(EccentricRotorDriveError::Cancelled));
    assert_eq!(
        motor
            .step(integrator(), initial, command(0.12), control(0.01), || {
                false
            })
            .unwrap(),
        accepted
    );
    assert!(
        rotor
            .current_spin_motor(drive(), RotorSpinDriveAxis::Stationary, || false)
            .is_err()
    );
    assert!(
        rotor
            .current_spin_motor(
                RotorCurrentDrive {
                    shaft_coordinate: 1,
                    ..drive()
                },
                RotorSpinDriveAxis::Shaft,
                || false
            )
            .is_err()
    );
    assert!(
        motor
            .step(
                integrator(),
                initial,
                command(0.12),
                control(f64::from_bits(1)),
                || false
            )
            .is_err()
    );
}

#[test]
fn g4_unresolved_native_angular_impulse_cannot_be_reported_as_motor_work() {
    let rotor = rotor([2.0, 0.0, 0.5, 0.0, 2.0, 0.0, 0.5, 0.0, 1.0], &[]);
    let large = 2.0_f64.powi(60);
    let initial = RigidBodyState::new(
        rotor.reference_pose(),
        Vec3::ZERO,
        Vec3::new(large, 0.0, 0.25 * large),
    )
    .unwrap();
    let motor = rotor
        .current_spin_motor(
            RotorCurrentDrive {
                maximum_terminal_voltage_v: 1e30,
                minimum_speed_rad_s: -1e30,
                maximum_speed_rad_s: 1e30,
                ..drive()
            },
            RotorSpinDriveAxis::Shaft,
            || false,
        )
        .unwrap();
    let mut limits = control(0.01);
    limits.electrical.maximum_source_input_w = 1e30;
    limits.electrical.maximum_source_return_w = 1e30;
    limits.electrical.absolute_power_tolerance_w = 1e10;
    // The scalar source resolves a small omega increment; the body's large
    // L_z cannot represent its 0.005 N m s impulse. The explicit native gate
    // must refuse even when a permissive electrical energy gate would pass.
    assert!(matches!(
        motor.step(integrator(), initial, command(1.0), limits, || false),
        Err(EccentricRotorDriveError::SourceMismatch { .. })
    ));
}
