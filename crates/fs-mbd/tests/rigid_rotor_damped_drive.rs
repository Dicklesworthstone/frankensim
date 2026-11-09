//! G1/G3/G4: electrical work and passive heat on one eccentric native body.
use fs_ga::{Mat3, Se3, Vec3 as GaVec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::rigid::{
    DampedEccentricCurrentMotor, DampedEccentricDriveControl, DampedEccentricDriveError,
    DampedEccentricRotor, EccentricRotorCurrentCommand, EccentricRotorDampingControl,
    EccentricRotorDriveControl, EccentricRotorDriveError, RadialBearing, RotorSpinDriveAxis,
};
use fs_mbd::rotordynamics::transient::drive::{
    RotorCurrentDrive, RotorDriveError, RotorDrivePowerControl,
};
use fs_mbd::{Gravity, RigidBodyIntegrator, RigidBodyState, Vec3};
use std::cell::Cell;

const MASS: f64 = 3.0;
const E: f64 = 0.25;
const POLAR: f64 = 0.8;

fn bearing() -> RadialBearing {
    RadialBearing {
        station_m: 0.0,
        stiffness_n_per_m: [[30.0, 5.0], [5.0, 18.0]],
        damping_ns_per_m: [[0.8, 0.1], [0.1, 0.5]],
    }
}
fn rotor(bearings: &[RadialBearing]) -> DampedEccentricRotor {
    DampedEccentricRotor::new(
        SpatialInertia::new(
            MASS,
            GaVec3::new(E, 0.0, 0.0),
            Mat3 {
                m: [1.2, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, POLAR],
            },
        )
        .unwrap(),
        Se3::identity(),
        bearings,
        || false,
    )
    .unwrap()
}
fn drive() -> RotorCurrentDrive {
    RotorCurrentDrive {
        shaft_coordinate: 0,
        torque_constant_nm_per_a: 1.0,
        winding_resistance_ohm: 0.4,
        maximum_current_a: 10.0,
        maximum_terminal_voltage_v: 100.0,
        minimum_speed_rad_s: -100.0,
        maximum_speed_rad_s: 100.0,
        motoring_efficiency: 0.8,
        regenerative_efficiency: 0.5,
    }
}
fn motor(rotor: &DampedEccentricRotor) -> DampedEccentricCurrentMotor<'_> {
    rotor.current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, || false).unwrap()
}
fn control(h: f64) -> DampedEccentricDriveControl {
    DampedEccentricDriveControl {
        drive: EccentricRotorDriveControl {
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
        },
        damping: EccentricRotorDampingControl {
            relative_tolerance: 1e-11,
            maximum_modal_step: 1.0,
            linear_impulse_tolerance_ns: 1e-12,
            angular_impulse_tolerance_nms: 1e-12,
            linear_speed_tolerance_m_s: 1e-11,
            angular_speed_tolerance_rad_s: 1e-11,
            source_energy_tolerance_j: 1e-10,
            total_energy_tolerance_j: 0.01,
        },
    }
}
fn command(current_a: f64, torque: f64) -> EccentricRotorCurrentCommand {
    EccentricRotorCurrentCommand { current_a, external_shaft_torque_nm: torque }
}
fn integrator() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn state(rotor: &DampedEccentricRotor, velocity: Vec3, speed: f64) -> RigidBodyState {
    let original = rotor.elastic_component();
    RigidBodyState::new(
        original.reference_pose(),
        velocity.scale(MASS),
        original.mass_properties().angular_momentum_body_checked(Vec3::new(0.0, 0.0, speed)).unwrap(),
    )
    .unwrap()
}
fn rates(rotor: &DampedEccentricRotor, state: RigidBodyState) -> [f64; 3] {
    let props = rotor.elastic_component().mass_properties();
    let v = state.center_of_mass_velocity_world(props).unwrap();
    let w = props.angular_velocity_body_checked(state.angular_momentum_body()).unwrap();
    [v.x, v.y, w.z]
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * (1.0 + b.abs()), "{a:e} != {b:e}");
}

// Independent continuous COM equations and separate physical heat/source integrals.
// State: X,Y,phi,Vx,Vy,omega,bearing_heat,draw,return,winding_heat,bridge_heat.
fn derivative(y: [f64; 11], current: f64, external: f64) -> [f64; 11] {
    let (s, c) = y[2].sin_cos();
    let q = [y[0] - E * c, y[1] - E * s];
    let u = [y[3] + E * y[5] * s, y[4] - E * y[5] * c];
    let force = [
        -30.0 * q[0] - 5.0 * q[1] - 0.8 * u[0] - 0.1 * u[1],
        -5.0 * q[0] - 18.0 * q[1] - 0.1 * u[0] - 0.5 * u[1],
    ];
    let terminal = (y[5] + 0.4 * current) * current;
    let source = if terminal >= 0.0 { terminal / 0.8 } else { terminal * 0.5 };
    [
        y[3], y[4], y[5], force[0] / MASS, force[1] / MASS,
        (current + external - E * (-s * force[0] + c * force[1])) / POLAR,
        u[0] * (0.8 * u[0] + 0.1 * u[1]) + u[1] * (0.1 * u[0] + 0.5 * u[1]),
        source.max(0.0), (-source).max(0.0), 0.4 * current * current, source - terminal,
    ]
}
fn reference(current: f64, external: f64) -> [f64; 11] {
    let mut y = [E, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let h = 2.0 / 8000.0;
    for _ in 0..8000 {
        let k1 = derivative(y, current, external);
        let k2 = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]), current, external);
        let k3 = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]), current, external);
        let k4 = derivative(core::array::from_fn(|i| y[i] + h * k3[i]), current, external);
        y = core::array::from_fn(|i| y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}
fn trajectory(steps: usize, oracle: [f64; 11]) -> (f64, f64, f64, f64) {
    let rotor = rotor(&[bearing()]);
    let motor = motor(&rotor);
    let mut state = state(&rotor, Vec3::ZERO, 3.0);
    let mut totals = [0.0; 5];
    let mut residual = 0.0;
    for _ in 0..steps {
        let receipt = motor.step(
            integrator(), state, command(0.4, -0.07), control(2.0 / steps as f64), || false,
        ).unwrap();
        for source in &receipt.current.sources {
            near(source.electrical.current_a, 0.4, 1e-14);
        }
        let energy = [
            receipt.bearing_heat_j, receipt.current.electrical_draw_j,
            receipt.current.electrical_return_j, receipt.current.winding_heat_j,
            receipt.current.bridge_heat_j,
        ];
        for i in 0..5 {
            assert!(energy[i] >= 0.0);
            totals[i] += energy[i];
        }
        residual += receipt.total_energy_residual_j;
        state = receipt.state_after;
    }
    let p = state.pose().position_world();
    let v = rates(&rotor, state);
    let axis = state.pose().orientation().rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    let error = (p.x - oracle[0]).abs() + (p.y - oracle[1]).abs()
        + (v[0] - oracle[3]).abs() + (v[1] - oracle[4]).abs() + (v[2] - oracle[5]).abs()
        + (axis.x - oracle[2].cos()).abs() + (axis.y - oracle[2].sin()).abs();
    near(totals[3], 2.0 * 0.4 * 0.4 * 0.4, 2e-13);
    let electrical_error: f64 = (1..5).map(|i| (totals[i] - oracle[i + 6]).abs()).sum();
    (error, (totals[0] - oracle[6]).abs(), electrical_error, residual.abs())
}

#[test]
fn g1_g3_powered_damped_acceleration_and_independent_energy_integrals_refine() {
    let rotor = rotor(&[bearing()]);
    let motor = motor(&rotor);
    let initial = state(&rotor, Vec3::ZERO, 3.0);
    let exact = derivative(
        [E, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.4, -0.07,
    );
    let small = motor.step(
        integrator(), initial, command(0.4, -0.07), control(1e-6), || false,
    ).unwrap();
    let v = rates(&rotor, small.state_after);
    near(v[0] / 1e-6, exact[3], 3e-5);
    near(v[1] / 1e-6, exact[4], 3e-5);
    near((v[2] - 3.0) / 1e-6, exact[5], 3e-5);
    let point = rotor.prepare_motion(integrator(), initial, 0.4 - 0.07, || false)
        .unwrap().point_motion(Vec3::new(-E, 0.0, 0.0)).unwrap();
    near(point.point_acceleration_world.x, exact[3] + E * 9.0, 3e-14);
    near(point.point_acceleration_world.y, exact[4] - E * exact[5], 3e-14);

    let oracle = reference(0.4, -0.07);
    let a = trajectory(50, oracle);
    let b = trajectory(100, oracle);
    let c = trajectory(200, oracle);
    assert!(a.0 / b.0 > 3.5 && b.0 / c.0 > 3.5, "{a:?} {b:?} {c:?}");
    assert!(a.1 / b.1 > 3.5 && b.1 / c.1 > 3.5, "{a:?} {b:?} {c:?}");
    assert!(a.2 / b.2 > 3.5 && b.2 / c.2 > 3.5, "{a:?} {b:?} {c:?}");
    assert!(c.0 < 5e-4 && c.1 < 1e-4 && c.2 < 1e-3 && c.3 < 1e-3, "{c:?}");
}

#[test]
fn g1_power_curtailment_and_regeneration_keep_three_physical_losses_separate() {
    let rotor = rotor(&[RadialBearing {
        stiffness_n_per_m: [[0.0; 2]; 2],
        damping_ns_per_m: [[0.0, 0.0], [0.0, 0.2]],
        ..bearing()
    }]);
    let motor = motor(&rotor);
    let initial = state(&rotor, Vec3::ZERO, 3.0);
    for requested in [2.0, -0.5] {
        let mut limits = control(0.05);
        limits.drive.electrical.maximum_source_input_w = 2.0;
        let step = motor.step(integrator(), initial, command(requested, 0.0), limits, || false).unwrap();
        let current = &step.current;
        let winding = 0.5 * 0.05 * 0.4 * current.sources.iter()
            .map(|s| s.electrical.current_a.powi(2)).sum::<f64>();
        near(current.winding_heat_j, winding, 3e-14);
        assert!(step.bearing_heat_j > 0.0 && current.winding_heat_j > 0.0 && current.bridge_heat_j > 0.0);
        if requested > 0.0 {
            assert!(current.sources.iter().all(|s| s.electrical.current_fraction < 1.0));
            assert!(current.electrical_draw_j > 0.0);
            near(current.electrical_return_j, 0.0, 0.0);
        } else {
            assert!(current.electrical_return_j > 0.0);
            near(current.electrical_draw_j, 0.0, 0.0);
        }
        let storage = initial.kinetic_energy_change_to(
            rotor.elastic_component().mass_properties(), step.state_after,
        ).unwrap();
        let balance = storage + step.bearing_heat_j + current.winding_heat_j
            + current.bridge_heat_j + current.electrical_return_j - current.electrical_draw_j;
        near(step.total_energy_residual_j, balance, 3e-14);
        near(balance, 0.0, 2e-12);
        for endpoint in step.electrical_endpoints {
            let voltage = endpoint.shaft_rate_rad_s + 0.4 * endpoint.current_a;
            near(endpoint.terminal_voltage_v, voltage, 3e-14);
            let power = voltage * endpoint.current_a;
            near(endpoint.source_power_w, if power >= 0.0 { power / 0.8 } else { power * 0.5 }, 3e-14);
            assert!(endpoint.source_power_w <= 2.0 + 1e-9);
        }
    }
}

#[test]
fn g4_passive_spin_up_requires_fresh_final_speed_voltage_and_both_power_gates() {
    let rotor = rotor(&[RadialBearing {
        stiffness_n_per_m: [[0.0; 2]; 2],
        damping_ns_per_m: [[0.0, 0.0], [0.0, 4.0]],
        ..bearing()
    }]);
    let original = rotor.elastic_component();
    let initial = state(&rotor, Vec3::new(0.0, 2.0, 0.0), 0.1);
    let model = RotorCurrentDrive {
        winding_resistance_ohm: 0.0,
        motoring_efficiency: 1.0,
        regenerative_efficiency: 1.0,
        ..drive()
    };
    for current in [0.1, -0.1] {
        let request = command(current, 0.0);
        let base = rotor.current_spin_motor(model, RotorSpinDriveAxis::Shaft, || false).unwrap()
            .step(integrator(), initial, request, control(0.1), || false).unwrap();
        let middle_speed = rates(&rotor, base.current.state_after)[2];
        let final_speed = rates(&rotor, base.state_after)[2];
        assert!(final_speed > middle_speed + 0.05);
        assert!(base.damping_sources[1].impulse.work.kinetic_energy_change_j < 0.0);
        let cap = 0.5 * (middle_speed + final_speed);
        for gate in 0..3 {
            let mut restricted = model;
            let mut limits = control(0.1);
            if gate == 0 {
                restricted.maximum_speed_rad_s = cap;
            } else if gate == 1 {
                restricted.maximum_terminal_voltage_v = cap;
            } else if current > 0.0 {
                limits.drive.electrical.maximum_source_input_w = cap * current;
            } else {
                limits.drive.electrical.maximum_source_return_w = -cap * current;
            }
            // Every existing middle-stage limit is genuinely satisfied. Only the
            // actual final passive source exposes the new physical incompatibility.
            original.current_spin_motor(restricted, RotorSpinDriveAxis::Shaft, || false).unwrap()
                .step(integrator(), base.damping_sources[0].impulse.state_after, request, limits.drive, || false).unwrap();
            let rejected = rotor.current_spin_motor(restricted, RotorSpinDriveAxis::Shaft, || false).unwrap()
                .step(integrator(), initial, request, limits, || false);
            match (gate, rejected) {
                (0, Err(DampedEccentricDriveError::Drive(EccentricRotorDriveError::Electrical(RotorDriveError::SpeedDomain { .. })))) => {}
                (1, Err(DampedEccentricDriveError::Drive(EccentricRotorDriveError::Electrical(RotorDriveError::VoltageCompliance { .. })))) => {}
                (2, Err(DampedEccentricDriveError::Drive(EccentricRotorDriveError::Electrical(RotorDriveError::PowerBalance { .. })))) => {}
                (_, result) => panic!("unexpected endpoint admission: {result:?}"),
            }
        }
    }
}

#[test]
fn g4_zero_limits_reduce_to_original_owners_and_late_failure_is_atomic() {
    let zero = rotor(&[RadialBearing { damping_ns_per_m: [[0.0; 2]; 2], ..bearing() }]);
    let initial = state(&zero, Vec3::ZERO, 3.0);
    let original = zero.elastic_component().current_spin_motor(drive(), RotorSpinDriveAxis::Shaft, || false).unwrap()
        .step(integrator(), initial, command(0.4, -0.07), control(0.01).drive, || false).unwrap();
    let combined = motor(&zero).step(integrator(), initial, command(0.4, -0.07), control(0.01), || false).unwrap();
    assert_eq!(combined.state_after, original.state_after);
    assert_eq!(combined.current, original);
    assert_eq!(combined.bearing_heat_j, 0.0);

    let rotor = rotor(&[bearing()]);
    let motor = motor(&rotor);
    let initial = state(&rotor, Vec3::ZERO, 3.0);
    let passive = rotor.step(integrator(), initial, -0.07, 0.01, control(0.01).damping, || false).unwrap();
    let combined = motor.step(integrator(), initial, command(0.0, -0.07), control(0.01), || false).unwrap();
    assert_eq!(combined.state_after, passive.state_after);
    assert_eq!(combined.bearing_heat_j, passive.dissipated_energy_j);
    assert_eq!(combined.current.winding_heat_j + combined.current.bridge_heat_j, 0.0);

    let polls = Cell::new(0);
    let expected = motor.step(integrator(), initial, command(0.4, -0.07), control(0.04), || {
        polls.set(polls.get() + 1); false
    }).unwrap();
    let last = polls.get();
    polls.set(0);
    assert_eq!(motor.step(integrator(), initial, command(0.4, -0.07), control(0.04), || {
        polls.set(polls.get() + 1); polls.get() == last
    }), Err(DampedEccentricDriveError::Cancelled));
    let mut tight = control(0.04);
    tight.damping.total_energy_tolerance_j = 1e-12;
    assert!(matches!(motor.step(integrator(), initial, command(0.4, -0.07), tight, || false),
        Err(DampedEccentricDriveError::EnergyBalance { .. })));
    assert_eq!(motor.step(integrator(), initial, command(0.4, -0.07), control(0.04), || false).unwrap(), expected);
}
