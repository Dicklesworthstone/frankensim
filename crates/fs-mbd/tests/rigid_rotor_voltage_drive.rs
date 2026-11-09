//! G1/G3/G4: continuing inductive current on one original eccentric rotor.
use fs_ga::{Mat3, Se3, Vec3 as GaVec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::rigid::{
    DampedEccentricRotor, EccentricRotorDampingControl, EccentricRotorVoltageCommand,
    EccentricRotorVoltageControl, EccentricRotorVoltageError, EccentricRotorWindingState,
    RadialBearing, RotorSpinDriveAxis,
};
use fs_mbd::rotordynamics::transient::drive::{
    RotorCurrentDrive, RotorDriveError, RotorDrivePowerControl,
    voltage::{RotorVoltageControl, RotorVoltageError, RotorWindingMotor},
};
use fs_mbd::{Gravity, RigidBodyIntegrator, RigidBodyState, Vec3};
use std::cell::Cell;

const MASS: f64 = 3.0;
const E: f64 = 0.25;
const POLAR: f64 = 0.8;
const INDUCTANCE: f64 = 0.6;
fn bearing() -> RadialBearing {
    RadialBearing {
        station_m: 0.0,
        stiffness_n_per_m: [[30.0, 5.0], [5.0, 18.0]],
        damping_ns_per_m: [[0.8, 0.1], [0.1, 0.5]],
    }
}
fn full(tensor: [f64; 9], bearings: &[RadialBearing]) -> DampedEccentricRotor {
    DampedEccentricRotor::new(
        SpatialInertia::new(MASS, GaVec3::new(E, 0.0, 0.0), Mat3 { m: tensor }).unwrap(),
        Se3::identity(), bearings, || false,
    ).unwrap()
}
fn rotor(bearings: &[RadialBearing]) -> DampedEccentricRotor {
    full([1.2, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, POLAR], bearings)
}
fn motor() -> RotorWindingMotor {
    RotorWindingMotor {
        drive: RotorCurrentDrive {
            shaft_coordinate: 0,
            torque_constant_nm_per_a: 1.0,
            winding_resistance_ohm: 0.4,
            maximum_current_a: 10.0,
            maximum_terminal_voltage_v: 20.0,
            minimum_speed_rad_s: -100.0,
            maximum_speed_rad_s: 100.0,
            motoring_efficiency: 0.8,
            regenerative_efficiency: 0.5,
        },
        inductance_h: INDUCTANCE,
    }
}
fn control(h: f64) -> EccentricRotorVoltageControl {
    EccentricRotorVoltageControl {
        duration_seconds: h,
        relative_tolerance: 1e-11,
        electrical: RotorVoltageControl {
            source: RotorDrivePowerControl {
                maximum_source_input_w: 100.0,
                maximum_source_return_w: 100.0,
                absolute_power_tolerance_w: 1e-9,
                absolute_voltage_tolerance_v: 1e-10,
            },
            absolute_current_tolerance_a: 1e-12,
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
        angular_impulse_tolerance_nms: 1e-12,
        shaft_rate_tolerance_rad_s: 1e-11,
    }
}
fn command(voltage: f64, external: f64) -> EccentricRotorVoltageCommand {
    EccentricRotorVoltageCommand { voltage_v: voltage, external_shaft_torque_nm: external }
}
fn integrator() -> RigidBodyIntegrator { RigidBodyIntegrator::new(Gravity::ZERO) }
fn state(rotor: &DampedEccentricRotor, speed: f64, current_a: f64) -> EccentricRotorWindingState {
    let original = rotor.elastic_component();
    EccentricRotorWindingState {
        body: RigidBodyState::new(
            original.reference_pose(), Vec3::ZERO,
            original.mass_properties().angular_momentum_body_checked(Vec3::new(0.0, 0.0, speed)).unwrap(),
        ).unwrap(),
        current_a,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn speed(rotor: &DampedEccentricRotor, body: RigidBodyState) -> f64 {
    rotor.elastic_component().mass_properties()
        .angular_velocity_body_checked(body.angular_momentum_body()).unwrap().z
}
// Independent two-equation elimination, not the production Schur implementation.
fn source_oracle(w: f64, i: f64, voltage: f64, h: f64, inertia: f64) -> (f64, f64, f64) {
    let midpoint = (2.0 * INDUCTANCE * i + h * (voltage - w))
        / (2.0 * INDUCTANCE + 0.4 * h + h * h / (2.0 * inertia));
    (w + h * midpoint / inertia, 2.0 * midpoint - i, midpoint)
}

#[test]
fn g1_short_circuit_stores_current_without_false_back_emf_voltage_refusal() {
    let rotor = rotor(&[]);
    let mut model = motor();
    model.drive.maximum_terminal_voltage_v = 0.0;
    let drive = rotor.voltage_spin_motor(model, RotorSpinDriveAxis::Shaft, || false).unwrap();
    let initial = state(&rotor, 3.0, 0.0);
    let mut limits = control(0.1);
    limits.electrical.source.maximum_source_input_w = 0.0;
    limits.electrical.source.maximum_source_return_w = 0.0;
    let step = drive.step(integrator(), initial, command(0.0, 0.0), limits, || false).unwrap();
    let first = source_oracle(3.0, 0.0, 0.0, 0.05, POLAR);
    let second = source_oracle(first.0, first.1, 0.0, 0.05, POLAR);
    near(step.sources[0].electrical.winding_current_after_a, first.1, 3e-14);
    near(step.state_after.current_a, second.1, 3e-14);
    near(speed(&rotor, step.state_after.body), second.0, 3e-14);
    near(step.magnetic_energy_change_j, 0.5 * INDUCTANCE * second.1 * second.1, 3e-14);
    near(step.winding_heat_j, 0.05 * 0.4 * (first.2 * first.2 + second.2 * second.2), 3e-14);
    assert!(step.state_after.current_a < 0.0 && step.magnetic_energy_change_j > 0.0);
    assert!(step.winding_heat_j > 0.0 && step.kinetic_energy_change_j < 0.0);
    assert_eq!(step.endpoint_terminal_voltage_v, [0.0; 2]);
    assert_eq!(step.electrical_draw_j + step.electrical_return_j + step.bridge_heat_j, 0.0);
    assert!(speed(&rotor, step.state_after.body) + 0.4 * step.state_after.current_a > 2.0);
    near(step.kinetic_energy_change_j + step.magnetic_energy_change_j + step.winding_heat_j, 0.0, 3e-13);
    let next = drive.step(integrator(), step.state_after, command(0.0, 0.0), limits, || false).unwrap();
    near(next.sources[0].electrical.magnetic_energy_before_j,
        0.5 * INDUCTANCE * step.state_after.current_a.powi(2), 3e-14);
    assert!(next.state_after.current_a < step.state_after.current_a);
}

#[test]
fn g1_nonprincipal_inductive_source_uses_original_mobility_and_tensor_work() {
    let rotor = full([2.0, 0.0, 0.5, 0.0, 2.0, 0.0, 0.5, 0.0, 1.0], &[]);
    let drive = rotor.voltage_spin_motor(motor(), RotorSpinDriveAxis::Shaft, || false).unwrap();
    near(drive.source_inertia_kg_m2(), 7.0 / 8.0, 3e-14);
    let initial = EccentricRotorWindingState {
        body: RigidBodyState::new(rotor.elastic_component().reference_pose(),
            Vec3::new(30.0, -60.0, 9.0), Vec3::new(0.4, -0.2, 0.7)).unwrap(),
        current_a: 0.2,
    };
    let step = drive.step(integrator(), initial, command(1.3, 0.0), control(0.04), || false).unwrap();
    let source = &step.sources[0];
    let oracle = source_oracle(1.2 / 1.75, 0.2, 1.3, 0.02, 7.0 / 8.0);
    let impulse = 0.02 * oracle.2;
    let after = source.impulse.state_after;
    assert_eq!(after.pose(), initial.body.pose());
    assert_eq!(after.linear_momentum_world(), initial.body.linear_momentum_world());
    near(after.angular_momentum_body().x, 0.4, 0.0);
    near(after.angular_momentum_body().y, -0.2, 0.0);
    near(after.angular_momentum_body().z, 0.7 + impulse, 3e-14);
    near(source.shaft_rate_after_rad_s, oracle.0, 3e-14);
    near(source.electrical.winding_current_after_a, oracle.1, 3e-14);
    let expected_work = (1.2 / 1.75) * impulse + 0.5 * (8.0 / 7.0) * impulse * impulse;
    near(source.impulse.work.kinetic_energy_change_j, expected_work, 3e-14);
    near(source.electrical.shaft_work_j, expected_work, 3e-14);
    near(source.winding_voltage_residual_v, 0.0, 3e-13);
    near(step.sources[1].electrical.magnetic_energy_before_j, 0.5 * INDUCTANCE * oracle.1.powi(2), 3e-14);
    near(step.sources[1].electrical.midpoint_current_a,
        0.5 * (oracle.1 + step.state_after.current_a), 3e-14);
    near(step.kinetic_energy_change_j, initial.body.kinetic_energy_change_to(
        rotor.elastic_component().mass_properties(), step.state_after.body).unwrap(), 3e-14);
}

// Continuous planar ODE: X,Y,phi,Vx,Vy,omega,I and five independent physical accounts.
fn derivative(y: [f64; 12]) -> [f64; 12] {
    let (s, c) = y[2].sin_cos();
    let q = [y[0] - E * c, y[1] - E * s];
    let u = [y[3] + E * y[5] * s, y[4] - E * y[5] * c];
    let force = [-30.0 * q[0] - 5.0 * q[1] - 0.8 * u[0] - 0.1 * u[1],
        -5.0 * q[0] - 18.0 * q[1] - 0.1 * u[0] - 0.5 * u[1]];
    let terminal = 2.2 * y[6];
    let source = if terminal >= 0.0 { terminal / 0.8 } else { terminal * 0.5 };
    [y[3], y[4], y[5], force[0] / MASS, force[1] / MASS,
        (y[6] + 0.05 - E * (-s * force[0] + c * force[1])) / POLAR,
        (2.2 - 0.4 * y[6] - y[5]) / INDUCTANCE,
        u[0] * (0.8 * u[0] + 0.1 * u[1]) + u[1] * (0.1 * u[0] + 0.5 * u[1]),
        source.max(0.0), (-source).max(0.0), 0.4 * y[6] * y[6], source - terminal]
}
fn reference() -> [f64; 12] {
    let mut y = [0.0; 12]; y[0] = E;
    let h = 2.0 / 8000.0;
    for _ in 0..8000 {
        let a = derivative(y);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = derivative(core::array::from_fn(|i| y[i] + h * c[i]));
        y = core::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}
fn trajectory(n: usize, exact: [f64; 12]) -> (f64, f64, f64) {
    let rotor = rotor(&[bearing()]);
    let drive = rotor.voltage_spin_motor(motor(), RotorSpinDriveAxis::Shaft, || false).unwrap();
    let mut current = state(&rotor, 0.0, 0.0);
    let mut totals = [0.0; 5];
    let mut defect = 0.0;
    for _ in 0..n {
        let step = drive.step(integrator(), current, command(2.2, 0.05), control(2.0 / n as f64), || false).unwrap();
        for source in &step.sources { near(source.electrical.applied_voltage_v, 2.2, 1e-14); }
        let values = [step.bearing_heat_j, step.electrical_draw_j, step.electrical_return_j,
            step.winding_heat_j, step.bridge_heat_j];
        for i in 0..5 { assert!(values[i] >= 0.0); totals[i] += values[i]; }
        defect += step.total_energy_residual_j;
        current = step.state_after;
    }
    let p = current.body.pose().position_world();
    let v = current.body.center_of_mass_velocity_world(rotor.elastic_component().mass_properties()).unwrap();
    let axis = current.body.pose().orientation().rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    let error = (p.x - exact[0]).abs() + (p.y - exact[1]).abs()
        + (v.x - exact[3]).abs() + (v.y - exact[4]).abs()
        + (speed(&rotor, current.body) - exact[5]).abs() + (current.current_a - exact[6]).abs()
        + (axis.x - exact[2].cos()).abs() + (axis.y - exact[2].sin()).abs();
    (error, (0..5).map(|i| (totals[i] - exact[i + 7]).abs()).sum(), defect.abs())
}

#[test]
fn g1_g3_startup_has_finite_current_rise_and_coupled_orbit_energy_refinement() {
    let rotor = rotor(&[bearing()]);
    let drive = rotor.voltage_spin_motor(motor(), RotorSpinDriveAxis::Shaft, || false).unwrap();
    let start = drive.step(integrator(), state(&rotor, 0.0, 0.0),
        command(2.2, 0.05), control(1e-6), || false).unwrap();
    near(start.state_after.current_a / 1e-6, 2.2 / INDUCTANCE, 3e-5);
    near(speed(&rotor, start.state_after.body) / 1e-6, 0.05 / POLAR, 3e-5);
    let exact = reference();
    let a = trajectory(50, exact);
    let b = trajectory(100, exact);
    let c = trajectory(200, exact);
    assert!(a.0 / b.0 > 3.5 && b.0 / c.0 > 3.5, "{a:?} {b:?} {c:?}");
    assert!(a.1 / b.1 > 3.5 && b.1 / c.1 > 3.5, "{a:?} {b:?} {c:?}");
    assert!(c.0 < 5e-4 && c.1 < 5e-4 && c.2 < 1e-3, "{c:?}");
}

#[test]
fn g1_voltage_power_curtailment_and_inductive_regeneration_preserve_source_storage() {
    let rotor = rotor(&[]);
    let drive = rotor.voltage_spin_motor(motor(), RotorSpinDriveAxis::Shaft, || false).unwrap();
    let mut limited = control(0.1);
    limited.electrical.source.maximum_source_input_w = 1.0;
    let start = drive.step(integrator(), state(&rotor, 0.0, 0.0),
        command(20.0, 0.0), limited, || false).unwrap();
    assert!(start.sources.iter().all(|s| s.electrical.voltage_fraction < 1.0));
    assert!(start.state_after.current_a > 0.0 && start.magnetic_energy_change_j > 0.0);
    for s in &start.sources {
        assert!(s.electrical.initial_source_power_w <= 1.0 + 1e-9);
        assert!(s.electrical.final_source_power_w <= 1.0 + 1e-9);
    }
    let mut model = motor(); model.drive.maximum_terminal_voltage_v = 1.0;
    let drive = rotor.voltage_spin_motor(model, RotorSpinDriveAxis::Shaft, || false).unwrap();
    let initial = state(&rotor, 3.0, -0.5);
    let regenerated = drive.step(integrator(), initial, command(1.0, 0.0), control(0.1), || false).unwrap();
    assert!(regenerated.electrical_return_j > 0.0 && regenerated.bridge_heat_j > 0.0);
    assert_eq!(regenerated.electrical_draw_j, 0.0);
    assert!(regenerated.winding_heat_j > 0.0);
    assert!(speed(&rotor, regenerated.state_after.body) + 0.4 * regenerated.state_after.current_a > 1.0);
    let balance = regenerated.kinetic_energy_change_j + regenerated.magnetic_energy_change_j
        + regenerated.winding_heat_j + regenerated.bridge_heat_j + regenerated.electrical_return_j;
    near(balance, 0.0, 3e-13);
    let bad = EccentricRotorWindingState { current_a: 10.1, ..initial };
    assert!(matches!(drive.step(integrator(), bad, command(1.0, 0.0), control(0.1), || false),
        Err(EccentricRotorVoltageError::Electrical(RotorVoltageError::CurrentCompliance { .. }))));
}

#[test]
fn g4_outer_damping_speed_and_final_cancellation_cannot_publish_partial_current() {
    let rotor = rotor(&[RadialBearing {
        stiffness_n_per_m: [[0.0; 2]; 2], damping_ns_per_m: [[0.0, 0.0], [0.0, 4.0]], ..bearing()
    }]);
    let mut initial = state(&rotor, 0.1, 0.0);
    initial.body = RigidBodyState::new(initial.body.pose(), Vec3::new(0.0, 2.0 * MASS, 0.0),
        initial.body.angular_momentum_body()).unwrap();
    let mut model = motor(); model.drive.maximum_terminal_voltage_v = 0.0;
    let drive = rotor.voltage_spin_motor(model, RotorSpinDriveAxis::Shaft, || false).unwrap();
    let polls = Cell::new(0);
    let expected = drive.step(integrator(), initial, command(0.0, 0.0), control(0.1), || {
        polls.set(polls.get() + 1); false
    }).unwrap();
    let middle = expected.sources[1].shaft_rate_after_rad_s;
    let final_speed = speed(&rotor, expected.state_after.body);
    assert!(final_speed > middle + 0.05);
    let cap = 0.5 * (middle + final_speed);
    assert!(expected.sources.iter().all(|s| s.shaft_rate_after_rad_s < cap));
    let mut constrained = model; constrained.drive.maximum_speed_rad_s = cap;
    assert!(matches!(rotor.voltage_spin_motor(constrained, RotorSpinDriveAxis::Shaft, || false).unwrap()
        .step(integrator(), initial, command(0.0, 0.0), control(0.1), || false),
        Err(EccentricRotorVoltageError::Electrical(RotorVoltageError::Drive(RotorDriveError::SpeedDomain { .. })))));
    let last = polls.get(); polls.set(0);
    assert_eq!(drive.step(integrator(), initial, command(0.0, 0.0), control(0.1), || {
        polls.set(polls.get() + 1); polls.get() == last
    }), Err(EccentricRotorVoltageError::Cancelled));
    assert_eq!(drive.step(integrator(), initial, command(0.0, 0.0), control(0.1), || false).unwrap(), expected);
}

#[test]
fn g4_source_cutback_cannot_discard_existing_magnetic_state() {
    let rotor = rotor(&[]);
    let mut model = motor();
    model.inductance_h = 0.01;
    model.drive.maximum_current_a = 1.0;
    let drive = rotor.voltage_spin_motor(model, RotorSpinDriveAxis::Shaft, || false).unwrap();
    let initial = state(&rotor, 3.0, 0.5);
    let expected = drive.step(integrator(), initial, command(4.0, 0.0), control(0.1), || false).unwrap();
    let mut disconnected = control(0.1);
    disconnected.electrical.source.maximum_source_input_w = 0.0;
    disconnected.electrical.source.maximum_source_return_w = 0.0;
    // With the original nonzero current, the needed positive voltage draws
    // unavailable power. Cutting it to zero would force current below -1 A;
    // the native source must refuse, never clip current or erase magnetic energy.
    assert!(matches!(drive.step(integrator(), initial, command(4.0, 0.0), disconnected, || false),
        Err(EccentricRotorVoltageError::Electrical(RotorVoltageError::CurrentCompliance { .. }))));
    assert_eq!(drive.step(integrator(), initial, command(4.0, 0.0), control(0.1), || false).unwrap(), expected);
}
