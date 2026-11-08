//! Finite-source power in the original full-coordinate rotor equations.
use fs_mbd::rotordynamics::RotorBearingModel;
use fs_mbd::rotordynamics::transient::RotorTransientState;
use fs_mbd::rotordynamics::transient::drive::voltage::battery::{
    BatteryRotorError, BatteryRotorState,
};
use fs_mbd::rotordynamics::transient::drive::voltage::{
    RotorVoltageControl, RotorWindingMotor, RotorWindingState,
};
use fs_mbd::rotordynamics::transient::drive::{RotorCurrentDrive, RotorDrivePowerControl};
use fs_phs::SolveBudget;
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::converter::battery_fed::BatteryFedBuckBoost;
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};

const BATTERY: BatteryParameters = BatteryParameters {
    series_resistance_ohm: 0.1,
    polarization_resistance_ohm: 0.2,
    polarization_capacitance_f: 0.5,
    bus_capacitance_f: 0.03,
};
const CONVERTER: BuckBoostParameters = BuckBoostParameters {
    inductance_h: 0.04,
    bus_capacitance_f: [0.04, 0.05],
    winding_resistance_ohm: 0.07,
    input_on_resistance_ohm: [0.02, 0.03],
    output_on_resistance_ohm: [0.02, 0.03],
};
const LEAD: f64 = 0.08;
const CONVERTER_R: f64 = 0.07 + 0.5 * 0.02 + 0.5 * 0.03 + 0.02;

fn source() -> BatteryFedBuckBoost {
    let battery = BatterySource::new(
        BatteryOcvCurve::new(80.0, &[(0.0, 12.0), (1.0, 12.0)]).unwrap(),
        BATTERY,
    )
    .unwrap();
    AveragedBuckBoost::new(
        CONVERTER,
        BridgeModulation {
            input_high_fraction: 0.5,
            output_high_fraction: 1.0,
        },
    )
    .unwrap()
    .with_battery_source(battery, LEAD)
    .unwrap()
}

fn motor() -> RotorWindingMotor {
    RotorWindingMotor {
        drive: RotorCurrentDrive {
            shaft_coordinate: 0,
            torque_constant_nm_per_a: 0.4,
            winding_resistance_ohm: 0.2,
            maximum_current_a: 100.0,
            maximum_terminal_voltage_v: 100.0,
            minimum_speed_rad_s: -100.0,
            maximum_speed_rad_s: 100.0,
            motoring_efficiency: 1.0,
            regenerative_efficiency: 1.0,
        },
        inductance_h: 0.02,
    }
}
fn electrical_controls() -> PowerLoadControls {
    PowerLoadControls {
        solve: SolveBudget {
            relative: 1e-13,
            absolute: 1e-13,
            max_iterations: 40,
        },
        backtracking_halvings: 8,
        absolute_energy_tolerance_j: 1e-9,
        relative_energy_tolerance: 1e-10,
    }
}
fn motor_controls() -> RotorVoltageControl {
    RotorVoltageControl {
        source: RotorDrivePowerControl {
            maximum_source_input_w: 10000.0,
            maximum_source_return_w: 10000.0,
            absolute_power_tolerance_w: 1e-7,
            absolute_voltage_tolerance_v: 1e-9,
        },
        absolute_current_tolerance_a: 1e-10,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn action(matrix: &[f64; 4], vector: &[f64; 2]) -> [f64; 2] {
    [
        matrix[0] * vector[0] + matrix[1] * vector[1],
        matrix[2] * vector[0] + matrix[3] * vector[1],
    ]
}
fn dot(a: &[f64; 2], b: &[f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn kinetic(mass: &[f64; 4], velocity: &[f64]) -> f64 {
    let v = [velocity[0], velocity[1]];
    0.5 * dot(&v, &action(mass, &v))
}

#[test]
fn finite_source_obeys_original_multicoordinate_motion_kirchhoff_and_work_equations() {
    let mass = [0.08, 0.01, 0.01, 0.12];
    // Active, cross-coupled damping and nonsymmetric stiffness must retain
    // their original signed work; neither is silently converted to heat.
    let damping = [-0.02, 0.003, -0.004, -0.01];
    let stiffness = [0.4, 0.2, -0.1, 0.3];
    let gyro = [0.0, 0.05, -0.05, 0.0];
    let model = RotorBearingModel::new(2, &mass, &damping, &stiffness, &gyro, || false).unwrap();
    let dt = 0.001;
    let native = model
        .at_speed(2.0, || false)
        .unwrap()
        .midpoint_stepper(dt, 1e-11, || false)
        .unwrap();
    let electrical = source();
    let mut state = BatteryRotorState {
        electrical: electrical.state(0.65, 0.1, 11.8, 0.5, [11.7, 5.8]).unwrap(),
        motor: RotorWindingState {
            rotor: RotorTransientState {
                displacement: vec![0.1, -0.2],
                velocity: vec![4.0, -1.5],
            },
            current_a: 0.1,
        },
    };
    let original = state.clone();
    let mut drive = native
        .voltage_motor(motor(), || false)
        .unwrap()
        .with_battery_converter(electrical, electrical_controls(), motor_controls(), || {
            false
        })
        .unwrap();
    let initial_energy = drive.electrical().system().hamiltonian(&state.electrical)
        + 0.5 * motor().inductance_h * state.motor.current_a.powi(2)
        + kinetic(&mass, &state.motor.rotor.velocity);
    let external = [-0.2, 0.03];
    let auxiliary = RegulatedPowerLoad::try_new(1.0, 10.0, 0.5, 100.0).unwrap();
    let mut accounted_work = 0.0;
    for iteration in 0..40 {
        let step = drive
            .step_controlled(&state, &external, Some(auxiliary), || false)
            .unwrap();
        let q: [f64; 2] = core::array::from_fn(|i| {
            0.5 * (state.motor.rotor.displacement[i]
                + step.motor.motion.state_after.displacement[i])
        });
        let v: [f64; 2] = core::array::from_fn(|i| {
            0.5 * (state.motor.rotor.velocity[i] + step.motor.motion.state_after.velocity[i])
        });
        let delta_v: [f64; 2] = core::array::from_fn(|i| {
            step.motor.motion.state_after.velocity[i] - state.motor.rotor.velocity[i]
        });
        let im = 0.5 * (state.motor.current_a + step.motor.winding_current_after_a);
        let midpoint: [f64; 6] =
            core::array::from_fn(|i| 0.5 * (state.electrical[i] + step.electrical_after[i]));
        let vp = midpoint[1] / BATTERY.polarization_capacitance_f;
        let vb = midpoint[2] / BATTERY.bus_capacitance_f;
        let vi = midpoint[4] / CONVERTER.bus_capacitance_f[0];
        let vo = midpoint[5] / CONVERTER.bus_capacitance_f[1];
        let ic = midpoint[3] / CONVERTER.inductance_h;
        let ib = (12.0 - vp - vb) / BATTERY.series_resistance_ohm;
        let lead = (vb - vi) / LEAD;
        // Independent SI KCL/KVL, including the actual mechanical winding.
        let derivative = [
            -ib,
            ib - vp / BATTERY.polarization_resistance_ohm,
            ib - lead,
            0.5 * vi - vo - CONVERTER_R * ic,
            lead - 0.5 * ic,
            ic - im - 1.0 / vo,
        ];
        for (i, derivative) in derivative.into_iter().enumerate() {
            near(
                step.electrical_after[i] - state.electrical[i],
                dt * derivative,
                2e-11,
            );
        }
        near(
            motor().inductance_h * (step.motor.winding_current_after_a - state.motor.current_a),
            dt * (vo - 0.2 * im - 0.4 * v[0]),
            2e-12,
        );
        let inertia = action(&mass, &delta_v);
        let c = action(&damping, &v);
        let k = action(&stiffness, &q);
        let g = action(&gyro, &v);
        for i in 0..2 {
            near(
                inertia[i] + dt * (c[i] + k[i] + 2.0 * g[i]),
                dt * (external[i] + if i == 0 { 0.4 * im } else { 0.0 }),
                2e-12,
            );
            near(
                step.motor.motion.state_after.displacement[i] - state.motor.rotor.displacement[i],
                dt * v[i],
                1e-14,
            );
        }
        near(
            step.motor.motion.velocity_force_work_j,
            dt * dot(&v, &c),
            2e-13,
        );
        near(
            step.motor.motion.stiffness_force_work_j,
            dt * dot(&v, &k),
            2e-13,
        );
        near(
            step.motor.motion.gyroscopic_work_j,
            2.0 * dt * dot(&v, &g),
            2e-13,
        );
        assert!(step.motor.motion.velocity_force_work_j < 0.0);
        assert!(step.heat_joules().into_iter().all(|heat| heat >= 0.0));
        near(step.electrical.external_work_j, -dt * vo * im, 2e-12);
        near(step.balance_residual_j, 0.0, 2e-10);
        assert_eq!(
            step.motor.applied_voltage_v,
            step.electrical.discrete_bus_voltage_v[2]
        );
        if iteration == 0 {
            assert!((vo - 5.8).abs() > 1e-6);
        }
        accounted_work += dt * dot(&v, &c)
            + dt * dot(&v, &k)
            + 2.0 * dt * dot(&v, &g)
            + step.heat_joules().iter().sum::<f64>()
            + dt
            - dt * dot(&v, &external);
        state = step.into_state();
    }
    let final_energy = drive.electrical().system().hamiltonian(&state.electrical)
        + 0.5 * motor().inductance_h * state.motor.current_a.powi(2)
        + kinetic(&mass, &state.motor.rotor.velocity);
    near(final_energy - initial_energy + accounted_work, 0.0, 2e-8);
    assert!(state.electrical[0] < original.electrical[0]);
    assert_ne!(
        state.motor.rotor.velocity[1],
        original.motor.rotor.velocity[1]
    );
}

#[test]
fn mechanical_generator_charges_chemical_storage_through_the_same_real_losses() {
    let model = RotorBearingModel::new(1, &[0.1], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let dt = 0.01;
    let native = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(dt, 1e-12, || false)
        .unwrap();
    let electrical = source();
    let vb = 12.0 + BATTERY.series_resistance_ohm + BATTERY.polarization_resistance_ohm;
    let vi = vb + LEAD;
    let vo = 0.5 * vi + 2.0 * CONVERTER_R;
    let speed = (vo + 0.4) / 0.4;
    let mut state = BatteryRotorState {
        electrical: electrical.state(0.5, -0.2, vb, -2.0, [vi, vo]).unwrap(),
        motor: RotorWindingState {
            rotor: RotorTransientState {
                displacement: vec![0.0],
                velocity: vec![speed],
            },
            current_a: -2.0,
        },
    };
    let initial = state.clone();
    let mut drive = native
        .voltage_motor(motor(), || false)
        .unwrap()
        .with_battery_converter(electrical, electrical_controls(), motor_controls(), || {
            false
        })
        .unwrap();
    for _ in 0..20 {
        let step = drive
            .step_controlled(&state, &[0.8], None, || false)
            .unwrap();
        near(step.electrical_after[0] - state.electrical[0], dt, 2e-12);
        near(
            step.electrical.battery.chemical_energy_change_j,
            12.0 * dt,
            2e-11,
        );
        near(step.motor.motion.state_after.velocity[0], speed, 2e-11);
        near(step.motor.winding_current_after_a, -2.0, 2e-11);
        near(step.electrical.lead_current_a, -1.0, 2e-11);
        near(step.motor.electrical_return_j, 2.0 * vo * dt, 2e-12);
        near(step.motor.winding_heat_j, 0.8 * dt, 2e-12);
        near(
            step.heat_joules().iter().sum::<f64>(),
            dt * (4.0 * CONVERTER_R + 0.1 + 0.2 + LEAD + 0.8),
            2e-12,
        );
        near(step.terminal_work_residual_j, 0.0, 2e-12);
        near(step.balance_residual_j, 0.0, 2e-10);
        assert_eq!(step.motor.bridge_heat_j, 0.0);
        state = step.into_state();
    }
    near(
        state.electrical[0] - initial.electrical[0],
        20.0 * dt,
        2e-11,
    );
}

#[test]
fn voltage_curtailment_exhaustion_and_cancellation_publish_neither_state() {
    let model = RotorBearingModel::new(1, &[0.1], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let native = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.001, 1e-12, || false)
        .unwrap();
    let input = BatteryRotorState {
        electrical: source().state(0.5, 0.0, 11.8, 0.5, [11.7, 5.8]).unwrap(),
        motor: RotorWindingState {
            rotor: RotorTransientState {
                displacement: vec![0.0],
                velocity: vec![0.0],
            },
            current_a: 0.0,
        },
    };
    let snapshot = input.clone();
    let mut limited_motor = motor();
    limited_motor.drive.maximum_terminal_voltage_v = 1.0;
    let mut limited = native
        .voltage_motor(limited_motor, || false)
        .unwrap()
        .with_battery_converter(source(), electrical_controls(), motor_controls(), || false)
        .unwrap();
    assert!(matches!(
        limited.step_controlled(&input, &[0.0], None, || false),
        Err(BatteryRotorError::VoltageCurtailment { .. })
    ));
    assert_eq!(input, snapshot);
    let mut drive = native
        .voltage_motor(motor(), || false)
        .unwrap()
        .with_battery_converter(source(), electrical_controls(), motor_controls(), || false)
        .unwrap();
    let mut empty = input.clone();
    empty.electrical[0] = 0.0;
    let empty_snapshot = empty.clone();
    assert!(matches!(
        drive.step_controlled(&empty, &[0.0], None, || false),
        Err(BatteryRotorError::Electrical(_))
    ));
    assert_eq!(empty, empty_snapshot);
    assert!(
        drive
            .step_controlled(&input, &[0.0], None, || true)
            .is_err()
    );
    assert_eq!(input, snapshot);
    let retry = drive
        .step_controlled(&input, &[0.0], None, || false)
        .unwrap();
    let mut fresh = native
        .voltage_motor(motor(), || false)
        .unwrap()
        .with_battery_converter(source(), electrical_controls(), motor_controls(), || false)
        .unwrap();
    let expected = fresh
        .step_controlled(&input, &[0.0], None, || false)
        .unwrap();
    assert_eq!(retry.electrical_after, expected.electrical_after);
    assert_eq!(retry.motor, expected.motor);
    let mut extra_bridge = motor();
    extra_bridge.drive.motoring_efficiency = 0.9;
    assert!(matches!(
        native
            .voltage_motor(extra_bridge, || false)
            .unwrap()
            .with_battery_converter(source(), electrical_controls(), motor_controls(), || false),
        Err(BatteryRotorError::Invalid(_))
    ));
}
