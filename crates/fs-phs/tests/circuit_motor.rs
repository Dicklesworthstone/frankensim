//! Analytical electrical/rotational exchange and actual finite-bus integration.
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget, PreparedExactZoh};
use fs_phs::{InterconnectionBudget, PreparedStepError, interconnect_mapped};

fn parameters() -> DcMotorParameters {
    DcMotorParameters {
        inductance_h: 0.5,
        winding_resistance_ohm: 0.4,
        rotor_inertia_kg_m2: 2.0,
        stator_inertia_kg_m2: 6.0,
        coupling_nm_per_a: 1.5,
        bearing_drag_nm_s_per_rad: 0.07,
    }
}

fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 6,
        max_omega_dt: 0.04,
    }
}

fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 6,
        max_dense_bytes: 1 << 20,
        max_multiply_terms: 1_000_000,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 4e-11 * scale,
        "{actual:.17e} != {expected:.17e} (scale={scale:e})",
    );
}

#[test]
fn voltage_startup_matches_independent_oscillator_with_stator_recoil() {
    let dt = 0.4;
    for k in [-1.5, 1.5] {
        for voltage in [-6.0, 6.0] {
            let p = DcMotorParameters {
                winding_resistance_ohm: 0.0,
                bearing_drag_nm_s_per_rad: 0.0,
                coupling_nm_per_a: k,
                ..parameters()
            };
            let reduced_inertia = p.rotor_inertia_kg_m2 * p.stator_inertia_kg_m2
                / (p.rotor_inertia_kg_m2 + p.stator_inertia_kg_m2);
            let frequency = k.abs() / (p.inductance_h * reduced_inertia).sqrt();
            let (sin, cos) = (frequency * dt).sin_cos();
            let current = voltage * sin / (p.inductance_h * frequency);
            let relative_speed = voltage / k * (1.0 - cos);
            let mean_current =
                voltage * (1.0 - cos) / (p.inductance_h * frequency * frequency * dt);
            let relative_angle = voltage / k * (dt - sin / frequency);
            let rotor_fraction =
                p.stator_inertia_kg_m2 / (p.rotor_inertia_kg_m2 + p.stator_inertia_kg_m2);
            let stator_fraction = rotor_fraction - 1.0;
            let mut prepared = DcMotor::new(p)
                .unwrap()
                .prepare_controlled(dt, budget(), heat_budget(), || false)
                .unwrap();
            let (mut next, mut mean) = ([0.0; 3], [0.0; 3]);
            let step = prepared
                .step_controlled(
                    &[0.0; 3],
                    &[voltage, 0.0, 0.0],
                    &mut next,
                    &mut mean,
                    || false,
                )
                .unwrap();
            close(next[0] / p.inductance_h, current, 10.0);
            close(
                next[1] / p.rotor_inertia_kg_m2,
                rotor_fraction * relative_speed,
                10.0,
            );
            close(
                next[2] / p.stator_inertia_kg_m2,
                stator_fraction * relative_speed,
                10.0,
            );
            close(next[1] + next[2], 0.0, 10.0);
            close(mean[0], mean_current, 10.0);
            close(
                step.shaft_angle_increment_rad[0],
                rotor_fraction * relative_angle,
                1.0,
            );
            close(
                step.shaft_angle_increment_rad[1],
                stator_fraction * relative_angle,
                1.0,
            );
            close(step.mean_back_emf_v, k * relative_angle / dt, 10.0);
            close(step.mean_electromagnetic_torque_nm, k * mean_current, 10.0);
            close(step.electrical_work_j, voltage * dt * mean_current, 10.0);
            close(step.native.delta_h, step.electrical_work_j, 10.0);
            assert_eq!(step.shaft_work_j, [0.0; 2]);
            assert_eq!(step.copper_heat_j, 0.0);
            assert_eq!(step.bearing_heat_j, 0.0);
            assert!(next[1] * next[2] < 0.0);
        }
    }
}

#[test]
fn all_four_quadrants_keep_reciprocal_torque_work_and_separate_heat() {
    let p = parameters();
    let dt = 0.37;
    for current in [-2.0, 2.0] {
        for relative_speed in [-5.0, 5.0] {
            let speeds = [relative_speed + 1.3, 1.3];
            let voltage = p.winding_resistance_ohm * current + p.coupling_nm_per_a * relative_speed;
            let torque =
                p.coupling_nm_per_a * current - p.bearing_drag_nm_s_per_rad * relative_speed;
            let input = [voltage, -torque, torque];
            let model = DcMotor::new(p).unwrap();
            let initial = model.state(current, speeds).unwrap();
            let mut prepared = model
                .prepare_controlled(dt, budget(), heat_budget(), || false)
                .unwrap();
            let (mut next, mut mean) = ([0.0; 3], [0.0; 3]);
            let step = prepared
                .step_controlled(&initial, &input, &mut next, &mut mean, || false)
                .unwrap();
            for i in 0..3 {
                close(next[i], initial[i], 30.0);
            }
            close(mean[0], current, 10.0);
            close(mean[1], speeds[0], 10.0);
            close(mean[2], speeds[1], 10.0);
            close(step.electrical_work_j, dt * voltage * current, 10.0);
            for i in 0..2 {
                close(step.shaft_work_j[i], dt * input[i + 1] * speeds[i], 10.0);
            }
            close(
                step.mean_electromagnetic_torque_nm,
                p.coupling_nm_per_a * current,
                10.0,
            );
            close(
                step.mean_bearing_torque_nm,
                -p.bearing_drag_nm_s_per_rad * relative_speed,
                1.0,
            );
            close(
                step.copper_heat_j,
                dt * p.winding_resistance_ohm * current * current,
                1.0,
            );
            close(
                step.bearing_heat_j,
                dt * p.bearing_drag_nm_s_per_rad * relative_speed * relative_speed,
                1.0,
            );
            close(
                step.native.dissipated,
                step.copper_heat_j + step.bearing_heat_j,
                10.0,
            );
            close(
                step.electrical_work_j + step.shaft_work_j.iter().sum::<f64>(),
                step.copper_heat_j + step.bearing_heat_j,
                10.0,
            );
            assert_eq!(step.electrical_work_j < 0.0, current * relative_speed < 0.0);
        }
    }
}

#[test]
fn shorted_generator_matches_damped_oscillator_and_two_distinct_heat_integrals() {
    // L = J_reduced = 1 and R/L = b/J_reduced = gamma. Eliminating only in
    // this oracle gives i=-w0 exp(-gamma*t)sin(k*t), w=w0 exp(-gamma*t)cos(k*t).
    let (gamma, k, w0, dt) = (0.3, 2.0, 4.0, 0.6_f64);
    let p = DcMotorParameters {
        inductance_h: 1.0,
        winding_resistance_ohm: gamma,
        rotor_inertia_kg_m2: 2.0,
        stator_inertia_kg_m2: 2.0,
        coupling_nm_per_a: k,
        bearing_drag_nm_s_per_rad: gamma,
    };
    let model = DcMotor::new(p).unwrap();
    let initial = model.state(0.0, [3.0, -1.0]).unwrap();
    let mut prepared = model
        .prepare_controlled(dt, budget(), heat_budget(), || false)
        .unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 3]);
    let step = prepared
        .step_controlled(&initial, &[0.0; 3], &mut next, &mut mean, || false)
        .unwrap();
    let decay = (-gamma * dt).exp();
    let (sin, cos) = (k * dt).sin_cos();
    let current = -w0 * decay * sin;
    let relative_speed = w0 * decay * cos;
    close(next[0], current, 10.0);
    close(next[1] / 2.0, 1.0 + relative_speed / 2.0, 10.0);
    close(next[2] / 2.0, 1.0 - relative_speed / 2.0, 10.0);
    close(next[1] + next[2], initial[1] + initial[2], 10.0);
    let integral_sin = (k - decay * (gamma * sin + k * cos)) / (gamma * gamma + k * k);
    let integral_cos = (gamma + decay * (-gamma * cos + k * sin)) / (gamma * gamma + k * k);
    close(mean[0], -w0 * integral_sin / dt, 10.0);
    close(mean[1] - mean[2], w0 * integral_cos / dt, 10.0);
    let (a, frequency) = (2.0 * gamma, 2.0 * k);
    let integral_decay = (1.0 - (-a * dt).exp()) / a;
    let integral_cosine = (a
        + (-a * dt).exp() * (-a * (frequency * dt).cos() + frequency * (frequency * dt).sin()))
        / (a * a + frequency * frequency);
    let copper = gamma * w0 * w0 * (integral_decay - integral_cosine) / 2.0;
    let bearing = gamma * w0 * w0 * (integral_decay + integral_cosine) / 2.0;
    close(step.copper_heat_j, copper, 10.0);
    close(step.bearing_heat_j, bearing, 10.0);
    close(step.native.dissipated, copper + bearing, 10.0);
    close(step.native.delta_h, -copper - bearing, 10.0);
    assert!(next[0] < 0.0 && next[1] < initial[1]);
    assert!((copper - bearing).abs() > 0.2);
    assert_eq!(step.electrical_work_j, 0.0);
}

#[test]
fn decoupled_winding_and_external_torques_keep_actual_absolute_work() {
    let p = DcMotorParameters {
        coupling_nm_per_a: 0.0,
        winding_resistance_ohm: 0.0,
        bearing_drag_nm_s_per_rad: 0.0,
        ..parameters()
    };
    let (dt, input) = (0.3, [4.0, 3.0, -0.5]);
    let model = DcMotor::new(p).unwrap();
    let initial = model.state(2.0, [1.0, -2.0]).unwrap();
    let mut prepared = model
        .prepare_controlled(dt, budget(), heat_budget(), || false)
        .unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 3]);
    let step = prepared
        .step_controlled(&initial, &input, &mut next, &mut mean, || false)
        .unwrap();
    for (i, storage) in [
        p.inductance_h,
        p.rotor_inertia_kg_m2,
        p.stator_inertia_kg_m2,
    ]
    .iter()
    .enumerate()
    {
        close(next[i], initial[i] + dt * input[i], 30.0);
        close(mean[i], (initial[i] + 0.5 * dt * input[i]) / storage, 30.0);
    }
    close(
        next[1] + next[2] - initial[1] - initial[2],
        dt * (input[1] + input[2]),
        10.0,
    );
    close(
        step.native.delta_h,
        step.electrical_work_j + step.shaft_work_j.iter().sum::<f64>(),
        10.0,
    );
    assert_eq!(step.mean_electromagnetic_torque_nm, 0.0);
    assert_eq!(step.mean_back_emf_v, 0.0);
}

#[test]
fn finite_converter_buses_drive_and_accept_regeneration_in_one_native_solve() {
    let converter_parameters = BuckBoostParameters {
        inductance_h: 0.1,
        bus_capacitance_f: [0.2, 0.3],
        winding_resistance_ohm: 0.0,
        input_on_resistance_ohm: [0.0; 2],
        output_on_resistance_ohm: [0.0; 2],
    };
    let p = DcMotorParameters {
        inductance_h: 0.05,
        rotor_inertia_kg_m2: 0.2,
        stator_inertia_kg_m2: 0.3,
        coupling_nm_per_a: 0.8,
        winding_resistance_ohm: 0.0,
        bearing_drag_nm_s_per_rad: 0.0,
    };
    for generating in [false, true] {
        let converter = AveragedBuckBoost::new(
            converter_parameters,
            BridgeModulation {
                input_high_fraction: 0.75,
                output_high_fraction: 1.0,
            },
        )
        .unwrap();
        let motor = DcMotor::new(p).unwrap();
        let bus_initial = converter.state(0.0, [4.0, 3.0]).unwrap();
        let rotor_speed = if generating { 10.0 } else { 0.0 };
        let motor_initial = motor
            .state(
                0.0,
                [
                    rotor_speed,
                    -p.rotor_inertia_kg_m2 / p.stator_inertia_kg_m2 * rotor_speed,
                ],
            )
            .unwrap();
        let initial = [
            bus_initial[0],
            bus_initial[1],
            bus_initial[2],
            motor_initial[0],
            motor_initial[1],
            motor_initial[2],
        ];
        let composite = interconnect_mapped(
            converter.into_system(),
            motor.into_system(),
            &[1],
            &[0],
            &[1.0],
            InterconnectionBudget {
                max_states: 6,
                max_total_ports: 5,
                max_dense_bytes: 1 << 20,
                max_coupling_multiply_terms: 1_000_000,
            },
        )
        .unwrap();
        assert_eq!(composite.state_dim(), 6);
        assert_eq!(composite.port_dim(), 3);
        let mut prepared = PreparedExactZoh::new_owned(composite, 0.02, budget()).unwrap();
        let (mut next, mut mean) = ([0.0; 6], [0.0; 3]);
        // Surviving ports: source current, rotor torque, stator torque. Zero
        // source current leaves genuinely finite stored energy in both buses.
        let step = prepared
            .step_into(&initial, &[0.0; 3], &mut next, &mut mean)
            .unwrap();
        let mechanical = |x: &[f64; 6]| {
            0.5 * (x[4] * x[4] / p.rotor_inertia_kg_m2 + x[5] * x[5] / p.stator_inertia_kg_m2)
        };
        let electrical = |x: &[f64; 6]| {
            0.5 * (x[0] * x[0] / converter_parameters.inductance_h
                + x[1] * x[1] / converter_parameters.bus_capacitance_f[0]
                + x[2] * x[2] / converter_parameters.bus_capacitance_f[1]
                + x[3] * x[3] / p.inductance_h)
        };
        let mechanical_change = mechanical(&next) - mechanical(&initial);
        let electrical_change = electrical(&next) - electrical(&initial);
        close(mechanical_change + electrical_change, 0.0, 30.0);
        close(next[4] + next[5], initial[4] + initial[5], 10.0);
        close(step.delta_h, 0.0, 30.0);
        assert_eq!(step.supplied, 0.0);
        assert_eq!(step.dissipated, 0.0);
        assert_eq!(next[3] < 0.0, generating);
        assert_eq!(mechanical_change < 0.0, generating);
        assert_eq!(electrical_change > 0.0, generating);
        assert_eq!(next[1] > initial[1], generating);
        assert_eq!(next[2] > initial[2], generating);
    }
}

#[test]
fn physical_admission_and_late_cancellation_publish_neither_output() {
    let p = parameters();
    for invalid in [
        DcMotorParameters {
            inductance_h: 0.0,
            ..p
        },
        DcMotorParameters {
            rotor_inertia_kg_m2: -1.0,
            ..p
        },
        DcMotorParameters {
            stator_inertia_kg_m2: f64::INFINITY,
            ..p
        },
        DcMotorParameters {
            winding_resistance_ohm: -0.1,
            ..p
        },
        DcMotorParameters {
            bearing_drag_nm_s_per_rad: -0.1,
            ..p
        },
        DcMotorParameters {
            coupling_nm_per_a: f64::NAN,
            ..p
        },
    ] {
        assert!(DcMotor::new(invalid).is_err());
    }
    let model = DcMotor::new(p).unwrap();
    assert!(model.state(0.0, [f64::NAN, 0.0]).is_err());
    let initial = model.state(0.8, [2.0, -0.5]).unwrap();
    let input = [6.0, -1.0, 0.25];
    let mut prepared = model
        .prepare_controlled(0.1, budget(), heat_budget(), || false)
        .unwrap();
    let (mut next, mut mean) = ([0.0; 3], [0.0; 3]);
    let mut polls = 0;
    prepared
        .step_controlled(&initial, &input, &mut next, &mut mean, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = (next, mean);
    let sentinel = ([17.0, 18.0, 19.0], [20.0, 21.0, 22.0]);
    (next, mean) = sentinel;
    let mut poll_number = 0;
    assert_eq!(
        prepared
            .step_controlled(&initial, &input, &mut next, &mut mean, || {
                poll_number += 1;
                poll_number == polls
            })
            .unwrap_err(),
        PreparedStepError::Cancelled
    );
    assert_eq!((next, mean), sentinel);
    assert!(
        prepared
            .step_controlled(
                &initial,
                &[f64::NAN, 0.0, 0.0],
                &mut next,
                &mut mean,
                || false
            )
            .is_err()
    );
    assert_eq!((next, mean), sentinel);
    prepared
        .step_controlled(&initial, &input, &mut next, &mut mean, || false)
        .unwrap();
    assert_eq!((next, mean), expected);
}
