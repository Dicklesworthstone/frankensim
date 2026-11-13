//! Finite mechanical power acts through the original full rotor dynamics.
use fs_mbd::rotordynamics::transient::{
    RotorPeakPowerControl, RotorPowerControl, RotorPowerError, RotorTransientState,
};
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::cell::Cell;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}

fn limits(motoring: f64, regeneration: f64) -> RotorPowerControl {
    RotorPowerControl {
        maximum_motoring_power_w: motoring,
        maximum_regenerative_power_w: regeneration,
        absolute_power_tolerance_w: 1e-11,
    }
}

fn state(velocity: f64) -> RotorTransientState {
    RotorTransientState {
        displacement: vec![0.0],
        velocity: vec![velocity],
    }
}

fn peak_limits(motoring: f64, regeneration: f64) -> RotorPeakPowerControl {
    RotorPeakPowerControl {
        maximum_motoring_power_w: motoring,
        maximum_regenerative_power_w: regeneration,
        absolute_power_tolerance_w: 1e-11,
    }
}

#[test]
fn peak_power_limits_stall_acceleration_instead_of_only_its_average() {
    let model = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.25, 1e-12, || false)
        .unwrap();
    let input = state(0.0);
    let peak = stepper
        .step_peak_power_limited(&input, &[8.0], peak_limits(2.0, 0.0), || false)
        .unwrap();
    // Torque=4, I=2, dt=1/4: endpoint power=2 W, mean power=1 W.
    near(peak.average.effort_fraction, 0.5, 1e-15);
    near(peak.initial_input_power_w, 0.0, 0.0);
    near(peak.final_input_power_w, 2.0, 1e-14);
    near(peak.average.average_input_power_w, 1.0, 1e-14);
    near(peak.average.motion.kinetic_energy_change_j, 0.25, 1e-14);
    assert_eq!(
        peak.average.motion,
        stepper.step(&input, &[4.0], || false).unwrap()
    );
    let full = stepper
        .step_peak_power_limited(&input, &[8.0], peak_limits(8.0, 0.0), || false)
        .unwrap();
    assert_eq!(full.average.effort_fraction, 1.0);
    assert_eq!(
        full.average.motion,
        stepper.step(&input, &[8.0], || false).unwrap()
    );
}

#[test]
fn zero_net_reversal_cannot_hide_positive_and_negative_power_peaks() {
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let input = state(-1.0);
    let average = stepper
        .step_power_limited(&input, &[2.0], limits(0.0, 0.0), || false)
        .unwrap();
    assert_eq!(average.effort_fraction, 1.0);
    assert_eq!(average.motion.input_work_j, 0.0);
    let peak = stepper
        .step_peak_power_limited(&input, &[2.0], peak_limits(0.0, 0.0), || false)
        .unwrap();
    assert_eq!(peak.average.effort_fraction, 0.0);
    // No effort: the rotor coasts through the step at its own speed.
    assert_eq!(
        peak.average.motion,
        stepper.step(&input, &[0.0], || false).unwrap()
    );
    assert_eq!(
        peak.average.motion.state_after,
        RotorTransientState {
            displacement: vec![-1.0],
            velocity: vec![-1.0],
        }
    );

    // Starting regeneration limits alpha before final motoring is considered.
    let peak = stepper
        .step_peak_power_limited(&input, &[4.0], peak_limits(0.25, 2.0), || false)
        .unwrap();
    let fraction = (1.0 + 2.0_f64.sqrt()) / 8.0;
    near(peak.average.effort_fraction, fraction, 2e-15);
    near(peak.initial_input_power_w, -4.0 * fraction, 1e-14);
    near(peak.final_input_power_w, 0.25, 1e-14);
    // Every affine-time power lies between the independently computed ends.
    for sample in 0..=20 {
        let t = f64::from(sample) / 20.0;
        let velocity = -1.0 + t * (peak.average.motion.state_after.velocity[0] + 1.0);
        let power = 4.0 * fraction * velocity;
        assert!((-2.0 - 1e-12..=0.25 + 1e-12).contains(&power));
    }
    let limited_regeneration = stepper
        .step_peak_power_limited(&input, &[4.0], peak_limits(100.0, 0.5), || false)
        .unwrap();
    near(limited_regeneration.average.effort_fraction, 0.125, 1e-15);
    near(limited_regeneration.initial_input_power_w, -0.5, 1e-15);
    near(limited_regeneration.final_input_power_w, -0.25, 1e-15);
}

#[test]
fn peak_power_retains_largest_disconnected_active_bearing_fraction() {
    let model = RotorBearingModel::new(1, &[1.0], &[-6.0], &[4.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let input = RotorTransientState {
        displacement: vec![0.75],
        velocity: vec![0.0],
    };
    // The full original equation gives v1=3-4*alpha. Therefore endpoint power
    // is 12*alpha-16*alpha². Its largest allowed root is beyond the forbidden
    // motoring interval, so a monotone search from zero would be wrong.
    let peak = stepper
        .step_peak_power_limited(&input, &[4.0], peak_limits(0.01, 0.25), || false)
        .unwrap();
    near(
        peak.average.effort_fraction,
        (12.0 + 160.0_f64.sqrt()) / 32.0,
        2e-15,
    );
    near(peak.initial_input_power_w, 0.0, 0.0);
    near(peak.final_input_power_w, -0.25, 3e-14);
    near(peak.average.motion.work_energy_residual_j, 0.0, 2e-13);
    let intermediate = stepper.step(&input, &[2.0], || false).unwrap();
    assert!(2.0 * intermediate.state_after.velocity[0] > 0.01);
}

#[test]
fn peak_power_is_scale_invariant_and_preserves_tiny_stall_effort() {
    for scale in [1e-180, 1.0, 1e180] {
        let model =
            RotorBearingModel::new(1, &[2.0 * scale], &[0.0], &[0.0], &[0.0], || false).unwrap();
        let stepper = model
            .at_speed(0.0, || false)
            .unwrap()
            .midpoint_stepper(0.25, 1e-12, || false)
            .unwrap();
        let peak = stepper
            .step_peak_power_limited(
                &state(0.0),
                &[8.0 * scale],
                RotorPeakPowerControl {
                    maximum_motoring_power_w: 2.0 * scale,
                    maximum_regenerative_power_w: 1e300,
                    absolute_power_tolerance_w: scale * 1e-12,
                },
                || false,
            )
            .unwrap();
        near(peak.average.effort_fraction, 0.5, 3e-15);
        near(peak.final_input_power_w / scale, 2.0, 3e-14);
        near(peak.average.motion.state_after.velocity[0], 0.5, 3e-14);
    }
    let model = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.25, 1e-12, || false)
        .unwrap();
    let peak = stepper
        .step_peak_power_limited(
            &state(0.0),
            &[8e150],
            RotorPeakPowerControl {
                maximum_motoring_power_w: 8e-300,
                maximum_regenerative_power_w: 0.0,
                absolute_power_tolerance_w: 1e-311,
            },
            || false,
        )
        .unwrap();
    near(peak.average.effort_fraction / 1e-300, 1.0, 2e-14);
    near(peak.final_input_power_w / 8e-300, 1.0, 2e-14);
}

#[test]
fn peak_power_refusal_and_cancellation_publish_no_motion() {
    let model = RotorBearingModel::new(1, &[1.3], &[0.07], &[2.1], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.3, 1e-12, || false)
        .unwrap();
    let input = state(0.37);
    let saved = input.clone();
    let control = peak_limits(0.43, 0.2);
    let calls = Cell::new(0);
    let accepted = stepper
        .step_peak_power_limited(&input, &[9.0], control, || {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    for stop in 1..=calls.get() {
        let count = Cell::new(0);
        assert_eq!(
            stepper.step_peak_power_limited(&input, &[9.0], control, || {
                count.set(count.get() + 1);
                count.get() == stop
            }),
            Err(RotorPowerError::Dynamics(RotorError::Cancelled))
        );
        assert_eq!(input, saved);
    }
    assert_eq!(
        stepper
            .step_peak_power_limited(&input, &[9.0], control, || false)
            .unwrap(),
        accepted
    );
    for invalid in [f64::NAN, f64::INFINITY, -1.0] {
        for bad in [
            RotorPeakPowerControl {
                maximum_motoring_power_w: invalid,
                ..control
            },
            RotorPeakPowerControl {
                maximum_regenerative_power_w: invalid,
                ..control
            },
            RotorPeakPowerControl {
                absolute_power_tolerance_w: invalid,
                ..control
            },
        ] {
            assert!(matches!(
                stepper.step_peak_power_limited(&input, &[9.0], bad, || false),
                Err(RotorPowerError::Invalid(_))
            ));
        }
    }
}

#[test]
fn finite_stall_torque_debits_actual_acceleration_work_and_preserves_full_command() {
    let model = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.25, 1e-12, || false)
        .unwrap();
    let input = state(0.0);
    let limited = stepper
        .step_power_limited(&input, &[8.0], limits(1.0, 0.0), || false)
        .unwrap();
    // I=2, dt=1/4, torque=4: omega=1/2, kinetic work=1/4 J.
    near(limited.effort_fraction, 0.5, 1e-15);
    near(limited.motion.state_after.velocity[0], 0.5, 1e-15);
    near(limited.motion.state_after.displacement[0], 0.0625, 1e-15);
    near(limited.average_input_power_w, 1.0, 1e-14);
    near(limited.motion.kinetic_energy_change_j, 0.25, 1e-14);
    near(limited.motion.input_work_j, 0.25, 1e-14);
    near(limited.motion.work_energy_residual_j, 0.0, 1e-14);

    let full = stepper
        .step_power_limited(&input, &[8.0], limits(4.0, 0.0), || false)
        .unwrap();
    assert_eq!(full.effort_fraction, 1.0);
    assert_eq!(full.motion, stepper.step(&input, &[8.0], || false).unwrap());
    let off = stepper
        .step_power_limited(&input, &[8.0], limits(0.0, 1e300), || false)
        .unwrap();
    assert_eq!(off.effort_fraction, 0.0);
    assert_eq!(off.motion.state_after, input);
    assert_eq!(off.average_input_power_w, 0.0);
    assert_eq!(
        stepper
            .step_power_limited(&input, &[8.0], limits(1.0, 0.0), || false)
            .unwrap(),
        limited
    );
}

#[test]
fn constant_power_acceleration_and_regeneration_match_energy_and_refine_position() {
    let model = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    for (initial_speed, torque, power) in [(1.0, 20.0, 3.0_f64), (2.0, -20.0, -1.0)] {
        let mut errors = Vec::new();
        for count in [8, 16, 32] {
            let h = 1.0 / f64::from(count);
            let stepper = model
                .at_speed(0.0, || false)
                .unwrap()
                .midpoint_stepper(h, 1e-12, || false)
                .unwrap();
            let mut current = state(initial_speed);
            let control = limits(power.max(0.0), (-power).max(0.0));
            let mut work = 0.0;
            for _ in 0..count {
                let step = stepper
                    .step_power_limited(&current, &[torque], control, || false)
                    .unwrap();
                assert!((0.0..1.0).contains(&step.effort_fraction));
                near(step.average_input_power_w, power, 1e-12);
                near(step.motion.work_energy_residual_j, 0.0, 1e-13);
                work += step.motion.input_work_j;
                current = step.motion.state_after;
            }
            // d(omega²)/dt=P for I=2; integrate that exact smooth velocity.
            let final_speed = (initial_speed * initial_speed + power).sqrt();
            let angle = 2.0 * (final_speed.powi(3) - initial_speed.powi(3)) / (3.0 * power);
            near(current.velocity[0], final_speed, 3e-13);
            near(work, power, 3e-12);
            near(
                current.velocity[0].powi(2) - initial_speed.powi(2),
                work,
                3e-12,
            );
            errors.push((current.displacement[0] - angle).abs());
        }
        for pair in errors.windows(2) {
            assert!((3.8..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
    }
}

#[test]
fn largest_fraction_survives_nonmonotone_work_reversal_and_active_bearings() {
    for damping in [0.0, -4.0] {
        let model =
            RotorBearingModel::new(1, &[1.0], &[damping], &[0.0], &[0.0], || false).unwrap();
        let stepper = model
            .at_speed(0.0, || false)
            .unwrap()
            .midpoint_stepper(1.0, 1e-12, || false)
            .unwrap();
        let input = state(-1.0);
        let (control, sign) = if damping == 0.0 {
            (limits(0.25, 0.01), 1.0)
        } else {
            (limits(0.01, 0.25), -1.0)
        };
        // P(alpha)= +/- (8 alpha² - 4 alpha). The lower-alpha interval is
        // disconnected from the largest admissible fraction by the other cap.
        let step = stepper
            .step_power_limited(&input, &[4.0], control, || false)
            .unwrap();
        near(step.effort_fraction, (4.0 + 24.0_f64.sqrt()) / 16.0, 2e-15);
        near(step.average_input_power_w, sign * 0.25, 2e-14);
        let larger = stepper
            .step(&input, &[4.0 * (step.effort_fraction + 1e-6)], || false)
            .unwrap();
        assert!(sign * larger.input_work_j > 0.25);
        let early = stepper.step(&input, &[1.0], || false).unwrap();
        assert!(sign * early.input_work_j < -0.01);

        let zero = stepper
            .step_power_limited(&input, &[4.0], limits(0.0, 0.0), || false)
            .unwrap();
        near(zero.effort_fraction, 0.5, 2e-15);
        near(zero.motion.state_after.velocity[0], 1.0, 2e-14);
        near(zero.average_input_power_w, 0.0, 2e-14);
    }
}

fn solve2(matrix: [f64; 4], rhs: [f64; 2]) -> [f64; 2] {
    let determinant = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    [
        (matrix[3] * rhs[0] - matrix[1] * rhs[1]) / determinant,
        (matrix[0] * rhs[1] - matrix[2] * rhs[0]) / determinant,
    ]
}

#[test]
fn full_inertia_gyro_and_cross_coupled_bearings_keep_original_power_and_units() {
    let mass = [2.0, 0.5, 0.5, 1.0];
    let damping = [0.2, 0.1, -0.4, 0.5];
    let stiffness = [3.0, 0.7, -0.1, 2.0];
    let gyro = [0.0, 0.3, -0.3, 0.0];
    let h = 0.05;
    let spin = 4.0;
    let q = [0.4, -0.2];
    let v = [0.3, -0.7];
    let force = [20.0, -5.0];
    let velocity_operator: [f64; 4] = core::array::from_fn(|i| damping[i] + spin * gyro[i]);
    let matrix: [f64; 4] = core::array::from_fn(|i| {
        mass[i] + 0.5 * h * velocity_operator[i] + 0.25 * h * h * stiffness[i]
    });
    let free_rhs: [f64; 2] = core::array::from_fn(|i| {
        -h * (0..2)
            .map(|j| {
                velocity_operator[2 * i + j] * v[j] + stiffness[2 * i + j] * (q[j] + 0.5 * h * v[j])
            })
            .sum::<f64>()
    });
    let free = solve2(matrix, free_rhs);
    let influence = solve2(matrix, force.map(|f| h * f));
    let a: f64 = (0..2).map(|i| 0.5 * force[i] * influence[i]).sum();
    let b: f64 = (0..2).map(|i| force[i] * (v[i] + 0.5 * free[i])).sum();
    let expected_fraction = 2.0 / (b + (b * b + 4.0 * a).sqrt());
    let end_a: f64 = (0..2).map(|i| force[i] * influence[i]).sum();
    let end_b: f64 = (0..2).map(|i| force[i] * (v[i] + free[i])).sum();
    let peak_fraction = 2.0 / (end_b + (end_b * end_b + 4.0 * end_a).sqrt());

    for scales in [[1.0, 1.0], [1e-6, 1e6]] {
        let map_matrix = |matrix: [f64; 4]| -> [f64; 4] {
            core::array::from_fn(|i| matrix[i] / scales[i / 2] / scales[i % 2])
        };
        let model = RotorBearingModel::new(
            2,
            &map_matrix(mass),
            &map_matrix(damping),
            &map_matrix(stiffness),
            &map_matrix(gyro),
            || false,
        )
        .unwrap();
        let stepper = model
            .at_speed(spin, || false)
            .unwrap()
            .midpoint_stepper(h, 1e-12, || false)
            .unwrap();
        let input = RotorTransientState {
            displacement: (0..2).map(|i| scales[i] * q[i]).collect(),
            velocity: (0..2).map(|i| scales[i] * v[i]).collect(),
        };
        let command: Vec<_> = (0..2).map(|i| force[i] / scales[i]).collect();
        let result = stepper
            .step_power_limited(&input, &command, limits(1.0, 0.75), || false)
            .unwrap();
        near(result.effort_fraction, expected_fraction, 2e-14);
        near(result.average_input_power_w, 1.0, 2e-13);
        near(result.motion.gyroscopic_work_j, 0.0, 1e-14);
        near(result.motion.work_energy_residual_j, 0.0, 1e-13);
        let peak = stepper
            .step_peak_power_limited(&input, &command, peak_limits(1.0, 0.75), || false)
            .unwrap();
        near(peak.average.effort_fraction, peak_fraction, 2e-14);
        near(peak.final_input_power_w, 1.0, 2e-13);
        near(
            peak.initial_input_power_w,
            peak_fraction * (force[0] * v[0] + force[1] * v[1]),
            2e-13,
        );
        near(peak.average.motion.work_energy_residual_j, 0.0, 1e-13);
        for i in 0..2 {
            let delta = free[i] + expected_fraction * influence[i];
            near(
                result.motion.state_after.velocity[i] / scales[i],
                v[i] + delta,
                2e-14,
            );
            near(
                result.motion.state_after.displacement[i] / scales[i],
                q[i] + h * (v[i] + 0.5 * delta),
                2e-14,
            );
            near(
                peak.average.motion.state_after.velocity[i] / scales[i],
                v[i] + free[i] + peak_fraction * influence[i],
                2e-14,
            );
        }
    }
}

#[test]
fn power_scaling_and_an_enormous_unused_reverse_allowance_do_not_erase_the_cap() {
    for scale in [1e-180, 1.0, 1e180] {
        let model =
            RotorBearingModel::new(1, &[2.0 * scale], &[0.0], &[0.0], &[0.0], || false).unwrap();
        let stepper = model
            .at_speed(0.0, || false)
            .unwrap()
            .midpoint_stepper(0.25, 1e-12, || false)
            .unwrap();
        let control = RotorPowerControl {
            maximum_motoring_power_w: scale,
            maximum_regenerative_power_w: 1e300,
            absolute_power_tolerance_w: scale * 1e-12,
        };
        let result = stepper
            .step_power_limited(&state(0.0), &[8.0 * scale], control, || false)
            .unwrap();
        near(result.effort_fraction, 0.5, 3e-15);
        near(result.average_input_power_w / scale, 1.0, 3e-14);
        near(result.motion.state_after.velocity[0], 0.5, 3e-14);
    }

    // The power ratio 1e-600 underflows, but its square root 1e-300 and
    // every actual force, motion and work quantity remain representable.
    let model = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(1.0, 1e-12, || false)
        .unwrap();
    let result = stepper
        .step_power_limited(
            &state(0.0),
            &[2e150],
            RotorPowerControl {
                maximum_motoring_power_w: 1e-300,
                maximum_regenerative_power_w: 0.0,
                absolute_power_tolerance_w: 1e-312,
            },
            || false,
        )
        .unwrap();
    near(result.effort_fraction / 1e-300, 1.0, 5e-14);
    near(result.motion.state_after.velocity[0] / 1e-150, 1.0, 5e-14);
    near(result.average_input_power_w / 1e-300, 1.0, 5e-14);
}

#[test]
fn cancellation_and_invalid_inputs_leave_original_state_and_command_unchanged() {
    let model = RotorBearingModel::new(1, &[2.0], &[0.3], &[1.7], &[0.0], || false).unwrap();
    let stepper = model
        .at_speed(0.0, || false)
        .unwrap()
        .midpoint_stepper(0.1, 1e-12, || false)
        .unwrap();
    let input = state(0.7);
    let saved = input.clone();
    let command = [9.0];
    let control = limits(0.5, 0.2);
    let calls = Cell::new(0);
    let accepted = stepper
        .step_power_limited(&input, &command, control, || {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    assert_ne!(accepted.motion.work_energy_residual_j, 0.0);
    assert!(matches!(
        stepper.step_power_limited(
            &input,
            &command,
            RotorPowerControl {
                absolute_power_tolerance_w: 0.0,
                ..control
            },
            || false,
        ),
        Err(RotorPowerError::PowerBalance { .. })
    ));
    for cancel_at in 1..=calls.get() {
        let seen = Cell::new(0);
        let result = stepper.step_power_limited(&input, &command, control, || {
            seen.set(seen.get() + 1);
            seen.get() == cancel_at
        });
        assert_eq!(
            result,
            Err(RotorPowerError::Dynamics(RotorError::Cancelled))
        );
        assert_eq!(input, saved);
    }
    assert_eq!(
        stepper
            .step_power_limited(&input, &command, control, || false)
            .unwrap(),
        accepted
    );
    for invalid in [f64::NAN, f64::INFINITY, -1.0] {
        for bad_control in [
            RotorPowerControl {
                maximum_motoring_power_w: invalid,
                ..control
            },
            RotorPowerControl {
                maximum_regenerative_power_w: invalid,
                ..control
            },
            RotorPowerControl {
                absolute_power_tolerance_w: invalid,
                ..control
            },
        ] {
            assert!(matches!(
                stepper.step_power_limited(&input, &command, bad_control, || false),
                Err(RotorPowerError::Invalid(_))
            ));
        }
    }
    assert!(
        stepper
            .step_power_limited(&input, &[], control, || false)
            .is_err()
    );
    assert!(
        stepper
            .step_power_limited(&input, &[f64::NAN], control, || false)
            .is_err()
    );
    assert_eq!(input, saved);
    assert_eq!(command, [9.0]);
}
