//! Independent decay/work oracles for the prepared stiff viscous subflow.

use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn diagonal(x: f64, y: f64, z: f64) -> Mat3 {
    Mat3 {
        m: [x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, z],
    }
}

fn rotation(value: Vec3) -> So3 {
    So3::exp(So3Tangent::new(value)).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn stiff_axial_relaxation_keeps_the_exponential_decay_and_physical_heat() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    let damping = diagonal(0.0, 0.0, 400.0);
    assert!(inertia.viscous_stepper(damping, 0.1).is_err());
    let stepper = inertia.viscous_exponential_stepper(damping, 0.1).unwrap();
    let result = stepper
        .step(
            So3::identity(),
            Vec3::new(0.0, 0.0, 6.0),
            HeldTorque::Body(Vec3::default()),
        )
        .unwrap();
    let expected = 6.0 * (-20.0_f64).exp();
    assert!(result.momentum_after_body.z > 0.0);
    close(result.momentum_after_body.z, expected, 2e-19);
    close(
        result.dissipated_energy_j,
        9.0 * (1.0 - (-40.0_f64).exp()),
        5e-14,
    );
    close(result.work_energy_residual_j, 0.0, 5e-14);
    assert!(result.angular_momentum_residual_space.norm() < 5e-14);

    // Complete stiff decay is an accepted limit rather than a spurious
    // midpoint sign flip. Attitude accuracy still requires step refinement.
    let saturated = inertia
        .viscous_exponential_stepper(diagonal(0.0, 0.0, 4e5), 0.1)
        .unwrap()
        .step(
            So3::identity(),
            Vec3::new(0.0, 0.0, 6.0),
            HeldTorque::Body(Vec3::default()),
        )
        .unwrap();
    close(saturated.momentum_after_body.z, 0.0, 5e-14);
    close(saturated.dissipated_energy_j, 9.0, 5e-14);
}

#[test]
fn stiff_driven_rotating_bath_integrates_source_work_and_slip_heat() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    let damping: f64 = 80.0;
    let drive: f64 = 1.4;
    let bath: f64 = 2.5;
    let initial: f64 = -0.7;
    let duration: f64 = 0.2;
    let lambda = damping / 2.0;
    let terminal = bath + drive / damping;
    let amplitude = initial - terminal;
    let decay = (-lambda * duration).exp();
    let expected = terminal + amplitude * decay;
    let angle = terminal * duration + amplitude * (1.0 - decay) / lambda;
    let slip_terminal = drive / damping;
    let heat = damping
        * (slip_terminal * slip_terminal * duration
            + 2.0 * slip_terminal * amplitude * (1.0 - decay) / lambda
            + amplitude * amplitude * (1.0 - decay * decay) / (2.0 * lambda));
    let result = inertia
        .viscous_exponential_stepper(diagonal(0.0, 0.0, damping), duration)
        .unwrap()
        .step_in_rotating_bath(
            So3::identity(),
            Vec3::new(0.0, 0.0, 2.0 * initial),
            HeldTorque::Body(Vec3::new(0.0, 0.0, drive)),
            Vec3::new(0.0, 0.0, bath),
        )
        .unwrap();
    close(result.motion.angular_velocity_after_body.z, expected, 2e-14);
    close(result.motion.drive_work_j, drive * angle, 2e-14);
    close(
        result.bath_work_j,
        bath * (2.0 * (expected - initial) - drive * duration),
        8e-14,
    );
    close(result.motion.dissipated_energy_j, heat, 8e-14);
    close(result.motion.work_energy_residual_j, 0.0, 8e-14);
}

#[test]
fn all_cross_axis_decay_modes_match_an_independent_closed_solution() {
    // Spherical inertia eliminates the gyroscopic momentum term, so the
    // original body momentum has a closed solution even when D is not diagonal.
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    let axes = rotation(Vec3::new(0.6, -0.4, 0.2));
    let damping = axes
        .matrix()
        .compose(diagonal(0.4, 20.0, 400.0))
        .compose(axes.matrix().transpose());
    let modes = Vec3::new(1.1, -0.8, 2.3);
    let initial = axes.rotate(modes).unwrap();
    let expected_modes = Vec3::new(
        modes.x * (-0.02_f64).exp(),
        modes.y * (-1.0_f64).exp(),
        modes.z * (-20.0_f64).exp(),
    );
    let expected = axes.rotate(expected_modes).unwrap();
    let result = inertia
        .viscous_exponential_stepper(damping, 0.1)
        .unwrap()
        .step(So3::identity(), initial, HeldTorque::Body(Vec3::default()))
        .unwrap();
    assert!((result.momentum_after_body - expected).norm() < 5e-13);
    let heat = (modes.dot(modes) - expected_modes.dot(expected_modes)) / 4.0;
    close(result.dissipated_energy_j, heat, 3e-13);
    close(result.work_energy_residual_j, 0.0, 3e-13);
}

#[test]
fn full_inertia_and_bath_are_covariant_and_phase_refines() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), So3::identity()).unwrap();
    let damping = Mat3 {
        m: [0.6, 0.1, -0.05, 0.1, 0.9, 0.2, -0.05, 0.2, 1.2],
    };
    let axes = rotation(Vec3::new(0.4, -0.3, 0.2));
    let rotated_inertia = InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), axes).unwrap();
    let rotated_damping = axes
        .matrix()
        .compose(damping)
        .compose(axes.matrix().transpose());
    let a = inertia.viscous_exponential_stepper(damping, 0.025).unwrap();
    let b = rotated_inertia
        .viscous_exponential_stepper(rotated_damping, 0.025)
        .unwrap();
    let initial_rotation = rotation(Vec3::new(-0.3, 0.2, 0.1));
    let initial_momentum = Vec3::new(1.2, -0.9, 2.8);
    let torque = HeldTorque::Space(Vec3::new(0.2, -0.1, 0.3));
    let bath = Vec3::new(-0.5, 0.7, 0.1);
    let original = a
        .step_in_rotating_bath(initial_rotation, initial_momentum, torque, bath)
        .unwrap();
    let changed = b
        .step_in_rotating_bath(
            initial_rotation.compose(axes.inverse()).unwrap(),
            axes.rotate(initial_momentum).unwrap(),
            torque,
            bath,
        )
        .unwrap();
    assert!(
        (changed.motion.momentum_after_body
            - axes.rotate(original.motion.momentum_after_body).unwrap())
        .norm()
            < 2e-12
    );
    close(
        changed.motion.dissipated_energy_j,
        original.motion.dissipated_energy_j,
        2e-13,
    );
    close(changed.bath_work_j, original.bath_work_j, 2e-13);
    close(changed.motion.work_energy_residual_j, 0.0, 2e-13);

    let spherical =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    let expected_angle = 3.0 * (1.0 - (-2.0_f64).exp()) / 2.0;
    let expected_rotation = rotation(Vec3::new(0.0, 0.0, expected_angle));
    let mut errors = Vec::new();
    for count in [10, 20, 40] {
        let method = spherical
            .viscous_exponential_stepper(diagonal(0.0, 0.0, 4.0), 1.0 / count as f64)
            .unwrap();
        let mut attitude = So3::identity();
        let mut momentum = Vec3::new(0.0, 0.0, 6.0);
        for _ in 0..count {
            let step = method
                .step(attitude, momentum, HeldTorque::Body(Vec3::default()))
                .unwrap();
            attitude = step.rotation_after;
            momentum = step.momentum_after_body;
        }
        close(momentum.z, 6.0 * (-2.0_f64).exp(), 2e-14);
        errors.push(
            attitude
                .body_minus(expected_rotation)
                .unwrap()
                .angular
                .norm(),
        );
    }
    assert!(errors[1] < 0.26 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.26 * errors[1], "{errors:?}");
}

#[test]
fn weak_drag_zero_drag_and_tensor_admission_keep_their_limits() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    let state = Vec3::new(0.0, 0.0, 6.0);
    let torque = HeldTorque::Body(Vec3::default());
    let weak = inertia
        .viscous_exponential_stepper(diagonal(0.0, 0.0, 2e-12), 0.1)
        .unwrap()
        .step(So3::identity(), state, torque)
        .unwrap();
    close(
        weak.dissipated_energy_j,
        9.0 * -(-2e-13_f64).exp_m1(),
        2e-26,
    );
    let zero = inertia
        .viscous_exponential_stepper(Mat3::zero(), 0.1)
        .unwrap()
        .step(So3::identity(), state, torque)
        .unwrap();
    let native = inertia
        .viscous_stepper(Mat3::zero(), 0.1)
        .unwrap()
        .step(So3::identity(), state, torque)
        .unwrap();
    assert_eq!(zero, native);
    assert!(
        inertia
            .viscous_exponential_stepper(diagonal(-1.0, 1.0, 1.0), 0.1)
            .is_err()
    );
    assert!(
        inertia
            .viscous_exponential_stepper(diagonal(1.0, 1.0, 1.0), 0.0)
            .is_err()
    );
}

#[test]
fn a_large_undamped_null_direction_does_not_hide_transverse_heat() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    let damping = Mat3 {
        m: [2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 4e20],
    };
    let initial = Vec3::new(1e20, -1e20, 2.0);
    let result = inertia
        .viscous_exponential_stepper(damping, 1e-20)
        .unwrap()
        .step(So3::identity(), initial, HeldTorque::Body(Vec3::default()))
        .unwrap();
    assert_eq!(result.momentum_after_body.x, initial.x);
    assert_eq!(result.momentum_after_body.y, initial.y);
    close(result.momentum_after_body.z, 2.0 * (-2.0_f64).exp(), 2e-14);
    close(result.dissipated_energy_j, 1.0 - (-4.0_f64).exp(), 2e-14);
    close(
        result.kinetic_energy_change_j,
        -result.dissipated_energy_j,
        2e-14,
    );
    close(result.work_energy_residual_j, 0.0, 2e-14);
}
