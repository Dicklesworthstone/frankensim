//! G0/G1/G3 full-inertia viscous torque: decay, work and material-frame covariance.

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
        "{actual:e} != {expected:e} (tol {tolerance:e})"
    );
}

#[test]
fn driven_viscous_spin_and_phase_refine_to_the_independent_exponential() {
    let inertia = InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity())
        .unwrap();
    let duration: f64 = 1.2;
    let initial_rate: f64 = -0.7;
    let drive: f64 = 1.4;
    let damping: f64 = 0.8;
    let lambda = damping / 2.0;
    let terminal = drive / damping;
    let amplitude = initial_rate - terminal;
    let exponential = (-lambda * duration).exp();
    let expected_rate = terminal + amplitude * exponential;
    let expected_angle =
        terminal * duration + amplitude * (1.0 - exponential) / lambda;
    let expected_heat = damping
        * (terminal * terminal * duration
            + 2.0 * terminal * amplitude * (1.0 - exponential) / lambda
            + amplitude * amplitude * (1.0 - exponential * exponential) / (2.0 * lambda));
    let initial_rotation = rotation(Vec3::new(0.3, -0.2, 0.4));
    let mut errors = Vec::new();
    for count in [12, 24, 48] {
        let stepper = inertia
            .viscous_stepper(diagonal(0.0, 0.0, damping), duration / count as f64)
            .unwrap();
        let mut attitude = initial_rotation;
        let mut momentum = Vec3::new(0.0, 0.0, 2.0 * initial_rate);
        let mut heat = 0.0;
        let mut work = 0.0;
        let mut drift = 0.0;
        for _ in 0..count {
            let step = stepper
                .step(
                    attitude,
                    momentum,
                    HeldTorque::Body(Vec3::new(0.0, 0.0, drive)),
                )
                .unwrap();
            assert!(step.dissipated_energy_j >= 0.0);
            close(step.work_energy_residual_j, 0.0, 5e-14);
            assert!(step.angular_momentum_residual_space.norm() < 5e-14);
            heat += step.dissipated_energy_j;
            work += step.drive_work_j;
            drift += step.drift_energy_defect_j;
            attitude = step.rotation_after;
            momentum = step.momentum_after_body;
        }
        let expected_rotation = initial_rotation
            .compose(rotation(Vec3::new(0.0, 0.0, expected_angle)))
            .unwrap();
        let rate = inertia.velocity(momentum).unwrap().z;
        let error = (rate - expected_rate).abs()
            + attitude.body_minus(expected_rotation).unwrap().angular.norm()
            + (heat - expected_heat).abs();
        errors.push(error);
        // Continuous work and phase approach the same independent integral.
        close(work, drive * expected_angle, 0.003 / (count * count) as f64 * 144.0);
        close(rate * rate - initial_rate * initial_rate, work - heat + drift, 2e-13);
    }
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
}

#[test]
fn passive_spin_decay_has_no_numerical_heat_from_gyroscopic_precession() {
    let inertia = InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity())
        .unwrap();
    let stepper = inertia.viscous_stepper(diagonal(0.0, 0.0, 4.0), 0.1).unwrap();
    let zero = HeldTorque::Body(Vec3::default());
    let mut attitude = So3::identity();
    let mut momentum = Vec3::new(0.0, 0.0, 6.0);
    let mut heat = 0.0;
    for _ in 0..10 {
        let step = stepper.step(attitude, momentum, zero).unwrap();
        assert!(step.momentum_after_body.z > 0.0);
        assert!(step.momentum_after_body.z < momentum.z);
        close(step.drive_work_j, 0.0, 0.0);
        close(step.drift_energy_defect_j, 0.0, 2e-14);
        heat += step.dissipated_energy_j;
        attitude = step.rotation_after;
        momentum = step.momentum_after_body;
    }
    let rate = momentum.z / 2.0;
    close(rate, 3.0 * (-2.0_f64).exp(), 0.001);
    close(heat, 9.0 - rate * rate, 2e-13);

    let asymmetric =
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), So3::identity()).unwrap();
    let free = asymmetric.viscous_stepper(Mat3::zero(), 0.1).unwrap();
    let step = free
        .step(So3::identity(), Vec3::new(1.2, -0.9, 2.8), zero)
        .unwrap();
    assert_eq!(step.dissipated_energy_j, 0.0);
    assert_eq!(step.drive_work_j, 0.0);
    close(step.kinetic_energy_change_j, step.drift_energy_defect_j, 1e-14);
}

#[test]
fn transverse_viscous_heat_survives_a_dominant_unchanged_spin() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    // The total axial energy has a nonzero DD product tail large enough to
    // absorb the transverse energy. A short interval still resolves rotation
    // and viscous relaxation: h omega_z = 0.5 and h D_x / I_x = 2.
    let initial = Vec3::new(2.0, 0.0, 1e20);
    let stepper = inertia
        .viscous_stepper(diagonal(4e20, 0.0, 0.0), 1e-20)
        .unwrap();
    let step = stepper
        .step(So3::identity(), initial, HeldTorque::Body(Vec3::default()))
        .unwrap();
    // Each implicit half-stage multiplies the transverse spin by 1/3.
    close(step.momentum_after_body.x, 2.0 / 9.0, 2e-15);
    assert_eq!(step.momentum_after_body.z, initial.z);
    close(step.kinetic_energy_change_j, -80.0 / 81.0, 2e-14);
    close(step.dissipated_energy_j, 80.0 / 81.0, 2e-14);
    assert_eq!(step.drift_energy_defect_j, 0.0);
    close(step.work_energy_residual_j, 0.0, 2e-14);
}

#[test]
fn full_inertia_and_cross_axis_drag_are_covariant_under_material_relabeling() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), So3::identity()).unwrap();
    let damping = Mat3 {
        m: [0.6, 0.1, -0.05, 0.1, 0.9, 0.2, -0.05, 0.2, 1.2],
    };
    let basis = rotation(Vec3::new(0.4, -0.3, 0.2));
    let changed_inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), basis).unwrap();
    let changed_damping = basis.matrix().compose(damping).compose(basis.matrix().transpose());
    let original = inertia.viscous_stepper(damping, 0.025).unwrap();
    let changed = changed_inertia.viscous_stepper(changed_damping, 0.025).unwrap();
    let mut attitude = rotation(Vec3::new(-0.3, 0.2, 0.1));
    let mut changed_attitude = attitude.compose(basis.inverse()).unwrap();
    let mut momentum = Vec3::new(1.2, -0.9, 2.8);
    let mut changed_momentum = basis.rotate(momentum).unwrap();
    let torque = HeldTorque::Space(Vec3::new(0.2, -0.1, 0.3));
    for _ in 0..16 {
        let a = original.step(attitude, momentum, torque).unwrap();
        let b = changed.step(changed_attitude, changed_momentum, torque).unwrap();
        assert!(
            (b.momentum_after_body - basis.rotate(a.momentum_after_body).unwrap()).norm()
                < 2e-12
        );
        assert!(
            b.rotation_after
                .body_minus(a.rotation_after.compose(basis.inverse()).unwrap())
                .unwrap()
                .angular
                .norm()
                < 2e-12
        );
        close(a.drive_work_j, b.drive_work_j, 2e-13);
        close(a.dissipated_energy_j, b.dissipated_energy_j, 2e-13);
        close(a.work_energy_residual_j, 0.0, 2e-13);
        close(b.work_energy_residual_j, 0.0, 2e-13);
        assert!(a.angular_momentum_residual_space.norm() < 2e-13);
        assert!(b.angular_momentum_residual_space.norm() < 2e-13);
        attitude = a.rotation_after;
        changed_attitude = b.rotation_after;
        momentum = a.momentum_after_body;
        changed_momentum = b.momentum_after_body;
    }
}

#[test]
fn zero_drag_retains_the_native_held_spatial_torque_trajectory() {
    let inertia = InertiaTensor3::from_principal(
        Vec3::new(2.0, 3.0, 4.0),
        rotation(Vec3::new(0.2, -0.3, 0.1)),
    )
    .unwrap();
    let attitude = rotation(Vec3::new(-0.1, 0.4, 0.2));
    let omega = Vec3::new(0.6, -0.3, 0.8);
    let momentum = inertia.momentum(omega).unwrap();
    let torque = HeldTorque::Space(Vec3::new(0.1, 0.2, -0.1));
    let duration = 0.03;
    let prepared = inertia.viscous_stepper(Mat3::zero(), duration).unwrap();
    let a = prepared.step(attitude, momentum, torque).unwrap();
    let (r, w) = inertia
        .torque_function_step(attitude, omega, 0.0, duration, |_, _| Ok(torque))
        .unwrap();
    assert!(a.rotation_after.body_minus(r).unwrap().angular.norm() < 2e-14);
    assert!((a.angular_velocity_after_body - w).norm() < 2e-14);
    assert_eq!(a.dissipated_energy_j, 0.0);
    if let HeldTorque::Space(value) = torque {
        assert!((a.drive_angular_impulse_space - value.scale(duration)).norm() < 2e-14);
    }
}

#[test]
fn invalid_drag_underresolved_relaxation_and_failed_steps_refuse_atomically() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(1.0, 1.0, 1.0), So3::identity()).unwrap();
    for bad in [
        diagonal(-0.1, 1.0, 1.0),
        diagonal(f64::NAN, 1.0, 1.0),
        Mat3 {
            m: [1.0, 0.5, 0.0, 0.1, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
        // A tolerance-stopped eigensolve can miss this indefinite zero row.
        Mat3 {
            m: [0.0, 1e-17, 0.0, 1e-17, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
        // All 2x2 principal minors are positive, but the determinant is negative.
        Mat3 {
            m: [1.0, -0.6, -0.6, -0.6, 1.0, -0.6, -0.6, -0.6, 1.0],
        },
    ] {
        assert!(inertia.viscous_stepper(bad, 0.01).is_err());
    }
    for duration in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(inertia.viscous_stepper(Mat3::zero(), duration).is_err());
    }
    assert!(inertia.viscous_stepper(diagonal(0.0, 0.0, 5.0), 1.0).is_err());
    assert!(inertia.viscous_stepper(diagonal(0.0, 0.0, 1.0), 0.1).is_ok());
    let prepared = inertia.viscous_stepper(diagonal(1.0, 1.0, 1.0), 0.1).unwrap();
    let rotation = So3::identity();
    let momentum = Vec3::new(0.2, -0.4, 0.7);
    assert!(
        prepared
            .step(rotation, momentum, HeldTorque::Body(Vec3::new(f64::NAN, 0.0, 0.0)))
            .is_err()
    );
    let first = prepared
        .step(rotation, momentum, HeldTorque::Body(Vec3::default()))
        .unwrap();
    let second = prepared
        .step(rotation, momentum, HeldTorque::Body(Vec3::default()))
        .unwrap();
    assert_eq!(first, second);
    assert_eq!(momentum, Vec3::new(0.2, -0.4, 0.7));
}
