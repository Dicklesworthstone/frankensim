//! G1/G3 prescribed rotating-bath torque: entrainment, source work and slip heat.

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
        "{actual:e} != {expected:e}"
    );
}

fn full_drag() -> Mat3 {
    Mat3 {
        m: [0.6, 0.1, -0.05, 0.1, 0.9, 0.2, -0.05, 0.2, 1.2],
    }
}

#[test]
fn g1_rotating_bath_entrains_reversing_spin_with_independent_work_and_heat() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    let duration: f64 = 1.2;
    let initial_rate: f64 = -0.5;
    let bath: f64 = 1.1;
    let drive: f64 = 0.3;
    let damping: f64 = 0.8;
    let lambda = damping / 2.0;
    let relative_terminal = drive / damping;
    let terminal = bath + relative_terminal;
    let amplitude = initial_rate - terminal;
    let exponential = (-lambda * duration).exp();
    let expected_rate = terminal + amplitude * exponential;
    let expected_angle = terminal * duration + amplitude * (1.0 - exponential) / lambda;
    let expected_heat = damping
        * (relative_terminal * relative_terminal * duration
            + 2.0 * relative_terminal * amplitude * (1.0 - exponential) / lambda
            + amplitude * amplitude * (1.0 - exponential * exponential) / (2.0 * lambda));
    let expected_bath_work = bath * (2.0 * (expected_rate - initial_rate) - drive * duration);
    assert!(expected_rate > 0.0 && expected_bath_work > 0.0 && expected_heat > 0.0);
    let mut errors = Vec::new();
    for count in [12, 24, 48] {
        let stepper = inertia
            .viscous_stepper(diagonal(0.0, 0.0, damping), duration / count as f64)
            .unwrap();
        let mut attitude = So3::identity();
        let mut momentum = Vec3::new(0.0, 0.0, 2.0 * initial_rate);
        let (mut work, mut bath_work, mut heat, mut drift) = (0.0, 0.0, 0.0, 0.0);
        for _ in 0..count {
            let step = stepper
                .step_in_rotating_bath(
                    attitude,
                    momentum,
                    HeldTorque::Body(Vec3::new(0.0, 0.0, drive)),
                    Vec3::new(0.0, 0.0, bath),
                )
                .unwrap();
            let motion = step.motion;
            close(motion.work_energy_residual_j, 0.0, 5e-14);
            assert!(motion.angular_momentum_residual_space.norm() < 5e-14);
            assert!(motion.dissipated_energy_j >= 0.0);
            close(
                step.bath_work_j,
                bath * motion.drag_angular_impulse_space.z,
                2e-14,
            );
            work += motion.drive_work_j;
            bath_work += step.bath_work_j;
            heat += motion.dissipated_energy_j;
            drift += motion.drift_energy_defect_j;
            attitude = motion.rotation_after;
            momentum = motion.momentum_after_body;
        }
        let rate = momentum.z / 2.0;
        let expected_rotation = rotation(Vec3::new(0.0, 0.0, expected_angle));
        errors.push(
            (rate - expected_rate).abs()
                + attitude
                    .body_minus(expected_rotation)
                    .unwrap()
                    .angular
                    .norm()
                + (bath_work - expected_bath_work).abs()
                + (heat - expected_heat).abs(),
        );
        close(
            work,
            drive * expected_angle,
            0.001 / (count * count) as f64 * 144.0,
        );
        close(
            rate * rate - initial_rate * initial_rate,
            work + bath_work - heat + drift,
            2e-13,
        );
    }
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
}

type Vector = [f64; 3];
type Reference = [f64; 10]; // quaternion, omega, drive work, bath work, heat

fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: Vector, b: Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn rotate(q: [f64; 4], v: Vector) -> Vector {
    let axis = [q[1], q[2], q[3]];
    let cross = cross(axis, v);
    let parallel = dot(axis, v);
    let scalar = q[0] * q[0] - dot(axis, axis);
    core::array::from_fn(|i| scalar * v[i] + 2.0 * axis[i] * parallel + 2.0 * q[0] * cross[i])
}

fn derivative(state: Reference) -> Reference {
    let q = [state[0], state[1], state[2], state[3]];
    let w = [state[4], state[5], state[6]];
    let inverse = [q[0], -q[1], -q[2], -q[3]];
    let drive = rotate(inverse, [0.2, -0.1, 0.3]);
    let bath = rotate(inverse, [-0.3, 0.5, 1.1]);
    let relative: Vector = core::array::from_fn(|i| w[i] - bath[i]);
    let d = full_drag().m;
    let drag: Vector = core::array::from_fn(|i| {
        -(d[3 * i] * relative[0] + d[3 * i + 1] * relative[1] + d[3 * i + 2] * relative[2])
    });
    let gyro = cross(w, [2.0 * w[0], 3.0 * w[1], 4.0 * w[2]]);
    let wdot: Vector =
        core::array::from_fn(|i| (drive[i] + drag[i] - gyro[i]) / [2.0, 3.0, 4.0][i]);
    let imaginary = cross([q[1], q[2], q[3]], w);
    [
        -0.5 * (q[1] * w[0] + q[2] * w[1] + q[3] * w[2]),
        0.5 * (q[0] * w[0] + imaginary[0]),
        0.5 * (q[0] * w[1] + imaginary[1]),
        0.5 * (q[0] * w[2] + imaginary[2]),
        wdot[0],
        wdot[1],
        wdot[2],
        dot(drive, w),
        dot(drag, bath),
        -dot(drag, relative),
    ]
}

fn reference() -> Reference {
    let mut state = [1.0, 0.0, 0.0, 0.0, 0.6, -0.3, 0.7, 0.0, 0.0, 0.0];
    let h = 0.8 / 8192.0;
    for _ in 0..8192 {
        let k1 = derivative(state);
        let k2 = derivative(core::array::from_fn(|i| state[i] + 0.5 * h * k1[i]));
        let k3 = derivative(core::array::from_fn(|i| state[i] + 0.5 * h * k2[i]));
        let k4 = derivative(core::array::from_fn(|i| state[i] + h * k3[i]));
        state = core::array::from_fn(|i| {
            state[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    state
}

#[test]
fn g1_full_tensor_rotating_bath_refines_against_independent_euler_rk4() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), So3::identity()).unwrap();
    let reference = reference();
    let expected_momentum = Vec3::new(2.0 * reference[4], 3.0 * reference[5], 4.0 * reference[6]);
    let q = [reference[0], reference[1], reference[2], reference[3]];
    let marker = rotate(q, [1.0, 0.0, 0.0]);
    let marker = Vec3::new(marker[0], marker[1], marker[2]);
    let mut errors = Vec::new();
    for count in [12, 24, 48] {
        let stepper = inertia
            .viscous_stepper(full_drag(), 0.8 / count as f64)
            .unwrap();
        let mut attitude = So3::identity();
        let mut momentum = Vec3::new(1.2, -0.9, 2.8);
        let (mut work, mut bath_work, mut heat) = (0.0, 0.0, 0.0);
        for _ in 0..count {
            let step = stepper
                .step_in_rotating_bath(
                    attitude,
                    momentum,
                    HeldTorque::Space(Vec3::new(0.2, -0.1, 0.3)),
                    Vec3::new(-0.3, 0.5, 1.1),
                )
                .unwrap();
            close(step.motion.work_energy_residual_j, 0.0, 3e-13);
            assert!(step.motion.angular_momentum_residual_space.norm() < 3e-13);
            work += step.motion.drive_work_j;
            bath_work += step.bath_work_j;
            heat += step.motion.dissipated_energy_j;
            attitude = step.motion.rotation_after;
            momentum = step.motion.momentum_after_body;
        }
        errors.push(
            (momentum - expected_momentum).norm()
                + (attitude.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap() - marker).norm()
                + (work - reference[7]).abs()
                + (bath_work - reference[8]).abs()
                + (heat - reference[9]).abs(),
        );
    }
    assert!(errors[1] < 0.27 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.27 * errors[1], "{errors:?}");
}

#[test]
fn g3_material_and_spatial_frame_changes_preserve_bath_work_and_body_motion() {
    let moments = Vec3::new(2.0, 3.0, 4.0);
    let basis = rotation(Vec3::new(0.4, -0.3, 0.2));
    let world = rotation(Vec3::new(-0.2, 0.5, 0.3));
    let inertia = InertiaTensor3::from_principal(moments, So3::identity()).unwrap();
    let changed_inertia = InertiaTensor3::from_principal(moments, basis).unwrap();
    let damping = full_drag();
    let changed_damping = basis
        .matrix()
        .compose(damping)
        .compose(basis.matrix().transpose());
    let original = inertia.viscous_stepper(damping, 0.025).unwrap();
    let changed = changed_inertia
        .viscous_stepper(changed_damping, 0.025)
        .unwrap();
    let mut attitude = rotation(Vec3::new(-0.3, 0.2, 0.1));
    let mut changed_attitude = world
        .compose(attitude)
        .unwrap()
        .compose(basis.inverse())
        .unwrap();
    let mut momentum = Vec3::new(1.2, -0.9, 2.8);
    let mut changed_momentum = basis.rotate(momentum).unwrap();
    let drive = Vec3::new(0.2, -0.1, 0.3);
    let bath = Vec3::new(-0.3, 0.5, 1.1);
    for _ in 0..16 {
        let a = original
            .step_in_rotating_bath(attitude, momentum, HeldTorque::Space(drive), bath)
            .unwrap();
        let b = changed
            .step_in_rotating_bath(
                changed_attitude,
                changed_momentum,
                HeldTorque::Space(world.rotate(drive).unwrap()),
                world.rotate(bath).unwrap(),
            )
            .unwrap();
        assert!(
            (b.motion.momentum_after_body - basis.rotate(a.motion.momentum_after_body).unwrap())
                .norm()
                < 2e-12
        );
        let expected_attitude = world
            .compose(a.motion.rotation_after)
            .unwrap()
            .compose(basis.inverse())
            .unwrap();
        assert!(
            b.motion
                .rotation_after
                .body_minus(expected_attitude)
                .unwrap()
                .angular
                .norm()
                < 2e-12
        );
        assert!(
            (b.motion.drag_angular_impulse_space
                - world.rotate(a.motion.drag_angular_impulse_space).unwrap())
            .norm()
                < 2e-13
        );
        close(a.bath_work_j, b.bath_work_j, 2e-13);
        close(a.motion.drive_work_j, b.motion.drive_work_j, 2e-13);
        close(
            a.motion.dissipated_energy_j,
            b.motion.dissipated_energy_j,
            2e-13,
        );
        close(
            a.bath_work_j,
            a.motion.drag_angular_impulse_space.dot(bath),
            2e-13,
        );
        assert!((a.motion.drive_angular_impulse_space - drive.scale(0.025)).norm() < 2e-13);
        attitude = a.motion.rotation_after;
        changed_attitude = b.motion.rotation_after;
        momentum = a.motion.momentum_after_body;
        changed_momentum = b.motion.momentum_after_body;
    }
}

#[test]
fn g0_corotation_extraction_stationary_parity_and_invalid_bath_refusal() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 2.0, 2.0), So3::identity()).unwrap();
    let stepper = inertia
        .viscous_stepper(diagonal(0.4, 0.8, 1.2), 0.1)
        .unwrap();
    let zero = HeldTorque::Body(Vec3::default());
    let corotating = stepper
        .step_in_rotating_bath(
            So3::identity(),
            Vec3::new(0.0, 0.0, 4.0),
            zero,
            Vec3::new(0.0, 0.0, 2.0),
        )
        .unwrap();
    close(corotating.motion.dissipated_energy_j, 0.0, 1e-28);
    close(corotating.bath_work_j, 0.0, 1e-14);
    close(corotating.motion.momentum_after_body.z, 4.0, 1e-14);
    for sign in [-1.0, 1.0] {
        let extraction = stepper
            .step_in_rotating_bath(
                So3::identity(),
                Vec3::new(0.0, 0.0, 6.0 * sign),
                zero,
                Vec3::new(0.0, 0.0, 1.2 * sign),
            )
            .unwrap();
        assert!(extraction.bath_work_j < 0.0);
        assert!(extraction.motion.kinetic_energy_change_j < 0.0);
        assert!(extraction.motion.dissipated_energy_j > 0.0);
        close(extraction.motion.work_energy_residual_j, 0.0, 1e-13);
    }
    let attitude = rotation(Vec3::new(0.3, -0.2, 0.1));
    let momentum = Vec3::new(0.2, -0.4, 0.7);
    let drive = HeldTorque::Space(Vec3::new(0.1, 0.2, -0.3));
    let stationary = stepper.step(attitude, momentum, drive).unwrap();
    let still_bath = stepper
        .step_in_rotating_bath(attitude, momentum, drive, Vec3::new(-0.0, 0.0, -0.0))
        .unwrap();
    assert_eq!(stationary, still_bath.motion);
    assert_eq!(still_bath.bath_work_j, 0.0);
    let free = inertia.viscous_stepper(Mat3::zero(), 0.1).unwrap();
    let without_drag = free.step(attitude, momentum, drive).unwrap();
    let moving = free
        .step_in_rotating_bath(attitude, momentum, drive, Vec3::new(2.0, -1.0, 0.5))
        .unwrap();
    assert_eq!(without_drag, moving.motion);
    assert_eq!(moving.bath_work_j, 0.0);
    // A disconnected bath cannot cause spurious work or overflow a frame
    // transform: only its finite input admission remains applicable.
    let disconnected = free
        .step_in_rotating_bath(
            attitude,
            momentum,
            drive,
            Vec3::new(f64::MAX, f64::MAX, f64::MAX),
        )
        .unwrap();
    assert_eq!(without_drag, disconnected.motion);
    assert_eq!(disconnected.bath_work_j, 0.0);
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            stepper
                .step_in_rotating_bath(attitude, momentum, drive, Vec3::new(value, 0.0, 0.0))
                .is_err()
        );
    }
    let bath = Vec3::new(0.2, 0.3, 0.4);
    let accepted = stepper
        .step_in_rotating_bath(attitude, momentum, drive, bath)
        .unwrap();
    assert_eq!(
        accepted,
        stepper
            .step_in_rotating_bath(attitude, momentum, drive, bath)
            .unwrap()
    );
    assert_eq!(momentum, Vec3::new(0.2, -0.4, 0.7));
}
