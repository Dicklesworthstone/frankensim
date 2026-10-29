use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3, driven_step, momentum_step};

fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}

fn vector_close(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() < tolerance, "{a:?} != {b:?}");
}

fn rotation_close(a: So3, b: So3, tolerance: f64) {
    assert!(a.body_minus(b).unwrap().angular.norm() < tolerance);
}

fn fixture() -> (InertiaTensor3, So3, Vec3) {
    let frame = rotation(Vec3::new(0.3, -0.7, 0.2));
    let principal = Vec3::new(1.2, 2.3, 3.1);
    let tensor = InertiaTensor3::from_principal(principal, frame).unwrap();
    (
        InertiaTensor3::new(tensor.matrix()).unwrap(),
        frame,
        principal,
    )
}

#[test]
fn full_inertia_reconstruction_inverse_and_extreme_scales() {
    let (tensor, frame, moments) = fixture();
    assert!(tensor.relative_reconstruction_error() < 4.0e-15);
    let omega = Vec3::new(0.8, -0.6, 0.2);
    let local = frame.inverse().rotate(omega).unwrap();
    let expected = frame
        .rotate(Vec3::new(
            moments.x * local.x,
            moments.y * local.y,
            moments.z * local.z,
        ))
        .unwrap();
    vector_close(tensor.momentum(omega).unwrap(), expected, 4.0e-15);
    vector_close(tensor.velocity(expected).unwrap(), omega, 4.0e-15);
    for scale in [1.0e-180, 1.0e180] {
        let scaled = InertiaTensor3::new(Mat3 {
            m: tensor.matrix().m.map(|x| x * scale),
        })
        .unwrap();
        vector_close(
            scaled.momentum(omega).unwrap().scale(1.0 / scale),
            expected,
            8.0e-15,
        );
        vector_close(
            scaled.velocity(expected.scale(scale)).unwrap(),
            omega,
            8.0e-15,
        );
    }
}

#[test]
fn diagonal_steps_retain_the_existing_arithmetic() {
    let principal = Vec3::new(3.0, 1.0, 2.0);
    let tensor = InertiaTensor3::new(Mat3 {
        m: [3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0],
    })
    .unwrap();
    let r = rotation(Vec3::new(0.2, 0.5, -0.1));
    let l = Vec3::new(0.8, 0.4, -0.7);
    let h = Vec3::new(-0.2, 0.3, 0.6);
    for dt in [0.01, -0.03, 0.5] {
        assert_eq!(
            tensor.momentum_step(r, l, h, dt).unwrap(),
            momentum_step(r, l, principal, h, dt).unwrap()
        );
        let torque = HeldTorque::Space(Vec3::new(0.1, -0.7, 0.3));
        assert_eq!(
            tensor
                .driven_step(r, l, h, h.scale(1.1), torque, dt)
                .unwrap(),
            driven_step(r, l, principal, h, h.scale(1.1), torque, dt).unwrap()
        );
    }
}

#[test]
fn cached_frame_matches_native_principal_motion_and_replays() {
    let (tensor, frame, moments) = fixture();
    let initial_rotation = rotation(Vec3::new(0.2, 0.1, -0.3));
    let initial_l = Vec3::new(0.8, -0.6, 0.7);
    let h = Vec3::new(0.3, 0.2, -0.1);
    let run = || {
        let mut state = (initial_rotation, initial_l);
        for _ in 0..400 {
            state = tensor.momentum_step(state.0, state.1, h, 0.005).unwrap();
        }
        state
    };
    let actual = run();
    assert_eq!(actual, run());
    let mut reference = (
        initial_rotation.compose(frame).unwrap(),
        frame.inverse().rotate(initial_l).unwrap(),
    );
    let h_principal = frame.inverse().rotate(h).unwrap();
    for _ in 0..400 {
        reference = momentum_step(reference.0, reference.1, moments, h_principal, 0.005).unwrap();
    }
    rotation_close(
        actual.0,
        reference.0.compose(frame.inverse()).unwrap(),
        2.0e-12,
    );
    vector_close(actual.1, frame.rotate(reference.1).unwrap(), 2.0e-12);
    vector_close(
        actual.0.rotate(actual.1).unwrap(),
        initial_rotation.rotate(initial_l).unwrap(),
        2.0e-12,
    );
    let reversed = tensor.momentum_step(actual.0, actual.1, h, -0.005).unwrap();
    let recovered = tensor
        .momentum_step(reversed.0, reversed.1, h, 0.005)
        .unwrap();
    rotation_close(recovered.0, actual.0, 5.0e-14);
    vector_close(recovered.1, actual.1, 5.0e-14);
}

#[test]
fn long_free_trajectory_does_not_repeat_principal_frame_roundoff() {
    let (tensor, _, _) = fixture();
    let h = Vec3::new(0.3, -0.2, 0.7);
    let initial = (
        rotation(Vec3::new(0.2, 0.1, -0.3)),
        Vec3::new(0.8, -0.6, 0.7),
    );
    let expected = initial.0.rotate(initial.1).unwrap();
    let mut state = initial;
    let mut maximum = 0.0_f64;
    for _ in 0..10_000 {
        state = tensor.momentum_step(state.0, state.1, h, 0.001).unwrap();
        maximum = maximum.max((state.0.rotate(state.1).unwrap() - expected).norm());
    }
    eprintln!("full-inertia 10000-step spatial-momentum drift: {maximum:e}");
    assert!(maximum < 2.0e-12);
}

#[test]
fn driven_tensor_closes_external_impulses_and_internal_recoil() {
    let (tensor, _, _) = fixture();
    let r0 = rotation(Vec3::new(0.2, -0.4, 0.1));
    let w0 = Vec3::new(0.3, 0.6, -0.8);
    let h0 = Vec3::new(0.2, -0.1, 0.3);
    let h1 = Vec3::new(-0.4, 0.5, 0.1);
    let torque = Vec3::new(0.4, -0.2, 0.7);
    let dt = 0.07;
    let l0 = r0.rotate(tensor.momentum(w0).unwrap() + h0).unwrap();
    for load in [
        HeldTorque::Space(torque),
        HeldTorque::Body(torque),
        HeldTorque::Space(Vec3::default()),
    ] {
        let (r1, w1) = tensor.driven_step(r0, w0, h0, h1, load, dt).unwrap();
        let l1 = r1.rotate(tensor.momentum(w1).unwrap() + h1).unwrap();
        let impulse = match load {
            HeldTorque::Space(value) => value.scale(dt),
            HeldTorque::Body(value) => {
                (r0.rotate(value).unwrap() + r1.rotate(value).unwrap()).scale(0.5 * dt)
            }
        };
        vector_close(l1 - l0, impulse, 2.0e-14);
        let reversed = tensor.driven_step(r1, w1, h1, h0, load, -dt).unwrap();
        rotation_close(reversed.0, r0, 3.0e-14);
        vector_close(reversed.1, w0, 3.0e-14);
    }
}

#[test]
fn world_frame_covariance_and_dense_euler_refinement() {
    let (tensor, _, _) = fixture();
    let r0 = rotation(Vec3::new(0.1, 0.2, -0.3));
    let world = rotation(Vec3::new(-0.4, 0.2, 0.5));
    let w0 = Vec3::new(0.6, -0.4, 0.8);
    let h = Vec3::new(0.3, -0.2, 0.7);
    let torque = Vec3::new(0.2, 0.4, -0.3);
    let a = tensor
        .driven_step(r0, w0, h, h, HeldTorque::Space(torque), 0.02)
        .unwrap();
    let b = tensor
        .driven_step(
            world.compose(r0).unwrap(),
            w0,
            h,
            h,
            HeldTorque::Space(world.rotate(torque).unwrap()),
            0.02,
        )
        .unwrap();
    rotation_close(b.0, world.compose(a.0).unwrap(), 2.0e-14);
    vector_close(a.1, b.1, 2.0e-14);

    // Independent dense inverse and ordinary RK4 Euler equation: no cached
    // eigensolver frame or native gyrostat drift is used by the reference.
    let m = tensor.matrix().m;
    let determinant = m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6]);
    let inverse = Mat3 {
        m: [
            m[4] * m[8] - m[5] * m[7],
            m[2] * m[7] - m[1] * m[8],
            m[1] * m[5] - m[2] * m[4],
            m[5] * m[6] - m[3] * m[8],
            m[0] * m[8] - m[2] * m[6],
            m[2] * m[3] - m[0] * m[5],
            m[3] * m[7] - m[4] * m[6],
            m[1] * m[6] - m[0] * m[7],
            m[0] * m[4] - m[1] * m[3],
        ]
        .map(|x| x / determinant),
    };
    let l0 = tensor.matrix().apply(w0) + h;
    let derivative = |l: Vec3| l.cross(inverse.apply(l - h));
    let reference_dt = 2.0 / 8192.0;
    let mut reference = l0;
    for _ in 0..8192 {
        let k1 = derivative(reference);
        let k2 = derivative(reference + k1.scale(reference_dt * 0.5));
        let k3 = derivative(reference + k2.scale(reference_dt * 0.5));
        let k4 = derivative(reference + k3.scale(reference_dt));
        reference = reference + (k1 + k2.scale(2.0) + k3.scale(2.0) + k4).scale(reference_dt / 6.0);
    }
    let mut errors = Vec::new();
    for steps in [50, 100, 200] {
        let mut state = (r0, l0);
        for _ in 0..steps {
            state = tensor
                .momentum_step(state.0, state.1, h, 2.0 / f64::from(steps))
                .unwrap();
        }
        errors.push((state.1 - reference).norm());
        vector_close(
            state.0.rotate(state.1).unwrap(),
            r0.rotate(l0).unwrap(),
            2.0e-12,
        );
    }
    eprintln!("full-inertia Euler errors: {errors:?}");
    assert!((3.8..4.2).contains(&(errors[0] / errors[1])));
    assert!((3.8..4.2).contains(&(errors[1] / errors[2])));
}

#[test]
fn invalid_tensors_and_unrepresentable_states_refuse() {
    for m in [
        [0.0; 9],
        [1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        [1.0, 0.2, 0.0, 0.3, 2.0, 0.0, 0.0, 0.0, 3.0],
        [1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 2.0],
        [f64::NAN, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0],
    ] {
        assert!(InertiaTensor3::new(Mat3 { m }).is_err());
    }
    let (tensor, _, _) = fixture();
    assert!(tensor.momentum(Vec3::new(f64::INFINITY, 0.0, 0.0)).is_err());
    assert!(tensor.velocity(Vec3::new(f64::NAN, 0.0, 0.0)).is_err());
    assert!(
        tensor
            .momentum_step(So3::identity(), Vec3::default(), Vec3::default(), 0.0)
            .is_err()
    );
    let tiny =
        InertiaTensor3::from_principal(Vec3::new(1.0e-310, 1.0, 1.0), So3::identity()).unwrap();
    assert!(tiny.velocity(Vec3::new(1.0, 0.0, 0.0)).is_err());
}
