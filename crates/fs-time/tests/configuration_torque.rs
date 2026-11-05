//! Physical regressions for live torque laws, independent of the splitting.

use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};
use fs_time::{LieStepError, rigid_body_torque_function_step};

fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}

fn close(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() < tolerance, "{a:?} != {b:?}");
}

fn angle_close(a: So3, b: So3, tolerance: f64) {
    assert!(
        a.body_minus(b).unwrap().angular.norm() < tolerance,
        "{a:?} != {b:?}"
    );
}

fn pendulum_torque(_: f64, r: So3) -> Result<HeldTorque, LieStepError> {
    let rod = r.rotate(Vec3::new(1.0, 0.0, 0.0))?;
    Ok(HeldTorque::Space(Vec3::new(0.0, 0.0, -2.0 * rod.y)))
}

// Independent scalar RK4 for theta'' = -sin(theta), using the physical
// inertia I_z=2 and torque -2 sin(theta), not a second copy of the split map.
fn pendulum_reference(duration: f64) -> [f64; 2] {
    let mut y = [1.1, 0.3];
    let h = duration / 20_000.0;
    let derivative = |y: [f64; 2]| [y[1], -y[0].sin()];
    for _ in 0..20_000 {
        let k1 = derivative(y);
        let k2 = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
        let k3 = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
        let k4 = derivative(core::array::from_fn(|i| y[i] + h * k3[i]));
        y = core::array::from_fn(|i| y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}

#[test]
fn nonlinear_torsion_has_second_order_attitude_and_spin() {
    let duration = 3.0;
    let expected = pendulum_reference(duration);
    let run = |n: usize| {
        let mut state = (rotation(Vec3::new(0.0, 0.0, 1.1)), Vec3::new(0.0, 0.0, 0.3));
        let h = duration / n as f64;
        for i in 0..n {
            state = rigid_body_torque_function_step(
                state.0,
                state.1,
                Vec3::new(1.0, 1.5, 2.0),
                i as f64 * h,
                h,
                pendulum_torque,
            )
            .unwrap();
        }
        let rod = state.0.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap();
        (rod.y.atan2(rod.x) - expected[0]).hypot(state.1.z - expected[1])
    };
    let coarse = run(60);
    let fine = run(120);
    assert!(fine < 2.0e-4, "fine error {fine:e}");
    assert!(coarse / fine > 3.8, "refinement {}", coarse / fine);
}

#[test]
fn changing_torque_is_sampled_even_when_initial_load_is_zero() {
    let h = 0.25;
    let mut samples = Vec::new();
    let accepted = rigid_body_torque_function_step(
        So3::identity(),
        Vec3::default(),
        Vec3::new(2.0, 2.0, 2.0),
        0.0,
        h,
        |t, r| {
            samples.push((t, r));
            Ok(HeldTorque::Space(Vec3::new(0.0, 0.0, 4.0 * t)))
        },
    )
    .unwrap();
    assert_eq!(samples.len(), 2);
    assert_eq!(samples[0], (0.0, So3::identity()));
    assert_eq!(samples[1], (h, accepted.0));
    close(accepted.1, Vec3::new(0.0, 0.0, h * h), 1.0e-15);
}

#[test]
fn prescribed_time_ramp_has_exact_impulse_and_second_order_phase() {
    let duration: f64 = 2.0;
    let exact_omega = 0.4 + (0.7 * duration + 0.15 * duration * duration) / 2.0;
    let exact_angle =
        0.2 + 0.4 * duration + 0.7 * duration * duration / 4.0 + 0.3 * duration.powi(3) / 12.0;
    let run = |n: usize| {
        let h = duration / n as f64;
        let mut state = (rotation(Vec3::new(0.0, 0.0, 0.2)), Vec3::new(0.0, 0.0, 0.4));
        for i in 0..n {
            state = rigid_body_torque_function_step(
                state.0,
                state.1,
                Vec3::new(2.0, 2.0, 2.0),
                i as f64 * h,
                h,
                |t, _| Ok(HeldTorque::Space(Vec3::new(0.0, 0.0, 0.7 + 0.3 * t))),
            )
            .unwrap();
        }
        close(state.1, Vec3::new(0.0, 0.0, exact_omega), 3.0e-14);
        let expected = rotation(Vec3::new(0.0, 0.0, exact_angle));
        state.0.body_minus(expected).unwrap().angular.norm()
    };
    let coarse = run(20);
    let fine = run(40);
    assert!(coarse / fine > 3.99, "refinement {}", coarse / fine);
    assert!(fine < 2.0e-4);
}

fn dense_inertia() -> InertiaTensor3 {
    InertiaTensor3::new(Mat3 {
        m: [1.5, 0.5, 0.0, 0.5, 1.5, 0.0, 0.0, 0.0, 3.0],
    })
    .unwrap()
}

fn magnetic_torque(r: So3, field: Vec3) -> Vec3 {
    r.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap().cross(field)
}

#[test]
fn full_inertia_closes_live_endpoint_impulse_in_both_frames() {
    let inertia = dense_inertia();
    let initial = (
        rotation(Vec3::new(0.3, -0.4, 0.2)),
        Vec3::new(0.5, 0.2, -0.7),
    );
    let field = Vec3::new(-0.2, 0.7, 0.4);
    let h = 0.12;
    let l0 = initial
        .0
        .rotate(inertia.momentum(initial.1).unwrap())
        .unwrap();
    let mut reference = None;
    for body_frame in [false, true] {
        let state = inertia
            .torque_function_step(initial.0, initial.1, 1.0, h, |_, r| {
                let torque = magnetic_torque(r, field);
                if body_frame {
                    Ok(HeldTorque::Body(r.inverse().rotate(torque)?))
                } else {
                    Ok(HeldTorque::Space(torque))
                }
            })
            .unwrap();
        let l1 = state.0.rotate(inertia.momentum(state.1).unwrap()).unwrap();
        let impulse =
            (magnetic_torque(initial.0, field) + magnetic_torque(state.0, field)).scale(0.5 * h);
        close(l1 - l0, impulse, 2.0e-14);
        if let Some(other) = reference {
            let other: (So3, Vec3) = other;
            angle_close(state.0, other.0, 2.0e-14);
            close(state.1, other.1, 2.0e-14);
        }
        reference = Some(state);
    }
}

#[test]
fn autonomous_tensor_load_reverses_and_respects_world_frame_covariance() {
    let inertia = dense_inertia();
    let initial = (
        rotation(Vec3::new(0.3, -0.4, 0.2)),
        Vec3::new(0.5, 0.2, -0.7),
    );
    let field = Vec3::new(-0.2, 0.7, 0.4);
    let h = 0.08;
    let force = |_: f64, r| Ok(HeldTorque::Space(magnetic_torque(r, field)));
    let state = inertia
        .torque_function_step(initial.0, initial.1, 0.0, h, force)
        .unwrap();
    assert_eq!(
        state,
        inertia
            .torque_function_step(initial.0, initial.1, 0.0, h, force)
            .unwrap()
    );
    let reversed = inertia
        .torque_function_step(state.0, state.1, h, -h, force)
        .unwrap();
    angle_close(reversed.0, initial.0, 3.0e-14);
    close(reversed.1, initial.1, 3.0e-14);
    let world = rotation(Vec3::new(-0.5, 0.1, 0.4));
    let rotated_field = world.rotate(field).unwrap();
    let transformed = inertia
        .torque_function_step(
            world.compose(initial.0).unwrap(),
            initial.1,
            0.0,
            h,
            |_, r| Ok(HeldTorque::Space(magnetic_torque(r, rotated_field))),
        )
        .unwrap();
    angle_close(transformed.0, world.compose(state.0).unwrap(), 3.0e-14);
    close(transformed.1, state.1, 3.0e-14);
}

#[test]
fn invalid_inputs_and_endpoint_callback_refusal_publish_no_state() {
    let initial = (So3::identity(), Vec3::new(0.2, 0.3, -0.4));
    for (time, h, inertia) in [
        (f64::NAN, 0.1, Vec3::new(1.0, 2.0, 3.0)),
        (0.0, 0.0, Vec3::new(1.0, 2.0, 3.0)),
        (f64::MAX, f64::MAX, Vec3::new(1.0, 2.0, 3.0)),
        (1.0e100, 0.1, Vec3::new(1.0, 2.0, 3.0)),
        (0.0, 0.1, Vec3::new(1.0, -2.0, 3.0)),
    ] {
        let mut calls = 0;
        let result =
            rigid_body_torque_function_step(initial.0, initial.1, inertia, time, h, |_, _| {
                calls += 1;
                Ok(HeldTorque::Body(Vec3::default()))
            });
        assert!(result.is_err());
        assert_eq!(calls, 0);
    }
    let refusal = LieStepError::InvalidParameter {
        context: "force field unavailable",
    };
    let mut calls = 0;
    let result = dense_inertia().torque_function_step(initial.0, initial.1, 0.0, 0.1, |_, _| {
        calls += 1;
        if calls == 2 {
            Err(refusal)
        } else {
            Ok(HeldTorque::Space(Vec3::new(0.1, 0.2, 0.3)))
        }
    });
    assert_eq!(result, Err(refusal));
    assert_eq!(calls, 2);
    assert_eq!(initial, (So3::identity(), Vec3::new(0.2, 0.3, -0.4)));
    assert!(
        rigid_body_torque_function_step(
            initial.0,
            initial.1,
            Vec3::new(1.0, 2.0, 3.0),
            0.0,
            0.1,
            |_, _| Ok(HeldTorque::Body(Vec3::new(f64::INFINITY, 0.0, 0.0))),
        )
        .is_err()
    );
}
