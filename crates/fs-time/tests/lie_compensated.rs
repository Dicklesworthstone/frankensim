//! Physical oracles for retained angular impulses in the free rigid flow.

use fs_ga::{Quat, So3, So3Tangent, Vec3};
use fs_math::dd::Dd;
use fs_time::lie::compensated::{rigid_body_momentum_step, rigid_body_momentum_step_in_frame};
use fs_time::lie::{rigid_body_spatial_torque_step, rigid_body_step, rigid_body_torque_step};

const LARGE: f64 = 9_007_199_254_740_992.0;

fn vector(values: [f64; 3]) -> Vec3 {
    Vec3::new(values[0], values[1], values[2])
}

fn values(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn attitude(v: [f64; 3]) -> So3 {
    So3::exp(So3Tangent::new(vector(v))).unwrap()
}

fn attitude_error(a: So3, b: So3) -> f64 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        .into_iter()
        .map(|probe| {
            let a = values(a.rotate(vector(probe)).unwrap());
            let b = values(b.rotate(vector(probe)).unwrap());
            (0..3).map(|d| (a[d] - b[d]).abs()).sum::<f64>()
        })
        .fold(0.0, f64::max)
}

// Independent scalar angle addition: neither the production DD trigonometry
// nor its Lie splitting participates in this one-radian phase oracle.
fn large_axis_attitude(axis: usize, sign: f64) -> So3 {
    let (base_sine, base_cosine) = (sign * 0.5 * LARGE).sin_cos();
    let (tail_sine, tail_cosine) = (sign * 0.5).sin_cos();
    let sine = base_sine * tail_cosine + base_cosine * tail_sine;
    let cosine = base_cosine * tail_cosine - base_sine * tail_sine;
    let mut xyz = [0.0; 3];
    xyz[axis] = sine;
    So3::try_from_quat(Quat {
        w: cosine,
        x: xyz[0],
        y: xyz[1],
        z: xyz[2],
    })
    .unwrap()
}

#[test]
fn spherical_midpoint_tail_changes_attitude_by_a_full_radian() {
    for sign in [-1.0, 1.0] {
        let momentum = [Dd::ZERO, Dd::ZERO, Dd::from_pair(sign * LARGE, sign)];
        let inertia = vector([1.0; 3]);
        let (next, retained) =
            rigid_body_momentum_step(So3::identity(), momentum, inertia, 1.0).unwrap();
        assert_eq!(retained, momentum, "spherical flow must leave L unchanged");
        assert!(attitude_error(next, large_axis_attitude(2, sign)) < 3.0e-15);

        let rounded = [Dd::ZERO, Dd::ZERO, Dd::from_f64(sign * LARGE)];
        let (rounded_attitude, _) =
            rigid_body_momentum_step(So3::identity(), rounded, inertia, 1.0).unwrap();
        assert!(attitude_error(next, rounded_attitude) > 0.9);
        assert_eq!(
            (retained[2] + Dd::from_f64(sign)).to_f64(),
            sign * (LARGE + 2.0)
        );
        assert_eq!((retained[2] - Dd::from_f64(sign)).to_f64(), sign * LARGE);
    }
}

#[test]
fn full_tensor_principal_axis_retains_the_same_resolved_phase() {
    let principal_to_body = attitude([0.37, 0.0, 0.0]);
    let inertia = vector([1.0, 1.5, 2.0]);
    // Rotating unequal y/z moments produces a genuine off-diagonal yz tensor,
    // while the x principal axis stays exact and has the same analytic oracle.
    for sign in [-1.0, 1.0] {
        let momentum = [Dd::from_pair(sign * LARGE, sign), Dd::ZERO, Dd::ZERO];
        let (next, retained) = rigid_body_momentum_step_in_frame(
            So3::identity(),
            momentum,
            inertia,
            principal_to_body,
            1.0,
        )
        .unwrap();
        assert!(attitude_error(next, large_axis_attitude(0, sign)) < 4.0e-15);
        for d in 0..3 {
            assert!((retained[d] - momentum[d]).to_f64().abs() < 2.0e-15);
        }
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn quaternion_product(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

fn inverse_inertia(moments: [f64; 3], frame: So3) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for axis in 0..3 {
        let mut unit = [0.0; 3];
        unit[axis] = 1.0;
        let direction = values(frame.rotate(vector(unit)).unwrap());
        for row in 0..3 {
            for col in 0..3 {
                result[row][col] += direction[row] * direction[col] / moments[axis];
            }
        }
    }
    result
}

fn apply(matrix: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * v[0] + row[1] * v[1] + row[2] * v[2])
}

// Independent continuous Euler equations and their first variation. This
// reference uses small-step classical RK4, with no Hamiltonian split, DD
// arithmetic, group exponential or production momentum transport.
fn derivative(state: [f64; 14], inverse: [[f64; 3]; 3]) -> [f64; 14] {
    let l = [state[0], state[1], state[2]];
    let q = [state[3], state[4], state[5], state[6]];
    let dl = [state[7], state[8], state[9]];
    let dq = [state[10], state[11], state[12], state[13]];
    let omega = apply(inverse, l);
    let domega = apply(inverse, dl);
    let ldot = cross(l, omega);
    let qdot = quaternion_product(q, [0.0, omega[0], omega[1], omega[2]]);
    let dldot1 = cross(dl, omega);
    let dldot2 = cross(l, domega);
    let dqdot1 = quaternion_product(dq, [0.0, omega[0], omega[1], omega[2]]);
    let dqdot2 = quaternion_product(q, [0.0, domega[0], domega[1], domega[2]]);
    let mut result = [0.0; 14];
    for d in 0..3 {
        result[d] = ldot[d];
        result[7 + d] = dldot1[d] + dldot2[d];
    }
    for d in 0..4 {
        result[3 + d] = 0.5 * qdot[d];
        result[10 + d] = 0.5 * (dqdot1[d] + dqdot2[d]);
    }
    result
}

fn rk4(mut state: [f64; 14], inverse: [[f64; 3]; 3], time: f64, steps: usize) -> [f64; 14] {
    let h = time / steps as f64;
    for _ in 0..steps {
        let k1 = derivative(state, inverse);
        let k2 = derivative(
            core::array::from_fn(|d| state[d] + 0.5 * h * k1[d]),
            inverse,
        );
        let k3 = derivative(
            core::array::from_fn(|d| state[d] + 0.5 * h * k2[d]),
            inverse,
        );
        let k4 = derivative(core::array::from_fn(|d| state[d] + h * k3[d]), inverse);
        for d in 0..14 {
            state[d] += h * (k1[d] + 2.0 * k2[d] + 2.0 * k3[d] + k4[d]) / 6.0;
        }
    }
    state
}

fn ode_state(momentum: [f64; 3], rotation: So3, variation: [f64; 3]) -> [f64; 14] {
    let q = rotation.as_quat();
    [
        momentum[0],
        momentum[1],
        momentum[2],
        q.w,
        q.x,
        q.y,
        q.z,
        variation[0],
        variation[1],
        variation[2],
        0.0,
        0.0,
        0.0,
        0.0,
    ]
}

fn ode_rotation(state: [f64; 14]) -> So3 {
    So3::try_from_quat(Quat {
        w: state[3],
        x: state[4],
        y: state[5],
        z: state[6],
    })
    .unwrap()
}

#[test]
fn asymmetric_full_tensor_converges_to_euler_equations_and_preserves_momentum() {
    let frame = attitude([0.3, -0.6, 0.15]);
    let initial_rotation = attitude([-0.2, 0.4, 0.1]);
    let moments = [2.0, 3.0, 4.0];
    let initial_momentum = [0.7, -0.4, 0.9];
    let reference = rk4(
        ode_state(initial_momentum, initial_rotation, [0.0; 3]),
        inverse_inertia(moments, frame),
        0.8,
        4096,
    );
    let spatial_initial = values(initial_rotation.rotate(vector(initial_momentum)).unwrap());
    let norm_initial = initial_momentum
        .map(Dd::from_f64)
        .map(|l| l * l)
        .into_iter()
        .fold(Dd::ZERO, |a, b| a + b);
    let mut errors = [0.0; 3];
    for (sample, steps) in [8, 16, 32].into_iter().enumerate() {
        let mut rotation = initial_rotation;
        let mut momentum = initial_momentum.map(Dd::from_f64);
        for _ in 0..steps {
            (rotation, momentum) = rigid_body_momentum_step_in_frame(
                rotation,
                momentum,
                vector(moments),
                frame,
                0.8 / steps as f64,
            )
            .unwrap();
        }
        errors[sample] = attitude_error(rotation, ode_rotation(reference));
        for d in 0..3 {
            errors[sample] += (momentum[d].to_f64() - reference[d]).abs();
        }
        let spatial = values(rotation.rotate(vector(momentum.map(Dd::to_f64))).unwrap());
        for d in 0..3 {
            assert!((spatial[d] - spatial_initial[d]).abs() < 2.0e-14);
        }
        let norm = momentum
            .map(|l| l * l)
            .into_iter()
            .fold(Dd::ZERO, |a, b| a + b);
        assert!((norm - norm_initial).to_f64().abs() < 5.0e-27);
    }
    for d in 0..2 {
        assert!(errors[d] / errors[d + 1] > 3.8, "errors: {errors:?}");
        assert!(errors[d] / errors[d + 1] < 4.2, "errors: {errors:?}");
    }
}

#[test]
fn off_axis_impulse_tail_drives_attitude_and_nonlinear_spin() {
    for frame in [So3::identity(), attitude([0.3, -0.6, 0.15])] {
        let scaled_moments = vector([2.0 * LARGE, 3.0 * LARGE, 4.0 * LARGE]);
        let normalized_momentum = [1.0, 1.25, -1.5];
        let tail = [1.0, 0.5, -0.25];
        let base = normalized_momentum.map(|l| Dd::from_f64(LARGE * l));
        let perturbed = core::array::from_fn(|d| base[d] + Dd::from_f64(tail[d]));
        let (increment, _) =
            rigid_body_momentum_step_in_frame(So3::identity(), base, scaled_moments, frame, 0.01)
                .unwrap();
        // Cancel the large base attitude so its sub-ULP variation is exposed
        // in ordinary f64 output components instead of rounded away at O(1).
        let initial_rotation = increment.inverse();
        let (base_rotation, base_momentum) =
            rigid_body_momentum_step_in_frame(initial_rotation, base, scaled_moments, frame, 0.01)
                .unwrap();
        let (next_rotation, next_momentum) = rigid_body_momentum_step_in_frame(
            initial_rotation,
            perturbed,
            scaled_moments,
            frame,
            0.01,
        )
        .unwrap();
        let reference = rk4(
            ode_state(normalized_momentum, initial_rotation, tail),
            inverse_inertia([2.0, 3.0, 4.0], frame),
            0.01,
            256,
        );
        let base_q = base_rotation.as_quat();
        let next_q = next_rotation.as_quat();
        for (d, (next, base)) in [next_q.x, next_q.y, next_q.z]
            .into_iter()
            .zip([base_q.x, base_q.y, base_q.z])
            .enumerate()
        {
            let response = (next - base) * LARGE;
            assert!(
                (response - reference[11 + d]).abs() < 8.0e-7,
                "axis {d}, attitude response {response}, reference {}",
                reference[11 + d]
            );
            let spin_response = (next_momentum[d] - base_momentum[d]).to_f64();
            assert!(
                (spin_response - reference[7 + d]).abs() < 8.0e-7,
                "axis {d}, spin response {spin_response}, reference {}",
                reference[7 + d]
            );
        }
    }
}

#[test]
fn compensated_full_tensor_flow_reverses_without_projecting_momentum() {
    let frame = attitude([0.3, -0.6, 0.15]);
    let initial_rotation = attitude([-0.2, 0.4, 0.1]);
    let inertia = vector([2.0, 3.0, 4.0]);
    let initial = [
        Dd::from_pair(0.7, 1.0e-17),
        Dd::from_pair(-0.4, -1.0e-17),
        Dd::from_pair(0.9, 2.0e-17),
    ];
    let (rotation, momentum) =
        rigid_body_momentum_step_in_frame(initial_rotation, initial, inertia, frame, 0.3).unwrap();
    let (rotation, momentum) =
        rigid_body_momentum_step_in_frame(rotation, momentum, inertia, frame, -0.3).unwrap();
    assert!(attitude_error(rotation, initial_rotation) < 2.0e-15);
    for d in 0..3 {
        assert!((momentum[d] - initial[d]).to_f64().abs() < 5.0e-29);
    }
}

#[test]
fn finite_large_quotient_and_subnormal_duration_keep_resolved_rotation() {
    let momentum = [Dd::from_f64(f64::MAX), Dd::ZERO, Dd::ZERO];
    let (rotation, retained) =
        rigid_body_momentum_step(So3::identity(), momentum, vector([3.0; 3]), 1.0).unwrap();
    assert_eq!(retained, momentum);
    assert!(rotation.as_quat().w.is_finite());

    let h = f64::from_bits(1);
    let (rotation, _) =
        rigid_body_momentum_step(So3::identity(), momentum, vector([1.0; 3]), h).unwrap();
    let expected_half_angle = (0.5 * f64::MAX) * h;
    assert!(expected_half_angle > 0.0);
    assert!((rotation.as_quat().x - expected_half_angle).abs() < 1.0e-30);
    assert!(rigid_body_momentum_step(So3::identity(), momentum, vector([1.0; 3]), 0.0).is_err());
    let invalid = [
        Dd {
            hi: 1.0,
            lo: f64::NAN,
        },
        Dd::ZERO,
        Dd::ZERO,
    ];
    assert!(rigid_body_momentum_step(So3::identity(), invalid, vector([1.0; 3]), 1.0).is_err());
}

#[test]
fn public_held_torque_owners_retain_midpoint_phase_and_complete_impulse() {
    for step in [rigid_body_torque_step, rigid_body_spatial_torque_step] {
        for sign in [-1.0, 1.0] {
            let omega = vector([0.0, 0.0, sign * LARGE]);
            let torque = vector([0.0, 0.0, sign * 2.0]);
            let inertia = vector([1.0; 3]);
            let (rotation, next_omega) =
                step(So3::identity(), omega, inertia, torque, 1.0).unwrap();
            assert_eq!(next_omega.z, sign * (LARGE + 2.0));
            assert!(attitude_error(rotation, large_axis_attitude(2, sign)) < 4.0e-15);
            let (returned_rotation, returned_omega) =
                step(rotation, next_omega, inertia, torque, -1.0).unwrap();
            assert_eq!(returned_omega, omega);
            assert!(attitude_error(returned_rotation, So3::identity()) < 3.0e-15);
        }
    }
}

#[test]
fn held_angular_impulse_survives_a_subnormal_duration() {
    let h = f64::from_bits(1);
    let torque = vector([0.0, 0.0, 1.0e308]);
    let expected = torque.z * h;
    assert!(expected > 0.0);
    for step in [rigid_body_torque_step, rigid_body_spatial_torque_step] {
        let (rotation, omega) = step(
            So3::identity(),
            vector([0.0; 3]),
            vector([1.0; 3]),
            torque,
            h,
        )
        .unwrap();
        assert_eq!(omega.z, expected);
        assert_eq!(rotation, So3::identity());
    }
}

#[test]
fn zero_held_torques_keep_exact_legacy_free_drift_and_refusals() {
    for (rotation, omega, inertia) in [
        (
            attitude([0.2, -0.3, 0.4]),
            vector([0.4, 0.7, 1.1]),
            vector([1.0, 2.0, 3.0]),
        ),
        (So3::identity(), vector([0.0, 0.0, 0.1]), vector([3.0; 3])),
    ] {
        for h in [0.03, -0.03, 0.0, f64::INFINITY, f64::NAN] {
            let free = rigid_body_step(rotation, omega, inertia, h);
            for zero in [vector([0.0; 3]), vector([-0.0, 0.0, -0.0])] {
                assert_eq!(
                    rigid_body_torque_step(rotation, omega, inertia, zero, h),
                    free
                );
                assert_eq!(
                    rigid_body_spatial_torque_step(rotation, omega, inertia, zero, h),
                    free
                );
            }
        }
    }
}

fn spatial_momentum(rotation: So3, omega: Vec3, inertia: Vec3) -> [f64; 3] {
    values(
        rotation
            .rotate(vector([
                inertia.x * omega.x,
                inertia.y * omega.y,
                inertia.z * omega.z,
            ]))
            .unwrap(),
    )
}

#[test]
fn held_torque_owners_preserve_their_distinct_endpoint_frame_balances() {
    let initial_rotation = attitude([0.2, -0.4, 0.3]);
    let omega = vector([0.5, -0.2, 0.7]);
    let inertia = vector([2.0, 3.0, 4.0]);
    let torque = vector([0.15, -0.25, 0.35]);
    let h = 0.3;
    let initial_momentum = spatial_momentum(initial_rotation, omega, inertia);
    let (body_rotation, body_omega) =
        rigid_body_torque_step(initial_rotation, omega, inertia, torque, h).unwrap();
    let body_momentum = spatial_momentum(body_rotation, body_omega, inertia);
    let body_start_torque = values(initial_rotation.rotate(torque).unwrap());
    let body_end_torque = values(body_rotation.rotate(torque).unwrap());
    let (space_rotation, space_omega) =
        rigid_body_spatial_torque_step(initial_rotation, omega, inertia, torque, h).unwrap();
    let space_momentum = spatial_momentum(space_rotation, space_omega, inertia);
    let spatial_torque = values(torque);
    for d in 0..3 {
        let body_impulse = 0.5 * h * (body_start_torque[d] + body_end_torque[d]);
        assert!((body_momentum[d] - initial_momentum[d] - body_impulse).abs() < 2.0e-15);
        assert!((space_momentum[d] - initial_momentum[d] - h * spatial_torque[d]).abs() < 2.0e-15);
    }
    assert!(attitude_error(body_rotation, space_rotation) > 1.0e-4);
    for (step, rotation, next_omega) in [
        (
            rigid_body_torque_step as fn(So3, Vec3, Vec3, Vec3, f64) -> _,
            body_rotation,
            body_omega,
        ),
        (rigid_body_spatial_torque_step, space_rotation, space_omega),
    ] {
        let (returned_rotation, returned_omega) =
            step(rotation, next_omega, inertia, torque, -h).unwrap();
        assert!(attitude_error(returned_rotation, initial_rotation) < 2.0e-15);
        for (actual, expected) in values(returned_omega).into_iter().zip(values(omega)) {
            assert!((actual - expected).abs() < 8.0e-16);
        }
    }
}
