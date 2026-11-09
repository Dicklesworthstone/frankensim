//! G0/G1/G3: native material-point acceleration reaches six-DOF perception.

use fs_flyer::adapter::{mbd_to_vec, state_to_mbd, vec_to_mbd};
use fs_flyer::perception::{CueSpec, PERCEPTION_HZ, RigidPointCueInput, perception_v1};
use fs_flyer::spine::{
    Loads, RigidBody, SixDofState, admit_and_canonicalize_state, step, tick_digest,
};
use fs_ga::{So3, So3Tangent, Vec3};
use fs_mbd::Vec3 as MbdVec3;

const G: f64 = 9.80665;

fn rest() -> SixDofState {
    SixDofState {
        pos_m: [0.0; 3],
        vel_mps: [0.0; 3],
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: [0.0; 3],
    }
}

fn close(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for axis in 0..3 {
        assert!(
            (actual[axis] - expected[axis]).abs() < tolerance,
            "axis {axis}: {actual:?} != {expected:?} within {tolerance}"
        );
    }
}

fn cue_input() -> RigidPointCueInput {
    RigidPointCueInput {
        point_body_m: [0.0; 3],
        gravity_world_mps2: [0.0, 0.0, G],
        pitch_relative_rad: 0.0,
        roll_relative_rad: 0.0,
        heave_reference_mps2: -G,
    }
}

#[test]
fn free_fall_and_supported_specific_force_keep_gravity_ownership() {
    let body = RigidBody::new(2.5, [2.0, 3.0, 4.0]).unwrap();
    let gravity = MbdVec3::new(0.0, 0.0, G);
    for (force, acceleration, specific) in [
        ([0.0, 0.0, body.mass_kg * G], [0.0, 0.0, G], [0.0; 3]),
        ([0.0; 3], [0.0; 3], [0.0, 0.0, -G]),
    ] {
        let point = body
            .prepare_point_motion(
                &rest(),
                Loads {
                    force_n: force,
                    moment_nm: [0.0; 3],
                },
            )
            .unwrap()
            .point_motion(MbdVec3::ZERO)
            .unwrap();
        close(
            mbd_to_vec(point.point_acceleration_world),
            acceleration,
            2e-14,
        );
        close(
            mbd_to_vec(point.specific_acceleration_world(gravity).unwrap()),
            specific,
            2e-14,
        );
    }
}

#[test]
fn off_center_point_retains_tangential_and_centripetal_acceleration() {
    let body = RigidBody::new(2.0, [2.0; 3]).unwrap();
    let state = SixDofState {
        pos_m: [1.0, 2.0, 3.0],
        vel_mps: [4.0, 5.0, 6.0],
        omega_body: [0.0, 0.0, 3.0],
        ..rest()
    };
    let loads = Loads {
        force_n: [2.0, -4.0, 6.0],
        moment_nm: [0.0, 0.0, 4.0],
    };
    let before = tick_digest(0, &state);
    let prepared = body.prepare_point_motion(&state, loads).unwrap();
    let point = prepared.point_motion(MbdVec3::new(2.0, 0.0, 0.0)).unwrap();
    close(
        mbd_to_vec(point.kinematics.point_world),
        [3.0, 2.0, 3.0],
        1e-14,
    );
    close(
        mbd_to_vec(point.kinematics.point_velocity_world),
        [4.0, 11.0, 6.0],
        1e-14,
    );
    close(
        mbd_to_vec(point.point_acceleration_world),
        [-17.0, 2.0, 3.0],
        1e-14,
    );
    close(
        mbd_to_vec(point.angular_acceleration_body),
        [0.0, 0.0, 2.0],
        1e-14,
    );
    assert_eq!(
        before,
        tick_digest(0, &state),
        "observation must not advance flight"
    );
    assert_eq!(prepared, body.prepare_point_motion(&state, loads).unwrap());
}

#[test]
fn physical_observation_is_invariant_to_the_quaternion_double_cover() {
    let body = RigidBody::new(2.0, [2.0, 3.0, 4.0]).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.7, -0.4, 1.1))).unwrap();
    let q = rotation.as_quat();
    let loads = Loads {
        force_n: [2.0, -4.0, 6.0],
        moment_nm: [0.2, -0.1, 0.4],
    };
    let arm = MbdVec3::new(0.8, -0.2, 0.6);
    let mut model = perception_v1(23);
    model.cues = [CueSpec {
        delay_ticks: 1,
        filter_tau_s: 0.04,
        gain: 1.0,
        remnant_sigma: 0.0,
    }; 6];
    for quaternion in [
        [q.w, q.x, q.y, q.z],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ] {
        // The second scale is inside the spine's admission but outside the
        // native replay transport's tighter norm tolerance. Physical
        // observation must accept every admitted representative without
        // weakening replay or changing the source state.
        for scale in [1.0, 1.0 + 2e-13] {
            let positive = SixDofState {
                pos_m: [1.0, 2.0, 3.0],
                vel_mps: [4.0, 5.0, 6.0],
                quat: quaternion.map(|value| value * scale),
                omega_body: [0.4, -0.3, 0.7],
            };
            let negative = SixDofState {
                quat: positive.quat.map(|value| -value),
                ..positive
            };
            let before = [tick_digest(0, &positive), tick_digest(0, &negative)];
            for state in [&positive, &negative] {
                admit_and_canonicalize_state(&body, state, 0.005).unwrap();
            }
            assert!(state_to_mbd(&negative).is_err());
            if scale != 1.0 {
                assert!(state_to_mbd(&positive).is_err());
            }
            let a = body.prepare_point_motion(&positive, loads).unwrap();
            let b = body.prepare_point_motion(&negative, loads).unwrap();
            assert_eq!(a, b);
            let point_a = a.point_motion(arm).unwrap();
            let point_b = b.point_motion(arm).unwrap();
            assert_eq!(point_a, point_b);
            assert_eq!(
                point_a
                    .specific_acceleration_world(MbdVec3::new(0.0, 0.0, G))
                    .unwrap(),
                point_b
                    .specific_acceleration_world(MbdVec3::new(0.0, 0.0, G))
                    .unwrap(),
            );
            let input = RigidPointCueInput {
                point_body_m: mbd_to_vec(arm),
                ..cue_input()
            };
            let mut state_a = model.init().unwrap();
            let mut state_b = model.init().unwrap();
            for _ in 0..3 {
                assert_eq!(
                    model
                        .step_with_rigid_motion(&mut state_a, &a, input)
                        .unwrap(),
                    model
                        .step_with_rigid_motion(&mut state_b, &b, input)
                        .unwrap(),
                );
            }
            assert_eq!(
                before,
                [tick_digest(0, &positive), tick_digest(0, &negative)],
            );
        }
    }
}

#[test]
fn euler_acceleration_includes_the_torque_free_gyroscopic_term() {
    let body = RigidBody::new(2.0, [2.0, 3.0, 4.0]).unwrap();
    let state = SixDofState {
        omega_body: [1.0, 2.0, 3.0],
        ..rest()
    };
    let point = body
        .prepare_point_motion(
            &state,
            Loads {
                force_n: [0.0; 3],
                moment_nm: [0.7, -0.2, 0.9],
            },
        )
        .unwrap()
        .point_motion(MbdVec3::ZERO)
        .unwrap();
    // Independent principal-axis Euler equations, not tau/I alone.
    close(
        mbd_to_vec(point.angular_acceleration_body),
        [(0.7 - 6.0) / 2.0, (-0.2 + 6.0) / 3.0, (0.9 - 2.0) / 4.0],
        2e-14,
    );
}

#[test]
fn full_tensor_point_motion_is_covariant_under_body_coordinate_changes() {
    let array = |v: Vec3| [v.x, v.y, v.z];
    let vector = |v: [f64; 3]| Vec3::new(v[0], v[1], v[2]);
    let basis = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.6, 0.2))).unwrap();
    let principal = RigidBody::new(2.0, [2.0, 3.0, 4.0]).unwrap();
    let tensor = basis
        .matrix()
        .compose(principal.inertia_tensor().matrix())
        .compose(basis.matrix().transpose());
    let full = RigidBody::from_inertia_tensor(2.0, tensor).unwrap();
    assert!(
        tensor.m[1].abs() > 0.05,
        "fixture must exercise products of inertia"
    );
    let mut principal_state = SixDofState {
        pos_m: [2.0, -1.0, 4.0],
        vel_mps: [5.0, 0.2, -0.7],
        omega_body: [0.8, -0.4, 1.2],
        ..rest()
    };
    let q = basis.as_quat();
    principal_state.quat = [q.w, q.x, q.y, q.z];
    let full_state = SixDofState {
        quat: rest().quat,
        omega_body: array(basis.rotate(vector(principal_state.omega_body)).unwrap()),
        ..principal_state
    };
    let torque = Vec3::new(-0.2, 0.3, 0.1);
    let arm = Vec3::new(0.7, -0.3, 0.9);
    let observe = |body: &RigidBody, state: &SixDofState, torque: Vec3, arm: Vec3| {
        body.prepare_point_motion(
            state,
            Loads {
                force_n: [0.4, 0.1, -0.2],
                moment_nm: array(torque),
            },
        )
        .unwrap()
        .point_motion(vec_to_mbd(array(arm)))
        .unwrap()
    };
    let a = observe(&principal, &principal_state, torque, arm);
    let b = observe(
        &full,
        &full_state,
        basis.rotate(torque).unwrap(),
        basis.rotate(arm).unwrap(),
    );
    for (x, y) in [
        (a.kinematics.point_world, b.kinematics.point_world),
        (
            a.kinematics.point_velocity_world,
            b.kinematics.point_velocity_world,
        ),
        (a.point_acceleration_world, b.point_acceleration_world),
        (a.angular_acceleration_world, b.angular_acceleration_world),
    ] {
        close(mbd_to_vec(x), mbd_to_vec(y), 2e-12);
    }
}

#[test]
fn native_point_acceleration_matches_the_spine_trajectory_derivative() {
    let body = RigidBody::new(2.0, [2.0, 3.0, 4.0]).unwrap();
    let start = SixDofState {
        pos_m: [1.0, 2.0, 3.0],
        vel_mps: [4.0, 5.0, 6.0],
        omega_body: [0.4, -0.3, 0.7],
        ..rest()
    };
    let loads = Loads {
        force_n: [2.0, -4.0, 6.0],
        moment_nm: [0.2, -0.1, 0.4],
    };
    let arm = MbdVec3::new(0.8, -0.2, 0.6);
    let original = body
        .prepare_point_motion(&start, loads)
        .unwrap()
        .point_motion(arm)
        .unwrap();
    let mut errors = Vec::new();
    for dt in [0.004, 0.002, 0.001] {
        let mut calls = 0;
        let end = step(&body, &start, 0.0, dt, |_, _| {
            calls += 1;
            loads
        })
        .unwrap();
        assert_eq!(calls, 2);
        let endpoint = body
            .prepare_point_motion(&end, loads)
            .unwrap()
            .point_motion(arm)
            .unwrap();
        let difference = endpoint
            .kinematics
            .point_velocity_world
            .sub(original.kinematics.point_velocity_world)
            .scale(1.0 / dt);
        errors.push(
            difference
                .sub(original.point_acceleration_world)
                .norm_squared()
                .sqrt(),
        );
    }
    assert!(errors[2] < 0.003, "finite-difference errors {errors:?}");
    for pair in errors.windows(2) {
        assert!(
            (1.9..2.1).contains(&(pair[0] / pair[1])),
            "errors {errors:?}"
        );
    }
}

#[test]
fn physical_specific_force_passes_through_the_original_delay_and_filter() {
    let body = RigidBody::new(2.0, [2.0; 3]).unwrap();
    let motion = body
        .prepare_point_motion(
            &rest(),
            Loads {
                force_n: [0.0, 0.0, 6.0],
                moment_nm: [0.0; 3],
            },
        )
        .unwrap();
    let mut model = perception_v1(19);
    model.cues = [CueSpec {
        delay_ticks: 2,
        filter_tau_s: 0.04,
        gain: 2.0,
        remnant_sigma: 0.0,
    }; 6];
    let mut state = model.init().unwrap();
    for _ in 0..2 {
        let output = model
            .step_with_rigid_motion(&mut state, &motion, cue_input())
            .unwrap();
        assert_eq!(output.values[2], 0.0, "cue cannot arrive before its delay");
    }
    let output = model
        .step_with_rigid_motion(&mut state, &motion, cue_input())
        .unwrap();
    let expected = 6.0 * (1.0 - (-1.0 / (PERCEPTION_HZ * 0.04)).exp());
    assert!((output.values[2] - expected).abs() < 2e-14, "{output:?}");

    // Existing checkpoint bytes retain the whole observer history; no new
    // hidden acceleration cache or different replay owner is introduced.
    let bytes = state.to_bytes();
    let (mut restored, used) = fs_flyer::perception::PerceptionState::from_bytes(&bytes).unwrap();
    assert_eq!(used, bytes.len());
    for _ in 0..12 {
        assert_eq!(
            model
                .step_with_rigid_motion(&mut state, &motion, cue_input())
                .unwrap(),
            model
                .step_with_rigid_motion(&mut restored, &motion, cue_input())
                .unwrap(),
        );
    }
}

#[test]
fn off_center_rotating_pilot_feels_heave_while_the_center_does_not() {
    let body = RigidBody::new(2.0, [2.0; 3]).unwrap();
    let motion = body
        .prepare_point_motion(
            &SixDofState {
                omega_body: [0.0, 2.0, 0.0],
                ..rest()
            },
            Loads {
                force_n: [0.0; 3],
                moment_nm: [0.0; 3],
            },
        )
        .unwrap();
    let mut model = perception_v1(7);
    model.cues = [CueSpec {
        delay_ticks: 0,
        filter_tau_s: 0.04,
        gain: 1.0,
        remnant_sigma: 0.0,
    }; 6];
    let mut at_center = model.init().unwrap();
    let mut off_center = model.init().unwrap();
    let center = model
        .step_with_rigid_motion(&mut at_center, &motion, cue_input())
        .unwrap();
    let offset = model
        .step_with_rigid_motion(
            &mut off_center,
            &motion,
            RigidPointCueInput {
                point_body_m: [0.0, 0.0, 1.0],
                ..cue_input()
            },
        )
        .unwrap();
    assert_eq!(center.values[2], 0.0);
    let filter = 1.0 - (-1.0 / (PERCEPTION_HZ * 0.04)).exp();
    assert!((offset.values[2] + 4.0 * filter).abs() < 2e-14);
    assert!((offset.values[1] - 2.0 * filter).abs() < 2e-14);
}

#[test]
fn specific_force_is_rotated_into_body_axes_before_the_heave_filter() {
    let body = RigidBody::new(2.0, [2.0; 3]).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(
        0.0,
        core::f64::consts::FRAC_PI_2,
        0.0,
    )))
    .unwrap();
    let q = rotation.as_quat();
    let motion = body
        .prepare_point_motion(
            &SixDofState {
                quat: [q.w, q.x, q.y, q.z],
                ..rest()
            },
            Loads {
                force_n: [6.0, 0.0, 2.0 * G],
                moment_nm: [0.0; 3],
            },
        )
        .unwrap();
    let mut model = perception_v1(11);
    model.cues = [CueSpec {
        delay_ticks: 0,
        filter_tau_s: 0.04,
        gain: 1.0,
        remnant_sigma: 0.0,
    }; 6];
    let mut state = model.init().unwrap();
    let output = model
        .step_with_rigid_motion(
            &mut state,
            &motion,
            RigidPointCueInput {
                heave_reference_mps2: 0.0,
                ..cue_input()
            },
        )
        .unwrap();
    // Body z points along world x. The world's vertical force is weight,
    // so only the horizontal 3 m/s² specific force reaches the heave cue.
    let expected = 3.0 * (1.0 - (-1.0 / (PERCEPTION_HZ * 0.04)).exp());
    assert!((output.values[2] - expected).abs() < 2e-14);
}

#[test]
fn invalid_observations_refuse_before_mutating_perception() {
    let body = RigidBody::new(2.0, [2.0; 3]).unwrap();
    let loads = Loads {
        force_n: [0.0; 3],
        moment_nm: [0.0; 3],
    };
    let motion = body.prepare_point_motion(&rest(), loads).unwrap();
    for quat in [
        [0.0; 4],
        [2.0, 0.0, 0.0, 0.0],
        [-2.0, 0.0, 0.0, 0.0],
        [1.0 + 1e-10, 0.0, 0.0, 0.0],
        [f64::NAN, 0.0, 0.0, 0.0],
        [0.0, f64::INFINITY, 0.0, 0.0],
    ] {
        let state = SixDofState { quat, ..rest() };
        let before = tick_digest(0, &state);
        assert_eq!(
            body.prepare_point_motion(&state, loads).unwrap_err().code,
            "orientation-outside-group",
        );
        assert_eq!(before, tick_digest(0, &state));
    }
    for bad in [
        Loads {
            force_n: [f64::NAN, 0.0, 0.0],
            ..loads
        },
        Loads {
            moment_nm: [0.0, f64::INFINITY, 0.0],
            ..loads
        },
    ] {
        assert!(body.prepare_point_motion(&rest(), bad).is_err());
    }
    let model = perception_v1(13);
    let mut state = model.init().unwrap();
    let before = state.clone();
    for input in [
        RigidPointCueInput {
            point_body_m: [f64::NAN, 0.0, 0.0],
            ..cue_input()
        },
        RigidPointCueInput {
            gravity_world_mps2: [0.0, f64::INFINITY, 0.0],
            ..cue_input()
        },
        RigidPointCueInput {
            heave_reference_mps2: f64::NAN,
            ..cue_input()
        },
    ] {
        assert!(
            model
                .step_with_rigid_motion(&mut state, &motion, input)
                .is_err()
        );
        assert_eq!(state, before);
    }
}
