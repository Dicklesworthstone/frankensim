//! G1/G3/G4 integration checks for continuing prescribed-base mechanics.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath, step_prescribed_base,
    step_prescribed_base_with_effort,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointLimits, JointModel, Link, SpatialInertia,
    apply_prescribed_base_impulses,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() <= tolerance, "{a:?} != {b:?}");
}

fn rotation(vector: Vec3) -> So3 {
    So3::exp(So3Tangent::new(vector)).unwrap()
}

fn scalar_root(limits: Option<JointLimits>) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "slider",
        None,
        Se3::identity(),
        JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), limits).unwrap(),
        SpatialInertia::new(
            2.0,
            Vec3::default(),
            Mat3 {
                m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            },
        )
        .unwrap(),
    )])
    .unwrap()
}

fn state(q: f64, velocity: f64) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 0.0,
        base: BaseState::stationary(Se3::identity()),
        joint_position: vec![q],
        joint_velocity: vec![velocity],
    }
}

#[test]
fn scalar_root_support_transmits_actuation_and_external_loads_once() {
    let model = scalar_root(None);
    let initial = state(0.2, 0.4);
    let snapshot = initial.clone();
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[6.0],
        external_body_wrench: &[Wrench::new(Vec3::default(), Vec3::new(2.0, 0.0, 0.0))],
        gravity_world: Vec3::default(),
    };
    let step = step_prescribed_base(
        &model,
        &initial,
        PrescribedBasePath::stationary(Se3::identity()),
        loads,
        0.1,
        || false,
    )
    .unwrap();
    // m=2, total force=8: a=4; the separate root actuator transmits 6 N.
    close(step.state_after.joint_position[0], 0.26, 2e-15);
    close(step.state_after.joint_velocity[0], 0.8, 2e-15);
    vector(
        step.support_wrench_body.force,
        Vec3::new(6.0, 0.0, 0.0),
        2e-14,
    );
    close(step.actuator_work, 0.36, 2e-14);
    close(step.external_work, 0.12, 2e-14);
    close(step.support_work, 0.0, 0.0);
    close(step.kinetic_energy_change, 0.48, 3e-14);
    close(step.kinetic_work_residual, 0.0, 3e-14);
    assert_eq!(initial, snapshot);
}

#[test]
fn moving_accelerating_support_retains_relative_and_support_work() {
    let model = scalar_root(None);
    for (acceleration, effort) in [(1.0, 6.0), (4.0, 0.0)] {
        let base_at = |t: f64| {
            BaseState::prescribed(
                Se3::from_parts(
                    So3::identity(),
                    Vec3::new(2.0 * t + 0.5 * acceleration * t * t, 0.0, 0.0),
                )
                .unwrap(),
                Twist::new(Vec3::default(), Vec3::new(2.0 + acceleration * t, 0.0, 0.0)),
                Twist::new(Vec3::default(), Vec3::new(acceleration, 0.0, 0.0)),
            )
        };
        let mut initial = state(0.1, 0.4);
        initial.base = base_at(0.0);
        let dt = 0.2;
        let relative_acceleration = effort / 2.0 - acceleration;
        let step = step_prescribed_base(
            &model,
            &initial,
            PrescribedBasePath {
                midpoint: base_at(dt / 2.0),
                endpoint: base_at(dt),
            },
            FreeFloatingMotionLoads {
                generalized_force: &[effort],
                external_body_wrench: &[Wrench::default()],
                gravity_world: Vec3::default(),
            },
            dt,
            || false,
        )
        .unwrap();
        let expected_rate = 0.4 + relative_acceleration * dt;
        close(step.state_after.joint_velocity[0], expected_rate, 3e-15);
        close(
            step.state_after.joint_position[0],
            0.1 + 0.4 * dt + 0.5 * relative_acceleration * dt * dt,
            3e-15,
        );
        close(step.support_wrench_body.force.x, effort, 3e-14);
        close(
            step.support_work,
            effort * (2.0 + acceleration * dt / 2.0) * dt,
            4e-14,
        );
        close(
            step.actuator_work,
            effort * (0.4 + relative_acceleration * dt / 2.0) * dt,
            4e-14,
        );
        let velocity_after = 2.0 + acceleration * dt + expected_rate;
        close(
            step.kinetic_energy_change,
            velocity_after.powi(2) - 2.4_f64.powi(2),
            5e-14,
        );
        close(step.kinetic_work_residual, 0.0, 5e-14);
        assert_eq!(step.state_after.base, base_at(dt));
    }
}

#[test]
fn fixed_root_full_tensor_and_offset_keep_support_force_and_couple() {
    let center = Vec3::new(0.2, 0.3, 0.1);
    let arm = Vec3::new(0.4, -0.1, 0.2);
    let inertia = SpatialInertia::new(
        2.0,
        center,
        Mat3 {
            m: [2.0, 0.2, 0.3, 0.2, 3.0, 0.1, 0.3, 0.1, 4.0],
        },
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "offset_fixed_body",
        None,
        Se3::from_parts(So3::identity(), arm).unwrap(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let omega = Vec3::new(0.0, 0.0, 1.3);
    let expected_force = omega.cross(omega.cross(arm + center)).scale(2.0);
    let expected_torque =
        (arm + center).cross(expected_force) + omega.cross(inertia.inertia_com().apply(omega));
    for origin in [Vec3::default(), Vec3::new(1e3, -2e3, 3e3)] {
        let at = |t: f64| {
            BaseState::prescribed(
                Se3::from_parts(rotation(omega.scale(t)), origin).unwrap(),
                Twist::new(omega, Vec3::default()),
                Twist::zero(),
            )
        };
        let initial = PrescribedBaseMotionState {
            time_s: 0.0,
            base: at(0.0),
            joint_position: vec![],
            joint_velocity: vec![],
        };
        let step = step_prescribed_base(
            &model,
            &initial,
            PrescribedBasePath {
                midpoint: at(0.025),
                endpoint: at(0.05),
            },
            FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &[Wrench::default()],
                gravity_world: Vec3::default(),
            },
            0.05,
            || false,
        )
        .unwrap();
        vector(step.support_wrench_body.force, expected_force, 8e-14);
        vector(step.support_wrench_body.torque, expected_torque, 8e-14);
        close(step.support_work, 0.0, 2e-14);
        close(step.kinetic_energy_change, 0.0, 1e-13);
        assert!(step.state_after.joint_position.is_empty());
    }
}

#[test]
fn native_impulse_continues_without_changing_configuration_or_clock_at_the_jump() {
    let model = scalar_root(None);
    let mut initial = state(0.3, 0.2);
    initial.time_s = 1.0;
    let impulse = apply_prescribed_base_impulses(
        &model,
        initial.base,
        &initial.joint_position,
        &initial.joint_velocity,
        &[0.8],
        &[Wrench::default()],
    )
    .unwrap();
    let after_jump = PrescribedBaseMotionState {
        base: impulse.base_after,
        joint_velocity: impulse.joint_velocity_after,
        ..initial.clone()
    };
    let step = step_prescribed_base(
        &model,
        &after_jump,
        PrescribedBasePath::stationary(Se3::identity()),
        FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default()],
            gravity_world: Vec3::default(),
        },
        0.25,
        || false,
    )
    .unwrap();
    close(after_jump.time_s, 1.0, 0.0);
    close(after_jump.joint_position[0], 0.3, 0.0);
    close(after_jump.joint_velocity[0], 0.6, 2e-15);
    close(step.state_after.time_s, 1.25, 0.0);
    close(step.state_after.joint_position[0], 0.45, 2e-15);
    close(step.kinetic_energy_change, 0.0, 2e-15);
    close(step.support_work, 0.0, 0.0);
}

#[test]
fn rigid_frame_changes_preserve_work_and_small_gravity_displacements() {
    let model = scalar_root(None);
    let base_state = state(0.2, 0.4);
    let gravity = Vec3::new(-1.0, -2.0, -3.0);
    let wrench = Wrench::new(Vec3::new(0.1, 0.2, 0.3), Vec3::new(1.0, 0.5, -0.2));
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[4.0],
        external_body_wrench: &[wrench],
        gravity_world: gravity,
    };
    let reference = step_prescribed_base(
        &model,
        &base_state,
        PrescribedBasePath::stationary(Se3::identity()),
        loads,
        0.1,
        || false,
    )
    .unwrap();
    let r = rotation(Vec3::new(0.2, -0.4, 0.3));
    // General huge translated rotations can be outside the original motor's
    // admission. Exercise frame rotation and huge pure translation separately.
    for (r, translation) in [
        (r, Vec3::new(2.0, -3.0, 0.7)),
        (So3::identity(), Vec3::new(2e15, -3e15, 1e15)),
    ] {
        let frame = Se3::from_parts(r, translation).unwrap();
        let mut transformed = base_state.clone();
        transformed.base = BaseState::stationary(frame);
        let actual = step_prescribed_base(
            &model,
            &transformed,
            PrescribedBasePath::stationary(frame),
            FreeFloatingMotionLoads {
                gravity_world: r.rotate(gravity).unwrap(),
                ..loads
            },
            0.1,
            || false,
        )
        .unwrap();
        close(
            actual.state_after.joint_position[0],
            reference.state_after.joint_position[0],
            3e-14,
        );
        close(
            actual.state_after.joint_velocity[0],
            reference.state_after.joint_velocity[0],
            3e-14,
        );
        vector(
            actual.support_wrench_body.force,
            reference.support_wrench_body.force,
            3e-13,
        );
        vector(
            actual.support_wrench_body.torque,
            reference.support_wrench_body.torque,
            3e-13,
        );
        close(actual.gravity_work, reference.gravity_work, 3e-13);
        close(actual.external_work, reference.external_work, 3e-13);
        close(actual.actuator_work, reference.actuator_work, 3e-13);
        close(
            actual.kinetic_energy_change,
            reference.kinetic_energy_change,
            3e-13,
        );
    }
}

#[test]
fn actual_stage_feedback_has_second_order_convergence_and_resets_held_efforts() {
    let model = scalar_root(None);
    let path = PrescribedBasePath::stationary(Se3::identity());
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default()],
        gravity_world: Vec3::default(),
    };
    let mut errors = Vec::new();
    for n in [40_u32, 80, 160] {
        let dt = 1.0 / f64::from(n);
        let mut current = state(1.0, 0.0);
        for _ in 0..n {
            let mut calls = 0;
            let next = step_prescribed_base_with_effort(
                &model,
                &current,
                path,
                loads,
                dt,
                |elapsed, stage, output| {
                    assert_eq!(&*output, &[0.0]);
                    close(stage.time_s, current.time_s + elapsed, 0.0);
                    close(elapsed, f64::from(calls) * dt / 2.0, 0.0);
                    calls += 1;
                    output[0] = -8.0 * stage.joint_position[0];
                    Ok(())
                },
                || false,
            )
            .unwrap();
            assert_eq!(calls, 2);
            current = next.state_after;
        }
        errors.push(
            (current.joint_position[0] - 2.0_f64.cos()).abs()
                + (current.joint_velocity[0] + 2.0 * 2.0_f64.sin()).abs(),
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 3.5 && pair[0] / pair[1] < 4.5,
            "{errors:?}"
        );
    }
}

#[test]
fn limits_invalid_base_data_and_last_cancellation_leave_borrowed_state_unchanged() {
    let limits = JointLimits::new(-1.0, 0.1, 3.0, 10.0).unwrap();
    let model = scalar_root(Some(limits));
    let initial = state(0.0, 0.8);
    let snapshot = initial.clone();
    let path = PrescribedBasePath::stationary(Se3::identity());
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default()],
        gravity_world: Vec3::default(),
    };
    assert!(matches!(
        step_prescribed_base(&model, &initial, path, loads, 0.2, || false),
        Err(PrescribedBaseMotionError::Dynamics(
            ArticulatedError::JointLimitViolation { .. }
        ))
    ));
    let mut bad_path = path;
    bad_path.endpoint.acceleration_body.linear.x = f64::NAN;
    let mut callbacks = 0;
    assert!(
        step_prescribed_base_with_effort(
            &model,
            &initial,
            bad_path,
            loads,
            0.01,
            |_, _, _| {
                callbacks += 1;
                Ok(())
            },
            || false,
        )
        .is_err()
    );
    assert_eq!(callbacks, 0);
    let mut polls = 0;
    step_prescribed_base(&model, &initial, path, loads, 0.01, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls] {
        let mut count = 0;
        assert_eq!(
            step_prescribed_base(&model, &initial, path, loads, 0.01, || {
                count += 1;
                count == stop
            }),
            Err(PrescribedBaseMotionError::Cancelled),
        );
    }
    let mut unresolvable = initial.clone();
    unresolvable.time_s = 1e20;
    assert!(matches!(
        step_prescribed_base(&model, &unresolvable, path, loads, 0.01, || false),
        Err(PrescribedBaseMotionError::Invalid(_))
    ));
    assert_eq!(initial, snapshot);
}
