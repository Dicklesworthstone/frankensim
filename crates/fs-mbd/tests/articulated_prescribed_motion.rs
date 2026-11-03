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
    vector(
        step.external_linear_impulse_world,
        Vec3::new(0.2, 0.0, 0.0),
        2e-15,
    );
    vector(
        step.support_linear_impulse_world,
        Vec3::new(0.6, 0.0, 0.0),
        2e-15,
    );
    vector(
        step.linear_momentum_change_world,
        Vec3::new(0.8, 0.0, 0.0),
        2e-15,
    );
    vector(step.linear_momentum_residual_world, Vec3::default(), 3e-15);
    vector(step.angular_momentum_residual_world, Vec3::default(), 0.0);
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
        vector(step.external_linear_impulse_world, Vec3::default(), 0.0);
        vector(
            step.support_linear_impulse_world,
            Vec3::new(effort * dt, 0.0, 0.0),
            2e-14,
        );
        vector(
            step.linear_momentum_change_world,
            Vec3::new(effort * dt, 0.0, 0.0),
            2e-14,
        );
        vector(step.linear_momentum_residual_world, Vec3::default(), 2e-14);
        vector(step.angular_momentum_residual_world, Vec3::default(), 0.0);
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
        let dt = 0.05;
        let midpoint = rotation(omega.scale(0.5 * dt));
        let endpoint = rotation(omega.scale(dt));
        let p0 = omega.cross(arm + center).scale(2.0);
        let dp = endpoint.rotate(p0).unwrap() - p0;
        let l0 = inertia.inertia_com().apply(omega) + (arm + center).cross(p0);
        let dl = endpoint.rotate(l0).unwrap() - l0 + origin.cross(dp);
        let applied_force = midpoint.rotate(expected_force).unwrap().scale(dt);
        let applied_moment =
            midpoint.rotate(expected_torque).unwrap().scale(dt) + origin.cross(applied_force);
        vector(step.external_linear_impulse_world, Vec3::default(), 0.0);
        vector(step.external_angular_impulse_world, Vec3::default(), 0.0);
        vector(step.support_linear_impulse_world, applied_force, 2e-13);
        vector(step.support_angular_impulse_world, applied_moment, 2e-10);
        vector(step.linear_momentum_change_world, dp, 2e-13);
        vector(step.angular_momentum_change_world, dl, 2e-10);
        vector(
            step.linear_momentum_residual_world,
            dp - applied_force,
            2e-13,
        );
        vector(
            step.angular_momentum_residual_world,
            dl - applied_moment,
            2e-10,
        );
        // Finite midpoint quadrature is not the exact integral of a rotating
        // force. A nonzero defect must remain visible, never supplied as heat.
        assert!(step.linear_momentum_residual_world.norm() > 1e-5);
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
        vector(
            actual.external_linear_impulse_world,
            r.rotate(reference.external_linear_impulse_world).unwrap(),
            3e-13,
        );
        vector(
            actual.support_linear_impulse_world,
            r.rotate(reference.support_linear_impulse_world).unwrap(),
            3e-13,
        );
        vector(
            actual.linear_momentum_change_world,
            r.rotate(reference.linear_momentum_change_world).unwrap(),
            3e-13,
        );
        if translation.norm() < 10.0 {
            vector(
                actual.external_angular_impulse_world,
                r.rotate(reference.external_angular_impulse_world).unwrap()
                    + translation.cross(actual.external_linear_impulse_world),
                3e-12,
            );
            vector(
                actual.support_angular_impulse_world,
                r.rotate(reference.support_angular_impulse_world).unwrap()
                    + translation.cross(actual.support_linear_impulse_world),
                3e-12,
            );
            vector(
                actual.angular_momentum_change_world,
                r.rotate(reference.angular_momentum_change_world).unwrap()
                    + translation.cross(actual.linear_momentum_change_world),
                3e-12,
            );
        }
    }
}

#[test]
fn stationary_support_balances_original_external_and_gravity_impulses() {
    let center = Vec3::new(0.25, -0.5, 0.125);
    let arm = Vec3::new(0.5, 0.25, -0.125);
    let model = ArticulatedModel::new(vec![Link::new(
        "fixed_offset",
        None,
        Se3::from_parts(So3::identity(), arm).unwrap(),
        JointModel::FIXED,
        SpatialInertia::new(2.0, center, Mat3::identity()).unwrap(),
    )])
    .unwrap();
    let frame = Se3::from_parts(
        rotation(Vec3::new(0.2, -0.3, 0.1)),
        Vec3::new(4.0, -3.0, 2.0),
    )
    .unwrap();
    let initial = PrescribedBaseMotionState {
        time_s: 0.0,
        base: BaseState::stationary(frame),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let body = Wrench::new(Vec3::new(0.25, -0.125, 0.5), Vec3::new(1.0, 2.0, -1.0));
    let gravity = Vec3::new(0.0, 0.0, -9.5);
    let dt = 0.125;
    let step = step_prescribed_base(
        &model,
        &initial,
        PrescribedBasePath::stationary(frame),
        FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &[body],
            gravity_world: gravity,
        },
        dt,
        || false,
    )
    .unwrap();
    let force = frame.rotation().rotate(body.force).unwrap();
    let couple = frame.rotation().rotate(body.torque).unwrap();
    let expected_linear = (force + gravity.scale(2.0)).scale(dt);
    let expected_angular = (couple
        + frame.transform_point(arm).unwrap().cross(force)
        + frame
            .transform_point(arm + center)
            .unwrap()
            .cross(gravity.scale(2.0)))
    .scale(dt);
    vector(step.external_linear_impulse_world, expected_linear, 2e-13);
    vector(step.external_angular_impulse_world, expected_angular, 2e-13);
    vector(
        step.support_linear_impulse_world,
        expected_linear.scale(-1.0),
        2e-13,
    );
    vector(
        step.support_angular_impulse_world,
        expected_angular.scale(-1.0),
        2e-13,
    );
    vector(step.linear_momentum_change_world, Vec3::default(), 0.0);
    vector(step.angular_momentum_change_world, Vec3::default(), 0.0);
    vector(step.linear_momentum_residual_world, Vec3::default(), 2e-13);
    vector(step.angular_momentum_residual_world, Vec3::default(), 2e-13);
    assert_eq!(step.support_work, 0.0);
}

fn three_branches(moving: bool) -> ArticulatedModel {
    let inertia = SpatialInertia::new(1.0, Vec3::default(), Mat3::identity()).unwrap();
    let mut links = vec![Link::new(
        "root",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )];
    for i in 0..3 {
        let joint = if moving {
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap()
        } else {
            JointModel::FIXED
        };
        links.push(Link::new(
            format!("child{i}"),
            Some(0),
            Se3::identity(),
            joint,
            inertia,
        ));
    }
    ArticulatedModel::new(links).unwrap()
}

#[test]
fn actual_momentum_change_survives_large_translation_and_rounded_endpoint_totals() {
    let model = three_branches(true);
    let large = (1u64 << 52) as f64;
    let origin = Vec3::new(0.0, (1u64 << 20) as f64, -((1u64 << 19) as f64));
    let frame = Se3::from_parts(So3::identity(), origin).unwrap();
    let initial = PrescribedBaseMotionState {
        time_s: 0.0,
        base: BaseState::stationary(frame),
        joint_position: vec![0.0; 3],
        joint_velocity: vec![large, 0.0, large],
    };
    let step = step_prescribed_base(
        &model,
        &initial,
        PrescribedBasePath::stationary(frame),
        FreeFloatingMotionLoads {
            generalized_force: &[6.0, 1.0, -6.0],
            external_body_wrench: &[Wrench::default(); 4],
            gravity_world: Vec3::default(),
        },
        1.0,
        || false,
    )
    .unwrap();
    // Each unit-mass velocity increment is represented. The total's one-unit
    // change is nevertheless below the resolution of the absolute 2^53 total.
    assert_eq!(
        step.momentum_before.linear_momentum_world,
        step.momentum_after.linear_momentum_world
    );
    let change = Vec3::new(1.0, 0.0, 0.0);
    vector(step.linear_momentum_change_world, change, 0.0);
    vector(
        step.angular_momentum_change_world,
        origin.cross(change),
        0.0,
    );
    vector(step.external_linear_impulse_world, Vec3::default(), 0.0);
    vector(step.external_angular_impulse_world, Vec3::default(), 0.0);
    // Internal actuation reaches the stationary support through native RNEA;
    // it is not counted again as a prescribed external body impulse.
    vector(step.support_linear_impulse_world, change, 0.0);
    vector(
        step.support_angular_impulse_world,
        origin.cross(change),
        0.0,
    );
    vector(step.linear_momentum_residual_world, Vec3::default(), 0.0);
    vector(step.angular_momentum_residual_world, Vec3::default(), 0.0);
    assert_eq!(step.kinetic_energy_change, 36.5);
    assert_eq!(step.actuator_work, 36.5);
    assert_eq!(step.support_work, 0.0);
}

#[test]
fn imposed_motion_keeps_half_ulp_kinetic_increment_separate_from_midpoint_source_work() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        SpatialInertia::new(1.0, Vec3::default(), Mat3::identity()).unwrap(),
    )])
    .unwrap();
    let large = (1u64 << 55) as f64;
    let velocity = Vec3::new(large, 0.0, large);
    let acceleration = Vec3::new(8.0, 1.0, -8.0);
    let at = |t: f64| {
        BaseState::prescribed(
            Se3::from_parts(
                So3::identity(),
                velocity.scale(t) + acceleration.scale(0.5 * t * t),
            )
            .unwrap(),
            Twist::new(Vec3::default(), velocity + acceleration.scale(t)),
            Twist::new(Vec3::default(), acceleration),
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
            midpoint: at(0.5),
            endpoint: at(1.0),
        },
        FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &[Wrench::default()],
            gravity_world: Vec3::default(),
        },
        1.0,
        || false,
    )
    .unwrap();
    // The endpoint velocity increments are all representable and v dot delta-v
    // cancels. Its quadratic increment is (8² + 1² + 8²)/2 = 64.5 J.
    assert_eq!(step.kinetic_energy_change, 64.5);
    vector(step.linear_momentum_change_world, acceleration, 0.0);
    vector(step.support_linear_impulse_world, acceleration, 0.0);
    vector(step.linear_momentum_residual_world, Vec3::default(), 0.0);
    // This supplied midpoint twist rounds the first half-ULP away. Keep that
    // actual source quadrature and report the resulting defect independently.
    assert_eq!(step.support_work, 32.5);
    assert_eq!(step.kinetic_work_residual, 32.0);
}

#[test]
fn source_impulse_retains_small_force_and_couple_between_opposite_large_loads() {
    let model = three_branches(false);
    let large = (1u64 << 55) as f64;
    let origin = Vec3::new((1u64 << 40) as f64, 0.0, 0.0);
    let frame = Se3::from_parts(So3::identity(), origin).unwrap();
    let initial = PrescribedBaseMotionState {
        time_s: 0.0,
        base: BaseState::stationary(frame),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let external = [
        Wrench::default(),
        Wrench::new(Vec3::default(), Vec3::new(0.0, large, 0.0)),
        Wrench::new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 1.0, 0.0)),
        Wrench::new(Vec3::default(), Vec3::new(0.0, -large, 0.0)),
    ];
    let dt = 0.125;
    let step = step_prescribed_base(
        &model,
        &initial,
        PrescribedBasePath::stationary(frame),
        FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &external,
            gravity_world: Vec3::default(),
        },
        dt,
        || false,
    )
    .unwrap();
    let applied = Vec3::new(0.0, dt, 0.0);
    let moment = origin.cross(applied) + Vec3::new(0.0, 0.0, dt);
    vector(step.external_linear_impulse_world, applied, 0.0);
    vector(step.external_angular_impulse_world, moment, 0.0);
    vector(step.linear_momentum_change_world, Vec3::default(), 0.0);
    vector(step.angular_momentum_change_world, Vec3::default(), 0.0);
    // The source is independent of whatever rounding the native reaction pass
    // retains. Its visible balance defect is not used to replace that source.
    vector(
        step.linear_momentum_residual_world,
        (applied + step.support_linear_impulse_world).scale(-1.0),
        0.0,
    );
    vector(
        step.angular_momentum_residual_world,
        (moment + step.support_angular_impulse_world).scale(-1.0),
        0.0,
    );
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
