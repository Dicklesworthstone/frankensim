//! G1: native mass actions on independently known coupled joint inertia.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel,
    Link, SpatialInertia, apply_free_floating_impulses, forward_dynamics, free_floating_momentum,
    stationary_base_joint_momentum, stationary_base_joint_velocity,
};

fn model() -> ArticulatedModel {
    let limits = JointLimits::new(-10., 10., 20., 1e-6).unwrap();
    let joint = JointModel::revolute(Vec3::new(0., 0., 1.), Some(limits)).unwrap();
    let inertia = |i| {
        SpatialInertia::new(
            2.,
            Vec3::default(),
            Mat3 {
                m: [i, 0., 0., 0., i, 0., 0., 0., i],
            },
        )
        .unwrap()
    };
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(1.),
        ),
        Link::new("first", Some(0), Se3::identity(), joint, inertia(2.)),
        Link::new("second", Some(1), Se3::identity(), joint, inertia(3.)),
    ])
    .unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} != {b}");
}

#[test]
fn native_recursions_apply_the_coupled_mass_law_without_effort_limit_confusion() {
    let model = model();
    let q = [0.4, -0.2];
    let v = [3., -1.];
    // T = 1/2 I1 v1^2 + 1/2 I2 (v1+v2)^2.
    // M = [[5,3],[3,3]], p = [12,6].
    let p = stationary_base_joint_momentum(&model, &q, &v).unwrap();
    close(p[0], 12.);
    close(p[1], 6.);
    let recovered = stationary_base_joint_velocity(&model, &q, &p).unwrap();
    close(recovered[0], 3.);
    close(recovered[1], -1.);
    // The same numbers as physical force MUST still violate actuator limits.
    assert!(matches!(
        forward_dynamics(
            &model,
            BaseState::stationary(Se3::identity()),
            &q,
            &[0.; 2],
            &p,
            Vec3::default(),
            &[Wrench::default(); 3]
        ),
        Err(ArticulatedError::JointEffortLimitViolation { .. })
    ));
}

#[test]
fn mass_actions_retain_position_rate_dimension_and_finite_admission() {
    let model = model();
    assert!(stationary_base_joint_momentum(&model, &[0.; 2], &[21., 0.]).is_err());
    assert!(stationary_base_joint_velocity(&model, &[11., 0.], &[0.; 2]).is_err());
    assert!(stationary_base_joint_velocity(&model, &[0.; 2], &[1.]).is_err());
    assert!(stationary_base_joint_velocity(&model, &[0.; 2], &[f64::NAN, 0.]).is_err());
    assert!(stationary_base_joint_velocity(&model, &[0.; 2], &[1e6, 0.]).is_err());
    assert_eq!(
        stationary_base_joint_velocity(&model, &[0.; 2], &[0.; 2]).unwrap(),
        vec![0.; 2]
    );
}

fn close_vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 2e-12,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn free_body_offset_impulse_matches_independent_com_and_inertia_solution() {
    let center = Vec3::new(0.2, -0.1, 0.3);
    let inertia = SpatialInertia::new(
        2.0,
        center,
        Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 2.0],
        },
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "free body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.7))).unwrap();
    let position = Vec3::new(1.0, -2.0, 0.5);
    let base = FreeFloatingBaseState::stationary(Se3::from_parts(rotation, position).unwrap());
    let arm = Vec3::new(0.5, 0.4, -0.2);
    let impulse = Vec3::new(2.0, -1.0, 0.5);
    let torque = (arm - center).cross(impulse);
    let expected_omega = Vec3::new(torque.x, torque.y / 1.5, torque.z / 2.0);
    let expected_origin_velocity = impulse.scale(0.5) - expected_omega.cross(center);
    let result = apply_free_floating_impulses(
        &model,
        base,
        &[],
        &[],
        &[],
        &[Wrench::new(arm.cross(impulse), impulse)],
    )
    .unwrap();
    assert_eq!(result.base_after.world_from_base, base.world_from_base);
    close_vector(result.base_after.twist_body.angular, expected_omega);
    close_vector(
        result.base_after.twist_body.linear,
        expected_origin_velocity,
    );
    let expected_linear = rotation.rotate(impulse).unwrap();
    let expected_angular =
        rotation.rotate(arm.cross(impulse)).unwrap() + position.cross(expected_linear);
    close_vector(result.momentum_after.linear_momentum_world, expected_linear);
    close_vector(
        result.momentum_after.angular_momentum_world,
        expected_angular,
    );
    close_vector(result.external_linear_impulse_world, expected_linear);
    close_vector(result.external_angular_impulse_world, expected_angular);
    let expected_energy = 0.25 * impulse.dot(impulse) + 0.5 * torque.dot(expected_omega);
    close(result.momentum_after.kinetic_energy, expected_energy);
    close(result.impulse_work, expected_energy);
    close(result.kinetic_work_residual, 0.0);
}

#[test]
fn free_base_and_joints_share_impulse_according_to_full_coupled_mass() {
    let model = model();
    let q = [0.4, -0.2];
    let v = [3.0, -1.0];
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::new(0.0, 0.0, 0.4), Vec3::new(0.2, -0.1, 0.3)),
    );
    let mut external = [Wrench::default(); 3];
    external[2] = Wrench::new(Vec3::new(0.0, 0.0, 0.6), Vec3::default());
    let result =
        apply_free_floating_impulses(&model, base, &q, &v, &[0.8, -0.3], &external).unwrap();
    // Independent generalized angular mass for [base spin, qdot1, qdot2]:
    // [[6,5,3],[5,5,3],[3,3,3]]. RHS = [0.6,1.4,0.3].
    // Its solution is delta-v = [-0.8,1.35,-0.45].
    close(result.base_after.twist_body.angular.z, -0.4);
    close(result.joint_velocity_after[0], 4.35);
    close(result.joint_velocity_after[1], -1.45);
    close_vector(result.base_after.twist_body.linear, base.twist_body.linear);
    close_vector(
        result.momentum_after.linear_momentum_world,
        result.momentum_before.linear_momentum_world,
    );
    close_vector(
        result.momentum_after.angular_momentum_world
            - result.momentum_before.angular_momentum_world,
        Vec3::new(0.0, 0.0, 0.6),
    );
    close(result.kinetic_work_residual, 0.0);
    // The momentum query must report the actual returned state's same totals.
    let actual =
        free_floating_momentum(&model, result.base_after, &q, &result.joint_velocity_after)
            .unwrap();
    assert_eq!(actual, result.momentum_after);
}

#[test]
fn internal_joint_impulses_recoil_base_without_creating_spatial_momentum() {
    let model = model();
    let base = FreeFloatingBaseState::stationary(Se3::identity());
    let result = apply_free_floating_impulses(
        &model,
        base,
        &[0.4, -0.2],
        &[0.0; 2],
        &[0.8, -0.3],
        &[Wrench::default(); 3],
    )
    .unwrap();
    close(result.base_after.twist_body.angular.z, -0.8);
    close(result.joint_velocity_after[0], 1.35);
    close(result.joint_velocity_after[1], -0.65);
    close_vector(result.momentum_after.linear_momentum_world, Vec3::default());
    close_vector(
        result.momentum_after.angular_momentum_world,
        Vec3::default(),
    );
    close_vector(result.external_angular_impulse_world, Vec3::default());
    // Absolute link spin increments are [-0.8,0.55,-0.1], not joint rates.
    let exact_energy = 0.5 * (0.8 * 0.8 + 2.0 * 0.55 * 0.55 + 3.0 * 0.1 * 0.1);
    close(result.momentum_after.kinetic_energy, exact_energy);
    close(result.impulse_work, exact_energy);
    close(result.kinetic_work_residual, 0.0);
}

#[test]
fn articulated_impulse_and_momentum_are_covariant_under_world_frame_changes() {
    let model = model();
    let pose = Se3::from_parts(So3::identity(), Vec3::new(1.0, -2.0, 0.5)).unwrap();
    let base = FreeFloatingBaseState::new(
        pose,
        Twist::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.2, 0.4, -0.1)),
    );
    let q = [0.4, -0.2];
    let v = [0.3, -0.2];
    let mut external = [Wrench::default(); 3];
    external[1] = Wrench::new(Vec3::new(0.2, -0.3, 0.1), Vec3::new(0.3, -0.1, 0.4));
    let original =
        apply_free_floating_impulses(&model, base, &q, &v, &[0.1, -0.05], &external).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(-0.3, 0.6, 0.2))).unwrap();
    let translation = Vec3::new(-0.3, 0.7, 0.2);
    let world_change = Se3::from_parts(rotation, translation).unwrap();
    let transformed_base =
        FreeFloatingBaseState::new(world_change.compose(pose).unwrap(), base.twist_body);
    let transformed =
        apply_free_floating_impulses(&model, transformed_base, &q, &v, &[0.1, -0.05], &external)
            .unwrap();
    assert_eq!(
        transformed.base_after.twist_body,
        original.base_after.twist_body
    );
    assert_eq!(
        transformed.joint_velocity_after,
        original.joint_velocity_after
    );
    close(transformed.impulse_work, original.impulse_work);
    close(
        transformed.momentum_after.kinetic_energy,
        original.momentum_after.kinetic_energy,
    );
    let rotated_linear = rotation
        .rotate(original.momentum_after.linear_momentum_world)
        .unwrap();
    close_vector(
        transformed.momentum_after.linear_momentum_world,
        rotated_linear,
    );
    close_vector(
        transformed.momentum_after.angular_momentum_world,
        rotation
            .rotate(original.momentum_after.angular_momentum_world)
            .unwrap()
            + translation.cross(rotated_linear),
    );
    close_vector(
        original.momentum_after.linear_momentum_world
            - original.momentum_before.linear_momentum_world,
        original.external_linear_impulse_world,
    );
    close_vector(
        original.momentum_after.angular_momentum_world
            - original.momentum_before.angular_momentum_world,
        original.external_angular_impulse_world,
    );
    close(original.kinetic_work_residual, 0.0);
}

#[test]
fn articulated_impulse_admission_is_atomic_and_does_not_clip_velocity() {
    let model = model();
    let base = FreeFloatingBaseState::stationary(Se3::identity());
    let q = [0.4, -0.2];
    let v = [0.3, -0.2];
    let zero = [Wrench::default(); 3];
    let no_impulse = apply_free_floating_impulses(&model, base, &q, &v, &[0.0; 2], &zero).unwrap();
    assert_eq!(no_impulse.base_after, base);
    assert_eq!(no_impulse.joint_velocity_after, v);
    assert_eq!(no_impulse.momentum_before, no_impulse.momentum_after);
    close(no_impulse.impulse_work, 0.0);
    assert!(matches!(
        apply_free_floating_impulses(&model, base, &q, &v, &[1e6, 0.0], &zero),
        Err(ArticulatedError::JointVelocityLimitViolation { .. })
    ));
    assert!(apply_free_floating_impulses(&model, base, &q, &v, &[f64::NAN, 0.0], &zero).is_err());
    assert!(apply_free_floating_impulses(&model, base, &q, &v, &[], &zero).is_err());
    assert!(apply_free_floating_impulses(&model, base, &q, &v, &[0.0; 2], &[]).is_err());
    let invalid = Wrench::new(Vec3::default(), Vec3::new(f64::INFINITY, 0.0, 0.0));
    assert!(apply_free_floating_impulses(&model, base, &q, &v, &[0.0; 2], &[invalid; 3]).is_err());
    assert_eq!(base.twist_body, Twist::zero());
    assert_eq!(q, [0.4, -0.2]);
    assert_eq!(v, [0.3, -0.2]);
}
