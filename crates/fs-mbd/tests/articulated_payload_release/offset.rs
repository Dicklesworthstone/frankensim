use super::*;

fn pose(axes: Vec3, translation: Vec3) -> Se3 {
    Se3::from_parts(rotation(axes), translation).unwrap()
}

fn point_velocity_world(pose: Se3, twist: Twist, point: Vec3) -> (Vec3, Vec3) {
    let angular = pose.rotation().rotate(twist.angular).unwrap();
    let velocity =
        pose.rotation().rotate(twist.linear).unwrap() + angular.cross(point - pose.translation());
    (angular, velocity)
}

fn offset_tree() -> (ArticulatedModel, FreeFloatingMotionState) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            pose(Vec3::new(0.2, -0.1, 0.15), Vec3::new(0.17, -0.04, 0.11)),
            JointModel::FIXED,
            inertia(
                2.3,
                Vec3::new(0.12, -0.21, 0.08),
                Vec3::new(1.2, 1.6, 2.),
                Vec3::new(0.2, -0.15, 0.3),
            ),
        ),
        Link::new(
            "rotor",
            Some(0),
            pose(Vec3::new(0.1, 0.25, -0.15), Vec3::new(0.6, 0.08, -0.11)),
            JointModel::revolute(Vec3::new(1., 2., -1.), None).unwrap(),
            inertia(
                1.1,
                Vec3::new(-0.08, 0.17, 0.04),
                Vec3::new(0.5, 0.7, 0.9),
                Vec3::new(0.3, -0.15, 0.3),
            ),
        ),
        Link::new(
            "screw",
            Some(1),
            pose(Vec3::new(0.11, -0.12, 0.18), Vec3::new(-0.2, 0.17, 0.33)),
            JointModel::helical(Vec3::new(2., -1., 1.), 0.13, None).unwrap(),
            inertia(
                0.4,
                Vec3::new(0.04, -0.12, 0.06),
                Vec3::new(0.2, 0.25, 0.3),
                Vec3::new(0.5, -0.15, 0.3),
            ),
        ),
    ])
    .unwrap();
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            pose(Vec3::new(0.31, -0.27, 0.19), Vec3::new(2., -1., 0.6)),
            Twist::new(Vec3::new(0.4, -0.3, 0.2), Vec3::new(-0.2, 0.6, 0.35)),
        ),
        joint_position: vec![0.2, 0.15],
        joint_velocity: vec![0.4, 0.3],
    };
    (model, state)
}

fn checked_release(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    relative: Se3,
    payload_inertia: SpatialInertia,
    principal: Vec3,
    principal_axes: Vec3,
    impulse: SeparationImpulse,
) -> PayloadRelease {
    let target = 2;
    let loaded = load(model, state, target, relative, payload_inertia, &[]);
    close(loaded.removed_relative_energy_j, 0.);
    let initial = attached(
        &loaded.model_after,
        &loaded.state_after,
        target,
        relative,
        payload_inertia,
    );
    let before_kin = kinematics(&loaded.model_after, &loaded.state_after);
    let released = release_payload(
        &loaded.model_after,
        &loaded.state_after,
        &loaded.attachment,
        Some(impulse),
        &[],
        controls(&[2e-10; 2], 20.),
        || false,
    )
    .unwrap();
    assert_eq!(released.model_after, *model);
    assert_eq!(
        released.state_after.base.world_from_base,
        state.base.world_from_base
    );
    assert_eq!(released.state_after.joint_position, state.joint_position);
    assert_eq!(released.payload_after.inertia, payload_inertia);
    assert_eq!(released.applied_impulse, Some(impulse));
    assert!(
        (released.state_after.base.twist_body.angular - loaded.state_after.base.twist_body.angular)
            .norm()
            > 1e-3
    );
    assert!(
        released
            .state_after
            .joint_velocity
            .iter()
            .zip(&loaded.state_after.joint_velocity)
            .any(|(new, old)| (new - old).abs() > 1e-3)
    );

    // Analytic COM inverse from the declared principal axes, independent of
    // the native six-dimensional spatial-inertia factorization used by release.
    let payload_from_link = relative.inverse().unwrap();
    let point = payload_from_link
        .transform_point(impulse.point_link_m)
        .unwrap();
    let force = payload_from_link
        .rotation()
        .rotate(impulse.impulse_on_payload_link_ns)
        .unwrap();
    let couple = payload_from_link
        .rotation()
        .rotate(impulse.couple_on_payload_link_nms)
        .unwrap();
    let torque_com = couple + (point - payload_inertia.center_of_mass()).cross(force);
    let axes = rotation(principal_axes).matrix();
    let principal_torque = axes.transpose().apply(torque_com);
    let delta_angular = axes.apply(Vec3::new(
        principal_torque.x / principal.x,
        principal_torque.y / principal.y,
        principal_torque.z / principal.z,
    ));
    let delta_linear = force.scale(1. / payload_inertia.mass())
        - delta_angular.cross(payload_inertia.center_of_mass());
    vector(
        released.payload_after.twist_body.angular - initial.twist_body.angular,
        delta_angular,
    );
    vector(
        released.payload_after.twist_body.linear - initial.twist_body.linear,
        delta_linear,
    );
    vector(released.impulse_on_payload_body.force, force);
    vector(
        released.impulse_on_payload_body.torque,
        point.cross(force) + couple,
    );
    vector(
        released.impulse_on_retained_link.force,
        impulse.impulse_on_payload_link_ns.scale(-1.),
    );
    vector(
        released.impulse_on_retained_link.torque,
        (impulse
            .point_link_m
            .cross(impulse.impulse_on_payload_link_ns)
            + impulse.couple_on_payload_link_nms)
            .scale(-1.),
    );

    // Compute the pair's work at one physical WORLD point. The point-force
    // contribution and free-couple contribution are evaluated independently;
    // their midpoint relative velocities remove a common translational boost.
    let after_kin = kinematics(&released.model_after, &released.state_after);
    let link_pose = before_kin.world_from_link[target];
    let world_point = link_pose.transform_point(impulse.point_link_m).unwrap();
    let (wp0, vp0) = point_velocity_world(initial.world_from_body, initial.twist_body, world_point);
    let (wp1, vp1) = point_velocity_world(
        released.payload_after.world_from_body,
        released.payload_after.twist_body,
        world_point,
    );
    let (wr0, vr0) = point_velocity_world(link_pose, before_kin.body_twist[target], world_point);
    let (wr1, vr1) = point_velocity_world(
        after_kin.world_from_link[target],
        after_kin.body_twist[target],
        world_point,
    );
    vector(vp0, vr0);
    vector(wp0, wr0);
    let force_world = link_pose
        .rotation()
        .rotate(impulse.impulse_on_payload_link_ns)
        .unwrap();
    let couple_world = link_pose
        .rotation()
        .rotate(impulse.couple_on_payload_link_nms)
        .unwrap();
    let work = force_world.dot((vp0 - vr0 + vp1 - vr1).scale(0.5))
        + couple_world.dot((wp0 - wr0 + wp1 - wr1).scale(0.5));
    close(released.supplied_work_j, work);
    assert!(work > 0.01);
    close(released.constraint_work_j, 0.);

    let before = physical(
        &loaded.model_after,
        &loaded.state_after,
        None,
        Vec3::default(),
    );
    let after = physical(
        &released.model_after,
        &released.state_after,
        Some(&released.payload_after),
        Vec3::default(),
    );
    vector(after.linear, before.linear);
    vector(after.angular, before.angular);
    let balance = &released.balance;
    vector(balance.momentum_before.linear_momentum_world, before.linear);
    vector(
        balance.momentum_before.angular_momentum_world,
        before.angular,
    );
    vector(balance.momentum_after.linear_momentum_world, after.linear);
    vector(balance.momentum_after.angular_momentum_world, after.angular);
    let before_common = physical(
        &loaded.model_after,
        &loaded.state_after,
        None,
        balance.accounting_boost_world,
    );
    let after_common = physical(
        &released.model_after,
        &released.state_after,
        Some(&released.payload_after),
        balance.accounting_boost_world,
    );
    close(balance.kinetic_before_j, before_common.kinetic);
    close(balance.kinetic_after_j, after_common.kinetic);
    close(after_common.kinetic - before_common.kinetic, work);
    close(balance.energy_residual_j, 0.);
    close(balance.reaction_work_j, 0.);
    assert!(balance.maximum_normalized_momentum_residual <= 1.);
    assert!(balance.multipliers.is_empty());
    for ((before, after), impulse) in balance
        .generalized_momentum_before
        .iter()
        .zip(&balance.generalized_momentum_after)
        .zip(&balance.generalized_constraint_impulse)
    {
        close(*after, *before);
        close(*impulse, 0.);
    }
    released
}

#[test]
fn offset_full_tensor_release_matches_com_recoil_and_shared_point_work_in_every_scene() {
    let (model, state) = offset_tree();
    let principal = Vec3::new(0.18, 0.24, 0.3);
    let principal_axes = Vec3::new(-0.2, 0.3, 0.15);
    let payload_inertia = inertia(0.7, Vec3::new(0.08, -0.06, 0.11), principal, principal_axes);
    let relative = pose(Vec3::new(0.25, -0.17, 0.31), Vec3::new(0.28, -0.19, 0.23));
    let impulse = SeparationImpulse {
        point_link_m: Vec3::new(-0.2, 0.3, 0.15),
        impulse_on_payload_link_ns: Vec3::new(0.35, -0.2, 0.1),
        couple_on_payload_link_nms: Vec3::new(0.04, -0.02, 0.03),
    };
    let reference = checked_release(
        &model,
        &state,
        relative,
        payload_inertia,
        principal,
        principal_axes,
        impulse,
    );
    let scene = pose(Vec3::new(-0.4, 0.2, 0.3), Vec3::new(4., -2.5, 1.7));
    for boost in [Vec3::default(), Vec3::new(3., -2., 4.)] {
        let mut transformed = state.clone();
        transformed.base.world_from_base = scene.compose(state.base.world_from_base).unwrap();
        let local_boost = transformed
            .base
            .world_from_base
            .rotation()
            .inverse()
            .rotate(boost)
            .unwrap();
        transformed.base.twist_body.linear = transformed.base.twist_body.linear + local_boost;
        let result = checked_release(
            &model,
            &transformed,
            relative,
            payload_inertia,
            principal,
            principal_axes,
            impulse,
        );
        vector(
            result.state_after.base.twist_body.angular,
            reference.state_after.base.twist_body.angular,
        );
        vector(
            result.state_after.base.twist_body.linear,
            reference.state_after.base.twist_body.linear + local_boost,
        );
        for (actual, expected) in result
            .state_after
            .joint_velocity
            .iter()
            .zip(&reference.state_after.joint_velocity)
        {
            close(*actual, *expected);
        }
        let expected_pose = scene
            .compose(reference.payload_after.world_from_body)
            .unwrap();
        vector(
            result.payload_after.world_from_body.translation(),
            expected_pose.translation(),
        );
        for (actual, expected) in result
            .payload_after
            .world_from_body
            .rotation()
            .matrix()
            .m
            .iter()
            .zip(expected_pose.rotation().matrix().m)
        {
            close(*actual, expected);
        }
        let payload_boost = expected_pose.rotation().inverse().rotate(boost).unwrap();
        vector(
            result.payload_after.twist_body.angular,
            reference.payload_after.twist_body.angular,
        );
        vector(
            result.payload_after.twist_body.linear,
            reference.payload_after.twist_body.linear + payload_boost,
        );
        vector(
            result.balance.accounting_boost_world,
            scene
                .rotation()
                .rotate(reference.balance.accounting_boost_world)
                .unwrap()
                + boost,
        );
        close(result.supplied_work_j, reference.supplied_work_j);
        close(
            result.balance.kinetic_before_j,
            reference.balance.kinetic_before_j,
        );
        close(
            result.balance.kinetic_after_j,
            reference.balance.kinetic_after_j,
        );
    }
}
