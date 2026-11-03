//! G0/G1/G3 evidence for smooth rigid material-point acceleration.

use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn properties() -> MassProperties {
    MassProperties::new(2.5, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap()
}

fn state_with_velocity(
    properties: MassProperties,
    pose: Pose,
    velocity_world: Vec3,
    angular_velocity_body: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        pose,
        velocity_world.scale(properties.mass()),
        properties
            .angular_momentum_body_checked(angular_velocity_body)
            .unwrap(),
    )
    .unwrap()
}

fn assert_vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    for (actual_component, expected_component) in [
        (actual.x, expected.x),
        (actual.y, expected.y),
        (actual.z, expected.z),
    ] {
        assert!(
            (actual_component - expected_component).abs()
                <= tolerance * expected_component.abs().max(1.0),
            "actual {actual:?}, expected {expected:?}, tolerance {tolerance}"
        );
    }
}

#[test]
fn applied_force_and_gravity_accelerate_com_and_nonrotating_attachments() {
    let properties = properties();
    let pose = Pose::new(Vec3::new(1.0, -2.0, 3.0), UnitQuaternion::IDENTITY).unwrap();
    let state = state_with_velocity(properties, pose, Vec3::new(0.5, -0.2, 0.1), Vec3::ZERO);
    let integrator = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, 0.0, -9.81)).unwrap());
    let prepared = integrator
        .prepare_point_motion(
            state,
            properties,
            Wrench {
                force_world: Vec3::new(5.0, -2.5, 7.5),
                torque_body: Vec3::ZERO,
            },
        )
        .unwrap();

    for arm in [Vec3::ZERO, Vec3::new(0.4, -0.3, 0.2)] {
        let motion = prepared.point_motion(arm).unwrap();
        assert_eq!(
            motion.kinematics,
            state.point_kinematics(properties, arm).unwrap()
        );
        assert_vector_close(
            motion.center_of_mass_acceleration_world,
            Vec3::new(2.0, -1.0, -6.81),
            1e-14,
        );
        assert_vector_close(
            motion.point_acceleration_world,
            Vec3::new(2.0, -1.0, -6.81),
            1e-14,
        );
        assert_eq!(motion.angular_acceleration_body, Vec3::ZERO);
        assert_eq!(motion.angular_acceleration_world, Vec3::ZERO);
    }
    assert_eq!(prepared.pose(), pose);
    assert_eq!(prepared.angular_velocity_body(), Vec3::ZERO);
}

#[test]
fn spinning_attachment_includes_centripetal_and_torque_accelerations() {
    let properties = properties();
    let state = state_with_velocity(
        properties,
        Pose::identity(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 3.0),
    );
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let arm = Vec3::new(2.0, 0.0, 0.0);
    for (torque_z, expected_alpha_z, expected_acceleration_y) in [(0.0, 0.0, 0.0), (8.0, 2.0, 4.0)]
    {
        let motion = integrator
            .prepare_point_motion(
                state,
                properties,
                Wrench {
                    force_world: Vec3::ZERO,
                    torque_body: Vec3::new(0.0, 0.0, torque_z),
                },
            )
            .unwrap()
            .point_motion(arm)
            .unwrap();
        assert_vector_close(
            motion.kinematics.point_velocity_world,
            Vec3::new(0.0, 6.0, 0.0),
            1e-14,
        );
        assert_vector_close(
            motion.angular_acceleration_world,
            Vec3::new(0.0, 0.0, expected_alpha_z),
            1e-14,
        );
        assert_vector_close(
            motion.point_acceleration_world,
            Vec3::new(-18.0, expected_acceleration_y, 0.0),
            1e-14,
        );
    }
}

#[test]
fn torque_free_asymmetric_body_has_gyroscopic_point_acceleration() {
    let properties = properties();
    let state =
        RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::new(2.0, 6.0, 12.0)).unwrap();
    let motion = RigidBodyIntegrator::new(Gravity::ZERO)
        .prepare_point_motion(state, properties, Wrench::ZERO)
        .unwrap()
        .point_motion(Vec3::new(0.4, -0.3, 0.2))
        .unwrap();

    // I = diag(2, 3, 4), omega = (1, 2, 3), omega x L = (6, -6, 2).
    // The explicit expected acceleration includes both Euler-top precession
    // and centripetal acceleration; zero applied torque does not mean alpha=0.
    assert_vector_close(
        motion.angular_acceleration_body,
        Vec3::new(-3.0, 2.0, -0.5),
        1e-13,
    );
    assert_vector_close(
        motion.point_acceleration_world,
        Vec3::new(-4.95, 5.4, -1.5),
        1e-13,
    );
}

#[test]
fn full_inertia_preserves_motion_when_body_coordinates_are_reexpressed() {
    let properties = properties();
    let position = Vec3::new(1.0, -2.0, 0.5);
    let state = state_with_velocity(
        properties,
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.2, 0.4, -0.1),
        Vec3::new(1.2, -0.7, 0.9),
    );
    let wrench = Wrench {
        force_world: Vec3::new(1.0, 2.0, -0.7),
        torque_body: Vec3::new(0.5, -0.4, 0.2),
    };
    let integrator = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, 0.0, -3.0)).unwrap());
    let arm = Vec3::new(0.4, -0.3, 0.2);
    let expected = integrator
        .prepare_point_motion(state, properties, wrench)
        .unwrap()
        .point_motion(arm)
        .unwrap();

    let principal_to_body =
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.7).unwrap();
    let [w, x, y, z] = principal_to_body.components();
    let body_to_world = UnitQuaternion::new(w, -x, -y, -z).unwrap();
    let reexpressed_properties = MassProperties::from_principal_frame(
        properties.mass(),
        Vec3::ZERO,
        properties.principal_inertia_body(),
        principal_to_body,
    )
    .unwrap();
    let reexpressed_state = RigidBodyState::new(
        Pose::new(position, body_to_world).unwrap(),
        state.linear_momentum_world(),
        principal_to_body.rotate_body_to_world(state.angular_momentum_body()),
    )
    .unwrap();
    let actual = integrator
        .prepare_point_motion(
            reexpressed_state,
            reexpressed_properties,
            Wrench {
                force_world: wrench.force_world,
                torque_body: principal_to_body.rotate_body_to_world(wrench.torque_body),
            },
        )
        .unwrap()
        .point_motion(principal_to_body.rotate_body_to_world(arm))
        .unwrap();

    assert_vector_close(
        actual.kinematics.point_world,
        expected.kinematics.point_world,
        2e-12,
    );
    assert_vector_close(
        actual.kinematics.point_velocity_world,
        expected.kinematics.point_velocity_world,
        2e-12,
    );
    assert_vector_close(
        actual.angular_acceleration_body,
        principal_to_body.rotate_body_to_world(expected.angular_acceleration_body),
        2e-12,
    );
    assert_vector_close(
        actual.angular_acceleration_world,
        expected.angular_acceleration_world,
        2e-12,
    );
    assert_vector_close(
        actual.point_acceleration_world,
        expected.point_acceleration_world,
        2e-12,
    );
}

#[test]
fn motion_is_covariant_under_an_inertial_frame_rotation_and_translation() {
    let properties = properties();
    let position = Vec3::new(1.0, -2.0, 0.5);
    let velocity = Vec3::new(0.2, 0.4, -0.1);
    let omega = Vec3::new(1.2, -0.7, 0.9);
    let gravity = Vec3::new(0.1, -0.2, -3.0);
    let force = Vec3::new(1.0, 2.0, -0.7);
    let torque = Vec3::new(0.5, -0.4, 0.2);
    let arm = Vec3::new(0.4, -0.3, 0.2);
    let state = state_with_velocity(
        properties,
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity,
        omega,
    );
    let expected = RigidBodyIntegrator::new(Gravity::new(gravity).unwrap())
        .prepare_point_motion(
            state,
            properties,
            Wrench {
                force_world: force,
                torque_body: torque,
            },
        )
        .unwrap()
        .point_motion(arm)
        .unwrap();

    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 3.0), 1.1).unwrap();
    let translation = Vec3::new(3.0, 2.0, -5.0);
    let rotated_state = state_with_velocity(
        properties,
        Pose::new(
            rotation.rotate_body_to_world(position).add(translation),
            rotation,
        )
        .unwrap(),
        rotation.rotate_body_to_world(velocity),
        omega,
    );
    let actual =
        RigidBodyIntegrator::new(Gravity::new(rotation.rotate_body_to_world(gravity)).unwrap())
            .prepare_point_motion(
                rotated_state,
                properties,
                Wrench {
                    force_world: rotation.rotate_body_to_world(force),
                    torque_body: torque,
                },
            )
            .unwrap()
            .point_motion(arm)
            .unwrap();

    assert_vector_close(
        actual.kinematics.point_world,
        rotation
            .rotate_body_to_world(expected.kinematics.point_world)
            .add(translation),
        2e-12,
    );
    for (actual_vector, expected_vector) in [
        (
            actual.kinematics.point_velocity_world,
            expected.kinematics.point_velocity_world,
        ),
        (
            actual.center_of_mass_acceleration_world,
            expected.center_of_mass_acceleration_world,
        ),
        (
            actual.angular_acceleration_world,
            expected.angular_acceleration_world,
        ),
        (
            actual.point_acceleration_world,
            expected.point_acceleration_world,
        ),
    ] {
        assert_vector_close(
            actual_vector,
            rotation.rotate_body_to_world(expected_vector),
            2e-12,
        );
    }
    assert_vector_close(
        actual.angular_acceleration_body,
        expected.angular_acceleration_body,
        1e-14,
    );
}

#[test]
fn native_accepted_steps_refine_to_the_reported_instantaneous_acceleration() {
    let properties = properties();
    let pose = Pose::new(
        Vec3::new(1.0, -2.0, 0.3),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.6).unwrap(),
    )
    .unwrap();
    let state = state_with_velocity(
        properties,
        pose,
        Vec3::new(0.2, 0.4, -0.1),
        Vec3::new(1.2, -0.7, 0.9),
    );
    let integrator = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, 0.0, -3.0)).unwrap());
    let wrench = Wrench {
        force_world: Vec3::new(1.0, 2.0, -0.7),
        torque_body: Vec3::new(0.5, -0.4, 0.2),
    };
    let arm = Vec3::new(0.4, -0.3, 0.2);
    let motion = integrator
        .prepare_point_motion(state, properties, wrench)
        .unwrap()
        .point_motion(arm)
        .unwrap();
    let velocity_before = state
        .point_kinematics(properties, arm)
        .unwrap()
        .point_velocity_world;
    let mut errors = [0.0; 3];
    for (index, duration) in [0.004, 0.002, 0.001].into_iter().enumerate() {
        let accepted = integrator
            .step(state, properties, wrench, duration)
            .unwrap();
        assert_eq!(accepted.state_before, state);
        let velocity_after = accepted
            .state_after
            .point_kinematics(properties, arm)
            .unwrap()
            .point_velocity_world;
        let secant_acceleration = velocity_after.sub(velocity_before).scale(duration.recip());
        errors[index] = secant_acceleration
            .sub(motion.point_acceleration_world)
            .norm_squared()
            .sqrt();
    }

    // A forward velocity secant converges at first order to the instantaneous
    // acceleration, using actual accepted steps from the native integrator.
    assert!(
        errors[1] < 0.65 * errors[0],
        "acceleration errors: {errors:?}"
    );
    assert!(
        errors[2] < 0.65 * errors[1],
        "acceleration errors: {errors:?}"
    );
    assert!(errors[2] < 0.01, "acceleration errors: {errors:?}");
}

#[test]
fn rejected_inputs_leave_a_prepared_motion_snapshot_reusable() {
    let properties = properties();
    let state = state_with_velocity(
        properties,
        Pose::identity(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 3.0),
    );
    let before = (state, properties);
    let integrator = RigidBodyIntegrator::new(Gravity::ZERO);
    let prepared = integrator
        .prepare_point_motion(state, properties, Wrench::ZERO)
        .unwrap();
    let arm = Vec3::new(2.0, 0.0, 0.0);
    let expected = prepared.point_motion(arm).unwrap();

    for wrench in [
        Wrench {
            force_world: Vec3::new(f64::NAN, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(0.0, f64::INFINITY, 0.0),
        },
    ] {
        assert!(
            integrator
                .prepare_point_motion(state, properties, wrench)
                .is_err()
        );
    }
    assert!(
        prepared
            .point_motion(Vec3::new(0.0, f64::NAN, 0.0))
            .is_err()
    );
    assert!(
        prepared
            .point_motion(Vec3::new(f64::MAX, 0.0, 0.0))
            .is_err()
    );

    let actual = prepared.point_motion(arm).unwrap();
    assert_eq!(actual.kinematics, expected.kinematics);
    assert_eq!(
        actual.point_acceleration_world,
        expected.point_acceleration_world
    );
    assert_eq!(
        actual.angular_acceleration_body,
        expected.angular_acceleration_body
    );
    assert_eq!(
        actual.angular_acceleration_world,
        expected.angular_acceleration_world
    );
    assert_eq!(
        actual.center_of_mass_acceleration_world,
        expected.center_of_mass_acceleration_world
    );
    assert_eq!((state, properties), before);
}
