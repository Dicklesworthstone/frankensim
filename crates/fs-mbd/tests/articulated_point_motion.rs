//! G1/G3 physical acceleration at original articulated material points.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link,
    SpatialInertia, forward_dynamics, forward_kinematics, free_floating_forward_dynamics,
    inverse_dynamics,
};

fn vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        (actual - expected).norm() <= tolerance,
        "{actual:?} != {expected:?}; error {}",
        (actual - expected).norm()
    );
}

fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}

fn inertia() -> SpatialInertia {
    SpatialInertia::new(
        2.0,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
        },
    )
    .unwrap()
}

fn slider() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "rotating arm",
            None,
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(),
        ),
        Link::new(
            "radial slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(),
        ),
    ])
    .unwrap()
}

fn fixed_body() -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(),
    )])
    .unwrap()
}

#[test]
fn rotating_slider_retains_coriolis_tangential_and_centripetal_acceleration() {
    let model = slider();
    let rot = rotation(Vec3::new(0.3, -0.2, 0.5));
    let base_velocity = Vec3::new(1.0, -2.0, 3.0);
    let base_acceleration = Vec3::new(0.2, 0.1, -0.4);
    let base = BaseState::prescribed(
        Se3::from_parts(rot, Vec3::new(2.0, -1.0, 3.0)).unwrap(),
        Twist::new(
            Vec3::default(),
            rot.inverse().rotate(base_velocity).unwrap(),
        ),
        Twist::new(
            Vec3::default(),
            rot.inverse().rotate(base_acceleration).unwrap(),
        ),
    );
    let gravity = Vec3::new(0.1, -9.81, 0.3);
    let q = [0.4, 1.3];
    let rate = [2.2, -0.7];
    let acceleration = [-0.3, 0.9];
    let point = Vec3::new(0.25, 0.15, -0.1);
    let kinematics = forward_kinematics(&model, base, &q, &rate).unwrap();
    let inverse = inverse_dynamics(
        &model,
        base,
        &q,
        &rate,
        &acceleration,
        gravity,
        &[Wrench::default(); 2],
    )
    .unwrap();
    let forward = forward_dynamics(
        &model,
        base,
        &q,
        &rate,
        &inverse.generalized_force,
        gravity,
        &[Wrench::default(); 2],
    )
    .unwrap();
    let from_inverse = inverse
        .prepare_point_motion(&kinematics, 1, gravity)
        .unwrap();
    let from_forward = forward
        .prepare_point_motion(&kinematics, 1, gravity)
        .unwrap();
    let theta = rotation(Vec3::new(0.0, 0.0, q[0]));
    let into_world = |v| rot.rotate(theta.rotate(v).unwrap()).unwrap();
    let radius = q[1] + point.x;
    // Direct Cartesian differentiation of R(theta) * [rho + x, y, z].
    let velocity = into_world(Vec3::new(
        rate[1] - rate[0] * point.y,
        rate[0] * radius,
        0.0,
    )) + base_velocity;
    let physical_acceleration = into_world(Vec3::new(
        acceleration[1] - acceleration[0] * point.y - rate[0] * rate[0] * radius,
        2.0 * rate[0] * rate[1] + acceleration[0] * radius - rate[0] * rate[0] * point.y,
        0.0,
    )) + base_acceleration;
    for prepared in [from_inverse, from_forward] {
        let motion = prepared.point_motion(point).unwrap();
        vector(motion.point_velocity_world, velocity, 2e-13);
        vector(
            motion.point_acceleration_world,
            physical_acceleration,
            2e-12,
        );
        vector(
            motion.angular_velocity_world,
            into_world(Vec3::new(0.0, 0.0, rate[0])),
            2e-13,
        );
        vector(
            motion.angular_acceleration_world,
            into_world(Vec3::new(0.0, 0.0, acceleration[0])),
            2e-13,
        );
    }
}

#[test]
fn supported_point_and_freely_falling_point_have_distinct_accelerometer_readings() {
    let model = fixed_body();
    let pose = Se3::from_parts(
        rotation(Vec3::new(0.2, 0.4, -0.1)),
        Vec3::new(1.0, 2.0, 3.0),
    )
    .unwrap();
    let base = BaseState::stationary(pose);
    let gravity = Vec3::new(0.0, 0.0, -9.81);
    let kinematics = forward_kinematics(&model, base, &[], &[]).unwrap();
    let supported = inverse_dynamics(&model, base, &[], &[], &[], gravity, &[Wrench::default()])
        .unwrap()
        .prepare_point_motion(&kinematics, 0, gravity)
        .unwrap()
        .point_motion(Vec3::new(0.2, -0.1, 0.3))
        .unwrap();
    vector(supported.point_acceleration_world, Vec3::default(), 3e-14);
    vector(
        supported.specific_acceleration_world(gravity).unwrap(),
        gravity.scale(-1.0),
        3e-14,
    );
    let falling = free_floating_forward_dynamics(
        &model,
        FreeFloatingBaseState::stationary(pose),
        &[],
        &[],
        &[],
        gravity,
        &[Wrench::default()],
    )
    .unwrap()
    .prepare_point_motion(&kinematics, 0)
    .unwrap()
    .point_motion(Vec3::new(0.2, -0.1, 0.3))
    .unwrap();
    vector(falling.point_acceleration_world, gravity, 3e-14);
    vector(
        falling.specific_acceleration_world(gravity).unwrap(),
        Vec3::default(),
        3e-14,
    );
}

#[test]
fn boosted_spinning_free_body_has_physical_origin_and_attachment_acceleration() {
    let model = fixed_body();
    let rot = rotation(Vec3::new(0.3, -0.2, 0.1));
    let pose = Se3::from_parts(rot, Vec3::new(3.0, -2.0, 1.0)).unwrap();
    let omega = Vec3::new(0.0, 0.0, 2.0);
    let gravity = Vec3::new(0.0, -9.81, 0.0);
    let point = Vec3::new(0.5, 0.3, -0.1);
    let expected = gravity + rot.rotate(Vec3::new(-2.0, -1.2, 0.0)).unwrap();
    for velocity in [Vec3::default(), Vec3::new(30.0, -40.0, 10.0)] {
        let twist = Twist::new(omega, velocity);
        let kinematics = forward_kinematics(
            &model,
            BaseState::prescribed(pose, twist, Twist::zero()),
            &[],
            &[],
        )
        .unwrap();
        let prepared = free_floating_forward_dynamics(
            &model,
            FreeFloatingBaseState::new(pose, twist),
            &[],
            &[],
            &[],
            gravity,
            &[Wrench::default()],
        )
        .unwrap()
        .prepare_point_motion(&kinematics, 0)
        .unwrap();
        vector(
            prepared
                .point_motion(Vec3::default())
                .unwrap()
                .point_acceleration_world,
            gravity,
            2e-13,
        );
        let motion = prepared.point_motion(point).unwrap();
        vector(motion.point_acceleration_world, expected, 2e-13);
        vector(
            motion.point_velocity_world,
            rot.rotate(velocity + Vec3::new(-0.6, 1.0, 0.0)).unwrap(),
            2e-13,
        );
        vector(motion.angular_acceleration_world, Vec3::default(), 2e-13);
    }
}

#[test]
fn point_acceleration_matches_independently_differentiated_forward_kinematics() {
    let model = slider();
    let base = BaseState::stationary(
        Se3::from_parts(
            rotation(Vec3::new(0.3, 0.4, 0.1)),
            Vec3::new(1.0, -2.0, 3.0),
        )
        .unwrap(),
    );
    let q = [0.6, 1.2];
    let rate = [0.7, -0.4];
    let acceleration = [-0.3, 0.5];
    let point = Vec3::new(0.2, -0.3, 0.1);
    let gravity = Vec3::new(0.0, -9.81, 0.0);
    let velocity = |time: f64| {
        let position = core::array::from_fn::<_, 2, _>(|i| {
            q[i] + time * rate[i] + 0.5 * time * time * acceleration[i]
        });
        let rates = core::array::from_fn::<_, 2, _>(|i| rate[i] + time * acceleration[i]);
        let kinematics = forward_kinematics(&model, base, &position, &rates).unwrap();
        let twist = kinematics.body_twist[1];
        kinematics.world_from_link[1]
            .rotation()
            .rotate(twist.linear + twist.angular.cross(point))
            .unwrap()
    };
    let kinematics = forward_kinematics(&model, base, &q, &rate).unwrap();
    let actual = inverse_dynamics(
        &model,
        base,
        &q,
        &rate,
        &acceleration,
        gravity,
        &[Wrench::default(); 2],
    )
    .unwrap()
    .prepare_point_motion(&kinematics, 1, gravity)
    .unwrap()
    .point_motion(point)
    .unwrap();
    let h = 2e-5;
    let difference = (velocity(h) - velocity(-h)).scale(0.5 / h);
    vector(actual.point_acceleration_world, difference, 2e-9);
}

#[test]
fn malformed_receipts_and_unrepresentable_point_motion_refuse_without_mutation() {
    let model = fixed_body();
    let base = BaseState::stationary(Se3::identity());
    let gravity = Vec3::default();
    let kinematics = forward_kinematics(&model, base, &[], &[]).unwrap();
    let result =
        inverse_dynamics(&model, base, &[], &[], &[], gravity, &[Wrench::default()]).unwrap();
    assert!(matches!(
        result.prepare_point_motion(&kinematics, usize::MAX, gravity),
        Err(ArticulatedError::InvalidLinkIndex { .. })
    ));
    let mut malformed = kinematics.clone();
    malformed.body_twist.clear();
    assert!(matches!(
        result.prepare_point_motion(&malformed, 0, gravity),
        Err(ArticulatedError::LengthMismatch { .. })
    ));
    let mut invalid_result = result.clone();
    invalid_result.body_acceleration[0].linear.x = f64::NAN;
    assert!(
        invalid_result
            .prepare_point_motion(&kinematics, 0, gravity)
            .is_err()
    );
    assert!(
        result
            .prepare_point_motion(&kinematics, 0, Vec3::new(f64::NAN, 0.0, 0.0))
            .is_err()
    );
    let prepared = result
        .prepare_point_motion(&kinematics, 0, gravity)
        .unwrap();
    assert!(
        prepared
            .point_motion(Vec3::new(f64::NAN, 0.0, 0.0))
            .is_err()
    );
    let mut excessive = kinematics.clone();
    excessive.body_twist[0].angular = Vec3::new(0.0, 0.0, 1e160);
    let prepared = result.prepare_point_motion(&excessive, 0, gravity).unwrap();
    assert!(prepared.point_motion(Vec3::new(1.0, 0.0, 0.0)).is_err());
    assert_eq!(
        kinematics,
        forward_kinematics(&model, base, &[], &[]).unwrap()
    );
    assert_eq!(result.body_acceleration, vec![Twist::zero()]);
}
