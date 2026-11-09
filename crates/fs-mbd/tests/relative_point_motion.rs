//! Native moving-coordinate queries, checked against Cartesian derivatives and
//! the actual equivalent articulated child, not a second transport helper.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Vec3 as GVec, Wrench as GWrench};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, JointModel, Link, SpatialInertia, forward_dynamics,
    forward_kinematics, inverse_dynamics,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion,
    Vec3, Wrench,
};

fn close(actual: [f64; 3], expected: [f64; 3]) {
    for (a, b) in actual.into_iter().zip(expected) {
        assert!((a - b).abs() < 5e-13 * b.abs().max(1.0), "{a:e} != {b:e}");
    }
}
fn components(v: Vec3) -> [f64; 3] { [v.x, v.y, v.z] }
fn gcomps(v: GVec) -> [f64; 3] { [v.x, v.y, v.z] }
fn rz(theta: f64, v: Vec3) -> Vec3 {
    let (s, c) = theta.sin_cos();
    Vec3::new(c * v.x - s * v.y, s * v.x + c * v.y, v.z)
}
fn prepared(theta: f64, omega: f64, alpha: f64) -> fs_mbd::PreparedRigidPointMotion {
    let properties = MassProperties::new(2.5, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let gravity = Vec3::new(0.0, -9.81, 0.0);
    let state = RigidBodyState::new(
        Pose::new(
            Vec3::new(2.0, -1.0, 0.5),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), theta).unwrap(),
        ).unwrap(),
        Vec3::new(7.0, -4.0, 2.0).scale(properties.mass()),
        Vec3::new(0.0, 0.0, 4.0 * omega),
    ).unwrap();
    RigidBodyIntegrator::new(Gravity::new(gravity).unwrap())
        .prepare_point_motion(state, properties, Wrench {
            force_world: Vec3::new(0.3, -0.2, 0.5).sub(gravity).scale(properties.mass()),
            torque_body: Vec3::new(0.0, 0.0, 4.0 * alpha),
        }).unwrap()
}

#[test]
fn native_velocity_and_acceleration_match_differentiated_rotating_coordinates() {
    let r = Vec3::new(0.8, -0.4, 0.2);
    let d = Vec3::new(-0.7, 0.9, -0.2);
    let dd = Vec3::new(0.4, -0.5, 0.1);
    for theta in [-1.1, 0.0, 0.6] {
        for omega in [-3.0, 0.0, 2.0] {
            let alpha = -0.3;
            let carrier = prepared(theta, omega, alpha);
            let point = carrier.relative_point_motion(r, d, dd).unwrap();
            // Differentiate x(t)=x_com(t)+Rz(theta(t))*r(t) componentwise.
            let v = rz(theta, Vec3::new(d.x - omega * r.y, d.y + omega * r.x, d.z))
                .add(Vec3::new(7.0, -4.0, 2.0));
            let a = rz(theta, Vec3::new(
                dd.x - alpha * r.y - omega * omega * r.x - 2.0 * omega * d.y,
                dd.y + alpha * r.x - omega * omega * r.y + 2.0 * omega * d.x,
                dd.z,
            )).add(Vec3::new(0.3, -0.2, 0.5));
            close(components(point.kinematics.point_velocity_world), components(v));
            let rigid = carrier.point_motion(r).unwrap().kinematics;
            close(components(point.kinematics.relative_velocity_world(rigid).unwrap()),
                components(rz(theta, d)));
            close(components(point.kinematics.velocity_relative_to_world(
                Vec3::new(7.0, -4.0, 2.0)).unwrap()),
                components(rz(theta, Vec3::new(
                    d.x - omega * r.y, d.y + omega * r.x, d.z))));
            close(components(point.point_acceleration_world), components(a));
            close(components(point.kinematics.point_world),
                components(Vec3::new(2.0, -1.0, 0.5).add(rz(theta, r))));
            assert_eq!(carrier.relative_point_motion(r, Vec3::ZERO, Vec3::ZERO).unwrap(),
                carrier.point_motion(r).unwrap());
        }
    }
}

fn slider() -> ArticulatedModel {
    let inertia = || SpatialInertia::new(2.0, GVec::default(), Mat3 {
        m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
    }).unwrap();
    ArticulatedModel::new(vec![
        Link::new("arm", None, Se3::identity(),
            JointModel::revolute(GVec::new(0.0, 0.0, 1.0), None).unwrap(), inertia()),
        Link::new("slider", Some(0), Se3::identity(),
            JointModel::prismatic(GVec::new(1.0, 0.0, 0.0), None).unwrap(), inertia()),
    ]).unwrap()
}

#[test]
fn relative_query_matches_real_slider_child_in_inverse_and_forward_dynamics() {
    let model = slider();
    let base = BaseState::stationary(Se3::from_parts(
        So3::exp(So3Tangent::new(GVec::new(0.3, -0.2, 0.5))).unwrap(),
        GVec::new(3.0, -1.0, 2.0),
    ).unwrap());
    let gravity = GVec::new(0.1, -9.81, 0.3);
    let q = [0.4, 1.3];
    let rate = [2.2, -0.7];
    let acceleration = [-0.3, 0.9];
    let arm = GVec::new(0.25, 0.15, -0.1);
    let kin = forward_kinematics(&model, base, &q, &rate).unwrap();
    let inverse = inverse_dynamics(&model, base, &q, &rate, &acceleration, gravity,
        &[GWrench::default(); 2]).unwrap();
    let forward = forward_dynamics(&model, base, &q, &rate, &inverse.generalized_force,
        gravity, &[GWrench::default(); 2]).unwrap();
    for (parent, child) in [
        (inverse.prepare_point_motion(&kin, 0, gravity).unwrap(),
         inverse.prepare_point_motion(&kin, 1, gravity).unwrap()),
        (forward.prepare_point_motion(&kin, 0, gravity).unwrap(),
         forward.prepare_point_motion(&kin, 1, gravity).unwrap()),
    ] {
        let actual = parent.relative_point_motion(
            GVec::new(q[1] + arm.x, arm.y, arm.z),
            GVec::new(rate[1], 0.0, 0.0),
            GVec::new(acceleration[1], 0.0, 0.0),
        ).unwrap();
        let expected = child.point_motion(arm).unwrap();
        close(gcomps(actual.point_world), gcomps(expected.point_world));
        close(gcomps(actual.point_velocity_world), gcomps(expected.point_velocity_world));
        close(gcomps(actual.point_acceleration_world), gcomps(expected.point_acceleration_world));
        assert_eq!(parent.relative_point_motion(arm, GVec::default(), GVec::default()).unwrap(),
            parent.point_motion(arm).unwrap());
        assert!(parent.relative_point_motion(arm, GVec::new(f64::NAN, 0.0, 0.0),
            GVec::default()).is_err());
        assert!(parent.relative_point_motion(arm, GVec::default(),
            GVec::new(0.0, f64::INFINITY, 0.0)).is_err());
        assert!(parent.relative_point_motion(GVec::default(), GVec::new(1e308, 0.0, 0.0),
            GVec::default()).is_err());
    }
}

#[test]
fn relative_motion_does_not_invent_centripetal_acceleration_from_slip() {
    // A stationary carrier with a moving point has exactly the authored local
    // velocity/acceleration, not a fictitious spin inferred from that velocity.
    let carrier = prepared(0.0, 0.0, 0.0);
    let point = carrier.relative_point_motion(Vec3::new(2.0, 0.0, 0.0),
        Vec3::new(0.0, 3.0, 0.0), Vec3::ZERO).unwrap();
    close(components(point.kinematics.point_velocity_world), [7.0, -1.0, 2.0]);
    close(components(point.point_acceleration_world), [0.3, -0.2, 0.5]);
}

#[test]
fn relative_queries_retain_coordinate_motion_under_common_translation() {
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    for boost in [0.0, 2.0_f64.powi(60)] {
        let translation = Vec3::new(boost, -boost, boost);
        let state =
            RigidBodyState::new(Pose::identity(), translation, Vec3::new(0.0, 0.0, 2.0)).unwrap();
        let carrier = RigidBodyIntegrator::new(Gravity::new(Vec3::ZERO).unwrap())
            .prepare_point_motion(state, properties, Wrench::default())
            .unwrap();
        let arm_a = Vec3::new(1.0, 0.0, 0.0);
        let a = carrier
            .relative_point_motion(arm_a, Vec3::new(0.0, 3.0, -0.5), Vec3::ZERO)
            .unwrap()
            .kinematics;
        let b = carrier
            .relative_point_motion(
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(4.0, 0.0, 0.25),
                Vec3::ZERO,
            )
            .unwrap()
            .kinematics;
        let rigid = carrier.point_motion(arm_a).unwrap().kinematics;
        // d(Rz(2t) * r(t))/dt: a has rotational velocity [0,2,0],
        // b has [-2,0,0], before either authored coordinate rate is added.
        assert_eq!(
            a.velocity_relative_to_world(translation).unwrap(),
            Vec3::new(0.0, 5.0, -0.5)
        );
        assert_eq!(
            b.velocity_relative_to_world(translation).unwrap(),
            Vec3::new(2.0, 0.0, 0.25)
        );
        assert_eq!(
            a.relative_velocity_world(b).unwrap(),
            Vec3::new(-2.0, 5.0, -0.75)
        );
        assert_eq!(
            b.relative_velocity_world(a).unwrap(),
            Vec3::new(2.0, -5.0, 0.75)
        );
        assert_eq!(
            a.relative_velocity_world(rigid).unwrap(),
            Vec3::new(0.0, 3.0, -0.5)
        );
        assert_eq!(a.relative_velocity_world(a).unwrap(), Vec3::ZERO);
        assert_eq!(a.center_of_mass_velocity_world, translation);
        assert_eq!(a.angular_velocity_world, Vec3::new(0.0, 0.0, 2.0));
        assert_eq!(rigid.relative_point_velocity_world, Vec3::ZERO);
        if boost != 0.0 {
            // Absolute f64 velocities cannot represent the small terms here.
            assert_eq!(a.point_velocity_world, translation);
            assert_eq!(b.point_velocity_world, translation);
        }
    }
}

#[test]
fn absolute_query_retains_spin_when_coordinate_motion_cancels_translation() {
    let boost = 2.0_f64.powi(60);
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let state = RigidBodyState::new(
        Pose::identity(),
        Vec3::new(boost, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    )
    .unwrap();
    let carrier = RigidBodyIntegrator::new(Gravity::new(Vec3::ZERO).unwrap())
        .prepare_point_motion(state, properties, Wrench::default())
        .unwrap();
    let point = carrier
        .relative_point_motion(
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(-boost, 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap()
        .kinematics;
    assert_eq!(point.point_velocity_world, Vec3::new(1.0, 0.0, 0.0));
    assert_eq!(
        point.velocity_relative_to_world(Vec3::ZERO).unwrap(),
        point.point_velocity_world
    );
}

#[test]
fn invalid_relative_inputs_and_overflow_refuse_without_changing_native_snapshot() {
    let carrier = prepared(0.0, 2.0, 0.0);
    let before = carrier;
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            carrier
                .relative_point_motion(Vec3::ZERO, Vec3::new(bad, 0.0, 0.0), Vec3::ZERO)
                .is_err()
        );
        assert!(
            carrier
                .relative_point_motion(Vec3::ZERO, Vec3::ZERO, Vec3::new(0.0, bad, 0.0))
                .is_err()
        );
    }
    assert!(
        carrier
            .relative_point_motion(Vec3::ZERO, Vec3::new(1e308, 0.0, 0.0), Vec3::ZERO)
            .is_err()
    );
    assert_eq!(carrier, before);
    assert!(carrier.point_motion(Vec3::ZERO).is_ok());
    let valid = carrier.point_motion(Vec3::ZERO).unwrap().kinematics;
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut invalid = valid;
        invalid.relative_point_velocity_world = Vec3::new(0.0, bad, 0.0);
        assert!(invalid.velocity_relative_to_world(Vec3::ZERO).is_err());
        assert!(invalid.relative_velocity_world(valid).is_err());
        assert!(valid.relative_velocity_world(invalid).is_err());
    }
}
