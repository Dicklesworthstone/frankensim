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
fn invalid_relative_inputs_and_overflow_refuse_without_changing_native_snapshot() {
    let carrier = prepared(0.0, 2.0, 0.0);
    let before = carrier;
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(carrier.relative_point_motion(Vec3::ZERO, Vec3::new(bad, 0.0, 0.0),
            Vec3::ZERO).is_err());
        assert!(carrier.relative_point_motion(Vec3::ZERO, Vec3::ZERO,
            Vec3::new(0.0, bad, 0.0)).is_err());
    }
    assert!(carrier.relative_point_motion(Vec3::ZERO, Vec3::new(1e308, 0.0, 0.0),
        Vec3::ZERO).is_err());
    assert_eq!(carrier, before);
    assert!(carrier.point_motion(Vec3::ZERO).is_ok());
}
