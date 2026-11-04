//! Native articulated motion reaches the original paired acoustic histories.

use fs_aeroac::AeroacError;
use fs_couple::pcm_wav::observation::moving::rigid::{
    BodyFrameCompactSourceStrength, RigidAcousticAttachment, RigidAcousticError,
    RigidAcousticHistory, RigidHistoryBudget,
};
use fs_couple::pcm_wav::observation::moving::rigid::articulated::ArticulatedAcousticError;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link,
    PreparedLinkPointMotion, SpatialInertia, forward_dynamics, forward_kinematics,
    free_floating_forward_dynamics, inverse_dynamics,
};

fn inertia() -> SpatialInertia {
    // Deliberately not centred on the link origin: the acoustic attachment is
    // a link coordinate, not an independent rigid body's COM-relative arm.
    SpatialInertia::new(
        2.0,
        Vec3::new(0.3, -0.2, 0.1),
        Mat3 { m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0] },
    ).unwrap()
}
fn slider() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new("arm", None, Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(), inertia()),
        Link::new("slider", Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(), inertia()),
    ]).unwrap()
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn attachment(point: Vec3) -> RigidAcousticAttachment {
    RigidAcousticAttachment::new(fs_mbd::Vec3::new(point.x, point.y, point.z)).unwrap()
}
fn components(v: Vec3) -> [f64; 3] { [v.x, v.y, v.z] }
fn near(actual: [f64; 3], expected: Vec3, tolerance: f64) {
    for (a, b) in actual.into_iter().zip(components(expected)) {
        assert!((a-b).abs() <= tolerance, "{a:e} != {b:e}");
    }
}
fn history(samples: usize) -> RigidAcousticHistory {
    RigidAcousticHistory::new(attachment(Vec3::new(0.25, 0.15, -0.1)),
        RigidHistoryBudget { maximum_samples: samples, maximum_payload_bytes: 262_144 }).unwrap()
}
fn quiet_motion() -> PreparedLinkPointMotion {
    let model = slider();
    let base = BaseState::stationary(Se3::identity());
    let kinematics = forward_kinematics(&model, base, &[0.4, 1.3], &[0.0; 2]).unwrap();
    inverse_dynamics(&model, base, &[0.4, 1.3], &[0.0; 2], &[0.0; 2],
        Vec3::default(), &[Wrench::default(); 2]).unwrap()
        .prepare_point_motion(&kinematics, 1, Vec3::default()).unwrap()
}

#[test]
fn rotating_slider_transports_native_point_motion_and_link_force_derivatives() {
    let model = slider();
    let rot = rotation(Vec3::new(0.3, -0.2, 0.5));
    let origin = Vec3::new(2.0, -1.0, 3.0);
    let base_acceleration = Vec3::new(0.2, 0.1, -0.4);
    let gravity = Vec3::new(0.1, -9.81, 0.3);
    let q = [0.4, 1.3];
    let rate = [2.2, -0.7];
    let acceleration = [-0.3, 0.9];
    let point = Vec3::new(0.25, 0.15, -0.1);
    let strength = BodyFrameCompactSourceStrength {
        volume_flow_m3_s: 0.002, volume_flow_rate_m3_s2: -0.003,
        force_on_fluid_body_n: [2.0, -1.0, 0.5],
        force_rate_on_fluid_body_n_s: [0.3, 0.2, -0.1],
    };
    // Independent Cartesian derivatives of R(theta)[rho+x,y,z].
    let theta = rotation(Vec3::new(0.0, 0.0, q[0]));
    let world = |v| rot.rotate(theta.rotate(v).unwrap()).unwrap();
    let radius = q[1] + point.x;
    let expected_position = origin + world(Vec3::new(radius, point.y, point.z));
    let expected_acceleration = base_acceleration + world(Vec3::new(
        acceleration[1] - acceleration[0]*point.y - rate[0]*rate[0]*radius,
        2.0*rate[0]*rate[1] + acceleration[0]*radius - rate[0]*rate[0]*point.y, 0.0));
    let expected_force = world(Vec3::new(2.0, -1.0, 0.5));
    let expected_force_rate = world(Vec3::new(0.3 + rate[0], 0.2 + 2.0*rate[0], -0.1));
    for boost in [Vec3::default(), Vec3::new(30.0, -40.0, 10.0)] {
        let base_velocity = Vec3::new(1.0, -2.0, 3.0) + boost;
        let base = BaseState::prescribed(Se3::from_parts(rot, origin).unwrap(),
            Twist::new(Vec3::default(), rot.inverse().rotate(base_velocity).unwrap()),
            Twist::new(Vec3::default(), rot.inverse().rotate(base_acceleration).unwrap()));
        let kinematics = forward_kinematics(&model, base, &q, &rate).unwrap();
        let inverse = inverse_dynamics(&model, base, &q, &rate, &acceleration,
            gravity, &[Wrench::default(); 2]).unwrap();
        let forward = forward_dynamics(&model, base, &q, &rate, &inverse.generalized_force,
            gravity, &[Wrench::default(); 2]).unwrap();
        let expected_velocity = base_velocity + world(Vec3::new(
            rate[1] - rate[0]*point.y, rate[0]*radius, 0.0));
        for motion in [inverse.prepare_point_motion(&kinematics, 1, gravity).unwrap(),
            forward.prepare_point_motion(&kinematics, 1, gravity).unwrap()] {
            let sample = attachment(point).sample_link(0.7, &motion, strength).unwrap();
            assert_eq!(sample.kinematics.source_time_s, sample.source.source_time_s);
            assert_eq!(sample.source.source_time_s, 0.7);
            assert_eq!(attachment(point).kinematics_link(0.7, &motion).unwrap(), sample.kinematics);
            near(sample.kinematics.kinematics.position_m, expected_position, 2e-12);
            near(sample.kinematics.kinematics.velocity_m_s, expected_velocity, 2e-12);
            near(sample.kinematics.kinematics.acceleration_m_s2, expected_acceleration, 3e-11);
            near(sample.source.strength.force_on_fluid_n, expected_force, 2e-12);
            near(sample.source.strength.force_rate_on_fluid_n_s, expected_force_rate, 2e-12);
            assert_eq!(sample.source.strength.volume_flow_m3_s, strength.volume_flow_m3_s);
            assert_eq!(sample.source.strength.volume_flow_rate_m3_s2, strength.volume_flow_rate_m3_s2);
        }
    }
}

#[test]
fn supported_and_falling_microphones_keep_physical_not_specific_acceleration() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body", None, Se3::identity(), JointModel::FIXED, inertia())]).unwrap();
    let pose = Se3::from_parts(rotation(Vec3::new(0.2, 0.4, -0.1)), Vec3::new(1.0, 2.0, 3.0)).unwrap();
    let base = BaseState::stationary(pose);
    let gravity = Vec3::new(0.0, 0.0, -9.81);
    let kin = forward_kinematics(&model, base, &[], &[]).unwrap();
    let supported = inverse_dynamics(&model, base, &[], &[], &[], gravity,
        &[Wrench::default()]).unwrap().prepare_point_motion(&kin, 0, gravity).unwrap();
    let falling = free_floating_forward_dynamics(&model,
        FreeFloatingBaseState::stationary(pose), &[], &[], &[], gravity,
        &[Wrench::default()]).unwrap().prepare_point_motion(&kin, 0).unwrap();
    for (motion, acceleration) in [(supported, Vec3::default()), (falling, gravity)] {
        let sample = attachment(Vec3::new(0.2, -0.1, 0.3))
            .sample_link(0.0, &motion, BodyFrameCompactSourceStrength::default()).unwrap();
        near(sample.kinematics.kinematics.acceleration_m_s2, acceleration, 3e-12);
        // Neither gravity nor the support wrench was invented as a fluid dipole.
        assert_eq!(sample.source.strength.force_on_fluid_n, [0.0; 3]);
        assert_eq!(sample.source.strength.force_rate_on_fluid_n_s, [0.0; 3]);
    }
}

#[test]
fn paired_link_history_rejects_late_cancellation_bad_source_clock_and_capacity_atomically() {
    let motion = quiet_motion();
    let saved_motion = motion;
    let strength = BodyFrameCompactSourceStrength::default();
    let mut h = history(2);
    let motion_pointer = h.kinematics().as_ptr();
    let source_pointer = h.sources().as_ptr();
    h.append_link(0.0, &motion, strength, &mut || false).unwrap();
    let before = (h.kinematics().to_vec(), h.sources().to_vec());
    let mut polls = 0;
    assert!(matches!(h.append_link(0.1, &motion, strength, &mut || {polls += 1; polls == 2}),
        Err(ArticulatedAcousticError::Observation(RigidAcousticError::Acoustic(
            AeroacError::Cancelled { .. })))));
    assert_eq!(polls, 2);
    for time in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(h.append_link(time, &motion, strength, &mut || false).is_err());
    }
    assert!(h.append_link(0.1, &motion, BodyFrameCompactSourceStrength {
        force_rate_on_fluid_body_n_s: [0.0, f64::NAN, 0.0], ..strength
    }, &mut || false).is_err());
    assert_eq!(h.kinematics(), before.0.as_slice());
    assert_eq!(h.sources(), before.1.as_slice());
    assert_eq!(h.remaining_capacity(), 1);
    h.append_link(0.1, &motion, strength, &mut || false).unwrap();
    assert!(h.append_link(0.2, &motion, strength, &mut || false).is_err());
    assert_eq!(h.len(), 2);
    assert_eq!(h.remaining_capacity(), 0);
    assert_eq!(motion, saved_motion);
    assert_eq!(motion_pointer, h.kinematics().as_ptr());
    assert_eq!(source_pointer, h.sources().as_ptr());
}

#[test]
fn overflowing_rotating_source_refuses_before_either_history_is_published() {
    let model = slider();
    let base = BaseState::stationary(Se3::identity());
    let kin = forward_kinematics(&model, base, &[0.0, 1.0], &[100.0, 0.0]).unwrap();
    let dynamics = inverse_dynamics(&model, base, &[0.0, 1.0], &[100.0, 0.0], &[0.0; 2],
        Vec3::default(), &[Wrench::default(); 2]).unwrap();
    let motion = dynamics.prepare_point_motion(&kin, 1, Vec3::default()).unwrap();
    let mut h = history(2);
    assert!(h.append_link(0.0, &motion, BodyFrameCompactSourceStrength {
        force_on_fluid_body_n: [1e308, 0.0, 0.0], ..Default::default()
    }, &mut || false).is_err());
    assert!(h.is_empty());
    assert!(h.kinematics().is_empty());
}

#[path = "articulated_acoustic_history/render.rs"]
mod render;
