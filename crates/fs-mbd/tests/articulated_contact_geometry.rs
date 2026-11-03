//! Original material-witness collinearity retains its transform input scales.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Vec3};
use fs_mbd::articulated::contact_geometry::{
    ContactCounterpart, ContactGeometryError, LinkContactPoint, prepare_contact_frame,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, JointModel, Link, SpatialInertia, forward_kinematics,
};

fn magnitude(value: Vec3) -> f64 {
    value.x.abs().max(value.y.abs()).max(value.z.abs())
}

#[test]
fn rotated_contact_at_world_origin_keeps_covariance_and_rejects_real_misalignment() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        SpatialInertia::new(
            1.0,
            Vec3::default(),
            Mat3 {
                m: [0.4, 0., 0., 0., 0.4, 0., 0., 0., 0.4],
            },
        )
        .unwrap(),
    )])
    .unwrap();
    let scene_rotations = [
        So3::identity(),
        So3::exp(So3Tangent {
            angular: Vec3::new(0.2, -0.3, 0.4),
        })
        .unwrap(),
    ];
    let mut exposes_cancelled_world_scale = false;
    for angle in [-0.375, 0.37, 0.83] {
        let shell = So3::exp(So3Tangent {
            angular: Vec3::new(0., 0., angle),
        })
        .unwrap();
        let first = LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: shell.inverse().rotate(Vec3::new(1., 0., 0.)).unwrap(),
        };
        let second = LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::new(-1., 0., 0.),
        };
        for rotation in scene_rotations {
            for origin in [Vec3::default(), Vec3::new(7., -2., 1.)] {
                let first_pose = Se3::from_parts(
                    rotation.compose(shell).unwrap(),
                    origin + rotation.rotate(Vec3::new(-1., 0., 0.)).unwrap(),
                )
                .unwrap();
                let second_pose = Se3::from_parts(
                    rotation,
                    origin + rotation.rotate(Vec3::new(1., 0., 0.)).unwrap(),
                )
                .unwrap();
                let kinematics = [
                    forward_kinematics(&model, BaseState::stationary(first_pose), &[], &[])
                        .unwrap(),
                    forward_kinematics(&model, BaseState::stationary(second_pose), &[], &[])
                        .unwrap(),
                ];
                let normal = rotation.rotate(Vec3::new(-1., 0., 0.)).unwrap();
                let tangent = rotation.rotate(Vec3::new(0., 1., 0.)).unwrap();
                let frame = prepare_contact_frame(
                    first,
                    ContactCounterpart::Body(second),
                    normal,
                    tangent,
                    17,
                    &kinematics,
                )
                .unwrap();
                let separation = frame.rows[0].world_point - frame.counterpart_point_world;
                let off_line = separation - frame.rows[0].axis_world.scale(frame.rows[0].gap);
                let old_scale = magnitude(frame.rows[0].world_point)
                    .max(magnitude(frame.counterpart_point_world))
                    .max(magnitude(separation));
                exposes_cancelled_world_scale |= origin == Vec3::default()
                    && magnitude(off_line) > 256.0 * f64::EPSILON * old_scale;
                assert!((frame.rows[0].world_point - origin).norm() < 1e-14);
                assert!((frame.counterpart_point_world - origin).norm() < 1e-14);
                assert!(frame.rows[0].gap.abs() < 1e-14);
                assert!((frame.rows[0].axis_world - normal).norm() < 1e-14);

                // The original source scale supplies a rounding allowance,
                // not an arbitrary geometric tolerance. A real 10-micron
                // material-point offset still refuses in every scene frame.
                let displaced = LinkContactPoint {
                    point_body: first.point_body
                        + first_pose
                            .rotation()
                            .inverse()
                            .rotate(tangent.scale(1e-5))
                            .unwrap(),
                    ..first
                };
                assert!(matches!(
                    prepare_contact_frame(
                        displaced,
                        ContactCounterpart::Body(second),
                        normal,
                        tangent,
                        17,
                        &kinematics,
                    ),
                    Err(ContactGeometryError::NonCollinearContact { index: 17 })
                ));
            }
        }
    }
    assert!(
        exposes_cancelled_world_scale,
        "fixture must expose real transform roundoff hidden by cancelled world points"
    );
}
