//! Original material-witness collinearity retains its transform input scales.

use fs_contact::multicontact::MultiImpactConfig;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::contact_geometry::{
    ArticulatedNormalContact, ContactCounterpart, ContactGeometryError, LinkContactPoint,
    PointDirection, prepare_contact_frame, relative_point_speed,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Kinematics, Link,
    SpatialInertia, forward_kinematics,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, resolve_articulated_contacts,
};

fn magnitude(value: Vec3) -> f64 {
    value.x.abs().max(value.y.abs()).max(value.z.abs())
}

/// G0/G3: a common translation must not round away the independently stored
/// angular motion before normal and tangential relative velocities are formed.
#[test]
fn contact_rows_retain_rotational_slip_under_common_translation() {
    let arm = Vec3::new(1.0, 1.0, 1.0);
    let point = LinkContactPoint {
        system_index: 0,
        link_index: 0,
        point_body: arm,
    };
    let other = LinkContactPoint {
        system_index: 1,
        ..point
    };
    let rotations = [
        So3::identity(),
        So3::exp(So3Tangent {
            angular: Vec3::new(0.2, -0.3, 0.4),
        })
        .unwrap(),
    ];
    // The independent rigid-point formula is omega cross arm = (-1, 2, -1).
    let angular = Vec3::new(1.0, 2.0, 3.0);
    for boost in [
        0.0,
        32.0,
        -32.0,
        2.0_f64.powi(56),
        -2.0_f64.powi(56),
        2.0_f64.powi(256),
    ] {
        for (rotation_index, rotation) in rotations.into_iter().enumerate() {
            // The large case isolates summation of exactly represented data;
            // it does not claim to undo roundoff in a prior frame conversion.
            if boost.abs() > 32.0 && rotation_index != 0 {
                continue;
            }
            let linear = Vec3::new(boost, boost, boost);
            let pose = Se3::from_parts(rotation, Vec3::default()).unwrap();
            let kinematics = [
                Kinematics {
                    world_from_link: vec![pose],
                    body_twist: vec![Twist::new(angular, linear)],
                },
                Kinematics {
                    world_from_link: vec![pose],
                    body_twist: vec![Twist::new(Vec3::default(), linear)],
                },
            ];
            for counterpart in [
                ContactCounterpart::Body(other),
                ContactCounterpart::PrescribedPlane {
                    point_world: pose.transform_point(arm).unwrap(),
                    velocity_world: rotation.rotate(linear).unwrap(),
                },
            ] {
                let frame = prepare_contact_frame(
                    point,
                    counterpart,
                    rotation.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap(),
                    rotation.rotate(Vec3::new(0.0, 1.0, 0.0)).unwrap(),
                    0,
                    &kinematics,
                )
                .unwrap();
                for (row, expected) in frame.rows.iter().zip([-1.0, 2.0, -1.0]) {
                    assert!(
                        (row.relative_normal_speed - expected).abs() < 1e-12,
                        "common translation {boost:e}, counterpart {counterpart:?}: {row:?}"
                    );
                }
            }
        }
    }
}

/// Cross-crate endpoint evaluation must reject malformed data before indexing
/// while retaining exactly the same signed relative velocity as preparation.
#[test]
fn checked_relative_point_speed_rejects_invalid_consumed_inputs() {
    let arm = Vec3::new(1.0, 1.0, 1.0);
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let row = PointDirection {
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: arm,
        },
        direction: Wrench::new(arm.cross(axis), axis),
    };
    let kinematics = [Kinematics {
        world_from_link: vec![Se3::identity()],
        body_twist: vec![Twist::new(
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(32.0, 0.0, 0.0),
        )],
    }];
    assert_eq!(
        relative_point_speed(row, None, 32.0, &kinematics).unwrap(),
        -1.0
    );

    for point in [
        LinkContactPoint {
            system_index: 1,
            ..row.point
        },
        LinkContactPoint {
            link_index: 1,
            ..row.point
        },
    ] {
        let invalid = PointDirection { point, ..row };
        assert!(matches!(
            relative_point_speed(invalid, None, 0.0, &kinematics),
            Err(ContactGeometryError::Invalid(_))
        ));
        assert!(matches!(
            relative_point_speed(row, Some(invalid), 0.0, &kinematics),
            Err(ContactGeometryError::Invalid(_))
        ));
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            relative_point_speed(row, None, bad, &kinematics),
            Err(ContactGeometryError::Invalid(_))
        ));
        for input in 0..4 {
            let mut invalid_row = row;
            let mut invalid_kinematics = kinematics.clone();
            match input {
                0 => invalid_row.direction.force.y = bad,
                1 => invalid_row.direction.torque.z = bad,
                2 => invalid_kinematics[0].body_twist[0].linear.z = bad,
                _ => invalid_kinematics[0].body_twist[0].angular.x = bad,
            }
            assert!(matches!(
                relative_point_speed(invalid_row, None, 0.0, &invalid_kinematics),
                Err(ContactGeometryError::Invalid(_))
            ));
        }
    }
}

fn impact_body() -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        SpatialInertia::new(
            1.0,
            Vec3::default(),
            Mat3 {
                m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            },
        )
        .unwrap(),
    )])
    .unwrap()
}

fn off_center_impact(
    model: &ArticulatedModel,
    rotation: So3,
    boost: f64,
) -> (ArticulatedContactSystem<'_>, ArticulatedNormalContact) {
    let arm = Vec3::new(0.0, 1.0, 0.0);
    let pose = Se3::from_parts(rotation, Vec3::default()).unwrap();
    let linear = Vec3::new(boost, 0.0, 0.0);
    (
        ArticulatedContactSystem {
            model,
            base: FreeFloatingBaseState::new(pose, Twist::new(Vec3::new(0.0, 0.0, 1.0), linear)),
            joint_position: &[],
            joint_velocity: &[],
        },
        ArticulatedNormalContact {
            point: LinkContactPoint {
                system_index: 0,
                link_index: 0,
                point_body: arm,
            },
            other: ContactCounterpart::PrescribedPlane {
                point_world: pose.transform_point(arm).unwrap(),
                velocity_world: rotation.rotate(linear).unwrap(),
            },
            normal_world: rotation.rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            activation_tolerance: 1e-8,
            restitution: 0.0,
        },
    )
}

/// G1/G3: W = 1/m + r²/I = 2, u = -1, e = 0, hence J = 1/2.
/// Both actual endpoint motion and the kinetic/work account are checked.
#[test]
fn off_center_impact_preserves_relative_motion_under_frame_and_boost_changes() {
    let model = impact_body();
    for rotation in [
        So3::identity(),
        So3::exp(So3Tangent {
            angular: Vec3::new(0.2, -0.3, 0.4),
        })
        .unwrap(),
    ] {
        for boost in [0.0, 32.0, -32.0] {
            let (system, contact) = off_center_impact(&model, rotation, boost);
            let result = resolve_articulated_contacts(
                &[system],
                &[contact],
                MultiImpactConfig::default(),
                || false,
            )
            .unwrap();
            for (actual, expected) in [
                (result.normal_impulses[0], 0.5),
                (result.relative_normal_speeds_after[0], 0.0),
                (result.dissipated_energy, 0.25),
                (result.plane_work, boost * 0.5),
                (result.energy_balance_residual, 0.0),
            ] {
                assert!(
                    (actual - expected).abs() < 1e-11,
                    "{actual:e} != {expected:e}"
                );
            }
            let after = result.responses[0].base_after.twist_body;
            assert!((after.linear - Vec3::new(boost + 0.5, 0.0, 0.0)).norm() < 1e-11);
            assert!((after.angular - Vec3::new(0.0, 0.0, 0.5)).norm() < 1e-11);
        }
    }
}

/// The approach is still present in the native angular state. If the required
/// linear endpoint increment cannot be represented, reject that event instead
/// of reporting a zero-impulse stationary contact.
#[test]
fn unrepresentable_boosted_impact_is_not_misclassified_as_rest() {
    let model = impact_body();
    let (system, contact) = off_center_impact(&model, So3::identity(), 2.0_f64.powi(56));
    let result =
        resolve_articulated_contacts(&[system], &[contact], MultiImpactConfig::default(), || {
            false
        });
    assert!(
        matches!(
            result,
            Err(ArticulatedContactError::Unresolved(
                "actual contact velocity failed admission"
            ))
        ),
        "unrepresentable endpoint must refuse: {result:?}"
    );
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
