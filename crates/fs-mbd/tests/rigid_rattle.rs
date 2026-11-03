use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::rigid_rattle::{PointPlaneRattle, PointPlaneRattleControls};
use fs_mbd::{MassProperties, Vec3 as MbdVec3};
use fs_time::se3::RattleProjection;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + actual.abs().max(expected.abs())),
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn full_tensor_point_plane_projection_uses_the_mass_metric() {
    let inertia = Mat3 {
        m: [2.0, 0.4, 0.0, 0.4, 1.6, 0.2, 0.0, 0.2, 1.2],
    };
    let properties = MassProperties::from_inertia_tensor(3.0, MbdVec3::ZERO, inertia).unwrap();
    let projection = PointPlaneRattle::new(
        properties,
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        PointPlaneRattleControls {
            position_tolerance_m: 2e-13,
            velocity_tolerance_mps: 2e-13,
            maximum_position_iterations: 16,
        },
    )
    .unwrap();

    let mut pose = Se3::from_parts(So3::identity(), Vec3::new(0.2, 0.0, 0.0)).unwrap();
    let initial_rotation = pose.rotation().matrix();
    let removed = projection.project_position(&mut pose).unwrap();
    close(removed, 0.2, 2e-12);
    assert!(projection.gap(pose).unwrap().abs() <= 2e-13);
    assert!((pose.translation().x - 0.2).abs() > 1e-4);
    assert!(
        pose.rotation()
            .matrix()
            .m
            .iter()
            .zip(initial_rotation.m)
            .any(|(a, b)| (a - b).abs() > 1e-5)
    );

    let mut twist = Twist::new(Vec3::new(0.2, -0.3, 0.7), Vec3::new(1.1, 0.4, -0.2));
    let before = twist;
    let speed = projection.normal_speed(pose, twist).unwrap();
    let removed_speed = projection.project_velocity(&pose, &mut twist).unwrap();
    close(removed_speed, speed.abs(), 2e-12);
    assert!(projection.normal_speed(pose, twist).unwrap().abs() <= 2e-13);

    let delta_w = twist.angular - before.angular;
    let delta_v = twist.linear - before.linear;
    let delta_l = properties
        .angular_momentum_body_checked(MbdVec3::new(delta_w.x, delta_w.y, delta_w.z))
        .unwrap();
    let n_body = pose
        .rotation()
        .inverse()
        .rotate(Vec3::new(1.0, 0.0, 0.0))
        .unwrap();
    let angular_row = Vec3::new(0.0, 1.0, 0.0).cross(n_body);
    let lambda = 3.0 * delta_v.dot(n_body);
    close(delta_l.x, lambda * angular_row.x, 2e-11);
    close(delta_l.y, lambda * angular_row.y, 2e-11);
    close(delta_l.z, lambda * angular_row.z, 2e-11);
    close(3.0 * delta_v.x, lambda * n_body.x, 2e-11);
    close(3.0 * delta_v.y, lambda * n_body.y, 2e-11);
    close(3.0 * delta_v.z, lambda * n_body.z, 2e-11);
    assert!(
        delta_w.x.abs() + delta_w.y.abs() > 1e-4,
        "products of inertia must couple the angular correction"
    );
}

#[test]
fn candidate_receipt_and_invalid_geometry_are_explicit() {
    let properties = MassProperties::new(2.0, MbdVec3::ZERO, MbdVec3::new(1.0, 1.2, 1.4)).unwrap();
    let projection = PointPlaneRattle::new(
        properties,
        Vec3::new(0.0, 0.5, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(2.0, 0.0, 0.0),
        PointPlaneRattleControls::default(),
    )
    .unwrap();
    close(
        projection.normal_mobility(Se3::identity()).unwrap(),
        0.5 + 0.25 / 1.4,
        2e-12,
    );
    let mut pose = Se3::from_parts(So3::identity(), Vec3::new(-0.03, 0.0, 0.0)).unwrap();
    let mut twist = Twist::new(Vec3::new(0.0, 0.0, 0.4), Vec3::new(-0.7, 0.0, 0.0));
    let receipt = projection.project_candidate(&mut pose, &mut twist).unwrap();
    assert!(receipt.position_violation_removed_m > 0.0);
    assert!(receipt.velocity_violation_removed_mps > 0.0);
    assert!(receipt.final_gap_m.abs() <= 1e-12);
    assert!(receipt.final_normal_speed_mps.abs() <= 1e-12);

    assert!(
        PointPlaneRattle::new(
            properties,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            PointPlaneRattleControls::default(),
        )
        .is_err()
    );
}

fn vector_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}

// G0/G3: a fixed plane cannot exert a tangential force or a moment about
// its normal. A configuration repair must not rotate those conserved momenta
// simply because a velocity was stored in the previous body's coordinates.
#[test]
fn candidate_preserves_plane_symmetries_and_applies_only_the_normal_impulse() {
    let inertia = Mat3 {
        m: [2.0, 0.4, 0.0, 0.4, 1.6, 0.2, 0.0, 0.2, 1.2],
    };
    let mass = 3.0;
    let properties = MassProperties::from_inertia_tensor(mass, MbdVec3::ZERO, inertia).unwrap();
    let point_body = Vec3::new(0.2, 0.5, -0.3);
    let original_pose = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.3, -0.4, 0.2))).unwrap(),
        Vec3::new(1.2, -0.5, 0.7),
    )
    .unwrap();
    let original_twist = Twist::new(Vec3::new(0.8, -0.4, 0.3), Vec3::new(1.3, -2.0, 0.7));
    let normal = Vec3::new(0.3, -0.4, 0.8).scale(1.0 / 0.89_f64.sqrt());
    let plane_point = original_pose.transform_point(point_body).unwrap() - normal.scale(0.15);
    let mut reference = None;
    for world_rotation in [
        So3::identity(),
        So3::exp(So3Tangent::new(Vec3::new(-0.7, 0.2, 0.4))).unwrap(),
    ] {
        let normal = world_rotation.rotate(normal).unwrap();
        let pose_before = Se3::from_parts(
            world_rotation.compose(original_pose.rotation()).unwrap(),
            world_rotation.rotate(original_pose.translation()).unwrap(),
        )
        .unwrap();
        let projection = PointPlaneRattle::new(
            properties,
            point_body,
            world_rotation.rotate(plane_point).unwrap(),
            normal,
            PointPlaneRattleControls::default(),
        )
        .unwrap();
        let p_before = pose_before
            .rotation()
            .rotate(original_twist.linear.scale(mass))
            .unwrap();
        let spin_before = pose_before
            .rotation()
            .rotate(inertia.apply(original_twist.angular))
            .unwrap();
        let mut pose = pose_before;
        let mut twist = original_twist;
        let receipt = projection.project_candidate(&mut pose, &mut twist).unwrap();
        assert!(receipt.position_violation_removed_m > 0.1);
        assert!(receipt.final_gap_m.abs() <= 1e-12);
        assert!(receipt.final_normal_speed_mps.abs() <= 1e-12);
        let p_after = pose.rotation().rotate(twist.linear.scale(mass)).unwrap();
        let spin_after = pose
            .rotation()
            .rotate(inertia.apply(twist.angular))
            .unwrap();
        let impulse = p_after - p_before;
        vector_close(impulse.cross(normal), Vec3::new(0.0, 0.0, 0.0), 2e-11);
        vector_close(
            spin_after - spin_before,
            pose.rotation().rotate(point_body).unwrap().cross(impulse),
            2e-11,
        );
        // The COM's position correction is itself normal, so its orbital
        // contribution cannot change angular momentum about the normal.
        vector_close(
            (pose.translation() - pose_before.translation()).cross(normal),
            Vec3::new(0.0, 0.0, 0.0),
            2e-12,
        );
        close(
            (spin_after + pose.translation().cross(p_after)).dot(normal),
            (spin_before + pose_before.translation().cross(p_before)).dot(normal),
            2e-11,
        );
        if let Some((reference_pose, reference_twist)) = reference {
            let reference_pose: Se3 = reference_pose;
            let reference_twist: Twist = reference_twist;
            vector_close(
                pose.translation(),
                world_rotation.rotate(reference_pose.translation()).unwrap(),
                2e-11,
            );
            let expected = world_rotation
                .compose(reference_pose.rotation())
                .unwrap()
                .matrix();
            for (actual, expected) in pose.rotation().matrix().m.into_iter().zip(expected.m) {
                close(actual, expected, 2e-11);
            }
            vector_close(twist.angular, reference_twist.angular, 2e-11);
            vector_close(twist.linear, reference_twist.linear, 2e-11);
        } else {
            reference = Some((pose, twist));
        }
    }
}

#[test]
fn center_of_mass_position_repair_does_not_manufacture_spin() {
    let properties = MassProperties::new(2.0, MbdVec3::ZERO, MbdVec3::new(1.0, 1.2, 1.4)).unwrap();
    let projection = PointPlaneRattle::new(
        properties,
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        PointPlaneRattleControls::default(),
    )
    .unwrap();
    let mut pose = Se3::from_parts(So3::identity(), Vec3::new(0.2, -0.4, 0.7)).unwrap();
    let mut twist = Twist::new(Vec3::new(0.2, -0.3, 0.4), Vec3::new(0.8, 1.3, -2.1));
    let spin_before = twist.angular;
    projection.project_candidate(&mut pose, &mut twist).unwrap();
    vector_close(pose.translation(), Vec3::new(0.0, -0.4, 0.7), 2e-12);
    vector_close(twist.angular, spin_before, 2e-12);
    vector_close(twist.linear, Vec3::new(0.0, 1.3, -2.1), 2e-12);
}

// G4: neither a failed configuration solve nor a late invalid twist may
// publish the earlier tentative position correction.
#[test]
fn refused_projection_leaves_original_configuration_and_twist_unchanged() {
    let properties = MassProperties::new(2.0, MbdVec3::ZERO, MbdVec3::new(1.0, 1.2, 1.4)).unwrap();
    let projection = PointPlaneRattle::new(
        properties,
        Vec3::new(0.0, 0.5, 0.2),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        PointPlaneRattleControls {
            position_tolerance_m: 1e-14,
            maximum_position_iterations: 1,
            ..PointPlaneRattleControls::default()
        },
    )
    .unwrap();
    let original = Se3::from_parts(So3::identity(), Vec3::new(0.4, -0.1, 0.3)).unwrap();
    let mut pose = original;
    assert!(projection.project_position(&mut pose).is_err());
    assert_eq!(pose, original);

    let projection = PointPlaneRattle::new(
        properties,
        Vec3::new(0.0, 0.5, 0.2),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        PointPlaneRattleControls::default(),
    )
    .unwrap();
    let mut twist = Twist::new(
        Vec3::new(f64::INFINITY, 0.0, 0.0),
        Vec3::new(0.2, 0.4, -0.3),
    );
    let twist_before = twist;
    assert!(projection.project_candidate(&mut pose, &mut twist).is_err());
    assert_eq!(pose, original);
    assert_eq!(twist, twist_before);

    // The finite candidate misses a deliberately zero velocity tolerance by
    // roundoff. It must not replace the accepted twist on that refusal either.
    let projection = PointPlaneRattle::new(
        properties,
        Vec3::new(0.2, 0.5, -0.3),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 2.0, 3.0),
        PointPlaneRattleControls {
            velocity_tolerance_mps: 0.0,
            ..PointPlaneRattleControls::default()
        },
    )
    .unwrap();
    let mut twist = Twist::new(Vec3::new(0.8, -0.4, 0.3), Vec3::new(1.3, -2.0, 0.7));
    let twist_before = twist;
    assert!(
        projection
            .project_velocity(&Se3::identity(), &mut twist)
            .is_err()
    );
    assert_eq!(twist, twist_before);
}
