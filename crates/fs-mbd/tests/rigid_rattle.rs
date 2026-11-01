use fs_ga::{Mat3, Se3, So3, Twist, Vec3};
use fs_mbd::rigid_rattle::{PointPlaneRattle, PointPlaneRattleControls};
use fs_mbd::{MassProperties, Vec3 as MbdVec3};
use fs_time::se3::RattleProjection;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * (1.0 + actual.abs().max(expected.abs())),
        "{actual:.16e} != {expected:.16e}");
}

#[test]
fn full_tensor_point_plane_projection_uses_the_mass_metric() {
    let inertia = Mat3 { m: [2.0, 0.4, 0.0, 0.4, 1.6, 0.2, 0.0, 0.2, 1.2] };
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
    ).unwrap();

    let mut pose = Se3::from_parts(So3::identity(), Vec3::new(0.2, 0.0, 0.0)).unwrap();
    let initial_rotation = pose.rotation().matrix();
    let removed = projection.project_position(&mut pose).unwrap();
    close(removed, 0.2, 2e-12);
    assert!(projection.gap(pose).unwrap().abs() <= 2e-13);
    assert!((pose.translation().x - 0.2).abs() > 1e-4);
    assert!(pose.rotation().matrix().m.iter().zip(initial_rotation.m).any(|(a, b)| (a - b).abs() > 1e-5));

    let mut twist = Twist::new(Vec3::new(0.2, -0.3, 0.7), Vec3::new(1.1, 0.4, -0.2));
    let before = twist;
    let speed = projection.normal_speed(pose, twist).unwrap();
    let removed_speed = projection.project_velocity(&pose, &mut twist).unwrap();
    close(removed_speed, speed.abs(), 2e-12);
    assert!(projection.normal_speed(pose, twist).unwrap().abs() <= 2e-13);

    let delta_w = twist.angular - before.angular;
    let delta_v = twist.linear - before.linear;
    let delta_l = properties.angular_momentum_body_checked(MbdVec3::new(delta_w.x, delta_w.y, delta_w.z)).unwrap();
    let n_body = pose.rotation().inverse().rotate(Vec3::new(1.0, 0.0, 0.0)).unwrap();
    let angular_row = Vec3::new(0.0, 1.0, 0.0).cross(n_body);
    let lambda = 3.0 * delta_v.dot(n_body);
    close(delta_l.x, lambda * angular_row.x, 2e-11);
    close(delta_l.y, lambda * angular_row.y, 2e-11);
    close(delta_l.z, lambda * angular_row.z, 2e-11);
    close(3.0 * delta_v.x, lambda * n_body.x, 2e-11);
    close(3.0 * delta_v.y, lambda * n_body.y, 2e-11);
    close(3.0 * delta_v.z, lambda * n_body.z, 2e-11);
    assert!(delta_w.x.abs() + delta_w.y.abs() > 1e-4, "products of inertia must couple the angular correction");
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
    ).unwrap();
    close(projection.normal_mobility(Se3::identity()).unwrap(), 0.5 + 0.25 / 1.4, 2e-12);
    let mut pose = Se3::from_parts(So3::identity(), Vec3::new(-0.03, 0.0, 0.0)).unwrap();
    let mut twist = Twist::new(Vec3::new(0.0, 0.0, 0.4), Vec3::new(-0.7, 0.0, 0.0));
    let receipt = projection.project_candidate(&mut pose, &mut twist).unwrap();
    assert!(receipt.position_violation_removed_m > 0.0);
    assert!(receipt.velocity_violation_removed_mps > 0.0);
    assert!(receipt.final_gap_m.abs() <= 1e-12);
    assert!(receipt.final_normal_speed_mps.abs() <= 1e-12);

    assert!(PointPlaneRattle::new(
        properties,
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        PointPlaneRattleControls::default(),
    ).is_err());
}
