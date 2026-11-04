//! Extended source field, its spatial tangent and actual source tidal torque.
use fs_mbd::gravity::{
    CentralGravity, GravityError, PointGravityField, QuadrupoleSourceExtent,
    QuadrupoleSourceGravity, RigidGravityField,
};
use fs_mbd::{MassProperties, Pose, UnitQuaternion, Vec3};

const FRAME: u64 = 41;
const G: f64 = 4.0;
const EXTENT: QuadrupoleSourceExtent = QuadrupoleSourceExtent {
    radius_m: 1.0,
    maximum_radius_ratio: 0.3,
};
fn axis(i: usize, value: f64) -> Vec3 {
    let mut a = [0.0; 3];
    a[i] = value;
    Vec3::new(a[0], a[1], a[2])
}
fn values(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn fixture() -> (Pose, MassProperties) {
    let axes = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.7).unwrap();
    let mass =
        MassProperties::from_principal_frame(2.5, Vec3::ZERO, Vec3::new(0.3, 0.4, 0.5), axes)
            .unwrap();
    let pose = Pose::new(
        Vec3::new(0.3, -0.8, 0.6),
        UnitQuaternion::from_axis_angle(Vec3::new(2.0, 1.0, -1.0), 0.4).unwrap(),
    )
    .unwrap();
    (pose, mass)
}

#[test]
fn extended_source_force_potential_tangent_and_spin_share_one_native_interaction() {
    let (pose, mass) = fixture();
    let field = QuadrupoleSourceGravity::new(FRAME, pose, mass, G, EXTENT).unwrap();
    for point in [Vec3::new(5.0, 3.0, -2.0), Vec3::new(-4.0, 2.0, 5.0)] {
        let sample = field.sample_point(FRAME, point).unwrap();
        // Existing finite-body owner, with the point target supplying unit mass,
        // supplies an independent native reciprocal force and source spin load.
        let source_load = CentralGravity::new(FRAME, point, G)
            .unwrap()
            .quadrupole(EXTENT.radius_m, EXTENT.maximum_radius_ratio)
            .unwrap()
            .load(FRAME, pose, mass)
            .unwrap();
        close(
            sample.specific_potential_j_kg,
            source_load.potential_energy_j,
            2e-15,
        );
        for i in 0..3 {
            close(
                values(sample.acceleration_world_m_s2)[i],
                -values(source_load.wrench.force_world)[i],
                2e-15,
            );
            let h = 1e-4;
            let plus = field.sample_point(FRAME, point.add(axis(i, h))).unwrap();
            let minus = field.sample_point(FRAME, point.sub(axis(i, h))).unwrap();
            close(
                -(plus.specific_potential_j_kg - minus.specific_potential_j_kg) / (2.0 * h),
                values(sample.acceleration_world_m_s2)[i],
                2e-10,
            );
            let derivative = values(
                plus.acceleration_world_m_s2
                    .sub(minus.acceleration_world_m_s2)
                    .scale(0.5 / h),
            );
            for j in 0..3 {
                close(
                    derivative[j],
                    sample.acceleration_gradient_s2.m[3 * j + i],
                    2e-10,
                );
            }
            let h = 1e-5;
            let energy = |angle: f64| {
                let rotated = Pose::new(
                    pose.position_world(),
                    pose.orientation().right_exp(axis(i, angle)).unwrap(),
                )
                .unwrap();
                QuadrupoleSourceGravity::new(FRAME, rotated, mass, G, EXTENT)
                    .unwrap()
                    .sample_point(FRAME, point)
                    .unwrap()
                    .specific_potential_j_kg
            };
            close(
                -(energy(h) - energy(-h)) / (2.0 * h),
                values(source_load.wrench.torque_body)[i],
                3e-11,
            );
            let reciprocal_spin = pose.orientation().rotate_world_to_body(
                point
                    .sub(pose.position_world())
                    .cross(sample.acceleration_world_m_s2)
                    .scale(-1.0),
            );
            close(
                -(energy(h) - energy(-h)) / (2.0 * h),
                values(reciprocal_spin)[i],
                3e-11,
            );
        }
        let m = sample.acceleration_gradient_s2.m;
        assert_eq!((m[1], m[2], m[5]), (m[3], m[6], m[7]));
        close(m[0] + m[4] + m[8], 0.0, 2e-16);
        let orbital = point
            .sub(pose.position_world())
            .cross(sample.acceleration_world_m_s2);
        let spin = pose
            .orientation()
            .rotate_body_to_world(source_load.wrench.torque_body);
        assert!(spin.x.hypot(spin.y).hypot(spin.z) > 1e-4);
        for value in values(orbital.add(spin)) {
            close(value, 0.0, 3e-16);
        }
    }
}

#[test]
fn isotropic_inertia_reduces_exactly_and_source_translation_preserves_relative_field() {
    let (pose, _) = fixture();
    let mass = MassProperties::new(2.5, Vec3::ZERO, Vec3::new(0.4, 0.4, 0.4)).unwrap();
    let field = QuadrupoleSourceGravity::new(FRAME, pose, mass, G, EXTENT).unwrap();
    let central = CentralGravity::new(FRAME, pose.position_world(), G * mass.mass()).unwrap();
    let point = Vec3::new(5.0, 3.0, -2.0);
    assert_eq!(
        field.sample_point(FRAME, point).unwrap(),
        central.sample_point(FRAME, point).unwrap()
    );
    let (pose, mass) = fixture();
    let shift = Vec3::new(10.0, -3.0, 7.0);
    let shifted = Pose::new(pose.position_world().add(shift), pose.orientation()).unwrap();
    let a = QuadrupoleSourceGravity::new(FRAME, pose, mass, G, EXTENT)
        .unwrap()
        .sample_point(FRAME, point)
        .unwrap();
    let b = QuadrupoleSourceGravity::new(FRAME, shifted, mass, G, EXTENT)
        .unwrap()
        .sample_point(FRAME, point.add(shift))
        .unwrap();
    close(a.specific_potential_j_kg, b.specific_potential_j_kg, 2e-15);
    for i in 0..3 {
        close(
            values(a.acceleration_world_m_s2)[i],
            values(b.acceleration_world_m_s2)[i],
            3e-16,
        );
    }
    for i in 0..9 {
        close(
            a.acceleration_gradient_s2.m[i],
            b.acceleration_gradient_s2.m[i],
            3e-16,
        );
    }
}

#[test]
fn extended_source_refuses_wrong_frame_incompatible_extent_and_near_field() {
    let (pose, mass) = fixture();
    let field = QuadrupoleSourceGravity::new(FRAME, pose, mass, G, EXTENT).unwrap();
    assert_eq!(
        field.sample_point(FRAME + 1, Vec3::new(5.0, 3.0, -2.0)),
        Err(GravityError::FrameMismatch)
    );
    assert_eq!(
        field.sample_point(FRAME, pose.position_world().add(Vec3::new(2.0, 0.0, 0.0))),
        Err(GravityError::OutsideExpansionDomain)
    );
    assert!(field.sample_point(FRAME, pose.position_world()).is_err());
    for extent in [
        QuadrupoleSourceExtent {
            radius_m: 0.001,
            ..EXTENT
        },
        QuadrupoleSourceExtent {
            maximum_radius_ratio: 1.0,
            ..EXTENT
        },
        QuadrupoleSourceExtent {
            radius_m: f64::INFINITY,
            ..EXTENT
        },
    ] {
        assert!(QuadrupoleSourceGravity::new(FRAME, pose, mass, G, extent).is_err());
    }
    assert!(QuadrupoleSourceGravity::new(FRAME, pose, mass, 0.0, EXTENT).is_err());
}
