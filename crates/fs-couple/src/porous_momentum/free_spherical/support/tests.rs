use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} within {tolerance:e}"
    );
}
fn close_v(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}
fn law() -> PrincipalPorousResistance {
    PrincipalPorousResistance {
        darcy_principal: [0.1, 0.4, 0.9],
        forchheimer: 0.7,
    }
}
fn ellipsoid() -> PorousBodySupport {
    PorousBodySupport::Ellipsoid {
        semiaxes_m: [2.0, 0.7, 0.3],
    }
}
fn drag(axes: [[f64; 3]; 3], law: PrincipalPorousResistance, slip: Vec3) -> Vec3 {
    let mut value = slip.scale(law.forchheimer * norm(slip));
    for (axis, coefficient) in axes.into_iter().zip(law.darcy_principal) {
        let direction = Vec3::new(axis[0], axis[1], axis[2]);
        value = value.add(direction.scale(coefficient * direction.dot(slip)));
    }
    value
}

#[test]
fn oriented_ellipsoid_uses_material_support_not_its_enclosing_sphere() {
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.71).unwrap();
    let support = PreparedSupport::new(ellipsoid(), frame);
    for point in [
        Vec3::ZERO,
        Vec3::new(1.5, 0.2, 0.1),
        Vec3::new(0.5, -0.4, -0.15),
    ] {
        assert!(support.contains(frame.rotate_body_to_world(point)));
    }
    for point in [
        Vec3::new(0.0, 0.8, 0.0),
        Vec3::new(0.0, 0.0, 0.4),
        Vec3::new(1.8, 0.4, 0.2),
    ] {
        assert!(norm(point) < ellipsoid().admit().unwrap());
        assert!(!support.contains(frame.rotate_body_to_world(point)));
    }
    let unrotated = PreparedSupport::new(ellipsoid(), UnitQuaternion::IDENTITY);
    assert!(unrotated.contains(Vec3::new(1.5, 0.0, 0.0)));
    assert!(!support.contains(Vec3::new(1.5, 0.0, 0.0)));
}

#[test]
fn material_midpoint_matches_independent_cayley_half_angle_and_work_axis() {
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 0.5), 1.1).unwrap();
    let omega = Vec3::new(0.4, -0.6, 0.7);
    let dt = 0.8;
    let magnitude = norm(omega);
    let half_angle = (0.5 * dt * magnitude).atan();
    let endpoint = frame
        .right_exp(omega.scale(2.0 * half_angle / magnitude))
        .unwrap();
    let expected = frame
        .right_exp(omega.scale(half_angle / magnitude))
        .unwrap();
    let midpoint = midpoint_orientation(frame, endpoint).unwrap();
    for point in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.4, -0.2, 1.0),
    ] {
        close_v(
            midpoint.rotate_body_to_world(point),
            expected.rotate_body_to_world(point),
            8e-16,
        );
    }
    close_v(
        midpoint.rotate_body_to_world(omega),
        frame.rotate_body_to_world(omega),
        6e-16,
    );
    close_v(
        midpoint.rotate_body_to_world(omega),
        endpoint.rotate_body_to_world(omega),
        6e-16,
    );
    assert_eq!(midpoint_orientation(frame, frame).unwrap(), frame);
}

#[test]
fn midpoint_crosses_quaternion_sign_cut_on_the_short_physical_path() {
    let axis = Vec3::new(0.0, 0.0, 1.0);
    let before = UnitQuaternion::from_axis_angle(axis, 170.0_f64.to_radians()).unwrap();
    let after = UnitQuaternion::from_axis_angle(axis, 190.0_f64.to_radians()).unwrap();
    let expected = UnitQuaternion::from_axis_angle(axis, core::f64::consts::PI).unwrap();
    let midpoint = midpoint_orientation(before, after).unwrap();
    let reverse = midpoint_orientation(after, before).unwrap();
    for point in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 1.0)] {
        close_v(
            midpoint.rotate_body_to_world(point),
            expected.rotate_body_to_world(point),
            4e-16,
        );
        close_v(
            midpoint.rotate_body_to_world(point),
            reverse.rotate_body_to_world(point),
            4e-16,
        );
    }
}

#[test]
fn rotated_anisotropic_drag_is_covariant_passive_and_not_a_scalar_substitute() {
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(-1.0, 2.0, 3.0), 0.9).unwrap();
    let prepared = PreparedSupport::new(ellipsoid(), frame);
    let body_slip = Vec3::new(0.3, -0.1, 0.2);
    let world_slip = frame.rotate_body_to_world(body_slip);
    let actual = drag(prepared.resistance_axes(law()), law(), world_slip);
    let expected_body = Vec3::new(
        law().darcy_principal[0] * body_slip.x,
        law().darcy_principal[1] * body_slip.y,
        law().darcy_principal[2] * body_slip.z,
    )
    .add(body_slip.scale(law().forchheimer * norm(body_slip)));
    close_v(actual, frame.rotate_body_to_world(expected_body), 3e-16);
    close(actual.dot(world_slip), expected_body.dot(body_slip), 2e-16);
    assert!(actual.dot(world_slip) > 0.0);
    assert!(norm(actual.cross(world_slip)) > 1e-3);
    for i in 0..3 {
        for j in 0..3 {
            let a = prepared.axes[i];
            let b = prepared.axes[j];
            close(
                Vec3::new(a[0], a[1], a[2]).dot(Vec3::new(b[0], b[1], b[2])),
                f64::from(i == j),
                8e-16,
            );
        }
    }
}

#[test]
fn sampled_anisotropic_exchange_pairs_full_force_torque_and_relative_work() {
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.6).unwrap();
    let support = PreparedSupport::new(ellipsoid(), orientation);
    let velocity = Vec3::new(0.2, -0.1, 0.05);
    let omega = Vec3::new(-0.3, 0.2, 0.6);
    let boost = Vec3::new(0.7, 0.5, -0.2);
    let mut impulse = Vec3::ZERO;
    let mut couple = Vec3::ZERO;
    let mut skeleton_work = 0.0;
    let mut fluid_work = 0.0;
    let mut loss = 0.0;
    for local in [
        Vec3::new(0.3, 0.2, 0.1),
        Vec3::new(-0.7, 0.1, -0.1),
        Vec3::new(0.4, -0.3, 0.15),
    ] {
        let arm = orientation.rotate_body_to_world(local);
        assert!(support.contains(arm));
        let skeleton = velocity.add(omega.cross(arm));
        let fluid = Vec3::new(0.04, 0.03, -0.02).add(local.scale(0.01));
        let slip = fluid.sub(skeleton);
        let j = drag(support.resistance_axes(law()), law(), slip);
        close_v(
            j,
            drag(
                support.resistance_axes(law()),
                law(),
                fluid.add(boost).sub(skeleton.add(boost)),
            ),
            4e-16,
        );
        impulse = impulse.add(j);
        couple = couple.add(arm.cross(j));
        skeleton_work += j.dot(skeleton);
        fluid_work -= j.dot(fluid);
        loss += j.dot(slip);
    }
    close(
        skeleton_work,
        velocity.dot(impulse) + omega.dot(couple),
        3e-16,
    );
    close(fluid_work + skeleton_work + loss, 0.0, 3e-16);
    assert!(loss > 0.0);
    assert!(norm(couple) > 1e-2);
}

#[test]
fn scalar_sphere_and_isotropic_tensor_keep_the_original_native_representation() {
    let shape = PorousBodySupport::Sphere { radius_m: 0.7 };
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 1.4).unwrap();
    let prepared = PreparedSupport::new(shape, frame);
    for point in [
        Vec3::ZERO,
        Vec3::new(0.4, 0.2, 0.1),
        Vec3::new(0.7, 0.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
    ] {
        assert_eq!(prepared.contains(point), norm(point) <= 0.7);
    }
    assert_eq!(
        prepared.resistance_axes(PrincipalPorousResistance::isotropic(0.4, 2.0)),
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
}

#[test]
fn support_and_resistance_refuse_unphysical_dimensions_without_inventing_mass() {
    assert_eq!(ellipsoid().admit().unwrap(), 2.0);
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            PorousBodySupport::Ellipsoid {
                semiaxes_m: [1.0, bad, 0.3]
            }
            .admit()
            .is_err()
        );
        assert!(PorousBodySupport::Sphere { radius_m: bad }.admit().is_err());
    }
    for bad in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            PrincipalPorousResistance {
                darcy_principal: [0.0, bad, 0.5],
                forchheimer: 0.0
            }
            .admit()
            .is_err()
        );
        assert!(
            PrincipalPorousResistance {
                darcy_principal: [0.0; 3],
                forchheimer: bad
            }
            .admit()
            .is_err()
        );
    }
    assert!(
        PrincipalPorousResistance {
            darcy_principal: [0.0, 1.0, 0.0],
            forchheimer: 0.0
        }
        .admit()
        .is_ok()
    );
    let thin = PreparedSupport::new(
        PorousBodySupport::Ellipsoid {
            semiaxes_m: [1.0, 1e-300, 0.3],
        },
        UnitQuaternion::IDENTITY,
    );
    assert!(thin.contains(Vec3::ZERO));
    assert!(!thin.contains(Vec3::new(0.0, 1e300, 0.0)));
}
