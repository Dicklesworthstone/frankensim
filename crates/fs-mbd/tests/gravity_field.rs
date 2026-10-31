//! G1/G3: the spatial tangent differentiates the same native gravity potential.
use fs_ga::Mat3;
use fs_mbd::gravity::{
    CentralGravity, GravityError, PointGravityField, PointGravitySample, RigidGravityField,
    UniformGravity,
};
use fs_mbd::{Gravity, MassProperties, Pose, UnitQuaternion, Vec3};

const FRAME: u64 = 29;

fn axis(i: usize, value: f64) -> Vec3 {
    let mut v = [0.0; 3];
    v[i] = value;
    Vec3::new(v[0], v[1], v[2])
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

fn apply(matrix: Mat3, direction: Vec3) -> Vec3 {
    let v = matrix.apply(fs_ga::Vec3::new(direction.x, direction.y, direction.z));
    Vec3::new(v.x, v.y, v.z)
}

#[test]
fn central_acceleration_gradient_matches_independent_spatial_differences() {
    let field = CentralGravity::new(FRAME, Vec3::new(0.4, -0.8, 1.1), 7.3).unwrap();
    for point in [Vec3::new(2.0, 1.3, -0.7), Vec3::new(-1.0, 2.4, 3.2)] {
        let sample = field.sample_point(FRAME, point).unwrap();
        assert_eq!(
            (
                sample.acceleration_world_m_s2,
                sample.specific_potential_j_kg
            ),
            field.sample(FRAME, point).unwrap()
        );
        let h = 1e-5;
        for j in 0..3 {
            let plus = field.sample(FRAME, point.add(axis(j, h))).unwrap();
            let minus = field.sample(FRAME, point.sub(axis(j, h))).unwrap();
            close(
                -(plus.1 - minus.1) / (2.0 * h),
                values(sample.acceleration_world_m_s2)[j],
                5e-10,
            );
            let derivative = values(plus.0.sub(minus.0).scale(0.5 / h));
            for i in 0..3 {
                close(
                    derivative[i],
                    sample.acceleration_gradient_s2.m[3 * i + j],
                    5e-10,
                );
            }
        }
        let m = sample.acceleration_gradient_s2.m;
        close(m[0] + m[4] + m[8], 0.0, 2e-15);
    }
}

#[test]
fn uniform_point_field_retains_potential_datum_and_native_mass_ownership() {
    let g = Gravity::new(Vec3::new(0.2, -0.7, -4.1)).unwrap();
    let datum = Vec3::new(1.2, -0.3, 2.4);
    let field = UniformGravity::new(FRAME, g, datum).unwrap();
    let point = Vec3::new(2.3, -1.0, 1.7);
    let sample = field.sample_point(FRAME, point).unwrap();
    assert_eq!(sample.acceleration_world_m_s2, g.acceleration_world());
    assert_eq!(sample.acceleration_gradient_s2.m, [0.0; 9]);
    assert_eq!(
        field
            .sample_point(FRAME, datum)
            .unwrap()
            .specific_potential_j_kg,
        0.0
    );
    let mass = MassProperties::new(3.5, Vec3::ZERO, Vec3::new(1.0, 1.2, 1.4)).unwrap();
    let load = field
        .load(
            FRAME,
            Pose::new(point, UnitQuaternion::IDENTITY).unwrap(),
            mass,
        )
        .unwrap();
    close(
        load.potential_energy_j,
        mass.mass() * sample.specific_potential_j_kg,
        1e-14,
    );
    assert_eq!(
        load.wrench.force_world,
        sample.acceleration_world_m_s2.scale(mass.mass())
    );
}

#[test]
fn central_spatial_derivative_rotates_with_the_declared_world_axes() {
    let center = Vec3::new(1.0, 2.0, -1.0);
    let point = Vec3::new(3.0, -1.0, 0.4);
    let shift = Vec3::new(-5.0, 3.0, 1.7);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.5), 0.7).unwrap();
    let field = CentralGravity::new(FRAME, center, 13.0).unwrap();
    let transformed = CentralGravity::new(
        FRAME + 1,
        rotation.rotate_body_to_world(center).add(shift),
        13.0,
    )
    .unwrap();
    let a = field.sample_point(FRAME, point).unwrap();
    let b = transformed
        .sample_point(FRAME + 1, rotation.rotate_body_to_world(point).add(shift))
        .unwrap();
    close(a.specific_potential_j_kg, b.specific_potential_j_kg, 3e-15);
    for i in 0..3 {
        let direction = axis(i, 1.0);
        let expected = rotation.rotate_body_to_world(apply(a.acceleration_gradient_s2, direction));
        let actual = apply(
            b.acceleration_gradient_s2,
            rotation.rotate_body_to_world(direction),
        );
        for (x, y) in values(expected).into_iter().zip(values(actual)) {
            close(x, y, 3e-15);
        }
    }
}

#[test]
fn central_gradient_avoids_radius_cubing_and_refuses_unrepresentable_tangents() {
    for (mu, radius, scale) in [(1e300, 1e150, 1e-150), (1e-300, 1e-150, 1e150)] {
        let field = CentralGravity::new(FRAME, Vec3::ZERO, mu).unwrap();
        let sample = field
            .sample_point(FRAME, Vec3::new(radius, 0.0, 0.0))
            .unwrap();
        close(sample.acceleration_gradient_s2.m[0] / scale, 2.0, 2e-15);
        close(sample.acceleration_gradient_s2.m[4] / scale, -1.0, 2e-15);
        close(sample.acceleration_gradient_s2.m[8] / scale, -1.0, 2e-15);
    }
    let field = CentralGravity::new(FRAME, Vec3::ZERO, 1e-300).unwrap();
    let point = Vec3::new(1e-210, 0.0, 0.0);
    assert!(field.sample(FRAME, point).is_ok());
    assert!(field.sample_point(FRAME, point).is_err());
    assert_eq!(
        field.sample_point(FRAME + 1, point),
        Err(GravityError::FrameMismatch)
    );
    assert!(field.sample_point(FRAME, Vec3::ZERO).is_err());
}

#[test]
fn point_field_admission_rejects_nonfinite_and_nonconservative_derivatives() {
    let mut sample = PointGravitySample {
        acceleration_world_m_s2: Vec3::ZERO,
        specific_potential_j_kg: 0.0,
        acceleration_gradient_s2: Mat3 { m: [0.0; 9] },
    };
    sample.acceleration_gradient_s2.m[1] = 1.0;
    assert!(sample.validate().is_err());
    sample.acceleration_gradient_s2.m[3] = 1.0;
    assert!(sample.validate().is_ok());
    sample.acceleration_world_m_s2.x = f64::NAN;
    assert!(sample.validate().is_err());
}
