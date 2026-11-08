//! G1/G3: exact radial source layers through the native gravity field owners.
use fs_mbd::gravity::{
    CentralGravity, GravityError, PointGravityField, RadialGravity, RadialGravityLayer,
    RigidGravityField, UniformSphereGravity,
};
use fs_mbd::{MassProperties, Pose, UnitQuaternion, Vec3};

const FRAME: u64 = 109;

fn layer(radius: f64, mass: f64) -> RadialGravityLayer {
    RadialGravityLayer {
        outer_radius_m: radius,
        mass_kg: mass,
    }
}

fn field(layers: &[RadialGravityLayer]) -> RadialGravity {
    RadialGravity::new(FRAME, Vec3::ZERO, 1.0, layers, layers.len()).unwrap()
}

fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn uniform_density_partition_recovers_the_exact_finite_sphere() {
    // Each shell's mass equals its difference of cubed radii: G*rho*4pi/3 = 1.
    let radial = field(&[layer(0.5, 0.125), layer(1.0, 0.875), layer(2.0, 7.0)]);
    let sphere = UniformSphereGravity::new(FRAME, Vec3::ZERO, 8.0, 2.0).unwrap();
    for r in [0.0, 0.1, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 100.0] {
        let point = Vec3::new(0.6 * r, 0.8 * r, 0.0);
        let actual = radial.sample_point(FRAME, point).unwrap();
        let expected = sphere.sample_point(FRAME, point).unwrap();
        assert!(
            norm(
                actual
                    .acceleration_world_m_s2
                    .sub(expected.acceleration_world_m_s2)
            ) < 3e-14
        );
        close(
            actual.specific_potential_j_kg,
            expected.specific_potential_j_kg,
            3e-14,
        );
        for (a, b) in actual
            .acceleration_gradient_s2
            .m
            .into_iter()
            .zip(expected.acceleration_gradient_s2.m)
        {
            close(a, b, 4e-14);
        }
    }
}

#[test]
fn hollow_shell_and_dense_core_follow_enclosed_mass_and_poisson_equation() {
    let hollow = field(&[layer(1.0, 0.0), layer(2.0, 7.0)]);
    for point in [
        Vec3::ZERO,
        Vec3::new(0.3, -0.4, 0.1),
        Vec3::new(0.99, 0.0, 0.0),
    ] {
        let sample = hollow.sample_point(FRAME, point).unwrap();
        assert_eq!(sample.acceleration_world_m_s2, Vec3::ZERO);
        assert_eq!(sample.acceleration_gradient_s2.m, [0.0; 9]);
        close(sample.specific_potential_j_kg, -4.5, 2e-14);
    }
    assert_eq!(
        hollow
            .specific_potential_increment(FRAME, Vec3::ZERO, Vec3::new(0.8, 0.0, 0.0))
            .unwrap(),
        Some(0.0)
    );
    let shell = hollow
        .sample_point(FRAME, Vec3::new(1.5, 0.0, 0.0))
        .unwrap();
    let enclosed = 1.5_f64.powi(3) - 1.0;
    close(
        shell.acceleration_world_m_s2.x,
        -enclosed / 1.5_f64.powi(2),
        2e-14,
    );
    close(
        shell.specific_potential_j_kg,
        -enclosed / 1.5 - 1.5 * (4.0 - 2.25),
        2e-14,
    );
    close(
        shell.acceleration_gradient_s2.m[0]
            + shell.acceleration_gradient_s2.m[4]
            + shell.acceleration_gradient_s2.m[8],
        -3.0,
        3e-14,
    );

    let dense = field(&[layer(1.0, 8.0), layer(2.0, 7.0)]);
    let core = dense.sample_point(FRAME, Vec3::new(0.5, 0.0, 0.0)).unwrap();
    close(core.acceleration_world_m_s2.x, -4.0, 2e-14);
    close(core.specific_potential_j_kg, -15.5, 3e-14);
    assert_eq!(
        core.acceleration_gradient_s2.m,
        [-8.0, 0.0, 0.0, 0.0, -8.0, 0.0, 0.0, 0.0, -8.0]
    );
    let interface = dense.sample_point(FRAME, Vec3::new(1.0, 0.0, 0.0)).unwrap();
    close(interface.acceleration_world_m_s2.x, -8.0, 3e-14);
    close(interface.specific_potential_j_kg, -12.5, 3e-14);
    close(interface.acceleration_gradient_s2.m[0], 13.0, 4e-14);
    close(interface.acceleration_gradient_s2.m[4], -8.0, 4e-14);
    let central = CentralGravity::new(FRAME, Vec3::ZERO, 15.0).unwrap();
    assert_eq!(
        dense
            .sample_point(FRAME, Vec3::new(4.0, -1.0, 2.0))
            .unwrap(),
        central
            .sample_point(FRAME, Vec3::new(4.0, -1.0, 2.0))
            .unwrap()
    );
}

#[test]
fn layered_force_and_tangent_differentiate_one_potential() {
    let model = field(&[layer(1.0, 8.0), layer(2.0, 7.0)]);
    let h = 2e-5;
    for p in [
        Vec3::new(0.3, -0.2, 0.4),
        Vec3::new(1.2, 0.4, -0.2),
        Vec3::new(2.5, -0.7, 0.9),
    ] {
        let sample = model.sample_point(FRAME, p).unwrap();
        let acceleration = [
            sample.acceleration_world_m_s2.x,
            sample.acceleration_world_m_s2.y,
            sample.acceleration_world_m_s2.z,
        ];
        for axis in 0..3 {
            let mut d = [0.0; 3];
            d[axis] = h;
            let d = Vec3::new(d[0], d[1], d[2]);
            let plus = model.sample_point(FRAME, p.add(d)).unwrap();
            let minus = model.sample_point(FRAME, p.sub(d)).unwrap();
            close(
                -(plus.specific_potential_j_kg - minus.specific_potential_j_kg) / (2.0 * h),
                acceleration[axis],
                2e-8,
            );
            let derivative = plus
                .acceleration_world_m_s2
                .sub(minus.acceleration_world_m_s2)
                .scale(0.5 / h);
            for (row, value) in [derivative.x, derivative.y, derivative.z]
                .into_iter()
                .enumerate()
            {
                close(
                    value,
                    sample.acceleration_gradient_s2.m[3 * row + axis],
                    3e-8,
                );
            }
        }
    }
}

#[test]
fn direct_work_survives_dominant_outer_shell_potential_and_native_mass_scaling() {
    let model = field(&[layer(1.0, 1.0), layer(2.0, 1e20)]);
    let transverse = 2.0_f64.powi(-27);
    let before = Vec3::new(0.5, 0.0, 0.0);
    let after = Vec3::new(0.5, transverse, 0.0);
    let expected = 0.5 * transverse * transverse;
    assert_eq!(
        model.sample(FRAME, before).unwrap().1,
        model.sample(FRAME, after).unwrap().1
    );
    close(
        model
            .specific_potential_increment(FRAME, before, after)
            .unwrap()
            .unwrap(),
        expected,
        expected * 2e-14,
    );
    let properties = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(0.1, 0.1, 0.1)).unwrap();
    let pose = |point| Pose::new(point, UnitQuaternion::IDENTITY).unwrap();
    close(
        model
            .potential_energy_increment(FRAME, pose(before), pose(after), properties)
            .unwrap()
            .unwrap(),
        3.0 * expected,
        expected * 8e-14,
    );
    let load = model.load(FRAME, pose(before), properties).unwrap();
    close(load.wrench.force_world.x, -1.5, 1e-14);
    assert_eq!(load.wrench.torque_body, Vec3::ZERO);
}

#[test]
fn work_crosses_density_interfaces_continuously_and_reverses() {
    let model = field(&[layer(1.0, 8.0), layer(2.0, 7.0)]);
    let point = |r| Vec3::new(r, 0.0, 0.0);
    let delta = |a, b| {
        model
            .specific_potential_increment(FRAME, point(a), point(b))
            .unwrap()
            .unwrap()
    };
    let all = delta(0.4, 4.0);
    close(all, 12.11, 4e-14);
    close(delta(4.0, 0.4), -all, 2e-14);
    let partition = [0.4, 1.0, 1.25, 2.0, 4.0];
    let sum = partition
        .windows(2)
        .map(|pair| delta(pair[0], pair[1]))
        .sum::<f64>();
    close(sum, all, 5e-14);
    let e = 2.0_f64.powi(-40);
    let expected =
        (8.0 * e - 4.0 * e * e) + 8.0 * e / (1.0 + e) + 0.5 * e * e * (3.0 + e) / (1.0 + e);
    close(delta(1.0 - e, 1.0 + e), expected, 5e-26);
}

#[test]
fn layer_budget_source_data_and_query_frame_are_admitted() {
    let layers = [layer(1.0, 1.0), layer(2.0, 7.0)];
    assert!(RadialGravity::new(FRAME, Vec3::ZERO, 1.0, &layers, 1).is_err());
    for invalid in [
        [layer(1.0, 1.0), layer(1.0, 1.0)],
        [layer(0.0, 1.0), layer(2.0, 1.0)],
        [layer(1.0, -1.0), layer(2.0, 1.0)],
        [layer(1.0, 1.0), layer(f64::NAN, 1.0)],
    ] {
        assert!(RadialGravity::new(FRAME, Vec3::ZERO, 1.0, &invalid, 2).is_err());
    }
    let model = field(&layers);
    assert!(matches!(
        model.sample_point(FRAME + 1, Vec3::ZERO),
        Err(GravityError::FrameMismatch)
    ));
    assert!(
        model
            .sample_point(FRAME, Vec3::new(f64::NAN, 0.0, 0.0))
            .is_err()
    );
    assert!(RadialGravity::new(FRAME, Vec3::ZERO, 1.0, &[], 0).is_err());
    let empty_mass = field(&[layer(1.0, 0.0), layer(2.0, 0.0)]);
    for radius in [0.0, 0.5, 1.0, 2.0, 10.0] {
        let sample = empty_mass
            .sample_point(FRAME, Vec3::new(radius, 0.0, 0.0))
            .unwrap();
        assert_eq!(sample.acceleration_world_m_s2, Vec3::ZERO);
        assert_eq!(sample.specific_potential_j_kg, 0.0);
        assert_eq!(sample.acceleration_gradient_s2.m, [0.0; 9]);
    }
}

#[test]
fn thin_shell_mass_fraction_is_not_lost_to_cubed_radius_cancellation() {
    let a: f64 = 1.0;
    let width = 2.0_f64.powi(-40);
    let b = a + width;
    let r = a + 0.25 * width;
    let model = field(&[layer(a, 0.0), layer(b, 1.0)]);
    let enclosed = ((r - a) / (b - a)) * ((r * r + r * a + a * a) / (b * b + b * a + a * a));
    let sample = model.sample_point(FRAME, Vec3::new(r, 0.0, 0.0)).unwrap();
    close(sample.acceleration_world_m_s2.x, -enclosed / (r * r), 2e-14);
    let s = a / b;
    let central_potential = -1.5 / b * (1.0 + s) / (1.0 + s + s * s);
    close(
        model.sample(FRAME, Vec3::ZERO).unwrap().1,
        central_potential,
        2e-14,
    );
    let k = 1.0 / ((b - a) * (b * b + b * a + a * a));
    close(
        sample.acceleration_gradient_s2.m[0],
        2.0 * enclosed / (r * r * r) - 3.0 * k,
        0.02,
    );
}
