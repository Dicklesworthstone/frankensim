//! G0/G1/G3: vorticity and velocity gradients of the native planar-wall field.

use core::f64::consts::{PI, TAU};
use fs_vpm::planar_wall::PlanarWall2;
use fs_vpm::{ParticleKinematics2, VortexParticle, induced_kinematics};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn sources() -> [VortexParticle; 3] {
    [
        VortexParticle::new([-0.4, 0.7], 1.3),
        VortexParticle::new([0.8, 0.3], -0.6),
        VortexParticle::new([0.1, 1.4], 2.2),
    ]
}

#[test]
fn same_native_velocity_and_explicit_image_kinematics() {
    let wall = PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap();
    let particles = sources();
    let mut explicit = particles.to_vec();
    explicit.extend(particles.iter().map(|source| {
        VortexParticle::new([source.pos[0], -source.pos[1]], -source.circulation)
    }));
    for point in [[0.2, 0.6], [-0.7, 1.1], [1.2, -0.4], [0.15, 0.0]] {
        let actual = wall.induced_kinematics(&particles, point, 0.2).unwrap();
        let expected = induced_kinematics(&explicit, point, 0.2).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(
            actual.velocity.map(f64::to_bits),
            wall.induced_velocity(&particles, point, 0.2)
                .unwrap()
                .map(f64::to_bits)
        );
        assert_eq!(actual.gradient[0][0] + actual.gradient[1][1], 0.0);
    }
}

#[test]
fn analytic_jacobian_matches_independent_native_velocity_differences() {
    let wall = PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap();
    let particles = sources();
    let h = 2e-6;
    for point in [[0.2, 0.6], [-0.7, 1.1], [1.2, -0.4], [0.15, 0.0]] {
        let sample = wall.induced_kinematics(&particles, point, 0.2).unwrap();
        for axis in 0..2 {
            let mut plus = point;
            let mut minus = point;
            plus[axis] += h;
            minus[axis] -= h;
            let up = wall.induced_velocity(&particles, plus, 0.2).unwrap();
            let um = wall.induced_velocity(&particles, minus, 0.2).unwrap();
            for component in 0..2 {
                near(
                    sample.gradient[component][axis],
                    (up[component] - um[component]) / (2.0 * h),
                    3e-8,
                );
            }
        }
        near(
            sample.vorticity,
            sample.gradient[1][0] - sample.gradient[0][1],
            2e-14,
        );
    }
}

#[test]
fn reflection_retains_polar_velocity_and_axial_vorticity() {
    let wall = PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap();
    let particles = sources();
    let upper = wall
        .induced_kinematics(&particles, [0.23, 0.6], 0.2)
        .unwrap();
    let lower = wall
        .induced_kinematics(&particles, [0.23, -0.6], 0.2)
        .unwrap();
    near(lower.velocity[0], upper.velocity[0], 2e-15);
    near(lower.velocity[1], -upper.velocity[1], 2e-15);
    for i in 0..2 {
        for j in 0..2 {
            let sign = if i == j { 1.0 } else { -1.0 };
            near(lower.gradient[i][j], sign * upper.gradient[i][j], 1e-14);
        }
    }
    near(lower.vorticity, -upper.vorticity, 1e-14);

    let boundary = wall
        .induced_kinematics(&particles, [0.23, 0.0], 0.2)
        .unwrap();
    near(boundary.velocity[1], 0.0, 2e-15);
    near(boundary.gradient[0][1], 0.0, 1e-14);
    near(boundary.gradient[1][0], 0.0, 1e-14);
    near(boundary.vorticity, 0.0, 1e-14);
}

#[test]
fn wall_rotation_and_translation_transform_the_full_field() {
    let particles = sources();
    let point = [0.2, 0.6];
    let wall = PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap();
    let sample = wall.induced_kinematics(&particles, point, 0.2).unwrap();
    let transform = |p: [f64; 2]| [4.0 - p[1], -3.0 + p[0]];
    let moved = particles.map(|source| {
        VortexParticle::new(transform(source.pos), source.circulation)
    });
    let rotated_wall = PlanarWall2::new(transform([0.0; 2]), [-1.0, 0.0]).unwrap();
    let rotated = rotated_wall
        .induced_kinematics(&moved, transform(point), 0.2)
        .unwrap();
    near(rotated.velocity[0], -sample.velocity[1], 5e-13);
    near(rotated.velocity[1], sample.velocity[0], 5e-13);
    let g = sample.gradient;
    let expected = [[g[1][1], -g[1][0]], [-g[0][1], g[0][0]]];
    for i in 0..2 {
        for j in 0..2 {
            near(rotated.gradient[i][j], expected[i][j], 5e-13);
        }
    }
    near(rotated.vorticity, sample.vorticity, 5e-13);
}

#[test]
fn core_center_and_thin_core_tail_keep_the_actual_vorticity() {
    let wall = PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap();
    let particle = VortexParticle::new([0.0, 1.0], TAU);
    let core = 0.25;
    let sample = wall
        .induced_kinematics(&[particle], particle.pos, core)
        .unwrap();
    let q = 4.0 + core * core;
    let expected = 2.0 / (core * core) - 2.0 * core * core / (q * q);
    near(sample.vorticity, expected, 2e-14);
    // The real source's center velocity vanishes; its image still transports it.
    assert!(sample.velocity[0] > 0.0);
    near(sample.velocity[1], 0.0, 0.0);

    let circulation = 1e300;
    let tiny = 1e-200;
    let tail = wall
        .induced_kinematics(
            &[VortexParticle::new([0.0, 1.0], circulation)],
            [1.0, 2.0],
            tiny,
        )
        .unwrap();
    assert_eq!(tiny * tiny, 0.0);
    assert_eq!(tail.gradient[1][0] - tail.gradient[0][1], 0.0);
    let expected_tail = ((circulation * tiny) * tiny) / PI * (0.25 - 0.01);
    assert!((tail.vorticity / expected_tail - 1.0).abs() < 3e-14);
}

#[test]
fn invalid_sources_singular_probes_and_unrepresentable_fields_refuse() {
    let wall = PlanarWall2::new([0.0; 2], [0.0, 1.0]).unwrap();
    let particle = VortexParticle::new([0.0, 1.0], 1.0);
    assert_eq!(
        wall.induced_kinematics(&[], [0.0; 2], 0.0).unwrap(),
        ParticleKinematics2::default()
    );
    for core in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(wall.induced_kinematics(&[], [0.0; 2], core).is_err());
    }
    assert!(
        wall.induced_kinematics(&[], [f64::NAN, 0.0], 0.1)
            .is_err()
    );
    for source in [
        VortexParticle::new([0.0; 2], 1.0),
        VortexParticle::new([0.0, -1.0], 1.0),
        VortexParticle::new([0.0, f64::INFINITY], 1.0),
        VortexParticle::new([0.0, 1.0], f64::NAN),
    ] {
        assert!(wall.induced_kinematics(&[source], [1.0, 1.0], 0.1).is_err());
    }
    let sources = [particle];
    let before = sources;
    assert!(wall.induced_kinematics(&sources, particle.pos, 0.0).is_err());
    assert!(wall.induced_kinematics(&sources, particle.pos, 1e-200).is_err());
    assert_eq!(sources, before);
    assert_eq!(
        wall.induced_kinematics(
            &[VortexParticle::new(particle.pos, 0.0)],
            particle.pos,
            0.0,
        )
        .unwrap(),
        ParticleKinematics2::default()
    );
}
