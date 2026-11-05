//! G0/G1/G3: analytic curl, independent derivatives and Stokes circulation.

use core::f64::consts::{PI, TAU};
use fs_vpm::{ParticleKinematics2, VortexParticle, induced_kinematics, induced_velocity};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn relative(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 3e-14,
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn positive_core_center_has_rotation_with_zero_velocity() {
    let particles = [VortexParticle::new([0.0; 2], TAU)];
    let sample = induced_kinematics(&particles, [0.0; 2], 0.25).unwrap();
    assert_eq!(sample.velocity, [0.0; 2]);
    near(sample.gradient[0][0], 0.0, 0.0);
    near(sample.gradient[1][1], 0.0, 0.0);
    near(sample.gradient[0][1], -16.0, 4e-15);
    near(sample.gradient[1][0], 16.0, 4e-15);
    near(sample.vorticity, 32.0, 8e-15);
}

#[test]
fn signed_cloud_jacobian_matches_independent_velocity_differences() {
    let particles = [
        VortexParticle::new([-0.5, 0.2], TAU),
        VortexParticle::new([0.7, -0.3], -1.3),
        VortexParticle::new([0.1, 0.0], 0.2),
    ];
    for point in [[0.3, 0.8], [-0.8, -0.4], [2.0, 0.1]] {
        let sample = induced_kinematics(&particles, point, 0.2).unwrap();
        assert_eq!(
            sample.velocity.map(f64::to_bits),
            induced_velocity(&particles, point, 0.2).map(f64::to_bits)
        );
        for axis in 0..2 {
            let mut p = point;
            let mut m = point;
            p[axis] += 1e-5;
            m[axis] -= 1e-5;
            let up = induced_velocity(&particles, p, 0.2);
            let um = induced_velocity(&particles, m, 0.2);
            for component in 0..2 {
                let derivative = (up[component] - um[component]) / 2e-5;
                near(sample.gradient[component][axis], derivative, 3e-9);
            }
        }
        assert_eq!(sample.gradient[0][0] + sample.gradient[1][1], 0.0);
        near(
            sample.vorticity,
            sample.gradient[1][0] - sample.gradient[0][1],
            2e-14,
        );
        assert_eq!(sample, induced_kinematics(&particles, point, 0.2).unwrap());
    }
}

#[test]
fn sampled_curl_recovers_finite_disk_stokes_circulation() {
    let particles = [VortexParticle::new([0.0; 2], 3.5)];
    let radius = 1.3;
    let core = 0.4;
    let intervals = 4096;
    let h = radius / intervals as f64;
    let mut integral = 0.0;
    for k in 0..=intervals {
        let r = k as f64 * h;
        let omega = induced_kinematics(&particles, [r, 0.0], core)
            .unwrap()
            .vorticity;
        let weight = if k == 0 || k == intervals {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * TAU * r * omega;
    }
    integral *= h / 3.0;
    let boundary = induced_velocity(&particles, [radius, 0.0], core);
    let enclosed = 3.5 * radius * radius / (radius * radius + core * core);
    near(integral, TAU * radius * boundary[1], 5e-12);
    near(integral, enclosed, 5e-12);
    assert!(integral < 3.5); // Algebraic-core vorticity has noncompact tails.
}

#[test]
fn point_vortex_is_irrotational_away_from_its_singular_center() {
    let particles = [VortexParticle::new([0.0; 2], TAU)];
    let sample = induced_kinematics(&particles, [2.0, 0.0], 0.0).unwrap();
    assert_eq!(sample.velocity, [0.0, 0.5]);
    near(sample.gradient[0][1], -0.25, 1e-16);
    near(sample.gradient[1][0], -0.25, 1e-16);
    assert_eq!(sample.vorticity, 0.0);
    assert!(induced_kinematics(&particles, [0.0; 2], 0.0).is_none());
    assert_eq!(
        induced_kinematics(&[VortexParticle::new([0.0; 2], 0.0)], [0.0; 2], 0.0),
        Some(ParticleKinematics2::default())
    );
}

#[test]
fn spatial_and_circulation_scaling_retains_kinematics_at_exceptional_ranges() {
    let source = VortexParticle::new([-0.25, 0.125], 1.3);
    let point = [0.5, 0.75];
    let core = 0.15;
    let reference = induced_kinematics(&[source], point, core).unwrap();
    // Length and circulation both scale by s: velocity stays fixed, whereas
    // the spatial Jacobian and vorticity scale by 1/s. Many squared lengths
    // and fourth powers here cannot be represented, but the fields can.
    for s in [1e-300, 1e-180, 1e-90, 1e-8, 1.0, 1e8, 1e90, 1e180, 1e300] {
        let particles = [VortexParticle::new(
            source.pos.map(|v| v * s),
            source.circulation * s,
        )];
        let p = point.map(|v| v * s);
        let sample = induced_kinematics(&particles, p, core * s).unwrap();
        assert_eq!(
            sample.velocity.map(f64::to_bits),
            induced_velocity(&particles, p, core * s).map(f64::to_bits)
        );
        for i in 0..2 {
            near(sample.velocity[i], reference.velocity[i], 2e-14);
            for j in 0..2 {
                near(sample.gradient[i][j] * s, reference.gradient[i][j], 2e-14);
            }
        }
        near(sample.vorticity * s, reference.vorticity, 2e-14);
    }
}

#[test]
fn rotation_and_translation_transform_the_full_gradient() {
    let source = VortexParticle::new([-0.25, 0.125], -2.0);
    let point = [0.5, 0.75];
    let original = induced_kinematics(&[source], point, 0.2).unwrap();
    let transform = |p: [f64; 2]| [3.0 - p[1], -4.0 + p[0]];
    let moved = induced_kinematics(
        &[VortexParticle::new(
            transform(source.pos),
            source.circulation,
        )],
        transform(point),
        0.2,
    )
    .unwrap();
    near(moved.velocity[0], -original.velocity[1], 1e-15);
    near(moved.velocity[1], original.velocity[0], 1e-15);
    let g = original.gradient;
    let expected = [[g[1][1], -g[1][0]], [-g[0][1], g[0][0]]];
    for i in 0..2 {
        for j in 0..2 {
            near(moved.gradient[i][j], expected[i][j], 2e-15);
        }
    }
    near(moved.vorticity, original.vorticity, 2e-15);
}

#[test]
fn tiny_core_curl_survives_square_underflow_and_jacobian_cancellation() {
    let circulation = 1e300;
    let core = 1e-200;
    let sample = induced_kinematics(
        &[VortexParticle::new([0.0; 2], circulation)],
        [1.0, 0.0],
        core,
    )
    .unwrap();
    assert_eq!(core * core, 0.0);
    assert_eq!(sample.gradient[1][0] - sample.gradient[0][1], 0.0);
    relative(sample.vorticity, ((circulation * core) * core) / PI);

    for (radius, gamma) in [
        (1e-180, 1e-300),
        (1e180, 1e300),
        (f64::MIN_POSITIVE, f64::from_bits(1)),
    ] {
        let center =
            induced_kinematics(&[VortexParticle::new([0.0; 2], gamma)], [0.0; 2], radius).unwrap();
        let expected = ((gamma / radius) / radius) / PI;
        relative(center.vorticity, expected);
        relative(center.gradient[0][1], -0.5 * expected);
        relative(center.gradient[1][0], 0.5 * expected);
        assert_eq!(center.velocity, [0.0; 2]);
    }
}

#[test]
fn invalid_singular_and_unrepresentable_fields_refuse_without_mutating_sources() {
    let particle = VortexParticle::new([0.0; 2], 1.0);
    assert_eq!(
        induced_kinematics(&[], [0.0; 2], 0.0),
        Some(ParticleKinematics2::default())
    );
    for core in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(induced_kinematics(&[], [0.0; 2], core).is_none());
    }
    assert!(induced_kinematics(&[], [f64::NAN, 0.0], 0.1).is_none());
    assert!(induced_kinematics(&[particle], [0.0; 2], 1e-200).is_none());
    for invalid in [
        VortexParticle::new([f64::INFINITY, 0.0], 1.0),
        VortexParticle::new([0.0; 2], f64::NAN),
    ] {
        assert!(induced_kinematics(&[invalid], [1.0, 0.0], 0.1).is_none());
    }
    assert!(
        induced_kinematics(
            &[VortexParticle::new([-f64::MAX, 0.0], 1.0)],
            [f64::MAX, 0.0],
            1.0,
        )
        .is_none()
    );
    let sources = [VortexParticle::new([0.0; 2], 1e308); 8];
    let before = sources;
    assert!(induced_kinematics(&sources, [0.0; 2], 1.0).is_none());
    assert_eq!(sources, before);
}
