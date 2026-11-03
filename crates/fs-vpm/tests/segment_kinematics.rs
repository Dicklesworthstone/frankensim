//! G0/G1: actual finite-core velocity derivatives, vorticity and length scaling.
//! The reference integrates the differentiated Biot–Savart kernel directly;
//! it does not differentiate the implementation's endpoint expression.

use core::f64::consts::PI;
use fs_vpm::filament3d::{
    SegmentKinematics, segment_kinematics_with_core, segment_velocity_with_core,
};

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn curl(g: [[f64; 3]; 3]) -> [f64; 3] {
    [g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1]]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn quadrature(p: [f64; 3], a: [f64; 3], b: [f64; 3], core2: f64) -> SegmentKinematics {
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let mut out = SegmentKinematics::default();
    let n = 4096;
    for k in 0..=n {
        let fraction = k as f64 / n as f64;
        let r = [
            p[0] - a[0] - fraction * d[0],
            p[1] - a[1] - fraction * d[1],
            p[2] - a[2] - fraction * d[2],
        ];
        let r2 = r.iter().map(|v| v * v).sum::<f64>() + core2;
        let inverse3 = 1.0 / (r2 * r2.sqrt());
        let inverse5 = inverse3 / r2;
        let dr = cross(d, r);
        let weight = if k == 0 || k == n {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let scale = weight / (3.0 * n as f64 * 4.0 * PI);
        for i in 0..3 {
            out.velocity[i] += scale * dr[i] * inverse3;
        }
        for j in 0..3 {
            let mut e = [0.0; 3];
            e[j] = 1.0;
            let de = cross(d, e);
            for i in 0..3 {
                out.gradient[i][j] += scale * (de[i] * inverse3 - 3.0 * dr[i] * r[j] * inverse5);
            }
        }
    }
    out
}

#[test]
fn exact_jacobian_matches_independent_integral_in_all_branches() {
    let a = [-0.7, 0.1, -0.3];
    let b = [0.6, 0.8, 1.4];
    for p in [[0.2, -0.3, 0.6], [4.2, 2.1, 3.0], [-3.0, -1.0, -4.0], a, b] {
        let actual = segment_kinematics_with_core(p, a, b, 0.12).unwrap();
        let reference = quadrature(p, a, b, 0.12);
        for i in 0..3 {
            near(actual.velocity[i], reference.velocity[i], 2e-12);
            near(
                actual.velocity[i],
                segment_velocity_with_core(p, a, b, 0.12)[i],
                2e-15,
            );
            for j in 0..3 {
                near(actual.gradient[i][j], reference.gradient[i][j], 3e-11);
            }
        }
        near(
            actual.gradient[0][0] + actual.gradient[1][1] + actual.gradient[2][2],
            0.0,
            2e-15,
        );
    }
}

#[test]
fn positive_core_axis_has_zero_velocity_but_finite_vorticity() {
    let radius = 0.2;
    let half_length = 1000.0;
    let sample = segment_kinematics_with_core(
        [0.0; 3],
        [0.0, 0.0, -half_length],
        [0.0, 0.0, half_length],
        radius * radius,
    )
    .unwrap();
    assert_eq!(sample.velocity, [0.0; 3]);
    let exact =
        half_length / (PI * radius * radius * (half_length * half_length + radius * radius).sqrt());
    near(curl(sample.gradient)[2], exact, 1e-13);
    near(curl(sample.gradient)[2], 1.0 / (PI * radius * radius), 2e-7);
    near(sample.gradient[0][1], -exact / 2.0, 1e-13);
    near(sample.gradient[1][0], exact / 2.0, 1e-13);
}

#[test]
fn long_filament_recovers_distributed_core_vorticity_and_circulation() {
    let core2 = 0.16;
    let radius = 0.6;
    let sample =
        segment_kinematics_with_core([radius, 0.0, 0.0], [0.0, 0.0, -1e5], [0.0, 0.0, 1e5], core2)
            .unwrap();
    let denominator = radius * radius + core2;
    near(sample.velocity[1], radius / (2.0 * PI * denominator), 1e-10);
    near(
        curl(sample.gradient)[2],
        core2 / (PI * denominator * denominator),
        1e-10,
    );
    // Stokes circulation through a finite radius: the algebraic core has tails.
    near(
        2.0 * PI * radius * sample.velocity[1],
        radius * radius / denominator,
        1e-10,
    );
}

#[test]
fn closed_zero_core_ring_has_no_off_filament_curl() {
    let ring = [
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
    ];
    for p in [[0.2, 0.3, 0.7], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]] {
        let mut g = [[0.0; 3]; 3];
        for k in 0..4 {
            let s = segment_kinematics_with_core(p, ring[k], ring[(k + 1) % 4], 0.0).unwrap();
            for i in 0..3 {
                for j in 0..3 {
                    g[i][j] += s.gradient[i][j];
                }
            }
        }
        for w in curl(g) {
            near(w, 0.0, 2e-15);
        }
    }
}

#[test]
fn length_scaling_and_reversal_preserve_physical_units() {
    let a = [-0.2, 0.1, -0.5];
    let b = [0.7, -0.4, 1.5];
    let p = [0.3, 0.5, 0.2];
    let reference = segment_kinematics_with_core(p, a, b, 0.07).unwrap();
    let reverse = segment_kinematics_with_core(p, b, a, 0.07).unwrap();
    for scale in [1e-8, 1e-3, 1.0, 1e3, 1e8] {
        let scaled = segment_kinematics_with_core(
            p.map(|x| x * scale),
            a.map(|x| x * scale),
            b.map(|x| x * scale),
            0.07 * scale * scale,
        )
        .unwrap();
        for i in 0..3 {
            near(scaled.velocity[i] * scale, reference.velocity[i], 5e-15);
            near(reverse.velocity[i], -reference.velocity[i], 5e-15);
            for j in 0..3 {
                near(
                    scaled.gradient[i][j] * scale * scale,
                    reference.gradient[i][j],
                    2e-14,
                );
                near(reverse.gradient[i][j], -reference.gradient[i][j], 2e-14);
            }
        }
    }
}

#[test]
fn thin_exterior_and_singular_axis_are_distinguished() {
    let a = [0.0, 0.0, 0.0];
    let b = [0.0, 0.0, 1.0];
    for core2 in [0.0, 1e-16] {
        let p = [1e-10, -2e-10, 30.0];
        let s = segment_kinematics_with_core(p, a, b, core2).unwrap();
        let q = quadrature(p, a, b, core2);
        for i in 0..3 {
            near(s.velocity[i], q.velocity[i], 1e-27);
            for j in 0..3 {
                near(s.gradient[i][j], q.gradient[i][j], 1e-18);
            }
        }
    }
    let outside = segment_kinematics_with_core([0.0, 0.0, 3.0], a, b, 0.0).unwrap();
    assert_eq!(outside.velocity, [0.0; 3]);
    assert!(outside.gradient[1][0] > 0.0);
    for z in [0.0, 0.5, 1.0] {
        assert!(segment_kinematics_with_core([0.0, 0.0, z], a, b, 0.0).is_none());
        assert!(segment_kinematics_with_core([0.0, 0.0, z], a, b, 0.01).is_some());
    }
    assert_eq!(
        segment_kinematics_with_core(a, b, b, 0.0),
        Some(SegmentKinematics::default())
    );
    for core2 in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(segment_kinematics_with_core(a, a, b, core2).is_none());
    }
    assert!(segment_kinematics_with_core(a, a, [1e200, 0.0, 0.0], 1.0).is_none());
    assert!(segment_kinematics_with_core(a, a, [1e-200, 0.0, 0.0], 1.0).is_none());
}
