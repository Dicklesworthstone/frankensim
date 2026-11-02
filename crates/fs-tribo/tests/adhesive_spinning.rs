//! One physical junction slip field supplies force, spin and dissipated work.
use core::f64::consts::PI;
use fs_tribo::adhesive::spinning::{AdhesiveSlidingPatch, AdhesiveSlidingPatchParameters};

fn patch(radius: f64, strength: f64) -> AdhesiveSlidingPatch {
    AdhesiveSlidingPatch::new(AdhesiveSlidingPatchParameters {
        contact_radius_m: radius,
        shear_strength_pa: strength,
    })
    .unwrap()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * scale,
        "{actual:e} != {expected:e}, scale {scale:e}"
    );
}

#[test]
fn translation_spin_and_detachment_have_their_actual_area_limits() {
    let a = 0.003;
    let tau = 2e6;
    let p = patch(a, tau);
    let capacity = tau * PI * a * a;
    close(p.sliding_capacity_n(), capacity, capacity);
    close(
        p.spinning_capacity_nm(),
        (2.0 / 3.0) * capacity * a,
        capacity * a,
    );
    let sliding = p.evaluate([0.3, 0.4], 0.0).unwrap();
    close(sliding.force_n[0], -0.6 * capacity, capacity);
    close(sliding.force_n[1], -0.8 * capacity, capacity);
    assert_eq!(sliding.spin_moment_nm, 0.0);
    close(sliding.dissipated_power_w, 0.5 * capacity, capacity);
    let spin = p.evaluate([0.0; 2], -2.0).unwrap();
    assert_eq!(spin.force_n, [0.0; 2]);
    close(
        spin.spin_moment_nm,
        (2.0 / 3.0) * capacity * a,
        capacity * a,
    );
    close(
        spin.dissipated_power_w,
        (4.0 / 3.0) * capacity * a,
        capacity * a,
    );
    assert!(p.evaluate([0.0; 2], 0.0).unwrap().at_rest);
    for p in [patch(0.0, tau), patch(a, 0.0)] {
        let w = p.evaluate([3.0, 4.0], 7.0).unwrap();
        assert_eq!(w.force_n, [0.0; 2]);
        assert_eq!(w.spin_moment_nm, 0.0);
        assert_eq!(w.dissipated_power_w, 0.0);
        assert!(!w.at_rest);
    }
}

// Independent positive surface integrals in polar coordinates centered on the
// instantaneous zero-slip point. Radial integration is exact. The remaining
// quarter-angle integrals use composite Simpson, with no elliptic functions,
// paired Taylor coefficients or sampled traction normalization from production.
fn surface_oracle(slip_to_spin: f64) -> (f64, f64) {
    let inside = slip_to_spin <= 1.0;
    let z = if inside {
        slip_to_spin
    } else {
        1.0 / slip_to_spin
    };
    let n = 8192;
    let h = 0.5 * PI / f64::from(n);
    let mut force = 0.0;
    let mut moment = 0.0;
    for i in 0..=n {
        let angle = f64::from(i) * h;
        let c = angle.cos();
        let root = ((1.0 - z) * (1.0 + z) + z * z * c * c).sqrt();
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        force += weight * c * c * root;
        moment += weight
            * if inside {
                root.powi(3)
            } else {
                if root == 0.0 { 0.0 } else { c.powi(4) / root }
            };
    }
    force *= (h / 3.0) * (4.0 / PI) * if inside { z } else { 1.0 };
    moment *= (h / 3.0) * (2.0 / PI) * if inside { 1.0 } else { z };
    (force, moment)
}

#[test]
fn combined_force_spin_and_power_match_independent_surface_integration() {
    let p = patch(0.003, 2e6);
    let force_scale = p.sliding_capacity_n();
    let moment_scale = p.spinning_capacity_nm();
    let omega = 2.0;
    for ratio in [
        1e-8, 0.05, 0.3, 0.49999, 0.5, 0.50001, 0.8, 0.9999, 1.0, 1.0001, 1.2, 1.9999, 2.0, 2.0001,
        10.0, 1e8,
    ] {
        let speed = ratio * omega * p.parameters().contact_radius_m;
        let w = p.evaluate([speed, 0.0], omega).unwrap();
        let (f, m) = surface_oracle(ratio);
        close(-w.force_n[0], force_scale * f, force_scale * f);
        close(-w.spin_moment_nm, moment_scale * m, moment_scale * m);
        let independent_power = force_scale * f * speed + moment_scale * m * omega;
        close(w.dissipated_power_w, independent_power, independent_power);
        close(w.power_residual_w, 0.0, independent_power);
        assert!(w.dissipated_power_w > 0.0);
    }
}

#[test]
fn slip_point_at_the_rim_is_finite_and_force_and_moment_are_coupled() {
    let p = patch(1.0, 1.0 / PI);
    let w = p.evaluate([1.0, 0.0], 1.0).unwrap();
    close(-w.force_n[0], 8.0 / (3.0 * PI), 1.0);
    close(-w.spin_moment_nm, 8.0 / (9.0 * PI), 1.0);
    assert!(-w.force_n[0] < p.sliding_capacity_n());
    assert!(-w.spin_moment_nm < p.spinning_capacity_nm());
    for ratio in [1.0 - 2.0 * f64::EPSILON, 1.0 + 2.0 * f64::EPSILON] {
        let adjacent = p.evaluate([ratio, 0.0], 1.0).unwrap();
        close(adjacent.force_n[0], w.force_n[0], 1.0);
        close(adjacent.spin_moment_nm, w.spin_moment_nm, 1.0);
    }
}

#[test]
fn frame_rotation_rate_homogeneity_and_sign_reversal_preserve_physical_wrench() {
    let p = patch(0.2, 2.0);
    let original = p.evaluate([0.03, 0.04], 0.5).unwrap();
    let rotated = p.evaluate([-0.04, 0.03], 0.5).unwrap();
    close(
        rotated.force_n[0],
        -original.force_n[1],
        p.sliding_capacity_n(),
    );
    close(
        rotated.force_n[1],
        original.force_n[0],
        p.sliding_capacity_n(),
    );
    close(
        rotated.spin_moment_nm,
        original.spin_moment_nm,
        p.spinning_capacity_nm(),
    );
    for scale in [1e-200, 1e-12, 1.0, 1e12, 1e200] {
        let w = p
            .evaluate([0.03 * scale, 0.04 * scale], 0.5 * scale)
            .unwrap();
        close(w.force_n[0], original.force_n[0], p.sliding_capacity_n());
        close(w.force_n[1], original.force_n[1], p.sliding_capacity_n());
        close(
            w.spin_moment_nm,
            original.spin_moment_nm,
            p.spinning_capacity_nm(),
        );
        close(
            w.dissipated_power_w,
            original.dissipated_power_w * scale,
            original.dissipated_power_w * scale,
        );
    }
    let reverse = p.evaluate([-0.03, -0.04], -0.5).unwrap();
    close(
        reverse.force_n[0],
        -original.force_n[0],
        p.sliding_capacity_n(),
    );
    close(
        reverse.force_n[1],
        -original.force_n[1],
        p.sliding_capacity_n(),
    );
    close(
        reverse.spin_moment_nm,
        -original.spin_moment_nm,
        p.spinning_capacity_nm(),
    );
    close(reverse.dissipated_power_w, original.dissipated_power_w, 1.0);
}

#[test]
fn every_generated_wrench_obeys_the_same_convex_dissipation_support() {
    let p = patch(0.2, 2.0);
    let motions = [
        [0.0, 0.0, 2.0],
        [0.3, 0.4, 0.0],
        [0.4, 0.0, 2.0],
        [-0.08, 0.06, -2.0],
        [1e-8, 0.0, 1.0],
        [1.0, 0.0, 1e-8],
    ];
    for v in motions {
        let w = p.evaluate([v[0], v[1]], v[2]).unwrap();
        for u in motions {
            let other = p.evaluate([u[0], u[1]], u[2]).unwrap();
            let support =
                -other.force_n[0] * v[0] - other.force_n[1] * v[1] - other.spin_moment_nm * v[2];
            assert!(support <= w.dissipated_power_w + 3e-14);
        }
    }
}

#[test]
fn invalid_geometry_and_nonfinite_motion_refuse_without_smoothing_tiny_rates() {
    for (a, strength) in [
        (-1.0, 1.0),
        (1.0, -1.0),
        (f64::NAN, 1.0),
        (1.0, f64::INFINITY),
        (1e-200, 1.0),
    ] {
        assert!(
            AdhesiveSlidingPatch::new(AdhesiveSlidingPatchParameters {
                contact_radius_m: a,
                shear_strength_pa: strength,
            })
            .is_err()
        );
    }
    let p = patch(1.0, 1.0);
    assert!(p.evaluate([f64::NAN, 0.0], 1.0).is_err());
    assert!(p.evaluate([0.0; 2], f64::INFINITY).is_err());
    let tiny = p.evaluate([1e-200, 0.0], 1.0).unwrap();
    assert!(tiny.force_n[0] < 0.0);
    close(tiny.force_n[0], -PI * 1e-200, PI * 1e-200);
    let tiny_spin = p.evaluate([1.0, 0.0], 1e-200).unwrap();
    assert!(tiny_spin.spin_moment_nm < 0.0);
    close(tiny_spin.spin_moment_nm, -0.25 * PI * 1e-200, PI * 1e-200);
}
