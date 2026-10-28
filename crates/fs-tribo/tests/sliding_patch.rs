//! G1/G3: exact Hertz sliding/spin integrals against independent spatial work.

use core::f64::consts::PI;

use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

fn patch(normal: f64, radius: f64, mu: f64) -> HertzSlidingPatch {
    HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: normal,
        contact_radius_m: radius,
        friction_coefficient: mu,
    })
    .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}

/// Independent midpoint surface quadrature; the production owner uses no grid.
/// r=a*sin(theta) turns p(r)*r*dr into 3P sin(theta) cos(theta)^2 dtheta.
/// Force, moment, and local-slip heat are accumulated independently.
fn spatial(
    parameters: HertzSlidingPatchParameters,
    slip: [f64; 2],
    spin: f64,
    radial: usize,
    angular: usize,
) -> [f64; 4] {
    let mut output = [0.0; 4];
    let dr = (PI / 2.0) / radial as f64;
    for i in 0..radial {
        let theta = (i as f64 + 0.5) * dr;
        let (sine, cosine) = theta.sin_cos();
        let radius = parameters.contact_radius_m * sine;
        let force = 3.0
            * parameters.normal_force_n
            * parameters.friction_coefficient
            * sine
            * cosine
            * cosine
            * dr
            / angular as f64;
        for j in 0..angular {
            let angle = 2.0 * PI * (j as f64 + 0.5) / angular as f64;
            let (s, c) = angle.sin_cos();
            let x = radius * c;
            let y = radius * s;
            let velocity = [slip[0] - spin * y, slip[1] + spin * x];
            let speed = velocity[0].hypot(velocity[1]);
            if speed > 0.0 {
                let local = velocity.map(|v| -force * v / speed);
                output[0] += local[0];
                output[1] += local[1];
                output[2] += x * local[1] - y * local[0];
                output[3] += force * speed;
            }
        }
    }
    output
}

#[test]
fn pure_translation_and_spin_use_pressure_first_moment_not_rms_radius() {
    let law = patch(12.0, 0.02, 0.4);
    let translate = law.evaluate([3.0, -4.0], 0.0).unwrap();
    close(translate.force_n[0], -2.88, 2e-15);
    close(translate.force_n[1], 3.84, 2e-15);
    assert_eq!(translate.spin_moment_nm, 0.0);
    close(translate.dissipated_power_w, 24.0, 8e-15);
    let expected = (3.0 * PI / 16.0) * 4.8 * 0.02;
    for spin in [-17.0, 17.0] {
        let result = law.evaluate([0.0; 2], spin).unwrap();
        assert_eq!(result.force_n, [0.0; 2]);
        close(result.spin_moment_nm, -expected.copysign(spin), 3e-17);
        close(result.dissipated_power_w, expected * spin.abs(), 3e-16);
        close(result.power_residual_w, 0.0, 3e-16);
    }
    assert!((expected - 0.4_f64.sqrt() * 4.8 * 0.02).abs() > 1e-3);
}

#[test]
fn sliding_and_spinning_share_one_friction_capacity() {
    let law = patch(1.0, 1.0, 1.0);
    let result = law.evaluate([1.0, 0.0], 1.0).unwrap();
    close(result.force_n[0], -9.0 * PI / 32.0, 4e-16);
    close(result.spin_moment_nm, -9.0 * PI / 128.0, 2e-16);
    assert!(result.force_n[0].abs() < law.sliding_capacity_n());
    assert!(result.spin_moment_nm.abs() < law.spinning_capacity_nm());
    let faster_translation = law.evaluate([10.0, 0.0], 1.0).unwrap();
    assert!(faster_translation.force_n[0].abs() > result.force_n[0].abs());
    assert!(faster_translation.spin_moment_nm.abs() < result.spin_moment_nm.abs());
}

#[test]
fn exact_wrench_matches_independent_hertz_surface_integrals() {
    let law = patch(3.2, 0.4, 0.7);
    for epsilon in [0.2, 0.75, 1.0, 1.5, 2.5, 20.0] {
        let slip = [0.6 * epsilon, -0.8 * epsilon];
        let result = law.evaluate(slip, 2.5).unwrap();
        let reference = spatial(law.parameters(), slip, 2.5, 384, 1024);
        close(result.force_n[0], reference[0], 2.0e-5);
        close(result.force_n[1], reference[1], 2.0e-5);
        close(result.spin_moment_nm, reference[2], 8.0e-6);
        close(result.dissipated_power_w, reference[3], 1.0e-4);
        close(result.power_residual_w, 0.0, 2e-14);
        assert!(result.dissipated_power_w > 0.0);
    }
}

#[test]
fn independent_surface_reference_has_a_live_three_rung_refinement() {
    let law = patch(1.0, 1.0, 1.0);
    let slip = [1.5, -2.0];
    let result = law.evaluate(slip, 1.0).unwrap();
    let exact = [
        result.force_n[0],
        result.force_n[1],
        result.spin_moment_nm,
        result.dissipated_power_w,
    ];
    let mut errors = Vec::new();
    for radial in [16, 32, 64] {
        let reference = spatial(law.parameters(), slip, 1.0, radial, 128);
        let error = reference
            .into_iter()
            .zip(exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push(error);
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order > 1.95 && order < 2.05,
            "radial errors {errors:?}, order {order}"
        );
    }
    // The angular axis is independently live near the interior stagnation point.
    let exact = law.evaluate([0.6, 0.0], 1.0).unwrap();
    let mut angular_errors = Vec::new();
    for angular in [16, 32, 64] {
        let reference = spatial(law.parameters(), [0.6, 0.0], 1.0, 512, angular);
        angular_errors.push((reference[0] - exact.force_n[0]).abs());
    }
    assert!(angular_errors[0] > angular_errors[1]);
    assert!(angular_errors[1] > angular_errors[2]);
    assert!(angular_errors[2] < angular_errors[0] / 8.0);
    eprintln!(
        "Hertz surface reference radial errors {errors:?}, angular errors {angular_errors:?}"
    );
}

#[test]
fn arbitrary_slip_directions_rotation_reflection_and_rate_scaling_are_covariant() {
    let law = patch(7.0, 0.12, 0.3);
    let slip = [0.3, -0.4];
    let spin = 4.0;
    let base = law.evaluate(slip, spin).unwrap();
    let angle: f64 = 0.713;
    let (s, c) = angle.sin_cos();
    let rotate = |v: [f64; 2]| [c * v[0] - s * v[1], s * v[0] + c * v[1]];
    let rotated = law.evaluate(rotate(slip), spin).unwrap();
    for (actual, expected) in rotated.force_n.into_iter().zip(rotate(base.force_n)) {
        close(actual, expected, 2e-15);
    }
    close(rotated.spin_moment_nm, base.spin_moment_nm, 2e-16);
    close(rotated.dissipated_power_w, base.dissipated_power_w, 2e-15);
    let reflected = law.evaluate([slip[0], -slip[1]], -spin).unwrap();
    assert_eq!(reflected.force_n, [base.force_n[0], -base.force_n[1]]);
    assert_eq!(reflected.spin_moment_nm, -base.spin_moment_nm);
    for scale in [1e-200, 1e-20, 1.0, 1e20, 1e200] {
        let scaled = law.evaluate(slip.map(|x| scale * x), scale * spin).unwrap();
        for (actual, expected) in scaled.force_n.into_iter().zip(base.force_n) {
            close(actual, expected, 2e-15);
        }
        close(scaled.spin_moment_nm, base.spin_moment_nm, 2e-16);
        close(
            scaled.dissipated_power_w / scale,
            base.dissipated_power_w,
            2e-15,
        );
    }
    assert_eq!(base, law.evaluate(slip, spin).unwrap());
}

#[test]
fn vanishing_spin_remains_resolved_without_large_ratio_cancellation() {
    let law = patch(1.0, 1.0, 1.0);
    for z in [1e-4, 1e-8, 1e-100, 1e-300] {
        let result = law.evaluate([1.0, 0.0], z).unwrap();
        // Independent leading moments of the Hertz pressure: <y^2>=a^2/5.
        close(-result.spin_moment_nm / z, 0.2, 2e-9);
        close(-result.force_n[0], 1.0 - z * z / 10.0, 4e-16);
    }
    for epsilon in [1.0, 2.0] {
        let before = law.evaluate([epsilon * (1.0 - 1e-9), 0.0], 1.0).unwrap();
        let after = law.evaluate([epsilon * (1.0 + 1e-9), 0.0], 1.0).unwrap();
        close(before.force_n[0], after.force_n[0], 2e-9);
        close(before.spin_moment_nm, after.spin_moment_nm, 2e-9);
    }
}

#[test]
fn rest_frictionless_limits_and_invalid_inputs_are_explicit() {
    let law = patch(1.0, 1.0, 0.5);
    let rest = law.evaluate([0.0; 2], 0.0).unwrap();
    assert!(rest.at_rest);
    assert_eq!(rest.force_n, [0.0; 2]);
    assert_eq!(rest.spin_moment_nm, 0.0);
    assert_eq!(rest.dissipated_power_w, 0.0);
    for no_friction in [patch(0.0, 1.0, 0.5), patch(1.0, 1.0, 0.0)] {
        let response = no_friction.evaluate([3.0, 4.0], 7.0).unwrap();
        assert!(!response.at_rest);
        assert_eq!(response.force_n, [0.0; 2]);
        assert_eq!(response.spin_moment_nm, 0.0);
        assert_eq!(response.dissipated_power_w, 0.0);
    }
    for (normal, radius, mu) in [
        (-1.0, 1.0, 0.5),
        (1.0, 0.0, 0.5),
        (1.0, 1.0, -0.5),
        (f64::INFINITY, 1.0, 0.5),
    ] {
        assert!(
            HertzSlidingPatch::new(HertzSlidingPatchParameters {
                normal_force_n: normal,
                contact_radius_m: radius,
                friction_coefficient: mu,
            })
            .is_err()
        );
    }
    assert!(law.evaluate([f64::NAN, 0.0], 0.0).is_err());
    assert!(law.evaluate([0.0; 2], f64::INFINITY).is_err());
    let enormous = patch(1e200, 1.0, 1.0);
    assert!(enormous.evaluate([1e200, 0.0], 0.0).is_err());
}
