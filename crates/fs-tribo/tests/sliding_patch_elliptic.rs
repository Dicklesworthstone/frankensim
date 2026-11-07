//! Original directional support, independently integrated over the physical
//! Hertz ellipse. All coefficients are authored fixtures, not material fits.

use core::f64::consts::{FRAC_PI_2, PI};
use fs_tribo::sliding_patch::elliptic::*;
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters, HertzSlidingWrench};
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};
use std::cell::Cell;

#[path = "sliding_patch_elliptic/resolved.rs"]
mod resolved;

fn direction(angle: f64) -> [f64; 2] {
    [angle.cos(), angle.sin()]
}
fn rotate(axis: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    [
        axis[0] * v[0] - axis[1] * v[1],
        axis[1] * v[0] + axis[0] * v[1],
    ]
}
fn point_law(mu: [f64; 2], axis: [f64; 2]) -> EllipticCoulomb {
    EllipticCoulomb::new(mu.map(|x| x + 0.3), mu, axis).unwrap()
}
fn model(
    normal: f64,
    axes: [f64; 2],
    angle: f64,
    mu: [f64; 2],
    texture: f64,
) -> EllipticHertzSlidingPatch {
    EllipticHertzSlidingPatch::new(
        EllipticHertzPatchParameters {
            normal_force_n: normal,
            semi_axes_m: axes,
            footprint_axis: direction(angle),
        },
        point_law(mu, direction(texture)),
    )
    .unwrap()
}
fn grid(angular_panels: usize) -> EllipticHertzQuadrature {
    EllipticHertzQuadrature { angular_panels }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}
fn values(r: HertzSlidingWrench) -> [f64; 4] {
    [
        r.force_n[0],
        r.force_n[1],
        r.spin_moment_nm,
        r.dissipated_power_w,
    ]
}

// An unrelated full two-dimensional surface rule. The radial substitution
// r=sin(alpha) removes the Hertz edge singularity; azimuth uses midpoint nodes.
// It evaluates local velocity/support/force/torque/power directly, with no
// zero-slip line coordinates or production radial moments.
fn surface_reference(
    model: &EllipticHertzSlidingPatch,
    slip: [f64; 2],
    spin: f64,
    radial_panels: usize,
    azimuths: usize,
) -> [f64; 4] {
    let p = model.parameters();
    let mu = model.friction_law().kinetic_coefficients();
    let texture = model.friction_law().initial_state().principal_axis();
    let rule = [
        (-0.8611363115940526, 0.34785484513745385),
        (-0.3399810435848563, 0.6521451548625461),
        (0.3399810435848563, 0.6521451548625461),
        (0.8611363115940526, 0.34785484513745385),
    ];
    let mut sum = [0.0; 4];
    let da = FRAC_PI_2 / radial_panels as f64;
    for i in 0..radial_panels {
        for (node, weight) in rule {
            let alpha = (i as f64 + 0.5 + 0.5 * node) * da;
            let (s, c) = alpha.sin_cos();
            let load = 3.0 * p.normal_force_n * s * c * c * (0.5 * da * weight) / azimuths as f64;
            for j in 0..azimuths {
                let theta = 2.0 * PI * (j as f64 + 0.5) / azimuths as f64;
                let (st, ct) = theta.sin_cos();
                let r = rotate(
                    p.footprint_axis,
                    [p.semi_axes_m[0] * s * ct, p.semi_axes_m[1] * s * st],
                );
                let v = [slip[0] - spin * r[1], slip[1] + spin * r[0]];
                let local = rotate([texture[0], -texture[1]], v);
                let norm = (mu[0] * local[0]).hypot(mu[1] * local[1]);
                if norm == 0.0 {
                    continue;
                }
                let force = rotate(
                    texture,
                    [
                        -load * mu[0] * mu[0] * local[0] / norm,
                        -load * mu[1] * mu[1] * local[1] / norm,
                    ],
                );
                sum[0] += force[0];
                sum[1] += force[1];
                sum[2] += r[0] * force[1] - r[1] * force[0];
                sum[3] += load * norm;
            }
        }
    }
    sum
}

#[test]
fn circular_isotropic_geometry_uses_the_original_exact_owner() {
    let directional = model(3.2, [0.4, 0.4], 0.71, [0.7, 0.7], -0.53);
    let circle = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 3.2,
        contact_radius_m: 0.4,
        friction_coefficient: 0.7,
    })
    .unwrap();
    for (u, w) in [
        ([0.0; 2], 0.0),
        ([0.0; 2], -3.0),
        ([0.7, -0.2], 0.0),
        ([0.2, 0.1], 2.5),
        ([1.7, -0.8], -2.5),
        ([1.0, 0.0], 1e-200),
    ] {
        for panels in [1, 256] {
            assert_eq!(
                directional.evaluate(u, w, grid(panels)).unwrap(),
                circle.evaluate(u, w).unwrap()
            );
        }
    }
    assert_eq!(directional.friction_law().static_coefficients(), [1.0; 2]);
}

#[test]
fn pure_translation_keeps_original_orthotropic_support_and_zero_centre_moment() {
    let patch = model(4.0, [0.9, 0.13], 0.24, [0.2, 0.8], 0.73);
    let law = *patch.friction_law();
    let state = law.initial_state();
    for u in [[0.4, -0.8], [-0.7, 0.2], [0.0; 2]] {
        let point = law
            .advance_kinetic(
                &state,
                &TangentialStepInput {
                    normal_n: 4.0,
                    slip_increment_m: u,
                    dt_s: 1.0,
                    contact_temperature_k: None,
                    nominal_pressure_pa: None,
                    compliance_m_per_n: 0.0,
                    solve: TangentialSolveConfig {
                        max_iterations: 128,
                        maximum_force_n: 1e6,
                        force_absolute_tolerance_n: 1e-12,
                        force_relative_tolerance: 1e-12,
                    },
                },
            )
            .unwrap();
        let result = patch.evaluate(u, 0.0, grid(1)).unwrap();
        for i in 0..2 {
            near(result.force_n[i], -point.traction_n[i], 2e-15);
        }
        assert_eq!(result.spin_moment_nm, 0.0);
        near(result.dissipated_power_w, point.dissipated_j, 2e-15);
        near(result.power_residual_w, 0.0, 3e-15);
    }
    assert_eq!(*patch.friction_law(), law);
    assert_eq!(state, law.initial_state());
}

#[test]
fn nonaligned_geometry_and_texture_match_full_surface_force_moment_and_heat() {
    let patch = model(3.2, [1.1, 0.4], 0.4, [0.2, 0.7], 0.9);
    for (u, w) in [([0.0; 2], 2.0), ([0.6, -0.3], 2.0), ([2.3, -1.7], -0.8)] {
        let actual = patch.evaluate(u, w, grid(128)).unwrap();
        let reference = surface_reference(&patch, u, w, 96, 2048);
        let tolerance = if u == [0.6, -0.3] { 3e-5 } else { 2e-10 };
        for (a, b) in values(actual).into_iter().zip(reference) {
            near(a, b, tolerance);
        }
        near(actual.power_residual_w, 0.0, 4e-14);
        assert!(actual.dissipated_power_w > 0.0);
        if u == [0.0; 2] {
            assert_eq!(actual.force_n, [0.0; 2]);
        }
    }
}

#[test]
fn angular_refinement_resolves_the_original_narrow_directional_response() {
    let patch = model(2.0, [1.3, 0.08], 0.17, [0.09, 0.8], 0.61);
    let u = [0.11, 0.02];
    let reference = values(patch.evaluate(u, 1.0, grid(256)).unwrap());
    let mut errors = Vec::new();
    // The combined footprint and texture anisotropy has a narrow angular
    // response. Resolve that feature before measuring the convergence rate;
    // four or eight panels do not provide a physical error certificate.
    for n in [32, 64, 128] {
        let result = patch.evaluate(u, 1.0, grid(n)).unwrap();
        let error = values(result)
            .into_iter()
            .zip(reference)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push(error);
        near(result.power_residual_w, 0.0, 3e-15);
    }
    eprintln!("Elliptical Hertz angular refinement: {errors:?}");
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < errors[0] / 8.0, "{errors:?}");
    assert!(errors[2] < 1e-7, "{errors:?}");
}

#[test]
fn rigid_rotation_axis_gauge_reversal_and_speed_scale_preserve_the_physical_wrench() {
    let n = 2.3;
    let axes = [0.8, 0.17];
    let shape = 0.23;
    let texture = 0.71;
    let mu = [0.3, 0.8];
    let patch = model(n, axes, shape, mu, texture);
    let u = [0.7, -0.3];
    let w = -2.2;
    let base = patch.evaluate(u, w, grid(64)).unwrap();
    let rotation = 0.83;
    let rotated = model(n, axes, shape + rotation, mu, texture + rotation)
        .evaluate(rotate(direction(rotation), u), w, grid(64))
        .unwrap();
    for (a, b) in rotated
        .force_n
        .into_iter()
        .zip(rotate(direction(rotation), base.force_n))
    {
        near(a, b, 3e-14);
    }
    near(rotated.spin_moment_nm, base.spin_moment_nm, 3e-14);
    near(rotated.dissipated_power_w, base.dissipated_power_w, 3e-14);
    for equivalent in [
        model(n, axes, shape + PI, mu, texture),
        model(n, [axes[1], axes[0]], shape + FRAC_PI_2, mu, texture),
        model(n, axes, shape, [mu[1], mu[0]], texture + FRAC_PI_2),
    ] {
        for (a, b) in values(equivalent.evaluate(u, w, grid(64)).unwrap())
            .into_iter()
            .zip(values(base))
        {
            near(a, b, 3e-14);
        }
    }
    let reverse = patch.evaluate(u.map(|x| -x), -w, grid(64)).unwrap();
    assert_eq!(reverse.force_n, base.force_n.map(|x| -x));
    assert_eq!(reverse.spin_moment_nm, -base.spin_moment_nm);
    assert_eq!(reverse.dissipated_power_w, base.dissipated_power_w);
    for scale in [1e-200, 1e-20, 1e20, 1e200] {
        let scaled = patch
            .evaluate(u.map(|x| x * scale), w * scale, grid(64))
            .unwrap();
        for (a, b) in scaled.force_n.into_iter().zip(base.force_n) {
            near(a, b, 3e-14);
        }
        near(scaled.spin_moment_nm, base.spin_moment_nm, 3e-14);
        near(
            scaled.dissipated_power_w / scale,
            base.dissipated_power_w,
            3e-14,
        );
    }
    assert_eq!(patch.evaluate(u, w, grid(64)).unwrap(), base);
}

#[test]
fn rank_one_local_friction_has_its_actual_hertz_first_moment() {
    let patch = model(7.0, [0.9, 0.13], 0.0, [0.4, 0.0], 0.0);
    let moment: f64 = 3.0 * 7.0 * 0.4 * 0.13 / 8.0;
    for spin in [-3.0, 3.0] {
        let result = patch.evaluate([0.0; 2], spin, grid(16)).unwrap();
        assert_eq!(result.force_n, [0.0; 2]);
        near(result.spin_moment_nm, -moment.copysign(spin), 3e-15);
        near(result.dissipated_power_w, moment * spin.abs(), 5e-15);
    }
    // Nonaligned texture requires splitting a true directional discontinuity.
    let angle = 0.61;
    let rotated_texture = model(7.0, [0.9, 0.13], 0.0, [0.4, 0.0], angle);
    let lever = (0.9 * angle.sin()).hypot(0.13 * angle.cos());
    let result = rotated_texture.evaluate([0.0; 2], 2.0, grid(8)).unwrap();
    near(result.spin_moment_nm, -3.0 * 7.0 * 0.4 * lever / 8.0, 2e-13);
    near(result.power_residual_w, 0.0, 2e-15);
}

#[test]
fn tiny_spin_retains_torque_from_the_original_support_derivative() {
    let patch = model(2.3, [1.3, 0.4], 0.0, [0.4, 0.7], 0.0);
    // At u=[1,0], dr_y/dv_y=mu_y^2/mu_x. Original Hertz <x^2>=a^2/5.
    let coefficient = 2.3 * (0.7 * 0.7 / 0.4) * (1.3 * 1.3 / 5.0);
    for spin in [1e-4, 1e-8, 1e-100, 1e-280] {
        let result = patch.evaluate([1.0, 0.0], spin, grid(32)).unwrap();
        near(-result.spin_moment_nm / spin, coefficient, 2e-8);
        near(result.force_n[0], -0.92, 2e-8);
        near(result.force_n[1], 0.0, 1e-15);
        near(result.power_residual_w, 0.0, 1e-15);
    }
}

#[test]
fn zero_slip_point_crossing_the_ellipse_rim_has_continuous_force_moment_and_heat() {
    let patch = model(3.2, [0.9, 0.13], 0.31, [0.2, 0.7], 0.76);
    // q=A^-1 J u is aligned to the first normalized ellipse axis.
    let r = rotate(direction(0.31), [0.9, 0.0]);
    let u = [r[1], -r[0]];
    let centre = values(patch.evaluate(u, 1.0, grid(128)).unwrap());
    for scale in [1.0 - 1e-9, 1.0 + 1e-9] {
        let result = patch
            .evaluate(u.map(|x| scale * x), 1.0, grid(128))
            .unwrap();
        for (a, b) in values(result).into_iter().zip(centre) {
            near(a, b, 5e-9);
        }
        near(result.power_residual_w, 0.0, 4e-15);
    }
}

#[test]
fn bounded_controls_source_data_rest_and_cancellation_are_explicit() {
    let patch = model(2.0, [0.9, 0.2], 0.1, [0.3, 0.8], 0.7);
    let rest = patch.evaluate([0.0; 2], 0.0, grid(1)).unwrap();
    assert!(rest.at_rest);
    assert_eq!(values(rest), [0.0; 4]);
    for empty in [
        model(0.0, [0.9, 0.2], 0.1, [0.3, 0.8], 0.7),
        model(2.0, [0.9, 0.2], 0.1, [0.0; 2], 0.7),
    ] {
        let result = empty.evaluate([1.0, 2.0], 3.0, grid(1)).unwrap();
        assert_eq!(values(result), [0.0; 4]);
        assert!(!result.at_rest);
    }
    for n in [0, 257, usize::MAX] {
        assert!(patch.evaluate([0.0; 2], 0.0, grid(n)).is_err());
    }
    assert!(patch.evaluate([f64::NAN, 0.0], 1.0, grid(1)).is_err());
    assert!(patch.evaluate([1e300, 0.0], 1e-300, grid(1)).is_err());
    assert!(patch.evaluate([0.0; 2], f64::INFINITY, grid(1)).is_err());
    for (normal, axes, axis) in [
        (-1.0, [1.0; 2], [1.0, 0.0]),
        (1.0, [0.0, 1.0], [1.0, 0.0]),
        (1.0, [1.0; 2], [0.0; 2]),
        (1.0, [1.0, f64::INFINITY], [1.0, 0.0]),
    ] {
        assert!(
            EllipticHertzSlidingPatch::new(
                EllipticHertzPatchParameters {
                    normal_force_n: normal,
                    semi_axes_m: axes,
                    footprint_axis: axis,
                },
                *patch.friction_law()
            )
            .is_err()
        );
    }
    let polls = Cell::new(0);
    let error = patch.evaluate_cancellable([0.2, 0.1], 2.0, grid(32), || {
        polls.set(polls.get() + 1);
        polls.get() >= 7
    });
    assert_eq!(error, Err(EllipticHertzError::Cancelled));
    assert_eq!(polls.get(), 7);
    let exact = model(2.0, [0.9; 2], 0.1, [0.3; 2], 0.7);
    let polls = Cell::new(0);
    assert_eq!(
        exact.evaluate_cancellable([1.0, 0.0], 0.0, grid(1), || {
            polls.set(polls.get() + 1);
            polls.get() == 2
        }),
        Err(EllipticHertzError::Cancelled)
    );
}
