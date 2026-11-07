//! Independent free-point force, trajectory and work solutions for directional
//! sustained friction. The original complete mass is two kg in every axis.

use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_solid::dynamics::sustained::elliptic::{
    SolidSustainedEllipticContact, SolidSustainedEllipticPath,
};

fn elliptic(h: f64, mu: [f64; 2], axis: [f64; 2]) -> SolidSustainedEllipticContact<'static> {
    let original = contact(h, mu[0].max(mu[1]), [0.0; 3]);
    SolidSustainedEllipticContact {
        id: original.id,
        velocity_weights: original.velocity_weights,
        gap_m: original.gap_m,
        prescribed_surface_velocity_mps: original.prescribed_surface_velocity_mps,
        friction: EllipticFriction {
            coefficients_mu: mu,
            principal_axis: axis,
        },
        scales: original.scales,
    }
}

// The native free-point trajectory has quadratic gap and affine velocity.
// Its exact interior gap minimum is checked. For a fixed force, the ellipse's
// normal cone is convex: support at both velocity endpoints therefore admits
// every intervening affine velocity, including a possible sticking endpoint.
fn admit(path: &SolidSustainedEllipticPath<'_, '_>) -> Result<&'static str, &'static str> {
    let h = path.duration_s;
    let acceleration: Vec<_> = path.total_force_n.iter().map(|f| f / 2.0).collect();
    for i in 0..3 {
        near(
            path.velocity_after_mps[i],
            path.before.velocity_mps()[i] + h * acceleration[i],
            1e-10,
        );
        near(
            path.displacement_change_m[i],
            h * path.before.velocity_mps()[i] + 0.5 * h * h * acceleration[i],
            1e-11,
        );
    }
    for (index, (c, f)) in path.contacts.iter().zip(path.forces_n).enumerate() {
        let v0: [f64; 3] = core::array::from_fn(|a| {
            dot(c.velocity_weights[a], path.before.velocity_mps())
                - c.prescribed_surface_velocity_mps[a]
        });
        let a: [f64; 3] = core::array::from_fn(|axis| dot(c.velocity_weights[axis], &acceleration));
        let gap = |t: f64| c.gap_m + v0[0] * t + 0.5 * a[0] * t * t;
        let mut minimum = gap(0.0).min(gap(h));
        if a[0] > 0.0 {
            let t = -v0[0] / a[0];
            if t > 0.0 && t < h {
                minimum = minimum.min(gap(t));
            }
        }
        if minimum < -1e-11 {
            return Err("interior penetration");
        }
        if f[0] > 1e-7 && (v0[0].abs() > 1e-9 || a[0].abs() > 1e-8 || c.gap_m > 1e-10) {
            return Err("loaded point does not dwell");
        }
        for axis in 0..3 {
            near(
                path.paths_m[index][axis],
                v0[axis] * h + 0.5 * a[axis] * h * h + if axis == 0 { c.gap_m } else { 0.0 },
                1e-11,
            );
        }
        let axis = c.friction.principal_axis;
        let norm = axis[0].hypot(axis[1]);
        let u = [axis[0] / norm, axis[1] / norm];
        let mu = c.friction.coefficients_mu;
        let fp = [u[0] * f[1] + u[1] * f[2], -u[1] * f[1] + u[0] * f[2]];
        let inside = if mu[0] == 0.0 && mu[1] == 0.0 {
            near(f[1].hypot(f[2]), 0.0, 1e-8);
            false
        } else {
            let gauge = (fp[0] / mu[0]).hypot(fp[1] / mu[1]);
            if gauge > f[0] + 1e-7 {
                return Err("physical force exceeds ellipse");
            }
            gauge < f[0] - 1e-7
        };
        for t in [0.0, h] {
            let v = [v0[1] + t * a[1], v0[2] + t * a[2]];
            let vp = [u[0] * v[0] + u[1] * v[1], -u[1] * v[0] + u[0] * v[1]];
            let support = f[0] * (mu[0] * vp[0]).hypot(mu[1] * vp[1]);
            let duality = support + f[1] * v[0] + f[2] * v[1];
            if duality.abs() > 1e-7 || (inside && v[0].hypot(v[1]) > 1e-9) {
                return Err("physical elliptic mode changes inside interval");
            }
        }
    }
    Ok("analytic continuous elliptic point path")
}

#[test]
fn principal_axes_have_different_breakaway_and_sliding_work() {
    with_cx(|cx, _| {
        let h = 0.125;
        for (velocity, expected_friction) in [
            ([2.0, 0.0, 0.0], [-12.0, 0.0]),
            ([0.0, 0.0, 2.0], [0.0, -4.0]),
        ] {
            let mut s = point(h, velocity, cx);
            let r = s
                .advance_sustained_elliptic_contacts(
                    &[elliptic(h, [0.6, 0.2], [1.0, 0.0])],
                    &[[20.0, 0.0, 0.0]],
                    &[0.0, -20.0, 0.0],
                    controls(),
                    cx,
                    admit,
                )
                .unwrap();
            near(r.forces_n[0][0], 20.0, 1e-8);
            let mut work = 0.0;
            for (axis, force) in [0, 2].into_iter().zip(expected_friction) {
                let displacement = h * velocity[axis] + h * h * force / 4.0;
                near(
                    s.velocity_mps()[axis],
                    velocity[axis] + h * force / 2.0,
                    1e-9,
                );
                near(s.displacement_m()[axis], displacement, 1e-10);
                near(r.forces_n[0][if axis == 0 { 1 } else { 2 }], force, 1e-8);
                work += force * displacement;
            }
            near(r.transformed.frictional_work_j, work, 1e-9);
            assert!(work < 0.0);
            energy(&r.transformed);
        }
        // A load inside the ellipse sticks. A six-newton load in the weak
        // direction exceeds its four-newton capacity despite being below the
        // twelve-newton capacity in the strong direction.
        for (load, expected_force, acceleration) in [
            ([9.0, 2.0], [-9.0, -2.0], [0.0, 0.0]),
            ([0.0, 6.0], [0.0, -4.0], [0.0, 1.0]),
        ] {
            let mut s = point(h, [0.0; 3], cx);
            let r = s
                .advance_sustained_elliptic_contacts(
                    &[elliptic(h, [0.6, 0.2], [1.0, 0.0])],
                    &[[0.0; 3]],
                    &[load[0], -20.0, load[1]],
                    controls(),
                    cx,
                    admit,
                )
                .unwrap();
            for k in 0..2 {
                near(r.forces_n[0][k + 1], expected_force[k], 1e-8);
                near(s.velocity_mps()[2 * k], h * acceleration[k], 1e-9);
                near(
                    s.displacement_m()[2 * k],
                    0.5 * h * h * acceleration[k],
                    1e-10,
                );
            }
            energy(&r.transformed);
        }
    });
}

#[test]
fn rotated_oblique_ellipse_uses_maximum_dissipation_not_opposite_velocity() {
    with_cx(|cx, _| {
        let h = 0.125;
        let velocity = [1.2, 0.0, 0.9];
        let mu: [f64; 2] = [0.8, 0.2];
        let u = [0.6, 0.8];
        let vp = [
            u[0] * velocity[0] + u[1] * velocity[2],
            -u[1] * velocity[0] + u[0] * velocity[2],
        ];
        let support = (mu[0] * vp[0]).hypot(mu[1] * vp[1]);
        let fp = [
            -20.0 * mu[0] * mu[0] * vp[0] / support,
            -20.0 * mu[1] * mu[1] * vp[1] / support,
        ];
        let expected = [u[0] * fp[0] - u[1] * fp[1], u[1] * fp[0] + u[0] * fp[1]];
        let mut s = point(h, velocity, cx);
        let r = s
            .advance_sustained_elliptic_contacts(
                &[elliptic(h, mu, [3.0, 4.0])],
                &[[20.0, 0.0, 0.0]],
                &[-expected[0], -20.0, -expected[1]],
                controls(),
                cx,
                admit,
            )
            .unwrap();
        for k in 0..2 {
            near(r.forces_n[0][k + 1], expected[k], 2e-8);
            near(s.velocity_mps()[2 * k], velocity[2 * k], 1e-9);
        }
        assert!((expected[0] * velocity[2] - expected[1] * velocity[0]).abs() > 1.0);
        near(
            r.transformed.frictional_dissipation_j,
            20.0 * support * h,
            1e-9,
        );
        energy(&r.transformed);
    });
}

#[test]
fn moving_surface_work_and_motion_are_covariant_under_tangent_rotation() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut results = Vec::new();
        for angle in [0.0_f64, 0.73] {
            let (sin, cos) = angle.sin_cos();
            let t1 = [cos, 0.0, sin];
            let t2 = [-sin, 0.0, cos];
            let mut c = elliptic(h, [0.6, 0.2], [cos, -sin]);
            c.velocity_weights = [&N, &t1, &t2];
            c.prescribed_surface_velocity_mps = [0.0, cos, -sin];
            let mut s = point(h, [0.0; 3], cx);
            let r = s
                .advance_sustained_elliptic_contacts(
                    &[c],
                    &[[20.0, 0.0, 0.0]],
                    &[0.0, -20.0, 0.0],
                    controls(),
                    cx,
                    admit,
                )
                .unwrap();
            near(r.forces_n[0][1], 12.0 * cos, 1e-8);
            near(r.forces_n[0][2], -12.0 * sin, 1e-8);
            near(s.velocity_mps()[0], 0.75, 1e-9);
            near(s.velocity_mps()[2], 0.0, 1e-9);
            near(s.displacement_m()[0], 0.046875, 1e-10);
            near(r.transformed.prescribed_surface_work_j, 1.5, 1e-9);
            near(r.transformed.frictional_dissipation_j, 0.9375, 1e-9);
            near(r.transformed.native.energy.kinetic_after_j, 0.5625, 1e-9);
            energy(&r.transformed);
            results.push((s.displacement_m().to_vec(), s.velocity_mps().to_vec()));
        }
        for (a, b) in results[0]
            .0
            .iter()
            .chain(&results[0].1)
            .zip(results[1].0.iter().chain(&results[1].1))
        {
            near(*a, *b, 1e-10);
        }
    });
}

#[test]
fn equal_coefficients_preserve_the_existing_isotropic_response() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut disk = point(h, [1.5, 0.0, 2.0], cx);
        let mut ellipse = point(h, [1.5, 0.0, 2.0], cx);
        let r = disk
            .advance_sustained_contacts(
                &[contact(h, 0.5, [0.0; 3])],
                &[[20.0, 0.0, 0.0]],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                admit_point,
            )
            .unwrap();
        let e = ellipse
            .advance_sustained_elliptic_contacts(
                &[elliptic(h, [0.5, 0.5], [3.0, 4.0])],
                &[[20.0, 0.0, 0.0]],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                admit,
            )
            .unwrap();
        assert_eq!(bits(&disk), bits(&ellipse));
        assert_eq!(r.solve.forces_n, e.forces_n);
        assert_eq!(r.solve.paths_m, e.paths_m);
        assert_eq!(
            r.solve.response_evaluations,
            e.transformed.solve.response_evaluations
        );
    });
}

#[test]
fn physical_path_refusal_singular_law_and_staging_budget_leave_native_state_unchanged() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut s = point(h, [2.0, 0.0, 0.0], cx);
        let before = bits(&s);
        let c = elliptic(h, [0.6, 0.2], [1.0, 0.0]);
        let refused = s.advance_sustained_elliptic_contacts(
            &[c],
            &[[20.0, 0.0, 0.0]],
            &[0.0, -20.0, 0.0],
            controls(),
            cx,
            |_| Err::<(), _>("unresolved contact feature"),
        );
        assert!(matches!(
            refused,
            Err(SolidSustainedError::Path("unresolved contact feature"))
        ));
        assert_eq!(before, bits(&s));
        let singular = elliptic(h, [0.6, 0.0], [1.0, 0.0]);
        assert!(matches!(
            s.advance_sustained_elliptic_contacts(
                &[singular],
                &[[20.0, 0.0, 0.0]],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                admit
            ),
            Err(SolidSustainedError::Invalid(_))
        ));
        let mut small = controls();
        small.maximum_workspace_bytes = 1;
        assert!(matches!(
            s.advance_sustained_elliptic_contacts(
                &[c],
                &[[20.0, 0.0, 0.0]],
                &[0.0, -20.0, 0.0],
                small,
                cx,
                admit
            ),
            Err(SolidSustainedError::Budget { .. })
        ));
        assert_eq!(before, bits(&s));
        assert_eq!(s.accepted_steps(), 0);
    });
}
