//! G1/G3/G4: independent Hertz moments, original uniform law and source refusal.
use super::*;
use fs_tribo::sliding_patch::elliptic::resolved::{
    ResolvedEllipticHertzPatch, ResolvedHertzError, ResolvedHertzQuadrature,
};
use fs_tribo::sliding_patch::rate_dependent::{HertzFrictionPatch, HertzFrictionQuadrature};
use fs_tribo::{FrictionLaw, TriboError};

fn controls(radial: usize, angular: usize) -> ResolvedHertzQuadrature {
    ResolvedHertzQuadrature {
        radial_panels: radial,
        angular_panels: angular,
        maximum_evaluations: 65_536,
    }
}

fn parameters() -> EllipticHertzPatchParameters {
    EllipticHertzPatchParameters {
        normal_force_n: 7.0,
        semi_axes_m: [0.3, 0.12],
        footprint_axis: direction(0.41),
    }
}

#[test]
fn affine_speed_drag_recovers_independent_elliptical_pressure_moments() {
    let p = parameters();
    let patch = ResolvedEllipticHertzPatch::new(p).unwrap();
    let beta = 0.7;
    for (slip, spin) in [
        ([0.8, -0.3], 0.0),
        ([0.0; 2], -1.7),
        ([0.07, -0.03], 2.3),
        ([0.8, -0.3], 2.3),
        ([0.8, -0.3], 1e-4),
    ] {
        let report = patch
            .evaluate(
                slip,
                spin,
                controls(16, 16),
                |point, _| {
                    let speed = point.slip_velocity_m_s[0].hypot(point.slip_velocity_m_s[1]);
                    EllipticCoulomb::new([10.0; 2], [beta * speed; 2], direction(0.8))
                },
                &mut || false,
            )
            .unwrap();
        // mu=beta|v| makes traction exactly -beta p v. Symmetry gives zero
        // first moment; both semiaxes enter <r^2>=(a^2+b^2)/5.
        let second = (p.semi_axes_m[0].powi(2) + p.semi_axes_m[1].powi(2)) / 5.0;
        for i in 0..2 {
            near(
                report.wrench.force_n[i],
                -beta * p.normal_force_n * slip[i],
                3e-10,
            );
        }
        near(
            report.wrench.spin_moment_nm,
            -beta * p.normal_force_n * spin * second,
            3e-10,
        );
        near(
            report.wrench.dissipated_power_w,
            beta * p.normal_force_n * (slip[0].powi(2) + slip[1].powi(2) + spin * spin * second),
            3e-10,
        );
        near(report.normal_force_residual_n, 0.0, 1e-10);
        near(report.wrench.power_residual_w, 0.0, 3e-15);
    }
}

#[test]
fn rays_cross_the_stationary_rim_without_a_change_in_physical_model() {
    let p = parameters();
    let law = point_law([0.22, 0.71], direction(1.2));
    let native = EllipticHertzSlidingPatch::new(p, law).unwrap();
    let resolved = ResolvedEllipticHertzPatch::new(p).unwrap();
    let spin = 1.7;
    for distance in [0.999, 1.0, 1.001, 2.0] {
        let zero_position = rotate(p.footprint_axis, [p.semi_axes_m[0] * distance, 0.0]);
        let slip = [spin * zero_position[1], -spin * zero_position[0]];
        let expected = native.evaluate(slip, spin, grid(128)).unwrap();
        let report = resolved
            .evaluate(
                slip,
                spin,
                controls(12, 64),
                |_, _| Ok::<_, ()>(law),
                &mut || false,
            )
            .unwrap();
        for (a, b) in values(report.wrench).into_iter().zip(values(expected)) {
            near(a, b, 2e-8);
        }
        near(report.normal_force_residual_n, 0.0, 3e-10);
    }
}

#[test]
fn resolved_uniform_texture_converges_to_the_original_analytic_radial_owner() {
    let p = parameters();
    let law = point_law([0.22, 0.71], direction(1.2));
    let native = EllipticHertzSlidingPatch::new(p, law).unwrap();
    let resolved = ResolvedEllipticHertzPatch::new(p).unwrap();
    let slip = [0.09, -0.04];
    let spin = 1.8;
    let expected = values(native.evaluate(slip, spin, grid(128)).unwrap());
    let mut errors = Vec::new();
    for panels in [4, 8, 16] {
        let report = resolved
            .evaluate(
                slip,
                spin,
                controls(panels, 2 * panels),
                |_, _| Ok::<_, ()>(law),
                &mut || false,
            )
            .unwrap();
        let error = values(report.wrench)
            .into_iter()
            .zip(expected)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
        near(report.wrench.power_residual_w, 0.0, 5e-15);
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1], "{errors:?}");
    assert!(errors[2] < 2e-6, "{errors:?}");
}

#[test]
fn local_stribeck_speed_matches_the_existing_circular_curve_owner() {
    let p = EllipticHertzPatchParameters {
        semi_axes_m: [0.2; 2],
        ..parameters()
    };
    let law = FrictionLaw::Stribeck {
        static_mu: 0.8,
        kinetic_mu: 0.3,
        characteristic_speed: 0.17,
        viscous_per_speed: 0.04,
    };
    let slip = [0.02, 0.045];
    let spin = -1.6;
    let expected = HertzFrictionPatch::new(p.normal_force_n, 0.2, law.clone())
        .unwrap()
        .evaluate(
            slip,
            spin,
            HertzFrictionQuadrature {
                radial_panels: 128,
                angular_panels: 128,
            },
        )
        .unwrap();
    let report = ResolvedEllipticHertzPatch::new(p)
        .unwrap()
        .evaluate(
            slip,
            spin,
            controls(24, 40),
            |point, _| {
                let speed = point.slip_velocity_m_s[0].hypot(point.slip_velocity_m_s[1]);
                // Independent direct Stribeck expression; the reference below
                // invokes the existing shared FrictionLaw owner instead.
                let mu = 0.3 + 0.5 * (-(speed / 0.17).powi(2)).exp() + 0.04 * speed;
                EllipticCoulomb::new([2.0; 2], [mu; 2], [1.0, 0.0])
            },
            &mut || false,
        )
        .unwrap();
    for (actual, expected) in values(report.wrench).into_iter().zip(values(expected)) {
        near(actual, expected, 5e-6);
    }
}

#[test]
fn original_local_axes_and_footprint_rotate_together_without_changing_heat() {
    let p = parameters();
    let slip = [0.04, 0.075];
    let spin = -1.3;
    let rotation = direction(0.73);
    let run = |p, slip, texture| {
        ResolvedEllipticHertzPatch::new(p)
            .unwrap()
            .evaluate(
                slip,
                spin,
                controls(4, 8),
                |point, _| {
                    let speed = point.slip_velocity_m_s[0].hypot(point.slip_velocity_m_s[1]);
                    EllipticCoulomb::new(
                        [1.0; 2],
                        [0.2 + 0.1 * speed, 0.5 / (1.0 + speed)],
                        texture,
                    )
                },
                &mut || false,
            )
            .unwrap()
    };
    let a = run(p, slip, direction(0.2));
    let b = run(
        EllipticHertzPatchParameters {
            footprint_axis: rotate(rotation, p.footprint_axis),
            ..p
        },
        rotate(rotation, slip),
        rotate(rotation, direction(0.2)),
    );
    let force = rotate(rotation, a.wrench.force_n);
    for i in 0..2 {
        near(b.wrench.force_n[i], force[i], 2e-14);
    }
    near(b.wrench.spin_moment_nm, a.wrench.spin_moment_nm, 2e-14);
    near(
        b.wrench.dissipated_power_w,
        a.wrench.dissipated_power_w,
        2e-14,
    );
    assert_eq!(a.evaluations, b.evaluations);
}

#[test]
fn controls_source_and_late_cancel_refuse_before_any_partial_wrench_is_returned() {
    let p = ResolvedEllipticHertzPatch::new(parameters()).unwrap();
    let calls = Cell::new(0);
    let sample = |_, _: &mut fn() -> bool| {
        calls.set(calls.get() + 1);
        Ok::<_, TriboError>(point_law([0.2, 0.4], [1.0, 0.0]))
    };
    let mut c = controls(4, 8);
    c.maximum_evaluations = 1;
    assert!(matches!(
        p.evaluate([0.1, 0.03], 2.0, c, sample, &mut (never as fn() -> bool)),
        Err(ResolvedHertzError::EvaluationBudget { .. })
    ));
    assert_eq!(calls.get(), 0);
    let source_error = p.evaluate(
        [0.1, 0.03],
        2.0,
        controls(2, 2),
        |_, _| {
            calls.set(calls.get() + 1);
            if calls.get() == 7 {
                Err("late local source")
            } else {
                Ok(point_law([0.2, 0.4], [1.0, 0.0]))
            }
        },
        &mut || false,
    );
    assert!(matches!(
        source_error,
        Err(ResolvedHertzError::Source("late local source"))
    ));
    calls.set(0);
    assert!(matches!(
        p.evaluate(
            [0.1, 0.03],
            2.0,
            controls(2, 2),
            |_, _| {
                calls.set(calls.get() + 1);
                Ok::<_, ()>(point_law([0.2, 0.4], [1.0, 0.0]))
            },
            &mut || calls.get() == 7
        ),
        Err(ResolvedHertzError::Cancelled)
    ));
    calls.set(0);
    let a = p
        .evaluate(
            [0.1, 0.03],
            2.0,
            controls(2, 2),
            sample,
            &mut (never as fn() -> bool),
        )
        .unwrap();
    calls.set(0);
    let b = p
        .evaluate(
            [0.1, 0.03],
            2.0,
            controls(2, 2),
            sample,
            &mut (never as fn() -> bool),
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(calls.get(), a.evaluations);
}

fn never() -> bool {
    false
}
