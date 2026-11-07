//! G1/G3/G4: actual directional hysteresis on the complete native force path.

use super::*;

#[test]
fn static_capture_and_nonproportional_kinetic_breakaway_follow_actual_motion() {
    with_cx(|cx, _| {
        let h = 0.125;
        let c = elliptic(h, [0.6, 0.2], [1.0, 0.0]);
        let kinetic = [[0.3, 0.05]];
        // The static ellipse admits (9,2) N, but the smaller kinetic ellipse
        // cannot hold it. An actually resting point retains static strength.
        let mut held = point(h, [0.0; 3], cx);
        let r = held
            .advance_sustained_static_kinetic_elliptic_contacts(
                &[c],
                &kinetic,
                &[[0.0; 3]],
                &[9.0, -20.0, 2.0],
                controls(),
                cx,
                admit,
            )
            .unwrap();
        assert_eq!(r.friction[0].coefficients_mu, [0.6, 0.2]);
        assert_eq!(r.transformed.friction_regime_resolves, 0);
        near(r.forces_n[0][1], -9.0, 1e-8);
        near(r.forces_n[0][2], -2.0, 1e-8);
        for v in held.velocity_mps() {
            near(*v, 0.0, 1e-9);
        }
        energy(&r.transformed);

        let mut sliding = point(h, [0.0; 3], cx);
        let r = sliding
            .advance_sustained_static_kinetic_elliptic_contacts(
                &[c],
                &kinetic,
                &[[0.0; 3]],
                &[0.0, -20.0, 6.0],
                controls(),
                cx,
                admit,
            )
            .unwrap();
        assert_eq!(r.friction[0].coefficients_mu, kinetic[0]);
        assert_eq!(r.transformed.friction_regime_resolves, 1);
        near(r.forces_n[0][2], -1.0, 1e-8);
        near(sliding.velocity_mps()[2], 2.5 * h, 1e-9);
        near(sliding.displacement_m()[2], 1.25 * h * h, 1e-10);
        near(r.transformed.frictional_dissipation_j, 1.25 * h * h, 1e-10);
        energy(&r.transformed);
        // Once moving, a load inside static strength but outside kinetic
        // strength must keep sliding at the smaller directional bound.
        let v = sliding.velocity_mps()[2];
        let r = sliding
            .advance_sustained_static_kinetic_elliptic_contacts(
                &[c],
                &kinetic,
                &r.forces_n,
                &[0.0, -20.0, 2.0],
                controls(),
                cx,
                admit,
            )
            .unwrap();
        assert_eq!(r.transformed.friction_regime_resolves, 0);
        assert_eq!(r.friction[0].coefficients_mu, kinetic[0]);
        near(sliding.velocity_mps()[2], v + 0.5 * h, 1e-9);
        near(r.forces_n[0][2], -1.0, 1e-8);
        energy(&r.transformed);
    });
}

#[test]
fn incoming_oblique_slip_uses_the_true_kinetic_ellipse_on_the_same_axes() {
    with_cx(|cx, _| {
        let h = 0.125;
        let velocity = [1.2, 0.0, 0.9];
        let mu: [f64; 2] = [0.3, 0.05];
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
        let force = [u[0] * fp[0] - u[1] * fp[1], u[1] * fp[0] + u[0] * fp[1]];
        let mut body = point(h, velocity, cx);
        let r = body
            .advance_sustained_static_kinetic_elliptic_contacts(
                &[elliptic(h, [0.8, 0.2], [3.0, 4.0])],
                &[mu],
                &[[20.0, 0.0, 0.0]],
                &[-force[0], -20.0, -force[1]],
                controls(),
                cx,
                admit,
            )
            .unwrap();
        assert_eq!(r.transformed.friction_regime_resolves, 0);
        assert_eq!(r.friction[0].principal_axis, [3.0, 4.0]);
        assert_eq!(r.friction[0].coefficients_mu, mu);
        for k in 0..2 {
            near(r.forces_n[0][k + 1], force[k], 2e-8);
            near(body.velocity_mps()[2 * k], velocity[2 * k], 1e-9);
        }
        assert!((force[0] * velocity[2] - force[1] * velocity[0]).abs() > 0.5);
        near(
            r.transformed.frictional_dissipation_j,
            20.0 * support * h,
            1e-9,
        );
        energy(&r.transformed);
    });
}

#[test]
fn isotropic_static_kinetic_pairs_recover_the_original_circular_path_bitwise() {
    with_cx(|cx, _| {
        let h = 0.125;
        for (velocity, fx) in [([0.0; 3], 9.0), ([0.0; 3], 16.0), ([2.0, 0.0, 0.0], 9.0)] {
            let mut circular = point(h, velocity, cx);
            let mut directional = point(h, velocity, cx);
            let load = [fx, -20.0, 0.0];
            let a = circular
                .advance_sustained_static_kinetic_contacts(
                    &[contact(h, 0.6, [0.0; 3])],
                    &[0.3],
                    &[[20.0, 0.0, 0.0]],
                    &load,
                    controls(),
                    cx,
                    admit_point,
                )
                .unwrap();
            let b = directional
                .advance_sustained_static_kinetic_elliptic_contacts(
                    &[elliptic(h, [0.6; 2], [3.0, 4.0])],
                    &[[0.3; 2]],
                    &[[20.0, 0.0, 0.0]],
                    &load,
                    controls(),
                    cx,
                    admit,
                )
                .unwrap();
            assert_eq!(bits(&directional), bits(&circular));
            assert_eq!(b.forces_n, a.solve.forces_n);
            assert_eq!(
                b.transformed.friction_regime_resolves,
                a.friction_regime_resolves
            );
            assert_eq!(
                b.transformed.solve.response_evaluations,
                a.solve.response_evaluations
            );
            assert_eq!(
                b.transformed.frictional_work_j.to_bits(),
                a.frictional_work_j.to_bits()
            );
        }
    });
}

#[test]
fn coupled_breakaway_resolves_normal_shear_and_the_other_contacts_static_reaction() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut mass = [0.0; 36];
        for i in 0..6 {
            mass[i * 6 + i] = 2.0;
        }
        mass[1] = 0.25;
        mass[6] = 0.25;
        mass[3] = 0.5;
        mass[18] = 0.5;
        let mut body = ExactLinearSolid::prepare(
            &assembly(&[0.0; 36], &mass, 6),
            h,
            &[0.0; 6],
            &[0.0; 6],
            native_controls(),
            cx,
        )
        .unwrap();
        let x0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let y0 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let z0 = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let x1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let y1 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let z1 = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let a = SolidSustainedEllipticContact {
            id: 10,
            velocity_weights: [&y0, &x0, &z0],
            ..elliptic(h, [0.6, 0.2], [1.0, 0.0])
        };
        let b = SolidSustainedEllipticContact {
            id: 20,
            velocity_weights: [&y1, &x1, &z1],
            ..elliptic(h, [0.5, 0.4], [1.0, 0.0])
        };
        // Normal/shear coupling gives N_A=20+.25*a; the native row
        // 2*a=14-.3*N_A fixes a independently. Keeping point B at rest
        // then needs f_B=.5*a-10, within its STATIC ellipse.
        let acceleration = 8.0 / 2.075;
        let normal = 20.0 + 0.25 * acceleration;
        let friction_a = -0.3 * normal;
        let friction_b = 0.5 * acceleration - 10.0;
        let calls = core::cell::Cell::new(0);
        let r = body
            .advance_sustained_static_kinetic_elliptic_contacts(
                &[b, a],
                &[[0.2, 0.1], [0.3, 0.05]],
                &[[20.0, 0.0, 0.0]; 2],
                &[14.0, -20.0, 0.0, 10.0, -20.0, 0.0],
                controls(),
                cx,
                |p| {
                    calls.set(calls.get() + 1);
                    assert_eq!(p.before.time_s(), 0.0);
                    near(p.velocity_after_mps[0], acceleration * h, 1e-8);
                    for axis in 1..6 {
                        near(p.velocity_after_mps[axis], 0.0, 1e-8);
                    }
                    assert!(friction_b.abs() < 0.5 * 20.0);
                    Ok::<_, &'static str>(())
                },
            )
            .unwrap();
        assert_eq!(calls.get(), 1);
        assert_eq!(r.transformed.contact_ids, [10, 20]);
        assert_eq!(r.transformed.friction_regime_resolves, 1);
        assert_eq!(r.friction[0].coefficients_mu, [0.3, 0.05]);
        assert_eq!(r.friction[1].coefficients_mu, [0.5, 0.4]);
        near(r.forces_n[0][0], normal, 1e-8);
        near(r.forces_n[0][1], friction_a, 1e-8);
        near(r.forces_n[1][1], friction_b, 1e-8);
        near(body.displacement_m()[0], 0.5 * acceleration * h * h, 1e-10);
        near(
            r.transformed.momentum_change_kg_m_s[1],
            0.25 * acceleration * h,
            1e-9,
        );
        near(
            r.transformed.momentum_change_kg_m_s[3],
            0.5 * acceleration * h,
            1e-9,
        );
        near(
            r.transformed.frictional_work_j,
            friction_a * 0.5 * acceleration * h * h,
            1e-9,
        );
        energy(&r.transformed);
    });
}

#[test]
fn regime_budgets_late_refusals_and_cancellation_preserve_the_entire_predecessor() {
    with_cx(|cx, gate| {
        let h = 0.125;
        let contacts = [elliptic(h, [0.6, 0.2], [1.0, 0.0])];
        let kinetic = [[0.3, 0.05]];
        let initial = [[20.0, 0.0, 0.0]];
        let load = [0.0, -20.0, 6.0];
        let mut reference = point(h, [0.0; 3], cx);
        let accepted = reference
            .advance_sustained_static_kinetic_elliptic_contacts(
                &contacts,
                &kinetic,
                &initial,
                &load,
                controls(),
                cx,
                admit,
            )
            .unwrap();
        let mut body = point(h, [0.0; 3], cx);
        let before = bits(&body);
        let mut spent = controls();
        spent.solve.maximum_response_evaluations =
            accepted.transformed.solve.response_evaluations - 1;
        assert!(
            body.advance_sustained_static_kinetic_elliptic_contacts(
                &contacts, &kinetic, &initial, &load, spent, cx, admit,
            )
            .is_err()
        );
        assert_eq!(bits(&body), before);
        for bad in [[[0.7, 0.05]], [[0.3, -0.1]], [[0.3, 0.0]]] {
            assert!(
                body.advance_sustained_static_kinetic_elliptic_contacts(
                    &contacts,
                    &bad,
                    &initial,
                    &load,
                    controls(),
                    cx,
                    admit,
                )
                .is_err()
            );
            assert_eq!(bits(&body), before);
        }
        assert!(matches!(
            body.advance_sustained_static_kinetic_elliptic_contacts(
                &contacts,
                &kinetic,
                &initial,
                &load,
                controls(),
                cx,
                |_| Err::<(), _>("changed geometry after breakaway"),
            ),
            Err(SolidSustainedError::Path(
                "changed geometry after breakaway"
            ))
        ));
        assert_eq!(bits(&body), before);
        assert!(matches!(
            body.advance_sustained_static_kinetic_elliptic_contacts(
                &contacts,
                &kinetic,
                &initial,
                &load,
                controls(),
                cx,
                |p| {
                    admit(p)?;
                    gate.request();
                    Ok::<_, &'static str>(())
                },
            ),
            Err(SolidSustainedError::Cancelled)
        ));
        assert_eq!(bits(&body), before);
        with_cx(|fresh, _| {
            let mut exact_budget = controls();
            exact_budget.solve.maximum_response_evaluations =
                accepted.transformed.solve.response_evaluations;
            body.advance_sustained_static_kinetic_elliptic_contacts(
                &contacts,
                &kinetic,
                &initial,
                &load,
                exact_budget,
                fresh,
                admit,
            )
            .unwrap();
            assert_eq!(bits(&body), bits(&reference));
        });
    });
}
