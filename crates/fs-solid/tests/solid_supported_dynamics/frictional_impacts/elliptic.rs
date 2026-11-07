//! G1/G3/G4 directional supported events, source transactions and native continuation.

#[path = "elliptic/source.rs"]
mod source;

use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_solid::dynamics::frictional_impact::elliptic::SolidEllipticContact;
use fs_solid::dynamics::frictional_impact::elliptic::source::{
    SolidEllipticSourceError, SolidImpactGeometry,
};
use fs_solid::dynamics::supported::SupportedSolidEllipticSourceError;

fn directional(
    c: SolidFrictionalContact<'_>,
    mu: [f64; 2],
    axis: [f64; 2],
) -> SolidEllipticContact<'_> {
    SolidEllipticContact {
        id: c.id,
        velocity_weights: c.velocity_weights,
        gap_m: c.gap_m,
        activation_tolerance_m: c.activation_tolerance_m,
        prescribed_surface_velocity_mps: c.prescribed_surface_velocity_mps,
        normal_law: c.normal_law,
        friction_coefficient: EllipticFriction {
            coefficients_mu: mu,
            principal_axis: axis,
        },
    }
}

fn geometry(c: SolidFrictionalContact<'_>) -> SolidImpactGeometry<'_> {
    SolidImpactGeometry {
        id: c.id,
        velocity_weights: c.velocity_weights,
        gap_m: c.gap_m,
        activation_tolerance_m: c.activation_tolerance_m,
        prescribed_surface_velocity_mps: c.prescribed_surface_velocity_mps,
        normal_law: c.normal_law,
        friction_coefficient: (),
    }
}

// Independent physical support-gradient oracle. No production cone transform.
fn support_impulse(pn: f64, slip: [f64; 2], mu: [f64; 2], axis: [f64; 2]) -> [f64; 3] {
    let [c, s] = axis;
    let v = [c * slip[0] + s * slip[1], -s * slip[0] + c * slip[1]];
    let denominator = (mu[0] * v[0]).hypot(mu[1] * v[1]);
    let p = [
        -pn * mu[0] * mu[0] * v[0] / denominator,
        -pn * mu[1] * mu[1] * v[1] / denominator,
    ];
    [pn, c * p[0] - s * p[1], s * p[0] + c * p[1]]
}

#[test]
fn isotropic_limit_keeps_supported_response_work_and_continuation_bit_identical() {
    with_cx(false, |cx| {
        for mu in [0.0, 0.25, 4.0] {
            let c = contact(mu, [0.1, -0.2, 0.3]);
            let e = directional(c, [mu; 2], [3e200, 4e200]);
            let initial = [-1.0, 2.4, 0.0, 0.4];
            let (mut circle, mut ellipse) = (
                make(initial, 1.0, controls(), cx),
                make(initial, 1.0, controls(), cx),
            );
            let a = circle
                .apply_frictional_impacts(&[c], impact_controls(), cx)
                .unwrap();
            let b = ellipse
                .apply_elliptic_impacts(&[e], impact_controls(), cx)
                .unwrap();
            assert_eq!(bits(&circle), bits(&ellipse));
            assert_eq!(a.contact.law.impulses, b.contact.law.impulses);
            assert_eq!(
                a.impulse.support_reaction_kg_m_s,
                b.impulse.support_reaction_kg_m_s
            );
            assert_eq!(a.impulse.energy, b.impulse.energy);
            assert_eq!(
                a.contact.relative_work_j.to_bits(),
                b.contact.relative_work_j.to_bits()
            );
            circle
                .advance_affine_supports(&[0.0], &[0.0; 4], &[0.0; 4], cx)
                .unwrap();
            ellipse
                .advance_affine_supports(&[0.0], &[0.0; 4], &[0.0; 4], cx)
                .unwrap();
            assert_eq!(bits(&circle), bits(&ellipse));
        }
    });
}

#[test]
fn oblique_material_friction_retains_direct_support_load_cross_mass_work_and_motion() {
    with_cx(false, |cx| {
        let (mu, axis, slip) = ([0.25, 0.6], [0.6, 0.8], [0.4, 0.3]);
        let p = support_impulse(3.0, slip, mu, axis);
        assert!((p[1] * slip[1] - p[2] * slip[0]).abs() > 0.1);
        let mut outputs = Vec::new();
        for boost in [0.0, 2.0] {
            let (surface, vc) = ([0.1 + boost, -0.2, 0.3 + boost], 0.4 + boost);
            let after = [
                surface[0] + 0.5,
                surface[1] + vc + slip[0],
                surface[2] + slip[1],
                vc,
            ];
            let delta = [p[0] / 2.0, p[1] / 3.0, p[2] / 4.0];
            let before = [
                after[0] - delta[0],
                after[1] - delta[1],
                after[2] - delta[2],
                vc,
            ];
            let mut body = make(before, 1.0, controls(), cx);
            let r = body
                .apply_elliptic_impacts(
                    &[directional(contact(0.0, surface), mu, axis)],
                    impact_controls(),
                    cx,
                )
                .unwrap();
            for i in 0..4 {
                near(body.velocity_mps()[i], after[i], 3e-10);
            }
            for i in 0..3 {
                near(r.contact.law.impulses[0][i], p[i], 3e-10);
            }
            let support_momentum = 0.35 * delta[0] - 0.2 * delta[1] + 0.1 * delta[2];
            let reaction = support_momentum + p[1]; // J_c = -p_t1.
            let relative_work = p[0] * (0.5 - 0.5 * delta[0])
                + p[1] * (slip[0] - 0.5 * delta[1])
                + p[2] * (slip[1] - 0.5 * delta[2]);
            let surface_work: f64 = p.iter().zip(surface).map(|(p, v)| p * v).sum();
            near(r.impulse.momentum_change_kg_m_s[3], support_momentum, 4e-10);
            near(r.impulse.support_reaction_kg_m_s[0], reaction, 4e-10);
            near(r.impulse.energy.support_work_j, reaction * vc, 4e-10);
            near(r.contact.relative_work_j, relative_work, 4e-10);
            near(r.contact.prescribed_surface_work_j, surface_work, 4e-10);
            near(
                r.impulse.energy.applied_work_j,
                relative_work + surface_work,
                5e-10,
            );
            near(
                kinetic(&after) - kinetic(&before),
                relative_work + surface_work + reaction * vc,
                5e-10,
            );
            assert!(r.contact.measured_residuals.sliding_impulse_residual < 1e-10);
            assert_eq!(body.displacement_m(), &U);
            assert_eq!(body.time_s(), 0.0);
            assert_eq!(body.accepted_steps(), 0);
            body.advance_affine_supports(&[0.0], &[0.0; 4], &[0.0; 4], cx)
                .unwrap();
            let (s, c) = 0.03_f64.sin_cos();
            for i in 0..3 {
                near(
                    body.displacement_m()[i],
                    U[3] + 0.01 * vc + (U[i] - U[3]) * c + (after[i] - vc) * s / 3.0,
                    4e-11,
                );
                near(
                    body.velocity_mps()[i],
                    vc - 3.0 * (U[i] - U[3]) * s + (after[i] - vc) * c,
                    4e-10,
                );
            }
            outputs.push(r);
        }
        near(
            outputs[0].contact.relative_work_j,
            outputs[1].contact.relative_work_j,
            5e-10,
        );
        near(
            outputs[1].contact.prescribed_surface_work_j
                - outputs[0].contact.prescribed_surface_work_j,
            2.0 * (p[0] + p[2]),
            5e-10,
        );
    });
}

#[test]
fn exact_free_material_axis_is_covariant_under_contact_frame_rotation() {
    with_cx(false, |cx| {
        let surface = [0.1, -0.2, 0.3];
        let before = [-0.9, 2.4, 1.1, 0.4];
        let rows = [ROWS[0], ROWS[2].map(|v| -v), ROWS[1]];
        let mut a = make(before, 1.0, controls(), cx);
        let mut b = make(before, 1.0, controls(), cx);
        let original = directional(contact(0.0, surface), [0.4, 0.0], [1.0, 0.0]);
        let rotated = SolidEllipticContact {
            velocity_weights: [&rows[0], &rows[1], &rows[2]],
            prescribed_surface_velocity_mps: [surface[0], -surface[2], surface[1]],
            friction_coefficient: EllipticFriction {
                coefficients_mu: [0.4, 0.0],
                principal_axis: [0.0, 1.0],
            },
            ..original
        };
        let ra = a
            .apply_elliptic_impacts(&[original], impact_controls(), cx)
            .unwrap();
        let rb = b
            .apply_elliptic_impacts(&[rotated], impact_controls(), cx)
            .unwrap();
        near(ra.contact.law.impulses[0][0], 3.0, 2e-10);
        near(ra.contact.law.impulses[0][1], -1.2, 2e-10);
        assert_eq!(ra.contact.law.impulses[0][2], 0.0);
        assert_eq!(a.velocity_mps()[2], before[2]);
        for i in 0..4 {
            near(a.velocity_mps()[i], b.velocity_mps()[i], 2e-10);
        }
        near(rb.contact.law.impulses[0][2], -1.2, 2e-10);
        assert_eq!(rb.contact.law.impulses[0][1], 0.0);
        near(
            ra.impulse.energy.support_work_j,
            rb.impulse.energy.support_work_j,
            2e-10,
        );
        near(
            ra.contact.relative_work_j,
            rb.contact.relative_work_j,
            2e-10,
        );
    });
}

#[test]
fn simultaneous_stick_slide_keeps_cross_contact_mobility_and_all_elastic_modes() {
    with_cx(false, |cx| {
        let mut phi = [[0.0; 6]; 6];
        let mut inverse = phi;
        for i in 0..6 {
            phi[i][i] = 1.0;
            inverse[i][i] = 1.0;
        }
        phi[0][1] = 0.2;
        phi[0][3] = 0.1;
        phi[3][4] = -0.15;
        inverse[0][1] = -0.2;
        inverse[0][3] = -0.1;
        inverse[0][4] = -0.015;
        inverse[3][4] = 0.15;
        let lambda = [0.0, 4.0, 9.0, 16.0, 25.0, 0.0];
        let cross_mass = [0.1, -0.2, 0.12, 0.04, -0.08, 0.03];
        let support_trace = [0.0, -1.0, 0.0, 0.0, 0.0, -0.5];
        let (mut mass, mut stiffness, mut mobility) = ([[0.0; 7]; 7], [[0.0; 7]; 7], [[0.0; 6]; 6]);
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    mass[i][j] += inverse[k][i] * inverse[k][j];
                    stiffness[i][j] += inverse[k][i] * lambda[k] * inverse[k][j];
                    mobility[i][j] += phi[i][k] * phi[j][k];
                }
                stiffness[i][6] -= stiffness[i][j];
                stiffness[6][i] -= stiffness[i][j];
                stiffness[6][6] += stiffness[i][j];
            }
            mass[i][6] = cross_mass[i];
            mass[6][i] = cross_mass[i];
        }
        mass[6][6] = 2.0;
        assert!(mobility[0][1] != 0.0 && mobility[0][3] != 0.0 && mobility[3][4] != 0.0);
        let (mut m, mut k) = (Coo::new(7, 7), Coo::new(7, 7));
        for i in 0..7 {
            for j in 0..7 {
                m.push(i, j, mass[i][j]);
                k.push(i, j, stiffness[i][j]);
            }
        }
        let pencil = TetElasticAssembly {
            mass: m.assemble(),
            stiffness: k.assemble(),
            free_dofs: (0..7).collect(),
            total_mass_kg: 9.0,
            minimum_scaled_jacobian: 1.0,
            element_volumes_m3: vec![1.0],
        };
        let surface = [0.1, -0.2, 0.3, -0.15, 0.05, -0.1];
        let outgoing = [0.0, 0.4, 0.3, 0.0, 0.0, 0.0];
        let p0 = support_impulse(1.2, [0.4, 0.3], [0.3, 0.6], [0.6, 0.8]);
        let p = [p0[0], p0[1], p0[2], 0.8, -0.07, 0.04];
        let delta: [f64; 6] = core::array::from_fn(|i| (0..6).map(|j| mobility[i][j] * p[j]).sum());
        let vc = 0.23;
        let after: [f64; 7] = core::array::from_fn(|i| {
            if i == 6 {
                vc
            } else {
                surface[i] + outgoing[i] - support_trace[i] * vc
            }
        });
        let before: [f64; 7] =
            core::array::from_fn(|i| if i == 6 { vc } else { after[i] - delta[i] });
        let u = [0.03, -0.04, 0.05, 0.06, -0.02, 0.01, 0.02];
        let basis: Vec<_> = phi.iter().flatten().copied().collect();
        let prepare = || {
            SupportedSpectralSolid::prepare_with_basis(
                &pencil,
                &[6],
                &lambda,
                &basis,
                0.01,
                &u,
                &before,
                &[0.0],
                controls(),
                cx,
            )
            .unwrap()
        };
        let (mut body, mut permuted) = (prepare(), prepare());
        let rows: [[f64; 7]; 6] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if j == 6 {
                    support_trace[i]
                } else if i == j {
                    1.0
                } else {
                    0.0
                }
            })
        });
        let contacts = core::array::from_fn::<_, 2, _>(|i| SolidEllipticContact {
            id: [40, 7][i],
            velocity_weights: [&rows[3 * i][..], &rows[3 * i + 1], &rows[3 * i + 2]],
            gap_m: 0.0,
            activation_tolerance_m: 1e-12,
            prescribed_surface_velocity_mps: [
                surface[3 * i],
                surface[3 * i + 1],
                surface[3 * i + 2],
            ],
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
            friction_coefficient: EllipticFriction {
                coefficients_mu: [[0.3, 0.6], [0.45, 0.2]][i],
                principal_axis: [[0.6, 0.8], [-0.8, 0.6]][i],
            },
        });
        let report = body
            .apply_elliptic_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        permuted
            .apply_elliptic_impacts(&[contacts[1], contacts[0]], impact_controls(), cx)
            .unwrap();
        assert_eq!(bits(&body), bits(&permuted));
        for i in 0..6 {
            near(body.velocity_mps()[i], after[i], 4e-9);
            near(report.contact.law.impulses[i / 3][i % 3], p[i], 4e-9);
        }
        let reaction: f64 = (0..6)
            .map(|i| cross_mass[i] * delta[i] - support_trace[i] * p[i])
            .sum();
        let relative_work: f64 = (0..6).map(|i| p[i] * (outgoing[i] - 0.5 * delta[i])).sum();
        let obstacle_work: f64 = (0..6).map(|i| p[i] * surface[i]).sum();
        near(report.impulse.support_reaction_kg_m_s[0], reaction, 4e-9);
        near(report.impulse.energy.support_work_j, reaction * vc, 4e-9);
        near(report.contact.relative_work_j, relative_work, 4e-9);
        near(
            report.impulse.kinetic_change_j,
            relative_work + obstacle_work + reaction * vc,
            5e-9,
        );
        body.advance_affine_supports(&[0.0], &[0.0; 7], &[0.0; 7], cx)
            .unwrap();
        let mut q = [0.0; 6];
        let mut dq = q;
        for i in 0..6 {
            let q0: f64 = (0..6).map(|j| inverse[i][j] * (u[j] - u[6])).sum();
            let v0: f64 = (0..6).map(|j| inverse[i][j] * (after[j] - vc)).sum();
            let w = lambda[i].sqrt();
            if w == 0.0 {
                q[i] = q0 + 0.01 * v0;
                dq[i] = v0;
            } else {
                let (s, c) = (0.01 * w).sin_cos();
                q[i] = q0 * c + v0 * s / w;
                dq[i] = -w * q0 * s + v0 * c;
            }
        }
        for i in 0..6 {
            near(
                body.displacement_m()[i],
                u[6] + 0.01 * vc + (0..6).map(|j| phi[i][j] * q[j]).sum::<f64>(),
                4e-10,
            );
            near(
                body.velocity_mps()[i],
                vc + (0..6).map(|j| phi[i][j] * dq[j]).sum::<f64>(),
                4e-9,
            );
        }
    });
}

#[derive(Debug, PartialEq)]
struct History {
    count: usize,
}

#[test]
fn source_callbacks_use_full_support_kinematics_and_publish_only_after_all_gates() {
    with_cx(false, |cx| {
        let initial = [-1.0, 2.4, 0.3, 0.4];
        let base = geometry(contact(0.0, [0.0, 0.1, 0.3]));
        let contacts = [
            SolidImpactGeometry { id: 40, ..base },
            SolidImpactGeometry { id: 7, ..base },
        ];
        let mut body = make(initial, 1.0, controls(), cx);
        let original = bits(&body);
        let mut histories = [History { count: 3 }, History { count: 5 }];
        let mut calls = Vec::new();
        let result = body.apply_elliptic_impacts_with(
            &contacts,
            &mut histories,
            impact_controls(),
            |old, q| {
                calls.push((q.contact_id, q.original_index));
                near(q.relative_velocity_mps[1], 1.9, 1e-14);
                near(q.relative_velocity_mps[2], 0.0, 1e-14);
                if q.contact_id == 40 {
                    return Err("second active source refuses");
                }
                Ok((
                    EllipticFriction {
                        coefficients_mu: [0.2, 0.4],
                        principal_axis: [1.0, 0.0],
                    },
                    History {
                        count: old.count + 1,
                    },
                ))
            },
            cx,
        );
        assert!(matches!(
            result,
            Err(SupportedSolidEllipticSourceError::Contact(
                SolidEllipticSourceError::Source {
                    contact_id: 40,
                    original_index: 0,
                    source: "second active source refuses"
                }
            ))
        ));
        assert_eq!(calls, [(7, 1), (40, 0)]);
        assert_eq!(histories, [History { count: 3 }, History { count: 5 }]);
        assert_eq!(bits(&body), original);
        let contacts = [
            contacts[0],
            SolidImpactGeometry {
                gap_m: 1.0,
                ..contacts[1]
            },
        ];
        let mut actual_calls = 0;
        let r = body
            .apply_elliptic_impacts_with(
                &contacts,
                &mut histories,
                impact_controls(),
                |old, q| {
                    actual_calls += 1;
                    assert_eq!(q.contact_id, 40);
                    assert_eq!(q.original_index, 0);
                    let mu =
                        0.1 + 0.1 * q.relative_velocity_mps[1].hypot(q.relative_velocity_mps[2]);
                    Ok::<_, &str>((
                        EllipticFriction {
                            coefficients_mu: [mu, 0.4],
                            principal_axis: [1.0, 0.0],
                        },
                        History {
                            count: old.count + 1,
                        },
                    ))
                },
                cx,
            )
            .unwrap();
        assert_eq!(actual_calls, 1);
        assert_eq!(histories, [History { count: 4 }, History { count: 5 }]);
        near(r.contact.law.impulses[0][1], -3.0 * 0.29, 3e-10);
        assert_eq!(r.contact.law.impulses[1], [0.0; 3]);
        assert_eq!(body.velocity_mps()[3], 0.4);
    });
}

#[test]
fn cancellation_and_late_native_refusal_preserve_directional_history_and_retry() {
    with_cx(false, |cx| {
        let initial = [-1.0, 2.4, 0.0, 0.4];
        let contacts = [geometry(contact(0.0, [0.0; 3]))];
        let mut body = make(initial, 1.0, controls(), cx);
        let before = bits(&body);
        let mut histories = [History { count: 5 }];
        let law = EllipticFriction {
            coefficients_mu: [0.25, 0.6],
            principal_axis: [1.0, 0.0],
        };
        let mut calls = 0;
        with_cx(true, |cancelled| {
            assert!(
                body.apply_elliptic_impacts_with(
                    &contacts,
                    &mut histories,
                    impact_controls(),
                    |old, _| {
                        calls += 1;
                        Ok::<_, &str>((
                            law,
                            History {
                                count: old.count + 1,
                            },
                        ))
                    },
                    cancelled
                )
                .is_err()
            );
        });
        assert_eq!(calls, 0);
        assert_eq!(bits(&body), before);
        assert_eq!(histories, [History { count: 5 }]);
        let mut limits = controls();
        limits.spectral.maximum_relative_operator_error = 0.02;
        limits.spectral.maximum_relative_state_error = 0.02;
        limits.spectral.maximum_work_error_j = 1.0;
        limits.maximum_work_error_j = 1e-16;
        let mut inaccurate = make(initial, 1.001, limits, cx);
        let inaccurate_before = bits(&inaccurate);
        let mut loose = impact_controls();
        loose.maximum_impulse_error_kg_m_s = 1.0;
        let result = inaccurate.apply_elliptic_impacts_with(
            &contacts,
            &mut histories,
            loose,
            |old, _| {
                calls += 1;
                Ok::<_, &str>((
                    law,
                    History {
                        count: old.count + 1,
                    },
                ))
            },
            cx,
        );
        assert!(matches!(
            result,
            Err(SupportedSolidEllipticSourceError::Solid(
                SupportedSolidError::Energy { .. }
            ))
        ));
        assert_eq!(calls, 1);
        assert_eq!(bits(&inaccurate), inaccurate_before);
        assert_eq!(histories, [History { count: 5 }]);
        let mut reference = make(initial, 1.0, controls(), cx);
        let mut expected = [History { count: 5 }];
        for (owner, history) in [(&mut body, &mut histories), (&mut reference, &mut expected)] {
            owner
                .apply_elliptic_impacts_with(
                    &contacts,
                    history,
                    impact_controls(),
                    |old, _| {
                        Ok::<_, &str>((
                            law,
                            History {
                                count: old.count + 1,
                            },
                        ))
                    },
                    cx,
                )
                .unwrap();
            owner
                .advance_affine_supports(&[0.0], &[0.0; 4], &[0.0; 4], cx)
                .unwrap();
        }
        assert_eq!(bits(&body), bits(&reference));
        assert_eq!(histories, expected);
    });
}
