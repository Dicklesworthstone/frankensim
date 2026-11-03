//! Directional multi-contact momentum through the retained native histories.
use super::*;
use fs_couple::inelastic_impact::multiple::elliptic::*;

fn audit(report: &InelasticMixedEllipticContactReport) {
    near(norm(report.linear_momentum_residual_kg_m_s), 0.0, 3e-10);
    near(norm(report.angular_momentum_residual_kg_m2_s), 0.0, 3e-10);
    near(report.work_energy_residual_j, 0.0, 3e-10);
    near(report.dissipation_residual_j, 0.0, 3e-10);
    near(report.endpoint_velocity_residual_m_s, 0.0, 3e-10);
    near(report.absolute_linear_residual_kg_m_s, 0.0, 3e-10);
    near(report.absolute_angular_residual_kg_m2_s, 0.0, 3e-10);
    near(report.absolute_work_residual_j, 0.0, 3e-10);
}
fn directional(
    coefficients_mu: [f64; 2],
    principal_axis: [f64; 2],
) -> [InelasticMixedEllipticContact; 2] {
    mixed_contacts(0.25, 0.0).map(|c| {
        c.with_friction(EllipticFriction {
            coefficients_mu,
            principal_axis,
        })
    })
}
fn near_bodies(a: &[RigidBodyState], b: &[RigidBodyState]) {
    for (a, b) in a.iter().zip(b) {
        assert_eq!(a.pose(), b.pose());
        near(
            norm(a.linear_momentum_world().sub(b.linear_momentum_world())),
            0.0,
            2e-9,
        );
        near(
            norm(a.angular_momentum_body().sub(b.angular_momentum_body())),
            0.0,
            2e-9,
        );
    }
}

pub(super) fn check_visco_directional(motion: &FiniteViscoTetMotion<'_>, cx: &Cx<'_>) {
    let isotropic = mixed_contacts(0.25, 0.0);
    let bodies: Vec<_> = isotropic
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let FiniteSolidContactCounterpart::Solid(point) = c.counterpart else {
                unreachable!()
            };
            let rigid = rigid_at(motion.point_kinematics(point, cx).unwrap().position_m, -0.8);
            RigidBodyState::new(
                rigid.pose(),
                Vec3::new(-1.6, if i == 0 { 1.6 } else { -0.8 }, 0.6),
                Vec3::ZERO,
            )
            .unwrap()
        })
        .collect();
    let properties = [properties(); 2];
    let circle = visco_mixed_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &isotropic,
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    // Equal coefficients preserve the original scalar coordinates even when
    // an arbitrary principal direction is supplied. This exercises the shared
    // native adapter, not a second diagonal-mobility implementation.
    let equal = visco_mixed_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &directional([0.25; 2], [3.0, 4.0]),
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    assert_eq!(equal.rigid_bodies, circle.rigid_bodies);
    assert_eq!(equal.solid.velocity_mps(), circle.solid.velocity_mps());
    assert_eq!(equal.report.contact.transformed, circle.report.contact);
    audit(&equal.report);
    let frictionless = visco_mixed_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &mixed_contacts(0.0, 0.0),
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    let zero = visco_mixed_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &directional([0.0; 2], [0.0, 1.0]),
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    assert_eq!(zero.rigid_bodies, frictionless.rigid_bodies);
    assert_eq!(zero.solid.velocity_mps(), frictionless.solid.velocity_mps());

    // The collapsed ellipse imposes exactly zero z traction. The original full
    // physical mobility still participates, including cross-point responses.
    let free_z = visco_mixed_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &directional([0.4, 0.0], [1.0, 0.0]),
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    let free_y = visco_mixed_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &directional([0.4, 0.0], [0.0, 1.0]),
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    for i in 0..2 {
        near(free_z.report.contact.impulses[i][2], 0.0, 1e-13);
        near(free_y.report.contact.impulses[i][1], 0.0, 1e-13);
        assert!(free_z.report.contact.impulses[i][1].abs() > 1e-4);
        assert!(free_y.report.contact.impulses[i][2].abs() > 1e-4);
    }
    for event in [&free_z, &free_y] {
        audit(&event.report);
        assert_eq!(
            event.solid.solid().element_states(),
            motion.solid().element_states()
        );
        assert_eq!(
            event.solid.solid().element_first_piola_pa(),
            motion.solid().element_first_piola_pa()
        );
        assert_eq!(
            event.solid.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        assert_eq!(
            event.solid.solid().stored_energy_j(),
            motion.solid().stored_energy_j()
        );
        assert_eq!(event.solid.solid().time_s(), motion.solid().time_s());
        assert_eq!(
            event.solid.reference_nodal_density_kg_m3(),
            motion.reference_nodal_density_kg_m3()
        );
    }
    assert!(motion.solid().element_states()[0].dissipated > 0.0);

    // Rotate the contact coordinates and texture together, retaining the same
    // PHYSICAL ellipse and complete native body response.
    let mut original = directional([0.2, 0.5], [0.6, 0.8]);
    original[1].friction_coefficient = EllipticFriction {
        coefficients_mu: [0.5, 0.2],
        principal_axis: [1.0, 0.0],
    };
    let rotated = original.map(|c| InelasticMixedEllipticContact {
        tangent_world: Vec3::new(0.0, 0.6, 0.8),
        friction_coefficient: EllipticFriction {
            principal_axis: [
                0.6 * c.friction_coefficient.principal_axis[0]
                    + 0.8 * c.friction_coefficient.principal_axis[1],
                -0.8 * c.friction_coefficient.principal_axis[0]
                    + 0.6 * c.friction_coefficient.principal_axis[1],
            ],
            ..c.friction_coefficient
        },
        ..c
    });
    let a = visco_mixed_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &original,
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    let b = visco_mixed_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &rotated,
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    near_bodies(&a.rigid_bodies, &b.rigid_bodies);
    for (a, b) in a.solid.velocity_mps().iter().zip(b.solid.velocity_mps()) {
        for axis in 0..3 {
            near(a[axis], b[axis], 2e-9);
        }
    }
    near(
        a.report.contact.dissipated_energy,
        b.report.contact.dissipated_energy,
        2e-9,
    );
    audit(&a.report);
    audit(&b.report);
    for (i, contact) in original.iter().enumerate() {
        let [c, s] = contact.friction_coefficient.principal_axis;
        let mu = contact.friction_coefficient.coefficients_mu;
        let v = a.report.contact.relative_velocities_after[i];
        let j = a.report.contact.impulses[i];
        let slip = [c * v[1] + s * v[2], -s * v[1] + c * v[2]];
        let traction = [c * j[1] + s * j[2], -s * j[1] + c * j[2]];
        let dual_norm = (mu[0] * slip[0]).hypot(mu[1] * slip[1]);
        assert!(dual_norm > 1e-4);
        for axis in 0..2 {
            near(
                traction[axis],
                -j[0] * mu[axis] * mu[axis] * slip[axis] / dual_norm,
                3e-10,
            );
        }
    }
    let mut reversed = original;
    reversed.reverse();
    let reordered = visco_mixed_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &reversed,
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    assert_eq!(a.rigid_bodies, reordered.rigid_bodies);
    assert_eq!(a.solid.velocity_mps(), reordered.solid.velocity_mps());

    // Refusal leaves the input histories and physical momenta available for an
    // exact valid retry, even after preparing the complete mobility.
    let mut invalid = original;
    invalid[1].friction_coefficient.coefficients_mu[0] = -0.1;
    assert!(
        visco_mixed_elliptic_contact_candidate(
            motion,
            &[],
            &bodies,
            &properties,
            &invalid,
            ORIGIN,
            mixed_controls(),
            cx,
        )
        .is_err()
    );
    let retry = visco_mixed_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &properties,
        &original,
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    assert_eq!(a.rigid_bodies, retry.rigid_bodies);
    assert_eq!(a.solid.velocity_mps(), retry.solid.velocity_mps());
    assert_eq!(a.report.contact, retry.report.contact);
}

pub(super) fn check_plastic_directional_supports(motion: &FinitePlasticTetMotion<'_>, cx: &Cx<'_>) {
    let contacts = directional([0.3, 0.15], [0.8, 0.6]);
    let bodies: Vec<_> = contacts
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let FiniteSolidContactCounterpart::Solid(point) = c.counterpart else {
                unreachable!()
            };
            let rigid = rigid_at(motion.point_kinematics(point, cx).unwrap().position_m, -0.8);
            RigidBodyState::new(
                rigid.pose(),
                Vec3::new(-1.6, if i == 0 { 0.6 } else { -0.8 }, 0.4),
                Vec3::ZERO,
            )
            .unwrap()
        })
        .collect();
    let held: Vec<_> = (0..12).filter(|&dof| dof != 3 && dof != 6).collect();
    let event = plastic_mixed_elliptic_contact_candidate(
        motion,
        &held,
        &bodies,
        &[properties(); 2],
        &contacts,
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    // A unique transformed impulse determines the physical ellipse image and
    // therefore moving-support work. The old conservative default would refuse
    // this valid driven-support event despite its uniqueness witness.
    assert!(matches!(
        event.report.contact.transformed.determinacy,
        FrictionalContactDeterminacy::UniqueImpulses
    ));
    assert!(event.report.support_work_j.abs() > 0.0);
    assert!(
        !event
            .report
            .solid_impulse
            .support_impulses_kg_m_s
            .is_empty()
    );
    for &dof in &held {
        assert_eq!(
            event.solid.velocity_mps()[dof / 3][dof % 3],
            motion.velocity_mps()[dof / 3][dof % 3]
        );
    }
    audit(&event.report);
    assert_eq!(
        event.solid.solid().element_states(),
        motion.solid().element_states()
    );
    assert_eq!(
        event.solid.solid().element_deformation(),
        motion.solid().element_deformation()
    );
    assert_eq!(
        event.solid.solid().displacement_m(),
        motion.solid().displacement_m()
    );
    assert_eq!(
        event.solid.solid().stored_energy_j(),
        motion.solid().stored_energy_j()
    );
    assert!(event.solid.solid().element_states()[0].alpha > 0.18);
    let zero = contacts[0].with_friction(EllipticFriction {
        coefficients_mu: [0.0; 2],
        principal_axis: [1.0, 0.0],
    });
    let duplicate = [zero, InelasticMixedEllipticContact { id: 99, ..zero }];
    let error = plastic_mixed_elliptic_contact_candidate(
        motion,
        &held,
        &bodies,
        &[properties(); 2],
        &duplicate,
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap_err();
    let InelasticMixedEllipticContactError::UnresolvedSupportWork(receipt) = error else {
        panic!("moving support needs physical impulse authority, got {error:?}");
    };
    assert!(matches!(
        receipt.transformed.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
}

pub(super) fn check_cancelled_visco(
    motion: &FiniteViscoTetMotion<'_>,
    rigid: RigidBodyState,
    cx: &Cx<'_>,
) {
    assert!(matches!(
        visco_mixed_elliptic_contact_candidate(
            motion,
            &[],
            &[rigid],
            &[properties()],
            &directional([0.2, 0.4], [1.0, 0.0])[..1],
            ORIGIN,
            mixed_controls(),
            cx,
        ),
        Err(InelasticMixedEllipticContactError::Cancelled)
    ));
}
