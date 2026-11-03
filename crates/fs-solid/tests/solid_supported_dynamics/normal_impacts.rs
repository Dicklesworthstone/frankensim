//! Independent full-mass support reactions and contact work at one event.

use super::*;
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_solid::dynamics::impact::{SolidNormalContact, SolidNormalImpactControls};
use fs_solid::dynamics::supported::{SupportedSolidError, SupportedSolidNormalImpactError};

fn impact_controls() -> SolidNormalImpactControls {
    let mut contact = MultiImpactConfig::default();
    contact.solver.gradient_tolerance = 1e-12;
    contact.energy_tolerance = 1e-12;
    SolidNormalImpactControls {
        contact,
        maximum_impulse_error_kg_m_s: 1e-11,
        maximum_contact_speed_error_mps: 1e-11,
        maximum_work_error_j: 1e-11,
        maximum_workspace_bytes: 1 << 20,
        maximum_response_multiply_terms: 1_000_000,
    }
}

fn normal(row: &[f64], surface_speed: f64, restitution: f64) -> SolidNormalContact<'_> {
    SolidNormalContact {
        normal_velocity_weights: row,
        gap_m: 0.0,
        activation_tolerance_m: 1e-12,
        prescribed_surface_speed_mps: surface_speed,
        restitution,
    }
}

#[test]
fn g1_supported_contact_keeps_cross_mass_and_moving_support_work_in_continuation() {
    with_cx(false, |cx| {
        let dt = 0.01;
        let u = [0.05, 0.02];
        let mut solid = prepare(dt, u, [-0.6, 0.4], 0.0, controls(), cx);
        let before_energy = energy(u, [-0.6, 0.4]);
        let report = solid
            .apply_normal_impacts(&[normal(&[1.0, -1.0], 0.0, 0.5)], impact_controls(), cx)
            .unwrap();

        // Relative speed -1 becomes +1/2. M_ff=2, so Delta v_f=1.5 and
        // p=3. Contact applies [3,-3]; the support supplies .35*1.5+3.
        near(solid.velocity_mps()[0], 0.9, 4e-13);
        assert_eq!(solid.velocity_mps()[1].to_bits(), 0.4_f64.to_bits());
        assert_eq!(solid.displacement_m(), &u);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(solid.accepted_steps(), 0);
        near(report.impulse.momentum_change_kg_m_s[0], 3.0, 1e-12);
        near(report.impulse.momentum_change_kg_m_s[1], 0.525, 1e-12);
        near(report.impulse.support_reaction_kg_m_s[0], 3.525, 1e-12);
        near(report.impulse.energy.applied_work_j, -0.75, 1e-12);
        near(report.impulse.energy.support_work_j, 1.41, 1e-12);
        near(report.impulse.kinetic_change_j, 0.66, 1e-12);
        near(report.contact.relative_work_j, -0.75, 1e-12);
        near(report.contact.dissipated_energy_j, 0.75, 1e-12);
        assert_eq!(report.contact.prescribed_surface_work_j, 0.0);
        near(
            energy(u, [solid.velocity_mps()[0], solid.velocity_mps()[1]]) - before_energy,
            report.impulse.energy.applied_work_j + report.impulse.energy.support_work_j,
            1e-12,
        );
        assert!(report.impulse.maximum_impulse_residual_kg_m_s < 1e-12);
        assert!(report.impulse.energy.support_work_j > 0.0);

        // The original elastic trajectory starts from the accepted impact,
        // retaining support velocity and the two separate cumulative works.
        let reference = Reference {
            support_u: u[1],
            support_v: 0.4,
            support_a: 0.0,
            jerk: 0.0,
            relative_u: u[0] - u[1],
            relative_v: 0.5,
            f0: 0.0,
            f_rate: 0.0,
            support_f0: 0.0,
            support_f_rate: 0.0,
        };
        solid
            .advance_affine_supports(&[0.0], &[0.0; 2], &[0.0; 2], cx)
            .unwrap();
        let expected = reference.state(dt);
        for i in 0..2 {
            near(solid.displacement_m()[i], expected.0[i], 2e-12);
            near(solid.velocity_mps()[i], expected.1[i], 2e-12);
        }
        near(solid.cumulative_applied_work_j(), -0.75, 1e-12);
        assert_eq!(solid.accepted_steps(), 1);
    });
}

#[test]
fn g3_prescribed_contact_speed_and_support_work_transform_with_the_same_boost() {
    with_cx(false, |cx| {
        let mut responses = Vec::new();
        for boost in [0.0, 3.0] {
            let before = [-0.75 + boost, -0.1 + boost];
            let mut solid = prepare(0.01, [0.0; 2], before, 0.0, controls(), cx);
            let surface_speed = 0.25 + boost;
            let report = solid
                .apply_normal_impacts(
                    &[normal(&[1.0, 0.0], surface_speed, 0.4)],
                    impact_controls(),
                    cx,
                )
                .unwrap();
            let impulse = report.impulse.momentum_change_kg_m_s[0];
            near(impulse, 2.8, 2e-12);
            near(report.impulse.support_reaction_kg_m_s[0], 0.49, 2e-12);
            near(solid.velocity_mps()[0] - surface_speed, 0.4, 1e-12);
            assert_eq!(solid.velocity_mps()[1].to_bits(), before[1].to_bits());
            let obstacle_work = impulse * surface_speed;
            let relative_work = report.impulse.energy.applied_work_j - obstacle_work;
            // Effective mass 2 and restitution .4 dissipate .5*2*(1-.4^2).
            near(relative_work, -0.84, 4e-12);
            near(report.contact.relative_work_j, -0.84, 4e-12);
            near(
                report.contact.prescribed_surface_work_j,
                obstacle_work,
                4e-12,
            );
            near(
                report.impulse.energy.applied_work_j + report.impulse.energy.support_work_j,
                report.impulse.kinetic_change_j,
                4e-12,
            );
            responses.push(report.impulse);
        }
        near(
            responses[1].energy.applied_work_j - responses[0].energy.applied_work_j,
            3.0 * 2.8,
            4e-12,
        );
        near(
            responses[1].energy.support_work_j - responses[0].energy.support_work_j,
            3.0 * 0.49,
            4e-12,
        );
    });
}

#[test]
fn g1_g4_redundant_contacts_share_one_native_response_and_locked_approach_refuses() {
    with_cx(false, |cx| {
        let mut solid = prepare(0.01, [0.05, 0.02], [-0.6, 0.4], 0.0, controls(), cx);
        let report = solid
            .apply_normal_impacts(
                &[
                    normal(&[1.0, -1.0], 0.0, 0.5),
                    normal(&[1.0, -1.0], 0.0, 0.5),
                ],
                impact_controls(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            report.contact.law.determinacy,
            MultiImpactDeterminacy::UniqueVelocity(_)
        ));
        near(solid.velocity_mps()[0], 0.9, 2e-12);
        near(report.impulse.support_reaction_kg_m_s[0], 3.525, 3e-12);
        let accepted = bits(&solid);

        // A contact attached only to the prescribed coordinate has zero
        // mobility. Its approaching obstacle cannot be resolved by free motion.
        assert!(
            solid
                .apply_normal_impacts(&[normal(&[0.0, 1.0], 1.0, 0.0)], impact_controls(), cx)
                .is_err()
        );
        assert_eq!(bits(&solid), accepted);
    });
}

#[test]
fn g4_supported_contact_refusals_preserve_both_owners_and_subsequent_motion() {
    with_cx(false, |cx| {
        let make = || prepare(0.01, [0.05, 0.02], [-0.6, 0.4], 0.0, controls(), cx);
        let (mut retry, mut reference) = (make(), make());
        let original = bits(&retry);
        let contacts = [normal(&[1.0, -1.0], 0.0, 0.5)];
        with_cx(true, |cancelled| {
            assert!(
                retry
                    .apply_normal_impacts(&contacts, impact_controls(), cancelled)
                    .is_err()
            );
        });
        for resource in 0..2 {
            let mut limited = impact_controls();
            if resource == 0 {
                limited.maximum_workspace_bytes = 1;
            } else {
                limited.maximum_response_multiply_terms = 1;
            }
            assert!(retry.apply_normal_impacts(&contacts, limited, cx).is_err());
            assert_eq!(bits(&retry), original);
        }
        retry
            .apply_normal_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        reference
            .apply_normal_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        retry
            .advance_affine_supports(&[0.0], &[0.1, 0.0], &[0.2, 0.0], cx)
            .unwrap();
        reference
            .advance_affine_supports(&[0.0], &[0.1, 0.0], &[0.2, 0.0], cx)
            .unwrap();
        assert_eq!(bits(&retry), bits(&reference));

        // Let the free basis admit a small mass-normalization error. The
        // contact solve and free impulse can succeed at their loose limits,
        // but the original full-system work gate must still reject publication.
        let mut strict = controls();
        strict.spectral.maximum_relative_operator_error = 0.02;
        strict.spectral.maximum_relative_state_error = 0.02;
        strict.spectral.maximum_work_error_j = 1.0;
        strict.maximum_work_error_j = 1e-16;
        let mut inaccurate = SupportedSpectralSolid::prepare_with_basis(
            &pencil(),
            &[8],
            &[9.0],
            &[1.001 / 2.0_f64.sqrt()],
            0.01,
            &[0.05, 0.02],
            &[-0.6, 0.4],
            &[0.0],
            strict,
            cx,
        )
        .unwrap();
        let before = bits(&inaccurate);
        let mut loose_impulse = impact_controls();
        loose_impulse.maximum_impulse_error_kg_m_s = 1.0;
        assert!(matches!(
            inaccurate.apply_normal_impacts(&contacts, loose_impulse, cx),
            Err(SupportedSolidNormalImpactError::Solid(
                SupportedSolidError::Energy { .. }
            ))
        ));
        assert_eq!(bits(&inaccurate), before);
    });
}
