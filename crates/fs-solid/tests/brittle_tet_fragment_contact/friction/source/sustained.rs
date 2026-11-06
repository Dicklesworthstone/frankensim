//! A source-held static ellipse carries physical shear across failed faces.

use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::friction::sustained::source::BrittleFaceSourcedSustainedFrictionError;

#[test]
fn original_static_sources_carry_shear_without_advancing_slip_or_publishing_failed_trials() {
    with_oblique_motion(&NODES, [0.0, 0.0, 1.0], [[0.0; 3]; 2], |motion, cx| {
        let law = law(motion.solid().materials()[0].material().material());
        let history = history(&law);
        let accepted_history = history.clone();
        let conditions = conditions();
        let vertex = BrittleVertexFrictionSource {
            law: &law,
            state: &history,
            conditions: &conditions,
            coefficient_role: SourceImpulseCoefficientRole::Static,
        };
        let binding = BrittleFaceFrictionSource {
            face_id: 951,
            tangent_direction: [1.0, 0.0, 0.0],
            vertices: [vertex; 3],
        };
        let acceleration = [3.0 / 8.0, 0.0, 1.0];
        let mut forces = Vec::new();
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                for axis in 0..3 {
                    forces.push((3 * node + axis, sign * acceleration[axis] / 4.0));
                }
            }
        }
        let h = 1.0 / 64.0;
        let controls = super::super::sustained::sustained_controls();
        let accepted = motion
            .advance_released_face_sustained_sourced_frictional(
                h,
                &[binding],
                &[],
                &forces,
                settings(),
                controls,
                cx,
            )
            .unwrap();
        assert_eq!(accepted.samples.len(), 3);
        for (i, sample) in accepted.samples.iter().enumerate() {
            assert_eq!(sample.row, accepted.step.rows[i]);
            assert_eq!(
                sample.sample.coefficient_role,
                SourceImpulseCoefficientRole::Static
            );
            assert_eq!(sample.sample.queried_speed_m_s, 0.0);
            close(sample.sample.coefficients_mu[0], 0.5, 1e-14);
            close(sample.sample.coefficients_mu[1], 0.6, 1e-14);
            // M=(I+ones)/20 for each unit-mass tet gives force a/8.
            // Shear ratio 3/8 exceeds kinetic mu=1/5 and fits static mu=1/2.
            close(accepted.step.forces_n[i][0], acceleration[2] / 8.0, 2e-8);
            close(accepted.step.forces_n[i][1], acceleration[0] / 8.0, 2e-8);
            close(accepted.step.forces_n[i][2], 0.0, 2e-8);
            let next = sample.sample.state.constitutive_state();
            let old = history.constitutive_state();
            assert_eq!(next.slip_m(), old.slip_m());
            assert_eq!(next.principal_axis(), old.principal_axis());
            assert_eq!(next.contact_temperature_k(), old.contact_temperature_k());
            for source in sample.sample.state.friction_sources() {
                assert_eq!(source.state_coordinate("speed"), Some(0.0));
                close(source.state_coordinate("path").unwrap(), 0.2525, 1e-14);
            }
        }
        let expected_work = 5.0 / 8.0 * acceleration.iter().map(|a| a * a).sum::<f64>() * h * h;
        close(
            accepted.step.energy.dead_load_mean_work_j,
            expected_work,
            1e-10,
        );
        close(
            accepted.step.energy.inertia.kinetic_change_j,
            expected_work,
            1e-10,
        );
        close(accepted.step.energy.total_balance_error_j, 0.0, 1e-10);
        close(accepted.step.normal_contact_work_j, 0.0, 1e-10);
        close(accepted.step.tangential_contact_work_j, 0.0, 1e-10);
        assert_eq!(history, accepted_history);
        assert_eq!(
            accepted.step.motion.accepted_steps(),
            motion.accepted_steps() + 1
        );

        let kinetic = BrittleFaceFrictionSource {
            vertices: [BrittleVertexFrictionSource {
                coefficient_role: SourceImpulseCoefficientRole::Kinetic,
                ..vertex
            }; 3],
            ..binding
        };
        assert!(matches!(
            motion.advance_released_face_sustained_sourced_frictional(
                h,
                &[kinetic],
                &[],
                &forces,
                settings(),
                controls,
                cx,
            ),
            Err(BrittleFaceSourcedSustainedFrictionError::Source(
                BrittleFaceSourcedFrictionError::Invalid(_)
            ))
        ));
        // Queries succeed, but this load requires unsupported finite sliding.
        let excessive: Vec<_> = forces
            .iter()
            .map(|&(dof, f)| (dof, if dof % 3 == 0 { f * 8.0 } else { f }))
            .collect();
        assert!(
            motion
                .advance_released_face_sustained_sourced_frictional(
                    h,
                    &[binding],
                    &[],
                    &excessive,
                    settings(),
                    controls,
                    cx,
                )
                .is_err()
        );
        assert_eq!(history, accepted_history);
        let retry = motion
            .advance_released_face_sustained_sourced_frictional(
                h,
                &[binding],
                &[],
                &forces,
                settings(),
                controls,
                cx,
            )
            .unwrap();
        assert_eq!(
            retry.step.motion.velocity_mps(),
            accepted.step.motion.velocity_mps()
        );
        assert_eq!(retry.samples[0].sample, accepted.samples[0].sample);
    });
}
