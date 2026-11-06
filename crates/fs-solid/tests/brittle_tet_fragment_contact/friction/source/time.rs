//! A complete sourced impact/static trajectory retains native source histories.

use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::friction::source::time::{
    BrittleFaceFrictionSourcePhase, BrittleFaceSourcedFrictionTimeError,
};
use fs_solid::brittle::solid::dynamics::fracture::contact::time::BrittleFaceTimeError;

#[test]
fn sourced_duration_captures_with_kinetic_claims_then_carries_static_shear_atomically() {
    let incoming = [SPEED / 8.0, 0.0, SPEED];
    let drift = [DRIFT[0], DRIFT[1], 0.0];
    let velocities = [
        core::array::from_fn(|c| drift[c] + incoming[c]),
        core::array::from_fn(|c| drift[c] - incoming[c]),
    ];
    with_oblique_motion_budget(&NODES, [0.0, 0.0, 1.0], velocities, 8, |motion, cx| {
        let law = law(motion.solid().materials()[0].material().material());
        let history = history(&law);
        let snapshot = history.clone();
        let conditions = conditions();
        let source = BrittleFaceFrictionSource {
            face_id: 951,
            tangent_direction: [1.0, 0.0, 0.0],
            vertices: [BrittleVertexFrictionSource {
                law: &law,
                state: &history,
                conditions: &conditions,
                coefficient_role: SourceImpulseCoefficientRole::Kinetic,
            }; 3],
        };
        let acceleration = [3.0 / 8.0, 0.0, 1.0];
        let mut forces = Vec::new();
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                for c in 0..3 {
                    let value = sign * acceleration[c] / 4.0;
                    forces.push((3 * node + c, [value; 2]));
                }
            }
        }
        let duration = 1.0 / 64.0;
        let controls = super::super::time::time_controls(duration / 2.0);
        let sustained = super::super::sustained::sustained_controls();
        let mut limited = controls;
        limited.maximum_intervals = 1;
        assert!(matches!(
            motion.advance_released_face_time_sourced_frictional(
                duration,
                &[],
                &forces,
                &[0.0],
                &[source],
                ORIGIN,
                settings(),
                limited,
                sustained,
                cx,
            ),
            Err(BrittleFaceSourcedFrictionTimeError::Time(
                BrittleFaceTimeError::WorkLimit {
                    intervals: 1,
                    impacts: 1,
                    ..
                }
            ))
        ));
        assert_eq!(history, snapshot);
        assert_eq!(motion.accepted_steps(), 1);
        let result = motion
            .advance_released_face_time_sourced_frictional(
                duration,
                &[],
                &forces,
                &[0.0],
                &[source],
                ORIGIN,
                settings(),
                controls,
                sustained,
                cx,
            )
            .unwrap();
        assert_eq!(result.step.impacts.len(), 1);
        assert_eq!(result.step.intervals.len(), 2);
        assert_eq!(result.step.sustained.len(), 2);
        assert_eq!(result.sources.len(), 3);
        assert_eq!(
            result.step.motion.accepted_steps(),
            motion.accepted_steps() + 2
        );
        assert_eq!(result.step.motion.time_s(), motion.time_s() + duration);

        for (stage, record) in result.sources.iter().enumerate() {
            assert_eq!(record.samples.len(), 3);
            assert_eq!(
                record.phase,
                if stage == 0 {
                    BrittleFaceFrictionSourcePhase::Impact
                } else {
                    BrittleFaceFrictionSourcePhase::StaticInterval
                }
            );
            if stage == 0 {
                assert_eq!(record.start_time_s, motion.time_s());
                assert_eq!(record.end_time_s, motion.time_s());
            } else {
                assert_eq!(
                    record.start_time_s,
                    result.step.sustained[stage - 1].start_time_s
                );
                assert_eq!(
                    record.end_time_s,
                    result.step.sustained[stage - 1].end_time_s
                );
            }
            for sample in &record.samples {
                let actual = &sample.sample;
                if stage == 0 {
                    assert_eq!(
                        actual.coefficient_role,
                        SourceImpulseCoefficientRole::Kinetic
                    );
                    close(actual.queried_speed_m_s, 2.0 * incoming[0], 1e-12);
                    close(actual.coefficients_mu[0], 0.2 + 0.2 * incoming[0], 1e-12);
                    close(actual.coefficients_mu[1], 0.12525, 1e-12);
                } else {
                    assert_eq!(
                        actual.coefficient_role,
                        SourceImpulseCoefficientRole::Static
                    );
                    close(actual.queried_speed_m_s, 0.0, 1e-9);
                    close(actual.coefficients_mu[0], 0.5, 1e-12);
                    close(actual.coefficients_mu[1], 0.6, 1e-12);
                }
                let old = history.constitutive_state();
                let next = actual.state.constitutive_state();
                assert_eq!(next.slip_m(), old.slip_m());
                assert_eq!(next.principal_axis(), old.principal_axis());
                assert_eq!(next.contact_temperature_k(), old.contact_temperature_k());
                for source in actual.state.friction_sources() {
                    close(source.state_coordinate("path").unwrap(), 0.2525, 1e-14);
                    close(
                        source.state_coordinate("speed").unwrap(),
                        actual.queried_speed_m_s,
                        1e-12,
                    );
                }
            }
        }
        for record in &result.step.sustained {
            for force in &record.forces_n {
                close(force[0], acceleration[2] / 8.0, 2e-8);
                close(force[1], acceleration[0] / 8.0, 2e-8);
                close(force[2], 0.0, 2e-8);
                // Required shear 3/8 exceeds kinetic ~1/5 but fits static1/2.
                assert!(force[1] > 0.3 * force[0]);
            }
        }
        let squared = |v: [f64; 3]| v.into_iter().map(|x| x * x).sum::<f64>();
        let final_internal = core::array::from_fn(|c| incoming[c] + acceleration[c] * duration);
        let loss = 3.0 / 8.0 * squared(incoming);
        let dead = 5.0 / 8.0 * (squared(final_internal) - squared(incoming));
        close(result.step.energy.contact_dissipation_j, loss, 2e-10);
        close(result.step.energy.dead_load_mean_work_j, dead, 2e-10);
        close(
            result.step.energy.kinetic_after_j,
            squared(drift) + 5.0 / 8.0 * squared(final_internal),
            2e-10,
        );
        close(result.step.energy.kinetic_change_j, dead - loss, 2e-10);
        close(result.step.energy.fracture_dissipation_j, 0.0, 2e-12);
        close(result.step.energy.normal_contact_work_j, 0.0, 2e-10);
        close(result.step.energy.tangential_contact_work_j, 0.0, 2e-10);
        close(result.step.energy.total_balance_error_j, 0.0, 2e-10);
        close(
            inventory(&result.step.motion).2,
            result.step.energy.kinetic_after_j,
            2e-12,
        );
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let apex = motion.solid().geometry().nodes_m()[node][2] != 0.0;
                for c in 0..3 {
                    close(
                        result.step.motion.velocity_mps()[node][c],
                        drift[c]
                            + if apex {
                                sign * 2.5 * final_internal[c]
                            } else {
                                0.0
                            },
                        2e-8,
                    );
                }
            }
        }
        assert_eq!(history, snapshot);
        let retry = motion
            .advance_released_face_time_sourced_frictional(
                duration,
                &[],
                &forces,
                &[0.0],
                &[source],
                ORIGIN,
                settings(),
                controls,
                sustained,
                cx,
            )
            .unwrap();
        assert_eq!(
            retry.step.motion.velocity_mps(),
            result.step.motion.velocity_mps()
        );
        for (left, right) in retry.sources.iter().zip(&result.sources) {
            for (a, b) in left.samples.iter().zip(&right.samples) {
                assert_eq!(a.row, b.row);
                assert_eq!(a.sample, b.sample);
            }
        }
    });
}
