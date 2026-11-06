//! A coupled arrival can release an old contact and continue on the new subset.

use super::super::time::time_controls as impact_time_controls;
use super::approach::{
    ARRIVAL_CORNER, SQUARE_NODES, SQUARE_TETS, arrival_gap, corner_motion, force_intervals,
    mixed_controls, square_topology,
};
use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::sustained::approach::BrittleFaceSustainedApproachError;
use fs_solid::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceTimeControls, BrittleFaceTimeError,
};

const RELEASED_CORNER: [f64; 3] = [1.0, 0.0, 0.0];

fn time_controls(duration: f64) -> BrittleFaceTimeControls {
    let (approach, sustained) = mixed_controls();
    let mut controls = impact_time_controls();
    controls.maximum_interval_s = duration;
    controls.maximum_trials = 2 * sustained.solve.maximum_response_evaluations;
    controls.approach = approach;
    controls.impact.activation_tolerance_m = approach.activation_tolerance_m;
    controls.impact.velocity_tolerance_mps = approach.velocity_tolerance_mps;
    controls
}

#[test]
fn coupled_arrival_impacts_release_the_old_corner_and_continue_with_the_new_force_subset() {
    with_cx(false, |cx| {
        let geometry =
            TetLoadGeometry3::new(&SQUARE_NODES, &SQUARE_TETS, budget(), 16, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = square_topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 16, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let gap = 1.0 / 512.0;
        let acceleration = 1.0;
        let duration = 3.0 / 128.0;
        let motion = corner_motion(&split.solid, gap, SPEED, 0.0, false, 4, cx);
        // Common affine transverse acceleration leaves normal contact
        // mechanics unchanged and independently probes global force means.
        let force = rigid_force_intervals(motion.solid(), [0.5, 1.5], acceleration);
        let result = motion
            .advance_released_face_time_with_sustained(
                duration,
                &[],
                &force,
                &[0.0, 0.0],
                ORIGIN,
                settings(),
                time_controls(duration),
                mixed_controls().1,
                cx,
            )
            .unwrap();
        assert_eq!(result.motion.time_s(), motion.time_s() + duration);
        assert_eq!(result.intervals.len(), 2);
        assert_eq!(result.sustained.len(), 2);
        assert_eq!(result.impacts.len(), 1);
        assert_eq!(result.motion.accepted_steps(), 2);
        assert_eq!(result.motion.remaining_steps(), 2);
        let normal_acceleration = 10.0 * acceleration / 7.0;
        let tau = gap / ((SPEED * SPEED + normal_acceleration * gap).sqrt() + SPEED);
        let remainder = duration - tau;
        let c = SPEED + normal_acceleration * tau;
        let b = 15.0 * acceleration * tau / 7.0;
        assert!(remainder > 0.0 && remainder < c / (5.0 * acceleration));
        let hit = &result.impacts[0];
        close(hit.location.time_s - motion.time_s(), tau, 1e-9);
        assert_eq!(hit.location.face_id, 952);
        assert_eq!(result.sustained[0].end_time_s, hit.location.time_s);
        assert_eq!(result.sustained[1].start_time_s, hit.location.time_s);
        assert_eq!(hit.rows.len(), 6);
        assert!(matches!(
            hit.contact.determinacy,
            MultiImpactDeterminacy::UniqueVelocity(_)
        ));
        let mut impulse = vec![0.0; target.nodes_m().len()];
        for (row, &value) in hit.rows.iter().zip(&hit.contact.normal_impulses) {
            impulse[row.nodes[1]] += value;
        }
        // Applying the assembled native mass to delta_v gives physical
        // impulses [c/35,0,c/35,3c/35] on upper corners [0,1,2,3].
        // Enforcing all four corners at rest would require a tensile impulse
        // at corner 1. The unilateral law must instead RELEASE that corner.
        for element in [1, 3] {
            for &node in &target.tetrahedra()[element] {
                let position = target.nodes_m()[node];
                let expected = if position[2] != 0.0 || position == RELEASED_CORNER {
                    0.0
                } else if position == ARRIVAL_CORNER {
                    3.0 * c / 35.0
                } else {
                    c / 35.0
                };
                close(impulse[node], expected, 1e-8);
            }
        }
        for (record, excluded) in result
            .sustained
            .iter()
            .zip([ARRIVAL_CORNER, RELEASED_CORNER])
        {
            assert_eq!(record.rows.len(), 5);
            let mut physical_force = vec![0.0; target.nodes_m().len()];
            for (row, force) in record.rows.iter().zip(&record.solve.forces_n) {
                assert_ne!(target.nodes_m()[row.nodes[0]], excluded);
                physical_force[row.nodes[1]] += force[0];
                assert_eq!(force[1], 0.0);
                assert_eq!(force[2], 0.0);
            }
            for element in [1, 3] {
                for &node in &target.tetrahedra()[element] {
                    let position = target.nodes_m()[node];
                    let expected = if position[2] != 0.0 || position == excluded {
                        0.0
                    } else if position == ARRIVAL_CORNER || position == RELEASED_CORNER {
                        acceleration / 7.0
                    } else {
                        3.0 * acceleration / 14.0
                    };
                    close(physical_force[node], expected, 1e-7);
                }
            }
            close(record.normal_contact_work_j, 0.0, 1e-10);
            assert!(record.absolute_contact_work_j <= 1e-10);
            close(record.contact_work_balance_error_j, 0.0, 1e-10);
        }
        let released_speed = -c / 7.0 + normal_acceleration * remainder;
        let apex_speed = b + 2.0 * c / 7.0 + 15.0 * acceleration * remainder / 7.0;
        let released_displacement =
            -c * remainder / 7.0 + 5.0 * acceleration * remainder * remainder / 7.0;
        let apex_displacement = 15.0 * acceleration * tau * tau / 14.0
            + (b + 2.0 * c / 7.0) * remainder
            + 15.0 * acceleration * remainder * remainder / 14.0;
        let mut vx = DRIFT[0];
        let mut dx = 0.0;
        for interval in &result.intervals {
            let start = interval.start_time_s - motion.time_s();
            let end = interval.end_time_s - motion.time_s();
            let mean_acceleration = 0.5 + 0.5 * (start + end) / duration;
            let next_vx = vx + (end - start) * mean_acceleration;
            dx += 0.5 * (end - start) * (vx + next_vx);
            vx = next_vx;
        }
        close(vx, DRIFT[0] + duration, 1e-14);
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element % 2 == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let position = target.nodes_m()[node];
                let (normal_speed, displacement) = if position == RELEASED_CORNER {
                    (released_speed, released_displacement)
                } else if position == ARRIVAL_CORNER {
                    (0.0, gap / 2.0)
                } else if position[2] != 0.0 {
                    (apex_speed, apex_displacement)
                } else {
                    (0.0, 0.0)
                };
                close(
                    result.motion.velocity_mps()[node][2],
                    sign * normal_speed,
                    1e-8,
                );
                close(
                    result.motion.solid().solution().displacement_m[node][2]
                        - motion.solid().solution().displacement_m[node][2],
                    sign * displacement,
                    1e-9,
                );
                close(result.motion.velocity_mps()[node][0], vx, 1e-10);
                close(
                    result.motion.solid().solution().displacement_m[node][0]
                        - motion.solid().solution().displacement_m[node][0],
                    dx,
                    1e-10,
                );
            }
        }
        assert!(minimum_gap(&result.motion) >= 0.0);
        assert!(arrival_gap(&result.motion) <= mixed_controls().0.activation_tolerance_m);
        let face = &result.motion.solid().released_faces()[0];
        let vertex = face.paired_nodes[0]
            .iter()
            .position(|&node| target.nodes_m()[node] == RELEASED_CORNER)
            .unwrap();
        let opened = result.motion.solid().solution().displacement_m[face.paired_nodes[1][vertex]]
            [2]
            - result.motion.solid().solution().displacement_m[face.paired_nodes[0][vertex]][2];
        close(opened, -2.0 * released_displacement, 1e-9);
        assert!(opened > 0.0);
        let normal_after = 0.1 * released_speed * released_speed
            + 0.1 * released_speed * apex_speed
            + 0.2 * apex_speed * apex_speed;
        let normal_change = normal_after - 0.1 * SPEED * SPEED;
        let transverse_work = 2.0 * (vx * vx - DRIFT[0] * DRIFT[0]); // total mass=4
        let loss = 3.0 * c * c / 35.0;
        let normal_work = 0.5 * acceleration * SPEED * tau
            + 10.0 * acceleration * acceleration * tau * tau / 7.0
            + acceleration * (b + 3.0 * c / 14.0) * remainder
            + 10.0 * acceleration * acceleration * remainder * remainder / 7.0;
        close(normal_work, normal_change + loss, 1e-14);
        close(result.energy.contact_dissipation_j, loss, 1e-9);
        close(
            result.energy.kinetic_change_j,
            normal_change + transverse_work,
            1e-9,
        );
        close(
            result.energy.dead_load_mean_work_j,
            normal_work + transverse_work,
            1e-9,
        );
        let before = inventory(&motion);
        let after = inventory(&result.motion);
        close(after.2 - before.2, normal_change + transverse_work, 1e-9);
        close(after.0[0] - before.0[0], 4.0 * duration, 1e-10);
        close(after.0[1], before.0[1], 1e-10);
        close(after.0[2], before.0[2], 1e-10);
        close(result.energy.normal_contact_work_j, 0.0, 1e-10);
        assert!(result.energy.absolute_normal_contact_work_j <= 1e-10);
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        assert_eq!(result.energy.support_mean_work_j, 0.0);
        assert_eq!(result.energy.support_impulse_work_j, 0.0);
        assert_eq!(
            result.motion.solid().fracture_dissipation_j(),
            motion.solid().fracture_dissipation_j()
        );
        assert_eq!(
            result.motion.solid().released_faces(),
            motion.solid().released_faces()
        );
        assert!(core::ptr::eq(
            result.motion.solid().materials().as_ptr(),
            motion.solid().materials().as_ptr()
        ));
    });
}

#[test]
fn late_mixed_interval_refusal_retains_the_original_and_valid_full_retries_match() {
    with_cx(false, |cx| {
        let geometry =
            TetLoadGeometry3::new(&SQUARE_NODES, &SQUARE_TETS, budget(), 16, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = square_topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 16, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let motion = corner_motion(&split.solid, 1.0 / 512.0, SPEED, 0.0, false, 4, cx);
        let force = force_intervals(motion.solid(), 1.0);
        let duration = 3.0 / 128.0;
        let original_u = motion.solid().solution().displacement_m.clone();
        let original_v = motion.velocity_mps().to_vec();
        let original_history = motion.solid().solution().element_states.clone();
        let limited = BrittleFaceTimeControls {
            maximum_intervals: 1,
            ..time_controls(duration)
        };
        assert!(matches!(
            motion.advance_released_face_time_with_sustained(
                duration,
                &[],
                &force,
                &[0.0, 0.0],
                ORIGIN,
                settings(),
                limited,
                mixed_controls().1,
                cx,
            ),
            Err(BrittleFaceTimeError::WorkLimit {
                intervals: 1,
                impacts: 1,
                ..
            })
        ));
        let limited = BrittleFaceTimeControls {
            maximum_trials: 1,
            ..time_controls(duration)
        };
        assert!(matches!(
            motion.advance_released_face_time_with_sustained(
                duration,
                &[],
                &force,
                &[0.0, 0.0],
                ORIGIN,
                settings(),
                limited,
                mixed_controls().1,
                cx,
            ),
            Err(BrittleFaceTimeError::SustainedApproach(
                BrittleFaceSustainedApproachError::ResponseBudget {
                    evaluations: 1,
                    maximum: 1
                }
            ))
        ));
        with_cx(true, |cancelled| {
            assert!(
                motion
                    .advance_released_face_time_with_sustained(
                        duration,
                        &[],
                        &force,
                        &[0.0, 0.0],
                        ORIGIN,
                        settings(),
                        time_controls(duration),
                        mixed_controls().1,
                        cancelled,
                    )
                    .is_err()
            );
        });
        assert_eq!(motion.solid().solution().displacement_m, original_u);
        assert_eq!(motion.velocity_mps(), original_v);
        assert_eq!(motion.solid().solution().element_states, original_history);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
        assert_eq!(motion.remaining_steps(), 4);
        let first = motion
            .advance_released_face_time_with_sustained(
                duration,
                &[],
                &force,
                &[0.0, 0.0],
                ORIGIN,
                settings(),
                time_controls(duration),
                mixed_controls().1,
                cx,
            )
            .unwrap();
        let second = motion
            .advance_released_face_time_with_sustained(
                duration,
                &[],
                &force,
                &[0.0, 0.0],
                ORIGIN,
                settings(),
                time_controls(duration),
                mixed_controls().1,
                cx,
            )
            .unwrap();
        assert_eq!(first.motion.time_s(), motion.time_s() + duration);
        assert_eq!(first.motion.velocity_mps(), second.motion.velocity_mps());
        assert_eq!(
            first.motion.solid().solution().displacement_m,
            second.motion.solid().solution().displacement_m
        );
        assert_eq!(
            first.motion.solid().solution().element_states,
            second.motion.solid().solution().element_states
        );
        assert_eq!(first.intervals, second.intervals);
        assert_eq!(first.energy, second.energy);
        assert_eq!(first.trials, second.trials);
        assert_eq!(first.sustained.len(), 2);
        for (a, b) in first.sustained.iter().zip(&second.sustained) {
            assert_eq!(a.rows, b.rows);
            assert_eq!(a.solve, b.solve);
            assert_eq!(a.native_energy, b.native_energy);
        }
        assert_eq!(
            first.impacts[0].contact.normal_impulses,
            second.impacts[0].contact.normal_impulses
        );
    });
}
