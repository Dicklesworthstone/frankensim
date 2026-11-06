//! Complete physical intervals retain native mass, force means and fracture loss.

use super::approach::{approach_controls, minimum_gap, open_motion, rigid_force_intervals};
use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::approach::BrittleFaceApproachError;
use fs_solid::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceTimeControls, BrittleFaceTimeError,
};

fn time_controls() -> BrittleFaceTimeControls {
    let approach = approach_controls();
    let mut impact = impact_controls();
    impact.activation_tolerance_m = approach.activation_tolerance_m;
    impact.velocity_tolerance_mps = approach.velocity_tolerance_mps;
    BrittleFaceTimeControls {
        maximum_interval_s: 1.0 / 32.0,
        maximum_intervals: 8,
        maximum_impacts: 4,
        maximum_trials: 160,
        approach,
        impact,
        total_energy_tolerance_j: 1e-10,
    }
}

#[test]
fn full_ballistic_duration_resolves_consistent_mass_impact_and_retains_the_exact_remainder() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let gap = 1.0 / 512.0;
        let duration = 1.0 / 32.0;
        let motion = open_motion(&split.solid, gap, SPEED, 3, cx);
        let initial = inventory(&motion);
        let e = 0.5;
        let result = motion
            .advance_released_face_time(
                duration,
                &[],
                &[],
                &[e],
                ORIGIN,
                settings(),
                time_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.motion.time_s(), motion.time_s() + duration);
        assert_eq!(result.intervals.len(), 2);
        assert_eq!(result.impacts.len(), 1);
        assert_eq!(result.motion.accepted_steps(), 2);
        assert_eq!(result.motion.remaining_steps(), 1);
        let hit = &result.impacts[0];
        let before_hit = hit.location.time_s - motion.time_s();
        let remainder = duration - before_hit;
        close(before_hit, gap / (2.0 * SPEED), 1e-9);
        assert_eq!(result.intervals[0].end_time_s, hit.location.time_s);
        assert_eq!(result.intervals[1].start_time_s, hit.location.time_s);
        assert_eq!(result.intervals[1].end_time_s, result.motion.time_s());
        assert_eq!(hit.rows.len(), 3);
        // M=(I+ones)/20 for each unit-mass tetrahedron. Its three paired
        // normal rows have W=40I-8ones, with uniform-mode eigenvalue 16.
        let p = (1.0 + e) * SPEED / 8.0;
        for &impulse in &hit.contact.normal_impulses {
            close(impulse, p, 1e-12);
        }
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let final_speed = if target.nodes_m()[node][2] == 0.0 {
                    -sign * e * SPEED
                } else {
                    sign * (2.5 + 1.5 * e) * SPEED
                };
                close(result.motion.velocity_mps()[node][2], final_speed, 1e-10);
                close(
                    result.motion.solid().solution().displacement_m[node][2],
                    motion.solid().solution().displacement_m[node][2]
                        + sign * SPEED * before_hit
                        + final_speed * remainder,
                    1e-10,
                );
                close(result.motion.velocity_mps()[node][0], DRIFT[0], 1e-12);
                close(result.motion.velocity_mps()[node][1], DRIFT[1], 1e-12);
            }
        }
        close(
            minimum_gap(&result.motion),
            hit.location.gap_m + 2.0 * e * SPEED * remainder,
            1e-10,
        );
        assert!(minimum_gap(&result.motion) > 0.0);
        let final_inventory = inventory(&result.motion);
        for c in 0..3 {
            close(final_inventory.0[c], initial.0[c], 1e-11);
            close(final_inventory.1[c], initial.1[c], 1e-11);
        }
        let loss = 3.0 * (1.0 - e * e) * SPEED * SPEED / 8.0;
        close(result.energy.contact_dissipation_j, loss, 1e-12);
        close(result.energy.kinetic_change_j, -loss, 1e-12);
        close(final_inventory.2, initial.2 - loss, 1e-12);
        close(result.energy.kinetic_after_j, final_inventory.2, 1e-12);
        close(result.energy.external_work_j(), 0.0, 1e-12);
        close(result.energy.total_balance_error_j, 0.0, 1e-12);
        assert!(result.energy.absolute_local_balance_error_j <= 1e-12);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
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
        assert_eq!(
            result.motion.reference_density_kg_m3(),
            motion.reference_density_kg_m3()
        );
        assert_eq!(
            result.motion.reference_mass_kg(),
            motion.reference_mass_kg()
        );
        for (before, after) in motion
            .solid()
            .solution()
            .element_states
            .iter()
            .zip(&result.motion.solid().solution().element_states)
        {
            assert!(after.failed());
            assert_eq!(before.damage(), after.damage());
            assert_eq!(before.dissipated_j_m3(), after.dissipated_j_m3());
            assert!(core::ptr::eq(
                before.cohesive_point().law(),
                after.cohesive_point().law()
            ));
        }
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
        close(minimum_gap(&motion), gap, 1e-15);
    });
}

#[test]
fn affine_transverse_force_retains_global_time_means_across_localization_and_remainder() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let motion = open_motion(&split.solid, 1.0 / 512.0, SPEED, 6, cx);
        let duration = 1.0 / 32.0;
        let force = rigid_force_intervals(motion.solid(), [0.5, 1.5], 0.0);
        let controls = BrittleFaceTimeControls {
            maximum_interval_s: 3.0 / 256.0,
            ..time_controls()
        };
        let result = motion
            .advance_released_face_time(
                duration,
                &[],
                &force,
                &[0.5],
                ORIGIN,
                settings(),
                controls,
                cx,
            )
            .unwrap();
        assert_eq!(result.impacts.len(), 1);
        assert_eq!(result.intervals.len(), 4);
        assert_eq!(result.motion.accepted_steps(), result.intervals.len());
        assert_eq!(result.motion.remaining_steps(), 2);
        assert_eq!(result.motion.time_s(), motion.time_s() + duration);
        let mut vx = DRIFT[0];
        let mut dx = 0.0;
        let mut work = 0.0;
        let mut previous = motion.time_s();
        for record in &result.intervals {
            assert_eq!(record.start_time_s, previous);
            let start = record.start_time_s - motion.time_s();
            let end = record.end_time_s - motion.time_s();
            let h = end - start;
            // a_x(t)=1/2+t/T over the ORIGINAL requested interval.
            // The native midpoint displacement uses this interval's mean.
            let mean = 0.5 + 0.5 * (start + end) / duration;
            let next_vx = vx + h * mean;
            let local_work = next_vx * next_vx - vx * vx; // total mass = 2
            dx += 0.5 * h * (vx + next_vx);
            work += local_work;
            close(record.energy.dead_load_mean_work_j, local_work, 1e-11);
            close(record.energy.inertia.kinetic_change_j, local_work, 1e-11);
            close(record.energy.support_mean_work_j, 0.0, 1e-12);
            close(record.energy.total_balance_error_j, 0.0, 1e-11);
            vx = next_vx;
            previous = record.end_time_s;
        }
        close(vx, DRIFT[0] + duration, 1e-14);
        for node in 0..target.nodes_m().len() {
            close(result.motion.velocity_mps()[node][0], vx, 1e-10);
            close(
                result.motion.solid().solution().displacement_m[node][0]
                    - motion.solid().solution().displacement_m[node][0],
                dx,
                1e-10,
            );
        }
        let initial = inventory(&motion);
        let final_inventory = inventory(&result.motion);
        close(final_inventory.0[0] - initial.0[0], 2.0 * duration, 1e-11);
        let loss = 3.0 * (1.0 - 0.5 * 0.5) * SPEED * SPEED / 8.0;
        close(result.energy.dead_load_mean_work_j, work, 1e-11);
        close(result.energy.contact_dissipation_j, loss, 1e-11);
        close(result.energy.kinetic_change_j, work - loss, 1e-11);
        close(final_inventory.2 - initial.2, work - loss, 1e-11);
        close(result.energy.total_balance_error_j, 0.0, 1e-11);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        assert!(minimum_gap(&result.motion) >= 0.0);
    });
}

#[test]
fn moving_affine_support_keeps_its_velocity_and_accounts_for_impulse_work_over_full_time() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let gap = 1.0 / 512.0;
        let duration = gap / (2.0 * SPEED) + 1.0 / 512.0;
        let motion = open_motion(&split.solid, gap, SPEED, 3, cx);
        let mut supports = Vec::new();
        for &node in &target.tetrahedra()[0] {
            for axis in 0..3 {
                supports.push((
                    3 * node + axis,
                    motion.solid().solution().displacement_m[node][axis]
                        + duration * motion.velocity_mps()[node][axis],
                ));
            }
        }
        supports.sort_unstable_by_key(|&(dof, _)| dof);
        let e = 0.5;
        let result = motion
            .advance_released_face_time(
                duration,
                &supports,
                &[],
                &[e],
                ORIGIN,
                settings(),
                time_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.impacts.len(), 1);
        assert_eq!(result.intervals.len(), 2);
        assert_eq!(result.motion.time_s(), motion.time_s() + duration);
        assert_eq!(result.motion.accepted_steps(), 2);
        for &(dof, endpoint) in &supports {
            close(
                result
                    .motion
                    .solid()
                    .solution()
                    .displacement_m
                    .as_flattened()[dof],
                endpoint,
                1e-12,
            );
            close(
                result.motion.velocity_mps().as_flattened()[dof],
                motion.velocity_mps().as_flattened()[dof],
                1e-11,
            );
        }
        for &node in &target.tetrahedra()[1] {
            let speed = if target.nodes_m()[node][2] == 0.0 {
                (1.0 + 2.0 * e) * SPEED
            } else {
                -(4.0 + 3.0 * e) * SPEED
            };
            close(result.motion.velocity_mps()[node][2], speed, 1e-10);
        }
        // With A held, W=20I-4ones. A support moving at +v supplies
        // sum(p*v); the contact owner's loss is a distinct positive account.
        let p = (1.0 + e) * SPEED / 4.0;
        let support_work = 3.0 * p * SPEED;
        let loss = 3.0 * (1.0 - e * e) * SPEED * SPEED / 4.0;
        for &impulse in &result.impacts[0].contact.normal_impulses {
            close(impulse, p, 1e-11);
        }
        close(result.energy.support_impulse_work_j, support_work, 1e-11);
        close(result.energy.support_mean_work_j, 0.0, 1e-11);
        close(result.energy.dead_load_mean_work_j, 0.0, 1e-12);
        close(result.energy.contact_dissipation_j, loss, 1e-11);
        close(result.energy.kinetic_change_j, support_work - loss, 1e-11);
        close(
            inventory(&result.motion).2 - inventory(&motion).2,
            support_work - loss,
            1e-11,
        );
        close(result.energy.total_balance_error_j, 0.0, 1e-11);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        assert!(minimum_gap(&result.motion) > 0.0);
        let mut inconsistent = supports.clone();
        inconsistent[0].1 += 1.0 / 1024.0;
        assert!(matches!(
            motion.advance_released_face_time(
                duration,
                &inconsistent,
                &[],
                &[e],
                ORIGIN,
                settings(),
                time_controls(),
                cx
            ),
            Err(BrittleFaceTimeError::Invalid(_))
        ));
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
    });
}

#[test]
fn late_budget_and_persistent_contact_refusals_leave_the_full_original_interval_retryable() {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let material = sources(&geometry);
        let solid = failed(&geometry, &material, cx);
        let cut = topology(
            &geometry,
            &[InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            }],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let gap = 1.0 / 512.0;
        let duration = 1.0 / 32.0;
        let motion = open_motion(&split.solid, gap, SPEED, 4, cx);
        let original_u = motion.solid().solution().displacement_m.clone();
        let original_v = motion.velocity_mps().to_vec();
        let original_history = motion.solid().solution().element_states.clone();
        let limited = BrittleFaceTimeControls {
            maximum_intervals: 1,
            ..time_controls()
        };
        assert!(matches!(
            motion.advance_released_face_time(
                duration,
                &[],
                &[],
                &[0.5],
                ORIGIN,
                settings(),
                limited,
                cx
            ),
            Err(BrittleFaceTimeError::WorkLimit {
                intervals: 1,
                impacts: 1,
                ..
            })
        ));
        let one_step = open_motion(&split.solid, gap, SPEED, 1, cx);
        assert!(matches!(
            one_step.advance_released_face_time(
                duration,
                &[],
                &[],
                &[0.5],
                ORIGIN,
                settings(),
                time_controls(),
                cx
            ),
            Err(BrittleFaceTimeError::Approach(
                BrittleFaceApproachError::Endpoint(BrittleTetDynamicError::StepBudget { .. })
            ))
        ));
        assert_eq!(one_step.accepted_steps(), 0);
        assert_eq!(one_step.remaining_steps(), 1);

        let resting = open_motion(&split.solid, 0.0, 0.0, 4, cx);
        let inward = rigid_force_intervals(resting.solid(), [0.0, 0.0], 4.0);
        let bounded = BrittleFaceTimeControls {
            maximum_trials: 8,
            ..time_controls()
        };
        // At zero gap and zero normal speed, g(h)=-4*h^2 for every
        // positive native interval. Impulses cannot supply a contact force.
        assert!(matches!(
            resting.advance_released_face_time(
                1.0 / 1024.0,
                &[],
                &inward,
                &[0.0],
                ORIGIN,
                settings(),
                bounded,
                cx
            ),
            Err(BrittleFaceTimeError::Approach(
                BrittleFaceApproachError::WorkLimit { trials: 8, .. }
            ))
        ));
        assert_eq!(resting.time_s(), 2.0);
        assert_eq!(resting.accepted_steps(), 0);
        assert_eq!(minimum_gap(&resting), 0.0);
        with_cx(true, |cancelled| {
            assert!(
                motion
                    .advance_released_face_time(
                        duration,
                        &[],
                        &[],
                        &[0.5],
                        ORIGIN,
                        settings(),
                        time_controls(),
                        cancelled
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
            .advance_released_face_time(
                duration,
                &[],
                &[],
                &[0.5],
                ORIGIN,
                settings(),
                time_controls(),
                cx,
            )
            .unwrap();
        let second = motion
            .advance_released_face_time(
                duration,
                &[],
                &[],
                &[0.5],
                ORIGIN,
                settings(),
                time_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(first.motion.time_s(), second.motion.time_s());
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
        assert_eq!(
            first.impacts[0].contact.normal_impulses,
            second.impacts[0].contact.normal_impulses
        );
    });
}
