//! G1/G3: complete thermal-fragment intervals retain source, impact and work history.

use super::approach::{
    approach_controls, balanced_force_intervals, minimum_gap, retained_thermal_matter,
    with_approaching_motion,
};
use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceTimeControls, BrittleFaceTimeError,
};

fn time_controls() -> BrittleFaceTimeControls {
    let approach = approach_controls();
    let mut impact = impact_controls();
    impact.activation_tolerance_m = approach.activation_tolerance_m;
    impact.velocity_tolerance_mps = approach.velocity_tolerance_mps;
    BrittleFaceTimeControls {
        maximum_interval_s: 2.0 * INTERVAL,
        maximum_intervals: 8,
        maximum_impacts: 4,
        maximum_trials: 160,
        approach,
        impact,
        total_energy_tolerance_j: 1e-10,
    }
}

#[test]
fn complete_thermal_duration_retries_after_late_limit_and_resolves_impact_plus_native_remainder() {
    with_approaching_motion(8, |motion, cx| {
        let duration = 2.0 * INTERVAL;
        let force = balanced_force_intervals(motion, [0.0; 2], 0.0);
        let initial = inventory(motion);
        let snapshot = motion.solid().solution().clone();
        let velocity = motion.velocity_mps().to_vec();
        let limited = BrittleFaceTimeControls {
            maximum_intervals: 1,
            ..time_controls()
        };
        // The prefix and impact can both finish internally, but the remainder
        // cannot be admitted. The entire requested interval must remain retryable.
        assert!(matches!(
            motion.advance_released_face_time(
                duration,
                &[],
                &force,
                &[0.5],
                ORIGIN,
                settings(),
                limited,
                cx,
            ),
            Err(BrittleFaceTimeError::WorkLimit {
                intervals: 1,
                impacts: 1,
                ..
            })
        ));
        assert_eq!(
            motion.solid().solution().displacement_m,
            snapshot.displacement_m
        );
        assert_eq!(
            motion.solid().solution().element_states,
            snapshot.element_states
        );
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(motion.accepted_steps(), 3);
        assert_eq!(motion.remaining_steps(), 5);

        let restitution = 0.5;
        let result = motion
            .advance_released_face_time(
                duration,
                &[],
                &force,
                &[restitution],
                ORIGIN,
                settings(),
                time_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.motion.time_s(), motion.time_s() + duration);
        assert_eq!(result.intervals.len(), 2);
        assert_eq!(result.impacts.len(), 1);
        assert_eq!(result.motion.accepted_steps(), 5);
        assert_eq!(result.motion.remaining_steps(), 3);
        assert!(result.trials <= time_controls().maximum_trials);
        let hit = &result.impacts[0];
        let before_hit = hit.location.time_s - motion.time_s();
        let remainder = duration - before_hit;
        close(before_hit, INTERVAL, 1e-9);
        assert_eq!(result.intervals[0].start_time_s, motion.time_s());
        assert_eq!(result.intervals[0].end_time_s, hit.location.time_s);
        assert_eq!(result.intervals[1].start_time_s, hit.location.time_s);
        assert_eq!(result.intervals[1].end_time_s, result.motion.time_s());
        assert_eq!(hit.rows.len(), 3);
        // Retained rho=6 gives one kg per tetrahedron, M=(I+ones)/20.
        // Thus W=40I-8ones and every normal impulse is (1+e)*S/8.
        let expected_impulse = (1.0 + restitution) * SPEED / 8.0;
        for &impulse in &hit.contact.normal_impulses {
            close(impulse, expected_impulse, 1e-12);
        }
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let final_speed = if motion.solid().geometry().nodes_m()[node][2] == 0.0 {
                    -sign * restitution * SPEED
                } else {
                    sign * (2.5 + 1.5 * restitution) * SPEED
                };
                close(result.motion.velocity_mps()[node][2], final_speed, 1e-10);
                close(
                    result.motion.solid().solution().displacement_m[node][2],
                    motion.solid().solution().displacement_m[node][2]
                        + sign * SPEED * before_hit
                        + final_speed * remainder,
                    1e-10,
                );
                close(result.motion.velocity_mps()[node][0], DRIFT[0], 1e-11);
                close(result.motion.velocity_mps()[node][1], DRIFT[1], 1e-11);
            }
        }
        close(
            minimum_gap(&result.motion),
            hit.location.gap_m + 2.0 * restitution * SPEED * remainder,
            1e-10,
        );
        assert!(minimum_gap(&result.motion) > 0.0);
        let final_inventory = inventory(&result.motion);
        for c in 0..3 {
            close(final_inventory.0[c], initial.0[c], 1e-10);
            close(final_inventory.1[c], initial.1[c], 1e-10);
        }
        let loss = 3.0 * (1.0 - restitution * restitution) * SPEED * SPEED / 8.0;
        close(result.energy.contact_dissipation_j, loss, 1e-11);
        close(result.energy.kinetic_change_j, -loss, 1e-11);
        close(final_inventory.2, initial.2 - loss, 1e-11);
        close(result.energy.kinetic_after_j, final_inventory.2, 1e-11);
        close(result.energy.external_work_j(), 0.0, 1e-11);
        close(result.energy.parameter_work_j, 0.0, 1e-12);
        close(result.energy.eigenstrain_work_j, 0.0, 1e-12);
        close(result.energy.material_work_j, 0.0, 1e-11);
        close(result.energy.fracture_dissipation_j, 0.0, 1e-12);
        close(result.energy.total_balance_error_j, 0.0, 1e-11);
        assert!(result.energy.absolute_local_balance_error_j <= 1e-10);
        for record in &result.intervals {
            close(record.energy.parameter_work_j, 0.0, 1e-12);
            close(record.energy.eigenstrain_work_j, 0.0, 1e-12);
            close(record.energy.fracture_dissipation_j, 0.0, 1e-12);
            close(record.energy.total_balance_error_j, 0.0, 1e-11);
            assert_eq!(record.energy.inertia.numerical_loss_j, 0.0);
        }
        retained_thermal_matter(motion, &result.motion);
        assert_eq!(
            motion.solid().solution().displacement_m,
            snapshot.displacement_m
        );
        assert_eq!(
            motion.solid().solution().element_states,
            snapshot.element_states
        );
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(motion.accepted_steps(), 3);
        assert_eq!(motion.time_s(), 3.0 + 3.0 * INTERVAL);
    });
}

#[test]
fn thermal_affine_force_means_use_original_global_time_across_contact_and_every_remainder() {
    with_approaching_motion(10, |motion, cx| {
        let duration = 2.0 * INTERVAL;
        let force = balanced_force_intervals(motion, [0.5, 1.5], 0.0);
        let controls = BrittleFaceTimeControls {
            maximum_interval_s: 0.75 * INTERVAL,
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
        assert_eq!(result.motion.accepted_steps(), 3 + result.intervals.len());
        assert_eq!(result.motion.remaining_steps(), 3);
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
            // Each native displacement uses its own midpoint velocity.
            let mean = 0.5 + 0.5 * (start + end) / duration;
            let next_vx = vx + h * mean;
            let local_work = next_vx * next_vx - vx * vx; // total mass=2 kg
            dx += 0.5 * h * (vx + next_vx);
            work += local_work;
            close(record.energy.dead_load_mean_work_j, local_work, 1e-10);
            close(record.energy.inertia.kinetic_change_j, local_work, 1e-10);
            close(record.energy.support_mean_work_j, 0.0, 1e-12);
            close(record.energy.mechanical_work_j, 0.0, 1e-10);
            close(record.energy.parameter_work_j, 0.0, 1e-12);
            close(record.energy.eigenstrain_work_j, 0.0, 1e-12);
            close(record.energy.total_balance_error_j, 0.0, 1e-10);
            vx = next_vx;
            previous = record.end_time_s;
        }
        close(vx, DRIFT[0] + duration, 1e-14);
        for node in 0..motion.solid().geometry().nodes_m().len() {
            close(result.motion.velocity_mps()[node][0], vx, 1e-10);
            close(
                result.motion.solid().solution().displacement_m[node][0]
                    - motion.solid().solution().displacement_m[node][0],
                dx,
                1e-10,
            );
        }
        let initial = inventory(motion);
        let final_inventory = inventory(&result.motion);
        close(final_inventory.0[0] - initial.0[0], 2.0 * duration, 1e-10);
        let loss = 3.0 * (1.0 - 0.5 * 0.5) * SPEED * SPEED / 8.0;
        close(result.energy.dead_load_mean_work_j, work, 1e-10);
        close(result.energy.support_impulse_work_j, 0.0, 1e-12);
        close(result.energy.contact_dissipation_j, loss, 1e-10);
        close(result.energy.kinetic_change_j, work - loss, 1e-10);
        close(final_inventory.2 - initial.2, work - loss, 1e-10);
        close(result.energy.parameter_work_j, 0.0, 1e-12);
        close(result.energy.eigenstrain_work_j, 0.0, 1e-12);
        close(result.energy.fracture_dissipation_j, 0.0, 1e-12);
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert!(minimum_gap(&result.motion) >= 0.0);
        retained_thermal_matter(motion, &result.motion);
        assert_eq!(motion.accepted_steps(), 3);
        assert_eq!(motion.time_s(), 3.0 + 3.0 * INTERVAL);
    });
}
