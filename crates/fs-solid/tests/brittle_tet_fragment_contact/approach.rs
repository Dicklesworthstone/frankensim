//! Native approaching prefixes use the accepted history and actual force mean.

use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::approach::{
    BrittleFaceApproachControls, BrittleFaceApproachError,
};

pub(super) fn approach_controls() -> BrittleFaceApproachControls {
    BrittleFaceApproachControls {
        work: midpoint_controls(),
        activation_tolerance_m: 1e-10,
        velocity_tolerance_mps: 1e-10,
        time_tolerance_s: 1e-10,
        maximum_trials: 80,
    }
}

pub(super) fn open_motion<'a>(
    solid: &BrittleTetSolid<'a>,
    gap: f64,
    speed: f64,
    maximum_steps: usize,
    cx: &Cx<'_>,
) -> BrittleTetMotion<'a> {
    let mut displacement = solid.solution().displacement_m.clone();
    let mut velocity = vec![[0.0; 3]; displacement.len()];
    for (element, tet) in solid.geometry().tetrahedra().iter().enumerate() {
        let sign = if element % 2 == 0 { 1.0 } else { -1.0 };
        for &node in tet {
            displacement[node][2] = solid.solution().displacement_m[node][2] - 0.5 * sign * gap;
            velocity[node] = [DRIFT[0], DRIFT[1], sign * speed];
        }
    }
    let opened = solid
        .advance(
            &prescribed(&displacement),
            &[],
            settings(),
            work_controls(),
            cx,
        )
        .unwrap()
        .solid;
    BrittleTetMotion::try_from_solid(&opened, &velocity, 2.0, maximum_steps, cx).unwrap()
}

pub(super) fn rigid_force_intervals(
    solid: &BrittleTetSolid<'_>,
    transverse_acceleration: [f64; 2],
    inward_acceleration: f64,
) -> Vec<(usize, [f64; 2])> {
    // For constant density, M times a uniform acceleration is rho*V/4 at
    // every corner. This integral does not replace M in the native solve.
    let mut force = vec![[0.0; 2]; 3 * solid.geometry().nodes_m().len()];
    for (element, tet) in solid.geometry().tetrahedra().iter().enumerate() {
        let weight = 6.0 * solid.geometry().element_volumes_m3()[element] / 4.0;
        let sign = if element % 2 == 0 { 1.0 } else { -1.0 };
        for &node in tet {
            for end in 0..2 {
                force[3 * node][end] += weight * transverse_acceleration[end];
                force[3 * node + 2][end] += weight * sign * inward_acceleration;
            }
        }
    }
    force
        .into_iter()
        .enumerate()
        .filter(|(_, values)| *values != [0.0; 2])
        .collect()
}

pub(super) fn minimum_gap(motion: &BrittleTetMotion<'_>) -> f64 {
    motion
        .solid()
        .released_faces()
        .iter()
        .flat_map(|face| {
            (0..3).map(move |vertex| {
                let a = motion.solid().solution().displacement_m[face.paired_nodes[0][vertex]];
                let b = motion.solid().solution().displacement_m[face.paired_nodes[1][vertex]];
                (0..3).map(|c| face.normal[c] * (b[c] - a[c])).sum::<f64>()
            })
        })
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn ballistic_overshoot_localizes_a_native_nonpenetrating_prefix_and_safe_duration_advances_once() {
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
        let motion = open_motion(&split.solid, gap, SPEED, 4, cx);
        let collision_time = gap / (2.0 * SPEED);
        let requested = 2.0 * collision_time;
        assert!(
            motion
                .advance_midpoint(requested, &[], &[], settings(), midpoint_controls(), cx)
                .is_err()
        );
        let result = motion
            .advance_toward_released_face_contact(
                requested,
                &[],
                &[],
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        assert!(result.trials <= approach_controls().maximum_trials);
        let location = result.contact.unwrap();
        let prefix = result.advance.unwrap();
        assert_eq!(location.face_id, 951);
        close(location.time_s, motion.time_s() + collision_time, 1e-9);
        assert_eq!(location.time_s, prefix.motion.time_s());
        assert!(location.bracket_end_time_s >= location.time_s);
        assert!(location.bracket_end_time_s - location.time_s <= 1.1e-10);
        assert!(
            location.gap_m >= 0.0 && location.gap_m <= approach_controls().activation_tolerance_m
        );
        assert!(minimum_gap(&prefix.motion) >= 0.0);
        close(location.relative_normal_speed_mps, -2.0 * SPEED, 1e-10);
        assert_eq!(prefix.motion.accepted_steps(), 1);
        assert_eq!(prefix.motion.remaining_steps(), 3);
        close(
            prefix.motion.kinetic_energy_j(),
            motion.kinetic_energy_j(),
            1e-12,
        );
        assert_eq!(
            prefix.motion.solid().fracture_dissipation_j(),
            motion.solid().fracture_dissipation_j()
        );
        close(prefix.energy.external_mean_work_j(), 0.0, 1e-12);
        close(prefix.energy.total_balance_error_j, 0.0, 1e-12);

        let before = motion
            .advance_toward_released_face_contact(
                0.5 * collision_time,
                &[],
                &[],
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        assert!(before.contact.is_none());
        let safe = before.advance.unwrap().motion;
        assert_eq!(safe.accepted_steps(), 1);
        close(safe.time_s(), motion.time_s() + 0.5 * collision_time, 1e-14);
        close(minimum_gap(&safe), 0.5 * gap, 1e-12);
        assert_eq!(motion.accepted_steps(), 0);
        assert_eq!(motion.time_s(), 2.0);
        close(minimum_gap(&motion), gap, 1e-15);
    });
}

#[test]
fn accelerated_contact_root_and_affine_transverse_load_use_the_actual_prefix_mean() {
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
        let motion = open_motion(&split.solid, gap, SPEED, 4, cx);
        let acceleration = 1.0;
        let duration = 1.0 / 32.0;
        let force = rigid_force_intervals(motion.solid(), [0.5, 1.5], acceleration);
        let result = motion
            .advance_toward_released_face_contact(
                duration,
                &[],
                &force,
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        let location = result.contact.unwrap();
        let prefix = result.advance.unwrap();
        let h = prefix.motion.time_s() - motion.time_s();
        // Constant inward acceleration gives this exact midpoint gap family.
        let expected_root = gap / ((SPEED * SPEED + acceleration * gap).sqrt() + SPEED);
        close(h, expected_root, 1e-9);
        let normal_speed = SPEED + acceleration * h;
        let mean_x_acceleration = 0.5 + 0.5 * h / duration;
        let expected_vx = DRIFT[0] + h * mean_x_acceleration;
        let expected_dx = DRIFT[0] * h + 0.5 * mean_x_acceleration * h * h;
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                close(
                    prefix.motion.velocity_mps()[node][2],
                    sign * normal_speed,
                    1e-9,
                );
                close(prefix.motion.velocity_mps()[node][0], expected_vx, 1e-10);
                close(
                    prefix.motion.solid().solution().displacement_m[node][0]
                        - motion.solid().solution().displacement_m[node][0],
                    expected_dx,
                    1e-10,
                );
            }
        }
        close(
            location.relative_normal_speed_mps,
            -2.0 * normal_speed,
            1e-9,
        );
        assert!(minimum_gap(&prefix.motion) >= 0.0);
        let expected_work = (expected_vx * expected_vx - DRIFT[0] * DRIFT[0])
            + (normal_speed * normal_speed - SPEED * SPEED);
        close(prefix.energy.dead_load_mean_work_j, expected_work, 1e-10);
        close(prefix.energy.inertia.kinetic_change_j, expected_work, 1e-10);
        close(
            prefix.motion.kinetic_energy_j(),
            inventory(&prefix.motion).2,
            1e-12,
        );
        close(prefix.energy.total_balance_error_j, 0.0, 1e-10);
        assert_eq!(prefix.energy.fracture_dissipation_j, 0.0);
    });
}

#[test]
fn initial_contact_strict_overlap_and_trial_limits_preserve_the_original_motion_for_retry() {
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
        let duration = 1.0 / 32.0;
        let closing = open_motion(&split.solid, 0.0, SPEED, 4, cx);
        let ready = closing
            .advance_toward_released_face_contact(
                duration,
                &[],
                &[],
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        assert!(ready.advance.is_none());
        let contact = ready.contact.unwrap();
        assert_eq!(contact.time_s, closing.time_s());
        assert_eq!(contact.bracket_end_time_s, closing.time_s());
        assert_eq!(closing.accepted_steps(), 0);
        let separating = open_motion(&split.solid, 0.0, -SPEED, 4, cx);
        let progress = separating
            .advance_toward_released_face_contact(
                duration,
                &[],
                &[],
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        assert!(progress.contact.is_none());
        close(
            minimum_gap(&progress.advance.unwrap().motion),
            2.0 * SPEED * duration,
            1e-12,
        );
        let penetrated = open_motion(&split.solid, -2.0_f64.powi(-45), SPEED, 4, cx);
        assert!(matches!(
            penetrated.advance_toward_released_face_contact(
                duration,
                &[],
                &[],
                settings(),
                approach_controls(),
                cx
            ),
            Err(BrittleFaceApproachError::Overlap { .. })
        ));

        let motion = open_motion(&split.solid, 1.0 / 512.0, SPEED, 4, cx);
        let saved_velocity = motion.velocity_mps().to_vec();
        let saved_displacement = motion.solid().solution().displacement_m.clone();
        let saved_history = motion.solid().solution().element_states.clone();
        let limited = BrittleFaceApproachControls {
            maximum_trials: 1,
            ..approach_controls()
        };
        assert!(matches!(
            motion.advance_toward_released_face_contact(
                duration,
                &[],
                &[],
                settings(),
                limited,
                cx
            ),
            Err(BrittleFaceApproachError::WorkLimit { .. })
        ));
        for invalid in [
            BrittleFaceApproachControls {
                activation_tolerance_m: 0.0,
                ..approach_controls()
            },
            BrittleFaceApproachControls {
                time_tolerance_s: f64::NAN,
                ..approach_controls()
            },
        ] {
            assert!(
                motion
                    .advance_toward_released_face_contact(
                        duration,
                        &[],
                        &[],
                        settings(),
                        invalid,
                        cx
                    )
                    .is_err()
            );
        }
        with_cx(true, |cancelled| {
            assert!(
                motion
                    .advance_toward_released_face_contact(
                        duration,
                        &[],
                        &[],
                        settings(),
                        approach_controls(),
                        cancelled
                    )
                    .is_err()
            );
        });
        assert_eq!(motion.velocity_mps(), saved_velocity);
        assert_eq!(motion.solid().solution().displacement_m, saved_displacement);
        assert_eq!(motion.solid().solution().element_states, saved_history);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
        let first = motion
            .advance_toward_released_face_contact(
                duration,
                &[],
                &[],
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        let second = motion
            .advance_toward_released_face_contact(
                duration,
                &[],
                &[],
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        let a = first.advance.unwrap().motion;
        let b = second.advance.unwrap().motion;
        assert_eq!(a.velocity_mps(), b.velocity_mps());
        assert_eq!(
            a.solid().solution().displacement_m,
            b.solid().solution().displacement_m
        );
        assert_eq!(
            a.solid().solution().element_states,
            b.solid().solution().element_states
        );
        assert_eq!(a.time_s(), b.time_s());
        assert_eq!(first.trials, second.trials);
    });
}
