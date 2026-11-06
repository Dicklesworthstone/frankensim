//! A complete physical duration can hold, impact, release and resume free motion.

use super::super::time::time_controls as impact_time_controls;
use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::time::BrittleFaceTimeControls;

fn time_controls(maximum_interval_s: f64) -> BrittleFaceTimeControls {
    let sustained = sustained_controls();
    let mut controls = impact_time_controls();
    controls.maximum_interval_s = maximum_interval_s;
    controls.maximum_trials = 2048;
    controls.approach.activation_tolerance_m = sustained.activation_tolerance_m;
    controls.approach.velocity_tolerance_mps = sustained.velocity_tolerance_mps;
    controls.impact.activation_tolerance_m = sustained.activation_tolerance_m;
    controls.impact.velocity_tolerance_mps = sustained.velocity_tolerance_mps;
    controls
}

#[test]
fn loaded_rest_advances_four_native_intervals_without_impulses_or_numerical_damping() {
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
        let motion = open_motion(&split.solid, 0.0, 0.0, 6, cx);
        let duration = 1.0 / 32.0;
        let interval = duration / 4.0;
        let acceleration = 1.0;
        let force: Vec<_> = body_forces(motion.solid(), acceleration)
            .into_iter()
            .map(|(dof, value)| (dof, [value, value]))
            .collect();
        let result = motion
            .advance_released_face_time_with_sustained(
                duration,
                &[],
                &force,
                &[0.0],
                ORIGIN,
                settings(),
                time_controls(interval),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.motion.time_s(), motion.time_s() + duration);
        assert_eq!(result.intervals.len(), 4);
        assert_eq!(result.sustained.len(), 4);
        assert!(result.impacts.is_empty());
        assert_eq!(result.motion.accepted_steps(), 4);
        assert_eq!(result.motion.remaining_steps(), 2);
        let response_evaluations: usize = result
            .sustained
            .iter()
            .map(|record| record.solve.response_evaluations)
            .sum();
        assert_eq!(result.trials, response_evaluations);
        for (i, record) in result.sustained.iter().enumerate() {
            assert_eq!(record.start_time_s, motion.time_s() + i as f64 * interval);
            assert_eq!(record.end_time_s, record.start_time_s + interval);
            for force in &record.solve.forces_n {
                close(force[0], acceleration / 8.0, 1e-8);
                assert_eq!(force[1], 0.0);
                assert_eq!(force[2], 0.0);
            }
            assert!(record.solve.paths_m.iter().all(|path| path[0] >= 0.0));
            assert!(record.absolute_contact_work_j <= 1e-10);
        }
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let apex = target.nodes_m()[node][2] != 0.0;
                let speed = if apex {
                    sign * 2.5 * acceleration * duration
                } else {
                    0.0
                };
                let displacement = if apex {
                    sign * 1.25 * acceleration * duration * duration
                } else {
                    0.0
                };
                close(result.motion.velocity_mps()[node][2], speed, 1e-8);
                close(
                    result.motion.solid().solution().displacement_m[node][2]
                        - motion.solid().solution().displacement_m[node][2],
                    displacement,
                    1e-10,
                );
            }
        }
        let work = 0.625 * acceleration * acceleration * duration * duration;
        close(result.energy.dead_load_mean_work_j, work, 1e-10);
        close(result.energy.kinetic_change_j, work, 1e-10);
        close(
            inventory(&result.motion).2 - inventory(&motion).2,
            work,
            1e-10,
        );
        close(result.energy.normal_contact_work_j, 0.0, 1e-10);
        assert!(result.energy.absolute_normal_contact_work_j <= 1e-10);
        assert_eq!(result.energy.contact_dissipation_j, 0.0);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert!(
            result
                .intervals
                .iter()
                .all(|record| record.energy.inertia.numerical_loss_j == 0.0)
        );
        assert!(minimum_gap(&result.motion) >= 0.0);
        assert_eq!(motion.accepted_steps(), 0);
        assert_eq!(motion.time_s(), 2.0);
    });
}

#[test]
fn accelerating_approach_resolves_one_inelastic_impact_then_holds_the_loaded_remainder() {
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
        let acceleration = 1.0;
        let motion = open_motion(&split.solid, gap, SPEED, 6, cx);
        let force: Vec<_> = body_forces(motion.solid(), acceleration)
            .into_iter()
            .map(|(dof, value)| (dof, [value, value]))
            .collect();
        let result = motion
            .advance_released_face_time_with_sustained(
                duration,
                &[],
                &force,
                &[0.0],
                ORIGIN,
                settings(),
                time_controls(duration),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.motion.time_s(), motion.time_s() + duration);
        assert_eq!(result.impacts.len(), 1);
        assert_eq!(result.sustained.len(), 1);
        assert_eq!(result.intervals.len(), 2);
        assert_eq!(result.motion.accepted_steps(), 2);
        assert_eq!(result.motion.remaining_steps(), 4);
        let impact_speed = (SPEED * SPEED + acceleration * gap).sqrt();
        let impact_time = gap / (impact_speed + SPEED);
        let hit = &result.impacts[0];
        close(hit.location.time_s - motion.time_s(), impact_time, 1e-9);
        assert_eq!(result.sustained[0].start_time_s, hit.location.time_s);
        assert_eq!(result.sustained[0].end_time_s, result.motion.time_s());
        for &impulse in &hit.contact.normal_impulses {
            close(impulse, impact_speed / 8.0, 1e-8);
        }
        for force in &result.sustained[0].solve.forces_n {
            close(force[0], acceleration / 8.0, 1e-7);
        }
        // The e=0 impulse leaves apex velocity 2.5*s. Sustained contact
        // then accelerates that native internal mode by 2.5*a. Its final
        // speed is 2.5*(v+a*T), not a stopped rigid fragment.
        let final_scale = SPEED + acceleration * duration;
        let post_impact_apex_displacement =
            1.25 * (final_scale * final_scale - impact_speed * impact_speed) / acceleration;
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let apex = target.nodes_m()[node][2] != 0.0;
                let speed = if apex { sign * 2.5 * final_scale } else { 0.0 };
                let displacement = sign
                    * (0.5 * gap
                        + if apex {
                            post_impact_apex_displacement
                        } else {
                            0.0
                        });
                close(result.motion.velocity_mps()[node][2], speed, 1e-8);
                close(
                    result.motion.solid().solution().displacement_m[node][2]
                        - motion.solid().solution().displacement_m[node][2],
                    displacement,
                    1e-9,
                );
            }
        }
        let loss = 3.0 * impact_speed * impact_speed / 8.0;
        let kinetic_change = 0.625 * final_scale * final_scale - SPEED * SPEED;
        close(result.energy.contact_dissipation_j, loss, 1e-9);
        close(result.energy.kinetic_change_j, kinetic_change, 1e-9);
        close(
            result.energy.dead_load_mean_work_j,
            kinetic_change + loss,
            1e-9,
        );
        close(
            inventory(&result.motion).2 - inventory(&motion).2,
            kinetic_change,
            1e-9,
        );
        close(result.energy.normal_contact_work_j, 0.0, 1e-10);
        assert!(result.energy.absolute_normal_contact_work_j <= 1e-10);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        assert_eq!(
            result.motion.solid().fracture_dissipation_j(),
            motion.solid().fracture_dissipation_j()
        );
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert!(minimum_gap(&result.motion) >= 0.0);
        assert!(minimum_gap(&result.motion) <= sustained_controls().activation_tolerance_m);
    });
}

#[test]
fn affine_load_reversal_releases_sustained_contact_and_finishes_in_free_flight() {
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
        let motion = open_motion(&split.solid, 0.0, 0.0, 6, cx);
        let duration = 1.0 / 32.0;
        let h = duration / 4.0;
        let force: Vec<_> = body_forces(motion.solid(), 1.0)
            .into_iter()
            .map(|(dof, value)| (dof, [value, -value]))
            .collect();
        let result = motion
            .advance_released_face_time_with_sustained(
                duration,
                &[],
                &force,
                &[0.0],
                ORIGIN,
                settings(),
                time_controls(h),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.motion.time_s(), motion.time_s() + duration);
        assert_eq!(result.intervals.len(), 4);
        assert_eq!(result.sustained.len(), 3);
        assert!(result.impacts.is_empty());
        assert_eq!(result.motion.accepted_steps(), 4);
        assert_eq!(result.motion.remaining_steps(), 2);
        // The original affine acceleration is 1-2*t/T. Its four native
        // means are [3/4,1/4,-1/4,-3/4]. The third interval releases;
        // the fourth starts with a positive gap and proceeds in free flight.
        for (record, expected_force) in result.sustained.iter().zip([0.75 / 8.0, 0.25 / 8.0, 0.0]) {
            for force in &record.solve.forces_n {
                close(force[0], expected_force, 1e-8);
            }
        }
        assert_eq!(result.sustained[2].end_time_s, motion.time_s() + 3.0 * h);
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let apex = target.nodes_m()[node][2] != 0.0;
                let speed = sign * if apex { 1.5 * h } else { -h };
                let displacement = sign
                    * if apex {
                        59.0 * h * h / 8.0
                    } else {
                        -0.75 * h * h
                    };
                close(result.motion.velocity_mps()[node][2], speed, 1e-8);
                close(
                    result.motion.solid().solution().displacement_m[node][2]
                        - motion.solid().solution().displacement_m[node][2],
                    displacement,
                    1e-10,
                );
            }
        }
        close(minimum_gap(&result.motion), 1.5 * h * h, 1e-10);
        let kinetic_change = 3.0 * h * h / 8.0;
        close(result.energy.kinetic_change_j, kinetic_change, 1e-10);
        close(result.energy.dead_load_mean_work_j, kinetic_change, 1e-10);
        close(
            inventory(&result.motion).2 - inventory(&motion).2,
            kinetic_change,
            1e-10,
        );
        assert_eq!(result.energy.contact_dissipation_j, 0.0);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        close(result.energy.normal_contact_work_j, 0.0, 1e-10);
        assert!(result.energy.absolute_normal_contact_work_j <= 1e-10);
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert!(result.energy.absolute_local_balance_error_j <= 1e-10);
    });
}
