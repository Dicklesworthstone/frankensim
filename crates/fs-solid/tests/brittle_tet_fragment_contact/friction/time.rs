//! Complete native time through directional impact, sticking, and release.
//! Mass, external work, and physical contact totals have independent oracles.

use super::super::approach::approach_controls;
use super::super::sustained::approach::{
    ARRIVAL_CORNER, SQUARE_NODES, SQUARE_TETS, corner_motion, square_topology,
};
use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::friction::time::BrittleFaceFrictionTimeError;
use fs_solid::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceTimeControls, BrittleFaceTimeError,
};

const DIRECTION: [f64; 3] = [0.25, -0.125, 1.0];
const COMMON: [f64; 3] = [DRIFT[0], DRIFT[1], 0.0];

pub(super) fn time_controls(
    maximum_interval_s: f64,
) -> BrittleFaceTimeControls<BrittleFaceFrictionControls> {
    let sustained = super::sustained::sustained_controls().contact;
    let mut approach = approach_controls();
    approach.activation_tolerance_m = sustained.activation_tolerance_m;
    approach.velocity_tolerance_mps = sustained.velocity_tolerance_mps;
    let mut impact = controls();
    impact.activation_tolerance_m = approach.activation_tolerance_m;
    impact.tangential_tolerance_m = sustained.tangential_tolerance_m;
    impact.velocity_tolerance_mps = approach.velocity_tolerance_mps;
    BrittleFaceTimeControls {
        maximum_interval_s,
        maximum_intervals: 8,
        maximum_impacts: 4,
        maximum_trials: 8192,
        approach,
        impact,
        total_energy_tolerance_j: 1e-10,
    }
}

pub(super) fn force_intervals(
    motion: &BrittleTetMotion<'_>,
    direction: [f64; 3],
    amplitude: [f64; 2],
) -> Vec<(usize, [f64; 2])> {
    let geometry = motion.solid().geometry();
    let mut loads = vec![[0.0; 2]; 3 * geometry.nodes_m().len()];
    for (element, tet) in geometry.tetrahedra().iter().enumerate() {
        let sign = if element % 2 == 0 { 1.0 } else { -1.0 };
        let weight = DENSITY_FOR_FIXTURE * geometry.element_volumes_m3()[element] / 4.0;
        for &node in tet {
            for c in 0..3 {
                for end in 0..2 {
                    loads[3 * node + c][end] += sign * weight * direction[c] * amplitude[end];
                }
            }
        }
    }
    loads
        .into_iter()
        .enumerate()
        .filter(|(_, f)| *f != [0.0; 2])
        .collect()
}

const DENSITY_FOR_FIXTURE: f64 = 6.0;

fn retained_history(
    before: &BrittleTetMotion<'_>,
    after: &BrittleTetMotion<'_>,
    duration: f64,
    steps: usize,
) {
    assert_eq!(after.time_s(), before.time_s() + duration);
    assert_eq!(after.accepted_steps(), before.accepted_steps() + steps);
    assert_eq!(after.remaining_steps() + steps, before.remaining_steps());
    assert!(core::ptr::eq(
        before.solid().materials().as_ptr(),
        after.solid().materials().as_ptr()
    ));
    assert!(core::ptr::eq(
        before.reference_density_kg_m3().as_ptr(),
        after.reference_density_kg_m3().as_ptr()
    ));
    assert_eq!(after.reference_mass_kg(), before.reference_mass_kg());
    assert_eq!(
        after.solid().released_faces(),
        before.solid().released_faces()
    );
    assert_eq!(
        after.solid().fracture_dissipation_j(),
        before.solid().fracture_dissipation_j()
    );
    close(after.solid().stored_energy_j(), 0.0, 1e-10);
    for (old, next) in before
        .solid()
        .solution()
        .element_states
        .iter()
        .zip(&after.solid().solution().element_states)
    {
        assert!(next.failed());
        assert_eq!(next.damage(), old.damage());
        assert_eq!(next.dissipated_j_m3(), old.dissipated_j_m3());
        assert!(core::ptr::eq(
            next.cohesive_point().law(),
            old.cohesive_point().law()
        ));
    }
}

#[test]
fn oblique_arrival_impacts_then_carries_static_shear_over_the_exact_remainder() {
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
        let mut opened = split.solid.solution().displacement_m.clone();
        let mut velocity = vec![[0.0; 3]; opened.len()];
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                for c in 0..3 {
                    // Initially open witnesses need a tangential offset so
                    // oblique approach actually reaches the same point.
                    opened[node][c] -= sign * gap * DIRECTION[c] / 2.0;
                    velocity[node][c] = COMMON[c] + sign * SPEED * DIRECTION[c];
                }
            }
        }
        let open = split
            .solid
            .advance(&prescribed(&opened), &[], settings(), work_controls(), cx)
            .unwrap()
            .solid;
        let motion = BrittleTetMotion::try_from_solid(&open, &velocity, 2.0, 6, cx).unwrap();
        let forces = force_intervals(&motion, DIRECTION, [1.0; 2]);
        let friction = [law([0.5, 0.25], [1.0, 0.0, 0.0], [1.0, 0.0])];
        let mut limited = time_controls(duration);
        limited.maximum_intervals = 1;
        let original = motion.solid().solution().clone();
        assert!(matches!(
            motion.advance_released_face_time_frictional(
                duration,
                &[],
                &forces,
                &[0.0],
                &friction,
                ORIGIN,
                settings(),
                limited,
                super::sustained::sustained_controls(),
                cx
            ),
            Err(BrittleFaceFrictionTimeError::Time(
                BrittleFaceTimeError::WorkLimit {
                    intervals: 1,
                    impacts: 1,
                    ..
                }
            ))
        ));
        assert_eq!(
            motion.solid().solution().displacement_m,
            original.displacement_m
        );
        assert_eq!(
            motion.solid().solution().element_states,
            original.element_states
        );
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
        assert_eq!(motion.remaining_steps(), 6);

        let result = motion
            .advance_released_face_time_frictional(
                duration,
                &[],
                &forces,
                &[0.0],
                &friction,
                ORIGIN,
                settings(),
                time_controls(duration),
                super::sustained::sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.intervals.len(), 2);
        assert_eq!(result.impacts.len(), 1);
        assert_eq!(result.sustained.len(), 1);
        retained_history(&motion, &result.motion, duration, 2);
        let c = (SPEED * SPEED + gap).sqrt();
        let tau = gap / (c + SPEED);
        let remainder = duration - tau;
        let hit = &result.impacts[0];
        close(hit.location.time_s - motion.time_s(), tau, 1e-9);
        assert_eq!(result.intervals[0].end_time_s, hit.location.time_s);
        assert_eq!(result.intervals[1].start_time_s, hit.location.time_s);
        assert_eq!(result.sustained[0].start_time_s, hit.location.time_s);
        assert_eq!(result.sustained[0].end_time_s, result.motion.time_s());
        for p in &hit.contact.impulses {
            for axis in 0..3 {
                close(
                    p[axis],
                    c * [DIRECTION[2], DIRECTION[0], DIRECTION[1]][axis] / 8.0,
                    2e-8,
                );
            }
        }
        for f in &result.sustained[0].forces_n {
            for axis in 0..3 {
                close(
                    f[axis],
                    [DIRECTION[2], DIRECTION[0], DIRECTION[1]][axis] / 8.0,
                    2e-7,
                );
            }
        }
        // M=(I+ones)/20: the inelastic event leaves native apex speed2.5*c.
        // Static forces then add2.5*a*remainder; no rigid-fragment reset occurs.
        let final_scale = SPEED + duration;
        let apex_travel = 2.5 * c * remainder + 1.25 * remainder * remainder;
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let apex = target.nodes_m()[node][2] != 0.0;
                for axis in 0..3 {
                    let relative_speed = if apex {
                        sign * 2.5 * final_scale * DIRECTION[axis]
                    } else {
                        0.0
                    };
                    let relative_travel =
                        sign * (gap / 2.0 + if apex { apex_travel } else { 0.0 }) * DIRECTION[axis];
                    close(
                        result.motion.velocity_mps()[node][axis],
                        COMMON[axis] + relative_speed,
                        2e-8,
                    );
                    close(
                        result.motion.solid().solution().displacement_m[node][axis]
                            - opened[node][axis],
                        COMMON[axis] * duration + relative_travel,
                        2e-9,
                    );
                }
            }
        }
        let norm2 = dot(DIRECTION, DIRECTION);
        let loss = 3.0 * norm2 * c * c / 8.0;
        let change = norm2 * (5.0 * final_scale * final_scale / 8.0 - SPEED * SPEED);
        let before = inventory(&motion);
        let after = inventory(&result.motion);
        close(after.2 - before.2, change, 2e-9);
        close(result.energy.kinetic_change_j, change, 2e-9);
        close(result.energy.contact_dissipation_j, loss, 2e-9);
        close(result.energy.dead_load_mean_work_j, change + loss, 2e-9);
        let torque_direction = cross([0.0, 0.0, 1.0], DIRECTION);
        for axis in 0..3 {
            close(after.0[axis], before.0[axis], 2e-9);
            // A/B centroids differ by -(1+1/128)*n/2 plus a vector parallel
            // to the applied load, whose cross product contributes no moment.
            close(
                after.1[axis] - before.1[axis],
                -(129.0 / 128.0) * duration * torque_direction[axis] / 2.0,
                2e-9,
            );
        }
        close(result.energy.normal_contact_work_j, 0.0, 1e-10);
        close(result.energy.tangential_contact_work_j, 0.0, 1e-10);
        assert!(result.energy.absolute_normal_contact_work_j <= 1e-10);
        assert!(result.energy.absolute_tangential_contact_work_j <= 1e-10);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert!(result.energy.absolute_local_balance_error_j <= 1e-10);
        let retry = motion
            .advance_released_face_time_frictional(
                duration,
                &[],
                &forces,
                &[0.0],
                &friction,
                ORIGIN,
                settings(),
                time_controls(duration),
                super::sustained::sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.motion.velocity_mps(), result.motion.velocity_mps());
        assert_eq!(
            retry.motion.solid().solution().displacement_m,
            result.motion.solid().solution().displacement_m
        );
        assert_eq!(retry.impacts[0].contact, result.impacts[0].contact);
        assert_eq!(retry.sustained[0].forces_n, result.sustained[0].forces_n);
        assert_eq!(retry.energy, result.energy);
    });
}

#[test]
fn affine_directional_loading_releases_static_friction_and_finishes_in_free_flight() {
    with_oblique_motion_budget(&NODES, [0.0, 0.0, 1.0], [COMMON; 2], 8, |motion, cx| {
        let duration = 1.0 / 32.0;
        let h = duration / 4.0;
        let forces = force_intervals(motion, DIRECTION, [1.0, -1.0]);
        let friction = [law([0.5, 0.25], [1.0, 0.0, 0.0], [1.0, 0.0])];
        let result = motion
            .advance_released_face_time_frictional(
                duration,
                &[],
                &forces,
                &[0.0],
                &friction,
                ORIGIN,
                settings(),
                time_controls(h),
                super::sustained::sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(result.intervals.len(), 4);
        assert_eq!(result.sustained.len(), 3);
        assert!(result.impacts.is_empty());
        retained_history(motion, &result.motion, duration, 4);
        for (record, amplitude) in result.sustained.iter().zip([0.75, 0.25, 0.0]) {
            for f in &record.forces_n {
                for c in 0..3 {
                    close(
                        f[c],
                        amplitude * [DIRECTION[2], DIRECTION[0], DIRECTION[1]][c] / 8.0,
                        2e-8,
                    );
                }
            }
            close(record.normal_contact_work_j, 0.0, 1e-10);
            close(record.tangential_contact_work_j, 0.0, 1e-10);
        }
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let apex = motion.solid().geometry().nodes_m()[node][2] != 0.0;
                // Interval means [3/4,1/4,-1/4,-3/4] first accelerate only
                // the native apex mode, then accelerate every released node.
                let speed = sign * if apex { 1.5 * h } else { -h };
                let travel = sign
                    * if apex {
                        59.0 * h * h / 8.0
                    } else {
                        -0.75 * h * h
                    };
                for c in 0..3 {
                    close(
                        result.motion.velocity_mps()[node][c],
                        COMMON[c] + speed * DIRECTION[c],
                        2e-8,
                    );
                    close(
                        result.motion.solid().solution().displacement_m[node][c]
                            - motion.solid().solution().displacement_m[node][c],
                        COMMON[c] * duration + travel * DIRECTION[c],
                        2e-10,
                    );
                }
            }
        }
        let face = &result.motion.solid().released_faces()[0];
        for vertex in 0..3 {
            for c in 0..3 {
                let separation = result.motion.solid().solution().displacement_m
                    [face.paired_nodes[1][vertex]][c]
                    - result.motion.solid().solution().displacement_m[face.paired_nodes[0][vertex]]
                        [c];
                close(separation, 1.5 * h * h * DIRECTION[c], 2e-10);
            }
        }
        let before = inventory(motion);
        let after = inventory(&result.motion);
        let work = 3.0 * dot(DIRECTION, DIRECTION) * h * h / 8.0;
        for c in 0..3 {
            close(after.0[c], before.0[c], 2e-9);
            close(after.1[c], before.1[c], 2e-9);
        }
        close(after.2 - before.2, work, 2e-10);
        close(result.energy.dead_load_mean_work_j, work, 2e-10);
        close(result.energy.kinetic_change_j, work, 2e-10);
        assert_eq!(result.energy.contact_dissipation_j, 0.0);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        close(result.energy.normal_contact_work_j, 0.0, 1e-10);
        close(result.energy.tangential_contact_work_j, 0.0, 1e-10);
        assert!(result.energy.absolute_normal_contact_work_j <= 1e-10);
        assert!(result.energy.absolute_tangential_contact_work_j <= 1e-10);
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert!(result.energy.absolute_local_balance_error_j <= 1e-10);
    });
}

#[test]
fn mixed_arrival_retains_shared_static_shear_while_a_different_corner_releases() {
    with_cx(false, |cx| {
        let geometry =
            TetLoadGeometry3::new(&SQUARE_NODES, &SQUARE_TETS, budget(), 16, cx).unwrap();
        let materials = sources(&geometry);
        let solid = failed(&geometry, &materials, cx);
        let cut = square_topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 16, cx).unwrap();
        let split = solid
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let gap = 1.0 / 512.0;
        let duration = 3.0 / 128.0;
        let motion = corner_motion(&split.solid, gap, SPEED, 0.0, false, 4, cx);
        let mut forces = force_intervals(&motion, [0.0, 0.0, 1.0], [1.0; 2]);
        let shear = [1.0 / 32.0, -1.0 / 64.0];
        let first_face = &motion.solid().released_faces()[0];
        for vertex in 0..3 {
            let position = target.nodes_m()[first_face.paired_nodes[0][vertex]];
            if position == [0.0; 3] || position == [1.0, 1.0, 0.0] {
                for side in 0..2 {
                    let sign = if side == 0 { 1.0 } else { -1.0 };
                    for c in 0..2 {
                        forces.push((
                            3 * first_face.paired_nodes[side][vertex] + c,
                            [sign * shear[c]; 2],
                        ));
                    }
                }
            }
        }
        let friction = [law([0.5, 0.25], [1.0, 0.0, 0.0], [1.0, 0.0]); 2];
        let mut sustained = super::sustained::sustained_controls();
        sustained.contact.solve.maximum_contacts = 6;
        sustained.contact.solve.maximum_response_evaluations = 4096;
        let result = motion
            .advance_released_face_time_frictional(
                duration,
                &[],
                &forces,
                &[0.0; 2],
                &friction,
                ORIGIN,
                settings(),
                time_controls(duration),
                sustained,
                cx,
            )
            .unwrap();
        assert_eq!(result.intervals.len(), 2);
        assert_eq!(result.sustained.len(), 2);
        assert_eq!(result.impacts.len(), 1);
        retained_history(&motion, &result.motion, duration, 2);
        let acceleration = 10.0 / 7.0;
        let tau = gap / ((SPEED * SPEED + acceleration * gap).sqrt() + SPEED);
        let remainder = duration - tau;
        let c = SPEED + acceleration * tau;
        let b = 15.0 * tau / 7.0;
        assert!(remainder > 0.0 && remainder < c / 5.0);
        let hit = &result.impacts[0];
        close(hit.location.time_s - motion.time_s(), tau, 1e-9);
        assert_eq!(hit.location.face_id, 952);
        assert_eq!(hit.rows.len(), 6);
        assert!(matches!(
            hit.contact.transformed.determinacy,
            FrictionalContactDeterminacy::UniqueVelocity(_)
        ));
        let mut impulse = vec![[0.0; 3]; target.nodes_m().len()];
        for (row, p) in hit.rows.iter().zip(&hit.contact.impulses) {
            impulse[row.nodes[1]][0] += p[1];
            impulse[row.nodes[1]][1] += p[2];
            impulse[row.nodes[1]][2] += p[0];
        }
        let released_corner = [1.0, 0.0, 0.0];
        for element in [1, 3] {
            for &node in &target.tetrahedra()[element] {
                let x = target.nodes_m()[node];
                let expected = if x[2] != 0.0 || x == released_corner {
                    0.0
                } else if x == ARRIVAL_CORNER {
                    3.0 * c / 35.0
                } else {
                    c / 35.0
                };
                close(impulse[node][0], 0.0, 2e-8);
                close(impulse[node][1], 0.0, 2e-8);
                close(impulse[node][2], expected, 2e-8);
            }
        }
        for (record, excluded) in result
            .sustained
            .iter()
            .zip([ARRIVAL_CORNER, released_corner])
        {
            assert_eq!(record.rows.len(), 5);
            let mut force = vec![[0.0; 3]; target.nodes_m().len()];
            for (row, f) in record.rows.iter().zip(&record.forces_n) {
                assert_ne!(target.nodes_m()[row.nodes[0]], excluded);
                force[row.nodes[1]][0] += f[1];
                force[row.nodes[1]][1] += f[2];
                force[row.nodes[1]][2] += f[0];
            }
            // Shared corners occur twice in the row list. Only these physical
            // totals, not a numerical split of duplicate row tractions, matter.
            for element in [1, 3] {
                for &node in &target.tetrahedra()[element] {
                    let x = target.nodes_m()[node];
                    let shared = x == [0.0; 3] || x == [1.0, 1.0, 0.0];
                    let normal = if x[2] != 0.0 || x == excluded {
                        0.0
                    } else if shared {
                        3.0 / 14.0
                    } else {
                        1.0 / 7.0
                    };
                    close(force[node][2], normal, 2e-7);
                    for axis in 0..2 {
                        close(
                            force[node][axis],
                            if shared { shear[axis] } else { 0.0 },
                            2e-7,
                        );
                    }
                }
            }
            close(record.normal_contact_work_j, 0.0, 1e-10);
            close(record.tangential_contact_work_j, 0.0, 1e-10);
        }
        let released_speed = -c / 7.0 + acceleration * remainder;
        let apex_speed = b + 2.0 * c / 7.0 + 15.0 * remainder / 7.0;
        let released_travel = -c * remainder / 7.0 + 5.0 * remainder * remainder / 7.0;
        let apex_travel = 15.0 * tau * tau / 14.0
            + (b + 2.0 * c / 7.0) * remainder
            + 15.0 * remainder * remainder / 14.0;
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element % 2 == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let x = target.nodes_m()[node];
                let (speed, travel) = if x == released_corner {
                    (released_speed, released_travel)
                } else if x == ARRIVAL_CORNER {
                    (0.0, gap / 2.0)
                } else if x[2] != 0.0 {
                    (apex_speed, apex_travel)
                } else {
                    (0.0, 0.0)
                };
                close(result.motion.velocity_mps()[node][2], sign * speed, 2e-8);
                close(
                    result.motion.solid().solution().displacement_m[node][2]
                        - motion.solid().solution().displacement_m[node][2],
                    sign * travel,
                    2e-9,
                );
                for axis in 0..2 {
                    close(result.motion.velocity_mps()[node][axis], DRIFT[axis], 2e-8);
                    close(
                        result.motion.solid().solution().displacement_m[node][axis]
                            - motion.solid().solution().displacement_m[node][axis],
                        DRIFT[axis] * duration,
                        2e-9,
                    );
                }
            }
        }
        let change = 0.1 * released_speed * released_speed
            + 0.1 * released_speed * apex_speed
            + 0.2 * apex_speed * apex_speed
            - 0.1 * SPEED * SPEED;
        let loss = 3.0 * c * c / 35.0;
        let work = 0.5 * SPEED * tau
            + 10.0 * tau * tau / 7.0
            + (b + 3.0 * c / 14.0) * remainder
            + 10.0 * remainder * remainder / 7.0;
        close(work, change + loss, 1e-14);
        let before = inventory(&motion);
        let after = inventory(&result.motion);
        for axis in 0..3 {
            close(after.0[axis], before.0[axis], 2e-9);
            close(after.1[axis], before.1[axis], 2e-9);
        }
        close(after.2 - before.2, change, 2e-9);
        close(result.energy.kinetic_change_j, change, 2e-9);
        close(result.energy.contact_dissipation_j, loss, 2e-9);
        close(result.energy.dead_load_mean_work_j, work, 2e-9);
        close(result.energy.normal_contact_work_j, 0.0, 1e-10);
        close(result.energy.tangential_contact_work_j, 0.0, 1e-10);
        assert!(result.energy.absolute_normal_contact_work_j <= 1e-10);
        assert!(result.energy.absolute_tangential_contact_work_j <= 1e-10);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert!(result.energy.absolute_local_balance_error_j <= 1e-10);
    });
}
