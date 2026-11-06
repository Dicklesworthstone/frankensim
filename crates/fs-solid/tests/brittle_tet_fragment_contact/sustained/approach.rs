//! Arrival at one corner while shared consistent-mass contact rows remain loaded.

use super::super::approach::approach_controls;
use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::approach::{
    BrittleFaceApproachControls, BrittleFaceApproachError,
};
use fs_solid::brittle::solid::dynamics::fracture::contact::sustained::BrittleFaceSustainedStep;
use fs_solid::brittle::solid::dynamics::fracture::contact::sustained::approach::BrittleFaceSustainedApproachError;

pub(in super::super) const SQUARE_NODES: [[f64; 3]; 6] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.5, 0.5, -1.0],
    [0.5, 0.5, 1.0],
];
pub(in super::super) const SQUARE_TETS: [[usize; 4]; 4] =
    [[0, 2, 1, 4], [0, 1, 2, 5], [0, 3, 2, 4], [0, 2, 3, 5]];
pub(in super::super) const ARRIVAL_CORNER: [f64; 3] = [0.0, 1.0, 0.0];

pub(in super::super) fn square_topology<'g>(
    geometry: &'g TetLoadGeometry3<'_>,
    cx: &Cx<'_>,
) -> TetFractureTopology<'g> {
    topology(
        geometry,
        &[
            InteriorSeam {
                id: 951,
                nodes: [0, 1, 2],
                side_a_element: 0,
            },
            InteriorSeam {
                id: 952,
                nodes: [0, 2, 3],
                side_a_element: 2,
            },
        ],
        cx,
    )
}

pub(super) fn mixed_controls() -> (BrittleFaceApproachControls, BrittleFaceSustainedControls) {
    let mut sustained = sustained_controls();
    sustained.solve.maximum_response_evaluations = 4096;
    let mut approach = approach_controls();
    approach.activation_tolerance_m = sustained.activation_tolerance_m;
    approach.velocity_tolerance_mps = sustained.velocity_tolerance_mps;
    (approach, sustained)
}

pub(in super::super) fn corner_motion<'a>(
    solid: &BrittleTetSolid<'a>,
    gap: f64,
    speed: f64,
    common_normal_speed: f64,
    held_lower: bool,
    maximum_steps: usize,
    cx: &Cx<'_>,
) -> BrittleTetMotion<'a> {
    let mut displacement = solid.solution().displacement_m.clone();
    let mut velocity = vec![[DRIFT[0], DRIFT[1], common_normal_speed]; displacement.len()];
    for (element, tet) in solid.geometry().tetrahedra().iter().enumerate() {
        let sign = if element % 2 == 0 { 1.0 } else { -1.0 };
        for &node in tet {
            if solid.geometry().nodes_m()[node] == ARRIVAL_CORNER {
                displacement[node][2] = solid.solution().displacement_m[node][2] - sign * gap / 2.0;
                velocity[node][2] = common_normal_speed
                    + if held_lower {
                        if sign > 0.0 { 0.0 } else { -2.0 * speed }
                    } else {
                        sign * speed
                    };
            }
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

pub(super) fn force_intervals(solid: &BrittleTetSolid<'_>, a: f64) -> Vec<(usize, [f64; 2])> {
    body_forces(solid, a)
        .into_iter()
        .map(|(dof, value)| (dof, [value, value]))
        .collect()
}

pub(super) fn arrival_gap(motion: &BrittleTetMotion<'_>) -> f64 {
    let face = motion
        .solid()
        .released_faces()
        .iter()
        .find(|face| face.id == 952)
        .unwrap();
    let vertex = face.paired_nodes[0]
        .iter()
        .position(|&node| motion.solid().geometry().nodes_m()[node] == ARRIVAL_CORNER)
        .unwrap();
    let a = motion.solid().solution().displacement_m[face.paired_nodes[0][vertex]];
    let b = motion.solid().solution().displacement_m[face.paired_nodes[1][vertex]];
    (0..3)
        .map(|axis| face.normal[axis] * (b[axis] - a[axis]))
        .sum()
}

fn assert_closed_reactions(step: &BrittleFaceSustainedStep<'_>, acceleration: f64) {
    // Shared vertices 0 and 2 each occur twice in the row list. Only their
    // physical force totals are determined; individual duplicate forces are not.
    assert_eq!(step.rows.len(), 5);
    let geometry = step.motion.solid().geometry();
    let mut force = vec![0.0; geometry.nodes_m().len()];
    for (row, values) in step.rows.iter().zip(&step.solve.forces_n) {
        assert_ne!(geometry.nodes_m()[row.nodes[0]], ARRIVAL_CORNER);
        assert_eq!(values[1], 0.0);
        assert_eq!(values[2], 0.0);
        force[row.nodes[1]] += values[0];
    }
    for &node in &geometry.tetrahedra()[1] {
        let position = geometry.nodes_m()[node];
        if position[2] != 0.0 {
            continue;
        }
        let expected = if position == [1.0, 0.0, 0.0] {
            acceleration / 7.0
        } else {
            3.0 * acceleration / 14.0
        };
        close(force[node], expected, 1e-8);
    }
}

#[test]
fn shared_mass_localizes_the_unforced_corner_while_existing_face_reactions_stay_loaded() {
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
        let motion = corner_motion(&split.solid, gap, SPEED, 0.0, false, 4, cx);
        let force = force_intervals(motion.solid(), acceleration);
        let (approach, sustained) = mixed_controls();
        let result = motion
            .advance_toward_released_face_contact_with_sustained(
                1.0 / 32.0,
                &[],
                &force,
                settings(),
                approach,
                sustained,
                cx,
            )
            .unwrap();
        let location = result.contact.unwrap();
        let prefix = result.advance.unwrap();
        let h = prefix.motion.time_s() - motion.time_s();
        // Each fragment consists of two m=1 tetrahedra. With its three
        // closed corners held by contact, the ORIGINAL free scalar mass on
        // [corner 3, apex] is [[1/10,1/20],[1/20,1/5]]. Applied loads are
        // [a/4,a/2], giving accelerations [10a/7,15a/7].
        let effective_acceleration = 10.0 * acceleration / 7.0;
        let collision_time = gap / ((SPEED * SPEED + effective_acceleration * gap).sqrt() + SPEED);
        close(h, collision_time, 1e-9);
        assert_eq!(location.face_id, 952);
        assert_eq!(location.time_s, prefix.motion.time_s());
        assert!(location.gap_m >= 0.0 && location.gap_m <= approach.activation_tolerance_m);
        assert!(location.bracket_end_time_s - location.time_s <= 1.1 * approach.time_tolerance_s);
        close(
            location.relative_normal_speed_mps,
            -2.0 * (SPEED + effective_acceleration * h),
            1e-8,
        );
        assert!(minimum_gap(&prefix.motion) >= 0.0);
        close(arrival_gap(&prefix.motion), location.gap_m, 1e-14);
        assert_closed_reactions(&prefix, acceleration);
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element % 2 == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let position = target.nodes_m()[node];
                let (speed, displacement) = if position == ARRIVAL_CORNER {
                    (
                        SPEED + effective_acceleration * h,
                        SPEED * h + 5.0 * acceleration * h * h / 7.0,
                    )
                } else if position[2] != 0.0 {
                    (
                        15.0 * acceleration * h / 7.0,
                        15.0 * acceleration * h * h / 14.0,
                    )
                } else {
                    (0.0, 0.0)
                };
                close(prefix.motion.velocity_mps()[node][2], sign * speed, 1e-8);
                close(
                    prefix.motion.solid().solution().displacement_m[node][2]
                        - motion.solid().solution().displacement_m[node][2],
                    sign * displacement,
                    1e-10,
                );
            }
        }
        let work =
            0.5 * acceleration * SPEED * h + 10.0 * acceleration * acceleration * h * h / 7.0;
        let before = inventory(&motion);
        let after = inventory(&prefix.motion);
        close(prefix.energy.dead_load_mean_work_j, work, 1e-10);
        close(prefix.energy.inertia.kinetic_change_j, work, 1e-10);
        close(after.2 - before.2, work, 1e-10);
        for axis in 0..3 {
            close(after.0[axis], before.0[axis], 1e-10);
            close(after.1[axis], before.1[axis], 1e-10);
        }
        close(prefix.normal_contact_work_j, 0.0, 1e-10);
        assert!(prefix.absolute_contact_work_j <= 1e-10);
        close(prefix.energy.total_balance_error_j, 0.0, 1e-10);
        assert_eq!(prefix.energy.fracture_dissipation_j, 0.0);
        assert_eq!(
            prefix.motion.solid().fracture_dissipation_j(),
            motion.solid().fracture_dissipation_j()
        );
        assert!(core::ptr::eq(
            prefix.motion.solid().materials().as_ptr(),
            motion.solid().materials().as_ptr()
        ));
        assert_eq!(prefix.motion.accepted_steps(), 1);
        assert_eq!(prefix.motion.remaining_steps(), 3);
        assert!(result.trials > 1 && result.trials <= approach.maximum_trials);
        assert!(result.response_evaluations >= result.trials);
        assert!(result.response_evaluations <= sustained.solve.maximum_response_evaluations);

        let safe = motion
            .advance_toward_released_face_contact_with_sustained(
                collision_time / 2.0,
                &[],
                &force,
                settings(),
                approach,
                sustained,
                cx,
            )
            .unwrap();
        assert!(safe.contact.is_none());
        let safe = safe.advance.unwrap();
        let h = safe.motion.time_s() - motion.time_s();
        close(
            arrival_gap(&safe.motion),
            gap - 2.0 * SPEED * h - effective_acceleration * h * h,
            1e-10,
        );
        assert_closed_reactions(&safe, acceleration);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
        close(arrival_gap(&motion), gap, 1e-15);
    });
}

#[test]
fn shared_moving_support_uses_the_actual_prefix_endpoint_and_retains_reaction_work() {
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
        let w = -1.0 / 32.0;
        let acceleration = 1.0;
        let duration = 1.0 / 32.0;
        let motion = corner_motion(&split.solid, gap, SPEED, w, true, 4, cx);
        let mut held: Vec<_> = [0, 2]
            .into_iter()
            .flat_map(|e| target.tetrahedra()[e])
            .collect();
        held.sort_unstable();
        held.dedup();
        let supports: Vec<_> = held
            .iter()
            .flat_map(|&node| (0..3).map(move |axis| (node, axis)))
            .map(|(node, axis)| {
                (
                    3 * node + axis,
                    motion.solid().solution().displacement_m[node][axis]
                        + duration * motion.velocity_mps()[node][axis],
                )
            })
            .collect();
        let force = force_intervals(motion.solid(), acceleration);
        let (approach, sustained) = mixed_controls();
        let result = motion
            .advance_toward_released_face_contact_with_sustained(
                duration,
                &supports,
                &force,
                settings(),
                approach,
                sustained,
                cx,
            )
            .unwrap();
        let location = result.contact.unwrap();
        let prefix = result.advance.unwrap();
        let h = prefix.motion.time_s() - motion.time_s();
        let effective_acceleration = 5.0 * acceleration / 7.0;
        let collision_time = gap / ((SPEED * SPEED + effective_acceleration * gap).sqrt() + SPEED);
        close(h, collision_time, 1e-9);
        assert_eq!(location.face_id, 952);
        assert_closed_reactions(&prefix, acceleration);
        for &node in &held {
            for axis in 0..3 {
                close(
                    prefix.motion.velocity_mps()[node][axis],
                    motion.velocity_mps()[node][axis],
                    1e-10,
                );
                close(
                    prefix.motion.solid().solution().displacement_m[node][axis]
                        - motion.solid().solution().displacement_m[node][axis],
                    h * motion.velocity_mps()[node][axis],
                    1e-10,
                );
            }
        }
        for element in [1, 3] {
            for &node in &target.tetrahedra()[element] {
                let position = target.nodes_m()[node];
                let speed = if position == ARRIVAL_CORNER {
                    w - 2.0 * SPEED - 10.0 * acceleration * h / 7.0
                } else if position[2] != 0.0 {
                    w - 15.0 * acceleration * h / 7.0
                } else {
                    w
                };
                close(prefix.motion.velocity_mps()[node][2], speed, 1e-8);
            }
        }
        // The held m=2 lower fragment has total external force +2a and
        // contact reaction -4a/7, so its support supplies -10a/7.
        let support_work = -10.0 * acceleration * w * h / 7.0;
        let dead_work =
            0.5 * acceleration * SPEED * h + 5.0 * acceleration * acceleration * h * h / 7.0;
        close(prefix.energy.support_mean_work_j, support_work, 1e-10);
        close(prefix.energy.dead_load_mean_work_j, dead_work, 1e-10);
        close(
            prefix.energy.inertia.kinetic_change_j,
            support_work + dead_work,
            1e-10,
        );
        let before = inventory(&motion);
        let after = inventory(&prefix.motion);
        close(after.2 - before.2, support_work + dead_work, 1e-10);
        close(
            after.0[2] - before.0[2],
            -10.0 * acceleration * h / 7.0,
            1e-10,
        );
        close(prefix.normal_contact_work_j, 0.0, 1e-10);
        assert!(prefix.absolute_contact_work_j <= 1e-10);
        close(prefix.contact_work_balance_error_j, 0.0, 1e-10);
        close(prefix.energy.total_balance_error_j, 0.0, 1e-10);
        assert_eq!(prefix.energy.fracture_dissipation_j, 0.0);
        assert!(minimum_gap(&prefix.motion) >= 0.0);
        assert_eq!(prefix.motion.accepted_steps(), 1);
        assert_eq!(prefix.motion.remaining_steps(), 3);
    });
}

#[test]
fn discarded_mixed_trials_and_response_budget_refuse_atomically_and_retry_identically() {
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
        let duration = 1.0 / 32.0;
        let (approach, sustained) = mixed_controls();
        let before_u = motion.solid().solution().displacement_m.clone();
        let before_v = motion.velocity_mps().to_vec();
        let before_history = motion.solid().solution().element_states.clone();
        let endpoint_limit = BrittleFaceApproachControls {
            maximum_trials: 1,
            ..approach
        };
        assert!(matches!(
            motion.advance_toward_released_face_contact_with_sustained(
                duration,
                &[],
                &force,
                settings(),
                endpoint_limit,
                sustained,
                cx,
            ),
            Err(BrittleFaceSustainedApproachError::Approach(
                BrittleFaceApproachError::WorkLimit { trials: 1, .. }
            ))
        ));
        let mut response_limit = sustained;
        response_limit.solve.maximum_response_evaluations = 1;
        assert!(matches!(
            motion.advance_toward_released_face_contact_with_sustained(
                duration,
                &[],
                &force,
                settings(),
                approach,
                response_limit,
                cx,
            ),
            Err(BrittleFaceSustainedApproachError::ResponseBudget {
                evaluations: 1,
                maximum: 1
            })
        ));
        let ready = corner_motion(&split.solid, 0.0, SPEED, 0.0, false, 4, cx)
            .advance_toward_released_face_contact_with_sustained(
                duration,
                &[],
                &force,
                settings(),
                approach,
                sustained,
                cx,
            )
            .unwrap();
        assert!(ready.advance.is_none());
        assert_eq!(ready.contact.unwrap().face_id, 952);
        assert_eq!(ready.trials, 0);
        assert_eq!(ready.response_evaluations, 0);
        let penetrated = corner_motion(&split.solid, -2.0_f64.powi(-45), SPEED, 0.0, false, 4, cx);
        assert!(matches!(
            penetrated.advance_toward_released_face_contact_with_sustained(
                duration,
                &[],
                &force,
                settings(),
                approach,
                sustained,
                cx,
            ),
            Err(BrittleFaceSustainedApproachError::Approach(
                BrittleFaceApproachError::Overlap { .. }
            ))
        ));
        with_cx(true, |cancelled| {
            assert!(
                motion
                    .advance_toward_released_face_contact_with_sustained(
                        duration,
                        &[],
                        &force,
                        settings(),
                        approach,
                        sustained,
                        cancelled,
                    )
                    .is_err()
            );
        });
        assert_eq!(motion.solid().solution().displacement_m, before_u);
        assert_eq!(motion.velocity_mps(), before_v);
        assert_eq!(motion.solid().solution().element_states, before_history);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
        assert_eq!(motion.remaining_steps(), 4);
        let first = motion
            .advance_toward_released_face_contact_with_sustained(
                duration,
                &[],
                &force,
                settings(),
                approach,
                sustained,
                cx,
            )
            .unwrap();
        let second = motion
            .advance_toward_released_face_contact_with_sustained(
                duration,
                &[],
                &force,
                settings(),
                approach,
                sustained,
                cx,
            )
            .unwrap();
        assert_eq!(first.contact, second.contact);
        assert_eq!(first.trials, second.trials);
        assert_eq!(first.response_evaluations, second.response_evaluations);
        let first = first.advance.unwrap();
        let second = second.advance.unwrap();
        assert_eq!(first.motion.velocity_mps(), second.motion.velocity_mps());
        assert_eq!(
            first.motion.solid().solution().displacement_m,
            second.motion.solid().solution().displacement_m
        );
        assert_eq!(
            first.motion.solid().solution().element_states,
            second.motion.solid().solution().element_states
        );
        assert_eq!(first.energy, second.energy);
        assert_eq!(first.native_energy, second.native_energy);
        assert_eq!(first.solve, second.solve);
    });
}
