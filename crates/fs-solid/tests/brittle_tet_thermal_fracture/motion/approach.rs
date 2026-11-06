//! G1/G3: thermal fragments localize contact with original matter and prestress.

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

/// Cooling, cutting, opening and reversal are actual accepted native steps.
/// No initializer can silently rebind inertia to the cold source's lower density.
pub(super) fn with_approaching_motion(
    maximum_steps: usize,
    run: impl FnOnce(&ThermalBrittleTetMotion<'_>, &Cx<'_>),
) {
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let reference = reference_sources(&card, &geometry);
        let cold = at_temperature(&reference, &card, 300.0);
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let initial =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 5], 3.0, maximum_steps, cx)
                .unwrap();
        let cooled = initial
            .advance(
                &cold,
                INTERVAL,
                &prescribed(&[[0.0; 3]; 5]),
                &[],
                settings(),
                backward_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let cut = topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = cooled
            .split_fracture(&cut, &target, fracture_controls(), 1e-12, cx)
            .unwrap()
            .motion;
        let mut opened = vec![[0.0; 3]; target.nodes_m().len()];
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                opened[node] = [
                    DRIFT[0] * INTERVAL,
                    DRIFT[1] * INTERVAL,
                    -sign * SPEED * INTERVAL,
                ];
            }
        }
        let outward = split
            .advance(
                &cold,
                INTERVAL,
                &prescribed(&opened),
                &[],
                settings(),
                backward_controls(),
                cx,
            )
            .unwrap()
            .motion;
        for u in &mut opened {
            for c in 0..2 {
                u[c] += DRIFT[c] * INTERVAL;
            }
        }
        let approaching = outward
            .advance_midpoint(
                &cold,
                INTERVAL,
                &prescribed(&opened),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap()
            .motion;
        assert_eq!(approaching.accepted_steps(), 3);
        assert_eq!(approaching.time_s(), 3.0 + 3.0 * INTERVAL);
        assert_eq!(approaching.reference_density_kg_m3(), &[DENSITY; 2]);
        assert!(core::ptr::eq(
            approaching.reference_density_kg_m3().as_ptr(),
            initial.reference_density_kg_m3().as_ptr()
        ));
        assert_eq!(
            approaching.solid().materials()[0]
                .band()
                .reference_density_kg_m3()
                .unwrap(),
            0.5 * DENSITY
        );
        close(approaching.reference_mass_kg(), 2.0, 1e-15);
        close(approaching.solid().stored_energy_j(), 0.04, 1e-13);
        close(approaching.solid().fracture_dissipation_j(), GC, 1e-14);
        run(&approaching, cx);
    });
}

/// Analytic loads for the two manufactured axis tetrahedra.
/// The cold transverse stress is diag(12,12,0) Pa. On each V=1/6 tetrahedron,
/// its nodal forces are (-2,-2,0) at the origin, +2x at the x corner, +2y at
/// the y corner, and zero at the apex. Balancing this real prestress keeps
/// rigid translation possible while retaining its 0.04 J of bulk storage.
/// Uniform acceleration adds (rho V/4)*a = a/4 at every corner.
pub(super) fn balanced_force_intervals(
    motion: &ThermalBrittleTetMotion<'_>,
    x_acceleration: [f64; 2],
    inward_acceleration: f64,
) -> Vec<(usize, [f64; 2])> {
    let geometry = motion.solid().geometry();
    let mut force = vec![[0.0; 2]; 3 * geometry.nodes_m().len()];
    for (element, tet) in geometry.tetrahedra().iter().enumerate() {
        let weight = DENSITY * geometry.element_volumes_m3()[element] / 4.0;
        let sign = if element == 0 { 1.0 } else { -1.0 };
        for &node in tet {
            let reference = geometry.nodes_m()[node];
            let transverse = if reference == [0.0; 3] {
                [-2.0, -2.0]
            } else if reference[0] == 1.0 {
                [2.0, 0.0]
            } else if reference[1] == 1.0 {
                [0.0, 2.0]
            } else {
                [0.0, 0.0]
            };
            for end in 0..2 {
                force[3 * node][end] += transverse[0] + weight * x_acceleration[end];
                force[3 * node + 1][end] += transverse[1];
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

pub(super) fn minimum_gap(motion: &ThermalBrittleTetMotion<'_>) -> f64 {
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

pub(super) fn retained_thermal_matter(
    before: &ThermalBrittleTetMotion<'_>,
    after: &ThermalBrittleTetMotion<'_>,
) {
    assert!(core::ptr::eq(
        before.solid().materials().as_ptr(),
        after.solid().materials().as_ptr()
    ));
    assert!(core::ptr::eq(
        before.reference_density_kg_m3().as_ptr(),
        after.reference_density_kg_m3().as_ptr()
    ));
    assert_eq!(after.reference_density_kg_m3(), &[DENSITY; 2]);
    close(after.reference_mass_kg(), 2.0, 1e-15);
    assert_eq!(
        after.solid().released_faces(),
        before.solid().released_faces()
    );
    close(after.solid().stored_energy_j(), 0.04, 1e-11);
    close(after.solid().fracture_dissipation_j(), GC, 1e-13);
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
        close(next.strain()[0], 0.01, 1e-12);
        close(next.strain()[1], 0.01, 1e-12);
    }
}

#[test]
fn thermal_contact_prefix_preserves_history_mass_and_retry_after_trial_or_cancellation_refusal() {
    with_approaching_motion(8, |motion, cx| {
        let duration = 2.0 * INTERVAL;
        let force = balanced_force_intervals(motion, [0.0; 2], 0.0);
        let snapshot = motion.solid().solution().clone();
        let velocity = motion.velocity_mps().to_vec();
        let limited = BrittleFaceApproachControls {
            maximum_trials: 1,
            ..approach_controls()
        };
        assert!(matches!(
            motion.advance_toward_released_face_contact(
                duration,
                &[],
                &force,
                settings(),
                limited,
                cx,
            ),
            Err(BrittleFaceApproachError::WorkLimit { trials: 1, .. })
        ));
        let gate = CancelGate::new_clock_free();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 599,
                    kernel_id: 599,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                motion
                    .advance_toward_released_face_contact(
                        duration,
                        &[],
                        &force,
                        settings(),
                        approach_controls(),
                        &cancelled,
                    )
                    .is_err()
            );
        });
        let located = motion
            .advance_toward_released_face_contact(
                duration,
                &[],
                &force,
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        let contact = located.contact.unwrap();
        let prefix = located.advance.unwrap();
        assert!(located.trials > 1 && located.trials <= approach_controls().maximum_trials);
        assert_eq!(contact.face_id, 598);
        close(contact.time_s, motion.time_s() + INTERVAL, 1e-9);
        assert_eq!(contact.time_s, prefix.motion.time_s());
        assert!(contact.bracket_end_time_s >= contact.time_s);
        assert!(contact.bracket_end_time_s - contact.time_s <= 1.1e-10);
        assert!(
            contact.gap_m >= 0.0 && contact.gap_m <= approach_controls().activation_tolerance_m
        );
        assert!(minimum_gap(&prefix.motion) >= 0.0);
        close(contact.relative_normal_speed_mps, -2.0 * SPEED, 1e-9);
        assert_eq!(prefix.motion.accepted_steps(), 4);
        assert_eq!(prefix.motion.remaining_steps(), 4);
        retained_thermal_matter(motion, &prefix.motion);
        close(prefix.energy.parameter_work_j, 0.0, 1e-12);
        close(prefix.energy.eigenstrain_work_j, 0.0, 1e-12);
        close(prefix.energy.external_mean_work_j(), 0.0, 1e-11);
        close(prefix.energy.total_balance_error_j, 0.0, 1e-11);
        close(
            prefix.motion.kinetic_energy_j(),
            motion.kinetic_energy_j(),
            1e-11,
        );
        let active = prefix
            .motion
            .advance_toward_released_face_contact(
                duration,
                &[],
                &force,
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        assert!(active.advance.is_none());
        assert_eq!(active.trials, 0);
        assert_eq!(active.contact.unwrap().time_s, prefix.motion.time_s());

        let safe = motion
            .advance_toward_released_face_contact(
                0.5 * INTERVAL,
                &[],
                &force,
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        assert!(safe.contact.is_none());
        assert_eq!(safe.trials, 1);
        let safe = safe.advance.unwrap().motion;
        close(minimum_gap(&safe), SPEED * INTERVAL, 1e-12);
        retained_thermal_matter(motion, &safe);
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
fn accelerated_thermal_contact_uses_the_native_prefix_mean_and_original_reference_density() {
    with_approaching_motion(6, |motion, cx| {
        let duration = 2.0 * INTERVAL;
        let acceleration = 1.0;
        let force = balanced_force_intervals(motion, [0.5, 1.5], acceleration);
        let located = motion
            .advance_toward_released_face_contact(
                duration,
                &[],
                &force,
                settings(),
                approach_controls(),
                cx,
            )
            .unwrap();
        let contact = located.contact.unwrap();
        let prefix = located.advance.unwrap();
        let h = prefix.motion.time_s() - motion.time_s();
        let gap = 2.0 * SPEED * INTERVAL;
        let root = gap / ((SPEED * SPEED + acceleration * gap).sqrt() + SPEED);
        close(h, root, 1e-9);
        let speed = SPEED + acceleration * h;
        // a_x(t)=1/2+t/duration. Velocity integrates its exact prefix mean;
        // displacement is the native midpoint expression, not exact cubic flight.
        let mean_x = 0.5 + 0.5 * h / duration;
        let vx = DRIFT[0] + h * mean_x;
        let dx = DRIFT[0] * h + 0.5 * mean_x * h * h;
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                close(prefix.motion.velocity_mps()[node][2], sign * speed, 1e-9);
                close(prefix.motion.velocity_mps()[node][0], vx, 1e-10);
                close(
                    prefix.motion.solid().solution().displacement_m[node][0]
                        - motion.solid().solution().displacement_m[node][0],
                    dx,
                    1e-10,
                );
            }
        }
        close(contact.relative_normal_speed_mps, -2.0 * speed, 1e-9);
        assert!(minimum_gap(&prefix.motion) >= 0.0);
        let work = vx * vx - DRIFT[0] * DRIFT[0] + speed * speed - SPEED * SPEED;
        close(prefix.energy.dead_load_mean_work_j, work, 1e-10);
        close(prefix.energy.inertia.kinetic_change_j, work, 1e-10);
        close(prefix.energy.total_balance_error_j, 0.0, 1e-10);
        close(prefix.energy.fracture_dissipation_j, 0.0, 1e-12);
        close(prefix.energy.parameter_work_j, 0.0, 1e-12);
        close(prefix.energy.eigenstrain_work_j, 0.0, 1e-12);
        close(
            prefix.motion.kinetic_energy_j(),
            inventory(&prefix.motion).2,
            1e-12,
        );
        retained_thermal_matter(motion, &prefix.motion);
        assert_eq!(prefix.motion.accepted_steps(), 4);
        assert_eq!(motion.accepted_steps(), 3);
    });
}
