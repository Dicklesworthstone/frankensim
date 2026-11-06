use super::*;

use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_solid::brittle::solid::dynamics::fracture::contact::{
    BrittleFaceImpactControls, SolidImpulseControls,
};
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::dynamics::{BrittleTetDynamicControls, BrittleTetDynamicError};
use fs_solid::brittle::solid::expansion::dynamics::ThermalBrittleTetMotion;
use fs_solver::CheckedCgConfig;

#[path = "motion/approach.rs"]
mod approach;
#[path = "motion/time.rs"]
mod time;

const SPEED: f64 = 1.0 / 16.0;
const DRIFT: [f64; 2] = [1.0 / 32.0, -1.0 / 64.0];
const INTERVAL: f64 = 1.0 / 64.0;
const ORIGIN: [f64; 3] = [0.125, -0.25, 0.5];

fn backward_controls() -> BrittleTetDynamicControls {
    BrittleTetDynamicControls {
        work: work_controls(),
        kinetic_work_tolerance_j: 1e-10,
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn midpoint_controls() -> BrittleTetMidpointControls {
    BrittleTetMidpointControls {
        material_energy_tolerance_j: 1e-10,
        equilibrium_work_tolerance_j: 1e-10,
        kinetic_work_tolerance_j: 1e-10,
        total_energy_tolerance_j: 1e-10,
        difference_step_strain: Some(1e-7),
    }
}

fn impact_controls() -> BrittleFaceImpactControls {
    let mut impact = MultiImpactConfig::default();
    impact.solver.gradient_tolerance = 1e-12;
    BrittleFaceImpactControls {
        impact,
        impulse: SolidImpulseControls {
            solve: CheckedCgConfig {
                tolerance: 1e-13,
                max_iterations: 64,
                max_corrections: 2,
            },
            maximum_workspace_components: 8192,
            linear_tolerance_kg_m_s: 1e-10,
            angular_tolerance_kg_m2_s: 1e-10,
            work_tolerance_j: 1e-10,
        },
        activation_tolerance_m: 1e-12,
        tangential_tolerance_m: 1e-12,
        velocity_tolerance_mps: 1e-9,
        energy_tolerance_j: 1e-10,
        maximum_contact_rows: 3,
    }
}

fn inventory(motion: &ThermalBrittleTetMotion<'_>) -> ([f64; 3], [f64; 3], f64) {
    // Positive degree-two tetrahedral quadrature is independent of native
    // consistent-mass assembly and retains all velocity cross terms.
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let geometry = motion.solid().geometry();
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    let mut kinetic = 0.0;
    for (element, tet) in geometry.tetrahedra().iter().enumerate() {
        let weight = motion.reference_density_kg_m3()[element]
            * geometry.element_volumes_m3()[element]
            / 4.0;
        for q in 0..4 {
            let mut x = [0.0; 3];
            let mut v = [0.0; 3];
            for (i, &node) in tet.iter().enumerate() {
                let shape = if i == q { a } else { b };
                for c in 0..3 {
                    x[c] += shape
                        * (geometry.nodes_m()[node][c]
                            + motion.solid().solution().displacement_m[node][c]
                            - ORIGIN[c]);
                    v[c] += shape * motion.velocity_mps()[node][c];
                }
            }
            let cross = [
                x[1] * v[2] - x[2] * v[1],
                x[2] * v[0] - x[0] * v[2],
                x[0] * v[1] - x[1] * v[0],
            ];
            for c in 0..3 {
                linear[c] += weight * v[c];
                angular[c] += weight * cross[c];
                kinetic += 0.5 * weight * v[c] * v[c];
            }
        }
    }
    (linear, angular, kinetic)
}

#[test]
fn thermal_motion_release_keeps_original_mass_and_continues_sources_with_both_opening_guards() {
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let reference = reference_sources(&card, &geometry);
        let cold = at_temperature(&reference, &card, 300.0);
        let warm = at_temperature(&reference, &card, 350.0);
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let velocity = NODES.map(|x| {
            [
                0.003 + 0.004 * x[0],
                -0.002 + 0.001 * x[1],
                0.001 + 0.0005 * x[2],
            ]
        });
        let initial =
            ThermalBrittleTetMotion::try_from_solid(&solid, &velocity, 2.0, 4, cx).unwrap();
        let endpoint = velocity.map(|v| v.map(|c| 0.01 * c));
        let cooled = initial
            .advance(
                &cold,
                0.01,
                &prescribed(&endpoint),
                &[],
                settings(),
                backward_controls(),
                cx,
            )
            .unwrap()
            .motion;
        assert!(
            cooled
                .solid()
                .solution()
                .element_states
                .iter()
                .all(|s| s.failed())
        );
        assert_eq!(
            cold[0].band().reference_density_kg_m3().unwrap(),
            0.5 * DENSITY
        );
        let before = inventory(&cooled);
        let cut = topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let handoff = cooled
            .split_fracture(&cut, &target, fracture_controls(), 1e-12, cx)
            .unwrap();
        let split = &handoff.motion;
        assert_eq!(handoff.inventory_before, handoff.inventory_after);
        assert_eq!(handoff.fragments.len(), 2);
        assert_eq!(handoff.energy.storage_defect_j, 0.0);
        assert_eq!(handoff.energy.fracture_history_defect_j, 0.0);
        close(handoff.kinetic_defect_j, 0.0, 1e-16);
        assert_eq!(split.reference_density_kg_m3(), &[DENSITY; 2]);
        assert!(core::ptr::eq(
            split.reference_density_kg_m3().as_ptr(),
            initial.reference_density_kg_m3().as_ptr()
        ));
        close(split.reference_mass_kg(), 2.0, 1e-15);
        assert_eq!(split.time_s(), cooled.time_s());
        assert_eq!(split.accepted_steps(), 1);
        assert_eq!(split.remaining_steps(), 3);
        assert!(core::ptr::eq(
            split.solid().materials().as_ptr(),
            cold.as_ptr()
        ));
        assert_eq!(
            split.solid().solution().element_states,
            cooled.solid().solution().element_states
        );
        assert_eq!(
            split.solid().solution().element_strain,
            cooled.solid().solution().element_strain
        );
        for (node, &original) in cut.original_nodes().iter().enumerate() {
            assert_eq!(split.velocity_mps()[node], cooled.velocity_mps()[original]);
            assert_eq!(
                split.solid().solution().displacement_m[node],
                cooled.solid().solution().displacement_m[original]
            );
        }
        let after = inventory(split);
        for c in 0..3 {
            close(after.0[c], before.0[c], 1e-14);
            close(after.1[c], before.1[c], 1e-14);
        }
        close(after.2, before.2, 1e-16);
        close(split.kinetic_energy_j(), after.2, 1e-16);

        let snapshot = split.solid().solution().clone();
        let delta_free = warm[0].free_linear_strain() - cold[0].free_linear_strain();
        let mut expanded = snapshot.displacement_m.clone();
        for (u, x) in expanded.iter_mut().zip(target.nodes_m()) {
            for c in 0..3 {
                u[c] += delta_free * x[c];
            }
        }
        let backward = split
            .advance(
                &warm,
                0.01,
                &prescribed(&expanded),
                &[],
                settings(),
                backward_controls(),
                cx,
            )
            .unwrap();
        let midpoint = split
            .advance_midpoint(
                &warm,
                0.01,
                &prescribed(&expanded),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap();
        for next in [&backward.motion, &midpoint.motion] {
            assert_eq!(next.reference_density_kg_m3(), &[DENSITY; 2]);
            close(next.reference_mass_kg(), 2.0, 1e-15);
            assert_eq!(next.accepted_steps(), 2);
            assert_eq!(
                next.solid().released_faces(),
                split.solid().released_faces()
            );
            assert!(core::ptr::eq(
                next.solid().materials().as_ptr(),
                warm.as_ptr()
            ));
            close(next.solid().fracture_dissipation_j(), GC, 1e-14);
            for (old, new) in snapshot
                .element_states
                .iter()
                .zip(&next.solid().solution().element_states)
            {
                for c in 0..6 {
                    close(new.strain()[c], old.strain()[c], 1e-15);
                }
                assert!(new.failed());
            }
        }
        close(backward.energy.total_balance_error_j, 0.0, 1e-10);
        close(midpoint.energy.total_balance_error_j, 0.0, 1e-10);
        assert_eq!(midpoint.energy.inertia.numerical_loss_j, 0.0);

        for &node in &target.tetrahedra()[1] {
            expanded[node][2] -= 1e-4;
        }
        assert!(matches!(
            split.advance(
                &warm,
                0.01,
                &prescribed(&expanded),
                &[],
                settings(),
                backward_controls(),
                cx,
            ),
            Err(BrittleTetDynamicError::Solid(
                BrittleTetError::ReleasedFaceClosure { face: 598, .. }
            ))
        ));
        assert!(matches!(
            split.advance_midpoint(
                &warm,
                0.01,
                &prescribed(&expanded),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            ),
            Err(BrittleTetDynamicError::Solid(
                BrittleTetError::ReleasedFaceClosure { face: 598, .. }
            ))
        ));
        assert_eq!(
            split.solid().solution().element_states,
            snapshot.element_states
        );
        assert_eq!(
            split.solid().solution().displacement_m,
            snapshot.displacement_m
        );
        assert_eq!(split.accepted_steps(), 1);
        assert_eq!(split.time_s(), cooled.time_s());
        assert!(core::ptr::eq(
            split.solid().materials().as_ptr(),
            cold.as_ptr()
        ));
    });
}

#[test]
fn thermally_fractured_body_impacts_with_retained_reference_mass_and_resumes_native_motion() {
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let reference = reference_sources(&card, &geometry);
        let cold = at_temperature(&reference, &card, 300.0);
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let initial =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 5], 3.0, 6, cx).unwrap();
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

        // Open, reverse, then arrive using actual native physical steps.
        // No motion initializer may replace the original rho=6 inventory with
        // the current source's rho=3, or reset the accepted-step count.
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
        let moving_out = split
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
        let mut reversed = opened.clone();
        for u in &mut reversed {
            for c in 0..2 {
                u[c] += DRIFT[c] * INTERVAL;
            }
        }
        let moving_in = moving_out
            .advance_midpoint(
                &cold,
                INTERVAL,
                &prescribed(&reversed),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let touching = vec![
            [3.0 * DRIFT[0] * INTERVAL, 3.0 * DRIFT[1] * INTERVAL, 0.0];
            target.nodes_m().len()
        ];
        let arrived = moving_in
            .advance_midpoint(
                &cold,
                INTERVAL,
                &prescribed(&touching),
                &[],
                settings(),
                midpoint_controls(),
                cx,
            )
            .unwrap()
            .motion;
        assert_eq!(arrived.accepted_steps(), 4);
        assert_eq!(arrived.remaining_steps(), 2);
        assert_eq!(arrived.reference_density_kg_m3(), &[DENSITY; 2]);
        close(arrived.reference_mass_kg(), 2.0, 1e-15);
        let before = inventory(&arrived);
        close(
            before.2,
            SPEED * SPEED + DRIFT[0].powi(2) + DRIFT[1].powi(2),
            1e-14,
        );
        close(arrived.solid().stored_energy_j(), 0.04, 1e-13);
        close(arrived.solid().fracture_dissipation_j(), GC, 1e-14);

        for restitution in [0.0, 0.5, 1.0] {
            let hit = arrived
                .impact_released_faces(&[restitution], &[], ORIGIN, impact_controls(), cx)
                .unwrap();
            assert_eq!(hit.rows.len(), 3);
            assert!(matches!(
                hit.contact.determinacy,
                MultiImpactDeterminacy::UniqueImpulses
            ));
            assert_eq!(
                hit.motion.solid().solution().displacement_m,
                arrived.solid().solution().displacement_m
            );
            assert_eq!(
                hit.motion.solid().solution().element_strain,
                arrived.solid().solution().element_strain
            );
            assert_eq!(
                hit.motion.solid().solution().element_states,
                arrived.solid().solution().element_states
            );
            assert_eq!(
                hit.motion.solid().solution().element_stress,
                arrived.solid().solution().element_stress
            );
            assert_eq!(
                hit.motion.solid().released_faces(),
                arrived.solid().released_faces()
            );
            assert!(core::ptr::eq(
                hit.motion.solid().materials().as_ptr(),
                cold.as_ptr()
            ));
            assert!(core::ptr::eq(
                hit.motion.reference_density_kg_m3().as_ptr(),
                initial.reference_density_kg_m3().as_ptr()
            ));
            assert_eq!(hit.motion.reference_density_kg_m3(), &[DENSITY; 2]);
            assert_eq!(hit.motion.reference_mass_kg(), arrived.reference_mass_kg());
            assert_eq!(hit.motion.time_s(), arrived.time_s());
            assert_eq!(hit.motion.accepted_steps(), 4);
            assert_eq!(hit.motion.remaining_steps(), 2);
            assert_eq!(
                hit.motion.solid().stored_energy_j(),
                arrived.solid().stored_energy_j()
            );
            assert_eq!(
                hit.motion.solid().fracture_dissipation_j(),
                arrived.solid().fracture_dissipation_j()
            );

            // Each tetrahedron has mass1, M=(I+ones)/20. The three paired
            // normal rows have W=40I-8ones; uniform eigenvalue16.
            let expected_impulse = (1.0 + restitution) * SPEED / 8.0;
            for (&impulse, &speed) in hit
                .contact
                .normal_impulses
                .iter()
                .zip(&hit.contact.relative_normal_speeds_after)
            {
                close(impulse, expected_impulse, 1e-12);
                close(speed, 2.0 * restitution * SPEED, 1e-11);
            }
            for (element, tet) in target.tetrahedra().iter().enumerate() {
                let sign = if element == 0 { 1.0 } else { -1.0 };
                for &node in tet {
                    let expected = if target.nodes_m()[node][2] == 0.0 {
                        -sign * restitution * SPEED
                    } else {
                        sign * (2.5 + 1.5 * restitution) * SPEED
                    };
                    close(hit.motion.velocity_mps()[node][2], expected, 1e-11);
                    close(hit.motion.velocity_mps()[node][0], DRIFT[0], 1e-12);
                    close(hit.motion.velocity_mps()[node][1], DRIFT[1], 1e-12);
                }
            }
            let after = inventory(&hit.motion);
            for c in 0..3 {
                close(after.0[c], before.0[c], 1e-12);
                close(after.1[c], before.1[c], 1e-12);
            }
            let loss = 3.0 * (1.0 - restitution * restitution) * SPEED * SPEED / 8.0;
            close(hit.contact.dissipated_energy, loss, 1e-12);
            close(after.2, before.2 - loss, 1e-12);
            close(hit.motion.kinetic_energy_j(), after.2, 1e-12);
            close(hit.impulse.kinetic_change_j, -loss, 1e-12);
            close(hit.relative_work_error_j, 0.0, 1e-12);
            close(hit.energy_balance_error_j, 0.0, 1e-12);
            close(hit.impulse.support_work_j, 0.0, 1e-12);
            close(hit.absolute_contact_couple_kg_m2_s, 0.0, 1e-12);

            let continued = hit
                .motion
                .advance_midpoint(
                    &cold,
                    1.0 / 4096.0,
                    &[],
                    &[],
                    settings(),
                    midpoint_controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(continued.motion.accepted_steps(), 5);
            assert_eq!(continued.motion.remaining_steps(), 1);
            assert_eq!(continued.motion.reference_density_kg_m3(), &[DENSITY; 2]);
            close(continued.motion.reference_mass_kg(), 2.0, 1e-15);
            close(continued.motion.solid().fracture_dissipation_j(), GC, 1e-14);
            close(
                continued.motion.kinetic_energy_j() + continued.motion.solid().stored_energy_j(),
                hit.motion.kinetic_energy_j() + hit.motion.solid().stored_energy_j(),
                1e-10,
            );
            close(continued.energy.total_balance_error_j, 0.0, 1e-10);
        }
    });
}
