//! Native warming creates pressure that carries static tangential load while
//! thermal work, original reference matter, and fracture history stay distinct.

use super::super::sustained::{sustained_controls as normal_controls, with_closed_thermal_motion};
use super::*;
use fs_contact::sustained::solve::SolveError;
use fs_solid::brittle::solid::dynamics::fracture::contact::sustained::BrittleFaceSustainedError;
use fs_solid::brittle::solid::dynamics::fracture::contact::time::BrittleFaceTimeControls;
use fs_solid::brittle::solid::expansion::dynamics::fracture::contact::friction::sustained::{
    BrittleFaceSustainedFrictionControls, BrittleFaceSustainedFrictionError,
};

fn sustained_controls() -> BrittleFaceSustainedFrictionControls {
    BrittleFaceSustainedFrictionControls {
        contact: normal_controls(),
        tangent_scale_n_per_m: 4096.0,
        slip_tolerance_m: 1e-10,
    }
}

#[test]
fn warming_under_static_shear_keeps_original_matter_and_signed_thermal_support_work() {
    with_closed_thermal_motion(true, 8, |stopped, warm, cx| {
        let velocity = [DRIFT[0], DRIFT[1], 0.0];
        let mut translated = stopped.solid().solution().displacement_m.clone();
        for u in &mut translated {
            for c in 0..3 {
                u[c] += velocity[c] * INTERVAL;
            }
        }
        // Accelerate through a real accepted native step, preserving the mass
        // originally bound at 400 K rather than rebuilding it from cold rho=3.
        let moving = stopped
            .advance(
                stopped.solid().materials(),
                INTERVAL,
                &prescribed(&translated),
                &[],
                settings(),
                backward_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let geometry = moving.solid().geometry();
        let mut supports = Vec::new();
        for &node in &geometry.tetrahedra()[0] {
            for c in 0..3 {
                supports.push((3 * node + c, translated[node][c] + velocity[c] * INTERVAL));
            }
        }
        for &node in &geometry.tetrahedra()[1] {
            if geometry.nodes_m()[node][2] != 0.0 {
                for c in 0..3 {
                    supports.push((3 * node + c, translated[node][c] + velocity[c] * INTERVAL));
                }
            }
        }
        supports.sort_unstable_by_key(|&(dof, _)| dof);
        let shear = [0.125, -0.0625];
        let mut forces = Vec::new();
        // Mean transverse stress is 9 Pa during 300->350 K warming:
        // V*B^T*sigma is (-1.5,-1.5), (1.5,0), (0,1.5) at B's face.
        // Balance that real bulk stress, then load the origin tangentially.
        // Thermal normal pressure at the origin must carry this extra shear.
        for &node in &moving.solid().released_faces()[0].paired_nodes[1] {
            let x = geometry.nodes_m()[node];
            let force = if x == [0.0; 3] {
                [-1.5 - shear[0], -1.5 - shear[1]]
            } else if x[0] == 1.0 {
                [1.5, 0.0]
            } else {
                [0.0, 1.5]
            };
            for c in 0..2 {
                forces.push((3 * node + c, force[c]));
            }
        }
        let friction = [BrittleFaceFriction {
            tangent_direction: [1.0, 0.0, 0.0],
            vertices: [EllipticFriction {
                coefficients_mu: [0.5, 0.25],
                principal_axis: [1.0, 0.0],
            }; 3],
        }];

        // An exhausted nonlinear response budget cannot publish the warm
        // source, an updated history, a new clock, or consumed physical time.
        let snapshot = moving.solid().solution().clone();
        let mut limited = sustained_controls();
        limited.contact.solve.maximum_response_evaluations = 1;
        let error = moving
            .advance_released_face_sustained_frictional_at_materials(
                warm,
                INTERVAL,
                &friction,
                &supports,
                &forces,
                settings(),
                limited,
                cx,
            )
            .unwrap_err();
        let BrittleFaceSustainedFrictionError::Sustained(BrittleFaceSustainedError::Solve(error)) =
            error
        else {
            panic!("expected native response budget refusal before source transfer");
        };
        assert!(matches!(
            *error,
            SolveError::ResponseBudget { evaluations: 1, .. }
        ));
        assert_eq!(
            moving.solid().solution().element_states,
            snapshot.element_states
        );
        assert_eq!(
            moving.solid().solution().displacement_m,
            snapshot.displacement_m
        );
        assert!(core::ptr::eq(
            moving.solid().materials().as_ptr(),
            stopped.solid().materials().as_ptr()
        ));
        assert_eq!(moving.accepted_steps(), 4);
        assert_eq!(moving.remaining_steps(), 4);
        assert_eq!(moving.time_s(), stopped.time_s() + INTERVAL);
        close(
            moving.solid().materials()[0].free_linear_strain(),
            -0.01,
            1e-15,
        );

        let heated = moving
            .advance_released_face_sustained_frictional_at_materials(
                warm,
                INTERVAL,
                &friction,
                &supports,
                &forces,
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(heated.rows.len(), 3);
        for ((row, force), path) in heated
            .rows
            .iter()
            .zip(&heated.forces_n)
            .zip(&heated.paths_m)
        {
            let expected = if geometry.nodes_m()[row.nodes[0]] == [0.0; 3] {
                [0.5, shear[0], shear[1]]
            } else {
                [0.0; 3]
            };
            for c in 0..3 {
                close(force[c], expected[c], 2e-8);
                close(path[c], 0.0, 1e-10);
            }
            assert!(path[0] >= 0.0);
        }
        for (node, u) in heated
            .motion
            .solid()
            .solution()
            .displacement_m
            .iter()
            .enumerate()
        {
            for c in 0..3 {
                close(heated.motion.velocity_mps()[node][c], velocity[c], 2e-8);
                close(u[c] - translated[node][c], velocity[c] * INTERVAL, 2e-10);
            }
        }
        for (state, stress) in heated
            .motion
            .solid()
            .solution()
            .element_states
            .iter()
            .zip(&heated.motion.solid().solution().element_stress)
        {
            for c in 0..3 {
                close(state.strain()[c], [0.005, 0.005, -0.005][c], 1e-10);
                close(stress[c], [6.0, 6.0, -6.0][c], 2e-8);
            }
        }

        // The source changes elastic storage by -.025 J. Equal/opposite
        // mechanical power passes through the sticking interface; none of it
        // is relabeled contact heat or hidden in the eigenstrain account.
        let work = (shear[0] * velocity[0] + shear[1] * velocity[1]) * INTERVAL;
        assert!(work > 0.0);
        close(heated.energy.dead_load_mean_work_j, -work, 2e-10);
        close(heated.energy.support_mean_work_j, work, 2e-10);
        close(heated.energy.external_mean_work_j(), 0.0, 2e-10);
        close(heated.energy.stored_before_j, 0.04, 1e-12);
        close(heated.energy.stored_after_j, 0.015, 2e-10);
        close(heated.energy.eigenstrain_work_j, -0.025, 2e-10);
        close(heated.energy.constitutive_strain_work_j, -0.025, 2e-10);
        close(heated.energy.mechanical_work_j, 0.0, 2e-10);
        close(heated.energy.parameter_work_j, 0.0, 1e-12);
        close(heated.energy.fracture_dissipation_j, 0.0, 1e-12);
        close(heated.energy.inertia.kinetic_change_j, 0.0, 2e-10);
        close(heated.energy.total_balance_error_j, 0.0, 2e-10);
        close(heated.normal_contact_work_j, 0.0, 2e-10);
        close(heated.tangential_contact_work_j, 0.0, 2e-10);
        assert!(heated.absolute_contact_work_j <= 2e-10);
        close(heated.coordinate_work_error_j, 0.0, 2e-10);
        close(heated.contact_work_balance_error_j, 0.0, 2e-10);
        close(heated.absolute_contact_couple_n_m, 0.0, 2e-10);
        let before = inventory(&moving);
        let after = inventory(&heated.motion);
        for c in 0..3 {
            close(after.0[c], before.0[c], 2e-9);
            close(after.1[c], before.1[c], 2e-9);
        }
        close(after.2, before.2, 2e-10);
        close(after.2, velocity.iter().map(|v| v * v).sum(), 2e-10);
        close(heated.motion.kinetic_energy_j(), after.2, 2e-12);
        assert_eq!(heated.energy.inertia.numerical_loss_j, 0.0);
        assert!(core::ptr::eq(
            heated.motion.solid().materials().as_ptr(),
            warm.as_ptr()
        ));
        assert!(core::ptr::eq(
            heated.motion.reference_density_kg_m3().as_ptr(),
            stopped.reference_density_kg_m3().as_ptr()
        ));
        assert_eq!(heated.motion.reference_density_kg_m3(), &[DENSITY; 2]);
        close(heated.motion.reference_mass_kg(), 2.0, 1e-15);
        for source in heated.motion.solid().materials() {
            close(source.band().reference_density_kg_m3().unwrap(), 4.5, 1e-14);
            close(source.free_linear_strain(), -0.005, 1e-15);
        }
        assert_eq!(
            heated.motion.solid().released_faces(),
            moving.solid().released_faces()
        );
        close(heated.motion.solid().fracture_dissipation_j(), GC, 1e-13);
        for (old, new) in moving
            .solid()
            .solution()
            .element_states
            .iter()
            .zip(&heated.motion.solid().solution().element_states)
        {
            assert!(new.failed());
            assert_eq!(new.damage(), old.damage());
            assert_eq!(new.dissipated_j_m3(), old.dissipated_j_m3());
        }
        assert_eq!(heated.motion.time_s(), moving.time_s() + INTERVAL);
        assert_eq!(heated.motion.accepted_steps(), 5);
        assert_eq!(heated.motion.remaining_steps(), 3);
        let retry = moving
            .advance_released_face_sustained_frictional_at_materials(
                warm,
                INTERVAL,
                &friction,
                &supports,
                &forces,
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.forces_n, heated.forces_n);
        assert_eq!(retry.motion.velocity_mps(), heated.motion.velocity_mps());
        assert_eq!(
            retry.motion.solid().solution().element_states,
            heated.motion.solid().solution().element_states
        );
        assert_eq!(retry.energy, heated.energy);

        // Continue the warmed, moving specimen through two complete native
        // intervals. At fixed350 K the mean bulk stress is now [6,6,-6], so
        // the origin pressure is1 N, not the preceding thermal ramp's0.5 N.
        let duration = 2.0 * INTERVAL;
        let held_end: Vec<_> = supports
            .iter()
            .map(|&(dof, _)| {
                (
                    dof,
                    heated
                        .motion
                        .solid()
                        .solution()
                        .displacement_m
                        .as_flattened()[dof]
                        + velocity[dof % 3] * duration,
                )
            })
            .collect();
        let mut held_forces = Vec::new();
        for &node in &heated.motion.solid().released_faces()[0].paired_nodes[1] {
            let x = geometry.nodes_m()[node];
            let force = if x == [0.0; 3] {
                [-1.0 - shear[0], -1.0 - shear[1]]
            } else if x[0] == 1.0 {
                [1.0, 0.0]
            } else {
                [0.0, 1.0]
            };
            for axis in 0..2 {
                held_forces.push((3 * node + axis, [force[axis]; 2]));
            }
        }
        let sustained = sustained_controls();
        let mut approach = super::super::approach::approach_controls();
        approach.activation_tolerance_m = sustained.contact.activation_tolerance_m;
        approach.velocity_tolerance_mps = sustained.contact.velocity_tolerance_mps;
        let mut impact = super::controls();
        impact.activation_tolerance_m = approach.activation_tolerance_m;
        impact.tangential_tolerance_m = sustained.contact.tangential_tolerance_m;
        impact.velocity_tolerance_mps = approach.velocity_tolerance_mps;
        let time = BrittleFaceTimeControls {
            maximum_interval_s: INTERVAL,
            maximum_intervals: 4,
            maximum_impacts: 2,
            maximum_trials: 2048,
            approach,
            impact,
            total_energy_tolerance_j: 1e-10,
        };
        let continued = heated
            .motion
            .advance_released_face_time_frictional(
                duration,
                &held_end,
                &held_forces,
                &[0.0],
                &friction,
                ORIGIN,
                settings(),
                time,
                sustained,
                cx,
            )
            .unwrap();
        assert_eq!(continued.intervals.len(), 2);
        assert_eq!(continued.sustained.len(), 2);
        assert!(continued.impacts.is_empty());
        for (index, (interval, contact)) in continued
            .intervals
            .iter()
            .zip(&continued.sustained)
            .enumerate()
        {
            assert_eq!(
                interval.start_time_s,
                heated.motion.time_s() + index as f64 * INTERVAL
            );
            assert_eq!(interval.end_time_s, interval.start_time_s + INTERVAL);
            assert_eq!(contact.start_time_s, interval.start_time_s);
            assert_eq!(contact.end_time_s, interval.end_time_s);
            for (row, force) in contact.rows.iter().zip(&contact.forces_n) {
                let expected = if geometry.nodes_m()[row.nodes[0]] == [0.0; 3] {
                    [1.0, shear[0], shear[1]]
                } else {
                    [0.0; 3]
                };
                for axis in 0..3 {
                    close(force[axis], expected[axis], 2e-8);
                }
            }
            close(interval.energy.dead_load_mean_work_j, -work, 2e-10);
            close(interval.energy.support_mean_work_j, work, 2e-10);
            close(interval.energy.parameter_work_j, 0.0, 1e-12);
            close(interval.energy.eigenstrain_work_j, 0.0, 1e-12);
            close(interval.energy.total_balance_error_j, 0.0, 1e-10);
            close(contact.normal_contact_work_j, 0.0, 1e-10);
            close(contact.tangential_contact_work_j, 0.0, 1e-10);
        }
        for (node, u) in continued
            .motion
            .solid()
            .solution()
            .displacement_m
            .iter()
            .enumerate()
        {
            for axis in 0..3 {
                close(
                    continued.motion.velocity_mps()[node][axis],
                    velocity[axis],
                    2e-8,
                );
                close(
                    u[axis] - heated.motion.solid().solution().displacement_m[node][axis],
                    velocity[axis] * duration,
                    2e-10,
                );
            }
        }
        close(continued.energy.stored_before_j, 0.015, 2e-10);
        close(continued.energy.stored_after_j, 0.015, 2e-10);
        close(continued.energy.kinetic_change_j, 0.0, 2e-10);
        close(continued.energy.dead_load_mean_work_j, -2.0 * work, 2e-10);
        close(continued.energy.support_mean_work_j, 2.0 * work, 2e-10);
        close(continued.energy.external_work_j(), 0.0, 2e-10);
        close(continued.energy.parameter_work_j, 0.0, 1e-12);
        close(continued.energy.eigenstrain_work_j, 0.0, 1e-12);
        assert_eq!(continued.energy.contact_dissipation_j, 0.0);
        assert_eq!(continued.energy.fracture_dissipation_j, 0.0);
        close(continued.energy.normal_contact_work_j, 0.0, 1e-10);
        close(continued.energy.tangential_contact_work_j, 0.0, 1e-10);
        assert!(continued.energy.absolute_normal_contact_work_j <= 1e-10);
        assert!(continued.energy.absolute_tangential_contact_work_j <= 1e-10);
        close(continued.energy.total_balance_error_j, 0.0, 1e-10);
        assert!(continued.energy.absolute_local_balance_error_j <= 1e-10);
        let final_inventory = inventory(&continued.motion);
        for axis in 0..3 {
            close(final_inventory.0[axis], after.0[axis], 2e-9);
            close(final_inventory.1[axis], after.1[axis], 2e-9);
        }
        close(final_inventory.2, after.2, 2e-10);
        close(
            continued.motion.kinetic_energy_j(),
            final_inventory.2,
            2e-12,
        );
        assert!(core::ptr::eq(
            continued.motion.solid().materials().as_ptr(),
            warm.as_ptr()
        ));
        assert!(core::ptr::eq(
            continued.motion.reference_density_kg_m3().as_ptr(),
            stopped.reference_density_kg_m3().as_ptr()
        ));
        assert_eq!(continued.motion.reference_density_kg_m3(), &[DENSITY; 2]);
        close(continued.motion.reference_mass_kg(), 2.0, 1e-15);
        close(continued.motion.solid().fracture_dissipation_j(), GC, 1e-13);
        assert_eq!(
            continued.motion.solid().released_faces(),
            heated.motion.solid().released_faces()
        );
        assert_eq!(continued.motion.time_s(), heated.motion.time_s() + duration);
        assert_eq!(continued.motion.accepted_steps(), 7);
        assert_eq!(continued.motion.remaining_steps(), 1);
        assert_eq!(heated.motion.accepted_steps(), 5);
    });
}
