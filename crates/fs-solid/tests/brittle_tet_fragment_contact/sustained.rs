//! Sustained normal forces use the native consistent mass and recovered compression.

use super::approach::{minimum_gap, open_motion, rigid_force_intervals};
use super::*;
use fs_contact::sustained::solve::{SolveControls, SolveError};
use fs_solid::brittle::solid::dynamics::fracture::contact::sustained::{
    BrittleFaceSustainedControls, BrittleFaceSustainedError,
};
use fs_solver::NewtonKrylovConfig;

#[path = "sustained/approach.rs"]
mod approach;
#[path = "sustained/time.rs"]
mod time;

fn sustained_controls() -> BrittleFaceSustainedControls {
    BrittleFaceSustainedControls {
        solve: SolveControls {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            maximum_response_evaluations: 512,
            difference_step_n: 1e-3,
            force_residual_tolerance_n: 1e-9,
            maximum_contacts: 12,
            ..SolveControls::default()
        },
        work: midpoint_controls(),
        // A conditioning scale in N/m, independent of material stiffness.
        normal_scale_n_per_m: 4096.0,
        activation_tolerance_m: 1e-10,
        tangential_tolerance_m: 1e-10,
        velocity_tolerance_mps: 1e-8,
        force_tolerance_n: 1e-8,
        complementarity_tolerance_j: 1e-10,
        contact_work_tolerance_j: 1e-10,
        force_moment_tolerance_n_m: 1e-10,
    }
}

fn body_forces(solid: &BrittleTetSolid<'_>, inward_acceleration: f64) -> Vec<(usize, f64)> {
    rigid_force_intervals(solid, [0.0, 0.0], inward_acceleration)
        .into_iter()
        .map(|(dof, values)| {
            assert_eq!(values[0], values[1]);
            (dof, values[0])
        })
        .collect()
}

#[test]
fn closed_failed_fragments_keep_consistent_mass_motion_and_recover_native_compression_storage() {
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
        let zero = vec![[0.0; 3]; target.nodes_m().len()];
        let unloaded = split
            .solid
            .advance(&prescribed(&zero), &[], settings(), work_controls(), cx)
            .unwrap()
            .solid;
        let h = 1.0 / 64.0;
        let acceleration = 1.0;
        for (base, compression) in [(&split.solid, false), (&unloaded, true)] {
            let motion = open_motion(base, 0.0, 0.0, 3, cx);
            let result = motion
                .advance_released_face_sustained(
                    h,
                    &[],
                    &body_forces(motion.solid(), acceleration),
                    settings(),
                    sustained_controls(),
                    cx,
                )
                .unwrap();
            // M=(I+ones)/20, E=1024, nu=0, V=1/6. Closed symmetric
            // faces stay at z=0. The two apex displacements are +/-d.
            // For compression, native mean stress supplies E*d/2, so
            // (1/(5*h^2)+E/12)*d=a/4. A rigid-body or mass-only response
            // cannot reproduce both this compliance and the stored energy.
            let young = if compression { 1024.0 } else { 0.0 };
            let denominator = 1.0 + 5.0 * young * h * h / 12.0;
            let d = 1.25 * acceleration * h * h / denominator;
            assert_eq!(result.rows.len(), 3);
            assert!(result.solve.nonlinear.converged);
            for (row, force) in result.rows.iter().zip(&result.solve.forces_n) {
                let position = target.nodes_m()[row.nodes[0]];
                let expected = if position == [0.0; 3] {
                    acceleration / 2.0 - 3.0 * acceleration / (8.0 * denominator)
                } else {
                    acceleration / 4.0 - acceleration / (8.0 * denominator)
                };
                close(force[0], expected, 1e-8);
                assert_eq!(force[1], 0.0);
                assert_eq!(force[2], 0.0);
            }
            for (element, tet) in target.tetrahedra().iter().enumerate() {
                let sign = if element == 0 { 1.0 } else { -1.0 };
                for &node in tet {
                    let displacement = if target.nodes_m()[node][2] == 0.0 {
                        0.0
                    } else {
                        sign * d
                    };
                    close(
                        result.motion.solid().solution().displacement_m[node][2]
                            - motion.solid().solution().displacement_m[node][2],
                        displacement,
                        1e-10,
                    );
                    close(
                        result.motion.velocity_mps()[node][2],
                        2.0 * displacement / h,
                        1e-8,
                    );
                }
                close(
                    result.motion.solid().solution().element_stress[element][2],
                    -young * d,
                    1e-8,
                );
                let before = &motion.solid().solution().element_states[element];
                let after = &result.motion.solid().solution().element_states[element];
                assert!(after.failed());
                assert_eq!(after.damage(), before.damage());
                assert_eq!(after.dissipated_j_m3(), before.dissipated_j_m3());
                assert!(core::ptr::eq(
                    after.cohesive_point().law(),
                    before.cohesive_point().law()
                ));
            }
            let storage = young * d * d / 6.0;
            let kinetic = 0.625 * (acceleration * h / denominator).powi(2);
            let work = 0.5 * acceleration * d;
            close(result.energy.stored_after_j, storage, 1e-10);
            close(result.energy.inertia.kinetic_change_j, kinetic, 1e-10);
            close(
                inventory(&result.motion).2 - inventory(&motion).2,
                kinetic,
                1e-10,
            );
            close(result.energy.dead_load_mean_work_j, work, 1e-10);
            close(
                result.energy.external_mean_work_j(),
                storage + kinetic,
                1e-10,
            );
            close(result.normal_contact_work_j, 0.0, 1e-10);
            assert!(result.absolute_contact_work_j <= 1e-10);
            close(result.contact_work_balance_error_j, 0.0, 1e-10);
            close(result.absolute_contact_couple_n_m, 0.0, 1e-10);
            close(result.energy.total_balance_error_j, 0.0, 1e-10);
            close(
                result.native_energy.dead_load_mean_work_j,
                result.energy.dead_load_mean_work_j + result.normal_contact_work_j,
                1e-10,
            );
            assert_eq!(result.energy.inertia.numerical_loss_j, 0.0);
            assert_eq!(result.energy.fracture_dissipation_j, 0.0);
            assert_eq!(
                result.motion.solid().fracture_dissipation_j(),
                motion.solid().fracture_dissipation_j()
            );
            assert!(minimum_gap(&result.motion) >= 0.0);
            assert!(minimum_gap(&result.motion) <= sustained_controls().activation_tolerance_m);
            assert_eq!(result.motion.time_s(), motion.time_s() + h);
            assert_eq!(result.motion.accepted_steps(), 1);
            assert_eq!(result.motion.remaining_steps(), 2);
            assert_eq!(motion.accepted_steps(), 0);
            assert_eq!(motion.time_s(), 2.0);
        }
    });
}

#[test]
fn outward_loading_releases_a_previously_loaded_face_without_damping_or_gluing() {
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
        let motion = open_motion(&split.solid, 0.0, 0.0, 3, cx);
        let h = 1.0 / 64.0;
        let acceleration = 1.0;
        let loaded = motion
            .advance_released_face_sustained(
                h,
                &[],
                &body_forces(motion.solid(), acceleration),
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert!(loaded.solve.forces_n.iter().all(|force| force[0] > 0.0));
        let released = loaded
            .motion
            .advance_released_face_sustained(
                h,
                &[],
                &body_forces(loaded.motion.solid(), -acceleration),
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert!(
            released
                .solve
                .forces_n
                .iter()
                .all(|force| *force == [0.0; 3])
        );
        close(minimum_gap(&released.motion), acceleration * h * h, 1e-10);
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let expected = if target.nodes_m()[node][2] == 0.0 {
                    -sign * acceleration * h
                } else {
                    sign * 1.5 * acceleration * h
                };
                close(released.motion.velocity_mps()[node][2], expected, 1e-8);
            }
        }
        // Previous apex velocity is +/-2.5*a*h. Outward force extracts
        // a^2*h^2/4 while the faces separate; this is external work,
        // without a contact loss or a numerical kinetic damping term.
        let work = -0.25 * acceleration * acceleration * h * h;
        close(released.energy.dead_load_mean_work_j, work, 1e-10);
        close(released.energy.inertia.kinetic_change_j, work, 1e-10);
        close(
            inventory(&released.motion).2 - inventory(&loaded.motion).2,
            work,
            1e-10,
        );
        assert_eq!(released.normal_contact_work_j, 0.0);
        assert_eq!(released.absolute_contact_work_j, 0.0);
        assert_eq!(released.energy.inertia.numerical_loss_j, 0.0);
        assert_eq!(released.energy.fracture_dissipation_j, 0.0);
        close(released.energy.total_balance_error_j, 0.0, 1e-10);
        assert_eq!(released.motion.time_s(), motion.time_s() + 2.0 * h);
        assert_eq!(released.motion.accepted_steps(), 2);
        assert_eq!(released.motion.remaining_steps(), 1);
        assert_eq!(
            released.motion.solid().released_faces(),
            motion.solid().released_faces()
        );
        assert_eq!(
            released.motion.solid().fracture_dissipation_j(),
            motion.solid().fracture_dissipation_j()
        );
    });
}

#[test]
fn moving_support_supplies_native_reaction_work_during_sustained_contact() {
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
        let common_speed = -1.0 / 32.0;
        let velocities = vec![[DRIFT[0], DRIFT[1], common_speed]; target.nodes_m().len()];
        let motion =
            BrittleTetMotion::try_from_solid(&split.solid, &velocities, 2.0, 3, cx).unwrap();
        let h = 1.0 / 64.0;
        let acceleration = 1.0;
        let mut supports = Vec::new();
        for &node in &target.tetrahedra()[0] {
            for axis in 0..3 {
                supports.push((
                    3 * node + axis,
                    motion.solid().solution().displacement_m[node][axis]
                        + h * motion.velocity_mps()[node][axis],
                ));
            }
        }
        supports.sort_unstable_by_key(|&(dof, _)| dof);
        let result = motion
            .advance_released_face_sustained(
                h,
                &supports,
                &body_forces(motion.solid(), acceleration),
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
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
                1e-12,
            );
        }
        for &node in &target.tetrahedra()[1] {
            let expected = if target.nodes_m()[node][2] == 0.0 {
                common_speed
            } else {
                common_speed - 2.5 * acceleration * h
            };
            close(result.motion.velocity_mps()[node][2], expected, 1e-8);
        }
        for force in &result.solve.forces_n {
            close(force[0], acceleration / 8.0, 1e-8);
        }
        // A's external force is +a and the three normal forces sum to
        // -3a/8, hence its prescribed motion needs reaction -5a/8.
        let support_work = -0.625 * acceleration * common_speed * h;
        let dead_work = 0.3125 * acceleration * acceleration * h * h;
        close(result.energy.support_mean_work_j, support_work, 1e-10);
        close(result.energy.dead_load_mean_work_j, dead_work, 1e-10);
        close(
            result.energy.inertia.kinetic_change_j,
            support_work + dead_work,
            1e-10,
        );
        close(
            inventory(&result.motion).2 - inventory(&motion).2,
            support_work + dead_work,
            1e-10,
        );
        close(result.normal_contact_work_j, 0.0, 1e-10);
        close(result.energy.total_balance_error_j, 0.0, 1e-10);
        assert_eq!(result.energy.fracture_dissipation_j, 0.0);
        assert_eq!(result.energy.inertia.numerical_loss_j, 0.0);
        assert!(minimum_gap(&result.motion) >= 0.0);
        assert_eq!(result.motion.time_s(), motion.time_s() + h);
        assert_eq!(result.motion.accepted_steps(), 1);
        assert_eq!(result.motion.remaining_steps(), 2);
    });
}

#[test]
fn sustained_response_budget_and_impact_mode_refusals_preserve_the_original_for_retry() {
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
        let motion = open_motion(&split.solid, 0.0, 0.0, 3, cx);
        let h = 1.0 / 64.0;
        let forces = body_forces(motion.solid(), 1.0);
        let original_u = motion.solid().solution().displacement_m.clone();
        let original_v = motion.velocity_mps().to_vec();
        let original_history = motion.solid().solution().element_states.clone();
        let mut limited = sustained_controls();
        limited.solve.maximum_response_evaluations = 1;
        let error = motion
            .advance_released_face_sustained(h, &[], &forces, settings(), limited, cx)
            .unwrap_err();
        let BrittleFaceSustainedError::Solve(error) = error else {
            panic!("expected the original contact response budget");
        };
        assert!(matches!(
            *error,
            SolveError::ResponseBudget { evaluations: 1, .. }
        ));

        let approaching = open_motion(&split.solid, 0.0, SPEED, 3, cx);
        assert!(matches!(
            approaching.advance_released_face_sustained(
                h,
                &[],
                &forces,
                settings(),
                sustained_controls(),
                cx
            ),
            Err(BrittleFaceSustainedError::ImpactRequired { .. })
        ));
        let penetrated = open_motion(&split.solid, -2.0_f64.powi(-45), 0.0, 3, cx);
        assert!(matches!(
            penetrated.advance_released_face_sustained(
                h,
                &[],
                &forces,
                settings(),
                sustained_controls(),
                cx
            ),
            Err(BrittleFaceSustainedError::Geometry { .. })
        ));
        let open = open_motion(&split.solid, 1.0 / 2048.0, 0.0, 3, cx);
        // A positive force held while this finite initial gap closes would
        // erase physical kinetic energy. It needs the existing impact route.
        assert!(matches!(
            open.advance_released_face_sustained(
                1.0 / 16.0,
                &[],
                &forces,
                settings(),
                sustained_controls(),
                cx
            ),
            Err(BrittleFaceSustainedError::Mode { .. })
        ));
        with_cx(true, |cancelled| {
            assert!(
                motion
                    .advance_released_face_sustained(
                        h,
                        &[],
                        &forces,
                        settings(),
                        sustained_controls(),
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
        assert_eq!(motion.remaining_steps(), 3);
        let first = motion
            .advance_released_face_sustained(h, &[], &forces, settings(), sustained_controls(), cx)
            .unwrap();
        let second = motion
            .advance_released_face_sustained(h, &[], &forces, settings(), sustained_controls(), cx)
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
        assert_eq!(first.solve, second.solve);
        assert_eq!(first.energy, second.energy);
        assert_eq!(first.native_energy, second.native_energy);
    });
}
