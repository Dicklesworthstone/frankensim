//! G1: exact closed-crack constraints with retained tetrahedral compression.

use super::*;
use fs_solid::brittle::solid::dynamics::modes::FrozenBrittleModesControls;
use fs_solid::brittle::solid::dynamics::modes::contact::FrozenBrittleContactControls;

fn modal_controls() -> FrozenBrittleContactControls {
    FrozenBrittleContactControls {
        assembly: FrozenBrittleModesControls {
            maximum_dofs: 24,
            maximum_staged_entries: 288,
        },
        maximum_free_dofs: 24,
        maximum_contact_rows: 3,
        maximum_dense_entries: 32_768,
        gap_tolerance_m: 1e-10,
        tangential_tolerance_m: 1e-10,
        velocity_tolerance_mps: 1e-8,
        minimum_normal_preload_n: 1e-8,
        rank_tolerance: 1e-12,
        constraint_tolerance: 1e-11,
        momentum_tolerance: 1e-9,
        spectrum_tolerance: 1e-9,
    }
}

fn with_loaded_fragments(run: impl FnOnce(&BrittleTetMotion<'_>, [[f64; 3]; 1], &Cx<'_>)) {
    with_cx(false, |cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let materials = sources(&geometry);
        let failed = failed(&geometry, &materials, cx);
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
        let split = failed
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let unloaded = split
            .solid
            .advance(
                &prescribed(&vec![[0.0; 3]; target.nodes_m().len()]),
                &[],
                settings(),
                work_controls(),
                cx,
            )
            .unwrap()
            .solid;
        let motion = open_motion(&unloaded, 0.0, 0.0, 3, cx);
        let loaded = motion
            .advance_released_face_sustained(
                1.0 / 64.0,
                &[],
                &body_forces(&unloaded, 1.0),
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        let mut preload = [[0.0; 3]];
        for (row, force) in loaded.rows.iter().zip(&loaded.solve.forces_n) {
            assert_eq!(row.face_id, loaded.motion.solid().released_faces()[0].id);
            preload[0][row.vertex] = force[0];
        }
        run(&loaded.motion, preload, cx);
    });
}

fn hold_tangent_dofs(motion: &BrittleTetMotion<'_>) -> Vec<usize> {
    (0..3 * motion.solid().geometry().nodes_m().len())
        .filter(|dof| dof % 3 != 2)
        .collect()
}

#[test]
fn g1_loaded_crack_modes_match_independent_consistent_mass_spectrum() {
    with_loaded_fragments(|motion, preload, cx| {
        let before = motion.solid().solution().clone();
        let fixed = hold_tangent_dofs(motion);
        let result = motion
            .frozen_contact_modes(&preload, &fixed, modal_controls(), cx)
            .unwrap();
        let pencil = result.pencil();
        assert_eq!(pencil.free_dofs().len(), 8);
        assert_eq!(result.constraints().len(), 3);
        assert_eq!(result.dimension(), 5);
        assert_eq!(result.modes().len(), 5);
        assert_eq!(pencil.element_damage(), &[1.0, 1.0]);
        close(pencil.reference_mass_kg(), 2.0, 1e-15);

        // Each unit tet has Kzz=(E/6)[1,-1;-1,1] between its
        // origin and apex, Mzz=(I+ones)/20. Identifying the three
        // paired face z-coordinates leaves three mechanisms and two
        // independent elastic eigenvalues: 5E/3 and 20E/3.
        let expected = [0.0, 0.0, 0.0, 5.0 * 1024.0 / 3.0, 20.0 * 1024.0 / 3.0];
        for (mode, lambda) in result.modes().iter().zip(expected) {
            close(mode.coordinates.lambda, lambda, 2e-8);
            assert_eq!(mode.displacement_shape.len(), 24);
            assert_eq!(mode.constraint_force_shape.len(), 3);
            for &dof in &fixed {
                assert_eq!(mode.displacement_shape[dof], 0.0);
            }
            let free = pencil.free_dofs();
            let mut mass_norm = 0.0;
            for (i, &dof) in free.iter().enumerate() {
                let mut stiffness_action = 0.0;
                let mut mass_action = 0.0;
                for (j, &other) in free.iter().enumerate() {
                    stiffness_action +=
                        pencil.stiffness_n_m().get(i, j) * mode.displacement_shape[other];
                    mass_action += pencil.mass_kg().get(i, j) * mode.displacement_shape[other];
                }
                mass_norm += mode.displacement_shape[dof] * mass_action;
                let force: f64 = result
                    .constraints()
                    .iter()
                    .zip(&mode.constraint_force_shape)
                    .map(|(constraint, value)| {
                        let sign = if dof / 3 == constraint.row.nodes[1] {
                            1.0
                        } else if dof / 3 == constraint.row.nodes[0] {
                            -1.0
                        } else {
                            0.0
                        };
                        sign * constraint.direction[dof % 3] * value
                    })
                    .sum();
                close(
                    stiffness_action - mode.coordinates.lambda * mass_action,
                    force,
                    2e-8,
                );
            }
            close(mass_norm, 1.0, 1e-11);
            for constraint in result.constraints() {
                let [a, b] = constraint.row.nodes;
                let normal_jump: f64 = (0..3)
                    .map(|axis| {
                        constraint.direction[axis]
                            * (mode.displacement_shape[3 * b + axis]
                                - mode.displacement_shape[3 * a + axis])
                    })
                    .sum();
                close(normal_jump, 0.0, 1e-11);
            }
            assert!(mode.constraint_relative_residual <= modal_controls().constraint_tolerance);
            assert!(mode.momentum_relative_residual <= modal_controls().momentum_tolerance);
        }
        // Preload selects a unilateral branch; it is not a penalty stiffness.
        let stronger = [preload[0].map(|force| 8.0 * force)];
        let repeated = motion
            .frozen_contact_modes(&stronger, &fixed, modal_controls(), cx)
            .unwrap();
        assert_eq!(
            result.projected_stiffness_n_m(),
            repeated.projected_stiffness_n_m()
        );
        assert_eq!(result.projected_mass_kg(), repeated.projected_mass_kg());
        assert_eq!(motion.solid().solution(), &before);
        assert_eq!(motion.accepted_steps(), 1);
    });
}

#[test]
fn g0_contact_modes_refuse_unloaded_touching_motion_and_excess_workspaces() {
    with_loaded_fragments(|motion, preload, cx| {
        let fixed = hold_tangent_dofs(motion);
        for invalid in [[[0.0; 3]], [[-1.0; 3]], [[f64::NAN; 3]]] {
            assert!(
                motion
                    .frozen_contact_modes(&invalid, &fixed, modal_controls(), cx)
                    .is_err()
            );
        }
        assert!(
            motion
                .frozen_contact_modes(&[], &fixed, modal_controls(), cx)
                .is_err()
        );
        assert!(
            motion
                .frozen_contact_modes(&preload, &[0, 0], modal_controls(), cx)
                .is_err()
        );
        for controls in [
            FrozenBrittleContactControls {
                maximum_free_dofs: 7,
                ..modal_controls()
            },
            FrozenBrittleContactControls {
                maximum_contact_rows: 2,
                ..modal_controls()
            },
            FrozenBrittleContactControls {
                maximum_dense_entries: 1,
                ..modal_controls()
            },
        ] {
            assert!(
                motion
                    .frozen_contact_modes(&preload, &fixed, controls, cx)
                    .is_err()
            );
        }
        let mut velocity = motion.velocity_mps().to_vec();
        let face = &motion.solid().released_faces()[0];
        velocity[face.paired_nodes[1][0]][2] += 1.0;
        let moving =
            BrittleTetMotion::try_from_solid(motion.solid(), &velocity, motion.time_s(), 1, cx)
                .unwrap();
        assert!(
            moving
                .frozen_contact_modes(&preload, &fixed, modal_controls(), cx)
                .is_err()
        );
        with_cx(true, |cancelled| {
            assert!(
                motion
                    .frozen_contact_modes(&preload, &fixed, modal_controls(), cancelled)
                    .is_err()
            );
        });
    });
}
