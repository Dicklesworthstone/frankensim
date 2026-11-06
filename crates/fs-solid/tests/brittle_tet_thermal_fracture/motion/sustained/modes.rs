//! G1: closed thermal cracks use mechanical strain and original reference inertia.

use super::*;
use fs_solid::brittle::solid::dynamics::fracture::contact::friction::{
    BrittleFaceFriction, EllipticFriction,
};
use fs_solid::brittle::solid::expansion::dynamics::modes::{
    FrozenBrittleContactControls, FrozenBrittleModesControls, FrozenBrittleStickingControls,
};

fn controls() -> FrozenBrittleContactControls {
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

#[test]
fn g1_thermal_contact_spectrum_retains_original_mass_and_compression_branch() {
    for compression in [false, true] {
        with_closed_thermal_motion(compression, 6, |motion, _warm, cx| {
            let loaded = motion
                .advance_released_face_sustained(
                    INTERVAL,
                    &[],
                    &body_forces(motion, 1.0),
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
            let before = loaded.motion.solid().solution().clone();
            let fixed: Vec<_> = (0..24).filter(|dof| dof % 3 != 2).collect();
            let result = loaded
                .motion
                .frozen_contact_modes(&preload, &fixed, controls(), cx)
                .unwrap();
            assert_eq!(result.dimension(), 5);
            assert_eq!(result.constraints().len(), 3);
            assert_eq!(result.pencil().element_damage(), &[1.0, 1.0]);
            close(result.pencil().reference_mass_kg(), 2.0, 1e-15);
            // The source now reports rho=3, but the original two tetrahedra
            // each retain one kg at rho=6. With Kzz=E/6 and consistent mass
            // (I+ones)/20, three tied face coordinates give 5E/3 and 20E/3.
            // Normal opening at complete failure has zero Kzz even though
            // the actual touching faces carry a positive normal preload.
            let young = if compression { 1200.0 } else { 0.0 };
            let expected = [0.0, 0.0, 0.0, 5.0 * young / 3.0, 20.0 * young / 3.0];
            for (mode, lambda) in result.modes().iter().zip(expected) {
                close(mode.coordinates.lambda, lambda, 2e-8);
                for &dof in &fixed {
                    assert_eq!(mode.displacement_shape[dof], 0.0);
                }
                let free = result.pencil().free_dofs();
                let mass_norm: f64 = free
                    .iter()
                    .enumerate()
                    .map(|(i, &a)| {
                        free.iter()
                            .enumerate()
                            .map(|(j, &b)| {
                                mode.displacement_shape[a]
                                    * result.pencil().mass_kg().get(i, j)
                                    * mode.displacement_shape[b]
                            })
                            .sum::<f64>()
                    })
                    .sum();
                close(mass_norm, 1.0, 1e-11);
                for constraint in result.constraints() {
                    let [a, b] = constraint.row.nodes;
                    let jump: f64 = (0..3)
                        .map(|axis| {
                            constraint.direction[axis]
                                * (mode.displacement_shape[3 * b + axis]
                                    - mode.displacement_shape[3 * a + axis])
                        })
                        .sum();
                    close(jump, 0.0, 1e-11);
                }
            }
            assert_eq!(loaded.motion.solid().solution(), &before);
            assert_eq!(loaded.motion.accepted_steps(), 4);
            assert_eq!(loaded.motion.reference_density_kg_m3(), &[DENSITY; 2]);
            assert!(
                loaded
                    .motion
                    .frozen_contact_modes(&[[0.0; 3]], &fixed, controls(), cx)
                    .is_err()
            );

            // The same thermally fractured body also supports a strictly
            // interior static ellipse, with physical constraints on all axes.
            let friction = [BrittleFaceFriction {
                tangent_direction: [3.0, 4.0, 0.0],
                vertices: [EllipticFriction {
                    coefficients_mu: [0.7, 0.25],
                    principal_axis: [4.0, -3.0],
                }; 3],
            }];
            let forces = [preload[0].map(|force| [force, 0.0, 0.0])];
            let sticking = loaded
                .motion
                .frozen_sticking_contact_modes(
                    &friction,
                    &forces,
                    &[],
                    FrozenBrittleStickingControls {
                        contact: FrozenBrittleContactControls {
                            maximum_contact_rows: 9,
                            maximum_dense_entries: 65_536,
                            ..controls()
                        },
                        minimum_friction_margin: 1e-6,
                    },
                    cx,
                )
                .unwrap();
            assert_eq!(sticking.dimension(), 15);
            assert_eq!(sticking.constraints().len(), 9);
            close(sticking.pencil().reference_mass_kg(), 2.0, 1e-15);
            let face = &loaded.motion.solid().released_faces()[0];
            for mode in sticking.modes() {
                let mut mass_norm = 0.0;
                for i in 0..24 {
                    for j in 0..24 {
                        mass_norm += mode.displacement_shape[i]
                            * sticking.pencil().mass_kg().get(i, j)
                            * mode.displacement_shape[j];
                    }
                }
                close(mass_norm, 1.0, 1e-11);
                for vertex in 0..3 {
                    for axis in 0..3 {
                        close(
                            mode.displacement_shape[3 * face.paired_nodes[0][vertex] + axis],
                            mode.displacement_shape[3 * face.paired_nodes[1][vertex] + axis],
                            1e-11,
                        );
                    }
                }
            }
            assert_eq!(loaded.motion.solid().solution(), &before);
            assert_eq!(loaded.motion.accepted_steps(), 4);
        });
    }
}
