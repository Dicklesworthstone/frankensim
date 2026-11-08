//! G0/G3: original failed material and consistent inertia under two-sided contact.
use super::*;
use fs_contact::sustained::solve::SolveControls;
use fs_solid::hyper3::plastic::ductile::fracture::contact::{
    FiniteDuctileFaceContactControls, FiniteDuctileFaceContactError,
};
use fs_solver::NewtonKrylovConfig;

fn contact_controls() -> FiniteDuctileFaceContactControls {
    let mut c = FiniteDuctileFaceContactControls::default();
    c.contact.solve = SolveControls {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-8,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        difference_step_n: 0.001,
        force_residual_tolerance_n: 1e-7,
        ..SolveControls::default()
    };
    c.contact.maximum_normal_work_defect_j = 1e-8;
    c.contact.work_balance_tolerance_j = 1e-8;
    c.normal_scale_n_per_m = 1e4;
    c.activation_tolerance_m = 1e-10;
    c.velocity_tolerance_mps = 1e-8;
    c.force_moment_tolerance_n_m = 1e-8;
    c
}

fn compressive_loads(
    face: fs_solid::hyper3::plastic::ductile::fracture::ReleasedDuctileFace,
    forces: [f64; 3],
) -> Vec<(usize, f64)> {
    let mut loads = Vec::new();
    for (vertex, force) in forces.into_iter().enumerate() {
        for axis in 0..3 {
            let value = force * face.normal[axis];
            if value != 0.0 {
                // External compression is opposed by the internal contact row.
                loads.push((3 * face.paired_nodes[0][vertex] + axis, value));
                loads.push((3 * face.paired_nodes[1][vertex] + axis, -value));
            }
        }
    }
    loads
}

#[test]
fn released_faces_carry_equal_opposite_normal_loads_without_resetting_failed_material() {
    with_cx(|gate, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let original_bands = bands(&geometry);
        let material = model(original_bands[0].band_width_m(), true);
        let materials = [&material, &material];
        let solid = failed_solid(&geometry, &materials, &original_bands, cx);
        let native =
            FiniteDuctileTetMotion::try_new(solid, &[12.0, 12.0], &[[0.0; 3]; 5], cx).unwrap();
        let topology = topology(&geometry, cx);
        let target = TetLoadGeometry3::new(
            topology.nodes_m(),
            topology.tetrahedra(),
            TetAssemblyBudget::standard(),
            8,
            cx,
        )
        .unwrap();
        let target_bands = bands(&target);
        let cut = native
            .split_fracture(&topology, &target, &target_bands, release_controls(), cx)
            .unwrap();
        let before = cut.motion.motion().solid().element_states().to_vec();
        let face = cut.motion.released_faces()[0];
        near(face.normal[0], 1.0);
        // These are exact closed-row solutions: external and contact nodal
        // forces cancel componentwise while the fully failed bulk stays native.
        let forces = [0.125, 0.25, 0.5];
        let loads = compressive_loads(face, forces);
        let held = cut
            .motion
            .advance_with_face_contact(
                &[],
                &loads,
                [0.0; 3],
                conditions(),
                Some(&forces),
                IncrementSettings::default(),
                motion_controls(),
                contact_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(held.contact.rows.len(), 3);
        for (i, row) in held.contact.rows.iter().enumerate() {
            assert_eq!(row.face_id, face.id);
            assert_eq!(row.vertex, i);
            assert_eq!(
                row.nodes,
                [face.paired_nodes[0][i], face.paired_nodes[1][i]]
            );
            near(held.contact.solve.forces_n[i][0], forces[i]);
            assert_eq!(held.contact.solve.forces_n[i][1..], [0.0, 0.0]);
            assert_eq!(held.contact.solve.paths_m[i], [0.0; 3]);
        }
        assert_eq!(held.step.motion.motion().solid().element_states(), before);
        assert_eq!(cut.motion.motion().solid().element_states(), before);
        assert_eq!(
            held.step.motion.released_faces(),
            cut.motion.released_faces()
        );
        assert_eq!(
            held.step.motion.motion().reference_mass_kg(),
            cut.motion.motion().reference_mass_kg(),
        );
        assert_eq!(
            held.step.motion.motion().solid().displacement_m(),
            cut.motion.motion().solid().displacement_m(),
        );
        near(held.step.material_energy.plastic_dissipation_j, 0.0);
        near(held.step.material_energy.damage_release_j, 0.0);
        near(held.step.inertia.kinetic_change_j, 0.0);
        near(held.step.inertia.numerical_loss_j, 0.0);
        near(held.contact.normal_contact_work_j, 0.0);
        near(held.contact.normal_work_defect_j, 0.0);
        near(held.contact.absolute_contact_couple_n_m, 0.0);
        near(held.contact.load_decomposition_error_j, 0.0);
        near(held.contact.work_balance_error_j, 0.0);
        assert_eq!(held.contact.minimum_gap_m, 0.0);
        assert!(held.contact.path_windows > 0);

        gate.request();
        let refused = held.step.motion.advance_with_face_contact(
            &[],
            &loads,
            [0.0; 3],
            conditions(),
            Some(&forces),
            IncrementSettings::default(),
            motion_controls(),
            contact_controls(),
            cx,
        );
        assert!(refused.is_err());
        assert_eq!(held.step.motion.motion().solid().element_states(), before);
    });
}

#[test]
fn separated_native_fragments_keep_mass_kinetic_work_and_zero_contact_force() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let original_bands = bands(&geometry);
        let material = model(original_bands[0].band_width_m(), true);
        let materials = [&material, &material];
        let native = FiniteDuctileTetMotion::try_new(
            failed_solid(&geometry, &materials, &original_bands, cx),
            &[12.0, 12.0],
            &[[0.0; 3]; 5],
            cx,
        )
        .unwrap();
        let topology = topology(&geometry, cx);
        let target = TetLoadGeometry3::new(
            topology.nodes_m(),
            topology.tetrahedra(),
            TetAssemblyBudget::standard(),
            8,
            cx,
        )
        .unwrap();
        let target_bands = bands(&target);
        let cut = native
            .split_fracture(&topology, &target, &target_bands, release_controls(), cx)
            .unwrap();
        let normal = cut.motion.released_faces()[0].normal;
        let acceleration = 0.03125;
        let mut field = vec![[0.0; 3]; target.nodes_m().len()];
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            let sign = if element == 0 { -1.0 } else { 1.0 };
            for &node in tet {
                field[node] = normal.map(|n| sign * acceleration * n);
            }
        }
        let loads = cut
            .motion
            .motion()
            .reference_acceleration_load_n(&field, cx)
            .unwrap();
        let separated = cut
            .motion
            .advance_with_face_contact(
                &[],
                &loads,
                [0.0; 3],
                conditions(),
                None,
                IncrementSettings::default(),
                motion_controls(),
                contact_controls(),
                cx,
            )
            .unwrap();
        let dt = conditions().duration_s;
        let expected_gap = 2.0 * acceleration * dt * dt;
        let total_mass = cut.motion.motion().reference_mass_kg();
        let expected_kinetic = 0.5 * total_mass * (acceleration * dt).powi(2);
        for (forces, path) in separated
            .contact
            .solve
            .forces_n
            .iter()
            .zip(&separated.contact.solve.paths_m)
        {
            assert_eq!(*forces, [0.0; 3]);
            near(path[0], expected_gap);
        }
        near(separated.step.inertia.kinetic_change_j, expected_kinetic);
        near(separated.step.inertia.numerical_loss_j, expected_kinetic);
        near(separated.contact.dead_load_work_j, 2.0 * expected_kinetic);
        near(separated.contact.normal_contact_work_j, 0.0);
        near(separated.contact.absolute_normal_work_defect_j, 0.0);
        near(separated.contact.work_balance_error_j, 0.0);
        near(separated.step.material_energy.damage_release_j, 0.0);
        assert!(
            separated
                .step
                .motion
                .motion()
                .solid()
                .element_states()
                .iter()
                .all(|s| s.is_fully_failed())
        );
        for element in 0..2 {
            near(
                separated.step.motion.motion().solid().element_states()[element]
                    .dissipated_damage_j_m3(),
                cut.motion.motion().solid().element_states()[element].dissipated_damage_j_m3(),
            );
        }
        for (node, velocity) in separated
            .step
            .motion
            .motion()
            .velocity_mps()
            .iter()
            .enumerate()
        {
            for axis in 0..3 {
                near(velocity[axis], dt * field[node][axis]);
            }
        }
        assert!(separated.contact.path_windows > 0);
        assert!(separated.contact.minimum_gap_m >= 0.0);
    });
}

#[test]
fn loaded_sliding_faces_refuse_instead_of_applying_an_artificial_force_couple() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let original_bands = bands(&geometry);
        let material = model(original_bands[0].band_width_m(), true);
        let materials = [&material, &material];
        let native = FiniteDuctileTetMotion::try_new(
            failed_solid(&geometry, &materials, &original_bands, cx),
            &[12.0, 12.0],
            &[[0.0; 3]; 5],
            cx,
        )
        .unwrap();
        let topology = topology(&geometry, cx);
        let target = TetLoadGeometry3::new(
            topology.nodes_m(),
            topology.tetrahedra(),
            TetAssemblyBudget::standard(),
            8,
            cx,
        )
        .unwrap();
        let target_bands = bands(&target);
        let cut = native
            .split_fracture(&topology, &target, &target_bands, release_controls(), cx)
            .unwrap();
        let mut displaced = cut.motion.motion().solid().displacement_m().to_vec();
        for &node in &target.tetrahedra()[1] {
            displaced[node][1] += 0.001;
        }
        let prescribed: Vec<_> = displaced
            .as_flattened()
            .iter()
            .copied()
            .enumerate()
            .collect();
        let sliding = cut
            .motion
            .advance(
                &prescribed,
                &[],
                [0.0; 3],
                conditions(),
                IncrementSettings::default(),
                motion_controls(),
                cx,
            )
            .unwrap();
        let old_u = sliding.motion.motion().solid().displacement_m().to_vec();
        let old_states = sliding.motion.motion().solid().element_states().to_vec();
        let face = sliding.motion.released_faces()[0];
        let forces = [0.125, 0.25, 0.5];
        let loads = compressive_loads(face, forces);
        let refused = sliding.motion.advance_with_face_contact(
            &[],
            &loads,
            [0.0; 3],
            conditions(),
            Some(&forces),
            IncrementSettings::default(),
            motion_controls(),
            contact_controls(),
            cx,
        );
        assert!(matches!(
            refused,
            Err(FiniteDuctileFaceContactError::Face {
                reason: "loaded pairs must remain matched, active and normally at rest",
                ..
            })
        ));
        assert_eq!(sliding.motion.motion().solid().displacement_m(), old_u);
        assert_eq!(sliding.motion.motion().solid().element_states(), old_states);
    });
}
