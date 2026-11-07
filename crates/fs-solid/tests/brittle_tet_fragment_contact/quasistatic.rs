//! Failed native volumes transmit compression across the released interface.
//! No independent fragment stiffness inverse exists in the tested normal mode.

use super::*;
use fs_solid::brittle::solid::fracture::contact::{
    BrittleFaceEquilibriumControls, BrittleFaceEquilibriumError,
};
use fs_solver::NewtonKrylovConfig;

#[path = "quasistatic/approach.rs"]
mod approach;

fn controls() -> BrittleFaceEquilibriumControls {
    BrittleFaceEquilibriumControls {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 0.0,
            forcing_maximum: 1e-8,
            ..NewtonKrylovConfig::default()
        },
        maximum_iterations: 24,
        maximum_operator_evaluations: 100,
        maximum_contact_rows: 3,
        maximum_workspace_entries: 65_536,
        displacement_scale_m: 1.0 / 128.0,
        force_scale_n: 1.0,
        force_tolerance_n: 1e-10,
        activation_tolerance_m: 1e-11,
        tangential_tolerance_m: 1e-11,
        complementarity_tolerance_j: 1e-11,
        contact_work_tolerance_j: 1e-11,
        force_moment_tolerance_n_m: 1e-11,
        work: work_controls(),
    }
}

fn supported_load(solid: &BrittleTetSolid<'_>) -> (Vec<(usize, f64)>, Vec<(usize, f64)>) {
    let face = solid.released_faces()[0];
    // The lower fragment is held; only the upper face's three normal DOFs
    // are free. Both native materials have already lost all tensile traction.
    let free: Vec<_> = face.paired_nodes[1]
        .iter()
        .map(|node| 3 * node + 2)
        .collect();
    let held = prescribed(&solid.solution().displacement_m)
        .into_iter()
        .filter(|(dof, _)| !free.contains(dof))
        .collect();
    let force = free
        .into_iter()
        .enumerate()
        .map(|(i, dof)| (dof, -((i + 1) as f64) / 8.0))
        .collect();
    (held, force)
}

#[test]
fn mixed_recontact_transmits_pressure_without_inverting_softened_fragments() {
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
        let original = split.solid.solution().clone();
        let history_loss = split.solid.fracture_dissipation_j();
        let (held, force) = supported_load(&split.solid);
        for (length_scale, force_scale, translation) in [
            (1.0 / 128.0, 1.0, 0.0),
            (1.0 / 64.0, 2.0, 0.0),
            (1.0 / 128.0, 1.0, 1.0 / 512.0),
        ] {
            let held: Vec<_> = held
                .iter()
                .map(|&(dof, value)| (dof, value + if dof % 3 == 2 { translation } else { 0.0 }))
                .collect();
            let mut numerical = controls();
            numerical.displacement_scale_m = length_scale;
            numerical.force_scale_n = force_scale;
            let step = split
                .solid
                .advance_released_face_equilibrium(&held, &force, &[0.5; 3], numerical, cx)
                .unwrap();
            assert_eq!(step.rows.len(), 3);
            assert!(step.nonlinear.converged);
            assert!(step.maximum_force_residual_n <= numerical.force_tolerance_n);
            assert!(step.operator_evaluations <= numerical.maximum_operator_evaluations);
            for (i, row) in step.rows.iter().enumerate() {
                close(row.normal_force_n, (i + 1) as f64 / 8.0, 1e-10);
                close(row.gap_m, 0.0, 1e-11);
                let lower = 3 * row.nodes[0] + 2;
                let index = held.iter().position(|&(dof, _)| dof == lower).unwrap();
                close(
                    step.native.solid.solution().increments[0].reactions_n[index],
                    row.normal_force_n,
                    1e-10,
                );
            }
            close(
                step.native.solid.fracture_dissipation_j(),
                history_loss,
                1e-12,
            );
            close(step.native.energy.fracture_dissipation_j, 0.0, 1e-12);
            close(step.native.energy.stored_after_j, 0.0, 1e-12);
            close(step.normal_contact_work_j, 0.0, 1e-11);
            close(step.absolute_contact_couple_n_m, 0.0, 1e-11);
            close(
                step.external_energy.dead_load_endpoint_work_j,
                -0.75 * translation,
                1e-10,
            );
            close(
                step.external_energy.support_endpoint_work_j,
                0.75 * translation,
                1e-10,
            );
            close(
                step.external_energy.external_endpoint_work_defect_j(),
                0.0,
                1e-10,
            );
            for (old, new) in original
                .displacement_m
                .iter()
                .zip(&step.native.solid.solution().displacement_m)
            {
                for axis in 0..3 {
                    close(
                        new[axis],
                        old[axis] + if axis == 2 { translation } else { 0.0 },
                        1e-10,
                    );
                }
            }
            assert_eq!(split.solid.solution(), &original);
        }
    });
}

#[test]
fn recontact_keeps_native_compression_and_exact_fracture_history() {
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
        let (mut held, force) = supported_load(&split.solid);
        let apex = target.tetrahedra()[1]
            .iter()
            .copied()
            .find(|&node| target.nodes_m()[node][2] == 1.0)
            .unwrap();
        held.iter_mut()
            .find(|(dof, _)| *dof == 3 * apex + 2)
            .unwrap()
            .1 = -1.0 / 1024.0;
        let step = split
            .solid
            .advance_released_face_equilibrium(&held, &force, &[0.5; 3], controls(), cx)
            .unwrap();
        // E=1024 Pa, epsilon_zz=-1/1024, nu=0: sigma_zz=-1 Pa.
        // On this unit simplex grad(N_origin)_z=-1 and V=1/6; the
        // original volume therefore adds +1/6 N at the origin's face node.
        close(
            step.native.solid.solution().element_stress[1][2],
            -1.0,
            1e-10,
        );
        for (index, row) in step.rows.iter().enumerate() {
            let extra = if target.nodes_m()[row.nodes[1]][0..2] == [0.0, 0.0] {
                1.0 / 6.0
            } else {
                0.0
            };
            close(row.normal_force_n, (index + 1) as f64 / 8.0 + extra, 1e-10);
            close(row.gap_m, 0.0, 1e-11);
        }
        close(step.native.solid.stored_energy_j(), 1.0 / 12_288.0, 1e-11);
        close(
            step.native.solid.fracture_dissipation_j(),
            split.solid.fracture_dissipation_j(),
            1e-12,
        );
        assert!(
            step.native
                .solid
                .solution()
                .element_states
                .iter()
                .all(|state| state.failed())
        );
        close(step.absolute_contact_work_j, 0.0, 1e-11);
        close(step.external_energy.equilibrium_work_residual_j, 0.0, 1e-10);
    });
}

#[test]
fn shared_crack_edges_report_aggregate_force_and_keep_reversed_face_witnesses() {
    with_cx(false, |cx| {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, -1.0],
            [0.5, 0.5, 1.0],
        ];
        let tets = [[0, 2, 1, 4], [0, 1, 2, 5], [0, 3, 2, 4], [0, 2, 3, 5]];
        let geometry = TetLoadGeometry3::new(&nodes, &tets, budget(), 16, cx).unwrap();
        let materials = sources(&geometry);
        let failed = failed(&geometry, &materials, cx);
        let cut = topology(
            &geometry,
            &[
                InteriorSeam {
                    id: 951,
                    nodes: [0, 1, 2],
                    side_a_element: 0,
                },
                // Reverse the side order: normal and paired nodes reverse together,
                // preserving the full positive contact covector on the shared edge.
                InteriorSeam {
                    id: 952,
                    nodes: [0, 2, 3],
                    side_a_element: 3,
                },
            ],
            cx,
        );
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 16, cx).unwrap();
        let split = failed
            .split_fracture(&cut, &target, fracture_controls(), cx)
            .unwrap();
        let mut free: Vec<_> = [1, 3]
            .into_iter()
            .flat_map(|element| target.tetrahedra()[element])
            .filter(|&node| target.nodes_m()[node][2] == 0.0)
            .map(|node| 3 * node + 2)
            .collect();
        free.sort_unstable();
        free.dedup();
        assert_eq!(free.len(), 4);
        let held: Vec<_> = prescribed(&split.solid.solution().displacement_m)
            .into_iter()
            .filter(|(dof, _)| !free.contains(dof))
            .collect();
        let force: Vec<_> = free.iter().map(|&dof| (dof, -1.0 / 8.0)).collect();
        let mut c = controls();
        c.maximum_contact_rows = 6;
        let step = split
            .solid
            .advance_released_face_equilibrium(&held, &force, &[0.5; 6], c, cx)
            .unwrap();
        assert_eq!(step.rows.len(), 4);
        assert_eq!(step.original_pairs.len(), 6);
        let mut counts = vec![0; 4];
        for original in &step.original_pairs {
            counts[original.aggregate_row] += 1;
            let row = step.rows[original.aggregate_row];
            assert!(
                (original.nodes == row.nodes && original.normal == row.normal)
                    || (original.nodes == [row.nodes[1], row.nodes[0]]
                        && original.normal == row.normal.map(|n| -n))
            );
        }
        counts.sort_unstable();
        assert_eq!(counts, [1, 1, 2, 2]);
        for row in &step.rows {
            close(row.normal_force_n, 1.0 / 8.0, 1e-10);
        }
        close(step.rows.iter().map(|r| r.normal_force_n).sum(), 0.5, 1e-10);
        close(
            step.native.solid.fracture_dissipation_j(),
            split.solid.fracture_dissipation_j(),
            1e-12,
        );
        c.maximum_contact_rows = 5;
        assert!(matches!(
            split
                .solid
                .advance_released_face_equilibrium(&held, &force, &[0.5; 6], c, cx),
            Err(BrittleFaceEquilibriumError::Budget("released contact rows"))
        ));
    });
}

#[test]
fn loaded_open_gap_is_refused_and_retry_retains_fracture_history() {
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
        let mut opening = split.solid.solution().displacement_m.clone();
        for &node in &target.tetrahedra()[1] {
            opening[node][2] += 1.0 / 1024.0;
        }
        let open = split
            .solid
            .advance(&prescribed(&opening), &[], settings(), work_controls(), cx)
            .unwrap()
            .solid;
        let before = open.solution().clone();
        let (held, force) = supported_load(&open);
        // An open unloaded interface carries exactly zero pressure. The native
        // accepted equilibrium needs no stiffness inverse or state projection.
        let unloaded = open
            .advance_released_face_equilibrium(&held, &[], &[0.0; 3], controls(), cx)
            .unwrap();
        for row in &unloaded.rows {
            assert_eq!(row.normal_force_n, 0.0);
            close(row.gap_m, 1.0 / 1024.0, 1e-12);
        }
        assert_eq!(
            unloaded.native.solid.solution().element_states,
            before.element_states
        );
        assert_eq!(
            unloaded.native.solid.solution().displacement_m,
            before.displacement_m
        );
        close(unloaded.native.energy.material_work_j, 0.0, 1e-12);
        assert!(matches!(
            open.advance_released_face_equilibrium(&held, &force, &[0.5; 3], controls(), cx),
            Err(BrittleFaceEquilibriumError::Geometry {
                reason: "loaded contact must touch at both endpoints",
                ..
            })
        ));
        assert_eq!(open.solution(), &before);

        let (held, force) = supported_load(&split.solid);
        let mut bounded = controls();
        bounded.maximum_contact_rows = 2;
        assert!(matches!(
            split
                .solid
                .advance_released_face_equilibrium(&held, &force, &[0.5; 3], bounded, cx),
            Err(BrittleFaceEquilibriumError::Budget("released contact rows"))
        ));
        let mut bounded = controls();
        bounded.newton.max_linear_cycles = usize::MAX;
        assert!(matches!(
            split
                .solid
                .advance_released_face_equilibrium(&held, &force, &[0.5; 3], bounded, cx),
            Err(BrittleFaceEquilibriumError::Budget(_))
        ));
        let mut bounded = controls();
        bounded.maximum_operator_evaluations = 1;
        assert!(matches!(
            split
                .solid
                .advance_released_face_equilibrium(&held, &force, &[0.5; 3], bounded, cx),
            Err(BrittleFaceEquilibriumError::Budget(
                "native force/tangent assemblies"
            ))
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                split.solid.advance_released_face_equilibrium(
                    &held,
                    &force,
                    &[0.5; 3],
                    controls(),
                    cancelled
                ),
                Err(BrittleFaceEquilibriumError::Solid(
                    fs_solid::brittle::solid::BrittleTetError::Cancelled
                ))
            ));
        });
        let retry = split
            .solid
            .advance_released_face_equilibrium(&held, &force, &[0.5; 3], controls(), cx)
            .unwrap();
        close(
            retry.native.solid.fracture_dissipation_j(),
            split.solid.fracture_dissipation_j(),
            1e-12,
        );
    });
}
