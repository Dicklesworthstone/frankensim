//! G1/G3 retained numerical gauge under material softening and native motion.

use super::*;
use fs_solid::plastic_shell::drilling::PlasticShellReferenceDrilling;
use fs_solid::plastic_shell::dynamics::{PlasticShellDensity, PlasticShellReferenceInertia};

fn triangle() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

#[test]
fn g1_reference_gauge_keeps_its_stiffness_and_work_when_the_material_softens() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let cold = section(IsotropicElastic::new(E, NU, 0.05).unwrap());
        let hot = section(IsotropicElastic::new(0.25 * E, NU, 0.05).unwrap());
        let zero: Vec<_> = (0..18).map(|d| (d, 0.0)).collect();
        let initial = problem(&mesh, &cold, &zero, &[], &[1.0], 0.37)
            .solve(cx)
            .unwrap();
        let fixed: Vec<_> = zero.iter().copied().filter(|(d, _)| *d != 5).collect();
        let weight = DRILLING_ALPHA * (E * H / (2.0 * (1.0 + NU))) * (0.5 / 3.0);
        let spin = 0.01;
        let loads = [(5, weight * spin)];
        let reference = {
            // This constructor retains only the mesh, not the temporary field
            // of section references used during a thermal trial.
            let fields = [&cold];
            let mut p = problem(&mesh, &cold, &fixed, &loads, &[1.0], 0.37);
            p.sections = PlasticShellSections::PerElementBorrowed(&fields);
            PlasticShellReferenceDrilling::prepare(&mesh, &p, cx).unwrap()
        };
        let first = problem(&mesh, &cold, &fixed, &loads, &[1.0], 0.37)
            .solve_from_with_reference_drilling(
                &initial.coordinates,
                &initial.section_states,
                &reference,
                cx,
            )
            .unwrap();
        near(first.solution.coordinates[0][5], spin, 1e-12);
        let expected = 0.5 * weight * spin * spin;
        near(first.drilling.final_j, expected, 1e-15);
        near(first.drilling.numerical_loss_j, expected, 1e-15);
        near(first.drilling.work_j, 2.0 * expected, 1e-15);
        assert_eq!(first.solution.section_states, initial.section_states);
        let fields = [&hot];
        let mut hot_problem = problem(&mesh, &hot, &fixed, &loads, &[1.0], 0.37);
        hot_problem.sections = PlasticShellSections::PerElementBorrowed(&fields);
        let hold = hot_problem
            .solve_from_with_reference_drilling(
                &first.solution.coordinates,
                &first.solution.section_states,
                &reference,
                cx,
            )
            .unwrap();
        near(hold.solution.coordinates[0][5], spin, 1e-12);
        near(hold.drilling.initial_j, expected, 1e-15);
        near(hold.drilling.final_j, expected, 1e-15);
        near(hold.drilling.work_j, 0.0, 1e-15);
        near(hold.drilling.numerical_loss_j, 0.0, 1e-15);
        assert_eq!(hold.solution.section_states, initial.section_states);

        // An explicitly fresh gauge uses the softened modulus. This observable
        // change would otherwise inject a temperature-dependent artificial mode.
        let refreshed = hot_problem
            .solve_from(
                &first.solution.coordinates,
                &first.solution.section_states,
                cx,
            )
            .unwrap();
        near(refreshed.coordinates[0][5], 4.0 * spin, 1e-12);
    });
}

#[test]
fn g1_reference_gauge_and_reference_mass_share_the_original_dynamic_solve() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let cold = section(IsotropicElastic::new(E, NU, 0.05).unwrap());
        let hot = section(IsotropicElastic::new(0.25 * E, NU, 0.05).unwrap());
        let prescribed: Vec<_> = (0..18)
            .map(|d| (d, if d == 5 { 0.008 } else { 0.0 }))
            .collect();
        let p = problem(&mesh, &cold, &prescribed, &[], &[1.0], 0.37);
        let gauge = p.prepare_reference_drilling(cx).unwrap();
        let inertia = {
            let fields = [&cold];
            let mut temporary = problem(&mesh, &cold, &prescribed, &[], &[1.0], 0.37);
            temporary.sections = PlasticShellSections::PerElementBorrowed(&fields);
            PlasticShellReferenceInertia::prepare(
                &mesh,
                &temporary,
                PlasticShellDensity::Uniform(2500.0),
                cx,
            )
            .unwrap()
        };
        let velocity = vec![[0.0; 6]; mesh.nodes.len()];
        let reference = p
            .advance_backward_euler_with_reference_drilling(
                None, &velocity, &inertia, &gauge, 0.1, cx,
            )
            .unwrap();
        let softened = problem(&mesh, &hot, &prescribed, &[], &[1.0], 0.37)
            .advance_backward_euler_with_reference_drilling(
                None, &velocity, &inertia, &gauge, 0.1, cx,
            )
            .unwrap();
        assert_eq!(
            reference.solution.coordinates,
            softened.solution.coordinates
        );
        assert_eq!(reference.velocity, softened.velocity);
        assert_eq!(reference.inertia, softened.inertia);
        assert_eq!(reference.drilling, softened.drilling);
        for (a, b) in reference.solution.increments[0]
            .reactions
            .iter()
            .zip(&softened.solution.increments[0].reactions)
        {
            near(*a, *b, 1e-13);
        }
        let weight = DRILLING_ALPHA * (E * H / (2.0 * (1.0 + NU))) * (0.5 / 3.0);
        near(
            softened.drilling.final_j,
            0.5 * weight * 0.008_f64.powi(2),
            1e-15,
        );
        near(softened.inertia.physical.kinetic_after_j, 0.0, 1e-15);
        assert!(softened.inertia.drilling_gauge.kinetic_after_j > 0.0);
    });
}

#[test]
fn g3_reference_binding_and_borrowed_sections_refuse_mismatches_atomically() {
    with_cx(|gate, cx| {
        let mesh = triangle();
        let section = section(IsotropicElastic::new(E, NU, 0.05).unwrap());
        let fixed: Vec<_> = (0..18).map(|d| (d, 0.0)).collect();
        let p = problem(&mesh, &section, &fixed, &[], &[1.0], 0.37);
        let reference = p.prepare_reference_drilling(cx).unwrap();
        let initial = p.solve(cx).unwrap();
        let saved = initial.clone();
        for (length, alpha, path) in [
            (0.74, DRILLING_ALPHA, vec![1.0]),
            (0.37, 2.0 * DRILLING_ALPHA, vec![1.0]),
            (0.37, DRILLING_ALPHA, vec![0.5, 1.0]),
        ] {
            let mut incompatible = problem(&mesh, &section, &fixed, &[], &path, length);
            incompatible.drilling_alpha = alpha;
            assert!(
                incompatible
                    .solve_from_with_reference_drilling(
                        &initial.coordinates,
                        &initial.section_states,
                        &reference,
                        cx,
                    )
                    .is_err()
            );
        }
        let mut changed = mesh.clone();
        changed.nodes[1][0] *= 1.1;
        assert!(
            problem(&changed, &section, &fixed, &[], &[1.0], 0.37)
                .solve_from_with_reference_drilling(
                    &initial.coordinates,
                    &initial.section_states,
                    &reference,
                    cx,
                )
                .is_err()
        );
        let mut layers = section.layers().to_vec();
        layers[0].bottom_z_m *= 2.0;
        layers[0].top_z_m *= 2.0;
        let thicker = LayeredPlaneStressSection::new(
            layers,
            ThicknessQuadrature::Gauss3,
            section.reduction_config(),
            SectionBudget {
                maximum_layers: 4,
                maximum_quadrature_points: 64,
            },
        )
        .unwrap();
        assert!(
            problem(&mesh, &thicker, &fixed, &[], &[1.0], 0.37)
                .solve_from_with_reference_drilling(
                    &initial.coordinates,
                    &initial.section_states,
                    &reference,
                    cx,
                )
                .is_err()
        );
        let empty = [];
        let mut missing = problem(&mesh, &section, &fixed, &[], &[1.0], 0.37);
        missing.sections = PlasticShellSections::PerElementBorrowed(&empty);
        assert!(
            missing
                .solve_from_with_reference_drilling(
                    &initial.coordinates,
                    &initial.section_states,
                    &reference,
                    cx,
                )
                .is_err()
        );
        gate.request();
        assert!(matches!(
            p.solve_from_with_reference_drilling(
                &initial.coordinates,
                &initial.section_states,
                &reference,
                cx,
            ),
            Err(SmallStrainTetError::Cancelled)
        ));
        assert_eq!(initial, saved);
    });
}
