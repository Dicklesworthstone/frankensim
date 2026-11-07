//! Native retained plastic shell motion, consistent inertia and support work.

use super::*;
use fs_plate::shell::{
    ReferenceDensityMoments, assemble_shell_sections_with_reference_density_moments,
};
use fs_solid::plastic_shell::dynamics::PlasticShellDensity;

fn triangle() -> ShellMesh {
    ShellMesh::new(
        vec![[0.2, -0.4, 0.8], [1.4, 0.3, 1.1], [-0.1, 0.8, 1.6]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|c| a[c] * b[c]).sum()
}

#[test]
fn g1_free_shell_rigid_motion_remains_ballistic_without_artificial_supports() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let material = section(metal());
        let p = problem(&mesh, &material, &[], &[], &[1.0], 0.4);
        let rho = 2500.0;
        let reference = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(rho), cx)
            .unwrap();
        let omega = [0.003, -0.002, 0.004];
        let translation = [0.01, -0.015, 0.02];
        let velocity: Vec<[f64; 6]> = mesh
            .nodes
            .iter()
            .map(|&x| {
                let spin = cross(omega, x);
                core::array::from_fn(|c| {
                    if c < 3 {
                        translation[c] + spin[c]
                    } else {
                        omega[c - 3]
                    }
                })
            })
            .collect();
        let g = mesh.facet(0).unwrap();
        let mass = rho * H * g.area_m2;
        let mut expected_kinetic = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                let gram = if i == j { 1.0 / 6.0 } else { 1.0 / 12.0 };
                expected_kinetic += 0.5
                    * mass
                    * gram
                    * (0..3).map(|c| velocity[i][c] * velocity[j][c]).sum::<f64>();
            }
        }
        expected_kinetic +=
            0.5 * mass * H * H / 12.0 * (dot3(omega, omega) - dot3(omega, g.frame[2]).powi(2));
        let dt = 0.03;
        let first = p
            .advance_backward_euler(None, &velocity, &reference, dt, cx)
            .unwrap();
        let second = p
            .advance_backward_euler(Some(&first.solution), &first.velocity, &reference, dt, cx)
            .unwrap();
        for (step, time) in [(&first, dt), (&second, 2.0 * dt)] {
            assert!(step.solution.increments[0].reactions.is_empty());
            assert_eq!(step.solution.increments[0].evolving_elements, 0);
            for i in 0..mesh.nodes.len() {
                for c in 0..6 {
                    near(
                        step.solution.coordinates[i][c],
                        time * velocity[i][c],
                        2e-10,
                    );
                    near(step.velocity[i][c], velocity[i][c], 2e-9);
                }
            }
            near(step.inertia.total.reference_mass_kg, mass, 2e-11);
            near(
                step.inertia.physical.kinetic_before_j,
                expected_kinetic,
                2e-10,
            );
            near(
                step.inertia.physical.kinetic_after_j,
                expected_kinetic,
                2e-10,
            );
            near(step.inertia.physical.numerical_loss_j, 0.0, 1e-15);
            near(step.inertia.drilling_gauge.kinetic_after_j, 0.0, 1e-16);
            near(step.inertia.total.inertial_work_j, 0.0, 2e-10);
        }
    });
}

#[test]
fn g1_free_coordinate_and_moving_supports_match_original_consistent_mass_equation() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let material = LayeredPlaneStressSection::new(
            vec![
                PlaneStressLayer {
                    bottom_z_m: -0.5 * H,
                    top_z_m: 0.0,
                    subdivisions: 2,
                    law: metal(),
                },
                PlaneStressLayer {
                    bottom_z_m: 0.0,
                    top_z_m: 0.5 * H,
                    subdivisions: 2,
                    law: metal(),
                },
            ],
            ThicknessQuadrature::Gauss3,
            section(metal()).reduction_config(),
            SectionBudget {
                maximum_layers: 2,
                maximum_quadrature_points: 64,
            },
        )
        .unwrap();
        for density in [
            vec![[700.0; 3], [2700.0; 3]],
            vec![[300.0, 900.0, 900.0], [2000.0, 2700.0, 3400.0]],
        ] {
            // Exact normalized moments at each corner of the two real layers.
            let moments = core::array::from_fn(|corner| {
                let lower = density[0][corner];
                let upper = density[1][corner];
                ReferenceDensityMoments::new(
                    (lower + upper) / 2.0,
                    (upper - lower) / 8.0,
                    (lower + upper) / 24.0,
                )
                .unwrap()
            });
            let dt = 0.07;
            let free = 14; // One global transverse translation remains genuinely free.
            let loads = [(free, 0.001), (0, -0.0003)];
            let prescribed: Vec<_> = (0..18)
                .filter(|&d| d != free)
                .map(|d| (d, 2e-6 * ((d * 7 % 9) as f64 - 4.0)))
                .collect();
            let velocity: Vec<[f64; 6]> = (0..3)
                .map(|i| core::array::from_fn(|c| 1e-4 * (((6 * i + c) * 5 % 13) as f64 - 6.0)))
                .collect();
            let original = assemble_shell_sections_with_reference_density_moments(
                &mesh,
                &[PlateSection::isotropic(E, NU, H, 1700.0).unwrap()],
                &[moments],
                &[],
                ShellSupport::Free,
            )
            .unwrap();
            let p = problem(&mesh, &material, &prescribed, &loads, &[1.0], 0.3);
            let field = [density.clone()];
            let reference = p
                .prepare_reference_inertia(PlasticShellDensity::PerElementCorners(&field), cx)
                .unwrap();
            assert_eq!(reference.reference_layer_densities(0), None);
            assert_eq!(
                reference
                    .reference_layer_corner_densities(0)
                    .unwrap()
                    .collect::<Vec<_>>(),
                density
            );
            // Independent exact simplex integration of rho(N) N_i N_j, including
            // each real layer's signed first and second thickness moments.
            let facet = mesh.facet(0).unwrap();
            let physical_kinetic = |densities: &[[f64; 3]], v: &[f64]| {
                let mut energy = 0.0;
                for (layer, rho) in densities.iter().enumerate() {
                    let z0 = (layer as f64 - 1.0) * H / 2.0;
                    let z1 = z0 + H / 2.0;
                    let dz = z1 - z0;
                    let first = (z1 * z1 - z0 * z0) / 2.0;
                    let second = (z1.powi(3) - z0.powi(3)) / 3.0;
                    for i in 0..3 {
                        let vi = core::array::from_fn(|c| v[6 * i + c]);
                        let di = cross(core::array::from_fn(|c| v[6 * i + 3 + c]), facet.frame[2]);
                        for j in 0..3 {
                            let vj = core::array::from_fn(|c| v[6 * j + c]);
                            let dj =
                                cross(core::array::from_fn(|c| v[6 * j + 3 + c]), facet.frame[2]);
                            let density_gram = facet.area_m2
                                * (rho.iter().sum::<f64>() + rho[i] + rho[j])
                                / if i == j { 30.0 } else { 60.0 };
                            energy += 0.5
                                * density_gram
                                * (dz * dot3(vi, vj)
                                    + first * (dot3(vi, dj) + dot3(di, vj))
                                    + second * dot3(di, dj));
                        }
                    }
                }
                energy
            };
            near(
                reference.reference_mass_kg(),
                facet.area_m2 * H * 1700.0,
                1e-11,
            );
            let constant_reference = p
                .prepare_reference_inertia(PlasticShellDensity::Layers(&[700.0, 2700.0]), cx)
                .unwrap();
            if density[0] == [700.0; 3] {
                assert_eq!(
                    reference.reference_mass_kg(),
                    constant_reference.reference_mass_kg()
                );
            } else {
                // All per-layer means match the constant case; the original mass
                // distribution and eccentric inertia nevertheless differ.
                let mean_inventory = [[700.0; 3], [2700.0; 3]];
                assert!(
                    (physical_kinetic(&density, velocity.as_flattened())
                        - physical_kinetic(&mean_inventory, velocity.as_flattened()))
                    .abs()
                        > 1e-10
                );
            }
            let mut old_momentum = vec![0.0; 18];
            original.m.spmv(velocity.as_flattened(), &mut old_momentum);
            let mut expected = vec![0.0; 18];
            for &(d, value) in &prescribed {
                expected[d] = value;
            }
            let effective = |i, j| original.k.get(i, j) + original.m.get(i, j) / (dt * dt);
            expected[free] = (loads[0].1 + old_momentum[free] / dt
                - (0..18)
                    .filter(|&d| d != free)
                    .map(|d| effective(free, d) * expected[d])
                    .sum::<f64>())
                / effective(free, free);
            let step = p
                .advance_backward_euler(None, &velocity, &reference, dt, cx)
                .unwrap();
            assert_eq!(step.solution.increments[0].evolving_elements, 0);
            for d in 0..18 {
                near(step.solution.coordinates[d / 6][d % 6], expected[d], 2e-12);
                near(step.velocity[d / 6][d % 6], expected[d] / dt, 3e-11);
            }
            for (&(d, _), &reaction) in prescribed
                .iter()
                .zip(&step.solution.increments[0].reactions)
            {
                let expected_reaction = (0..18).map(|j| effective(d, j) * expected[j]).sum::<f64>()
                    - old_momentum[d] / dt
                    - loads
                        .iter()
                        .filter(|&&(row, _)| row == d)
                        .map(|&(_, f)| f)
                        .sum::<f64>();
                near(reaction, expected_reaction, 3e-10);
            }
            let next_velocity: Vec<_> = expected.iter().map(|q| q / dt).collect();
            let mut next_momentum = vec![0.0; 18];
            original.m.spmv(&next_velocity, &mut next_momentum);
            let expected_before = 0.5
                * velocity
                    .as_flattened()
                    .iter()
                    .zip(&old_momentum)
                    .map(|(v, p)| v * p)
                    .sum::<f64>();
            let expected_after = 0.5
                * next_velocity
                    .iter()
                    .zip(&next_momentum)
                    .map(|(v, p)| v * p)
                    .sum::<f64>();
            near(step.inertia.total.kinetic_before_j, expected_before, 2e-12);
            near(step.inertia.total.kinetic_after_j, expected_after, 2e-12);
            near(
                step.inertia.physical.kinetic_before_j,
                physical_kinetic(&density, velocity.as_flattened()),
                2e-12,
            );
            near(
                step.inertia.physical.kinetic_after_j,
                physical_kinetic(&density, &next_velocity),
                2e-12,
            );
            assert!(step.inertia.total.numerical_loss_j > 0.0);
        }
    });
}

#[test]
fn g3_dynamic_release_retains_plastic_history_independently_of_rotation_scaling() {
    with_cx(|_, cx| {
        let (mesh, fixed, loads) = cantilever(2, 1);
        let material = section(metal());
        let loaded = problem(&mesh, &material, &fixed, &loads, &[0.5, 1.0], 1.0)
            .solve(cx)
            .unwrap();
        assert!(
            loaded
                .section_states
                .iter()
                .flat_map(|s| &s.points)
                .any(|p| p.history.alpha > 0.0)
        );
        let saved = loaded.clone();
        let zero_velocity = vec![[0.0; 6]; mesh.nodes.len()];
        let baseline_problem = problem(&mesh, &material, &fixed, &[], &[1.0], 1.0);
        let baseline_mass = baseline_problem
            .prepare_reference_inertia(PlasticShellDensity::Uniform(3.0), cx)
            .unwrap();
        let baseline = baseline_problem
            .advance_backward_euler(Some(&loaded), &zero_velocity, &baseline_mass, 0.5, cx)
            .unwrap();
        assert!(baseline.velocity.iter().flatten().any(|v| v.abs() > 1e-6));
        for length in [0.1, 10.0] {
            let p = problem(&mesh, &material, &fixed, &[], &[1.0], length);
            let reference = p
                .prepare_reference_inertia(PlasticShellDensity::Uniform(3.0), cx)
                .unwrap();
            let step = p
                .advance_backward_euler(Some(&loaded), &zero_velocity, &reference, 0.5, cx)
                .unwrap();
            for (a, b) in baseline
                .solution
                .coordinates
                .iter()
                .flatten()
                .zip(step.solution.coordinates.iter().flatten())
            {
                near(*a, *b, 3e-9);
            }
            for (a, b) in baseline
                .velocity
                .iter()
                .flatten()
                .zip(step.velocity.iter().flatten())
            {
                near(*a, *b, 6e-9);
            }
            for ((a, b), old) in baseline
                .solution
                .section_states
                .iter()
                .zip(&step.solution.section_states)
                .zip(&loaded.section_states)
            {
                for ((a, b), old) in a.points.iter().zip(&b.points).zip(&old.points) {
                    near(a.history.alpha, b.history.alpha, 2e-9);
                    assert!(b.history.alpha >= old.history.alpha);
                    for c in 0..6 {
                        near(
                            a.history.plastic_strain[c],
                            b.history.plastic_strain[c],
                            2e-9,
                        );
                    }
                }
            }
            near(
                step.inertia.physical.kinetic_after_j,
                baseline.inertia.physical.kinetic_after_j,
                2e-11,
            );
            near(
                step.inertia.drilling_gauge.kinetic_after_j,
                baseline.inertia.drilling_gauge.kinetic_after_j,
                2e-11,
            );
        }
        assert_eq!(loaded, saved);
    });
}

#[test]
fn g0_retained_mass_refusals_and_cancellation_leave_accepted_motion_unchanged() {
    with_cx(|gate, cx| {
        let mesh = triangle();
        let material = section(metal());
        let p = problem(&mesh, &material, &[], &[], &[1.0], 1.0);
        let reference = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(1200.0), cx)
            .unwrap();
        let velocity = vec![[0.001, -0.002, 0.003, 0.0, 0.0, 0.0]; mesh.nodes.len()];
        let accepted = p
            .advance_backward_euler(None, &velocity, &reference, 0.02, cx)
            .unwrap();
        let saved = accepted.solution.clone();
        let saved_velocity = accepted.velocity.clone();
        let mut shifted_nodes = mesh.nodes.clone();
        shifted_nodes[1][0] += 0.1;
        let other_mesh = ShellMesh::new(shifted_nodes, mesh.tris.clone()).unwrap();
        assert!(
            problem(&other_mesh, &material, &[], &[], &[1.0], 1.0)
                .advance_backward_euler(None, &velocity, &reference, 0.02, cx)
                .is_err()
        );
        let thicker = LayeredPlaneStressSection::new(
            vec![PlaneStressLayer {
                bottom_z_m: -0.7 * H,
                top_z_m: 0.7 * H,
                subdivisions: 4,
                law: metal(),
            }],
            ThicknessQuadrature::Gauss3,
            material.reduction_config(),
            SectionBudget {
                maximum_layers: 4,
                maximum_quadrature_points: 64,
            },
        )
        .unwrap();
        assert!(
            problem(&mesh, &thicker, &[], &[], &[1.0], 1.0)
                .advance_backward_euler(None, &velocity, &reference, 0.02, cx)
                .is_err()
        );
        assert!(
            p.prepare_reference_inertia(PlasticShellDensity::Uniform(0.0), cx)
                .is_err()
        );
        assert!(
            p.prepare_reference_inertia(PlasticShellDensity::Layers(&[]), cx)
                .is_err()
        );
        assert!(
            p.prepare_reference_inertia(PlasticShellDensity::PerElement(&[]), cx)
                .is_err()
        );
        for field in [
            vec![],
            vec![vec![]],
            vec![vec![[1200.0; 3]], vec![[1200.0; 3]]],
            vec![vec![[1200.0; 3], [1200.0; 3]]],
        ] {
            assert!(
                p.prepare_reference_inertia(PlasticShellDensity::PerElementCorners(&field), cx)
                    .is_err()
            );
        }
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let field = [vec![[1200.0, 1200.0, invalid]]];
            assert!(
                p.prepare_reference_inertia(PlasticShellDensity::PerElementCorners(&field), cx)
                    .is_err()
            );
        }
        for dt in [0.0, -0.1, f64::NAN] {
            assert!(
                p.advance_backward_euler(
                    Some(&accepted.solution),
                    &accepted.velocity,
                    &reference,
                    dt,
                    cx
                )
                .is_err()
            );
        }
        assert!(
            p.advance_backward_euler(
                Some(&accepted.solution),
                &accepted.velocity[..2],
                &reference,
                0.02,
                cx
            )
            .is_err()
        );
        assert!(
            problem(&mesh, &material, &[], &[], &[0.5], 1.0)
                .advance_backward_euler(
                    Some(&accepted.solution),
                    &accepted.velocity,
                    &reference,
                    0.02,
                    cx
                )
                .is_err()
        );
        let excessive: Vec<_> = (0..18)
            .map(|d| (d, if d == 3 { 0.2 } else { 0.0 }))
            .collect();
        assert!(
            problem(&mesh, &material, &excessive, &[], &[1.0], 1.0)
                .advance_backward_euler(
                    Some(&accepted.solution),
                    &accepted.velocity,
                    &reference,
                    0.02,
                    cx
                )
                .is_err()
        );
        assert_eq!(accepted.solution, saved);
        assert_eq!(accepted.velocity, saved_velocity);
        gate.request();
        assert!(matches!(
            p.advance_backward_euler(
                Some(&accepted.solution),
                &accepted.velocity,
                &reference,
                0.02,
                cx
            ),
            Err(SmallStrainTetError::Cancelled)
        ));
        assert_eq!(accepted.solution, saved);
        assert_eq!(accepted.velocity, saved_velocity);
    });
}
