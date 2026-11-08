//! Physical inventory in the original loaded equilibrium and vibration owner.
use super::super::tests::reference_newton;
use super::*;
use fs_material::sheet::SheetMatrix;
use fs_sparse::Csr;

fn fixture() -> (ShellMesh, Vec<PlateSection>, Vec<usize>, Vec<ShellBeam>) {
    let mut nodes = Vec::new();
    let mut triangles = Vec::new();
    let mut fixed = Vec::new();
    for j in 0..3 {
        for i in 0..3 {
            nodes.push([i as f64 * 0.25, j as f64 * 0.25, 0.0]);
            if i != 1 || j != 1 {
                fixed.push(3 * j + i);
            }
        }
    }
    for j in 0..2 {
        for i in 0..2 {
            let a = 3 * j + i;
            triangles.push([a, a + 1, a + 4]);
            triangles.push([a, a + 4, a + 3]);
        }
    }
    let mesh = ShellMesh::new(nodes, triangles).unwrap();
    let sections = mesh
        .tris
        .iter()
        .map(|_| {
            PlateSection::orthotropic_plane_stress_at_angle(
                9e9, 0.6e9, 0.3, 0.4e9, 0.004, 430.0, 0.0,
            )
            .unwrap()
        })
        .collect();
    let beams = vec![ShellBeam {
        nodes: [3, 4],
        offsets_m: [[0.0, 0.0, 0.03]; 2],
        section_z: [0.0, 0.0, 1.0],
        section: crate::shell::stiffened::BeamSection {
            young_pa: 8e9,
            shear_pa: 0.8e9,
            density_kg_m3: 450.0,
            area_m2: 0.0001,
            iy_m4: 1e-8,
            iz_m4: 3e-9,
            torsion_m4: 4e-9,
        },
    }];
    (mesh, sections, fixed, beams)
}

fn inventory(mesh: &ShellMesh, scale: f64) -> Vec<[f64; 3]> {
    (0..mesh.tris.len())
        .map(|i| {
            [
                scale * (200.0 + 10.0 * i as f64),
                scale * (650.0 - 5.0 * i as f64),
                scale * (450.0 + 12.0 * i as f64),
            ]
        })
        .collect()
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

fn same_matrix(a: &Csr, b: &Csr) {
    assert_eq!(a.nrows(), b.nrows());
    for i in 0..a.nrows() {
        assert_eq!(a.row(i), b.row(i));
    }
}

fn quadratic(a: &Csr, v: &[f64]) -> f64 {
    let mut product = vec![0.0; v.len()];
    a.spmv(v, &mut product);
    v.iter().zip(product).map(|(x, y)| x * y).sum()
}

#[test]
fn loaded_reinforcement_retains_original_graded_mass_without_changing_dead_load_equilibrium() {
    let (mesh, sections, fixed, beams) = fixture();
    let saved = mesh.clone();
    let density = inventory(&mesh, 1.0);
    let beam_density = [[100.0, 900.0]];
    let mut load = vec![0.0; 6 * mesh.nodes.len()];
    load[26] = -10.0;
    let ordinary = super::super::equilibrate_stiffened_shell(
        &mesh,
        &sections,
        &fixed,
        ShellSupport::Clamped,
        &beams,
        &load,
        PreloadOptions::default(),
        &mut reference_newton,
    )
    .unwrap();
    for mass in [BeamMass::Lumped, BeamMass::Consistent] {
        let original = assemble_stiffened_shell_with_reference_density(
            &mesh,
            &sections,
            &density,
            &fixed,
            ShellSupport::Clamped,
            &beams,
            &beam_density,
            mass,
        )
        .unwrap();
        let loaded = equilibrate_stiffened_shell_with_reference_density(
            &mesh,
            &sections,
            &density,
            &fixed,
            ShellSupport::Clamped,
            &beams,
            &beam_density,
            mass,
            &load,
            &[],
            PreloadOptions::default(),
            &mut reference_newton,
        )
        .unwrap();
        same_matrix(&loaded.model.m, &original.m);
        same_matrix(&loaded.model.k, &ordinary.model.k);
        for (&a, &b) in loaded.displacement.iter().zip(&ordinary.displacement) {
            near(a, b, 2e-11);
        }
        near(loaded.stored_energy_j, ordinary.stored_energy_j, 2e-11);
        assert!(loaded.displacement[26] < 0.0);
        assert!(loaded.residual_force_n < 1e-6);
        assert!(loaded.iterations > 0);
        assert!(loaded.model.m.get(0, 0) != ordinary.model.m.get(0, 0));
        // Preload changes the native tangent while reference inertia survives.
        let z = loaded.model.dof_map[26].unwrap();
        assert!((loaded.model.k.get(z, z) - original.k.get(z, z)).abs() > 1e-3);
        for &node in &fixed {
            assert!(
                loaded.displacement[6 * node..6 * node + 6]
                    .iter()
                    .all(|v| *v == 0.0)
            );
        }
    }
    assert_eq!(mesh, saved);
}

#[test]
fn material_sections_and_tethers_keep_density_scaling_in_the_loaded_vibration_pencil() {
    let (mesh, plates, fixed, beams) = fixture();
    let sections: Vec<_> = plates
        .into_iter()
        .map(|plate| {
            ShellSection::new(
                plate,
                SheetMatrix::orthotropic(80e6, 12e6, 2e6, 4e6).unwrap(),
            )
            .unwrap()
        })
        .collect();
    let tethers = [ShellTether {
        nodes: mesh.tris[0],
        weights: [0.0, 0.0, 1.0],
        arm_m: [0.0; 3],
        anchor_m: [0.75, 0.25, -0.04],
        reference_tension_n: 10.0,
        axial_stiffness_n_per_m: 2000.0,
    }];
    let mut load = vec![0.0; 6 * mesh.nodes.len()];
    load[26] = -3.0;
    let mut results = Vec::new();
    for scale in [1.0, 2.0] {
        let density = inventory(&mesh, scale);
        let beam_density = [[100.0 * scale, 900.0 * scale]];
        let original = assemble_stiffened_shell_materials_with_reference_density(
            &mesh,
            &sections,
            &density,
            &fixed,
            ShellSupport::Clamped,
            &beams,
            &beam_density,
            BeamMass::Consistent,
        )
        .unwrap();
        let loaded = equilibrate_shell_materials_with_reference_density(
            &mesh,
            &sections,
            &density,
            &fixed,
            ShellSupport::Clamped,
            &beams,
            &beam_density,
            BeamMass::Consistent,
            &load,
            &tethers,
            PreloadOptions::default(),
            &mut reference_newton,
        )
        .unwrap();
        same_matrix(&loaded.model.m, &original.m);
        assert_eq!(loaded.tether_responses.len(), 1);
        assert!(loaded.tether_responses[0].tension_n > 0.0);
        results.push(loaded);
    }
    let (a, b) = (&results[0], &results[1]);
    same_matrix(&a.model.k, &b.model.k);
    assert_eq!(a.displacement, b.displacement);
    for i in 0..a.model.free {
        for j in 0..a.model.free {
            near(b.model.m.get(i, j), 2.0 * a.model.m.get(i, j), 2e-14);
        }
    }
    // Every Rayleigh quotient, and hence the full generalized spectrum,
    // scales inversely with inventory at fixed equilibrium and stiffness.
    for selector in 0..3 {
        let v: Vec<_> = (0..a.model.free)
            .map(|i| {
                if i % 3 == selector {
                    0.7 + i as f64 * 0.13
                } else {
                    -0.19
                }
            })
            .collect();
        let stiffness = quadratic(&a.model.k, &v);
        let q_a = stiffness / quadratic(&a.model.m, &v);
        let q_b = stiffness / quadratic(&b.model.m, &v);
        assert!(q_a > 0.0 && q_b > 0.0);
        near(q_b, 0.5 * q_a, 3e-14);
    }
}

#[test]
fn invalid_inventory_and_false_solver_convergence_cannot_publish_a_loaded_pencil() {
    let (mesh, sections, fixed, beams) = fixture();
    let saved = mesh.clone();
    let good = inventory(&mesh, 1.0);
    let good_beam = vec![[100.0, 900.0]];
    let mut load = vec![0.0; 6 * mesh.nodes.len()];
    load[26] = -10.0;
    let mut missing = good.clone();
    missing.pop();
    let mut nonfinite = good.clone();
    nonfinite[2][1] = f64::NAN;
    let mut calls = 0;
    for (panel, beam) in [
        (missing, good_beam.clone()),
        (nonfinite, good_beam.clone()),
        (good.clone(), vec![]),
        (good.clone(), vec![[100.0, -1.0]]),
    ] {
        assert!(
            equilibrate_stiffened_shell_with_reference_density(
                &mesh,
                &sections,
                &panel,
                &fixed,
                ShellSupport::Clamped,
                &beams,
                &beam,
                BeamMass::Consistent,
                &load,
                &[],
                PreloadOptions::default(),
                &mut |_, _, _| {
                    calls += 1;
                    Err("invalid inventory reached solver".into())
                },
            )
            .is_err()
        );
    }
    assert_eq!(calls, 0);
    assert!(
        equilibrate_stiffened_shell_with_reference_density(
            &mesh,
            &sections,
            &good,
            &fixed,
            ShellSupport::Clamped,
            &beams,
            &good_beam,
            BeamMass::Consistent,
            &load,
            &[],
            PreloadOptions::default(),
            &mut |problem, _, _| Ok((vec![0.0; problem.dimension()], 0)),
        )
        .is_err()
    );
    assert_eq!(mesh, saved);
    assert_eq!(good, inventory(&mesh, 1.0));
}
