//! G1: source thermal strain reaches actual thin-solid equilibrium and recovery.
use super::*;
use fs_material::hyper::incompressible_plane_stress::incompressible_plane_stress;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_solid::hyper2d::PlaneMaterial;
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch};

fn settings() -> NewtonSettings {
    NewtonSettings {
        tol: 2e-9,
        max_iters: 40,
        load_steps: 1,
        max_backtracks: 30,
    }
}

#[test]
fn source_anisotropic_expansion_drives_full_and_exact_thin_equilibrium() {
    with_cx(|_, cx| {
        let source = card(STRETCHES, "T", false);
        let material = bind_expansion(&source, Q, 310.0);
        let reference = ExpansionThermalReference::try_new(
            &source,
            path(&source, false, "T", point(300.0)),
            310.0,
            frame(Q),
        )
        .unwrap();
        let g = *reference.reference_deformation().deformation_gradient();
        // The same oriented source supplies both off-diagonal in-plane terms.
        assert!(g[1].abs() > 0.01 && g[3].abs() > 0.01);
        assert!(material.native_hyperelastic().is_none());
        let mesh = Mesh2::triangles(1.0, 0.7, 2, 2);
        let left = |_: f64, y: f64| [g[1] * y, (g[4] - 1.0) * y];
        let problem = HyperProblem {
            mesh: &mesh,
            material: &material,
            dirichlet: vec![(Patch::Left, &left)],
            traction: vec![],
            settings: settings(),
        };
        let local = PlaneStressSettings::default();
        for exact in [false, true] {
            let (displacement, report) = if exact {
                problem.solve_incompressible_plane_stress(local, cx)
            } else {
                problem.solve_plane_stress(local, cx)
            }
            .unwrap();
            assert!(report.histories.iter().all(|h| h.last().unwrap() < &2e-9));
            for (u, [x, y]) in displacement.iter().zip(&mesh.nodes) {
                near(u[0], (g[0] - 1.0) * x + g[1] * y, 3e-8);
                near(u[1], g[3] * x + (g[4] - 1.0) * y, 3e-8);
            }
            let flat: Vec<f64> = displacement.iter().flatten().copied().collect();
            if exact {
                for point in problem
                    .incompressible_plane_stress_states(&flat, local, cx)
                    .unwrap()
                {
                    near(
                        det(&point.deformation_gradient),
                        reference.geometry().relative_volume(),
                        3e-12,
                    );
                    near(point.deformation_gradient[8], g[8], 3e-8);
                    for stress in point.full_piola {
                        near(stress, 0.0, 3e-7);
                    }
                }
            } else {
                for point in problem.plane_stress_states(&flat, local, cx).unwrap() {
                    near(point.deformation_gradient[8], g[8], 3e-8);
                    near(point.energy_density_pa, 0.0, 3e-10);
                    for stress in point.piola {
                        near(stress, 0.0, 3e-7);
                    }
                }
            }
        }
        near(material.reference_density_kg_m3(), 6.0, 0.0);
    });
}

#[test]
fn source_domain_and_large_thermal_volume_survive_thin_recovery() {
    with_cx(|_, cx| {
        // Existing fixture: G=4I with the original elastic domain [1/3,3].
        // Total F=4I is outside the bare law but Fe=I is stress-free.
        let source = super::super::thermal::thermal_card(0.3, 1.0, false);
        let material = super::super::thermal::prepared(&source);
        let total = [4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 4.0];
        assert!(material.law().law().tangent(&total).is_err());
        assert!(material.native_hyperelastic().is_none());
        let mesh = Mesh2::triangles(1.0, 0.7, 1, 1);
        let problem = HyperProblem {
            mesh: &mesh,
            material: &material,
            dirichlet: vec![],
            traction: vec![],
            settings: settings(),
        };
        let flat: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|[x, y]| [3.0 * x, 3.0 * y])
            .collect();
        let local = PlaneStressSettings::default();
        for state in problem.plane_stress_states(&flat, local, cx).unwrap() {
            near(state.deformation_gradient[8], 4.0, 3e-11);
            for stress in state.piola {
                near(stress, 0.0, 3e-10);
            }
        }
        let exact = incompressible_plane_stress(&material, &[4.0, 0.0, 0.0, 4.0], local).unwrap();
        near(det(&exact.deformation_gradient), 64.0, 3e-12);
        for stress in exact.full_piola {
            near(stress, 0.0, 3e-10);
        }
        // This request has elastic in-plane stretches 1/4 and must still fail.
        assert!(incompressible_plane_stress(&material, &[1.0, 0.0, 0.0, 1.0], local).is_err());
        near(material.reference_density_kg_m3(), 6.0, 0.0);
    });
}

#[test]
fn source_exact_reduction_tangent_matches_force_differences() {
    let source = card(STRETCHES, "T", false);
    let material = bind_expansion(&source, Q, 310.0);
    let a = [1.18, 0.08, 0.04, 0.93];
    let local = PlaneStressSettings::default();
    let center = incompressible_plane_stress(&material, &a, local).unwrap();
    let h = 2e-6;
    for column in 0..4 {
        let mut plus = a;
        let mut minus = a;
        plus[column] += h;
        minus[column] -= h;
        let p = incompressible_plane_stress(&material, &plus, local).unwrap();
        let m = incompressible_plane_stress(&material, &minus, local).unwrap();
        near(
            (p.energy_density_pa - m.energy_density_pa) / (2.0 * h),
            center.piola[column],
            3e-7,
        );
        for row in 0..4 {
            near(
                (p.piola[row] - m.piola[row]) / (2.0 * h),
                center.tangent[row][column],
                3e-6,
            );
        }
    }
}
