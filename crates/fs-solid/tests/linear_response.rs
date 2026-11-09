//! G1/G3: native plane elasticity storage, work, stresses and support forces.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::{Formulation, LinearProblem, Mesh2, Patch, PlaneKind, SolidError};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 721,
                    kernel_id: 721,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn problem(mesh: &Mesh2, plane: PlaneKind, formulation: Formulation) -> LinearProblem<'_> {
    LinearProblem {
        mesh,
        youngs: 2400.0,
        poisson: 0.3,
        plane,
        formulation,
        body_force: None,
        dirichlet: vec![],
        traction: vec![],
        symmetry: vec![],
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * scale.max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_uniaxial_plane_equilibria_retain_thickness_stress_storage_and_reactions() {
    with_cx(|cx, _| {
        for mesh in [
            Mesh2::quads(2.0, 1.0, 3, 2),
            Mesh2::triangles(2.0, 1.0, 3, 2),
        ] {
            for plane in [PlaneKind::Stress, PlaneKind::Strain] {
                for formulation in [Formulation::Standard, Formulation::BBar] {
                    for scale in [1.0, 3.0] {
                        let traction = |_: f64, _: f64| [12.0, 0.0];
                        let mut p = problem(&mesh, plane, formulation);
                        p.traction.push((Patch::Right, &traction));
                        p.symmetry = vec![(Patch::Left, 0), (Patch::Bottom, 1)];
                        let scales = vec![scale; mesh.elems.len()];
                        let u = p.solve_with_element_scales(&scales).unwrap();
                        let response = p.response_with_element_scales(&u, &scales, cx).unwrap();
                        let axial = 12.0 / (p.youngs * scale);
                        let expected = match plane {
                            PlaneKind::Stress => [axial, -p.poisson * axial, -p.poisson * axial],
                            PlaneKind::Strain => [
                                (1.0 - p.poisson.powi(2)) * axial,
                                -p.poisson * (1.0 + p.poisson) * axial,
                                0.0,
                            ],
                        };
                        for point in &response.quadrature {
                            for i in 0..3 {
                                close(point.strain[i], expected[i], 0.01);
                            }
                            close(point.stress_pa[0], 12.0, 12.0);
                            close(point.stress_pa[1], 0.0, 12.0);
                            close(
                                point.stress_pa[2],
                                if plane == PlaneKind::Strain { 3.6 } else { 0.0 },
                                12.0,
                            );
                            close(point.stress_pa[3], 0.0, 12.0);
                        }
                        close(response.stored_energy_j_m, 12.0 * expected[0], 0.1);
                        let reactions: [f64; 2] = std::array::from_fn(|component| {
                            response
                                .support_reactions_n_m
                                .iter()
                                .filter(|(dof, _)| dof % 2 == component)
                                .map(|(_, r)| r)
                                .sum()
                        });
                        close(reactions[0], -12.0, 12.0);
                        close(reactions[1], 0.0, 12.0);
                        for (node, r) in response.nodal_residual_n_m.iter().enumerate() {
                            for (component, force) in r.iter().enumerate() {
                                if !response
                                    .support_reactions_n_m
                                    .iter()
                                    .any(|(d, _)| *d == 2 * node + component)
                                {
                                    close(*force, 0.0, 12.0);
                                }
                            }
                        }
                    }
                }
            }
        }
    });
}

#[test]
fn g3_nonaffine_bbar_storage_and_force_share_the_original_stiffness() {
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let (a, b) = (0.02, -0.01);
    let u: Vec<_> = mesh
        .nodes
        .iter()
        .map(|p| [a * p[0] * p[1], b * p[0] * p[1]])
        .collect();
    with_cx(|cx, _| {
        for plane in [PlaneKind::Stress, PlaneKind::Strain] {
            for formulation in [Formulation::Standard, Formulation::BBar] {
                let mut p = problem(&mesh, plane, formulation);
                p.poisson = 0.4999;
                let response = p.response(&u, cx).unwrap();
                let (lambda, mu) = fs_solid::linear::lame(p.youngs, p.poisson, plane);
                let integrated_trace_squared = match formulation {
                    Formulation::Standard => (a * a + b * b) / 3.0 + a * b / 2.0,
                    Formulation::BBar => ((a + b) / 2.0).powi(2),
                };
                let expected_energy =
                    0.5 * (lambda + mu) * integrated_trace_squared + mu * (a * a + b * b) / 3.0;
                close(response.stored_energy_j_m, expected_energy, expected_energy);
                if formulation == Formulation::BBar {
                    for point in &response.quadrature {
                        close(point.strain[0] + point.strain[1], (a + b) / 2.0, 0.01);
                        let raw = point.compatible_in_plane_strain;
                        close(raw[0], a * point.reference_position_m[1], 0.02);
                        close(raw[1], b * point.reference_position_m[0], 0.02);
                    }
                }
                // Existing public group operator uses the original element
                // stiffness, independent of this stress/virtual-work observer.
                let map: Vec<_> = (0..8).map(Some).collect();
                let matrix = fs_solid::stability::group_stiffness(&p, &map, &|_| true).to_dense();
                for dof in 0..8 {
                    let force: f64 = (0..8).map(|j| matrix[8 * dof + j] * u[j / 2][j % 2]).sum();
                    let recovered = response.nodal_residual_n_m[dof / 2][dof % 2];
                    close(recovered, force, 1.0 + force.abs());
                    let mut plus = u.clone();
                    let mut minus = u.clone();
                    let h = 1e-6;
                    plus[dof / 2][dof % 2] += h;
                    minus[dof / 2][dof % 2] -= h;
                    let derivative = (p.response(&plus, cx).unwrap().stored_energy_j_m
                        - p.response(&minus, cx).unwrap().stored_energy_j_m)
                        / (2.0 * h);
                    close(derivative, recovered, 1.0 + recovered.abs());
                }
            }
        }
    });
}

#[test]
fn g1_support_reactions_include_body_and_edge_loads_on_prescribed_rows() {
    with_cx(|cx, _| {
        for mesh in [
            Mesh2::quads(1.0, 1.0, 1, 1),
            Mesh2::triangles(1.0, 1.0, 1, 1),
        ] {
            let zero = |_: f64, _: f64| [0.0; 2];
            let body = |_: f64, _: f64| [3.0, -4.0];
            let edge = |_: f64, _: f64| [10.0, 8.0];
            let mut p = problem(&mesh, PlaneKind::Stress, Formulation::BBar);
            p.body_force = Some(&body);
            p.traction = vec![(Patch::Right, &edge)];
            p.dirichlet = vec![(Patch::Left, &zero), (Patch::Right, &zero)];
            let u = p.solve().unwrap();
            let response = p.response(&u, cx).unwrap();
            assert_eq!(response.stored_energy_j_m, 0.0);
            let mut resultant = [0.0; 2];
            let mut moment = 0.0;
            for &(dof, reaction) in &response.support_reactions_n_m {
                resultant[dof % 2] += reaction;
                let [x, y] = mesh.nodes[dof / 2];
                moment += if dof % 2 == 0 {
                    -y * reaction
                } else {
                    x * reaction
                };
            }
            close(resultant[0], -13.0, 13.0);
            close(resultant[1], -4.0, 4.0);
            close(moment, 0.5, 1.0);
        }
    });
}

#[test]
fn g3_shear_and_translation_retain_physical_work_on_skew_elements() {
    let mesh = Mesh2::mapped_quads(2, 2, &|s, t| [2.0 * s + 0.3 * t, t + 0.2 * s]);
    let gamma = 0.002;
    with_cx(|cx, _| {
        for plane in [PlaneKind::Stress, PlaneKind::Strain] {
            for formulation in [Formulation::Standard, Formulation::BBar] {
                let p = problem(&mesh, plane, formulation);
                let mu = p.youngs / (2.0 * (1.0 + p.poisson));
                let u: Vec<_> = mesh
                    .nodes
                    .iter()
                    .map(|p| [0.5 * gamma * p[1], 0.5 * gamma * p[0]])
                    .collect();
                let response = p.response(&u, cx).unwrap();
                for point in &response.quadrature {
                    close(point.strain[3], gamma / 2.0, gamma);
                    close(point.stress_pa[3], mu * gamma, mu * gamma);
                    close(
                        point.stored_energy_density_j_m3,
                        0.5 * mu * gamma * gamma,
                        0.01,
                    );
                }
                close(
                    response.stored_energy_j_m,
                    1.94 * 0.5 * mu * gamma * gamma,
                    0.01,
                );
                let rigid = vec![[2.0_f64.powi(40), -2.0_f64.powi(42)]; mesh.nodes.len()];
                let translated = p.response(&rigid, cx).unwrap();
                assert_eq!(translated.stored_energy_j_m, 0.0);
                assert!(
                    translated
                        .nodal_residual_n_m
                        .iter()
                        .flatten()
                        .all(|&v| v == 0.0)
                );
                assert!(
                    translated
                        .quadrature
                        .iter()
                        .all(|q| q.strain == [0.0; 6] && q.stress_pa == [0.0; 6])
                );
            }
        }
    });
}

#[test]
fn g3_heterogeneous_scales_change_original_stress_and_energy_without_averaging() {
    let mesh = Mesh2::quads(2.0, 1.0, 2, 1);
    let u: Vec<_> = mesh
        .nodes
        .iter()
        .map(|p| [0.001 * p[0], -0.0003 * p[1]])
        .collect();
    with_cx(|cx, _| {
        let p = problem(&mesh, PlaneKind::Stress, Formulation::Standard);
        let base = p.response(&u, cx).unwrap();
        let scaled = p.response_with_element_scales(&u, &[1.0, 4.0], cx).unwrap();
        for (a, b) in base.quadrature.iter().zip(&scaled.quadrature) {
            let scale = if a.element == 0 { 1.0 } else { 4.0 };
            assert_eq!(a.strain, b.strain);
            for i in 0..6 {
                close(b.stress_pa[i], scale * a.stress_pa[i], 10.0);
            }
            close(
                b.stored_energy_density_j_m3,
                scale * a.stored_energy_density_j_m3,
                0.01,
            );
        }
        close(scaled.stored_energy_j_m, 2.5 * base.stored_energy_j_m, 0.01);
    });
}

#[test]
fn g0_solving_and_recovery_share_material_geometry_and_load_refusals() {
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let u = vec![[0.0; 2]; 4];
    with_cx(|cx, gate| {
        for (youngs, poisson) in [(-1.0, 0.3), (1.0, 0.6), (f64::INFINITY, 0.0)] {
            let mut p = problem(&mesh, PlaneKind::Stress, Formulation::Standard);
            p.youngs = youngs;
            p.poisson = poisson;
            assert!(matches!(p.solve(), Err(SolidError::InvalidInput { .. })));
            assert!(matches!(
                p.response(&u, cx),
                Err(SolidError::InvalidInput { .. })
            ));
        }
        let mut broken = mesh.clone();
        broken.elems[0][0] = usize::MAX;
        let p = problem(&broken, PlaneKind::Strain, Formulation::BBar);
        assert!(p.solve().is_err());
        assert!(p.response(&u, cx).is_err());
        let mut p = problem(&mesh, PlaneKind::Stress, Formulation::Standard);
        assert!(p.response(&u[..3], cx).is_err());
        assert!(p.response_with_element_scales(&u, &[0.0], cx).is_err());
        p.symmetry.push((Patch::Left, 2));
        assert!(p.solve().is_err());
        assert!(p.response(&u, cx).is_err());
        p.symmetry.clear();
        let zero = |_: f64, _: f64| [0.0; 2];
        let invalid = |_: f64, _: f64| [f64::NAN, 0.0];
        p.dirichlet = vec![(Patch::Left, &zero), (Patch::Right, &zero)];
        p.body_force = Some(&invalid);
        assert!(p.solve().is_err());
        assert!(p.response(&u, cx).is_err());
        p.body_force = None;
        p.traction.push((Patch::Left, &invalid));
        assert!(p.solve().is_err());
        assert!(p.response(&u, cx).is_err());
        p.traction.clear();
        gate.request();
        assert!(matches!(p.response(&u, cx), Err(SolidError::Cancelled)));
    });
}
