//! G1/G3: the actual plane Newton operator includes consistent physical inertia.
use super::*;
use crate::inertia::ReferenceMass;
use crate::{NewtonSettings, Patch};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 47,
                kernel_id: 61,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn settings() -> NewtonSettings {
    NewtonSettings {
        tol: 2e-9,
        max_iters: 40,
        load_steps: 1,
        max_backtracks: 30,
    }
}

fn law(scale: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 0.8 * scale,
            c01: 0.2 * scale,
            kappa: 25.0 * scale,
        },
        5.0,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn physical_dynamic_force_potential_tangent_and_support_elimination_share_one_operator() {
    with_cx(|cx| {
        for mesh in [
            Mesh2::triangles(1.1, 0.7, 1, 1),
            Mesh2::mapped_quads(1, 1, &|x, y| [1.1 * x + 0.1 * x * y, 0.7 * y + 0.03 * x]),
        ] {
            let material = law(1.0);
            let thickness = 0.13;
            let density = vec![2.7; mesh.elems.len()];
            let mass = ReferenceMass::from_plane(&mesh, thickness, &density, 1 << 20, cx).unwrap();
            let dt = 0.04;
            let old_u: Vec<f64> = mesh
                .nodes
                .iter()
                .flat_map(|[x, y]| [0.03 * x + 0.01 * x * y, -0.02 * y])
                .collect();
            let old_v: Vec<f64> = mesh
                .nodes
                .iter()
                .flat_map(|[x, y]| [0.2 + 0.1 * x, -0.15 + 0.04 * y])
                .collect();
            let u: Vec<f64> = old_u
                .iter()
                .zip(&old_v)
                .enumerate()
                .map(|(i, (u, v))| u + dt * v + 0.002 * (i as f64 + 1.0))
                .collect();
            let inertia = BackwardEulerInertia::new(&mass, &old_u, &old_v, dt, cx).unwrap();
            let traction = |_: f64, _: f64| [0.3, -0.1];
            let problem = HyperProblem {
                mesh: &mesh,
                material: &material,
                dirichlet: vec![],
                traction: vec![(Patch::Right, &traction)],
                settings: settings(),
            };
            let nodal = [(0, 0.04), (0, 0.06), (u.len() - 1, -0.03)];
            for exact in [false, true] {
                let local = PlaneStressSettings::default();
                let prepared =
                    Prepared::dynamic(&problem, local, exact, thickness, &inertia, &nodal, cx)
                        .unwrap();
                let (force, tangent) = prepared
                    .assemble(&problem, &u, 1.0, &BTreeMap::new(), cx)
                    .unwrap();
                let cheap = prepared.force_residual(&problem, &u, 1.0, cx).unwrap();
                let (bare_force, bare_tangent) = if exact {
                    problem.incompressible_plane_stress_residual_and_tangent(&u, 1.0, local, cx)
                } else {
                    problem.plane_stress_residual_and_tangent(&u, 1.0, local, cx)
                }
                .unwrap();
                let mut expected: Vec<f64> = bare_force.iter().map(|f| thickness * f).collect();
                inertia.add_force(&u, &mut expected, cx).unwrap();
                for &(dof, value) in &nodal {
                    expected[dof] -= value;
                }
                for row in 0..u.len() {
                    near(force[row], expected[row], 2e-12);
                    near(cheap[row], force[row], 2e-12);
                    for column in 0..u.len() {
                        near(
                            tangent.get(row, column),
                            thickness * bare_tangent.get(row, column)
                                + mass.matrix().get(row, column) / (dt * dt),
                            3e-12,
                        );
                    }
                }
                // Independent centered differences exercise both external-work
                // terms and the kinetic potential in the actual Newton owner.
                let h = 2e-6;
                for column in 0..u.len() {
                    let mut plus = u.clone();
                    let mut minus = u.clone();
                    plus[column] += h;
                    minus[column] -= h;
                    let ep = prepared.potential_energy(&problem, &plus, 1.0, cx).unwrap();
                    let em = prepared
                        .potential_energy(&problem, &minus, 1.0, cx)
                        .unwrap();
                    near((ep - em) / (2.0 * h), force[column], 2e-7);
                    let fp = prepared.force_residual(&problem, &plus, 1.0, cx).unwrap();
                    let fm = prepared.force_residual(&problem, &minus, 1.0, cx).unwrap();
                    for row in 0..u.len() {
                        near(
                            (fp[row] - fm[row]) / (2.0 * h),
                            tangent.get(row, column),
                            2e-6,
                        );
                    }
                }
                let fixed = BTreeMap::from([(0, u[0]), (u.len() - 1, u[u.len() - 1])]);
                let (supported, eliminated) =
                    prepared.assemble(&problem, &u, 1.0, &fixed, cx).unwrap();
                for row in 0..u.len() {
                    near(
                        supported[row],
                        if fixed.contains_key(&row) {
                            0.0
                        } else {
                            force[row]
                        },
                        2e-12,
                    );
                    for column in 0..u.len() {
                        let expected = if fixed.contains_key(&row) || fixed.contains_key(&column) {
                            if row == column { 1.0 } else { 0.0 }
                        } else {
                            tangent.get(row, column)
                        };
                        near(eliminated.get(row, column), expected, 2e-12);
                    }
                }
                assert!(
                    Prepared::dynamic(&problem, local, exact, 0.0, &inertia, &nodal, cx).is_err()
                );
                assert!(
                    Prepared::dynamic(
                        &problem,
                        local,
                        exact,
                        thickness,
                        &inertia,
                        &[(u.len(), 1.0)],
                        cx
                    )
                    .is_err()
                );
                assert!(
                    Prepared::dynamic(
                        &problem,
                        local,
                        exact,
                        thickness,
                        &inertia,
                        &[(0, f64::NAN)],
                        cx
                    )
                    .is_err()
                );
            }
        }
    });
}

#[test]
fn unsupported_rigid_translation_and_mass_weighted_acceleration_survive_force_unit_changes() {
    with_cx(|cx| {
        for mesh in [
            Mesh2::triangles(1.0, 0.7, 1, 1),
            Mesh2::quads(1.0, 0.7, 2, 1),
        ] {
            let initial = vec![[3.0, -4.0]; mesh.nodes.len()];
            let preserved = initial.clone();
            let old_u: Vec<f64> = initial.iter().flatten().copied().collect();
            let old_v: Vec<f64> = (0..mesh.nodes.len()).flat_map(|_| [0.4, -0.3]).collect();
            let dt = 0.02;
            for scale in [1e-12, 1.0, 1e12] {
                let material = law(scale);
                let density: Vec<f64> = (0..mesh.elems.len())
                    .map(|i| scale * (2.0 + i as f64))
                    .collect();
                let thickness = 0.12;
                let mass =
                    ReferenceMass::from_plane(&mesh, thickness, &density, 1 << 20, cx).unwrap();
                let inertia = BackwardEulerInertia::new(&mass, &old_u, &old_v, dt, cx).unwrap();
                let problem = HyperProblem {
                    mesh: &mesh,
                    material: &material,
                    dirichlet: vec![],
                    traction: vec![],
                    settings: settings(),
                };
                for acceleration in [[0.0, 0.0], [0.7, -0.4]] {
                    let acceleration_field: Vec<f64> =
                        (0..mesh.nodes.len()).flat_map(|_| acceleration).collect();
                    let force = mass.action(&acceleration_field, cx).unwrap();
                    let nodal: Vec<(usize, f64)> = force
                        .iter()
                        .enumerate()
                        .flat_map(|(dof, force)| [(dof, 0.4 * force), (dof, 0.6 * force)])
                        .collect();
                    for exact in [false, true] {
                        let (u, report) = problem
                            .solve_plane_stress_dynamic_from(
                                &initial,
                                PlaneStressSettings::default(),
                                exact,
                                thickness,
                                &inertia,
                                &nodal,
                                cx,
                            )
                            .unwrap();
                        for value in &u {
                            near(value[0], 3.0 + dt * 0.4 + dt * dt * acceleration[0], 2e-10);
                            near(value[1], -4.0 - dt * 0.3 + dt * dt * acceleration[1], 2e-10);
                        }
                        assert!(report.histories.last().unwrap().last().unwrap() < &settings().tol);
                        let flat: Vec<f64> = u.iter().flatten().copied().collect();
                        let residual = problem
                            .plane_stress_dynamic_unconstrained_residual(
                                &flat,
                                1.0,
                                PlaneStressSettings::default(),
                                exact,
                                thickness,
                                &inertia,
                                &nodal,
                                cx,
                            )
                            .unwrap();
                        let residual_norm = residual.iter().fold(0.0_f64, |a, b| a.hypot(*b));
                        assert!(residual_norm / report.force_scales[0] < settings().tol);
                    }
                }
            }
            assert_eq!(initial, preserved);
        }
    });
}
