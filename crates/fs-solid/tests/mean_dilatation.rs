//! G0/G1 finite-strain Q1 mean-dilatation checks. These compare independently
//! derived continuum/element identities and actual compression solves; they do
//! not establish arbitrary-mesh stability or a fully incompressible limit.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 49,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn material(mu: f64, lambda: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 8.0).unwrap()
}

fn problem<'a>(mesh: &'a Mesh2, law: &'a Hyperelastic) -> HyperProblem<'a> {
    HyperProblem {
        mesh,
        material: law,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings {
            tol: 1e-8,
            max_iters: 60,
            load_steps: 8,
            max_backtracks: 30,
        },
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn g0_condensed_energy_force_and_hessian_include_cross_quadrature_terms() {
    with_cx(|cx, _| {
        let mesh = Mesh2::mapped_quads(2, 1, &|x, y| [x + 0.18 * x * y, y + 0.07 * x]);
        for law in [
            material(2.0, 30.0),
            Hyperelastic::new(
                HyperelasticModel::MooneyRivlin {
                    c10: 0.7,
                    c01: 0.3,
                    kappa: 50.0,
                },
                8.0,
            )
            .unwrap(),
            Hyperelastic::new(
                HyperelasticModel::Ogden {
                    terms: vec![OgdenTerm {
                        mu_pa: 2.0,
                        alpha: 2.0,
                    }],
                    volumetric: Volumetric::LogSquared { kappa: 40.0 },
                },
                8.0,
            )
            .unwrap(),
            Hyperelastic::new(
                HyperelasticModel::HolzapfelGasserOgden {
                    c10: 0.8,
                    k1: 0.3,
                    k2: 2.0,
                    dispersion: 0.0,
                    fibres: vec![[0.0, 1.0, 0.0]],
                    invariants: FibreInvariants::Full,
                    volumetric: Volumetric::Quadratic { kappa: 40.0 },
                },
                8.0,
            )
            .unwrap(),
        ] {
            let p = problem(&mesh, &law);
            let u: Vec<f64> = mesh
                .nodes
                .iter()
                .flat_map(|[x, y]| [-0.14 * x + 0.03 * x * y, 0.08 * y - 0.02 * x * y])
                .collect();
            let direction: Vec<f64> = (0..u.len())
                .map(|i| ((i * 7 % 11) as f64 - 5.0) / 13.0)
                .collect();
            let (r, k) = p.mean_dilatation_residual_and_tangent(&u, 0.0, cx).unwrap();
            let h = 2e-6;
            let plus: Vec<f64> = u.iter().zip(&direction).map(|(a, b)| a + h * b).collect();
            let minus: Vec<f64> = u.iter().zip(&direction).map(|(a, b)| a - h * b).collect();
            let rp = p
                .mean_dilatation_residual_and_tangent(&plus, 0.0, cx)
                .unwrap()
                .0;
            let rm = p
                .mean_dilatation_residual_and_tangent(&minus, 0.0, cx)
                .unwrap()
                .0;
            let ep = p.mean_dilatation_potential_energy(&plus, 0.0, cx).unwrap();
            let em = p.mean_dilatation_potential_energy(&minus, 0.0, cx).unwrap();
            let derivative = dot(&r, &direction);
            assert!(((ep - em) / (2.0 * h) - derivative).abs() < 2e-7 * (1.0 + derivative.abs()));
            for i in 0..u.len() {
                let exact: f64 = (0..u.len()).map(|j| k.get(i, j) * direction[j]).sum();
                assert!(((rp[i] - rm[i]) / (2.0 * h) - exact).abs() < 2e-6 * (1.0 + exact.abs()));
                for j in 0..u.len() {
                    assert_eq!(k.get(i, j), k.get(j, i));
                }
            }
        }
    });
}

#[test]
fn g1_homogeneous_finite_compression_matches_analytic_energy_and_reactions() {
    with_cx(|cx, _| {
        let (a, b, mu, lambda) = (0.78_f64, 1.17_f64, 2.0_f64, 35.0_f64);
        let j = a * b;
        let expected_energy =
            0.5 * mu * (a * a + b * b - 2.0) - mu * j.ln() + 0.5 * lambda * j.ln().powi(2);
        let px = mu * (a - a.recip()) + lambda * j.ln() / a;
        let py = mu * (b - b.recip()) + lambda * j.ln() / b;
        let law = material(mu, lambda);
        for n in [2, 3] {
            let mesh = Mesh2::quads(1.0, 1.0, n, n);
            let imposed = |x: f64, y: f64| [(a - 1.0) * x, (b - 1.0) * y];
            let mut p = problem(&mesh, &law);
            p.dirichlet = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top]
                .into_iter()
                .map(|patch| (patch, &imposed as _))
                .collect();
            let (solution, report) = p.solve_mean_dilatation(cx).unwrap();
            assert_eq!(report.histories.len(), p.settings.load_steps);
            for (point, u) in mesh.nodes.iter().zip(&solution) {
                let exact = imposed(point[0], point[1]);
                assert!((u[0] - exact[0]).abs() < 3e-8 && (u[1] - exact[1]).abs() < 3e-8);
            }
            let flat: Vec<f64> = solution.iter().flatten().copied().collect();
            let energy = p.mean_dilatation_potential_energy(&flat, 1.0, cx).unwrap();
            assert!((energy - expected_energy).abs() < 2e-8);
            let (reaction, _) = p
                .mean_dilatation_residual_and_tangent(&flat, 1.0, cx)
                .unwrap();
            let rx: f64 = mesh
                .patch_nodes(Patch::Right)
                .iter()
                .map(|&i| reaction[2 * i])
                .sum();
            let ry: f64 = mesh
                .patch_nodes(Patch::Top)
                .iter()
                .map(|&i| reaction[2 * i + 1])
                .sum();
            assert!((rx - px).abs() < 3e-7 && (ry - py).abs() < 3e-7);
            for state in p.mean_dilatation_states(&flat, cx).unwrap() {
                assert!((state.mean_volume_ratio - j).abs() < 3e-8);
                assert!((state.volumetric_cauchy_stress - lambda * j.ln() / j).abs() < 2e-6);
            }
        }
    });
}

#[test]
fn g0_mean_constraint_removes_the_parasitic_bulk_stiffness_of_nonuniform_dilatation() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        // v_x = x(y - 1/2), v_y = 0 has zero mean divergence.
        // Independent integration gives v^T K_iso v = mu/2, while the
        // unprojected log-J penalty adds lambda * integral(y-1/2)^2 = lambda/12.
        let direction: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|[x, y]| [x * (y - 0.5), 0.0])
            .collect();
        let zero = vec![0.0; direction.len()];
        let mu = 2.0;
        for lambda in [10.0, 10_000.0, 100_000_000.0] {
            let law = material(mu, lambda);
            let p = problem(&mesh, &law);
            let (_, projected) = p
                .mean_dilatation_residual_and_tangent(&zero, 0.0, cx)
                .unwrap();
            let (_, standard) = p.residual_and_tangent(&zero, 0.0).unwrap();
            let quadratic = |k: &fs_sparse::Csr| {
                (0..direction.len())
                    .map(|i| {
                        (0..direction.len())
                            .map(|j| direction[i] * k.get(i, j) * direction[j])
                            .sum::<f64>()
                    })
                    .sum::<f64>()
            };
            assert!((quadratic(&projected) - 0.5 * mu).abs() < 2e-8);
            assert!(
                (quadratic(&standard) - (0.5 * mu + lambda / 12.0)).abs() < 2e-8 * (1.0 + lambda)
            );
            // Finite compressed background: Fxx=.8+.25(y-.5), Fyy=1.25.
            // Its actual mean J is one, despite locally distinct volumes.
            let compressed: Vec<f64> = mesh
                .nodes
                .iter()
                .flat_map(|[x, y]| [-0.2 * x + 0.25 * x * (y - 0.5), 0.25 * y])
                .collect();
            let state = p.mean_dilatation_states(&compressed, cx).unwrap()[0];
            assert!((state.mean_volume_ratio - 1.0).abs() < 4e-15);
            assert!(
                state.minimum_quadrature_volume_ratio < 0.95
                    && state.maximum_quadrature_volume_ratio > 1.05
            );
            let gauss = [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()];
            let mut independent = 0.0;
            for x in gauss {
                for y in gauss {
                    let fxx = 0.8 + 0.25 * (y - 0.5);
                    let fxy = 0.25 * x;
                    let det = fxx * 1.25;
                    independent += 0.25
                        * (0.5 * mu * (fxx * fxx + fxy * fxy + 1.25_f64.powi(2) - 2.0)
                            - mu * det.ln());
                }
            }
            assert!(
                (p.mean_dilatation_potential_energy(&compressed, 0.0, cx)
                    .unwrap()
                    - independent)
                    .abs()
                    < 2e-12
            );
        }
    });
}

#[test]
fn g1_nonuniform_compression_responds_without_the_standard_q1_bulk_constraint() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(2.0, 1.0, 4, 2);
        let law = material(1.0, 10_000.0);
        let clamp = |_: f64, _: f64| [0.0, 0.0];
        let compress = |_: f64, _: f64| [-0.02, 0.0];
        let mut p = problem(&mesh, &law);
        p.dirichlet = vec![(Patch::Left, &clamp), (Patch::Right, &compress)];
        let projected = p.solve_mean_dilatation(cx).unwrap().0;
        let standard = p.solve().unwrap().0;
        let flat = |u: &Vec<[f64; 2]>| u.iter().flatten().copied().collect::<Vec<_>>();
        let pf = flat(&projected);
        let sf = flat(&standard);
        let pr = p
            .mean_dilatation_residual_and_tangent(&pf, 1.0, cx)
            .unwrap()
            .0;
        let sr = p.residual_and_tangent(&sf, 1.0).unwrap().0;
        let reaction = |r: &[f64]| {
            mesh.patch_nodes(Patch::Right)
                .iter()
                .map(|&a| -r[2 * a])
                .sum::<f64>()
        };
        let projected_force = reaction(&pr);
        let standard_force = reaction(&sr);
        assert!(projected_force > 0.0);
        // Independent reference: analytic NH Piola/cofactor derivatives,
        // explicit 2x2 integration and a dense Newton solve of this mesh.
        // These are discrete-model values, not a continuum accuracy claim.
        assert!((projected_force - 0.048_903_180_828_632_375).abs() < 1e-6);
        assert!((standard_force - 6.394_411_794_327_845).abs() < 1e-6);
        assert!(
            standard_force > 4.0 * projected_force,
            "standard={standard_force}, mean-dilatation={projected_force}"
        );
        for i in mesh.patch_nodes(Patch::Top) {
            if mesh.nodes[i][0] > 0.0 && mesh.nodes[i][0] < 2.0 {
                assert!(projected[i][1] > 0.0);
            }
        }
        for state in p.mean_dilatation_states(&pf, cx).unwrap() {
            assert!((state.mean_volume_ratio - 1.0).abs() < 5e-4);
        }
    });
}

#[test]
fn g0_objective_under_superposed_rigid_motion_and_refuses_hidden_local_failure() {
    with_cx(|cx, gate| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let law = material(2.0, 40.0);
        let p = problem(&mesh, &law);
        let u: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|[x, y]| [-0.18 * x + 0.04 * x * y, 0.15 * y - 0.02 * x * y])
            .collect();
        let (s, c) = 0.73_f64.sin_cos();
        let rotated: Vec<f64> = mesh
            .nodes
            .iter()
            .enumerate()
            .flat_map(|(i, [x, y])| {
                let px = x + u[2 * i];
                let py = y + u[2 * i + 1];
                [c * px - s * py - x + 2.0, s * px + c * py - y - 1.0]
            })
            .collect();
        let energy = p.mean_dilatation_potential_energy(&u, 0.0, cx).unwrap();
        let rotated_energy = p
            .mean_dilatation_potential_energy(&rotated, 0.0, cx)
            .unwrap();
        assert!((energy - rotated_energy).abs() < 2e-12);
        let f = p
            .mean_dilatation_residual_and_tangent(&u, 0.0, cx)
            .unwrap()
            .0;
        let rf = p
            .mean_dilatation_residual_and_tangent(&rotated, 0.0, cx)
            .unwrap()
            .0;
        for i in 0..mesh.nodes.len() {
            assert!((rf[2 * i] - (c * f[2 * i] - s * f[2 * i + 1])).abs() < 2e-11);
            assert!((rf[2 * i + 1] - (s * f[2 * i] + c * f[2 * i + 1])).abs() < 2e-11);
        }
        let mut folded = vec![0.0; 8];
        folded[6] = -1.5; // mean area remains positive; upper corner is reversed.
        assert!(matches!(
            p.mean_dilatation_potential_energy(&folded, 0.0, cx),
            Err(SolidError::MaterialRefused { .. })
        ));
        assert!(matches!(
            p.mean_dilatation_states(&folded, cx),
            Err(SolidError::MaterialRefused { .. })
        ));
        assert!(matches!(
            p.mean_dilatation_residual_and_tangent(&u[..6], 0.0, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let triangles = Mesh2::triangles(1.0, 1.0, 1, 1);
        assert!(matches!(
            problem(&triangles, &law).solve_mean_dilatation(cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let limited = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 40.0,
            },
            1.2,
        )
        .unwrap();
        let shear: Vec<f64> = mesh.nodes.iter().flat_map(|[_, y]| [*y, 0.0]).collect();
        assert!(matches!(
            problem(&mesh, &limited).mean_dilatation_states(&shear, cx),
            Err(SolidError::MaterialRefused { .. })
        ));
        let cancel_during_load = |_: f64, _: f64| {
            gate.request();
            [0.0, 0.0]
        };
        let mut cancelled_problem = problem(&mesh, &law);
        cancelled_problem.traction = vec![(Patch::Right, &cancel_during_load)];
        assert!(matches!(
            cancelled_problem.mean_dilatation_residual_and_tangent(&u, 1.0, cx),
            Err(SolidError::Cancelled)
        ));
        assert_eq!(
            p.solve_mean_dilatation(cx).unwrap_err(),
            SolidError::Cancelled
        );
    });
}
