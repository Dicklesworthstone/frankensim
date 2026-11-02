//! G0/G1/G3: actual thin-solid equilibrium and independent work derivatives.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{FibreInvariants, Hyperelastic, HyperelasticModel, Volumetric};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 49,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn nh(mu: f64, lambda: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap()
}

fn problem<'a>(mesh: &'a Mesh2, material: &'a Hyperelastic) -> HyperProblem<'a> {
    HyperProblem {
        mesh,
        material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings {
            tol: 2e-9,
            max_iters: 40,
            load_steps: 5,
            max_backtracks: 30,
        },
    }
}

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * (1.0 + b.abs()), "{a:e} != {b:e}");
}

#[test]
fn g1_dead_traction_solve_recovers_finite_lateral_and_thickness_contraction() {
    with_cx(|cx, _| {
        let (mu, lambda, axial): (f64, f64, f64) = (2.0, 40.0, 1.3);
        // Independent uniaxial NH reference: Pyy=Pzz=0 with equal lateral
        // stretches t gives mu(t²-1)+lambda ln(axial*t²)=0.
        let mut lo = 1.0 / axial.sqrt();
        let mut hi = 1.0;
        for _ in 0..70 {
            let t = 0.5 * (lo + hi);
            if mu * (t * t - 1.0) + lambda * (axial * t * t).ln() > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
        }
        let lateral = 0.5 * (lo + hi);
        let nominal = mu * (axial - lateral * lateral / axial);
        let law = nh(mu, lambda);
        let left = |_: f64, y: f64| [0.0, (lateral - 1.0) * y];
        let pull = |_: f64, _: f64| [nominal, 0.0];
        let local = PlaneStressSettings::default();
        for mesh in [
            Mesh2::triangles(2.0, 1.0, 2, 2),
            Mesh2::quads(2.0, 1.0, 2, 2),
        ] {
            let mut p = problem(&mesh, &law);
            p.dirichlet = vec![(Patch::Left, &left)];
            p.traction = vec![(Patch::Right, &pull)];
            let (u, report) = p.solve_plane_stress(local, cx).unwrap();
            assert_eq!(report.histories.len(), p.settings.load_steps);
            assert!(
                report
                    .histories
                    .iter()
                    .all(|h| h.last().unwrap() < &p.settings.tol)
            );
            for (v, [x, y]) in u.iter().zip(&mesh.nodes) {
                close(v[0], (axial - 1.0) * x, 2e-8);
                close(v[1], (lateral - 1.0) * y, 2e-8);
            }
            let flat: Vec<f64> = u.iter().flatten().copied().collect();
            for state in p.plane_stress_states(&flat, local, cx).unwrap() {
                close(state.deformation_gradient[8], lateral, 2e-8);
                assert!(state.deformation_gradient[8] < 0.99);
                close(state.piola[0], nominal, 2e-8);
                close(state.piola[3], 0.0, 2e-8);
                assert!(state.relative_transverse_residual <= local.relative_tolerance);
            }
            let (residual, _) = p
                .plane_stress_residual_and_tangent(&flat, 1.0, local, cx)
                .unwrap();
            let reaction: f64 = mesh
                .patch_nodes(Patch::Left)
                .iter()
                .map(|&i| residual[2 * i])
                .sum();
            close(reaction, -nominal, 2e-8);
        }
    });
}

#[test]
fn g0_skew_elements_energy_force_and_tangent_share_the_relaxed_material() {
    with_cx(|cx, _| {
        let local = PlaneStressSettings::default();
        let fibre = Hyperelastic::new(
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 1.0,
                k1: 1.8,
                k2: 1.1,
                dispersion: 0.05,
                fibres: vec![[0.8, 0.0, 0.6]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 20.0 },
            },
            4.0,
        )
        .unwrap();
        for law in [nh(2.0, 40.0), fibre] {
            for mesh in [
                Mesh2::triangles(1.0, 0.7, 1, 1),
                Mesh2::mapped_quads(2, 1, &|x, y| [x + 0.15 * x * y + 0.1 * y, y + 0.1 * x]),
            ] {
                let traction = |x: f64, y: f64| [0.1 + x * 0.04, -0.03 * y];
                let mut p = problem(&mesh, &law);
                p.traction = vec![(Patch::Right, &traction)];
                let u: Vec<f64> = mesh
                    .nodes
                    .iter()
                    .flat_map(|[x, y]| [0.14 * x + 0.03 * x * y, -0.04 * y + 0.02 * x * y])
                    .collect();
                let (force, tangent) = p
                    .plane_stress_residual_and_tangent(&u, 0.7, local, cx)
                    .unwrap();
                let h = 2e-6;
                for col in 0..u.len() {
                    let mut plus = u.clone();
                    let mut minus = u.clone();
                    plus[col] += h;
                    minus[col] -= h;
                    let ep = p
                        .plane_stress_potential_energy(&plus, 0.7, local, cx)
                        .unwrap();
                    let em = p
                        .plane_stress_potential_energy(&minus, 0.7, local, cx)
                        .unwrap();
                    close((ep - em) / (2.0 * h), force[col], 3e-7);
                    let rp = p
                        .plane_stress_residual_and_tangent(&plus, 0.7, local, cx)
                        .unwrap()
                        .0;
                    let rm = p
                        .plane_stress_residual_and_tangent(&minus, 0.7, local, cx)
                        .unwrap()
                        .0;
                    for row in 0..u.len() {
                        close((rp[row] - rm[row]) / (2.0 * h), tangent.get(row, col), 5e-6);
                        close(tangent.get(row, col), tangent.get(col, row), 5e-11);
                    }
                }
            }
        }
    });
}

#[test]
fn g3_rigid_motion_preserves_energy_and_rotates_internal_forces() {
    with_cx(|cx, _| {
        let local = PlaneStressSettings::default();
        let mesh = Mesh2::mapped_quads(2, 1, &|x, y| [x + 0.12 * y, y + 0.07 * x * y]);
        let law = nh(2.0, 40.0);
        let p = problem(&mesh, &law);
        let u: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|[x, y]| [0.18 * x + 0.02 * x * y, -0.04 * y])
            .collect();
        let (s, c) = 0.64_f64.sin_cos();
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
        close(
            p.plane_stress_potential_energy(&u, 0.0, local, cx).unwrap(),
            p.plane_stress_potential_energy(&rotated, 0.0, local, cx)
                .unwrap(),
            2e-11,
        );
        let (force, _) = p
            .plane_stress_residual_and_tangent(&u, 0.0, local, cx)
            .unwrap();
        let (rf, _) = p
            .plane_stress_residual_and_tangent(&rotated, 0.0, local, cx)
            .unwrap();
        for i in 0..mesh.nodes.len() {
            close(rf[2 * i], c * force[2 * i] - s * force[2 * i + 1], 2e-10);
            close(
                rf[2 * i + 1],
                s * force[2 * i] + c * force[2 * i + 1],
                2e-10,
            );
        }
    });
}

#[test]
fn g0_rejects_folded_domains_bad_controls_and_cancellation_without_partial_output() {
    with_cx(|cx, gate| {
        let local = PlaneStressSettings::default();
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let law = nh(2.0, 40.0);
        let p = problem(&mesh, &law);
        let mut folded = vec![0.0; 8];
        folded[6] = -1.5;
        assert!(matches!(
            p.plane_stress_states(&folded, local, cx),
            Err(SolidError::MaterialRefused { .. })
        ));
        assert!(matches!(
            p.plane_stress_states(&[0.0; 6], local, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let mut limited = law.clone();
        limited.stretch_limit = 1.1;
        let stretch: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|[x, _]| [0.3 * x, 0.0])
            .collect();
        assert!(matches!(
            problem(&mesh, &limited).plane_stress_states(&stretch, local, cx),
            Err(SolidError::MaterialRefused { .. })
        ));
        let mut bad = problem(&mesh, &law);
        bad.settings.load_steps = 0;
        assert!(matches!(
            bad.solve_plane_stress(local, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let mut invalid_mesh = mesh.clone();
        invalid_mesh.elems[0][0] = usize::MAX;
        assert!(matches!(
            problem(&invalid_mesh, &law).solve_plane_stress(local, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let cancel_during_load = |_: f64, _: f64| {
            gate.request();
            [0.0, 0.0]
        };
        let mut cancelled = problem(&mesh, &law);
        cancelled.traction = vec![(Patch::Right, &cancel_during_load)];
        assert!(matches!(
            cancelled.plane_stress_potential_energy(&[0.0; 8], 1.0, local, cx),
            Err(SolidError::Cancelled)
        ));
        assert!(matches!(
            p.solve_plane_stress(local, cx),
            Err(SolidError::Cancelled)
        ));
    });
}
