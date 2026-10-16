//! Three-dimensional finite-strain hyperelasticity on P1 tetrahedra.
//!
//! G0: the assembled tangent is the exact derivative of the internal force
//! (central-difference directional check); at the reference state it equals
//! the independently written small-strain `linear3` stiffness.
//! G1: the affine patch test is reproduced exactly at finite strain
//! (interior nodes and stored energy).
//! G2: compressible Neo-Hookean uniaxial stretch to lambda = 1.6 against the
//! closed-form lateral stretch and nominal stress.
//! G3: a rigid rotation of the boundary stores no energy and needs no force.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::linear3::{TetLinearElasticProblem, TetMaterialField};
use fs_solid::{
    HyperTetError, HyperTetProblem, HyperTetSettings, TetAssemblyBudget, TetElasticMaterial,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Conforming box mesh: every cube split into six tetrahedra around its
/// main diagonal (Kuhn), positively oriented.
fn box_mesh(n: [usize; 3], size: [f64; 3]) -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let id = |i: usize, j: usize, k: usize| (k * (n[1] + 1) + j) * (n[0] + 1) + i;
    let mut nodes = Vec::new();
    for k in 0..=n[2] {
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                nodes.push([
                    size[0] * i as f64 / n[0] as f64,
                    size[1] * j as f64 / n[1] as f64,
                    size[2] * k as f64 / n[2] as f64,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                for perm in [
                    [0, 1, 2],
                    [0, 2, 1],
                    [1, 0, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [2, 1, 0],
                ] {
                    let mut c = [i, j, k];
                    let mut tet = [id(c[0], c[1], c[2]); 4];
                    for (step, axis) in perm.iter().enumerate() {
                        c[*axis] += 1;
                        tet[step + 1] = id(c[0], c[1], c[2]);
                    }
                    let x = tet.map(|v| nodes[v]);
                    let e = |a: usize| [0, 1, 2].map(|d| x[a][d] - x[0][d]);
                    let (a, b, cc) = (e(1), e(2), e(3));
                    let det = (a[1] * b[2] - a[2] * b[1]) * cc[0]
                        + (a[2] * b[0] - a[0] * b[2]) * cc[1]
                        + (a[0] * b[1] - a[1] * b[0]) * cc[2];
                    if det < 0.0 {
                        tet.swap(2, 3);
                    }
                    tets.push(tet);
                }
            }
        }
    }
    (nodes, tets)
}

fn neo_hookean(mu: f64, lambda: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 5.0).unwrap()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    material: &'a Hyperelastic,
    prescribed: &'a [(usize, f64)],
    forces: &'a [(usize, f64)],
    load_steps: usize,
) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        material,
        prescribed_m: prescribed,
        nodal_forces_n: forces,
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings {
            load_steps,
            ..HyperTetSettings::default()
        },
    }
}

fn on_boundary(x: [f64; 3], size: [f64; 3]) -> bool {
    (0..3).any(|d| x[d].abs() < 1e-12 || (x[d] - size[d]).abs() < 1e-12)
}

#[test]
fn affine_patch_test_is_exact_at_finite_strain() {
    let size = [1.0, 1.0, 1.0];
    let (nodes, tets) = box_mesh([3, 3, 3], size);
    let material = neo_hookean(1.0, 2.0);
    let f0 = [1.2, 0.1, 0.0, 0.0, 0.9, 0.05, 0.02, 0.0, 1.1];
    let affine = |x: [f64; 3], i: usize| {
        (0..3)
            .map(|j| (f0[3 * i + j] - f64::from(u8::from(i == j))) * x[j])
            .sum::<f64>()
    };
    let mut prescribed = Vec::new();
    for (node, &x) in nodes.iter().enumerate() {
        if on_boundary(x, size) {
            for i in 0..3 {
                prescribed.push((3 * node + i, affine(x, i)));
            }
        }
    }
    let p = problem(&nodes, &tets, &material, &prescribed, &[], 3);
    let solution = with_cx(|cx| p.solve(cx)).unwrap();
    let mut worst = 0.0f64;
    for (node, &x) in nodes.iter().enumerate() {
        for i in 0..3 {
            worst = worst.max((solution.displacement_m[node][i] - affine(x, i)).abs());
        }
    }
    let exact_energy = material.energy(&f0);
    println!(
        "patch: worst interior error {worst:e}, energy {} vs {exact_energy}",
        solution.strain_energy_j
    );
    assert!(worst < 1e-10, "affine field not reproduced: {worst}");
    assert!((solution.strain_energy_j - exact_energy).abs() < 1e-12 * exact_energy);
    // Global equilibrium: the support forces balance.
    for i in 0..3 {
        let total: f64 = solution
            .reactions_n
            .iter()
            .filter(|(dof, _)| dof % 3 == i)
            .map(|(_, r)| r)
            .sum();
        assert!(
            total.abs() < 1e-10,
            "component {i} reaction resultant {total}"
        );
    }
}

/// Lateral stretch `t` of compressible Neo-Hookean uniaxial stress:
/// `mu (t^2 - 1) + lambda ln(stretch t^2) = 0`.
fn lateral_stretch(mu: f64, lambda: f64, stretch: f64) -> f64 {
    let mut t = 1.0f64;
    for _ in 0..60 {
        let g = mu * (t * t - 1.0) + lambda * (stretch * t * t).ln();
        let dg = 2.0 * mu * t + 2.0 * lambda / t;
        t -= g / dg;
    }
    t
}

#[test]
fn neo_hookean_uniaxial_stretch_matches_closed_form() {
    let size = [2.0, 1.0, 1.0];
    let (nodes, tets) = box_mesh([4, 2, 2], size);
    let (mu, lambda) = (1.0, 1.5);
    let material = neo_hookean(mu, lambda);
    let stretch = 1.6;
    let mut prescribed = Vec::new();
    for (node, &x) in nodes.iter().enumerate() {
        if x[0].abs() < 1e-12 {
            prescribed.push((3 * node, 0.0));
        }
        if (x[0] - size[0]).abs() < 1e-12 {
            prescribed.push((3 * node, (stretch - 1.0) * size[0]));
        }
        if x[1].abs() < 1e-12 {
            prescribed.push((3 * node + 1, 0.0));
        }
        if x[2].abs() < 1e-12 {
            prescribed.push((3 * node + 2, 0.0));
        }
    }
    let p = problem(&nodes, &tets, &material, &prescribed, &[], 4);
    let solution = with_cx(|cx| p.solve(cx)).unwrap();
    let t = lateral_stretch(mu, lambda, stretch);
    let j = stretch * t * t;
    let nominal = mu * (stretch - 1.0 / stretch) + lambda * j.ln() / stretch;
    let axial: f64 = prescribed
        .iter()
        .zip(&solution.reactions_n)
        .filter(|((dof, value), _)| dof % 3 == 0 && *value > 0.0)
        .map(|(_, (_, r))| r)
        .sum();
    let corner = nodes
        .iter()
        .position(|x| x.iter().zip(&size).all(|(a, b)| (a - b).abs() < 1e-12))
        .unwrap();
    println!(
        "uniaxial lambda={stretch}: lateral {:.12} vs {t:.12}, force {axial:.12} vs {nominal:.12}, steps {:?}",
        1.0 + solution.displacement_m[corner][1],
        solution
            .steps
            .iter()
            .map(|s| s.residual_history.len())
            .collect::<Vec<_>>()
    );
    assert!((1.0 + solution.displacement_m[corner][1] - t).abs() < 1e-9);
    assert!((1.0 + solution.displacement_m[corner][2] - t).abs() < 1e-9);
    assert!(
        (axial - nominal).abs() < 1e-9 * nominal,
        "{axial} vs {nominal}"
    );
    // Newton's local quadratic rate: few iterations per load step.
    for step in &solution.steps {
        assert!(step.residual_history.len() <= 8, "{step:?}");
        assert_eq!(step.indefinite_fallbacks, 0);
    }
}

#[test]
fn tangent_is_the_derivative_of_the_internal_force() {
    let (nodes, tets) = box_mesh([2, 2, 2], [1.0, 1.0, 1.0]);
    for material in [
        neo_hookean(1.0, 3.0),
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 0.4,
                c01: 0.1,
                kappa: 5.0,
            },
            5.0,
        )
        .unwrap(),
    ] {
        let p = problem(&nodes, &tets, &material, &[], &[], 1);
        // A smooth, finite deformation and a pseudo-random direction.
        let u: Vec<f64> = nodes
            .iter()
            .flat_map(|x| {
                [
                    0.2 * x[1] * x[2],
                    0.1 * x[0] * x[0] - 0.05 * x[2],
                    0.15 * x[0] * x[1],
                ]
            })
            .collect();
        let d: Vec<f64> = (0..u.len())
            .map(|i| ((i * 7919 % 101) as f64 / 101.0) - 0.5)
            .collect();
        let (k, plus, minus) = with_cx(|cx| {
            let eps = 1e-6;
            let up: Vec<f64> = u.iter().zip(&d).map(|(a, b)| a + eps * b).collect();
            let um: Vec<f64> = u.iter().zip(&d).map(|(a, b)| a - eps * b).collect();
            (
                p.tangent_matrix(&u, cx).unwrap(),
                p.internal_force(&up, cx).unwrap(),
                p.internal_force(&um, cx).unwrap(),
            )
        });
        let mut kd = vec![0.0; u.len()];
        k.spmv(&d, &mut kd);
        let fd: Vec<f64> = plus
            .iter()
            .zip(&minus)
            .map(|(a, b)| (a - b) / 2e-6)
            .collect();
        let err = kd
            .iter()
            .zip(&fd)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        let scale = kd.iter().map(|a| a * a).sum::<f64>().sqrt();
        println!(
            "{:?}: |K d - FD| / |K d| = {:e}",
            material.model,
            err / scale
        );
        assert!(err / scale < 1e-6);
    }
}

#[test]
fn reference_tangent_equals_the_small_strain_stiffness() {
    let (nodes, tets) = box_mesh([2, 1, 2], [1.0, 0.5, 1.0]);
    let (mu, lambda) = (0.8, 1.7);
    let material = neo_hookean(mu, lambda);
    let young = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
    let poisson = lambda / (2.0 * (lambda + mu));
    let linear = TetElasticMaterial::try_new(1.0, young, poisson, [7; 32]).unwrap();
    let (k_hyper, k_linear) = with_cx(|cx| {
        let hyper = problem(&nodes, &tets, &material, &[], &[], 1);
        let reference = TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            materials: TetMaterialField::Uniform(&linear),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
        };
        (
            hyper
                .tangent_matrix(&vec![0.0; 3 * nodes.len()], cx)
                .unwrap(),
            reference.assemble(cx).unwrap().stiffness,
        )
    });
    let (mut worst, mut scale) = (0.0f64, 0.0f64);
    for r in 0..k_linear.nrows() {
        for c in 0..k_linear.ncols() {
            worst = worst.max((k_hyper.get(r, c) - k_linear.get(r, c)).abs());
            scale = scale.max(k_linear.get(r, c).abs());
        }
    }
    println!("reference tangent vs linear3: max |dK| = {worst:e} (scale {scale:e})");
    assert!(worst < 1e-12 * scale);
}

#[test]
fn rigid_rotation_stores_no_energy_and_needs_no_force() {
    let size = [1.0, 1.0, 1.0];
    let (nodes, tets) = box_mesh([2, 2, 2], size);
    let material = neo_hookean(1.0, 1.0);
    let (c, s) = (0.6f64, 0.8f64); // 53 degrees about z
    let rotation = [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0];
    let mut prescribed = Vec::new();
    for (node, &x) in nodes.iter().enumerate() {
        if on_boundary(x, size) {
            for i in 0..3 {
                let rx: f64 = (0..3).map(|j| rotation[3 * i + j] * x[j]).sum();
                prescribed.push((3 * node + i, rx - x[i]));
            }
        }
    }
    let p = problem(&nodes, &tets, &material, &prescribed, &[], 4);
    let solution = with_cx(|cx| p.solve(cx)).unwrap();
    let largest = solution
        .reactions_n
        .iter()
        .map(|(_, r)| r.abs())
        .fold(0.0, f64::max);
    println!(
        "rotation: energy {:e}, largest reaction {largest:e}",
        solution.strain_energy_j
    );
    assert!(solution.strain_energy_j.abs() < 1e-12);
    assert!(largest < 1e-10);
    assert!((solution.min_det_f - 1.0).abs() < 1e-12);
}

#[test]
fn refusals_are_structured() {
    let (nodes, tets) = box_mesh([1, 1, 1], [1.0, 1.0, 1.0]);
    let material = neo_hookean(1.0, 1.0);
    // Inverted element.
    let mut flipped = tets.clone();
    flipped[0].swap(0, 1);
    let err =
        with_cx(|cx| problem(&nodes, &flipped, &material, &[], &[], 1).solve(cx)).unwrap_err();
    assert!(
        matches!(err, HyperTetError::DegenerateElement { element: 0, .. }),
        "{err}"
    );
    // Repeated prescribed DOF.
    let twice = [(0, 0.0), (0, 0.1)];
    let err =
        with_cx(|cx| problem(&nodes, &tets, &material, &twice, &[], 1).solve(cx)).unwrap_err();
    assert!(matches!(err, HyperTetError::InvalidInput { .. }), "{err}");
    // A tripped context publishes nothing.
    let gate = CancelGate::new();
    gate.request();
    let err = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0)];
        problem(&nodes, &tets, &material, &fixed, &[(21, 0.1)], 1).solve(&cx)
    });
    assert_eq!(err.unwrap_err(), HyperTetError::Cancelled);
}
