//! G1/G3: finite elastic equilibrium must survive common material/load scaling.
//! The references below are independent homogeneous Neo-Hookean solutions.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 49,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

#[derive(Clone, Copy, Debug)]
enum Formulation {
    PlaneStrain,
    MeanDilatation,
    PlaneStress,
    IncompressiblePlaneStress,
}

fn exact(formulation: Formulation) -> (f64, f64, f64, f64) {
    let (mu, lambda, axial): (f64, f64, f64) = (2.0, 5.0, 1.25);
    match formulation {
        Formulation::PlaneStrain | Formulation::MeanDilatation => {
            let lateral = 0.9;
            let logarithm = (axial * lateral).ln();
            (
                axial,
                lateral,
                mu * (axial - axial.recip()) + lambda * logarithm / axial,
                mu * (lateral - lateral.recip()) + lambda * logarithm / lateral,
            )
        }
        Formulation::PlaneStress => {
            // Equal lateral/thickness stretches satisfy Pyy=Pzz=0.
            let mut lo = axial.sqrt().recip();
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
            (
                axial,
                lateral,
                mu * (axial - lateral * lateral / axial),
                0.0,
            )
        }
        Formulation::IncompressiblePlaneStress => (
            axial,
            axial.sqrt().recip(),
            mu * (axial - axial.powi(-2)),
            0.0,
        ),
    }
}

fn check_equilibrium(formulation: Formulation, stress_scale: f64, length: f64, cx: &Cx<'_>) {
    let mesh = Mesh2::quads(2.0 * length, length, 2, 1);
    let material = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0 * stress_scale,
            lambda: 5.0 * stress_scale,
        },
        4.0,
    )
    .unwrap();
    let (axial, lateral, px, py) = exact(formulation);
    let left = |_: f64, y: f64| [0.0, (lateral - 1.0) * y];
    let right = |_: f64, _: f64| [stress_scale * px, 0.0];
    let top = |_: f64, _: f64| [0.0, stress_scale * py];
    let bottom = |_: f64, _: f64| [0.0, -stress_scale * py];
    let problem = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![(Patch::Left, &left)],
        traction: vec![
            (Patch::Right, &right),
            (Patch::Top, &top),
            (Patch::Bottom, &bottom),
        ],
        settings: NewtonSettings {
            tol: 2e-9,
            max_iters: 40,
            load_steps: 1,
            max_backtracks: 30,
        },
    };
    let local = PlaneStressSettings::default();
    let result = match formulation {
        Formulation::PlaneStrain => problem.solve(),
        Formulation::MeanDilatation => problem.solve_mean_dilatation(cx),
        Formulation::PlaneStress => problem.solve_plane_stress(local, cx),
        Formulation::IncompressiblePlaneStress => {
            problem.solve_incompressible_plane_stress(local, cx)
        }
    };
    let (motion, report) = result.unwrap_or_else(|error| {
        panic!("{formulation:?}, stress={stress_scale:e}, length={length:e}: {error}")
    });
    assert_eq!(report.histories.len(), 1);
    assert_eq!(report.force_scales.len(), 1);
    assert!(report.force_scales[0] > 0.0);
    assert!(
        report.histories[0].len() > 1,
        "a loaded body must actually move"
    );
    assert!(*report.histories[0].last().unwrap() < problem.settings.tol);
    for (actual, point) in motion.iter().zip(&mesh.nodes) {
        let expected = [(axial - 1.0) * point[0], (lateral - 1.0) * point[1]];
        for i in 0..2 {
            assert!(
                (actual[i] - expected[i]).abs() / length < 2e-8,
                "{formulation:?}, stress={stress_scale:e}, length={length:e}: {actual:?} != {expected:?}"
            );
        }
    }
    // Check physical reactions/free forces independently of the reported norm.
    let flat: Vec<f64> = motion.iter().flatten().copied().collect();
    let residual = match formulation {
        Formulation::PlaneStrain => problem.unconstrained_residual(&flat, 1.0),
        Formulation::MeanDilatation => problem
            .mean_dilatation_residual_and_tangent(&flat, 1.0, cx)
            .map(|(r, _)| r),
        Formulation::PlaneStress => problem
            .plane_stress_residual_and_tangent(&flat, 1.0, local, cx)
            .map(|(r, _)| r),
        Formulation::IncompressiblePlaneStress => problem
            .incompressible_plane_stress_residual_and_tangent(&flat, 1.0, local, cx)
            .map(|(r, _)| r),
    }
    .unwrap();
    let supports = mesh.patch_nodes(Patch::Left);
    let reaction: f64 = supports.iter().map(|&node| residual[2 * node]).sum();
    assert!((reaction / (stress_scale * length) + px).abs() < 3e-8);
    for (node, force) in residual.chunks_exact(2).enumerate() {
        if !supports.contains(&node) {
            for value in force {
                assert!(value.abs() / (stress_scale * length) < 3e-8);
            }
        }
    }
}

#[test]
fn g3_all_plane_formulations_preserve_equilibrium_under_material_and_load_rescaling() {
    with_cx(|cx| {
        for formulation in [
            Formulation::PlaneStrain,
            Formulation::MeanDilatation,
            Formulation::PlaneStress,
            Formulation::IncompressiblePlaneStress,
        ] {
            for stress_scale in [1e-12, 1.0, 1e12] {
                check_equilibrium(formulation, stress_scale, 1.0, cx);
            }
        }
    });
}

#[test]
fn g3_loaded_elasticity_preserves_strain_on_small_and_large_specimens() {
    with_cx(|cx| {
        for length in [1e-8, 1.0, 1e8] {
            check_equilibrium(Formulation::PlaneStrain, 1e-6, length, cx);
            check_equilibrium(Formulation::PlaneStress, 1e-6, length, cx);
        }
    });
}

#[test]
fn g0_zero_load_and_invalid_controls_cannot_manufacture_convergence() {
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let material = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 5.0,
        },
        4.0,
    )
    .unwrap();
    let clamp = |_: f64, _: f64| [0.0, 0.0];
    let mut problem = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![(Patch::Left, &clamp)],
        traction: vec![],
        settings: NewtonSettings::default(),
    };
    let (motion, report) = problem.solve().unwrap();
    assert!(motion.iter().flatten().all(|&value| value == 0.0));
    assert_eq!(report.histories, vec![vec![0.0]]);
    assert_eq!(report.force_scales, vec![0.0]);
    for tolerance in [0.0, -1.0, 1.0, f64::INFINITY, f64::NAN] {
        problem.settings.tol = tolerance;
        assert!(matches!(
            problem.solve(),
            Err(SolidError::InvalidInput { .. })
        ));
    }
    problem.settings = NewtonSettings {
        load_steps: 0,
        ..NewtonSettings::default()
    };
    assert!(matches!(
        problem.solve(),
        Err(SolidError::InvalidInput { .. })
    ));
}
