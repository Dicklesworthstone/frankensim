//! Native reduced-plane paths retain support coupling and physical work.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_solid::continuation::{ArcSettings, PathState};
use fs_solid::hyper2d::HyperelasticMaterials;
use fs_solid::hyper2d::continuation::{PlanePathFormulation as Form, PlanePathScales};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};

const H: f64 = 0.02;
const A: f64 = 0.2;
const TRACTION: f64 = 0.1;
const SCALES: PlanePathScales = PlanePathScales {
    displacement_m: 0.1,
    force_n: 0.04,
};

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 743,
                kernel_id: 319,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate);
    });
}

fn law() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 40.0,
        },
        4.0,
    )
    .unwrap()
}
fn ramp(x: f64, y: f64) -> [f64; 2] {
    [A * x, A * y]
}
fn traction(_: f64, _: f64) -> [f64; 2] {
    [TRACTION, 0.0]
}
fn problem<'a>(mesh: &'a Mesh2, material: &'a dyn HyperelasticMaterials) -> HyperProblem<'a> {
    HyperProblem {
        mesh,
        material,
        dirichlet: vec![
            (Patch::Left, &ramp),
            (Patch::Right, &ramp),
            (Patch::Bottom, &ramp),
            (Patch::Top, &ramp),
        ],
        // Deliberately applied only to prescribed rows: the free interior is
        // driven solely by K_fc*g, not by a nonzero free dead-load vector.
        traction: vec![(Patch::Right, &traction)],
        settings: NewtonSettings::default(),
    }
}
fn settings() -> ArcSettings {
    ArcSettings {
        ds: 0.04,
        ds_max: 0.08,
        tol: 1e-9,
        ..ArcSettings::default()
    }
}
fn cases() -> Vec<(Mesh2, Form)> {
    let local = PlaneStressSettings::default();
    let mut cases = vec![];
    for form in [
        Form::PlaneStress(local),
        Form::IncompressiblePlaneStress(local),
    ] {
        cases.push((Mesh2::triangles(1.0, 1.0, 2, 2), form));
        cases.push((Mesh2::quads(1.0, 1.0, 2, 2), form));
    }
    cases.push((Mesh2::quads(1.0, 1.0, 2, 2), Form::MeanDilatation));
    cases
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn same_state(a: &PathState, b: &PathState) {
    assert_eq!(a.u, b.u);
    assert_eq!(a.lambda, b.lambda);
    assert_eq!(a.tangent, b.tangent);
    assert_eq!(a.ds, b.ds);
    assert_eq!(a.step, b.step);
    assert_eq!(a.events, b.events);
    assert_eq!(a.trace, b.trace);
    assert_eq!(a.last_definite, b.last_definite);
    assert_eq!(a.pending_switch, b.pending_switch);
}

/// Independent affine biaxial NH stored-energy and nominal-stress formulas.
/// Finite plane stress solves only its scalar Pzz=0 reference equation here.
fn oracle(form: Form, stretch: f64) -> (f64, f64) {
    let (mu, lame) = (2.0, 40.0);
    match form {
        Form::PlaneStress(_) => {
            let (mut lo, mut hi) = (1.0 / (stretch * stretch), 1.0);
            for _ in 0..80 {
                let t = 0.5 * (lo + hi);
                if mu * (t * t - 1.0) + lame * (stretch * stretch * t).ln() > 0.0 {
                    hi = t;
                } else {
                    lo = t;
                }
            }
            let t = 0.5 * (lo + hi);
            let log_j = (stretch * stretch * t).ln();
            (
                0.5 * mu * (2.0 * stretch * stretch + t * t - 3.0) - mu * log_j
                    + 0.5 * lame * log_j * log_j,
                mu * (stretch - t * t / stretch),
            )
        }
        Form::IncompressiblePlaneStress(_) => (
            0.5 * mu * (2.0 * stretch * stretch + stretch.powi(-4) - 3.0),
            mu * (stretch - stretch.powi(-5)),
        ),
        Form::MeanDilatation => (
            mu * (stretch * stretch - 1.0) - 2.0 * mu * stretch.ln()
                + 2.0 * lame * stretch.ln().powi(2),
            mu * stretch + (2.0 * lame * stretch.ln() - mu) / stretch,
        ),
    }
}

#[test]
fn support_driven_native_paths_match_affine_energy_reactions_and_exact_resume() {
    with_cx(|cx, _| {
        let material = law();
        for (mesh, form) in cases() {
            let p = problem(&mesh, &material);
            let owner = p.prepare_continuation(form, H, SCALES, 2, cx).unwrap();
            assert_eq!(owner.free_dofs(), &[8, 9]);
            assert!(
                owner
                    .free_dofs()
                    .iter()
                    .all(|&dof| owner.reference_force_n()[dof] == 0.0)
            );
            let control = settings();
            let initial = owner.start(&control, cx).unwrap();
            let straight = owner.advance(&initial, &control, 14, cx).unwrap();
            let half = owner.advance(&initial, &control, 7, cx).unwrap();
            let resumed = owner.advance(&half, &control, 7, cx).unwrap();
            same_state(&straight, &resumed);
            assert!(straight.lambda > 0.4);
            let response = owner.response(&straight, cx).unwrap();
            for (actual, &[x, y]) in response.displacement_m.iter().zip(&mesh.nodes) {
                close(actual[0], straight.lambda * A * x, 2e-9);
                close(actual[1], straight.lambda * A * y, 2e-9);
            }
            let (energy, piola) = oracle(form, 1.0 + A * straight.lambda);
            close(response.strain_energy_j, H * energy, 3e-9);
            close(
                response.potential_energy_j,
                H * (energy - TRACTION * A * straight.lambda * straight.lambda),
                3e-9,
            );
            let right = mesh.patch_nodes(Patch::Right);
            let reaction = response
                .reactions_n
                .iter()
                .filter(|(dof, _)| dof % 2 == 0 && right.contains(&(dof / 2)))
                .map(|(_, force)| force)
                .sum::<f64>();
            close(reaction, H * (piola - TRACTION * straight.lambda), 5e-9);
            assert!(
                response
                    .free_force_residual_n
                    .iter()
                    .fold(0.0_f64, |n, value| n.hypot(*value))
                    <= 1.1 * SCALES.force_n * control.tol
            );
            let from_physical = owner
                .start_from(&response.displacement_m, straight.lambda, &control, cx)
                .unwrap();
            close(from_physical.lambda, straight.lambda, 0.0);
            for (a, b) in from_physical.u.iter().zip(&straight.u) {
                close(*a, *b, 2e-15);
            }
        }
    });
}

#[test]
fn complete_reduced_tangent_and_support_parameter_column_are_actual_derivatives() {
    with_cx(|cx, _| {
        let material = law();
        for (mesh, form) in cases() {
            let p = problem(&mesh, &material);
            let owner = p.prepare_continuation(form, H, SCALES, 2, cx).unwrap();
            let (q, lambda, eps) = ([0.64, 0.58], 0.6, 2e-6);
            let tangent = owner.tangent(&q, lambda, cx).unwrap();
            let load = owner.load_vector_at(&q, lambda, cx).unwrap();
            assert!(load.iter().any(|v| v.abs() > 0.01));
            for column in 0..2 {
                let (mut plus, mut minus) = (q, q);
                plus[column] += eps;
                minus[column] -= eps;
                let rp = owner.residual(&plus, lambda, cx).unwrap();
                let rm = owner.residual(&minus, lambda, cx).unwrap();
                for row in 0..2 {
                    close(
                        tangent.get(row, column),
                        (rp[row] - rm[row]) / (2.0 * eps),
                        8e-6,
                    );
                }
            }
            let rp = owner.residual(&q, lambda + eps, cx).unwrap();
            let rm = owner.residual(&q, lambda - eps, cx).unwrap();
            for row in 0..2 {
                close(load[row], -(rp[row] - rm[row]) / (2.0 * eps), 8e-6);
            }
        }
    });
}

#[test]
fn held_reference_laws_retain_native_force_and_energy_and_frozen_history_refuses() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
        let original = law();
        let references: Vec<_> = mesh
            .elems
            .iter()
            .enumerate()
            .map(|(i, _)| {
                ReferenceDeformation::try_new([
                    1.0 + 0.01 * i as f64,
                    0.02,
                    0.03,
                    0.0,
                    0.98,
                    -0.02,
                    0.0,
                    0.0,
                    1.04,
                ])
                .unwrap()
            })
            .collect();
        let laws: Vec<_> = references
            .iter()
            .map(|&g| ReferencedHyperelastic::new(&original, g))
            .collect();
        let p = problem(&mesh, &laws);
        for exact in [false, true] {
            let local = PlaneStressSettings::default();
            let form = if exact {
                Form::IncompressiblePlaneStress(local)
            } else {
                Form::PlaneStress(local)
            };
            let owner = p.prepare_continuation(form, H, SCALES, 2, cx).unwrap();
            // Observation is intentionally off equilibrium; every original law
            // must survive restriction without averaging away the reference field.
            let mut state = PathState::start(2, &settings());
            state.u = vec![0.34, 0.29];
            state.lambda = 0.3;
            let response = owner.response(&state, cx).unwrap();
            let flat = response.displacement_m.as_flattened();
            let (force, _) = if exact {
                p.incompressible_plane_stress_residual_and_tangent(flat, 0.3, local, cx)
            } else {
                p.plane_stress_residual_and_tangent(flat, 0.3, local, cx)
            }
            .unwrap();
            let energy = if exact {
                p.incompressible_plane_stress_potential_energy(flat, 0.3, local, cx)
            } else {
                p.plane_stress_potential_energy(flat, 0.3, local, cx)
            }
            .unwrap();
            close(response.potential_energy_j, H * energy, 2e-12);
            for (i, &dof) in owner.free_dofs().iter().enumerate() {
                close(response.free_force_residual_n[i], H * force[dof], 2e-12);
            }
            for &(dof, reaction) in &response.reactions_n {
                close(reaction, H * force[dof], 2e-12);
            }
        }
        let visco = ViscoHyperelastic::new(
            law(),
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 3.0 },
                tau_s: 0.6,
            }],
        )
        .unwrap();
        let accepted = visco.initial_state();
        let frozen = visco.freeze(&accepted, 0.25).unwrap();
        let p = problem(&mesh, &frozen);
        assert!(
            p.prepare_continuation(
                Form::PlaneStress(PlaneStressSettings::default()),
                H,
                SCALES,
                2,
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn invalid_resume_domain_and_cancellation_preserve_the_accepted_path() {
    with_cx(|cx, gate| {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
        let material = law();
        let p = problem(&mesh, &material);
        let form = Form::PlaneStress(PlaneStressSettings::default());
        assert!(p.prepare_continuation(form, H, SCALES, 1, cx).is_err());
        let owner = p.prepare_continuation(form, H, SCALES, 2, cx).unwrap();
        let control = settings();
        let initial = owner.start(&control, cx).unwrap();
        let accepted = owner.advance(&initial, &control, 3, cx).unwrap();
        let saved = accepted.clone();
        let mut edited = accepted.clone();
        edited.u[0] += 0.01;
        assert!(owner.advance(&edited, &control, 1, cx).is_err());
        let mut physical = owner.response(&accepted, cx).unwrap().displacement_m;
        physical[0][0] = 1e-8;
        assert!(
            owner
                .start_from(&physical, accepted.lambda, &control, cx)
                .is_err()
        );
        // Collapse/fold is rejected by the original native element/domain owner.
        assert!(owner.residual(&[20.0, 20.0], accepted.lambda, cx).is_err());
        gate.request();
        assert!(matches!(
            owner.advance(&accepted, &control, 1, cx),
            Err(SolidError::Cancelled)
        ));
        same_state(&accepted, &saved);
    });
}
