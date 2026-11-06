//! G0/G1/G3: the native plane continuum retains prescribed strain and history.

use super::*;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_solid::hyper2d::HyperelasticMaterials;
use fs_sparse::Csr;

#[derive(Clone, Copy, Debug)]
enum Form {
    Strain,
    Stress,
    Exact,
}

const FORMS: [Form; 3] = [Form::Strain, Form::Stress, Form::Exact];

fn sheet<'a>(mesh: &'a Mesh2, material: &'a dyn HyperelasticMaterials) -> HyperProblem<'a> {
    HyperProblem {
        mesh,
        material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings {
            tol: 2e-9,
            max_iters: 40,
            load_steps: 1,
            max_backtracks: 30,
        },
    }
}

fn probe(p: &HyperProblem<'_>, u: &[f64], form: Form, cx: &Cx<'_>) -> (Vec<f64>, Csr) {
    let settings = PlaneStressSettings::default();
    match form {
        Form::Strain => p.residual_and_tangent(u, 0.0),
        Form::Stress => p.plane_stress_residual_and_tangent(u, 0.0, settings, cx),
        Form::Exact => p.incompressible_plane_stress_residual_and_tangent(u, 0.0, settings, cx),
    }
    .unwrap()
}

fn potential(p: &HyperProblem<'_>, u: &[f64], form: Form, cx: &Cx<'_>) -> f64 {
    let settings = PlaneStressSettings::default();
    match form {
        Form::Strain => p.potential_energy(u, 0.0).unwrap(),
        Form::Stress => p
            .plane_stress_potential_energy(u, 0.0, settings, cx)
            .unwrap(),
        Form::Exact => p
            .incompressible_plane_stress_potential_energy(u, 0.0, settings, cx)
            .unwrap(),
    }
}

fn solve(p: &HyperProblem<'_>, form: Form, cx: &Cx<'_>) -> Vec<[f64; 2]> {
    let settings = PlaneStressSettings::default();
    let (u, report) = match form {
        Form::Strain => p.solve(),
        Form::Stress => p.solve_plane_stress(settings, cx),
        Form::Exact => p.solve_incompressible_plane_stress(settings, cx),
    }
    .unwrap();
    assert!(
        report
            .histories
            .iter()
            .all(|h| *h.last().unwrap() <= p.settings.tol)
    );
    u
}

fn layer(y: f64) -> usize {
    usize::from(y > 0.5)
}

fn free_motion(x: f64, y: f64) -> [f64; 2] {
    let low = y.min(0.5);
    let high = (y - 0.5).max(0.0);
    [
        0.08 * x + 0.06 * low - 0.03 * high,
        0.04 * low - 0.02 * high,
    ]
}

#[test]
fn g1_native_equilibrium_recovers_compatible_heterogeneous_reference_strain() {
    with_cx(|cx, _| {
        let original = nh(2.0, 40.0);
        let settings = PlaneStressSettings::default();
        for mesh in [
            Mesh2::triangles(1.0, 1.0, 2, 2),
            Mesh2::quads(1.0, 1.0, 2, 2),
        ] {
            for form in FORMS {
                let references: Vec<_> = mesh
                    .elems
                    .iter()
                    .map(|nodes| {
                        let y = nodes.iter().map(|&n| mesh.nodes[n][1]).sum::<f64>()
                            / nodes.len() as f64;
                        let upper = layer(y) == 1;
                        let (shear, stretch) = if upper { (-0.03, 0.98) } else { (0.06, 1.04) };
                        let (xz, yz, zz) = if matches!(form, Form::Strain) {
                            (0.0, 0.0, 1.0)
                        } else if upper {
                            (-0.04, 0.03, 0.97)
                        } else {
                            (0.05, -0.02, 1.06)
                        };
                        ReferenceDeformation::try_new([
                            1.08, shear, xz, 0.0, stretch, yz, 0.0, 0.0, zz,
                        ])
                        .unwrap()
                    })
                    .collect();
                let laws: Vec<_> = references
                    .iter()
                    .map(|&reference| ReferencedHyperelastic::new(&original, reference))
                    .collect();
                let left = |x: f64, y: f64| free_motion(x, y);
                let mut p = sheet(&mesh, &laws);
                p.dirichlet = vec![(Patch::Left, &left)];
                let u = solve(&p, form, cx);
                for (actual, &[x, y]) in u.iter().zip(&mesh.nodes) {
                    let expected = free_motion(x, y);
                    close(actual[0], expected[0], 3e-8);
                    close(actual[1], expected[1], 3e-8);
                }
                let flat: Vec<_> = u.iter().flatten().copied().collect();
                let residual = probe(&p, &flat, form, cx).0;
                assert!(
                    residual.iter().all(|f| f.abs() < 2e-8),
                    "{form:?}: {residual:?}"
                );
                close(potential(&p, &flat, form, cx), 0.0, 2e-9);
                if !matches!(form, Form::Strain) {
                    let gradients: Vec<_> = match form {
                        Form::Stress => p
                            .plane_stress_states(&flat, settings, cx)
                            .unwrap()
                            .into_iter()
                            .map(|s| s.deformation_gradient)
                            .collect(),
                        Form::Exact => p
                            .incompressible_plane_stress_states(&flat, settings, cx)
                            .unwrap()
                            .into_iter()
                            .map(|s| s.deformation_gradient)
                            .collect(),
                        Form::Strain => unreachable!(),
                    };
                    let points = if mesh.elems[0].len() == 3 { 1 } else { 4 };
                    for (i, f) in gradients.iter().enumerate() {
                        for (actual, expected) in
                            f.iter().zip(references[i / points].deformation_gradient())
                        {
                            close(*actual, *expected, 4e-8);
                        }
                    }
                }
            }
        }
    });
}

#[test]
fn g0_general_reference_assembly_preserves_work_tangent_and_objectivity() {
    with_cx(|cx, _| {
        let original = Hyperelastic::new(
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
        let reference =
            ReferenceDeformation::try_new([1.07, 0.06, 0.08, -0.03, 0.97, 0.04, 0.02, -0.05, 1.10])
                .unwrap();
        let law = ReferencedHyperelastic::new(&original, reference);
        let mesh = Mesh2::mapped_quads(1, 1, &|x, y| [x + 0.12 * y, y + 0.07 * x]);
        let p = sheet(&mesh, &law);
        let u: Vec<_> = mesh
            .nodes
            .iter()
            .flat_map(|[x, y]| [0.14 * x + 0.05 * y, -0.02 * x - 0.03 * y])
            .collect();
        for form in FORMS {
            let (force, tangent) = probe(&p, &u, form, cx);
            let h = 2e-6;
            for col in 0..u.len() {
                let mut plus = u.clone();
                let mut minus = u.clone();
                plus[col] += h;
                minus[col] -= h;
                close(
                    (potential(&p, &plus, form, cx) - potential(&p, &minus, form, cx)) / (2.0 * h),
                    force[col],
                    5e-7,
                );
                let rp = probe(&p, &plus, form, cx).0;
                let rm = probe(&p, &minus, form, cx).0;
                for row in 0..u.len() {
                    close((rp[row] - rm[row]) / (2.0 * h), tangent.get(row, col), 8e-6);
                    close(tangent.get(row, col), tangent.get(col, row), 8e-11);
                }
            }
            let (s, c) = 0.43_f64.sin_cos();
            let rotated: Vec<_> = mesh
                .nodes
                .iter()
                .enumerate()
                .flat_map(|(i, [x, y])| {
                    let a = x + u[2 * i];
                    let b = y + u[2 * i + 1];
                    [c * a - s * b - x + 2.0, s * a + c * b - y - 1.0]
                })
                .collect();
            close(
                potential(&p, &rotated, form, cx),
                potential(&p, &u, form, cx),
                4e-10,
            );
            let rf = probe(&p, &rotated, form, cx).0;
            for i in 0..mesh.nodes.len() {
                close(rf[2 * i], c * force[2 * i] - s * force[2 * i + 1], 5e-9);
                close(rf[2 * i + 1], s * force[2 * i] + c * force[2 * i + 1], 5e-9);
            }
        }
        assert!(matches!(
            p.solve_mean_dilatation(cx),
            Err(SolidError::InvalidInput { .. })
        ));
    });
}

#[test]
fn g0_relaxed_state_can_be_admissible_when_identity_is_outside_elastic_domain() {
    with_cx(|cx, _| {
        let original = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 40.0,
            },
            1.15,
        )
        .unwrap();
        let reference =
            ReferenceDeformation::try_new([2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0]).unwrap();
        let law = ReferencedHyperelastic::new(&original, reference);
        let mesh = Mesh2::triangles(1.0, 0.7, 1, 1);
        let p = sheet(&mesh, &law);
        let expanded: Vec<_> = mesh.nodes.iter().flat_map(|[x, y]| [*x, *y]).collect();
        let settings = PlaneStressSettings::default();
        for state in p.plane_stress_states(&expanded, settings, cx).unwrap() {
            close(state.deformation_gradient[8], 2.0, 2e-10);
            assert!(state.piola.iter().all(|p| p.abs() < 1e-9));
        }
        for state in p
            .incompressible_plane_stress_states(&expanded, settings, cx)
            .unwrap()
        {
            close(state.deformation_gradient[8], 2.0, 2e-10);
            assert!(state.piola.iter().all(|p| p.abs() < 1e-9));
        }
        let identity = vec![0.0; 2 * mesh.nodes.len()];
        assert!(p.plane_stress_states(&identity, settings, cx).is_err());
        assert!(
            p.incompressible_plane_stress_states(&identity, settings, cx)
                .is_err()
        );
        // Plane strain constrains total thickness to one; it cannot claim the
        // same stress-free state merely because in-plane motion is admitted.
        assert!(p.residual_and_tangent(&expanded, 0.0).is_err());
    });
}

#[test]
fn g1_frozen_prony_sheet_equilibrium_retains_history_and_incremental_potential() {
    with_cx(|cx, _| {
        let equilibrium = Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 1.0,
                c01: 0.2,
                kappa: 30.0,
            },
            4.0,
        )
        .unwrap();
        let law = ViscoHyperelastic::new(
            equilibrium,
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 3.0 },
                tau_s: 0.6,
            }],
        )
        .unwrap();
        let previous_f = [1.18, 0.08, 0.05, 0.0, 0.95, 0.02, 0.0, 0.0, 1.01];
        let accepted = law.step(&law.initial_state(), &previous_f, 0.15).unwrap().0;
        let before = accepted.clone();
        let dt = 0.25;
        let frozen = law.freeze(&accepted, dt).unwrap();
        let in_plane = [1.12, 0.04, 0.01, 0.96];
        let affine = |x: f64, y: f64| [0.12 * x + 0.04 * y, 0.01 * x - 0.04 * y];
        let settings = PlaneStressSettings::default();
        for mesh in [
            Mesh2::triangles(1.0, 0.7, 2, 2),
            Mesh2::quads(1.0, 0.7, 2, 2),
        ] {
            let mut p = sheet(&mesh, &frozen);
            p.dirichlet = vec![
                (Patch::Left, &affine),
                (Patch::Right, &affine),
                (Patch::Bottom, &affine),
                (Patch::Top, &affine),
            ];
            let u = solve(&p, Form::Stress, cx);
            for (actual, &[x, y]) in u.iter().zip(&mesh.nodes) {
                let expected = affine(x, y);
                close(actual[0], expected[0], 3e-8);
                close(actual[1], expected[1], 3e-8);
            }
            let flat: Vec<_> = u.iter().flatten().copied().collect();
            let states = p.plane_stress_states(&flat, settings, cx).unwrap();
            for state in &states {
                for (i, j) in [0, 1, 3, 4].into_iter().enumerate() {
                    close(state.deformation_gradient[j], in_plane[i], 3e-8);
                }
                let candidate = law
                    .step_response(&accepted, &state.deformation_gradient, dt)
                    .unwrap();
                for (i, j) in [0, 1, 3, 4].into_iter().enumerate() {
                    close(state.piola[i], candidate.first_piola_pa[j], 3e-9);
                }
                close(
                    state.energy_density_pa,
                    candidate.incremental_potential_pa,
                    3e-10,
                );
                assert!(candidate.dissipated_pa > 0.0);
                assert!((candidate.incremental_potential_pa - candidate.state.stored).abs() > 1e-6);
                let backbone = law.equilibrium.piola(&state.deformation_gradient).unwrap();
                assert!(
                    [0, 1, 3, 4]
                        .into_iter()
                        .any(|i| (backbone[i] - candidate.first_piola_pa[i]).abs() > 1e-4)
                );
            }
            let density = frozen.energy_density(&states[0].deformation_gradient);
            close(
                p.plane_stress_potential_energy(&flat, 0.0, settings, cx)
                    .unwrap(),
                0.7 * density,
                3e-9,
            );
            assert!(matches!(
                p.solve_mean_dilatation(cx),
                Err(SolidError::InvalidInput { .. })
            ));
        }
        assert_eq!(accepted, before);
    });
}

#[test]
fn g1_continuation_uses_admitted_deformed_geometry_without_resetting_reference() {
    with_cx(|cx, _| {
        let original = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 40.0,
            },
            1.15,
        )
        .unwrap();
        let settings = PlaneStressSettings::default();
        for mesh in [
            Mesh2::triangles(1.0, 0.7, 2, 2),
            Mesh2::quads(1.0, 0.7, 2, 2),
        ] {
            for form in FORMS {
                let thickness = if matches!(form, Form::Strain) {
                    1.0
                } else {
                    4.1
                };
                let reference = ReferenceDeformation::try_new([
                    4.1, 0.0, 0.0, 0.0, 4.1, 0.0, 0.0, 0.0, thickness,
                ])
                .unwrap();
                let material = ReferencedHyperelastic::new(&original, reference);
                let left = |_: f64, y: f64| [0.0, 3.1 * y];
                let mut p = sheet(&mesh, &material);
                p.dirichlet = vec![(Patch::Left, &left)];
                // The previous body's in-plane stretch is four. A new modest
                // change of natural strain is admissible there; restarting at
                // identity lies far outside this same elastic material domain.
                let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| [3.0 * x, 3.0 * y]).collect();
                let before = initial.clone();
                let zero_start = match form {
                    Form::Strain => p.solve(),
                    Form::Stress => p.solve_plane_stress(settings, cx),
                    Form::Exact => p.solve_incompressible_plane_stress(settings, cx),
                };
                assert!(
                    zero_start.is_err(),
                    "{form:?} must admit actual initial geometry"
                );
                let run = |start: &[[f64; 2]]| match form {
                    Form::Strain => p.solve_from(start),
                    Form::Stress => p.solve_plane_stress_from(start, settings, cx),
                    Form::Exact => p.solve_incompressible_plane_stress_from(start, settings, cx),
                };
                let (u, report) = run(&initial).unwrap();
                assert_eq!(report.histories.len(), 1);
                assert!(*report.histories[0].last().unwrap() < p.settings.tol);
                for (actual, &[x, y]) in u.iter().zip(&mesh.nodes) {
                    close(actual[0], 3.1 * x, 4e-8);
                    close(actual[1], 3.1 * y, 4e-8);
                }
                let flat: Vec<_> = u.iter().flatten().copied().collect();
                close(potential(&p, &flat, form, cx), 0.0, 2e-9);
                assert_eq!(initial, before);
                assert!(matches!(
                    run(&initial[..initial.len() - 1]),
                    Err(SolidError::InvalidInput { .. })
                ));
                let mut invalid = initial.clone();
                invalid[0][0] = f64::NAN;
                assert!(matches!(
                    run(&invalid),
                    Err(SolidError::InvalidInput { .. })
                ));
                assert!(matches!(
                    run(&vec![[0.0; 2]; mesh.nodes.len()]),
                    Err(SolidError::MaterialRefused { .. })
                ));
                assert_eq!(initial, before);
            }
        }
    });
}
