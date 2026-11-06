//! G1/G3: element-owned hyperelastic laws retain equilibrium and admissibility.
//! A layered affine reference has zero traction across its material interface.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper2d::{HyperelasticMaterials, NewtonReport};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};
use fs_sparse::Csr;

#[path = "hyper_material_field/points.rs"]
mod points;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 51,
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
    Incompressible,
}

const FORMULATIONS: [Formulation; 4] = [
    Formulation::PlaneStrain,
    Formulation::MeanDilatation,
    Formulation::PlaneStress,
    Formulation::Incompressible,
];

impl Formulation {
    fn solve(
        self,
        p: &HyperProblem<'_>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        let local = PlaneStressSettings::default();
        match self {
            Self::PlaneStrain => p.solve(),
            Self::MeanDilatation => p.solve_mean_dilatation(cx),
            Self::PlaneStress => p.solve_plane_stress(local, cx),
            Self::Incompressible => p.solve_incompressible_plane_stress(local, cx),
        }
    }

    fn probe(
        self,
        p: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Csr), SolidError> {
        let local = PlaneStressSettings::default();
        match self {
            Self::PlaneStrain => p.residual_and_tangent(u, load),
            Self::MeanDilatation => p.mean_dilatation_residual_and_tangent(u, load, cx),
            Self::PlaneStress => p.plane_stress_residual_and_tangent(u, load, local, cx),
            Self::Incompressible => {
                p.incompressible_plane_stress_residual_and_tangent(u, load, local, cx)
            }
        }
    }

    fn stored_energy(self, p: &HyperProblem<'_>, u: &[f64], cx: &Cx<'_>) -> Option<f64> {
        // Zero load removes dead-load work, without changing the supplied motion.
        let local = PlaneStressSettings::default();
        match self {
            Self::PlaneStrain => p.potential_energy(u, 0.0),
            Self::MeanDilatation => p.mean_dilatation_potential_energy(u, 0.0, cx).ok(),
            Self::PlaneStress => p.plane_stress_potential_energy(u, 0.0, local, cx).ok(),
            Self::Incompressible => p
                .incompressible_plane_stress_potential_energy(u, 0.0, local, cx)
                .ok(),
        }
    }
}

fn law(mu: f64) -> Hyperelastic {
    Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda: 0.0 }, 4.0).unwrap()
}

fn problem<'a>(mesh: &'a Mesh2, material: &'a dyn HyperelasticMaterials) -> HyperProblem<'a> {
    HyperProblem {
        mesh,
        material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings {
            tol: 2e-9,
            max_iters: 50,
            load_steps: 1,
            max_backtracks: 30,
        },
    }
}

fn element_mu(mesh: &Mesh2, index: usize) -> f64 {
    let y = 0.25
        * mesh.elems[index]
            .iter()
            .map(|&n| mesh.nodes[n][1])
            .sum::<f64>();
    if y < 0.5 { 1.0 } else { 8.0 }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 5e-8 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_layered_affine_equilibrium_reactions_and_energy_use_each_original_law() {
    with_cx(|cx| {
        let mesh = Mesh2::quads(2.0, 1.0, 2, 2);
        let laws: Vec<_> = (0..mesh.elems.len())
            .map(|i| law(element_mu(&mesh, i)))
            .collect();
        for formulation in FORMULATIONS {
            let axial: f64 = 1.25;
            let (lateral, px, energy_density) =
                if matches!(formulation, Formulation::Incompressible) {
                    (
                        axial.sqrt().recip(),
                        axial - axial.powi(-2),
                        0.5 * (axial * axial + 2.0 / axial - 3.0),
                    )
                } else {
                    (
                        1.0,
                        axial - axial.recip(),
                        0.5 * (axial * axial - 1.0) - axial.ln(),
                    )
                };
            let left = |_: f64, y: f64| [0.0, (lateral - 1.0) * y];
            let right = |_: f64, y: f64| [(if y < 0.5 { 1.0 } else { 8.0 }) * px, 0.0];
            let mut p = problem(&mesh, &laws);
            p.dirichlet = vec![(Patch::Left, &left)];
            p.traction = vec![(Patch::Right, &right)];
            let (motion, report) = formulation.solve(&p, cx).unwrap();
            assert_eq!(report.histories.len(), 1);
            assert!(
                report.histories[0].len() > 1,
                "{formulation:?} must actually solve"
            );
            assert!(*report.histories[0].last().unwrap() < p.settings.tol);
            for (u, point) in motion.iter().zip(&mesh.nodes) {
                close(u[0], (axial - 1.0) * point[0]);
                close(u[1], (lateral - 1.0) * point[1]);
            }
            let flat: Vec<f64> = motion.iter().flatten().copied().collect();
            let (forces, _) = formulation.probe(&p, &flat, 1.0, cx).unwrap();
            let supports = mesh.patch_nodes(Patch::Left);
            let reaction: f64 = supports.iter().map(|&node| forces[2 * node]).sum();
            close(reaction, -0.5 * (1.0 + 8.0) * px);
            for (node, force) in forces.chunks_exact(2).enumerate() {
                if !supports.contains(&node) {
                    close(force[0], 0.0);
                    close(force[1], 0.0);
                }
            }
            // Four area-1/2 elements: two of each original shear modulus.
            let weighted_area: f64 = (0..mesh.elems.len())
                .map(|i| 0.5 * element_mu(&mesh, i))
                .sum();
            close(
                formulation.stored_energy(&p, &flat, cx).unwrap(),
                weighted_area * energy_density,
            );
            let local = PlaneStressSettings::default();
            let recovered: Vec<[f64; 4]> = match formulation {
                Formulation::PlaneStress => p
                    .plane_stress_states(&flat, local, cx)
                    .unwrap()
                    .into_iter()
                    .map(|state| state.piola)
                    .collect(),
                Formulation::Incompressible => p
                    .incompressible_plane_stress_states(&flat, local, cx)
                    .unwrap()
                    .into_iter()
                    .map(|state| state.piola)
                    .collect(),
                _ => vec![],
            };
            if !recovered.is_empty() {
                assert_eq!(recovered.len(), 4 * mesh.elems.len());
                for (element, states) in recovered.chunks_exact(4).enumerate() {
                    for piola in states {
                        close(piola[0], element_mu(&mesh, element) * px);
                        for &value in &piola[1..] {
                            close(value, 0.0);
                        }
                    }
                }
            }
        }
    });
}

#[test]
fn g3_uniform_and_repeated_vector_array_or_borrowed_slice_are_equivalent() {
    with_cx(|cx| {
        let mesh = Mesh2::quads(2.0, 1.0, 2, 2);
        let uniform = law(3.0);
        let repeated = vec![uniform.clone(); mesh.elems.len()];
        let array = [
            uniform.clone(),
            uniform.clone(),
            uniform.clone(),
            uniform.clone(),
        ];
        let borrowed = repeated.as_slice();
        let baseline = problem(&mesh, &uniform);
        let mapped = [
            problem(&mesh, &repeated),
            problem(&mesh, &array),
            problem(&mesh, &borrowed),
        ];
        let motion: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|&[x, y]| [0.08 * x + 0.015 * x * y, -0.03 * y + 0.01 * x * y])
            .collect();
        for formulation in FORMULATIONS {
            let (force, tangent) = formulation.probe(&baseline, &motion, 0.0, cx).unwrap();
            let energy = formulation.stored_energy(&baseline, &motion, cx).unwrap();
            for p in &mapped {
                let (other_force, other_tangent) = formulation.probe(p, &motion, 0.0, cx).unwrap();
                assert_eq!(other_force, force, "{formulation:?}");
                assert_eq!(other_tangent.nrows(), tangent.nrows());
                for row in 0..tangent.nrows() {
                    assert_eq!(
                        other_tangent.row(row),
                        tangent.row(row),
                        "{formulation:?}, row {row}"
                    );
                }
                assert_eq!(formulation.stored_energy(p, &motion, cx).unwrap(), energy);
            }
        }
    });
}

#[test]
fn g0_field_count_and_late_material_refusals_cover_all_formulations() {
    with_cx(|cx| {
        let mesh = Mesh2::quads(2.0, 1.0, 2, 2);
        let zero = vec![0.0; 2 * mesh.node_count()];
        for count in [0, mesh.elems.len() - 1, mesh.elems.len() + 1] {
            let laws = vec![law(1.0); count];
            let p = problem(&mesh, &laws);
            for formulation in FORMULATIONS {
                assert!(matches!(
                    formulation.solve(&p, cx),
                    Err(SolidError::InvalidInput { .. })
                ));
                assert!(matches!(
                    formulation.probe(&p, &zero, 0.0, cx),
                    Err(SolidError::InvalidInput { .. })
                ));
                assert!(formulation.stored_energy(&p, &zero, cx).is_none());
            }
        }
        let mut laws = vec![law(1.0); mesh.elems.len()];
        laws.last_mut().unwrap().model = HyperelasticModel::NeoHookean {
            mu: -1.0,
            lambda: 0.0,
        };
        let p = problem(&mesh, &laws);
        for formulation in FORMULATIONS {
            assert!(matches!(
                formulation.solve(&p, cx),
                Err(SolidError::MaterialRefused { .. })
            ));
            assert!(matches!(
                formulation.probe(&p, &zero, 0.0, cx),
                Err(SolidError::MaterialRefused { .. })
            ));
            assert!(formulation.stored_energy(&p, &zero, cx).is_none());
        }
        // The first law remains valid: only the final element has a tighter domain.
        let mut laws = vec![law(1.0); mesh.elems.len()];
        laws.last_mut().unwrap().stretch_limit = 1.1;
        let p = problem(&mesh, &laws);
        let motion: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|&[x, _]| [0.25 * x, 0.0])
            .collect();
        for formulation in FORMULATIONS {
            assert!(matches!(
                formulation.probe(&p, &motion, 0.0, cx),
                Err(SolidError::MaterialRefused { .. })
            ));
            assert!(formulation.stored_energy(&p, &motion, cx).is_none());
        }
        let local = PlaneStressSettings::default();
        assert!(p.unconstrained_residual(&motion, 0.0).is_err());
        assert!(p.mean_dilatation_states(&motion, cx).is_err());
        assert!(p.plane_stress_states(&motion, local, cx).is_err());
        assert!(
            p.incompressible_plane_stress_states(&motion, local, cx)
                .is_err()
        );
    });
}
