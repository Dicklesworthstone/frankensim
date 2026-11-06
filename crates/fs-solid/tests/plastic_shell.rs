//! G0/G1/G3 native CST/DKT elasticity, plastic bending, springback and restart.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::{
    LayeredPlaneStressSection, PlaneStressLayer, SectionBudget, ThicknessQuadrature,
};
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_plate::PlateSection;
use fs_plate::shell::{DRILLING_ALPHA, ShellMesh, ShellSupport, assemble_shell};
use fs_solid::plastic_shell::{
    PlasticShellBudget, PlasticShellSections, SmallStrainShellProblem, SmallStrainShellSolution,
};
use fs_solid::{IncrementSettings, SmallStrainTetError};

const E: f64 = 1000.0;
const NU: f64 = 0.3;
const H: f64 = 0.1;
const MOMENT: f64 = 0.00045;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 83,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "{actual:e} != {expected:e}; allowance {allowance:e}"
    );
}

fn metal() -> J2Plasticity {
    J2Plasticity::new_mixed(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        1.0,
        100.0,
        50.0,
    )
    .unwrap()
}

fn section<L: SmallStrainLaw>(law: L) -> LayeredPlaneStressSection<L> {
    LayeredPlaneStressSection::new(
        vec![PlaneStressLayer {
            bottom_z_m: -0.5 * H,
            top_z_m: 0.5 * H,
            subdivisions: 4,
            law,
        }],
        ThicknessQuadrature::Gauss3,
        ReductionConfig {
            max_iterations: 48,
            stress_tolerance_pa: 1e-12,
            strain_limit: 0.05,
        },
        SectionBudget {
            maximum_layers: 4,
            maximum_quadrature_points: 64,
        },
    )
    .unwrap()
}

fn problem<'a, L: SmallStrainLaw>(
    mesh: &'a ShellMesh,
    section: &'a LayeredPlaneStressSection<L>,
    prescribed: &'a [(usize, f64)],
    loads: &'a [(usize, f64)],
    path: &'a [f64],
    coordinate_length_m: f64,
) -> SmallStrainShellProblem<'a, L> {
    SmallStrainShellProblem {
        mesh,
        sections: PlasticShellSections::Uniform(section),
        coordinate_length_m,
        rotation_limit_rad: 0.05,
        drilling_alpha: DRILLING_ALPHA,
        prescribed,
        nodal_loads: loads,
        load_path: path,
        budget: PlasticShellBudget::default(),
        settings: IncrementSettings {
            max_newton_iterations: 40,
            ..IncrementSettings::default()
        },
    }
}

#[test]
fn g1_all_prescribed_elastic_shell_matches_original_full_matrix_in_oblique_frames() {
    with_cx(|_, cx| {
        for nodes in [
            vec![[0.1, -0.2, 0.0], [1.3, 0.1, 0.0], [0.4, 0.9, 0.0]],
            vec![[0.1, -0.2, 0.3], [1.3, 0.1, 0.4], [0.4, 0.9, 0.7]],
        ] {
            let mesh = ShellMesh::new(nodes, vec![[0, 1, 2]]).unwrap();
            let material = section(IsotropicElastic::new(E, NU, 0.05).unwrap());
            let original = assemble_shell(
                &mesh,
                &PlateSection::isotropic(E, NU, H, 2700.0).unwrap(),
                &[],
                ShellSupport::Free,
            )
            .unwrap();
            let coordinates: Vec<f64> = (0..18)
                .map(|d| 0.00017 * ((17 * d % 11) as f64 - 5.0))
                .collect();
            let prescribed: Vec<_> = coordinates.iter().copied().enumerate().collect();
            let mut expected = vec![0.0; 18];
            original.k.spmv(&coordinates, &mut expected);
            let solution = problem(&mesh, &material, &prescribed, &[], &[1.0], 0.37)
                .solve(cx)
                .unwrap();
            for (actual, reference) in solution.increments[0].reactions.iter().zip(expected) {
                near(*actual, reference, 2e-12);
            }
            for (actual, reference) in solution.coordinates.iter().flatten().zip(coordinates) {
                near(*actual, reference, 4e-19);
            }
            assert_eq!(solution.increments[0].evolving_elements, 0);
            assert_eq!(solution.section_states.len(), 3);
            for state in &solution.section_states {
                material.validate_state(state).unwrap();
            }
        }
    });
}

#[test]
fn g1_eccentric_section_follows_physical_director_sign_and_membrane_bending_work() {
    with_cx(|_, cx| {
        let mesh = ShellMesh::new(
            vec![[0.0, 0.0, 0.0], [0.6, 0.0, 0.0], [0.2, 0.4, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        // Reference surface at the physical bottom; positive z is along the
        // facet normal. This makes B nonzero and pins the curvature sign.
        let material = LayeredPlaneStressSection::new(
            vec![PlaneStressLayer {
                bottom_z_m: 0.0,
                top_z_m: H,
                subdivisions: 3,
                law: IsotropicElastic::new(E, NU, 0.05).unwrap(),
            }],
            ThicknessQuadrature::Gauss3,
            ReductionConfig {
                max_iterations: 32,
                stress_tolerance_pa: 1e-12,
                strain_limit: 0.05,
            },
            SectionBudget {
                maximum_layers: 1,
                maximum_quadrature_points: 32,
            },
        )
        .unwrap();
        let membrane = [0.001, -0.0007, 0.0008];
        let (a, b, c) = (0.012, 0.004, -0.006);
        let curvature = [-a, -c, -2.0 * b];
        let mut coordinates = Vec::new();
        for &[x, y, _] in &mesh.nodes {
            coordinates.extend([
                membrane[0] * x + 0.5 * membrane[2] * y,
                membrane[1] * y + 0.5 * membrane[2] * x,
                0.5 * a * x * x + b * x * y + 0.5 * c * y * y,
                b * x + c * y,
                -a * x - b * y,
                0.0,
            ]);
        }
        let prescribed: Vec<_> = coordinates.iter().copied().enumerate().collect();
        let solution = problem(&mesh, &material, &prescribed, &[], &[1.0], 1.0)
            .solve(cx)
            .unwrap();
        let diagonal = E / (1.0 - NU * NU);
        let q = [
            [diagonal, NU * diagonal, 0.0],
            [NU * diagonal, diagonal, 0.0],
            [0.0, 0.0, E / (2.0 * (1.0 + NU))],
        ];
        let mut resultants = [0.0; 6];
        for i in 0..3 {
            for j in 0..3 {
                resultants[i] += q[i][j] * (H * membrane[j] + 0.5 * H * H * curvature[j]);
                resultants[i + 3] +=
                    q[i][j] * (0.5 * H * H * membrane[j] + H * H * H * curvature[j] / 3.0);
            }
        }
        let generalized = [
            membrane[0],
            membrane[1],
            membrane[2],
            curvature[0],
            curvature[1],
            curvature[2],
        ];
        for (state, actual) in solution
            .section_states
            .iter()
            .zip(&solution.section_resultants)
        {
            for i in 0..6 {
                near(state.generalized_strain[i], generalized[i], 3e-17);
                near(actual[i], resultants[i], 2e-15);
            }
            for (point, layer_state) in material.points().iter().zip(&state.points) {
                near(layer_state.strain[0], membrane[0] - point.z_m * a, 3e-18);
                near(layer_state.strain[1], membrane[1] - point.z_m * c, 3e-18);
                near(
                    layer_state.strain[3],
                    0.5 * membrane[2] - point.z_m * b,
                    3e-18,
                );
            }
        }
        let external_work: f64 = solution.increments[0]
            .reactions
            .iter()
            .zip(coordinates)
            .map(|(reaction, displacement)| reaction * displacement)
            .sum();
        let material_work: f64 = resultants
            .iter()
            .zip(generalized)
            .map(|(force, strain)| force * strain)
            .sum::<f64>()
            * 0.12;
        near(external_work, material_work, 2e-17);
    });
}

fn cantilever(nx: usize, ny: usize) -> (ShellMesh, Vec<(usize, f64)>, Vec<(usize, f64)>) {
    let mut nodes = Vec::new();
    for i in 0..=nx {
        for j in 0..=ny {
            nodes.push([
                i as f64 / nx as f64,
                0.2 * (j as f64 / ny as f64 - 0.5),
                0.0,
            ]);
        }
    }
    let node = |i, j| i * (ny + 1) + j;
    let mut tris = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            tris.push([node(i, j), node(i + 1, j), node(i + 1, j + 1)]);
            tris.push([node(i, j), node(i + 1, j + 1), node(i, j + 1)]);
        }
    }
    let prescribed = (0..=ny)
        .flat_map(|j| (0..6).map(move |c| (6 * node(0, j) + c, 0.0)))
        .collect();
    let loads = (0..=ny)
        .map(|j| {
            let weight = if j == 0 || j == ny { 0.5 } else { 1.0 };
            (6 * node(nx, j) + 4, MOMENT * weight / ny as f64)
        })
        .collect();
    (ShellMesh::new(nodes, tris).unwrap(), prescribed, loads)
}

fn mean_tip_rotation(mesh: &ShellMesh, solution: &SmallStrainShellSolution<J2State>) -> f64 {
    let values: Vec<_> = mesh
        .nodes
        .iter()
        .zip(&solution.coordinates)
        .filter(|(p, _)| p[0] == 1.0)
        .map(|(_, q)| q[4])
        .collect();
    values.iter().sum::<f64>() / values.len() as f64
}

fn reaction_wrench(mesh: &ShellMesh, fixed: &[(usize, f64)], reactions: &[f64]) -> [f64; 6] {
    let mut wrench = [0.0; 6];
    for (&(dof, _), &value) in fixed.iter().zip(reactions) {
        let component = dof % 6;
        if component >= 3 {
            wrench[component] += value;
        } else {
            wrench[component] += value;
            let mut force = [0.0; 3];
            force[component] = value;
            let r = mesh.nodes[dof / 6];
            wrench[3] += r[1] * force[2] - r[2] * force[1];
            wrench[4] += r[2] * force[0] - r[0] * force[2];
            wrench[5] += r[0] * force[1] - r[1] * force[0];
        }
    }
    wrench
}

#[test]
fn g1_free_cantilever_yields_through_thickness_and_springs_back_to_permanent_shape() {
    with_cx(|_, cx| {
        let (mesh, fixed, loads) = cantilever(4, 2);
        let material = section(metal());
        let loaded = problem(
            &mesh,
            &material,
            &fixed,
            &loads,
            &[0.25, 0.5, 0.75, 1.0],
            1.0,
        )
        .solve(cx)
        .unwrap();
        assert!(
            loaded
                .increments
                .iter()
                .any(|step| step.evolving_elements > 0)
        );
        assert!(
            loaded
                .increments
                .iter()
                .any(|step| step.residual_history.len() > 2)
        );
        let all_points: Vec<_> = loaded
            .section_states
            .iter()
            .flat_map(|s| &s.points)
            .collect();
        assert!(all_points.iter().any(|p| p.history.alpha > 1e-5));
        assert!(all_points.iter().any(|p| p.history.alpha == 0.0));
        let peak = mean_tip_rotation(&mesh, &loaded);
        assert!(peak > 0.0 && peak < 0.05, "peak rotation {peak}");
        let reaction_moment: f64 = fixed
            .iter()
            .zip(&loaded.increments.last().unwrap().reactions)
            .filter(|((dof, _), _)| dof % 6 == 4)
            .map(|(_, value)| value)
            .sum();
        near(reaction_moment, -MOMENT, 2e-10);
        let unloaded = problem(&mesh, &material, &fixed, &loads, &[0.5, 0.0], 1.0)
            .solve_from(&loaded.coordinates, &loaded.section_states, cx)
            .unwrap();
        let residual = mean_tip_rotation(&mesh, &unloaded);
        assert!(
            residual > 0.03 * peak && residual < 0.9 * peak,
            "springback peak={peak}, residual={residual}"
        );
        assert!(unloaded.coordinates.iter().any(|q| q[2].abs() > 1e-4));
        for (before, after) in loaded.section_states.iter().zip(&unloaded.section_states) {
            for (p, q) in before.points.iter().zip(&after.points) {
                assert!(q.history.alpha >= p.history.alpha);
            }
        }
        let last = unloaded.increments.last().unwrap();
        assert!(last.residual_history.last().unwrap() <= &1e-10);
        // Residual clamp tractions may self-equilibrate. It is their complete
        // resultant wrench, not each individual support reaction, that vanishes.
        for reaction in reaction_wrench(&mesh, &fixed, &last.reactions) {
            near(reaction, 0.0, 2e-9);
        }
    });
}

#[test]
fn g3_numerical_rotation_length_does_not_change_plastic_cantilever_response() {
    with_cx(|_, cx| {
        let (mesh, fixed, loads) = cantilever(2, 1);
        let material = section(metal());
        let path = [0.5, 1.0, 0.0];
        let reference = problem(&mesh, &material, &fixed, &loads, &path, 1.0)
            .solve(cx)
            .unwrap();
        for length in [0.1, 10.0] {
            let candidate = problem(&mesh, &material, &fixed, &loads, &path, length)
                .solve(cx)
                .unwrap();
            for (a, b) in reference
                .coordinates
                .iter()
                .flatten()
                .zip(candidate.coordinates.iter().flatten())
            {
                near(*a, *b, 5e-9);
            }
            for (a, b) in reference
                .section_states
                .iter()
                .zip(&candidate.section_states)
            {
                for (a, b) in a.points.iter().zip(&b.points) {
                    near(a.history.alpha, b.history.alpha, 2e-9);
                    for i in 0..6 {
                        near(
                            a.history.plastic_strain[i],
                            b.history.plastic_strain[i],
                            2e-9,
                        );
                    }
                }
            }
            for (a, b) in reference.increments.iter().zip(&candidate.increments) {
                for (left, right) in a.reactions.iter().zip(&b.reactions) {
                    near(*left, *right, 3e-10);
                }
            }
        }
    });
}

#[test]
fn g0_checkpoint_resume_late_refusal_and_cancellation_keep_original_state() {
    with_cx(|gate, cx| {
        let (mesh, fixed, loads) = cantilever(2, 1);
        let material = section(metal());
        let loaded = problem(&mesh, &material, &fixed, &loads, &[0.5, 1.0], 1.0)
            .solve(cx)
            .unwrap();
        let saved = loaded.clone();
        let resumed = problem(&mesh, &material, &fixed, &loads, &[0.0], 1.0)
            .solve_from(&loaded.coordinates, &loaded.section_states, cx)
            .unwrap();
        let continuous = problem(&mesh, &material, &fixed, &loads, &[0.5, 1.0, 0.0], 1.0)
            .solve(cx)
            .unwrap();
        assert_eq!(resumed.coordinates, continuous.coordinates);
        assert_eq!(resumed.section_states, continuous.section_states);
        let refused = problem(&mesh, &material, &fixed, &loads, &[0.0, 1e6], 1.0).solve_from(
            &loaded.coordinates,
            &loaded.section_states,
            cx,
        );
        assert!(
            refused.is_err(),
            "a late large-rotation proposal must refuse"
        );
        assert_eq!(loaded, saved);
        let mut altered = loaded.coordinates.clone();
        altered.last_mut().unwrap()[4] += 1e-4;
        assert!(
            problem(&mesh, &material, &fixed, &loads, &[0.0], 1.0)
                .solve_from(&altered, &loaded.section_states, cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            problem(&mesh, &material, &fixed, &loads, &[0.0], 1.0).solve_from(
                &loaded.coordinates,
                &loaded.section_states,
                cx
            ),
            Err(SmallStrainTetError::Cancelled)
        ));
        assert_eq!(loaded, saved);
    });
}
