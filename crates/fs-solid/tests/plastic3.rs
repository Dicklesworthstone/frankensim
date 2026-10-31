//! Incremental small-strain J2 plasticity on P1 tetrahedra.
//!
//! G1: homogeneous uniaxial tension with linear isotropic hardening follows
//! the closed-form stress-strain path through yield, and elastic unloading
//! to zero force leaves exactly the plastic strain and the plastic
//! (isochoric) lateral contraction.
//! G0: the assembled algorithmic tangent is the derivative of the internal
//! force at a plastic state (central differences), and the elastic law's
//! tangent equals the independent `linear3` stiffness.
//! G3: a cantilever under increasing tip displacement yields progressively
//! (monotone plastic zone, softening secant stiffness).

use fs_alloc::{ArenaConfig, ArenaPool};
#[path = "plastic3/dynamics.rs"]
mod dynamics;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::linear3::{TetLinearElasticProblem, TetMaterialField};
use fs_solid::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, TetAssemblyBudget,
    TetElasticMaterial,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 5,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

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

const E: f64 = 1000.0;
const NU: f64 = 0.3;
const SIGMA_Y: f64 = 1.0;
const H: f64 = 100.0;

fn j2() -> J2Plasticity {
    J2Plasticity::new(IsotropicElastic::new(E, NU, 0.05).unwrap(), SIGMA_Y, H).unwrap()
}

#[test]
fn uniaxial_hardening_path_and_unloading_match_closed_form() {
    let size = [2.0, 1.0, 1.0];
    let (nodes, tets) = box_mesh([4, 2, 2], size);
    let law = j2();
    // Reference prescribed field: unit axial strain on the far face, rollers
    // on the three symmetry planes.
    let mut prescribed = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        if x[0].abs() < 1e-12 {
            prescribed.push((3 * node, 0.0));
        }
        if (x[0] - size[0]).abs() < 1e-12 {
            prescribed.push((3 * node, size[0]));
        }
        if x[1].abs() < 1e-12 {
            prescribed.push((3 * node + 1, 0.0));
        }
        if x[2].abs() < 1e-12 {
            prescribed.push((3 * node + 2, 0.0));
        }
    }
    let stress = |eps: f64| {
        if E * eps <= SIGMA_Y {
            E * eps
        } else {
            E * (SIGMA_Y + H * eps) / (E + H)
        }
    };
    let peak = 0.006;
    let plastic = peak - stress(peak) / E;
    let path = [0.0005, 0.001, 0.002, 0.004, peak, plastic];
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::Uniform(&law),
        prescribed_m: &prescribed,
        nodal_forces_n: &[],
        load_path: &path,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    let solution = with_cx(|cx| problem.solve(cx)).unwrap();
    let expected = [
        stress(path[0]),
        stress(path[1]),
        stress(path[2]),
        stress(path[3]),
        stress(peak),
        0.0,
    ];
    for ((record, &eps), &expected) in solution.increments.iter().zip(&path).zip(&expected) {
        let axial: f64 = prescribed
            .iter()
            .zip(&record.reactions_n)
            .filter(|((dof, value), _)| dof % 3 == 0 && *value > 0.0)
            .map(|(_, r)| r)
            .sum();
        println!(
            "eps {eps:.7}: force {axial:.12} expected {expected:.12} iterations {} evolving {}",
            record.residual_history.len(),
            record.evolving_elements
        );
        assert!(
            (axial - expected).abs() < 1e-9,
            "eps {eps}: {axial} vs {expected}"
        );
        assert!(record.residual_history.len() <= 6, "{record:?}");
    }
    // Permanent state after unloading: every element carries alpha = eps_p,
    // and the free lateral faces contracted plastically by eps_p / 2.
    for state in &solution.element_states {
        let J2State { alpha, .. } = state;
        // Newton's 1e-10 force gate bounds the committed state to ~1e-9.
        assert!((alpha - plastic).abs() < 1e-8 * plastic, "alpha {alpha} vs {plastic}");
    }
    let corner = nodes
        .iter()
        .position(|x| x.iter().zip(&size).all(|(a, b)| (a - b).abs() < 1e-12))
        .unwrap();
    let lateral = solution.displacement_m[corner][1] / size[1];
    println!(
        "residual lateral strain {lateral:.12} vs {:.12}",
        -plastic / 2.0
    );
    assert!((lateral + plastic / 2.0).abs() < 1e-8 * plastic);
}

fn problem_for<'a, L: SmallStrainLaw>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    law: &'a L,
) -> SmallStrainTetProblem<'a, L> {
    SmallStrainTetProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: &[],
        nodal_forces_n: &[],
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}

#[test]
fn algorithmic_tangent_is_the_derivative_at_a_plastic_state() {
    let (nodes, tets) = box_mesh([3, 2, 2], [1.5, 1.0, 1.0]);
    let law = j2();
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::Uniform(&law),
        prescribed_m: &[],
        nodal_forces_n: &[],
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    // A committed state with some plastic history, then a non-uniform
    // displacement that drives every element further into flow.
    let states: Vec<J2State> = (0..tets.len())
        .map(|e| J2State {
            plastic_strain: [
                0.0005,
                -0.00025,
                -0.00025,
                0.0001 * (e % 3) as f64,
                0.0,
                0.0,
            ],
            alpha: 0.0005,
        })
        .collect();
    let u: Vec<f64> = nodes
        .iter()
        .flat_map(|x| {
            [
                0.004 * x[0] + 0.002 * x[1] * x[2],
                -0.001 * x[1],
                0.003 * x[0] * x[2],
            ]
        })
        .collect();
    let d: Vec<f64> = (0..u.len())
        .map(|i| ((i * 6151 % 97) as f64 / 97.0) - 0.5)
        .collect();
    let (k, plus, minus) = with_cx(|cx| {
        let eps = 1e-9;
        let up: Vec<f64> = u.iter().zip(&d).map(|(a, b)| a + eps * b).collect();
        let um: Vec<f64> = u.iter().zip(&d).map(|(a, b)| a - eps * b).collect();
        (
            problem.tangent_matrix(&u, &states, cx).unwrap(),
            problem.internal_force(&up, &states, cx).unwrap(),
            problem.internal_force(&um, &states, cx).unwrap(),
        )
    });
    // Plastic flow is active: the algorithmic tangent differs from the
    // elastic one at this state.
    let elastic = IsotropicElastic::new(E, NU, 0.05).unwrap();
    let k_elastic = with_cx(|cx| {
        SmallStrainTetProblem {
            laws: SmallStrainLawField::Uniform(&elastic),
            ..problem_for(&nodes, &tets, &elastic)
        }
        .tangent_matrix(&u, &vec![elastic.initial_state(); tets.len()], cx)
        .unwrap()
    });
    let mut gap = 0.0f64;
    for r in 0..k.nrows() {
        for c in 0..k.ncols() {
            gap = gap.max((k.get(r, c) - k_elastic.get(r, c)).abs());
        }
    }
    let mut kd = vec![0.0; u.len()];
    k.spmv(&d, &mut kd);
    let fd: Vec<f64> = plus
        .iter()
        .zip(&minus)
        .map(|(a, b)| (a - b) / 2e-9)
        .collect();
    let err = kd
        .iter()
        .zip(&fd)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt();
    let scale = kd.iter().map(|a| a * a).sum::<f64>().sqrt();
    println!(
        "max |K_alg - K_el| = {gap:e}; |K d - FD| / |K d| = {:e}",
        err / scale
    );
    assert!(gap > 1.0, "no plastic flow at the probe state");
    assert!(err / scale < 1e-5);
}

#[test]
fn elastic_law_tangent_equals_the_linear3_stiffness() {
    let (nodes, tets) = box_mesh([2, 2, 1], [1.0, 1.0, 0.5]);
    let law = IsotropicElastic::new(E, NU, 0.05).unwrap();
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::Uniform(&law),
        prescribed_m: &[],
        nodal_forces_n: &[],
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    let reference = TetElasticMaterial::try_new(1.0, E, NU, [9; 32]).unwrap();
    let (k, k_ref) = with_cx(|cx| {
        let states = vec![law.initial_state(); tets.len()];
        (
            problem
                .tangent_matrix(&vec![0.0; 3 * nodes.len()], &states, cx)
                .unwrap(),
            TetLinearElasticProblem {
                nodes_m: &nodes,
                tetrahedra: &tets,
                materials: TetMaterialField::Uniform(&reference),
                fixed_dofs: &[],
                budget: TetAssemblyBudget::standard(),
            }
            .assemble(cx)
            .unwrap()
            .stiffness,
        )
    });
    let (mut worst, mut scale) = (0.0f64, 0.0f64);
    for r in 0..k_ref.nrows() {
        for c in 0..k_ref.ncols() {
            worst = worst.max((k.get(r, c) - k_ref.get(r, c)).abs());
            scale = scale.max(k_ref.get(r, c).abs());
        }
    }
    println!("elastic tangent vs linear3: max |dK| = {worst:e} (scale {scale:e})");
    assert!(worst < 1e-12 * scale);
}

#[test]
fn cantilever_yields_progressively() {
    let size = [4.0, 1.0, 1.0];
    let (nodes, tets) = box_mesh([8, 2, 2], size);
    let law = j2();
    let mut prescribed = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        if x[0].abs() < 1e-12 {
            for i in 0..3 {
                prescribed.push((3 * node + i, 0.0));
            }
        }
        if (x[0] - size[0]).abs() < 1e-12 {
            prescribed.push((3 * node + 2, -1.0));
        }
    }
    // Tip displacement path (m): elastic, first yield, spreading plasticity.
    let path = [0.002, 0.004, 0.008, 0.012, 0.016, 0.020];
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::Uniform(&law),
        prescribed_m: &prescribed,
        nodal_forces_n: &[],
        load_path: &path,
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    let solution = with_cx(|cx| problem.solve(cx)).unwrap();
    let mut forces = Vec::new();
    let mut yielded = Vec::new();
    for record in &solution.increments {
        let tip: f64 = prescribed
            .iter()
            .zip(&record.reactions_n)
            .filter(|((dof, value), _)| dof % 3 == 2 && *value < 0.0)
            .map(|(_, r)| -r)
            .sum();
        forces.push(tip);
        println!(
            "tip {:.3}: force {tip:.6}, newly evolving {}, iterations {}",
            record.load_factor,
            record.evolving_elements,
            record.residual_history.len()
        );
    }
    for state in &solution.element_states {
        yielded.push(state.alpha > 0.0);
    }
    let secant: Vec<f64> = forces.iter().zip(&path).map(|(f, d)| f / d).collect();
    // Elastic first increment, then yielding: the secant stiffness falls and
    // the plastic zone is non-empty at the end.
    assert_eq!(solution.increments[0].evolving_elements, 0);
    assert!(secant.windows(2).skip(1).all(|w| w[1] < w[0]), "{secant:?}");
    assert!(yielded.iter().filter(|&&y| y).count() > 0);
    for record in &solution.increments {
        assert!(record.residual_history.len() <= 12, "{record:?}");
    }
}

fn continuation_prescribed(nodes: &[[f64; 3]], length: f64, drive_end: bool) -> Vec<(usize, f64)> {
    let mut prescribed = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        if x[0] == 0.0 {
            prescribed.push((3 * node, 0.0));
        }
        if drive_end && x[0] == length {
            prescribed.push((3 * node, length));
        }
        if x[1] == 0.0 {
            prescribed.push((3 * node + 1, 0.0));
        }
        if x[2] == 0.0 {
            prescribed.push((3 * node + 2, 0.0));
        }
    }
    prescribed
}

#[test]
fn split_loading_reversal_and_reloading_preserve_the_native_path() {
    let (nodes, tets) = box_mesh([2, 1, 1], [2.0, 1.0, 1.0]);
    let law = j2();
    let prescribed = continuation_prescribed(&nodes, 2.0, true);
    let path = [0.002, 0.004, 0.006, 0.0045, -0.002, 0.007];
    let problem = SmallStrainTetProblem {
        prescribed_m: &prescribed,
        load_path: &path,
        ..problem_for(&nodes, &tets, &law)
    };
    let whole = with_cx(|cx| problem.solve(cx)).unwrap();
    let prefix_problem = SmallStrainTetProblem {
        load_path: &path[..3],
        ..problem
    };
    let loaded = with_cx(|cx| prefix_problem.solve(cx)).unwrap();
    let accepted = loaded.clone();
    assert!(loaded.element_states.iter().all(|state| state.alpha > 0.0));
    let remainder = SmallStrainTetProblem {
        load_path: &path[3..],
        ..prefix_problem
    };
    let resumed = with_cx(|cx| {
        remainder.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.element_states,
            cx,
        )
    })
    .unwrap();
    assert_eq!(resumed.displacement_m, whole.displacement_m);
    assert_eq!(resumed.element_strain, whole.element_strain);
    assert_eq!(resumed.element_stress, whole.element_stress);
    assert_eq!(resumed.element_states, whole.element_states);
    let mut joined = loaded.increments.clone();
    joined.extend(resumed.increments.clone());
    assert_eq!(joined, whole.increments);
    assert_eq!(loaded, accepted);
}

#[test]
fn releasing_driven_constraints_retains_closed_form_permanent_strain() {
    let (nodes, tets) = box_mesh([2, 1, 1], [2.0, 1.0, 1.0]);
    let law = j2();
    let driven = continuation_prescribed(&nodes, 2.0, true);
    let peak = 0.006;
    let peak_path = [peak];
    let problem = SmallStrainTetProblem {
        prescribed_m: &driven,
        load_path: &peak_path,
        ..problem_for(&nodes, &tets, &law)
    };
    let loaded = with_cx(|cx| problem.solve(cx)).unwrap();
    let accepted = loaded.clone();
    let rollers = continuation_prescribed(&nodes, 2.0, false);
    let release = SmallStrainTetProblem {
        prescribed_m: &rollers,
        ..problem_for(&nodes, &tets, &law)
    };
    let unloaded = with_cx(|cx| {
        release.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.element_states,
            cx,
        )
    })
    .unwrap();
    // Independent uniaxial energetic-hardening solution; the driven end is
    // now free, so it returns to plastic strain rather than zero displacement.
    let peak_stress = E * (SIGMA_Y + H * peak) / (E + H);
    let plastic = peak - peak_stress / E;
    for (x, u) in nodes.iter().zip(&unloaded.displacement_m) {
        let expected = [plastic * x[0], -0.5 * plastic * x[1], -0.5 * plastic * x[2]];
        for (&actual, reference) in u.iter().zip(expected) {
            assert!((actual - reference).abs() < 1e-9, "{actual} != {reference}");
        }
    }
    for (old, next) in loaded.element_states.iter().zip(&unloaded.element_states) {
        assert!((next.alpha - plastic).abs() < 1e-9);
        assert!((next.alpha - old.alpha).abs() < 1e-9);
    }
    assert!(
        unloaded
            .element_stress
            .iter()
            .flatten()
            .all(|s| s.abs() < 1e-8)
    );
    assert_eq!(loaded, accepted);
}

#[test]
fn changed_element_yield_data_continue_history_at_fixed_displacement() {
    let (nodes, tets) = box_mesh([2, 1, 1], [2.0, 1.0, 1.0]);
    let law = j2();
    let driven = continuation_prescribed(&nodes, 2.0, true);
    let problem = SmallStrainTetProblem {
        prescribed_m: &driven,
        load_path: &[0.006],
        ..problem_for(&nodes, &tets, &law)
    };
    let loaded = with_cx(|cx| problem.solve(cx)).unwrap();
    let endpoint: Vec<_> = loaded
        .displacement_m
        .iter()
        .flatten()
        .copied()
        .enumerate()
        .collect();
    let softened = J2Plasticity::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        0.5 * SIGMA_Y,
        H,
    )
    .unwrap();
    let laws: Vec<_> = (0..tets.len())
        .map(|element| {
            if element % 2 == 0 {
                softened.clone()
            } else {
                law.clone()
            }
        })
        .collect();
    let changed = SmallStrainTetProblem {
        laws: SmallStrainLawField::PerElement(&laws),
        prescribed_m: &endpoint,
        ..problem_for(&nodes, &tets, &law)
    };
    let next = with_cx(|cx| {
        changed.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.element_states,
            cx,
        )
    })
    .unwrap();
    assert_eq!(next.displacement_m, loaded.displacement_m);
    // At fixed total strain, reducing only yield releases deviatoric elastic
    // storage. The radial-return consistency equation gives this delta alpha.
    let shear = E / (2.0 * (1.0 + NU));
    let delta_alpha = 0.5 * SIGMA_Y / (3.0 * shear + H);
    for (element, (old, state)) in loaded
        .element_states
        .iter()
        .zip(&next.element_states)
        .enumerate()
    {
        let expected = old.alpha + if element % 2 == 0 { delta_alpha } else { 0.0 };
        assert!((state.alpha - expected).abs() < 1e-9);
        if element % 2 == 0 {
            assert!(state.alpha > old.alpha);
        }
    }
}

#[test]
fn continuation_refusal_does_not_consume_the_accepted_history() {
    let (nodes, tets) = box_mesh([2, 1, 1], [2.0, 1.0, 1.0]);
    let law = j2();
    let prescribed = continuation_prescribed(&nodes, 2.0, true);
    let problem = SmallStrainTetProblem {
        prescribed_m: &prescribed,
        load_path: &[0.006],
        ..problem_for(&nodes, &tets, &law)
    };
    let loaded = with_cx(|cx| problem.solve(cx)).unwrap();
    let accepted = loaded.clone();
    let failing = SmallStrainTetProblem {
        load_path: &[0.006, 0.012],
        settings: IncrementSettings {
            max_newton_iterations: 1,
            ..IncrementSettings::default()
        },
        ..problem
    };
    let result = with_cx(|cx| {
        failing.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.element_states,
            cx,
        )
    });
    assert!(matches!(
        result,
        Err(fs_solid::SmallStrainTetError::NewtonStalled { increment: 1, .. })
    ));
    assert!(matches!(
        with_cx(|cx| failing.solve_from(
            loaded.displacement_m[..nodes.len() - 1].as_flattened(),
            &loaded.element_states,
            cx,
        )),
        Err(fs_solid::SmallStrainTetError::InvalidInput { .. })
    ));
    assert!(matches!(
        with_cx(|cx| failing.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.element_states[..tets.len() - 1],
            cx,
        )),
        Err(fs_solid::SmallStrainTetError::InvalidInput { .. })
    ));
    let mut nonfinite = loaded.displacement_m.clone();
    nonfinite[0][0] = f64::NAN;
    assert!(matches!(
        with_cx(|cx| failing.solve_from(nonfinite.as_flattened(), &loaded.element_states, cx,)),
        Err(fs_solid::SmallStrainTetError::InvalidInput { .. })
    ));
    assert_eq!(loaded, accepted);
    let retry = SmallStrainTetProblem {
        settings: IncrementSettings::default(),
        ..failing
    };
    assert!(
        with_cx(|cx| {
            retry.solve_from(
                loaded.displacement_m.as_flattened(),
                &loaded.element_states,
                cx,
            )
        })
        .is_ok()
    );
    assert_eq!(loaded, accepted);
}

#[test]
fn overstress_solution_publishes_one_physical_update_per_increment() {
    use fs_material::plastic::rate::LinearOverstressJ2;

    let nodes = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let tets = [[0, 1, 2, 3]];
    let dt = 0.1;
    let viscosity = 100.0;
    let law = LinearOverstressJ2::new(j2(), viscosity, dt).unwrap();
    let strain = [0.006, 0.0, 0.0, 0.0, 0.0, 0.0];
    let displacement: Vec<_> = nodes.iter().flat_map(|x| [strain[0] * x[0], 0.0, 0.0]).collect();
    let prescribed: Vec<_> = displacement.iter().copied().enumerate().collect();
    let problem = SmallStrainTetProblem {
        prescribed_m: &prescribed,
        ..problem_for(&nodes, &tets, &law)
    };
    let mut state = law.initial_state();
    let mut old_u = vec![0.0; 12];
    let shear = E / (2.0 * (1.0 + NU));
    for step in 0..3 {
        let old_alpha = state.alpha;
        let old_plastic = state.plastic_strain;
        let result = with_cx(|cx| problem.solve_from(&old_u, &[state.clone()], cx)).unwrap();
        let next = &result.element_states[0];
        // Independently reconstruct the actual endpoint C:(eps-eps_p),
        // which must agree with both reported stress and support reactions.
        let elastic_strain = std::array::from_fn(|i| strain[i] - next.plastic_strain[i]);
        let expected = law.backbone().elastic.stress(&elastic_strain, &());
        for (actual, expected) in result.element_stress[0].iter().zip(expected) {
            assert!((actual - expected).abs() < 2e-12, "step {step}: {actual} != {expected}");
        }
        assert!((result.increments[0].reactions_n[3] - expected[0] / 6.0).abs() < 2e-12);
        // Scalar backward-Euler Perzyna equation, including further stress
        // relaxation at fixed displacement on the second and third intervals.
        let trial_q = 2.0 * shear * (strain[0] - 1.5 * old_plastic[0]);
        let delta = (trial_q - SIGMA_Y - H * old_alpha).max(0.0)
            / (3.0 * shear + H + viscosity / dt);
        assert!((next.alpha - old_alpha - delta).abs() < 2e-15);
        assert!(delta > 0.0);
        state = next.clone();
        old_u = result.displacement_m.as_flattened().to_vec();
    }
}
