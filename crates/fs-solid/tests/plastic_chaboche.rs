//! Nonproportional recovery must reach actual continuum equilibrium without
//! symmetrizing the constitutive tangent. Material oracles live in fs-material;
//! these tests independently apply uniform physical surface tractions.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2};
use fs_material::plastic::reduction::{PlaneStressReduction, ReducedState, ReductionConfig};
use fs_material::plastic::voce::VoceJ2;
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem, SmallStrainTetProblem, TetAssemblyBudget,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 9,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}
fn law() -> ChabocheJ2<2> {
    ChabocheJ2::new(
        VoceJ2::new(
            J2Plasticity::new_mixed(
                IsotropicElastic::new(1200.0, 0.3, 0.05).unwrap(),
                1.0,
                30.0,
                20.0,
            )
            .unwrap(),
            0.6,
            80.0,
        )
        .unwrap(),
        [
            BackstressBranch::new(300.0, 250.0).unwrap(),
            BackstressBranch::new(80.0, 30.0).unwrap(),
        ],
    )
    .unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.14e} != {b:.14e}, allowance {tol}"
    );
}
fn tet_displacement(e: Voigt) -> [f64; 12] {
    [
        0.0,
        0.0,
        0.0,
        e[0],
        0.0,
        0.0,
        2.0 * e[3],
        e[1],
        0.0,
        2.0 * e[5],
        2.0 * e[4],
        e[2],
    ]
}

#[test]
fn nonproportional_tet_uses_full_nonsymmetric_tangent_and_reaches_force_equilibrium() {
    let law = law();
    assert!(!law.admissibility().tangent_symmetric);
    let nodes = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let old_strain = [0.007, -0.002, -0.001, 0.0015, 0.0, 0.0];
    let target = [0.003, 0.006, -0.002, -0.004, 0.001, 0.002];
    let old = law
        .response(&old_strain, &law.initial_state())
        .unwrap()
        .state;
    let expected = law.response(&target, &old).unwrap();
    let s = expected.stress_pa;
    // Gradients are unit Cartesian vectors and volume=1/6. These are the
    // consistent nodal forces from the physical uniform boundary traction.
    let loads = [
        (3, s[0] / 6.0),
        (6, s[3] / 6.0),
        (7, s[1] / 6.0),
        (9, s[5] / 6.0),
        (10, s[4] / 6.0),
        (11, s[2] / 6.0),
    ];
    let supports = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
    let mut problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &[[0, 1, 2, 3]],
        laws: SmallStrainLawField::Uniform(&law),
        prescribed_m: &supports,
        nodal_forces_n: &loads,
        load_path: &[1.0],
        budget: TetAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    };
    let incoming = [old.clone()];
    let u = tet_displacement(target);
    with_cx(|cx| {
        let matrix = problem.tangent_matrix(&u, &incoming, cx).unwrap();
        let mut asymmetry = 0.0_f64;
        for &(i, _) in &loads {
            for &(j, _) in &loads {
                asymmetry = asymmetry.max((matrix.get(i, j) - matrix.get(j, i)).abs());
            }
        }
        assert!(
            asymmetry > 0.1,
            "fixture must exercise a nonsymmetric free tangent"
        );
        for &(column, _) in &loads {
            let h = 1e-7;
            let (mut plus, mut minus) = (u, u);
            plus[column] += h;
            minus[column] -= h;
            let fp = problem.internal_force(&plus, &incoming, cx).unwrap();
            let fm = problem.internal_force(&minus, &incoming, cx).unwrap();
            for row in 0..12 {
                close(
                    matrix.get(row, column),
                    (fp[row] - fm[row]) / (2.0 * h),
                    2e-5,
                );
            }
        }
        let result = problem
            .solve_from(&tet_displacement(old_strain), &incoming, cx)
            .unwrap();
        for (a, b) in result.displacement_m.iter().flatten().zip(u) {
            close(*a, b, 2e-9);
        }
        for i in 0..6 {
            close(result.element_stress[0][i], s[i], 2e-8);
        }
        assert!(result.increments[0].residual_history.len() > 1);
        assert_eq!(result.increments[0].indefinite_fallbacks, 0);
        assert!(result.element_states[0].plastic.alpha > old.plastic.alpha);
        // Exhausting a linear budget publishes neither a displacement nor history.
        problem.settings.max_linear_iterations = 0;
        assert!(
            problem
                .solve_from(&tet_displacement(old_strain), &incoming, cx)
                .is_err()
        );
        assert_eq!(incoming[0], old);
    });
}

#[test]
fn nonproportional_plane_stress_and_bbar_strain_preserve_the_general_jacobian() {
    let law = law();
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let config = ReductionConfig {
        max_iterations: 40,
        stress_tolerance_pa: 1e-11,
        strain_limit: 0.04,
    };
    let reduction = PlaneStressReduction::new(&law, config).unwrap();
    let old_in_plane = [0.007, -0.002, 0.0015];
    let target = [0.003, 0.006, -0.004];
    let displacement = |e: [f64; 3]| -> Vec<f64> {
        mesh.nodes
            .iter()
            .flat_map(|x| [e[0] * x[0] + 2.0 * e[2] * x[1], e[1] * x[1]])
            .collect()
    };
    for (plane, formulation) in [
        (PlaneKind::Stress, Formulation::Standard),
        (PlaneKind::Strain, Formulation::BBar),
    ] {
        let (old, stress) = if plane == PlaneKind::Stress {
            let old = reduction
                .evaluate(old_in_plane, &reduction.initial_state())
                .unwrap()
                .state;
            let stress = reduction.evaluate(target, &old).unwrap().stress_pa;
            (old, stress)
        } else {
            let old_strain = [
                old_in_plane[0],
                old_in_plane[1],
                0.0,
                old_in_plane[2],
                0.0,
                0.0,
            ];
            let history = law
                .response(&old_strain, &law.initial_state())
                .unwrap()
                .state;
            let s = law
                .response(&[target[0], target[1], 0.0, target[2], 0.0, 0.0], &history)
                .unwrap()
                .stress_pa;
            (
                ReducedState {
                    strain: old_strain,
                    history,
                },
                [s[0], s[1], s[3]],
            )
        };
        let [sx, sy, tau] = stress;
        let tractions = [
            -sx - tau,
            -tau - sy,
            sx - tau,
            tau - sy,
            -sx + tau,
            -tau + sy,
            sx + tau,
            tau + sy,
        ];
        let loads: Vec<_> = tractions
            .iter()
            .enumerate()
            .map(|(i, f)| (i, 0.5 * f))
            .collect();
        let problem = SmallStrainPlaneProblem {
            mesh: &mesh,
            laws: SmallStrainLawField::Uniform(&law),
            plane,
            formulation,
            thickness_m: 1.0,
            reduction: config,
            rotation_limit_rad: 0.04,
            prescribed_m: &[(0, 0.0), (1, 0.0), (3, 0.0)],
            nodal_forces_n: &loads,
            load_path: &[1.0],
            budget: PlaneAssemblyBudget::default(),
            settings: IncrementSettings::default(),
        };
        let incoming = vec![old.clone(); 4];
        let u = displacement(target);
        with_cx(|cx| {
            let matrix = problem.tangent_matrix(&u, &incoming, cx).unwrap();
            let mut asymmetry = 0.0_f64;
            for i in [2, 4, 5, 6, 7] {
                for j in [2, 4, 5, 6, 7] {
                    asymmetry = asymmetry.max((matrix.get(i, j) - matrix.get(j, i)).abs());
                }
            }
            assert!(
                asymmetry > 0.1,
                "{plane:?} must retain a nonsymmetric tangent"
            );
            for column in [2, 4, 5, 6, 7] {
                let h = 1e-7;
                let (mut plus, mut minus) = (u.clone(), u.clone());
                plus[column] += h;
                minus[column] -= h;
                let fp = problem.internal_force(&plus, &incoming, cx).unwrap();
                let fm = problem.internal_force(&minus, &incoming, cx).unwrap();
                for row in 0..8 {
                    close(
                        matrix.get(row, column),
                        (fp[row] - fm[row]) / (2.0 * h),
                        3e-5,
                    );
                }
            }
            let result = problem
                .solve_from(&displacement(old_in_plane), &incoming, cx)
                .unwrap();
            for (a, b) in result.displacement_m.iter().flatten().zip(&u) {
                close(*a, *b, 2e-9);
            }
            assert_eq!(incoming, vec![old.clone(); 4]);
        });
    }
}
