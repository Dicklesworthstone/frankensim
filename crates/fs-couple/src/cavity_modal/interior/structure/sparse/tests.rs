use super::super::time::{StructuralInteriorTimeError, StructuralInteriorTimeLimits};
use super::super::{ConformingTetFace, StructuralInterior};
use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_feec::fixtures::single_tet;
use fs_modal::harmonic::HarmonicResponseControls;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn matrix(n: usize, data: &[f64]) -> Csr {
    let mut result = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            if data[i * n + j] != 0.0 {
                result.push(i, j, data[i * n + j]);
            }
        }
    }
    result.assemble()
}
fn dense_controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 128,
        maximum_payload_bytes: 16 * 1024 * 1024,
        maximum_relative_residual: 1e-10,
    }
}

#[test]
fn full_heterogeneous_pressure_trace_refines_to_existing_exact_coupled_owner() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 2, medium).unwrap();
    let interior = HeterogeneousHelmholtzFe::new(
        &geometry,
        5,
        &mut |_, p: [f64; 3]| {
            Ok::<_, AcousticFeError>(Medium {
                density_kg_m3: 1.1 + 0.2 * p[0],
                sound_speed_m_s: 3.0 + 0.4 * p[1],
            })
        },
        &mut || false,
    )
    .unwrap();
    let mass = matrix(2, &[2.0, 0.4, 0.4, 1.2]);
    let stiffness = matrix(2, &[70.0, -13.0, -13.0, 40.0]);
    let damping = matrix(2, &[0.7, 0.1, 0.1, 0.5]);
    let traces = [
        vec![([0, 1, 2], [0.7, -0.3, 0.9])],
        vec![([2, 0, 1], [1.1, -0.2, 0.8])],
    ];
    let walls = [ResistiveBoundary {
        faces: vec![[0, 1, 3]],
        admittance_m_pa_s: 0.02,
    }];
    let matrices = || StructuralMatrices {
        mass: &mass,
        stiffness: &stiffness,
        damping: Some(&damping),
    };
    with_cx(|cx| {
        let original =
            StructuralInterior::new(&interior, matrices(), &traces, &walls, dense_controls(), cx)
                .unwrap();
        let total = 0.08;
        let exact = original
            .prepare_exact_time(total, StructuralInteriorTimeLimits::default(), cx)
            .unwrap();
        let u = [0.02, -0.01];
        let v = [0.01, 0.03];
        let psi = (0..interior.ndof())
            .map(|i| 1e-4 * (i as f64 - 2.0))
            .collect::<Vec<_>>();
        let p = (0..interior.ndof())
            .map(|i| 0.001 * (i as f64 + 1.0))
            .collect::<Vec<_>>();
        let native_source = [([1, 2, 3], [0.001, -0.002, 0.003])];
        let source = geometry.boundary_linear_load(&native_source).unwrap();
        let force = [0.3, -0.2];
        let initial = exact.state(0.0, &u, &v, &psi, &p, cx).unwrap();
        let reference = exact.step(&initial, &force, &source, cx).unwrap();
        let mut errors = Vec::new();
        for count in [8, 16, 32] {
            let sparse = PreparedSparseStructuralInterior::new(
                &interior,
                matrices(),
                &traces,
                &walls,
                total / count as f64,
                SparseStructuralTimeLimits::default(),
                cx,
            )
            .unwrap();
            assert_eq!(sparse.pressure_dofs(), geometry.ndof());
            assert_eq!(
                sparse.project_inward_velocity(&native_source, cx).unwrap(),
                source
            );
            for a in 0..interior.ndof() {
                for s in 0..2 {
                    assert_eq!(sparse.coupling().get(a, s), original.coupling()[a * 2 + s]);
                }
            }
            let mut state = sparse.state(0.0, &u, &v, &psi, &p, cx).unwrap();
            for _ in 0..count {
                let step = sparse.step(&state, &force, &source, cx).unwrap();
                assert!(step.work.structural_dissipation_j >= 0.0);
                assert!(step.work.wall_dissipation_j >= 0.0);
                assert!(step.work.structural_balance_residual_j().abs() < 1e-11);
                assert!(step.work.acoustic_balance_residual_j().abs() < 1e-11);
                state = step.state;
            }
            let mut error = 0.0;
            for (a, b) in [
                (
                    state.structural_displacement(),
                    reference.state.structural_displacement(),
                ),
                (
                    state.structural_velocity(),
                    reference.state.structural_velocity(),
                ),
                (
                    state.pressure_potential_pa_s(),
                    reference.state.pressure_potential_pa_s(),
                ),
                (state.pressure_pa(), reference.state.pressure_pa()),
            ] {
                error += a.iter().zip(b).map(|(&x, &y)| (x - y).abs()).sum::<f64>();
            }
            errors.push(error);
        }
        assert!(
            errors[0] / errors[1] > 3.6 && errors[1] / errors[2] > 3.6,
            "{errors:?}"
        );
    });
}

#[test]
fn original_tetrahedral_solid_and_high_order_pressure_run_beyond_dense_time_cap() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 4, medium).unwrap();
    let interior = HeterogeneousHelmholtzFe::new(
        &geometry,
        6,
        &mut |_, _| Ok::<_, AcousticFeError>(medium),
        &mut || false,
    )
    .unwrap();
    let solid_points = [points[0], points[1], points[2], [0.0, 0.0, -1.0]];
    let tetrahedra = [[0, 1, 2, 3]];
    let material = TetElasticMaterial::try_new(6.0, 120.0, 0.25, [7; 32]).unwrap();
    let problem = TetLinearElasticProblem {
        nodes_m: &solid_points,
        tetrahedra: &tetrahedra,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &[],
        budget: TetAssemblyBudget::standard(),
    };
    with_cx(|cx| {
        let assembly = problem.assemble(cx).unwrap();
        let interface = [ConformingTetFace {
            fluid: [0, 1, 2],
            solid: [0, 1, 2],
        }];
        let solid = || ConformingTetStructure {
            assembly: &assembly,
            nodes_m: &solid_points,
            tetrahedra: &tetrahedra,
            interface: &interface,
        };
        let original = StructuralInterior::from_tet_elastic(
            &interior,
            solid(),
            None,
            &[],
            dense_controls(),
            cx,
        )
        .unwrap();
        assert!(assembly.free_dofs.len() + interior.ndof() > 42);
        assert!(matches!(
            original.prepare_exact_time(0.0005, StructuralInteriorTimeLimits::default(), cx),
            Err(StructuralInteriorTimeError::Limit { .. })
        ));
        let sparse = PreparedSparseStructuralInterior::from_tet_elastic(
            &interior,
            solid(),
            None,
            &[],
            0.0005,
            SparseStructuralTimeLimits::default(),
            cx,
        )
        .unwrap();
        let ns = assembly.free_dofs.len();
        let na = interior.ndof();
        assert_eq!(sparse.structural_dofs(), ns);
        assert_eq!(sparse.pressure_dofs(), na);
        for a in 0..na {
            for s in 0..ns {
                assert_eq!(sparse.coupling().get(a, s), original.coupling()[a * ns + s]);
            }
        }
        let mut force = vec![0.0; ns];
        force[2] = 0.2;
        let source = vec![0.0; na];
        let mut state = sparse.zero_state();
        let mut supplied = 0.0;
        for _ in 0..8 {
            let next = sparse.step(&state, &force, &source, cx).unwrap();
            supplied += next.work.structural_input_work_j;
            assert!(next.work.balance_residual_j().abs() < 1e-11);
            state = next.state;
        }
        assert!(state.pressure_pa().iter().any(|v| v.abs() > 1e-8));
        assert!((sparse.energy(&state, cx).unwrap().total_j() - supplied).abs() < 1e-11);
    });
}

#[test]
fn sparse_public_source_admission_keeps_reciprocal_scope_explicit() {
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = HeterogeneousHelmholtzFe::new(
        &geometry,
        3,
        &mut |_, _| Ok::<_, AcousticFeError>(medium),
        &mut || false,
    )
    .unwrap();
    let mass = Csr::identity(2);
    let stiffness = Csr::identity(2);
    let damping = matrix(2, &[0.1, 0.2, -0.2, 0.1]);
    let traces = [vec![([0, 1, 2], [0.0, 0.0, 1.0])], Vec::new()];
    with_cx(|cx| {
        assert!(
            PreparedSparseStructuralInterior::new(
                &interior,
                StructuralMatrices {
                    mass: &mass,
                    stiffness: &stiffness,
                    damping: Some(&damping)
                },
                &traces,
                &[],
                0.01,
                SparseStructuralTimeLimits::default(),
                cx
            )
            .is_err()
        );
        let mut limits = SparseStructuralTimeLimits::default();
        limits.maximum_trace_terms = 1;
        let walls = [ResistiveBoundary {
            faces: vec![[0, 1, 3]],
            admittance_m_pa_s: 0.1,
        }];
        assert!(matches!(
            PreparedSparseStructuralInterior::new(
                &interior,
                StructuralMatrices {
                    mass: &mass,
                    stiffness: &stiffness,
                    damping: None
                },
                &traces,
                &walls,
                0.01,
                limits,
                cx
            ),
            Err(SparseStructuralInteriorError::Time(
                SparseStructuralTimeError::Limit { .. }
            ))
        ));
        assert!(
            PreparedSparseStructuralInterior::new(
                &interior,
                StructuralMatrices {
                    mass: &mass,
                    stiffness: &stiffness,
                    damping: None
                },
                &[],
                &[],
                0.01,
                SparseStructuralTimeLimits::default(),
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn sparse_original_maxwell_histories_refine_to_existing_exact_coupled_owner() {
    use super::super::time::relaxation::StructuralInteriorRelaxationTimeLimits;
    let (complex, points) = single_tet();
    let medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 3.0,
    };
    let geometry = HelmholtzFe::new(&complex, &points, 1, medium).unwrap();
    let interior = HeterogeneousHelmholtzFe::new(
        &geometry,
        3,
        &mut |_, _| Ok::<_, AcousticFeError>(medium),
        &mut || false,
    )
    .unwrap();
    let mass = matrix(2, &[2.0, 0.4, 0.4, 1.2]);
    let stiffness = matrix(2, &[70.0, -13.0, -13.0, 40.0]);
    let damping = matrix(2, &[0.7, 0.1, 0.1, 0.5]);
    let traces = [
        vec![([0, 1, 2], [0.7, -0.3, 0.9])],
        vec![([2, 0, 1], [1.1, -0.2, 0.8])],
    ];
    let walls = [ResistiveBoundary {
        faces: vec![[0, 1, 3]],
        admittance_m_pa_s: 0.02,
    }];
    let matrices = || StructuralMatrices {
        mass: &mass,
        stiffness: &stiffness,
        damping: Some(&damping),
    };
    let n = 2 + interior.ndof();
    let mut solid = vec![0.0; n];
    solid[..2].copy_from_slice(&[0.7, -0.3]);
    let mut acoustic = vec![0.0; n];
    for (i, b) in acoustic[2..].iter_mut().enumerate() {
        *b = 0.02 * (i as f64 - 1.0);
    }
    let branches = [
        MaxwellDampingBranch {
            strain: solid,
            stiffness: 2.0,
            relaxation_time_s: 0.01,
        },
        MaxwellDampingBranch {
            strain: acoustic,
            stiffness: 0.8,
            relaxation_time_s: 0.02,
        },
        MaxwellDampingBranch {
            strain: vec![0.0; n],
            stiffness: 0.4,
            relaxation_time_s: 0.03,
        },
    ];
    let stresses = [0.03, 0.001, -0.04];
    let u = [0.02, -0.01];
    let v = [0.01, 0.03];
    let psi = vec![0.001; interior.ndof()];
    let p = (0..interior.ndof())
        .map(|i| 0.001 * (i + 1) as f64)
        .collect::<Vec<_>>();
    let source = geometry
        .boundary_linear_load(&[([1, 2, 3], [0.001, -0.002, 0.003])])
        .unwrap();
    let force = [0.02, -0.03];
    with_cx(|cx| {
        let original =
            StructuralInterior::new(&interior, matrices(), &traces, &walls, dense_controls(), cx)
                .unwrap();
        let total = 0.04;
        let exact = original
            .prepare_relaxation_exact_time(
                &branches,
                total,
                StructuralInteriorRelaxationTimeLimits::default(),
                cx,
            )
            .unwrap();
        let initial = exact.state(0.0, &u, &v, &psi, &p, &stresses, cx).unwrap();
        let reference = exact.step(&initial, &force, &source, cx).unwrap();
        let mut errors = Vec::new();
        for count in [8, 16, 32] {
            let sparse = PreparedSparseStructuralInterior::new_relaxing(
                &interior,
                matrices(),
                &traces,
                &walls,
                &branches,
                total / count as f64,
                SparseStructuralTimeLimits::default(),
                SparseStructuralRelaxationLimits::default(),
                cx,
            )
            .unwrap();
            assert_eq!(sparse.branch_count(), branches.len());
            let mut state = sparse
                .state_with_stresses(0.0, &u, &v, &psi, &p, &stresses, cx)
                .unwrap();
            for _ in 0..count {
                let step = sparse.step(&state, &force, &source, cx).unwrap();
                assert_eq!(step.memory.len(), branches.len());
                assert!(step.memory[2].dissipation_j > 0.0);
                assert_eq!(step.memory[2].structural_work_j, 0.0);
                assert_eq!(step.memory[2].acoustic_work_j, 0.0);
                state = step.state;
            }
            let mut error = 0.0;
            for (a, b) in [
                (
                    state.structural_displacement(),
                    reference.state.structural_displacement(),
                ),
                (
                    state.structural_velocity(),
                    reference.state.structural_velocity(),
                ),
                (
                    state.pressure_potential_pa_s(),
                    reference.state.pressure_potential_pa_s(),
                ),
                (state.pressure_pa(), reference.state.pressure_pa()),
                (state.branch_stresses(), reference.state.branch_stresses()),
            ] {
                error += a.iter().zip(b).map(|(&x, &y)| (x - y).abs()).sum::<f64>();
            }
            errors.push(error);
        }
        assert!(
            errors[0] / errors[1] > 3.5 && errors[1] / errors[2] > 3.5,
            "{errors:?}"
        );
    });
}
