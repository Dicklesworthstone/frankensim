//! Analytical initial-strain checks through the original native solid kernels.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::{eigenstrain::*, *};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}
const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

#[test]
fn generic_and_thermal_sources_share_the_same_numeric_mechanical_path() {
    with_cx(|cx| {
        let material = TetElasticMaterial::try_new(700., 1e6, 0.25, [1; 32]).unwrap();
        let p = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &[0, 1, 2, 4, 5, 8],
            budget: TetAssemblyBudget::standard(),
        };
        let eps = [0.001, 0.002, -0.0005, 0., 0., 0.];
        let generic =
            TetEigenstrainState::try_new(eps, [1; 32], ContentHash([2; 32]), 0.02).unwrap();
        let thermal =
            TetThermalStrainState::try_new(310., 300., eps, [1; 32], ContentHash([3; 32]), 0.02)
                .unwrap();
        let g = TetEigenstrainField::Uniform(&generic);
        let t = TetThermalStrainField::Uniform(&thermal);
        let gl = p.assemble_eigenstrain_load(g, cx).unwrap();
        let tl = p.assemble_thermal_load(t, cx).unwrap();
        assert_eq!(gl.full_equivalent_force_n, tl.full_equivalent_force_n);
        assert_ne!(gl.identity, tl.identity);
        let gs = p
            .solve_eigenstrain_displacement(g, TetStaticSolveConfig::default(), cx)
            .unwrap();
        let ts = p
            .solve_thermal_displacement(t, TetStaticSolveConfig::default(), cx)
            .unwrap();
        assert_eq!(gs.displacement_m(), ts.displacement_m());
        for (x, u) in NODES.iter().zip(gs.displacement_m()) {
            for k in 0..3 {
                near(u[k], eps[k] * x[k], 1e-13);
            }
        }
        let stress = p
            .recover_eigenstrain_stress(&gs, g, ContentHash([4; 32]), 0.02, cx)
            .unwrap();
        for s in stress.elements()[0].stress_mandel_pa() {
            near(*s, 0., 1e-8);
        }
        let mesh = p.update_geometry_from_displacement(&gs, 0.1, cx).unwrap();
        near(
            mesh.element_volumes_m3()[0],
            (1. + eps[0]) * (1. + eps[1]) * (1. + eps[2]) / 6.,
            1e-14,
        );
        // Existing modal/stability owner receives the same physical wet mass.
        let tangent = p
            .assemble_eigenstrain_tangent(&gs, g, ContentHash([4; 32]), 0.02, cx)
            .unwrap();
        near(tangent.pencil().total_mass_kg, 700. / 6., 1e-12);
    });
}

#[test]
fn constrained_swelling_recovers_compression_and_energy_without_fictitious_temperature() {
    with_cx(|cx| {
        let material = TetElasticMaterial::try_new(700., 1000., 0.25, [1; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let p = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
        };
        let strain = TetEigenstrainState::try_new(
            [0.01, 0.01, 0.01, 0., 0., 0.],
            [1; 32],
            ContentHash([2; 32]),
            0.02,
        )
        .unwrap();
        let field = TetEigenstrainField::Uniform(&strain);
        let solution = p
            .solve_eigenstrain_displacement(field, TetStaticSolveConfig::default(), cx)
            .unwrap();
        assert_eq!(solution.displacement_m(), &[[0.; 3]; 4]);
        let result = p
            .recover_eigenstrain_stress(&solution, field, ContentHash([4; 32]), 0.02, cx)
            .unwrap();
        for &s in &result.elements()[0].stress_mandel_pa()[..3] {
            near(s, -20., 1e-12);
        }
        near(result.elastic_energy_j(), 0.05, 1e-14);
        for k in 0..3 {
            near(
                result.nodal_internal_force_n().iter().map(|v| v[k]).sum(),
                0.,
                1e-12,
            );
        }
        let foreign =
            TetEigenstrainState::try_new([0.; 6], [8; 32], ContentHash([2; 32]), 0.02).unwrap();
        assert!(matches!(
            p.assemble_eigenstrain_load(TetEigenstrainField::Uniform(&foreign), cx),
            Err(TetElasticError::EigenstrainElasticStateMismatch { element: 0 })
        ));
        let changed =
            TetEigenstrainState::try_new([0.005; 6], [1; 32], ContentHash([2; 32]), 0.02).unwrap();
        assert!(matches!(
            p.recover_eigenstrain_stress(
                &solution,
                TetEigenstrainField::Uniform(&changed),
                ContentHash([4; 32]),
                0.02,
                cx
            ),
            Err(TetElasticError::StressRecoveryInputMismatch)
        ));
    });
}

#[test]
fn principal_axes_rotate_and_unsupported_states_refuse() {
    let material = TetElasticMaterial::try_new(700., 1000., 0.25, [1; 32]).unwrap();
    let q = [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]];
    let strain = TetEigenstrainState::try_from_principal_strains(
        &material,
        [0.01, 0.02, -0.03],
        q,
        ContentHash([2; 32]),
        ContentHash([3; 32]),
        ContentHash([4; 32]),
        0.04,
    )
    .unwrap();
    assert_eq!(
        strain.free_strain_mandel(),
        &[0.02, 0.01, -0.03, 0., 0., 0.]
    );
    assert!(
        TetEigenstrainState::try_from_principal_strains(
            &material,
            [0.01; 3],
            q,
            ContentHash([2; 32]),
            ContentHash([2; 32]),
            ContentHash([4; 32]),
            0.04
        )
        .is_err()
    );
    assert!(TetEigenstrainState::try_new([0.1; 6], [1; 32], ContentHash([4; 32]), 0.04).is_err());
    let gate = CancelGate::new();
    gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let p = TetLinearElasticProblem {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &[0, 1, 2, 4, 5, 8],
            budget: TetAssemblyBudget::standard(),
        };
        assert!(matches!(
            p.solve_eigenstrain_displacement(
                TetEigenstrainField::Uniform(&strain),
                TetStaticSolveConfig::default(),
                &cx
            ),
            Err(TetElasticError::Cancelled)
        ));
    });
}
