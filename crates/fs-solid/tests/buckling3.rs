//! Physical reference loads enter the existing 3-D initial-stress owner.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::SolidError;
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem, TetMaterialField,
    TetStaticSolveConfig, TetThermalStrainField, TetThermalStrainState,
};
use fs_solid::prestress::PrestressedTetAssembly;

fn with_cx<R>(cancelled: bool, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn strain(material: &TetElasticMaterial, value: [f64; 6]) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        310.0,
        300.0,
        value,
        material.material_state_identity,
        ContentHash([32; 32]),
        0.02,
    )
    .unwrap()
}

fn tets(a: f64, cx: &Cx<'_>) -> PrestressedTetAssembly {
    let nodes = [
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [0., 0., -1.],
    ];
    let cells = [[0, 1, 2, 3], [0, 2, 1, 4]];
    let fixed: Vec<_> = (0..15).filter(|&dof| dof != 4).collect();
    let material = TetElasticMaterial::try_new(6.0, 1000.0, 0.25, [31; 32]).unwrap();
    let thermal = strain(&material, [a, a, a, 0., 0., 0.]);
    let field = TetThermalStrainField::Uniform(&thermal);
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &cells,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
    };
    let solution = problem
        .solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
        .unwrap();
    problem
        .assemble_thermal_tangent(&solution, field, ContentHash([33; 32]), 0.02, cx)
        .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-9 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn source_bound_tetrahedra_return_total_load_and_original_xyz_mode() {
    with_cx(false, |cx| {
        let reference = tets(0.01, cx);
        let mass = reference.pencil().mass.to_dense();
        let result = reference.proportional_buckling(1, 100, cx).unwrap();
        // Both tets are hydrostatically restrained: sigma = -3 K a I.
        // Only node 1's y displacement is free: K_e = V G, K_g = -V 3 K a.
        // lambda = G / (3 K a) = 20, not 19 from double-counting K_g.
        close(result.result().loads[0], 20.0);
        assert_eq!(result.reference_assembly_identity(), reference.identity());
        assert_eq!(result.result().dof_map.len(), 15);
        let mode = result.nodal_mode(0).unwrap();
        assert_eq!(mode.len(), 5);
        for dof in 0..15 {
            if dof != 4 {
                assert_eq!(mode[dof / 3][dof % 3], 0.0);
            }
        }
        close(
            -reference.geometric_stiffness().get(0, 0) * mode[1][1].powi(2),
            1.0,
        );
        assert_eq!(reference.pencil().mass.to_dense(), mass);
        assert!(result.result().relative_residuals[0] <= 1e-8);
        assert!(result.nodal_mode(1).is_err());
    });
}

#[test]
fn hexahedral_buckling_retains_the_accepted_strain_formulation() {
    with_cx(false, |cx| {
        let nodes = [
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 0., 1.],
            [1., 1., 1.],
            [0., 1., 1.],
        ];
        let cells = [[0, 1, 2, 3, 4, 5, 6, 7]];
        let fixed: Vec<_> = (0..24).filter(|&dof| dof != 1).collect();
        let material = TetElasticMaterial::try_new(6.0, 1000.0, 0.0, [31; 32]).unwrap();
        let thermal = strain(&material, [0.01, 0., 0., 0., 0., 0.]);
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let problem = Hex8ElasticProblem {
                nodes_m: &nodes,
                hexahedra: &cells,
                materials: TetMaterialField::Uniform(&material),
                formulation,
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                stiffness_points_per_dir: 2,
            };
            let reference = problem
                .solve_thermal_tangent(
                    TetThermalStrainField::Uniform(&thermal),
                    TetStaticSolveConfig::default(),
                    ContentHash([33; 32]),
                    0.02,
                    cx,
                )
                .unwrap();
            let result = reference.proportional_buckling(1, 100, cx).unwrap();
            // Unit-cube Q1 integrals: integral (dN/dx)^2 = 1/9,
            // integral dN/dy = -1/4. sigma_xx = -E a throughout.
            let expected = match formulation {
                Formulation::Standard => 200.0,
                Formulation::BBar => 89.0 / (48.0 * 0.01),
            };
            close(result.result().loads[0], expected);
            close(reference.geometric_stiffness().get(0, 0), -10.0 / 9.0);
            assert_eq!(result.reference_assembly_identity(), reference.identity());
            assert_eq!(result.nodal_mode(0).unwrap().len(), 8);
        }
    });
}

#[test]
fn tension_zero_reference_and_cancellation_do_not_fabricate_a_mode() {
    with_cx(false, |cx| {
        for a in [0.0, -0.01] {
            let result = tets(a, cx).proportional_buckling(1, 100, cx).unwrap();
            assert!(result.result().loads.is_empty());
            assert!(result.nodal_mode(0).is_err());
        }
        let reference = tets(0.01, cx);
        with_cx(true, |cancelled| {
            assert!(matches!(
                reference.proportional_buckling(1, 100, cancelled),
                Err(SolidError::Cancelled)
            ));
        });
    });
}
