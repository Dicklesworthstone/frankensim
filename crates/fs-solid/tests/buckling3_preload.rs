//! Fixed thermal prestress and incremental support loading have distinct factors.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::SolidError;
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::{
    QuadraticPrestressedAssembly, TetAssemblyBudget, TetElasticMaterial, TetLinearElasticProblem,
    TetMaterialField, TetStaticSolveConfig, TetThermalStrainField, TetThermalStrainState,
};
use fs_solid::prestress::PrestressedTetAssembly;

const FRAME: ContentHash = ContentHash([133; 32]);

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
                seed: 131,
                kernel_id: 131,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn law(youngs: f64, poisson: f64, tag: u8) -> TetElasticMaterial {
    TetElasticMaterial::try_new(6.0, youngs, poisson, [tag; 32]).unwrap()
}

fn thermal(material: &TetElasticMaterial, free: [f64; 6]) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        300.0 + 1000.0 * free[0],
        300.0,
        free,
        material.material_state_identity,
        ContentHash([132; 32]),
        0.04,
    )
    .unwrap()
}

fn tet_state(
    material: &TetElasticMaterial,
    expansion: f64,
    compression: f64,
    origin: [f64; 3],
    free_dof: usize,
    cx: &Cx<'_>,
) -> PrestressedTetAssembly {
    let nodes = [
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [0., 0., -1.],
    ]
    .map(|p| core::array::from_fn(|i| p[i] + origin[i]));
    let cells = [[0, 1, 2, 3], [0, 2, 1, 4]];
    let fixed: Vec<_> = (0..15).filter(|&dof| dof != free_dof).collect();
    let state = thermal(material, [expansion, expansion, expansion, 0., 0., 0.]);
    let field = TetThermalStrainField::Uniform(&state);
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &cells,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
    };
    if compression == 0.0 {
        let solution = problem
            .solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
            .unwrap();
        problem
            .assemble_thermal_tangent(&solution, field, FRAME, 0.04, cx)
            .unwrap()
    } else {
        // Node 1's x support gives u_x=-compression*x in both tetrahedra.
        // The sole incremental/free y DOF sees no transverse applied force.
        let prescribed = [(3, -compression)];
        let solution = problem
            .solve_thermal_displacement_with_prescribed_dofs(
                field,
                &[],
                &prescribed,
                0.04,
                TetStaticSolveConfig::default(),
                cx,
            )
            .unwrap();
        problem
            .assemble_thermal_tangent_with_prescribed_dofs(
                &solution,
                field,
                &[],
                &prescribed,
                0.04,
                FRAME,
                0.04,
                cx,
            )
            .unwrap()
    }
}

fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() < 2e-9 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_fixed_thermal_preload_is_not_scaled_with_incremental_support_compression() {
    with_cx(false, |cx| {
        let material = law(1000.0, 0.25, 131);
        for expansion in [-0.01, 0.0, 0.01] {
            let base = tet_state(&material, expansion, 0.0, [0.0; 3], 4, cx);
            for compression in [0.005, 0.01] {
                let reference = tet_state(&material, expansion, compression, [0.0; 3], 4, cx);
                let result = base.incremental_buckling(&reference, 1, 100, cx).unwrap();
                // V=1/3; G=400, 3K=2000, C_xxxx=1200 Pa. The fixed
                // hydrostatic stress is -3K*a, while delta_sigma_xx=-C_xxxx*d.
                // V cancels: lambda=(G-3K*a)/(C_xxxx*d).
                let expected = (400.0 - 2000.0 * expansion) / (1200.0 * compression);
                close(result.result().loads[0], expected);
                assert_eq!(result.base_assembly_identity(), Some(base.identity()));
                assert_eq!(result.reference_assembly_identity(), reference.identity());
                assert_eq!(result.result().dof_map.len(), 15);
                let mode = result.nodal_mode(0).unwrap();
                for dof in 0..15 {
                    if dof != 4 {
                        assert_eq!(mode[dof / 3][dof % 3], 0.0);
                    }
                }
                let increment = reference.geometric_stiffness().get(0, 0)
                    - base.geometric_stiffness().get(0, 0);
                close(-increment * mode[1][1].powi(2), 1.0);
                assert!(result.result().relative_residuals[0] <= 1e-8);
            }
        }
        let base = tet_state(&material, 0.01, 0.0, [0.0; 3], 4, cx);
        assert_eq!(
            base.proportional_buckling(1, 100, cx)
                .unwrap()
                .base_assembly_identity(),
            None
        );
    });
}

fn hex_state(
    expansion: f64,
    compression: f64,
    formulation: Formulation,
    cx: &Cx<'_>,
) -> QuadraticPrestressedAssembly {
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
    let material = law(1000.0, 0.0, 131);
    let state = thermal(&material, [expansion, 0., 0., 0., 0., 0.]);
    let problem = Hex8ElasticProblem {
        nodes_m: &nodes,
        hexahedra: &cells,
        materials: TetMaterialField::Uniform(&material),
        formulation,
        fixed_dofs: &fixed,
        budget: TetAssemblyBudget::standard(),
        stiffness_points_per_dir: 2,
    };
    let prescribed: Vec<_> = nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 1.0)
        .map(|(node, _)| (3 * node, -compression))
        .collect();
    problem
        .solve_thermal_tangent_with_prescribed_dofs(
            TetThermalStrainField::Uniform(&state),
            &[],
            &prescribed,
            0.04,
            TetStaticSolveConfig::default(),
            FRAME,
            0.04,
            cx,
        )
        .unwrap()
}

#[test]
fn g1_hex8_preload_retains_standard_or_bbar_elastic_work() {
    with_cx(false, |cx| {
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let base = hex_state(0.01, 0.0, formulation, cx);
            let reference = hex_state(0.01, 0.005, formulation, cx);
            let result = base.incremental_buckling(&reference, 1, 100, cx).unwrap();
            // Unit-cube integral (dN/dx)^2 = 1/9. Ke*9/E is 2 for compatible
            // strain and 89/48 for the volume-projected B-bar column.
            let stiffness = match formulation {
                Formulation::Standard => 2.0,
                Formulation::BBar => 89.0 / 48.0,
            };
            close(result.result().loads[0], (stiffness - 0.01) / 0.005);
            assert_eq!(result.base_assembly_identity(), Some(base.identity()));
            assert_eq!(result.reference_assembly_identity(), reference.identity());
            assert_eq!(result.nodal_mode(0).unwrap().len(), 8);
        }
    });
}

#[test]
fn g0_pairing_requires_original_geometry_material_and_incremental_supports() {
    with_cx(false, |cx| {
        let material = law(1000.0, 0.25, 131);
        let base = tet_state(&material, 0.01, 0.0, [0.0; 3], 4, cx);
        let moved = tet_state(&material, 0.01, 0.005, [2.0, 0.0, 0.0], 4, cx);
        // Translation preserves the complete reduced elastic matrix. Matrix
        // equality alone therefore cannot admit original-geometry compatibility.
        assert_eq!(base.elastic_stiffness(), moved.elastic_stiffness());
        assert!(matches!(
            base.incremental_buckling(&moved, 1, 100, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        for other in [law(1000.0, 0.25, 132), law(1100.0, 0.25, 131)] {
            let changed = tet_state(&other, 0.01, 0.005, [0.0; 3], 4, cx);
            assert!(matches!(
                base.incremental_buckling(&changed, 1, 100, cx),
                Err(SolidError::InvalidInput { .. })
            ));
        }
        let supports = tet_state(&material, 0.01, 0.0, [0.0; 3], 3, cx);
        assert!(matches!(
            base.incremental_buckling(&supports, 1, 100, cx),
            Err(SolidError::InvalidInput { .. })
        ));
        let compatible = hex_state(0.01, 0.0, Formulation::Standard, cx);
        let projected = hex_state(0.01, 0.005, Formulation::BBar, cx);
        assert!(matches!(
            compatible.incremental_buckling(&projected, 1, 100, cx),
            Err(SolidError::InvalidInput { .. })
        ));
    });
}

#[test]
fn g0_zero_or_tensile_increment_does_not_hide_an_unstable_base() {
    with_cx(false, |cx| {
        let material = law(1000.0, 0.25, 131);
        let base = tet_state(&material, 0.01, 0.0, [0.0; 3], 4, cx);
        assert!(
            base.incremental_buckling(&base, 1, 100, cx)
                .unwrap()
                .result()
                .loads
                .is_empty()
        );
        let tension = tet_state(&material, 0.01, -0.005, [0.0; 3], 4, cx);
        assert!(
            base.incremental_buckling(&tension, 1, 100, cx)
                .unwrap()
                .result()
                .loads
                .is_empty()
        );
        let stiff_bulk = law(1000.0, 0.4995, 131);
        let unstable = tet_state(&stiff_bulk, 0.001, 0.0, [0.0; 3], 4, cx);
        assert!(unstable.pencil().stiffness.get(0, 0) < 0.0);
        assert!(matches!(
            unstable.incremental_buckling(&unstable, 1, 100, cx),
            Err(SolidError::SolveFailed { .. })
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                base.incremental_buckling(&tension, 1, 100, cancelled),
                Err(SolidError::Cancelled)
            ));
        });
    });
}
