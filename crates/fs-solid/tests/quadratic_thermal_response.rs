//! Physical thermal equilibrium and stress/energy/reaction observations on
//! quadratic elements, checked against compatible and fully restrained states.
#[allow(dead_code, unused_imports)]
#[path = "hex/support.rs"]
mod hex_support;
#[allow(dead_code)]
#[path = "tet10/support.rs"]
mod tet_support;

use fs_blake3::ContentHash;
use fs_solid::linear3::tet10::QuadraticTetProblem;
use fs_solid::linear3::{
    QuadraticThermalResponse, TetAssemblyBudget, TetElasticError, TetElasticMaterial,
    TetMaterialField, TetStaticSolveConfig, TetThermalStrainField, TetThermalStrainState,
};
use tet_support::{block, material, promote, with_cx};

const FRAME: ContentHash = ContentHash([84; 32]);

fn thermal(law: &TetElasticMaterial, strain: [f64; 6]) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        400.0,
        300.0,
        strain,
        law.material_state_identity,
        ContentHash([85; 32]),
        0.01,
    )
    .unwrap()
}

fn supports(nodes: &[[f64; 3]]) -> Vec<usize> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, p)| {
            (0..3)
                .filter(move |&axis| p[axis] == 0.0)
                .map(move |axis| 3 * node + axis)
        })
        .collect()
}

fn check_expansion(
    nodes: &[[f64; 3]],
    response: &QuadraticThermalResponse,
    expected: impl Fn([f64; 3]) -> [f64; 3],
) {
    assert!(response.true_relative_residual() <= 1e-10);
    assert_eq!(response.mesh_frame(), FRAME);
    for (&p, u) in nodes.iter().zip(response.displacement_m()) {
        for (actual, exact) in u.iter().zip(expected(p)) {
            assert!((actual - exact).abs() < 1e-11, "{actual} != {exact}");
        }
    }
    assert!(
        response
            .elements()
            .iter()
            .all(|e| { e.mean_stress_mandel_pa().iter().all(|v| v.abs() < 1e-5) })
    );
    assert!(response.elastic_energy_j() < 1e-14);
    assert!(
        response
            .nodal_internal_force_n()
            .iter()
            .flatten()
            .all(|v| v.abs() < 1e-5)
    );
}

#[test]
fn g1_quadratic_equilibrium_recovers_free_expansion_and_zero_physical_stress() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let state = thermal(&law, [0.001, 0.001, 0.001, 0.0, 0.0, 0.0]);
        let field = TetThermalStrainField::Uniform(&state);
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let fixed = supports(&nodes);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        let response = problem
            .solve_thermal_response(field, TetStaticSolveConfig::default(), FRAME, 0.01, cx)
            .unwrap();
        check_expansion(&nodes, &response, |p| p.map(|x| 0.001 * x));
        assert_eq!(
            response,
            problem
                .solve_thermal_response(field, TetStaticSolveConfig::default(), FRAME, 0.01, cx,)
                .unwrap(),
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed = supports(&mesh.nodes_m);
            let problem = hex_support::problem(&mesh, &law, &fixed);
            let response = problem
                .solve_thermal_response(field, TetStaticSolveConfig::default(), FRAME, 0.01, cx)
                .unwrap();
            check_expansion(&mesh.nodes_m, &response, |p| p.map(|x| 0.001 * x));
            assert_eq!(
                response.thermal_load_identity(),
                problem.assemble_thermal_load(field, cx).unwrap().identity,
            );
        }
    });
}

fn check_clamped(nodes: &[[f64; 3]], response: &QuadraticThermalResponse) {
    assert_eq!(response.iterations(), 0);
    assert_eq!(response.true_relative_residual(), 0.0);
    assert!(
        response
            .displacement_m()
            .iter()
            .flatten()
            .all(|&u| u == 0.0)
    );
    assert!((response.elastic_energy_j() - 6.0).abs() < 1e-10);
    assert!(
        (response
            .elements()
            .iter()
            .map(|e| e.quadrature_volume_m3())
            .sum::<f64>()
            - 1.0)
            .abs()
            < 1e-12
    );
    for element in response.elements() {
        for i in 0..6 {
            let exact = if i < 3 { -4000.0 } else { 0.0 };
            assert!((element.mean_stress_mandel_pa()[i] - exact).abs() < 1e-8);
        }
        assert!((element.elastic_energy_j() - 6.0 * element.quadrature_volume_m3()).abs() < 1e-10);
    }
    // The x=1 support pushes toward the body: physical reaction is -4000 N.
    let reaction: f64 = nodes
        .iter()
        .zip(response.nodal_internal_force_n())
        .filter(|(p, _)| p[0] == 1.0)
        .map(|(_, f)| f[0])
        .sum();
    assert!((reaction + 4000.0).abs() < 1e-8, "{reaction}");
    for axis in 0..3 {
        let net: f64 = response
            .nodal_internal_force_n()
            .iter()
            .map(|f| f[axis])
            .sum();
        assert!(net.abs() < 1e-8);
    }
}

#[test]
fn g1_fully_restrained_quadratic_heating_recovers_stress_energy_and_reactions() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let state = thermal(&law, [0.001, 0.001, 0.001, 0.0, 0.0, 0.0]);
        let field = TetThermalStrainField::Uniform(&state);
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let fixed: Vec<_> = (0..3 * nodes.len()).collect();
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        check_clamped(
            &nodes,
            &problem
                .solve_thermal_response(field, TetStaticSolveConfig::default(), FRAME, 0.01, cx)
                .unwrap(),
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed: Vec<_> = (0..3 * mesh.nodes_m.len()).collect();
            check_clamped(
                &mesh.nodes_m,
                &hex_support::problem(&mesh, &law, &fixed)
                    .solve_thermal_response(field, TetStaticSolveConfig::default(), FRAME, 0.01, cx)
                    .unwrap(),
            );
        }
    });
}

#[test]
fn g1_spatial_eigenstrain_drives_compatible_nonuniform_expansion() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let strain = |x| {
            thermal(
                &law,
                [if x < 0.5 { 0.001 } else { 0.002 }, 0.0, 0.0, 0.0, 0.0, 0.0],
            )
        };
        let expected = |p: [f64; 3]| {
            [
                0.001 * p[0].min(0.5) + 0.002 * (p[0] - 0.5).max(0.0),
                0.0,
                0.0,
            ]
        };
        let (nodes, elements) = promote(&block([2, 1, 1], [1.0; 3]), cx);
        let states: Vec<_> = elements
            .iter()
            .map(|e| strain(e[..4].iter().map(|&node| nodes[node][0]).sum::<f64>() / 4.0))
            .collect();
        let fixed = supports(&nodes);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        check_expansion(
            &nodes,
            &problem
                .solve_thermal_response(
                    TetThermalStrainField::PerElement(&states),
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.01,
                    cx,
                )
                .unwrap(),
            expected,
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([2, 1, 1], [1.0; 3], family, cx);
            let states: Vec<_> = match &mesh.elements {
                hex_support::HexElements::Hex20(elements) => elements
                    .iter()
                    .map(|e| {
                        strain(
                            e[..8]
                                .iter()
                                .map(|&node| mesh.nodes_m[node][0])
                                .sum::<f64>()
                                / 8.0,
                        )
                    })
                    .collect(),
                hex_support::HexElements::Hex27(elements) => elements
                    .iter()
                    .map(|e| {
                        strain(
                            e[..8]
                                .iter()
                                .map(|&node| mesh.nodes_m[node][0])
                                .sum::<f64>()
                                / 8.0,
                        )
                    })
                    .collect(),
            };
            let fixed = supports(&mesh.nodes_m);
            check_expansion(
                &mesh.nodes_m,
                &hex_support::problem(&mesh, &law, &fixed)
                    .solve_thermal_response(
                        TetThermalStrainField::PerElement(&states),
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.01,
                        cx,
                    )
                    .unwrap(),
                expected,
            );
        }
    });
}

#[test]
fn g0_quadratic_response_refuses_unrestrained_bodies_and_invalid_controls() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let state = thermal(&law, [0.001, 0.001, 0.001, 0.0, 0.0, 0.0]);
        let field = TetThermalStrainField::Uniform(&state);
        let mesh = hex_support::mesh([1; 3], [1.0; 3], hex_support::HexFamily::Lagrange27, cx);
        assert!(matches!(
            hex_support::problem(&mesh, &law, &[]).solve_thermal_response(
                field,
                TetStaticSolveConfig::default(),
                FRAME,
                0.01,
                cx,
            ),
            Err(TetElasticError::UnconstrainedRigidModes { .. })
        ));
        let fixed: Vec<_> = (0..3 * mesh.nodes_m.len()).collect();
        let problem = hex_support::problem(&mesh, &law, &fixed);
        for (frame, ceiling) in [
            (ContentHash([0; 32]), 0.01),
            (FRAME, 0.0),
            (FRAME, f64::NAN),
            (FRAME, 0.0005),
        ] {
            assert!(matches!(
                problem.solve_thermal_response(
                    field,
                    TetStaticSolveConfig::default(),
                    frame,
                    ceiling,
                    cx
                ),
                Err(TetElasticError::InvalidStressObservation { .. })
            ));
        }
        assert!(matches!(
            problem.solve_thermal_response(
                field,
                TetStaticSolveConfig {
                    maximum_iterations: 0,
                    ..TetStaticSolveConfig::default()
                },
                FRAME,
                0.01,
                cx,
            ),
            Err(TetElasticError::InvalidStaticSolveConfig { .. })
        ));
    });
}
