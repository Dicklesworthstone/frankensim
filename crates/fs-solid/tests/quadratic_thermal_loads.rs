//! Mixed thermal eigenstrain and dead traction on the existing quadratic solids.
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
    TetMaterialField, TetStaticSolveConfig, TetThermalLoad, TetThermalStrainField,
    TetThermalStrainState,
};
use tet_support::{block, material, promote, with_cx};

const FRAME: ContentHash = ContentHash([86; 32]);

fn thermal(law: &TetElasticMaterial, xx: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        400.0,
        300.0,
        [xx, 0.0, 0.0, 0.0, 0.0, 0.0],
        law.material_state_identity,
        ContentHash([87; 32]),
        0.01,
    )
    .unwrap()
}

fn clamp_left(nodes: &[[f64; 3]]) -> Vec<usize> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0)
        .flat_map(|(node, _)| [3 * node, 3 * node + 1, 3 * node + 2])
        .collect()
}

// With nu=0, auxiliary epsilon_xx=.001 gives uniform sigma_xx=2000 Pa.
// Restrict its consistent equivalent load to x=1 to apply that dead traction.
// All source element families must independently integrate a 2000 N resultant.
fn end_traction(nodes: &[[f64; 3]], load: TetThermalLoad) -> Vec<(usize, f64)> {
    let forces: Vec<_> = nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 1.0)
        .map(|(node, _)| (3 * node, load.full_equivalent_force_n[3 * node]))
        .collect();
    assert!((forces.iter().map(|(_, force)| force).sum::<f64>() - 2000.0).abs() < 1e-8);
    forces
}

fn check_uniaxial(nodes: &[[f64; 3]], response: &QuadraticThermalResponse) {
    assert!(response.true_relative_residual() <= 1e-10);
    for (node, &[x, _, _]) in nodes.iter().enumerate() {
        let u = response.displacement_m()[node];
        assert!((u[0] - 0.00125 * x).abs() < 1e-11);
        assert!(u[1].abs() < 1e-11 && u[2].abs() < 1e-11);
        if x != 0.0 {
            assert!(
                response.nodal_force_balance_n()[node]
                    .iter()
                    .all(|v| v.abs() < 1e-5)
            );
        }
    }
    for element in response.elements() {
        assert!((element.mean_stress_mandel_pa()[0] - 2000.0).abs() < 1e-5);
        assert!(
            element.mean_stress_mandel_pa()[1..]
                .iter()
                .all(|v| v.abs() < 1e-5)
        );
        assert!((element.mean_elastic_strain_mandel()[0] - 0.001).abs() < 1e-11);
    }
    assert!((response.elastic_energy_j() - 1.0).abs() < 1e-8);
    let support: f64 = nodes
        .iter()
        .zip(response.nodal_force_balance_n())
        .filter(|(p, _)| p[0] == 0.0)
        .map(|(_, force)| force[0])
        .sum();
    assert!((support + 2000.0).abs() < 1e-5);
}

#[test]
fn g1_dead_traction_and_thermal_expansion_share_quadratic_equilibrium() {
    with_cx(|cx| {
        let law = material(2e6, 0.0);
        let actual = thermal(&law, 0.00025);
        let auxiliary = thermal(&law, 0.001);
        let field = TetThermalStrainField::Uniform(&actual);
        macro_rules! check {
            ($problem:expr, $nodes:expr) => {{
                let problem = $problem;
                let nodes: &[[f64; 3]] = $nodes;
                let forces = end_traction(
                    nodes,
                    problem
                        .assemble_thermal_load(TetThermalStrainField::Uniform(&auxiliary), cx)
                        .unwrap(),
                );
                let pure = problem.assemble_thermal_load(field, cx).unwrap();
                let combined = problem
                    .assemble_thermal_load_with_nodal_forces(field, &forces, cx)
                    .unwrap();
                let mut expected = pure.full_equivalent_force_n.clone();
                for &(dof, force) in &forces {
                    expected[dof] += force;
                }
                assert_eq!(combined.full_equivalent_force_n, expected);
                let response = problem
                    .solve_thermal_response_with_nodal_forces(
                        field,
                        &forces,
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.01,
                        cx,
                    )
                    .unwrap();
                assert_eq!(response.thermal_load_identity(), pure.identity);
                assert_eq!(response.load_identity(), combined.identity);
                assert_ne!(response.load_identity(), response.thermal_load_identity());
                check_uniaxial(nodes, &response);
                for &(dof, force) in &forces {
                    let internal = response.nodal_internal_force_n()[dof / 3][dof % 3];
                    let balance = response.nodal_force_balance_n()[dof / 3][dof % 3];
                    assert!((internal - balance - force).abs() < 1e-10);
                }
            }};
        }
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let fixed = clamp_left(&nodes);
        check!(
            QuadraticTetProblem {
                nodes_m: &nodes,
                elements: &elements,
                materials: TetMaterialField::Uniform(&law),
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                quadrature: tet_support::Q,
            },
            &nodes
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed = clamp_left(&mesh.nodes_m);
            check!(hex_support::problem(&mesh, &law, &fixed), &mesh.nodes_m);
        }
    });
}

fn check_fixed(pure: &QuadraticThermalResponse, loaded: &QuadraticThermalResponse) {
    assert_eq!(loaded.iterations(), 0);
    assert_eq!(loaded.displacement_m(), pure.displacement_m());
    assert_eq!(loaded.elements(), pure.elements());
    assert_eq!(
        loaded.nodal_internal_force_n(),
        pure.nodal_internal_force_n()
    );
    assert_eq!(loaded.elastic_energy_j(), pure.elastic_energy_j());
    assert_eq!(loaded.thermal_load_identity(), pure.thermal_load_identity());
    assert_ne!(loaded.load_identity(), pure.load_identity());
    for element in loaded.elements() {
        assert!((element.mean_stress_mandel_pa()[0] + 500.0).abs() < 1e-8);
    }
    let reaction: f64 = loaded
        .nodal_force_balance_n()
        .iter()
        .map(|force| force[0])
        .sum();
    assert!((reaction + 17.0).abs() < 1e-8);
}

#[test]
fn g0_g1_fixed_loads_change_reactions_and_malformed_force_lists_refuse() {
    with_cx(|cx| {
        let law = material(2e6, 0.0);
        let state = thermal(&law, 0.00025);
        let field = TetThermalStrainField::Uniform(&state);
        macro_rules! check {
            ($problem:expr, $nodes:expr) => {{
                let problem = $problem;
                let nodes: &[[f64; 3]] = $nodes;
                let pure = problem
                    .solve_thermal_response(field, TetStaticSolveConfig::default(), FRAME, 0.01, cx)
                    .unwrap();
                let loaded = problem
                    .solve_thermal_response_with_nodal_forces(
                        field,
                        &[(0, 17.0)],
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.01,
                        cx,
                    )
                    .unwrap();
                check_fixed(&pure, &loaded);
                for forces in [
                    vec![(1, 1.0), (0, 2.0)],
                    vec![(0, 1.0), (0, 2.0)],
                    vec![(3 * nodes.len(), 1.0)],
                    vec![(0, f64::NAN)],
                ] {
                    assert!(matches!(
                        problem.assemble_thermal_load_with_nodal_forces(field, &forces, cx),
                        Err(TetElasticError::InvalidNodalForce { .. })
                    ));
                    assert!(matches!(
                        problem.solve_thermal_response_with_nodal_forces(
                            field,
                            &forces,
                            TetStaticSolveConfig::default(),
                            FRAME,
                            0.01,
                            cx,
                        ),
                        Err(TetElasticError::InvalidNodalForce { .. })
                    ));
                }
            }};
        }
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let fixed: Vec<_> = (0..3 * nodes.len()).collect();
        check!(
            QuadraticTetProblem {
                nodes_m: &nodes,
                elements: &elements,
                materials: TetMaterialField::Uniform(&law),
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                quadrature: tet_support::Q,
            },
            &nodes
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed: Vec<_> = (0..3 * mesh.nodes_m.len()).collect();
            check!(hex_support::problem(&mesh, &law, &fixed), &mesh.nodes_m);
        }
    });
}
