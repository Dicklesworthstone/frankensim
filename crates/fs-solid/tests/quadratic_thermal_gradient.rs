//! G1 thermal-gradient bending, stored energy and signed vibration stiffness.
//! The analytic fields are independent of the assembly/quadrature arithmetic.
#[allow(dead_code, unused_imports)]
#[path = "hex/support.rs"]
mod hex_support;
#[allow(dead_code)]
#[path = "tet10/support.rs"]
mod tet_support;

use fs_blake3::ContentHash;
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_solid::linear3::hex::{HexConnectivity, HexElasticProblem, hex20_point, hex27_point};
use fs_solid::linear3::tet10::{QuadraticTetProblem, tet10_point};
use fs_solid::linear3::{
    QuadraticPrestressedAssembly, QuadraticThermalResponse, QuadraticThermalStrainField,
    TetAssemblyBudget, TetElasticError, TetElasticMaterial, TetMaterialField, TetStaticSolveConfig,
    TetThermalStrainField, TetThermalStrainState,
};
use tet_support::{block, material, promote, quadratic_form, with_cx};

const FRAME: ContentHash = ContentHash([184; 32]);
const GRADIENT_PER_M: f64 = 0.001;

fn thermal(law: &TetElasticMaterial, z: f64, reference_z: f64) -> TetThermalStrainState {
    // Fixed anisotropic expansion alpha_x = 1e-5/K and dT/dz = 100 K/m,
    // independent of the elastic symmetry. No quadrature rule enters this law.
    TetThermalStrainState::try_new(
        300.0 + 100.0 * (z - reference_z),
        300.0,
        [GRADIENT_PER_M * (z - reference_z), 0.0, 0.0, 0.0, 0.0, 0.0],
        law.material_state_identity,
        ContentHash([185; 32]),
        0.01,
    )
    .unwrap()
}

fn tet_samples(
    problem: &QuadraticTetProblem<'_>,
    law: &TetElasticMaterial,
    reference_z: f64,
) -> Vec<TetThermalStrainState> {
    let rule = duffy_quadrature(problem.quadrature.stiffness_points_per_dir);
    problem
        .elements
        .iter()
        .enumerate()
        .flat_map(|(e, nodes)| {
            rule.iter().enumerate().map(move |(i, &(r, _))| {
                let p = tet10_point(nodes.map(|n| problem.nodes_m[n]), r, 1e-12, e, i).unwrap();
                thermal(law, p.position_m[2], reference_z)
            })
        })
        .collect()
}

fn hex_samples(
    problem: &HexElasticProblem<'_>,
    law: &TetElasticMaterial,
    reference_z: f64,
) -> Vec<TetThermalStrainState> {
    let rule = problem.quadrature.points().unwrap();
    let positions: Vec<_> = match problem.elements {
        HexConnectivity::Hex20(elements) => elements
            .iter()
            .enumerate()
            .flat_map(|(e, nodes)| {
                rule.iter().enumerate().map(move |(i, &(r, _))| {
                    hex20_point(nodes.map(|n| problem.nodes_m[n]), r, 1e-12, e, i)
                        .unwrap()
                        .position_m
                })
            })
            .collect(),
        HexConnectivity::Hex27(elements) => elements
            .iter()
            .enumerate()
            .flat_map(|(e, nodes)| {
                rule.iter().enumerate().map(move |(i, &(r, _))| {
                    hex27_point(nodes.map(|n| problem.nodes_m[n]), r, 1e-12, e, i)
                        .unwrap()
                        .position_m
                })
            })
            .collect(),
    };
    positions
        .iter()
        .map(|p| thermal(law, p[2], reference_z))
        .collect()
}

fn sampled(states: &[TetThermalStrainState]) -> QuadraticThermalStrainField<'_> {
    QuadraticThermalStrainField::AtQuadrature {
        states,
        points_per_element: 27,
    }
}

fn fixed_face(nodes: &[[f64; 3]], axis: usize) -> Vec<usize> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[axis] == 0.0)
        .flat_map(|(node, _)| [3 * node, 3 * node + 1, 3 * node + 2])
        .collect()
}

fn close(actual: f64, exact: f64, tolerance: f64) {
    assert!(
        (actual - exact).abs() <= tolerance,
        "actual={actual:e}, exact={exact:e}, tolerance={tolerance:e}"
    );
}

fn check_bending(nodes: &[[f64; 3]], result: &QuadraticThermalResponse) {
    for (&[x, _, z], displacement) in nodes.iter().zip(result.displacement_m()) {
        let expected = [GRADIENT_PER_M * x * z, 0.0, -0.5 * GRADIENT_PER_M * x * x];
        for (&actual, exact) in displacement.iter().zip(expected) {
            close(actual, exact, 2e-11);
        }
    }
    assert!(result.elastic_energy_j() < 1e-14);
    assert!(
        result
            .elements()
            .iter()
            .all(|e| e.mean_stress_mandel_pa().iter().all(|s| s.abs() < 1e-5))
    );
    assert!(result.true_relative_residual() < 1e-10);
}

#[test]
fn g1_through_thickness_heating_produces_exact_compatible_quadratic_bending() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let fixed = fixed_face(&nodes, 0);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        let samples = tet_samples(&problem, &law, 0.0);
        let response = problem
            .solve_thermal_response(
                sampled(&samples),
                TetStaticSolveConfig::default(),
                FRAME,
                0.01,
                cx,
            )
            .unwrap();
        check_bending(&nodes, &response);
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed = fixed_face(&mesh.nodes_m, 0);
            let problem = hex_support::problem(&mesh, &law, &fixed);
            let samples = hex_samples(&problem, &law, 0.0);
            let response = problem
                .solve_thermal_response(
                    sampled(&samples),
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.01,
                    cx,
                )
                .unwrap();
            check_bending(&mesh.nodes_m, &response);
        }
    });
}

fn check_restrained(result: &QuadraticThermalResponse, law: &TetElasticMaterial) {
    // 1/2 C_xxxx b^2 integral_0^1 (z-1/2)^2 dz = C_xxxx b^2/24.
    let expected = law.stiffness_mandel_pa()[0][0] * GRADIENT_PER_M.powi(2) / 24.0;
    close(result.elastic_energy_j(), expected, 1e-12);
    assert_eq!(result.iterations(), 0);
    let mean_force_density: f64 = result
        .elements()
        .iter()
        .map(|e| e.quadrature_volume_m3() * e.mean_stress_mandel_pa()[0])
        .sum();
    close(mean_force_density, 0.0, 1e-9);
}

#[test]
fn g1_zero_mean_thermal_gradient_retains_restrained_elastic_energy() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
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
        let samples = tet_samples(&problem, &law, 0.5);
        check_restrained(
            &problem
                .solve_thermal_response(
                    sampled(&samples),
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.01,
                    cx,
                )
                .unwrap(),
            &law,
        );
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed: Vec<_> = (0..3 * mesh.nodes_m.len()).collect();
            let problem = hex_support::problem(&mesh, &law, &fixed);
            let samples = hex_samples(&problem, &law, 0.5);
            check_restrained(
                &problem
                    .solve_thermal_response(
                        sampled(&samples),
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.01,
                        cx,
                    )
                    .unwrap(),
                &law,
            );
            let mean = thermal(&law, 0.5, 0.5);
            assert_eq!(
                problem
                    .solve_thermal_response(
                        TetThermalStrainField::Uniform(&mean),
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.01,
                        cx
                    )
                    .unwrap()
                    .elastic_energy_j(),
                0.0
            );
        }
    });
}

fn check_prestress(
    nodes: &[[f64; 3]],
    result: &QuadraticPrestressedAssembly,
    law: &TetElasticMaterial,
) {
    // Trial displacement v=(z^2,0,0) vanishes on the fixed z=0 face.
    // integral grad(v)' sigma grad(v) = -4 C_zzxx b integral z^2(z-1/2) dz.
    let mode: Vec<_> = result
        .pencil()
        .free_dofs
        .iter()
        .map(|&dof| {
            if dof % 3 == 0 {
                nodes[dof / 3][2].powi(2)
            } else {
                0.0
            }
        })
        .collect();
    let expected = -law.stiffness_mandel_pa()[2][0] * GRADIENT_PER_M / 3.0;
    close(
        quadratic_form(result.geometric_stiffness(), &mode),
        expected,
        1e-8,
    );
    assert!(
        result
            .response()
            .displacement_m()
            .iter()
            .flatten()
            .all(|x| *x == 0.0)
    );
}

#[test]
fn g1_vibration_tangent_retains_gradient_stress_even_when_mean_stress_is_zero() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let fixed = fixed_face(&nodes, 2);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        let samples = tet_samples(&problem, &law, 0.5);
        let thermal_load = problem
            .assemble_thermal_load(sampled(&samples), cx)
            .unwrap();
        let loads: Vec<_> = thermal_load
            .full_equivalent_force_n
            .iter()
            .enumerate()
            .map(|(dof, force)| (dof, -force))
            .collect();
        let tangent = problem
            .solve_thermal_tangent_with_nodal_forces(
                sampled(&samples),
                &loads,
                TetStaticSolveConfig::default(),
                FRAME,
                0.01,
                cx,
            )
            .unwrap();
        check_prestress(&nodes, &tangent, &law);
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let fixed = fixed_face(&mesh.nodes_m, 2);
            let problem = hex_support::problem(&mesh, &law, &fixed);
            let samples = hex_samples(&problem, &law, 0.5);
            let thermal_load = problem
                .assemble_thermal_load(sampled(&samples), cx)
                .unwrap();
            let loads: Vec<_> = thermal_load
                .full_equivalent_force_n
                .iter()
                .enumerate()
                .map(|(dof, force)| (dof, -force))
                .collect();
            let tangent = problem
                .solve_thermal_tangent_with_nodal_forces(
                    sampled(&samples),
                    &loads,
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.01,
                    cx,
                )
                .unwrap();
            check_prestress(&mesh.nodes_m, &tangent, &law);
        }
    });
}

#[test]
fn g0_quadrature_count_and_late_source_mismatch_refuse_without_changing_input() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let foreign = TetElasticMaterial::try_new(1.0, 2e6, 0.25, [186; 32]).unwrap();
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let problem = hex_support::problem(&mesh, &law, &[]);
            let mut samples = hex_samples(&problem, &law, 0.0);
            let valid = problem
                .assemble_thermal_load(sampled(&samples), cx)
                .unwrap();
            assert!(matches!(
                problem.assemble_thermal_load(
                    QuadraticThermalStrainField::AtQuadrature {
                        states: &samples,
                        points_per_element: 8
                    },
                    cx
                ),
                Err(TetElasticError::InvalidThermalStrain { .. })
            ));
            assert!(matches!(
                problem.assemble_thermal_load(sampled(&samples[..26]), cx),
                Err(TetElasticError::ThermalStrainCountMismatch {
                    expected: 27,
                    got: 26
                })
            ));
            let saved = samples[26];
            samples[26] = thermal(&foreign, 0.9, 0.0);
            assert!(matches!(
                problem.assemble_thermal_load(sampled(&samples), cx),
                Err(TetElasticError::ThermalElasticStateMismatch { element: 0 })
            ));
            samples[26] = saved;
            assert_eq!(
                valid,
                problem
                    .assemble_thermal_load(sampled(&samples), cx)
                    .unwrap()
            );
        }
    });
}

#[test]
fn g3_constant_quadrature_samples_reproduce_the_original_elementwise_load() {
    with_cx(|cx| {
        let law = material(2e6, 0.25);
        let state = thermal(&law, 0.7, 0.0);
        let (nodes, elements) = promote(&block([1; 3], [1.0; 3]), cx);
        let problem = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
            quadrature: tet_support::Q,
        };
        let pointwise = problem
            .assemble_thermal_load(sampled(&vec![state; 27 * elements.len()]), cx)
            .unwrap();
        let elementwise = problem
            .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
            .unwrap();
        assert_eq!(
            pointwise.full_equivalent_force_n,
            elementwise.full_equivalent_force_n
        );
        assert_eq!(pointwise.reduced_force_n, elementwise.reduced_force_n);
        for family in hex_support::FAMILIES {
            let mesh = hex_support::mesh([1; 3], [1.0; 3], family, cx);
            let problem = hex_support::problem(&mesh, &law, &[]);
            let pointwise = problem
                .assemble_thermal_load(sampled(&[state; 27]), cx)
                .unwrap();
            let elementwise = problem
                .assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
                .unwrap();
            assert_eq!(
                pointwise.full_equivalent_force_n,
                elementwise.full_equivalent_force_n
            );
            assert_eq!(pointwise.reduced_force_n, elementwise.reduced_force_n);
        }
    });
}
