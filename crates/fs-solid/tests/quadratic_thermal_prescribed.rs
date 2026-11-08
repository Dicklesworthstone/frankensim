//! Thermal/support superposition and prestress use one physical equilibrium.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex::{
    HEX20_REFERENCE_NODES, HEX27_REFERENCE_NODES, HexConnectivity, HexElasticProblem, HexFamily,
    HexQuadrature,
};
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::{
    QuadraticPrestressedAssembly, QuadraticThermalStrainField, TetAssemblyBudget,
    TetElasticAssembly, TetElasticError, TetElasticMaterial, TetMaterialField,
    TetStaticSolveConfig, TetThermalStrainField, TetThermalStrainState,
};

const FRAME: ContentHash = ContentHash([121; 32]);
const YOUNG: f64 = 2e6;
const FREE: f64 = 0.001;
const FORCE: f64 = 37.0;

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 121,
                kernel_id: 121,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ));
    });
}

fn config() -> TetStaticSolveConfig {
    TetStaticSolveConfig {
        relative_residual_tolerance: 1e-12,
        maximum_iterations: 2000,
    }
}

fn material(youngs: f64, tag: u8) -> TetElasticMaterial {
    TetElasticMaterial::try_new(1.0, youngs, 0.0, [tag; 32]).unwrap()
}

fn thermal(law: &TetElasticMaterial, axial: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        350.0,
        300.0,
        [axial, 0.0, 0.0, 0.0, 0.0, 0.0],
        law.material_state_identity,
        ContentHash([122; 32]),
        0.02,
    )
    .unwrap()
}

fn supports(nodes: &[[f64; 3]]) -> Vec<usize> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0 || p[0] == 1.0)
        .flat_map(|(node, _)| [3 * node, 3 * node + 1, 3 * node + 2])
        .collect()
}

fn prescribed(nodes: &[[f64; 3]], axial: f64) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 1.0)
        .map(|(node, _)| (3 * node, axial))
        .collect()
}

fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite());
    assert!(
        (actual - expected).abs() <= 2e-8 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn geometric_work(nodes: &[[f64; 3]], result: &QuadraticPrestressedAssembly) -> f64 {
    // Quadratic families represent v_y=x(1-x); the two HEX8 layers represent
    // its piecewise-linear interpolant. Both vanish on the prescribed ends.
    let probe: Vec<_> = result
        .pencil()
        .free_dofs
        .iter()
        .map(|&dof| {
            let x = nodes[dof / 3][0];
            if dof % 3 == 1 { x * (1.0 - x) } else { 0.0 }
        })
        .collect();
    let mut force = vec![0.0; probe.len()];
    result.geometric_stiffness().spmv(&probe, &mut force);
    probe.iter().zip(force).map(|(u, f)| u * f).sum()
}

fn check_family(
    nodes: &[[f64; 3]],
    original: &TetElasticAssembly,
    gradient_energy: f64,
    solve: impl Fn(&[(usize, f64)], &[(usize, f64)]) -> QuadraticPrestressedAssembly,
) {
    let force_dof = 3 * nodes.iter().position(|p| p[0] == 1.0).unwrap();
    for axial in [-FREE, FREE, 3.0 * FREE] {
        let supports = prescribed(nodes, axial);
        let result = solve(&[(force_dof, FORCE)], &supports);
        let response = result.response();
        let sigma = YOUNG * (axial - FREE);
        for (x, u) in nodes.iter().zip(response.displacement_m()) {
            close(u[0], axial * x[0]);
            close(u[1], 0.0);
            close(u[2], 0.0);
        }
        for element in response.elements() {
            for point in element.stress_points() {
                close(point.stress_mandel_pa()[0], sigma);
                for &s in &point.stress_mandel_pa()[1..] {
                    close(s, 0.0);
                }
            }
        }
        close(
            response.elastic_energy_j(),
            0.5 * YOUNG * (axial - FREE).powi(2),
        );
        let face_force = |x| {
            nodes
                .iter()
                .zip(response.nodal_force_balance_n())
                .filter(|(p, _)| p[0] == x)
                .map(|(_, force)| force[0])
                .sum::<f64>()
        };
        close(face_force(0.0), -sigma);
        close(face_force(1.0), sigma - FORCE);
        for &dof in &original.free_dofs {
            close(response.nodal_force_balance_n()[dof / 3][dof % 3], 0.0);
        }
        close(geometric_work(nodes, &result), sigma * gradient_energy);
        assert_eq!(result.pencil().free_dofs, original.free_dofs);
        assert_eq!(result.pencil().mass.to_dense(), original.mass.to_dense());
        assert_eq!(
            result.elastic_stiffness().to_dense(),
            original.stiffness.to_dense()
        );
        assert!(response.true_relative_residual() <= config().relative_residual_tolerance);
        assert_eq!(
            result.identity(),
            solve(&[(force_dof, FORCE)], &supports).identity()
        );
    }
}

#[test]
fn g1_prescribed_thermal_superposition_drives_all_element_prestress_pencils() {
    with_cx(false, |cx| {
        let law = material(YOUNG, 121);
        let state = thermal(&law, FREE);
        let field = TetThermalStrainField::Uniform(&state);
        let corners: [[f64; 3]; 8] =
            core::array::from_fn(|i| HEX27_REFERENCE_NODES[i].map(|v| (f64::from(v) + 1.0) * 0.5));
        let tets = [
            [0, 1, 2, 6],
            [0, 2, 3, 6],
            [0, 3, 7, 6],
            [0, 7, 4, 6],
            [0, 4, 5, 6],
            [0, 5, 1, 6],
        ];
        let (nodes, elements) = promote_p1_to_p2(&corners, &tets, None, cx).unwrap();
        let fixed = supports(&nodes);
        let tet = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: TetQuadrature {
                stiffness_points_per_dir: 3,
                mass_points_per_dir: 4,
            },
        };
        check_family(&nodes, &tet.assemble(cx).unwrap(), 1.0 / 3.0, |loads, g| {
            tet.solve_thermal_tangent_with_prescribed_dofs(
                field,
                loads,
                g,
                0.02,
                config(),
                FRAME,
                0.02,
                cx,
            )
            .unwrap()
        });

        let nodes20 = HEX20_REFERENCE_NODES.map(|p| p.map(|v| (f64::from(v) + 1.0) * 0.5));
        let nodes27 = HEX27_REFERENCE_NODES.map(|p| p.map(|v| (f64::from(v) + 1.0) * 0.5));
        let elements20 = [core::array::from_fn(|i| i)];
        let elements27 = [core::array::from_fn(|i| i)];
        for (family, nodes, elements) in [
            (
                HexFamily::Serendipity20,
                nodes20.as_slice(),
                HexConnectivity::Hex20(&elements20),
            ),
            (
                HexFamily::Lagrange27,
                nodes27.as_slice(),
                HexConnectivity::Hex27(&elements27),
            ),
        ] {
            let fixed = supports(nodes);
            let hex = HexElasticProblem {
                nodes_m: nodes,
                family,
                elements,
                materials: TetMaterialField::Uniform(&law),
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                quadrature: HexQuadrature { points_per_dir: 3 },
            };
            check_family(nodes, &hex.assemble(cx).unwrap(), 1.0 / 3.0, |loads, g| {
                hex.solve_thermal_tangent_with_prescribed_dofs(
                    field,
                    loads,
                    g,
                    0.02,
                    config(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap()
            });
        }

        let nodes: [[f64; 3]; 12] =
            core::array::from_fn(|i| [(i / 4) as f64 * 0.5, (i % 4 / 2) as f64, (i % 2) as f64]);
        let cells = [[0, 4, 6, 2, 1, 5, 7, 3], [4, 8, 10, 6, 5, 9, 11, 7]];
        let fixed = supports(&nodes);
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let hex = Hex8ElasticProblem {
                nodes_m: &nodes,
                hexahedra: &cells,
                materials: TetMaterialField::Uniform(&law),
                formulation,
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                stiffness_points_per_dir: 2,
            };
            check_family(&nodes, &hex.assemble(cx).unwrap(), 0.25, |loads, g| {
                hex.solve_thermal_tangent_with_prescribed_dofs(
                    field,
                    loads,
                    g,
                    0.02,
                    config(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap()
            });
        }
    });
}

#[test]
fn g1_prescribed_lift_retains_local_thermal_stiffness_and_reference_mass() {
    with_cx(false, |cx| {
        let reference = material(YOUNG, 121);
        let nodes = HEX27_REFERENCE_NODES.map(|p| p.map(|v| (f64::from(v) + 1.0) * 0.5));
        let cells = [core::array::from_fn(|i| i)];
        let fixed = supports(&nodes);
        let quadrature = HexQuadrature { points_per_dir: 3 };
        let rule = quadrature.points().unwrap();
        let axial = 2.0 * FREE;
        let sigma = YOUNG * FREE;
        let local: Vec<_> = rule
            .iter()
            .enumerate()
            .map(|(i, (r, _))| {
                let x = (r[0] + 1.0) * 0.5;
                material(YOUNG * (1.0 + 0.5 * x), (i + 1) as u8)
            })
            .collect();
        let states: Vec<_> = local
            .iter()
            .map(|law| thermal(law, axial - sigma / law.stiffness_mandel_pa()[0][0]))
            .collect();
        let field = QuadraticThermalStrainField::ThermomechanicalAtQuadrature {
            states: &states,
            materials: &local,
            points_per_element: rule.len(),
        };
        let problem = HexElasticProblem {
            nodes_m: &nodes,
            family: HexFamily::Lagrange27,
            elements: HexConnectivity::Hex27(&cells),
            materials: TetMaterialField::Uniform(&reference),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature,
        };
        let result = problem
            .solve_thermal_tangent_with_prescribed_dofs(
                field,
                &[],
                &prescribed(&nodes, axial),
                0.02,
                config(),
                FRAME,
                0.02,
                cx,
            )
            .unwrap();
        for (p, u) in nodes.iter().zip(result.response().displacement_m()) {
            close(u[0], axial * p[0]);
            close(u[1], 0.0);
            close(u[2], 0.0);
        }
        for point in result.response().elements()[0].stress_points() {
            close(point.stress_mandel_pa()[0], sigma);
            for &s in &point.stress_mandel_pa()[1..] {
                close(s, 0.0);
            }
        }
        // Exact original quadrature storage for sigma constant and E(x) varying.
        let energy: f64 = local
            .iter()
            .zip(&rule)
            .map(|(law, (_, weight))| {
                0.5 * sigma * sigma / law.stiffness_mandel_pa()[0][0] * weight / 8.0
            })
            .sum();
        close(result.response().elastic_energy_j(), energy);
        close(geometric_work(&nodes, &result), sigma / 3.0);
        assert_eq!(
            result.pencil().mass.to_dense(),
            problem.assemble(cx).unwrap().mass.to_dense()
        );
    });
}

#[test]
fn g1_warped_anisotropic_bricks_retain_stress_free_imposed_support_motion() {
    with_cx(false, |cx| {
        let diagonal = [11.0, 14.0, 9.0, 6.0, 5.0, 8.0];
        let coupling = [0.8, -0.4, 0.3, 0.7, -0.6, 0.2];
        let stiffness = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                1e6 * ((if i == j { diagonal[i] } else { 0.0 }) + coupling[i] * coupling[j])
            })
        });
        let law = TetElasticMaterial::try_new_mandel(1.0, stiffness, [123; 32]).unwrap();
        let gradient = [
            [0.001, 0.0002, -0.0001],
            [0.0002, -0.0004, -0.00015],
            [-0.0001, -0.00015, 0.0007],
        ];
        let s = core::f64::consts::SQRT_2;
        let free = [
            0.001,
            -0.0004,
            0.0007,
            s * 0.0002,
            -s * 0.00015,
            -s * 0.0001,
        ];
        let state = TetThermalStrainState::try_new(
            350.0,
            300.0,
            free,
            law.material_state_identity,
            ContentHash([124; 32]),
            0.02,
        )
        .unwrap();
        let nodes: [[f64; 3]; 8] = core::array::from_fn(|i| {
            let [x, y, z] = HEX27_REFERENCE_NODES[i].map(|v| (f64::from(v) + 1.0) * 0.5);
            [x + 0.12 * y * z, y + 0.08 * x * z, z + 0.06 * x * y]
        });
        let expected: [[f64; 3]; 8] = nodes.map(|p| {
            core::array::from_fn(|i| {
                [0.003, -0.002, 0.001][i] + (0..3).map(|j| gradient[i][j] * p[j]).sum::<f64>()
            })
        });
        let cells = [[0, 1, 2, 3, 4, 5, 6, 7]];
        let fixed = [0, 1, 2, 4, 5, 11];
        let prescribed = fixed.map(|dof| (dof, expected[dof / 3][dof % 3]));
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let problem = Hex8ElasticProblem {
                nodes_m: &nodes,
                hexahedra: &cells,
                materials: TetMaterialField::Uniform(&law),
                formulation,
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                stiffness_points_per_dir: 2,
            };
            let result = problem
                .solve_thermal_response_with_prescribed_dofs(
                    TetThermalStrainField::Uniform(&state),
                    &[],
                    &prescribed,
                    0.02,
                    config(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap();
            for (actual, exact) in result
                .displacement_m()
                .iter()
                .flatten()
                .zip(expected.iter().flatten())
            {
                assert!((actual - exact).abs() < 2e-11);
            }
            for point in result.elements()[0].stress_points() {
                assert!(point.stress_mandel_pa().iter().all(|s| s.abs() < 2e-5));
            }
            assert!(result.elastic_energy_j() < 1e-12);
        }
    });
}

#[test]
fn g0_prescriptions_refuse_wrong_supports_and_preserve_fully_fixed_response() {
    with_cx(false, |cx| {
        let law = material(YOUNG, 121);
        let state = thermal(&law, FREE);
        let field = TetThermalStrainField::Uniform(&state);
        let nodes = HEX27_REFERENCE_NODES.map(|p| p.map(|v| (f64::from(v) + 1.0) * 0.5));
        let cells = [core::array::from_fn(|i| i)];
        let fixed = supports(&nodes);
        let problem = HexElasticProblem {
            nodes_m: &nodes,
            family: HexFamily::Lagrange27,
            elements: HexConnectivity::Hex27(&cells),
            materials: TetMaterialField::Uniform(&law),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
            quadrature: HexQuadrature { points_per_dir: 3 },
        };
        let free_dof = (0..3 * nodes.len())
            .find(|dof| !fixed.contains(dof))
            .unwrap();
        for invalid in [
            vec![(free_dof, 0.001)],
            vec![(fixed[1], 0.001), (fixed[0], 0.001)],
            vec![(fixed[0], 0.001), (fixed[0], 0.001)],
            vec![(fixed[0], f64::NAN)],
            vec![(fixed[0], 0.03)],
        ] {
            assert!(matches!(
                problem.solve_thermal_response_with_prescribed_dofs(
                    field,
                    &[],
                    &invalid,
                    0.02,
                    config(),
                    FRAME,
                    0.02,
                    cx
                ),
                Err(TetElasticError::InvalidStaticSolveConfig { .. })
            ));
        }
        let empty = problem
            .solve_thermal_response_with_prescribed_dofs(
                field,
                &[],
                &[],
                0.02,
                config(),
                FRAME,
                0.02,
                cx,
            )
            .unwrap();
        let old = problem
            .solve_thermal_response(field, config(), FRAME, 0.02, cx)
            .unwrap();
        assert_eq!(empty.displacement_m(), old.displacement_m());
        assert_eq!(empty.nodal_force_balance_n(), old.nodal_force_balance_n());
        assert_eq!(empty.thermal_load_identity(), old.thermal_load_identity());

        let fixed: Vec<_> = (0..3 * nodes.len()).collect();
        let full = HexElasticProblem {
            fixed_dofs: &fixed,
            ..problem
        };
        let prescribed: Vec<_> = nodes
            .iter()
            .enumerate()
            .map(|(i, p)| (3 * i, 2.0 * FREE * p[0]))
            .collect();
        let result = full
            .solve_thermal_response_with_prescribed_dofs(
                field,
                &[],
                &prescribed,
                0.02,
                config(),
                FRAME,
                0.02,
                cx,
            )
            .unwrap();
        assert_eq!(result.iterations(), 0);
        close(
            result.elements()[0].mean_stress_mandel_pa()[0],
            YOUNG * FREE,
        );
        close(result.elastic_energy_j(), 0.5 * YOUNG * FREE * FREE);
        assert!(matches!(
            full.solve_thermal_tangent_with_prescribed_dofs(
                field,
                &[],
                &prescribed,
                0.02,
                config(),
                FRAME,
                0.02,
                cx
            ),
            Err(TetElasticError::NoFreeDofs)
        ));
    });
}
