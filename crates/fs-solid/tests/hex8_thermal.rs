//! G0/G1/G3/G4 source-bound HEX8 thermal equilibrium and physical response.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::{ExactLinearSolid, ExactLinearSolidControls};
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticError, TetElasticMaterial, TetMaterialField, TetStaticSolveConfig,
    TetThermalStrainField, TetThermalStrainState,
};

const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const CELLS: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const FIXED: [usize; 6] = [0, 1, 2, 4, 5, 11];
const FREE: [f64; 6] = [0.001, 0.0006, -0.0003, 0.0002, -0.0001, 0.00015];
const FRAME: ContentHash = ContentHash([93; 32]);
const FORMS: [Formulation; 2] = [Formulation::Standard, Formulation::BBar];

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 93,
                kernel_id: 93,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn material(tag: u8) -> TetElasticMaterial {
    let diagonal = [11.0, 14.0, 9.0, 6.0, 5.0, 8.0];
    let a = [0.8, -0.4, 0.3, 0.7, -0.6, 0.2];
    let c = std::array::from_fn(|i| {
        std::array::from_fn(|j| 1e6 * ((if i == j { diagonal[i] } else { 0.0 }) + a[i] * a[j]))
    });
    TetElasticMaterial::try_new_mandel(1500.0, c, [tag; 32]).unwrap()
}

fn thermal(law: &TetElasticMaterial, strain: [f64; 6], temperature: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        temperature,
        300.0,
        strain,
        law.material_state_identity,
        ContentHash([94; 32]),
        0.02,
    )
    .unwrap()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    law: &'a TetElasticMaterial,
    fixed: &'a [usize],
    formulation: Formulation,
) -> Hex8ElasticProblem<'a> {
    Hex8ElasticProblem {
        nodes_m: nodes,
        hexahedra: &CELLS,
        materials: TetMaterialField::Uniform(law),
        formulation,
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        stiffness_points_per_dir: 2,
    }
}

fn config() -> TetStaticSolveConfig {
    TetStaticSolveConfig {
        relative_residual_tolerance: 1e-11,
        maximum_iterations: 2000,
    }
}
fn near(a: f64, b: f64, allowance: f64) {
    assert!(
        (a - b).abs() <= allowance,
        "actual={a:e}, expected={b:e}, allowance={allowance:e}"
    );
}
fn dot(a: [f64; 6], b: [f64; 6]) -> f64 {
    (0..6).map(|i| a[i] * b[i]).sum()
}
fn stress(law: &TetElasticMaterial, e: [f64; 6]) -> [f64; 6] {
    std::array::from_fn(|i| (0..6).map(|j| law.stiffness_mandel_pa()[i][j] * e[j]).sum())
}
fn tensor(s: [f64; 6]) -> [[f64; 3]; 3] {
    let h = std::f64::consts::FRAC_1_SQRT_2;
    [
        [s[0], h * s[3], h * s[5]],
        [h * s[3], s[1], h * s[4]],
        [h * s[5], h * s[4], s[2]],
    ]
}
fn free_displacement([x, y, z]: [f64; 3], e: [f64; 6]) -> [f64; 3] {
    let s = std::f64::consts::SQRT_2;
    [
        e[0] * x + s * (e[3] * y + e[5] * z),
        e[1] * y + s * e[4] * z,
        e[2] * z,
    ]
}

#[test]
fn g1_g3_anisotropic_expansion_is_stress_free_on_the_warped_reference_geometry() {
    with_cx(false, |cx| {
        let law = material(93);
        let state = thermal(&law, FREE, 400.0);
        let nodes = NODES.map(|[r, s, t]| [r + 0.12 * s * t, s + 0.08 * r * t, t + 0.06 * r * s]);
        for form in FORMS {
            let p = problem(&nodes, &law, &FIXED, form);
            let field = TetThermalStrainField::Uniform(&state);
            let result = p
                .solve_thermal_response(field, config(), FRAME, 0.02, cx)
                .unwrap();
            assert_eq!(result.formulation(), form);
            assert_eq!(
                result.thermal_load_identity(),
                p.assemble_thermal_load(field, cx).unwrap().identity
            );
            for (x, u) in nodes.iter().zip(result.displacement_m()) {
                for (actual, expected) in u.iter().zip(free_displacement(*x, FREE)) {
                    near(*actual, expected, 3e-11);
                }
            }
            assert!(result.elastic_energy_j() < 1e-12);
            assert!(
                result.elements()[0]
                    .mean_stress_mandel_pa()
                    .iter()
                    .all(|s| s.abs() < 3e-4)
            );
            assert!(
                result
                    .nodal_force_balance_n()
                    .iter()
                    .flatten()
                    .all(|f| f.abs() < 3e-4)
            );
        }
    });
}

#[test]
fn g1_heterogeneous_thermal_restraint_keeps_each_source_and_its_physical_reaction() {
    with_cx(false, |cx| {
        let index = |x, y, z| 4 * x + 2 * y + z;
        let mut nodes = Vec::new();
        for x in 0..3 {
            for y in 0..2 {
                for z in 0..2 {
                    nodes.push([x as f64, y as f64, z as f64]);
                }
            }
        }
        let cells: Vec<_> = (0..2)
            .map(|x| {
                [
                    index(x, 0, 0),
                    index(x + 1, 0, 0),
                    index(x + 1, 1, 0),
                    index(x, 1, 0),
                    index(x, 0, 1),
                    index(x + 1, 0, 1),
                    index(x + 1, 1, 1),
                    index(x, 1, 1),
                ]
            })
            .collect();
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d % 3 != 0 || nodes[d / 3][0] != 1.0)
            .collect();
        let moduli = [8e6, 20e6];
        let laws: Vec<_> = moduli
            .iter()
            .enumerate()
            .map(|(tag, &e)| {
                let diagonal = [e, 12e6, 12e6, 6e6, 6e6, 6e6];
                let c = std::array::from_fn(|i| {
                    std::array::from_fn(|j| if i == j { diagonal[i] } else { 0.0 })
                });
                TetElasticMaterial::try_new_mandel(1500.0, c, [95 + tag as u8; 32]).unwrap()
            })
            .collect();
        let strains = [0.002, 0.0005];
        let states: Vec<_> = laws
            .iter()
            .zip(strains)
            .enumerate()
            .map(|(i, (law, e))| thermal(law, [e, 0.0, 0.0, 0.0, 0.0, 0.0], [400.0, 325.0][i]))
            .collect();
        // Independent series-bar solution: both ends fixed, equal axial stress.
        let sigma = -(strains[0] + strains[1]) / (1.0 / moduli[0] + 1.0 / moduli[1]);
        let interface = strains[0] + sigma / moduli[0];
        let energy = 0.5 * sigma * sigma * (1.0 / moduli[0] + 1.0 / moduli[1]);
        for form in FORMS {
            let mut p = problem(&nodes, &laws[0], &fixed, form);
            p.hexahedra = &cells;
            p.materials = TetMaterialField::PerElement(&laws);
            let response = p
                .solve_thermal_response(
                    TetThermalStrainField::PerElement(&states),
                    config(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap();
            for (x, u) in nodes.iter().zip(response.displacement_m()) {
                near(u[0], if x[0] == 1.0 { interface } else { 0.0 }, 2e-12);
                near(u[1], 0.0, 1e-15);
                near(u[2], 0.0, 1e-15);
            }
            for (element, &e) in response.elements().iter().zip(&moduli) {
                near(element.mean_stress_mandel_pa()[0], sigma, 2e-6);
                near(element.mean_elastic_strain_mandel()[0], sigma / e, 2e-12);
            }
            near(response.elastic_energy_j(), energy, 2e-8);
            for side in [0.0, 2.0] {
                let reaction: f64 = nodes
                    .iter()
                    .zip(response.nodal_force_balance_n())
                    .filter(|(x, _)| x[0] == side)
                    .map(|(_, f)| f[0])
                    .sum();
                near(reaction, if side == 0.0 { -sigma } else { sigma }, 2e-6);
            }
        }
    });
}

/// Exact polynomial oracle for elastic displacement u_x=a*x*y on [0,1]^3.
/// The three strain coefficients multiply 1, x, y. Integrals of a trilinear
/// shape derivative and its first moments are elementary tensor products;
/// no production geometry, B matrix, quadrature or linear solver is reused.
fn manufactured(
    law: &TetElasticMaterial,
    form: Formulation,
    a: f64,
) -> ([[f64; 3]; 8], f64, [f64; 6]) {
    let h = std::f64::consts::FRAC_1_SQRT_2;
    let coefficients = if form == Formulation::BBar {
        [
            [a / 6.0, a / 6.0, a / 6.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, a * h, 0.0, 0.0],
            [2.0 * a / 3.0, -a / 3.0, -a / 3.0, 0.0, 0.0, 0.0],
        ]
    } else {
        [
            [0.0; 6],
            [0.0, 0.0, 0.0, a * h, 0.0, 0.0],
            [a, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    };
    let s = coefficients.map(|e| stress(law, e));
    let tensors = s.map(tensor);
    let mean = std::array::from_fn(|i| s[0][i] + 0.5 * (s[1][i] + s[2][i]));
    let energy = 0.5
        * (dot(coefficients[0], s[0])
            + dot(coefficients[0], s[1])
            + dot(coefficients[0], s[2])
            + (dot(coefficients[1], s[1]) + dot(coefficients[2], s[2])) / 3.0
            + dot(coefficients[1], s[2]) / 2.0);
    let force = NODES.map(|x| {
        let moments = std::array::from_fn::<_, 3, _>(|axis| {
            let sign = 2.0 * x[axis] - 1.0;
            [
                sign / 4.0,
                if axis == 0 {
                    sign / 8.0
                } else {
                    sign * (1.0 + x[0]) / 12.0
                },
                if axis == 1 {
                    sign / 8.0
                } else {
                    sign * (1.0 + x[1]) / 12.0
                },
            ]
        });
        std::array::from_fn(|c| {
            let mut f: f64 = (0..3)
                .map(|j| {
                    (0..3)
                        .map(|q| moments[j][q] * tensors[q][c][j])
                        .sum::<f64>()
                })
                .sum();
            if form == Formulation::BBar {
                let trace_mean = mean[0] + mean[1] + mean[2];
                let trace_gradient: f64 = (0..3)
                    .map(|q| moments[c][q] * (s[q][0] + s[q][1] + s[q][2]))
                    .sum();
                f += (moments[c][0] * trace_mean - trace_gradient) / 3.0;
            }
            f
        })
    });
    (force, energy, mean)
}

fn free_forces(forces: &[[f64; 3]; 8]) -> Vec<(usize, f64)> {
    (0..24)
        .filter(|d| !FIXED.contains(d))
        .map(|d| (d, forces[d / 3][d % 3]))
        .collect()
}

#[test]
fn g1_projected_equilibrium_stress_force_and_energy_match_independent_work() {
    with_cx(false, |cx| {
        let law = material(93);
        let state = thermal(&law, FREE, 400.0);
        let amplitude = 0.002;
        for form in FORMS {
            let p = problem(&NODES, &law, &FIXED, form);
            let (force, energy, mean_stress) = manufactured(&law, form, amplitude);
            let response = p
                .solve_thermal_response_with_nodal_forces(
                    TetThermalStrainField::Uniform(&state),
                    &free_forces(&force),
                    config(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap();
            for (node, &x) in NODES.iter().enumerate() {
                let mut expected = free_displacement(x, FREE);
                expected[0] += amplitude * x[0] * x[1];
                for c in 0..3 {
                    near(response.displacement_m()[node][c], expected[c], 4e-11);
                    near(
                        response.nodal_internal_force_n()[node][c],
                        force[node][c],
                        5e-6,
                    );
                    near(
                        response.nodal_force_balance_n()[node][c],
                        if FIXED.contains(&(3 * node + c)) {
                            force[node][c]
                        } else {
                            0.0
                        },
                        5e-6,
                    );
                }
            }
            near(response.elastic_energy_j(), energy, 3e-8);
            for i in 0..6 {
                near(
                    response.elements()[0].mean_stress_mandel_pa()[i],
                    mean_stress[i],
                    5e-6,
                );
            }
        }
    });
}

#[test]
fn g1_fully_clamped_heating_keeps_physical_stress_and_fixed_dof_loads() {
    with_cx(false, |cx| {
        let law = material(93);
        let state = thermal(&law, FREE, 400.0);
        let fixed: Vec<_> = (0..24).collect();
        let sigma = stress(&law, FREE).map(|s| -s);
        let physical = tensor(sigma);
        for form in FORMS {
            let p = problem(&NODES, &law, &fixed, form);
            assert!(matches!(p.assemble(cx), Err(TetElasticError::NoFreeDofs)));
            let response = p
                .solve_thermal_response_with_nodal_forces(
                    TetThermalStrainField::Uniform(&state),
                    &[(18, 17.0)],
                    config(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap();
            assert_eq!(response.iterations(), 0);
            assert_eq!(response.true_relative_residual(), 0.0);
            assert!(
                response
                    .displacement_m()
                    .iter()
                    .flatten()
                    .all(|v| *v == 0.0)
            );
            near(response.elastic_energy_j(), -0.5 * dot(FREE, sigma), 2e-9);
            for i in 0..6 {
                near(
                    response.elements()[0].mean_stress_mandel_pa()[i],
                    sigma[i],
                    2e-6,
                );
            }
            for (node, x) in NODES.iter().enumerate() {
                for c in 0..3 {
                    let expected: f64 = (0..3)
                        .map(|j| physical[c][j] * (2.0 * x[j] - 1.0) / 4.0)
                        .sum();
                    near(response.nodal_internal_force_n()[node][c], expected, 2e-6);
                    near(
                        response.nodal_force_balance_n()[node][c],
                        expected - if 3 * node + c == 18 { 17.0 } else { 0.0 },
                        2e-6,
                    );
                }
            }
        }
    });
}

#[test]
fn g1_native_motion_starts_at_the_same_thermal_equilibrium_after_load_release() {
    with_cx(false, |cx| {
        let law = material(93);
        let state = thermal(&law, FREE, 400.0);
        let fixed: Vec<_> = (0..24).filter(|d| *d != 18).collect();
        let p = problem(&NODES, &law, &fixed, Formulation::BBar);
        let field = TetThermalStrainField::Uniform(&state);
        let response = p
            .solve_thermal_response_with_nodal_forces(
                field,
                &[(18, 3.0)],
                config(),
                FRAME,
                0.02,
                cx,
            )
            .unwrap();
        let load = p.assemble_thermal_load(field, cx).unwrap();
        let pencil = p.assemble(cx).unwrap();
        let k = pencil.stiffness.to_dense()[0];
        let m = pencil.mass.to_dense()[0];
        let omega = (k / m).sqrt();
        let phase = 0.2_f64;
        let equilibrium = load.reduced_force_n[0] / k;
        let initial = response.displacement_m()[6][0];
        near(initial, equilibrium + 3.0 / k, 1e-13);
        let mut motion = ExactLinearSolid::prepare(
            &pencil,
            phase / omega,
            &[initial],
            &[0.0],
            ExactLinearSolidControls {
                maximum_free_dofs: 1,
                maximum_preparation_bytes: 64 * 1024 * 1024,
                maximum_exponential_multiply_terms: 100_000_000,
                maximum_relative_operator_error: 1e-12,
                maximum_work_error_j: 1e-8,
                maximum_steps: 1,
            },
            cx,
        )
        .unwrap();
        motion.advance(&load.reduced_force_n, cx).unwrap();
        near(
            motion.displacement_m()[0],
            equilibrium + (initial - equilibrium) * phase.cos(),
            1e-11,
        );
        near(
            motion.velocity_mps()[0],
            -omega * (initial - equilibrium) * phase.sin(),
            1e-9,
        );
    });
}

#[test]
fn g0_g4_source_support_and_raw_strain_admission_are_not_hidden_by_projection() {
    with_cx(false, |cx| {
        let law = material(93);
        let state = thermal(&law, [0.0; 6], 300.0);
        let p = problem(&NODES, &law, &FIXED, Formulation::BBar);
        let other = thermal(&material(97), FREE, 400.0);
        assert!(matches!(
            p.assemble_thermal_load(TetThermalStrainField::Uniform(&other), cx),
            Err(TetElasticError::ThermalElasticStateMismatch { .. })
        ));
        assert!(matches!(
            p.assemble_thermal_load(TetThermalStrainField::PerElement(&[]), cx),
            Err(TetElasticError::ThermalStrainCountMismatch { .. })
        ));
        let free = problem(&NODES, &law, &[], Formulation::BBar);
        assert!(matches!(
            free.solve_thermal_response(
                TetThermalStrainField::Uniform(&state),
                config(),
                FRAME,
                0.02,
                cx
            ),
            Err(TetElasticError::UnconstrainedRigidModes { .. })
        ));
        let (force, _, _) = manufactured(&law, Formulation::BBar, 0.01);
        // The projected xx strain peaks near 0.00693, but the compatible
        // strain reaches 0.00789 at a stiffness point. This ceiling must fail.
        assert!(matches!(
            p.solve_thermal_response_with_nodal_forces(
                TetThermalStrainField::Uniform(&state),
                &free_forces(&force),
                config(),
                FRAME,
                0.0073,
                cx
            ),
            Err(TetElasticError::InvalidStressObservation { .. })
        ));
    });
    with_cx(true, |cx| {
        let law = material(93);
        let state = thermal(&law, FREE, 400.0);
        assert!(matches!(
            problem(&NODES, &law, &FIXED, Formulation::BBar).solve_thermal_response(
                TetThermalStrainField::Uniform(&state),
                config(),
                FRAME,
                0.02,
                cx
            ),
            Err(TetElasticError::Cancelled)
        ));
    });
}
