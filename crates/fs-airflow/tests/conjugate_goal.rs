//! G1/G3: actual P1 solid/air solves with an independent effective-conductance
//! adjoint oracle, downstream feedback and nonlinear remainder scaling.
//! G4: invalid binding, nonconverged reference, budget and cancellation refusal.

use fs_airflow::conjugate::goal::{
    CoupledGoalConfig, CoupledGoalError, InterfaceSolveConfig, compare_discrete_goal,
};
use fs_airflow::conjugate::{AirPath, AirSegment, ConjugateConfig, SolidRegionState};
use fs_airflow::graph::thermal::solve_conjugate_branches_iqn;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::assemble::{assemble_operator, full_residual};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductionSolution, ConductivityModel, ConductivityTable,
    InitialGuess, ScalarField, SolveConfig, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_couple::iqn_ils::IqnIlsConfig;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 907,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn config() -> CoupledGoalConfig {
    let mut solid = SolveConfig::default().linear;
    solid.tolerance = 1e-8;
    CoupledGoalConfig {
        solid,
        primal: ConjugateConfig {
            temperature_tolerance_k: 1e-10,
            max_iterations: 40,
            ..ConjugateConfig::default()
        },
        interface: InterfaceSolveConfig {
            max_iterations: 24,
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            relaxation: 1.0,
        },
        acceleration: IqnIlsConfig::default(),
    }
}
fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([3, 2, 2], [0.1, 1.0, 1.0]);
    ConductionMesh::new(complex, positions).unwrap()
}
fn right_path() -> Vec<AirPath> {
    vec![
        AirPath::new(
            290.0,
            0.01,
            1000.0,
            vec![AirSegment::new("right", 1.0, 80.0).unwrap()],
        )
        .unwrap(),
    ]
}
fn boundary(mesh: &ConductionMesh, paths: &[AirPath], references: &[f64]) -> ThermalBoundary {
    let mut builder = ThermalBoundaryBuilder::new(mesh);
    let mut index = 0;
    for path in paths {
        for segment in path.segments() {
            let x = if segment.region() == "left" { 0.0 } else { 0.1 };
            builder = builder
                .region(
                    segment.region(),
                    |face| on_box_face(face.centroid[0], x),
                    ThermalBc::robin(segment.htc_w_per_m2_k(), references[index]).unwrap(),
                )
                .unwrap();
            index += 1;
        }
    }
    builder.adiabatic_remainder().finish().unwrap()
}
fn solved(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    material: &ConductivityModel,
    source: &ScalarField,
    paths: &[AirPath],
) -> (ThermalBoundary, ConductionSolution) {
    let mut retained = None;
    solve_conjugate_branches_iqn(
        cx,
        paths,
        &config().primal,
        IqnIlsConfig::default(),
        |cx, references| {
            let boundary = boundary(mesh, paths, references);
            let mut config = SolveConfig::default();
            config.initial = InitialGuess::Uniform(300.0);
            config.linear.tolerance = 1e-12;
            config.stop.residual_rtol = 1e-12;
            config.stop.residual_atol = 0.0;
            config.stop.step_atol = 0.0;
            let solution = fs_conduction::solve(
                cx,
                ConductionProblem {
                    mesh,
                    boundary: &boundary,
                    material,
                    element_materials: None,
                    source,
                },
                config,
            )
            .unwrap();
            let states = paths
                .iter()
                .flat_map(|path| path.regions())
                .map(|name| {
                    SolidRegionState::from_robin_flux(
                        solution
                            .report
                            .robin_fluxes
                            .iter()
                            .find(|flux| flux.region == name)
                            .unwrap(),
                    )
                })
                .collect();
            retained = Some((boundary, solution));
            Ok(states)
        },
    )
    .unwrap();
    retained.unwrap()
}
fn wall_weights(mesh: &ConductionMesh, x: f64) -> Vec<f64> {
    let mut weights = vec![0.0; mesh.vertex_count()];
    for face in mesh.boundary() {
        if on_box_face(face.centroid[0], x) {
            for vertex in face.vertices {
                weights[vertex as usize] += face.area / 3.0;
            }
        }
    }
    weights
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn coupled_residual_contributions_match_effective_conductance_not_frozen_air() {
    with_gate(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let paths = right_path();
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(100.0);
        let (boundary, reference) = solved(cx, &mesh, &material, &source, &paths);
        let weights = wall_weights(&mesh, 0.1);
        let approximate: Vec<f64> = reference
            .temperature
            .iter()
            .zip(mesh.positions())
            .map(|(t, p)| t + 0.7 + 0.3 * p[0] / 0.1 + 0.2 * p[1] * p[2])
            .collect();
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let report = compare_discrete_goal(
            cx,
            problem,
            None,
            &paths,
            config(),
            &reference.temperature,
            &approximate,
            &weights,
        )
        .unwrap();
        let effective = 10.0 * (-(-8.0_f64).exp_m1());
        let actual_wall = dot(&weights, &reference.temperature);
        assert!((actual_wall - (290.0 + 10.0 / effective)).abs() < 1e-7);
        // For this wall-mean goal, the physical adjoint is identically
        // 1/(C epsilon) throughout the insulated, source-driven solid.
        let residual = |field: &[f64]| {
            let air = paths[0].march(&[dot(&weights, field)]).unwrap();
            let bc = boundary_for_path(&mesh, &paths, &air.reference_temperatures_k());
            let system = assemble_operator(cx, &mesh, &bc, &material, &source, field).unwrap();
            full_residual(&system, field)
        };
        let r = residual(&reference.temperature);
        let a = residual(&approximate);
        for ((contribution, r), a) in report.goal.nodal_contributions.iter().zip(r).zip(a) {
            assert!(
                (contribution - (r - a) / effective).abs() < 2e-7,
                "{contribution} vs {}",
                (r - a) / effective
            );
        }
        assert!(report.interface_iterations > 1);
        assert!(report.goal.dual_relative_residual < config().solid.tolerance);
        assert!(report.goal.primal_relative_residual < config().solid.tolerance);
        assert!(report.goal.linearization_remainder.abs() < 1e-7);
        let frozen = fs_conduction::adjoint::compare_discrete_goal(
            cx,
            problem,
            None,
            config().solid,
            &reference.temperature,
            &approximate,
            &weights,
        )
        .unwrap();
        assert!(
            frozen
                .nodal_contributions
                .iter()
                .zip(&report.goal.nodal_contributions)
                .any(|(a, b)| (a - b).abs() > 1e-3),
            "fixture must distinguish the physical adjoint from frozen air"
        );
        for scale in [-3.0, 1e-150, 1e150] {
            let scaled: Vec<_> = weights.iter().map(|w| w * scale).collect();
            let next = compare_discrete_goal(
                cx,
                problem,
                None,
                &paths,
                config(),
                &reference.temperature,
                &approximate,
                &scaled,
            )
            .unwrap();
            assert!(
                (next.goal.signed_residual_change / scale - report.goal.signed_residual_change)
                    .abs()
                    < 1e-7
            );
            assert_eq!(next.interface_iterations, report.interface_iterations);
        }
    });
}
// Avoid shadowing the boundary variable in physical-oracle closures.
fn boundary_for_path(
    mesh: &ConductionMesh,
    paths: &[AirPath],
    references: &[f64],
) -> ThermalBoundary {
    boundary(mesh, paths, references)
}

#[test]
fn nonlinear_coupled_goal_retains_a_quadratic_material_remainder() {
    with_gate(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let paths = right_path();
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(250.0, 1.0), (600.0, 36.0)]).unwrap(),
        );
        let source = ScalarField::Uniform(1000.0);
        let (boundary, reference) = solved(cx, &mesh, &material, &source, &paths);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let weights = vec![1.0 / mesh.vertex_count() as f64; mesh.vertex_count()];
        let mut remainders = Vec::new();
        for scale in [1.0, 0.5, 0.25] {
            let approximate: Vec<_> = reference
                .temperature
                .iter()
                .zip(mesh.positions())
                .map(|(t, p)| t + scale * (2.0 + 20.0 * p[0]))
                .collect();
            let report = compare_discrete_goal(
                cx,
                problem,
                None,
                &paths,
                config(),
                &reference.temperature,
                &approximate,
                &weights,
            )
            .unwrap();
            assert!(report.goal.uses_nonlinear_jacobian);
            assert!(report.goal.linearization_remainder.abs() > 1e-5);
            assert!(
                (report.goal.signed_residual_change + report.goal.linearization_remainder
                    - report.goal.signed_goal_change)
                    .abs()
                    < 1e-12
            );
            remainders.push(report.goal.linearization_remainder);
        }
        for pair in remainders.windows(2) {
            assert!((pair[0] / pair[1] - 4.0).abs() < 2e-4, "{remainders:?}");
        }
    });
}

#[test]
fn serial_heating_and_independent_branches_each_close_their_actual_transpose() {
    with_gate(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(100.0);
        for paths in [
            vec![
                AirPath::new(
                    290.0,
                    0.01,
                    1000.0,
                    vec![
                        AirSegment::new("left", 1.0, 20.0).unwrap(),
                        AirSegment::new("right", 1.0, 60.0).unwrap(),
                    ],
                )
                .unwrap(),
            ],
            vec![
                AirPath::new(
                    330.0,
                    0.01,
                    1000.0,
                    vec![AirSegment::new("left", 1.0, 20.0).unwrap()],
                )
                .unwrap(),
                AirPath::new(
                    290.0,
                    0.02,
                    1000.0,
                    vec![AirSegment::new("right", 1.0, 60.0).unwrap()],
                )
                .unwrap(),
            ],
        ] {
            let (boundary, reference) = solved(cx, &mesh, &material, &source, &paths);
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let mut weights = wall_weights(&mesh, 0.0);
            for (w, right) in weights.iter_mut().zip(wall_weights(&mesh, 0.1)) {
                *w -= 0.7 * right;
            }
            let approximate: Vec<_> = reference
                .temperature
                .iter()
                .zip(mesh.positions())
                .map(|(t, p)| t + 1.0 - 15.0 * p[0])
                .collect();
            let report = compare_discrete_goal(
                cx,
                problem,
                None,
                &paths,
                config(),
                &reference.temperature,
                &approximate,
                &weights,
            )
            .unwrap();
            assert!(report.goal.dual_relative_residual < config().solid.tolerance);
            assert!(report.goal.linearization_remainder.abs() < 1e-7);
            let difference: Vec<_> = reference
                .temperature
                .iter()
                .zip(&approximate)
                .map(|(a, b)| a - b)
                .collect();
            assert!((report.goal.signed_residual_change - dot(&weights, &difference)).abs() < 1e-7);
        }
    });
}

#[test]
fn invalid_binding_reference_budget_and_cancellation_refuse_without_a_result() {
    with_gate(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let paths = right_path();
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(100.0);
        let (boundary, reference) = solved(cx, &mesh, &material, &source, &paths);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let weights = wall_weights(&mesh, 0.1);
        let run = |paths: &[AirPath], config, reference: &[f64]| {
            compare_discrete_goal(
                cx, problem, None, paths, config, reference, reference, &weights,
            )
        };
        let mut short = config();
        short.interface.max_iterations = 1;
        assert!(matches!(
            run(&paths, short, &reference.temperature),
            Err(CoupledGoalError::DidNotConverge { .. })
        ));
        let mut bad_budget = config();
        bad_budget.interface.max_iterations = 0;
        assert!(matches!(
            run(&paths, bad_budget, &reference.temperature),
            Err(CoupledGoalError::InvalidInput(_))
        ));
        let bad_paths = vec![
            AirPath::new(
                290.0,
                0.01,
                1000.0,
                vec![AirSegment::new("right", 1.0, 81.0).unwrap()],
            )
            .unwrap(),
        ];
        assert!(matches!(
            run(&bad_paths, config(), &reference.temperature),
            Err(CoupledGoalError::InvalidInput(_))
        ));
        let wrong_inlet = vec![
            AirPath::new(
                280.0,
                0.01,
                1000.0,
                vec![AirSegment::new("right", 1.0, 80.0).unwrap()],
            )
            .unwrap(),
        ];
        assert!(matches!(
            run(&wrong_inlet, config(), &reference.temperature),
            Err(CoupledGoalError::Air(_))
        ));
        let shifted: Vec<_> = reference.temperature.iter().map(|t| t + 0.1).collect();
        assert!(matches!(
            run(&paths, config(), &shifted),
            Err(CoupledGoalError::Solid(_))
        ));
        let gate = CancelGate::new_clock_free();
        gate.request();
        with_gate(&gate, |cx| {
            assert!(matches!(
                compare_discrete_goal(
                    cx,
                    problem,
                    None,
                    &paths,
                    config(),
                    &reference.temperature,
                    &reference.temperature,
                    &weights
                ),
                Err(CoupledGoalError::Interrupted)
            ))
        });
        let zero = compare_discrete_goal(
            cx,
            problem,
            None,
            &paths,
            config(),
            &reference.temperature,
            &reference.temperature,
            &vec![0.0; weights.len()],
        )
        .unwrap();
        assert_eq!(zero.goal.dual_iterations, 0);
        assert_eq!(zero.goal.signed_residual_change, 0.0);
    });
}
