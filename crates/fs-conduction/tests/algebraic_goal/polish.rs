use super::*;
use fs_conduction::adjoint::{LinearGoalSolveConfig, LinearGoalStop, analyze_linear_maximum};
use fs_conduction::{InitialGuess, SolveConfig, polish_linear_maximum};

fn loose_config() -> SolveConfig {
    let mut cfg = SolveConfig::default();
    cfg.stop.residual_rtol = 1.0;
    cfg.initial = InitialGuess::Uniform(300.0);
    cfg
}

fn control(iterations: usize) -> LinearGoalSolveConfig {
    LinearGoalSolveConfig {
        absolute_tolerance: 1e-8,
        max_primal_iterations: iterations,
        check_every: 1,
        max_defect_corrections: 2,
    }
}

#[test]
fn goal_polish_changes_a_residual_accepted_field_and_refreshes_its_energy_report() {
    let f = Fixture::new(true, false);
    with_cx(|cx| {
        let original = fs_conduction::solve(cx, f.problem(), loose_config()).unwrap();
        assert!(original.report.final_residual <= original.report.residual_threshold);
        let region: Vec<_> = (0..f.mesh.vertex_count()).collect();
        let result = polish_linear_maximum(
            cx,
            f.problem(),
            None,
            linear(),
            &original,
            &region,
            config(),
            control(64),
        )
        .unwrap();
        assert!(result.initial_bound_k.unwrap() > control(64).absolute_tolerance);
        assert!(result.goal_met && result.candidate_accepted);
        assert_eq!(result.stop, LinearGoalStop::GoalTolerance);
        assert!(result.primal_iterations > 0);
        assert!(result.physical_gate_refusal.is_none());
        let solved = &result.solution;
        assert_ne!(solved.temperature, original.temperature);
        assert!(solved.report.final_residual <= original.report.residual_threshold);
        assert!(solved.report.energy.closure_w.abs() < original.report.energy.closure_w.abs());
        assert_eq!(solved.report.iterations, original.report.iterations);
        assert_eq!(
            solved.report.residual_history,
            original.report.residual_history
        );
        assert_eq!(
            solved.report.linear, original.report.linear,
            "no invented Krylov evidence"
        );
        let checked = analyze_linear_maximum(
            cx,
            f.problem(),
            None,
            linear(),
            &solved.temperature,
            &region,
            config(),
        )
        .unwrap();
        assert_eq!(result.analysis, checked);
        let exact = oracle(&f, cx);
        for (&actual, &expected) in solved.temperature.iter().zip(&exact) {
            assert!((actual - expected).abs() < 1e-9);
        }
        assert!(original.temperature.iter().all(|&value| value == 300.0));
    });
}

#[test]
fn corrected_contact_and_robin_reports_follow_the_new_temperature_field() {
    with_cx(|cx| {
        contact_model(1.0, |problem, interfaces| {
            // The side-wall anchor gives this fixture strict unscaled
            // comparison dominance while contact and both Robin faces stay
            // active. The original all-Robin slab has no verified bound.
            let boundary = ThermalBoundaryBuilder::new(problem.mesh)
                .region(
                    "hot",
                    |face| on_box_face(face.centroid[0], 0.0),
                    ThermalBc::robin(20.0, 330.0).unwrap(),
                )
                .unwrap()
                .region(
                    "cold",
                    |face| on_box_face(face.centroid[0], 2.0),
                    ThermalBc::robin(40.0, 300.0).unwrap(),
                )
                .unwrap()
                .region(
                    "anchor",
                    |face| on_box_face(face.centroid[1], 0.0),
                    ThermalBc::dirichlet(315.0).unwrap(),
                )
                .unwrap()
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let problem = ConductionProblem {
                boundary: &boundary,
                ..problem
            };
            let mut cfg = loose_config();
            cfg.initial = InitialGuess::Uniform(315.0);
            let original =
                fs_conduction::solve_with_interfaces(cx, problem, interfaces, cfg).unwrap();
            let region: Vec<_> = (0..problem.mesh.vertex_count()).collect();
            let result = polish_linear_maximum(
                cx,
                problem,
                Some(interfaces),
                linear(),
                &original,
                &region,
                config(),
                control(256),
            )
            .unwrap();
            assert!(result.goal_met && result.candidate_accepted, "{result:?}");
            let solved = result.solution;
            assert_eq!(
                solved.report.interface_fluxes,
                interfaces.fluxes(&solved.temperature).unwrap()
            );
            assert_ne!(
                solved.report.interface_fluxes,
                original.report.interface_fluxes
            );
            assert_ne!(solved.report.robin_fluxes, original.report.robin_fluxes);
            let system = fs_conduction::assemble::assemble_operator_with_interfaces(
                cx,
                problem.mesh,
                problem.boundary,
                problem.material,
                problem.source,
                &original.temperature,
                interfaces,
            )
            .unwrap();
            let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count()).unwrap();
            let (matrix, rhs) = reduce(&system, &dofs);
            let exact = dofs.scatter(&dense_solve(&matrix, &rhs));
            for (&actual, &expected) in solved.temperature.iter().zip(&exact) {
                assert!((actual - expected).abs() < 1e-8);
            }
            // Independent surface integration of the dense oracle; the
            // side-wall heat transfer makes a 1-D series formula inapplicable.
            for (name, plane, htc, reference) in
                [("hot", 0.0, 20.0, 330.0), ("cold", 2.0, 40.0, 300.0)]
            {
                let heat: f64 = problem
                    .mesh
                    .boundary()
                    .iter()
                    .filter(|face| on_box_face(face.centroid[0], plane))
                    .map(|face| {
                        let delta = face
                            .vertices
                            .iter()
                            .map(|&vertex| exact[vertex as usize] - reference)
                            .sum::<f64>()
                            / 3.0;
                        htc * face.area * delta
                    })
                    .sum();
                let row = solved
                    .report
                    .robin_fluxes
                    .iter()
                    .find(|row| row.region == name)
                    .unwrap();
                assert!((row.heat_rate_w - heat).abs() < 1e-6);
            }
            assert!(solved.report.energy.closure_w.abs() < 1e-6);
        })
    });
}

#[test]
fn zero_work_and_missing_inverse_preserve_the_original_physical_field() {
    with_cx(|cx| {
        for (all_fixed, iterations, stability, stop) in [
            (true, 0, 1024, LinearGoalStop::IterationBudget),
            (false, 64, 0, LinearGoalStop::BoundUnavailable),
        ] {
            let f = Fixture::new(all_fixed, false);
            let original = fs_conduction::solve(cx, f.problem(), loose_config()).unwrap();
            let region: Vec<_> = (0..f.mesh.vertex_count()).collect();
            let mut analysis = config();
            analysis.max_stability_iterations = stability;
            let result = polish_linear_maximum(
                cx,
                f.problem(),
                None,
                linear(),
                &original,
                &region,
                analysis,
                control(iterations),
            )
            .unwrap();
            assert_eq!(result.stop, stop);
            assert_eq!(result.primal_iterations, 0);
            assert!(!result.goal_met && !result.candidate_accepted);
            assert_eq!(result.solution, original);
        }
    });
}

#[test]
fn false_baseline_residual_invalid_scale_nonlinearity_and_cancellation_refuse() {
    let mut f = Fixture::new(true, false);
    let original = with_cx(|cx| fs_conduction::solve(cx, f.problem(), loose_config()).unwrap());
    let region: Vec<_> = (0..f.mesh.vertex_count()).collect();
    with_cx(|cx| {
        let mut bad = original.clone();
        bad.report.residual_threshold = original.report.final_residual * 0.5;
        bad.report.final_residual = 0.0; // A forged summary cannot replace reassembly.
        assert!(matches!(
            polish_linear_maximum(
                cx,
                f.problem(),
                None,
                linear(),
                &bad,
                &region,
                config(),
                control(64)
            ),
            Err(ConductionError::NotConverged { .. })
        ));
        bad.report.residual_threshold = f64::INFINITY;
        assert!(
            polish_linear_maximum(
                cx,
                f.problem(),
                None,
                linear(),
                &bad,
                &region,
                config(),
                control(64)
            )
            .is_err()
        );
    });
    with_cx(|cx| {
        let mut invalid_scale = original.clone();
        invalid_scale.report.energy.scale_w = f64::NAN;
        assert!(
            polish_linear_maximum(
                cx,
                f.problem(),
                None,
                linear(),
                &invalid_scale,
                &region,
                config(),
                control(64)
            )
            .is_err()
        );
    });
    with_cancelled_cx(|cx| {
        assert!(matches!(
            polish_linear_maximum(
                cx,
                f.problem(),
                None,
                linear(),
                &original,
                &region,
                config(),
                control(64)
            ),
            Err(ConductionError::Cancelled { .. })
        ));
    });
    f.material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.0), (400.0, 2.0)]).unwrap(),
    );
    with_cx(|cx| {
        assert!(
            polish_linear_maximum(
                cx,
                f.problem(),
                None,
                linear(),
                &original,
                &region,
                config(),
                control(64)
            )
            .is_err()
        );
    });
}
