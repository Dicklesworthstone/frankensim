//! A real consistent-contact matrix defeats scaled comparison dominance but
//! admits a residual-verified approximate inverse under the same work cap.

use super::*;
use fs_conduction::adjoint::{LinearGoalSolveConfig, LinearGoalStop};

#[test]
fn contact_maximum_uses_checked_inverse_columns_and_encloses_dense_solution() {
    with_cx(|cx| {
        contact_model(0.1, |problem, interfaces| {
            let initial = vec![315.0; problem.mesh.vertex_count()];
            let region: Vec<_> = (0..initial.len()).collect();
            let weights = vec![0.0; initial.len()];
            let comparison = LinearGoalAnalyzer::new(
                cx,
                problem,
                Some(interfaces),
                linear(),
                &initial,
                &weights,
                config(),
            )
            .unwrap();
            assert!(
                comparison
                    .analyze(cx, &initial)
                    .unwrap()
                    .enclosure
                    .inverse_infinity_upper()
                    .is_none()
            );
            let analyzer = LinearGoalAnalyzer::new_for_maximum(
                cx,
                problem,
                Some(interfaces),
                linear(),
                &initial,
                config(),
            )
            .unwrap();
            assert_eq!(
                analyzer.inverse_columns().unwrap().len(),
                analyzer.dofs().n()
            );
            let checked = analyzer.analyze_maximum(cx, &initial, &region).unwrap();
            assert!(
                checked.linear_analysis().stability_iterations <= config().max_stability_iterations
            );
            let system = fs_conduction::assemble::assemble_operator_with_interfaces(
                cx,
                problem.mesh,
                problem.boundary,
                problem.material,
                problem.source,
                &initial,
                interfaces,
            )
            .unwrap();
            let (matrix, rhs) = reduce(&system, analyzer.dofs());
            let exact = analyzer.dofs().scatter(&dense_solve(&matrix, &rhs));
            let maximum = exact.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let [lower, upper] = checked.interval_k().unwrap();
            assert!(lower <= maximum && maximum <= upper);
            let result = analyzer
                .solve_maximum_to_goal(
                    cx,
                    &initial,
                    &region,
                    LinearGoalSolveConfig {
                        absolute_tolerance: 1e-8,
                        max_primal_iterations: 256,
                        check_every: 1,
                        max_defect_corrections: 2,
                    },
                )
                .unwrap();
            assert_eq!(result.stop, LinearGoalStop::GoalTolerance);
            let [lower, upper] = result.analysis.interval_k().unwrap();
            assert!(lower <= maximum && maximum <= upper);
            assert!(
                result
                    .temperature
                    .iter()
                    .zip(&exact)
                    .all(|(a, b)| (a - b).abs() < 1e-8)
            );
            assert!(initial.iter().all(|value| *value == 315.0));
        })
    });
}

#[test]
fn inverse_proposals_share_work_and_respect_structural_admission() {
    with_cx(|cx| {
        contact_model(0.1, |problem, interfaces| {
            let initial = vec![315.0; problem.mesh.vertex_count()];
            let region: Vec<_> = (0..initial.len()).collect();
            let system = fs_conduction::assemble::assemble_operator_with_interfaces(
                cx,
                problem.mesh,
                problem.boundary,
                problem.material,
                problem.source,
                &initial,
                interfaces,
            )
            .unwrap();
            let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count()).unwrap();
            let (matrix, _) = reduce(&system, &dofs);
            for (iterations, entries) in [(0, 65536), (3, 65536), (1024, matrix.nnz())] {
                let mut budget = config();
                budget.max_stability_iterations = iterations;
                budget.residual_limits.max_nonzeros = entries;
                let analyzer = LinearGoalAnalyzer::new_for_maximum(
                    cx,
                    problem,
                    Some(interfaces),
                    linear(),
                    &initial,
                    budget,
                )
                .unwrap();
                assert!(analyzer.inverse_columns().is_none());
                let checked = analyzer.analyze_maximum(cx, &initial, &region).unwrap();
                assert!(checked.algebraic_half_width_k().is_none());
                assert!(checked.linear_analysis().stability_iterations <= iterations);
            }
        })
    });
}
