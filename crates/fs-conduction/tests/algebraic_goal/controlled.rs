//! Goal-controlled use of the actual assembled thermal analyzer.
use super::*;
use fs_conduction::adjoint::{LinearGoalSolveConfig, LinearGoalStop};

fn policy() -> LinearGoalSolveConfig {
    LinearGoalSolveConfig {
        absolute_tolerance: 1e-8, max_primal_iterations: 128,
        check_every: 1, max_defect_corrections: 2,
    }
}

#[test]
fn loose_goal_skips_primal_work_while_strict_goal_solves_the_physical_field() {
    let f = Fixture::new(true, false);
    let initial = f.initial();
    with_cx(|cx| {
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &f.weights(), config()).unwrap();
        let mut loose = policy(); loose.absolute_tolerance = 1.0;
        let skipped = analyzer.solve_to_goal(cx, &initial, loose).unwrap();
        assert_eq!(skipped.stop, LinearGoalStop::GoalTolerance);
        assert_eq!(skipped.primal_iterations, 0);
        assert_eq!(skipped.goal_checks, 1);
        assert_eq!(skipped.temperature, initial);
        // The field was NOT declared converged just because its goal passed.
        assert!(skipped.analysis.enclosure.primal_residual_infinity_upper() > 0.01);
        let solved = analyzer.solve_to_goal(cx, &initial, policy()).unwrap();
        assert_eq!(solved.stop, LinearGoalStop::GoalTolerance);
        assert_eq!(solved.primal_iterations, 1);
        assert_eq!(solved.goal_checks, 2);
        assert!(solved.analysis.meets_absolute_tolerance(policy().absolute_tolerance));
        assert_eq!(solved.analysis, analyzer.analyze(cx, &solved.temperature).unwrap());
        let exact = oracle(&f, cx);
        for (&a,&b) in solved.temperature.iter().zip(&exact) { assert!((a-b).abs() < 1e-10); }
        for &v in analyzer.dofs().fixed() { assert_eq!(solved.temperature[v].to_bits(), initial[v].to_bits()); }
    });
}

#[test]
fn no_inverse_and_zero_iteration_budget_never_claim_the_goal() {
    let f = Fixture::new(false, false);
    with_cx(|cx| {
        let initial = f.initial();
        let mut without_scaling = config(); without_scaling.max_stability_iterations = 0;
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &f.weights(), without_scaling).unwrap();
        let missing = analyzer.solve_to_goal(cx, &initial, policy()).unwrap();
        assert_eq!(missing.stop, LinearGoalStop::BoundUnavailable);
        assert_eq!(missing.primal_iterations, 0);
        assert!(!missing.analysis.meets_absolute_tolerance(policy().absolute_tolerance));
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &f.weights(), config()).unwrap();
        let mut zero = policy(); zero.max_primal_iterations = 0;
        let exhausted = analyzer.solve_to_goal(cx, &initial, zero).unwrap();
        assert_eq!(exhausted.stop, LinearGoalStop::IterationBudget);
        assert_eq!(exhausted.primal_iterations, 0);
        assert_eq!(exhausted.temperature, initial);
    });
}

#[test]
fn primal_work_and_corrections_obey_the_single_shared_budget() {
    let f = Fixture::new(false, false);
    with_cx(|cx| {
        let initial = f.initial();
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &f.weights(), config()).unwrap();
        for budget in [1, 2, 5] {
            let mut limited = policy(); limited.max_primal_iterations = budget;
            limited.max_defect_corrections = 1;
            limited.absolute_tolerance = 1e-40;
            let result = analyzer.solve_to_goal(cx, &initial, limited).unwrap();
            assert_ne!(result.stop, LinearGoalStop::GoalTolerance);
            assert!(result.primal_iterations <= budget);
            assert!(result.defect_corrections <= 1);
            assert!(!result.analysis.meets_absolute_tolerance(limited.absolute_tolerance));
            assert_eq!(result.analysis, analyzer.analyze(cx, &result.temperature).unwrap());
        }
    });
}

#[test]
fn rounded_residual_zero_is_not_an_arbitrarily_precise_goal_certificate() {
    let f = Fixture::new(true, false);
    with_cx(|cx| {
        let initial = f.initial();
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &f.weights(), config()).unwrap();
        let mut impossible = policy(); impossible.absolute_tolerance = 1e-40;
        let result = analyzer.solve_to_goal(cx, &initial, impossible).unwrap();
        assert_ne!(result.stop, LinearGoalStop::GoalTolerance);
        assert!(result.analysis.enclosure.goal_error().unwrap().magnitude_upper() > impossible.absolute_tolerance);
        assert!(result.primal_iterations <= impossible.max_primal_iterations);
        assert_eq!(result.analysis, analyzer.analyze(cx, &result.temperature).unwrap());
    });
}

#[test]
fn signed_extreme_goal_scales_preserve_the_solved_temperature() {
    let f = Fixture::new(true, false);
    with_cx(|cx| {
        let initial = f.initial();
        let weights = f.weights();
        for factor in [1e-100, -2.0, 1e100] {
            let scaled: Vec<_> = weights.iter().map(|w| w*factor).collect();
            let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &scaled, config()).unwrap();
            let mut cfg = policy(); cfg.absolute_tolerance *= factor.abs();
            let result = analyzer.solve_to_goal(cx, &initial, cfg).unwrap();
            assert_eq!(result.stop, LinearGoalStop::GoalTolerance);
            assert_eq!(result.primal_iterations, 1);
            assert!(result.analysis.meets_absolute_tolerance(cfg.absolute_tolerance));
            let exact = oracle(&f, cx);
            for (&a,&b) in result.temperature.iter().zip(&exact) { assert!((a-b).abs() < 1e-10); }
        }
    });
}

#[test]
fn cancellation_after_a_successful_candidate_check_prevents_publication() {
    let f = Fixture::new(true, false);
    let initial = f.initial();
    let analyzer = with_cx(|cx| LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &f.weights(), config()).unwrap());
    support::with_gate(|gate,cx| {
        let mut saw_solved_candidate = false;
        let result = analyzer.solve_to_goal_observed(cx, &initial, policy(), |iteration, report| {
            if iteration > 0 {
                saw_solved_candidate = report.meets_absolute_tolerance(policy().absolute_tolerance);
                gate.request();
            }
        });
        assert!(saw_solved_candidate);
        assert!(matches!(result, Err(ConductionError::Cancelled { .. })));
    });
    let first = with_cx(|cx| analyzer.solve_to_goal(cx, &initial, policy()).unwrap());
    let again = with_cx(|cx| analyzer.solve_to_goal(cx, &initial, policy()).unwrap());
    assert_eq!(first, again);
    assert!(initial.iter().all(|&v| v == 300.0));
}

#[test]
fn malformed_goal_policy_is_refused_before_the_observer() {
    let f = Fixture::new(true, false);
    with_cx(|cx| {
        let initial = f.initial();
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &f.weights(), config()).unwrap();
        for (tolerance, check_every) in [(0.0,1), (-1.0,1), (f64::NAN,1), (f64::INFINITY,1), (1e-8,0), (1e-8,33)] {
            let mut cfg = policy(); cfg.absolute_tolerance = tolerance; cfg.check_every = check_every;
            let mut observed = false;
            assert!(matches!(analyzer.solve_to_goal_observed(cx, &initial, cfg, |_,_| observed = true), Err(ConductionError::Config { .. })));
            assert!(!observed);
        }
    });
}

#[test]
fn a_budget_limited_inexact_dual_retains_its_error_instead_of_refusing_analysis() {
    let f = Fixture::new(false, false);
    with_cx(|cx| {
        let initial = f.initial();
        let mut cheap = linear(); cheap.max_iterations = 1; cheap.tolerance = f64::EPSILON;
        let mut weights = f.weights();
        for (i,w) in weights.iter_mut().enumerate() { *w *= (i as f64 + 1.0) / 7.0; }
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, cheap, &initial, &weights, config()).unwrap();
        let report = analyzer.analyze(cx, &initial).unwrap();
        assert!(report.dual_iterations <= 1);
        assert!(report.dual_relative_residual.is_finite());
        assert!(report.enclosure.dual_error_upper().unwrap() > 0.0);
        let exact = oracle(&f, cx);
        let error: f64 = weights.iter().zip(exact.iter().zip(&initial)).map(|(w,(a,b))| w*(a-b)).sum();
        let bound = report.enclosure.goal_error().unwrap();
        assert!(bound.lower() <= error && error <= bound.upper());
        let result = analyzer.solve_to_goal(cx, &initial, policy()).unwrap();
        assert_eq!(result.stop, LinearGoalStop::GoalTolerance);
        assert!(result.analysis.meets_absolute_tolerance(policy().absolute_tolerance));
        assert_eq!(result.analysis.dual_iterations, report.dual_iterations);
    });
}
