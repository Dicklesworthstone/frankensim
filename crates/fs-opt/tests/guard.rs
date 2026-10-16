//! Battery for the Goodhart guard (addendum Proposal D). Covers the policy
//! engine (no-steps→provisional, all-pass→cleared→honored, veto→failed+finding,
//! unavailable-step→provisional-never-cleared), the amended "converged AND
//! guard-cleared" contract, determinism, and the concrete δ-perturbation step
//! (smooth optimum passes; found-better and sharp-crack exploits are vetoed;
//! non-finite fails closed).

use fs_opt::{
    DeltaPerturbationStep, DescentReport, DescentStop, Endpoint, EscalationKind, EscalationStep,
    GoodhartGuard, GuardStatus, StepOutcome, converged_and_guard_cleared,
};

/// A fixed-outcome step of a chosen kind, for exercising the aggregator.
struct Stub(EscalationKind, StepOutcome);
impl EscalationStep for Stub {
    fn kind(&self) -> EscalationKind {
        self.0
    }
    fn evaluate(&self, _: &Endpoint) -> StepOutcome {
        self.1.clone()
    }
}

fn endpoint() -> Endpoint {
    Endpoint::new("bracket-opt-endpoint", vec![0.0, 0.0], 0.0)
}

fn all_pass_guard() -> GoodhartGuard {
    let mut g = GoodhartGuard::new();
    for k in EscalationKind::ORDER {
        g = g.with_step(Box::new(Stub(k, StepOutcome::Passed)));
    }
    g
}

#[test]
fn no_steps_registered_is_provisional_not_honored() {
    // The honest baseline: nothing could be checked → provisional, never honored.
    let report = GoodhartGuard::new().evaluate(&endpoint());
    assert_eq!(report.status, GuardStatus::Provisional);
    assert!(!report.is_honored());
    assert_eq!(report.steps.len(), 4);
    assert!(report.steps.iter().all(|s| s.outcome.is_not_performed()));
    assert!(report.findings.is_empty());
}

#[test]
fn all_steps_pass_is_cleared_and_honored() {
    let report = all_pass_guard().evaluate(&endpoint());
    assert_eq!(report.status, GuardStatus::Cleared);
    assert!(report.is_honored(), "{}", report.diagnosis());
    assert!(report.findings.is_empty());
}

#[test]
fn any_veto_is_failed_with_finding() {
    // three pass, cross-representation vetoes.
    let g = GoodhartGuard::new()
        .with_step(Box::new(Stub(
            EscalationKind::RungKPlus1,
            StepOutcome::Passed,
        )))
        .with_step(Box::new(Stub(
            EscalationKind::CrossRepresentation,
            StepOutcome::Vetoed {
                reason: "SDF and mesh paths disagree beyond tolerance".to_string(),
            },
        )))
        .with_step(Box::new(Stub(
            EscalationKind::DeltaPerturbation,
            StepOutcome::Passed,
        )))
        .with_step(Box::new(Stub(
            EscalationKind::EstimatorIndependence,
            StepOutcome::Passed,
        )));
    let report = g.evaluate(&endpoint());
    assert_eq!(report.status, GuardStatus::Failed);
    assert!(!report.is_honored());
    assert_eq!(report.findings.len(), 1);
    assert_eq!(report.findings[0].step, EscalationKind::CrossRepresentation);
    // the finding is treasure: a filable tombstone/bug-report line.
    let s = report.findings[0].summary();
    assert!(s.contains("cross-representation") && s.contains("bracket-opt-endpoint"));
}

#[test]
fn unavailable_step_keeps_provisional_never_cleared() {
    // three steps pass, estimator-independence is NOT registered.
    // The hardening rule: NEVER "cleared" on a skipped check.
    let g = GoodhartGuard::new()
        .with_step(Box::new(Stub(
            EscalationKind::RungKPlus1,
            StepOutcome::Passed,
        )))
        .with_step(Box::new(Stub(
            EscalationKind::CrossRepresentation,
            StepOutcome::Passed,
        )))
        .with_step(Box::new(Stub(
            EscalationKind::DeltaPerturbation,
            StepOutcome::Passed,
        )));
    let report = g.evaluate(&endpoint());
    assert_eq!(report.status, GuardStatus::Provisional);
    assert!(!report.is_honored());
    // the skipped step is recorded as NotPerformed, not silently passed.
    let est = report
        .steps
        .iter()
        .find(|s| s.kind == EscalationKind::EstimatorIndependence)
        .unwrap();
    assert!(est.outcome.is_not_performed());
}

#[test]
fn steps_run_in_fixed_order() {
    let report = GoodhartGuard::new().evaluate(&endpoint());
    let kinds: Vec<EscalationKind> = report.steps.iter().map(|s| s.kind).collect();
    assert_eq!(kinds, EscalationKind::ORDER.to_vec());
}

#[test]
fn amended_contract_requires_both_converged_and_cleared() {
    let cleared = all_pass_guard().evaluate(&endpoint());
    let provisional = GoodhartGuard::new().evaluate(&endpoint());
    assert!(converged_and_guard_cleared(true, &cleared));
    // converged but not guard-cleared → NOT honored.
    assert!(!converged_and_guard_cleared(true, &provisional));
    // guard-cleared but the optimizer did not converge → NOT honored.
    assert!(!converged_and_guard_cleared(false, &cleared));
}

#[test]
fn evaluation_is_deterministic() {
    let g = all_pass_guard();
    let ep = endpoint();
    assert_eq!(g.evaluate(&ep), g.evaluate(&ep));
    let empty = GoodhartGuard::new();
    assert_eq!(empty.evaluate(&ep), empty.evaluate(&ep));
}

#[test]
fn first_registered_step_of_a_kind_is_used() {
    // register DeltaPerturbation twice: Passed first, Vetoed second.
    // The FIRST wins → not vetoed (deterministic registration).
    let g = GoodhartGuard::new()
        .with_step(Box::new(Stub(
            EscalationKind::DeltaPerturbation,
            StepOutcome::Passed,
        )))
        .with_step(Box::new(Stub(
            EscalationKind::DeltaPerturbation,
            StepOutcome::Vetoed {
                reason: "should not be reached".to_string(),
            },
        )));
    let report = g.evaluate(&endpoint());
    let delta = report
        .steps
        .iter()
        .find(|s| s.kind == EscalationKind::DeltaPerturbation)
        .unwrap();
    assert!(delta.outcome.is_pass());
    assert!(report.findings.is_empty());
}

#[test]
fn from_descent_builds_endpoint() {
    let report = DescentReport {
        x: vec![1.5, -2.0],
        f0: 9.0,
        f_final: 0.25,
        evals: 42,
        steps_taken: 10,
        stop: DescentStop::StepLimit,
        budget_stopped: false,
        work_upper_bound: 1_024,
        workspace_upper_bound_bytes: 4_096,
    };
    let ep = Endpoint::from_descent("study-node-7", &report);
    assert_eq!(ep.design, vec![1.5, -2.0]);
    assert_eq!(ep.objective.to_bits(), 0.25f64.to_bits());
    assert_eq!(ep.label, "study-node-7");
}

// ---- the concrete δ-perturbation step ------------------------------------

fn delta_step<F: Fn(&[f64]) -> f64 + 'static>(f: F) -> Box<dyn EscalationStep> {
    // better_tol tiny (any real improvement is a veto); sharpness_tol=1.0.
    Box::new(DeltaPerturbationStep::new(0.1, 1e-6, 1.0, f))
}

fn assert_veto_reason(outcome: StepOutcome, expected: &str) {
    let StepOutcome::Vetoed { reason } = outcome else {
        panic!("expected a fail-closed veto, got {outcome:?}");
    };
    assert!(
        reason.contains(expected),
        "veto reason {reason:?} must identify {expected:?}"
    );
}

#[test]
fn delta_perturbation_passes_a_smooth_optimum() {
    // f(x) = Σ x_i²  has a smooth minimum at 0; perturbing rises gently.
    let step = delta_step(|x: &[f64]| x.iter().map(|v| v * v).sum());
    let out = step.evaluate(&Endpoint::new("smooth", vec![0.0, 0.0], 0.0));
    assert!(out.is_pass(), "smooth optimum must pass: {out:?}");
}

#[test]
fn delta_perturbation_accepts_inclusive_tolerance_boundaries() {
    let step = DeltaPerturbationStep::new(0.25, 0.25, 0.5, |x: &[f64]| {
        if x[0].is_sign_positive() { -0.25 } else { 0.5 }
    });
    let out = step.evaluate(&Endpoint::new("inclusive-boundaries", vec![0.0], 0.0));
    assert!(
        out.is_pass(),
        "equality with either tolerance must pass: {out:?}"
    );

    for zero in [0.0, -0.0] {
        let flat = DeltaPerturbationStep::new(0.25, zero, zero, |_: &[f64]| 0.0);
        let out = flat.evaluate(&Endpoint::new("zero-tolerances", vec![0.0], 0.0));
        assert!(out.is_pass(), "flat probes satisfy zero tolerance: {out:?}");
    }
}

#[test]
fn delta_perturbation_vetoes_a_found_better_point() {
    // the endpoint claims objective 0 at x=0, but honest re-eval nearby is
    // LOWER — the endpoint was not a true optimum (a discretization artifact).
    let step = delta_step(|x: &[f64]| if x[0].abs() < 1e-9 { 0.0 } else { -1.0 });
    let out = step.evaluate(&Endpoint::new("found-better", vec![0.0], 0.0));
    assert!(out.is_veto(), "a nearby better point must veto: {out:?}");
}

#[test]
fn delta_perturbation_vetoes_a_sharp_crack() {
    // the endpoint sits on a downward spike: honest re-eval nearby jumps UP.
    let step = delta_step(|x: &[f64]| if x[0].abs() < 1e-9 { 0.0 } else { 1.0e6 });
    let out = step.evaluate(&Endpoint::new("crack", vec![0.0], 0.0));
    assert!(out.is_veto(), "an optimum in a crack must veto: {out:?}");
}

#[test]
fn delta_perturbation_fails_closed_on_nonfinite() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let step = delta_step(|x: &[f64]| x.iter().sum());
        assert_veto_reason(
            step.evaluate(&Endpoint::new("nonfinite-endpoint", vec![0.0], value)),
            "endpoint objective",
        );

        let step = delta_step(move |_: &[f64]| value);
        assert_veto_reason(
            step.evaluate(&Endpoint::new("nonfinite-probe", vec![0.0], 0.0)),
            "objective is non-finite",
        );
    }
}

#[test]
fn delta_perturbation_refuses_malformed_numeric_policy() {
    for (delta, better_tol, sharpness_tol, expected) in [
        (f64::NAN, 0.0, 1.0, "delta"),
        (0.0, 0.0, 1.0, "delta"),
        (-0.0, 0.0, 1.0, "delta"),
        (-0.1, 0.0, 1.0, "delta"),
        (f64::INFINITY, 0.0, 1.0, "delta"),
        (f64::NEG_INFINITY, 0.0, 1.0, "delta"),
        (0.1, f64::NAN, 1.0, "better_tol"),
        (0.1, -1.0, 1.0, "better_tol"),
        (0.1, f64::INFINITY, 1.0, "better_tol"),
        (0.1, f64::NEG_INFINITY, 1.0, "better_tol"),
        (0.1, 0.0, f64::NAN, "sharpness_tol"),
        (0.1, 0.0, f64::INFINITY, "sharpness_tol"),
        (0.1, 0.0, f64::NEG_INFINITY, "sharpness_tol"),
        (0.1, 0.0, -1.0, "sharpness_tol"),
    ] {
        let step = DeltaPerturbationStep::new(delta, better_tol, sharpness_tol, |_: &[f64]| 0.0);
        assert_veto_reason(
            step.evaluate(&Endpoint::new("malformed-policy", vec![0.0], 0.0)),
            expected,
        );
    }
}

#[test]
fn delta_perturbation_refuses_nonfinite_design_even_when_objective_ignores_it() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let step = delta_step(|_: &[f64]| 0.0);
        assert_veto_reason(
            step.evaluate(&Endpoint::new("nonfinite-design", vec![value], 0.0)),
            "design coord 0",
        );
    }
}

#[test]
fn delta_perturbation_refuses_nonfinite_probe_and_difference() {
    let overflowing_probe =
        DeltaPerturbationStep::new(f64::MAX, f64::MAX, f64::MAX, |_: &[f64]| 0.0);
    assert_veto_reason(
        overflowing_probe.evaluate(&Endpoint::new("probe-overflow", vec![f64::MAX], 0.0)),
        "makes coord 0 non-finite",
    );
    assert_veto_reason(
        overflowing_probe.evaluate(&Endpoint::new(
            "negative-probe-overflow",
            vec![-f64::MAX],
            0.0,
        )),
        "makes coord 0 non-finite",
    );

    let overflowing_difference =
        DeltaPerturbationStep::new(0.1, f64::MAX, f64::MAX, |_: &[f64]| f64::MAX);
    assert_veto_reason(
        overflowing_difference.evaluate(&Endpoint::new(
            "difference-overflow",
            vec![0.0],
            -f64::MAX,
        )),
        "objective difference is non-finite",
    );
}

#[test]
fn delta_perturbation_refuses_a_radius_absorbed_by_coordinate_scale() {
    let objective_invoked = std::cell::Cell::new(false);
    let step = DeltaPerturbationStep::new(1.0, 0.0, 0.0, |_: &[f64]| {
        objective_invoked.set(true);
        0.0
    });
    assert_veto_reason(
        step.evaluate(&Endpoint::new("rounded-away", vec![f64::MAX], 0.0)),
        "does not change coord 0",
    );
    assert!(
        !objective_invoked.get(),
        "an unchanged probe must be refused before objective evaluation"
    );
}

#[test]
fn empty_design_is_vacuously_robust() {
    let step = delta_step(|_: &[f64]| 0.0);
    let out = step.evaluate(&Endpoint::new("scalarless", vec![], 0.0));
    assert!(out.is_pass());
}

#[test]
fn realistic_v0_delta_only_is_provisional() {
    // The honest v0 reality: only δ-perturbation machinery exists. It passes a
    // smooth optimum, but the other three steps are NotPerformed → the endpoint
    // is PROVISIONAL, not honored. The guard never over-claims.
    let g = GoodhartGuard::new().with_step(delta_step(|x: &[f64]| x.iter().map(|v| v * v).sum()));
    let report = g.evaluate(&Endpoint::new("v0", vec![0.0, 0.0], 0.0));
    assert_eq!(report.status, GuardStatus::Provisional);
    assert!(!report.is_honored());
    let delta = report
        .steps
        .iter()
        .find(|s| s.kind == EscalationKind::DeltaPerturbation)
        .unwrap();
    assert!(delta.outcome.is_pass());
    assert_eq!(
        report
            .steps
            .iter()
            .filter(|s| s.outcome.is_not_performed())
            .count(),
        3
    );
}
