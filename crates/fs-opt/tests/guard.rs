//! Battery for the Goodhart guard (addendum Proposal D). Covers the policy
//! engine (no-steps→provisional, all-pass→cleared→honored, veto→failed+finding,
//! unavailable-step→provisional-never-cleared), the amended "converged AND
//! guard-cleared" contract, determinism, and the concrete δ-perturbation step
//! (smooth optimum passes; found-better and sharp-crack exploits are vetoed;
//! non-finite fails closed).

use fs_opt::{
    DeltaPerturbationStep, DescentReport, Endpoint, EscalationKind, EscalationStep, GoodhartGuard,
    GuardStatus, StepOutcome, converged_and_guard_cleared,
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
        budget_stopped: false,
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

#[test]
fn delta_perturbation_passes_a_smooth_optimum() {
    // f(x) = Σ x_i²  has a smooth minimum at 0; perturbing rises gently.
    let step = delta_step(|x: &[f64]| x.iter().map(|v| v * v).sum());
    let out = step.evaluate(&Endpoint::new("smooth", vec![0.0, 0.0], 0.0));
    assert!(out.is_pass(), "smooth optimum must pass: {out:?}");
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
    // non-finite endpoint objective → veto.
    let step = delta_step(|x: &[f64]| x.iter().sum());
    let out = step.evaluate(&Endpoint::new("nan-endpoint", vec![0.0], f64::NAN));
    assert!(out.is_veto());
    // non-finite objective under perturbation → veto.
    let step2 = delta_step(|x: &[f64]| {
        if x[0].abs() < 1e-9 {
            0.0
        } else {
            f64::INFINITY
        }
    });
    let out2 = step2.evaluate(&Endpoint::new("nan-probe", vec![0.0], 0.0));
    assert!(out2.is_veto());
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
