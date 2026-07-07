//! The Goodhart guard (plan addendum, Proposal D): treat every optimizer
//! endpoint as an ADVERSARIAL EXAMPLE by default.
//!
//! Adjoint optimizers converge to nonphysical exploits of the discretization,
//! and an agent operator is, in the limit, a very creative adjoint optimizer:
//! it finds the cracks in the model and calls them design. A human designer
//! wants tools that trust them; an agent operator needs tools that DON'T. So
//! the system distrusts its outputs most at exactly the points the optimizer
//! is proudest of (design principle P5): a converged design's certificate is
//! **provisional** until an escalation ladder clears it.
//!
//! This module is the POLICY ENGINE (Phase-0 scaffold). It runs a FIXED
//! four-step escalation in order — rung k+1, cross-representation re-solve,
//! δ-perturbation robustness, estimator independence — where each step is a
//! pluggable [`EscalationStep`]. Steps whose machinery does not yet exist
//! honestly report [`StepOutcome::NotPerformed`], and the endpoint stays
//! [`GuardStatus::Provisional`] — never reported "cleared" on a skipped
//! check. A step that VETOES makes the endpoint [`GuardStatus::Failed`]; the
//! veto is TREASURE — captured as a [`GuardFinding`] the caller files as a
//! tombstone + estimator/discretization bug report (this module is L4 and
//! does not call HELM/the ledger itself — it produces the findings).
//!
//! One concrete step ships today: [`DeltaPerturbationStep`], which needs no
//! external machinery — it re-evaluates the objective at small deterministic
//! perturbations of the design and vetoes an endpoint that is not robust (an
//! optimum living in a crack is sharp; a real optimum is smooth at the
//! certificate's scale). The other three steps are injected by callers once
//! the fidelity ladder, Rep Router, and multiple estimator families exist.
//!
//! Determinism: the ladder runs steps in fixed order and the perturbation
//! probes are the fixed ±δ coordinate directions (no RNG), so a replayed
//! guard evaluation reproduces the same report.

use crate::eval::DescentReport;

/// The four escalation steps, in the fixed order the guard runs them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EscalationKind {
    /// Re-solve at the next fidelity rung (needs the fidelity-ladder registry).
    RungKPlus1,
    /// Re-route through a DIFFERENT representation via the Rep Router and
    /// re-solve — cross-representation agreement is an independent check only
    /// a multi-representation system can perform.
    CrossRepresentation,
    /// The optimum must persist, within bounds, under small perturbations of
    /// the design. Shippable today (see [`DeltaPerturbationStep`]).
    DeltaPerturbation,
    /// Re-verify with a DIFFERENT a-posteriori estimator family.
    EstimatorIndependence,
}

impl EscalationKind {
    /// The fixed escalation order.
    pub const ORDER: [EscalationKind; 4] = [
        EscalationKind::RungKPlus1,
        EscalationKind::CrossRepresentation,
        EscalationKind::DeltaPerturbation,
        EscalationKind::EstimatorIndependence,
    ];

    /// Stable machine-readable name (for structured logging).
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            EscalationKind::RungKPlus1 => "rung-k+1",
            EscalationKind::CrossRepresentation => "cross-representation",
            EscalationKind::DeltaPerturbation => "delta-perturbation",
            EscalationKind::EstimatorIndependence => "estimator-independence",
        }
    }
}

/// The verdict of one escalation step.
#[derive(Debug, Clone, PartialEq)]
pub enum StepOutcome {
    /// The step ran and the endpoint survived it.
    Passed,
    /// The step ran and REJECTED the endpoint (a caught exploit).
    Vetoed {
        /// Why the endpoint was rejected.
        reason: String,
    },
    /// The step could not run (its machinery is not registered). The
    /// certificate stays provisional — this is NOT a pass.
    NotPerformed {
        /// Why the step could not run.
        reason: String,
    },
}

impl StepOutcome {
    /// Did the endpoint survive this step?
    #[must_use]
    pub fn is_pass(&self) -> bool {
        matches!(self, StepOutcome::Passed)
    }

    /// Did this step reject the endpoint?
    #[must_use]
    pub fn is_veto(&self) -> bool {
        matches!(self, StepOutcome::Vetoed { .. })
    }

    /// Was this step skipped for lack of machinery?
    #[must_use]
    pub fn is_not_performed(&self) -> bool {
        matches!(self, StepOutcome::NotPerformed { .. })
    }
}

/// A pluggable escalation step. Implementations bring the actual machinery
/// (a rung solve, a Router re-route, an estimator family); the guard only
/// sequences them and interprets their verdicts.
pub trait EscalationStep {
    /// Which escalation slot this step fills.
    fn kind(&self) -> EscalationKind;
    /// Judge an endpoint.
    fn evaluate(&self, endpoint: &Endpoint) -> StepOutcome;
}

/// An optimizer endpoint under scrutiny: a converged design and its objective
/// value, plus a provenance label. Decoupled from any particular optimizer.
#[derive(Debug, Clone, PartialEq)]
pub struct Endpoint {
    /// Provenance label (design id / study node) for logging + findings.
    pub label: String,
    /// The converged design point.
    pub design: Vec<f64>,
    /// The objective value the optimizer reported at the endpoint.
    pub objective: f64,
}

impl Endpoint {
    /// A labelled endpoint.
    #[must_use]
    pub fn new(label: impl Into<String>, design: Vec<f64>, objective: f64) -> Endpoint {
        Endpoint {
            label: label.into(),
            design,
            objective,
        }
    }

    /// The endpoint of an fs-opt descent (`x`, `f_final`).
    #[must_use]
    pub fn from_descent(label: impl Into<String>, report: &DescentReport) -> Endpoint {
        Endpoint::new(label, report.x.clone(), report.f_final)
    }
}

/// The overall guard verdict for an endpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GuardStatus {
    /// Every one of the four steps ran and PASSED — the certificate may be
    /// honored (this is the only honored state).
    Cleared,
    /// No step vetoed, but at least one could not run — the certificate stays
    /// provisional (honest incompleteness, never a false clear).
    Provisional,
    /// At least one step VETOED — a caught exploit.
    Failed,
}

/// The captured record of a veto — "treasure" the caller files as a tombstone
/// (Proposal E) AND an estimator/discretization bug report (Proposal 6).
#[derive(Debug, Clone, PartialEq)]
pub struct GuardFinding {
    /// The offending endpoint's label.
    pub endpoint: String,
    /// Which step caught it.
    pub step: EscalationKind,
    /// The step's reason.
    pub reason: String,
}

impl GuardFinding {
    /// A one-line summary for a tombstone / bug-report row.
    #[must_use]
    pub fn summary(&self) -> String {
        format!(
            "goodhart-veto[{}] on endpoint '{}': {}",
            self.step.name(),
            self.endpoint,
            self.reason
        )
    }
}

/// The full guard report: the per-step verdicts in fixed order, the overall
/// status, and any findings.
#[derive(Debug, Clone, PartialEq)]
pub struct GuardReport {
    /// The endpoint's label.
    pub endpoint: String,
    /// The overall verdict.
    pub status: GuardStatus,
    /// Every step's verdict, in [`EscalationKind::ORDER`].
    pub steps: Vec<StepReport>,
    /// Findings produced by vetoes (empty unless `status == Failed`).
    pub findings: Vec<GuardFinding>,
}

/// One step's verdict inside a [`GuardReport`].
#[derive(Debug, Clone, PartialEq)]
pub struct StepReport {
    /// Which step.
    pub kind: EscalationKind,
    /// Its outcome.
    pub outcome: StepOutcome,
}

impl GuardReport {
    /// May the endpoint's certificate be honored? True ONLY when every step
    /// passed (`Cleared`). Provisional and Failed are both un-honored.
    #[must_use]
    pub fn is_honored(&self) -> bool {
        self.status == GuardStatus::Cleared
    }

    /// A one-line structured diagnosis for logging (never printed to stdout by
    /// library code — the caller decides).
    #[must_use]
    pub fn diagnosis(&self) -> String {
        let performed = self.steps.iter().filter(|s| s.outcome.is_pass()).count();
        let skipped = self
            .steps
            .iter()
            .filter(|s| s.outcome.is_not_performed())
            .count();
        format!(
            "guard[{}] endpoint '{}': {performed}/4 passed, {skipped} not-performed, {} veto(s)",
            match self.status {
                GuardStatus::Cleared => "cleared",
                GuardStatus::Provisional => "provisional",
                GuardStatus::Failed => "failed",
            },
            self.endpoint,
            self.findings.len()
        )
    }
}

/// The Goodhart guard: a registry of escalation steps that sequences them
/// over an endpoint. Register at most one step per [`EscalationKind`]; an
/// unregistered kind is reported [`StepOutcome::NotPerformed`].
#[derive(Default)]
pub struct GoodhartGuard {
    steps: Vec<Box<dyn EscalationStep>>,
}

impl GoodhartGuard {
    /// A guard with no steps registered — every endpoint comes back
    /// `Provisional` (nothing could be checked). The honest baseline.
    #[must_use]
    pub fn new() -> GoodhartGuard {
        GoodhartGuard { steps: Vec::new() }
    }

    /// Register an escalation step (builder style). A later registration of
    /// the same kind shadows an earlier one is NOT allowed silently — the
    /// FIRST registered step of a kind is the one used, so registration order
    /// is deterministic and explicit.
    #[must_use]
    pub fn with_step(mut self, step: Box<dyn EscalationStep>) -> GoodhartGuard {
        self.steps.push(step);
        self
    }

    /// Evaluate an endpoint through the fixed four-step ladder. Deterministic:
    /// steps run in [`EscalationKind::ORDER`]; the first registered step of
    /// each kind is used; unregistered kinds are `NotPerformed`.
    #[must_use]
    pub fn evaluate(&self, endpoint: &Endpoint) -> GuardReport {
        let mut steps = Vec::with_capacity(4);
        let mut findings = Vec::new();
        let mut any_veto = false;
        let mut any_skipped = false;
        for kind in EscalationKind::ORDER {
            let outcome = match self.steps.iter().find(|s| s.kind() == kind) {
                Some(step) => step.evaluate(endpoint),
                None => StepOutcome::NotPerformed {
                    reason: format!("no {} capability registered", kind.name()),
                },
            };
            match &outcome {
                StepOutcome::Vetoed { reason } => {
                    any_veto = true;
                    findings.push(GuardFinding {
                        endpoint: endpoint.label.clone(),
                        step: kind,
                        reason: reason.clone(),
                    });
                }
                StepOutcome::NotPerformed { .. } => any_skipped = true,
                StepOutcome::Passed => {}
            }
            steps.push(StepReport { kind, outcome });
        }
        let status = if any_veto {
            GuardStatus::Failed
        } else if any_skipped {
            GuardStatus::Provisional
        } else {
            GuardStatus::Cleared
        };
        GuardReport {
            endpoint: endpoint.label.clone(),
            status,
            steps,
            findings,
        }
    }
}

/// The amended optimization contract: "converged" is redefined as "converged
/// AND guard-cleared". An endpoint may be honored only if the optimizer
/// actually converged AND the guard cleared every step.
#[must_use]
pub fn converged_and_guard_cleared(descent_converged: bool, report: &GuardReport) -> bool {
    descent_converged && report.is_honored()
}

/// δ-PERTURBATION robustness — the escalation step that needs no external
/// machinery. It re-evaluates a supplied objective at small deterministic
/// perturbations (`±δ` along each coordinate) of the endpoint and vetoes an
/// endpoint that is not robust: an optimum living in a discretization crack is
/// sharp (the honest re-evaluation jumps), and a point that a perturbation
/// beats is not a true optimum. A smooth true optimum passes.
pub struct DeltaPerturbationStep<F> {
    /// The perturbation radius.
    pub delta: f64,
    /// How much LOWER a perturbed objective may be before the endpoint is
    /// judged "not actually a minimum" (a found-better exploit).
    pub better_tol: f64,
    /// How much HIGHER a perturbed objective may rise before the endpoint is
    /// judged "sharp / living in a crack".
    pub sharpness_tol: f64,
    /// The objective, re-evaluated HONESTLY (whatever an exploit gamed at the
    /// endpoint, this closure evaluates the real value at nearby designs).
    objective: F,
}

impl<F> DeltaPerturbationStep<F>
where
    F: Fn(&[f64]) -> f64,
{
    /// A δ-perturbation step with explicit tolerances and an objective.
    pub fn new(delta: f64, better_tol: f64, sharpness_tol: f64, objective: F) -> Self {
        DeltaPerturbationStep {
            delta,
            better_tol,
            sharpness_tol,
            objective,
        }
    }
}

impl<F> EscalationStep for DeltaPerturbationStep<F>
where
    F: Fn(&[f64]) -> f64,
{
    fn kind(&self) -> EscalationKind {
        EscalationKind::DeltaPerturbation
    }

    fn evaluate(&self, endpoint: &Endpoint) -> StepOutcome {
        // Fail closed on a non-finite endpoint objective.
        if !endpoint.objective.is_finite() {
            return StepOutcome::Vetoed {
                reason: "endpoint objective is non-finite".to_string(),
            };
        }
        // A 0-dim design has nothing to perturb — vacuously robust.
        if endpoint.design.is_empty() {
            return StepOutcome::Passed;
        }
        let f0 = endpoint.objective;
        let mut worst_lower = 0.0_f64; // most-negative (perturbation beat the endpoint)
        let mut worst_higher = 0.0_f64; // most-positive (endpoint sits on a spike)
        for k in 0..endpoint.design.len() {
            for signed in [self.delta, -self.delta] {
                let mut probe = endpoint.design.clone();
                probe[k] += signed;
                let fp = (self.objective)(&probe);
                if !fp.is_finite() {
                    return StepOutcome::Vetoed {
                        reason: format!(
                            "objective is non-finite under a δ={} perturbation of coord {k}",
                            self.delta
                        ),
                    };
                }
                let diff = fp - f0;
                worst_lower = worst_lower.min(diff);
                worst_higher = worst_higher.max(diff);
            }
        }
        if worst_lower < -self.better_tol {
            return StepOutcome::Vetoed {
                reason: format!(
                    "a δ-perturbation found a better objective ({worst_lower:.3e} below the endpoint): \
                     the endpoint is not a true optimum (likely a discretization artifact)"
                ),
            };
        }
        if worst_higher > self.sharpness_tol {
            return StepOutcome::Vetoed {
                reason: format!(
                    "the objective rises sharply ({worst_higher:.3e}) under a δ-perturbation: \
                     the optimum lives in a crack, not a smooth basin"
                ),
            };
        }
        StepOutcome::Passed
    }
}
