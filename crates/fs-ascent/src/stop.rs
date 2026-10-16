//! The stopping-condition algebra (Appendix C): `until (any …)` /
//! `until (all …)` combinators over gradient norms, objective drops,
//! evaluation budgets, and stagnation — budget-aware, resumable, and
//! attributable (a stop REASON travels with every report, so "why did
//! the study end" is data, not archaeology).

/// A single stopping criterion.
#[derive(Debug, Clone)]
pub enum StopRule {
    /// ‖g‖∞ below the threshold.
    GradNorm(f64),
    /// Objective at or below the target value.
    ObjectiveBelow(f64),
    /// Total function/gradient evaluations at or above the budget.
    Budget(usize),
    /// No relative objective improvement above `rel` for `window`
    /// consecutive iterations.
    Stall {
        /// Relative-improvement floor.
        rel: f64,
        /// Consecutive-iteration window.
        window: usize,
    },
    /// Satisfied when ANY child is.
    Any(Vec<StopRule>),
    /// Satisfied when ALL children are.
    All(Vec<StopRule>),
}

/// Why an optimization stopped (attribution, not archaeology).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StopReason {
    /// Gradient-norm criterion met.
    GradNorm,
    /// Objective target reached.
    ObjectiveBelow,
    /// Evaluation budget exhausted.
    Budget,
    /// Stagnation window triggered.
    Stall,
    /// Composite (first satisfied child of an Any / the All itself).
    Composite,
    /// The iteration cap of the driving loop (no rule fired).
    IterationCap,
}

/// Observable state a rule is checked against.
#[derive(Debug, Clone, Copy)]
pub struct StopObservation<'h> {
    /// Current ‖g‖∞.
    pub grad_norm: f64,
    /// Current objective.
    pub objective: f64,
    /// Evaluations spent so far.
    pub evals: usize,
    /// Objective history (most recent last).
    pub history: &'h [f64],
}

impl StopRule {
    /// Check the rule; `Some(reason)` when satisfied.
    #[must_use]
    pub fn check(&self, obs: &StopObservation<'_>) -> Option<StopReason> {
        match self {
            StopRule::GradNorm(t) => (obs.grad_norm <= *t).then_some(StopReason::GradNorm),
            StopRule::ObjectiveBelow(t) => {
                (obs.objective <= *t).then_some(StopReason::ObjectiveBelow)
            }
            StopRule::Budget(b) => (obs.evals >= *b).then_some(StopReason::Budget),
            StopRule::Stall { rel, window } => {
                if obs.history.len() < window + 1 {
                    return None;
                }
                let now = *obs.history.last().expect("nonempty");
                let then = obs.history[obs.history.len() - 1 - window];
                let improved = (then - now) > rel * then.abs().max(1e-30);
                (!improved).then_some(StopReason::Stall)
            }
            StopRule::Any(rules) => rules.iter().find_map(|r| r.check(obs)),
            StopRule::All(rules) => rules
                .iter()
                .all(|r| r.check(obs).is_some())
                .then_some(StopReason::Composite),
        }
    }
}
