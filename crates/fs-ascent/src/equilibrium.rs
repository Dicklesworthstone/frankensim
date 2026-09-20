//! Bounded physical inverse design using the existing small-dense SQP engine.
//!
//! `fs-couple` owns every preload solve and implicit equilibrium adjoint; this
//! module supplies their results, physical constraints and parameter bounds to SQP. It
//! does not differentiate solver iterations, relax contact admission, or infer
//! material identifiability from a small optimization residual.

use crate::sqp::{SqpError, SqpRunReport, SqpSample, SqpState, SqpStop};
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::{
    DesignControl, DesignError, DesignEvaluation, DesignWork, EquilibriumDesign,
    constraints::ConstraintSense,
};
use fs_exec::CancelGate;

/// Original physics errors remain inspectable through `SqpError::Evaluation`.
/// Cancellation, including cancellation inside physics, has one attribution.
pub type EquilibriumStudyError = SqpError<DesignError>;

fn poll(gate: &CancelGate) -> Result<(), EquilibriumStudyError> {
    if gate.is_requested() { Err(SqpError::Cancelled) } else { Ok(()) }
}

fn normalize_error(error: EquilibriumStudyError) -> EquilibriumStudyError {
    match error {
        SqpError::Evaluation(DesignError::Cancelled) => SqpError::Cancelled,
        other => other,
    }
}

fn admit(problem: &EquilibriumDesign, maximum_kkt_dimension: usize)
    -> Result<(), EquilibriumStudyError>
{
    let n = problem.variables().len();
    if n == 0 || n.checked_mul(3).and_then(|d| d.checked_add(problem.constraints().len()))
        .is_none_or(|d| d > maximum_kkt_dimension) {
        return Err(SqpError::Invalid("decisions, box faces and physical constraints exceed the dense KKT cap"));
    }
    for variable in problem.variables() {
        // All constraint residuals must be representable in decision units.
        // Do not silently rescale the caller's physical design coordinates.
        let width = (variable.maximum - variable.minimum) / variable.scale;
        if !width.is_finite() || width <= 0.0 {
            return Err(SqpError::Invalid("physical bound width is not representable in decision units"));
        }
    }
    Ok(())
}

fn sample(problem: &EquilibriumDesign, evaluation: &DesignEvaluation) -> SqpSample {
    let n = problem.variables().len();
    let mut ci = Vec::with_capacity(2 * n);
    let mut ji = vec![0.0; 2 * n * n];
    for (i, (variable, value)) in problem.variables().iter()
        .zip(&evaluation.physical_parameters).enumerate()
    {
        // p = reference + scale*x. Dividing residuals by scale gives exact
        // +/-1 Jacobian rows and multipliers in dimensionless decision units.
        // Keep the signed residuals; clipping them destroys complementarity.
        ci.push((variable.minimum - value) / variable.scale);
        ci.push((value - variable.maximum) / variable.scale);
        ji[(2 * i) * n + i] = -1.0;
        ji[(2 * i + 1) * n + i] = 1.0;
    }
    let mut ce = Vec::new();
    let mut je = Vec::new();
    // Box inequalities remain first. Physical rows retain declaration order
    // within each equality/inequality family, including signed feasible slack.
    for (constraint, row) in problem.constraints().iter().zip(&evaluation.constraints) {
        if constraint.sense == ConstraintSense::Equal {
            ce.push(row.residual);
            je.extend_from_slice(&row.gradient);
        } else {
            ci.push(row.residual);
            ji.extend_from_slice(&row.gradient);
        }
    }
    SqpSample { f: evaluation.value, gradient: evaluation.gradient.clone(), ce, ci, je, ji }

}

/// A solve bound to one immutable physical problem and caller-owned work ledger.
///
/// `DesignControl` limits all physical attempts, even failed initialization or
/// failed line searches. It is borrowed exclusively, not copied into a fresh
/// budget on resume. `SqpState` separately counts this study's callback attempts.
/// Neither mutable optimizer state nor a replacement problem can be supplied.
///
/// Domain bounds are actual SQP inequalities. This is a local box-constrained
/// stationary-design search with optional physical response equalities/inequalities,
/// not a global optimum or identifiability certificate.
/// The oracle's fixed modal bases and contact activity-margin restrictions still
/// apply. Dense QP/BFGS phases are bounded but not internally cancellable.
pub struct EquilibriumStudy<'problem, 'work> {
    problem: &'problem EquilibriumDesign,
    control: &'work mut DesignControl,
    state: SqpState,
    accepted: DesignEvaluation,
    last_domain_rejection: Option<DesignError>,
}

impl<'problem, 'work> EquilibriumStudy<'problem, 'work> {
    /// Evaluate a complete initial load-case family once, then initialize SQP.
    /// The dense cap includes decisions, both box faces and all physical rows. An invalid
    /// initial physical point is an error, not an artificial penalty sample.
    /// The caller's work ledger remains charged even if initialization fails.
    pub fn new(
        problem: &'problem EquilibriumDesign,
        point: &[f64],
        control: &'work mut DesignControl,
        maximum_kkt_dimension: usize,
        gate: &CancelGate,
    ) -> Result<Self, EquilibriumStudyError> {
        poll(gate)?;
        admit(problem, maximum_kkt_dimension)?;
        if point.len() != problem.variables().len() {
            return Err(SqpError::Shape { field: "physical design point",
                expected: problem.variables().len(), actual: point.len() });
        }
        let mut initial = None;
        let state = SqpState::try_new(point, maximum_kkt_dimension, &mut |x| {
            let evaluation = problem.evaluate(x, control, gate)?;
            let result = sample(problem, &evaluation);
            initial = Some(evaluation);
            Ok(Some(result))
        }, None).map_err(normalize_error)?;
        poll(gate)?;
        let accepted = initial.ok_or(SqpError::Invalid("missing complete initial physical sample"))?;
        Ok(Self { problem, control, state, accepted, last_domain_rejection: None })
    }

    /// Accepted optimizer point, derivatives, history and cumulative callbacks.
    #[must_use]
    pub fn optimizer(&self) -> &SqpState { &self.state }

    /// Complete independent-case evidence for exactly `optimizer().point()`.
    /// Never the most recent rejected trial, and never a partial family.
    #[must_use]
    pub fn accepted(&self) -> &DesignEvaluation { &self.accepted }

    /// Caller-owned physical work, including failed/rejected attempts and any
    /// work that preceded this study. This is not a FLOP or elapsed-time count.
    #[must_use]
    pub fn work(&self) -> DesignWork { self.control.work() }

    /// Most recent original out-of-bounds trial error. All other physics errors
    /// propagate immediately rather than masquerading as poor objective values.
    #[must_use]
    pub fn last_domain_rejection(&self) -> Option<&DesignError> {
        self.last_domain_rejection.as_ref()
    }

    /// Raise physical allowances without resetting accumulated work. The SQP
    /// callback ceiling supplied to `run` is independent and study-relative.
    pub fn extend_physics_budget(&mut self, maximum_evaluations: usize, maximum_case_solves: usize)
        -> Result<(), DesignError>
    {
        self.control.extend(maximum_evaluations, maximum_case_solves)
    }

    /// Run a bounded number of additional accepted steps. The callback ceiling
    /// includes initialization, prior segments and failed/rejected trials.
    ///
    /// The existing SQP engine owns QP steps, BFGS, merit backtracking and KKT
    /// certification. Out-of-bounds trial values are explicitly unavailable;
    /// they are not clamped, assigned a fabricated gradient, or accepted. This
    /// includes floating-point overshoots at a physical bound. Original primal,
    /// contact-margin, adjoint, cancellation and physical-budget errors propagate.
    ///
    /// A failed search retains the last accepted physical and optimizer states
    /// and spent work. Cancellation is polled before/after each accepted step
    /// and inside the oracle. Resume restarts an interrupted search; splitting
    /// only between accepted steps does not introduce extra evaluations.
    pub fn run(
        &mut self,
        tolerance: f64,
        additional_iterations: usize,
        maximum_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<SqpRunReport, EquilibriumStudyError> {
        poll(gate)?;
        // A zero-step call obtains stop attribution without an oracle callback.
        let mut report = self.advance(tolerance, 0, maximum_evaluations, gate)?;
        for _ in 0..additional_iterations {
            if report.stop != SqpStop::IterationLimit { break; }
            poll(gate)?;
            report = self.advance(tolerance, 1, maximum_evaluations, gate)?;
        }
        poll(gate)?;
        Ok(report)
    }

    fn advance(&mut self, tolerance: f64, steps: usize, maximum_evaluations: usize, gate: &CancelGate)
        -> Result<SqpRunReport, EquilibriumStudyError>
    {
        let problem = self.problem;
        let control = &mut *self.control;
        let last_rejection = &mut self.last_domain_rejection;
        let mut candidate = None;
        let outcome = self.state.try_run(&mut |point| {
            match problem.evaluate(point, control, gate) {
                Ok(evaluation) => {
                    let result = sample(problem, &evaluation);
                    candidate = Some((point.to_vec(), evaluation));
                    Ok(Some(result))
                }
                Err(error @ DesignError::OutsideBounds { .. }) => {
                    *last_rejection = Some(error);
                    Ok(None)
                }
                Err(error) => Err(error),
            }
        }, tolerance, steps, maximum_evaluations, None);
        // At most one accepted step can occur in this call. A later rejected
        // sample must never replace its predecessor's physical evidence.
        if let Some((point, evaluation)) = candidate {
            if point.as_slice() == self.state.point() { self.accepted = evaluation; }
        }
        let report = outcome.map_err(normalize_error)?;
        poll(gate)?;
        Ok(report)
    }
}

#[cfg(test)]
#[path = "equilibrium/constraints_tests.rs"]
mod constraints_tests;

/// Finite-scenario worst-case physical optimization through the same SQP engine.
pub mod scenarios;
