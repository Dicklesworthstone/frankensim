//! Powell--Hestenes--Rockafellar scalarization on the same live reverse tape.
use super::{ProblemEvaluation, ReverseProblemError, capacity, poll};
use crate::{ConstraintKind, Problem};
use fs_exec::Cx;
use fs_qty::Dims;

fn invalid(index: Option<usize>, what: &'static str) -> ReverseProblemError {
    ReverseProblemError::Augmentation { index, what }
}
fn finite(index: usize, value: f64) -> Result<f64, ReverseProblemError> {
    if value.is_finite() { Ok(value) }
    else { Err(ReverseProblemError::NonFiniteAugmentation { index, bits: value.to_bits() }) }
}

/// A frozen augmented objective for one multiplier/penalty iteration. Both the
/// signed weighted objective and every constraint must be dimensionless: put
/// physical normalization in the IR, not hidden unit-dependent penalty weights.
/// Equalities retain signed multipliers; inequality multipliers must be >= 0.
///
/// Uses f + lambda*c + mu*c^2/2 for equalities and the PHR expression
/// f + (max(0,lambda+mu*c)^2-lambda^2)/(2*mu) for c <= 0. Its Hessian action
/// includes constraint curvature AND J^T mu_active J. At an inequality's
/// shifted zero it uses the inactive-side generalized Hessian (zero penalty
/// curvature); there is no classical twice-differentiable claim at that kink.
/// No dense Jacobian or parameter Hessian is formed. This value does not run
/// an outer multiplier loop, certify feasibility, or imply global convergence.
#[derive(Debug, Clone)]
pub struct AugmentedLagrangian<'p> {
    problem: &'p Problem,
    multipliers: Vec<f64>,
    penalty: f64,
}
struct Terms { value: f64, weights: Vec<f64>, active: Vec<bool> }

impl<'p> AugmentedLagrangian<'p> {
    /// Validate identity, units, penalty and multiplier signs before any solve.
    pub fn new(problem: &'p Problem, multipliers: &[f64], penalty: f64)
        -> Result<Self, ReverseProblemError>
    {
        if !penalty.is_finite() || penalty <= 0.0 {
            return Err(invalid(None, "penalty must be finite and positive"));
        }
        if multipliers.len() != problem.constraints().len() {
            return Err(ReverseProblemError::SeedCount { expected: problem.constraints().len(), actual: multipliers.len() });
        }
        for objective in problem.objectives() {
            if problem.node_dims(objective.node)? != Dims::NONE {
                return Err(invalid(None, "normalize objectives to dimensionless values in the IR"));
            }
        }
        for (i, (constraint, &weight)) in problem.constraints().iter().zip(multipliers).enumerate() {
            if !weight.is_finite() { return Err(ReverseProblemError::SeedNonFinite { index: i, bits: weight.to_bits() }); }
            if constraint.kind == ConstraintKind::LeZero && weight < 0.0 {
                return Err(invalid(Some(i), "inequality multiplier must be nonnegative"));
            }
            if problem.node_dims(constraint.node)? != Dims::NONE {
                return Err(invalid(Some(i), "normalize constraints to dimensionless residuals in the IR"));
            }
        }
        let mut retained = capacity(multipliers.len(), "augmented/multipliers")?;
        retained.extend_from_slice(multipliers);
        Ok(Self { problem, multipliers: retained, penalty })
    }
    /// Sealed problem to which this multiplier vector belongs.
    pub fn problem(&self) -> &'p Problem { self.problem }
    /// Frozen multipliers in original constraint declaration order.
    pub fn multipliers(&self) -> &[f64] { &self.multipliers }
    /// Positive, fixed penalty during this subproblem.
    pub fn penalty(&self) -> f64 { self.penalty }

    fn terms(&self, tape: &ProblemEvaluation<'_, '_>, cx: Option<&Cx<'_>>) -> Result<Terms, ReverseProblemError> {
        poll(cx, 0)?;
        if !core::ptr::eq(self.problem, tape.owner.problem) {
            return Err(invalid(None, "augmentation and tape belong to different problems"));
        }
        let mut value = tape.objective_value();
        let mut weights = capacity(self.multipliers.len(), "augmented/weights")?;
        let mut active = capacity(self.multipliers.len(), "augmented/active")?;
        for (i, ((constraint, &lambda), &c)) in self.problem.constraints().iter()
            .zip(&self.multipliers).zip(tape.constraint_values()).enumerate()
        {
            poll(cx, i)?;
            let shift = finite(i, lambda + self.penalty * c)?;
            let on = constraint.kind == ConstraintKind::EqZero || shift > 0.0;
            let term = if on {
                // Avoid subtracting two nearly equal squares in the PHR form.
                c * (lambda + (0.5 * self.penalty) * c)
            } else { -(0.5 * lambda / self.penalty) * lambda };
            value = finite(i, value + finite(i, term)?)?;
            weights.push(if on { shift } else { 0.0 });
            active.push(on);
        }
        poll(cx, 0)?;
        Ok(Terms { value, weights, active })
    }
    /// Augmented value and total ambient gradient; one shared reverse sweep.
    pub fn value_gradient(&self, tape: &ProblemEvaluation<'_, '_>, cx: Option<&Cx<'_>>)
        -> Result<(f64, Vec<Vec<f64>>), ReverseProblemError>
    {
        let terms = self.terms(tape, cx)?;
        let gradient = match cx {
            Some(cx) => tape.lagrangian_gradient_cancellable(&terms.weights, cx)?,
            None => tape.lagrangian_gradient(&terms.weights)?,
        };
        Ok((terms.value, gradient))
    }
    /// Complete Hessian action (generalized inactive-side action at a PHR kink).
    /// A directional forward sweep and reverse pullback add J^T mu_active J to
    /// the existing Lagrangian HVP. Bound physical first/second derivatives are
    /// reused under their original contracts, without numerical differencing.
    pub fn hessian_vector_product(&self, tape: &ProblemEvaluation<'_, '_>,
        direction: &[Vec<f64>], cx: Option<&Cx<'_>>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let terms = self.terms(tape, cx)?;
        let mut jv = tape.selected_constraint_direction(&terms.active, direction, cx)?;
        for (i, value) in jv.iter_mut().enumerate() {
            poll(cx, i)?;
            *value = if terms.active[i] { finite(i, self.penalty * *value)? } else { 0.0 };
        }
        let mut result = tape.lagrangian_hessian_vector_product(&terms.weights, direction, cx)?;
        let additional = match cx {
            Some(cx) => tape.constraint_pullback_cancellable(&jv, cx)?,
            None => tape.constraint_pullback(&jv)?,
        };
        for (variable, (block, other)) in result.iter_mut().zip(additional).enumerate() {
            for (i, (value, term)) in block.iter_mut().zip(other).enumerate() {
                poll(cx, i)?;
                *value = finite(variable, *value + term)?;
            }
        }
        poll(cx, 0)?;
        Ok(result)
    }
    /// Projected multiplier update at this tape, without changing the model.
    /// The outer algorithm decides whether its subproblem is sufficiently solved.
    pub fn updated_multipliers(&self, tape: &ProblemEvaluation<'_, '_>, cx: Option<&Cx<'_>>)
        -> Result<Vec<f64>, ReverseProblemError> {
        Ok(self.terms(tape, cx)?.weights)
    }
}

impl ProblemEvaluation<'_, '_> {
    /// Matrix-free J_constraint*v in declaration order. First-order physics
    /// bindings suffice. Objective-only subgraphs are not differentiated.
    pub fn constraint_jacobian_vector_product(&self, direction: &[Vec<f64>], cx: Option<&Cx<'_>>)
        -> Result<Vec<f64>, ReverseProblemError> {
        let mut active = capacity(self.constraint_values().len(), "constraint-direction/active")?;
        for i in 0..self.constraint_values().len() { poll(cx, i)?; active.push(true); }
        self.selected_constraint_direction(&active, direction, cx)
    }
    fn selected_constraint_direction(&self, active: &[bool], direction: &[Vec<f64>], cx: Option<&Cx<'_>>)
        -> Result<Vec<f64>, ReverseProblemError> {
        let mut selected = capacity(self.owner.objective_weights.len(), "constraint-direction/roots")?;
        for i in 0..self.owner.objective_count { poll(cx, i)?; selected.push(false); }
        selected.extend_from_slice(active);
        let values = self.tape.selected_jacobian_vector_product(Some(&selected), direction, cx)?;
        let mut result = capacity(active.len(), "constraint-direction/result")?;
        for (i, &v) in values[self.owner.objective_count..].iter().enumerate() { poll(cx, i)?; result.push(v); }
        Ok(result)
    }
}
