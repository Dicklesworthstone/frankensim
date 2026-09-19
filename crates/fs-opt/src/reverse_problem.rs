//! Weighted objective and constraint-Jacobian products over one shared tape.

use crate::reverse::{HessianError, ReverseError, ReverseEvaluation, ReverseLimits, ReverseProgram};
use crate::{OptError, Problem, Sense};
use fs_exec::Cx;

/// Refusal from the solver-facing problem derivative adapter.
#[derive(Debug, Clone, PartialEq)]
pub enum ReverseProblemError {
    /// Preserve the reverse evaluator's complete refusal.
    Reverse(ReverseError),
    /// Original second-order directional derivative refusal.
    Hessian(HessianError),
    /// A scalarized optimization objective must exist.
    NoObjectives,
    /// The adapter does not execute chance, multi-fidelity or bilevel metadata.
    UnsupportedProblemTags,
    /// Weighted scalarization became non-finite.
    NonFiniteObjective {
        /// Objective-list index at which accumulation failed.
        index: usize,
        /// Exact non-finite accumulated value.
        bits: u64,
    },
    /// Constraint multipliers must match the declared constraint count.
    SeedCount {
        /// Number of declared constraints.
        expected: usize,
        /// Number of supplied multipliers.
        actual: usize,
    },
    /// A non-finite constraint multiplier was supplied.
    SeedNonFinite {
        /// Index in constraint declaration order, not combined root order.
        index: usize,
        /// Exact invalid multiplier bits.
        bits: u64,
    },
    /// Combined root-count arithmetic cannot be represented.
    SizeOverflow,
}

impl From<HessianError> for ReverseProblemError {
    fn from(error: HessianError) -> Self { Self::Hessian(error) }
}

impl From<ReverseError> for ReverseProblemError {
    fn from(error: ReverseError) -> Self {
        Self::Reverse(error)
    }
}

impl From<OptError> for ReverseProblemError {
    fn from(error: OptError) -> Self {
        Self::Reverse(ReverseError::Evaluation(error))
    }
}

impl core::fmt::Display for ReverseProblemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Reverse(error) => write!(f, "{error}"),
            Self::Hessian(error) => write!(f, "{error}"),
            Self::NoObjectives => write!(f, "reverse problem requires at least one objective"),
            Self::UnsupportedProblemTags => write!(
                f, "reverse problem cannot execute chance, multi-fidelity or bilevel tags"
            ),
            Self::NonFiniteObjective { index, bits } => {
                write!(f, "weighted objective {index} is non-finite ({bits:#018x})")
            }
            Self::SeedCount { expected, actual } => {
                write!(f, "expected {expected} constraint multipliers, received {actual}")
            }
            Self::SeedNonFinite { index, bits } => {
                write!(f, "constraint multiplier {index} is non-finite ({bits:#018x})")
            }
            Self::SizeOverflow => write!(f, "reverse problem root count overflow"),
        }
    }
}

impl std::error::Error for ReverseProblemError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Reverse(error) => Some(error),
            Self::Hessian(error) => Some(error),
            _ => None,
        }
    }
}

fn capacity<T>(len: usize, resource: &'static str) -> Result<Vec<T>, ReverseProblemError> {
    let mut out = Vec::new();
    out.try_reserve_exact(len).map_err(|_| OptError::RuntimeAllocationRefused {
        path: resource,
        node: None,
        variable: None,
        elements: len as u64,
        element_bytes: core::mem::size_of::<T>() as u64,
    })?;
    Ok(out)
}

fn poll(cx: Option<&Cx<'_>>, index: usize) -> Result<(), ReverseProblemError> {
    if index % 256 == 0 {
        if let Some(cx) = cx {
            cx.checkpoint().map_err(|_| OptError::Cancelled)?;
        }
    }
    Ok(())
}

fn filled(
    len: usize,
    resource: &'static str,
    cx: Option<&Cx<'_>>,
) -> Result<Vec<f64>, ReverseProblemError> {
    let mut out = capacity(len, resource)?;
    for i in 0..len {
        poll(cx, i)?;
        out.push(0.0);
    }
    Ok(out)
}

/// Solver-facing derivative adapter for a deterministic, untagged IR problem.
///
/// Objective senses and weights come from the sealed Problem. Constraints
/// retain declaration order and their original zero/less-than-zero residuals;
/// this adapter does not turn constraints into penalties or drop inequalities.
/// Metadata-only chance/bilevel/multi-fidelity tags refuse rather than silently
/// changing the optimization problem. Evaluation budgets belong to the caller's
/// solve loop, just as with eval and ReverseProgram; this is not a solve loop.
#[derive(Debug)]
pub struct ReverseProblem<'p> {
    problem: &'p Problem,
    program: ReverseProgram<'p>,
    objective_weights: Vec<f64>,
    objective_count: usize,
}

impl<'p> ReverseProblem<'p> {
    /// Compile all objectives and constraints into one shared reverse program.
    pub fn new(problem: &'p Problem, limits: ReverseLimits) -> Result<Self, ReverseProblemError> {
        Self::compile(problem, limits, None)
    }

    /// Observe cancellation before/after compilation. The underlying reverse
    /// compiler is one non-interruptible phase; this wrapper does not change it.
    pub fn new_cancellable(
        problem: &'p Problem,
        limits: ReverseLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, ReverseProblemError> {
        Self::compile(problem, limits, Some(cx))
    }

    fn compile(
        problem: &'p Problem,
        limits: ReverseLimits,
        cx: Option<&Cx<'_>>,
    ) -> Result<Self, ReverseProblemError> {
        poll(cx, 0)?;
        if problem.objectives().is_empty() {
            return Err(ReverseProblemError::NoObjectives);
        }
        if !problem.tags().is_empty() {
            return Err(ReverseProblemError::UnsupportedProblemTags);
        }
        let count = problem.objectives().len()
            .checked_add(problem.constraints().len())
            .ok_or(ReverseProblemError::SizeOverflow)?;
        if count > limits.max_nodes {
            return Err(OptError::CapExceeded {
                what: "reverse problem roots",
                count: count as u64,
                cap: limits.max_nodes as u64,
            }.into());
        }
        let mut roots = capacity(count, "reverse-problem/roots")?;
        let mut objective_weights = filled(count, "reverse-problem/objective-weights", cx)?;
        for (i, objective) in problem.objectives().iter().enumerate() {
            poll(cx, i)?;
            roots.push(objective.node);
            objective_weights[i] = match objective.sense {
                Sense::Minimize => objective.weight,
                Sense::Maximize => -objective.weight,
            };
        }
        for (i, constraint) in problem.constraints().iter().enumerate() {
            poll(cx, i)?;
            roots.push(constraint.node);
        }
        let program = ReverseProgram::new(problem, &roots, limits)?;
        poll(cx, 0)?;
        Ok(Self {
            problem,
            program,
            objective_weights,
            objective_count: problem.objectives().len(),
        })
    }

    /// The immutable Problem, including constraint kinds in declaration order.
    #[must_use]
    pub fn problem(&self) -> &Problem {
        self.problem
    }

    /// Evaluate all objectives and constraint residuals in one forward sweep.
    pub fn evaluate<'a>(
        &'a self,
        bindings: &[Vec<f64>],
    ) -> Result<ProblemEvaluation<'a, 'p>, ReverseProblemError> {
        self.forward(bindings, None)
    }

    /// Cancellation-aware evaluation and scalarization.
    pub fn evaluate_cancellable<'a>(
        &'a self,
        bindings: &[Vec<f64>],
        cx: &Cx<'_>,
    ) -> Result<ProblemEvaluation<'a, 'p>, ReverseProblemError> {
        self.forward(bindings, Some(cx))
    }

    fn forward<'a>(
        &'a self,
        bindings: &[Vec<f64>],
        cx: Option<&Cx<'_>>,
    ) -> Result<ProblemEvaluation<'a, 'p>, ReverseProblemError> {
        let tape = self.program.evaluate(bindings, cx)?;
        let mut objective = 0.0;
        for i in 0..self.objective_count {
            poll(cx, i)?;
            objective += self.objective_weights[i] * tape.values()[i];
            if !objective.is_finite() {
                return Err(ReverseProblemError::NonFiniteObjective {
                    index: i, bits: objective.to_bits(),
                });
            }
        }
        poll(cx, 0)?;
        Ok(ProblemEvaluation { owner: self, tape, objective })
    }
}

/// A problem valuation and its reusable objective/constraint derivative tape.
/// The tape owns its primal point values; later edits to the input bindings do
/// not change this valuation or its derivatives.
#[derive(Debug)]
pub struct ProblemEvaluation<'a, 'p> {
    owner: &'a ReverseProblem<'p>,
    tape: ReverseEvaluation<'a, 'p>,
    objective: f64,
}

impl ProblemEvaluation<'_, '_> {
    /// Sum of weight*f for minimization and -weight*f for maximization.
    #[must_use]
    pub fn objective_value(&self) -> f64 {
        self.objective
    }

    /// Raw individual objective values, before applying sense or weight.
    #[must_use]
    pub fn objective_values(&self) -> &[f64] {
        &self.tape.values()[..self.owner.objective_count]
    }

    /// Original constraint residuals, in declaration order (no clipping).
    #[must_use]
    pub fn constraint_values(&self) -> &[f64] {
        &self.tape.values()[self.owner.objective_count..]
    }

    /// Ambient gradient of the signed, weighted objective in one reverse sweep.
    pub fn objective_gradient(&self) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        Ok(self.tape.pullback(&self.owner.objective_weights, None)?)
    }

    /// Cancellation-aware objective gradient.
    pub fn objective_gradient_cancellable(
        &self, cx: &Cx<'_>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        Ok(self.tape.pullback(&self.owner.objective_weights, Some(cx))?)
    }

    /// Weighted objective gradient in authoritative manifold parameters.
    pub fn objective_parameter_gradient(&self) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        Ok(self.tape.parameter_pullback(&self.owner.objective_weights, None)?)
    }

    /// Cancellation-aware objective gradient in manifold parameters.
    pub fn objective_parameter_gradient_cancellable(
        &self, cx: &Cx<'_>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        Ok(self.tape.parameter_pullback(&self.owner.objective_weights, Some(cx))?)
    }

    /// Matrix-free constraint Jacobian-transpose product in ambient coordinates.
    /// One multiplier is required per declared constraint, regardless of kind.
    pub fn constraint_pullback(
        &self, multipliers: &[f64],
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, false, None)?;
        Ok(self.tape.pullback(&seeds, None)?)
    }

    /// Cancellation-aware constraint Jacobian-transpose product.
    pub fn constraint_pullback_cancellable(
        &self, multipliers: &[f64], cx: &Cx<'_>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, false, Some(cx))?;
        Ok(self.tape.pullback(&seeds, Some(cx))?)
    }

    /// Constraint Jacobian-transpose product in manifold parameters.
    pub fn constraint_parameter_pullback(
        &self, multipliers: &[f64],
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, false, None)?;
        Ok(self.tape.parameter_pullback(&seeds, None)?)
    }

    /// Cancellation-aware constraint product in manifold parameters.
    pub fn constraint_parameter_pullback_cancellable(
        &self, multipliers: &[f64], cx: &Cx<'_>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, false, Some(cx))?;
        Ok(self.tape.parameter_pullback(&seeds, Some(cx))?)
    }

    /// Gradient of the weighted objective plus multiplier*constraint residuals.
    /// All contributions share ONE reverse sweep, including shared roots.
    /// Multiplier sign feasibility is the constrained solver's responsibility.
    pub fn lagrangian_gradient(
        &self, multipliers: &[f64],
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, true, None)?;
        Ok(self.tape.pullback(&seeds, None)?)
    }

    /// Cancellation-aware Lagrangian gradient in ambient coordinates.
    pub fn lagrangian_gradient_cancellable(
        &self, multipliers: &[f64], cx: &Cx<'_>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, true, Some(cx))?;
        Ok(self.tape.pullback(&seeds, Some(cx))?)
    }

    /// Lagrangian gradient in authoritative manifold parameters.
    pub fn lagrangian_parameter_gradient(
        &self, multipliers: &[f64],
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, true, None)?;
        Ok(self.tape.parameter_pullback(&seeds, None)?)
    }

    /// Cancellation-aware Lagrangian gradient in manifold parameters.
    pub fn lagrangian_parameter_gradient_cancellable(
        &self, multipliers: &[f64], cx: &Cx<'_>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, true, Some(cx))?;
        Ok(self.tape.parameter_pullback(&seeds, Some(cx))?)
    }

    fn seeds(
        &self, multipliers: &[f64], objective: bool, cx: Option<&Cx<'_>>,
    ) -> Result<Vec<f64>, ReverseProblemError> {
        poll(cx, 0)?;
        let expected = self.constraint_values().len();
        if multipliers.len() != expected {
            return Err(ReverseProblemError::SeedCount { expected, actual: multipliers.len() });
        }
        // Validate before allocation. Error indices remain in constraint order.
        for (i, &weight) in multipliers.iter().enumerate() {
            poll(cx, i)?;
            if !weight.is_finite() {
                return Err(ReverseProblemError::SeedNonFinite { index: i, bits: weight.to_bits() });
            }
        }
        let mut seeds = filled(
            self.owner.objective_weights.len(), "reverse-problem/pullback-seeds", cx,
        )?;
        if objective {
            for (i, &weight) in self.owner.objective_weights.iter().enumerate() {
                poll(cx, i)?;
                seeds[i] = weight;
            }
        }
        for (i, &weight) in multipliers.iter().enumerate() {
            poll(cx, i)?;
            seeds[self.owner.objective_count + i] = weight;
        }
        poll(cx, 0)?;
        Ok(seeds)
    }
}

impl ProblemEvaluation<'_, '_> {
    /// Hessian of the signed weighted objective applied to an ambient direction.
    /// Reuses the accepted primal tape; no new objective evaluation, coordinate
    /// perturbations, or dense Hessian. `None` explicitly disables cancellation.
    /// Direction blocks use point storage, not manifold parameters. This is NOT
    /// a Riemannian Hessian and has no interval/error-certificate interpretation.
    pub fn objective_hessian_vector_product(
        &self, direction: &[Vec<f64>], cx: Option<&Cx<'_>>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        Ok(self.tape.hessian_vector_product(&self.owner.objective_weights, direction, cx)?)
    }

    /// Hessian action of objective + sum(multiplier * constraint), holding all
    /// multipliers fixed. Constraints keep declaration order and raw residuals.
    /// Primal values are shared with the objective and constraint derivatives.
    pub fn lagrangian_hessian_vector_product(
        &self, multipliers: &[f64], direction: &[Vec<f64>], cx: Option<&Cx<'_>>,
    ) -> Result<Vec<Vec<f64>>, ReverseProblemError> {
        let seeds = self.seeds(multipliers, true, cx)?;
        Ok(self.tape.hessian_vector_product(&seeds, direction, cx)?)
    }
}
