//! Matrix-free Newton--Krylov studies over the live smooth optimization IR.
//! The existing trust-region/Steihaug kernel remains the numerical authority.

use crate::runner::Packing;
use crate::stop::{StopObservation, StopReason, StopRule};
use crate::trust::{TrustRegionProgress, TrustRegionReport, TrustRegionState};
use fs_exec::Cx;
use fs_opt::reverse::{HessianError, ReverseError};
use fs_opt::{Manifold, OptError, ProblemEvaluation, ReverseProblem, ReverseProblemError};
use std::cell::{Cell, RefCell};
use std::sync::Arc;

/// Structural, evaluator or cancellation refusal; accepted state remains usable.
#[derive(Debug, Clone, PartialEq)]
pub enum ReverseNewtonError {
    /// An unconstrained trust region cannot silently drop constraints.
    ConstraintsUnsupported,
    /// Ambient Hessians do not implement a Riemannian trust region.
    NonEuclidean { variable: usize },
    /// Packed input does not bind all variable point coordinates.
    PackedPointLength { expected: usize, actual: usize },
    /// A stop rule has an invalid scalar, window, or empty composite.
    InvalidRule,
    /// The original first- or second-order evaluator refusal.
    Evaluation(ReverseProblemError),
    /// Cancellation was observed before publishing the staged iteration.
    Cancelled,
    /// Attempted-work accounting cannot be represented.
    CounterOverflow,
}
impl From<ReverseProblemError> for ReverseNewtonError {
    fn from(error: ReverseProblemError) -> Self {
        match error {
            ReverseProblemError::Reverse(ReverseError::Evaluation(OptError::Cancelled))
            | ReverseProblemError::Hessian(HessianError::Reverse(ReverseError::Evaluation(OptError::Cancelled))) => Self::Cancelled,
            other => Self::Evaluation(other),
        }
    }
}
impl core::fmt::Display for ReverseNewtonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ConstraintsUnsupported => write!(f, "reverse Newton requires an unconstrained problem"),
            Self::NonEuclidean { variable } => write!(f, "reverse Newton variable {variable} is not Euclidean"),
            Self::PackedPointLength { expected, actual } => write!(f, "reverse Newton needs {expected} coordinates, received {actual}"),
            Self::InvalidRule => write!(f, "invalid reverse Newton stopping rule"),
            Self::Evaluation(error) => write!(f, "{error}"),
            Self::Cancelled => write!(f, "reverse Newton cancelled; accepted checkpoint retained"),
            Self::CounterOverflow => write!(f, "reverse Newton work counter overflow"),
        }
    }
}
impl std::error::Error for ReverseNewtonError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Evaluation(error) => Some(error), _ => None }
    }
}

/// A Hessian-work ceiling is distinct from objective work or convergence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReverseNewtonStop {
    /// User stop rule, objective budget, numerical stall or outer-iteration cap.
    Stopped(StopReason),
    /// The next Hessian product would exceed the caller's cumulative ceiling.
    HessianBudget,
}

/// Accepted solution and actual attempted work, including rolled-back searches.
#[derive(Debug, Clone)]
pub struct ReverseNewtonReport {
    /// Explicit stopping attribution.
    pub stop: ReverseNewtonStop,
    /// Objective evaluations and Hessian products include failed attempts.
    pub solution: TrustRegionReport,
}

fn poll(cx: Option<&Cx<'_>>) -> Result<(), ReverseNewtonError> {
    if let Some(cx) = cx { cx.checkpoint().map_err(|_| ReverseNewtonError::Cancelled)?; }
    Ok(())
}

fn budget(rule: &StopRule) -> Result<usize, ReverseNewtonError> {
    match rule {
        StopRule::GradNorm(t) if t.is_finite() && *t >= 0.0 => Ok(usize::MAX),
        StopRule::ObjectiveBelow(t) if t.is_finite() => Ok(usize::MAX),
        StopRule::Stall { rel, window } if rel.is_finite() && *rel >= 0.0 && *window > 0 && *window < usize::MAX => Ok(usize::MAX),
        StopRule::Budget(cap) => Ok(*cap),
        StopRule::Any(rules) | StopRule::All(rules) if !rules.is_empty() => {
            let mut cap = usize::MAX;
            for rule in rules { cap = cap.min(budget(rule)?); }
            Ok(cap)
        }
        _ => Err(ReverseNewtonError::InvalidRule),
    }
}

fn flatten(blocks: Vec<Vec<f64>>, n: usize) -> Result<Vec<f64>, ReverseProblemError> {
    let mut packed = Vec::new();
    packed.try_reserve_exact(n).map_err(|_| OptError::RuntimeAllocationRefused {
        path: "reverse-newton/packed-derivative", node: None, variable: None,
        elements: n as u64, element_bytes: 8,
    })?;
    for block in blocks { packed.extend(block); }
    Ok(packed)
}

fn sample<'a, 'p>(
    oracle: &'a ReverseProblem<'p>, packing: &Packing, point: &[f64], cx: Option<&Cx<'_>>,
) -> Result<(Arc<ProblemEvaluation<'a, 'p>>, Vec<f64>), ReverseProblemError> {
    let bindings = packing.unpack(point);
    let tape = match cx {
        Some(cx) => oracle.evaluate_cancellable(&bindings, cx)?,
        None => oracle.evaluate(&bindings)?,
    };
    let gradient = match cx {
        Some(cx) => tape.objective_gradient_cancellable(cx)?,
        None => tape.objective_gradient()?,
    };
    Ok((Arc::new(tape), flatten(gradient, packing.dim)?))
}

fn domain_trial(error: &ReverseProblemError) -> bool {
    matches!(error, ReverseProblemError::NonFiniteObjective { .. }
        | ReverseProblemError::Reverse(ReverseError::NonFiniteAdjoint { .. })
        | ReverseProblemError::Reverse(ReverseError::Evaluation(OptError::EvalNonFinite { .. })))
}

/// Resumable unconstrained Euclidean Newton study with exact chain-rule HVPs.
/// The compiled oracle and accepted primal tape are shared across checkpoints;
/// radius, point, gradient, history and actual attempted-work counts are retained.
/// The immutable borrowed oracle prevents changing problem meaning on resume.
///
/// One existing trust-region iteration is staged on a cloned state. A fatal
/// callback error, cancellation or insufficient Hessian allowance discards that
/// staged iteration, not its spent work. No infallible callback placeholder can
/// be published. Completed rejected trials retain their radius contraction.
/// The legacy trust kernel and all other study drivers remain unchanged.
#[derive(Debug, Clone)]
pub struct ReverseNewtonStudy<'a, 'p> {
    oracle: &'a ReverseProblem<'p>,
    packing: Packing,
    tape: Arc<ProblemEvaluation<'a, 'p>>,
    state: TrustRegionState,
    evaluation_limit: usize,
    evals: usize,
    hv_evals: usize,
    rejected: usize,
    last_rejection: Option<ReverseProblemError>,
    numerical_stop: Option<StopReason>,
}

impl<'a, 'p> ReverseNewtonStudy<'a, 'p> {
    /// Bind the problem and evaluate the initial point/gradient exactly once.
    /// The problem's positive/unlimited objective budget includes this sample.
    pub fn new(oracle: &'a ReverseProblem<'p>, point: &[f64], cx: Option<&Cx<'_>>)
        -> Result<Self, ReverseNewtonError>
    {
        poll(cx)?;
        if !oracle.problem().constraints().is_empty() { return Err(ReverseNewtonError::ConstraintsUnsupported); }
        for (variable, declaration) in oracle.problem().vars().iter().enumerate() {
            if !matches!(declaration.manifold, Manifold::Rn { .. }) {
                return Err(ReverseNewtonError::NonEuclidean { variable });
            }
        }
        let packing = Packing::new(oracle.problem());
        if point.len() != packing.dim {
            return Err(ReverseNewtonError::PackedPointLength { expected: packing.dim, actual: point.len() });
        }
        let evaluation_limit = oracle.problem().budget().limit.maximum().map_or(
            usize::MAX, |cap| usize::try_from(cap.get()).unwrap_or(usize::MAX));
        let (tape, gradient) = sample(oracle, &packing, point, cx)?;
        let state = TrustRegionState::new(point, &mut |_| (tape.objective_value(), gradient.clone()));
        poll(cx)?;
        Ok(Self { oracle, packing, tape, state, evaluation_limit, evals: 1, hv_evals: 0,
            rejected: 0, last_rejection: None, numerical_stop: None })
    }

    /// Snapshot without graph evaluation; counts cover rolled-back work too.
    pub fn snapshot(&self) -> TrustRegionReport {
        let mut report = self.state.report();
        report.evals = self.evals; report.hv_evals = self.hv_evals;
        report
    }
    /// Trust radius after the last complete outer iteration.
    pub fn radius(&self) -> f64 { self.state.radius() }
    /// Accepted gradient, in packed Euclidean coordinates.
    pub fn gradient(&self) -> &[f64] { self.state.gradient() }
    /// Initial objective and each completed outer iteration, including rejects.
    pub fn history(&self) -> &[f64] { self.state.history() }
    /// Non-finite objective/gradient trial count (not model-ratio rejections).
    pub fn domain_rejections(&self) -> usize { self.rejected }
    /// Latest original non-finite trial refusal.
    pub fn last_rejection(&self) -> Option<&ReverseProblemError> { self.last_rejection.as_ref() }

    fn report(&self, stop: ReverseNewtonStop) -> ReverseNewtonReport {
        ReverseNewtonReport { stop, solution: self.snapshot() }
    }

    /// Run additional complete trust-region iterations. All Budget leaves and
    /// the problem ceiling limit objective attempts; `max_hessian_products` is
    /// a separate cumulative HVP ceiling, checked before each Krylov/model call.
    /// A final permitted call may finish its iteration. Mid-Krylov exhaustion
    /// preserves the prior radius/iterate; increasing the limit restarts that
    /// solve and may repeat products. No spent work is refunded.
    ///
    /// Finite-domain trial failures contract the radius. Initial or Hessian
    /// failures are typed errors, never fabricated convergence. Cancellation
    /// runs inside AD sweeps and between callbacks. Packing, state cloning and
    /// legacy vector algebra remain whole phases; no wall-clock bound is claimed.
    pub fn run(&mut self, rule: &StopRule, additional_iters: usize,
        max_hessian_products: usize, cx: Option<&Cx<'_>>)
        -> Result<ReverseNewtonReport, ReverseNewtonError>
    {
        let objective_cap = self.evaluation_limit.min(budget(rule)?);
        for completed in 0..=additional_iters {
            poll(cx)?;
            if self.evals >= objective_cap {
                return Ok(self.report(ReverseNewtonStop::Stopped(StopReason::Budget)));
            }
            let current = self.snapshot();
            let observation = StopObservation { grad_norm: current.grad_norm, objective: current.f,
                evals: self.evals, history: self.state.history() };
            if let Some(reason) = rule.check(&observation).or_else(|| self.numerical_stop.clone()) {
                return Ok(self.report(ReverseNewtonStop::Stopped(reason)));
            }
            if completed == additional_iters {
                return Ok(self.report(ReverseNewtonStop::Stopped(StopReason::IterationCap)));
            }
            if self.hv_evals >= max_hessian_products {
                return Ok(self.report(ReverseNewtonStop::HessianBudget));
            }
            let mut candidate = self.state.clone();
            let error = RefCell::new(None);
            let hessian_limit = Cell::new(false);
            let mut trial = None;
            let packing = &self.packing; let oracle = self.oracle; let accepted = &self.tape;
            let evals = &mut self.evals; let hv_evals = &mut self.hv_evals;
            let rejected = &mut self.rejected; let last_rejection = &mut self.last_rejection;
            let progress = {
                let mut fg = |point: &[f64]| {
                    let refusal = || (f64::NAN, vec![f64::NAN; packing.dim]);
                    if error.borrow().is_some() || hessian_limit.get() { return refusal(); }
                    let Some(next) = evals.checked_add(1) else {
                        *error.borrow_mut() = Some(ReverseNewtonError::CounterOverflow); return refusal();
                    };
                    *evals = next;
                    match sample(oracle, packing, point, cx) {
                        Ok((tape, gradient)) => {
                            let f = tape.objective_value(); trial = Some((point.to_vec(), tape));
                            (f, gradient)
                        }
                        Err(source) if domain_trial(&source) => {
                            *rejected += 1; *last_rejection = Some(source); refusal()
                        }
                        Err(source) => { *error.borrow_mut() = Some(source.into()); refusal() }
                    }
                };
                let mut hv = |_: &[f64], direction: &[f64]| {
                    let refusal = || vec![f64::NAN; packing.dim];
                    if error.borrow().is_some() || hessian_limit.get() { return refusal(); }
                    if *hv_evals >= max_hessian_products { hessian_limit.set(true); return refusal(); }
                    *hv_evals += 1;
                    let result = accepted.objective_hessian_vector_product(&packing.unpack(direction), cx)
                        .and_then(|blocks| flatten(blocks, packing.dim));
                    match result {
                        Ok(product) => product,
                        Err(source) => { *error.borrow_mut() = Some(source.into()); refusal() }
                    }
                };
                candidate.run(&mut fg, &mut hv, &StopRule::GradNorm(0.0), 1).progress
            };
            // The infallible kernel only saw temporary sentinels. Refuse the
            // entire staged state before it can expose any placeholder result.
            if let Some(error) = error.into_inner() { return Err(error); }
            if hessian_limit.get() { return Ok(self.report(ReverseNewtonStop::HessianBudget)); }
            poll(cx)?;
            if let Some((point, tape)) = trial {
                if candidate.report().x == point { self.tape = tape; }
            }
            self.numerical_stop = match progress {
                TrustRegionProgress::Stopped(StopReason::IterationCap) => None,
                TrustRegionProgress::Stopped(reason) => Some(reason),
                TrustRegionProgress::Paused => return Err(ReverseNewtonError::Cancelled),
            };
            self.state = candidate;
        }
        unreachable!("inclusive loop returns at its final boundary")
    }
}
