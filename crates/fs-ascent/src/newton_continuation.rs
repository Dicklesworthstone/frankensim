//! Warm-started continuation through immutable live-IR Newton problems.
//!
//! Each stage keeps its own physics bindings and declared model. Only accepted
//! coordinates cross a stage boundary; tapes, multipliers, merit history and
//! trust radii do not. The existing augmented Newton solver owns every step.

use crate::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions,
    AugmentedNewtonReport, AugmentedNewtonStop, AugmentedNewtonStudy, ReverseNewtonError};
use crate::trust::KrylovOptions;
use fs_exec::Cx;
use fs_opt::{AugmentedLagrangian, Manifold, ReverseProblem};

/// One immutable continuation model, with an explicit KKT tolerance and policy.
/// All stages must use the same variable names, dimensions and packed layout.
/// The caller must also keep their coordinate meaning fixed (for example, use
/// the SAME LogisticBox); that fact cannot be inferred from a provider callback.
#[derive(Debug, Clone, Copy)]
pub struct NewtonStage<'a, 'p> {
    pub oracle: &'a ReverseProblem<'p>,
    pub tolerance: f64,
    pub options: AugmentedNewtonOptions,
    pub krylov: KrylovOptions,
}

/// Work across every stage, including failed entry attempts and discarded steps.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ContinuationWork {
    /// Includes one conservative objective reservation per stage-entry attempt.
    pub objective_evals: usize,
    pub hessian_products: usize,
    pub multiplier_attempts: usize,
    pub newton_iterations: usize,
    pub stage_initializations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ContinuationError {
    InvalidInput(&'static str),
    /// Zero-based stage at which the original failure occurred. On a failed
    /// entry this may be one beyond active_stage(); the old checkpoint remains.
    Stage { stage: usize, source: ReverseNewtonError },
    CounterOverflow,
}
impl std::fmt::Display for ContinuationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput(s) => write!(f, "invalid Newton continuation: {s}"),
            Self::Stage { stage, source } => write!(f, "Newton continuation stage {stage}: {source}"),
            Self::CounterOverflow => write!(f, "Newton continuation counter overflow"),
        }
    }
}
impl std::error::Error for ContinuationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Stage { source, .. } => Some(source), _ => None }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContinuationStop {
    /// EVERY stage met its own numerical KKT tolerance. Not global optimality.
    Complete,
    /// This call used its additional Newton-iteration allowance.
    IterationCap,
    /// Entering the first/next model requires an unfunded objective reservation.
    EntryBudget,
    /// Preserve the precise optimizer stop; never silently skip a stalled stage.
    StageStopped(AugmentedNewtonStop),
}

#[derive(Debug, Clone, PartialEq)]
pub struct ContinuationReport {
    pub stop: ContinuationStop,
    /// Zero-based current stage (zero while the first stage is uninitialized).
    pub stage: usize,
    /// Current model's original objective, constraints and KKT diagnostics.
    /// None means no stage has successfully initialized; no value is fabricated.
    pub solution: Option<AugmentedNewtonReport>,
    pub work: ContinuationWork,
}

/// In-process resumable continuation, with absolute optimizer-work ceilings.
///
/// Construction performs admission, not physics. run(0, ...) neither evaluates
/// a point nor enters a new stage. A stage advances only after numerical KKT
/// convergence. A failed or cancelled entry retains the prior accepted study;
/// the attempted entry is still charged. Increasing the budget retries it.
///
/// Physics providers own their separate native work limits. Share a provider
/// resource handle across stages (e.g. SimpNewtonResources); this generic driver
/// cannot count hidden PDE solves or silently reset provider budgets. Cloning
/// preserves optimizer continuation state, but does not clone provider caches.
/// No self-referential ownership, new linear solver, process serialization or
/// whole-process memory bound is introduced.
#[derive(Debug, Clone)]
pub struct NewtonContinuation<'a, 'p> {
    stages: Vec<NewtonStage<'a, 'p>>,
    initial: Vec<f64>,
    index: usize,
    active: Option<AugmentedNewtonStudy<'a, 'p>>,
    completed: Vec<AugmentedNewtonReport>,
    // Retired work excludes each stage's initial evaluation, which is charged
    // through entries BEFORE construction. Active counters are read live, even
    // after an error/unwind, rather than reconciled only after a successful run.
    retired: ContinuationWork,
    entries: usize,
}

fn poll(stage: usize, cx: Option<&Cx<'_>>) -> Result<(), ContinuationError> {
    if let Some(cx) = cx {
        cx.checkpoint().map_err(|_| ContinuationError::Stage {
            stage, source: ReverseNewtonError::Cancelled,
        })?;
    }
    Ok(())
}
fn local_work(study: &AugmentedNewtonStudy<'_, '_>) -> ContinuationWork {
    let current = study.inner().snapshot();
    ContinuationWork { objective_evals: current.evals, hessian_products: current.hv_evals,
        multiplier_attempts: study.multiplier_attempts(), newton_iterations: study.newton_iterations(),
        stage_initializations: 0 }
}

impl<'a, 'p> NewtonContinuation<'a, 'p> {
    pub fn new(stages: &[NewtonStage<'a, 'p>], point: &[f64], cx: Option<&Cx<'_>>)
        -> Result<Self, ContinuationError>
    {
        let invalid = ContinuationError::InvalidInput;
        let first = stages.first().ok_or(invalid("at least one stage is required"))?;
        let variables = first.oracle.problem().variables();
        let mut dimension = 0_usize;
        for (i, v) in variables.iter().enumerate() {
            if i % 256 == 0 { poll(0, cx)?; }
            let Manifold::Rn { dim } = &v.manifold else { return Err(invalid("Euclidean variables required")); };
            dimension = dimension.checked_add(*dim as usize).ok_or(ContinuationError::CounterOverflow)?;
        }
        if dimension == 0 || point.len() != dimension || point.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite initial point must match the nonempty packed layout"));
        }
        for (i, stage) in stages.iter().enumerate() {
            poll(i, cx)?;
            if stage.oracle.problem().variables() != variables {
                return Err(invalid("stage variable names, dimensions and layouts must agree"));
            }
            if !stage.tolerance.is_finite() || stage.tolerance <= 0.0 || 0.1 * stage.tolerance == 0.0 {
                return Err(invalid("each stage needs a finite positive resolved KKT tolerance"));
            }
            // Admit dimensionless objective/constraints now, before any stage
            // can spend physical work. This constructs no evaluation tape.
            let zeros = vec![0.0; stage.oracle.problem().constraints().len()];
            AugmentedLagrangian::new(stage.oracle.problem(), &zeros, 1.0)
                .map_err(|e| ContinuationError::Stage { stage: i, source: e.into() })?;
        }
        Ok(Self { stages: stages.to_vec(), initial: point.to_vec(), index: 0,
            active: None, completed: Vec::with_capacity(stages.len().saturating_sub(1)),
            retired: ContinuationWork::default(), entries: 0 })
    }

    pub fn active_stage(&self) -> usize { self.index }
    pub fn active(&self) -> Option<&AugmentedNewtonStudy<'a, 'p>> { self.active.as_ref() }
    /// Completed models before the current stage, in declared order. Their
    /// objectives are model-specific and must not be read as one descent trace.
    pub fn completed_stages(&self) -> &[AugmentedNewtonReport] { &self.completed }

    pub fn work(&self) -> ContinuationWork {
        let mut work = self.retired;
        work.objective_evals += self.entries;
        work.stage_initializations = self.entries;
        if let Some(active) = &self.active {
            let live = local_work(active);
            work.objective_evals += live.objective_evals - 1;
            work.hessian_products += live.hessian_products;
            work.multiplier_attempts += live.multiplier_attempts;
            work.newton_iterations += live.newton_iterations;
        }
        work
    }

    fn local_budget(&self, budget: AugmentedNewtonBudget) -> AugmentedNewtonBudget {
        let total = self.work();
        let local = self.active.as_ref().map(local_work).unwrap_or_default();
        // Subtraction, rather than local + remaining, remains representable at
        // usize::MAX. Retired attempts cannot buy another stage a fresh budget.
        AugmentedNewtonBudget {
            objective_evals: budget.objective_evals.saturating_sub(total.objective_evals - local.objective_evals),
            hessian_products: budget.hessian_products.saturating_sub(total.hessian_products - local.hessian_products),
            multiplier_updates: budget.multiplier_updates.saturating_sub(total.multiplier_attempts - local.multiplier_attempts),
        }
    }

    fn report(&self, stop: ContinuationStop, solution: Option<AugmentedNewtonReport>) -> ContinuationReport {
        ContinuationReport { stop, stage: self.index, solution, work: self.work() }
    }

    fn enter(&mut self, next: usize, previous: Option<AugmentedNewtonReport>, cx: Option<&Cx<'_>>)
        -> Result<(), ContinuationError>
    {
        poll(next, cx)?;
        let point = previous.as_ref().map_or(self.initial.as_slice(), |r| r.x.as_slice());
        self.entries = self.entries.checked_add(1).ok_or(ContinuationError::CounterOverflow)?;
        let stage = self.stages[next];
        let mut candidate = AugmentedNewtonStudy::new(stage.oracle, point, stage.options, cx)
            .map_err(|source| ContinuationError::Stage { stage: next, source })?;
        candidate.set_krylov_options(stage.krylov);
        poll(next, cx)?;
        if let Some(old) = &self.active {
            let work = local_work(old);
            self.retired.objective_evals += work.objective_evals - 1;
            self.retired.hessian_products += work.hessian_products;
            self.retired.multiplier_attempts += work.multiplier_attempts;
            self.retired.newton_iterations += work.newton_iterations;
        }
        if let Some(previous) = previous { self.completed.push(previous); }
        self.active = Some(candidate);
        self.index = next;
        Ok(())
    }

    /// Run up to additional_iterations accepted/rejected Newton outer steps
    /// across the schedule, using ONE absolute budget for all stage work.
    /// A later problem's independent IR budget still applies locally. Native
    /// provider limits remain separate. No stage is skipped because it stalls,
    /// exhausts a resource, or fails preparation. All prior successful stages
    /// remain inspectable, and the last accepted study remains available on Err.
    pub fn run(&mut self, additional_iterations: usize, budget: AugmentedNewtonBudget, cx: Option<&Cx<'_>>)
        -> Result<ContinuationReport, ContinuationError>
    {
        let target = self.work().newton_iterations.checked_add(additional_iterations)
            .ok_or(ContinuationError::CounterOverflow)?;
        loop {
            poll(self.index, cx)?;
            if self.active.is_none() {
                if additional_iterations == 0 { return Ok(self.report(ContinuationStop::IterationCap, None)); }
                if self.work().objective_evals >= budget.objective_evals {
                    return Ok(self.report(ContinuationStop::EntryBudget, None));
                }
                self.enter(0, None, cx)?;
            }
            let stage = self.stages[self.index];
            let local_budget = self.local_budget(budget);
            let current = self.active.as_mut().expect("initialized stage")
                .run(stage.tolerance, 0, local_budget, cx)
                .map_err(|source| ContinuationError::Stage { stage: self.index, source })?;
            if current.stop == AugmentedNewtonStop::Converged && self.index + 1 == self.stages.len() {
                return Ok(self.report(ContinuationStop::Complete, Some(current)));
            }
            if self.work().newton_iterations == target {
                return Ok(self.report(ContinuationStop::IterationCap, Some(current)));
            }
            if current.stop == AugmentedNewtonStop::Converged {
                if self.work().objective_evals >= budget.objective_evals {
                    return Ok(self.report(ContinuationStop::EntryBudget, Some(current)));
                }
                self.enter(self.index + 1, Some(current), cx)?;
                continue;
            }
            let progress = self.active.as_mut().expect("initialized stage")
                .run(stage.tolerance, 1, local_budget, cx)
                .map_err(|source| ContinuationError::Stage { stage: self.index, source })?;
            match &progress.stop {
                AugmentedNewtonStop::Converged | AugmentedNewtonStop::IterationCap => {}
                stop => return Ok(self.report(ContinuationStop::StageStopped(stop.clone()), Some(progress))),
            }
        }
    }
}
