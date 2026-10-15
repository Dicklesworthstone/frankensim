//! Resumable constrained Newton using the shared physical derivative/solver path.
use super::{ReverseNewtonError, ReverseNewtonStudy, ReverseNewtonStop, poll};
use crate::{KktResidual, StopReason, StopRule};
use crate::trust::KrylovOptions;
use fs_exec::Cx;
use fs_opt::{AugmentedLagrangian, ConstraintKind, ReverseProblem, ReverseProblemError};

fn invalid(what: &'static str) -> ReverseNewtonError {
    ReverseProblemError::Augmentation { index: None, what }.into()
}
fn finite(value: f64) -> Result<f64, ReverseNewtonError> {
    if value.is_finite() { Ok(value) } else {
        Err(ReverseProblemError::NonFiniteAugmentation { index: 0, bits: value.to_bits() }.into())
    }
}

/// Fixed penalty-update policy. Constraints and objective are normalized in IR.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AugmentedNewtonOptions {
    initial_penalty: f64,
    penalty_growth: f64,
    maximum_penalty: f64,
    required_reduction: f64,
}
impl Default for AugmentedNewtonOptions {
    fn default() -> Self {
        Self { initial_penalty: 10.0, penalty_growth: 10.0,
            maximum_penalty: 1e12, required_reduction: 0.5 }
    }
}
impl AugmentedNewtonOptions {
    /// Increase penalty when feasibility has not decreased by `reduction`.
    /// Require finite 0 < initial <= maximum, growth > 1, 0 < reduction < 1.
    pub fn new(initial: f64, growth: f64, maximum: f64, reduction: f64)
        -> Result<Self, ReverseNewtonError>
    {
        if ![initial, growth, maximum, reduction].iter().all(|x| x.is_finite())
            || initial <= 0.0 || growth <= 1.0 || maximum < initial
            || reduction <= 0.0 || reduction >= 1.0 {
            return Err(invalid("invalid augmented Newton penalty policy"));
        }
        Ok(Self { initial_penalty: initial, penalty_growth: growth,
            maximum_penalty: maximum, required_reduction: reduction })
    }
}

/// Absolute cumulative work ceilings, not fresh allowances on each run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AugmentedNewtonBudget {
    /// Primal point/gradient attempts, including the initial point and failures.
    pub objective_evals: usize,
    /// Complete augmented Hessian-action attempts, including discarded searches.
    /// Each includes a Lagrangian HVP and matrix-free J^T mu_active J work.
    pub hessian_products: usize,
    /// Attempted multiplier/penalty retargets, precharged before tape reweighting.
    /// These reverse sweeps need no primal solve; failures receive no refund.
    pub multiplier_updates: usize,
}

/// KKT convergence is separate from budget, penalty and inner-solver stops.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AugmentedNewtonStop {
    /// All four numerical KKT residuals for the reported point/multipliers pass.
    Converged,
    /// This call completed its permitted additional Newton iterations.
    IterationCap,
    /// The next objective attempt is not funded.
    ObjectiveBudget,
    /// The next complete augmented Hessian action is not funded.
    HessianBudget,
    /// A multiplier/penalty retarget is not funded.
    MultiplierBudget,
    /// Further required penalty growth exceeds policy. NOT an infeasibility proof.
    PenaltyLimit,
    /// Inner solve stalled before its gradient target; not KKT convergence.
    InnerStopped(StopReason),
}

/// Original objective/residuals and their candidate KKT multipliers.
#[derive(Debug, Clone, PartialEq)]
pub struct AugmentedNewtonReport {
    pub stop: AugmentedNewtonStop,
    pub x: Vec<f64>,
    /// Signed weighted original objective, not the augmented merit.
    pub objective: f64,
    /// Original residuals in constraint declaration order.
    pub constraints: Vec<f64>,
    /// Projected lambda+mu*c at this point. The retained augmented gradient is
    /// exactly the Lagrangian gradient for THESE multipliers, not necessarily
    /// for `frozen_multipliers()`. These values define the reported KKT pair.
    pub multipliers: Vec<f64>,
    /// Numerical stationarity/feasibility/dual/complementarity residuals only;
    /// no derivative-error, global-optimality or continuum certificate.
    pub kkt: KktResidual,
    pub penalty: f64,
    pub newton_iterations: usize,
    pub multiplier_updates: usize,
    pub multiplier_attempts: usize,
    pub objective_evals: usize,
    pub hessian_products: usize,
}

/// PHR outer multiplier iterations with the existing Euclidean Newton kernel.
/// Inner steps may be infeasible. Only the four KKT checks imply `Converged`;
/// an inner gradient stop or exhausted work budget does not. The implementation
/// uses exact IR/physical Hessian actions, not dense SQP matrices or numerical
/// derivatives. Derivatives at PHR shifted-zero kinks use the documented
/// inactive-side generalized action.
///
/// Cloning preserves the current inner state, fixed merit, accepted tape and
/// counters. Retargeting reweights that tape without a primal solve and resets
/// merit-specific radius/history while retaining cumulative work and Krylov
/// policy. Physical-provider caches/budgets remain provider-owned and shared.
/// No process-restart serialization, global memory cap or wall-clock guarantee.
#[derive(Debug, Clone)]
pub struct AugmentedNewtonStudy<'a, 'p> {
    inner: ReverseNewtonStudy<'a, 'p>,
    options: AugmentedNewtonOptions,
    previous_feasibility: f64,
    newton_iterations: usize,
    multiplier_updates: usize,
    multiplier_attempts: usize,
}
impl<'a, 'p> AugmentedNewtonStudy<'a, 'p> {
    /// Initialize with zero multipliers and one physical point evaluation.
    /// The immutable IR problem's own objective ceiling still applies.
    pub fn new(oracle: &'a ReverseProblem<'p>, point: &[f64], options: AugmentedNewtonOptions,
        cx: Option<&Cx<'_>>,
    ) -> Result<Self, ReverseNewtonError> {
        poll(cx)?;
        let multipliers = vec![0.0; oracle.problem().constraints().len()];
        let model = AugmentedLagrangian::new(oracle.problem(), &multipliers, options.initial_penalty)?;
        let inner = ReverseNewtonStudy::new_augmented(oracle, point, model, cx)?;
        Ok(Self { inner, options, previous_feasibility: f64::INFINITY,
            newton_iterations: 0, multiplier_updates: 0, multiplier_attempts: 0 })
    }
    /// Read-only inner checkpoint. `snapshot().f` is the current augmented merit.
    pub fn inner(&self) -> &ReverseNewtonStudy<'a, 'p> { &self.inner }
    /// Frozen multipliers currently used by the inner augmented subproblem.
    pub fn frozen_multipliers(&self) -> &[f64] {
        self.inner.augmentation().expect("augmented constructor").multipliers()
    }
    /// Work already spent on successful or failed merit retargets.
    pub fn multiplier_attempts(&self) -> usize { self.multiplier_attempts }
    /// Completed Newton iterations across merit retargets, also readable after failure.
    pub fn newton_iterations(&self) -> usize { self.newton_iterations }
    /// Configure approximate inner Newton steps; does not weaken KKT criteria.
    pub fn set_krylov_options(&mut self, options: KrylovOptions) { self.inner.set_krylov_options(options); }

    fn report(&self, stop: AugmentedNewtonStop, cx: Option<&Cx<'_>>)
        -> Result<AugmentedNewtonReport, ReverseNewtonError>
    {
        poll(cx)?;
        let model = self.inner.augmentation().expect("augmented constructor");
        let multipliers = model.updated_multipliers(&self.inner.tape, cx)?;
        let mut kkt = KktResidual { stationarity: 0.0, feasibility: 0.0,
            dual_feasibility: 0.0, complementarity: 0.0 };
        for (i, &g) in self.inner.gradient().iter().enumerate() {
            if i % 256 == 0 { poll(cx)?; }
            kkt.stationarity = kkt.stationarity.max(finite(g)?.abs());
        }
        for (i, ((constraint, &c), &lambda)) in self.inner.oracle.problem().constraints().iter()
            .zip(self.inner.constraint_values()).zip(&multipliers).enumerate()
        {
            if i % 256 == 0 { poll(cx)?; }
            finite(c)?;
            match constraint.kind {
                ConstraintKind::EqZero => kkt.feasibility = kkt.feasibility.max(c.abs()),
                ConstraintKind::LeZero => {
                    kkt.feasibility = kkt.feasibility.max(c.max(0.0));
                    kkt.dual_feasibility = kkt.dual_feasibility.max((-lambda).max(0.0));
                    kkt.complementarity = kkt.complementarity.max(finite(lambda * c)?.abs());
                }
            }
        }
        let current = self.inner.snapshot();
        poll(cx)?;
        Ok(AugmentedNewtonReport { stop, x: current.x, objective: self.inner.objective_value(),
            constraints: self.inner.constraint_values().to_vec(), multipliers, kkt,
            penalty: model.penalty(), newton_iterations: self.newton_iterations,
            multiplier_updates: self.multiplier_updates, multiplier_attempts: self.multiplier_attempts,
            objective_evals: current.evals, hessian_products: current.hv_evals })
    }

    /// Complete up to `additional_iterations` Newton iterations across outer
    /// multiplier phases. Ceiling increases resume the retained checkpoint;
    /// no admitted callback work is refunded. `tolerance` is a positive finite
    /// absolute tolerance for ALL four dimensionless KKT residuals. A phase is
    /// retargeted only when its augmented gradient norm is <= tolerance/10.
    /// Zero additional iterations performs no solve or multiplier update.
    ///
    /// A cancelled or failed Newton iteration keeps its prior accepted state.
    /// A cancelled retarget keeps the old merit/tape/radius but retains its
    /// attempted-update charge. Inner numerical stalls are returned explicitly.
    pub fn run(&mut self, tolerance: f64, additional_iterations: usize,
        budget: AugmentedNewtonBudget, cx: Option<&Cx<'_>>,
    ) -> Result<AugmentedNewtonReport, ReverseNewtonError> {
        let inner_tolerance = 0.1 * tolerance;
        if !tolerance.is_finite() || tolerance <= 0.0 || inner_tolerance == 0.0 {
            return Err(invalid("KKT tolerance must be finite, positive and represent tolerance/10"));
        }
        let target = self.newton_iterations.checked_add(additional_iterations)
            .ok_or(ReverseNewtonError::CounterOverflow)?;
        let rule = StopRule::Any(vec![StopRule::GradNorm(inner_tolerance), StopRule::Budget(budget.objective_evals)]);
        loop {
            let mut report = self.report(AugmentedNewtonStop::IterationCap, cx)?;
            if report.kkt.within_tolerance(tolerance) {
                report.stop = AugmentedNewtonStop::Converged; return Ok(report);
            }
            if self.newton_iterations == target { return Ok(report); }
            if self.inner.snapshot().grad_norm <= inner_tolerance {
                if self.multiplier_attempts >= budget.multiplier_updates {
                    report.stop = AugmentedNewtonStop::MultiplierBudget; return Ok(report);
                }
                let mut penalty = report.penalty;
                if report.kkt.feasibility > self.options.required_reduction * self.previous_feasibility {
                    if penalty >= self.options.maximum_penalty {
                        report.stop = AugmentedNewtonStop::PenaltyLimit; return Ok(report);
                    }
                    penalty = (penalty * self.options.penalty_growth).min(self.options.maximum_penalty);
                }
                // Precharge before allocation or reverse reweighting; a failed
                // retarget cannot buy free repeated derivative work on resume.
                self.multiplier_attempts = self.multiplier_attempts.checked_add(1)
                    .ok_or(ReverseNewtonError::CounterOverflow)?;
                let next_updates = self.multiplier_updates.checked_add(1)
                    .ok_or(ReverseNewtonError::CounterOverflow)?;
                let model = AugmentedLagrangian::new(self.inner.oracle.problem(), &report.multipliers, penalty)?;
                self.inner.retarget_augmented(model, cx)?;
                self.previous_feasibility = report.kkt.feasibility;
                self.multiplier_updates = next_updates;
                continue;
            }
            let before = self.inner.snapshot().iters;
            let next_iterations = self.newton_iterations.checked_add(1)
                .ok_or(ReverseNewtonError::CounterOverflow)?;
            let result = self.inner.run(&rule, 1, budget.hessian_products, cx);
            let made_step = self.inner.snapshot().iters > before;
            if made_step { self.newton_iterations = next_iterations; }
            let progress = result?;
            match progress.stop {
                ReverseNewtonStop::HessianBudget => return self.report(AugmentedNewtonStop::HessianBudget, cx),
                ReverseNewtonStop::PreconditionerBudget => unreachable!("outer driver uses unpreconditioned entry point"),
                ReverseNewtonStop::Stopped(StopReason::Budget) => return self.report(AugmentedNewtonStop::ObjectiveBudget, cx),
                ReverseNewtonStop::Stopped(StopReason::GradNorm | StopReason::IterationCap) => {
                    if !made_step && self.inner.snapshot().grad_norm > inner_tolerance {
                        return self.report(AugmentedNewtonStop::InnerStopped(StopReason::IterationCap), cx);
                    }
                }
                ReverseNewtonStop::Stopped(reason) => return self.report(AugmentedNewtonStop::InnerStopped(reason), cx),
            }
        }
    }
}
