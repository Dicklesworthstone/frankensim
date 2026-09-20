//! Linear-storage PHR augmented Lagrangian for a box and ONE smooth inequality.
//!
//! Minimize f(x), l <= x <= u, c(x) <= 0. Bounds are projected exactly, not
//! represented by 2n dense Jacobian rows. A spectral projected-gradient inner
//! step uses monotone Armijo backtracking; multiplier updates happen only after
//! the current box-constrained PHR subproblem meets its inner stationarity gate.
//! Storage is O(n), including trial vectors. No dense Hessian/KKT matrix, inner
//! Krylov solve, numerical derivative, or logistic change of coordinates.
//!
//! PHR is the same inequality model used by `auglag`; this specialized path
//! adds exact box handling and fallible resumable steps, not general sparse SQP.
//! Related method: Gomes-Ruggiero et al., UNICAMP report 2000/41 (ALSPG).
//! This implementation uses a monotone search, not that paper's complete method.
//! Residuals are numerical first-order diagnostics, not global optimality proofs.
use std::ops::ControlFlow;
use crate::auglag::KktResidual;

/// Coherent value/first derivative sample; c <= 0 means feasible.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectedAlSample {
    pub objective: f64,
    pub gradient: Vec<f64>,
    pub constraint: f64,
    pub constraint_gradient: Vec<f64>,
}
impl ProjectedAlSample {
    fn validate<E>(&self, n: usize) -> Result<(), ProjectedAlError<E>> {
        if self.gradient.len() != n || self.constraint_gradient.len() != n {
            return Err(ProjectedAlError::Invalid("sample derivative dimension changed"));
        }
        if !self.objective.is_finite() || !self.constraint.is_finite()
            || self.gradient.iter().chain(&self.constraint_gradient).any(|v| !v.is_finite()) {
            return Err(ProjectedAlError::Invalid("nonfinite sample"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy)]
pub struct ProjectedAlOptions {
    pub max_dimension: usize,
    /// Cumulative callback attempts, including initialization and rejected trials.
    pub max_evaluations: usize,
    pub max_multiplier_updates: usize,
    pub max_backtracks: usize,
    pub tolerance: f64,
    pub initial_penalty: f64,
    pub max_penalty: f64,
    /// Initial inner projected-gradient tolerance, tightened at multiplier updates.
    pub inner_tolerance: f64,
    pub min_spectral_step: f64,
    pub max_spectral_step: f64,
}
impl Default for ProjectedAlOptions {
    fn default() -> Self {
        Self { max_dimension: 1_000_000, max_evaluations: 4096,
            max_multiplier_updates: 100, max_backtracks: 40, tolerance: 1e-6,
            initial_penalty: 10.0, max_penalty: 1e10, inner_tolerance: 1e-2,
            min_spectral_step: 1e-8, max_spectral_step: 1e8 }
    }
}
impl ProjectedAlOptions {
    /// Admission can be run before constructing expensive physical callbacks.
    pub fn validate<E>(&self, n: usize) -> Result<(), ProjectedAlError<E>> {
        if n == 0 || n > self.max_dimension || self.max_evaluations == 0
            || self.max_multiplier_updates == 0 || self.max_backtracks == 0 || self.max_backtracks > 64
            || !self.tolerance.is_finite() || self.tolerance <= 0.0 || self.tolerance >= 1.0
            || !self.initial_penalty.is_finite() || self.initial_penalty <= 0.0
            || !self.max_penalty.is_finite() || self.max_penalty < self.initial_penalty
            || !self.inner_tolerance.is_finite() || self.inner_tolerance < self.tolerance
            || !self.min_spectral_step.is_finite() || self.min_spectral_step <= 0.0
            || !self.max_spectral_step.is_finite() || self.max_spectral_step < self.min_spectral_step {
            return Err(ProjectedAlError::Invalid("invalid projected AL policy or dimension budget"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq)]
pub enum ProjectedAlError<E> {
    Evaluation(E),
    Invalid(&'static str),
    Cancelled,
}
impl<E: std::fmt::Display> std::fmt::Display for ProjectedAlError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Evaluation(e) => write!(f, "projected AL evaluation failed: {e}"),
            Self::Invalid(s) => write!(f, "projected AL refused: {s}"),
            Self::Cancelled => write!(f, "projected AL cancelled"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ProjectedAlError<E> {}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectedAlStop {
    Converged,
    IterationLimit,
    EvaluationLimit,
    MultiplierLimit,
    PenaltyLimit,
    /// No admitted Armijo step; NOT a nonlinear infeasibility certificate.
    Stalled,
}
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ProjectedAlWork {
    pub iterations: usize,
    pub evaluations: usize,
    pub multiplier_updates: usize,
    pub rejected_trials: usize,
}
#[derive(Debug, Clone)]
pub struct ProjectedAlReport {
    pub stop: ProjectedAlStop,
    pub objective: f64,
    pub constraint: f64,
    pub multiplier: f64,
    pub penalty: f64,
    /// Stationarity is ||x-P_box(x-grad L)||_inf, not ||grad L||_inf.
    /// This accounts for bound normals without constructing dense bound rows.
    pub kkt: KktResidual,
    pub work: ProjectedAlWork,
}
fn poll<E>(work: ProjectedAlWork, cp: &mut impl FnMut(ProjectedAlWork) -> ControlFlow<()>)
    -> Result<(), ProjectedAlError<E>> {
    if cp(work).is_break() { Err(ProjectedAlError::Cancelled) } else { Ok(()) }
}
fn finite<E>(v: f64) -> Result<f64, ProjectedAlError<E>> {
    if v.is_finite() { Ok(v) } else { Err(ProjectedAlError::Invalid("nonfinite optimization arithmetic")) }
}
// Project without overflowing x-alpha*g when the answer is a finite endpoint.
fn project(x: f64, g: f64, alpha: f64, lo: f64, hi: f64) -> f64 {
    if g > 0.0 && alpha >= (x-lo)/g { lo }
    else if g < 0.0 && alpha >= (hi-x)/(-g) { hi }
    else { (x-alpha*g).clamp(lo,hi) }
}
fn phr<E>(c: f64, multiplier: f64, penalty: f64) -> Result<(f64,f64), ProjectedAlError<E>> {
    let shifted = finite(multiplier + penalty*c)?;
    let value = if shifted > 0.0 { multiplier*c + (0.5*penalty*c)*c }
        else { -0.5*(multiplier/penalty)*multiplier };
    Ok((finite(value)?, shifted.max(0.0)))
}

/// Accepted checkpoint with bounded O(n) storage; clone is an in-memory restart.
/// Callbacks must keep identical problem meaning across run segments. Errors
/// preserve the last accepted point/sample and spent evaluation counters. Dual
/// updates are complete checkpoint transitions and are not refunded either.
/// Nonlinear inequality feasibility is NOT promised at every accepted inner step.
#[derive(Debug, Clone)]
pub struct ProjectedAlState {
    x: Vec<f64>, lower: Vec<f64>, upper: Vec<f64>, sample: ProjectedAlSample,
    options: ProjectedAlOptions, work: ProjectedAlWork,
    multiplier: f64, penalty: f64, spectral_step: f64,
    inner_tolerance: f64, previous_violation: f64,
}
impl ProjectedAlState {
    pub fn try_new<E>(x: &[f64], lower: &[f64], upper: &[f64], options: ProjectedAlOptions,
        evaluate: &mut impl FnMut(&[f64]) -> Result<Option<ProjectedAlSample>, E>,
        mut checkpoint: impl FnMut(ProjectedAlWork) -> ControlFlow<()>) -> Result<Self, ProjectedAlError<E>> {
        options.validate(x.len())?;
        if lower.len() != x.len() || upper.len() != x.len()
            || x.iter().zip(lower).zip(upper).any(|((&v,&lo),&hi)|
                !v.is_finite() || !lo.is_finite() || !hi.is_finite() || !(hi-lo).is_finite()
                || lo > hi || v < lo || v > hi) {
            return Err(ProjectedAlError::Invalid("finite box and already-in-box start required"));
        }
        let mut work = ProjectedAlWork::default(); poll(work,&mut checkpoint)?;
        work.evaluations = 1;
        let result = evaluate(x).map_err(ProjectedAlError::Evaluation);
        // Count the attempted callback even when it returns an error.
        poll(work,&mut checkpoint)?;
        let sample = result?.ok_or(ProjectedAlError::Invalid("initial sample unavailable"))?;
        sample.validate(x.len())?;
        let state = Self { x: x.to_vec(), lower: lower.to_vec(), upper: upper.to_vec(), sample,
            options, work, multiplier: 0.0, penalty: options.initial_penalty,
            spectral_step: 1.0_f64.clamp(options.min_spectral_step,options.max_spectral_step),
            inner_tolerance: options.inner_tolerance, previous_violation: f64::INFINITY };
        state.report::<E>(ProjectedAlStop::IterationLimit)?;
        Ok(state)
    }
    #[must_use] pub fn point(&self) -> &[f64] { &self.x }
    #[must_use] pub fn sample(&self) -> &ProjectedAlSample { &self.sample }
    #[must_use] pub const fn work(&self) -> ProjectedAlWork { self.work }
    fn gradient<E>(&self, sample: &ProjectedAlSample, multiplier: f64) -> Result<Vec<f64>,ProjectedAlError<E>> {
        sample.gradient.iter().zip(&sample.constraint_gradient).map(|(g,v)| finite(g+multiplier*v)).collect()
    }
    fn projected_norm(&self, g: &[f64]) -> f64 {
        self.x.iter().zip(g).zip(&self.lower).zip(&self.upper)
            .map(|(((&x,&g),&lo),&hi)| (x-project(x,g,1.0,lo,hi)).abs()).fold(0.0,f64::max)
    }
    fn report<E>(&self, stop: ProjectedAlStop) -> Result<ProjectedAlReport,ProjectedAlError<E>> {
        let gradient = self.gradient(&self.sample,self.multiplier)?;
        Ok(ProjectedAlReport { stop, objective: self.sample.objective, constraint: self.sample.constraint,
            multiplier: self.multiplier, penalty: self.penalty, work: self.work,
            kkt: KktResidual { stationarity: self.projected_norm(&gradient),
                feasibility: self.sample.constraint.max(0.0), dual_feasibility: 0.0,
                complementarity: finite((self.multiplier*self.sample.constraint).abs())? } })
    }
    /// Run at most `additional_steps` accepted primal updates. Zero steps report
    /// cached state without evaluating or changing multipliers. Chunking preserves
    /// the spectral step and outer schedule. A cancelled line search restarts
    /// from the same accepted point; repeated probes still consume the budget.
    pub fn try_run<E>(&mut self, additional_steps: usize,
        evaluate: &mut impl FnMut(&[f64]) -> Result<Option<ProjectedAlSample>, E>,
        mut checkpoint: impl FnMut(ProjectedAlWork) -> ControlFlow<()>) -> Result<ProjectedAlReport,ProjectedAlError<E>> {
        let mut completed = 0;
        loop {
            poll(self.work,&mut checkpoint)?;
            let mut report = self.report(ProjectedAlStop::IterationLimit)?;
            if report.kkt.within_tolerance(self.options.tolerance) {
                report.stop = ProjectedAlStop::Converged; return Ok(report);
            }
            if completed == additional_steps { return Ok(report); }
            let (base_term, effective) = phr(self.sample.constraint,self.multiplier,self.penalty)?;
            let gradient = self.gradient(&self.sample,effective)?;
            if self.projected_norm(&gradient) <= self.inner_tolerance {
                if self.work.multiplier_updates == self.options.max_multiplier_updates {
                    return self.report(ProjectedAlStop::MultiplierLimit);
                }
                let violation = self.sample.constraint.max(-self.multiplier/self.penalty).abs();
                let stalled = violation > 0.5*self.previous_violation;
                if stalled && self.penalty == self.options.max_penalty {
                    return self.report(ProjectedAlStop::PenaltyLimit);
                }
                let next_penalty = if stalled { (self.penalty*10.0).min(self.options.max_penalty) } else { self.penalty };
                finite::<E>((effective*self.sample.constraint).abs())?;
                poll(self.work,&mut checkpoint)?;
                self.multiplier = effective;
                self.penalty = next_penalty;
                self.previous_violation = violation;
                self.inner_tolerance = (0.2*self.inner_tolerance).max(0.1*self.options.tolerance);
                self.spectral_step = 1.0_f64.clamp(self.options.min_spectral_step,self.options.max_spectral_step);
                self.work.multiplier_updates += 1;
                continue;
            }
            let projected: Vec<f64> = self.x.iter().zip(&gradient).zip(&self.lower).zip(&self.upper)
                .map(|(((&x,&g),&lo),&hi)| project(x,g,self.spectral_step,lo,hi)).collect();
            let direction: Vec<f64> = projected.iter().zip(&self.x).map(|(p,x)| p-x).collect();
            let slope: f64 = finite(gradient.iter().zip(&direction).map(|(g,d)| g*d).sum())?;
            if slope >= 0.0 { return self.report(ProjectedAlStop::Stalled); }
            let mut alpha = 1.0;
            let mut accepted = None;
            for _ in 0..self.options.max_backtracks {
                poll(self.work,&mut checkpoint)?;
                if self.work.evaluations == self.options.max_evaluations {
                    return self.report(ProjectedAlStop::EvaluationLimit);
                }
                let trial: Vec<f64> = self.x.iter().zip(&projected).zip(&self.lower).zip(&self.upper)
                    .map(|(((&x,&p),&lo),&hi)| if alpha == 1.0 {p} else {(x+alpha*(p-x)).clamp(lo,hi)}).collect();
                if trial == self.x { return self.report(ProjectedAlStop::Stalled); }
                self.work.evaluations += 1;
                let result = evaluate(&trial).map_err(ProjectedAlError::Evaluation);
                poll(self.work,&mut checkpoint)?;
                if let Some(sample) = result? {
                    sample.validate(self.x.len())?;
                    let (term, next_effective) = phr(sample.constraint,self.multiplier,self.penalty)?;
                    let delta = finite((sample.objective-self.sample.objective)+(term-base_term))?;
                    if delta < 0.0 && delta <= 1e-4*alpha*slope {
                        let next_gradient = self.gradient(&sample,next_effective)?;
                        let mut ss = 0.0; let mut sy = 0.0;
                        for (((&x,&next),&g),&ng) in self.x.iter().zip(&trial).zip(&gradient).zip(&next_gradient) {
                            let s = next-x; ss += s*s; sy += s*(ng-g);
                        }
                        let bb = if sy.is_finite() && sy > 0.0 && ss.is_finite() { ss/sy } else {1.0};
                        let bb = if bb.is_finite() { bb } else { self.options.max_spectral_step };
                        // The new sample must also yield a finite ORIGINAL KKT
                        // report before it can replace any accepted state.
                        self.gradient(&sample,self.multiplier)?;
                        finite::<E>((self.multiplier*sample.constraint).abs())?;
                        poll(self.work,&mut checkpoint)?;
                        accepted = Some((trial,sample,bb.clamp(self.options.min_spectral_step,self.options.max_spectral_step)));
                        break;
                    }
                }
                self.work.rejected_trials += 1;
                alpha *= 0.5;
            }
            let Some((x,sample,step)) = accepted else { return self.report(ProjectedAlStop::Stalled); };
            self.x = x; self.sample = sample; self.spectral_step = step;
            self.work.iterations += 1; completed += 1;
        }
    }
}
