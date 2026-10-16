//! Small-dimension, volume/bound-constrained response design using fs-ascent SQP.
//! No optimizer, KKT solver, BFGS update or merit search is reimplemented here.
use super::*;
use fs_ascent::sqp::{SqpError, SqpRunReport, SqpSample, SqpState, SqpStop};

#[derive(Debug, Clone, Copy)]
pub struct ResponseDesignOptions3 {
    pub response: ResponseOptions3,
    /// Numerical projected-volume inequality V(rho)-volume_cap <= 0.
    pub volume_cap: f64,
    /// Raw density bounds are [density_floor, 1].
    pub density_floor: f64,
    /// Numerical KKT gate in the declared objective/constraint coordinates.
    pub tolerance: f64,
    /// Cumulative SQP sample attempts, including unavailable/rejected trials.
    pub max_evaluations: usize,
    /// Decisions plus ALL constraints; admitted before any physical evaluation.
    /// Hard ceiling 1024 keeps this explicitly in the existing dense SQP regime.
    pub max_kkt_dimension: usize,
}
impl Default for ResponseDesignOptions3 {
    fn default() -> Self {
        Self { response: ResponseOptions3::default(), volume_cap: 0.5, density_floor: 1e-3,
            tolerance: 1e-7, max_evaluations: 256, max_kkt_dimension: 256 }
    }
}
/// One accepted SQP point. Acceptance is merit-based and is NOT an assertion
/// that every intermediate point meets the nonlinear material-volume cap.
#[derive(Debug, Clone)]
pub struct ResponseDesignIteration3 {
    pub iteration: usize,
    pub objective: f64,
    pub volume_fraction: f64,
    pub constraint_violation: f64,
}
fn violation(sample: &SqpSample) -> f64 { sample.ci.iter().map(|v| v.max(0.0)).fold(0.0, f64::max) }
fn row(e: &ResponseEvaluation3, options: ResponseDesignOptions3) -> SqpSample {
    let n = e.rho.len();
    let mut ci = Vec::with_capacity(1+2*n); ci.push(e.volume_fraction-options.volume_cap);
    ci.extend(e.rho.iter().map(|r| options.density_floor-r));
    ci.extend(e.rho.iter().map(|r| r-1.0));
    let mut ji = vec![0.0; (1+2*n)*n]; ji[..n].copy_from_slice(&e.volume_gradient);
    for i in 0..n { ji[(1+i)*n+i] = -1.0; ji[(1+n+i)*n+i] = 1.0; }
    SqpSample { f: e.objective, gradient: e.gradient.clone(), ce: Vec::new(), je: Vec::new(), ci, ji }
}
fn evaluate<O: AdaptiveSdf3Elasticity>(study: &mut CutDensityStudy3<O>, cases: &[ResponseCase3<'_>],
    options: ResponseDesignOptions3, rho: &[f64], last: &mut Option<ResponseEvaluation3>,
    control: &mut SolveControl<'_>) -> Result<Option<SqpSample>, ResponseError3> {
    control.checkpoint("response-sqp-evaluate")?;
    // Only physically unavailable density-box probes become line-search
    // rejections. Original physics, preconditioner and budget failures propagate.
    if rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) { return Ok(None); }
    let evaluation = study.evaluate_responses(rho, cases, options.response, control)?;
    let sample = row(&evaluation, options);
    *last = Some(evaluation);
    Ok(Some(sample))
}

/// Resumable optimizer bound to one geometry, experiment family and work ledger.
/// Targets, motions and numerical options cannot change under cached BFGS/KKT
/// state. Only read-only physical access is exposed while the session is live.
/// Trial evaluations restore scales; ONLY a completely accepted SQP point is
/// installed. An interrupted line search can repeat work on resume, never
/// refunds it, and never replaces accepted fields with a rejected candidate.
///
/// Nonlinear constraints may be violated at intermediate accepted points;
/// inspect `constraint_violation()` and the returned numerical KKT residuals.
/// This is small dense constrained design, not large-scale sparse SQP, a global
/// optimum, physical material identification, or continuum certification.
pub struct ResponseDesignStudy3<'a, 'callback, O: AdaptiveSdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    cases: &'a [ResponseCase3<'a>],
    control: &'a mut SolveControl<'callback>,
    options: ResponseDesignOptions3,
    state: SqpState,
    accepted: ResponseEvaluation3,
    history: Vec<ResponseDesignIteration3>,
}
impl<'a, 'callback, O: AdaptiveSdf3Elasticity> ResponseDesignStudy3<'a, 'callback, O> {
    pub fn new(study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseCase3<'a>], rho0: &[f64],
        options: ResponseDesignOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, SqpError<ResponseError3>> {
        let n = study.cells();
        let dim = n.checked_mul(3).and_then(|v| v.checked_add(1));
        if n == 0 || rho0.len() != n || options.max_kkt_dimension > 1024
            || dim.is_none_or(|v| v > options.max_kkt_dimension) || options.max_evaluations == 0
            || !options.volume_cap.is_finite() || options.volume_cap <= 0.0 || options.volume_cap > 1.0
            || !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
            || !options.tolerance.is_finite() || options.tolerance <= 0.0 || options.tolerance >= 1.0 {
            return Err(SqpError::Invalid("invalid response SQP options or dense KKT dimension"));
        }
        admit(study, cases, options.response).map_err(SqpError::Evaluation)?;
        let mut last = None;
        let state = SqpState::try_new(rho0, options.max_kkt_dimension,
            &mut |rho| evaluate(study, cases, options, rho, &mut last, control), None)?;
        control.checkpoint("response-sqp-initialize").map_err(|e| SqpError::Evaluation(e.into()))?;
        let accepted = last.expect("initial SQP sample has complete response evidence");
        study.operator.set_scales(&accepted.scales).expect("evaluated scales are admitted");
        let history = vec![ResponseDesignIteration3 { iteration: 0, objective: accepted.objective,
            volume_fraction: accepted.volume_fraction, constraint_violation: violation(state.sample()) }];
        Ok(Self { study, cases, control, options, state, accepted, history })
    }
    #[must_use] pub fn accepted(&self) -> &ResponseEvaluation3 { &self.accepted }
    /// Physical model with scales matching the accepted fields, never last trial.
    #[must_use] pub fn study(&self) -> &CutDensityStudy3<O> { self.study }
    #[must_use] pub fn point(&self) -> &[f64] { self.state.point() }
    #[must_use] pub fn sample(&self) -> &SqpSample { self.state.sample() }
    #[must_use] pub fn iterations(&self) -> usize { self.state.iterations() }
    #[must_use] pub fn evaluations(&self) -> usize { self.state.evaluations() }
    #[must_use] pub fn history(&self) -> &[ResponseDesignIteration3] { &self.history }
    #[must_use] pub fn work(&self) -> SolveWork { self.control.work() }
    #[must_use] pub fn constraint_violation(&self) -> f64 { violation(self.state.sample()) }

    /// Continue the SAME SQP state for at most this many further accepted steps.
    /// Dense QP/BFGS phases cannot be interrupted internally. Poll before each
    /// bounded step and all physics/adjoint/setup work; synchronize accepted
    /// evidence before a post-step cancellation is returned. A zero-step call
    /// reports cached numerical state without repeating any physical solve.
    pub fn run(&mut self, additional_steps: usize) -> Result<SqpRunReport, SqpError<ResponseError3>> {
        let mut remaining = additional_steps;
        loop {
            self.control.checkpoint("response-sqp-step").map_err(|e| SqpError::Evaluation(e.into()))?;
            let before = self.state.iterations();
            let mut last = None;
            let outcome = {
                let study = &mut *self.study; let control = &mut *self.control;
                let cases = self.cases; let options = self.options;
                self.state.try_run(&mut |rho| evaluate(study, cases, options, rho, &mut last, control),
                    options.tolerance, remaining.min(1), options.max_evaluations, None)
            };
            if self.state.iterations() > before {
                let accepted = last.expect("accepted SQP step has an evaluated response");
                assert_eq!(accepted.rho, self.state.point(), "accepted physics and SQP point must match");
                self.study.operator.set_scales(&accepted.scales).expect("accepted scales are admitted");
                self.accepted = accepted;
                self.history.push(ResponseDesignIteration3 { iteration: self.state.iterations(), objective: self.accepted.objective,
                    volume_fraction: self.accepted.volume_fraction, constraint_violation: violation(self.state.sample()) });
                self.control.checkpoint("response-sqp-accepted").map_err(|e| SqpError::Evaluation(e.into()))?;
            }
            let report = outcome?;
            if report.stop != SqpStop::IterationLimit || remaining <= 1 { return Ok(report); }
            remaining -= 1;
        }
    }
}
