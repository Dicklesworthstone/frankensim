//! Response fitting without a dense KKT system or one bound row per density.
//! The physical evaluator, prescribed-motion adjoints and filter stay unchanged.
use super::*;
use std::cell::RefCell;
use fs_ascent::projected_al::{ProjectedAlError,ProjectedAlOptions,ProjectedAlReport,ProjectedAlSample,ProjectedAlState,ProjectedAlStop,ProjectedAlWork};

#[derive(Debug, Clone, Copy)]
pub struct ProjectedResponseOptions3 {
    pub response: ResponseOptions3,
    pub volume_cap: f64,
    pub density_floor: f64,
    /// Optimizer f and gradient are divided by this declared positive scale.
    /// Physical ResponseEvaluation3 objectives and responses are never rescaled.
    pub objective_scale: f64,
    pub optimizer: ProjectedAlOptions,
}
impl Default for ProjectedResponseOptions3 {
    fn default() -> Self {
        Self { response: ResponseOptions3::default(), volume_cap: 0.5,
            density_floor: 1e-3, objective_scale: 1.0, optimizer: ProjectedAlOptions::default() }
    }
}
#[derive(Debug, Clone)]
pub struct ProjectedResponseIteration3 {
    pub iteration: usize,
    pub objective: f64,
    pub volume_fraction: f64,
    /// Accepted PHR steps can be infeasible for the nonlinear volume constraint.
    pub constraint_violation: f64,
}
fn sample<O: AdaptiveSdf3Elasticity>(study: &mut CutDensityStudy3<O>, cases: &[ResponseCase3<'_>],
    rho: &[f64], options: ProjectedResponseOptions3, last: &mut Option<ResponseEvaluation3>,
    control: &mut SolveControl<'_>) -> Result<Option<ProjectedAlSample>,ResponseError3> {
    control.checkpoint("response-projected-evaluate")?;
    let result = study.evaluate_responses(rho,cases,options.response,control)?;
    let sample = ProjectedAlSample { objective: result.objective/options.objective_scale,
        gradient: result.gradient.iter().map(|g|g/options.objective_scale).collect(),
        constraint: result.volume_fraction-options.volume_cap, constraint_gradient: result.volume_gradient.clone() };
    *last = Some(result);
    Ok(Some(sample))
}
fn row(e: &ResponseEvaluation3, iteration: usize, cap: f64) -> ProjectedResponseIteration3 {
    ProjectedResponseIteration3 { iteration, objective: e.objective,
        volume_fraction: e.volume_fraction, constraint_violation: (e.volume_fraction-cap).max(0.0) }
}

/// Resumable large-design path for a single material-volume inequality and
/// exact raw-density bounds. Optimization state uses O(n) vectors; scalar
/// history is bounded by the cumulative callback allowance. Physics fields and
/// preconditioner memory retain their own existing geometry/solve limits.
///
/// Every trial reuses evaluate_responses, including CURRENT-density lifting,
/// full Nitsche/ghost derivatives and shared primal/adjoint preparation. Only
/// completely accepted optimizer points install physical scales. An error keeps
/// accepted fields and spent work; there is no finite-difference or fixed-load
/// shortcut. PHR descent need not decrease the original objective, and nonlinear
/// volume feasibility is assessed separately. This is not general sparse SQP.
pub struct ProjectedResponseStudy3<'a,'callback,O: AdaptiveSdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseCase3<'a>],
    control: &'a mut SolveControl<'callback>, options: ProjectedResponseOptions3,
    state: ProjectedAlState, accepted: ResponseEvaluation3, history: Vec<ProjectedResponseIteration3>,
}
impl<'a,'callback,O: AdaptiveSdf3Elasticity> ProjectedResponseStudy3<'a,'callback,O> {
    pub fn new(study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseCase3<'a>], rho: &[f64],
        options: ProjectedResponseOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self,ProjectedAlError<ResponseError3>> {
        let n = study.cells(); options.optimizer.validate(n)?;
        if rho.len()!=n || !options.volume_cap.is_finite() || options.volume_cap<=0.0 || options.volume_cap>1.0
            || !options.density_floor.is_finite() || options.density_floor<=0.0 || options.density_floor>=1.0
            || !options.objective_scale.is_finite() || options.objective_scale<=0.0
            || rho.iter().any(|r|!r.is_finite() || *r<options.density_floor || *r>1.0) {
            return Err(ProjectedAlError::Invalid("invalid projected response policy or starting design"));
        }
        admit(study,cases,options.response).map_err(ProjectedAlError::Evaluation)?;
        let mut last = None;
        let state = {
            // Callbacks are synchronous and never overlap. Each RefCell borrow
            // ends before the optimizer invokes its next poll or evaluation.
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(rho,&vec![options.density_floor;n],&vec![1.0;n],options.optimizer,
                &mut |x| sample(study,cases,x,options,&mut last,&mut **ledger.borrow_mut()),
                |_| if ledger.borrow_mut().checkpoint("response-projected-control").is_ok() {ControlFlow::Continue(())} else {ControlFlow::Break(())})?
        };
        control.checkpoint("response-projected-initialize").map_err(|e|ProjectedAlError::Evaluation(e.into()))?;
        let accepted = last.expect("accepted initial sample has complete physical evidence");
        study.operator.set_scales(&accepted.scales).expect("evaluated scales are admitted");
        let history = vec![row(&accepted,0,options.volume_cap)];
        Ok(Self { study,cases,control,options,state,accepted,history })
    }
    #[must_use] pub fn accepted(&self) -> &ResponseEvaluation3 { &self.accepted }
    #[must_use] pub fn study(&self) -> &CutDensityStudy3<O> { self.study }
    #[must_use] pub fn point(&self) -> &[f64] { self.state.point() }
    #[must_use] pub fn evaluations(&self) -> usize { self.state.work().evaluations }
    #[must_use] pub fn history(&self) -> &[ProjectedResponseIteration3] { &self.history }
    #[must_use] pub fn optimizer_work(&self) -> ProjectedAlWork { self.state.work() }
    #[must_use] pub fn work(&self) -> SolveWork { self.control.work() }
    #[must_use] pub fn constraint_violation(&self) -> f64 { self.state.sample().constraint.max(0.0) }

    /// Continue the SAME multiplier/spectral state. Synchronize a committed
    /// physical step even if a later poll reports cancellation. Zero steps use
    /// cached derivatives and do not re-solve. No history or work is refunded.
    pub fn run(&mut self, additional_steps: usize) -> Result<ProjectedAlReport,ProjectedAlError<ResponseError3>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations; let mut last = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study; let cases = self.cases; let options = self.options;
                self.state.try_run(remaining.min(1),
                    &mut |x| sample(study,cases,x,options,&mut last,&mut **ledger.borrow_mut()),
                    |_| if ledger.borrow_mut().checkpoint("response-projected-control").is_ok() {ControlFlow::Continue(())} else {ControlFlow::Break(())})
            };
            if self.state.work().iterations > before {
                let accepted = last.expect("accepted optimizer step has complete physical evidence");
                assert_eq!(accepted.rho,self.state.point(),"accepted physics and projected point must match");
                self.study.operator.set_scales(&accepted.scales).expect("accepted scales are admitted");
                self.accepted = accepted;
                self.history.push(row(&self.accepted,self.state.work().iterations,self.options.volume_cap));
                self.control.checkpoint("response-projected-accepted").map_err(|e|ProjectedAlError::Evaluation(e.into()))?;
            }
            let report = result?;
            if report.stop!=ProjectedAlStop::IterationLimit || remaining<=1 {return Ok(report);}
            remaining-=1;
        }
    }
}
