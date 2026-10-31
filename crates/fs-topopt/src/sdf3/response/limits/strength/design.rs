//! Minimum material with force/displacement requirements AND sampled strength.
//! Use the original projected AL engine; no dense Jacobian, KKT solve or second
//! optimizer algorithm. Intermediate accepted AL points need not be feasible.
use super::*;
use std::{cell::RefCell, ops::ControlFlow};
use crate::sdf3::{PhysicalRegion3, design::monotone_stress_branch};
use fs_ascent::projected_al::{ProjectedAlCheckpoint, ProjectedAlError, ProjectedAlOptions,
    ProjectedAlReport, ProjectedAlSample, ProjectedAlState, ProjectedAlStop, ProjectedAlWork};

type DesignError = ProjectedAlError<ResponseStressError3>;

#[derive(Debug, Clone, Copy)]
pub struct ResponseStressDesignOptions3 {
    pub limits: ResponseStressOptions3,
    /// Exact raw-density bound. The projected floor must also remain above the
    /// same ersatz stress turnover admitted by the existing stress-design owner.
    pub density_floor: f64,
    pub optimizer: ProjectedAlOptions,
}
impl Default for ResponseStressDesignOptions3 {
    fn default() -> Self {
        Self { limits: ResponseStressOptions3::default(), density_floor: 0.1,
            optimizer: ProjectedAlOptions::default() }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct ResponseStressDesignIteration3 {
    pub iteration: usize,
    pub volume_fraction: f64,
    pub constraint: f64,
    pub worst_response_residual: f64,
    pub stress_aggregate: f64,
    pub sampled_physical_max: f64,
    pub feasible: bool,
}
fn row(e: &ResponseStressEvaluation3, iteration: usize, tolerance: f64) -> ResponseStressDesignIteration3 {
    ResponseStressDesignIteration3 { iteration, volume_fraction: e.response.volume_fraction,
        constraint: e.constraint, worst_response_residual: e.response.worst_residual,
        stress_aggregate: e.stress.aggregate, sampled_physical_max: e.stress.sampled_physical_max,
        feasible: e.feasible(tolerance) }
}

/// Search state without cached physical fields. Consumers must bind this to the
/// original geometry, material law, complete experiment family and policy.
#[derive(Debug, Clone)]
pub struct ResponseStressDesignCheckpoint3 {
    pub optimizer: ProjectedAlCheckpoint,
    /// Complete prepared body map and native source attribution, when attached.
    pub material_body_load_fingerprint: Option<fs_blake3::ContentHash>,
    pub best_feasible_density: Option<Vec<f64>>,
    pub history: Vec<ResponseStressDesignIteration3>,
    pub restoration_evaluations: usize,
}
impl ResponseStressDesignCheckpoint3 {
    #[must_use]
    pub fn restoration_cost(&self) -> usize {
        1 + usize::from(self.best_feasible_density.as_ref()
            .is_some_and(|rho| rho != &self.optimizer.point))
    }
}
fn sample<O: AdaptiveSdf3Elasticity>(study: &mut CutDensityStudy3<O>,
    cases: &[ResponseLimitCase3<'_>], rho: &[f64], options: ResponseStressDesignOptions3,
    last: &mut Option<ResponseStressEvaluation3>, control: &mut SolveControl<'_>)
    -> Result<Option<ProjectedAlSample>, ResponseStressError3> {
    control.checkpoint("response-stress-design-evaluate")?;
    let e = study.evaluate_response_stress_limits(rho, cases, options.limits, control)?;
    // Identical stress-model admission to the existing single-cap optimizer.
    // Protected voids are fixed model data, not a free ersatz foldback branch.
    if e.stress.projected_rho.iter().enumerate().any(|(i, r)| {
        let fixed_void = *r == 0.0 && study.physical_regions()
            .is_some_and(|labels| labels[i] == PhysicalRegion3::Void);
        !fixed_void && !monotone_stress_branch(*r, study.params(), options.limits.stress.relaxation_power)
    }) {
        return Err(ResponseStressError3::Invalid("joint design entered the ersatz stress turnover"));
    }
    let sample = ProjectedAlSample { objective: e.response.volume_fraction,
        gradient: e.response.volume_gradient.clone(), constraint: e.constraint,
        constraint_gradient: e.constraint_gradient.clone() };
    *last = Some(e);
    Ok(Some(sample))
}

/// In-memory continuation for a fixed geometry, material law and experiment
/// family. A single joint constraint differentiates ALL original requirements.
/// Bound projections and AL steps are delegated to fs-ascent. Returned physics
/// is checked separately for original response and sampled-stress feasibility.
///
/// Every probe restores material; only an accepted transition installs scales.
/// The least-volume feasible incumbent remains separate from the latest AL
/// point and its KKT report. No global optimum, continuum safety or actuator
/// energy is inferred. This owner does not add disk restart or adaptive transfer;
/// the response-only checkpoint and DWR APIs do not encode this joint objective.
pub struct ResponseStressDesignStudy3<'a, 'callback, O: AdaptiveSdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    cases: &'a [ResponseLimitCase3<'a>],
    control: &'a mut SolveControl<'callback>,
    options: ResponseStressDesignOptions3,
    state: ProjectedAlState,
    accepted: ResponseStressEvaluation3,
    best_feasible: Option<ResponseStressEvaluation3>,
    history: Vec<ResponseStressDesignIteration3>,
    restoration_evaluations: usize,
}
impl<'a, 'callback, O: AdaptiveSdf3Elasticity> ResponseStressDesignStudy3<'a, 'callback, O> {
    pub fn new(study: &'a mut CutDensityStudy3<O>, cases: &'a [ResponseLimitCase3<'a>],
        rho: &[f64], options: ResponseStressDesignOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, DesignError> {
        let n = study.cells();
        options.optimizer.validate(n)?;
        if !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
            || rho.len() != n || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
            return Err(ProjectedAlError::Invalid("invalid joint design bounds or initial density"));
        }
        admit(study, rho, cases, options.limits).map_err(ProjectedAlError::Evaluation)?;
        let p = study.params();
        let floor = crate::filter::heaviside(options.density_floor, p.beta, p.eta);
        if !monotone_stress_branch(floor, p, options.limits.stress.relaxation_power) {
            return Err(ProjectedAlError::Invalid("joint density floor permits ersatz stress foldback"));
        }
        let mut last = None;
        let state = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(rho, &vec![options.density_floor;n], &vec![1.0;n], options.optimizer,
                &mut |r| sample(study, cases, r, options, &mut last, &mut **ledger.borrow_mut()),
                |_| if ledger.borrow_mut().checkpoint("response-stress-design-control").is_ok() {
                    ControlFlow::Continue(())
                } else { ControlFlow::Break(()) })?
        };
        control.checkpoint("response-stress-design-initialize")
            .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
        let accepted = last.expect("initial joint evaluation contains complete evidence");
        study.operator.set_scales(&accepted.stress.scales).expect("evaluated scales remain valid");
        let first = row(&accepted, 0, options.optimizer.tolerance);
        let best_feasible = first.feasible.then(|| accepted.clone());
        Ok(Self { study, cases, control, options, state, accepted, best_feasible,
            history:vec![first], restoration_evaluations:0 })
    }
    #[must_use] pub fn accepted(&self) -> &ResponseStressEvaluation3 { &self.accepted }
    #[must_use] pub fn best_feasible(&self) -> Option<&ResponseStressEvaluation3> { self.best_feasible.as_ref() }
    #[must_use] pub fn study(&self) -> &CutDensityStudy3<O> { self.study }
    #[must_use] pub fn history(&self) -> &[ResponseStressDesignIteration3] { &self.history }
    #[must_use] pub fn point(&self) -> &[f64] { self.state.point() }
    #[must_use] pub fn work(&self) -> SolveWork { self.control.work() }
    #[must_use] pub fn optimizer_work(&self) -> ProjectedAlWork { self.state.work() }
    #[must_use] pub fn feasible(&self) -> bool { self.accepted.feasible(self.options.optimizer.tolerance) }

    #[must_use]
    pub fn checkpoint(&self) -> ResponseStressDesignCheckpoint3 {
        ResponseStressDesignCheckpoint3 {
            optimizer:self.state.checkpoint(),
            material_body_load_fingerprint:self.study.material_body_load_fingerprint(),
            best_feasible_density:self.best_feasible.as_ref().map(|e|e.response.rho.clone()),
            history:self.history.clone(),
            restoration_evaluations:self.restoration_evaluations,
        }
    }

    /// Restore the projected-AL search only after freshly re-solving its retained
    /// physical endpoints under the caller-supplied original fixed-grid problem.
    pub fn restore(study:&'a mut CutDensityStudy3<O>, cases:&'a [ResponseLimitCase3<'a>],
        checkpoint:ResponseStressDesignCheckpoint3, options:ResponseStressDesignOptions3,
        control:&'a mut SolveControl<'callback>) -> Result<Self,DesignError> {
        if checkpoint.material_body_load_fingerprint != study.material_body_load_fingerprint() {
            return Err(ProjectedAlError::Invalid("joint restart changed its material body load"));
        }
        let n=study.cells();
        let cost=checkpoint.restoration_cost();
        let restoration_evaluations=checkpoint.restoration_evaluations.checked_add(cost)
            .ok_or(ProjectedAlError::Invalid("joint restoration counter overflow"))?;
        let mut optimizer=checkpoint.optimizer;
        let previous_evaluations=optimizer.work.evaluations;
        optimizer.work.evaluations=previous_evaluations.checked_add(cost)
            .ok_or(ProjectedAlError::Invalid("joint restoration evaluation overflow"))?;
        let restored=ProjectedAlState::try_restore(optimizer,
            &vec![options.density_floor;n],&vec![1.0;n],options.optimizer)?;
        let history=checkpoint.history;
        if checkpoint.restoration_evaluations>=previous_evaluations
            || history.len()!=restored.work().iterations+1
            || history.iter().enumerate().any(|(i,r)| {
                r.iteration!=i || !r.volume_fraction.is_finite()
                    || !(0.0..=1.0).contains(&r.volume_fraction)
                    || !r.constraint.is_finite() || !r.worst_response_residual.is_finite()
                    || !r.stress_aggregate.is_finite() || r.stress_aggregate<0.0
                    || !r.sampled_physical_max.is_finite() || r.sampled_physical_max<0.0
            }) {
            return Err(ProjectedAlError::Invalid("invalid joint checkpoint history or work"));
        }
        let best_row=history.iter().filter(|r|r.feasible).fold(None,
            |best:Option<&ResponseStressDesignIteration3>,r| {
                if best.is_none_or(|b|r.volume_fraction<b.volume_fraction){Some(r)}else{best}
            });
        if best_row.is_some()!=checkpoint.best_feasible_density.is_some() {
            return Err(ProjectedAlError::Invalid("joint checkpoint lost its feasible incumbent"));
        }
        let previous_scales=study.operator.scales().to_vec();
        let mut result=Self::new(study,cases,restored.point(),options,control)?;
        let rebuilt=(|| {
            if result.state.sample()!=restored.sample()
                || history.last()!=Some(&row(&result.accepted,restored.work().iterations,
                    options.optimizer.tolerance)) {
                return Err(ProjectedAlError::Invalid("restored joint point changed its physical response"));
            }
            let best=match checkpoint.best_feasible_density {
                None=>None,
                Some(rho) if rho==result.accepted.response.rho=>Some(result.accepted.clone()),
                Some(rho)=>{
                    if rho.len()!=n || rho.iter().any(|r|!r.is_finite()
                        || !(options.density_floor..=1.0).contains(r)) {
                        return Err(ProjectedAlError::Invalid("invalid retained joint feasible density"));
                    }
                    let mut evaluated=None;
                    sample(result.study,cases,&rho,options,&mut evaluated,result.control)
                        .map_err(ProjectedAlError::Evaluation)?;
                    evaluated
                }
            };
            if let (Some(e),Some(expected))=(&best,best_row) {
                if row(e,expected.iteration,options.optimizer.tolerance)!=*expected {
                    return Err(ProjectedAlError::Invalid("retained joint feasible design changed its physical response"));
                }
            }
            Ok(best)
        })();
        let best=match rebuilt {
            Ok(best)=>best,
            Err(error)=>{
                result.study.operator.set_scales(&previous_scales)
                    .expect("incoming joint scales remain valid");
                return Err(error);
            }
        };
        result.state=restored;
        result.best_feasible=best;
        result.history=history;
        result.restoration_evaluations=restoration_evaluations;
        result.study.operator.set_scales(&result.accepted.stress.scales)
            .expect("restored accepted joint scales remain valid");
        Ok(result)
    }

    /// Continue without resetting multipliers, spectral steps or spent work.
    /// Late cancellation still installs a completed accepted transition before
    /// returning its error. Zero steps inspect cached state without PDE work.
    pub fn run(&mut self, additional_steps: usize) -> Result<ProjectedAlReport, DesignError> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study;
                let cases = self.cases;
                let options = self.options;
                self.state.try_run(remaining.min(1),
                    &mut |r| sample(study, cases, r, options, &mut last, &mut **ledger.borrow_mut()),
                    |_| if ledger.borrow_mut().checkpoint("response-stress-design-control").is_ok() {
                        ControlFlow::Continue(())
                    } else { ControlFlow::Break(()) })
            };
            if self.state.work().iterations > before {
                let accepted = last.expect("accepted joint point contains complete evidence");
                assert_eq!(accepted.response.rho, self.state.point(), "joint physical fields must match optimizer state");
                self.study.operator.set_scales(&accepted.stress.scales).expect("accepted scales remain valid");
                self.accepted = accepted;
                let next = row(&self.accepted, self.state.work().iterations, self.options.optimizer.tolerance);
                if next.feasible && self.best_feasible.as_ref()
                    .is_none_or(|b| self.accepted.response.volume_fraction < b.response.volume_fraction) {
                    self.best_feasible = Some(self.accepted.clone());
                }
                self.history.push(next);
                self.control.checkpoint("response-stress-design-accepted")
                    .map_err(|e| ProjectedAlError::Evaluation(e.into()))?;
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 { return Ok(report); }
            remaining -= 1;
        }
    }
    /// Explicitly install the feasible incumbent before consuming this session.
    /// Its derivatives/fields are retained; the latest AL KKT report does not
    /// become a stationarity statement about this potentially older design.
    pub fn into_best_feasible(self) -> Result<ResponseStressEvaluation3, DesignError> {
        let best = self.best_feasible.ok_or(ProjectedAlError::Invalid("no joint feasible design was accepted"))?;
        self.study.operator.set_scales(&best.stress.scales).expect("retained scales remain valid");
        Ok(best)
    }
}
