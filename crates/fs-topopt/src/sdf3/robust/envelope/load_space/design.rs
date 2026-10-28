//! Optimize the continuous load set with the existing projected AL solver.
pub mod limits;
mod resume;
pub use resume::LoadSpaceDesignCheckpoint3;
use super::*;
use super::super::optimize::EnvelopeDesignOptions3;
use crate::sdf3::robust::{RobustDesignOptions3, continuation::controlled_sdf3_robust_feasible_start};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport, ProjectedAlSample,
    ProjectedAlState, ProjectedAlStop, ProjectedAlWork};
use std::{cell::RefCell, ops::ControlFlow};

#[derive(Debug, Clone)]
pub struct LoadSpaceDesignPoint3 { pub rho: Vec<f64>, pub evaluation: LoadSpaceEvaluation3 }
#[derive(Debug, Clone, PartialEq)]
pub struct LoadSpaceDesignIteration3 {
    pub iteration: usize,
    pub smooth_maximum: f64,
    pub maximum: f64,
    pub dilated_volume: f64,
    pub feasible: bool,
}
fn row(p: &LoadSpaceDesignPoint3, iteration: usize, options: EnvelopeDesignOptions3,
    cap: Option<&limits::LoadSpaceCap3>) -> LoadSpaceDesignIteration3 {
    LoadSpaceDesignIteration3 { iteration, smooth_maximum: p.evaluation.spectrum.smooth_maximum,
        maximum: p.evaluation.spectrum.maximum, dilated_volume: p.evaluation.physical.dilated_volume_fraction,
        feasible: cap.map_or_else(
            || p.evaluation.physical.dilated_volume_fraction/options.volume_limit-1.0 <= options.optimizer.tolerance,
            |cap| cap.admits(&p.evaluation, options.optimizer.tolerance)) }
}
fn sample<O: Sdf3Elasticity>(study: &mut CutDensityStudy3<O>, rho: &[f64], basis: &[&[f64]],
    policy: &LoadSpace3, options: EnvelopeDesignOptions3, minimum_material: bool,
    last: &mut Option<LoadSpaceDesignPoint3>,
    control: &mut SolveControl<'_>) -> Result<Option<ProjectedAlSample>, EvaluationStop> {
    let evaluation = study.evaluate_robust_load_space(rho, basis, options.robustness, policy, control)?;
    let sample = if minimum_material {
        ProjectedAlSample { objective: evaluation.physical.dilated_volume_fraction,
            gradient: evaluation.physical.dilated_volume_gradient.clone(),
            constraint: evaluation.spectrum.smooth_maximum-1.0,
            constraint_gradient: evaluation.gradient.clone() }
    } else {
        ProjectedAlSample { objective: evaluation.spectrum.smooth_maximum,
        gradient: evaluation.gradient.clone(),
        constraint: evaluation.physical.dilated_volume_fraction/options.volume_limit-1.0,
        constraint_gradient: evaluation.physical.dilated_volume_gradient.iter().map(|g| g/options.volume_limit).collect() }
    };
    *last = Some(LoadSpaceDesignPoint3 { rho: rho.to_vec(), evaluation });
    Ok(Some(sample))
}

fn admit<O: Sdf3Elasticity>(study: &CutDensityStudy3<O>, basis: &[&[f64]], rho: &[f64],
    policy: &LoadSpace3, options: EnvelopeDesignOptions3) -> Result<(), ProjectedAlError<EvaluationStop>> {
    options.optimizer.validate::<EvaluationStop>(study.cells())?;
    if basis.len() != policy.dimension() || !options.robustness.valid(study.params())
        || basis.iter().any(|f| f.len() != study.operator().n() || f.iter().any(|v| !v.is_finite()))
        || !options.volume_limit.is_finite() || options.volume_limit <= 0.0 || options.volume_limit > 1.0
        || !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
        || rho.len() != study.cells() || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0)
        || study.prescribed_solid_fraction() > options.volume_limit*(1.0+options.optimizer.tolerance) {
        return Err(ProjectedAlError::Invalid("invalid load-space design policy, basis or raw start"));
    }
    Ok(())
}

/// Fixed-grid, resumable optimization of f=Fq, ||q||<=1. Basis amplitudes,
/// correlations, reference, temperature and all material parameters stay fixed.
/// The smooth spectral maximum is minimized under the DILATED material cap.
/// This is not a finite-scenario maximum, CVaR, or a probability model.
///
/// KKT diagnostics belong to the accepted AL point; the best material-feasible
/// point may be earlier. Segments retain all AL state. Recovery is explicitly
/// separate and cannot replace the state or its history. All physics shares
/// the original control; incoming nominal operator scales remain unchanged.
pub struct LoadSpaceDesignStudy3<'a, 'callback, O: Sdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    basis: &'a [&'a [f64]],
    control: &'a mut SolveControl<'callback>,
    policy: LoadSpace3,
    options: EnvelopeDesignOptions3,
    cap: Option<limits::LoadSpaceCap3>,
    check: LoadSpaceGradientCheck3,
    state: ProjectedAlState,
    accepted: LoadSpaceDesignPoint3,
    best: Option<LoadSpaceDesignPoint3>,
    history: Vec<LoadSpaceDesignIteration3>,
    restoration_evaluations: usize,
}
impl<'a, 'callback, O: Sdf3Elasticity> LoadSpaceDesignStudy3<'a, 'callback, O> {
    pub fn new(study: &'a mut CutDensityStudy3<O>, basis: &'a [&'a [f64]], rho: &[f64],
        policy: LoadSpace3, options: EnvelopeDesignOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<EvaluationStop>> {
        Self::initialize(study, basis, rho, policy, options, None, control)
    }

    // Both public formulations share one numerical loop. Only objective,
    // constraint and feasible-incumbent ordering differ; the original public
    // constructor always selects the unchanged volume-constrained formulation.
    #[allow(clippy::too_many_arguments)]
    fn initialize(study: &'a mut CutDensityStudy3<O>, basis: &'a [&'a [f64]], rho: &[f64],
        policy: LoadSpace3, options: EnvelopeDesignOptions3, cap: Option<limits::LoadSpaceCap3>,
        control: &'a mut SolveControl<'callback>) -> Result<Self, ProjectedAlError<EvaluationStop>> {
        admit(study, basis, rho, &policy, options)?;
        let check = controlled_load_space_gradient_check(study, rho, basis, options.robustness,
            &policy, options.gradient_check, control).map_err(ProjectedAlError::Evaluation)?;
        if !check.numerical.passed() { return Err(ProjectedAlError::Invalid("load-space gradient gate failed")); }
        let mut last = None;
        let n = study.cells();
        let state = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(rho, &vec![options.density_floor; n], &vec![1.0; n], options.optimizer,
                &mut |rho| sample(study, rho, basis, &policy, options, cap.is_some(), &mut last, &mut **ledger.borrow_mut()),
                |_| if ledger.borrow_mut().checkpoint("sdf3-load-space-design-control").is_ok() {
                    ControlFlow::Continue(())
                } else { ControlFlow::Break(()) })?
        };
        let accepted = last.expect("admitted sample retains its complete physical fields");
        let first = row(&accepted, 0, options, cap.as_ref());
        let best = first.feasible.then(|| accepted.clone());
        Ok(Self { study, basis, control, policy, options, cap, check, state, accepted, best,
            history: vec![first], restoration_evaluations: 0 })
    }
    #[must_use] pub fn accepted(&self) -> &LoadSpaceDesignPoint3 { &self.accepted }
    #[must_use] pub fn best_feasible(&self) -> Option<&LoadSpaceDesignPoint3> { self.best.as_ref() }
    #[must_use] pub fn history(&self) -> &[LoadSpaceDesignIteration3] { &self.history }
    #[must_use] pub fn gradient_check(&self) -> &LoadSpaceGradientCheck3 { &self.check }
    #[must_use] pub fn work(&self) -> SolveWork { self.control.work() }
    #[must_use] pub fn optimizer_work(&self) -> ProjectedAlWork { self.state.work() }

    /// Fresh strict-volume recovery of the accepted raw controls. It can improve
    /// a stale feasible incumbent after a short AL run, but improvement is NOT
    /// guaranteed. The caller compares returned and accepted feasible results.
    pub fn recover_feasible(&mut self) -> Result<LoadSpaceDesignPoint3, EvaluationStop> {
        if self.cap.is_some() { return Err(failure("sdf3-load-space-volume-recovery-mode")); }
        self.control.checkpoint("sdf3-load-space-recovery")?;
        let rho = controlled_sdf3_robust_feasible_start(self.study, &self.accepted.rho,
            RobustDesignOptions3 { robustness: self.options.robustness, volume_limit: self.options.volume_limit,
                density_floor: self.options.density_floor, compliance_scale: 1.0,
                optimizer: self.options.optimizer }, self.control)?;
        let evaluation = self.study.evaluate_robust_load_space(&rho, self.basis,
            self.options.robustness, &self.policy, self.control)?;
        if evaluation.physical.dilated_volume_fraction > self.options.volume_limit {
            return Err(failure("sdf3-load-space-recovery-volume"));
        }
        Ok(LoadSpaceDesignPoint3 { rho, evaluation })
    }

    pub fn run(&mut self, additional_steps: usize) -> Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study; let basis = self.basis;
                let policy = &self.policy; let options = self.options; let minimum_material = self.cap.is_some();
                self.state.try_run(remaining.min(1),
                    &mut |rho| sample(study, rho, basis, policy, options, minimum_material, &mut last, &mut **ledger.borrow_mut()),
                    |_| if ledger.borrow_mut().checkpoint("sdf3-load-space-design-control").is_ok() {
                        ControlFlow::Continue(())
                    } else { ControlFlow::Break(()) })
            };
            if self.state.work().iterations > before {
                let accepted = last.expect("accepted load-space update retains complete fields");
                assert_eq!(accepted.rho, self.state.point(), "load-space fields detached from AL state");
                self.accepted = accepted;
                let next = row(&self.accepted, self.state.work().iterations, self.options, self.cap.as_ref());
                if next.feasible && self.best.as_ref().is_none_or(|best| if self.cap.is_some() {
                    next.dilated_volume < best.evaluation.physical.dilated_volume_fraction
                } else { next.smooth_maximum < best.evaluation.spectrum.smooth_maximum }) {
                    self.best = Some(self.accepted.clone());
                }
                self.history.push(next);
                self.control.checkpoint("sdf3-load-space-design-accepted").map_err(ProjectedAlError::Evaluation)?;
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 { return Ok(report); }
            remaining -= 1;
        }
    }
}

#[cfg(test)]
#[path = "design_tests.rs"]
mod tests;
