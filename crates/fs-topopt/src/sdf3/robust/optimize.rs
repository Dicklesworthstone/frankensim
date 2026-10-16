//! Fixed-geometry robust compliance design using the existing O(n)-storage AL.
use super::{RobustEvaluation3, RobustOptions3};
use super::super::{CutDensityStudy3, Sdf3Elasticity, assert_loads};
use crate::pipeline::LoadCase;
use crate::{EvaluationStop, SolveControl, SolveWork};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions, ProjectedAlReport,
    ProjectedAlSample, ProjectedAlState, ProjectedAlStop, ProjectedAlWork};
use std::{cell::RefCell, ops::ControlFlow};

#[derive(Debug, Clone, Copy)]
pub struct RobustDesignOptions3 {
    pub robustness: RobustOptions3,
    /// Material cap on the DILATED field, as a fraction of the cut domain.
    pub volume_limit: f64,
    /// Exact raw-control box [density_floor, 1], not a physical material mask.
    pub density_floor: f64,
    /// Explicit positive objective scale: minimize eroded compliance / scale.
    /// This changes conditioning and KKT coordinates, not the physical result.
    /// A baseline compliance is often useful; never recomputed during a run.
    pub compliance_scale: f64,
    /// The inequality is dilated_volume / volume_limit - 1. Feasibility and
    /// KKT tolerances are expressed in these normalized problem coordinates.
    pub optimizer: ProjectedAlOptions,
}
impl Default for RobustDesignOptions3 {
    fn default() -> Self {
        Self { robustness: RobustOptions3::default(), volume_limit: 0.5,
            density_floor: 1e-3, compliance_scale: 1.0,
            optimizer: ProjectedAlOptions::default() }
    }
}

/// Raw design and all of its matching, independently solved physical fields.
#[derive(Debug, Clone)]
pub struct RobustDesignPoint3 {
    pub rho: Vec<f64>,
    pub evaluation: RobustEvaluation3,
}
#[derive(Debug, Clone, PartialEq)]
pub struct RobustDesignIteration3 {
    pub iteration: usize,
    /// Unscaled, weighted independently solved compliance of eroded material.
    pub eroded_compliance: f64,
    pub eroded_volume_fraction: f64,
    pub nominal_volume_fraction: f64,
    pub dilated_volume_fraction: f64,
    pub constraint_violation: f64,
    /// Satisfies the dilated constraint at the explicit optimizer tolerance.
    pub feasible: bool,
}
fn row(e: &RobustEvaluation3, iteration: usize, options: RobustDesignOptions3) -> RobustDesignIteration3 {
    let constraint_violation = (e.dilated_volume_fraction / options.volume_limit - 1.0).max(0.0);
    RobustDesignIteration3 { iteration, eroded_compliance: e.eroded.objective.compliance,
        eroded_volume_fraction: e.eroded.volume_fraction,
        nominal_volume_fraction: e.nominal_volume_fraction,
        dilated_volume_fraction: e.dilated_volume_fraction, constraint_violation,
        feasible: constraint_violation <= options.optimizer.tolerance }
}
fn sample<O: Sdf3Elasticity>(study: &mut CutDensityStudy3<O>, loads: &[LoadCase<'_>],
    rho: &[f64], options: RobustDesignOptions3, last: &mut Option<RobustDesignPoint3>,
    control: &mut SolveControl<'_>) -> Result<Option<ProjectedAlSample>, EvaluationStop> {
    let evaluation = study.evaluate_robust(rho, loads, options.robustness, control)?;
    let sample = ProjectedAlSample {
        objective: evaluation.eroded.objective.compliance / options.compliance_scale,
        gradient: evaluation.eroded.objective.gradient.iter().map(|g| g / options.compliance_scale).collect(),
        constraint: evaluation.dilated_volume_fraction / options.volume_limit - 1.0,
        constraint_gradient: evaluation.dilated_volume_gradient.iter().map(|g| g / options.volume_limit).collect(),
    };
    *last = Some(RobustDesignPoint3 { rho: rho.to_vec(), evaluation });
    Ok(Some(sample))
}

/// Resumable, exclusively bound raw-SDF robust design. No new optimizer, dense
/// KKT matrix, numerical-gradient optimization or surrogate physics is used.
/// Geometry, independent loads, regions, robustness and material model stay fixed.
///
/// Accepted AL steps need not be feasible or decrease original compliance.
/// `best_feasible()` separately retains the lowest-compliance ACCEPTED design
/// meeting the declared dilated-volume tolerance. Rejected or cancelled trial
/// fields never replace it. Both physics and optimizer work remain charged.
///
/// Every evaluation restores incoming scales. The borrowed study therefore
/// keeps its original scales even on success; use the returned point's eroded
/// fields, not those scales, for the optimized solution. A returned run error
/// preserves the session's latest accepted state for inspection and continuation.
/// This solves the discrete threshold model, not a continuum safety guarantee.
pub struct RobustDesignStudy3<'a, 'callback, O: Sdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    loads: &'a [LoadCase<'a>],
    control: &'a mut SolveControl<'callback>,
    options: RobustDesignOptions3,
    state: ProjectedAlState,
    accepted: RobustDesignPoint3,
    best_feasible: Option<RobustDesignPoint3>,
    history: Vec<RobustDesignIteration3>,
}
impl<'a, 'callback, O: Sdf3Elasticity> RobustDesignStudy3<'a, 'callback, O> {
    /// The start may violate the material cap, but must satisfy the raw box.
    /// Model/load shape follows the evaluator's panic-before-mutation contract;
    /// invalid design policy or protected-solid infeasibility returns Invalid.
    pub fn new(study: &'a mut CutDensityStudy3<O>, loads: &'a [LoadCase<'a>], rho: &[f64],
        options: RobustDesignOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<EvaluationStop>> {
        let n = study.cells();
        options.optimizer.validate(n)?;
        if !options.robustness.valid(study.params())
            || !options.volume_limit.is_finite() || options.volume_limit <= 0.0 || options.volume_limit > 1.0
            || !options.compliance_scale.is_finite() || options.compliance_scale <= 0.0
            || !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
            || rho.len() != n || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
            return Err(ProjectedAlError::Invalid("invalid robust design policy or initial density"));
        }
        assert_loads(study, loads);
        if study.prescribed_solid_fraction() > options.volume_limit * (1.0 + options.optimizer.tolerance) {
            return Err(ProjectedAlError::Invalid("prescribed solid already exceeds the dilated material cap"));
        }
        let mut last = None;
        let state = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(rho, &vec![options.density_floor; n], &vec![1.0; n], options.optimizer,
                &mut |rho| sample(study, loads, rho, options, &mut last, &mut **ledger.borrow_mut()),
                |_| if ledger.borrow_mut().checkpoint("sdf3-robust-design-control").is_ok() {
                    ControlFlow::Continue(())
                } else { ControlFlow::Break(()) })?
        };
        control.checkpoint("sdf3-robust-design-initialize").map_err(ProjectedAlError::Evaluation)?;
        let accepted = last.expect("initial optimizer sample contains complete robust fields");
        let first = row(&accepted.evaluation, 0, options);
        let best_feasible = first.feasible.then(|| accepted.clone());
        Ok(Self { study, loads, control, options, state, accepted, best_feasible, history: vec![first] })
    }
    #[must_use] pub fn options(&self) -> RobustDesignOptions3 { self.options }
    #[must_use] pub fn accepted(&self) -> &RobustDesignPoint3 { &self.accepted }
    #[must_use] pub fn best_feasible(&self) -> Option<&RobustDesignPoint3> { self.best_feasible.as_ref() }
    #[must_use] pub fn point(&self) -> &[f64] { self.state.point() }
    #[must_use] pub fn history(&self) -> &[RobustDesignIteration3] { &self.history }
    #[must_use] pub fn work(&self) -> SolveWork { self.control.work() }
    #[must_use] pub fn optimizer_work(&self) -> ProjectedAlWork { self.state.work() }
    #[must_use] pub fn feasible(&self) -> bool {
        self.state.sample().constraint.max(0.0) <= self.options.optimizer.tolerance
    }

    /// At most `additional_steps` accepted primal steps, sharing the original
    /// evaluation/Krylov budgets across segments. Zero steps uses cached fields
    /// and does not update the multipliers. Never restart AL between segments.
    pub fn run(&mut self, additional_steps: usize) -> Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study;
                let loads = self.loads;
                let options = self.options;
                self.state.try_run(remaining.min(1),
                    &mut |rho| sample(study, loads, rho, options, &mut last, &mut **ledger.borrow_mut()),
                    |_| if ledger.borrow_mut().checkpoint("sdf3-robust-design-control").is_ok() {
                        ControlFlow::Continue(())
                    } else { ControlFlow::Break(()) })
            };
            // try_run may stop at its next poll AFTER accepting one point.
            // Pair that state with its fields before propagating any such error.
            if self.state.work().iterations > before {
                let accepted = last.expect("accepted optimizer point contains complete robust fields");
                assert_eq!(accepted.rho, self.state.point(), "robust fields must match the accepted optimizer point");
                self.accepted = accepted;
                let next = row(&self.accepted.evaluation, self.state.work().iterations, self.options);
                if next.feasible && self.best_feasible.as_ref().is_none_or(|best|
                    self.accepted.evaluation.eroded.objective.compliance < best.evaluation.eroded.objective.compliance) {
                    self.best_feasible = Some(self.accepted.clone());
                }
                self.history.push(next);
                self.control.checkpoint("sdf3-robust-design-accepted").map_err(ProjectedAlError::Evaluation)?;
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 { return Ok(report); }
            remaining -= 1;
        }
    }
}
