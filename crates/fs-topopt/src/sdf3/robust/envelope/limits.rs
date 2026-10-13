//! Minimum dilated material subject to ALL declared eroded compliance caps.
//! The smooth maximum is a conservative numerical surrogate, not an average.
pub mod start;
pub mod refinement;
pub mod refit;
use super::{EnvelopeEvaluation3, EnvelopeGradientCheck3, LoadEnvelope3,
    RobustOptions3, controlled_envelope_gradient_check};
use super::optimize::EnvelopeDesignPoint3;
use crate::sdf3::{CutDensityStudy3, Sdf3Elasticity};
use crate::{EvaluationStop, GradientCheckOptions, SolveControl, SolveWork};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions, ProjectedAlReport,
    ProjectedAlSample, ProjectedAlState, ProjectedAlStop, ProjectedAlWork};
use std::{cell::RefCell, ops::ControlFlow};

/// Immutable physical UPPER limits, in original load order. Unlike envelope
/// references used only for scaling, these are actual discrete requirements.
/// `temperature*ln(case_count)` is additional conservative smoothing, NOT a
/// tolerance for violating a cap. It must be < 1 to admit a strictly feasible
/// smooth problem even at zero compliance. No probability interpretation.
#[derive(Debug, Clone, PartialEq)]
pub struct ComplianceCaps3 { envelope: LoadEnvelope3 }
impl ComplianceCaps3 {
    pub fn new(caps: Vec<f64>, temperature: f64, max_load_cases: usize)
        -> Result<Self, &'static str> {
        let envelope = LoadEnvelope3::new(caps, temperature, max_load_cases)?;
        if envelope.bias_bound() >= 1.0 {
            return Err("compliance-cap smoothing consumes the entire admissible margin");
        }
        Ok(Self { envelope })
    }
    #[must_use] pub fn caps(&self) -> &[f64] { self.envelope.references() }
    #[must_use] pub fn envelope(&self) -> &LoadEnvelope3 { &self.envelope }
    /// Compare physical values directly: division can round a violating ratio
    /// to 1. Neither normalized AL tolerances nor smoothing loosen these caps.
    #[must_use]
    pub fn satisfied_by(&self, compliances: &[f64]) -> bool {
        compliances.len() == self.caps().len() && compliances.iter().zip(self.caps())
            .all(|(&c, &cap)| c.is_finite() && c >= 0.0 && c <= cap)
    }
    fn admits(&self, e: &EnvelopeEvaluation3, tolerance: f64) -> bool {
        e.policy == self.envelope && e.envelope.value.is_finite()
            && e.envelope.value - 1.0 <= tolerance
            && self.satisfied_by(&e.physical.eroded.objective.case_compliances)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ComplianceLimitedOptions3 {
    pub robustness: RobustOptions3,
    /// Exact raw-control lower bound, not a physical solid/void mask.
    pub density_floor: f64,
    /// Checks both the envelope constraint and dilated-material objective.
    pub gradient_check: GradientCheckOptions,
    /// f = dilated volume fraction; c = envelope(C_i/cap_i) - 1.
    /// Its feasibility tolerance is not permission to exceed a physical cap.
    pub optimizer: ProjectedAlOptions,
}
impl Default for ComplianceLimitedOptions3 {
    fn default() -> Self {
        Self { robustness: RobustOptions3::default(), density_floor: 1e-3,
            gradient_check: GradientCheckOptions::default(), optimizer: ProjectedAlOptions::default() }
    }
}

fn admit<O: Sdf3Elasticity>(study: &CutDensityStudy3<O>, forces: &[&[f64]], rho: &[f64],
    caps: &ComplianceCaps3, options: ComplianceLimitedOptions3)
    -> Result<(), ProjectedAlError<EvaluationStop>> {
    options.optimizer.validate::<EvaluationStop>(study.cells())?;
    if forces.len() != caps.caps().len() || !options.robustness.valid(study.params())
        || forces.iter().any(|f| f.len() != study.operator().n() || f.iter().any(|v| !v.is_finite()))
        || !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
        || rho.len() != study.cells() || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
        return Err(ProjectedAlError::Invalid("invalid compliance-limited design policy, forces or raw start"));
    }
    Ok(())
}
#[derive(Debug, Clone, PartialEq)]
pub struct ComplianceLimitedIteration3 {
    pub iteration: usize,
    /// Actual objective, NOT the envelope or nominal material fraction.
    pub dilated_volume: f64,
    pub envelope: f64,
    pub maximum_ratio: f64,
    pub limiting_case: usize,
    pub smooth_violation: f64,
    /// Each physical eroded compliance is <= its own declared cap AND the
    /// smooth surrogate is feasible at the optimizer's numerical tolerance.
    pub feasible: bool,
}
fn row(point: &EnvelopeDesignPoint3, iteration: usize, caps: &ComplianceCaps3,
    options: ComplianceLimitedOptions3) -> ComplianceLimitedIteration3 {
    let e = &point.evaluation;
    ComplianceLimitedIteration3 { iteration, dilated_volume: e.physical.dilated_volume_fraction,
        envelope: e.envelope.value, maximum_ratio: e.envelope.maximum,
        limiting_case: e.envelope.worst_case, smooth_violation: (e.envelope.value - 1.0).max(0.0),
        feasible: caps.admits(e, options.optimizer.tolerance) }
}
fn sample<O: Sdf3Elasticity>(study: &mut CutDensityStudy3<O>, forces: &[&[f64]], rho: &[f64],
    caps: &ComplianceCaps3, options: ComplianceLimitedOptions3,
    last: &mut Option<EnvelopeDesignPoint3>, control: &mut SolveControl<'_>)
    -> Result<Option<ProjectedAlSample>, EvaluationStop> {
    let evaluation = study.evaluate_robust_envelope(rho, forces, options.robustness, caps.envelope(), control)?;
    let sample = ProjectedAlSample {
        objective: evaluation.physical.dilated_volume_fraction,
        gradient: evaluation.physical.dilated_volume_gradient.clone(),
        constraint: evaluation.envelope.value - 1.0,
        constraint_gradient: evaluation.gradient.clone(),
    };
    *last = Some(EnvelopeDesignPoint3 { rho: rho.to_vec(), evaluation });
    Ok(Some(sample))
}

/// Fixed-grid, in-memory minimum-material design using the EXISTING projected
/// AL, eroded equilibrium solves and exact discrete envelope derivatives.
/// Every declared independent dead load is a constraint case; none is dropped
/// by averaging or probability weights. Geometry/supports/model/caps stay fixed.
///
/// Only accepted points satisfying ALL physical caps may become best_feasible.
/// It minimizes DILATED material, not envelope compliance. The final AL point
/// may differ and may violate a cap; KKT diagnostics describe that final point.
/// A converged AL status is therefore not itself a usable-design certificate.
/// No global optimum, continuum load guarantee or geometric tolerance claim.
/// Incoming nominal operator scales remain unchanged on returned outcomes.
pub struct ComplianceLimitedDesign3<'a, 'callback, O: Sdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    forces: &'a [&'a [f64]],
    control: &'a mut SolveControl<'callback>,
    caps: ComplianceCaps3,
    options: ComplianceLimitedOptions3,
    check: EnvelopeGradientCheck3,
    state: ProjectedAlState,
    accepted: EnvelopeDesignPoint3,
    best: Option<EnvelopeDesignPoint3>,
    history: Vec<ComplianceLimitedIteration3>,
}
impl<'a, 'callback, O: Sdf3Elasticity> ComplianceLimitedDesign3<'a, 'callback, O> {
    /// The start must be in the raw box; it need not yet meet the caps.
    /// No claim of infeasibility follows from having no feasible incumbent.
    /// Shape/policy refusals precede solves; gradient probes and initialization
    /// use the caller's existing cumulative physical work budget.
    pub fn new(study: &'a mut CutDensityStudy3<O>, forces: &'a [&'a [f64]], rho: &[f64],
        caps: ComplianceCaps3, options: ComplianceLimitedOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<EvaluationStop>> {
        admit(study, forces, rho, &caps, options)?;
        let n = study.cells();
        let check = controlled_envelope_gradient_check(study, rho, forces, options.robustness,
            caps.envelope(), options.gradient_check, control).map_err(ProjectedAlError::Evaluation)?;
        if !check.numerical.passed() {
            return Err(ProjectedAlError::Invalid("compliance-limited baseline gradient gate failed"));
        }
        let mut last = None;
        let state = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(rho, &vec![options.density_floor; n], &vec![1.0; n], options.optimizer,
                &mut |rho| sample(study, forces, rho, &caps, options, &mut last, &mut **ledger.borrow_mut()),
                |_| if ledger.borrow_mut().checkpoint("sdf3-compliance-limited-control").is_ok() {
                    ControlFlow::Continue(())
                } else { ControlFlow::Break(()) })?
        };
        let accepted = last.expect("initial sample has complete physical fields");
        let first = row(&accepted, 0, &caps, options);
        let best = first.feasible.then(|| accepted.clone());
        Ok(Self { study, forces, control, caps, options, check, state, accepted, best, history: vec![first] })
    }
    #[must_use] pub fn caps(&self) -> &ComplianceCaps3 { &self.caps }
    #[must_use] pub fn options(&self) -> ComplianceLimitedOptions3 { self.options }
    #[must_use] pub fn gradient_check(&self) -> &EnvelopeGradientCheck3 { &self.check }
    #[must_use] pub fn accepted(&self) -> &EnvelopeDesignPoint3 { &self.accepted }
    #[must_use] pub fn best_feasible(&self) -> Option<&EnvelopeDesignPoint3> { self.best.as_ref() }
    #[must_use] pub fn history(&self) -> &[ComplianceLimitedIteration3] { &self.history }
    #[must_use] pub fn work(&self) -> SolveWork { self.control.work() }
    #[must_use] pub fn optimizer_work(&self) -> ProjectedAlWork { self.state.work() }

    /// Freshly re-solve the retained raw design and recheck every physical cap.
    /// This is a new numerical solve, NOT an independent solver or continuum
    /// certificate. It charges the physical budget, not optimization callback
    /// counters, and never mutates the trajectory or cached feasible incumbent.
    pub fn recheck_best(&mut self) -> Result<Option<EnvelopeDesignPoint3>, EvaluationStop> {
        let Some(best) = &self.best else { return Ok(None); };
        self.control.checkpoint("sdf3-compliance-limited-recheck")?;
        let evaluation = self.study.evaluate_robust_envelope(&best.rho, self.forces,
            self.options.robustness, self.caps.envelope(), self.control)?;
        if !self.caps.admits(&evaluation, self.options.optimizer.tolerance) {
            return Err(EvaluationStop::Breakdown { stage: "sdf3-compliance-limited-recheck" });
        }
        Ok(Some(EnvelopeDesignPoint3 { rho: best.rho.clone(), evaluation }))
    }

    /// Continue the same AL state for at most additional_steps accepted primal
    /// updates; zero steps uses cached fields. Trials and numerical failures
    /// remain charged. An accepted transition is paired with its fields before
    /// propagating a late cancellation; rejected trials never replace the best.
    pub fn run(&mut self, additional_steps: usize) -> Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study;
                let forces = self.forces;
                let caps = &self.caps;
                let options = self.options;
                self.state.try_run(remaining.min(1),
                    &mut |rho| sample(study, forces, rho, caps, options, &mut last, &mut **ledger.borrow_mut()),
                    |_| if ledger.borrow_mut().checkpoint("sdf3-compliance-limited-control").is_ok() {
                        ControlFlow::Continue(())
                    } else { ControlFlow::Break(()) })
            };
            if self.state.work().iterations > before {
                let accepted = last.expect("accepted material step has complete fields");
                assert_eq!(accepted.rho, self.state.point(), "compliance-limited fields detached from AL state");
                self.accepted = accepted;
                let next = row(&self.accepted, self.state.work().iterations, &self.caps, self.options);
                if next.feasible && self.best.as_ref().is_none_or(|best|
                    next.dilated_volume < best.evaluation.physical.dilated_volume_fraction) {
                    self.best = Some(self.accepted.clone());
                }
                self.history.push(next);
                self.control.checkpoint("sdf3-compliance-limited-accepted").map_err(ProjectedAlError::Evaluation)?;
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 { return Ok(report); }
            remaining -= 1;
        }
    }
}

#[cfg(test)]
#[path = "limits_tests.rs"]
mod tests;
