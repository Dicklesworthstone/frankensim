//! Fixed-grid worst-load design using the existing linear-storage projected AL.
use super::{EnvelopeEvaluation3, EnvelopeGradientCheck3, LoadEnvelope3,
    RobustOptions3, controlled_envelope_gradient_check};
use crate::sdf3::{CutDensityStudy3, Sdf3Elasticity};
use crate::{EvaluationStop, GradientCheckOptions, SolveControl, SolveWork};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions, ProjectedAlReport,
    ProjectedAlSample, ProjectedAlState, ProjectedAlStop, ProjectedAlWork};
use std::{cell::RefCell, ops::ControlFlow};

#[derive(Debug, Clone, Copy)]
pub struct EnvelopeDesignOptions3 {
    pub robustness: RobustOptions3,
    /// Material cap is on the DILATED field, as a fraction of cut volume.
    pub volume_limit: f64,
    pub density_floor: f64,
    /// Required baseline envelope AND dilated-volume numerical gradient gate.
    pub gradient_check: GradientCheckOptions,
    /// Objective is the dimensionless envelope; c = dilated_volume/cap - 1.
    pub optimizer: ProjectedAlOptions,
}
impl Default for EnvelopeDesignOptions3 {
    fn default() -> Self {
        Self { robustness: RobustOptions3::default(), volume_limit: 0.5, density_floor: 1e-3,
            gradient_check: GradientCheckOptions::default(), optimizer: ProjectedAlOptions::default() }
    }
}
#[derive(Debug, Clone)]
pub struct EnvelopeDesignPoint3 { pub rho: Vec<f64>, pub evaluation: EnvelopeEvaluation3 }
#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopeDesignIteration3 {
    pub iteration: usize,
    pub envelope: f64,
    pub maximum: f64,
    pub worst_case: usize,
    pub dilated_volume: f64,
    /// Only material feasibility; references are not extra hard constraints.
    pub feasible: bool,
}
fn row(point: &EnvelopeDesignPoint3, iteration: usize, options: EnvelopeDesignOptions3) -> EnvelopeDesignIteration3 {
    let e = &point.evaluation;
    EnvelopeDesignIteration3 { iteration, envelope: e.envelope.value, maximum: e.envelope.maximum,
        worst_case: e.envelope.worst_case, dilated_volume: e.physical.dilated_volume_fraction,
        feasible: (e.physical.dilated_volume_fraction/options.volume_limit - 1.0).max(0.0) <= options.optimizer.tolerance }
}
fn sample<O: Sdf3Elasticity>(study: &mut CutDensityStudy3<O>, forces: &[&[f64]],
    rho: &[f64], policy: &LoadEnvelope3, options: EnvelopeDesignOptions3,
    last: &mut Option<EnvelopeDesignPoint3>, control: &mut SolveControl<'_>)
    -> Result<Option<ProjectedAlSample>, EvaluationStop> {
    let evaluation = study.evaluate_robust_envelope(rho, forces, options.robustness, policy, control)?;
    let sample = ProjectedAlSample { objective: evaluation.envelope.value,
        gradient: evaluation.gradient.clone(),
        constraint: evaluation.physical.dilated_volume_fraction/options.volume_limit - 1.0,
        constraint_gradient: evaluation.physical.dilated_volume_gradient.iter().map(|g| g/options.volume_limit).collect() };
    *last = Some(EnvelopeDesignPoint3 { rho: rho.to_vec(), evaluation });
    Ok(Some(sample))
}

/// In-memory resumable worst-load study. All physics, baseline gradient probes
/// and rejected trials share the caller's original budget. The AL state, loads,
/// model, references and temperature remain bound across run segments.
///
/// `best_feasible` minimizes the SMOOTH ENVELOPE among accepted material-feasible
/// designs, not the uniform mean or necessarily the unsmoothed maximum. Both
/// scores are reported. KKT diagnostics belong to the final AL iterate, not to
/// an earlier incumbent. A finite list of loads is not distributional robustness.
/// Incoming nominal scales remain unchanged; read the returned eroded fields.
pub struct EnvelopeDesignStudy3<'a, 'callback, O: Sdf3Elasticity> {
    study: &'a mut CutDensityStudy3<O>,
    forces: &'a [&'a [f64]],
    control: &'a mut SolveControl<'callback>,
    policy: LoadEnvelope3,
    options: EnvelopeDesignOptions3,
    check: EnvelopeGradientCheck3,
    state: ProjectedAlState,
    accepted: EnvelopeDesignPoint3,
    best: Option<EnvelopeDesignPoint3>,
    history: Vec<EnvelopeDesignIteration3>,
}
impl<'a, 'callback, O: Sdf3Elasticity> EnvelopeDesignStudy3<'a, 'callback, O> {
    /// The raw start may exceed the material cap; AL feasibility is not assumed.
    /// Invalid policy returns before solves; modeling/load errors retain the
    /// evaluator's panic-before-mutation contract. A failed gradient gate refuses
    /// optimization; its numerical probes still consume physical work.
    pub fn new(study: &'a mut CutDensityStudy3<O>, forces: &'a [&'a [f64]], rho: &[f64],
        policy: LoadEnvelope3, options: EnvelopeDesignOptions3, control: &'a mut SolveControl<'callback>)
        -> Result<Self, ProjectedAlError<EvaluationStop>> {
        let n = study.cells();
        options.optimizer.validate::<EvaluationStop>(n)?;
        if forces.len() != policy.references().len() || !options.robustness.valid(study.params())
            || !options.volume_limit.is_finite() || options.volume_limit <= 0.0 || options.volume_limit > 1.0
            || !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
            || rho.len() != n || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
            return Err(ProjectedAlError::Invalid("invalid envelope design policy, cases or initial density"));
        }
        if study.prescribed_solid_fraction() > options.volume_limit*(1.0 + options.optimizer.tolerance) {
            return Err(ProjectedAlError::Invalid("protected solid exceeds the dilated material cap"));
        }
        let check = controlled_envelope_gradient_check(study, rho, forces, options.robustness,
            &policy, options.gradient_check, control).map_err(ProjectedAlError::Evaluation)?;
        if !check.numerical.passed() { return Err(ProjectedAlError::Invalid("envelope baseline gradient gate failed")); }
        let mut last = None;
        let state = {
            let ledger = RefCell::new(&mut *control);
            ProjectedAlState::try_new(rho, &vec![options.density_floor; n], &vec![1.0; n], options.optimizer,
                &mut |rho| sample(study, forces, rho, &policy, options, &mut last, &mut **ledger.borrow_mut()),
                |_| if ledger.borrow_mut().checkpoint("sdf3-envelope-design-control").is_ok() {
                    ControlFlow::Continue(())
                } else { ControlFlow::Break(()) })?
        };
        let accepted = last.expect("admitted initial sample has complete fields");
        let first = row(&accepted, 0, options);
        let best = first.feasible.then(|| accepted.clone());
        Ok(Self { study, forces, control, policy, options, check, state, accepted, best, history: vec![first] })
    }
    #[must_use] pub fn policy(&self) -> &LoadEnvelope3 { &self.policy }
    #[must_use] pub fn options(&self) -> EnvelopeDesignOptions3 { self.options }
    #[must_use] pub fn gradient_check(&self) -> &EnvelopeGradientCheck3 { &self.check }
    #[must_use] pub fn accepted(&self) -> &EnvelopeDesignPoint3 { &self.accepted }
    #[must_use] pub fn best_feasible(&self) -> Option<&EnvelopeDesignPoint3> { self.best.as_ref() }
    #[must_use] pub fn history(&self) -> &[EnvelopeDesignIteration3] { &self.history }
    #[must_use] pub fn work(&self) -> SolveWork { self.control.work() }
    #[must_use] pub fn optimizer_work(&self) -> ProjectedAlWork { self.state.work() }

    /// Recover a separately solved, material-feasible candidate from the current
    /// accepted AL point. A bounded AL run can have improved but infeasible
    /// iterates and retain only its initial feasible design. This explicit
    /// postprocessing does NOT change AL multipliers, history, callback counts,
    /// accepted fields or best_feasible; it is not an optimizer convergence step.
    ///
    /// Reuse the robust DILATED-volume restoration toward the declared raw floor
    /// (at most 64 bisections plus two endpoint filter evaluations), then freshly
    /// solve every eroded load and its complete envelope derivative. Respect
    /// physical regions and the strict material cap without increasing tolerance.
    /// A rejected floor is failure of this path, not global infeasibility.
    ///
    /// All recovery work consumes the original physical budget, separately from
    /// AL callback counts. Errors leave cached optimizer results and nominal
    /// physics intact. Compare the returned envelope with best_feasible before
    /// choosing a result: recovery guarantees neither descent nor optimality.
    pub fn recover_feasible(&mut self) -> Result<EnvelopeDesignPoint3, EvaluationStop> {
        use crate::sdf3::robust::{RobustDesignOptions3,
            continuation::controlled_sdf3_robust_feasible_start};
        self.control.checkpoint("sdf3-envelope-recovery-start")?;
        let restoration = RobustDesignOptions3 {
            robustness: self.options.robustness, volume_limit: self.options.volume_limit,
            density_floor: self.options.density_floor, compliance_scale: 1.0,
            optimizer: self.options.optimizer,
        };
        let rho = controlled_sdf3_robust_feasible_start(self.study, &self.accepted.rho,
            restoration, self.control)?;
        let evaluation = self.study.evaluate_robust_envelope(&rho, self.forces,
            self.options.robustness, &self.policy, self.control)?;
        let volume = evaluation.physical.dilated_volume_fraction;
        if !volume.is_finite() || volume > self.options.volume_limit {
            return Err(EvaluationStop::Breakdown { stage: "sdf3-envelope-recovery-volume" });
        }
        self.control.checkpoint("sdf3-envelope-recovery-publish")?;
        Ok(EnvelopeDesignPoint3 { rho, evaluation })
    }

    /// Up to additional_steps accepted updates, preserving AL state and work.
    /// Zero steps uses cached fields. Rejected trial fields never replace an
    /// incumbent; errors preserve any step already accepted before a late poll.
    pub fn run(&mut self, additional_steps: usize) -> Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let result = {
                let ledger = RefCell::new(&mut *self.control);
                let study = &mut *self.study;
                let forces = self.forces;
                let options = self.options;
                let policy = &self.policy;
                self.state.try_run(remaining.min(1),
                    &mut |rho| sample(study, forces, rho, policy, options, &mut last, &mut **ledger.borrow_mut()),
                    |_| if ledger.borrow_mut().checkpoint("sdf3-envelope-design-control").is_ok() {
                        ControlFlow::Continue(())
                    } else { ControlFlow::Break(()) })
            };
            if self.state.work().iterations > before {
                let accepted = last.expect("accepted envelope step has complete fields");
                assert_eq!(accepted.rho, self.state.point(), "envelope fields detached from AL state");
                self.accepted = accepted;
                let next = row(&self.accepted, self.state.work().iterations, self.options);
                if next.feasible && self.best.as_ref().is_none_or(|best|
                    self.accepted.evaluation.envelope.value < best.evaluation.envelope.value) {
                    self.best = Some(self.accepted.clone());
                }
                self.history.push(next);
                self.control.checkpoint("sdf3-envelope-design-accepted").map_err(ProjectedAlError::Evaluation)?;
            }
            let report = result?;
            if report.stop != ProjectedAlStop::IterationLimit || remaining <= 1 { return Ok(report); }
            remaining -= 1;
        }
    }
}

#[cfg(test)]
#[path = "optimize_tests.rs"]
mod tests;
