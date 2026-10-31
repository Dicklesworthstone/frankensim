//! Refit a complete simultaneous-load requirement on an admitted finer space.
//! Transfer RAW controls; inherited-stiffness assessment is a different experiment.
use std::ops::ControlFlow;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport, ProjectedAlStop};
use fs_cutfem::elastic3::{adaptive::enrichment::AdaptiveTransfer3, surface::ReferenceLoad3};
use super::{LoadSpaceCap3, LoadSpaceLimitedDesign3};
use super::super::{LoadSpaceDesignIteration3, LoadSpaceDesignPoint3};
use super::start::{LoadSpaceStartReport3, LoadSpaceStartStop3, controlled_load_space_feasible_start};
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3};
use crate::sdf3_goal::GoalRefinementError3;
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};

// Same policy units and error owners as finite-case refitting. The constraint
// is NOT the finite-case envelope: every evaluation uses the complete spectrum.
pub use crate::sdf3::robust::envelope::limits::refit::{
    ComplianceCapRefitOptions3 as LoadSpaceCapRefitOptions3,
    ComplianceCapRefitError3 as LoadSpaceCapRefitError3,
};

/// Uninstalled candidate, with all fields bound to its own grid and fixed cap.
/// Restoration seeds, accepted optimizer prefixes and fresh checks are separate.
pub struct LoadSpaceCapRefit3<O: AdaptiveSdf3Elasticity> {
    pub study: CutDensityStudy3<O>,
    pub params: SimpParams,
    pub cap: LoadSpaceCap3,
    /// Supplied initial controls, or checked parent-inherited controls on a refit.
    pub inherited_rho: Vec<f64>,
    pub restoration: LoadSpaceStartReport3,
    pub history: Vec<LoadSpaceDesignIteration3>,
    pub best: Option<LoadSpaceDesignPoint3>,
    /// None means restoration did not admit a complete optimizer baseline.
    pub outcome: Option<Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>>>,
    /// Fresh COMPLETE-basis solve; never just the previously worst direction.
    pub recheck: Option<Result<LoadSpaceDesignPoint3, EvaluationStop>>,
    pub options: LoadSpaceCapRefitOptions3,
    pub work: SolveWork,
}
impl<O: AdaptiveSdf3Elasticity> LoadSpaceCapRefit3<O> {
    /// Bounded optimization completed and a fresh solve met the full load-set
    /// requirement. Not stationarity, cross-grid descent or continuum accuracy.
    #[must_use]
    pub fn completed_and_rechecked(&self) -> bool {
        matches!(&self.outcome, Some(Ok(r))
            if matches!(r.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit))
            && self.recheck.as_ref().is_some_and(|r| r.as_ref().is_ok_and(|p|
                self.cap.admits(&p.evaluation, self.options.design.optimizer.tolerance)))
    }
}

fn invalid(message: &'static str) -> LoadSpaceCapRefitError3 {
    GoalRefinementError3::Invalid(message).into()
}
fn poll(control: &mut SolveControl<'_>, stage: &'static str, stop: &mut Option<EvaluationStop>)
    -> ControlFlow<()> {
    match control.checkpoint(stage) {
        Ok(()) => ControlFlow::Continue(()),
        Err(error) => { *stop = Some(error); ControlFlow::Break(()) }
    }
}
fn admit(cells: usize, rho: &[f64], params: SimpParams, laws: &[ReferenceLoad3<'_>],
    cap: &LoadSpaceCap3, options: LoadSpaceCapRefitOptions3) -> Result<(), LoadSpaceCapRefitError3> {
    let design = options.design;
    design.optimizer.validate::<EvaluationStop>(cells)?;
    design.robustness.assert_valid(params);
    design.gradient_check.assert_valid();
    if laws.len() != cap.policy().dimension() || laws.len() > options.max_load_cases
        || options.max_transfer_terms == 0 || !design.density_floor.is_finite()
        || design.density_floor <= 0.0 || design.density_floor >= 1.0
        || rho.len() != cells || rho.iter().any(|r| !r.is_finite() || *r < design.density_floor || *r > 1.0)
        || options.start.max_bisections > 64 || !options.start.margin.is_finite()
        || options.start.margin < 0.0 || options.start.margin >= 1.0-cap.policy().bias_bound() {
        return Err(invalid("invalid combined-load refit inputs, raw bounds or resource policy"));
    }
    Ok(())
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Build a new minimum-material problem on the supplied finer operator.
    /// Parent-inherit raw controls and prescribed physical labels, preserve the
    /// source filter radius, reintegrate every reference load, restore the SAME
    /// combined-load cap, check gradients and run a fresh constrained optimizer.
    ///
    /// Load amplitudes/correlations/order and the cap/temperature remain fixed.
    /// No old optimizer state, force vector, gradient, displacement, or worst
    /// direction is transferred. Even a dependent basis is assessed as authored.
    ///
    /// The source is immutable and this API never installs the candidate. The
    /// caller must preserve SDF/material/support/integration identity and own the
    /// geometry budget; checked transfer additionally enforces compatibility.
    /// Preparation errors cannot damage the source. After restoration, a late
    /// stop retains complete seeds/accepted fields. All work shares control.
    #[allow(clippy::too_many_arguments)]
    pub fn refit_load_space_cap<F: AdaptiveSdf3Elasticity>(
        &self, operator: F, rho: &[f64], params: SimpParams, laws: &[ReferenceLoad3<'_>],
        cap: &LoadSpaceCap3, options: LoadSpaceCapRefitOptions3, control: &mut SolveControl<'_>,
    ) -> Result<LoadSpaceCapRefit3<F>, LoadSpaceCapRefitError3> {
        self.require_homogeneous_load_space().map_err(invalid)?;
        admit(self.cells(), rho, params, laws, cap, options)?;
        options.design.optimizer.validate::<EvaluationStop>(operator.cells())?;
        control.checkpoint("sdf3-load-space-refit-start")?;
        let mut stop = None;
        let mapped = AdaptiveTransfer3::new(self.operator().adaptive(), operator.adaptive(),
            options.max_transfer_terms, || poll(control, "sdf3-load-space-refit-transfer", &mut stop));
        if let Some(error) = stop { return Err(error.into()); }
        let transfer = mapped.map_err(GoalRefinementError3::from)?;
        if !operator.adaptive().leaves().iter().zip(transfer.parents()).any(|(leaf, &parent)|
            leaf.level() > self.operator().adaptive().leaves()[parent].level()) {
            return Err(invalid("combined-load refit requires actual active-space enrichment"));
        }
        let mut inherited = Vec::with_capacity(operator.cells());
        for &parent in transfer.parents() {
            control.checkpoint("sdf3-load-space-refit-transfer")?;
            inherited.push(rho[parent]);
        }
        let regions = self.physical_regions().map(|regions|
            crate::sdf3::regions::inherit_physical_regions(regions, transfer.parents(), control)).transpose()?;
        drop(transfer);
        let mut study = CutDensityStudy3::new(operator, self.filter_radius, params);
        if let Some(regions) = regions { study = study.with_physical_regions(regions)?; }
        study.solve_reference_load_space_cap(&inherited, laws, cap, options, control)
    }

    /// Initial-grid solve using EXACTLY the same reference-load integration,
    /// full-spectrum restoration, gradient gate, optimizer and recheck as refits.
    /// Consume the working study and return it paired with any admitted prefix.
    /// Resource exhaustion never becomes a completed/rechecked result.
    pub fn solve_reference_load_space_cap(
        self, incoming: &[f64], laws: &[ReferenceLoad3<'_>], cap: &LoadSpaceCap3,
        options: LoadSpaceCapRefitOptions3, control: &mut SolveControl<'_>,
    ) -> Result<LoadSpaceCapRefit3<O>, LoadSpaceCapRefitError3> {
        self.require_homogeneous_load_space().map_err(invalid)?;
        let params = self.params();
        admit(self.cells(), incoming, params, laws, cap, options)?;
        let mut study = self;
        let mut forces = Vec::with_capacity(laws.len());
        for &law in laws {
            control.checkpoint("sdf3-load-space-refit-load")?;
            let mut stop = None;
            let assembled = study.operator().adaptive().reference_load(law,
                || poll(control, "sdf3-load-space-refit-load", &mut stop));
            if let Some(error) = stop { return Err(error.into()); }
            forces.push(assembled.map_err(GoalRefinementError3::from)?);
        }
        let basis: Vec<_> = forces.iter().map(Vec::as_slice).collect();
        let restoration = controlled_load_space_feasible_start(&mut study, &basis, incoming,
            cap, options.design, options.start, control)?;
        let (mut history, mut best, mut outcome, mut recheck) = (Vec::new(), None, None, None);
        if matches!(restoration.stop, LoadSpaceStartStop3::IncomingSatisfied | LoadSpaceStartStop3::Restored) {
            if let Some(seed) = &restoration.point {
                match LoadSpaceLimitedDesign3::new(&mut study, &basis, &seed.rho, cap.clone(), options.design, control) {
                    Err(error) => outcome = Some(Err(error)),
                    Ok(mut session) => {
                        let result = session.run(options.accepted_steps);
                        history = session.history().to_vec();
                        best = session.best_feasible().cloned();
                        if matches!(&result, Ok(r) if matches!(r.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit)) {
                            recheck = Some(session.recheck_best().and_then(|point| point.ok_or(
                                EvaluationStop::Breakdown { stage: "sdf3-load-space-refit-no-feasible-point" })));
                        }
                        outcome = Some(result);
                    }
                }
            }
        }
        Ok(LoadSpaceCapRefit3 { study, params, cap: cap.clone(), inherited_rho: incoming.to_vec(),
            restoration, history, best, outcome, recheck, options, work: control.work() })
    }
}

#[cfg(test)]
#[path = "refit_tests.rs"]
mod tests;
