//! A new cap-constrained design on a checked finer space, never an in-place
//! replacement of a source result. Caps, load order and smoothing stay fixed.
use std::ops::ControlFlow;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport, ProjectedAlStop};
use fs_cutfem::elastic3::{
    adaptive::enrichment::AdaptiveTransfer3, surface::ReferenceLoad3,
};
use super::{ComplianceCaps3, ComplianceLimitedDesign3, ComplianceLimitedIteration3,
    ComplianceLimitedOptions3, EnvelopeDesignPoint3};
use super::start::{ComplianceStartOptions3, ComplianceStartReport3, ComplianceStartStop3,
    controlled_compliance_feasible_start};
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3};
use crate::sdf3_goal::GoalRefinementError3;
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};

#[derive(Debug, Clone, Copy)]
pub struct ComplianceCapRefitOptions3 {
    pub design: ComplianceLimitedOptions3,
    pub start: ComplianceStartOptions3,
    pub accepted_steps: usize,
    pub max_transfer_terms: usize,
    pub max_load_cases: usize,
}
impl Default for ComplianceCapRefitOptions3 {
    fn default() -> Self {
        Self { design: ComplianceLimitedOptions3::default(), start: ComplianceStartOptions3::default(),
            accepted_steps: 100, max_transfer_terms: 2_000_000, max_load_cases: 64 }
    }
}
#[derive(Debug, Clone)]
pub enum ComplianceCapRefitError3 {
    Goal(GoalRefinementError3),
    Design(ProjectedAlError<EvaluationStop>),
}
impl From<GoalRefinementError3> for ComplianceCapRefitError3 {
    fn from(e: GoalRefinementError3) -> Self { Self::Goal(e) }
}
impl From<EvaluationStop> for ComplianceCapRefitError3 {
    fn from(e: EvaluationStop) -> Self { Self::Goal(e.into()) }
}
impl From<ProjectedAlError<EvaluationStop>> for ComplianceCapRefitError3 {
    fn from(e: ProjectedAlError<EvaluationStop>) -> Self { Self::Design(e) }
}
impl std::fmt::Display for ComplianceCapRefitError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "cap-constrained refit: {self:?}") }
}
impl std::error::Error for ComplianceCapRefitError3 {}

/// Candidate fields belong to study, NEVER to the source coarse grid.
/// A restoration seed is not an optimized design. A final AL status is not
/// itself a cap-feasibility, continuum, or cross-grid-descent certificate.
pub struct ComplianceCapRefit3<O: AdaptiveSdf3Elasticity> {
    pub study: CutDensityStudy3<O>,
    pub params: SimpParams,
    pub caps: ComplianceCaps3,
    /// Raw input on an initial grid, or checked parent-inherited input on a refit.
    pub inherited_rho: Vec<f64>,
    pub restoration: ComplianceStartReport3,
    pub history: Vec<ComplianceLimitedIteration3>,
    pub best: Option<EnvelopeDesignPoint3>,
    /// None means restoration stopped before optimizer initialization.
    pub outcome: Option<Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>>>,
    /// Fresh best-design re-solve after a completed optimizer allowance.
    /// Late failure retains the preceding complete best point above.
    pub recheck: Option<Result<EnvelopeDesignPoint3, EvaluationStop>>,
    pub options: ComplianceCapRefitOptions3,
    pub work: SolveWork,
}
impl<O: AdaptiveSdf3Elasticity> ComplianceCapRefit3<O> {
    /// A bounded optimizer allowance completed and every cap passed a fresh
    /// re-solve. This is NOT a stationarity or continuum guarantee.
    #[must_use]
    pub fn completed_and_rechecked(&self) -> bool {
        matches!(&self.outcome, Some(Ok(r))
            if matches!(r.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit))
            && self.recheck.as_ref().is_some_and(|r| r.as_ref().is_ok_and(|p|
                self.caps.admits(&p.evaluation, self.options.design.optimizer.tolerance)))
    }
}

fn poll(control: &mut SolveControl<'_>, stage: &'static str, stop: &mut Option<EvaluationStop>)
    -> ControlFlow<()> {
    match control.checkpoint(stage) {
        Ok(()) => ControlFlow::Continue(()),
        Err(error) => { *stop = Some(error); ControlFlow::Break(()) }
    }
}
fn invalid(message: &'static str) -> ComplianceCapRefitError3 {
    GoalRefinementError3::Invalid(message).into()
}

fn admit_grid(cells: usize, rho: &[f64], params: SimpParams, laws: &[ReferenceLoad3<'_>],
    caps: &ComplianceCaps3, options: ComplianceCapRefitOptions3) -> Result<(), ComplianceCapRefitError3> {
    let design = options.design;
    design.optimizer.validate::<EvaluationStop>(cells)?;
    design.robustness.assert_valid(params);
    design.gradient_check.assert_valid();
    if laws.len() != caps.caps().len() || laws.len() > options.max_load_cases
        || options.max_transfer_terms == 0 || !design.density_floor.is_finite()
        || design.density_floor <= 0.0 || design.density_floor >= 1.0
        || rho.len() != cells || rho.iter().any(|r| !r.is_finite() || *r < design.density_floor || *r > 1.0)
        || options.start.max_bisections > 64 || !options.start.margin.is_finite()
        || options.start.margin < 0.0 || options.start.margin >= 1.0 - caps.envelope().bias_bound() {
        return Err(invalid("invalid cap refit inputs, raw bounds or resource policy"));
    }
    Ok(())
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Parent-inherit RAW controls and authored solid/void regions, rebuild
    /// filtering and numerical volume, reintegrate every reference law, recover
    /// a cap-feasible seed and minimize dilated material on the candidate.
    /// The previous optimizer, displacement fields and gradients are not reused.
    ///
    /// This consumes a genuinely finer operator supplied by the caller, who
    /// owns unchanged SDF/material/support/integration identity and its geometry
    /// build budget. The sparse transfer additionally checks compatibility.
    /// Physical caps and temperature are immutable across this operation.
    ///
    /// The source is borrowed immutably. No candidate is installed by this API,
    /// including after a successful run. Failed initialization preserves any
    /// already solved restoration seed; later errors retain the optimizer's
    /// accepted prefix. All solves and the final recheck share control.
    #[allow(clippy::too_many_arguments)]
    pub fn refit_compliance_caps<F: AdaptiveSdf3Elasticity>(
        &self, operator: F, rho: &[f64], params: SimpParams, laws: &[ReferenceLoad3<'_>],
        caps: &ComplianceCaps3, options: ComplianceCapRefitOptions3, control: &mut SolveControl<'_>,
    ) -> Result<ComplianceCapRefit3<F>, ComplianceCapRefitError3> {
        if self.has_material_body_load() {
            return Err(invalid("material body-load refits require a load prepared on the candidate geometry"));
        }
        admit_grid(self.cells(), rho, params, laws, caps, options)?;
        options.design.optimizer.validate::<EvaluationStop>(operator.cells())?;
        control.checkpoint("sdf3-cap-refit-start")?;
        let mut stop = None;
        let mapped = AdaptiveTransfer3::new(self.operator().adaptive(), operator.adaptive(),
            options.max_transfer_terms, || poll(control, "sdf3-cap-refit-transfer", &mut stop));
        if let Some(error) = stop { return Err(error.into()); }
        let transfer = mapped.map_err(GoalRefinementError3::from)?;
        if !operator.adaptive().leaves().iter().zip(transfer.parents()).any(|(leaf, &parent)|
            leaf.level() > self.operator().adaptive().leaves()[parent].level()) {
            return Err(invalid("cap refit requires actual active-space enrichment"));
        }
        let mut inherited_rho = Vec::with_capacity(operator.cells());
        for &parent in transfer.parents() {
            control.checkpoint("sdf3-cap-refit-transfer")?;
            inherited_rho.push(rho[parent]);
        }
        let regions = self.physical_regions().map(|regions|
            crate::sdf3::regions::inherit_physical_regions(regions, transfer.parents(), control)).transpose()?;
        drop(transfer);
        let mut study = CutDensityStudy3::new(operator, self.filter_radius, params);
        if let Some(regions) = regions { study = study.with_physical_regions(regions)?; }
        study.solve_reference_compliance_caps(&inherited_rho, laws, caps, options, control)
    }

    /// Solve the initial grid with the SAME reference-load integration,
    /// feasible-start recovery, gradient gate, optimizer and fresh recheck as a
    /// refined candidate. No adaptive transfer or geometry rebuilding occurs.
    ///
    /// This consumes the working study. The returned owner retains it together
    /// with the restoration seed and accepted prefix even when later numerical
    /// work stops. Shape/policy and load-assembly errors precede any admitted
    /// design. `inherited_rho` records the supplied initial controls on this
    /// route; for refits it records their parent-inherited counterpart.
    /// A completed allowance is not an accuracy or stationarity certificate.
    pub fn solve_reference_compliance_caps(
        self, incoming: &[f64], laws: &[ReferenceLoad3<'_>], caps: &ComplianceCaps3,
        options: ComplianceCapRefitOptions3, control: &mut SolveControl<'_>,
    ) -> Result<ComplianceCapRefit3<O>, ComplianceCapRefitError3> {
        let params = self.params();
        admit_grid(self.cells(), incoming, params, laws, caps, options)?;
        let design = options.design;
        let inherited_rho = incoming.to_vec();
        let mut study = self;
        let mut forces = Vec::with_capacity(laws.len());
        for &law in laws {
            control.checkpoint("sdf3-cap-refit-load")?;
            let mut stop = None;
            let result = study.operator().adaptive().reference_load(law,
                || poll(control, "sdf3-cap-refit-load", &mut stop));
            if let Some(error) = stop { return Err(error.into()); }
            forces.push(result.map_err(GoalRefinementError3::from)?);
        }
        let borrowed: Vec<_> = forces.iter().map(Vec::as_slice).collect();
        let restoration = controlled_compliance_feasible_start(&mut study, &borrowed, &inherited_rho,
            caps, design, options.start, control)?;
        let (mut history, mut best, mut outcome, mut recheck) = (Vec::new(), None, None, None);
        if matches!(restoration.stop, ComplianceStartStop3::IncomingSatisfied | ComplianceStartStop3::Restored) {
            if let Some(seed) = &restoration.point {
                match ComplianceLimitedDesign3::new(&mut study, &borrowed, &seed.rho, caps.clone(), design, control) {
                    Err(error) => outcome = Some(Err(error)),
                    Ok(mut session) => {
                        let result = session.run(options.accepted_steps);
                        history = session.history().to_vec();
                        best = session.best_feasible().cloned();
                        if matches!(&result, Ok(r) if matches!(r.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit)) {
                            recheck = Some(session.recheck_best().and_then(|point| point.ok_or(
                                EvaluationStop::Breakdown { stage: "sdf3-cap-refit-no-feasible-point" })));
                        }
                        outcome = Some(result);
                    }
                }
            }
        }
        // No fallible publication step may discard the already solved prefix.
        Ok(ComplianceCapRefit3 { study, params, caps: caps.clone(), inherited_rho, restoration,
            history, best, outcome, recheck, options, work: control.work() })
    }
}
