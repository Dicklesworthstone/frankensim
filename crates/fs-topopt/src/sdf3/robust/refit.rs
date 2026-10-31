//! Transfer raw controls to a checked refined space, then solve a NEW robust
//! design problem. The immutable source and its solved fields are never edited.
use std::ops::ControlFlow;
use fs_ascent::projected_al::ProjectedAlError;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use super::continuation::{RobustContinuationOptions3, RobustContinuationReport3,
    controlled_sdf3_robust_continuation};
use super::super::{AdaptiveSdf3Elasticity, CutDensityStudy3};
use crate::pipeline::LoadCase;
use crate::sdf3_goal::{GoalReferenceLoad3, GoalRefinementError3};
use crate::{EvaluationStop, SimpParams, SolveControl};

/// One fresh optimizer per material model; physical work is never reset.
#[derive(Debug, Clone, Copy)]
pub struct RobustRefitOptions3 {
    pub continuation: RobustContinuationOptions3,
    /// Total sparse Q1 transfer coefficients, not coefficients per node.
    pub max_transfer_terms: usize,
    pub max_load_cases: usize,
}
impl Default for RobustRefitOptions3 {
    fn default() -> Self {
        Self { continuation: RobustContinuationOptions3::default(),
            max_transfer_terms: 2_000_000, max_load_cases: 64 }
    }
}

#[derive(Debug, Clone)]
pub enum RobustRefitError3 {
    Goal(GoalRefinementError3),
    Design(ProjectedAlError<EvaluationStop>),
}
impl From<GoalRefinementError3> for RobustRefitError3 {
    fn from(e: GoalRefinementError3) -> Self { Self::Goal(e) }
}
impl From<EvaluationStop> for RobustRefitError3 {
    fn from(e: EvaluationStop) -> Self { Self::Goal(e.into()) }
}
impl From<ProjectedAlError<EvaluationStop>> for RobustRefitError3 {
    fn from(e: ProjectedAlError<EvaluationStop>) -> Self { Self::Design(e) }
}
impl std::fmt::Display for RobustRefitError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "robust reference refit stopped: {self:?}")
    }
}
impl std::error::Error for RobustRefitError3 {}

/// Uninstalled candidate. A stopped report is NOT permission to replace the
/// source. `report.last`, when present, belongs only to this candidate grid.
pub struct RobustReferenceRefit3<O: AdaptiveSdf3Elasticity> {
    pub study: CutDensityStudy3<O>,
    /// Parent-inherited RAW values before fine filtering or cap restoration.
    pub inherited_rho: Vec<f64>,
    /// Includes the restored baseline, numerical gradient gate, accepted
    /// trajectory and best feasible point. A late stop preserves its prefix.
    pub report: RobustContinuationReport3,
}

fn invalid(message: &'static str) -> RobustRefitError3 {
    GoalRefinementError3::Invalid(message).into()
}
pub(super) fn poll(control: &mut SolveControl<'_>, stage: &'static str) -> ControlFlow<()> {
    if control.checkpoint(stage).is_ok() { ControlFlow::Continue(()) }
    else { ControlFlow::Break(()) }
}

// Complete admission before user load callbacks or nodal load allocations.
fn admit(cells: usize, rho: &[f64], loads: &[GoalReferenceLoad3<'_>],
    schedule: &[SimpParams], options: RobustRefitOptions3) -> Result<(), RobustRefitError3> {
    let d = options.continuation.design;
    d.optimizer.validate::<EvaluationStop>(cells)?;
    options.continuation.gradient.assert_valid();
    if schedule.is_empty() || options.max_transfer_terms == 0 || loads.is_empty()
        || loads.len() > options.max_load_cases || !loads.iter().any(|l| l.weight > 0.0)
        || loads.iter().any(|l| !l.weight.is_finite() || l.weight < 0.0)
        || !d.volume_limit.is_finite() || d.volume_limit <= 0.0 || d.volume_limit > 1.0
        || !d.compliance_scale.is_finite() || d.compliance_scale <= 0.0
        || !d.density_floor.is_finite() || d.density_floor <= 0.0 || d.density_floor >= 1.0
        || rho.len() != cells || rho.iter().any(|r| !r.is_finite() || *r < d.density_floor || *r > 1.0) {
        return Err(invalid("invalid robust refit policy, raw controls or reference loads"));
    }
    // Modeling panics, consistent with the existing evaluator, precede work.
    for &params in schedule { d.robustness.assert_valid(params); }
    Ok(())
}

pub(super) fn reference_design<O: AdaptiveSdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, loads: &[GoalReferenceLoad3<'_>], rho: &[f64],
    schedule: &[SimpParams], options: RobustRefitOptions3, control: &mut SolveControl<'_>,
) -> Result<RobustContinuationReport3, RobustRefitError3> {
    admit(study.cells(), rho, loads, schedule, options)?;
    let mut forces = Vec::with_capacity(loads.len());
    for law in loads {
        control.checkpoint("sdf3-robust-reference-load")?;
        forces.push(study.operator().adaptive().reference_load(law.load,
            || poll(control, "sdf3-robust-reference-load")).map_err(GoalRefinementError3::from)?);
    }
    let nodal: Vec<_> = forces.iter().zip(loads)
        .map(|(force, law)| LoadCase { force, weight: law.weight }).collect();
    Ok(controlled_sdf3_robust_continuation(study, &nodal, rho, schedule,
        options.continuation, control)?)
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Reintegrate the SAME independent reference laws and reoptimize on an
    /// explicitly constructed refined operator. No nodal force interpolation,
    /// prolonged displacement, old gradient or optimizer multiplier is reused.
    ///
    /// The source filter radius and physical solid/void labels are preserved.
    /// Only raw controls transfer through AdaptiveTransfer3's checked parents;
    /// any positive declared raw floor is respected (not silently raised to
    /// the nominal OC helper's 1e-3 floor). The candidate gets its own filter,
    /// numerical cut volumes, dilated-cap restoration and eroded equilibrium.
    ///
    /// `params` explicitly defines the new model, e.g. a continuation report's
    /// final params while the source retains its incoming nominal params.
    /// Input is raw controls, NOT a purported certified or cached solved point.
    /// The caller owns unchanged SDF, base material, supports and pure load laws;
    /// transfer additionally checks bounds/material/support compatibility.
    ///
    /// The returned study is never installed into the source. Read `report.stop`
    /// and `report.last` before publication; initialization failures and later
    /// stops cannot damage the original. All work shares the supplied control.
    /// There is no cross-grid descent, continuum or length-scale guarantee.
    #[allow(clippy::too_many_arguments)]
    pub fn refit_robust_reference<F: AdaptiveSdf3Elasticity>(
        &self, operator: F, rho: &[f64], params: SimpParams,
        loads: &[GoalReferenceLoad3<'_>], options: RobustRefitOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<RobustReferenceRefit3<F>, RobustRefitError3> {
        self.require_material_body_recipe("sdf3-robust-refit-requires-material-acceleration")?;
        admit(self.cells(), rho, loads, &[params], options)?;
        options.continuation.design.optimizer.validate::<EvaluationStop>(operator.cells())?;
        control.checkpoint("sdf3-robust-refit-start")?;
        let (inherited_rho, regions) = {
            let coarse = self.operator().adaptive();
            let fine = operator.adaptive();
            let transfer = AdaptiveTransfer3::new(coarse, fine, options.max_transfer_terms,
                || poll(control, "sdf3-robust-refit-transfer")).map_err(GoalRefinementError3::from)?;
            if !fine.leaves().iter().zip(transfer.parents())
                .any(|(leaf, &p)| leaf.level() > coarse.leaves()[p].level()) {
                return Err(invalid("robust refit requires actual active-space enrichment"));
            }
            let mut inherited = Vec::with_capacity(fine.cells());
            for &parent in transfer.parents() {
                control.checkpoint("sdf3-robust-refit-transfer")?;
                inherited.push(rho[parent]);
            }
            let regions = self.physical_regions().map(|regions|
                super::super::regions::inherit_physical_regions(regions, transfer.parents(), control)).transpose()?;
            (inherited, regions)
        };
        let mut study = CutDensityStudy3::new(operator, self.filter_radius, params);
        self.inherit_material_acceleration(&mut study, control)?;
        if let Some(regions) = regions { study = study.with_physical_regions(regions)?; }
        let report = reference_design(&mut study, loads, &inherited_rho, &[params], options, control)?;
        // Do not discard an admitted candidate prefix with a post-run poll.
        // Publication belongs to the caller and needs its own checkpoint.
        Ok(RobustReferenceRefit3 { study, inherited_rho, report })
    }
}

#[cfg(test)]
#[path = "refit_tests.rs"]
mod tests;
