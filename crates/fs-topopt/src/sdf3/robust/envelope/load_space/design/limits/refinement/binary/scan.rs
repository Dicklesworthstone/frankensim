//! Select a binary endpoint only after the full load set passes on BOTH grids.
use std::ops::ControlFlow;
use super::{binary_error, BinaryLoadSpaceCapRefinement3};
use super::super::{invalid, LoadSpaceCap3};
use crate::sdf3::{AdaptiveSdf3Elasticity, CutDensityStudy3};
use crate::sdf3::robust::envelope::load_space::design::limits::binary::{
    BinaryLoadSpaceScanOptions3, BinaryLoadSpaceScanStop3,
};
use crate::sdf3_goal::{GoalRefinementError3, GoalRefinementOptions3};
use crate::{SimpParams, SolveControl, SolveWork};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, surface::ReferenceLoad3};

#[derive(Debug, Clone, Copy)]
pub struct BinaryRefinedScanOptions3 {
    pub selection: BinaryLoadSpaceScanOptions3,
    pub enrichment: GoalRefinementOptions3,
    /// Full normalized Gram change, not merely a maximum or diagonal change.
    pub change_tolerance: f64,
}
impl Default for BinaryRefinedScanOptions3 {
    fn default() -> Self {
        Self { selection: BinaryLoadSpaceScanOptions3::default(),
            enrichment: GoalRefinementOptions3::default(), change_tolerance: 0.02 }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct BinaryRefinedCandidate3 {
    pub index: usize,
    pub threshold: f64,
    /// Coarse numerical material measure; NOT a two-grid volume-error estimate.
    pub dilated_volume: f64,
    pub coarse_worst_compliance: f64,
    pub fine_worst_compliance: f64,
    pub normalized_gram_change: f64,
    pub two_grid_passed: bool,
    pub feasible: bool,
}
#[derive(Debug, Clone)]
pub struct BinaryRefinedScan3 {
    pub cap: LoadSpaceCap3,
    pub params: SimpParams,
    pub options: BinaryRefinedScanOptions3,
    /// Only complete two-grid experiments, in authored order.
    pub candidates: Vec<BinaryRefinedCandidate3>,
    /// Least coarse dilated volume among endpoints passing ALL requested gates.
    /// Contains masks and independently recomputed worst fields on both grids.
    pub best: Option<BinaryLoadSpaceCapRefinement3>,
    pub selected_index: Option<usize>,
    pub attempted_candidates: usize,
    pub stop: BinaryLoadSpaceScanStop3,
    pub error: Option<GoalRefinementError3>,
    pub work: SolveWork,
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Search an authored threshold schedule with the existing binary two-grid
    /// evaluator. Coarse feasibility alone NEVER admits a candidate. Passing
    /// both physical caps and the full-Gram change criterion is also required,
    /// together with the independently authored coarse dilated-volume limit.
    ///
    /// The builder supplies the SAME enriched background for each candidate,
    /// honors its checkpoint, and owns ONE cumulative geometry budget including
    /// failed builds. Active leaf identities are checked for equality; SDF,
    /// material, integration policy and support identity remain its obligation.
    /// The existing checked transfer refuses an unchanged or incompatible space.
    ///
    /// Each threshold is a fresh full two-grid experiment (including its source
    /// filter); no claim of one filter/build for this scan. Fine stiffness is
    /// inherited from that candidate's binary eroded mask, never refiltered raw
    /// controls. No new solver, spectrum algorithm or direction sampling.
    ///
    /// Late stops retain completed candidates and their best checked endpoint;
    /// unsuccessful candidates are not silently skipped. Equal-volume ties keep
    /// the first entry. Only the current/best full field families are retained.
    /// A complete scan with best=None is not success or global infeasibility.
    /// All solves share control and the borrowed study keeps its incoming state.
    /// This does not certify zero-void behavior, continuum error or manufacture.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_binary_load_space_cap_enrichment(&mut self, rho: &[f64],
        laws: &[ReferenceLoad3<'_>], params: SimpParams, cap: &LoadSpaceCap3,
        thresholds: &[f64], options: BinaryRefinedScanOptions3, control: &mut SolveControl<'_>,
        mut build: impl FnMut(&mut dyn FnMut() -> ControlFlow<()>) -> Result<AdaptiveElasticity3, GoalRefinementError3>,
    ) -> Result<BinaryRefinedScan3, GoalRefinementError3> {
        self.require_homogeneous_load_space().map_err(invalid)?;
        options.selection.validate(params, thresholds).map_err(binary_error)?;
        options.selection.robustness.assert_valid(params);
        let goal = options.enrichment;
        if !options.change_tolerance.is_finite() || options.change_tolerance < 0.0
            || laws.len() != cap.policy().dimension() || laws.len() > options.selection.max_load_cases
            || laws.len() > goal.max_load_cases || goal.max_transfer_terms == 0
            || ![goal.numerical.residual_tolerance,goal.numerical.identity_tolerance]
                .iter().all(|t| t.is_finite() && *t > 0.0 && *t < 1.0)
            || rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(invalid("invalid refined binary schedule, load basis, raw controls or policy"));
        }
        let mut report = BinaryRefinedScan3 { cap:cap.clone(), params, options,
            candidates:Vec::new(), best:None, selected_index:None, attempted_candidates:0,
            stop:BinaryLoadSpaceScanStop3::ScheduleComplete, error:None, work:control.work() };
        let attempt = (|| -> Result<(), GoalRefinementError3> {
            let mut fine_leaves: Option<Vec<fs_cutfem::octree3::Octant3>> = None;
            for (index,&threshold) in thresholds.iter().enumerate() {
                control.checkpoint("sdf3-binary-refined-scan-candidate")?;
                report.attempted_candidates += 1;
                let mut stopped = None;
                let built = build(&mut || match control.checkpoint("sdf3-binary-refined-scan-build") {
                    Ok(()) => ControlFlow::Continue(()),
                    Err(e) => { stopped = Some(e); ControlFlow::Break(()) }
                });
                if let Some(e) = stopped { return Err(e.into()); }
                control.checkpoint("sdf3-binary-refined-scan-build")?;
                let fine = built?;
                if let Some(leaves) = &fine_leaves {
                    if fine.leaves() != leaves.as_slice() {
                        return Err(invalid("refined binary scan builder changed its active background"));
                    }
                } else { fine_leaves = Some(fine.leaves().to_vec()); }
                let checked = self.estimate_binary_load_space_cap_enrichment(fine,laws,rho,params,cap,
                    options.selection.endpoint(threshold),goal,control)?;
                let two_grid_passed = checked.stable_and_feasible(options.selection.smooth_tolerance,options.change_tolerance);
                let volume = checked.physical.dilated.volume_fraction;
                let feasible = two_grid_passed && volume <= options.selection.dilated_volume_limit;
                report.candidates.push(BinaryRefinedCandidate3 { index,threshold,dilated_volume:volume,
                    coarse_worst_compliance:checked.coarse.worst_compliance,
                    fine_worst_compliance:checked.fine.worst_compliance,
                    normalized_gram_change:checked.normalized_gram_change,two_grid_passed,feasible });
                if feasible && report.best.as_ref().is_none_or(|best| volume < best.physical.dilated.volume_fraction) {
                    report.best = Some(checked); report.selected_index = Some(index);
                }
                // Publish the complete paired result before observing a late stop.
                control.checkpoint("sdf3-binary-refined-scan-accepted")?;
            }
            Ok(())
        })();
        if let Err(e) = attempt {
            report.stop = BinaryLoadSpaceScanStop3::EvaluationStopped; report.error = Some(e);
        }
        report.work = control.work();
        Ok(report)
    }

    /// Exhaustively search every distinct three-field binary mask state generated
    /// by the retained gray design, then require the SAME two-grid gates as the
    /// authored-schedule scan. A bounded source-projection prepass derives the
    /// complete threshold set; if it exceeds max_candidates the call refuses
    /// before any enriched geometry build or reference-load callback.
    ///
    /// The prepass itself is charged to control and temporarily installs the
    /// explicit nominal params only while deriving projected fields. Candidate
    /// evaluation still uses the existing fresh two-grid evaluator, so no
    /// stiffness, displacement, spectrum or feasibility is inferred from the
    /// prepass. This intentionally favors correctness over filter/build reuse.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_all_binary_load_space_cap_enrichment(&mut self, rho: &[f64],
        laws: &[ReferenceLoad3<'_>], params: SimpParams, cap: &LoadSpaceCap3,
        options: BinaryRefinedScanOptions3, control: &mut SolveControl<'_>,
        build: impl FnMut(&mut dyn FnMut() -> ControlFlow<()>) -> Result<AdaptiveElasticity3, GoalRefinementError3>,
    ) -> Result<BinaryRefinedScan3, GoalRefinementError3> {
        self.require_homogeneous_load_space().map_err(invalid)?;
        options.selection.validate(params, &[1.0]).map_err(binary_error)?;
        options.selection.robustness.assert_valid(params);
        let goal = options.enrichment;
        if !options.change_tolerance.is_finite() || options.change_tolerance < 0.0
            || laws.len() != cap.policy().dimension() || laws.len() > options.selection.max_load_cases
            || laws.len() > goal.max_load_cases || goal.max_transfer_terms == 0
            || ![goal.numerical.residual_tolerance,goal.numerical.identity_tolerance]
                .iter().all(|t| t.is_finite() && *t > 0.0 && *t < 1.0)
            || rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(invalid("invalid exhaustive refined binary model, load basis, raw controls or policy"));
        }
        let incoming = self.params;
        self.params = params;
        let source = self.binary_cap_source(rho, options.selection.robustness, control);
        self.params = incoming;
        let source = source?;
        let thresholds = crate::sdf3::robust::envelope::limits::binary::critical_thresholds(
            &source, options.selection.max_candidates).map_err(binary_error)?;
        self.scan_binary_load_space_cap_enrichment(
            rho, laws, params, cap, &thresholds, options, control, build)
    }

}

#[cfg(test)]
#[path = "scan_tests.rs"]
mod tests;
