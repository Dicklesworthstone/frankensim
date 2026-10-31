//! Eroded goal assessment -> selective background refinement -> robust refit.
//! The probe inherits physical stiffness; the optimizer inherits RAW controls.
use std::ops::ControlFlow;
use fs_cutfem::octree3::{Octree3, OctreeError3};
use fs_dwr::elasticity3::GoalMarking3;
use super::{RobustDesignPoint3, continuation::{RobustContinuationReport3, RobustContinuationStop3}};
use super::refit::{RobustRefitError3, RobustRefitOptions3, poll, reference_design};
use super::super::{AdaptiveSdf3Elasticity, CutDensityStudy3};
use crate::sdf3_goal::{ComplianceRefinement3, GoalReferenceLoad3, GoalRefinementError3, GoalRefinementOptions3};
use crate::{EvaluationStop, SimpParams, SolveControl, SolveWork};

#[derive(Debug, Clone, Copy)]
pub struct RobustAdaptiveOptions3 {
    pub refit: RobustRefitOptions3,
    pub enrichment: GoalRefinementOptions3,
    pub max_refinements: usize,
    pub marking_fraction: f64,
    pub max_marks: usize,
    /// Optional maximum relative coarse/fine compliance difference for EACH
    /// positive-weight load. No cancellation across loads can pass this gate.
    /// This is a two-grid numerical observation, not a continuum certificate.
    /// When present, the final allowed grid is probed before a cap stop.
    pub eroded_relative_tolerance: Option<f64>,
}
impl Default for RobustAdaptiveOptions3 {
    fn default() -> Self {
        Self { refit: RobustRefitOptions3::default(), enrichment: GoalRefinementOptions3::default(),
            max_refinements: 2, marking_fraction: 0.5, max_marks: 2, eroded_relative_tolerance: None }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobustAdaptiveStop3 {
    /// Requested count reached; not a numerical accuracy or convergence pass.
    RefinementLimit,
    /// Every positive-weight eroded compliance passed on the CURRENT grid.
    ErodedComplianceTolerance,
    /// Initial or proposed material optimization did not finish its allowance.
    DesignStopped,
    /// Empty residual marks do not establish accuracy.
    NoRefinementSignal,
    RefinementStopped,
}
#[derive(Debug, Clone)]
pub enum RobustAdaptiveError3 {
    Background(OctreeError3),
    Refit(RobustRefitError3),
}
impl From<RobustRefitError3> for RobustAdaptiveError3 {
    fn from(e: RobustRefitError3) -> Self { Self::Refit(e) }
}
impl From<GoalRefinementError3> for RobustAdaptiveError3 {
    fn from(e: GoalRefinementError3) -> Self { Self::Refit(e.into()) }
}
impl From<EvaluationStop> for RobustAdaptiveError3 {
    fn from(e: EvaluationStop) -> Self { Self::Refit(e.into()) }
}
impl From<OctreeError3> for RobustAdaptiveError3 {
    fn from(e: OctreeError3) -> Self {
        if e == OctreeError3::Cancelled { EvaluationStop::Cancelled.into() }
        else { Self::Background(e) }
    }
}
impl std::fmt::Display for RobustAdaptiveError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "robust adaptation stopped: {self:?}") }
}
impl std::error::Error for RobustAdaptiveError3 {}

#[derive(Debug, Clone)]
pub struct RobustGridDesign3 {
    pub background_cells: usize,
    pub active_cells: usize,
    /// Initial full schedule, or one fresh final-model refit on a later grid.
    pub design: RobustContinuationReport3,
}
#[derive(Debug, Clone)]
pub struct RobustGridProbe3 {
    /// Index in grids: the design ASSESSED, never a proposed fine-grid design.
    pub source_grid: usize,
    pub estimate: ComplianceRefinement3,
    pub maximum_relative_gap: f64,
    pub marking: Option<GoalMarking3>,
    pub target_background_cells: Option<usize>,
    /// Unfiltered parent-inherited controls, before candidate cap restoration.
    pub inherited_rho: Option<Vec<f64>>,
    pub installed: bool,
}
#[derive(Debug, Clone)]
pub struct RobustAdaptiveReport3 {
    /// Only installed grids; initial interrupted continuation retains its prefix.
    pub grids: Vec<RobustGridDesign3>,
    pub probes: Vec<RobustGridProbe3>,
    /// Rejected candidate prefix, if fitting reached that stage. These fields
    /// are NEVER associated with the installed study or the current background.
    pub rejected_design: Option<RobustContinuationReport3>,
    pub stop: RobustAdaptiveStop3,
    pub error: Option<RobustAdaptiveError3>,
    pub options: RobustAdaptiveOptions3,
    pub work: SolveWork,
}
impl RobustAdaptiveReport3 {
    #[must_use] pub fn last(&self) -> Option<&RobustDesignPoint3> {
        self.grids.last().and_then(|grid| grid.design.last.as_ref())
    }
}
fn check_tree<O: AdaptiveSdf3Elasticity>(operator: &O, tree: &Octree3) -> Result<(), GoalRefinementError3> {
    if operator.adaptive().leaves().iter().all(|cell| tree.leaves().contains(cell)) { Ok(()) }
    else { Err(GoalRefinementError3::Invalid("robust operator cells differ from the proposed background")) }
}
fn build_operator<O: AdaptiveSdf3Elasticity>(tree: &Octree3, control: &mut SolveControl<'_>,
    stage: &'static str,
    build: &mut impl FnMut(&Octree3, &mut dyn FnMut() -> ControlFlow<()>) -> Result<O, GoalRefinementError3>,
) -> Result<O, GoalRefinementError3> {
    control.checkpoint(stage)?;
    let mut stopped = false;
    let result = build(tree, &mut || {
        let flow = poll(control, stage); stopped |= flow.is_break(); flow
    });
    if stopped { return Err(EvaluationStop::Cancelled.into()); }
    control.checkpoint(stage)?;
    let operator = result?;
    check_tree(&operator, tree)?;
    Ok(operator)
}
fn relative_gap(a: f64, b: f64) -> f64 {
    let scale = a.abs().max(b.abs());
    if !a.is_finite() || !b.is_finite() { f64::INFINITY }
    else if scale == 0.0 { 0.0 } else { (a / scale - b / scale).abs() }
}
fn maximum_gap(estimate: &ComplianceRefinement3) -> f64 {
    estimate.cases.iter().zip(&estimate.weights).filter(|(_, w)| **w > 0.0)
        .map(|(case, _)| relative_gap(case.coarse_value, case.fine_value)).fold(0.0, f64::max)
}

/// Run the initial material schedule, then repeatedly assess and reoptimize the
/// final model on SELECTIVELY refined backgrounds. Existing CutFEM, DWR, raw
/// transfer, gradient audit and robust AL implementations do the numerical work.
///
/// The builder must preserve domain/SDF/material/support/integration identity
/// and enforce ONE cumulative geometry budget across all builds, including
/// failures. It receives cooperative cancellation; callbacks are indivisible.
/// Octree depth/leaf caps remain on the tree. All physical solves share control.
///
/// A globally enriched probe measures the current ERODED physical stiffness.
/// Its Dörfler marks refine the current original tree, NOT the probe tree. The
/// candidate refit rebuilds filtering, reference loads and dilated feasibility.
/// Its optimizer is fresh; no cross-grid compliance descent is asserted.
///
/// The tree, study and matching feasible report change together ONLY after a
/// completed candidate and a final checkpoint. Failed proposals preserve the
/// previous grid. Post-install cancellation preserves the complete new grid.
/// Operator scales remain nominal working state; use report.last()'s eroded
/// fields. A two-grid compliance pass does not assess dilated-volume error,
/// minimum feature size, continuum safety, or optimization stationarity.
#[allow(clippy::too_many_arguments)]
pub fn controlled_adaptive_robust_design<O: AdaptiveSdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, tree: &mut Octree3, loads: &[GoalReferenceLoad3<'_>],
    rho0: &[f64], schedule: &[SimpParams], options: RobustAdaptiveOptions3,
    control: &mut SolveControl<'_>,
    mut build: impl FnMut(&Octree3, &mut dyn FnMut() -> ControlFlow<()>) -> Result<O, GoalRefinementError3>,
) -> Result<RobustAdaptiveReport3, RobustAdaptiveError3> {
    study.require_material_body_recipe("sdf3-robust-adapt-requires-material-acceleration")?;
    if !options.marking_fraction.is_finite() || options.marking_fraction <= 0.0
        || options.marking_fraction > 1.0 || options.max_marks == 0
        || options.eroded_relative_tolerance.is_some_and(|t| !t.is_finite() || t <= 0.0 || t >= 1.0)
        || options.enrichment.max_transfer_terms == 0 || loads.len() > options.enrichment.max_load_cases
        || ![options.enrichment.numerical.residual_tolerance, options.enrichment.numerical.identity_tolerance]
            .iter().all(|t| t.is_finite() && *t > 0.0 && *t < 1.0) {
        return Err(GoalRefinementError3::Invalid("invalid robust adaptive policy").into());
    }
    check_tree(study.operator(), tree)?;
    let initial = reference_design(study, loads, rho0, schedule, options.refit, control)?;
    let complete = initial.stop == RobustContinuationStop3::ScheduleComplete && initial.last.is_some();
    let mut report = RobustAdaptiveReport3 {
        grids: vec![RobustGridDesign3 { background_cells: tree.leaves().len(), active_cells: study.cells(), design: initial }],
        probes: Vec::new(), rejected_design: None,
        stop: if complete { RobustAdaptiveStop3::RefinementLimit } else { RobustAdaptiveStop3::DesignStopped },
        error: None, options, work: control.work(),
    };
    if !complete { return Ok(report); }
    let mut installed = 0;
    loop {
        if installed == options.max_refinements && options.eroded_relative_tolerance.is_none() { break; }
        let attempt = (|| -> Result<bool, RobustAdaptiveError3> {
            let source_grid = report.grids.len() - 1;
            let current = &report.grids[source_grid].design;
            let params = current.params;
            let raw = &current.last.as_ref().expect("admitted grid has a design").rho;
            let probe_tree = tree.refined(study.operator().adaptive().leaves(),
                || poll(control, "sdf3-robust-adaptive-probe-tree"))?;
            let probe = build_operator(&probe_tree, control, "sdf3-robust-adaptive-probe-build", &mut build)?;
            let assessment = study.estimate_robust_reference_compliance_enrichment(probe.into_adaptive(), loads,
                raw, params, options.refit.continuation.design.robustness, options.enrichment, control)?;
            let gap = maximum_gap(&assessment.estimate);
            report.probes.push(RobustGridProbe3 { source_grid, estimate: assessment.estimate,
                maximum_relative_gap: gap, marking: None, target_background_cells: None,
                inherited_rho: None, installed: false });
            if options.eroded_relative_tolerance.is_some_and(|t| gap <= t) {
                report.stop = RobustAdaptiveStop3::ErodedComplianceTolerance;
                return Ok(false);
            }
            if installed == options.max_refinements { return Ok(false); }
            let diagnostic = report.probes.last_mut().expect("completed probe");
            let marks = diagnostic.estimate.mark(options.marking_fraction, options.max_marks,
                || poll(control, "sdf3-robust-adaptive-mark")).map_err(GoalRefinementError3::from)?;
            let empty = marks.marked.is_empty();
            diagnostic.marking = Some(marks);
            if empty { report.stop = RobustAdaptiveStop3::NoRefinementSignal; return Ok(false); }
            let next_tree = tree.refined(&diagnostic.marking.as_ref().expect("retained marks").marked,
                || poll(control, "sdf3-robust-adaptive-selected-tree"))?;
            diagnostic.target_background_cells = Some(next_tree.leaves().len());
            let operator = build_operator(&next_tree, control, "sdf3-robust-adaptive-candidate-build", &mut build)?;
            let candidate = study.refit_robust_reference(operator, raw, params, loads, options.refit, control)?;
            diagnostic.inherited_rho = Some(candidate.inherited_rho);
            if candidate.report.stop != RobustContinuationStop3::ScheduleComplete || candidate.report.last.is_none() {
                report.rejected_design = Some(candidate.report);
                report.stop = RobustAdaptiveStop3::DesignStopped;
                return Ok(false);
            }
            if let Err(error) = control.checkpoint("sdf3-robust-adaptive-promote") {
                report.rejected_design = Some(candidate.report);
                return Err(error.into());
            }
            // No fallible work between installing geometry and its matching result.
            *study = candidate.study;
            *tree = next_tree;
            diagnostic.installed = true;
            report.grids.push(RobustGridDesign3 { background_cells: tree.leaves().len(),
                active_cells: study.cells(), design: candidate.report });
            installed += 1;
            control.checkpoint("sdf3-robust-adaptive-promoted")?;
            Ok(true)
        })();
        match attempt {
            Ok(true) => {}
            Ok(false) => break,
            Err(error) => {
                report.stop = RobustAdaptiveStop3::RefinementStopped;
                report.error = Some(error);
                break;
            }
        }
    }
    report.work = control.work();
    Ok(report)
}

#[cfg(test)]
#[path = "adaptive_tests.rs"]
mod tests;
