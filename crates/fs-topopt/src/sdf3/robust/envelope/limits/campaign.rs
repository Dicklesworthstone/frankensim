//! Resumable cap-constrained refinement, with no replay of completed probes.
use std::ops::ControlFlow;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::{Octree3, OctreeError3};
use super::refit::{ComplianceCapRefit3, ComplianceCapRefitError3};
use super::refinement::{ComplianceCapCheck3, ComplianceCapRefinement3};
use crate::sdf3::AdaptiveSdf3Elasticity;
use crate::sdf3_goal::{GoalRefinementError3, GoalRefinementOptions3};
use crate::{EvaluationStop, SolveControl, SolveWork};

#[derive(Debug, Clone, Copy)]
pub struct ComplianceCapCampaignOptions3 {
    /// Successful grid replacements over ALL run segments.
    pub max_refits: usize,
    /// Candidate build/refit attempts, including failures. Never refunded.
    pub max_candidate_attempts: usize,
    pub marking_fraction: f64,
    pub max_marks: usize,
    /// Maximum absolute coarse/fine change divided by EACH authored cap.
    /// Zero requests exact two-grid agreement, not a continuum certificate.
    pub change_tolerance: f64,
    pub enrichment: GoalRefinementOptions3,
}
impl Default for ComplianceCapCampaignOptions3 {
    fn default() -> Self {
        Self { max_refits: 4, max_candidate_attempts: 8, marking_fraction: 0.5,
            max_marks: 2, change_tolerance: 0.02, enrichment: GoalRefinementOptions3::default() }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComplianceCapCampaignStop3 {
    Paused,
    /// Current material passes every cap on BOTH grids and the change gate.
    StableAndFeasible,
    /// Current material has been probed but the requested check still fails.
    RefinementLimit,
    CandidateAttemptLimit,
    NoRefinementSignal,
    /// Inspect rejected(). Restoration, optimizer or final recheck stopped.
    CandidateStopped,
}
#[derive(Debug, Clone)]
pub enum ComplianceCapCampaignError3 {
    Goal(GoalRefinementError3),
    Refit(ComplianceCapRefitError3),
    Background(OctreeError3),
}
impl From<GoalRefinementError3> for ComplianceCapCampaignError3 {
    fn from(e: GoalRefinementError3) -> Self { Self::Goal(e) }
}
impl From<ComplianceCapRefitError3> for ComplianceCapCampaignError3 {
    fn from(e: ComplianceCapRefitError3) -> Self { Self::Refit(e) }
}
impl From<EvaluationStop> for ComplianceCapCampaignError3 {
    fn from(e: EvaluationStop) -> Self { Self::Goal(e.into()) }
}
impl From<OctreeError3> for ComplianceCapCampaignError3 {
    fn from(e: OctreeError3) -> Self {
        if e == OctreeError3::Cancelled { EvaluationStop::Cancelled.into() }
        else { Self::Background(e) }
    }
}
impl std::fmt::Display for ComplianceCapCampaignError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "compliance-cap campaign: {self:?}") }
}
impl std::error::Error for ComplianceCapCampaignError3 {}

/// Small retained summary; old operators and full displacement families are
/// released on promotion. The current owner retains the complete solved fields.
#[derive(Debug, Clone, PartialEq)]
pub struct ComplianceCapGrid3 {
    pub active_cells: usize,
    pub background_cells: usize,
    pub dilated_volume: f64,
    pub compliances: Vec<f64>,
    pub work: SolveWork,
}
#[derive(Debug, Clone, PartialEq)]
pub struct ComplianceCapProbe3 {
    pub source_grid: usize,
    pub checks: Vec<ComplianceCapCheck3>,
    pub coarse_envelope: f64,
    pub fine_envelope: f64,
    pub passed: bool,
    pub work: SolveWork,
}
fn grid<O: AdaptiveSdf3Elasticity>(current: &ComplianceCapRefit3<O>, tree: &Octree3) -> ComplianceCapGrid3 {
    let e = &current.recheck.as_ref().expect("completed solve").as_ref().expect("rechecked point").evaluation.physical;
    ComplianceCapGrid3 { active_cells: current.study.cells(), background_cells: tree.leaves().len(),
        dilated_volume: e.dilated_volume_fraction, compliances: e.eroded.objective.case_compliances.clone(), work: current.work }
}
fn check_tree<O: AdaptiveSdf3Elasticity>(op: &O, tree: &Octree3) -> Result<(), GoalRefinementError3> {
    if op.adaptive().leaves().iter().all(|cell| tree.leaves().contains(cell)) { Ok(()) }
    else { Err(GoalRefinementError3::Invalid("cap campaign operator does not belong to its background")) }
}
fn build_checked<O: AdaptiveSdf3Elasticity>(tree: &Octree3, control: &mut SolveControl<'_>, stage: &'static str,
    build: &mut impl FnMut(&Octree3, &mut dyn FnMut() -> ControlFlow<()>) -> Result<O, GoalRefinementError3>,
) -> Result<O, GoalRefinementError3> {
    control.checkpoint(stage)?;
    let mut stopped = None;
    let result = build(tree, &mut || match control.checkpoint(stage) {
        Ok(()) => ControlFlow::Continue(()),
        Err(e) => { stopped = Some(e); ControlFlow::Break(()) }
    });
    if let Some(e) = stopped { return Err(e.into()); }
    control.checkpoint(stage)?;
    let op = result?; check_tree(&op, tree)?; Ok(op)
}

/// Exclusive in-memory campaign: caps, model, laws, builder and the physical
/// budget remain bound across segments. Start with solve_reference_compliance_caps.
/// The builder must preserve SDF/material/support/integration identity, honor
/// its checkpoint and retain ONE cumulative geometry allowance including failures.
///
/// Each probe inherits physical eroded stiffness; each candidate refit instead
/// transfers raw controls, rebuilds its filter and minimizes material under the
/// unchanged caps. Probe success never reuses evidence from an earlier design.
/// A successful candidate is held until the publication checkpoint, so retrying
/// a cancelled publication need not rebuild or re-solve it. Failed candidates
/// are separately inspectable, not installed; their retries consume another
/// attempt and share the original physical budget.
///
/// Run segments pause between refits, not inside a candidate's AL solve. This
/// does not provide cross-process checkpoints, continuum load guarantees,
/// volume discretization bounds, stationarity, or cross-grid descent.
pub struct ComplianceCapCampaign3<'a, 'callback, O: AdaptiveSdf3Elasticity, B> {
    current: ComplianceCapRefit3<O>,
    tree: Octree3,
    laws: &'a [ReferenceLoad3<'a>],
    control: &'a mut SolveControl<'callback>,
    build: B,
    options: ComplianceCapCampaignOptions3,
    grids: Vec<ComplianceCapGrid3>,
    probes: Vec<ComplianceCapProbe3>,
    checked: Option<ComplianceCapRefinement3>,
    pending: Option<(Octree3, ComplianceCapRefit3<O>)>,
    rejected: Option<(Octree3, ComplianceCapRefit3<O>)>,
    attempts: usize,
    terminal: Option<ComplianceCapCampaignStop3>,
}
impl<'a, 'callback, O: AdaptiveSdf3Elasticity, B> ComplianceCapCampaign3<'a, 'callback, O, B>
where B: FnMut(&Octree3, &mut dyn FnMut() -> ControlFlow<()>) -> Result<O, GoalRefinementError3> {
    pub fn new(current: ComplianceCapRefit3<O>, tree: Octree3, laws: &'a [ReferenceLoad3<'a>],
        options: ComplianceCapCampaignOptions3, control: &'a mut SolveControl<'callback>, build: B,
    ) -> Result<Self, ComplianceCapCampaignError3> {
        current.study.require_material_body_recipe("sdf3-cap-campaign-requires-material-acceleration")?;
        if !current.completed_and_rechecked() || current.work != control.work()
            || laws.len() != current.caps.caps().len() || laws.len() > options.enrichment.max_load_cases
            || !options.change_tolerance.is_finite() || options.change_tolerance < 0.0
            || !options.marking_fraction.is_finite() || options.marking_fraction <= 0.0 || options.marking_fraction > 1.0
            || options.max_marks == 0 || options.enrichment.max_transfer_terms == 0
            || ![options.enrichment.numerical.residual_tolerance, options.enrichment.numerical.identity_tolerance]
                .iter().all(|t| t.is_finite() && *t > 0.0 && *t < 1.0) {
            return Err(GoalRefinementError3::Invalid("cap campaign needs a completed rechecked design, original work and valid policy").into());
        }
        check_tree(current.study.operator(), &tree)?;
        let p = current.study.params(); let q = current.params;
        let point = current.recheck.as_ref().unwrap().as_ref().unwrap();
        if [p.penal, p.beta, p.eta, p.e_min] != [q.penal, q.beta, q.eta, q.e_min]
            || point.rho.len() != current.study.cells()
            || point.rho.iter().any(|r| !r.is_finite() || *r < current.options.design.density_floor || *r > 1.0)
            || point.evaluation.physical.robustness != current.options.design.robustness {
            return Err(GoalRefinementError3::Invalid("cap campaign model or raw design does not match its recheck").into());
        }
        let grids = vec![grid(&current, &tree)];
        Ok(Self { current, tree, laws, control, build, options, grids, probes: Vec::new(),
            checked: None, pending: None, rejected: None, attempts: 0, terminal: None })
    }
    pub fn current(&self) -> &ComplianceCapRefit3<O> { &self.current }
    pub fn tree(&self) -> &Octree3 { &self.tree }
    pub fn grids(&self) -> &[ComplianceCapGrid3] { &self.grids }
    pub fn probes(&self) -> &[ComplianceCapProbe3] { &self.probes }
    pub fn current_probe(&self) -> Option<&ComplianceCapRefinement3> { self.checked.as_ref() }
    pub fn pending(&self) -> Option<&(Octree3, ComplianceCapRefit3<O>)> { self.pending.as_ref() }
    pub fn rejected(&self) -> Option<&(Octree3, ComplianceCapRefit3<O>)> { self.rejected.as_ref() }
    pub fn work(&self) -> SolveWork { self.control.work() }
    pub fn candidate_attempts(&self) -> usize { self.attempts }

    /// Up to additional_refits promotions. Zero is a true no-op. A positive
    /// segment always assesses its latest grid before pausing or reaching a
    /// count limit. Completed probes are cached across run calls. Terminal
    /// outcomes are also cached; only StableAndFeasible is a numerical pass.
    pub fn run(&mut self, additional_refits: usize) -> Result<ComplianceCapCampaignStop3, ComplianceCapCampaignError3> {
        use ComplianceCapCampaignStop3 as Stop;
        if let Some(stop) = self.terminal { return Ok(stop); }
        if additional_refits == 0 { return Ok(Stop::Paused); }
        let mut promoted = 0;
        loop {
            if self.checked.is_none() {
                let tree = self.tree.refined(self.current.study.operator().adaptive().leaves(), || {
                    if self.control.checkpoint("sdf3-cap-campaign-probe-tree").is_ok() { ControlFlow::Continue(()) }
                    else { ControlFlow::Break(()) }
                })?;
                let op = build_checked(&tree, self.control, "sdf3-cap-campaign-probe-build", &mut self.build)?;
                let point = self.current.recheck.as_ref().unwrap().as_ref().unwrap();
                let checked = self.current.study.estimate_compliance_cap_enrichment(op.into_adaptive(), self.laws,
                    &point.rho, self.current.params, self.current.options.design.robustness,
                    &self.current.caps, self.options.enrichment, self.control)?;
                let passed = checked.stable_and_feasible(self.current.options.design.optimizer.tolerance, self.options.change_tolerance);
                self.probes.push(ComplianceCapProbe3 { source_grid: self.grids.len()-1,
                    checks: checked.checks.clone(), coarse_envelope: checked.coarse_envelope,
                    fine_envelope: checked.fine_envelope, passed, work: checked.work });
                self.checked = Some(checked);
            }
            let stop = if self.probes.last().unwrap().passed { Some(Stop::StableAndFeasible) }
                else if self.grids.len()-1 == self.options.max_refits { Some(Stop::RefinementLimit) }
                else { None };
            if let Some(stop) = stop { self.terminal = Some(stop); return Ok(stop); }
            if promoted == additional_refits { return Ok(Stop::Paused); }
            if self.pending.is_none() {
                if self.attempts == self.options.max_candidate_attempts {
                    self.terminal = Some(Stop::CandidateAttemptLimit); return Ok(Stop::CandidateAttemptLimit);
                }
                let marks = self.checked.as_ref().unwrap().mark(self.options.marking_fraction, self.options.max_marks, || {
                    if self.control.checkpoint("sdf3-cap-campaign-mark").is_ok() { ControlFlow::Continue(()) }
                    else { ControlFlow::Break(()) }
                }).map_err(GoalRefinementError3::from)?;
                if marks.marked.is_empty() { self.terminal = Some(Stop::NoRefinementSignal); return Ok(Stop::NoRefinementSignal); }
                let tree = self.tree.refined(&marks.marked, || {
                    if self.control.checkpoint("sdf3-cap-campaign-select-tree").is_ok() { ControlFlow::Continue(()) }
                    else { ControlFlow::Break(()) }
                })?;
                self.attempts += 1; // Reserve before the builder; failures are charged.
                let op = build_checked(&tree, self.control, "sdf3-cap-campaign-candidate-build", &mut self.build)?;
                let raw = &self.current.recheck.as_ref().unwrap().as_ref().unwrap().rho;
                let candidate = self.current.study.refit_compliance_caps(op, raw, self.current.params,
                    self.laws, &self.current.caps, self.current.options, self.control)?;
                if !candidate.completed_and_rechecked() {
                    self.rejected = Some((tree, candidate)); return Ok(Stop::CandidateStopped);
                }
                self.pending = Some((tree, candidate));
            }
            self.control.checkpoint("sdf3-cap-campaign-promote")?;
            let (tree, candidate) = self.pending.take().expect("complete rechecked candidate");
            // No fallible computation between installing geometry and its fields.
            self.current = candidate; self.tree = tree;
            self.grids.push(grid(&self.current, &self.tree));
            self.checked = None; self.rejected = None;
            promoted += 1;
            self.control.checkpoint("sdf3-cap-campaign-promoted")?;
        }
    }
}
