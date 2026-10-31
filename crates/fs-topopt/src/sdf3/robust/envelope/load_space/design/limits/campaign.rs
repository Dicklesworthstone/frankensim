//! Repeated full-load-set assessment and refitting without replaying completed work.
use std::ops::ControlFlow;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use super::refit::LoadSpaceCapRefit3;
use super::refinement::LoadSpaceCapRefinement3;
use crate::sdf3::AdaptiveSdf3Elasticity;
use crate::sdf3_goal::{GoalRefinementError3, GoalRefinementOptions3};
use crate::{EvaluationStop, SolveControl, SolveWork};
pub use crate::sdf3::robust::envelope::limits::campaign::ComplianceCapCampaignError3
    as LoadSpaceCapCampaignError3;

#[derive(Debug, Clone, Copy)]
pub struct LoadSpaceCapCampaignOptions3 {
    /// Successful replacements across ALL run segments, not per call.
    pub max_refits: usize,
    /// Candidate build/refit attempts, including failures; never refunded.
    pub max_candidate_attempts: usize,
    /// Maximum normalized Frobenius change of the WHOLE compliance Gram.
    /// Zero requests exact two-grid agreement, not a continuum certificate.
    pub change_tolerance: f64,
    pub enrichment: GoalRefinementOptions3,
}
impl Default for LoadSpaceCapCampaignOptions3 {
    fn default() -> Self {
        Self { max_refits: 2, max_candidate_attempts: 4, change_tolerance: 0.02,
            enrichment: GoalRefinementOptions3::default() }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoadSpaceCapCampaignStop3 {
    Paused,
    /// The CURRENT material passed both physical caps and the full-Gram gate.
    StableAndFeasible,
    /// The current material was assessed, but still fails the requested check.
    RefinementLimit,
    CandidateAttemptLimit,
    /// Inspect rejected(): restoration, optimization or final recheck stopped.
    CandidateStopped,
}
#[derive(Debug, Clone, PartialEq)]
pub struct LoadSpaceCapGrid3 {
    pub active_cells: usize,
    pub background_cells: usize,
    pub dilated_volume: f64,
    pub worst_compliance: f64,
    pub work: SolveWork,
}
#[derive(Debug, Clone, PartialEq)]
pub struct LoadSpaceCapProbe3 {
    pub source_grid: usize,
    pub coarse_worst_compliance: f64,
    pub fine_worst_compliance: f64,
    pub normalized_gram_change: f64,
    pub passed: bool,
    pub work: SolveWork,
}
fn grid<O: AdaptiveSdf3Elasticity>(current: &LoadSpaceCapRefit3<O>, tree: &Octree3) -> LoadSpaceCapGrid3 {
    let e = &current.recheck.as_ref().expect("completed solve").as_ref().expect("rechecked point").evaluation;
    LoadSpaceCapGrid3 { active_cells: current.study.cells(), background_cells: tree.leaves().len(),
        dilated_volume: e.physical.dilated_volume_fraction, worst_compliance: e.worst_compliance, work: current.work }
}
fn check_tree<O: AdaptiveSdf3Elasticity>(op: &O, tree: &Octree3) -> Result<(), GoalRefinementError3> {
    if op.adaptive().leaves().iter().all(|cell| tree.leaves().contains(cell)) { Ok(()) }
    else { Err(GoalRefinementError3::Invalid("load-space campaign operator differs from its background")) }
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

/// An in-memory adaptive design session for the SAME authored load ellipsoid.
/// Start with solve_reference_load_space_cap. All parameters, caps, reference
/// laws, builder and the original physical work context stay bound to this owner.
/// The builder must preserve physical identity, honor its checkpoint and enforce
/// ONE cumulative geometry budget including failed builds; the tree caps growth.
///
/// Every active source cell is refined. Basis-wise DWR indicators do not localize
/// the entire load-set error, so this API does not silently use them as marks.
/// The diagnostic inherits physical eroded stiffness; the candidate instead
/// inherits raw controls, refilters and freshly minimizes material. Completed
/// probes and prepared candidates are retained across run calls. Failed candidates
/// remain separately inspectable and retries consume another candidate attempt.
/// Tree, operator and rechecked fields are installed together after a checkpoint.
///
/// Segments pause between grid refits, not within an AL solve. This is not disk
/// checkpointing, cross-grid descent, a continuum bound, or a stationarity claim.
/// Old full operators/fields are released on promotion; summaries retain history.
pub struct LoadSpaceCapCampaign3<'a, 'callback, O: AdaptiveSdf3Elasticity, B> {
    current: LoadSpaceCapRefit3<O>,
    tree: Octree3,
    laws: &'a [ReferenceLoad3<'a>],
    control: &'a mut SolveControl<'callback>,
    build: B,
    options: LoadSpaceCapCampaignOptions3,
    grids: Vec<LoadSpaceCapGrid3>,
    probes: Vec<LoadSpaceCapProbe3>,
    checked: Option<(Octree3, LoadSpaceCapRefinement3)>,
    pending: Option<(Octree3, LoadSpaceCapRefit3<O>)>,
    rejected: Option<(Octree3, LoadSpaceCapRefit3<O>)>,
    attempts: usize,
    terminal: Option<LoadSpaceCapCampaignStop3>,
}
impl<'a, 'callback, O: AdaptiveSdf3Elasticity, B> LoadSpaceCapCampaign3<'a, 'callback, O, B>
where B: FnMut(&Octree3, &mut dyn FnMut() -> ControlFlow<()>) -> Result<O, GoalRefinementError3> {
    pub fn new(current: LoadSpaceCapRefit3<O>, tree: Octree3, laws: &'a [ReferenceLoad3<'a>],
        options: LoadSpaceCapCampaignOptions3, control: &'a mut SolveControl<'callback>, build: B,
    ) -> Result<Self, LoadSpaceCapCampaignError3> {
        current.study.require_homogeneous_load_space().map_err(GoalRefinementError3::Invalid)?;
        if !current.completed_and_rechecked() || current.work != control.work()
            || laws.len() != current.cap.policy().dimension() || laws.len() > options.enrichment.max_load_cases
            || !options.change_tolerance.is_finite() || options.change_tolerance < 0.0
            || options.enrichment.max_transfer_terms == 0
            || ![options.enrichment.numerical.residual_tolerance, options.enrichment.numerical.identity_tolerance]
                .iter().all(|t| t.is_finite() && *t > 0.0 && *t < 1.0) {
            return Err(GoalRefinementError3::Invalid("load-space campaign needs a rechecked design, original work and valid policy").into());
        }
        check_tree(current.study.operator(), &tree)?;
        let p = current.study.params(); let q = current.params;
        let point = current.recheck.as_ref().unwrap().as_ref().unwrap();
        if [p.penal, p.beta, p.eta, p.e_min] != [q.penal, q.beta, q.eta, q.e_min]
            || point.rho.len() != current.study.cells()
            || point.rho.iter().any(|r| !r.is_finite() || *r < current.options.design.density_floor || *r > 1.0)
            || point.evaluation.physical.robustness != current.options.design.robustness {
            return Err(GoalRefinementError3::Invalid("load-space campaign model does not match its recheck").into());
        }
        let grids = vec![grid(&current, &tree)];
        Ok(Self { current, tree, laws, control, build, options, grids, probes: Vec::new(),
            checked: None, pending: None, rejected: None, attempts: 0, terminal: None })
    }
    pub fn current(&self) -> &LoadSpaceCapRefit3<O> { &self.current }
    pub fn tree(&self) -> &Octree3 { &self.tree }
    pub fn grids(&self) -> &[LoadSpaceCapGrid3] { &self.grids }
    pub fn probes(&self) -> &[LoadSpaceCapProbe3] { &self.probes }
    pub fn current_probe(&self) -> Option<&LoadSpaceCapRefinement3> { self.checked.as_ref().map(|(_, p)| p) }
    pub fn pending(&self) -> Option<&(Octree3, LoadSpaceCapRefit3<O>)> { self.pending.as_ref() }
    pub fn rejected(&self) -> Option<&(Octree3, LoadSpaceCapRefit3<O>)> { self.rejected.as_ref() }
    pub fn work(&self) -> SolveWork { self.control.work() }
    pub fn candidate_attempts(&self) -> usize { self.attempts }
    pub fn into_current(self) -> (Octree3, LoadSpaceCapRefit3<O>) { (self.tree, self.current) }

    /// At most additional_refits promotions; zero is a no-op. A positive segment
    /// assesses its final installed grid before pausing or reaching a count cap.
    /// Only StableAndFeasible is a two-grid pass. Cached terminal outcomes and
    /// completed probes consume no additional solves. A stopped promotion keeps
    /// its completed candidate for retry; a stopped refit keeps the old grid.
    pub fn run(&mut self, additional_refits: usize) -> Result<LoadSpaceCapCampaignStop3, LoadSpaceCapCampaignError3> {
        use LoadSpaceCapCampaignStop3 as Stop;
        if let Some(stop) = self.terminal { return Ok(stop); }
        if additional_refits == 0 { return Ok(Stop::Paused); }
        let mut promoted = 0;
        loop {
            if self.checked.is_none() {
                let tree = self.tree.refined(self.current.study.operator().adaptive().leaves(), || {
                    if self.control.checkpoint("sdf3-load-space-campaign-probe-tree").is_ok() { ControlFlow::Continue(()) }
                    else { ControlFlow::Break(()) }
                })?;
                let op = build_checked(&tree, self.control, "sdf3-load-space-campaign-probe-build", &mut self.build)?;
                let point = self.current.recheck.as_ref().unwrap().as_ref().unwrap();
                let checked = self.current.study.estimate_load_space_cap_enrichment(op.into_adaptive(), self.laws,
                    &point.rho, self.current.params, self.current.options.design.robustness,
                    &self.current.cap, self.options.enrichment, self.control)?;
                let passed = checked.stable_and_feasible(self.current.options.design.optimizer.tolerance, self.options.change_tolerance);
                self.probes.push(LoadSpaceCapProbe3 { source_grid: self.grids.len()-1,
                    coarse_worst_compliance: checked.coarse.worst_compliance,
                    fine_worst_compliance: checked.fine.worst_compliance,
                    normalized_gram_change: checked.normalized_gram_change, passed, work: checked.work });
                self.checked = Some((tree, checked));
            }
            let stop = if self.probes.last().unwrap().passed { Some(Stop::StableAndFeasible) }
                else if self.grids.len()-1 >= self.options.max_refits { Some(Stop::RefinementLimit) }
                else { None };
            if let Some(stop) = stop { self.terminal = Some(stop); return Ok(stop); }
            if promoted == additional_refits { return Ok(Stop::Paused); }
            if self.pending.is_none() {
                if self.attempts >= self.options.max_candidate_attempts {
                    self.terminal = Some(Stop::CandidateAttemptLimit); return Ok(Stop::CandidateAttemptLimit);
                }
                let tree = self.checked.as_ref().expect("complete probe").0.clone();
                self.attempts += 1; // Charge before building; failures are not refunded.
                let op = build_checked(&tree, self.control, "sdf3-load-space-campaign-candidate-build", &mut self.build)?;
                let raw = &self.current.recheck.as_ref().unwrap().as_ref().unwrap().rho;
                let candidate = self.current.study.refit_load_space_cap(op, raw, self.current.params,
                    self.laws, &self.current.cap, self.current.options, self.control)?;
                if !candidate.completed_and_rechecked() {
                    self.rejected = Some((tree, candidate)); return Ok(Stop::CandidateStopped);
                }
                self.pending = Some((tree, candidate));
            }
            self.control.checkpoint("sdf3-load-space-campaign-promote")?;
            let (tree, candidate) = self.pending.take().expect("complete rechecked candidate");
            self.current = candidate; self.tree = tree;
            self.grids.push(grid(&self.current, &self.tree));
            self.checked = None; self.rejected = None;
            promoted += 1;
            self.control.checkpoint("sdf3-load-space-campaign-promoted")?;
        }
    }
}

#[cfg(test)]
#[path = "campaign_tests.rs"]
mod tests;
