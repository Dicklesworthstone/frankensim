//! Automatic native-grid assessment and reoptimization for frequency design.
//! The same material pencil, optimizer and immutable physical requirements are
//! retained across probes, candidate solves and atomic grid promotions.

use super::design::{Sdf3EigenfrequencyGoal3, admit, normalized};
use super::refinement::{
    Sdf3EigenfrequencyAdaptiveError3, Sdf3EigenfrequencyRefinement3,
    Sdf3EigenfrequencyRefinementOptions3,
};
use super::refit::Sdf3EigenfrequencyRefit3;
use crate::sdf3::{AdaptiveSdf3Elasticity, Sdf3Elasticity};
use crate::sdf3::response::limits::design::MinimumMaterialObjective3;
use crate::sdf3_goal::GoalRefinementError3;
use crate::{SolveControl, SolveControlCheckpoint, SolveWork};
use fs_cutfem::octree3::Octree3;
use std::ops::ControlFlow;

/// The same physical builder serves diagnostic and optimization grids.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sdf3EigenfrequencyBuild3 {
    Probe,
    Refit,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Sdf3EigenfrequencyMarking3 {
    UniformActive,
    /// Coarse-mode residuals weighted by projection onto the fine selected
    /// subspace. This is a diagnostic indicator, not a certified error bound.
    /// Unresolved physical inventory, or empty marks, uses uniform refinement.
    ModalResidual { theta: f64, max_marks: usize },
}

#[derive(Debug, Clone, Copy)]
pub struct Sdf3EigenfrequencyCampaignOptions3 {
    /// Successful promotions across all run segments.
    pub max_refits: usize,
    /// Candidate attempts, including failed geometry or optimization work.
    pub max_candidate_attempts: usize,
    /// Maximum selected eigenvalue change divided by the original reference.
    pub eigenvalue_tolerance: f64,
    /// Change in the actual normalized frequency or inventory objective.
    pub objective_tolerance: f64,
    /// Absolute change in numerical physical volume fraction.
    pub volume_tolerance: f64,
    /// |m_f-m_c|/max(m_f,m_c); physical kg use linear occupancy.
    pub mass_relative_tolerance: f64,
    /// Fine-mass distance between whole selected modal subspaces.
    pub subspace_tolerance: f64,
    pub marking: Sdf3EigenfrequencyMarking3,
    pub enrichment: Sdf3EigenfrequencyRefinementOptions3,
}
impl Default for Sdf3EigenfrequencyCampaignOptions3 {
    fn default() -> Self {
        Self {
            max_refits: 2,
            max_candidate_attempts: 4,
            eigenvalue_tolerance: 0.02,
            objective_tolerance: 0.02,
            volume_tolerance: 1e-4,
            mass_relative_tolerance: 0.01,
            subspace_tolerance: 0.05,
            marking: Sdf3EigenfrequencyMarking3::UniformActive,
            enrichment: Sdf3EigenfrequencyRefinementOptions3::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sdf3EigenfrequencyCampaignStop3 {
    Paused,
    /// Both native grids pass original physical bounds and every change gate.
    StableAndFeasible,
    RefinementLimit,
    CandidateAttemptLimit,
    /// The failed candidate and its completed optimizer prefix are retained.
    CandidateStopped,
}

/// Compact accepted-grid history; superseded operators and full fields drop.
#[derive(Debug, Clone, PartialEq)]
pub struct Sdf3EigenfrequencyGridSummary3 {
    pub active_cells: usize,
    pub background_cells: usize,
    pub objective: f64,
    pub aggregate: f64,
    pub lowest_eigenvalue: f64,
    pub volume_fraction: f64,
    pub physical_mass_kg: f64,
    pub work: SolveWork,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sdf3EigenfrequencyProbe3 {
    pub source_grid: usize,
    pub coarse_aggregate: f64,
    pub fine_aggregate: f64,
    pub eigenvalue_change: f64,
    pub objective_change: f64,
    pub volume_change: f64,
    pub mass_change_kg: f64,
    pub mass_relative_change: f64,
    pub subspace_change: f64,
    pub passed: bool,
    pub work: SolveWork,
}

fn invalid(message: &'static str) -> Sdf3EigenfrequencyAdaptiveError3 {
    GoalRefinementError3::Invalid(message).into()
}

// Pure geometry builders use the existing refinement error owner. Preserve
// the original control failure through their narrower cancellation callback.
fn controlled<T>(
    control: &mut SolveControl<'_>,
    stage: &'static str,
    action: impl FnOnce(&mut dyn FnMut() -> ControlFlow<()>) -> Result<T, GoalRefinementError3>,
) -> Result<T, Sdf3EigenfrequencyAdaptiveError3> {
    control.checkpoint(stage)?;
    let mut stopped = None;
    let result = action(&mut || match control.checkpoint(stage) {
        Ok(()) => ControlFlow::Continue(()),
        Err(error) => {
            stopped = Some(error);
            ControlFlow::Break(())
        }
    });
    if let Some(error) = stopped { return Err(error.into()); }
    let value = result?;
    control.checkpoint(stage)?;
    Ok(value)
}

fn check_tree<O: AdaptiveSdf3Elasticity>(
    op: &O, tree: &Octree3, control: &mut SolveControl<'_>,
) -> Result<(), Sdf3EigenfrequencyAdaptiveError3> {
    for (i, cell) in op.adaptive().leaves().iter().enumerate() {
        if i % 64 == 0 { control.checkpoint("sdf3-modal-campaign-background")?; }
        if !tree.leaves().contains(cell) {
            return Err(invalid("modal operator differs from its declared background"));
        }
    }
    Ok(())
}

fn grid<O: AdaptiveSdf3Elasticity>(
    current: &Sdf3EigenfrequencyRefit3<O>, tree: &Octree3,
) -> Result<Sdf3EigenfrequencyGridSummary3, Sdf3EigenfrequencyAdaptiveError3> {
    let point = current.recheck.as_ref().and_then(|r| r.as_ref().ok())
        .ok_or_else(|| invalid("modal grid has no freshly rechecked point"))?;
    let e = &point.evaluation;
    let objective = match current.goal {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { .. } =>
            normalized(-e.objective.aggregate, current.options.design.eigenvalue_reference)?,
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::VolumeFraction, volume_scale,
        } => normalized(e.volume_fraction, volume_scale)?,
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { scale_kg }, ..
        } => normalized(e.physical_mass_kg, scale_kg)?,
    };
    let lowest_eigenvalue = e.objective.modal.values.first().copied()
        .ok_or_else(|| invalid("modal grid has no eigenvalue"))?;
    Ok(Sdf3EigenfrequencyGridSummary3 {
        active_cells: current.study.cells(), background_cells: tree.leaves().len(),
        objective, aggregate: e.objective.aggregate, lowest_eigenvalue,
        volume_fraction: e.volume_fraction, physical_mass_kg: e.physical_mass_kg,
        work: current.work,
    })
}

fn mass_relative_change(report: &Sdf3EigenfrequencyRefinement3) -> Option<f64> {
    let a = report.coarse.physical_mass_kg;
    let b = report.fine.physical_mass_kg;
    if !a.is_finite() || !b.is_finite() || a < 0.0 || b < 0.0 { return None; }
    let scale = a.max(b);
    let value = if scale == 0.0 { 0.0 } else { (b - a).abs() / scale };
    value.is_finite().then_some(value)
}
fn volume_and_mass_resolved(
    report: &Sdf3EigenfrequencyRefinement3, options: Sdf3EigenfrequencyCampaignOptions3,
) -> bool {
    report.volume_change.is_finite() && report.volume_change >= 0.0
        && report.volume_change <= options.volume_tolerance
        && mass_relative_change(report).is_some_and(|v| v <= options.mass_relative_tolerance)
}
fn resolved(
    report: &Sdf3EigenfrequencyRefinement3, options: Sdf3EigenfrequencyCampaignOptions3,
) -> bool {
    report.physically_feasible()
        && volume_and_mass_resolved(report, options)
        && report.eigenvalue_change.is_finite() && report.eigenvalue_change >= 0.0
        && report.eigenvalue_change <= options.eigenvalue_tolerance
        && report.objective_change().is_some_and(|v|
            v.is_finite() && v >= 0.0 && v <= options.objective_tolerance)
        && report.subspace_change.is_finite() && report.subspace_change >= 0.0
        && report.subspace_change <= options.subspace_tolerance
}


/// Owned pause state for rebinding an adaptive campaign to another cancellation
/// callback without rebuilding validated physics or renewing numerical work.
/// This is an in-memory handoff, not a byte-serialization or cross-process
/// identity claim. The state is opaque: only a live campaign can create it.
pub struct Sdf3EigenfrequencyCampaignState3<O: AdaptiveSdf3Elasticity> {
    current: Sdf3EigenfrequencyRefit3<O>,
    tree: Octree3,
    options: Sdf3EigenfrequencyCampaignOptions3,
    grids: Vec<Sdf3EigenfrequencyGridSummary3>,
    probes: Vec<Sdf3EigenfrequencyProbe3>,
    checked: Option<(Octree3, Sdf3EigenfrequencyRefinement3)>,
    pending: Option<(Octree3, Sdf3EigenfrequencyRefit3<O>)>,
    rejected: Option<(Octree3, Sdf3EigenfrequencyRefit3<O>)>,
    attempts: usize,
    terminal: Option<Sdf3EigenfrequencyCampaignStop3>,
    control: SolveControlCheckpoint,
}

impl<O: AdaptiveSdf3Elasticity> Sdf3EigenfrequencyCampaignState3<O> {
    #[must_use]
    pub fn control_checkpoint(&self) -> SolveControlCheckpoint { self.control }
    #[must_use]
    pub fn current(&self) -> &Sdf3EigenfrequencyRefit3<O> { &self.current }
    #[must_use]
    pub fn tree(&self) -> &Octree3 { &self.tree }
    #[must_use]
    pub fn pending(&self) -> Option<&(Octree3, Sdf3EigenfrequencyRefit3<O>)> {
        self.pending.as_ref()
    }
    #[must_use]
    pub fn candidate_attempts(&self) -> usize { self.attempts }
}

/// In-memory adaptive design with native material-preserving transfer.
///
/// Probes copy physical occupancy and stiffness onto a uniformly enriched
/// active space and independently solve both complete selected modal groups.
/// Candidates inherit RAW controls, refilter on their own grid, prepare their
/// own native inertia and rerun the same constrained design owner. A candidate
/// is installed only after the original physical requirements pass a fresh
/// modal solve. No monotone frequency-restoration or cross-grid descent is
/// assumed: changing occupancy changes stiffness AND inertia.
///
/// The builder must retain physical domain, supports and source material
/// assignments and honor one cumulative geometry budget across all calls.
/// Native transfer checks actual enrichment and each source material. The
/// original eigenvalue reference, aggregation policy and goals remain fixed.
/// Every installed grid receives a fresh full probe before a normal pause
/// or refinement-limit return.
///
/// Whole-subspace comparison tolerates rotations within a selected eigenspace;
/// it is numerical resolution evidence, not a continuum error or completeness
/// certificate. Pending rechecked candidates survive cancellation immediately
/// before promotion. Failed candidates retain their own completed solve prefix.
pub struct Sdf3EigenfrequencyCampaign3<'a, 'callback, O: AdaptiveSdf3Elasticity, B> {
    current: Sdf3EigenfrequencyRefit3<O>,
    tree: Octree3,
    options: Sdf3EigenfrequencyCampaignOptions3,
    control: &'a mut SolveControl<'callback>,
    build: B,
    grids: Vec<Sdf3EigenfrequencyGridSummary3>,
    probes: Vec<Sdf3EigenfrequencyProbe3>,
    checked: Option<(Octree3, Sdf3EigenfrequencyRefinement3)>,
    pending: Option<(Octree3, Sdf3EigenfrequencyRefit3<O>)>,
    rejected: Option<(Octree3, Sdf3EigenfrequencyRefit3<O>)>,
    attempts: usize,
    terminal: Option<Sdf3EigenfrequencyCampaignStop3>,
}

impl<'a, 'callback, O: AdaptiveSdf3Elasticity, B>
    Sdf3EigenfrequencyCampaign3<'a, 'callback, O, B>
where
    B: FnMut(
        &Octree3, Sdf3EigenfrequencyBuild3, &mut dyn FnMut() -> ControlFlow<()>,
    ) -> Result<O, GoalRefinementError3>,
{
    pub fn new(
        current: Sdf3EigenfrequencyRefit3<O>, tree: Octree3,
        options: Sdf3EigenfrequencyCampaignOptions3,
        control: &'a mut SolveControl<'callback>, build: B,
    ) -> Result<Self, Sdf3EigenfrequencyAdaptiveError3> {
        let marking_valid = match options.marking {
            Sdf3EigenfrequencyMarking3::UniformActive => true,
            Sdf3EigenfrequencyMarking3::ModalResidual { theta, max_marks } =>
                theta.is_finite() && theta > 0.0 && theta <= 1.0 && max_marks > 0,
        };
        if !current.completed_and_rechecked() || current.work != control.work()
            || !marking_valid || options.enrichment.max_transfer_terms == 0
            || options.enrichment.max_mass_terms == 0
            || !options.enrichment.orthogonality_tolerance.is_finite()
            || options.enrichment.orthogonality_tolerance <= 0.0
            || options.enrichment.orthogonality_tolerance >= 1.0
            || ![options.eigenvalue_tolerance, options.objective_tolerance,
                 options.volume_tolerance, options.mass_relative_tolerance,
                 options.subspace_tolerance].iter().all(|v| v.is_finite() && *v >= 0.0)
        {
            return Err(invalid("modal campaign needs a rechecked design, original work and valid options"));
        }
        current.study.operator().validate_material_mass(&current.mass)?;
        current.study.require_material_body_recipe("sdf3-modal-campaign-material-recipe")?;
        let point = current.recheck.as_ref().unwrap().as_ref().unwrap();
        admit(&current.study, &point.rho, current.goal, current.options.design)?;
        check_tree(current.study.operator(), &tree, control)?;
        let grids = vec![grid(&current, &tree)?];
        Ok(Self {
            current, tree, options, control, build, grids, probes: Vec::new(),
            checked: None, pending: None, rejected: None, attempts: 0, terminal: None,
        })
    }

    /// Detach all validated grid/probe/candidate state from the borrowed
    /// callback. The returned control receipt must be used to construct the
    /// next SolveControl; otherwise resume refuses before any physics work.
    pub fn detach(self) -> Sdf3EigenfrequencyCampaignState3<O> {
        let control = self.control.checkpoint_state();
        Sdf3EigenfrequencyCampaignState3 {
            current: self.current,
            tree: self.tree,
            options: self.options,
            grids: self.grids,
            probes: self.probes,
            checked: self.checked,
            pending: self.pending,
            rejected: self.rejected,
            attempts: self.attempts,
            terminal: self.terminal,
            control,
        }
    }

    /// Rebind a detached campaign to a new callback/control owner. Completed
    /// probes and freshly rechecked pending candidates are reused exactly; no
    /// geometry, modal solve, gradient gate or optimizer step is replayed.
    pub fn resume(
        state: Sdf3EigenfrequencyCampaignState3<O>,
        control: &'a mut SolveControl<'callback>,
        build: B,
    ) -> Result<Self, Sdf3EigenfrequencyAdaptiveError3> {
        if state.control != control.checkpoint_state() {
            return Err(invalid("modal campaign resume would change budget or cumulative work"));
        }
        state.current.study.operator().validate_material_mass(&state.current.mass)?;
        state.current.study.require_material_body_recipe("sdf3-modal-campaign-material-recipe")?;
        let point = state.current.recheck.as_ref().and_then(|r| r.as_ref().ok())
            .ok_or_else(|| invalid("modal campaign resume lacks a rechecked current design"))?;
        admit(&state.current.study, &point.rho, state.current.goal, state.current.options.design)?;
        check_tree(state.current.study.operator(), &state.tree, control)?;
        if let Some((tree, candidate)) = &state.pending {
            if !candidate.completed_and_rechecked() {
                return Err(invalid("modal campaign pending candidate lost its physical recheck"));
            }
            check_tree(candidate.study.operator(), tree, control)?;
        }
        Ok(Self {
            current: state.current,
            tree: state.tree,
            options: state.options,
            control,
            build,
            grids: state.grids,
            probes: state.probes,
            checked: state.checked,
            pending: state.pending,
            rejected: state.rejected,
            attempts: state.attempts,
            terminal: state.terminal,
        })
    }

    pub fn current(&self) -> &Sdf3EigenfrequencyRefit3<O> { &self.current }
    pub fn tree(&self) -> &Octree3 { &self.tree }
    pub fn grids(&self) -> &[Sdf3EigenfrequencyGridSummary3] { &self.grids }
    pub fn probes(&self) -> &[Sdf3EigenfrequencyProbe3] { &self.probes }
    pub fn current_probe(&self) -> Option<&Sdf3EigenfrequencyRefinement3> {
        self.checked.as_ref().map(|(_, report)| report)
    }
    pub fn pending(&self) -> Option<&(Octree3, Sdf3EigenfrequencyRefit3<O>)> {
        self.pending.as_ref()
    }
    pub fn rejected(&self) -> Option<&(Octree3, Sdf3EigenfrequencyRefit3<O>)> {
        self.rejected.as_ref()
    }
    pub fn work(&self) -> SolveWork { self.control.work() }
    pub fn candidate_attempts(&self) -> usize { self.attempts }
    pub fn into_current(self) -> (Octree3, Sdf3EigenfrequencyRefit3<O>) {
        (self.tree, self.current)
    }

    fn build_checked(
        &mut self, tree: &Octree3, kind: Sdf3EigenfrequencyBuild3,
    ) -> Result<O, Sdf3EigenfrequencyAdaptiveError3> {
        let stage = match kind {
            Sdf3EigenfrequencyBuild3::Probe => "sdf3-modal-campaign-probe-build",
            Sdf3EigenfrequencyBuild3::Refit => "sdf3-modal-campaign-refit-build",
        };
        let build = &mut self.build;
        let operator = controlled(self.control, stage, |poll| build(tree, kind, poll))?;
        check_tree(&operator, tree, self.control)?;
        Ok(operator)
    }

    fn candidate_tree(&mut self) -> Result<Octree3, Sdf3EigenfrequencyAdaptiveError3> {
        let (uniform, report) = self.checked.as_ref().expect("complete modal probe");
        let Sdf3EigenfrequencyMarking3::ModalResidual { theta, max_marks } = self.options.marking
        else { return Ok(uniform.clone()); };
        // Residual indicators do not localize unresolved quadrature inventory.
        // Minimum-material design also needs its actual inventory objective
        // resolved before modal marking can select a smaller candidate set.
        if !volume_and_mass_resolved(report, self.options)
            || (matches!(self.current.goal, Sdf3EigenfrequencyGoal3::MinimumMaterial { .. })
                && !report.objective_change().is_some_and(|v| v <= self.options.objective_tolerance))
        {
            return Ok(uniform.clone());
        }
        let marks = controlled(self.control, "sdf3-modal-campaign-mark", |poll| {
            report.mark(theta, max_marks, poll).map_err(GoalRefinementError3::from)
        })?;
        if marks.marked.is_empty() { return Ok(uniform.clone()); }
        controlled(self.control, "sdf3-modal-campaign-candidate-tree", |poll| {
            self.tree.refined(&marks.marked, poll).map_err(Into::into)
        })
    }

    /// Run at most additional_refits successful grid replacements. A zero
    /// allowance still assesses the current design. Finishing a positive
    /// segment always probes its final installed grid before returning.
    pub fn run(&mut self, additional_refits: usize)
        -> Result<Sdf3EigenfrequencyCampaignStop3, Sdf3EigenfrequencyAdaptiveError3> {
        use Sdf3EigenfrequencyCampaignStop3 as Stop;
        if let Some(stop) = self.terminal { return Ok(stop); }
        let mut promoted = 0;
        loop {
            if self.checked.is_none() {
                let tree = controlled(self.control, "sdf3-modal-campaign-probe-tree", |poll| {
                    self.tree.refined(self.current.study.operator().adaptive().leaves(), poll)
                        .map_err(Into::into)
                })?;
                let mut operator = self.build_checked(&tree, Sdf3EigenfrequencyBuild3::Probe)?;
                let point = self.current.recheck.as_ref().unwrap().as_ref().unwrap();
                let report = self.current.study.estimate_eigenfrequency_enrichment(
                    &self.current.mass, &mut operator, &point.rho, self.current.goal,
                    self.current.options.design, self.options.enrichment, self.control,
                )?;
                let objective_change = report.objective_change()
                    .ok_or_else(|| invalid("unrepresentable modal objective grid change"))?;
                let mass_relative_change = mass_relative_change(&report)
                    .ok_or_else(|| invalid("unrepresentable physical mass grid change"))?;
                let passed = resolved(&report, self.options);
                self.probes.push(Sdf3EigenfrequencyProbe3 {
                    source_grid: self.grids.len() - 1,
                    coarse_aggregate: report.coarse.aggregate,
                    fine_aggregate: report.fine.aggregate,
                    eigenvalue_change: report.eigenvalue_change, objective_change,
                    volume_change: report.volume_change, mass_change_kg: report.mass_change_kg,
                    mass_relative_change, subspace_change: report.subspace_change,
                    passed, work: report.work,
                });
                self.checked = Some((tree, report));
            }
            let stop = if self.probes.last().expect("complete modal probe").passed {
                Some(Stop::StableAndFeasible)
            } else if self.grids.len() - 1 >= self.options.max_refits {
                Some(Stop::RefinementLimit)
            } else { None };
            if let Some(stop) = stop {
                self.terminal = Some(stop);
                return Ok(stop);
            }
            if promoted == additional_refits { return Ok(Stop::Paused); }
            if self.pending.is_none() {
                if self.attempts >= self.options.max_candidate_attempts {
                    self.terminal = Some(Stop::CandidateAttemptLimit);
                    return Ok(Stop::CandidateAttemptLimit);
                }
                self.attempts += 1;
                let tree = self.candidate_tree()?;
                let operator = self.build_checked(&tree, Sdf3EigenfrequencyBuild3::Refit)?;
                let rho = &self.current.recheck.as_ref().unwrap().as_ref().unwrap().rho;
                let candidate = self.current.study.refit_eigenfrequency_design(
                    operator, rho, self.current.goal, self.current.options, self.control,
                )?;
                if !candidate.completed_and_rechecked() {
                    self.rejected = Some((tree, candidate));
                    return Ok(Stop::CandidateStopped);
                }
                self.pending = Some((tree, candidate));
            }
            let (tree, candidate) = self.pending.as_ref().expect("complete modal candidate");
            let summary = grid(candidate, tree)?;
            self.control.checkpoint("sdf3-modal-campaign-promote")?;
            let (tree, candidate) = self.pending.take().expect("complete modal candidate");
            self.current = candidate;
            self.tree = tree;
            self.grids.push(summary);
            self.checked = None;
            self.rejected = None;
            promoted += 1;
            self.control.checkpoint("sdf3-modal-campaign-promoted")?;
        }
    }
}
