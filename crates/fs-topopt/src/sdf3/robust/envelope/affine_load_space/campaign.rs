//! Automatic native refinement and reoptimization of an affine material design.
use super::design::{AffineLoadSpaceDesignGoal3, admit_model};
use super::grid::normalize;
use super::refinement::{AffineLoadSpaceRefinement3, AffineLoadSpaceRefinementOptions3};
use super::refit::{AffineLoadSpaceRefit3, AffineLoadSpaceRefitError3};
use super::response::reference::{self as response_reference, ReferenceAffineResponseLimits3};
use crate::sdf3::AdaptiveSdf3Elasticity;
use crate::sdf3::response::limits::design::MinimumMaterialObjective3;
use crate::sdf3_goal::GoalRefinementError3;
use crate::{SolveControl, SolveWork};
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use std::ops::ControlFlow;

pub use crate::sdf3::robust::envelope::limits::campaign::ComplianceCapCampaignError3 as AffineLoadSpaceCampaignError3;

/// The same physical builder supplies both diagnostic and optimization grids.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AffineLoadSpaceBuild3 {
    Probe,
    Refit,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AffineLoadSpaceMarking3 {
    /// Refine every active source cell; this is the default for a full load set.
    UniformActive,
    /// Optional diagnostic heuristic from two frozen maximizing directions.
    /// These indicators do not localize the complete affine model. Unresolved
    /// volume/inventory, or empty marks, always falls back to all active cells.
    FrozenWitness { theta: f64, max_marks: usize },
}

#[derive(Debug, Clone, Copy)]
pub struct AffineLoadSpaceCampaignOptions3 {
    /// Successful grid replacements over ALL calls to run.
    pub max_refits: usize,
    /// Candidate construction/refit attempts, including failures; never refunded.
    pub max_candidate_attempts: usize,
    /// Absolute normalized change of the full G, h and c model on the unit ball.
    pub change_tolerance: f64,
    /// Maximum normalized whole-model change of EACH original response limit.
    pub response_change_tolerance: f64,
    /// Absolute change in dilated numerical volume fraction.
    pub volume_tolerance: f64,
    /// Absolute inventory change divided by its declared objective scale.
    pub objective_tolerance: f64,
    pub marking: AffineLoadSpaceMarking3,
    pub enrichment: AffineLoadSpaceRefinementOptions3,
}
impl Default for AffineLoadSpaceCampaignOptions3 {
    fn default() -> Self {
        Self {
            max_refits: 2,
            max_candidate_attempts: 4,
            change_tolerance: 0.02,
            response_change_tolerance: 0.02,
            volume_tolerance: 1e-4,
            objective_tolerance: 1e-4,
            marking: AffineLoadSpaceMarking3::UniformActive,
            enrichment: AffineLoadSpaceRefinementOptions3::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AffineLoadSpaceCampaignStop3 {
    Paused,
    /// Current design passes both grids' physical bounds and all change gates.
    StableAndFeasible,
    /// The last installed design has been probed but remains unresolved.
    RefinementLimit,
    CandidateAttemptLimit,
    /// The attempted grid and its completed prefix remain in rejected().
    CandidateStopped,
}

/// Compact history; superseded full operators and field families are released.
#[derive(Debug, Clone, PartialEq)]
pub struct AffineLoadSpaceGridSummary3 {
    pub active_cells: usize,
    pub background_cells: usize,
    pub dilated_volume_fraction: f64,
    pub dilated_mass_kg: Option<f64>,
    pub objective: f64,
    pub worst_compliance: f64,
    pub response_constraint: Option<f64>,
    pub worst_response_residual: Option<f64>,
    pub work: SolveWork,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AffineLoadSpaceProbe3 {
    pub source_grid: usize,
    pub coarse_worst_compliance: f64,
    pub fine_worst_compliance: f64,
    pub normalized_model_change: f64,
    pub response_model_change: Option<f64>,
    pub volume_change: f64,
    pub objective_change: f64,
    pub mass_change_kg: Option<f64>,
    pub passed: bool,
    pub work: SolveWork,
}

fn invalid(message: &'static str) -> AffineLoadSpaceCampaignError3 {
    GoalRefinementError3::Invalid(message).into()
}

fn inventory_goal(
    goal: AffineLoadSpaceDesignGoal3,
) -> Result<(MinimumMaterialObjective3, f64), AffineLoadSpaceCampaignError3> {
    match goal {
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective,
            volume_scale,
        } => Ok((objective, volume_scale)),
        AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => Err(invalid(
            "automatic affine refinement requires a minimum-material compliance cap",
        )),
    }
}

fn grid<O: AdaptiveSdf3Elasticity>(
    current: &AffineLoadSpaceRefit3<O>,
    tree: &Octree3,
) -> Result<AffineLoadSpaceGridSummary3, AffineLoadSpaceCampaignError3> {
    let e = &current
        .recheck
        .as_ref()
        .expect("completed solve")
        .as_ref()
        .expect("rechecked point")
        .evaluation;
    let (objective, volume_scale) = inventory_goal(current.goal)?;
    let value = match objective {
        MinimumMaterialObjective3::VolumeFraction => normalize(
            e.physical.dilated_volume_fraction,
            volume_scale,
            "sdf3-affine-campaign-inventory",
        )?,
        MinimumMaterialObjective3::MassKg { scale_kg } => normalize(
            e.physical
                .dilated_mass_kg
                .ok_or_else(|| invalid("missing affine campaign kilograms"))?,
            scale_kg,
            "sdf3-affine-campaign-inventory",
        )?,
    };
    if value < 0.0 {
        return Err(invalid("negative affine campaign inventory"));
    }
    Ok(AffineLoadSpaceGridSummary3 {
        active_cells: current.study.cells(),
        background_cells: tree.leaves().len(),
        dilated_volume_fraction: e.physical.dilated_volume_fraction,
        dilated_mass_kg: e.physical.dilated_mass_kg,
        objective: value,
        worst_compliance: e.worst_compliance,
        response_constraint: e.responses.as_ref().map(|response| response.constraint),
        worst_response_residual: e.responses.as_ref().map(|response| response.worst_residual),
        work: current.work,
    })
}

fn check_tree<O: AdaptiveSdf3Elasticity>(
    op: &O,
    tree: &Octree3,
    control: &mut SolveControl<'_>,
) -> Result<(), AffineLoadSpaceCampaignError3> {
    for (i, cell) in op.adaptive().leaves().iter().enumerate() {
        if i % 64 == 0 {
            control.checkpoint("sdf3-affine-campaign-background")?;
        }
        if !tree.leaves().contains(cell) {
            return Err(invalid(
                "affine campaign operator differs from its background",
            ));
        }
    }
    Ok(())
}

// Preserve the original physical budget/cancellation reason through all native
// callbacks. A native generic Cancelled must not erase an exhausted solve budget.
fn controlled<T>(
    control: &mut SolveControl<'_>,
    stage: &'static str,
    action: impl FnOnce(&mut dyn FnMut() -> ControlFlow<()>) -> Result<T, AffineLoadSpaceCampaignError3>,
) -> Result<T, AffineLoadSpaceCampaignError3> {
    control.checkpoint(stage)?;
    let mut stopped = None;
    let result = action(&mut || match control.checkpoint(stage) {
        Ok(()) => ControlFlow::Continue(()),
        Err(error) => {
            stopped = Some(error);
            ControlFlow::Break(())
        }
    });
    if let Some(error) = stopped {
        return Err(error.into());
    }
    let value = result?;
    control.checkpoint(stage)?;
    Ok(value)
}

fn build_checked<O: AdaptiveSdf3Elasticity>(
    tree: &Octree3,
    kind: AffineLoadSpaceBuild3,
    control: &mut SolveControl<'_>,
    build: &mut impl FnMut(
        &Octree3,
        AffineLoadSpaceBuild3,
        &mut dyn FnMut() -> ControlFlow<()>,
    ) -> Result<O, GoalRefinementError3>,
) -> Result<O, AffineLoadSpaceCampaignError3> {
    let stage = match kind {
        AffineLoadSpaceBuild3::Probe => "sdf3-affine-campaign-probe-build",
        AffineLoadSpaceBuild3::Refit => "sdf3-affine-campaign-refit-build",
    };
    let operator = controlled(control, stage, |poll| {
        build(tree, kind, poll).map_err(Into::into)
    })?;
    check_tree(&operator, tree, control)?;
    Ok(operator)
}

fn inventory_resolved(
    report: &AffineLoadSpaceRefinement3,
    options: AffineLoadSpaceCampaignOptions3,
) -> bool {
    report.inventory.as_ref().is_some_and(|inventory| {
        inventory.volume_change.is_finite()
            && inventory.volume_change >= 0.0
            && inventory.volume_change <= options.volume_tolerance
            && inventory.objective_change.is_finite()
            && inventory.objective_change >= 0.0
            && inventory.objective_change <= options.objective_tolerance
    })
}

/// In-memory adaptive minimum-material design for one fixed affine load model.
/// Start with a reference affine solve, including its response family when used.
/// The load/response laws and builder preserve physical domain, supports, complete
/// native materials, amplitudes/correlations and purity throughout this owner.
/// The builder honors
/// its checkpoint and ONE cumulative geometry budget, including failed builds.
/// Native transfer checks material identity and genuine active-space enrichment.
///
/// Full probes inherit eroded AND dilated physical fields on uniformly refined
/// active cells. Candidates inherit RAW controls, rebuild native self-weight,
/// refilter on their own grid and run the same constrained optimizer. A candidate
/// is installed only after a fresh complete load-set recheck passes the physical
/// cap and every original response bound. No monotone self-weight restoration
/// or cross-grid descent is assumed.
/// Every installed grid receives a fresh full probe before a normal pause/limit.
///
/// Completed probes and pending candidates survive run segments. Failed refits
/// retain their own baseline, accepted prefix, outcome and recheck separately.
/// Promotion atomically replaces tree/operator/rechecked fields after a callback;
/// old full fields are released. This is not serialized checkpointing, global
/// design optimality, stationarity, or a continuum feasibility certificate.
pub struct AffineLoadSpaceCampaign3<'a, 'callback, O: AdaptiveSdf3Elasticity, B> {
    current: AffineLoadSpaceRefit3<O>,
    tree: Octree3,
    dead: ReferenceLoad3<'a>,
    basis: &'a [ReferenceLoad3<'a>],
    responses: Option<ReferenceAffineResponseLimits3<'a>>,
    options: AffineLoadSpaceCampaignOptions3,
    control: &'a mut SolveControl<'callback>,
    build: B,
    grids: Vec<AffineLoadSpaceGridSummary3>,
    probes: Vec<AffineLoadSpaceProbe3>,
    checked: Option<(Octree3, AffineLoadSpaceRefinement3)>,
    pending: Option<(Octree3, AffineLoadSpaceRefit3<O>)>,
    rejected: Option<(Octree3, AffineLoadSpaceRefit3<O>)>,
    attempts: usize,
    terminal: Option<AffineLoadSpaceCampaignStop3>,
}

impl<'a, 'callback, O: AdaptiveSdf3Elasticity, B> AffineLoadSpaceCampaign3<'a, 'callback, O, B>
where
    B: FnMut(
        &Octree3,
        AffineLoadSpaceBuild3,
        &mut dyn FnMut() -> ControlFlow<()>,
    ) -> Result<O, GoalRefinementError3>,
{
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        current: AffineLoadSpaceRefit3<O>,
        tree: Octree3,
        dead: ReferenceLoad3<'a>,
        basis: &'a [ReferenceLoad3<'a>],
        options: AffineLoadSpaceCampaignOptions3,
        control: &'a mut SolveControl<'callback>,
        build: B,
    ) -> Result<Self, AffineLoadSpaceCampaignError3> {
        Self::new_optional(current, tree, dead, basis, None, options, control, build)
    }

    /// Preserve one ordered family of reference response laws and physical
    /// bounds across every probe and candidate. Displacement observations are
    /// reintegrated on each actual grid; reaction modes stay native and live.
    /// The caller keeps their physical identity fixed, as for dead/basis laws.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_response_limits(
        current: AffineLoadSpaceRefit3<O>,
        tree: Octree3,
        dead: ReferenceLoad3<'a>,
        basis: &'a [ReferenceLoad3<'a>],
        responses: ReferenceAffineResponseLimits3<'a>,
        options: AffineLoadSpaceCampaignOptions3,
        control: &'a mut SolveControl<'callback>,
        build: B,
    ) -> Result<Self, AffineLoadSpaceCampaignError3> {
        Self::new_optional(
            current,
            tree,
            dead,
            basis,
            Some(responses),
            options,
            control,
            build,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_optional(
        current: AffineLoadSpaceRefit3<O>,
        tree: Octree3,
        dead: ReferenceLoad3<'a>,
        basis: &'a [ReferenceLoad3<'a>],
        responses: Option<ReferenceAffineResponseLimits3<'a>>,
        options: AffineLoadSpaceCampaignOptions3,
        control: &'a mut SolveControl<'callback>,
        build: B,
    ) -> Result<Self, AffineLoadSpaceCampaignError3> {
        let (objective, volume_scale) = inventory_goal(current.goal)?;
        let audit = current.gradient_check.as_ref().ok_or_else(|| {
            invalid("affine campaign requires the original successful gradient audit")
        })?;
        let (audited_objective, audited_volume_scale) = inventory_goal(audit.goal)?;
        if !audit.passed()
            || objective != audited_objective
            || (matches!(objective, MinimumMaterialObjective3::VolumeFraction)
                && volume_scale != audited_volume_scale)
        {
            return Err(invalid(
                "affine campaign goal differs from its optimized objective",
            ));
        }
        let constraint_policy = super::constraints::AffineLoadSpaceConstraintPolicy3 {
            responses: responses.map(|family| family.options),
            stress: None,
            sharpness: responses.map_or(40.0, |family| family.options.sharpness),
        };
        if audit.constraints != constraint_policy {
            return Err(invalid(
                "affine campaign must retain its optimized constraint policy",
            ));
        }
        if let Some(family) = responses {
            response_reference::admit(current.study.operator(), family)?;
        }
        current
            .study
            .require_material_body_recipe("sdf3-affine-campaign-material-recipe")?;
        let marking_valid = match options.marking {
            AffineLoadSpaceMarking3::UniformActive => true,
            AffineLoadSpaceMarking3::FrozenWitness { theta, max_marks } => {
                theta.is_finite() && theta > 0.0 && theta <= 1.0 && max_marks > 0
            }
        };
        if !current.completed_and_rechecked()
            || current.work != control.work()
            || !marking_valid
            || current.options.max_transfer_terms == 0
            || options.enrichment.max_transfer_terms == 0
            || ![
                options.change_tolerance,
                options.response_change_tolerance,
                options.volume_tolerance,
                options.objective_tolerance,
            ]
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0)
            || ![
                options.enrichment.numerical.residual_tolerance,
                options.enrichment.numerical.identity_tolerance,
            ]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0 && *v < 1.0)
        {
            return Err(invalid(
                "affine campaign needs a rechecked design, original work and valid options",
            ));
        }
        let point = current.recheck.as_ref().unwrap().as_ref().unwrap();
        let response_binding = match (responses, &point.evaluation.responses) {
            (None, None) => true,
            (Some(family), Some(report)) => {
                family.options == report.options
                    && family.limits.len() == report.witnesses.len()
                    && family
                        .limits
                        .iter()
                        .zip(&report.witnesses)
                        .all(|(limit, witness)| {
                            limit.bound == witness.bound
                                && limit.sense == witness.sense
                                && limit.scale == witness.scale
                        })
            }
            _ => false,
        };
        if !response_binding {
            return Err(invalid(
                "affine campaign original response bounds differ from its recheck",
            ));
        }
        let p = current.study.params();
        let q = current.params;
        if [p.penal, p.beta, p.eta, p.e_min] != [q.penal, q.beta, q.eta, q.e_min]
            || point.evaluation.policy != current.policy
            || point.evaluation.robustness != current.options.design.robustness
            || point.evaluation.physical.rho != point.rho
        {
            return Err(invalid("affine campaign model does not match its recheck"));
        }
        admit_model(
            &current.study,
            &point.rho,
            basis.len(),
            &current.policy,
            current.goal,
            current.options.design,
        )
        .map_err(AffineLoadSpaceRefitError3::from)?;
        check_tree(current.study.operator(), &tree, control)?;
        let grids = vec![grid(&current, &tree)?];
        Ok(Self {
            current,
            tree,
            dead,
            basis,
            responses,
            options,
            control,
            build,
            grids,
            probes: Vec::new(),
            checked: None,
            pending: None,
            rejected: None,
            attempts: 0,
            terminal: None,
        })
    }

    pub fn current(&self) -> &AffineLoadSpaceRefit3<O> {
        &self.current
    }
    pub fn tree(&self) -> &Octree3 {
        &self.tree
    }
    pub fn grids(&self) -> &[AffineLoadSpaceGridSummary3] {
        &self.grids
    }
    pub fn probes(&self) -> &[AffineLoadSpaceProbe3] {
        &self.probes
    }
    pub fn current_probe(&self) -> Option<&AffineLoadSpaceRefinement3> {
        self.checked.as_ref().map(|(_, report)| report)
    }
    pub fn pending(&self) -> Option<&(Octree3, AffineLoadSpaceRefit3<O>)> {
        self.pending.as_ref()
    }
    pub fn rejected(&self) -> Option<&(Octree3, AffineLoadSpaceRefit3<O>)> {
        self.rejected.as_ref()
    }
    pub fn work(&self) -> SolveWork {
        self.control.work()
    }
    pub fn candidate_attempts(&self) -> usize {
        self.attempts
    }
    pub fn into_current(self) -> (Octree3, AffineLoadSpaceRefit3<O>) {
        (self.tree, self.current)
    }

    fn candidate_tree(&mut self) -> Result<Octree3, AffineLoadSpaceCampaignError3> {
        let (uniform, report) = self.checked.as_ref().expect("complete full probe");
        let AffineLoadSpaceMarking3::FrozenWitness { theta, max_marks } = self.options.marking
        else {
            return Ok(uniform.clone());
        };
        if !inventory_resolved(report, self.options)
            || report.responses.as_ref().is_some_and(|responses| {
                !responses.stable_and_feasible(self.options.response_change_tolerance)
            })
        {
            return Ok(uniform.clone());
        }
        let marks = controlled(self.control, "sdf3-affine-campaign-mark", |poll| {
            report
                .mark(theta, max_marks, poll)
                .map_err(GoalRefinementError3::from)
                .map_err(Into::into)
        })?;
        if marks.marked.is_empty() {
            return Ok(uniform.clone());
        }
        controlled(
            self.control,
            "sdf3-affine-campaign-candidate-tree",
            |poll| self.tree.refined(&marks.marked, poll).map_err(Into::into),
        )
    }

    /// At most additional_refits successful promotions. Zero assesses the
    /// current grid without building a candidate. A positive segment also
    /// assesses its final installed design before pausing or reaching a cap.
    /// Only StableAndFeasible passes all numerical two-grid checks. Failures
    /// preserve the current design; completed pending promotions are reusable.
    pub fn run(
        &mut self,
        additional_refits: usize,
    ) -> Result<AffineLoadSpaceCampaignStop3, AffineLoadSpaceCampaignError3> {
        use AffineLoadSpaceCampaignStop3 as Stop;
        if let Some(stop) = self.terminal {
            return Ok(stop);
        }
        let mut promoted = 0;
        loop {
            if self.checked.is_none() {
                let tree = controlled(self.control, "sdf3-affine-campaign-probe-tree", |poll| {
                    self.tree
                        .refined(self.current.study.operator().adaptive().leaves(), poll)
                        .map_err(Into::into)
                })?;
                let operator = build_checked(
                    &tree,
                    AffineLoadSpaceBuild3::Probe,
                    self.control,
                    &mut self.build,
                )?;
                let point = self.current.recheck.as_ref().unwrap().as_ref().unwrap();
                let (objective, volume_scale) = inventory_goal(self.current.goal)?;
                let report = if let Some(family) = self.responses {
                    self.current
                        .study
                        .estimate_affine_load_space_enrichment_with_responses(
                            operator,
                            &point.rho,
                            self.dead,
                            self.basis,
                            self.current.options.design.robustness,
                            &self.current.policy,
                            objective,
                            volume_scale,
                            family,
                            self.options.enrichment,
                            self.control,
                        )?
                } else {
                    self.current
                        .study
                        .estimate_affine_load_space_enrichment_with_inventory(
                            operator,
                            &point.rho,
                            self.dead,
                            self.basis,
                            self.current.options.design.robustness,
                            &self.current.policy,
                            objective,
                            volume_scale,
                            self.options.enrichment,
                            self.control,
                        )?
                };
                let inventory = report.inventory.as_ref().expect("inventory requested");
                let passed = report.stable_and_feasible_with_response_tolerance(
                    self.options.change_tolerance,
                    self.options.response_change_tolerance,
                ) && inventory_resolved(&report, self.options);
                self.probes.push(AffineLoadSpaceProbe3 {
                    source_grid: self.grids.len() - 1,
                    coarse_worst_compliance: report.coarse.worst_compliance,
                    fine_worst_compliance: report.fine.worst_compliance,
                    normalized_model_change: report.normalized_model_change,
                    response_model_change: report
                        .responses
                        .as_ref()
                        .map(|r| r.maximum_model_change),
                    volume_change: inventory.volume_change,
                    objective_change: inventory.objective_change,
                    mass_change_kg: inventory.mass_change_kg,
                    passed,
                    work: report.work,
                });
                self.checked = Some((tree, report));
            }
            let stop = if self.probes.last().expect("complete probe").passed {
                Some(Stop::StableAndFeasible)
            } else if self.grids.len() - 1 >= self.options.max_refits {
                Some(Stop::RefinementLimit)
            } else {
                None
            };
            if let Some(stop) = stop {
                self.terminal = Some(stop);
                return Ok(stop);
            }
            if promoted == additional_refits {
                return Ok(Stop::Paused);
            }
            if self.pending.is_none() {
                if self.attempts >= self.options.max_candidate_attempts {
                    self.terminal = Some(Stop::CandidateAttemptLimit);
                    return Ok(Stop::CandidateAttemptLimit);
                }
                self.attempts += 1; // Charge before tree/build/refit; never refund.
                let tree = self.candidate_tree()?;
                let operator = build_checked(
                    &tree,
                    AffineLoadSpaceBuild3::Refit,
                    self.control,
                    &mut self.build,
                )?;
                let rho = &self.current.recheck.as_ref().unwrap().as_ref().unwrap().rho;
                let candidate = if let Some(family) = self.responses {
                    self.current
                        .study
                        .refit_reference_affine_load_space_with_responses(
                            operator,
                            rho,
                            self.dead,
                            self.basis,
                            &self.current.policy,
                            self.current.goal,
                            self.current.options,
                            family,
                            self.control,
                        )?
                } else {
                    self.current.study.refit_reference_affine_load_space(
                        operator,
                        rho,
                        self.dead,
                        self.basis,
                        &self.current.policy,
                        self.current.goal,
                        self.current.options,
                        self.control,
                    )?
                };
                if !candidate.completed_and_rechecked() {
                    self.rejected = Some((tree, candidate));
                    return Ok(Stop::CandidateStopped);
                }
                self.pending = Some((tree, candidate));
            }
            let (tree, candidate) = self.pending.as_ref().expect("complete rechecked candidate");
            let summary = grid(candidate, tree)?;
            self.control.checkpoint("sdf3-affine-campaign-promote")?;
            let (tree, candidate) = self.pending.take().expect("complete rechecked candidate");
            self.current = candidate;
            self.tree = tree;
            self.grids.push(summary);
            self.checked = None;
            self.rejected = None;
            promoted += 1;
            self.control.checkpoint("sdf3-affine-campaign-promoted")?;
        }
    }
}
