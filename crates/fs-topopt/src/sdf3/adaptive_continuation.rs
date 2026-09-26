//! DWR-guided background refinement BETWEEN gradient-checked density stages.
//! The implicit boundary remains unchanged. An enriched probe inherits physical
//! stiffness; an optimization proposal instead refilters inherited RAW density.
//! These are different experiments and their compliances are not compared.

use std::ops::ControlFlow;

use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::octree3::{Octree3, OctreeError3};
use fs_dwr::elasticity3::GoalMarking3;

use super::continuation::controlled_gradient_checked_sdf3_continuation;
use super::{AdaptiveSdf3Elasticity, CutDensityStudy3};
use crate::pipeline::LoadCase;
use crate::sdf3_goal::{
    ComplianceRefinement3, GoalReferenceLoad3, GoalRefinementError3, GoalRefinementOptions3,
};
use crate::{
    ContinuationTermination, EvaluationStop, GradientCheckOptions, MultiLoadContinuationReport,
    MultiLoadOcOptions, SimpParams, SolveControl,
};

/// Policies for the existing numerical stages, enrichment and deterministic
/// Dörfler marking. Octree level/leaf caps remain attached to the input tree.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveContinuationOptions3 {
    pub optimization: MultiLoadOcOptions,
    pub gradient: GradientCheckOptions,
    pub enrichment: GoalRefinementOptions3,
    /// Fraction of summed absolute local goal residual contributions to seek.
    pub marking_fraction: f64,
    /// A smaller cap can produce a valid mark set without meeting the fraction.
    pub max_marks: usize,
}

impl Default for AdaptiveContinuationOptions3 {
    fn default() -> Self {
        Self {
            optimization: MultiLoadOcOptions::default(),
            gradient: GradientCheckOptions::default(),
            enrichment: GoalRefinementOptions3::default(),
            marking_fraction: 0.5,
            max_marks: 2,
        }
    }
}

/// A refinement refusal retains the previous accepted tree, model and fields.
#[derive(Debug, Clone, PartialEq)]
pub enum AdaptiveContinuationError3 {
    Background(OctreeError3),
    Goal(GoalRefinementError3),
    /// No positive local residual signal; this does not establish accuracy.
    NoRefinementSignal,
}

impl From<OctreeError3> for AdaptiveContinuationError3 {
    fn from(error: OctreeError3) -> Self {
        if error == OctreeError3::Cancelled {
            Self::Goal(EvaluationStop::Cancelled.into())
        } else {
            Self::Background(error)
        }
    }
}
impl From<GoalRefinementError3> for AdaptiveContinuationError3 {
    fn from(error: GoalRefinementError3) -> Self {
        Self::Goal(error)
    }
}
impl From<EvaluationStop> for AdaptiveContinuationError3 {
    fn from(error: EvaluationStop) -> Self {
        Self::Goal(error.into())
    }
}

/// A complete goal experiment on the preceding accepted design. The probe
/// inherits that design's stiffness; the proposed stage restores a new filtered
/// design and model. `installed` distinguishes the two even on interruption.
#[derive(Debug, Clone)]
pub struct ContinuationRefinement3 {
    pub destination_stage: usize,
    pub source_background_cells: usize,
    pub target_background_cells: usize,
    pub source_active_cells: usize,
    pub target_active_cells: usize,
    pub estimate: ComplianceRefinement3,
    pub marking: GoalMarking3,
    /// True only after the complete proposed optimization stage is accepted.
    pub installed: bool,
}

/// The retained continuation prefix uses the same reports as fixed-background
/// continuation. Refinement failures set `EvaluationStopped` and retain their
/// typed reason separately; no refinement creates a continuum error bound.
#[derive(Debug, Clone)]
pub struct AdaptiveContinuationReport3 {
    pub continuation: MultiLoadContinuationReport,
    pub refinements: Vec<ContinuationRefinement3>,
    pub refinement_error: Option<AdaptiveContinuationError3>,
}

fn poll(control: &mut SolveControl<'_>, stage: &'static str) -> ControlFlow<()> {
    match control.checkpoint(stage) {
        Ok(()) => ControlFlow::Continue(()),
        Err(_) => ControlFlow::Break(()),
    }
}

fn check_tree<O: AdaptiveSdf3Elasticity>(
    operator: &O,
    tree: &Octree3,
) -> Result<(), GoalRefinementError3> {
    if operator
        .adaptive()
        .leaves()
        .iter()
        .all(|cell| tree.leaves().contains(cell))
    {
        Ok(())
    } else {
        Err(GoalRefinementError3::Invalid(
            "operator active cells differ from the supplied background",
        ))
    }
}

fn build_operator<O>(
    tree: &Octree3,
    control: &mut SolveControl<'_>,
    stage: &'static str,
    build: &mut impl FnMut(
        &Octree3,
        &mut dyn FnMut() -> ControlFlow<()>,
    ) -> Result<O, GoalRefinementError3>,
) -> Result<O, GoalRefinementError3> {
    control.checkpoint(stage)?;
    let mut stopped = false;
    let result = build(tree, &mut || {
        let flow = poll(control, stage);
        stopped |= flow.is_break();
        flow
    });
    if stopped {
        return Err(EvaluationStop::Cancelled.into());
    }
    control.checkpoint(stage)?;
    result
}

fn reference_stage<O: AdaptiveSdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    loads: &[GoalReferenceLoad3<'_>],
    rho: &[f64],
    params: SimpParams,
    options: AdaptiveContinuationOptions3,
    control: &mut SolveControl<'_>,
) -> Result<MultiLoadContinuationReport, GoalRefinementError3> {
    let mut forces = Vec::with_capacity(loads.len());
    for load in loads {
        control.checkpoint("sdf3-adaptive-reference-load")?;
        forces.push(
            study
                .operator()
                .adaptive()
                .reference_load(load.load, || poll(control, "sdf3-adaptive-reference-load"))?,
        );
    }
    let nodal: Vec<_> = loads
        .iter()
        .zip(&forces)
        .map(|(load, force)| LoadCase {
            force,
            weight: load.weight,
        })
        .collect();
    Ok(controlled_gradient_checked_sdf3_continuation(
        study,
        &nodal,
        rho,
        &[params],
        options.optimization,
        options.gradient,
        control,
    ))
}

/// Run a material/projection schedule with actual compliance-DWR background
/// refinement before EVERY stage after the first. The builder receives the
/// proposed balanced tree and a cooperative checkpoint; it must build the SAME
/// implicit domain, reference material, clamp and integration policy each time.
/// It may return a bare, Jacobi or hierarchy-backed operator. The retained
/// filter radius is reused automatically. No SDF surface meshing occurs.
///
/// Reference loads are reintegrated on each grid. The probe inherits PHYSICAL
/// stiffness scales, solves the real enriched compliance family and marks its
/// original cells by weighted absolute residual contributions. The candidate
/// inherits only RAW density and receives its own volume restoration, runtime
/// compliance/volume gradient gate and independently solved baseline.
///
/// Every proposed tree and study stays separate until its complete stage and
/// final checkpoint succeed. Cancellation, caps, a failed gradient gate or
/// numerical failure in a refined proposal leaves the previous accepted tree,
/// scales, model, design and fields intact. Stage zero uses the fixed-grid
/// driver's accepted-prefix semantics. All work, including rejected proposals,
/// remains charged to the supplied control. Geometry/quadrature limits and
/// cumulative geometry work belong to the builder; callbacks are indivisible.
///
/// This is goal-guided BACKGROUND adaptivity, not a moving implicit boundary,
/// continuum certificate, mesh-independent optimum or cross-model descent.
#[allow(clippy::too_many_arguments)]
pub fn controlled_adaptive_sdf3_continuation<O: AdaptiveSdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    tree: &mut Octree3,
    loads: &[GoalReferenceLoad3<'_>],
    rho0: &[f64],
    schedule: &[SimpParams],
    options: AdaptiveContinuationOptions3,
    control: &mut SolveControl<'_>,
    build: impl FnMut(
        &Octree3,
        &mut dyn FnMut() -> ControlFlow<()>,
    ) -> Result<O, GoalRefinementError3>,
) -> AdaptiveContinuationReport3 {
    match controlled_adaptive_sdf3_continuation_observed(
        study, tree, loads, rho0, schedule, options, control, build,
        |_, _, _| Ok::<(), std::convert::Infallible>(()),
    ) {
        Ok(report) => report,
        Err(never) => match never {},
    }
}

/// Run the same numerical driver with a fallible accepted-stage observer.
///
/// The observer receives the installed study, complete background and matching
/// report AFTER each whole stage succeeds, and BEFORE any work on the next
/// stage. It is never called for an unassessed baseline or a rejected proposal.
/// The report's work counters include all work through that boundary. A caller
/// can commit a durable checkpoint here without cloning the finite-element
/// operator or publishing fields from a different mesh.
///
/// An observer error returns immediately and unchanged: no later geometry,
/// solve, or observer is run. Persistence is the caller's responsibility; an
/// observer failure does not roll back the already accepted in-memory study.
/// Physics/budget/cancellation stops still return `Ok(report)` with their typed
/// stop cause. A prefix of the schedule may be supplied for bounded execution;
/// the caller must distinguish completing that prefix from the full schedule.
#[allow(clippy::too_many_arguments)]
pub fn controlled_adaptive_sdf3_continuation_observed<O: AdaptiveSdf3Elasticity, E>(
    study: &mut CutDensityStudy3<O>,
    tree: &mut Octree3,
    loads: &[GoalReferenceLoad3<'_>],
    rho0: &[f64],
    schedule: &[SimpParams],
    options: AdaptiveContinuationOptions3,
    control: &mut SolveControl<'_>,
    mut build: impl FnMut(
        &Octree3,
        &mut dyn FnMut() -> ControlFlow<()>,
    ) -> Result<O, GoalRefinementError3>,
    mut accepted: impl FnMut(
        &CutDensityStudy3<O>,
        &Octree3,
        &AdaptiveContinuationReport3,
    ) -> Result<(), E>,
) -> Result<AdaptiveContinuationReport3, E> {
    assert!(
        !schedule.is_empty(),
        "continuation requires at least one stage"
    );
    for params in schedule {
        params.assert_valid();
    }
    options.gradient.assert_valid();
    assert!(
        options.marking_fraction.is_finite()
            && options.marking_fraction > 0.0
            && options.marking_fraction <= 1.0
            && options.max_marks > 0,
        "invalid adaptive marking policy"
    );
    assert!(
        !loads.is_empty()
            && loads.len() <= options.enrichment.max_load_cases
            && loads.iter().any(|load| load.weight > 0.0)
            && loads
                .iter()
                .all(|load| load.weight.is_finite() && load.weight >= 0.0),
        "invalid reference load family"
    );
    let mut report = AdaptiveContinuationReport3 {
        continuation: MultiLoadContinuationReport {
            last: None,
            params: study.params(),
            stages: Vec::new(),
            termination: ContinuationTermination::ScheduleComplete,
            stopped_stage: None,
            evaluation_stop: None,
            rejected_gradient_check: None,
            work: control.work(),
        },
        refinements: Vec::new(),
        refinement_error: None,
    };
    let mut rho = rho0.to_vec();
    for (stage, &params) in schedule.iter().enumerate() {
        let result = (|| -> Result<MultiLoadContinuationReport, AdaptiveContinuationError3> {
            control.checkpoint("sdf3-adaptive-stage")?;
            check_tree(study.operator(), tree)?;
            if stage == 0 {
                return Ok(reference_stage(
                    study, loads, &rho, params, options, control,
                )?);
            }
            let accepted = report
                .continuation
                .last
                .as_ref()
                .expect("preceding stage solved");
            let source = study.operator().adaptive();
            let probe_tree = tree.refined(source.leaves(), || {
                poll(control, "sdf3-adaptive-enrichment-tree")
            })?;
            let probe = build_operator(
                &probe_tree,
                control,
                "sdf3-adaptive-enrichment-geometry",
                &mut build,
            )?;
            check_tree(&probe, &probe_tree)?;
            let estimate = study.estimate_reference_compliance_enrichment(
                probe.into_adaptive(),
                loads,
                &accepted.displacements,
                options.enrichment,
                control,
            )?;
            let marking = estimate
                .mark(options.marking_fraction, options.max_marks, || {
                    poll(control, "sdf3-adaptive-mark")
                })
                .map_err(GoalRefinementError3::from)?;
            if marking.marked.is_empty() {
                return Err(AdaptiveContinuationError3::NoRefinementSignal);
            }
            let next_tree = tree.refined(&marking.marked, || {
                poll(control, "sdf3-adaptive-marked-tree")
            })?;
            let operator = build_operator(
                &next_tree,
                control,
                "sdf3-adaptive-candidate-geometry",
                &mut build,
            )?;
            check_tree(&operator, &next_tree)?;
            let transfer = AdaptiveTransfer3::new(
                source,
                operator.adaptive(),
                options.enrichment.max_transfer_terms,
                || poll(control, "sdf3-adaptive-transfer"),
            )
            .map_err(GoalRefinementError3::from)?;
            if !operator
                .adaptive()
                .leaves()
                .iter()
                .zip(transfer.parents())
                .any(|(leaf, &parent)| leaf.level() > source.leaves()[parent].level())
            {
                return Err(GoalRefinementError3::Invalid(
                    "marked proposal did not enrich active space",
                )
                .into());
            }
            let mut seed = Vec::with_capacity(operator.cells());
            for &parent in transfer.parents() {
                control.checkpoint("sdf3-adaptive-raw-transfer")?;
                seed.push(accepted.rho[parent]);
            }
            report.refinements.push(ContinuationRefinement3 {
                destination_stage: stage,
                source_background_cells: tree.leaves().len(),
                target_background_cells: next_tree.leaves().len(),
                source_active_cells: study.cells(),
                target_active_cells: operator.cells(),
                estimate,
                marking,
                installed: false,
            });
            let mut candidate =
                CutDensityStudy3::new(operator, study.filter_radius, study.params());
            let next = reference_stage(&mut candidate, loads, &seed, params, options, control)?;
            if next.termination == ContinuationTermination::ScheduleComplete && next.last.is_some()
            {
                control.checkpoint("sdf3-adaptive-accept")?;
                *study = candidate;
                *tree = next_tree;
                report
                    .refinements
                    .last_mut()
                    .expect("completed proposal")
                    .installed = true;
            }
            Ok(next)
        })();
        match result {
            Ok(mut next) => {
                let complete = next.termination == ContinuationTermination::ScheduleComplete;
                if stage == 0 || complete {
                    if let Some(last) = next.last.take() {
                        rho.clone_from(&last.rho);
                        report.continuation.params = next.params;
                        report.continuation.last = Some(last);
                        for mut row in next.stages {
                            row.stage = stage;
                            report.continuation.stages.push(row);
                        }
                    }
                }
                if !complete {
                    report.continuation.termination = next.termination;
                    report.continuation.stopped_stage = Some(stage);
                    report.continuation.evaluation_stop = next.evaluation_stop;
                    report.continuation.rejected_gradient_check = next.rejected_gradient_check;
                    break;
                }
                report.continuation.work = control.work();
                accepted(study, tree, &report)?;
            }
            Err(error) => {
                report.continuation.termination = ContinuationTermination::EvaluationStopped;
                report.continuation.stopped_stage = Some(stage);
                if let AdaptiveContinuationError3::Goal(GoalRefinementError3::Evaluation(stop)) =
                    &error
                {
                    report.continuation.evaluation_stop = Some(stop.clone());
                }
                report.refinement_error = Some(error);
                break;
            }
        }
    }
    report.continuation.work = control.work();
    Ok(report)
}

#[cfg(test)]
#[path = "adaptive_continuation/checkpoint_tests.rs"]
mod checkpoint_tests;
