//! Native minimum-material and fixed-volume worst-compliance affine design.
//! Fixed reference preload and native self-weight share the common offset,
//! separately from the uncertain external-load basis.
use super::*;
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlReport, ProjectedAlStop, ProjectedAlWork,
};
use fs_topopt::sdf3::response::limits::design::MinimumMaterialObjective3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::{
    AffineLoadSpace3,
    constraints::{AffineLoadSpaceConstraintPolicy3, AffineLoadSpaceConstraints3},
    design::{
        AffineLoadSpaceDesign3, AffineLoadSpaceDesignCheckpoint3, AffineLoadSpaceDesignGoal3,
        AffineLoadSpaceDesignIteration3, AffineLoadSpaceDesignOptions3,
        AffineLoadSpaceDesignPoint3, AffineLoadSpaceGradientCheck3,
    },
    response::{AffineResponseLimits3, AffineResponseOptions3},
    stress::AffineStressLimit3,
};
#[path = "affine_load_space/output.rs"]
mod output;
#[path = "affine_load_space/resume.rs"]
mod resume;
#[path = "affine_load_space/response_binding.rs"]
mod response_binding;
#[path = "affine_load_space/grid_assessment.rs"]
pub(super) mod grid_assessment;
#[path = "affine_load_space/grid_refinement.rs"]
pub(super) mod grid_refinement;

const SCOPE: &str = "Estimated linear-elastic design for f=f_dead+B*rho_eroded+Fq, ||q||_2<=1. The declared goal either minimizes dilated material under a worst-compliance cap, or minimizes worst compliance under a dilated-volume cap. In the latter mode the compliance reference is a numerical normalization scale, not a physical cap. The common load combines optional fixed reference body/surface preload with native material self-weight at eroded occupancy. Fixed preload is not scaled by occupancy; pressure is inward-positive on the reference normal, never a follower load. Omitted preload is zero. Force columns retain physical amplitudes and correlations, not probabilities. A bounded quadratic-ball solve selects the complete affine compliance witness; its primal/dual gap and residuals are numerical diagnostics, not certified bounds. Eroded material supplies stiffness and self-weight; dilated material supplies inventory. With self-weight, erosion is not asserted to be the worst manufacturing realization. When strength is declared, each retained bulk quadrature point has its own adverse load direction; an unweighted p-norm bounds both relaxed and physical sampled stresses over the load ball in exact arithmetic. The declared allowable is not inferred from stiffness. Signed displacement and support-reaction limits use each inequality's own adverse load over the same full affine ball. Displacement means use each grid's reference cut volume; support forces are numerical tractions on the solid and moments are signed about the authored reference origin. The goal's physical cap, responses and strength share one conservative joint constraint and the same physical primal fields. Initial full-chain objective and constraint gradients are independently checked. Hard or unresolved load maximizers cannot supply a smooth gradient. Only accepted points passing every declared physical cap and the complete joint constraint enter the feasible archive, ranked by the actual objective; KKT values describe the last accepted point. Accepted results are retained before more optimization. An optional active-cell uniform enrichment re-solves the entire affine load family and assesses every original response and sampled strength while inheriting physical occupancy. Separate objective, inventory and physical change gates apply on both grids. An optional single bounded native refit optimizes inherited raw controls with the same goal on an unresolved enriched grid, with separate optimizer allowances and cumulative physical work. Its source design stays selected until a complete fresh candidate recheck passes, after which the promoted design receives its own finer-grid assessment. Failed assessment preserves its accepted checkpoint and charged work for a bounded retry. Cross-process resume restores accepted optimizer state under the same source and executable, verifies the native physical-problem and goal binding, rebuilds geometry and re-solves accepted and distinct feasible incumbent endpoints against the original work allowances. Earlier optimizer steps and the initial gradient experiment are not replayed. Failed restoration and uncheckpointed work lost after a crash cannot be durably charged; previous results remain intact. Individual kernels and ledger writes are indivisible; wall and cancellation limits are cooperative. Memory is an admission envelope. No continuum stress maximum, buckling, continuum-error, manufacturing-tolerance, global design-optimality, safety or experimental-validation claim.";

// The original authored mode retains its exact numerical objective. Every
// live solve, audit, restore and refinement instead uses Config::goal.
pub(super) const GOAL: AffineLoadSpaceDesignGoal3 = AffineLoadSpaceDesignGoal3::MinimumMaterial {
    objective: MinimumMaterialObjective3::VolumeFraction,
    volume_scale: 1.0,
};

#[derive(Clone, Copy, Debug)]
pub(super) struct Config {
    pub goal: AffineLoadSpaceDesignGoal3,
    /// Compliance cap for minimum material; numerical scale for worst compliance.
    pub cap_j: f64,
    pub options: AffineLoadSpaceDesignOptions3,
    pub stress: Option<AffineStressLimit3>,
    pub response_options: Option<AffineResponseOptions3>,
    pub joint_sharpness: f64,
}
impl Config {
    pub(super) fn constraints<'a>(
        self,
        responses: Option<AffineResponseLimits3<'a>>,
    ) -> AffineLoadSpaceConstraints3<'a> {
        AffineLoadSpaceConstraints3 {
            responses,
            stress: self.stress,
            sharpness: self.joint_sharpness,
        }
    }

    pub(super) fn constraint_policy(self) -> AffineLoadSpaceConstraintPolicy3 {
        AffineLoadSpaceConstraintPolicy3 {
            responses: self.response_options,
            stress: self.stress,
            sharpness: self.joint_sharpness,
        }
    }

    fn policy(self, dimension: usize) -> Result<AffineLoadSpace3> {
        AffineLoadSpace3::new(dimension, self.cap_j).map_err(|e| fail("cli-study-sdf3-input", e))
    }

    pub(super) fn validate(self, spec: &Spec) -> Result<()> {
        self.policy(spec.loads.len())?;
        match self.goal {
            AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit }
                if volume_limit.is_finite() && volume_limit > 0.0 && volume_limit <= 1.0 => {}
            AffineLoadSpaceDesignGoal3::MinimumMaterial {
                objective: MinimumMaterialObjective3::VolumeFraction, volume_scale: 1.0,
            } => {}
            _ => return Err(fail("cli-study-sdf3-input",
                "affine design requires minimum dilated volume or worst compliance with a finite dilated-volume limit in (0,1]")),
        }
        if self.response_options != spec.affine_responses.as_ref().map(|family| family.options) {
            return Err(fail("cli-study-sdf3-input", "affine response policy must match its authored physical observations"));
        }
        if let Some(family) = &spec.affine_responses {
            family.validate(spec)?;
        }
        if let Some(policy) = spec.affine_grid_assessment {
            policy.validate(spec)?;
        }
        if let Some(policy) = spec.affine_grid_refinement {
            policy.validate(spec)?;
        }
        self.options
            .optimizer
            .validate::<EvaluationStop>(spec.leaves)
            .map_err(|e| fail("cli-study-sdf3-input", e.to_string()))?;
        if spec.material_acceleration.is_none()
            || (spec.affine_grid_assessment.is_none() && spec.max_level != spec.level)
            || spec.schedule.len() != 1
            || spec.loads.iter().any(|(_, weight)| *weight != 1.0)
            || !(1..=1024).contains(&spec.updates)
            || !(0.001..=0.9).contains(&self.options.density_floor)
            || !(self.options.density_floor..=1.0).contains(&spec.density)
            || !(1e-12..=1e15).contains(&self.cap_j)
            || !(0.0..0.5).contains(&self.options.robustness.threshold_offset)
            || !(1..=10_000).contains(&self.options.optimizer.max_evaluations)
            || !(1e-10..=1e-2).contains(&self.options.optimizer.tolerance)
            || !(1e-6..=1e6).contains(&self.joint_sharpness)
        {
            return Err(fail(
                "cli-study-sdf3-input",
                "affine load-space design requires material-acceleration, an admitted grid and fixed model, unit basis weights, bounded updates/evaluations, positive compliance cap or reference scale, 0<=threshold-offset<0.5 and initial density above its floor",
            ));
        }
        if let Some(stress) = self.stress {
            let params = spec.schedule[0];
            let q = stress.options.relaxation_power;
            let floor = fs_topopt::heaviside(self.options.density_floor, params.beta,
                params.eta + self.options.robustness.threshold_offset);
            let turnover = q * params.e_min / ((params.penal - q) * (1.0 - params.e_min));
            if !(1e-12..=1e15).contains(&stress.limit)
                || !q.is_finite() || q < 1.0 || q >= params.penal
                || !(2.0..=64.0).contains(&stress.options.aggregation_power)
                || !(1..=1_000_000).contains(&stress.options.max_points)
                || spec.loads.len().checked_add(1) != Some(stress.options.max_cases)
                || !floor.is_finite() || floor <= 0.0
                || !(fs_math::det::pow(floor, params.penal) > turnover)
            {
                return Err(fail("cli-study-sdf3-input",
                    "affine-strength requires a finite positive stress cap, 1<=q<penal, 2<=p<=64, 1..=1000000 total stress point-fields covering the common offset and every basis column, and an eroded density floor above ersatz stress foldback"));
            }
        }
        Ok(())
    }

    /// Physical cap checks do not inherit the AL stopping tolerance.
    fn feasible(self, point: &AffineLoadSpaceDesignPoint3) -> bool {
        let e = &point.evaluation;
        self.policy(e.policy.dimension()).is_ok_and(|policy|
            e.robustness == self.options.robustness
                && e.constraints == self.constraint_policy()
                && point.physically_feasible(self.goal, &policy))
    }
}

type NumericalResult = std::result::Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>>;

struct State {
    accepted: AffineLoadSpaceDesignPoint3,
    best: Option<AffineLoadSpaceDesignPoint3>,
    history: Vec<AffineLoadSpaceDesignIteration3>,
    gradient: AffineLoadSpaceGradientCheck3,
    report: Option<ProjectedAlReport>,
    error: Option<ProjectedAlError<EvaluationStop>>,
    optimizer: ProjectedAlWork,
    spent: Spent,
    status: &'static str,
    optimization_status: &'static str,
    assessment: Option<grid_assessment::Outcome>,
    refinement: Option<grid_refinement::Context>,
    refinement_recheck: Option<std::result::Result<AffineLoadSpaceDesignPoint3, EvaluationStop>>,
    refinement_failure: Option<String>,
    checkpoint: AffineLoadSpaceDesignCheckpoint3,
}
impl State {
    fn capture(
        design: &AffineLoadSpaceDesign3<'_, '_, AdaptiveSolveSpace3>,
        result: NumericalResult,
    ) -> Self {
        Self {
            accepted: design.accepted().clone(),
            best: design.best_feasible().cloned(),
            history: design.history().to_vec(),
            gradient: design.gradient_check().clone(),
            report: result.as_ref().ok().cloned(),
            error: result.err(),
            optimizer: design.optimizer_work(),
            spent: Spent::default(),
            status: "checkpointed",
            optimization_status: "checkpointed",
            assessment: None,
            refinement: None,
            refinement_recheck: None,
            refinement_failure: None,
            checkpoint: design.checkpoint(),
        }
    }
    fn selected(&self) -> &AffineLoadSpaceDesignPoint3 {
        self.best.as_ref().unwrap_or(&self.accepted)
    }
    fn terminal(&self, config: Config, gate: &CancelGate, expired: bool) -> &'static str {
        if gate.is_requested() {
            return "cancelled";
        }
        if expired {
            return "budget-exhausted";
        }
        if let Some(error) = &self.error {
            return match error {
                ProjectedAlError::Cancelled
                | ProjectedAlError::Evaluation(EvaluationStop::Cancelled) => "cancelled",
                ProjectedAlError::Evaluation(
                    EvaluationStop::LinearBudget { .. } | EvaluationStop::TotalBudget { .. },
                ) => "budget-exhausted",
                _ => "numerical-failure",
            };
        }
        match self.report.as_ref().map(|r| r.stop) {
            Some(ProjectedAlStop::Converged) if config.feasible(&self.accepted) => "completed",
            Some(ProjectedAlStop::Stalled | ProjectedAlStop::PenaltyLimit) => "no-feasible-descent",
            _ => "budget-exhausted",
        }
    }
}

struct Computation {
    study: CutDensityStudy3<AdaptiveSolveSpace3>,
    state: State,
    tree: Octree3,
}

fn numerical_failure(
    error: ProjectedAlError<EvaluationStop>,
    gate: &CancelGate,
    expired: bool,
) -> Failure {
    let budget = matches!(
        &error,
        ProjectedAlError::Evaluation(
            EvaluationStop::LinearBudget { .. } | EvaluationStop::TotalBudget { .. }
        )
    );
    let cancelled = matches!(
        &error,
        ProjectedAlError::Cancelled | ProjectedAlError::Evaluation(EvaluationStop::Cancelled)
    );
    Failure {
        code: "cli-study-sdf3-affine-load-space",
        message: error.to_string(),
        exit: if gate.is_requested() {
            exit::CANCELLED
        } else if expired || budget {
            exit::BUDGET
        } else if cancelled {
            exit::CANCELLED
        } else {
            exit::REFUSED
        },
    }
}

fn compute_observed(
    spec: &Spec,
    gate: &CancelGate,
    limit: usize,
    recovery: Option<&resume::Recovery>,
    preparation_wall_s: f64,
    observed: impl FnMut(&CutDensityStudy3<AdaptiveSolveSpace3>, &State) -> Result<()>,
) -> Result<Computation> {
    compute_on_grid(spec, gate, limit, recovery, preparation_wall_s, None, observed)
}

#[allow(clippy::too_many_arguments)]
fn compute_on_grid(
    spec: &Spec,
    gate: &CancelGate,
    limit: usize,
    recovery: Option<&resume::Recovery>,
    preparation_wall_s: f64,
    prepared: Option<grid_refinement::Prepared>,
    mut observed: impl FnMut(&CutDensityStudy3<AdaptiveSolveSpace3>, &State) -> Result<()>,
) -> Result<Computation> {
    let config = spec.affine_load_space.ok_or_else(|| {
        fail(
            "cli-study-sdf3-input",
            "missing affine load-space optimizer",
        )
    })?;
    config.validate(spec)?;
    let policy = config.policy(spec.loads.len())?;
    let context = prepared.as_ref().map(|p| p.context.clone())
        .or_else(|| recovery.and_then(|r| r.refinement.clone()));
    let prior = prepared.as_ref().map(|p| p.spent)
        .unwrap_or_else(|| recovery.map_or(Spent::default(), |r| r.spent)).add(
        SolveWork::default(),
        QuadratureWork3::default(),
        preparation_wall_s,
    )?;
    prior.require_remaining(spec)?;
    let start = Instant::now();
    let poll = || {
        if gate.is_requested() || prior.wall_s + start.elapsed().as_secs_f64() >= spec.wall_s {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    if poll().is_break() {
        return Err(Failure {
            code: "cli-study-sdf3-stopped",
            message: "affine load-space study stopped before geometry".into(),
            exit: if gate.is_requested() {
                exit::CANCELLED
            } else {
                exit::BUDGET
            },
        });
    }
    let tree = match &prepared {
        Some(prepared) => prepared.tree.clone(),
        None => grid_refinement::tree(spec, context.as_ref(), &poll)
            .map_err(|e| grid_refinement::numerical(e, gate, poll().is_break()))?,
    };
    let legacy = Domain {
        height: spec.height,
        curvature: spec.curvature,
    };
    let physical = PhysicalDomain {
        bounds: spec.bounds,
        height: spec.height,
        curvature: spec.curvature,
    };
    let domain: &dyn CutSdf3 = if spec.physical { &physical } else { &legacy };
    let bounds = HexCell::try_new(spec.bounds.0, spec.bounds.1)
        .map_err(|e| fail("cli-study-sdf3-domain", e.to_string()))?;
    let material = loading::material(spec)?;
    let mut cp = |_| poll();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            max_boxes: spec.boxes - prior.geometry.boxes,
            max_points: spec.points - prior.geometry.points,
            ..Default::default()
        },
        &mut cp,
    )
    .map_err(|e| fail("cli-study-sdf3-geometry", e.to_string()))?;
    let (mut study, raw) = if let Some(prepared) = prepared {
        (prepared.study, prepared.inherited_rho)
    } else {
        let operator = (|| -> std::result::Result<AdaptiveSolveSpace3, GoalRefinementError3> {
            let mut build = |grid: &Octree3, checkpoint: &mut dyn FnMut() -> ControlFlow<()>| {
                if checkpoint().is_break() {
                    return Err(EvaluationStop::Cancelled.into());
                }
                loading::build_operator(spec, bounds, grid, domain, &material, &mut quadrature)
                    .map_err(GoalRefinementError3::from)
            };
            let corrections = grid_refinement::corrections(spec, &tree, &mut build, &mut || poll())?;
            let operator = build(&tree, &mut || poll())?;
            let refs: Vec<_> = corrections.iter().collect();
            spec.solver.wrap(operator, &refs, || poll())
        })().map_err(|e| Failure {
            code: "cli-study-sdf3-geometry",
            message: e.to_string(),
            exit: if gate.is_requested() {
                exit::CANCELLED
            } else if poll().is_break() || geometry_budget(&e) || solver::setup_budget(&e) {
                exit::BUDGET
            } else {
                exit::REFUSED
            },
        })?;
        let study = regions::bind(
            CutDensityStudy3::new(operator, spec.radius, spec.schedule[0]), spec,
        )?;
        let study = loading::bind_material_acceleration(study, spec, gate, || poll())?;
        let raw = vec![spec.density; study.cells()];
        (study, raw)
    };
    let geometry = quadrature.work();
    let mut cp = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: spec.linear - prior.linear.linear_iterations,
            per_solve_iterations: spec.per_solve,
        },
        &mut cp,
    );
    let mut state = loading::with_affine_laws(spec, |dead_law, laws| -> Result<State> {
        let assemble = |law| study.operator().elasticity().reference_load(law, &poll)
            .map_err(|e| Failure {
                code: "cli-study-sdf3-affine-load-space",
                message: e.to_string(),
                exit: if gate.is_requested() {
                    exit::CANCELLED
                } else if poll().is_break() {
                    exit::BUDGET
                } else {
                    exit::REFUSED
                },
            });
        let dead = assemble(dead_law)?;
        let forces = laws.iter().map(|law| assemble(law.load))
            .collect::<Result<Vec<_>>>()?;
        let basis: Vec<_> = forces.iter().map(Vec::as_slice).collect();
        // The changing material body force is assembled by the native owner
        // once at eroded occupancy. Fixed preload and self-weight must never
        // be baked into each uncertain basis column.
        let response_fields = response_binding::assemble(spec, study.operator().elasticity(), &poll)
            .map_err(|e| Failure {
                code: "cli-study-sdf3-affine-response",
                message: e.to_string(),
                exit: if gate.is_requested() {
                    exit::CANCELLED
                } else if poll().is_break() {
                    exit::BUDGET
                } else {
                    exit::REFUSED
                },
            })?;
        response_binding::with_limits(response_fields.as_ref(), |responses| {
            let mut design = match recovery {
                None => AffineLoadSpaceDesign3::new_with_constraints(
                    &mut study,
                    &dead,
                    &basis,
                    &raw,
                    policy,
                    config.goal,
                    config.options,
                    config.constraints(responses),
                    &mut control,
                ),
                Some(recovery) => AffineLoadSpaceDesign3::restore_with_constraints(
                    &mut study,
                    &dead,
                    &basis,
                    &mut recovery.checkpoint.clone(),
                    policy,
                    config.goal,
                    config.options,
                    config.constraints(responses),
                    &mut control,
                ),
            }
            .map_err(|e| numerical_failure(e, gate, poll().is_break()))?;
            let result = design.run(0);
            let mut state = State::capture(&design, result);
            state.refinement = context.clone();
            if let Some(recovery) = recovery.filter(|r| r.optimizer_finished(spec)) {
                if let (Some(report), Some(stop)) = (&mut state.report, recovery.optimizer_stop) {
                    // Restoring an endpoint is not a new search that can
                    // replace its original Stalled/PenaltyLimit termination.
                    report.stop = stop;
                }
            }
            let updates = if recovery.is_some_and(|r| r.optimizer_finished(spec)) { 0 }
                else { limit.min(spec.updates - design.optimizer_work().iterations) };
            for _ in 0..updates {
                if state.error.is_some()
                    || state
                        .report
                        .as_ref()
                        .is_some_and(|r| r.stop != ProjectedAlStop::IterationLimit)
                {
                    break;
                }
                state.spent = prior.add(design.work(), geometry, start.elapsed().as_secs_f64())?;
                observed(design.study(), &state)?;
                let result = design.run(1);
                state = State::capture(&design, result);
                state.refinement = context.clone();
            }
            match context.as_ref().map(|c| c.phase) {
                Some(grid_refinement::Phase::Candidate) if state.error.is_none()
                    && state.report.as_ref().is_some_and(|r| r.stop == ProjectedAlStop::Converged
                        || (r.stop == ProjectedAlStop::IterationLimit
                            && state.optimizer.iterations == spec.updates)) => {
                    // Save the complete accepted candidate before its independent
                    // physical recheck. A failed check never installs its fields.
                    state.spent = prior.add(design.work(), geometry, start.elapsed().as_secs_f64())?;
                    state.optimization_status = recovery.filter(|r| r.optimizer_finished(spec))
                        .map_or_else(|| state.terminal(config, gate, state.spent.wall_s >= spec.wall_s),
                            |r| r.optimization_status);
                    observed(design.study(), &state)?;
                    state.refinement_recheck = Some(design.recheck_best().and_then(|point| point.ok_or(
                        EvaluationStop::Breakdown { stage: "cli-sdf3-affine-refinement-no-feasible-candidate" }
                    )));
                    state.checkpoint = design.checkpoint();
                }
                Some(grid_refinement::Phase::Promoted) => {
                    // Restore independently solved both native retained endpoints;
                    // keep that evidence without replaying an optimizer experiment.
                    if !config.feasible(state.selected()) {
                        return Err(fail("cli-study-sdf3-refinement",
                            "restored promoted design failed its complete physical caps"));
                    }
                    state.refinement_recheck = Some(Ok(state.selected().clone()));
                }
                _ => {}
            }
            Ok(state)
        })
    })?;
    state.spent = prior.add(control.work(), geometry, start.elapsed().as_secs_f64())?;
    state.optimization_status = recovery.filter(|r| r.optimizer_finished(spec))
        .map_or_else(|| state.terminal(config, gate, state.spent.wall_s >= spec.wall_s),
            |r| r.optimization_status);
    state.status = state.optimization_status;
    if let Some(Err(error)) = &state.refinement_recheck {
        state.status = if gate.is_requested() { "cancelled" }
            else if state.spent.wall_s >= spec.wall_s
                || matches!(error, EvaluationStop::LinearBudget { .. } | EvaluationStop::TotalBudget { .. })
            { "budget-exhausted" }
            else if matches!(error, EvaluationStop::Cancelled) { "cancelled" }
            else { "numerical-failure" };
        state.refinement_failure = Some(error.to_string());
    }
    if let Some(assessment) = spec.affine_grid_assessment.filter(|_| context.as_ref()
        .is_none_or(|c| c.phase == grid_refinement::Phase::Promoted)) {
        // Retain the complete accepted AL slice before any enriched physics.
        // Pending assessment is not a claim of full study completion.
        state.status = "checkpointed";
        observed(&study, &state)?;
        grid_assessment::attach(assessment, spec, &tree, &mut study, &mut state,
            &mut quadrature, &mut control, gate, &poll);
        state.spent = prior.add(control.work(), quadrature.work(), start.elapsed().as_secs_f64())?;
    }
    Ok(Computation { study, state, tree })
}

pub(super) fn drive(
    spec: &Spec,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
) -> Result<Outcome> {
    if spec.affine_grid_refinement.is_some() {
        grid_refinement::drive(spec, ledger, cap, gate, None)
    } else {
        resume::drive(spec, ledger, cap, gate, None)
    }
}

pub(super) fn resume(
    ledger: &Ledger,
    old: &Loaded,
    cap: Option<usize>,
    gate: &CancelGate,
) -> Result<Outcome> {
    let source = linked(ledger, &old.value, "source", "study-source")?;
    let source = std::str::from_utf8(&source)
        .map_err(|e| fail("cli-study-sdf3-affine-load-space-checkpoint", e.to_string()))?;
    let spec = spec::parse(source)?;
    if spec.affine_grid_refinement.is_some() {
        grid_refinement::drive(&spec, ledger, cap, gate, Some(old))
    } else {
        resume::drive(&spec, ledger, cap, gate, Some(old))
    }
}

#[cfg(test)]
#[path = "affine_load_space/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "affine_load_space/strength_tests.rs"]
mod strength_tests;

#[cfg(test)]
#[path = "affine_load_space/response_tests.rs"]
mod response_tests;

#[cfg(test)]
#[path = "affine_load_space/grid_assessment_tests.rs"]
mod grid_assessment_tests;

#[cfg(test)]
#[path = "affine_load_space/grid_refinement_tests.rs"]
mod grid_refinement_tests;

#[cfg(test)]
#[path = "affine_load_space/compliance_tests.rs"]
mod compliance_tests;

#[cfg(test)]
#[path = "affine_load_space/dead_load_tests.rs"]
mod dead_load_tests;
