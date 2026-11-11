//! One native enriched-grid assessment of the selected fixed-grid design.
//! The optimizer remains on its original grid and keeps its exact checkpoint.
use super::*;
use fs_solver::LinearOp;
use fs_topopt::sdf3::robust::envelope::affine_load_space::constraints::ReferenceAffineLoadSpaceConstraints3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::refinement::{
    AffineLoadSpaceRefinement3, AffineLoadSpaceRefinementOptions3,
};

#[derive(Debug, Clone, Copy)]
pub(crate) struct Policy {
    pub max_transfer_terms: usize,
    pub compliance_tolerance: f64,
    pub response_tolerance: f64,
    pub stress_tolerance: f64,
    pub inventory_tolerance: f64,
}
impl Policy {
    /// The optimized functional and the separately retained inventory have
    /// different units and change gates in a fixed-volume stiffness design.
    pub(super) fn objective_tolerance(self, goal: AffineLoadSpaceDesignGoal3) -> f64 {
        match goal {
            AffineLoadSpaceDesignGoal3::WorstCompliance { .. } => self.compliance_tolerance,
            AffineLoadSpaceDesignGoal3::MinimumMaterial { .. } => self.inventory_tolerance,
        }
    }

    pub(crate) fn validate(self, spec: &Spec) -> Result<()> {
        let fine_leaves = spec
            .level
            .checked_add(1)
            .and_then(|level| level.checked_mul(3))
            .and_then(|shift| 1usize.checked_shl(shift))
            .unwrap_or(usize::MAX);
        if spec.affine_load_space.is_none()
            || !(1..=2).contains(&spec.level)
            || spec.max_level <= spec.level
            || spec.leaves < fine_leaves
            || !(1..=4_000_000).contains(&self.max_transfer_terms)
            || [
                self.compliance_tolerance,
                self.response_tolerance,
                self.stress_tolerance,
                self.inventory_tolerance,
            ]
            .iter()
            .any(|v| !(0.0..=1e6).contains(v))
        {
            return Err(fail(
                "cli-study-sdf3-input",
                "affine-grid-assessment requires the affine optimizer, room for one uniform enrichment in maximum-level/leaves, 1..=4000000 transfer terms, and finite change tolerances in [0,1e6]",
            ));
        }
        Ok(())
    }

    pub(crate) fn memory_envelope(self, spec: &Spec) -> usize {
        // A simultaneous fine operator, two-grid transfer/field workspace,
        // two additional sampled-strength/response buffers and fine ladder.
        let stress = spec
            .affine_load_space
            .and_then(|c| c.stress)
            .map_or(0, |s| s.options.max_points);
        (64usize * 1024 * 1024)
            .saturating_add(65_536usize.saturating_mul(spec.leaves))
            .saturating_add(128usize.saturating_mul(self.max_transfer_terms))
            .saturating_add(
                spec.material_acceleration
                    .map_or(0, |law| law.memory_envelope()),
            )
            .saturating_add(2048usize.saturating_mul(stress))
            .saturating_add(
                spec.affine_responses
                    .as_ref()
                    .map_or(0, |f| f.memory_envelope().saturating_mul(2)),
            )
            .saturating_add(spec.solver.memory_envelope(spec.level.saturating_add(1)))
    }
}

#[derive(Debug)]
pub(super) struct Outcome {
    pub report: Option<AffineLoadSpaceRefinement3>,
    pub reason: Option<String>,
    pub status: &'static str,
    pub coarse_cells: usize,
    pub fine_cells: usize,
    pub coarse_dofs: usize,
    pub fine_dofs: usize,
    pub coarse_mean_volumes: Vec<Option<f64>>,
    pub fine_mean_volumes: Vec<Option<f64>>,
    /// Assessment-only increments; State.spent retains the cumulative totals.
    pub work: SolveWork,
    pub geometry: QuadratureWork3,
}
impl Outcome {
    pub(super) fn passed(&self, policy: Policy, goal: AffineLoadSpaceDesignGoal3) -> bool {
        self.status == "checked"
            && self.report.as_ref().is_some_and(|r| {
                r.goal == Some(goal)
                    && r.stable_and_feasible_with_constraint_tolerances(
                        policy.compliance_tolerance,
                        policy.response_tolerance,
                        policy.stress_tolerance,
                    )
                    && r.inventory.as_ref().is_some_and(|i| {
                        i.volume_change.is_finite()
                            && i.volume_change >= 0.0
                            && i.volume_change <= policy.inventory_tolerance
                    })
                    && r.objective_change()
                        .is_some_and(|v| v <= policy.objective_tolerance(goal))
            })
    }
}

fn budget(error: &GoalRefinementError3) -> bool {
    geometry_budget(error)
        || solver::setup_budget(error)
        || matches!(
            error,
            GoalRefinementError3::Evaluation(
                EvaluationStop::LinearBudget { .. }
                    | EvaluationStop::TotalBudget { .. }
                    | EvaluationStop::Breakdown {
                        stage: "sdf3-material-body-load-budget"
                            | "sdf3-affine-stress-point-budget"
                    }
            ) | GoalRefinementError3::Background(
                OctreeError3::LeafBudget | OctreeError3::LevelBudget
            ) | GoalRefinementError3::Physics(ElasticityError3::Invalid(
                "transfer coefficient budget exhausted"
                    | "material body-load term allowance exhausted"
            ))
        )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    policy: Policy,
    spec: &Spec,
    tree: &Octree3,
    study: &mut CutDensityStudy3<AdaptiveSolveSpace3>,
    rho: &[f64],
    quadrature: &mut QuadratureControl3<'_>,
    control: &mut SolveControl<'_>,
    gate: &CancelGate,
    mut poll: impl FnMut() -> ControlFlow<()>,
) -> Outcome {
    let before = control.work();
    let geometry = quadrature.work();
    let mut out = Outcome {
        report: None,
        reason: None,
        status: "numerical-failure",
        coarse_cells: study.cells(),
        fine_cells: 0,
        coarse_dofs: study.operator().elasticity().n(),
        fine_dofs: 0,
        coarse_mean_volumes: Vec::new(),
        fine_mean_volumes: Vec::new(),
        work: SolveWork::default(),
        geometry: QuadratureWork3::default(),
    };
    let result = (|| -> std::result::Result<AffineLoadSpaceRefinement3, GoalRefinementError3> {
        control.checkpoint("cli-sdf3-affine-grid-assessment")?;
        let config = spec
            .affine_load_space
            .ok_or(GoalRefinementError3::Invalid("missing affine optimizer"))?;
        let load_policy = AffineLoadSpace3::new(spec.loads.len(), config.cap_j)
            .map_err(|_| GoalRefinementError3::Invalid("invalid affine load policy"))?;
        let coarse_responses =
            response_binding::prepare_reference(spec, study.operator().elasticity(), &mut poll)?;
        out.coarse_mean_volumes = coarse_responses
            .as_ref()
            .map_or_else(Vec::new, response_binding::ReferencePrepared::mean_volumes);
        let fine_tree = tree.refined(study.operator().elasticity().leaves(), &mut poll)?;
        let bounds = HexCell::try_new(spec.bounds.0, spec.bounds.1)
            .map_err(|_| GoalRefinementError3::Invalid("invalid assessment domain bounds"))?;
        let material = loading::material(spec)
            .map_err(|_| GoalRefinementError3::Invalid("invalid assessment material"))?;
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
        let mut build = |grid: &Octree3, checkpoint: &mut dyn FnMut() -> ControlFlow<()>| {
            if checkpoint().is_break() {
                return Err(EvaluationStop::Cancelled.into());
            }
            loading::build_operator(spec, bounds, grid, domain, &material, quadrature)
                .map_err(GoalRefinementError3::from)
        };
        let corrections = grid_refinement::corrections(spec, tree, &mut build, &mut poll)?;
        let operator = build(&fine_tree, &mut poll)?;
        out.fine_cells = operator.cells();
        out.fine_dofs = operator.n();
        let fine_responses = response_binding::prepare_reference(spec, &operator, &mut poll)?;
        out.fine_mean_volumes = fine_responses
            .as_ref()
            .map_or_else(Vec::new, response_binding::ReferencePrepared::mean_volumes);
        let enriched = {
            let mut refs = Vec::new();
            if matches!(spec.solver, solver::Policy::Multilevel { .. }) {
                refs.push(study.operator().elasticity());
                refs.extend(corrections.iter());
            }
            spec.solver.wrap(operator, &refs, &mut poll)?
        };
        loading::with_affine_laws(spec, |dead, laws| {
            let basis: Vec<_> = laws.iter().map(|law| law.load).collect();
            response_binding::with_reference_limits(coarse_responses.as_ref(), |coarse| {
                response_binding::with_reference_limits(fine_responses.as_ref(), |fine| {
                    study.estimate_affine_load_space_enrichment_for_goal_with_grid_constraints(
                        enriched,
                        rho,
                        dead,
                        &basis,
                        config.options.robustness,
                        &load_policy,
                        config.goal,
                        ReferenceAffineLoadSpaceConstraints3 {
                            responses: coarse,
                            stress: config.stress,
                            sharpness: config.joint_sharpness,
                        },
                        ReferenceAffineLoadSpaceConstraints3 {
                            responses: fine,
                            stress: config.stress,
                            sharpness: config.joint_sharpness,
                        },
                        AffineLoadSpaceRefinementOptions3 {
                            max_transfer_terms: policy.max_transfer_terms,
                            ..Default::default()
                        },
                        control,
                    )
                })
            })
        })
    })();
    let after = control.work();
    let after_geometry = quadrature.work();
    out.work = SolveWork {
        linear_iterations: after.linear_iterations - before.linear_iterations,
        linear_solves: after.linear_solves - before.linear_solves,
        preconditioner_operator_applications: after.preconditioner_operator_applications
            - before.preconditioner_operator_applications,
        preconditioner_galerkin_products: after.preconditioner_galerkin_products
            - before.preconditioner_galerkin_products,
    };
    out.geometry = QuadratureWork3 {
        boxes: after_geometry.boxes - geometry.boxes,
        points: after_geometry.points - geometry.points,
        field_evaluations: after_geometry.field_evaluations - geometry.field_evaluations,
    };
    match result {
        Ok(report) => {
            out.report = Some(report);
            out.status = "checked";
        }
        Err(error) => {
            out.status = if gate.is_requested() {
                "cancelled"
            } else if poll().is_break() || budget(&error) {
                "budget-exhausted"
            } else if matches!(
                error,
                GoalRefinementError3::Evaluation(EvaluationStop::Cancelled)
                    | GoalRefinementError3::Physics(ElasticityError3::Cancelled)
                    | GoalRefinementError3::Background(OctreeError3::Cancelled)
            ) {
                "cancelled"
            } else {
                "numerical-failure"
            };
            out.reason = Some(error.to_string());
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
pub(super) fn attach(
    policy: Policy,
    spec: &Spec,
    tree: &Octree3,
    study: &mut CutDensityStudy3<AdaptiveSolveSpace3>,
    state: &mut State,
    quadrature: &mut QuadratureControl3<'_>,
    control: &mut SolveControl<'_>,
    gate: &CancelGate,
    poll: impl FnMut() -> ControlFlow<()>,
) {
    let out = run(
        policy,
        spec,
        tree,
        study,
        &state.selected().rho,
        quadrature,
        control,
        gate,
        poll,
    );
    state.status = if out.status != "checked" {
        out.status
    } else if state.optimization_status == "completed"
        && !spec.affine_load_space.is_some_and(|config| out.passed(policy, config.goal))
    {
        "grid-unresolved"
    } else {
        state.optimization_status
    };
    state.assessment = Some(out);
}
