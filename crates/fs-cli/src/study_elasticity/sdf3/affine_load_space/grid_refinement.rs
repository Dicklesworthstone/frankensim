//! One bounded redesign on a genuinely enriched native grid.
//! Candidate AL state belongs to that grid; the installed source stays public
//! until a fresh complete physical recheck authorizes promotion.
use super::*;
use fs_topopt::sdf3::robust::envelope::affine_load_space::constraints::ReferenceAffineLoadSpaceConstraints3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::refit::{
    AffineLoadSpaceRefitError3, AffineLoadSpaceRefitOptions3,
};

#[path = "grid_refinement/drive.rs"]
mod driver;
pub(super) use driver::drive;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Policy {
    pub max_updates: usize,
    pub max_evaluations: usize,
}
impl Policy {
    pub(crate) fn validate(self, spec: &Spec) -> Result<()> {
        let required = spec
            .level
            .checked_add(2)
            .and_then(|l| l.checked_mul(3))
            .and_then(|shift| 1usize.checked_shl(shift))
            .unwrap_or(usize::MAX);
        if spec.affine_load_space.is_none()
            || spec.affine_grid_assessment.is_none()
            || spec.max_level < spec.level.saturating_add(2)
            || spec.leaves < required
            || !(1..=1024).contains(&self.max_updates)
            || !(2..=10_000).contains(&self.max_evaluations)
        {
            return Err(fail(
                "cli-study-sdf3-input",
                "affine-grid-refinement requires affine-grid-assessment, room for one redesign and its final probe, 1..=1024 candidate updates and 2..=10000 candidate optimizer evaluations",
            ));
        }
        Ok(())
    }

    pub(crate) fn memory_envelope(self, spec: &Spec) -> usize {
        let stress = spec
            .affine_load_space
            .and_then(|c| c.stress)
            .map_or(0, |s| s.options.max_points);
        (64usize * 1024 * 1024)
            .saturating_add(65_536usize.saturating_mul(spec.leaves))
            .saturating_add(
                64usize
                    .saturating_mul(spec.leaves)
                    .saturating_mul(spec.leaves),
            )
            .saturating_add(1024usize.saturating_mul(stress))
            .saturating_add(
                spec.affine_responses
                    .as_ref()
                    .map_or(0, |f| f.memory_envelope()),
            )
            .saturating_add(
                spec.material_acceleration
                    .map_or(0, |l| l.memory_envelope()),
            )
            .saturating_add(spec.solver.memory_envelope(spec.level.saturating_add(2)))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Phase {
    Preparing,
    Candidate,
    Promoted,
}
impl Phase {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Preparing => "preparing",
            Self::Candidate => "candidate",
            Self::Promoted => "promoted",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Context {
    pub source_receipt: ContentHash,
    /// Sorted indices into the authored initial uniform background tree.
    /// Replaying these topological marks needs no source physics or optimization.
    pub source_marks: Vec<usize>,
    pub phase: Phase,
    pub attempt: usize,
}
impl Context {
    pub(crate) fn new(source: &Computation, source_receipt: ContentHash) -> Result<Self> {
        let background: Vec<_> = source.tree.leaves().iter().copied().collect();
        let mut source_marks = source
            .study
            .operator()
            .elasticity()
            .leaves()
            .iter()
            .map(|leaf| {
                background.binary_search(leaf).map_err(|_| {
                    fail(
                        "cli-study-sdf3-refinement",
                        "source active leaf is absent from its background",
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?;
        source_marks.sort_unstable();
        if source_marks.is_empty() || source_marks.windows(2).any(|w| w[0] >= w[1]) {
            return Err(fail(
                "cli-study-sdf3-refinement",
                "empty or repeated source refinement marks",
            ));
        }
        Ok(Self {
            source_receipt,
            source_marks,
            phase: Phase::Preparing,
            attempt: 1,
        })
    }
    pub(crate) fn working_generation(&self) -> usize {
        usize::from(self.phase != Phase::Preparing)
    }
    pub(crate) fn uses_source_selection(&self) -> bool {
        self.phase != Phase::Promoted
    }
    pub(crate) fn json(&self) -> String {
        format!(
            "{{\"source_receipt\":{:?},\"source_marks\":{:?},\"phase\":{:?},\"attempt\":{}}}",
            self.source_receipt.to_hex(),
            self.source_marks,
            self.phase.as_str(),
            self.attempt
        )
    }
    pub(crate) fn read(spec: &Spec, receipt: &JsonValue) -> Result<Option<Self>> {
        let Some(value) = receipt.get("grid_refinement") else {
            return Ok(None);
        };
        if spec.affine_grid_refinement.is_none() {
            return Err(fail(
                "cli-study-sdf3-refinement",
                "unrequested retained grid refinement",
            ));
        }
        let bad = || {
            fail(
                "cli-study-sdf3-refinement",
                "invalid retained refinement context",
            )
        };
        let source_receipt = value
            .str_field("source_receipt")
            .and_then(ContentHash::from_hex)
            .ok_or_else(bad)?;
        let phase = match value.str_field("phase") {
            Some("preparing") => Phase::Preparing,
            Some("candidate") => Phase::Candidate,
            Some("promoted") => Phase::Promoted,
            _ => return Err(bad()),
        };
        let count = 8usize.pow(spec.level);
        let values = value
            .get("source_marks")
            .and_then(JsonValue::as_array)
            .ok_or_else(bad)?;
        if values.is_empty() || values.len() > count || integer(value, "attempt")? != 1 {
            return Err(bad());
        }
        let source_marks = values
            .iter()
            .map(|v| {
                v.as_f64()
                    .filter(|n| n.is_finite() && *n >= 0.0 && *n < count as f64 && n.fract() == 0.0)
                    .map(|n| n as usize)
                    .ok_or_else(bad)
            })
            .collect::<Result<Vec<_>>>()?;
        if source_marks.windows(2).any(|w| w[0] >= w[1]) {
            return Err(bad());
        }
        Ok(Some(Self {
            source_receipt,
            source_marks,
            phase,
            attempt: 1,
        }))
    }
}

pub(super) fn effective_spec(original: &Spec, context: Option<&Context>) -> Result<Spec> {
    let mut spec = original.clone();
    if let Some(context) = context {
        let policy = spec.affine_grid_refinement.ok_or_else(|| {
            fail(
                "cli-study-sdf3-refinement",
                "refinement context has no authored policy",
            )
        })?;
        if context.working_generation() != 0 {
            spec.updates = policy.max_updates;
            spec.affine_load_space
                .as_mut()
                .ok_or_else(|| fail("cli-study-sdf3-refinement", "missing affine optimizer"))?
                .options
                .optimizer
                .max_evaluations = policy.max_evaluations;
        }
    }
    Ok(spec)
}

pub(super) fn tree(
    spec: &Spec,
    context: Option<&Context>,
    mut poll: impl FnMut() -> ControlFlow<()>,
) -> std::result::Result<Octree3, GoalRefinementError3> {
    let base = Octree3::uniform(spec.level as u8, spec.max_level as u8, spec.leaves)?;
    let Some(context) = context.filter(|c| c.working_generation() != 0) else {
        return Ok(base);
    };
    let leaves: Vec<_> = base.leaves().iter().copied().collect();
    let marks = context
        .source_marks
        .iter()
        .map(|&i| {
            leaves.get(i).copied().ok_or(GoalRefinementError3::Invalid(
                "refinement mark exceeds initial background",
            ))
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    if poll().is_break() {
        return Err(EvaluationStop::Cancelled.into());
    }
    base.refined(&marks, poll).map_err(Into::into)
}

/// Build the complete nearest-first uniform correction ladder below the
/// actual working grid. Initial-grid specs still own all physical patches.
pub(super) fn corrections(
    spec: &Spec,
    tree: &Octree3,
    build: &mut impl FnMut(
        &Octree3,
        &mut dyn FnMut() -> ControlFlow<()>,
    ) -> std::result::Result<AdaptiveElasticity3, GoalRefinementError3>,
    poll: &mut dyn FnMut() -> ControlFlow<()>,
) -> std::result::Result<Vec<AdaptiveElasticity3>, GoalRefinementError3> {
    let solver::Policy::Multilevel { coarsest_level, .. } = spec.solver else {
        return Ok(Vec::new());
    };
    let level = tree
        .leaves()
        .iter()
        .map(|c| u32::from(c.level()))
        .max()
        .ok_or(GoalRefinementError3::Invalid("empty working background"))?;
    let mut result = Vec::new();
    for level in (coarsest_level..level).rev() {
        if poll().is_break() {
            return Err(EvaluationStop::Cancelled.into());
        }
        let grid = Octree3::uniform(level as u8, spec.max_level as u8, spec.leaves)?;
        result.push(build(&grid, poll)?);
    }
    Ok(result)
}

pub(super) struct Prepared {
    pub study: CutDensityStudy3<AdaptiveSolveSpace3>,
    pub inherited_rho: Vec<f64>,
    pub tree: Octree3,
    pub context: Context,
    pub spent: Spent,
}
pub(super) struct Preparation {
    pub result: Result<Prepared>,
    pub spent: Spent,
}

pub(super) fn numerical(error: GoalRefinementError3, gate: &CancelGate, expired: bool) -> Failure {
    let budget = geometry_budget(&error)
        || solver::setup_budget(&error)
        || matches!(
            &error,
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
            ))
        );
    let cancelled = matches!(
        &error,
        GoalRefinementError3::Evaluation(EvaluationStop::Cancelled)
    );
    Failure {
        code: "cli-study-sdf3-refinement",
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

pub(super) fn prepare(
    spec: &Spec,
    source: &Computation,
    mut context: Context,
    gate: &CancelGate,
) -> Preparation {
    let prior = source.state.spent;
    let started = Instant::now();
    let poll = || {
        if gate.is_requested() || prior.wall_s + started.elapsed().as_secs_f64() >= spec.wall_s {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    let early = |error| Preparation {
        result: Err(error),
        spent: prior,
    };
    if let Err(e) = prior.require_remaining(spec) {
        return early(e);
    }
    let effective = match effective_spec(
        spec,
        Some(&Context {
            phase: Phase::Candidate,
            ..context.clone()
        }),
    ) {
        Ok(s) => s,
        Err(e) => return early(e),
    };
    context.phase = Phase::Candidate;
    let mut qp = |_| poll();
    let mut quadrature = match QuadratureControl3::new(
        QuadratureOptions3 {
            max_boxes: spec.boxes - prior.geometry.boxes,
            max_points: spec.points - prior.geometry.points,
            ..Default::default()
        },
        &mut qp,
    ) {
        Ok(q) => q,
        Err(e) => return early(fail("cli-study-sdf3-refinement", e.to_string())),
    };
    let mut cp = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: spec.linear - prior.linear.linear_iterations,
            per_solve_iterations: spec.per_solve,
        },
        &mut cp,
    );
    let result = (|| {
        let tree = tree(spec, Some(&context), &poll)?;
        let bounds = HexCell::try_new(spec.bounds.0, spec.bounds.1)
            .map_err(|_| GoalRefinementError3::Invalid("invalid refinement domain"))?;
        let material = loading::material(spec)
            .map_err(|_| GoalRefinementError3::Invalid("invalid refinement material"))?;
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
        let mut build = |tree: &Octree3, checkpoint: &mut dyn FnMut() -> ControlFlow<()>| {
            if checkpoint().is_break() {
                return Err(EvaluationStop::Cancelled.into());
            }
            loading::build_operator(spec, bounds, tree, domain, &material, &mut quadrature)
                .map_err(Into::into)
        };
        let coarse = corrections(spec, &tree, &mut build, &mut || poll())?;
        let operator = build(&tree, &mut || poll())?;
        let observations = response_binding::prepare_reference(spec, &operator, &poll)
            .map_err(GoalRefinementError3::from)?;
        let references: Vec<_> = coarse.iter().collect();
        let operator = spec.solver.wrap(operator, &references, &poll)?;
        let config = effective
            .affine_load_space
            .ok_or(GoalRefinementError3::Invalid("missing affine design"))?;
        let policy = AffineLoadSpace3::new(spec.loads.len(), config.cap_j)
            .map_err(|_| GoalRefinementError3::Invalid("invalid affine load family"))?;
        let assessment = spec
            .affine_grid_assessment
            .ok_or(GoalRefinementError3::Invalid("missing assessment policy"))?;
        let seed = response_binding::with_reference_limits(observations.as_ref(), |responses| {
            source
                .study
                .prepare_affine_load_space_refit_with_constraints(
                    operator,
                    &source.state.selected().rho,
                    spec.loads.len(),
                    &policy,
                    config.goal,
                    AffineLoadSpaceRefitOptions3 {
                        design: config.options,
                        accepted_steps: effective.updates,
                        max_transfer_terms: assessment.max_transfer_terms,
                    },
                    ReferenceAffineLoadSpaceConstraints3 {
                        responses,
                        stress: config.stress,
                        sharpness: config.joint_sharpness,
                    },
                    &mut control,
                )
        })?;
        Ok::<_, AffineLoadSpaceRefitError3>((tree, seed))
    })();
    let spent = match prior.add(
        control.work(),
        quadrature.work(),
        started.elapsed().as_secs_f64(),
    ) {
        Ok(s) => s,
        Err(e) => return early(e),
    };
    let result = result
        .map(|(tree, seed)| Prepared {
            study: seed.study,
            inherited_rho: seed.inherited_rho,
            tree,
            context,
            spent,
        })
        .map_err(|e| match e {
            AffineLoadSpaceRefitError3::Goal(e) => numerical(e, gate, poll().is_break()),
            AffineLoadSpaceRefitError3::Design(e) => numerical_failure(e, gate, poll().is_break()),
        });
    Preparation { result, spent }
}

pub(super) fn assess(spec: &Spec, gate: &CancelGate, run: &mut Computation) -> Result<()> {
    let policy = spec.affine_grid_assessment.ok_or_else(|| {
        fail(
            "cli-study-sdf3-refinement",
            "missing final assessment policy",
        )
    })?;
    let prior = run.state.spent;
    prior.require_remaining(spec)?;
    let started = Instant::now();
    let poll = || {
        if gate.is_requested() || prior.wall_s + started.elapsed().as_secs_f64() >= spec.wall_s {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    let mut qp = |_| poll();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            max_boxes: spec.boxes - prior.geometry.boxes,
            max_points: spec.points - prior.geometry.points,
            ..Default::default()
        },
        &mut qp,
    )
    .map_err(|e| fail("cli-study-sdf3-refinement", e.to_string()))?;
    let mut cp = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: spec.linear - prior.linear.linear_iterations,
            per_solve_iterations: spec.per_solve,
        },
        &mut cp,
    );
    grid_assessment::attach(
        policy,
        spec,
        &run.tree,
        &mut run.study,
        &mut run.state,
        &mut quadrature,
        &mut control,
        gate,
        poll,
    );
    run.state.spent = prior.add(
        control.work(),
        quadrature.work(),
        started.elapsed().as_secs_f64(),
    )?;
    Ok(())
}
