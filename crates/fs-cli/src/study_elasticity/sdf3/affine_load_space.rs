//! Minimum dilated material for f = B*rho_eroded + F*q, ||q|| <= 1.
//! Native self-weight is a common offset, never another uncertain basis column.
use super::*;
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlReport, ProjectedAlStop, ProjectedAlWork,
};
use fs_topopt::sdf3::response::limits::design::MinimumMaterialObjective3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::{
    AffineLoadSpace3,
    design::{
        AffineLoadSpaceDesign3, AffineLoadSpaceDesignCheckpoint3, AffineLoadSpaceDesignGoal3,
        AffineLoadSpaceDesignIteration3, AffineLoadSpaceDesignOptions3,
        AffineLoadSpaceDesignPoint3, AffineLoadSpaceGradientCheck3,
    },
};
#[path = "affine_load_space/output.rs"]
mod output;
#[path = "affine_load_space/resume.rs"]
mod resume;

const SCOPE: &str = "Estimated fixed-grid linear-elastic minimum-material design for f=B*rho_eroded+Fq, ||q||_2<=1. The common load is native material self-weight at eroded occupancy; no separate dead load is authored. The force columns retain physical amplitudes and correlations, not probabilities. A bounded quadratic-ball solve selects the complete affine compliance witness; its primal/dual gap and residuals are numerical diagnostics, not certified bounds. Eroded material supplies stiffness and self-weight; dilated material supplies the minimized volume. With self-weight, erosion is not asserted to be the worst manufacturing realization. Initial full-chain objective and constraint gradients are independently checked. Hard or unresolved load maximizers cannot supply a smooth gradient. Only physically cap-feasible accepted points enter the least-material archive; KKT values describe the last accepted point. Accepted results are retained before more optimization. Cross-process resume restores accepted optimizer state under the same source and executable, verifies the native physical-problem binding, rebuilds geometry and re-solves the accepted and distinct feasible incumbent endpoints against the original work allowances. Earlier optimizer steps and the initial gradient experiment are not replayed. Failed restoration and uncheckpointed work lost after a crash cannot be durably charged; previous results remain intact. Individual kernels and ledger writes are indivisible; wall and cancellation limits are cooperative. Memory is an admission envelope. No stress, buckling, continuum-error, manufacturing-tolerance, global design-optimality, safety or experimental-validation claim.";

const GOAL: AffineLoadSpaceDesignGoal3 = AffineLoadSpaceDesignGoal3::MinimumMaterial {
    objective: MinimumMaterialObjective3::VolumeFraction,
    volume_scale: 1.0,
};

#[derive(Clone, Copy, Debug)]
pub(super) struct Config {
    pub cap_j: f64,
    pub options: AffineLoadSpaceDesignOptions3,
}
impl Config {
    fn policy(self, dimension: usize) -> Result<AffineLoadSpace3> {
        AffineLoadSpace3::new(dimension, self.cap_j).map_err(|e| fail("cli-study-sdf3-input", e))
    }

    pub(super) fn validate(self, spec: &Spec) -> Result<()> {
        self.policy(spec.loads.len())?;
        self.options
            .optimizer
            .validate::<EvaluationStop>(spec.leaves)
            .map_err(|e| fail("cli-study-sdf3-input", e.to_string()))?;
        if spec.material_acceleration.is_none()
            || spec.max_level != spec.level
            || spec.schedule.len() != 1
            || spec.loads.iter().any(|(_, weight)| *weight != 1.0)
            || !(1..=1024).contains(&spec.updates)
            || !(0.001..=0.9).contains(&self.options.density_floor)
            || !(self.options.density_floor..=1.0).contains(&spec.density)
            || !(1e-12..=1e15).contains(&self.cap_j)
            || !(0.0..0.5).contains(&self.options.robustness.threshold_offset)
            || !(1..=10_000).contains(&self.options.optimizer.max_evaluations)
            || !(1e-10..=1e-2).contains(&self.options.optimizer.tolerance)
        {
            return Err(fail(
                "cli-study-sdf3-input",
                "affine-load-space-limited-simp requires material-acceleration, a fixed grid and model, unit basis weights, bounded updates/evaluations, positive compliance cap, 0<=threshold-offset<0.5 and initial density above its floor",
            ));
        }
        Ok(())
    }

    /// Physical cap checks do not inherit the AL stopping tolerance.
    fn feasible(self, point: &AffineLoadSpaceDesignPoint3) -> bool {
        let e = &point.evaluation;
        e.policy.reference() == self.cap_j
            && e.maximum.is_finite()
            && (0.0..=1.0).contains(&e.maximum)
            && e.worst_compliance.is_finite()
            && (0.0..=self.cap_j).contains(&e.worst_compliance)
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
    let prior = recovery.map_or(Spent::default(), |r| r.spent).add(
        SolveWork::default(),
        QuadratureWork3::default(),
        preparation_wall_s,
    )?;
    prior.validate(spec)?;
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
    let tree = Octree3::uniform(spec.level as u8, spec.max_level as u8, spec.leaves)
        .map_err(|e| fail("cli-study-sdf3-background", e.to_string()))?;
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
    let prepared = (|| -> std::result::Result<AdaptiveSolveSpace3, GoalRefinementError3> {
        let mut build = |grid: &Octree3, checkpoint: &mut dyn FnMut() -> ControlFlow<()>| {
            if checkpoint().is_break() {
                return Err(EvaluationStop::Cancelled.into());
            }
            loading::build_operator(spec, bounds, grid, domain, &material, &mut quadrature)
                .map_err(GoalRefinementError3::from)
        };
        let corrections = solver::correction_spaces(spec, &mut build, &mut || poll())?;
        let operator = build(&tree, &mut || poll())?;
        let refs: Vec<_> = corrections.iter().collect();
        spec.solver.wrap(operator, &refs, || poll())
    })();
    let operator = prepared.map_err(|e| Failure {
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
    let geometry = quadrature.work();
    let study = regions::bind(
        CutDensityStudy3::new(operator, spec.radius, spec.schedule[0]),
        spec,
    )?;
    let mut study = loading::bind_material_acceleration(study, spec, gate, || poll())?;
    let raw = vec![spec.density; study.cells()];
    let mut cp = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: spec.linear - prior.linear.linear_iterations,
            per_solve_iterations: spec.per_solve,
        },
        &mut cp,
    );
    let mut state = loading::with_laws(spec, |laws| -> Result<State> {
        let forces = laws
            .iter()
            .map(|law| {
                study
                    .operator()
                    .elasticity()
                    .reference_load(law.load, &poll)
            })
            .collect::<std::result::Result<Vec<_>, _>>()
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
            })?;
        let basis: Vec<_> = forces.iter().map(Vec::as_slice).collect();
        // The changing material body force is assembled by the native owner
        // once at eroded occupancy. It must not be baked into each basis load.
        let dead = vec![0.0; forces[0].len()];
        let mut design = match recovery {
            None => AffineLoadSpaceDesign3::new(
                &mut study,
                &dead,
                &basis,
                &raw,
                policy,
                GOAL,
                config.options,
                &mut control,
            ),
            Some(recovery) => AffineLoadSpaceDesign3::restore(
                &mut study,
                &dead,
                &basis,
                &mut recovery.checkpoint.clone(),
                policy,
                GOAL,
                config.options,
                &mut control,
            ),
        }
        .map_err(|e| numerical_failure(e, gate, poll().is_break()))?;
        let result = design.run(0);
        let mut state = State::capture(&design, result);
        let updates = limit.min(spec.updates - design.optimizer_work().iterations);
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
        }
        Ok(state)
    })?;
    state.spent = prior.add(control.work(), geometry, start.elapsed().as_secs_f64())?;
    state.status = state.terminal(config, gate, state.spent.wall_s >= spec.wall_s);
    Ok(Computation { study, state })
}

pub(super) fn drive(
    spec: &Spec,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
) -> Result<Outcome> {
    resume::drive(spec, ledger, cap, gate, None)
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
    resume::drive(&spec, ledger, cap, gate, Some(old))
}

#[cfg(test)]
#[path = "affine_load_space/tests.rs"]
mod tests;
