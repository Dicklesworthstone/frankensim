//! Native minimum-material design for the continuous load set f = Fq, |q| <= 1.
//! Geometry, supports, projection thresholds and the load basis remain fixed.
use super::*;
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlReport, ProjectedAlStop, ProjectedAlWork,
};
use fs_topopt::sdf3::robust::envelope::load_space::{
    LoadSpaceGradientCheck3,
    design::{
        LoadSpaceDesignIteration3, LoadSpaceDesignPoint3,
        limits::{LoadSpaceCap3, LoadSpaceLimitedDesign3, LoadSpaceLimitedOptions3},
    },
};

#[path = "load_space/output.rs"]
mod output;

const SCOPE: &str = "Estimated fixed-grid linear-elastic minimum-material design over every simultaneous force f=Fq with ||q||_2<=1. Authored basis forces include their amplitudes and correlations; they are not probabilities or mutually exclusive scenarios. The worst compliance is computed from the complete cross-compliance spectrum on eroded material. A smooth spectral upper surrogate drives the existing projected augmented-Lagrangian optimizer; physical feasibility additionally checks the unrelaxed compliance cap. The objective is dilated material fraction, with nominal and eroded fields retained separately. Projection-threshold robustness is a discrete density model, not a geometric manufacturing tolerance. Initial spectral and material gradients are checked by independently solved directional differences. Only accepted feasible points compete for the least-material incumbent; KKT values refer to the last accepted iterate. Accepted outputs are durable before more optimization. Partial runs remain exportable but this mode does not restore optimizer state across processes. No stress, buckling, continuum error, mesh-independent optimum, safety or experimental-validation claim. Memory is an admission envelope and wall/cancellation checks are cooperative.";

#[derive(Clone, Copy, Debug)]
pub(super) struct Config {
    pub cap_j: f64,
    pub temperature: f64,
    pub options: LoadSpaceLimitedOptions3,
}
impl Config {
    fn cap(self, dimension: usize) -> Result<LoadSpaceCap3> {
        LoadSpaceCap3::new(dimension, self.cap_j, self.temperature)
            .map_err(|e| fail("cli-study-sdf3-input", e))
    }
    pub(super) fn validate(self, spec: &Spec) -> Result<()> {
        self.cap(spec.loads.len())?;
        self.options
            .optimizer
            .validate::<EvaluationStop>(spec.leaves)
            .map_err(|e| fail("cli-study-sdf3-input", e.to_string()))?;
        if spec.max_level != spec.level
            || !(1..=1024).contains(&spec.updates)
            || !(0.001..=0.9).contains(&self.options.density_floor)
            || !(self.options.density_floor..=1.0).contains(&spec.density)
            || !(1e-12..=1e15).contains(&self.cap_j)
            || !(1e-6..=1.0).contains(&self.temperature)
            || !(0.0..0.5).contains(&self.options.robustness.threshold_offset)
            || !(1..=10_000).contains(&self.options.optimizer.max_evaluations)
            || !(1e-10..=1e-2).contains(&self.options.optimizer.tolerance)
        {
            return Err(fail(
                "cli-study-sdf3-input",
                "load-space-limited-simp requires a fixed grid, bounded updates/evaluations, a positive compliance cap and spectral temperature, 0<=threshold-offset<0.5, and an initial density above its declared floor",
            ));
        }
        Ok(())
    }
}

type NumericalResult = std::result::Result<ProjectedAlReport, ProjectedAlError<EvaluationStop>>;

struct State {
    accepted: LoadSpaceDesignPoint3,
    best: Option<LoadSpaceDesignPoint3>,
    history: Vec<LoadSpaceDesignIteration3>,
    gradient: LoadSpaceGradientCheck3,
    report: Option<ProjectedAlReport>,
    error: Option<ProjectedAlError<EvaluationStop>>,
    optimizer: ProjectedAlWork,
    spent: Spent,
    status: &'static str,
}
impl State {
    fn capture(
        design: &LoadSpaceLimitedDesign3<'_, '_, AdaptiveSolveSpace3>,
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
        }
    }
    fn selected(&self) -> &LoadSpaceDesignPoint3 {
        self.best.as_ref().unwrap_or(&self.accepted)
    }
    fn terminal(&self, gate: &CancelGate, expired: bool) -> &'static str {
        if gate.is_requested() {
            return "cancelled";
        }
        if expired {
            return "budget-exhausted";
        }
        if let Some(error) = &self.error {
            return match error {
                ProjectedAlError::Evaluation(
                    EvaluationStop::LinearBudget { .. } | EvaluationStop::TotalBudget { .. },
                ) => "budget-exhausted",
                _ => "numerical-failure",
            };
        }
        match self.report.as_ref().map(|r| r.stop) {
            Some(ProjectedAlStop::Converged) if self.history.last().is_some_and(|r| r.feasible) => {
                "completed"
            }
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
    Failure {
        code: "cli-study-sdf3-load-space",
        message: error.to_string(),
        exit: if gate.is_requested() {
            exit::CANCELLED
        } else if expired || budget {
            exit::BUDGET
        } else {
            exit::REFUSED
        },
    }
}

fn compute_observed(
    spec: &Spec,
    gate: &CancelGate,
    limit: usize,
    mut observed: impl FnMut(&CutDensityStudy3<AdaptiveSolveSpace3>, &State) -> Result<()>,
) -> Result<Computation> {
    let config = spec
        .load_space
        .ok_or_else(|| fail("cli-study-sdf3-input", "missing load-space optimizer"))?;
    let cap = config.cap(spec.loads.len())?;
    let start = Instant::now();
    let poll = || {
        if gate.is_requested() || start.elapsed().as_secs_f64() >= spec.wall_s {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    if poll().is_break() {
        return Err(Failure {
            code: "cli-study-sdf3-stopped",
            message: "load-space study stopped before geometry".into(),
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
    let material = IsotropicElastic::new(spec.youngs, spec.poisson, 1.0)
        .map_err(|e| fail("cli-study-sdf3-material", e.to_string()))?;
    let mut cp = |_| poll();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            max_boxes: spec.boxes,
            max_points: spec.points,
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
    let mut study = regions::bind(
        CutDensityStudy3::new(operator, spec.radius, spec.schedule[0]),
        spec,
    )?;
    let raw = vec![spec.density; study.cells()];
    let mut cp = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: spec.linear,
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
                code: "cli-study-sdf3-load-space",
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
        let mut design = LoadSpaceLimitedDesign3::new(
            &mut study,
            &basis,
            &raw,
            cap,
            config.options,
            &mut control,
        )
        .map_err(|e| numerical_failure(e, gate, poll().is_break()))?;
        let result = design.run(0);
        let mut state = State::capture(&design, result);
        let updates = limit.min(spec.updates);
        for _ in 0..updates {
            if state.error.is_some()
                || state
                    .report
                    .as_ref()
                    .is_some_and(|r| r.stop != ProjectedAlStop::IterationLimit)
            {
                break;
            }
            state.spent =
                Spent::default().add(design.work(), geometry, start.elapsed().as_secs_f64())?;
            observed(design.study(), &state)?;
            let result = design.run(1);
            state = State::capture(&design, result);
        }
        Ok(state)
    })?;
    state.spent = Spent::default().add(control.work(), geometry, start.elapsed().as_secs_f64())?;
    state.status = state.terminal(gate, state.spent.wall_s >= spec.wall_s);
    Ok(Computation { study, state })
}

pub(super) fn drive(
    spec: &Spec,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
) -> Result<Outcome> {
    let limit = cap.unwrap_or(spec.updates).min(spec.updates);
    let mut previous = None;
    let run = compute_observed(spec, gate, limit, |study, state| {
        let out = output::persist(spec, ledger, study, state, previous, limit)?;
        previous = out
            .pointer
            .strip_prefix("study-")
            .and_then(ContentHash::from_hex);
        // A progress pointer names already durable fields, even after a crash.
        use std::io::Write as _;
        let _ = writeln!(
            std::io::stderr().lock(),
            "{{\"schema\":\"frankensim.cli.sdf3-load-space-progress.v1\",\"run_id\":{},\"iterations_completed\":{}}}",
            quoted(&out.pointer),
            state.optimizer.iterations,
        );
        Ok(())
    })?;
    output::persist(spec, ledger, &run.study, &run.state, previous, limit)
}

#[cfg(test)]
#[path = "load_space/tests.rs"]
mod tests;
