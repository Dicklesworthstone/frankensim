//! Canonical minimum-volume studies using the existing exact-discrete stress
//! adjoint and projected augmented-Lagrangian optimizer. The domain is fixed.
use super::*;
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlReport, ProjectedAlStop, ProjectedAlWork,
};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::design::{StressDesignCheckpoint3, StressDesignIteration3, StressDesignOptions3, StressDesignStudy3};
use fs_topopt::sdf3::stress::{StressError3, StressEvaluation3, StressMeasure3};

#[path = "stress/output.rs"]
mod output;
#[path = "stress/resume.rs"]
mod resume;

const SCOPE: &str = "Estimated fixed-background 3-D linear-elastic SIMP minimum-volume design. The stress_measure declaration identifies the constrained functional in Pa: the original normalized average or an explicitly selected unweighted sampled-peak bound of both relaxed and physical stresses. Only the latter bounds retained numerical samples, including every declared case; neither bounds continuum maximum stress. Independent body, reference-pressure and traction cases remain separate. The initial stress and volume adjoints pass bounded directional finite differences before optimization. Authored solid and zero-density regions are enforced in the physical map before every solve and pullback, including endpoint restoration; the raw density floor continues to guard optimizable material against ersatz stress foldback. Accepted augmented-Lagrangian steps can be infeasible; the least-volume accepted feasible design is retained separately and selected for export when available. KKT residuals describe the last accepted iterate, not an earlier exported incumbent. Only numerical convergence with a feasible current iterate reports completed. Every accepted update is durable before further optimization. Cross-process resume restores the complete accepted optimizer state under the same source and executable, rebuilding geometry and re-solving the accepted and distinct feasible incumbent endpoints against the original allowances. Earlier optimizer steps are not replayed. Failed restoration and work after a process crash cannot be durably charged; previous checkpoints remain intact. Report and package export all retained results without another solve. No continuum safety, global optimum, manufacturing, adaptivity or physical-validation claim is made. Memory is an admission envelope, not measured RSS. Wall time and cancellation are cooperative; quadrature, solver, stress-cell and optimizer boundaries poll, while individual kernels and ledger I/O are indivisible.";
const EVALUATION_CAP: &str = "stress evaluation allowance exhausted";
type StressResult<T> = std::result::Result<T, ProjectedAlError<StressError3>>;

#[derive(Clone, Debug)]
struct Probe {
    direction: &'static str,
    active: usize,
    stress_analytic: f64,
    stress_difference: f64,
    stress_relative_error: f64,
    volume_analytic: f64,
    volume_difference: f64,
    volume_relative_error: f64,
}

#[derive(Clone, Default)]
struct Audit {
    probes: Vec<Probe>,
    evaluations: usize,
    passed: bool,
}

struct State {
    accepted: Option<StressEvaluation3>,
    best: Option<StressEvaluation3>,
    history: Vec<StressDesignIteration3>,
    audit: Audit,
    report: Option<ProjectedAlReport>,
    error: Option<ProjectedAlError<StressError3>>,
    optimizer_work: ProjectedAlWork,
    spent: Spent,
    status: &'static str,
    checkpoint: Option<StressDesignCheckpoint3>,
}
impl State {
    fn new() -> Self {
        Self {
            accepted: None,
            best: None,
            history: Vec::new(),
            audit: Audit::default(),
            report: None,
            error: None,
            optimizer_work: ProjectedAlWork::default(),
            spent: Spent::default(),
            status: "checkpointed",
            checkpoint: None,
        }
    }
    fn iterations(&self) -> usize {
        self.optimizer_work.iterations
    }
    fn selected(&self) -> Option<&StressEvaluation3> {
        self.best.as_ref().or(self.accepted.as_ref())
    }
    fn capture(
        &mut self,
        design: &StressDesignStudy3<'_, '_, AdaptiveSolveSpace3>,
        result: StressResult<ProjectedAlReport>,
    ) {
        self.accepted = Some(design.accepted().clone());
        self.best = design.best_feasible().cloned();
        self.history = design.history().to_vec();
        self.optimizer_work = design.optimizer_work();
        self.checkpoint = Some(design.checkpoint());
        self.report = result.as_ref().ok().cloned();
        self.error = result.err();
    }
}

struct Computation {
    study: CutDensityStudy3<AdaptiveSolveSpace3>,
    state: State,
}

fn exhausted() -> ProjectedAlError<StressError3> {
    ProjectedAlError::Evaluation(StressError3::Invalid(EVALUATION_CAP))
}

/// Two complementary heterogeneous second-order directions exercise both
/// pullbacks. Probes respect the declared density box and charge the same
/// evaluation/Krylov allowances as optimization; no partial audit passes.
fn gradient_gate(
    study: &mut CutDensityStudy3<AdaptiveSolveSpace3>,
    loads: &[LoadCase<'_>],
    rho: &[f64],
    options: StressDesignOptions3,
    measure: StressMeasure3,
    audit: &mut Audit,
    control: &mut SolveControl<'_>,
) -> StressResult<()> {
    let h = 1e-4;
    let mut evaluate = |point: &[f64]| {
        if audit.evaluations == options.optimizer.max_evaluations {
            return Err(exhausted());
        }
        audit.evaluations += 1;
        study
            .evaluate_stress_with_measure(point, loads, options.stress, measure, control)
            .map_err(ProjectedAlError::Evaluation)
    };
    let baseline = evaluate(rho)?;
    for increasing in [true, false] {
        let direction: Vec<f64> = rho
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                let d = if increasing {
                    if i % 2 == 0 { 1.0 } else { 0.5 }
                } else if i % 2 == 0 {
                    -0.5
                } else {
                    -1.0
                };
                if (options.density_floor..=1.0).contains(&(r + 2.0 * h * d)) {
                    d
                } else {
                    0.0
                }
            })
            .collect();
        let active = direction.iter().filter(|&&d| d != 0.0).count();
        if active == 0 {
            continue;
        }
        let one: Vec<_> = rho.iter().zip(&direction).map(|(r, d)| r + h * d).collect();
        let two: Vec<_> = rho
            .iter()
            .zip(&direction)
            .map(|(r, d)| r + 2.0 * h * d)
            .collect();
        let one = evaluate(&one)?;
        let two = evaluate(&two)?;
        let derivative = |g: &[f64]| g.iter().zip(&direction).map(|(g, d)| g * d).sum::<f64>();
        let difference = |base: f64, a: f64, b: f64| (2.0 * (a - base) - 0.5 * (b - base)) / h;
        let relative = |a: f64, b: f64| {
            if !a.is_finite() || !b.is_finite() {
                f64::INFINITY
            } else {
                let scale = a.abs().max(b.abs());
                if scale == 0.0 {
                    0.0
                } else {
                    (a / scale - b / scale).abs()
                }
            }
        };
        let sa = derivative(&baseline.gradient);
        let sd = difference(baseline.aggregate, one.aggregate, two.aggregate);
        let va = derivative(&baseline.volume_gradient);
        let vd = difference(
            baseline.volume_fraction,
            one.volume_fraction,
            two.volume_fraction,
        );
        let probe = Probe {
            direction: if increasing {
                "increasing"
            } else {
                "decreasing"
            },
            active,
            stress_analytic: sa,
            stress_difference: sd,
            stress_relative_error: relative(sa, sd),
            volume_analytic: va,
            volume_difference: vd,
            volume_relative_error: relative(va, vd),
        };
        if !probe.stress_relative_error.is_finite() || !probe.volume_relative_error.is_finite() {
            return Err(ProjectedAlError::Invalid(
                "nonfinite stress gradient comparison",
            ));
        }
        audit.probes.push(probe);
    }
    audit.passed = !audit.probes.is_empty()
        && audit
            .probes
            .iter()
            .all(|p| p.stress_relative_error <= 5e-4 && p.volume_relative_error <= 5e-4);
    if audit.passed {
        Ok(())
    } else {
        Err(ProjectedAlError::Invalid(
            "stress/volume directional gradient gate failed",
        ))
    }
}

fn status(state: &State, gate: &CancelGate, expired: bool) -> &'static str {
    if gate.is_requested() {
        return "cancelled";
    }
    if expired {
        return "budget-exhausted";
    }
    if let Some(error) = &state.error {
        return match error {
            ProjectedAlError::Evaluation(
                StressError3::PointBudget
                | StressError3::Invalid(EVALUATION_CAP)
                | StressError3::Evaluation(
                    EvaluationStop::LinearBudget { .. } | EvaluationStop::TotalBudget { .. },
                ),
            ) => "budget-exhausted",
            _ => "numerical-failure",
        };
    }
    match state.report.as_ref().map(|r| r.stop) {
        Some(ProjectedAlStop::Converged) if state.history.last().is_some_and(|r| r.feasible) => {
            "completed"
        }
        Some(ProjectedAlStop::Stalled | ProjectedAlStop::PenaltyLimit) => "no-feasible-descent",
        _ => "budget-exhausted",
    }
}

fn restoration_failure(error: ProjectedAlError<StressError3>, gate: &CancelGate, expired: bool) -> Failure {
    let budget = matches!(&error, ProjectedAlError::Evaluation(
        StressError3::PointBudget | StressError3::Invalid(EVALUATION_CAP)
            | StressError3::Evaluation(EvaluationStop::LinearBudget { .. } | EvaluationStop::TotalBudget { .. })));
    Failure {
        code: "cli-study-sdf3-stress-restore", message: error.to_string(),
        exit: if gate.is_requested() { exit::CANCELLED }
            else if expired || budget { exit::BUDGET } else { exit::REFUSED },
    }
}

fn compute_observed(
    spec: &Spec,
    gate: &CancelGate,
    limit: usize,
    recovery: Option<&resume::Recovery>,
    preparation_wall_s: f64,
    mut accepted: impl FnMut(&CutDensityStudy3<AdaptiveSolveSpace3>, &State) -> Result<()>,
) -> Result<Computation> {
    let options = spec
        .stress
        .ok_or_else(|| fail("cli-study-sdf3-input", "missing stress optimizer"))?;
    let prior = recovery.map_or(Spent::default(), |r| r.spent)
        .add(SolveWork::default(), QuadratureWork3::default(), preparation_wall_s)?;
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
            message: "stress study stopped before geometry".into(),
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
            max_boxes: spec.boxes - prior.geometry.boxes,
            max_points: spec.points - prior.geometry.points,
            ..Default::default()
        },
        &mut cp,
    )
    .map_err(|e| fail("cli-study-sdf3-geometry", e.to_string()))?;
    let prepared = (|| -> std::result::Result<AdaptiveSolveSpace3, GoalRefinementError3> {
        let mut build = |grid: &Octree3, checkpoint: &mut dyn FnMut() -> ControlFlow<()>| {
            if checkpoint().is_break() { return Err(EvaluationStop::Cancelled.into()); }
            loading::build_operator(spec, bounds, grid, domain, &material, &mut quadrature)
                .map_err(GoalRefinementError3::from)
        };
        // Geometry-only correction spaces use the same original quadrature
        // allowance and support law, including during direct-state recovery.
        let corrections = solver::correction_spaces(spec, &mut build, &mut || poll())?;
        let operator = build(&tree, &mut || poll())?;
        let correction_refs: Vec<_> = corrections.iter().collect();
        spec.solver.wrap(operator, &correction_refs, || poll())
    })();
    let operator = prepared.map_err(|e| Failure {
        code: "cli-study-sdf3-geometry",
        message: e.to_string(),
        exit: if gate.is_requested() { exit::CANCELLED }
            else if poll().is_break() || geometry_budget(&e) || solver::setup_budget(&e) {
                exit::BUDGET
            } else { exit::REFUSED },
    })?;
    let geometry = quadrature.work();
    // Bind the physical map before gradient admission, every candidate solve,
    // and accepted/incumbent endpoint re-evaluation on optimizer restoration.
    let mut study = regions::bind(
        CutDensityStudy3::new(operator, spec.radius, spec.schedule[0]),
        spec,
    )?;
    let raw = vec![spec.density; study.cells()];
    let mut cp = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: spec.linear - prior.linear.linear_iterations,
            per_solve_iterations: spec.per_solve,
        },
        &mut cp,
    );
    let mut state = State::new();
    if let Some(recovery) = recovery {
        state.audit = recovery.audit.clone();
    }
    loading::with_laws(spec, |laws| -> Result<()> {
        let forces = laws
            .iter()
            .map(|law| {
                study
                    .operator()
                    .elasticity()
                    .reference_load(law.load, &poll)
                    .map_err(|e| ProjectedAlError::Evaluation(StressError3::from(e)))
            })
            .collect::<StressResult<Vec<_>>>();
        let forces = match forces {
            Ok(f) => f,
            Err(e) => {
                if recovery.is_some() {
                    return Err(restoration_failure(e, gate, poll().is_break()));
                }
                state.error = Some(e);
                return Ok(());
            }
        };
        let loads: Vec<_> = laws
            .iter()
            .zip(&forces)
            .map(|(law, force)| LoadCase {
                force,
                weight: law.weight,
            })
            .collect();
        if recovery.is_none() && let Err(e) = gradient_gate(
            &mut study,
            &loads,
            &raw,
            options,
            spec.stress_measure,
            &mut state.audit,
            &mut control,
        ) {
            state.error = Some(e);
            return Ok(());
        }
        let mut remaining_options = options;
        remaining_options.optimizer.max_evaluations -= state.audit.evaluations;
        if remaining_options.optimizer.max_evaluations == 0 {
            state.error = Some(exhausted());
            return Ok(());
        }
        // Initialization is one admitted callback attempt even if its solve
        // fails before a complete optimizer session can be returned.
        state.optimizer_work.evaluations = 1;
        let initialized = match recovery {
            None => StressDesignStudy3::new_with_measure(
                &mut study, &loads, &raw, remaining_options, spec.stress_measure, &mut control,
            ),
            Some(recovery) => StressDesignStudy3::restore_with_measure(
                &mut study, &loads, recovery.checkpoint.clone(), remaining_options, spec.stress_measure, &mut control,
            ),
        };
        let mut design = match initialized {
            Ok(design) => design,
            Err(e) => {
                if recovery.is_some() {
                    return Err(restoration_failure(e, gate, poll().is_break()));
                }
                state.error = Some(e);
                return Ok(());
            }
        };
        let result = design.run(0);
        state.capture(&design, result);
        if state.error.is_some()
            || state
                .report
                .as_ref()
                .is_some_and(|r| r.stop == ProjectedAlStop::Converged)
        {
            return Ok(());
        }
        state.spent =
            prior.add(design.work(), geometry, start.elapsed().as_secs_f64())?;
        accepted(design.study(), &state)?;
        let updates = limit.min(spec.updates - design.optimizer_work().iterations);
        for step in 0..updates {
            let result = design.run(1);
            state.capture(&design, result);
            if state.error.is_some()
                || state
                    .report
                    .as_ref()
                    .is_some_and(|r| r.stop != ProjectedAlStop::IterationLimit)
            {
                break;
            }
            if step + 1 < updates {
                state.spent =
                    prior.add(design.work(), geometry, start.elapsed().as_secs_f64())?;
                accepted(design.study(), &state)?;
            }
        }
        Ok(())
    })?;
    state.spent = prior.add(control.work(), geometry, start.elapsed().as_secs_f64())?;
    state.status = status(&state, gate, state.spent.wall_s >= spec.wall_s);
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

pub(super) fn resume(ledger: &Ledger, old: &Loaded, cap: Option<usize>, gate: &CancelGate) -> Result<Outcome> {
    let source = linked(ledger, &old.value, "source", "study-source")?;
    let source = std::str::from_utf8(&source)
        .map_err(|e| fail("cli-study-sdf3-stress-checkpoint", e.to_string()))?;
    let spec = spec::parse(source)?;
    resume::drive(&spec, ledger, cap, gate, Some(old))
}

#[cfg(test)]
#[path = "stress/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "stress/peak_tests.rs"]
mod peak_tests;
