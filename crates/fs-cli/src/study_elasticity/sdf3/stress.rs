//! Canonical minimum-volume studies using the existing exact-discrete stress
//! adjoint and projected augmented-Lagrangian optimizer. The domain is fixed.
use super::*;
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlReport, ProjectedAlStop, ProjectedAlWork,
};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::design::{StressDesignIteration3, StressDesignOptions3, StressDesignStudy3};
use fs_topopt::sdf3::stress::{StressError3, StressEvaluation3};

#[path = "stress/output.rs"]
mod output;

const SCOPE: &str = "Estimated fixed-background 3-D linear-elastic SIMP minimum-volume design. The constraint is a normalized volume-and-load-weighted qp von Mises aggregate in Pa; it is not a limit on sampled or continuum maximum stress. Independent body, reference-pressure and traction cases remain separate. The initial stress and volume adjoints pass bounded directional finite differences before optimization. Accepted augmented-Lagrangian steps can be infeasible; the least-volume accepted feasible design is retained separately and selected for export when available. KKT residuals describe the last accepted iterate, not an earlier exported incumbent. Only numerical convergence with a feasible current iterate reports completed. Every accepted update is durable before further optimization. Cross-process optimizer resume is unsupported; report and package export all retained results without another solve. No continuum safety, global optimum, manufacturing, adaptivity or physical-validation claim is made. Memory is an admission envelope, not measured RSS. Wall time and cancellation are cooperative; quadrature, solver, stress-cell and optimizer boundaries poll, while individual kernels and ledger I/O are indivisible.";
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

#[derive(Default)]
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
            .evaluate_stress(point, loads, options.stress, control)
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

fn compute_observed(
    spec: &Spec,
    gate: &CancelGate,
    limit: usize,
    mut accepted: impl FnMut(&CutDensityStudy3<AdaptiveSolveSpace3>, &State) -> Result<()>,
) -> Result<Computation> {
    let options = spec
        .stress
        .ok_or_else(|| fail("cli-study-sdf3-input", "missing stress optimizer"))?;
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
            max_boxes: spec.boxes,
            max_points: spec.points,
            ..Default::default()
        },
        &mut cp,
    )
    .map_err(|e| fail("cli-study-sdf3-geometry", e.to_string()))?;
    let operator = loading::build_operator(spec, bounds, &tree, domain, &material, &mut quadrature)
        .map_err(|e| Failure {
            code: "cli-study-sdf3-geometry",
            message: e.to_string(),
            exit: if gate.is_requested() {
                exit::CANCELLED
            } else if poll().is_break()
                || matches!(
                    e,
                    ElasticityError3::Quadrature(
                        QuadratureError3::BoxBudget | QuadratureError3::PointBudget
                    )
                )
            {
                exit::BUDGET
            } else {
                exit::REFUSED
            },
        })?;
    let geometry = quadrature.work();
    let mut study = CutDensityStudy3::new(
        AdaptiveSolveSpace3::jacobi(operator, 100_000_000),
        spec.radius,
        spec.schedule[0],
    );
    let raw = vec![spec.density; study.cells()];
    let mut cp = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: spec.linear,
            per_solve_iterations: spec.per_solve,
        },
        &mut cp,
    );
    let mut state = State::new();
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
        if let Err(e) = gradient_gate(
            &mut study,
            &loads,
            &raw,
            options,
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
        let mut design = match StressDesignStudy3::new(
            &mut study,
            &loads,
            &raw,
            remaining_options,
            &mut control,
        ) {
            Ok(design) => design,
            Err(e) => {
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
            Spent::default().add(design.work(), geometry, start.elapsed().as_secs_f64())?;
        accepted(design.study(), &state)?;
        for step in 0..limit.min(spec.updates) {
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
            if step + 1 < limit.min(spec.updates) {
                state.spent =
                    Spent::default().add(design.work(), geometry, start.elapsed().as_secs_f64())?;
                accepted(design.study(), &state)?;
            }
        }
        Ok(())
    })?;
    state.spent = Spent::default().add(control.work(), geometry, start.elapsed().as_secs_f64())?;
    state.status = status(&state, gate, start.elapsed().as_secs_f64() >= spec.wall_s);
    Ok(Computation { study, state })
}

pub(super) fn drive(
    spec: &Spec,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
) -> Result<Outcome> {
    let limit = cap.unwrap_or(spec.updates).min(spec.updates);
    let mut predecessor = None;
    let result = compute_observed(spec,gate,limit,|study,state| {
        let out = output::persist(spec,ledger,study,state,predecessor,limit)?;
        predecessor = out.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        use std::io::Write as _;
        let _ = writeln!(std::io::stderr().lock(),
            "{{\"schema\":\"frankensim.cli.sdf3-stress-progress.v1\",\"run_id\":{},\"iterations_completed\":{}}}",
            quoted(&out.pointer), state.iterations());
        Ok(())
    }).map_err(|mut e| {
        if let Some(previous) = predecessor { e.message.push_str(&format!("; last durable result: study-{}",previous.to_hex())); }
        e
    })?;
    output::persist(
        spec,
        ledger,
        &result.study,
        &result.state,
        predecessor,
        limit,
    )
}

#[cfg(test)]
#[path = "stress/tests.rs"]
mod tests;
