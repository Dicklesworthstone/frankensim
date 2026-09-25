//! Real 3-D implicit-domain SIMP study through the canonical CLI/ledger surface.
//! The background adapts; the declared implicit boundary remains fixed.
use super::*;
use std::ops::ControlFlow;
use std::cell::Cell;

use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{
    AdaptivePreconditionError3, AdaptiveSolveSpace3,
};
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::{ElasticityError3, ElasticityOptions3};
use fs_cutfem::octree3::{Octree3, OctreeError3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureError3, QuadratureOptions3, QuadratureWork3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::op::two_level::{TwoLevelBudget, TwoLevelError};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::adaptive_continuation::{
    AdaptiveContinuationError3, AdaptiveContinuationOptions3, AdaptiveContinuationReport3,
    controlled_adaptive_sdf3_continuation_observed,
};
use fs_topopt::sdf3_goal::{
    GoalPreconditioner3, GoalReferenceLoad3, GoalRefinementError3, GoalRefinementOptions3,
};
use fs_topopt::{
    ContinuationTermination, EvaluationStop, MultiLoadOcOptions, SimpParams, SolveBudget,
    SolveControl, SolveWork,
};

#[path = "sdf3/checkpoint.rs"]
mod checkpoint;
use checkpoint::Spent;

#[path = "sdf3/output.rs"]
mod output;
#[path = "sdf3/spec.rs"]
mod spec;
use spec::Spec;

const SCOPE: &str = "Estimated 3-D linear-elastic SIMP compliance on a fixed raw implicit height field. The octree background is refined by numerical two-grid goal residuals; raw densities are transferred, physical-volume feasibility restored, and compliance/volume directional derivatives checked before every stage. These are discrete numerical comparisons, not continuum error enclosures, moving-boundary optimization, mesh-independent optima, manufacturing guarantees or physical validation. Compliance descent applies within each stage. Memory admission is an envelope, not measured RSS. Cancellation and wall time are checked at cooperative quadrature/solver boundaries; individual kernels and ledger I/O are indivisible. Every completed stage is durably retained before more physics. Resume verifies the retained stage prefix by replay under the same executable; replay spends the original wall, Krylov and geometry allowances. This is not direct optimizer-state restoration. Failed replay or a process crash preserves previous checkpoints but cannot durably charge later uncheckpointed work. The final ledger write is outside its own recorded wall measurement.";

struct Domain {
    height: f64,
    curvature: f64,
}
impl CutSdf3 for Domain {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - self.height - self.curvature * p[0] * (1.0 - p[0])
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        Interval::new(lo[2], hi[2])
            - Interval::new(self.height, self.height)
            - Interval::new(self.curvature, self.curvature) * x * (Interval::new(1.0, 1.0) - x)
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        match axis {
            HeightAxis::X => {
                Interval::new(-self.curvature, -self.curvature)
                    * (Interval::new(1.0, 1.0)
                        - Interval::new(2.0, 2.0) * Interval::new(lo[0], hi[0]))
            }
            HeightAxis::Y => Interval::new(0.0, 0.0),
            HeightAxis::Z => Interval::new(1.0, 1.0),
        }
    }
}

struct Computation {
    study: CutDensityStudy3<AdaptiveSolveSpace3>,
    report: AdaptiveContinuationReport3,
    tree: Octree3,
    spent: Spent,
    completed_stages: usize,
    requested_stages: usize,
    status: &'static str,
}

/// A borrowed accepted state avoids copying the operator just to persist it.
struct View<'a> {
    study: &'a CutDensityStudy3<AdaptiveSolveSpace3>,
    report: &'a AdaptiveContinuationReport3,
    tree: &'a Octree3,
    spent: Spent,
    completed_stages: usize,
    requested_stages: usize,
    status: &'static str,
}
impl Computation {
    fn view(&self) -> View<'_> {
        View {
            study: &self.study, report: &self.report, tree: &self.tree,
            spent: self.spent, completed_stages: self.completed_stages,
            requested_stages: self.requested_stages, status: self.status,
        }
    }
}

fn geometry_budget(error: &GoalRefinementError3) -> bool {
    matches!(
        error,
        GoalRefinementError3::Physics(ElasticityError3::Quadrature(
            QuadratureError3::BoxBudget | QuadratureError3::PointBudget
        ))
    )
}

fn refinement_budget(error: &AdaptiveContinuationError3) -> bool {
    match error {
        AdaptiveContinuationError3::Background(
            OctreeError3::LeafBudget | OctreeError3::LevelBudget,
        ) => true,
        AdaptiveContinuationError3::Goal(GoalRefinementError3::Preconditioner(
            AdaptivePreconditionError3::Coarse(TwoLevelError::Budget(_)),
        )) => true,
        AdaptiveContinuationError3::Goal(error) => geometry_budget(error),
        _ => false,
    }
}

#[cfg(test)]
fn compute(spec: &Spec, gate: &CancelGate) -> Result<Computation> {
    compute_observed(spec, gate, Spent::default(), spec.schedule.len(), |_| Ok(()))
}

fn compute_observed(
    spec: &Spec,
    gate: &CancelGate,
    prior: Spent,
    requested_stages: usize,
    mut accepted: impl FnMut(&View<'_>) -> Result<()>,
) -> Result<Computation> {
    if requested_stages == 0 || requested_stages > spec.schedule.len() {
        return Err(fail("cli-study-sdf3-stage-budget", "invalid stage prefix"));
    }
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
            code: if gate.is_requested() {
                "cli-study-sdf3-cancelled"
            } else {
                "cli-study-sdf3-wall-budget"
            },
            message: "study stopped before geometry work".into(),
            exit: if gate.is_requested() {
                exit::CANCELLED
            } else {
                exit::BUDGET
            },
        });
    }
    let mut tree = Octree3::uniform(spec.level as u8, spec.max_level as u8, spec.leaves)
        .map_err(|e| fail("cli-study-sdf3-background", e.to_string()))?;
    let domain = Domain {
        height: spec.height,
        curvature: spec.curvature,
    };
    let bounds = HexCell::try_new([0.0; 3], [1.0; 3])
        .map_err(|e| fail("cli-study-sdf3-domain", e.to_string()))?;
    let material = IsotropicElastic::new(spec.youngs, spec.poisson, 1.0)
        .map_err(|e| fail("cli-study-sdf3-material", e.to_string()))?;
    let geometry = Cell::new(QuadratureWork3::default());
    let mut build = |tree: &Octree3,
                     checkpoint: &mut dyn FnMut() -> ControlFlow<()>|
     -> std::result::Result<AdaptiveSolveSpace3, GoalRefinementError3> {
        let limits = QuadratureOptions3 {
            max_boxes: spec.boxes - prior.geometry.boxes - geometry.get().boxes,
            max_points: spec.points - prior.geometry.points - geometry.get().points,
            ..Default::default()
        };
        let mut gate = |_| {
            if poll().is_break() {
                ControlFlow::Break(())
            } else {
                checkpoint()
            }
        };
        let mut quadrature =
            QuadratureControl3::new(limits, &mut gate).map_err(ElasticityError3::from)?;
        let result = AdaptiveElasticity3::build(
            bounds,
            tree,
            &domain,
            &material,
            &|p| p[0] == 0.0,
            ElasticityOptions3 {
                max_cells: spec.leaves,
                max_dofs: 50_000,
                ..Default::default()
            },
            &mut quadrature,
        );
        let spent = quadrature.work();
        let old = geometry.get();
        let add = |a: usize, b: usize| a.checked_add(b)
            .ok_or(GoalRefinementError3::Invalid("geometry work counter overflow"));
        geometry.set(QuadratureWork3 {
            boxes: add(old.boxes, spent.boxes)?,
            points: add(old.points, spent.points)?,
            field_evaluations: add(old.field_evaluations, spent.field_evaluations)?,
        });
        Ok(AdaptiveSolveSpace3::jacobi(result?, 100_000_000))
    };
    let operator = build(&tree, &mut || ControlFlow::Continue(()))
        .map_err(|e| Failure {
            code: "cli-study-sdf3-geometry",
            exit: if gate.is_requested() { exit::CANCELLED }
                else if geometry_budget(&e) || prior.wall_s + start.elapsed().as_secs_f64() >= spec.wall_s { exit::BUDGET }
                else { exit::REFUSED },
            message: format!("initial 3-D geometry could not complete: {e}; no displacement or optimized design is available"),
        })?;
    let mut study = CutDensityStudy3::new(operator, spec.radius, spec.schedule[0]);
    let raw = vec![spec.density; study.cells()];
    // Independent constant reference body densities. Reintegrating these laws
    // on each background preserves physical forcing without nodal transfer.
    let laws: Vec<_> = spec
        .loads
        .iter()
        .map(|(force, _)| move |_: [f64; 3]| *force)
        .collect();
    let loads: Vec<_> = laws
        .iter()
        .zip(&spec.loads)
        .map(|(law, (_, weight))| GoalReferenceLoad3 {
            load: ReferenceLoad3::body(law),
            weight: *weight,
        })
        .collect();
    let mut checkpoint = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: spec.linear - prior.linear.linear_iterations,
            per_solve_iterations: spec.per_solve,
        },
        &mut checkpoint,
    );
    let mut completed_stages = 0;
    let report = controlled_adaptive_sdf3_continuation_observed(
        &mut study,
        &mut tree,
        &loads,
        &raw,
        &spec.schedule[..requested_stages],
        AdaptiveContinuationOptions3 {
            optimization: MultiLoadOcOptions {
                volume_fraction: spec.volume,
                move_limit: spec.move_limit,
                max_iterations: spec.updates,
                ..Default::default()
            },
            enrichment: GoalRefinementOptions3 {
                max_load_cases: 4,
                preconditioner: GoalPreconditioner3::TwoLevel {
                    budget: TwoLevelBudget {
                        max_fine_dofs: 50_000,
                        ..Default::default()
                    },
                    max_diagonal_contributions: 100_000_000,
                },
                ..Default::default()
            },
            marking_fraction: spec.marking,
            max_marks: spec.max_marks,
            ..Default::default()
        },
        &mut control,
        &mut build,
        |study, tree, report| {
            completed_stages = report.continuation.stages.len();
            let status = if completed_stages == spec.schedule.len() { "completed" }
                else if completed_stages == requested_stages { "budget-exhausted" }
                else { "checkpointed" };
            accepted(&View {
                study, tree, report,
                spent: prior.add(report.continuation.work, geometry.get(), start.elapsed().as_secs_f64())?,
                completed_stages, requested_stages, status,
            })
        },
    )?;
    let result = &report.continuation;
    let status = if result.termination == ContinuationTermination::ScheduleComplete {
        if completed_stages == spec.schedule.len() { "completed" } else { "budget-exhausted" }
    } else if gate.is_requested() {
        "cancelled"
    } else if prior.wall_s + start.elapsed().as_secs_f64() >= spec.wall_s
        || matches!(
            result.evaluation_stop,
            Some(EvaluationStop::LinearBudget { .. } | EvaluationStop::TotalBudget { .. })
        )
        || report
            .refinement_error
            .as_ref()
            .is_some_and(refinement_budget)
    {
        "budget-exhausted"
    } else {
        "numerical-failure"
    };
    Ok(Computation {
        spent: prior.add(report.continuation.work, geometry.get(), start.elapsed().as_secs_f64())?,
        study, report, tree, completed_stages, requested_stages, status,
    })
}

pub(super) fn study(
    source: &str,
    ledger_path: &Path,
    override_text: Option<&str>,
    gate: &CancelGate,
) -> Result<Outcome> {
    let cap = budget(override_text)?;
    let spec = spec::parse(source)?;
    let ledger = Ledger::open(
        ledger_path
            .to_str()
            .ok_or_else(|| fail("cli-study-sdf3-ledger", "ledger path is not UTF-8"))?,
    )?;
    checkpoint::drive(&spec, &ledger, cap, gate, None, checkpoint::announce)
}

pub(super) fn resume(
    ledger: &Ledger,
    old: &Loaded,
    override_text: Option<&str>,
    gate: &CancelGate,
) -> Result<Outcome> {
    let cap = budget(override_text)?;
    let source = linked(ledger, &old.value, "source", "study-source")?;
    let source = std::str::from_utf8(&source)
        .map_err(|error| fail("cli-study-sdf3-receipt", error.to_string()))?;
    let spec = spec::parse(source)?;
    checkpoint::drive(&spec, ledger, cap, gate, Some(old), checkpoint::announce)
}

#[cfg(test)]
#[path = "sdf3/tests.rs"]
mod tests;
