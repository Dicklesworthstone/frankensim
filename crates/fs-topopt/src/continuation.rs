//! Volume-feasible continuation of independently loaded SIMP studies.
//!
//! Each stage re-solves equilibrium under its own material/projection model.
//! Compliance descent is checked WITHIN a stage, never across different SIMP
//! models. A sharper projection can violate the old material budget; transition
//! restoration therefore evaluates physical volume before admitting a baseline.
//! All transition, filter, equilibrium and rejected-trial work shares one
//! `SolveControl`. A returned stop retains only a completely solved design.
//!
//! This is a fixed-mesh density method, not the raw-SDF/CutFEM marquee, a
//! stationarity certificate, or a proof of mesh-independent manufacturability.

use std::ops::ControlFlow;

use crate::control::{EvaluationStop, SolveBudget, SolveControl, SolveWork};
use crate::elasticity::DensityElasticity;
use crate::gradient_check::{
    GradientCheckOptions, MultiLoadGradientCheck, controlled_multi_load_gradient_check,
};
use crate::multi_load::{
    MultiLoadOcIteration, MultiLoadOcOptions, MultiLoadOcReport, MultiLoadOcTermination,
    controlled_multi_load_optimality_criteria,
};
use crate::oc::assert_valid_oc_inputs;
use crate::pipeline::{DesignPipeline, LoadCase, SimpParams, assert_valid_load_cases};

const DENSITY_FLOOR: f64 = 1e-3;

/// Why a continuation study ended. No variant asserts optimality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuationTermination {
    /// Every scheduled model received its bounded optimization stage.
    /// Individual stages may have stopped without finding an acceptable step.
    ScheduleComplete,
    /// Cancellation, exhausted linear work, or numerical failure.
    EvaluationStopped,
    /// Even the tested density-floor design exceeded the material cap.
    /// This is a failed restoration, not a global infeasibility certificate.
    FeasibilityRestorationFailed,
    /// A complete stage-start numerical gradient audit exceeded its tolerance.
    GradientCheckFailed,
}

/// History for one stage with at least one complete equilibrium evaluation.
#[derive(Debug, Clone)]
pub struct ContinuationStageReport {
    /// Zero-based index in the requested schedule.
    pub stage: usize,
    /// Exact model used for every row in this stage.
    pub params: SimpParams,
    /// The incoming raw design's volume under THIS stage's projection.
    pub incoming_volume_fraction: f64,
    /// Retained fraction of the incoming design's distance from the density
    /// floor. One means no restoration; zero means the floor design.
    pub restoration_scale: f64,
    /// Optional numerical audit at this stage's restored starting design.
    /// None means no runtime gradient gate was requested. It does not describe
    /// later accepted iterates or certify individual gradient components.
    pub gradient_check: Option<MultiLoadGradientCheck>,
    /// Aligned accepted-design history, including the new stage's baseline.
    pub history: Vec<MultiLoadOcIteration>,
    /// The fixed-stage driver's actual stopping reason.
    pub termination: MultiLoadOcTermination,
}

/// Accepted continuation prefix and the last complete physical state.
#[derive(Debug, Clone)]
pub struct MultiLoadContinuationReport {
    /// Last stage's accepted design/fields. None means no baseline was solved;
    /// the caller's incoming parameters and operator then remain unchanged.
    /// Its work counter is the snapshot at that stage's return; `work` below
    /// also includes any later, rejected stage preparation.
    pub last: Option<MultiLoadOcReport>,
    /// Parameters matching `last`, or the original parameters when it is None.
    pub params: SimpParams,
    /// Only stages having a complete baseline, in schedule order.
    pub stages: Vec<ContinuationStageReport>,
    /// Overall stopping reason, independent of individual stage convergence.
    pub termination: ContinuationTermination,
    /// Stage that could not finish; None after the whole schedule was consumed.
    pub stopped_stage: Option<usize>,
    /// Detailed numerical/cancellation stop, when applicable.
    pub evaluation_stop: Option<EvaluationStop>,
    /// Failed numerical audit of the rejected stage, never an accepted design.
    pub rejected_gradient_check: Option<MultiLoadGradientCheck>,
    /// All consumed work, including failed transitions and rejected trials.
    pub work: SolveWork,
}

struct StageStart {
    rho: Vec<f64>,
    incoming_volume: f64,
    scale: f64,
}

fn projected_volume(
    pipeline: &DesignPipeline,
    rho: &[f64],
    cell_vol: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<f64, EvaluationStop> {
    let (_, projected, _) = pipeline.try_forward(rho, control)?;
    let total: f64 = cell_vol.iter().sum();
    let value: f64 = projected
        .iter()
        .zip(cell_vol)
        .map(|(r, v)| r * (v / total))
        .sum();
    if !value.is_finite() {
        return Err(EvaluationStop::Breakdown {
            stage: "continuation-volume",
        });
    }
    Ok(value)
}

fn restore_start(
    pipeline: &DesignPipeline,
    rho: &[f64],
    cell_vol: &[f64],
    options: MultiLoadOcOptions,
    control: &mut SolveControl<'_>,
) -> Result<Option<StageStart>, EvaluationStop> {
    let cap = options.volume_fraction + options.volume_tolerance;
    let incoming_volume = projected_volume(pipeline, rho, cell_vol, control)?;
    if incoming_volume <= cap {
        return Ok(Some(StageStart {
            rho: rho.to_vec(),
            incoming_volume,
            scale: 1.0,
        }));
    }
    let mut feasible = vec![DENSITY_FLOOR; rho.len()];
    if projected_volume(pipeline, &feasible, cell_vol, control)? > cap {
        return Ok(None);
    }
    let (mut low, mut high) = (0.0, 1.0);
    // Retain an actually evaluated feasible endpoint, not a multiplier root
    // rounded toward infeasibility. This remains safe even if finite-element
    // filtering does not give a globally monotone volume along the segment.
    for _ in 0..64 {
        control.checkpoint("continuation-volume-search")?;
        let scale = low + 0.5 * (high - low);
        if scale <= low || scale >= high {
            break;
        }
        let trial: Vec<f64> = rho
            .iter()
            .map(|r| (DENSITY_FLOOR + scale * (r - DENSITY_FLOOR)).clamp(DENSITY_FLOOR, 1.0))
            .collect();
        if projected_volume(pipeline, &trial, cell_vol, control)? <= cap {
            low = scale;
            feasible = trial;
        } else {
            high = scale;
        }
    }
    Ok(Some(StageStart {
        rho: feasible,
        incoming_volume,
        scale: low,
    }))
}

/// Run a schedule with default linear-work limits and a cancellation hook.
/// `options.max_iterations` is the accepted-update allowance for EACH stage.
#[allow(clippy::too_many_arguments)]
pub fn multi_load_continuation(
    pipeline: &mut DesignPipeline,
    elasticity: &mut DensityElasticity,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    cell_vol: &[f64],
    schedule: &[SimpParams],
    options: MultiLoadOcOptions,
    mut checkpoint: impl FnMut() -> ControlFlow<()>,
) -> MultiLoadContinuationReport {
    let mut callback = |_| checkpoint();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    controlled_multi_load_continuation(
        pipeline, elasticity, loads, rho0, cell_vol, schedule, options, &mut control,
    )
}

/// Run a continuation schedule under one cumulative work/cancellation budget.
///
/// The entire schedule and model inputs are admitted before mutating the
/// operator. Each transition first tests the incoming design under the NEW
/// projection. If necessary, a bounded search toward the density floor retains
/// a physically volume-feasible start. Restoration may increase compliance;
/// stage history records a fresh equilibrium baseline rather than claiming
/// descent across different objectives. Within a stage the existing OC
/// feasibility and non-increasing-compliance acceptance rules are unchanged.
///
/// Every returned stop restores parameters AND moduli to the last accepted
/// state, including a partial stage with a fully solved baseline. Failed
/// preparation of a later stage does not replace earlier accepted fields.
/// Invalid modeling inputs retain the component APIs' panic contract.
#[allow(clippy::too_many_arguments)]
pub fn controlled_multi_load_continuation(
    pipeline: &mut DesignPipeline,
    elasticity: &mut DensityElasticity,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    cell_vol: &[f64],
    schedule: &[SimpParams],
    options: MultiLoadOcOptions,
    control: &mut SolveControl<'_>,
) -> MultiLoadContinuationReport {
    run_continuation(pipeline, elasticity, loads, rho0, cell_vol, schedule, options, None, control)
}

/// Continuation with a real compliance AND physical-volume gradient gate at
/// every stage's restored starting design. A failed finite-difference audit
/// stops before the new stage can replace the previously accepted model.
/// Numerical checks, restoration and optimization all consume the same budget.
/// A successful directional audit is not a continuum or optimality certificate.
#[allow(clippy::too_many_arguments)]
pub fn controlled_gradient_checked_multi_load_continuation(
    pipeline: &mut DesignPipeline,
    elasticity: &mut DensityElasticity,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    cell_vol: &[f64],
    schedule: &[SimpParams],
    options: MultiLoadOcOptions,
    gradient_options: GradientCheckOptions,
    control: &mut SolveControl<'_>,
) -> MultiLoadContinuationReport {
    gradient_options.assert_valid();
    run_continuation(pipeline, elasticity, loads, rho0, cell_vol, schedule, options,
        Some(gradient_options), control)
}

#[allow(clippy::too_many_arguments)]
fn run_continuation(
    pipeline: &mut DesignPipeline,
    elasticity: &mut DensityElasticity,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    cell_vol: &[f64],
    schedule: &[SimpParams],
    options: MultiLoadOcOptions,
    gradient_options: Option<GradientCheckOptions>,
    control: &mut SolveControl<'_>,
) -> MultiLoadContinuationReport {
    assert!(!schedule.is_empty(), "continuation requires at least one stage");
    pipeline.params.assert_valid();
    for params in schedule {
        params.assert_valid();
    }
    assert_valid_load_cases(elasticity, loads);
    assert_valid_oc_inputs(
        elasticity, loads[0].force, rho0, cell_vol,
        options.volume_fraction, options.move_limit,
    );
    assert!(rho0.iter().all(|r| *r >= DENSITY_FLOOR),
        "OC starting densities must be at least 1e-3");
    assert!(options.change_tolerance.is_finite() && options.change_tolerance >= 0.0,
        "design-change tolerance must be finite and nonnegative");
    assert!(options.volume_tolerance.is_finite() && options.volume_tolerance >= 0.0
        && options.volume_tolerance < options.volume_fraction,
        "volume tolerance must be finite, nonnegative, and smaller than the cap");
    assert!(options.max_backtracks <= 64, "at most 64 OC backtracks are admitted");

    let mut accepted_moduli = elasticity.moduli.clone();
    let mut rho = rho0.to_vec();
    let mut report = MultiLoadContinuationReport {
        last: None,
        params: pipeline.params,
        stages: Vec::new(),
        termination: ContinuationTermination::ScheduleComplete,
        stopped_stage: None,
        evaluation_stop: None,
        rejected_gradient_check: None,
        work: control.work(),
    };
    for (stage, &params) in schedule.iter().enumerate() {
        let start = (|| {
            control.checkpoint("continuation-stage")?;
            pipeline.params = params;
            restore_start(pipeline, &rho, cell_vol, options, control)
        })();
        let start = match start {
            Ok(Some(start)) => start,
            Ok(None) => {
                report.termination = ContinuationTermination::FeasibilityRestorationFailed;
                report.stopped_stage = Some(stage);
                break;
            }
            Err(stop) => {
                report.termination = ContinuationTermination::EvaluationStopped;
                report.stopped_stage = Some(stage);
                report.evaluation_stop = Some(stop);
                break;
            }
        };
        let gradient_check = if let Some(check_options) = gradient_options {
            match controlled_multi_load_gradient_check(
                pipeline, elasticity, loads, &start.rho, cell_vol, check_options, control,
            ) {
                Ok(check) if check.passed() => Some(check),
                Ok(check) => {
                    report.termination = ContinuationTermination::GradientCheckFailed;
                    report.stopped_stage = Some(stage);
                    report.rejected_gradient_check = Some(check);
                    break;
                }
                Err(stop) => {
                    report.termination = ContinuationTermination::EvaluationStopped;
                    report.stopped_stage = Some(stage);
                    report.evaluation_stop = Some(stop);
                    break;
                }
            }
        } else { None };
        let run = controlled_multi_load_optimality_criteria(
            pipeline, elasticity, loads, &start.rho, cell_vol, options, control,
        );
        let stopped = matches!(run.termination,
            MultiLoadOcTermination::Cancelled | MultiLoadOcTermination::LinearBudget
                | MultiLoadOcTermination::NumericalFailure);
        if stopped {
            report.termination = ContinuationTermination::EvaluationStopped;
            report.stopped_stage = Some(stage);
            report.evaluation_stop = run.evaluation_stop.clone();
        }
        if !run.history.is_empty() {
            report.params = params;
            accepted_moduli = elasticity.moduli.clone();
            rho.clone_from(&run.rho);
            report.stages.push(ContinuationStageReport {
                stage,
                params,
                incoming_volume_fraction: start.incoming_volume,
                restoration_scale: start.scale,
                gradient_check,
                history: run.history.clone(),
                termination: run.termination,
            });
            report.last = Some(run);
        }
        if stopped {
            break;
        }
    }
    pipeline.params = report.params;
    elasticity.moduli = accepted_moduli;
    report.work = control.work();
    report
}
