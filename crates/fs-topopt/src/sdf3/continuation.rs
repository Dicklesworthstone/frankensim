//! SIMP/projection continuation on the existing 3-D CutFEM density pipeline.
//!
//! Geometry, quadrature, supports and loads stay fixed across stages. Each new
//! material model restores actual projected-volume feasibility before its own
//! baseline solve; compliance comparisons only apply WITHIN that model. The
//! optional gradient gate uses the same finite-difference experiment as the
//! tetrahedral pipeline, evaluated through real independent CutFEM solves.

use super::{
    CutDensityStudy3, Sdf3Elasticity, assert_loads, assert_oc_inputs,
    controlled_sdf3_optimality_criteria,
};
use crate::continuation::{
    ContinuationStageReport, ContinuationTermination, MultiLoadContinuationReport,
};
use crate::gradient_check::{
    GradientCheckOptions, GradientSample, MultiLoadGradientCheck, audit_design_map,
};
use crate::pipeline::LoadCase;
use crate::{EvaluationStop, MultiLoadOcOptions, MultiLoadOcTermination, SimpParams, SolveControl};

/// Audit compliance and physical-volume derivatives through the complete
/// graph-filter/projection/SIMP/CutFEM chain. At most five independent load
/// families are solved at a single model and design. The heterogeneous,
/// bound-aware directions collectively touch every density coordinate.
///
/// Incoming scales are restored on success and every returned numerical stop.
/// All probe work counts against `control`; cancellation never publishes a
/// partial check. Modeling-input panics follow the existing density APIs and
/// occur before operator mutation. A passing check is directional numerical
/// evidence, not a continuum error bound or a componentwise gradient proof.
pub fn controlled_sdf3_gradient_check<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    loads: &[LoadCase<'_>],
    rho: &[f64],
    options: GradientCheckOptions,
    control: &mut SolveControl<'_>,
) -> Result<MultiLoadGradientCheck, EvaluationStop> {
    options.assert_valid();
    study.params.assert_valid();
    assert_loads(study, loads);
    assert_eq!(
        rho.len(),
        study.cells(),
        "one raw density per active cut cell required"
    );
    assert!(
        !rho.is_empty() && rho.iter().all(|r| r.is_finite() && (0.0..=1.0).contains(r)),
        "gradient-check densities must be nonempty, finite, and lie in [0, 1]"
    );
    let previous = study.operator.scales().to_vec();
    let result = audit_design_map(study.params, rho, options, control, |point, _, control| {
        let evaluated = study.evaluate(point, loads, control)?;
        Ok(GradientSample {
            compliance: evaluated.objective.compliance,
            compliance_gradient: evaluated.objective.gradient,
            volume: evaluated.volume_fraction,
            volume_gradient: evaluated.volume_gradient,
        })
    });
    study
        .operator
        .set_scales(&previous)
        .expect("previous scales are valid");
    result
}

/// Execute a finite continuation schedule on Cartesian or adaptive CutFEM.
/// `options.max_iterations` is the accepted-update allowance PER stage.
/// Projection/material changes receive a fresh volume-feasible equilibrium;
/// reports never compare compliance across different material models.
///
/// This borrows the existing geometry-bound pipeline and changes only density
/// parameters/scales: no geometry reconstruction, meshing, force transfer or
/// quadrature is performed. Geometry refinement is a separate between-study
/// operation. A stop restores both parameters and stiffness scales matching
/// the last fully solved accepted report, or the incoming state if none exists.
/// Invalid modeling inputs panic before any model or operator mutation.
pub fn controlled_sdf3_continuation<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    schedule: &[SimpParams],
    options: MultiLoadOcOptions,
    control: &mut SolveControl<'_>,
) -> MultiLoadContinuationReport {
    run(study, loads, rho0, schedule, options, None, control)
}

/// The same continuation with a required numerical compliance AND volume
/// gradient audit at every restored stage baseline. A failed gate retains its
/// measured discrepancies separately and cannot replace the previous model.
/// Volume restoration, all probes, state solves and rejected optimizer trials
/// share one cumulative budget. A completed schedule is not convergence.
pub fn controlled_gradient_checked_sdf3_continuation<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    schedule: &[SimpParams],
    options: MultiLoadOcOptions,
    gradient_options: GradientCheckOptions,
    control: &mut SolveControl<'_>,
) -> MultiLoadContinuationReport {
    gradient_options.assert_valid();
    run(
        study,
        loads,
        rho0,
        schedule,
        options,
        Some(gradient_options),
        control,
    )
}

fn run<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>,
    loads: &[LoadCase<'_>],
    rho0: &[f64],
    schedule: &[SimpParams],
    options: MultiLoadOcOptions,
    gradient_options: Option<GradientCheckOptions>,
    control: &mut SolveControl<'_>,
) -> MultiLoadContinuationReport {
    assert!(
        !schedule.is_empty(),
        "continuation requires at least one stage"
    );
    study.params.assert_valid();
    for params in schedule {
        params.assert_valid();
    }
    assert_oc_inputs(study, loads, rho0, options);

    let mut accepted_scales = study.operator.scales().to_vec();
    let mut rho = rho0.to_vec();
    let mut report = MultiLoadContinuationReport {
        last: None,
        params: study.params,
        stages: Vec::new(),
        termination: ContinuationTermination::ScheduleComplete,
        stopped_stage: None,
        evaluation_stop: None,
        rejected_gradient_check: None,
        work: control.work(),
    };
    for (stage, &params) in schedule.iter().enumerate() {
        let start = (|| {
            control.checkpoint("sdf3-continuation-stage")?;
            study.params = params;
            study.restore_start_with_floor(
                &rho,
                options.volume_fraction,
                options.volume_tolerance,
                1e-3,
                control,
            )
        })();
        let start = match start {
            Ok(start) => start,
            Err(stop) => {
                report.termination = if matches!(
                    stop,
                    EvaluationStop::Breakdown {
                        stage: "sdf3-volume-restoration"
                    }
                ) {
                    ContinuationTermination::FeasibilityRestorationFailed
                } else {
                    ContinuationTermination::EvaluationStopped
                };
                report.stopped_stage = Some(stage);
                report.evaluation_stop = Some(stop);
                break;
            }
        };
        let gradient_check = if let Some(check_options) = gradient_options {
            match controlled_sdf3_gradient_check(study, loads, &start.rho, check_options, control) {
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
        } else {
            None
        };
        let run = controlled_sdf3_optimality_criteria(study, loads, &start.rho, options, control);
        let stopped = matches!(
            run.termination,
            MultiLoadOcTermination::Cancelled
                | MultiLoadOcTermination::LinearBudget
                | MultiLoadOcTermination::NumericalFailure
        );
        if stopped {
            report.termination = ContinuationTermination::EvaluationStopped;
            report.stopped_stage = Some(stage);
            report.evaluation_stop = run.evaluation_stop.clone();
        }
        if !run.history.is_empty() {
            report.params = params;
            accepted_scales = study.operator.scales().to_vec();
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
    study.params = report.params;
    study
        .operator
        .set_scales(&accepted_scales)
        .expect("accepted scales are valid");
    report.work = control.work();
    report
}
