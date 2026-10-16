//! Gradient-checked continuation of the three-field material model.
//!
//! Only raw controls cross a stage boundary. Dilated volume is restored under
//! the NEW projection before a new eroded equilibrium and gradient audit. Each
//! stage starts a fresh AL problem; old multipliers, gradients and compliance
//! values are never reused under a different material model.

use super::{RobustDesignIteration3, RobustDesignOptions3, RobustDesignPoint3,
    RobustDesignStudy3, RobustGradientCheck3, controlled_sdf3_robust_gradient_check};
use super::super::{CutDensityStudy3, Sdf3Elasticity, assert_loads};
use crate::pipeline::LoadCase;
use crate::{EvaluationStop, GradientCheckOptions, SimpParams, SolveControl, SolveWork};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlReport, ProjectedAlStop};

type DesignError = ProjectedAlError<EvaluationStop>;

#[derive(Debug, Clone, Copy)]
pub struct RobustContinuationOptions3 {
    /// Fixed uncertainty, raw box, material cap and objective scale. AL limits
    /// apply PER stage, because each new material model is a fresh problem.
    /// All stages, probes and restorations share the supplied physical budget.
    pub design: RobustDesignOptions3,
    pub steps_per_stage: usize,
    /// Required numerical eroded-compliance AND dilated-volume audit at each
    /// restored baseline; it is not a continuum or componentwise certificate.
    pub gradient: GradientCheckOptions,
}
impl Default for RobustContinuationOptions3 {
    fn default() -> Self {
        Self { design: RobustDesignOptions3::default(), steps_per_stage: 100,
            gradient: GradientCheckOptions::default() }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobustContinuationStop3 {
    /// All models received their bounded stage. This is not convergence.
    ScheduleComplete,
    /// The tested raw-floor design exceeded the new dilated cap. This does
    /// not prove that every admissible design is infeasible.
    RestorationFailed,
    GradientCheckFailed,
    EvaluationStopped,
    /// Includes stalls and resource limits; none establishes infeasibility.
    OptimizerStopped(ProjectedAlStop),
}

/// One admitted model, with a freshly solved and gradient-checked baseline.
#[derive(Debug, Clone)]
pub struct RobustContinuationStage3 {
    pub stage: usize,
    pub params: SimpParams,
    pub incoming_dilated_volume_fraction: f64,
    /// Fraction retained along floor + scale * (incoming - floor).
    pub restoration_scale: f64,
    pub starting_rho: Vec<f64>,
    pub gradient_check: RobustGradientCheck3,
    /// Actual accepted AL trajectory, including potentially infeasible steps.
    pub history: Vec<RobustDesignIteration3>,
    /// Index in history of the BEST FEASIBLE accepted point. Only its raw
    /// controls initialize the next stage; it need not be the final AL iterate.
    pub selected_iteration: usize,
    /// Applies to the final AL iterate, not necessarily the selected point.
    /// Late interruptions retain this stage's complete accepted prefix.
    pub outcome: Result<ProjectedAlReport, DesignError>,
}

#[derive(Debug, Clone)]
pub struct RobustContinuationReport3 {
    pub stages: Vec<RobustContinuationStage3>,
    /// Feasible incumbent of the latest admitted model, not a cross-model
    /// minimum. None means no stage baseline was admitted.
    pub last: Option<RobustDesignPoint3>,
    /// Model matching last; the incoming model when last is None.
    pub params: SimpParams,
    pub options: RobustContinuationOptions3,
    pub stop: RobustContinuationStop3,
    pub stopped_stage: Option<usize>,
    pub error: Option<DesignError>,
    pub rejected_gradient_check: Option<RobustGradientCheck3>,
    /// Includes failed transitions and all gradient probes, not just AL work.
    pub work: SolveWork,
}

/// Restore actual DILATED-volume feasibility without running elasticity.
/// Uses the same filter, bisection, raw floor and physical-region overrides as
/// nominal restoration, with the explicit dilated projection threshold. The
/// strict cap leaves the optimizer's normalized tolerance for numerical slack.
/// Neither the incoming nominal parameters nor stiffness scales are changed.
/// Invalid modeling inputs panic before mutation, like feasible_start.
pub fn controlled_sdf3_robust_feasible_start<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, incoming: &[f64], options: RobustDesignOptions3,
    control: &mut SolveControl<'_>,
) -> Result<Vec<f64>, EvaluationStop> {
    admit_design(study, incoming, options).expect("invalid robust restoration policy");
    restore(study, incoming, options, control).map(|start| start.rho)
}

fn admit_design<O: Sdf3Elasticity>(study: &CutDensityStudy3<O>, rho: &[f64],
    options: RobustDesignOptions3) -> Result<(), DesignError> {
    study.params().assert_valid();
    options.optimizer.validate(study.cells())?;
    if !options.robustness.valid(study.params())
        || !options.volume_limit.is_finite() || options.volume_limit <= 0.0 || options.volume_limit > 1.0
        || !options.compliance_scale.is_finite() || options.compliance_scale <= 0.0
        || !options.density_floor.is_finite() || options.density_floor <= 0.0 || options.density_floor >= 1.0
        || rho.len() != study.cells()
        || rho.iter().any(|r| !r.is_finite() || *r < options.density_floor || *r > 1.0) {
        return Err(ProjectedAlError::Invalid("invalid robust continuation policy or raw design"));
    }
    Ok(())
}

fn restore<O: Sdf3Elasticity>(study: &mut CutDensityStudy3<O>, rho: &[f64],
    options: RobustDesignOptions3, control: &mut SolveControl<'_>)
    -> Result<super::super::refine::RestoredDensity3, EvaluationStop> {
    let nominal = study.params;
    study.params.eta -= options.robustness.threshold_offset;
    let result = study.restore_start_with_floor(rho, options.volume_limit, 0.0,
        options.density_floor, control);
    study.params = nominal;
    result
}

/// Fixed-geometry continuation on Cartesian or conforming adaptive CutFEM.
/// Independent dead loads, physical regions and filter radius remain fixed.
/// Every stage restores the NEW dilated cap, audits its exact numerical map,
/// and initializes the existing robust optimizer from scratch. The selected
/// within-stage feasible incumbent is the only state passed to the next model.
///
/// The borrowed study keeps its INCOMING nominal model and stiffness scales on
/// success and every returned stop. Read the returned eroded fields instead of
/// interpreting those nominal scales as the robust solution. Failed transition
/// preparation cannot replace the preceding admitted model or physical fields.
/// No geometry rebuild, certified length-scale claim or cross-model descent.
pub fn controlled_sdf3_robust_continuation<O: Sdf3Elasticity>(
    study: &mut CutDensityStudy3<O>, loads: &[LoadCase<'_>], rho0: &[f64],
    schedule: &[SimpParams], options: RobustContinuationOptions3,
    control: &mut SolveControl<'_>,
) -> Result<RobustContinuationReport3, DesignError> {
    admit_design(study, rho0, options.design)?;
    if schedule.is_empty() { return Err(ProjectedAlError::Invalid("empty robust continuation schedule")); }
    // Validate the ENTIRE schedule before changing even the first model.
    for &params in schedule { options.design.robustness.assert_valid(params); }
    options.gradient.assert_valid();
    assert_loads(study, loads);
    let incoming_params = study.params;
    let result = run(study, loads, rho0, schedule, options, control);
    study.params = incoming_params;
    Ok(result)
}

fn run<O: Sdf3Elasticity>(study: &mut CutDensityStudy3<O>, loads: &[LoadCase<'_>],
    rho0: &[f64], schedule: &[SimpParams], options: RobustContinuationOptions3,
    control: &mut SolveControl<'_>) -> RobustContinuationReport3 {
    let mut report = RobustContinuationReport3 {
        stages: Vec::new(), last: None, params: study.params, options,
        stop: RobustContinuationStop3::ScheduleComplete, stopped_stage: None,
        error: None, rejected_gradient_check: None, work: control.work(),
    };
    let mut rho = rho0.to_vec();
    for (stage, &params) in schedule.iter().enumerate() {
        let attempt = (|| -> Result<(), DesignError> {
            control.checkpoint("sdf3-robust-continuation-stage").map_err(ProjectedAlError::Evaluation)?;
            study.params = params;
            let start = match restore(study, &rho, options.design, control) {
                Ok(start) => start,
                Err(error) => {
                    if matches!(error, EvaluationStop::Breakdown { stage: "sdf3-volume-restoration" }) {
                        report.stop = RobustContinuationStop3::RestorationFailed;
                    }
                    return Err(ProjectedAlError::Evaluation(error));
                }
            };
            let check = controlled_sdf3_robust_gradient_check(study, loads, &start.rho,
                options.design.robustness, options.gradient, control).map_err(ProjectedAlError::Evaluation)?;
            if !check.numerical.passed() {
                report.stop = RobustContinuationStop3::GradientCheckFailed;
                report.rejected_gradient_check = Some(check);
                return Ok(());
            }
            let mut session = RobustDesignStudy3::new(study, loads, &start.rho, options.design, control)?;
            // Restoration is verified AGAIN by the actual robust baseline map.
            if session.best_feasible().is_none() {
                report.stop = RobustContinuationStop3::RestorationFailed;
                return Err(ProjectedAlError::Invalid("restored robust baseline violates the dilated cap"));
            }
            let outcome = session.run(options.steps_per_stage);
            let selected = session.best_feasible().expect("feasible baseline retained").clone();
            let history = session.history().to_vec();
            // Match the session's strict-improvement incumbent rule, including
            // the earliest point in the event of equal objective values.
            let mut selected_iteration = 0;
            for (index, row) in history.iter().enumerate().skip(1) {
                if row.feasible && row.eroded_compliance < history[selected_iteration].eroded_compliance {
                    selected_iteration = index;
                }
            }
            report.work = session.work();
            drop(session);
            match &outcome {
                Err(error) => { report.error = Some(error.clone()); report.stop = RobustContinuationStop3::EvaluationStopped; }
                Ok(done) if !matches!(done.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit) => {
                    report.stop = RobustContinuationStop3::OptimizerStopped(done.stop);
                }
                _ => {}
            }
            rho.clone_from(&selected.rho);
            report.params = params;
            report.last = Some(selected);
            report.stages.push(RobustContinuationStage3 {
                stage, params, incoming_dilated_volume_fraction: start.incoming_volume,
                restoration_scale: start.scale, starting_rho: start.rho,
                gradient_check: check, history, selected_iteration, outcome,
            });
            Ok(())
        })();
        if let Err(error) = attempt {
            report.error = Some(error);
            if report.stop == RobustContinuationStop3::ScheduleComplete {
                report.stop = RobustContinuationStop3::EvaluationStopped;
            }
        }
        if report.stop != RobustContinuationStop3::ScheduleComplete {
            report.stopped_stage = Some(stage);
            break;
        }
    }
    report.work = control.work();
    report
}
