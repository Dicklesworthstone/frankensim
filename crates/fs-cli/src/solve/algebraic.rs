//! Algebraic error of the published region maximum. Only the complete fixed
//! linear solid model enters the outward analyzer; a frozen coupled or
//! nonlinear operator cannot certify the original equations.

use super::{
    BTreeMap, EvidenceWork, PropagatedTerm, QoiRegionTraceError, RungSolved, SolveRefusal,
    canonical_f64, conduction_error, json_string, trace_qoi_region_vertices,
};

mod coupled;

/// Each mesh solve uses a tenth of the requested primary-QoI accuracy for
/// algebraic error. This is an allocation, not a discretization observation.
const ALGEBRAIC_ACCURACY_FRACTION: f64 = 0.1;

#[derive(Default)]
pub(super) struct MaximumEvidence {
    pub(super) term: Option<PropagatedTerm>,
    pub(super) control_json: Option<String>,
    pub(super) primal_iterations: usize,
}

fn unavailable(status: &str, reason: String, estimate_fallback: bool) -> MaximumEvidence {
    MaximumEvidence {
        control_json: Some(format!(
            "{{\"status\":{},\"goal_met\":false,\"requested_tolerance_k\":null,\
             \"primal_iterations\":0,\"reason\":{}}}",
            json_string(status),
            json_string(&reason),
        )),
        term: (!estimate_fallback).then_some(PropagatedTerm::Unmeasured { reason }),
        primal_iterations: 0,
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn maximum_term(
    cx: &fs_exec::Cx<'_>,
    solved: &mut RungSolved,
    region: &str,
    region_ids: &BTreeMap<String, u32>,
    memory_bytes: u64,
    accuracy_rel: f64,
    reference_k: Option<f64>,
    work: EvidenceWork<'_>,
) -> Result<MaximumEvidence, SolveRefusal> {
    let cancelled = || {
        conduction_error(
            "cli-solve-cancelled",
            "the published-maximum algebraic analysis was cancelled",
            "rerun the solve",
        )
    };
    let gap = |reason: String| Ok(unavailable("analysis-unavailable", reason, false));
    let Some(data) = &solved.adjoint_data else {
        return gap("the published rung retained no final-state operator".to_string());
    };
    // The existing whole-model tolerance comparison remains available for
    // radiation. Air feedback is assessed below as the complete coupled model.
    if data.radiating_boundary.is_some() {
        return Ok(unavailable("unsupported-model",
            "maximum-goal corrections do not cover radiation; the whole-model tolerance estimate remains available".to_string(), true));
    }
    for element in 0..solved.labels.len() {
        if element % 1024 == 0 {
            cx.checkpoint().map_err(|_| cancelled())?;
        }
        match data.materials.model_for(element) {
            Ok(model) if model.is_temperature_dependent() => return Ok(unavailable(
                "unsupported-model", "maximum-goal corrections require temperature-independent conductivity; the nonlinear tolerance estimate remains available".to_string(), true)),
            Ok(_) => {}
            Err(error) => return gap(format!("algebraic material binding refused: {error}")),
        }
    }
    let Some(&region_id) = region_ids.get(region) else {
        return gap("the temperature-maximum region has no element label".to_string());
    };
    let (vertices, _) = trace_qoi_region_vertices(
        &solved.labels,
        &solved.mesh.complex().tets,
        solved.mesh.vertex_count(),
        region_id,
        work,
    )
    .map_err(|error| match error {
        QoiRegionTraceError::Cancelled { .. } => cancelled(),
        QoiRegionTraceError::UnitCountOverflow => conduction_error(
            "cli-solve-conduction-algebraic",
            "the algebraic region trace exceeded its representable work count",
            "rerun with a smaller admitted mesh",
        ),
    })?;
    let problem = fs_conduction::ConductionProblem {
        mesh: &solved.mesh,
        boundary: &data.boundary,
        material: &data.fallback,
        element_materials: Some(&data.materials),
        source: &data.source,
    };
    // Explicit structural caps constrain the extra matrix/vector work. They
    // are not a claim about allocator peak memory, which the producer does
    // not measure. The existing Cx checkpoints bound cancellation latency.
    let config = fs_conduction::adjoint::LinearGoalAnalysisConfig {
        residual_limits: fs_solver::goal::GoalResidualLimits {
            max_rows: usize::try_from(memory_bytes / 256).unwrap_or(usize::MAX),
            max_nonzeros: usize::try_from(memory_bytes / 64).unwrap_or(usize::MAX),
        },
        max_stability_iterations: data.linear.max_iterations,
    };
    // Freeze the same physical scale as the adaptive accuracy rule, before
    // corrections. Never use the residual's dimensionless tolerance as K.
    let nominal_k = vertices
        .iter()
        .map(|&v| solved.solution.temperature[v])
        .fold(f64::NEG_INFINITY, f64::max);
    let scale_k = reference_k.map_or(nominal_k.abs(), |reference| (nominal_k - reference).abs());
    let requested_k = ALGEBRAIC_ACCURACY_FRACTION * accuracy_rel * scale_k;
    if !data.air_paths.is_empty() {
        return coupled::maximum_evidence(
            cx, problem, data.interfaces.as_ref(), &data.air_paths, data.linear,
            &solved.solution.temperature, &vertices, memory_bytes, config, requested_k,
        );
    }
    let (analysis, control_json, primal_iterations, control_summary) = if requested_k.is_finite()
        && requested_k > 0.0
    {
        let control = fs_conduction::adjoint::LinearGoalSolveConfig {
            absolute_tolerance: requested_k,
            max_primal_iterations: data.linear.max_iterations,
            check_every: 8,
            max_defect_corrections: 2,
        };
        let polished = match fs_conduction::polish_linear_maximum(
            cx,
            problem,
            data.interfaces.as_ref(),
            data.linear,
            &solved.solution,
            &vertices,
            config,
            control,
        ) {
            Ok(polished) => polished,
            Err(fs_conduction::ConductionError::Cancelled { .. }) => return Err(cancelled()),
            // Normal unresolved goals return typed evidence with actual work.
            // An invalid baseline or arithmetic failure may occur after work
            // starts; do not publish a successful receipt inventing zero work.
            Err(error) => {
                return Err(conduction_error(
                    "cli-solve-conduction-algebraic",
                    format!("maximum-goal correction refused: {error}"),
                    "inspect the reported physical model or numerical refusal before rerunning",
                ));
            }
        };
        let receipt = control_receipt(
            &polished,
            accuracy_rel,
            nominal_k,
            reference_k,
            scale_k,
            control.max_primal_iterations,
        )?;
        let iterations = polished.primal_iterations;
        let summary = format!(
            "requested algebraic allowance {} K from 10 percent of requested accuracy; \
             accepted-field goal met {}, correction stop {:?}, {} correction PCG iterations, \
             candidate accepted {}, physical gate refused {}",
            requested_k,
            polished.goal_met,
            polished.stop,
            iterations,
            polished.candidate_accepted,
            polished.physical_gate_refusal.is_some(),
        );
        // This field has independently passed the original physical residual
        // gate, and all energy/contact/Robin outputs have been recomputed.
        solved.solution = polished.solution;
        (polished.analysis, receipt, iterations, summary)
    } else {
        let analysis = match fs_conduction::adjoint::analyze_linear_maximum(
            cx,
            problem,
            data.interfaces.as_ref(),
            data.linear,
            &solved.solution.temperature,
            &vertices,
            config,
        ) {
            Ok(analysis) => analysis,
            Err(fs_conduction::ConductionError::Cancelled { .. }) => return Err(cancelled()),
            Err(error) => {
                return gap(format!(
                    "published-maximum algebraic analysis refused: {error}"
                ));
            }
        };
        // Zero rise or an unrepresentable allocation admits no positive K
        // target. Preserve the measured bound and disclose the missing goal.
        (
            analysis,
            format!(
                "{{\"status\":\"unallocated-accuracy-scale\",\"goal_met\":false,\
             \"requested_tolerance_k\":null,\"primal_iterations\":0,\"reason\":{}}}",
                json_string(
                    "the requested temperature-rise accuracy has no finite positive algebraic allocation; no absolute-temperature floor is invented"
                ),
            ),
            0,
            "no finite positive algebraic accuracy allocation; no goal completion claimed"
                .to_string(),
        )
    };
    Ok(MaximumEvidence {
        term: Some(analysis_term(&analysis, vertices.len(), &control_summary)),
        control_json: Some(control_json),
        primal_iterations,
    })
}

fn analysis_term(
    analysis: &fs_conduction::adjoint::LinearMaximumAnalysis,
    vertices: usize,
    control_summary: &str,
) -> PropagatedTerm {
    let linear_analysis = analysis.linear_analysis();
    let enclosure = &linear_analysis.enclosure;
    let Some(half_width_k) = analysis.algebraic_half_width_k() else {
        return PropagatedTerm::Unmeasured {
            reason: format!(
                "the stored linear operator has no finite verified maximum-error bound; \
             inverse status {:?}, primal residual upper {:e}; no tolerance \
             comparison is substituted for missing inverse evidence; {control_summary}",
                enclosure.status(),
                enclosure.primal_residual_infinity_upper(),
            ),
        };
    };
    PropagatedTerm::Measured {
        half_width_k,
        method: "outward-linear-maximum-enclosure",
        detail: format!(
            "published mesh and temperature field; fixed linear conductivity, Robin \
             references and matching contact; {} region vertices ({} free), nominal \
             maximum {:e} K, enclosed maximum {:?} K; outward residual infinity \
             upper {:e}, verified inverse infinity upper {:?}, stability proposal \
             relative residual {:?}, stability iterations {}; \
             ||A^-1||_infinity * ||b-A*T||_infinity bounds the full region maximum, \
             including a change of hottest node and exact prescribed values; \
             stored floating-point system only, not assembly, discretization or physical \
             uncertainty; residual evaluation rounding is already included in this term; \
             {control_summary}; Estimated",
            vertices,
            analysis.free_vertices(),
            analysis.nominal_k(),
            analysis.interval_k(),
            enclosure.primal_residual_infinity_upper(),
            enclosure.inverse_infinity_upper(),
            linear_analysis.stability_relative_residual,
            linear_analysis.stability_iterations,
        ),
        vertices: Vec::new(),
    }
}

fn number(value: f64) -> Result<String, SolveRefusal> {
    canonical_f64(value).ok_or_else(|| {
        conduction_error(
            "cli-solve-conduction-algebraic-nonfinite",
            "a maximum-goal control value is not finite",
            "inspect the declared accuracy and physical temperature scale",
        )
    })
}

fn optional_number(value: Option<f64>) -> Result<String, SolveRefusal> {
    value
        .map(number)
        .transpose()
        .map(|value| value.unwrap_or_else(|| "null".to_string()))
}

fn stop_tag(stop: fs_conduction::adjoint::LinearGoalStop) -> &'static str {
    use fs_conduction::adjoint::LinearGoalStop;
    match stop {
        LinearGoalStop::GoalTolerance => "goal-tolerance",
        LinearGoalStop::IterationBudget => "iteration-budget",
        LinearGoalStop::BoundUnavailable => "bound-unavailable",
        LinearGoalStop::NoProgress => "no-progress",
        LinearGoalStop::DefectCorrectionBudget => "defect-correction-budget",
    }
}

fn control_receipt(
    result: &fs_conduction::LinearMaximumPolish,
    accuracy_rel: f64,
    initial_maximum_k: f64,
    reference_k: Option<f64>,
    scale_k: f64,
    max_primal_iterations: usize,
) -> Result<String, SolveRefusal> {
    let status = if result.goal_met {
        "goal-tolerance"
    } else if result.physical_gate_refusal.is_some() {
        "physical-residual-gate"
    } else {
        stop_tag(result.stop)
    };
    let physical_refusal = match &result.physical_gate_refusal {
        Some(refusal) => format!(
            "{{\"candidate_residual_w\":{},\"residual_threshold_w\":{}}}",
            number(refusal.candidate_residual_w)?,
            number(refusal.residual_threshold_w)?,
        ),
        None => "null".to_string(),
    };
    Ok(format!(
        "{{\"status\":{},\"goal_met\":{},\"candidate_accepted\":{},\
         \"accuracy_rel\":{},\"allocation_fraction\":{},\"tolerance_basis\":{},\
         \"initial_maximum_k\":{},\"reference_k\":{},\"scale_k\":{},\
         \"requested_tolerance_k\":{},\"initial_bound_k\":{},\"final_bound_k\":{},\
         \"primal_iterations\":{},\"max_primal_iterations\":{},\
         \"defect_corrections\":{},\"max_defect_corrections\":2,\"goal_checks\":{},\
         \"candidate_stop\":{},\"physical_gate_refusal\":{},\"scope\":{}}}",
        json_string(status),
        result.goal_met,
        result.candidate_accepted,
        number(accuracy_rel)?,
        number(ALGEBRAIC_ACCURACY_FRACTION)?,
        json_string(if reference_k.is_some() {
            "initial-temperature-rise"
        } else {
            "initial-absolute-temperature"
        }),
        number(initial_maximum_k)?,
        optional_number(reference_k)?,
        number(scale_k)?,
        number(result.requested_tolerance_k)?,
        optional_number(result.initial_bound_k)?,
        optional_number(result.analysis.algebraic_half_width_k())?,
        result.primal_iterations,
        max_primal_iterations,
        result.defect_corrections,
        result.goal_checks,
        json_string(stop_tag(result.stop)),
        physical_refusal,
        json_string(
            "10 percent of the requested QoI accuracy on the frozen pre-correction temperature scale; not 10 percent of a measured discretization error; bounded corrections after the physical solve, with independent residual and energy/flux re-evaluation; stored linear system only"
        ),
    ))
}
