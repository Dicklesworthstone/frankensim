//! Retarget a fixed-linear rung without resetting its admitted correction work.
//! The caller must re-probe a changed coarse/fine pair before finalizing its
//! measured discretization allowance.

use super::{
    MaximumEvidence, PropagatedTerm, RungSolved, SolveRefusal, analysis_term, conduction_error,
    json_string, number, optional_number, stop_tag,
};
use fs_conduction::LinearMaximumPolish;
use fs_conduction::adjoint::{
    LinearGoalAnalysisConfig, LinearGoalSolveConfig, LinearGoalStop, LinearMaximumAnalysis,
};

const MAX_RETARGET_CALLS: usize = 4;

/// Work belongs to the rung, including attempts whose candidate was rejected.
#[derive(Clone)]
pub(in crate::solve) struct LinearWork {
    initial_bound_k: Option<f64>,
    primal_limit: usize,
    stability_iterations: usize,
    defect_corrections: usize,
    goal_checks: usize,
    retarget_calls: usize,
    last_target_k: Option<f64>,
    last_stop: Option<LinearGoalStop>,
    candidate_accepted: bool,
    physical_gate_refused: bool,
}

impl LinearWork {
    pub(in crate::solve) fn from_analysis(analysis: &LinearMaximumAnalysis, limit: usize) -> Self {
        Self {
            initial_bound_k: analysis.algebraic_half_width_k(),
            primal_limit: limit,
            stability_iterations: analysis.linear_analysis().stability_iterations,
            defect_corrections: 0,
            goal_checks: 0,
            retarget_calls: 0,
            last_target_k: None,
            last_stop: None,
            candidate_accepted: false,
            physical_gate_refused: false,
        }
    }

    pub(in crate::solve) fn from_polish(result: &LinearMaximumPolish, limit: usize) -> Self {
        Self {
            initial_bound_k: result.initial_bound_k,
            defect_corrections: result.defect_corrections,
            goal_checks: result.goal_checks,
            last_target_k: Some(result.requested_tolerance_k),
            last_stop: Some(result.stop),
            candidate_accepted: result.candidate_accepted,
            physical_gate_refused: result.physical_gate_refusal.is_some(),
            ..Self::from_analysis(&result.analysis, limit)
        }
    }
}

pub(in crate::solve) fn maximum_bound(evidence: &MaximumEvidence) -> Option<f64> {
    match &evidence.term {
        Some(PropagatedTerm::Measured { half_width_k, .. })
            if half_width_k.is_finite() && *half_width_k >= 0.0 =>
        {
            Some(*half_width_k)
        }
        _ => None,
    }
}

fn refused(what: impl Into<String>) -> SolveRefusal {
    conduction_error(
        "cli-solve-conduction-adaptive-algebraic",
        what,
        "inspect the measured discretization scale and bounded physical correction",
    )
}

fn poll(cx: &fs_exec::Cx<'_>) -> Result<(), SolveRefusal> {
    cx.checkpoint().map_err(|_| {
        conduction_error(
            "cli-solve-cancelled",
            "adaptive maximum correction cancelled",
            "rerun the solve",
        )
    })
}

pub(in crate::solve) fn retarget_linear_maximum(
    cx: &fs_exec::Cx<'_>,
    solved: &mut RungSolved,
    vertices: &[usize],
    absolute_k: f64,
    memory_bytes: u64,
) -> Result<bool, SolveRefusal> {
    poll(cx)?;
    if !(absolute_k.is_finite() && absolute_k > 0.0) {
        return Err(refused(
            "adaptive correction requires a finite positive absolute target",
        ));
    }
    let Some(mut totals) = solved.algebraic.linear_work.clone() else {
        return Ok(false);
    };
    let remaining = totals
        .primal_limit
        .saturating_sub(solved.algebraic.primal_iterations);
    if remaining == 0 || totals.retarget_calls == MAX_RETARGET_CALLS {
        return Ok(false);
    }
    let data = solved
        .adjoint_data
        .as_ref()
        .ok_or_else(|| refused("adaptive correction lost its retained linear operator"))?;
    let problem = fs_conduction::ConductionProblem {
        mesh: &solved.mesh,
        boundary: &data.boundary,
        material: &data.fallback,
        element_materials: Some(&data.materials),
        source: &data.source,
    };
    let config = LinearGoalAnalysisConfig {
        residual_limits: fs_solver::goal::GoalResidualLimits {
            max_rows: usize::try_from(memory_bytes / 256).unwrap_or(usize::MAX),
            max_nonzeros: usize::try_from(memory_bytes / 64).unwrap_or(usize::MAX),
        },
        // Re-creating this analyzer must retain its stability proposal. The
        // initial call plus four retargets explicitly cap the total at 5*limit.
        max_stability_iterations: totals.primal_limit,
    };
    let result = fs_conduction::polish_linear_maximum(
        cx,
        problem,
        data.interfaces.as_ref(),
        data.linear,
        &solved.solution,
        vertices,
        config,
        LinearGoalSolveConfig {
            absolute_tolerance: absolute_k,
            max_primal_iterations: remaining,
            check_every: 8,
            max_defect_corrections: 2,
        },
    )
    .map_err(|error| match error {
        fs_conduction::ConductionError::Cancelled { .. } => conduction_error(
            "cli-solve-cancelled",
            "adaptive maximum correction cancelled",
            "rerun the solve",
        ),
        other => refused(format!("adaptive physical correction refused: {other}")),
    })?;
    let add = |a: usize, b: usize| {
        a.checked_add(b)
            .ok_or_else(|| refused("adaptive correction work count overflow"))
    };
    let primal_iterations = add(solved.algebraic.primal_iterations, result.primal_iterations)?;
    totals.stability_iterations = add(
        totals.stability_iterations,
        result.analysis.linear_analysis().stability_iterations,
    )?;
    totals.defect_corrections = add(totals.defect_corrections, result.defect_corrections)?;
    totals.goal_checks = add(totals.goal_checks, result.goal_checks)?;
    totals.retarget_calls += 1;
    totals.last_target_k = Some(absolute_k);
    totals.last_stop = Some(result.stop);
    totals.candidate_accepted |= result.candidate_accepted;
    totals.physical_gate_refused = result.physical_gate_refusal.is_some();
    let changed = result.candidate_accepted;
    let old_bound = maximum_bound(&solved.algebraic);
    let new_bound = result.analysis.algebraic_half_width_k();
    let keep_old = !changed && old_bound.is_some_and(|old| new_bound.is_none_or(|new| old <= new));
    poll(cx)?;
    if !keep_old {
        solved.algebraic.term = Some(analysis_term(
            &result.analysis,
            vertices.len(),
            "adaptive retarget with independent physical residual and energy/flux re-evaluation; final allowance is determined only after re-probing the pair",
        ));
    }
    solved.solution = result.solution;
    solved.algebraic.primal_iterations = primal_iterations;
    solved.algebraic.linear_work = Some(totals);
    Ok(changed)
}

/// Publish only the allowance recomputed from the final pair. A previous inner
/// target may differ because corrections changed the measured discrepancy.
pub(in crate::solve) fn finish_discretization_balance(
    solved: &mut RungSolved,
    discretization_k: f64,
    allowance_k: f64,
) -> Result<(), SolveRefusal> {
    if !(discretization_k.is_finite()
        && discretization_k >= 0.0
        && allowance_k.is_finite()
        && allowance_k >= 0.0)
        || allowance_k > fs_math::next_down(discretization_k / 10.0).max(0.0)
    {
        return Err(refused(
            "adaptive algebraic allowance exceeds one tenth of its finite measured term",
        ));
    }
    let Some(totals) = &solved.algebraic.linear_work else {
        return Ok(());
    };
    let bound = maximum_bound(&solved.algebraic);
    let met = bound.is_some_and(|bound| bound <= allowance_k);
    let preparation_limit = totals
        .primal_limit
        .checked_mul(MAX_RETARGET_CALLS + 1)
        .ok_or_else(|| refused("adaptive preparation budget overflow"))?;
    let candidate_stop = totals
        .last_stop
        .map_or_else(|| "null".to_string(), |stop| json_string(stop_tag(stop)));
    solved.algebraic.control_json = Some(format!(
        "{{\"status\":{},\"goal_met\":{},\"correction_supported\":true,\"candidate_accepted\":{},\
         \"tolerance_basis\":\"measured-adaptive-discretization\",\"allocation_fraction\":0.1,\
         \"discretization_half_width_k\":{},\"requested_tolerance_k\":{},\"last_correction_target_k\":{},\
         \"initial_bound_k\":{},\"final_bound_k\":{},\"primal_iterations\":{},\"max_primal_iterations\":{},\
         \"stability_iterations\":{},\"max_stability_iterations\":{},\"retarget_calls\":{},\"max_retarget_calls\":4,\
         \"defect_corrections\":{},\"max_defect_corrections_per_call\":2,\"goal_checks\":{},\
         \"candidate_stop\":{},\"physical_gate_refused\":{},\"scope\":{}}}",
        json_string(if met {
            "discretization-allowance-met"
        } else {
            "discretization-allowance-unresolved"
        }),
        met,
        totals.candidate_accepted,
        number(discretization_k)?,
        number(allowance_k)?,
        optional_number(totals.last_target_k)?,
        optional_number(totals.initial_bound_k)?,
        optional_number(bound)?,
        solved.algebraic.primal_iterations,
        totals.primal_limit,
        totals.stability_iterations,
        preparation_limit,
        totals.retarget_calls,
        totals.defect_corrections,
        totals.goal_checks,
        candidate_stop,
        totals.physical_gate_refused,
        json_string(
            "one tenth of the final re-probed Estimated adaptive discretization term; shared per-rung primal cap and at most four retarget calls; actual stability, defect and driver-check work is cumulative; each changed field independently satisfies the original physical residual gate; no continuum, coupled or nonlinear error bound"
        ),
    ));
    Ok(())
}
