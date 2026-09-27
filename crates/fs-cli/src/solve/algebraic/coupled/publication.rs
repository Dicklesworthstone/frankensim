//! One publication boundary for the corrected field, boundary and air receipt.
//! No caller-visible state changes until physical admission AND serialization
//! finish. Rejected candidates contribute work, never their bound to an old field.

use fs_airflow::conjugate::{ConjugateConfig, goal::maximum::{
    SpectralMaximumControl,
    physical::{PhysicalAirMaximumPolish, PhysicalCoolingGates, polish_linear_maximum_with_spectral},
}};
use fs_conduction::{ConductionSolution, ThermalBoundary};
use fs_conduction::adjoint::LinearGoalSolveConfig;
use crate::solve::RungSolved;
use super::{AirPath, ConductionProblem, CoupledGoalError, Cx, LinearConfig,
    LinearGoalAnalysisConfig, MaximumEvidence, PropagatedTerm, SolveRefusal,
    ThermalInterfaces, METHOD, cancelled, conduction_error, json_string,
    maximum_evidence, number, optional_number, policy, spectral_policy};

mod receipt;

const SCHEMA: &str = "frankensim.cli.coupled-maximum-publication.v1";
const SCOPE: &str = "full stored linear solid/air correction with independent production-law physical admission; field, Robin boundary and air receipt are published together; assembly, continuum, nonlinear/radiation, flow uncertainty and experimental validation remain outside this bound";

struct Replacement {
    solution: ConductionSolution,
    boundary: ThermalBoundary,
    conjugate: String,
}

struct Publication {
    evidence: MaximumEvidence,
    replacement: Option<Replacement>,
}

fn error(what: impl Into<String>) -> SolveRefusal {
    conduction_error("cli-solve-cooling-publication", what,
        "inspect the retained model, numerical budget and physical cooling gates; no partial replacement was published")
}

/// Called before the rung feeds mesh comparison, roundoff, QoI or export.
/// Every fallible operation, including receipt reconstruction, precedes mutation.
pub(in crate::solve::algebraic) fn polish_rung(
    cx: &Cx<'_>, solved: &mut RungSolved, vertices: &[usize], memory_bytes: u64,
    solid_config: LinearGoalAnalysisConfig, requested_k: f64,
) -> Result<MaximumEvidence, SolveRefusal> {
    cx.checkpoint().map_err(|_| cancelled())?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| error("missing retained operator"))?;
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &data.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    if !(requested_k.is_finite() && requested_k > 0.0) {
        return maximum_evidence(cx, problem, data.interfaces.as_ref(), &data.air_paths,
            data.linear, &solved.solution.temperature, vertices, memory_bytes, solid_config, requested_k);
    }
    let original_fragment = solved.conjugate_fragment.as_deref()
        .ok_or_else(|| error("missing original air-exchange receipt"))?;
    let publication = prepare(cx, problem, data.interfaces.as_ref(), &data.air_paths,
        data.linear, &solved.solution, original_fragment, vertices, memory_bytes,
        solid_config, requested_k)?;
    cx.checkpoint().map_err(|_| cancelled())?;
    // This mutable borrow is acquired only after all work has succeeded. Moving
    // the three owned values has no callback, allocation or cancellation seam.
    let data = solved.adjoint_data.as_mut().ok_or_else(|| error("retained operator disappeared"))?;
    let Publication { evidence, replacement } = publication;
    if let Some(replacement) = replacement {
        solved.solution = replacement.solution;
        data.boundary = replacement.boundary;
        solved.conjugate_fragment = Some(replacement.conjugate);
    }
    Ok(evidence)
}

#[allow(clippy::too_many_arguments)]
fn prepare(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    paths: &[AirPath], linear: LinearConfig, original: &ConductionSolution,
    original_fragment: &str, vertices: &[usize], memory_bytes: u64,
    solid_config: LinearGoalAnalysisConfig, requested_k: f64,
) -> Result<Publication, SolveRefusal> {
    let feedback = policy(memory_bytes, solid_config, linear);
    let spectral_limits = spectral_policy(memory_bytes, feedback);
    let exchange = ConjugateConfig::default();
    let gates = PhysicalCoolingGates {
        // Additional physical acceptance, NOT inferred from the temperature goal.
        energy_relative_tolerance: 1e-6,
        reference_tolerance_k: exchange.temperature_tolerance_k,
        balance_tolerance_w: exchange.balance_tolerance_w,
        balance_relative_tolerance: exchange.balance_relative_tolerance,
    };
    let control = LinearGoalSolveConfig { absolute_tolerance: requested_k,
        max_primal_iterations: linear.max_iterations, check_every: 8, max_defect_corrections: 2 };
    // A dimensioned proposal only: integrated Robin conductance per vertex.
    // Neither this heuristic nor a factorization success grants inverse authority.
    let mut conductance = f64::INFINITY;
    for flux in &original.report.robin_fluxes {
        cx.checkpoint().map_err(|_| cancelled())?;
        let value = flux.mean_htc_w_per_m2_k * flux.area_m2;
        if value.is_finite() && value > 0.0 { conductance = conductance.min(value); }
    }
    let initial_shift = conductance / problem.mesh.vertex_count().max(1) as f64;
    if !(initial_shift.is_finite() && initial_shift > 0.0) {
        return Err(error("no finite positive Robin conductance for a spectral shift proposal"));
    }
    let result = polish_linear_maximum_with_spectral(cx, problem, interfaces, paths,
        linear, original, vertices, solid_config, feedback, control,
        SpectralMaximumControl { initial_shift, limits: spectral_limits }, gates)
        .map_err(|failure| match failure {
            CoupledGoalError::Interrupted
            | CoupledGoalError::Solid(fs_conduction::ConductionError::Cancelled { .. })
            | CoupledGoalError::Air(fs_airflow::AirflowError::Cancelled { .. }) => cancelled(),
            other => error(format!("coupled physical correction refused: {other}")),
        })?;
    let mut publication = project(cx, result, original_fragment, paths, requested_k, linear.max_iterations,
        gates, spectral_limits)?;
    let mut references = Vec::new();
    for path in paths {
        for segment in path.segments() {
            cx.checkpoint().map_err(|_| cancelled())?;
            let index = problem.boundary.region_names().iter().position(|name| name == segment.region())
                .ok_or_else(|| error("missing original coupling reference"))?;
            let value = match &problem.boundary.conditions()[index] {
                fs_conduction::ThermalBc::Robin { t_ref: fs_conduction::ScalarField::Uniform(value), .. } => *value,
                _ => return Err(error("original coupling reference is not uniform Robin")),
            };
            references.push(format!("{{\"region\":{},\"reference_k\":{}}}",
                json_string(segment.region()), number(value)?));
        }
    }
    if let Some(control) = &mut publication.evidence.control_json {
        control.pop();
        control.push_str(&format!(",\"analysis_reference_origin\":[{}]}}", references.join(",")));
    }
    cx.checkpoint().map_err(|_| cancelled())?;
    Ok(publication)
}

#[allow(clippy::too_many_arguments)]
fn project(
    cx: &Cx<'_>, mut result: PhysicalAirMaximumPolish, original_fragment: &str,
    paths: &[AirPath], requested_k: f64, primal_limit: usize, gates: PhysicalCoolingGates,
    spectral_limits: fs_solver::goal::inverse::spectral::SpectralInverseLimits,
) -> Result<Publication, SolveRefusal> {
    let correction = &result.correction.solution.solid;
    if let Some(accepted) = &result.accepted {
        if accepted.solid.temperature.len() != correction.temperature.len()
            || accepted.solid.temperature.iter().zip(&correction.temperature)
                .any(|(a, b)| a.to_bits() != b.to_bits())
        { return Err(error("accepted physical field and numerical goal field differ")); }
    }
    let analysis = if result.accepted.is_some() { &correction.analysis }
        else { &result.correction.initial_analysis };
    let bound = analysis.algebraic_half_width_k();
    let initial_bound = result.correction.initial_analysis.algebraic_half_width_k();
    let goal_met = analysis.meets_absolute_tolerance(requested_k);
    let changed = result.accepted.as_ref().is_some_and(|accepted| accepted.temperature_changed);
    let physical_accepted = result.accepted.is_some();
    let refusal = result.physical_refusal.as_ref().map(|value| json_string(&format!("{value:?}")))
        .unwrap_or_else(|| "null".into());
    let interval = match analysis.interval_k() {
        Some([lo, hi]) => format!("[{},{}]", number(lo)?, number(hi)?),
        None => "null".into(),
    };
    let status = if !physical_accepted { "physical-gate-refused" }
        else if bound.is_none() { "coupled-bound-unavailable" }
        else if goal_met { "coupled-goal-tolerance" } else { "coupled-goal-unresolved" };
    let preparation = match &result.correction.preparation {
        Some(p) => format!(
            "{{\"stop\":{},\"work_entries\":{},\"peak_storage_entries\":{},\"shift_attempts\":{},\"shift\":{},\"defect_upper\":{},\"coercivity_lower\":{},\"max_work_entries\":{},\"max_storage_entries\":{}}}",
            json_string(&format!("{:?}", p.stop)), p.work_entries, p.peak_storage_entries,
            p.shift_attempts, optional_number(p.shift)?, optional_number(p.defect_upper)?,
            optional_number(p.coercivity_lower)?, spectral_limits.max_work_entries,
            spectral_limits.max_storage_entries),
        None => "null".into(),
    };
    let control_json = format!(
        "{{\"schema\":{},\"mode\":\"physical-goal-correction\",\"status\":{},\"correction_supported\":true,\"physical_accepted\":{},\"candidate_accepted\":{},\"goal_met\":{},\"requested_tolerance_k\":{},\"initial_bound_k\":{},\"candidate_bound_k\":{},\"final_bound_k\":{},\"maximum_interval_k\":{},\"primal_iterations\":{},\"max_primal_iterations\":{},\"defect_corrections\":{},\"goal_checks\":{},\"preliminary_goal_checks\":1,\"response_iterations\":{},\"physical_checks\":{},\"physical_rejections\":{},\"candidate_stop\":{},\"physical_refusal\":{},\"spectral_preparation\":{},\"physical_gates\":{{\"energy_relative\":{},\"reference_k\":{},\"balance_absolute_w\":{},\"balance_relative\":{}}},\"scope\":{}}}",
        json_string(SCHEMA), json_string(status), physical_accepted, changed, goal_met,
        number(requested_k)?, optional_number(initial_bound)?,
        optional_number(correction.analysis.algebraic_half_width_k())?, optional_number(bound)?,
        interval, correction.primal_iterations, primal_limit, correction.defect_corrections,
        correction.goal_checks, correction.analysis.response_iterations(), result.physical_checks,
        result.physical_rejections, json_string(super::super::stop_tag(correction.stop)), refusal,
        preparation, number(gates.energy_relative_tolerance)?, number(gates.reference_tolerance_k)?,
        number(gates.balance_tolerance_w)?, number(gates.balance_relative_tolerance)?, json_string(SCOPE));
    let term = match bound {
        Some(half_width_k) => PropagatedTerm::Measured { half_width_k, method: METHOD,
            detail: format!("{SCHEMA}; published field: {}; nominal {} K; maximum interval {:?}; full coupled residual upper {}; physical checks {}; {SCOPE}",
                if physical_accepted { "accepted correction" } else { "unchanged original" },
                analysis.nominal_k(), analysis.interval_k(), analysis.coupled().residual_infinity_upper(),
                result.physical_checks), vertices: Vec::new() },
        None => PropagatedTerm::Unmeasured {
            reason: "the published field has no finite coupled maximum enclosure; candidate-only bounds are not applied to an unchanged field".into(),
        },
    };
    let primal_iterations = correction.primal_iterations;
    // Preserve the original coupling origin in the replay evidence. Replacing
    // a physical Robin reference changes the rounded RHS, not the modeled air law.
    let mut control_json = control_json;
    control_json.pop();
    control_json.push_str(&format!(
        ",\"bound_target\":\"pre-correction-stored-affine-system\",\"air_paths\":{},\"ports\":{},\"max_response_iterations\":{},\"coupled_residual_infinity_upper\":{},\"solid_inverse_infinity_upper\":{},\"coupled_inverse_infinity_upper\":{},\"bound_status\":{}}}",
        paths.len(), paths.iter().map(|path| path.segments().len()).sum::<usize>(),
        primal_limit, number(analysis.coupled().residual_infinity_upper())?,
        optional_number(analysis.coupled().solid_inverse_infinity_upper())?,
        optional_number(analysis.coupled().coupled_inverse_infinity_upper())?,
        json_string(&format!("{:?}", analysis.coupled().status()))));
    let evidence = MaximumEvidence { term: Some(term), control_json: Some(control_json),
        primal_iterations, linear_work: None };
    let replacement = match result.accepted.take() {
        Some(accepted) => {
            let conjugate = receipt::rebuild(cx, original_fragment, paths, &accepted, primal_iterations)?;
            Some(Replacement { solution: accepted.solid, boundary: accepted.boundary, conjugate })
        }
        None => None,
    };
    cx.checkpoint().map_err(|_| cancelled())?;
    Ok(Publication { evidence, replacement })
}

#[cfg(test)]
mod tests;
