//! Product projection of the actual solid/air residual, not a frozen-solid
//! solve. This path assesses the published field without changing its physical
//! temperature, fluxes, air march, mesh-ladder row or energy balance.

use fs_airflow::conjugate::{AirPath, goal::CoupledGoalError};
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, RobinFeedbackAnalysisConfig};
use fs_conduction::{ConductionProblem, LinearConfig, ThermalInterfaces};
use fs_exec::Cx;
use fs_solver::goal::feedback::FeedbackResidualLimits;

use super::{MaximumEvidence, PropagatedTerm, SolveRefusal, conduction_error, json_string,
    number, optional_number, unavailable};

const SCHEMA: &str = "frankensim.cli.coupled-maximum-evidence.v1";
const METHOD: &str = "outward-coupled-linear-maximum-enclosure";
const SCOPE: &str = "assessment of the published stored linear solid/air system; all downstream air references vary with the walls, while geometry, flow, material coefficients and transport properties are fixed; residual and response-solve rounding are included; no correction, coefficient-assembly, continuum, nonlinear/radiation or physical-validation claim; Estimated";

/// These caps cover the additional coupled analysis, not total allocator peak.
/// Response iterations are shared by every port, never multiplied per branch.
fn policy(memory_bytes: u64, solid: LinearGoalAnalysisConfig, linear: LinearConfig)
    -> RobinFeedbackAnalysisConfig
{
    let count = |divisor| usize::try_from(memory_bytes / divisor).unwrap_or(usize::MAX);
    RobinFeedbackAnalysisConfig {
        residual: FeedbackResidualLimits {
            solid: solid.residual_limits,
            max_ports: 64,
            max_transfer_nonzeros: count(128),
            max_response_entries: count(64),
            max_verification_entries: count(8),
        },
        max_response_iterations: linear.max_iterations,
        max_lowering_entries: count(16),
    }
}

fn cancelled() -> SolveRefusal {
    conduction_error("cli-solve-cancelled", "coupled maximum analysis was cancelled",
        "rerun the solve; no partial coupled evidence was published")
}

#[allow(clippy::too_many_arguments)]
pub(super) fn maximum_evidence(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    paths: &[AirPath], linear: LinearConfig, temperature: &[f64], vertices: &[usize],
    memory_bytes: u64, solid_config: LinearGoalAnalysisConfig, requested_k: f64,
) -> Result<MaximumEvidence, SolveRefusal> {
    cx.checkpoint().map_err(|_| cancelled())?;
    let config = policy(memory_bytes, solid_config, linear);
    let analyzer = match fs_airflow::conjugate::goal::maximum::prepare_linear_maximum(
        cx, problem, interfaces, paths, linear, temperature, solid_config, config,
    ) {
        Ok(analyzer) => analyzer,
        Err(CoupledGoalError::Interrupted
            | CoupledGoalError::Solid(fs_conduction::ConductionError::Cancelled { .. })) => {
            return Err(cancelled());
        }
        Err(error) => {
            // Preparation may have spent response work before refusing. An
            // ordinary input/range/limit refusal is not a zero-work success.
            // It is a typed gap, never a frozen-solid or tolerance fallback.
            cx.checkpoint().map_err(|_| cancelled())?;
            return Ok(unavailable("coupled-analysis-unavailable",
                format!("{SCHEMA}: full solid/air maximum preparation refused: {error}; no frozen-solid bound or tolerance comparison is substituted"), false));
        }
    };
    let analysis = match analyzer.analyze_maximum(cx, temperature, vertices) {
        Ok(analysis) => analysis,
        Err(fs_conduction::ConductionError::Cancelled { .. }) => return Err(cancelled()),
        Err(error) => {
            cx.checkpoint().map_err(|_| cancelled())?;
            return Ok(unavailable("coupled-analysis-unavailable",
                format!("{SCHEMA}: full solid/air maximum evaluation refused: {error}; no frozen-solid bound or tolerance comparison is substituted"), false));
        }
    };
    let coupled = analysis.coupled();
    let requested = (requested_k.is_finite() && requested_k > 0.0).then_some(requested_k);
    let goal_met = requested.is_some_and(|tolerance| analysis.meets_absolute_tolerance(tolerance));
    let bound = analysis.algebraic_half_width_k();
    let interval = match analysis.interval_k() {
        Some([lo, hi]) => format!("[{},{}]", number(lo)?, number(hi)?),
        None => "null".to_string(),
    };
    let status = if bound.is_none() { "coupled-bound-unavailable" }
        else if goal_met { "coupled-goal-tolerance" }
        else { "coupled-goal-unresolved" };
    let detail = format!(
        "{SCHEMA}; published mesh and temperature field; {} independent air paths, {} ordered Robin ports, {} selected vertices ({} free); nominal {:e} K, interval {:?} K; full coupled residual infinity upper {:e}, solid inverse {:?}, feedback gain {:?}, coupled inverse {:?}, inverse route {:?}, port Schur inverse {:?}, disposition {:?}; {} response iterations shared across all ports (cap {}); {} response residuals checked; no new primal solve or field mutation; requested allowance {:?} K, goal met {}; {SCOPE}",
        paths.len(), analyzer.ports().len(), vertices.len(), analysis.free_vertices(),
        analysis.nominal_k(), analysis.interval_k(), coupled.residual_infinity_upper(),
        coupled.solid_inverse_infinity_upper(), coupled.gain_infinity_upper(),
        coupled.coupled_inverse_infinity_upper(), coupled.inverse_method(),
        coupled.schur_inverse_infinity_upper(), coupled.status(), analysis.response_iterations(),
        config.max_response_iterations, coupled.response_residual_infinity_upper().len(),
        requested, goal_met,
    );
    let term = match bound {
        Some(half_width_k) => PropagatedTerm::Measured {
            half_width_k, method: METHOD, detail, vertices: Vec::new(),
        },
        None => PropagatedTerm::Unmeasured { reason: format!(
            "no finite coupled maximum bound was established; a failed sufficient check is not a singularity claim; no frozen-solid or tolerance fallback; {detail}"
        ) },
    };
    let inverse_method = coupled.inverse_method()
        .map(|method| json_string(method.tag())).unwrap_or_else(|| "null".to_string());
    let control_json = format!(
        "{{\"schema\":{},\"status\":{},\"mode\":\"assessment-only\",\"goal_met\":{},\"candidate_accepted\":false,\"requested_tolerance_k\":{},\"initial_bound_k\":{},\"final_bound_k\":{},\"maximum_interval_k\":{},\"primal_iterations\":0,\"response_iterations\":{},\"max_response_iterations\":{},\"max_stability_iterations\":{},\"air_paths\":{},\"ports\":{},\"coupled_residual_infinity_upper\":{},\"solid_inverse_infinity_upper\":{},\"feedback_gain_infinity_upper\":{},\"coupled_inverse_infinity_upper\":{},\"inverse_method\":{},\"schur_inverse_infinity_upper\":{},\"bound_status\":{},\"scope\":{}}}",
        json_string(SCHEMA), json_string(status), goal_met, optional_number(requested)?,
        optional_number(bound)?, optional_number(bound)?, interval,
        analysis.response_iterations(), config.max_response_iterations,
        solid_config.max_stability_iterations, paths.len(), analyzer.ports().len(),
        number(coupled.residual_infinity_upper())?,
        optional_number(coupled.solid_inverse_infinity_upper())?,
        optional_number(coupled.gain_infinity_upper())?,
        optional_number(coupled.coupled_inverse_infinity_upper())?,
        inverse_method,
        optional_number(coupled.schur_inverse_infinity_upper())?,
        json_string(&format!("{:?}", coupled.status())), json_string(SCOPE),
    );
    cx.checkpoint().map_err(|_| cancelled())?;
    Ok(MaximumEvidence { term: Some(term), control_json: Some(control_json), primal_iterations: 0 })
}

#[cfg(test)]
mod tests;
