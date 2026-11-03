//! Reuse the source-backed contact step without advancing trial history.

use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContacts};

use super::*;

pub(super) mod enthalpy;

pub(super) const MAX_RESPONSES: usize = 64;
pub(super) const NO_CLAIM: &str = "Source-backed matching contact queries each pinned resistance law at the arithmetic mean of its two independent face-temperature means. Pressure, finish, medium, geometry and claim selection remain fixed. Every counted contact trial uses the same immutable old field; the accepted unchanged endpoint passes freshly queried complete nodal joule-residual and separate storage/heat gates. Published contact heat and reservoir reactions describe that physical endpoint, while inner Krylov residuals describe its counted numerical response. Contact heat is internal transfer and is not added to external input. Scalar and continuous linear source curves share the complete capacity/bulk/contact Newton solve, with generalized primal source tangents at curve knots and inward source endpoints; every physical trial still passes its pinned source query. Surface laws retain the bounded contact fixed point. Neither algorithm establishes arbitrary-law convergence or uniqueness. This capacity path does not supply total-enthalpy storage; declared total-enthalpy trajectories use the separate shared contact/enthalpy owner with complete Newton tangents for scalar/continuous-curve sources and a bounded physical fixed point for source surfaces. Radiation, derived natural/airflow convection and trajectory derivatives are not supplied for temperature-varying contact by this capacity adapter.";

pub(super) const fn method(coupled_newton: bool) -> &'static str {
    if coupled_newton {
        "source-backed-contact-newton"
    } else {
        "source-backed-contact-fixed-point"
    }
}

pub(super) const fn max_responses(coupled_newton: bool) -> usize {
    if coupled_newton { 1 } else { MAX_RESPONSES }
}

pub(super) struct Step {
    pub(super) conduction: StepSolution,
    pub(super) receipt: String,
    pub(super) responses: usize,
    pub(super) krylov_iterations: usize,
    pub(super) nonlinear_updates: usize,
    pub(super) nonlinear_backtracks: usize,
}

/// The full physical residual is in joules, independent of the temperature
/// offset or an inner relative Krylov residual. The energy audit remains a
/// separate gate; the factor leaves room for summed and rounded residuals.
pub(super) fn residual_tolerance_j(energy_tolerance_j: f64, vertices: usize) -> Result<f64, SolveRefusal> {
    let tolerance = 0.1 * energy_tolerance_j
        / fs_math::det::sqrt(vertices.max(1) as f64);
    if tolerance.is_finite() && tolerance > 0.0 {
        Ok(tolerance)
    } else {
        Err(bad("temperature-contact transient residual allowance is not positive and finite"))
    }
}

pub(super) fn face_states(
    cx: &Cx<'_>,
    contacts: &fs_conduction::interface::temperature::TemperatureContactEvaluation,
    deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    let mut states = Vec::with_capacity(contacts.states().len());
    for state in contacts.states() {
        poll(cx, deadline)?;
        states.push(format!(concat!(
            "{{\"interface\":{},\"side_a_face\":{},\"side_b_face\":{},",
            "\"side_a_temperature_k\":{},\"side_b_temperature_k\":{},\"query_temperature_k\":{},",
            "\"resistance_m2_k_w\":{},\"interface_card\":{},\"selected_claim\":{},\"usage_receipt\":{}}}"
        ),
            json_string(&state.interface), state.pair.side_a, state.pair.side_b,
            number(state.side_a_temperature_k)?, number(state.side_b_temperature_k)?,
            number(state.interface_temperature_k)?, number(state.resistance.value_m2_k_per_w())?,
            json_string(&state.resistance.card_identity().to_hex()),
            json_string(&state.resistance.receipt().selected.0.to_hex()),
            json_string(&state.resistance.receipt().content_hash().to_hex()),
        ));
    }
    Ok(states.join(","))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, coupled_newton: bool, max_face_pairs: usize,
    old: &[f64], dt_s: f64,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    let tolerance = residual_tolerance_j(config.energy_tolerance_j, problem.mesh.vertex_count())?;
    let coupling = ContactIterationConfig {
        max_iterations: max_responses(coupled_newton),
        relaxation: 0.5,
        residual_tolerance_j: tolerance,
    };
    let accepted = if coupled_newton {
        let nonlinear = nonlinear.ok_or_else(|| bad("coupled contact Newton lost its nonlinear controls"))?;
        contacts.advance_capacity_newton(
            cx, engine, problem, old, dt_s, config, nonlinear, coupling, max_face_pairs,
        )
    } else {
        contacts.advance_capacity(cx, engine, problem, old, dt_s, config, nonlinear, coupling)
    }.map_err(|error| match error {
        cancelled @ ConductionError::Cancelled { .. } => lower(cancelled),
        other => bad(format!("source-backed contact step refused: {other}")),
    })?;
    poll(cx, deadline)?;

    let states = face_states(cx, &accepted.physical_contacts, deadline)?;
    let receipt = format!(concat!(
        "{{\"method\":{},\"thermal_responses\":{},",
        "\"krylov_iterations\":{},\"nonlinear_updates\":{},\"nonlinear_backtracks\":{},",
        "\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},",
        "\"physical_energy_residual_j\":{},\"energy_tolerance_j\":{},\"face_states\":[{}]}}"
    ),
        json_string(method(coupled_newton)), accepted.iterations,
        accepted.krylov_iterations, accepted.nonlinear_updates,
        accepted.nonlinear_backtracks, number(accepted.physical.residual_norm_j)?,
        number(tolerance)?, number(accepted.physical.energy_residual_j)?,
        number(config.energy_tolerance_j)?, states,
    );
    let mut conduction = accepted.conduction;
    // The counted inner solve and physical audit describe the same unchanged
    // field. Publish freshly queried contact heat and physical reservoirs;
    // its inner Krylov residual remains explicitly numerical solve evidence.
    conduction.contact_fluxes = accepted.physical_contacts.frozen_interfaces()
        .fluxes(&conduction.temperature).map_err(lower)?;
    conduction.robin_fluxes = accepted.physical.convective_robin_fluxes;
    conduction.source_w = accepted.physical.source_w;
    conduction.neumann_out_w = accepted.physical.neumann_out_w;
    conduction.robin_out_w = accepted.physical.convective_out_w;
    conduction.dirichlet_in_w = accepted.physical.dirichlet_in_w;
    conduction.stored_energy_change_j = accepted.physical.stored_energy_change_j;
    conduction.energy_residual_j = accepted.physical.energy_residual_j;
    poll(cx, deadline)?;
    Ok(Step {
        conduction, receipt, responses: accepted.iterations,
        krylov_iterations: accepted.krylov_iterations,
        nonlinear_updates: accepted.nonlinear_updates,
        nonlinear_backtracks: accepted.nonlinear_backtracks,
    })
}
