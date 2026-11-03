//! Reuse the source-backed contact step without advancing trial history.

use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContacts};

use super::*;

pub(super) const MAX_RESPONSES: usize = 64;
pub(super) const NO_CLAIM: &str = "Source-backed matching contact queries each pinned resistance law at the arithmetic mean of its two independent face-temperature means. Pressure, finish, medium, geometry and claim selection remain fixed. Every counted contact trial uses the same immutable old field; the accepted unchanged endpoint passes freshly queried complete nodal joule-residual and separate storage/heat gates. Published contact heat and reservoir reactions describe that physical endpoint, while inner Krylov residuals describe its counted numerical response. Contact heat is internal transfer and is not added to external input. This bounded fixed point does not establish arbitrary-law convergence or uniqueness. Total-enthalpy storage, radiation, derived natural/airflow convection and trajectory derivatives are not supplied for temperature-varying contact by this native adapter.";

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

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, old: &[f64], dt_s: f64,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    let tolerance = residual_tolerance_j(config.energy_tolerance_j, problem.mesh.vertex_count())?;
    let accepted = contacts.advance_capacity(
        cx, engine, problem, old, dt_s, config, nonlinear,
        ContactIterationConfig {
            max_iterations: MAX_RESPONSES,
            relaxation: 0.5,
            residual_tolerance_j: tolerance,
        },
    ).map_err(lower)?;
    poll(cx, deadline)?;

    let mut states = Vec::with_capacity(accepted.physical_contacts.states().len());
    for state in accepted.physical_contacts.states() {
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
    let receipt = format!(concat!(
        "{{\"method\":\"source-backed-contact-fixed-point\",\"thermal_responses\":{},",
        "\"krylov_iterations\":{},\"nonlinear_updates\":{},\"nonlinear_backtracks\":{},",
        "\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},",
        "\"physical_energy_residual_j\":{},\"energy_tolerance_j\":{},\"face_states\":[{}]}}"
    ),
        accepted.iterations, accepted.krylov_iterations, accepted.nonlinear_updates,
        accepted.nonlinear_backtracks, number(accepted.physical.residual_norm_j)?,
        number(tolerance)?, number(accepted.physical.energy_residual_j)?,
        number(config.energy_tolerance_j)?, states.join(","),
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
