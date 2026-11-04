//! Reuse the source-backed contact step without advancing trial history.

use fs_conduction::interface::temperature::{
    ContactIterationConfig, TemperatureContactEvaluation, TemperatureContacts,
};
use fs_conduction::transient::backward_euler::EndpointEvaluation;

use super::*;

pub(super) mod enthalpy;

pub(super) const MAX_RESPONSES: usize = 64;
pub(super) const NO_CLAIM: &str = "Source-backed matching contact queries each pinned resistance law at the arithmetic mean of its two independent face-temperature means. Pressure, finish, medium, geometry and claim selection remain fixed. Every counted contact trial uses the same immutable old field; the accepted unchanged endpoint passes freshly queried complete nodal joule-residual and separate storage/heat gates. Published contact heat and reservoir reactions describe that physical endpoint, while inner Krylov residuals describe its counted numerical response. Contact heat is internal transfer and is not added to external input. Scalar and continuous linear source curves share the complete capacity/bulk/contact Newton solve, with generalized primal source tangents at curve knots and inward source endpoints; every physical trial still passes its pinned source query. Surface laws retain the bounded contact fixed point. Neither algorithm establishes arbitrary-law convergence or uniqueness. This capacity path does not supply total-enthalpy storage; declared total-enthalpy trajectories use the separate shared contact/enthalpy owner with complete Newton tangents for scalar/continuous-curve sources and a bounded physical fixed point for source surfaces. Natural-convection trials use these same live contact responses and the original old temperature. Acceptance re-queries contact and evaluates the complete unchanged endpoint under its updated physical convection coefficient, with all nested responses counted. Transported-air callbacks use these same contact-aware capacity responses. Natural convection can enclose complete air solves on distinct targets, followed by a fresh combined contact/natural/air audit. Area-mean ambient gray or closed-enclosure radiation encloses these same contact responses inside each prescribed or derived-boundary trial; predictor and accepted physical audits include the actual radiative law. Only forced convection heats air; ambient and off-path reservoir heat retain separate accounts. Contact and closed reflected exchange remain internal, while authored Neumann heat and prescribed-node storage reactions remain external. Radiation callbacks and all nested contact responses, Newton updates, backtracks and Krylov iterations are counted separately. Optional final/sampled-peak capacity adjoints retain the complete smooth source-contact tangent for prescribed or scheduled thermal boundaries, with no replacement forward solve. Open-enclosure/pointwise radiation and contact-history derivatives with derived convection, radiation or material controls remain unsupported by this capacity adapter.";

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

#[derive(Clone, Copy)]
pub(super) struct Config {
    pub(super) coupled_newton: bool,
    pub(super) max_face_pairs: usize,
    pub(super) residual_tolerance_j: f64,
}

impl Config {
    pub(super) fn max_responses(self) -> usize {
        max_responses(self.coupled_newton)
    }
}

pub(super) struct Step {
    pub(super) conduction: StepSolution,
    pub(super) receipt: String,
    pub(super) responses: usize,
    pub(super) krylov_iterations: usize,
    pub(super) nonlinear_updates: usize,
    pub(super) nonlinear_backtracks: usize,
}

pub(super) fn admit_radiation(
    lowering: Option<&radiation::LoweredRadiation>,
    model: Option<&fs_conduction::radiation::GrayDiffuseEnclosure>,
) -> Result<(), SolveRefusal> {
    if lowering.is_some_and(|row| row.is_pointwise() || row.open_ambient_temperature_k().is_some())
    {
        return Err(bad("live capacity contact requires area-mean ambient or closed-enclosure radiation; pointwise and open-enclosure radiation are unsupported"));
    }
    if lowering.is_some_and(radiation::LoweredRadiation::is_enclosure) != model.is_some() {
        return Err(bad("live capacity contact lost its declared closed-enclosure radiation model"));
    }
    Ok(())
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

/// Re-query the pinned laws at the unchanged candidate, then audit that field
/// under the caller's actual boundary. This never advances thermal history.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, old: &[f64], candidate: &[f64],
    dt_s: f64, config: StepConfig,
    lowering: Option<&radiation::LoweredRadiation>,
    model: Option<&fs_conduction::radiation::GrayDiffuseEnclosure>,
    deadline: Option<(Instant, f64)>,
) -> Result<(EndpointEvaluation, TemperatureContactEvaluation), SolveRefusal> {
    poll(cx, deadline)?;
    admit_radiation(lowering, model)?;
    if old.len() != problem.mesh.vertex_count() || candidate.len() != old.len() {
        return Err(bad("contact capacity audit requires complete original and candidate fields"));
    }
    let actual = contacts.evaluate(cx, problem.boundary, candidate).map_err(lower)?;
    poll(cx, deadline)?;
    let physical = radiative::evaluate(
        cx, engine, problem, Some(actual.frozen_interfaces()), old, candidate, dt_s,
        config, lowering, model, deadline,
    )?;
    poll(cx, deadline)?;
    Ok((physical, actual))
}

/// Publish physical heat at the unchanged counted field. Original convection
/// stays separate in its trace rows; total outward Robin heat includes radiation.
pub(super) fn publish_physical(
    conduction: &mut StepSolution,
    physical: &EndpointEvaluation,
    contacts: &TemperatureContactEvaluation,
) -> Result<(), SolveRefusal> {
    conduction.contact_fluxes = contacts.frozen_interfaces()
        .fluxes(&conduction.temperature).map_err(lower)?;
    conduction.robin_fluxes = physical.convective_robin_fluxes.clone();
    conduction.source_w = physical.source_w;
    conduction.neumann_out_w = physical.neumann_out_w;
    conduction.robin_out_w = physical.convective_out_w + physical.radiation_out_w;
    conduction.dirichlet_in_w = physical.dirichlet_in_w;
    conduction.stored_energy_change_j = physical.stored_energy_change_j;
    conduction.energy_residual_j = physical.energy_residual_j;
    Ok(())
}

/// The accepted boundary owner supplies its fresh physical audit and all
/// counted inner responses; numerical Krylov evidence retains its own meaning.
#[allow(clippy::too_many_arguments)]
pub(super) fn step_receipt(
    cx: &Cx<'_>, controls: Config, physical: &EndpointEvaluation,
    contacts: &TemperatureContactEvaluation, responses: usize, krylov_iterations: usize,
    nonlinear_updates: usize, nonlinear_backtracks: usize, energy_tolerance_j: f64,
    deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    poll(cx, deadline)?;
    let states = face_states(cx, contacts, deadline)?;
    let receipt = format!(concat!(
        "{{\"method\":{},\"thermal_responses\":{},",
        "\"krylov_iterations\":{},\"nonlinear_updates\":{},\"nonlinear_backtracks\":{},",
        "\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},",
        "\"physical_energy_residual_j\":{},\"energy_tolerance_j\":{},\"face_states\":[{}]}}"
    ),
        json_string(method(controls.coupled_newton)), responses,
        krylov_iterations, nonlinear_updates,
        nonlinear_backtracks, number(physical.residual_norm_j)?,
        number(controls.residual_tolerance_j)?, number(physical.energy_residual_j)?,
        number(energy_tolerance_j)?, states,
    );
    poll(cx, deadline)?;
    Ok(receipt)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, controls: Config,
    old: &[f64], dt_s: f64,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    let coupling = ContactIterationConfig {
        max_iterations: controls.max_responses(),
        relaxation: 0.5,
        residual_tolerance_j: controls.residual_tolerance_j,
    };
    let accepted = if controls.coupled_newton {
        let nonlinear = nonlinear.ok_or_else(|| bad("coupled contact Newton lost its nonlinear controls"))?;
        contacts.advance_capacity_newton(
            cx, engine, problem, old, dt_s, config, nonlinear, coupling, controls.max_face_pairs,
        )
    } else {
        contacts.advance_capacity(cx, engine, problem, old, dt_s, config, nonlinear, coupling)
    }.map_err(|error| match error {
        cancelled @ ConductionError::Cancelled { .. } => lower(cancelled),
        other => bad(format!("source-backed contact step refused: {other}")),
    })?;
    let receipt = step_receipt(
        cx, controls, &accepted.physical, &accepted.physical_contacts, accepted.iterations,
        accepted.krylov_iterations, accepted.nonlinear_updates, accepted.nonlinear_backtracks,
        config.energy_tolerance_j, deadline,
    )?;
    let mut conduction = accepted.conduction;
    // The counted inner solve and physical audit describe the same unchanged
    // field. Publish freshly queried contact heat and physical reservoirs;
    // its inner Krylov residual remains explicitly numerical solve evidence.
    publish_physical(&mut conduction, &accepted.physical, &accepted.physical_contacts)?;
    poll(cx, deadline)?;
    Ok(Step {
        conduction, receipt, responses: accepted.iterations,
        krylov_iterations: accepted.krylov_iterations,
        nonlinear_updates: accepted.nonlinear_updates,
        nonlinear_backtracks: accepted.nonlinear_backtracks,
    })
}
