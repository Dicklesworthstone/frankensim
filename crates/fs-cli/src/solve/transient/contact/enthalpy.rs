//! Total enthalpy uses the native complete contact Newton or physical fixed point.

use fs_conduction::interface::temperature::{
    ContactIterationConfig, TemperatureContactEvaluation, TemperatureContacts,
};
use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use fs_conduction::transient::enthalpy::{
    EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
    EnthalpyEndpointEvaluation,
    heterogeneous::HeterogeneousEnthalpyBackwardEuler,
};
use fs_conduction::ConductionProblem;
use fs_exec::Cx;
use std::time::Instant;

use crate::solve::{SolveRefusal, json_string, radiation};
use crate::solve::transient::{bad, number, poll, lower, enthalpy::{
    failed, advance_pointwise_contact, evaluate_pointwise_endpoint, pointwise_contact_receipt,
}};
use super::{face_states, residual_tolerance_j};

/// Scalar and curve laws retain the complete contact Newton tangent. Source
/// surfaces use the existing bounded fixed point, without requiring a unique
/// derivative where adjacent source triangles meet.
#[derive(Clone, Copy)]
pub(in crate::solve::transient) struct Config {
    pub(in crate::solve::transient) native: ContactEnthalpyNewtonConfig,
    pub(in crate::solve::transient) coupled_newton: bool,
}

impl Config {
    pub(in crate::solve::transient) fn max_responses(self) -> usize {
        super::max_responses(self.coupled_newton)
    }

    pub(in crate::solve::transient) fn coupling(self) -> ContactIterationConfig {
        ContactIterationConfig {
            max_iterations: self.max_responses(),
            relaxation: 0.5,
            residual_tolerance_j: self.native.residual_tolerance_j,
        }
    }

    fn method(self) -> &'static str {
        if self.coupled_newton {
            "source-backed-contact-enthalpy-newton"
        } else {
            "source-backed-contact-enthalpy-fixed-point"
        }
    }

    pub(in crate::solve::transient) fn max_residual_evaluations(self) -> Result<usize, SolveRefusal> {
        self.native.max_residual_evaluations.checked_mul(self.max_responses())
            .ok_or_else(|| bad("contact enthalpy residual allowance overflows"))
    }
}

pub(in crate::solve::transient) const NO_CLAIM: &str = "Source-backed matching contact resistance from scalar or continuous linear source curves and primary specific enthalpy are solved in one native Newton/FGMRES system with the complete contact-temperature tangent and each regional phase chart's temperature response. Temperature-source surfaces instead use the existing bounded relaxed contact fixed point; every inner enthalpy solve keeps the original old state, and acceptance re-queries the source at the unchanged endpoint and checks the complete physical residual and energy balance. Response, Newton and Krylov counts include discarded contact trials. A successful fixed point does not establish arbitrary-law convergence or uniqueness. Every residual and rejected trial retains the same accepted old enthalpy. The unchanged accepted field passes a freshly queried full nodal joule-residual gate and a separate storage/heat energy gate. Contact flux is internal transfer; prescribed reservoir reactions include boundary-node storage once. All retained contact rows describe accepted physical endpoints under the originally pinned material claims. Both time grids retain every workload and prescribed-boundary switch. Explicit pointwise ambient, closed-enclosure and open-enclosure radiation use the original local-emission and reflected-irradiation owners inside every inner enthalpy Newton residual and tangent, including every source-surface contact response. Their actual endpoint heat shares the independent physical acceptance gates; closed exchange is internal and open-environment heat is counted once. Natural convection encloses these contact-aware enthalpy responses, re-queries the original card at the candidate wall temperature, and rechecks the unchanged endpoint with both its fresh contact law and updated physical convection coefficient. Quasi-steady transported-air callbacks use the same contact-aware responses, including nested natural convection on separate targets. The unchanged endpoint passes fresh contact, updated-air-reference, physical residual and actual-air energy gates. Only forced convection heats transported air; natural reservoirs and radiation retain separate heat accounts. All nested response, Newton and Krylov work is counted. Pressure, geometry, source support, view factors and reference mass are fixed. Patch-mean radiation, phase-dependent contact resistance, contact trajectory derivatives, arbitrary-law convergence and continuum bounds are not inferred.";

pub(in crate::solve::transient) struct Step {
    pub(in crate::solve::transient) conduction: EnthalpyStepSolution,
    pub(in crate::solve::transient) physical: EnthalpyEndpointEvaluation,
    pub(in crate::solve::transient) receipt: String,
    pub(in crate::solve::transient) radiation_receipt: Option<String>,
    pub(in crate::solve::transient) responses: usize,
    pub(in crate::solve::transient) nonlinear_updates: usize,
    pub(in crate::solve::transient) krylov_iterations: usize,
}

/// Retain both grids' physical face rows and bounded contact tangent scratch
/// before either history is allocated. Include temporary serialization copies;
/// this is declared-memory admission, not allocator overhead or peak RSS.
pub(in crate::solve::transient) fn storage_bytes(
    face_pairs: usize,
    total_steps: usize,
) -> Result<usize, SolveRefusal> {
    if face_pairs == 0 {
        return Err(bad("temperature-dependent enthalpy contact has no matching face pairs"));
    }
    face_pairs.checked_mul(total_steps)
        .and_then(|rows| rows.checked_mul(4096))
        .and_then(|rows| face_pairs.checked_mul(32768).and_then(|scratch| rows.checked_add(scratch)))
        .ok_or_else(|| bad("enthalpy contact face-state and tangent allowance overflows"))
}

pub(in crate::solve::transient) fn config(
    step: EnthalpyStepConfig,
    face_pairs: usize,
    max_backtracks: usize,
    vertices: usize,
    coupled_newton: bool,
) -> Result<Config, SolveRefusal> {
    // Initial residual plus each update's current residual and every admitted
    // full-step/backtrack candidate. The native owner enforces this callback cap.
    let max_residual_evaluations = max_backtracks.checked_add(2)
        .and_then(|evaluations| evaluations.checked_mul(step.max_newton_iterations))
        .and_then(|evaluations| evaluations.checked_add(1))
        .ok_or_else(|| bad("contact enthalpy residual evaluation allowance overflows"))?;
    Ok(Config {
        native: ContactEnthalpyNewtonConfig {
            step,
            max_residual_evaluations,
            max_face_pairs: face_pairs,
            residual_tolerance_j: residual_tolerance_j(step.energy_tolerance_j, vertices)?,
        },
        coupled_newton,
    })
}

/// Re-query the original contact laws at an unchanged primary enthalpy field,
/// then audit it under the caller's actual physical boundary. The original
/// per-vertex charts resolve temperature directly, including latent plateaus;
/// no temperature-to-enthalpy inversion or replacement PDE solve occurs.
/// The returned contact rows belong to this exact physical evaluation.
#[allow(clippy::too_many_arguments)]
pub(in crate::solve::transient) fn evaluate(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>,
    old_h: &[f64],
    candidate_h: &[f64],
    dt_s: f64,
    config: EnthalpyStepConfig,
    radiation: Option<&radiation::LoweredRadiation>,
    enclosure: Option<&fs_conduction::radiation::GrayDiffuseEnclosure>,
    deadline: Option<(Instant, f64)>,
) -> Result<(EnthalpyEndpointEvaluation, TemperatureContactEvaluation), SolveRefusal> {
    poll(cx, deadline)?;
    let vertices = problem.mesh.vertex_count();
    if old_h.len() != vertices || candidate_h.len() != vertices {
        return Err(bad("contact enthalpy audit needs complete old and candidate fields on the original mesh"));
    }
    if radiation.is_none() && enclosure.is_some() {
        return Err(bad("contact enthalpy audit lost its declared pointwise radiation model"));
    }
    let mut temperature = Vec::new();
    temperature.try_reserve_exact(vertices)
        .map_err(|_| bad("contact enthalpy audit temperature allocation refused"))?;
    for (vertex, &enthalpy) in candidate_h.iter().enumerate() {
        if vertex % 256 == 0 { poll(cx, deadline)?; }
        let chart = engine.phase_curve_at(vertex)
            .ok_or_else(|| bad("contact enthalpy audit lost an original vertex phase chart"))?;
        let state = chart.state_at_specific_enthalpy(enthalpy)
            .map_err(|source| failed(EnthalpyError::Phase { vertex, source }))?;
        temperature.push(state.temperature_k());
    }
    poll(cx, deadline)?;
    let actual = contacts.evaluate(cx, problem.boundary, &temperature).map_err(lower)?;
    poll(cx, deadline)?;
    let physical = match radiation {
        Some(lowering) => evaluate_pointwise_endpoint(
            cx, engine, problem, Some(actual.frozen_interfaces()), old_h, candidate_h,
            dt_s, config, lowering, enclosure,
        )?,
        None => engine.evaluate_endpoint(
            cx, problem, Some(actual.frozen_interfaces()), old_h, candidate_h,
            dt_s, config, None,
        ).map_err(failed)?,
    };
    poll(cx, deadline)?;
    Ok((physical, actual))
}

/// Serialize the same audited physical field and contact queries, with the
/// caller's complete counted work. Outer boundary iteration may replace the
/// applied-boundary audit without losing its accepted source-state evidence.
#[allow(clippy::too_many_arguments)]
pub(in crate::solve::transient) fn step_receipt(
    cx: &Cx<'_>,
    config: Config,
    physical: &EnthalpyEndpointEvaluation,
    contacts: &TemperatureContactEvaluation,
    responses: usize,
    nonlinear_updates: usize,
    krylov_iterations: usize,
    deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    poll(cx, deadline)?;
    let states = face_states(cx, contacts, deadline)?;
    let receipt = format!(concat!(
        "{{\"method\":{},\"thermal_responses\":{},",
        "\"krylov_iterations\":{},\"nonlinear_updates\":{},",
        "\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},",
        "\"physical_energy_residual_j\":{},\"energy_tolerance_j\":{},\"face_states\":[{}]}}"
    ),
        json_string(config.method()), responses, krylov_iterations, nonlinear_updates,
        number(physical.residual_norm_j)?, number(config.native.residual_tolerance_j)?,
        number(physical.energy_residual_j)?, number(config.native.step.energy_tolerance_j)?,
        states,
    );
    poll(cx, deadline)?;
    Ok(receipt)
}

#[allow(clippy::too_many_arguments)]
pub(in crate::solve::transient) fn advance(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>,
    old_h: &[f64],
    dt_s: f64,
    config: Config,
    radiation: Option<&radiation::LoweredRadiation>,
    enclosure: Option<&fs_conduction::radiation::GrayDiffuseEnclosure>,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    let accepted = match (config.coupled_newton, radiation) {
        (_, Some(lowering)) => advance_pointwise_contact(
            cx, engine, problem, contacts, old_h, dt_s, config, lowering, enclosure,
        )?,
        (true, None) => contacts.advance_enthalpy_newton(
            cx, engine, problem, old_h, dt_s, config.native,
        ).map_err(failed)?,
        (false, None) => contacts.advance_enthalpy(
            cx, engine, problem, old_h, dt_s, config.native.step,
            config.coupling(),
        ).map_err(failed)?,
    };
    poll(cx, deadline)?;
    let receipt = step_receipt(
        cx, config, &accepted.physical, &accepted.physical_contacts,
        accepted.iterations, accepted.nonlinear_iterations, accepted.krylov_iterations,
        deadline,
    )?;
    let radiation_receipt = radiation.map(|lowering| pointwise_contact_receipt(
        &accepted.conduction, &accepted.physical, lowering, config.native.residual_tolerance_j,
        accepted.iterations, accepted.nonlinear_iterations, accepted.krylov_iterations,
    )).transpose()?;
    let mut conduction = accepted.conduction;
    // Publish the fresh accepted physical operator and heat reactions. Newton
    // history remains the original counted numerical solve; storage is never
    // advanced a second time to manufacture endpoint evidence.
    conduction.contact_fluxes = accepted.physical_contacts.frozen_interfaces()
        .fluxes(&conduction.temperature).map_err(lower)?;
    conduction.robin_fluxes = accepted.physical.convective_robin_fluxes.clone();
    conduction.source_w = accepted.physical.source_w;
    conduction.neumann_out_w = accepted.physical.neumann_out_w;
    conduction.robin_out_w = accepted.physical.convective_out_w;
    conduction.dirichlet_in_w = accepted.physical.dirichlet_in_w;
    conduction.stored_energy_change_j = accepted.physical.stored_energy_change_j;
    conduction.energy_residual_j = accepted.physical.energy_residual_j;
    poll(cx, deadline)?;
    Ok(Step {
        conduction, physical: accepted.physical, receipt, radiation_receipt,
        responses: accepted.iterations, nonlinear_updates: accepted.nonlinear_iterations,
        krylov_iterations: accepted.krylov_iterations,
    })
}

pub(in crate::solve::transient) fn receipt(
    config: Config,
    boundary_responses: usize,
    coarse_responses: usize,
    fine_responses: usize,
    nonlinear_updates: usize,
    krylov_iterations: usize,
) -> Result<String, SolveRefusal> {
    let max_responses = config.max_responses().checked_mul(boundary_responses)
        .ok_or_else(|| bad("coupled contact enthalpy response allowance overflows"))?;
    let max_residual_evaluations = config.max_residual_evaluations()?
        .checked_mul(boundary_responses)
        .ok_or_else(|| bad("coupled contact enthalpy residual allowance overflows"))?;
    Ok(format!(concat!(
        ",\"temperature_contact\":{{\"method\":{},",
        "\"face_pairs\":{},\"max_thermal_responses_per_step\":{},\"max_residual_evaluations_per_step\":{},",
        "\"coarse_thermal_responses\":{},\"fine_thermal_responses\":{},\"total_thermal_responses\":{},",
        "\"total_nonlinear_updates\":{},\"total_krylov_iterations\":{},\"physical_residual_tolerance_j\":{},",
        "\"history\":\"unchanged-physical-old-specific-enthalpy-in-every-contact-and-newton-trial\",",
        "\"scope\":{}}}"
    ),
        json_string(config.method()), config.native.max_face_pairs, max_responses,
        max_residual_evaluations, coarse_responses, fine_responses,
        coarse_responses.checked_add(fine_responses)
            .ok_or_else(|| bad("contact enthalpy response count overflows"))?,
        nonlinear_updates, krylov_iterations,
        number(config.native.residual_tolerance_j)?, json_string(NO_CLAIM),
    ))
}
