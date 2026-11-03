//! Total enthalpy and source-backed contact share one native Newton endpoint.

use fs_conduction::interface::temperature::TemperatureContacts;
use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use fs_conduction::transient::enthalpy::{
    EnthalpyStepConfig, EnthalpyStepSolution,
    EnthalpyEndpointEvaluation,
    heterogeneous::HeterogeneousEnthalpyBackwardEuler,
};
use fs_conduction::ConductionProblem;
use fs_exec::Cx;
use std::time::Instant;

use crate::solve::{SolveRefusal, json_string};
use crate::solve::transient::{bad, number, poll, lower, enthalpy::failed};
use super::{face_states, residual_tolerance_j};

pub(in crate::solve::transient) const NO_CLAIM: &str = "Source-backed matching contact resistance from scalar or continuous linear source curves and primary specific enthalpy are solved in one native Newton/FGMRES system with the complete contact-temperature tangent and each regional phase chart's temperature response. Every residual and rejected trial retains the same accepted old enthalpy. The unchanged accepted field passes a freshly queried full nodal joule-residual gate and a separate storage/heat energy gate. Contact flux is internal transfer; prescribed reservoir reactions include boundary-node storage once. All retained contact rows describe accepted physical endpoints under the originally pinned material claims. Both time grids retain every workload and prescribed-boundary switch. Pressure, geometry, source support and reference mass are fixed. Surface-law contact fixed points, radiation, derived natural/airflow feedback, phase-dependent contact resistance and contact trajectory derivatives are not inferred.";

pub(in crate::solve::transient) struct Step {
    pub(in crate::solve::transient) conduction: EnthalpyStepSolution,
    pub(in crate::solve::transient) physical: EnthalpyEndpointEvaluation,
    pub(in crate::solve::transient) receipt: String,
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
) -> Result<ContactEnthalpyNewtonConfig, SolveRefusal> {
    // Initial residual plus each update's current residual and every admitted
    // full-step/backtrack candidate. The native owner enforces this callback cap.
    let max_residual_evaluations = max_backtracks.checked_add(2)
        .and_then(|evaluations| evaluations.checked_mul(step.max_newton_iterations))
        .and_then(|evaluations| evaluations.checked_add(1))
        .ok_or_else(|| bad("contact enthalpy residual evaluation allowance overflows"))?;
    Ok(ContactEnthalpyNewtonConfig {
        step,
        max_residual_evaluations,
        max_face_pairs: face_pairs,
        residual_tolerance_j: residual_tolerance_j(step.energy_tolerance_j, vertices)?,
    })
}

#[allow(clippy::too_many_arguments)]
pub(in crate::solve::transient) fn advance(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>,
    old_h: &[f64],
    dt_s: f64,
    config: ContactEnthalpyNewtonConfig,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    let accepted = contacts.advance_enthalpy_newton(
        cx, engine, problem, old_h, dt_s, config,
    ).map_err(failed)?;
    poll(cx, deadline)?;
    let states = face_states(cx, &accepted.physical_contacts, deadline)?;
    let receipt = format!(concat!(
        "{{\"method\":\"source-backed-contact-enthalpy-newton\",\"thermal_responses\":{},",
        "\"krylov_iterations\":{},\"nonlinear_updates\":{},",
        "\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},",
        "\"physical_energy_residual_j\":{},\"energy_tolerance_j\":{},\"face_states\":[{}]}}"
    ),
        accepted.iterations, accepted.krylov_iterations, accepted.nonlinear_iterations,
        number(accepted.physical.residual_norm_j)?, number(config.residual_tolerance_j)?,
        number(accepted.physical.energy_residual_j)?, number(config.step.energy_tolerance_j)?,
        states,
    );
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
    Ok(Step { conduction, physical: accepted.physical, receipt })
}

pub(in crate::solve::transient) fn receipt(
    config: ContactEnthalpyNewtonConfig,
    coarse_steps: usize,
    fine_steps: usize,
    nonlinear_updates: usize,
) -> Result<String, SolveRefusal> {
    Ok(format!(concat!(
        ",\"temperature_contact\":{{\"method\":\"source-backed-contact-enthalpy-newton\",",
        "\"face_pairs\":{},\"max_thermal_responses_per_step\":1,\"max_residual_evaluations_per_step\":{},",
        "\"coarse_thermal_responses\":{},\"fine_thermal_responses\":{},\"total_thermal_responses\":{},",
        "\"total_nonlinear_updates\":{},\"physical_residual_tolerance_j\":{},",
        "\"history\":\"unchanged-physical-old-specific-enthalpy-in-every-newton-trial\",",
        "\"scope\":{}}}"
    ),
        config.max_face_pairs, config.max_residual_evaluations,
        coarse_steps, fine_steps, coarse_steps + fine_steps, nonlinear_updates,
        number(config.residual_tolerance_j)?, json_string(NO_CLAIM),
    ))
}
