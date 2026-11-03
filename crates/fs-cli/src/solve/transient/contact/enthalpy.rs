//! Total enthalpy uses the native complete contact Newton or physical fixed point.

use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContacts};
use fs_conduction::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use fs_conduction::transient::enthalpy::{
    EnthalpyStepConfig, EnthalpyStepSolution,
    EnthalpyEndpointEvaluation,
    heterogeneous::HeterogeneousEnthalpyBackwardEuler,
};
use fs_conduction::ConductionProblem;
use fs_exec::Cx;
use std::time::Instant;

use crate::solve::{SolveRefusal, json_string, radiation};
use crate::solve::transient::{bad, number, poll, lower, enthalpy::{
    failed, advance_pointwise_contact, pointwise_contact_receipt,
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

pub(in crate::solve::transient) const NO_CLAIM: &str = "Source-backed matching contact resistance from scalar or continuous linear source curves and primary specific enthalpy are solved in one native Newton/FGMRES system with the complete contact-temperature tangent and each regional phase chart's temperature response. Temperature-source surfaces instead use the existing bounded relaxed contact fixed point; every inner enthalpy solve keeps the original old state, and acceptance re-queries the source at the unchanged endpoint and checks the complete physical residual and energy balance. Response, Newton and Krylov counts include discarded contact trials. A successful fixed point does not establish arbitrary-law convergence or uniqueness. Every residual and rejected trial retains the same accepted old enthalpy. The unchanged accepted field passes a freshly queried full nodal joule-residual gate and a separate storage/heat energy gate. Contact flux is internal transfer; prescribed reservoir reactions include boundary-node storage once. All retained contact rows describe accepted physical endpoints under the originally pinned material claims. Both time grids retain every workload and prescribed-boundary switch. Explicit pointwise ambient, closed-enclosure and open-enclosure radiation use the original local-emission and reflected-irradiation owners inside that same Newton residual and tangent. Their actual endpoint heat shares the independent physical acceptance gates; closed exchange is internal and open-environment heat is counted once. Pressure, geometry, source support, view factors and reference mass are fixed. Surface-law contacts combined with radiation, patch-mean radiation, derived natural/airflow feedback, phase-dependent contact resistance, contact trajectory derivatives, arbitrary-law convergence and continuum bounds are not inferred.";

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
        (true, Some(lowering)) => advance_pointwise_contact(
            cx, engine, problem, contacts, old_h, dt_s, config.native, lowering, enclosure,
        )?,
        (true, None) => contacts.advance_enthalpy_newton(
            cx, engine, problem, old_h, dt_s, config.native,
        ).map_err(failed)?,
        (false, None) => contacts.advance_enthalpy(
            cx, engine, problem, old_h, dt_s, config.native.step,
            ContactIterationConfig {
                max_iterations: config.max_responses(),
                relaxation: 0.5,
                residual_tolerance_j: config.native.residual_tolerance_j,
            },
        ).map_err(failed)?,
        (false, Some(_)) => return Err(bad(
            "surface-law enthalpy contact cannot use the pointwise-radiation Newton path",
        )),
    };
    poll(cx, deadline)?;
    let states = face_states(cx, &accepted.physical_contacts, deadline)?;
    let receipt = format!(concat!(
        "{{\"method\":{},\"thermal_responses\":{},",
        "\"krylov_iterations\":{},\"nonlinear_updates\":{},",
        "\"physical_residual_norm_j\":{},\"physical_residual_tolerance_j\":{},",
        "\"physical_energy_residual_j\":{},\"energy_tolerance_j\":{},\"face_states\":[{}]}}"
    ),
        json_string(config.method()), accepted.iterations,
        accepted.krylov_iterations, accepted.nonlinear_iterations,
        number(accepted.physical.residual_norm_j)?, number(config.native.residual_tolerance_j)?,
        number(accepted.physical.energy_residual_j)?, number(config.native.step.energy_tolerance_j)?,
        states,
    );
    let radiation_receipt = radiation.map(|lowering| pointwise_contact_receipt(
        &accepted.conduction, &accepted.physical, lowering, config.native.residual_tolerance_j,
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
    coarse_responses: usize,
    fine_responses: usize,
    nonlinear_updates: usize,
    krylov_iterations: usize,
) -> Result<String, SolveRefusal> {
    Ok(format!(concat!(
        ",\"temperature_contact\":{{\"method\":{},",
        "\"face_pairs\":{},\"max_thermal_responses_per_step\":{},\"max_residual_evaluations_per_step\":{},",
        "\"coarse_thermal_responses\":{},\"fine_thermal_responses\":{},\"total_thermal_responses\":{},",
        "\"total_nonlinear_updates\":{},\"total_krylov_iterations\":{},\"physical_residual_tolerance_j\":{},",
        "\"history\":\"unchanged-physical-old-specific-enthalpy-in-every-contact-and-newton-trial\",",
        "\"scope\":{}}}"
    ),
        json_string(config.method()), config.native.max_face_pairs, config.max_responses(),
        config.max_residual_evaluations()?, coarse_responses, fine_responses,
        coarse_responses.checked_add(fine_responses)
            .ok_or_else(|| bad("contact enthalpy response count overflows"))?,
        nonlinear_updates, krylov_iterations,
        number(config.native.residual_tolerance_j)?, json_string(NO_CLAIM),
    ))
}
