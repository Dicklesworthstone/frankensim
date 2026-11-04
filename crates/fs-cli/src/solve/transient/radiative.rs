//! Select the physical radiation owner once; all capacity coupling loops share
//! its solve, work evidence, and fresh endpoint audit.

use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContacts};
use fs_conduction::transient::backward_euler::EndpointEvaluation;
use super::*;

#[derive(Clone, Copy)]
pub(super) struct PhysicalEnergy {
    pub(super) neumann_out_w: f64,
    pub(super) robin_out_w: f64,
    pub(super) dirichlet_in_w: f64,
    pub(super) net_input_w: f64,
    pub(super) energy_residual_j: f64,
}

impl PhysicalEnergy {
    pub(super) fn inner(step: &StepSolution) -> Self {
        Self {
            neumann_out_w: step.neumann_out_w,
            robin_out_w: step.robin_out_w,
            dirichlet_in_w: step.dirichlet_in_w,
            net_input_w: step.source_w + step.dirichlet_in_w
                - step.neumann_out_w - step.robin_out_w,
            energy_residual_j: step.energy_residual_j,
        }
    }

    pub(super) fn physical(endpoint: &EndpointEvaluation) -> Self {
        let robin_out_w = endpoint.convective_out_w + endpoint.radiation_out_w;
        Self {
            neumann_out_w: endpoint.neumann_out_w,
            robin_out_w,
            dirichlet_in_w: endpoint.dirichlet_in_w,
            net_input_w: endpoint.source_w + endpoint.dirichlet_in_w
                - endpoint.neumann_out_w - robin_out_w,
            energy_residual_j: endpoint.energy_residual_j,
        }
    }
}

pub(super) struct Step {
    pub(super) conduction: StepSolution,
    pub(super) nonlinear_evidence: String,
    pub(super) radiation_evidence: String,
    pub(super) contact_evidence: Option<String>,
    pub(super) radiation: RadiativeEndpoint,
    /// Radiation callbacks, distinct from their nested contact responses.
    pub(super) radiation_trials: usize,
    pub(super) solid_solves: usize,
    pub(super) nonlinear_updates: usize,
    pub(super) nonlinear_backtracks: usize,
    pub(super) krylov_iterations: usize,
    pub(super) physical_residual_tolerance_j: f64,
    pub(super) physical_energy: PhysicalEnergy,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    temperature_contacts: Option<&TemperatureContacts<'_, '_>>,
    contact_config: Option<contact::Config>, old: &[f64], dt_s: f64,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
    lowering: &radiation::LoweredRadiation, model: Option<&enclosure::Model>,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    if temperature_contacts.is_some() != contact_config.is_some()
        || (temperature_contacts.is_some() && interfaces.is_some())
    {
        return Err(bad("capacity radiation requires matching live contacts and controls without duplicate frozen interfaces"));
    }
    if lowering.is_enclosure() != model.is_some()
        || lowering.open_ambient_temperature_k().is_some() != model.is_some_and(enclosure::Model::is_open)
    {
        return Err(bad("capacity radiation iteration lost its exact enclosure model"));
    }
    if temperature_contacts.is_some() {
        contact::admit_radiation(Some(lowering), model)?;
    }
    if let Some(contacts) = temperature_contacts.filter(|_| model.is_none()) {
        let controls = contact_config
            .ok_or_else(|| bad("capacity radiation lost its live contact controls"))?;
        let coupling = ContactIterationConfig {
            max_iterations: controls.max_responses(), relaxation: 0.5,
            residual_tolerance_j: controls.residual_tolerance_j,
        };
        let solved = if controls.coupled_newton {
            contacts.advance_capacity_newton_with_ambient_radiation(
                cx, engine, problem, old, dt_s, config,
                nonlinear.ok_or_else(|| bad("radiative contact Newton lost its nonlinear controls"))?,
                coupling, controls.max_face_pairs, &lowering.patches, lowering.config,
            )
        } else {
            contacts.advance_capacity_with_ambient_radiation(
                cx, engine, problem, old, dt_s, config, nonlinear,
                coupling, &lowering.patches, lowering.config,
            )
        }.map_err(lower)?;
        poll(cx, deadline)?;
        let physical = &solved.physical;
        let contact_evidence = contact::step_receipt(
            cx, contact::Config {
                residual_tolerance_j: solved.physical_residual_tolerance_j, ..controls
            }, physical, &solved.physical_contacts, solved.thermal_responses,
            solved.krylov_iterations, solved.nonlinear_updates, solved.nonlinear_backtracks,
            config.energy_tolerance_j, deadline,
        )?;
        let radiation_evidence = radiation_row(
            &solved.radiation, physical.convective_out_w, physical.residual_norm_j,
            solved.physical_residual_tolerance_j, physical.energy_residual_j,
            physical.dirichlet_in_w,
        )?;
        let physical_energy = PhysicalEnergy::physical(physical);
        let mut conduction = solved.conduction;
        contact::publish_physical(&mut conduction, physical, &solved.physical_contacts)?;
        return Ok(Step {
            conduction, nonlinear_evidence: "null".to_string(), radiation_evidence,
            contact_evidence: Some(contact_evidence),
            radiation_trials: solved.radiation.iterations,
            solid_solves: solved.thermal_responses,
            nonlinear_updates: solved.nonlinear_updates,
            nonlinear_backtracks: solved.nonlinear_backtracks,
            krylov_iterations: solved.krylov_iterations,
            physical_residual_tolerance_j: solved.physical_residual_tolerance_j,
            physical_energy,
            radiation: RadiativeEndpoint::Ambient {
                // The retained numerical radiation report still describes its
                // applied secants; downstream physical heat uses this audit.
                combined_boundary: physical.combined_boundary.clone(),
                convective_robin_fluxes: physical.convective_robin_fluxes.clone(),
                convective_out_w: physical.convective_out_w,
                solid_solves: solved.thermal_responses,
                report: solved.radiation,
            },
        });
    }
    if let Some(model) = model {
        let solved = enclosure::advance(
            cx, engine, problem, interfaces, temperature_contacts, contact_config,
            old, dt_s, config, nonlinear,
            model, lowering.config, deadline,
        )?;
        let nonlinear_evidence = solved.nonlinear.as_ref().map(nonlinear_row)
            .transpose()?.unwrap_or_else(|| "null".to_string());
        let report = &solved.endpoint.report;
        let physical_energy = PhysicalEnergy {
            neumann_out_w: solved.physical_neumann_out_w,
            robin_out_w: solved.endpoint.convective_out_w
                + model.surface_heat_w(&report.radiosity)?,
            dirichlet_in_w: solved.physical_dirichlet_in_w,
            net_input_w: solved.net_input_w,
            energy_residual_j: solved.endpoint.physical_energy_residual_j,
        };
        return Ok(Step {
            conduction: solved.conduction,
            nonlinear_evidence,
            radiation_evidence: solved.row,
            contact_evidence: solved.contact_evidence,
            radiation_trials: report.iterations,
            solid_solves: solved.solid_solves,
            nonlinear_updates: solved.nonlinear_updates,
            nonlinear_backtracks: solved.nonlinear_backtracks,
            krylov_iterations: report.krylov_iterations,
            physical_residual_tolerance_j: solved.endpoint.physical_residual_tolerance_j,
            physical_energy,
            radiation: RadiativeEndpoint::Enclosure {
                endpoint: solved.endpoint, solid_solves: solved.solid_solves,
                nonlinear_backtracks: solved.nonlinear_backtracks,
            },
        });
    }
    poll(cx, deadline)?;
    let solved = engine.advance_prescribed_with_ambient_radiation(
        cx, problem, interfaces, old, dt_s, config, nonlinear,
        &lowering.patches, lowering.config,
    ).map_err(lower)?;
    poll(cx, deadline)?;
    let nonlinear_evidence = solved.nonlinear.as_ref().map(nonlinear_row)
        .transpose()?.unwrap_or_else(|| "null".to_string());
    let radiation_evidence = radiation_row(
        &solved.radiation, solved.convective_out_w, solved.physical_residual_norm_j,
        solved.physical_residual_tolerance_j, solved.physical_energy_residual_j,
        solved.physical_dirichlet_in_w,
    )?;
    let robin_out_w = solved.convective_out_w + solved.radiation.nonlinear_radiation_out_w;
    let physical_energy = PhysicalEnergy {
        neumann_out_w: solved.conduction.neumann_out_w,
        robin_out_w,
        dirichlet_in_w: solved.physical_dirichlet_in_w,
        net_input_w: solved.conduction.source_w + solved.physical_dirichlet_in_w
            - solved.conduction.neumann_out_w - robin_out_w,
        energy_residual_j: solved.physical_energy_residual_j,
    };
    Ok(Step {
        conduction: solved.conduction,
        nonlinear_evidence,
        radiation_evidence,
        contact_evidence: None,
        radiation_trials: solved.radiation.iterations,
        solid_solves: solved.radiation.iterations,
        nonlinear_updates: if nonlinear.is_some() { solved.radiation.solid_iterations } else { 0 },
        nonlinear_backtracks: solved.nonlinear_backtracks,
        krylov_iterations: solved.radiation.krylov_iterations,
        physical_residual_tolerance_j: solved.physical_residual_tolerance_j,
        physical_energy,
        radiation: RadiativeEndpoint::Ambient {
            combined_boundary: solved.combined_boundary,
            convective_robin_fluxes: solved.convective_robin_fluxes,
            convective_out_w: solved.convective_out_w,
            solid_solves: solved.radiation.iterations,
            report: solved.radiation,
        },
    })
}

/// A surrounding natural/air solve changes convection at the same unchanged
/// field. Keep its last counted radiation response while publishing the
/// fresh complete physical boundary, residual and separated convective heat.
pub(super) fn refresh(
    cx: &Cx<'_>,
    endpoint: &mut RadiativeEndpoint,
    model: Option<&enclosure::Model>,
    physical: &EndpointEvaluation,
    physical_residual_tolerance_j: f64,
    deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    poll(cx, deadline)?;
    match endpoint {
        RadiativeEndpoint::Ambient {
            combined_boundary, convective_robin_fluxes, convective_out_w, report, ..
        } if model.is_none() => {
            *combined_boundary = physical.combined_boundary.clone();
            *convective_robin_fluxes = physical.convective_robin_fluxes.clone();
            *convective_out_w = physical.convective_out_w;
            radiation_row(
                report, physical.convective_out_w, physical.residual_norm_j,
                physical_residual_tolerance_j, physical.energy_residual_j, physical.dirichlet_in_w,
            )
        }
        RadiativeEndpoint::Enclosure { endpoint, solid_solves, nonlinear_backtracks } => {
            let model = model.ok_or_else(|| bad("capacity radiation refresh lost its enclosure"))?;
            enclosure::refresh(
                cx, endpoint, model, physical, physical_residual_tolerance_j,
                *solid_solves, *nonlinear_backtracks, deadline,
            )
        }
        _ => Err(bad("capacity radiation refresh has inconsistent ambient/enclosure ownership")),
    }
}

/// No primal step is taken here. Candidate and old temperature retain their
/// identities while the actual convection and reflected exchange are audited.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>, old: &[f64], candidate: &[f64],
    dt_s: f64, config: StepConfig, lowering: Option<&radiation::LoweredRadiation>,
    model: Option<&enclosure::Model>, deadline: Option<(Instant, f64)>,
) -> Result<EndpointEvaluation, SolveRefusal> {
    poll(cx, deadline)?;
    if lowering.is_some_and(radiation::LoweredRadiation::is_enclosure) != model.is_some() {
        return Err(bad("capacity physical boundary audit lost its declared radiation model"));
    }
    let physical = if let Some(model) = model {
        let lowering = lowering.ok_or_else(|| bad("enclosure endpoint lost its radiation policy"))?;
        match model {
            enclosure::Model::Closed(model) => engine.evaluate_prescribed_enclosure_endpoint(
                cx, problem, interfaces, old, candidate, dt_s, config, model, lowering.config),
            enclosure::Model::Open(model, ambient) => engine.evaluate_prescribed_open_enclosure_endpoint(
                cx, problem, interfaces, old, candidate, dt_s, config, model, *ambient, lowering.config),
        }.map_err(lower)?.physical
    } else {
        engine.evaluate_prescribed_endpoint(
            cx, problem, interfaces, old, candidate, dt_s, config,
            lowering.map(|row| row.patches.as_slice()),
        ).map_err(lower)?
    };
    poll(cx, deadline)?;
    Ok(physical)
}

