//! Native capacity trajectories use the shared reflected-radiation endpoint.
//! Returned inner solve data retain their applied boundary; fresh physical heat
//! is explicit so outer convection owners never advance or relabel history.

use fs_conduction::radiation::{GrayDiffuseEnclosure, OpenDiffuseEnclosure, RadiationSurface, RadiosityReport};
use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContacts};
use fs_conduction::transient::backward_euler::EndpointEvaluation;
use super::*;

/// One geometric binding per spatial rung, shared by both time grids and all
/// natural/air/contact trials. The black reservoir never becomes a solid patch.
pub(super) enum Model {
    Closed(GrayDiffuseEnclosure),
    Open(OpenDiffuseEnclosure, f64),
}

impl Model {
    pub(super) fn bind(cx: &Cx<'_>, problem: ConductionProblem<'_>, lowering: &radiation::LoweredRadiation)
        -> Result<Option<Self>, ConductionError> {
        if let Some(ambient) = lowering.open_ambient_temperature_k() {
            lowering.bind_open_enclosure(cx, problem).map(|model| model.map(|model| Self::Open(model, ambient)))
        } else {
            lowering.bind_enclosure(cx, problem).map(|model| model.map(Self::Closed))
        }
    }

    pub(super) fn surfaces(&self) -> &[RadiationSurface] {
        match self { Self::Closed(model) => model.surfaces(), Self::Open(model, _) => model.surfaces() }
    }

    pub(super) fn is_open(&self) -> bool { matches!(self, Self::Open(..)) }

    pub(super) fn closed(&self) -> Result<&GrayDiffuseEnclosure, SolveRefusal> {
        match self {
            Self::Closed(model) => Ok(model),
            Self::Open(..) => Err(bad("open area-mean capacity radiation does not expose trajectory adjoints")),
        }
    }

    pub(super) fn surface_heat_w(&self, report: &RadiosityReport) -> Result<f64, SolveRefusal> {
        let heat = report.net_outward_heat_w.get(..self.surfaces().len())
            .ok_or_else(|| bad("capacity enclosure lost its physical-surface heat rows"))?;
        finite(heat.iter().sum())
    }

    pub(super) fn into_solution(self, conduction: ConductionSolution, endpoint: radiation::EnclosureEndpoint)
        -> radiation::SolidSolution {
        match self {
            Self::Closed(model) => radiation::SolidSolution::from_enclosure_endpoint(conduction, model, endpoint),
            Self::Open(model, ambient) => radiation::SolidSolution::from_open_enclosure_endpoint(
                conduction, model, ambient, endpoint),
        }
    }
}

pub(super) struct Step {
    pub(super) conduction: StepSolution,
    pub(super) nonlinear: Option<NonlinearStepSolution>,
    pub(super) contact_evidence: Option<String>,
    pub(super) solid_solves: usize,
    pub(super) nonlinear_updates: usize,
    pub(super) nonlinear_backtracks: usize,
    pub(super) endpoint: radiation::EnclosureEndpoint,
    pub(super) net_input_w: f64,
    pub(super) physical_neumann_out_w: f64,
    pub(super) physical_dirichlet_in_w: f64,
    pub(super) row: String,
}

fn finite(value: f64) -> Result<f64, SolveRefusal> {
    if value.is_finite() { Ok(value) }
    else { Err(bad("capacity enclosure arithmetic is not finite")) }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    temperature_contacts: Option<&TemperatureContacts<'_, '_>>,
    contact_config: Option<contact::Config>, old: &[f64], dt_s: f64,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
    model: &Model, radiation_config: fs_conduction::AmbientRadiationConfig,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    if temperature_contacts.is_some() != contact_config.is_some()
        || (temperature_contacts.is_some() && interfaces.is_some())
    {
        return Err(bad("capacity enclosure requires matching live contacts and controls without duplicate frozen interfaces"));
    }
    if let Some(contacts) = temperature_contacts {
        return advance_contact(
            cx, engine, problem, contacts,
            contact_config.ok_or_else(|| bad("capacity enclosure lost its contact controls"))?,
            old, dt_s, config, nonlinear, model, radiation_config, deadline,
        );
    }
    let solved = match model {
        Model::Closed(model) => engine.advance_prescribed_with_enclosure_radiation(
            cx, problem, interfaces, old, dt_s, config, nonlinear, model, radiation_config),
        Model::Open(model, ambient) => engine.advance_prescribed_with_open_enclosure_radiation(
            cx, problem, interfaces, old, dt_s, config, nonlinear, model, *ambient, radiation_config),
    }.map_err(lower)?;
    poll(cx, deadline)?;
    let report = &solved.radiation;
    let fresh = &report.radiosity;
    let applied = finite(report.applied_heat_w.iter().sum())?;
    let decomposition = finite(
        solved.conduction.robin_out_w + solved.conduction.neumann_out_w
            - solved.convective_out_w - solved.prescribed_neumann_out_w - applied,
    )?;
    if decomposition.abs() > radiation_config.balance_tolerance_w {
        return Err(bad("capacity enclosure failed its original boundary decomposition allowance"));
    }
    let net_input_w = finite(
        solved.conduction.source_w + solved.physical_dirichlet_in_w
            - solved.prescribed_neumann_out_w - solved.convective_out_w
            - model.surface_heat_w(fresh)?,
    )?;
    let solid_solves = report.iterations;
    let nonlinear_updates = if nonlinear.is_some() { report.solid_iterations } else { 0 };
    let endpoint = radiation::EnclosureEndpoint {
        combined_boundary: solved.combined_boundary,
        convective_robin_fluxes: solved.convective_robin_fluxes,
        convective_out_w: solved.convective_out_w,
        report: solved.radiation,
        dt_s,
        time_integration: "backward-euler-capacity",
        physical_residual_norm_j: solved.physical_residual_norm_j,
        physical_residual_tolerance_j: solved.physical_residual_tolerance_j,
        physical_energy_residual_j: solved.physical_energy_residual_j,
        decomposition_residual_w: decomposition,
    };
    let row = step_row(
        cx, model, &endpoint, solid_solves, solved.nonlinear_backtracks,
        solved.prescribed_neumann_out_w, solved.physical_dirichlet_in_w, deadline,
    )?;
    Ok(Step {
        conduction: solved.conduction,
        nonlinear: solved.nonlinear,
        contact_evidence: None,
        solid_solves,
        nonlinear_updates,
        nonlinear_backtracks: solved.nonlinear_backtracks,
        endpoint,
        net_input_w,
        physical_neumann_out_w: solved.prescribed_neumann_out_w,
        physical_dirichlet_in_w: solved.physical_dirichlet_in_w,
        row,
    })
}

#[allow(clippy::too_many_arguments)]
fn advance_contact(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, controls: contact::Config,
    old: &[f64], dt_s: f64, config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
    model: &Model, radiation_config: fs_conduction::AmbientRadiationConfig,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    let coupling = ContactIterationConfig {
        max_iterations: controls.max_responses(), relaxation: 0.5,
        residual_tolerance_j: controls.residual_tolerance_j,
    };
    let solved = match (controls.coupled_newton, model) {
        (true, Model::Closed(model)) => contacts.advance_capacity_newton_with_enclosure_radiation(
            cx, engine, problem, old, dt_s, config,
            nonlinear.ok_or_else(|| bad("enclosure contact Newton lost its nonlinear controls"))?,
            coupling, controls.max_face_pairs, model, radiation_config,
        ),
        (false, Model::Closed(model)) => contacts.advance_capacity_with_enclosure_radiation(
            cx, engine, problem, old, dt_s, config, nonlinear, coupling, model, radiation_config,
        ),
        (true, Model::Open(model, ambient)) => contacts.advance_capacity_newton_with_open_enclosure_radiation(
            cx, engine, problem, old, dt_s, config,
            nonlinear.ok_or_else(|| bad("open enclosure contact Newton lost its nonlinear controls"))?,
            coupling, controls.max_face_pairs, model, *ambient, radiation_config,
        ),
        (false, Model::Open(model, ambient)) => contacts.advance_capacity_with_open_enclosure_radiation(
            cx, engine, problem, old, dt_s, config, nonlinear, coupling, model, *ambient, radiation_config,
        ),
    }.map_err(lower)?;
    poll(cx, deadline)?;
    let physical = &solved.physical;
    // Inspect the counted applied response before publishing fresh physical
    // heat. Authored Neumann and convection remain distinct from reflection.
    let applied = finite(solved.radiation.applied_heat_w.iter().sum())?;
    let decomposition = finite(
        solved.conduction.robin_out_w + solved.conduction.neumann_out_w
            - physical.convective_out_w - physical.neumann_out_w - applied,
    )?;
    if decomposition.abs() > radiation_config.balance_tolerance_w {
        return Err(bad("capacity contact/enclosure failed its original boundary decomposition allowance"));
    }
    let contact_evidence = contact::step_receipt(
        cx, contact::Config {
            residual_tolerance_j: solved.physical_residual_tolerance_j, ..controls
        }, physical, &solved.physical_contacts, solved.thermal_responses,
        solved.krylov_iterations, solved.nonlinear_updates, solved.nonlinear_backtracks,
        config.energy_tolerance_j, deadline,
    )?;
    let net_input_w = finite(physical.source_w + physical.dirichlet_in_w
        - physical.neumann_out_w - physical.convective_out_w - physical.radiation_out_w)?;
    let mut conduction = solved.conduction;
    contact::publish_physical(&mut conduction, physical, &solved.physical_contacts)?;
    let endpoint = radiation::EnclosureEndpoint {
        combined_boundary: physical.combined_boundary.clone(),
        convective_robin_fluxes: physical.convective_robin_fluxes.clone(),
        convective_out_w: physical.convective_out_w,
        report: solved.radiation,
        dt_s, time_integration: "backward-euler-capacity",
        physical_residual_norm_j: physical.residual_norm_j,
        physical_residual_tolerance_j: solved.physical_residual_tolerance_j,
        physical_energy_residual_j: physical.energy_residual_j,
        decomposition_residual_w: decomposition,
    };
    let row = step_row(
        cx, model, &endpoint, solved.thermal_responses, solved.nonlinear_backtracks,
        physical.neumann_out_w, physical.dirichlet_in_w, deadline,
    )?;
    Ok(Step {
        conduction, nonlinear: None, contact_evidence: Some(contact_evidence),
        solid_solves: solved.thermal_responses, nonlinear_updates: solved.nonlinear_updates,
        nonlinear_backtracks: solved.nonlinear_backtracks, endpoint, net_input_w,
        physical_neumann_out_w: physical.neumann_out_w,
        physical_dirichlet_in_w: physical.dirichlet_in_w, row,
    })
}

#[allow(clippy::too_many_arguments)]
fn step_row(
    cx: &Cx<'_>, model: &Model, endpoint: &radiation::EnclosureEndpoint,
    solid_solves: usize, nonlinear_backtracks: usize,
    prescribed_neumann_out_w: f64, physical_dirichlet_in_w: f64,
    deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    let report = &endpoint.report;
    let fresh = &report.radiosity;
    let count = model.surfaces().len();
    let rows = count + usize::from(model.is_open());
    if fresh.surface_temperatures_k.len() != rows || fresh.net_outward_heat_w.len() != rows
        || report.applied_heat_w.len() != count
    {
        return Err(bad("capacity enclosure evidence lost its ordered physical patches"));
    }
    let mut patches = Vec::with_capacity(count);
    for (index, surface) in model.surfaces().iter().enumerate() {
        poll(cx, deadline)?;
        patches.push(format!(
            "{{\"region\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{}}}",
            json_string(surface.name()), number(surface.area_m2())?,
            number(fresh.surface_temperatures_k[index])?,
            number(report.applied_heat_w[index])?, number(fresh.net_outward_heat_w[index])?,
        ));
    }
    let (name, bath) = match model {
        Model::Closed(_) => ("closed-gray-diffuse-enclosure", String::new()),
        Model::Open(_, ambient) => ("open-gray-diffuse-enclosure", format!(
            ",\"ambient_temperature_k\":{},\"ambient_absorbed_w\":{},\"enclosure_energy_closure_w\":{}",
            number(*ambient)?, number(-fresh.net_outward_heat_w[count])?, number(fresh.enclosure_energy_closure_w)?)),
    };
    Ok(format!(
        "{{\"model\":\"{name}\",\"storage\":\"volumetric-heat-capacity\",\"solid_solves\":{},\"radiation_trials\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"nonlinear_backtracks\":{},\"max_heat_mismatch_w\":{},\"applied_radiation_out_w\":{},\"nonlinear_radiation_out_w\":{},\"convective_out_w\":{},\"prescribed_neumann_out_w\":{},\"physical_dirichlet_in_w\":{},\"radiosity_residual_w_m2\":{},\"physical_residual_j\":{},\"physical_residual_threshold_j\":{},\"physical_energy_residual_j\":{},\"decomposition_residual_w\":{},\"patches\":[{}]{bath}}}",
        solid_solves, report.iterations, report.solid_iterations, report.krylov_iterations,
        nonlinear_backtracks, number(report.max_heat_mismatch_w)?,
        number(finite(report.applied_heat_w.iter().sum())?)?,
        number(model.surface_heat_w(fresh)?)?, number(endpoint.convective_out_w)?,
        number(prescribed_neumann_out_w)?, number(physical_dirichlet_in_w)?,
        number(fresh.linear_residual_max_w_m2)?, number(endpoint.physical_residual_norm_j)?,
        number(endpoint.physical_residual_tolerance_j)?, number(endpoint.physical_energy_residual_j)?,
        number(endpoint.decomposition_residual_w)?, patches.join(","),
    ))
}

/// The accepted temperature is unchanged, so its exact patch means and
/// radiosity remain the counted enclosure result. Refresh only the complete
/// physical boundary and original convection/reactions audited by its caller.
#[allow(clippy::too_many_arguments)]
pub(super) fn refresh(
    cx: &Cx<'_>, endpoint: &mut radiation::EnclosureEndpoint,
    model: &Model, physical: &EndpointEvaluation,
    physical_residual_tolerance_j: f64, solid_solves: usize, nonlinear_backtracks: usize,
    deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    poll(cx, deadline)?;
    endpoint.combined_boundary = physical.combined_boundary.clone();
    endpoint.convective_robin_fluxes = physical.convective_robin_fluxes.clone();
    endpoint.convective_out_w = physical.convective_out_w;
    endpoint.physical_residual_norm_j = physical.residual_norm_j;
    endpoint.physical_residual_tolerance_j = physical_residual_tolerance_j;
    endpoint.physical_energy_residual_j = physical.energy_residual_j;
    step_row(cx, model, endpoint, solid_solves, nonlinear_backtracks,
        physical.neumann_out_w, physical.dirichlet_in_w, deadline)
}

