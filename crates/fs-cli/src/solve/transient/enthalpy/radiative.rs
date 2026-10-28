//! One native boundary interface for both shared radiation owners.
//! Coupling loops keep original convection separate from internal exchange.
use fs_conduction::radiation::GrayDiffuseEnclosure;
use fs_conduction::transient::enthalpy::EnthalpyEndpointEvaluation;

use super::*;

pub(super) struct Step {
    pub(super) endpoint: EnthalpyStepSolution,
    pub(super) ambient: Option<RadiativeEndpoint>,
    pub(super) enclosure: Option<radiation::EnclosureEndpoint>,
    pub(super) evidence: String,
    pub(super) solid_solves: usize,
    pub(super) updates: usize,
    pub(super) krylov_iterations: usize,
    pub(super) net_input_w: f64,
    pub(super) residual_norm_j: f64,
    pub(super) residual_tolerance_j: f64,
    pub(super) energy_residual_j: f64,
    pub(super) dirichlet_in_w: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    old: &[f64],
    dt_s: f64,
    config: EnthalpyStepConfig,
    lowering: &radiation::LoweredRadiation,
    model: Option<&GrayDiffuseEnclosure>,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    if lowering.is_enclosure() != model.is_some() {
        return Err(bad("radiation boundary iteration lost its exact enclosure model"));
    }
    if let Some(model) = model {
        let solved = enclosure::advance(
            cx, engine, problem, interfaces, old, dt_s, config,
            model, lowering.config, deadline,
        )?;
        Ok(Step {
            solid_solves: solved.endpoint.report.iterations,
            updates: solved.endpoint.report.solid_iterations,
            krylov_iterations: solved.endpoint.report.krylov_iterations,
            net_input_w: solved.net_input_w,
            residual_norm_j: solved.endpoint.physical_residual_norm_j,
            residual_tolerance_j: solved.endpoint.physical_residual_tolerance_j,
            energy_residual_j: solved.endpoint.physical_energy_residual_j,
            dirichlet_in_w: solved.physical_dirichlet_in_w,
            endpoint: solved.conduction,
            ambient: None,
            enclosure: Some(solved.endpoint),
            evidence: solved.row,
        })
    } else {
        let solved = engine.advance_with_ambient_radiation(
            cx, problem, interfaces, old, dt_s, config,
            &lowering.patches, lowering.config,
        ).map_err(failed)?;
        poll(cx, deadline)?;
        let evidence = radiation_step_row(&solved)?;
        Ok(Step {
            solid_solves: solved.radiation.iterations,
            updates: solved.radiation.solid_iterations,
            krylov_iterations: solved.radiation.krylov_iterations,
            net_input_w: finite(solved.conduction.source_w + solved.physical_dirichlet_in_w
                - solved.conduction.neumann_out_w - solved.convective_out_w
                - solved.radiation.nonlinear_radiation_out_w)?,
            residual_norm_j: solved.physical_residual_norm_j,
            residual_tolerance_j: solved.physical_residual_tolerance_j,
            energy_residual_j: solved.physical_energy_residual_j,
            dirichlet_in_w: solved.physical_dirichlet_in_w,
            endpoint: solved.conduction,
            ambient: Some(RadiativeEndpoint {
                combined_boundary: solved.combined_boundary,
                convective_robin_fluxes: solved.convective_robin_fluxes,
                convective_out_w: solved.convective_out_w,
                report: solved.radiation,
            }),
            enclosure: None,
            evidence,
        })
    }
}

/// Re-evaluate the unchanged candidate with the outer owner's actual physical
/// boundary. No nonlinear solve or update to accepted history occurs here.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    old: &[f64],
    candidate: &[f64],
    dt_s: f64,
    config: EnthalpyStepConfig,
    radiation: Option<&radiation::LoweredRadiation>,
    model: Option<&GrayDiffuseEnclosure>,
) -> Result<EnthalpyEndpointEvaluation, SolveRefusal> {
    match (radiation, model) {
        (Some(lowering), Some(model)) if lowering.is_enclosure() => {
            engine.evaluate_endpoint_with_enclosure(
                cx, problem, interfaces, old, candidate, dt_s, config,
                model, lowering.config,
            ).map_err(failed)
        }
        (Some(lowering), None) if !lowering.is_enclosure() => {
            engine.evaluate_endpoint(
                cx, problem, interfaces, old, candidate, dt_s, config,
                Some(&lowering.patches),
            ).map_err(failed)
        }
        (None, None) => engine.evaluate_endpoint(
            cx, problem, interfaces, old, candidate, dt_s, config, None,
        ).map_err(failed),
        _ => Err(bad("physical boundary audit lost its declared radiation model")),
    }
}

pub(super) fn fluxes<'a>(
    endpoint: &'a EnthalpyStepSolution,
    ambient: Option<&'a RadiativeEndpoint>,
    enclosure: Option<&'a radiation::EnclosureEndpoint>,
) -> &'a [fs_conduction::RobinFlux] {
    if let Some(enclosure) = enclosure {
        &enclosure.convective_robin_fluxes
    } else if let Some(ambient) = ambient {
        &ambient.convective_robin_fluxes
    } else {
        &endpoint.robin_fluxes
    }
}

pub(super) fn convective_out_w(
    endpoint: &EnthalpyStepSolution,
    ambient: Option<&RadiativeEndpoint>,
    enclosure: Option<&radiation::EnclosureEndpoint>,
) -> f64 {
    enclosure.map_or_else(
        || ambient.map_or(endpoint.robin_out_w, |row| row.convective_out_w),
        |row| row.convective_out_w,
    )
}
