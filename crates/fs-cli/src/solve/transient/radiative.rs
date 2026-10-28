//! Select the physical radiation owner once; all capacity coupling loops share
//! its solve, work evidence, and fresh endpoint audit.

use fs_conduction::radiation::GrayDiffuseEnclosure;
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
    pub(super) radiation: RadiativeEndpoint,
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
    interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
    lowering: &radiation::LoweredRadiation, model: Option<&GrayDiffuseEnclosure>,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    if lowering.is_enclosure() != model.is_some() {
        return Err(bad("capacity radiation iteration lost its exact enclosure model"));
    }
    if let Some(model) = model {
        let solved = enclosure::advance(
            cx, engine, problem, interfaces, old, dt_s, config, nonlinear,
            model, lowering.config, deadline,
        )?;
        let nonlinear_evidence = solved.nonlinear.as_ref().map(nonlinear_row)
            .transpose()?.unwrap_or_else(|| "null".to_string());
        let report = &solved.endpoint.report;
        let physical_energy = PhysicalEnergy {
            neumann_out_w: solved.physical_neumann_out_w,
            robin_out_w: solved.endpoint.convective_out_w
                + report.radiosity.enclosure_energy_closure_w,
            dirichlet_in_w: solved.physical_dirichlet_in_w,
            net_input_w: solved.net_input_w,
            energy_residual_j: solved.endpoint.physical_energy_residual_j,
        };
        return Ok(Step {
            conduction: solved.conduction,
            nonlinear_evidence,
            radiation_evidence: solved.row,
            solid_solves: report.iterations,
            nonlinear_updates: if nonlinear.is_some() { report.solid_iterations } else { 0 },
            nonlinear_backtracks: solved.nonlinear_backtracks,
            krylov_iterations: report.krylov_iterations,
            physical_residual_tolerance_j: solved.endpoint.physical_residual_tolerance_j,
            physical_energy,
            radiation: RadiativeEndpoint::Enclosure(solved.endpoint),
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
    let radiation_evidence = radiation_row(&solved)?;
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
            report: solved.radiation,
        },
    })
}

/// No primal step is taken here. Candidate and old temperature retain their
/// identities while the actual convection and reflected exchange are audited.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>, old: &[f64], candidate: &[f64],
    dt_s: f64, config: StepConfig, lowering: Option<&radiation::LoweredRadiation>,
    model: Option<&GrayDiffuseEnclosure>, deadline: Option<(Instant, f64)>,
) -> Result<EndpointEvaluation, SolveRefusal> {
    poll(cx, deadline)?;
    if lowering.is_some_and(radiation::LoweredRadiation::is_enclosure) != model.is_some() {
        return Err(bad("capacity physical boundary audit lost its declared radiation model"));
    }
    let physical = if let Some(model) = model {
        let lowering = lowering.ok_or_else(|| bad("enclosure endpoint lost its radiation policy"))?;
        engine.evaluate_prescribed_enclosure_endpoint(
            cx, problem, interfaces, old, candidate, dt_s, config, model, lowering.config,
        ).map_err(lower)?.physical
    } else {
        engine.evaluate_prescribed_endpoint(
            cx, problem, interfaces, old, candidate, dt_s, config,
            lowering.map(|row| row.patches.as_slice()),
        ).map_err(lower)?
    };
    poll(cx, deadline)?;
    Ok(physical)
}

