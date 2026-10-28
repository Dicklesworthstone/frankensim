//! Native capacity trajectories use the shared reflected-radiation endpoint.
//! Returned inner solve data retain their applied boundary; fresh physical heat
//! is explicit so outer convection owners never advance or relabel history.

use fs_conduction::radiation::GrayDiffuseEnclosure;
use super::*;

pub(super) struct Step {
    pub(super) conduction: StepSolution,
    pub(super) nonlinear: Option<NonlinearStepSolution>,
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
    interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
    model: &GrayDiffuseEnclosure, radiation_config: fs_conduction::AmbientRadiationConfig,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    let solved = engine.advance_prescribed_with_enclosure_radiation(
        cx, problem, interfaces, old, dt_s, config, nonlinear, model, radiation_config,
    ).map_err(lower)?;
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
            - fresh.enclosure_energy_closure_w,
    )?;
    let mut patches = Vec::with_capacity(model.surfaces().len());
    for (index, surface) in model.surfaces().iter().enumerate() {
        poll(cx, deadline)?;
        patches.push(format!(
            "{{\"region\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{}}}",
            json_string(surface.name()), number(surface.area_m2())?,
            number(fresh.surface_temperatures_k[index])?,
            number(report.applied_heat_w[index])?, number(fresh.net_outward_heat_w[index])?,
        ));
    }
    let row = format!(
        "{{\"model\":\"closed-gray-diffuse-enclosure\",\"storage\":\"volumetric-heat-capacity\",\"solid_solves\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"nonlinear_backtracks\":{},\"max_heat_mismatch_w\":{},\"applied_radiation_out_w\":{},\"nonlinear_radiation_out_w\":{},\"convective_out_w\":{},\"prescribed_neumann_out_w\":{},\"physical_dirichlet_in_w\":{},\"radiosity_residual_w_m2\":{},\"physical_residual_j\":{},\"physical_residual_threshold_j\":{},\"physical_energy_residual_j\":{},\"decomposition_residual_w\":{},\"patches\":[{}]}}",
        report.iterations, report.solid_iterations, report.krylov_iterations,
        solved.nonlinear_backtracks, number(report.max_heat_mismatch_w)?, number(applied)?,
        number(fresh.enclosure_energy_closure_w)?, number(solved.convective_out_w)?,
        number(solved.prescribed_neumann_out_w)?, number(solved.physical_dirichlet_in_w)?,
        number(fresh.linear_residual_max_w_m2)?, number(solved.physical_residual_norm_j)?,
        number(solved.physical_residual_tolerance_j)?, number(solved.physical_energy_residual_j)?,
        number(decomposition)?, patches.join(","),
    );
    Ok(Step {
        conduction: solved.conduction,
        nonlinear: solved.nonlinear,
        nonlinear_backtracks: solved.nonlinear_backtracks,
        endpoint: radiation::EnclosureEndpoint {
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
        },
        net_input_w,
        physical_neumann_out_w: solved.prescribed_neumann_out_w,
        physical_dirichlet_in_w: solved.physical_dirichlet_in_w,
        row,
    })
}

