//! Native time marching delegates reflected exchange to the shared enthalpy owner.
//! No iteration updates history; only the accepted complete endpoint advances it.
use fs_conduction::radiation::GrayDiffuseEnclosure;

use super::*;

pub(super) struct Step {
    pub(super) conduction: EnthalpyStepSolution,
    pub(super) endpoint: radiation::EnclosureEndpoint,
    pub(super) net_input_w: f64,
    pub(super) physical_neumann_out_w: f64,
    pub(super) physical_dirichlet_in_w: f64,
    pub(super) row: String,
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
    model: &GrayDiffuseEnclosure,
    radiation_config: fs_conduction::AmbientRadiationConfig,
    deadline: Option<(Instant, f64)>,
) -> Result<Step, SolveRefusal> {
    poll(cx, deadline)?;
    let solved = engine.advance_with_enclosure_radiation(
        cx, problem, interfaces, old, dt_s, config, model, radiation_config,
    ).map_err(failed)?;
    poll(cx, deadline)?;
    let report = &solved.radiation;
    let fresh = &report.radiosity;
    let applied = finite(report.applied_heat_w.iter().sum())?;
    let decomposition = finite(
        solved.conduction.robin_out_w + solved.conduction.neumann_out_w
            - solved.convective_out_w - solved.prescribed_neumann_out_w - applied,
    )?;
    // Applied exchange may be carried by either Robin or Neumann rows.
    // Only the original prescribed flux belongs to the external heat ledger.
    if decomposition.abs() > radiation_config.balance_tolerance_w {
        return Err(bad("transient enclosure failed its original boundary decomposition allowance"));
    }
    let net_input_w = finite(
        solved.conduction.source_w + solved.physical_dirichlet_in_w
            - solved.prescribed_neumann_out_w - solved.convective_out_w
            - fresh.enclosure_energy_closure_w,
    )?;
    let mut patches = Vec::with_capacity(model.surfaces().len());
    for (i, surface) in model.surfaces().iter().enumerate() {
        poll(cx, deadline)?;
        patches.push(format!(
            "{{\"region\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{}}}",
            json_string(surface.name()), number(surface.area_m2())?,
            number(fresh.surface_temperatures_k[i])?,
            number(report.applied_heat_w[i])?, number(fresh.net_outward_heat_w[i])?,
        ));
    }
    let row = format!(
        "{{\"model\":\"closed-gray-diffuse-enclosure\",\"solid_solves\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"max_heat_mismatch_w\":{},\"applied_radiation_out_w\":{},\"nonlinear_radiation_out_w\":{},\"convective_out_w\":{},\"prescribed_neumann_out_w\":{},\"physical_dirichlet_in_w\":{},\"radiosity_residual_w_m2\":{},\"physical_residual_j\":{},\"physical_residual_threshold_j\":{},\"physical_energy_residual_j\":{},\"decomposition_residual_w\":{},\"patches\":[{}]}}",
        report.iterations, report.solid_iterations, report.krylov_iterations,
        number(report.max_heat_mismatch_w)?, number(applied)?,
        number(fresh.enclosure_energy_closure_w)?, number(solved.convective_out_w)?,
        number(solved.prescribed_neumann_out_w)?, number(solved.physical_dirichlet_in_w)?,
        number(fresh.linear_residual_max_w_m2)?, number(solved.physical_residual_norm_j)?,
        number(solved.physical_residual_tolerance_j)?, number(solved.physical_energy_residual_j)?,
        number(decomposition)?, patches.join(","),
    );
    Ok(Step {
        conduction: solved.conduction,
        endpoint: radiation::EnclosureEndpoint {
            combined_boundary: solved.combined_boundary,
            convective_robin_fluxes: solved.convective_robin_fluxes,
            convective_out_w: solved.convective_out_w,
            report: solved.radiation,
            time_integration: "backward-euler-total-enthalpy",
            dt_s,
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

