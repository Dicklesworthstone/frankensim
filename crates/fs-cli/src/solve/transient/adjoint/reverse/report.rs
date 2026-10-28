//! Deterministic derivative rows for the retained conduction receipt.
use super::{Cx, Derivatives, Instant, Objective, SolveRefusal, StepConfig, bad, finite, poll};
use super::super::super::{json_string, number};

pub(super) fn render(
    cx: &Cx<'_>, result: &Derivatives, ends: &[f64], objective: Objective,
    config: StepConfig, deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    let mut totals = vec![0.0; result.region_names.len()];
    let mut intervals = Vec::new();
    let mut start = 0.0;
    for (i, &end) in ends.iter().enumerate() {
        poll(cx, deadline)?;
        let mut rows = Vec::new();
        for (region, name) in result.region_names.iter().enumerate() {
            let derivative = result.source_controls[i][region];
            totals[region] = finite(totals[region] + derivative)?;
            rows.push(format!("{{\"region\":{},\"dtemperature_ddelivered_watts_k_per_w\":{}}}",
                json_string(name), number(derivative)?));
        }
        intervals.push(format!("{{\"step\":{},\"start_s\":{},\"end_s\":{},\"dtemperature_dsource_multiplier_k\":{},\"regional_power\":[{}]}}",
            i + 1, number(start)?, number(end)?, number(result.source_multipliers[i])?, rows.join(",")));
        start = end;
    }
    let regional = result.region_names.iter().zip(&totals).map(|(name, &value)| Ok(format!(
        "{{\"region\":{},\"dtemperature_duniform_added_watts_k_per_w\":{}}}", json_string(name), number(value)?)))
        .collect::<Result<Vec<_>, SolveRefusal>>()?;
    let robin = result.ports.iter().enumerate().map(|(i, name)| Ok(format!(
        "{{\"region\":{},\"dtemperature_dreference_k_per_k\":{},\"dtemperature_dlog_htc_k\":{}}}",
        json_string(name), number(result.references[i])?, number(result.coefficients[i])?)))
        .collect::<Result<Vec<_>, SolveRefusal>>()?;
    let mut initial = Vec::new();
    let mut initial_uniform = 0.0;
    for (i, &value) in result.initial.iter().enumerate() {
        if i % 256 == 0 { poll(cx, deadline)?; }
        initial.push(number(value)?);
        initial_uniform = finite(initial_uniform + value)?;
    }
    let selected = &result.selected;
    let time = if selected.state == 0 { 0.0 } else { ends[selected.state - 1] };
    let iteration_limit = ends.len().checked_mul(config.linear.max_iterations)
        .and_then(|v| v.checked_mul(result.max_dual_solves_per_step))
        .ok_or_else(|| bad("reverse iteration allowance overflows"))?;
    let reconstruction_limit = if result.revalidated_steps > 0
        || !result.radiation.is_empty() || !result.convection.is_empty() {
        0
    } else {
        iteration_limit
    };
    // Old request names retain the exact old scope and JSON bytes. The new
    // project-bound names explicitly select the extra material coordinates.
    let material_scope = objective.material_controls().then(|| SCOPE.replace(
        "Surface-power, prescribed-temperature, conductivity, contact, geometry and time-grid derivatives are not supplied.",
        MATERIAL_SCOPE,
    ));
    poll(cx, deadline)?;
    Ok(format!("{{\"schema\":\"fs-cli-native-transient-adjoint-v1\",\"output\":{},\"method\":\"discrete-backward-euler-reverse\",\"value_k\":{},\"time_s\":{},\"state_index\":{},\"active_vertex\":{},\"active_ties\":{},\"unique_sampled_maximum\":{},\"dtemperature_duniform_initial_k_per_k\":{},\"dtemperature_dinitial_temperatures\":[{}],\"dtemperature_dcapacity_multiplier_k\":{},\"regional_power\":[{}],\"robin\":[{}],\"intervals\":[{}],\"reconstructed_solid_endpoints\":{},\"revalidated_solid_endpoints\":{},\"reconstruction_krylov_iterations\":{},\"adjoint_krylov_iterations\":{},\"maximum_dual_relative_residual\":{},\"max_reverse_steps\":{},\"max_krylov_iterations_per_linear_solve\":{},\"max_reconstruction_krylov_iterations\":{},\"max_adjoint_krylov_iterations\":{},\"additional_storage_estimate_bytes\":{},\"authority\":\"Estimated\",\"scope\":{}{}{}{}}}",
        json_string(objective.name()), number(selected.value)?, number(time)?, selected.state,
        selected.vertex, selected.ties, selected.ties == 1, number(initial_uniform)?, initial.join(","),
        number(result.capacity)?, regional.join(","), robin.join(","), intervals.join(","), result.reconstructed_steps, result.revalidated_steps,
        result.primal_iterations, result.dual_iterations, number(result.worst_dual)?, ends.len(),
        config.linear.max_iterations, reconstruction_limit, iteration_limit,
        result.additional_storage_bytes, json_string(material_scope.as_deref().unwrap_or(SCOPE)),
        result.radiation, result.convection, result.material_controls))
}

const MATERIAL_SCOPE: &str = "Independent regional conductivity multipliers scale each region's full tensor/curve at one throughout the selected trajectory. Named matching-interface resistance multipliers scale every retained face resistance together at one, preserving mapped values, both temperature traces and consistent P1 integration. Both controls contract the same complete watt-coordinate multiplier without another timestep factor or solve. They retain distinct material/interface identities and have zero future influence at an initial selected peak. They are not absolute tensor-entry derivatives, new material claims or uncertainty propagation. Surface-power, prescribed-temperature, geometry, constitutive-knot and time-grid derivatives are not supplied.";

const SCOPE: &str = "Selected final or earliest sampled-peak nodal branch on the accepted fine time grid, including the initial state for sampled peaks. Ties have no unique maximum derivative. Future controls have zero influence on a selected earlier state. Complete storage and smooth endpoint k(T) tangent with fixed matching contact, geometry, prescribed temperatures and material laws, including declared natural-convection, streamwise quasi-steady airflow and ambient or closed-enclosure reflected-radiation mean feedback. Regional controls add absolute delivered watts through the same P1 projector, including zero-baseline regions; sum fine-interval derivatives to differentiate a declared workload segment. Robin reference/log(h) controls act on prescribed convection throughout. Natural reservoirs, air branches, independent fan banks and radiation surfaces have physical ambient, inlet, absolute speed, emissivity and reservoir controls. Closed enclosure surfaces retain their explicit view-factor ordering and shared reflected-radiation transpose; their controls are independent absolute emissivities, while original convection remains separate. The capacity control scales all regional capacities together. Surface-power, prescribed-temperature, conductivity, contact, geometry and time-grid derivatives are not supplied. Every counted reconstruction reproduces retained temperatures bit for bit. Coupled cooling endpoints instead revalidate the unchanged field under its original physical residual target and energy gate, with no primal reconstruction. No physical step, energy total or field is replaced. Natural ambient controls include signed wall and film-density coefficient dependence at fixed pressure. Air references include complete upstream and local wall feedback, including their radiative weighting. Each supported fan-bank speed includes hydraulic flow, air heat capacity and local forced-convection dependence; a refused bank retains its reason without discarding available controls. Not a continuous-time maximum, error-bound, uncertainty or physical-validation certificate.";
