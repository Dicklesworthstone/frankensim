//! Opt-in nominal derivatives on the FINAL native conduction mesh and field.
//!
//! A report output named `temperature-max-adjoint` selects the smallest-index
//! hottest vertex of the declared temperature-max region. At a tie this is a
//! selected nodal functional, NOT a unique derivative of a relocating maximum.
//! It is useful as a frozen mean control even when an active branch changes.
//! `temperature-max-contact-adjoint` selects the same complete derivative and
//! adds named contact-resistance controls; choose one report, not both.
//! `temperature-max-boundary-adjoint` instead adds prescribed-temperature
//! controls through the complete lift, without changing legacy report bytes.
//! `temperature-max-contact-boundary-adjoint` requests BOTH control families
//! from that same solve/dual; it does not request two independent analyses.
//! `temperature-volume-mean-adjoint` instead differentiates the spatial volume
//! mean over that same declared requirement region, including both control
//! families. The maximum requirement and its uncertainty budget are unchanged.
//!
//! Reuse production adjoints, contact operators and complete boundary feedback.
//! Only contractions of native input laws live here; no perturbed primal or
//! new solver. Nonlinear material and radiation derivatives are NOT frozen.
//! Independent series/parallel fan-bank controls compose hydraulic, capacity,
//! convection and thermal feedback. Regional material-law scales reuse that
//! same complete dual. Unsupported or nonsmooth controls remain explicit rows.

use std::collections::BTreeMap;
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc};
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalAnalyzer, RobinFeedbackAnalysisConfig};
use fs_exec::Cx;
use super::{EvidenceWork, ProjectSpec, RungSolved, SolveRefusal, canonical_f64,
    conduction_error, json_string, temperature_maximum_region, trace_qoi_region_vertices};

mod air_inlet;
mod contractions;
mod contact_controls;
mod coupled_thermal;
mod fan_speed;
mod material_controls;
mod natural_feedback;
mod nonlinear_solid;
mod prescribed_controls;
mod radiative_feedback;
mod surface_power;
mod volume_mean;

const OUTPUT: &str = "temperature-max-adjoint";
const COMBINED_OUTPUT: &str = "temperature-max-contact-boundary-adjoint";
const MAX_PARAMETERS: usize = 256;
const SCOPE: &str = "Estimated derivative of a selected hottest nodal temperature on the final accepted native mesh. Fixed geometry and matching contact; the hydraulic operating point is fixed except for the explicitly admitted fan-speed control. Power rows differentiate declared pre-duty watts: volume sources retain regional nodal mixing, while surface sources retain their inward P1 face load and actual patch-area normalization. Smooth heterogeneous k(T) uses the full nonsymmetric material Jacobian; material slope discontinuities and validity endpoints refuse. Constant-conductivity solids and nonradiating linear solid/air models retain their existing linear analysis. Natural-convection and radiation modes include smooth k(T) and complete area-mean feedback. Coupled nonlinear-solid/air and radiation/air modes retain the full material, consistent radiative secant, weighted reference and stream-wise air Jacobian together. Only original convective heat enters the air law. Every nonlinear mode checks the complete constitutive primal residual at the unchanged field, not just the last frozen outer iterate, and supplies no inverse or gradient-error certificate. Emissivity is fixed at the selected card query; its row is a local coefficient partial, not a card-selection or temperature-dependent-emissivity derivative. Not a unique maximum derivative at a tie, a continuum/shape derivative, experimental validation, or a parameter-uncertainty bound. Independent air-inlet derivatives include upstream segments at fixed flow and transport properties. Independent single/series/parallel fan-bank speed derivatives use native quadratic vent/leakage losses, each member's hydraulic response and fan affinity, include all branch capacity and smooth-card convection changes, and differentiate each absolute speed ratio, not its logarithm or a common system speed. This is a nominal local model derivative, not a derivative of interval root-bracket endpoints or pressure-tolerance uncertainty. Conductivity-multiplier rows differentiate a shared dimensionless scale s on each region's whole effective conductivity tensor/curve K(T), evaluated at s=1; they retain the complete thermal dual and prescribed-temperature lift. They are not absolute scalar-conductivity, tensor-entry, card-selection, heat-capacity or material-uncertainty derivatives, and do not modify or validate the selected material claim. Nonsmooth fan-curve knots, nonunique parallel fan inverses, card regime boundaries, natural convection combined with airflow and Dirichlet-temperature derivatives are not supplied.";

fn bad(message: impl Into<String>) -> SolveRefusal {
    conduction_error("cli-solve-nominal-adjoint", message,
        "request a report named temperature-max-adjoint for an admitted thermal solve; inspect unsupported parameter rows")
}
fn poll(cx: &Cx<'_>) -> Result<(), SolveRefusal> {
    cx.checkpoint().map_err(|_| conduction_error("cli-solve-cancelled",
        "native nominal adjoint was cancelled", "resume the accepted pipeline prefix"))
}
fn number(value: f64) -> Result<String, SolveRefusal> {
    canonical_f64(value).ok_or_else(|| bad("nonfinite native adjoint or contraction"))
}
fn finite(value: f64) -> Result<f64, SolveRefusal> {
    if value.is_finite() { Ok(value) } else { Err(bad("nonfinite native adjoint arithmetic")) }
}
fn lower(error: fs_conduction::ConductionError) -> SolveRefusal {
    match error {
        fs_conduction::ConductionError::Cancelled { .. } => conduction_error(
            "cli-solve-cancelled", "nominal adjoint interrupted", "resume the accepted pipeline prefix"),
        error => bad(format!("production nominal adjoint refused: {error}")),
    }
}
fn zeros(n: usize) -> Result<Vec<f64>, SolveRefusal> {
    let mut out = Vec::new();
    out.try_reserve_exact(n).map_err(|_| bad("nominal adjoint allocation refused"))?;
    out.resize(n, 0.0);
    Ok(out)
}

/// An unsupported state law invalidates ALL contractions, not only its own
/// parameter row. Never fall back to a frozen Robin power derivative.
fn admit_state_laws(spec: &ProjectSpec) -> Result<(), SolveRefusal> {
    let Some(setup) = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref()) else { return Ok(()); };
    let mut natural = false;
    let mut airflow = false;
    for row in &setup.boundaries {
        match &row.condition {
            fs_project::ThermalBoundaryCondition::NaturalConvection { correlation, .. } => {
                natural = true;
                if correlation != natural_feedback::CARD {
                    return Err(bad("natural convection requires a differentiated Churchill-Chu card; frozen Robin derivatives are not admitted"));
                }
            }
            fs_project::ThermalBoundaryCondition::AirflowConvection { .. } => airflow = true,
            _ => {}
        }
    }
    if natural && airflow { return Err(bad("a combined natural/airflow adjoint is not supplied")); }
    Ok(())
}

pub(super) fn requested(spec: &ProjectSpec) -> Result<bool, SolveRefusal> {
    let rows = spec.outputs.as_deref().unwrap_or(&[]);
    let mut found = false;
    for row in rows.iter().filter(|row| row.name == OUTPUT || row.name == contact_controls::OUTPUT
        || row.name == prescribed_controls::OUTPUT || row.name == COMBINED_OUTPUT
        || row.name == volume_mean::OUTPUT) {
        if found || row.kind != "report" {
            return Err(bad("choose exactly one report: temperature-max-adjoint, temperature-max-contact-adjoint, temperature-max-boundary-adjoint, temperature-max-contact-boundary-adjoint or temperature-volume-mean-adjoint"));
        }
        if row.name == volume_mean::OUTPUT && crate::SOLVE_DRIVER_VERSION < 48 {
            return Err(bad("the volume-mean adjoint requires solve driver 48 so older cached reports cannot substitute for this objective"));
        }
        if row.name == volume_mean::OUTPUT && row.region.is_some() {
            return Err(bad("the volume-mean adjoint uses the existing temperature-max requirement region; :region on an output is a surface selector, not a volume selector"));
        }
        found = true;
    }
    if found && temperature_maximum_region(spec).is_none() {
        return Err(bad("the adjoint needs the existing declared temperature-max requirement and region"));
    }
    if found { admit_state_laws(spec)?; }
    Ok(found)
}

fn row(target: &str, entity: &str, ordinal: usize, unit: &str, value: f64) -> Result<String, SolveRefusal> {
    Ok(format!("{{\"target\":{},\"entity\":{},\"ordinal\":{},\"parameter_unit\":{},\"derivative_unit\":{},\"derivative\":{}}}",
        json_string(target), json_string(entity), ordinal, json_string(unit),
        json_string(&format!("K/({unit})")), number(value)?))
}
fn unsupported(target: &str, entity: &str, reason: &str) -> String {
    format!("{{\"target\":{},\"entity\":{},\"reason\":{}}}",
        json_string(target), json_string(entity), json_string(reason))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn extract(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved,
    ids: &BTreeMap<String, u32>, audited: &fs_mesh::AuditedLabeledTetComplex,
    work: EvidenceWork<'_>,
) -> Result<String, SolveRefusal> {
    poll(cx)?;
    admit_state_laws(spec)?;
    let combined = spec.outputs.as_deref().unwrap_or(&[]).iter()
        .any(|row| row.name == COMBINED_OUTPUT);
    let volume_mean_requested = spec.outputs.as_deref().unwrap_or(&[]).iter()
        .any(|row| row.name == volume_mean::OUTPUT);
    let contact_requested = volume_mean_requested || combined || contact_controls::requested(spec);
    let boundary_requested = volume_mean_requested || combined || prescribed_controls::requested(spec);
    let output = if volume_mean_requested { volume_mean::OUTPUT }
        else if combined { COMBINED_OUTPUT }
        else if contact_requested { contact_controls::OUTPUT }
        else if boundary_requested { prescribed_controls::OUTPUT } else { OUTPUT };
    let region = temperature_maximum_region(spec).ok_or_else(|| bad("missing maximum region"))?;
    let region_id = *ids.get(region).ok_or_else(|| bad("maximum region has no mesh label"))?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("final operator was not retained"))?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing conduction setup"))?;
    if data.radiating_boundary.is_some() != setup.radiation.is_some() {
        return Err(bad("radiation intent differs from the retained native operator"));
    }
    let n = solved.mesh.vertex_count();
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let count = |divisor| usize::try_from(memory / divisor).unwrap_or(usize::MAX);
    let fan_count = spec.cooling.as_ref().and_then(|c| c.fan_system.as_ref())
        .map_or(0, |system| system.banks.len());
    let radiation_count = setup.radiation.as_ref().map_or(0, |r| r.surfaces.len());
    let material_count = spec.materials.as_deref().unwrap_or(&[]).len();
    let contact_count = if contact_requested { spec.interface_cards.as_deref().unwrap_or(&[]).len() } else { 0 };
    let parameter_count = spec.power.as_deref().unwrap_or(&[]).len()
        .checked_add(setup.boundaries.len().saturating_mul(2))
        .and_then(|n| n.checked_add(radiation_count.saturating_mul(2)))
        .and_then(|n| n.checked_add(fan_count))
        .and_then(|n| n.checked_add(material_count))
        .and_then(|n| n.checked_add(contact_count))
        .ok_or_else(|| bad("adjoint parameter count overflow"))?;
    // Logical extra vectors/records, not a total-allocator/RSS promise. The
    // numerical owners independently enforce their matrix and iteration caps.
    if n > count(128) || solved.mesh.element_count() > count(64)
        || parameter_count > MAX_PARAMETERS {
        return Err(bad("nominal adjoint exceeds declared memory or 256-parameter work envelope"));
    }
    let temperature = &solved.solution.temperature;
    let (weights, functional_fields) = if volume_mean_requested {
        let goal = volume_mean::prepare(cx, &solved.mesh, &solved.labels, region_id, temperature)?;
        let fields = format!(concat!("\"functional\":\"region-volume-mean-temperature\",\"region\":{},",
            "\"selected_vertex\":null,\"value_k\":{},\"tied_maximum_vertices\":null,\"runner_up_gap_k\":null,",
            "\"region_volume_m3\":{},\"region_elements\":{},\"region_label\":{},"),
            json_string(region), number(goal.value_k)?, number(goal.volume_m3)?, goal.elements, region_id);
        (goal.weights, fields)
    } else {
        let (vertices, _) = trace_qoi_region_vertices(&solved.labels, &solved.mesh.complex().tets,
            n, region_id, work).map_err(|error| {
                if work.is_requested() { conduction_error("cli-solve-cancelled",
                    "nominal region trace interrupted", "resume the accepted pipeline prefix") }
                else { bad(format!("nominal region trace refused: {error:?}")) }
            })?;
        let mut selected = *vertices.first().ok_or_else(|| bad("empty maximum region"))?;
        for (i, &v) in vertices.iter().enumerate() {
            if i % 512 == 0 { poll(cx)?; }
            finite(temperature[v])?;
            if temperature[v] > temperature[selected]
                || (temperature[v] == temperature[selected] && v < selected) { selected = v; }
        }
        let mut tied = 0;
        let mut second = None::<f64>;
        for (i, &v) in vertices.iter().enumerate() {
            if i % 512 == 0 { poll(cx)?; }
            if temperature[v] == temperature[selected] { tied += 1; }
            if v != selected { second = Some(second.map_or(temperature[v], |x| x.max(temperature[v]))); }
        }
        let mut weights = zeros(n)?;
        weights[selected] = 1.0;
        let gap = second.map(|value| number(temperature[selected] - value)).transpose()?
            .unwrap_or_else(|| "null".into());
        let fields = format!(concat!("\"functional\":\"selected-nodal-temperature\",\"region\":{},\"selected_vertex\":{},",
            "\"value_k\":{},\"tied_maximum_vertices\":{},\"runner_up_gap_k\":{},"),
            json_string(region), selected, number(temperature[selected])?, tied, gap);
        (weights, fields)
    };
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &data.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let config = LinearGoalAnalysisConfig {
        residual_limits: fs_solver::goal::GoalResidualLimits { max_rows: count(256), max_nonzeros: count(64) },
        max_stability_iterations: data.linear.max_iterations,
    };
    let has_natural = setup.boundaries.iter().any(|b| matches!(
        b.condition, fs_project::ThermalBoundaryCondition::NaturalConvection { .. }));
    let mut natural_ambient = BTreeMap::new();
    let mut radiation_rows = Vec::new();
    let mut inlet_rows = Vec::new();
    let (lambda, residual, dual_iterations, stability_iterations, response_iterations, mode) =
    if !data.air_paths.is_empty()
        && (setup.radiation.is_some() || nonlinear_solid::needed(cx, problem)?) {
        let result = coupled_thermal::pullback(cx, spec, solved, &weights)?;
        inlet_rows = result.air_rows;
        radiation_rows = result.radiation_rows;
        let mode = if setup.radiation.is_some() { "radiation-full-air-feedback" }
            else { "nonlinear-solid-full-air-feedback" };
        (result.gradient.nodal_load, result.gradient.relative_residual, result.gradient.iterations,
            None, None, mode)
    } else if setup.radiation.is_some() {
        let result = radiative_feedback::pullback(cx, spec, solved, &weights)?;
        natural_ambient = result.natural_ambient;
        radiation_rows = result.radiation_rows;
        (result.gradient.nodal_load, result.gradient.relative_residual, result.gradient.iterations,
            None, None, "radiation-full-wall-feedback")
    } else if has_natural {
        let (gradient, ambient) = natural_feedback::pullback(cx, spec, solved, &weights)?;
        natural_ambient = ambient;
        (gradient.nodal_load, gradient.relative_residual, gradient.iterations,
            None, None, "natural-convection-full-wall-feedback")
    } else if data.air_paths.is_empty() && nonlinear_solid::needed(cx, problem)? {
        let gradient = nonlinear_solid::pullback(cx, problem, data.interfaces.as_ref(),
            data.linear, temperature, &weights)?;
        (gradient.nodal_load, gradient.relative_residual, gradient.iterations,
            None, None, "nonlinear-solid-full-material-feedback")
    } else if data.air_paths.is_empty() {
        let analyzer = LinearGoalAnalyzer::new(cx, problem, data.interfaces.as_ref(), data.linear,
            temperature, &weights, config).map_err(lower)?;
        let analysis = analyzer.analyze(cx, temperature).map_err(lower)?;
        if !analysis.dual_relative_residual.is_finite()
            || analysis.dual_relative_residual >= data.linear.tolerance {
            return Err(bad("nominal dual exhausted its work before meeting the true-residual gate"));
        }
        let mut lambda = zeros(n)?;
        for (i, (&vertex, &value)) in analyzer.dofs().free().iter().zip(analyzer.free_dual()).enumerate() {
            if i % 512 == 0 { poll(cx)?; }
            lambda[vertex] = finite(value)?;
        }
        (lambda, analysis.dual_relative_residual, analysis.dual_iterations,
            Some(analysis.stability_iterations), Some(0), "linear-solid")
    } else {
        let feedback = RobinFeedbackAnalysisConfig {
            residual: fs_solver::goal::feedback::FeedbackResidualLimits {
                solid: config.residual_limits, max_ports: 64,
                max_transfer_nonzeros: count(128), max_response_entries: count(64),
                max_verification_entries: count(8),
            },
            max_response_iterations: data.linear.max_iterations, max_lowering_entries: count(16),
        };
        let analyzer = fs_airflow::conjugate::goal::maximum::prepare_linear_maximum(cx,
            problem, data.interfaces.as_ref(), &data.air_paths, data.linear, temperature, config, feedback)
            .map_err(|error| {
                if work.is_requested() { conduction_error("cli-solve-cancelled",
                    "nominal coupled adjoint interrupted", "resume the accepted pipeline prefix") }
                else { bad(format!("complete solid/air adjoint preparation refused: {error}")) }
            })?;
        // This owner refuses missing complete-system inverse evidence. Never
        // replace a failed coupled derivative with a frozen-solid derivative.
        // The coupled pullback caps restart storage at 32 vectors; narrowing
        // restart does not increase the shared iteration budget or tolerance.
        let linear = fs_conduction::LinearConfig { restart: data.linear.restart.clamp(1, 32), ..data.linear };
        let gradient = analyzer.pullback_affine_controls(cx, temperature, &weights, linear)
            .map_err(lower)?;
        let port_names: Vec<_> = analyzer.ports().iter().map(|port| port.name.as_str()).collect();
        inlet_rows = air_inlet::rows(cx, setup, &data.air_paths, &port_names, &gradient.references)?;
        (gradient.nodal_load, gradient.relative_residual, gradient.iterations,
            None, Some(analyzer.response_iterations()), "linear-solid-full-air-feedback")
    };
    let source_bar = regional_source_pullback(cx, solved, &lambda)?;
    let volumes: BTreeMap<u32, f64> = audited.witness().per_region_auditor.iter()
        .map(|(region, volume)| (region.0, *volume)).collect();
    // Use the SAME explicit Surface/Region classification as the primal.
    // A surface has no volume ID; it is a Neumann load on retained face slots.
    let surface_sources = super::surface_heat(spec)?;
    let mut rows = radiation_rows;
    rows.extend(inlet_rows);
    rows.extend(material_controls::rows(cx, spec, solved, ids, &lambda, count(64))?);
    if contact_requested {
        rows.extend(contact_controls::rows(cx, spec, solved, &lambda, count(64))?);
    }
    let mut missing = Vec::new();
    for (ordinal, power) in spec.power.as_deref().unwrap_or(&[]).iter().enumerate() {
        poll(cx)?;
        let value = if surface_sources.names.contains(&power.region) {
            surface_power::pullback(cx, solved, &power.region, power.duty, &lambda)?
        } else {
            let id = *ids.get(&power.region).ok_or_else(|| bad("power region has no label"))?;
            let volume = *volumes.get(&id).ok_or_else(|| bad("power region has no audited volume"))?;
            if !(volume.is_finite() && volume > 0.0) { return Err(bad("invalid audited source volume")); }
            finite(*source_bar.get(&id).ok_or_else(|| bad("source region has no element contribution"))? * (power.duty / volume))?
        };
        rows.push(row("power", &power.region, ordinal, "W", value)?);
    }
    let boundary_bar = boundary_pullback(cx, solved, &lambda)?;
    for (ordinal, declared) in setup.boundaries.iter().enumerate() {
        poll(cx)?;
        use fs_project::spec::ThermalBoundaryCondition as B;
        let value = boundary_bar.get(&declared.target).ok_or_else(|| bad("boundary derivative has no native trace"))?;
        match &declared.condition {
            B::Convection { .. } => {
                rows.push(row("convection-coefficient", &declared.target, ordinal, "W/m^2/K", value[0])?);
                rows.push(row("convection-temperature", &declared.target, ordinal, "K", value[1])?);
            }
            B::HeatFlux { .. } => rows.push(row("heat-flux", &declared.target, ordinal, "W/m^2", value[2])?),
            B::FixedTemperature { .. } => {
                if !boundary_requested { missing.push(unsupported("fixed-temperature", &declared.target,
                    "prescribed values require their complete lift derivative")); }
            }
            B::NaturalConvection { .. } => rows.push(row("natural-convection-ambient", &declared.target,
                ordinal, "K", *natural_ambient.get(&declared.target)
                    .ok_or_else(|| bad("natural boundary has no complete ambient derivative"))?)?),
            // Already emitted once per branch, not once per segment.
            B::AirflowConvection { .. } => {}
        }
    }
    fan_speed::append(cx, spec, solved, &lambda, &mut rows, &mut missing)?;
    if boundary_requested {
        prescribed_controls::append(cx, spec, solved, &weights, &lambda, &mut rows, &mut missing, count(32))?;
    }
    let scope = if volume_mean_requested { format!("{} {} {}",
            SCOPE.replace("Estimated derivative of a selected hottest nodal temperature on the final accepted native mesh.", volume_mean::SCOPE)
                .replace("Not a unique maximum derivative at a tie, a continuum/shape derivative", "Not a continuum/shape derivative")
                .replace(" and Dirichlet-temperature derivatives", ""),
            contact_controls::SCOPE,
            prescribed_controls::SCOPE.replace("direct selected-node objective term", "direct volume-mean objective weights on prescribed nodes")) }
        else if combined { format!("{} {} {}",
            SCOPE.replace(" and Dirichlet-temperature derivatives", ""),
            contact_controls::SCOPE, prescribed_controls::SCOPE) }
        else if contact_requested { format!("{SCOPE} {}", contact_controls::SCOPE) }
        else if boundary_requested { format!("{} {}",
            SCOPE.replace(" and Dirichlet-temperature derivatives", ""), prescribed_controls::SCOPE) }
        else { SCOPE.to_string() };
    poll(cx)?;
    Ok(format!(concat!("{{\"schema\":\"frankensim.cli.nominal-adjoint.v1\",\"output\":{},",
        "{}\"mode\":{},",
        "\"true_relative_residual\":{},\"dual_iterations\":{},\"stability_iterations\":{},\"response_iterations\":{},",
        "\"parameters\":[{}],\"unsupported\":[{}],\"authority\":\"Estimated\",\"scope\":{}}}"),
        json_string(output), functional_fields, json_string(mode),
        number(residual)?, dual_iterations,
        stability_iterations.map_or_else(|| "null".into(), |n| n.to_string()),
        response_iterations.map_or_else(|| "null".into(), |n| n.to_string()),
        rows.join(","), missing.join(","), json_string(&scope)))
}

/// Transpose BOTH the consistent source mass and native regional nodal mixing.
/// Using only V/4*lambda is wrong at shared regional nodes. Normalization uses
/// the exact audited region volumes upstream, not a replacement mesh sum.
fn regional_source_pullback(cx: &Cx<'_>, solved: &RungSolved, lambda: &[f64])
    -> Result<BTreeMap<u32, f64>, SolveRefusal>
{
    let mesh = &solved.mesh;
    let mut mass_bar = zeros(mesh.vertex_count())?;
    let mut lumped = zeros(mesh.vertex_count())?;
    for (e, tet) in mesh.complex().tets.iter().enumerate() {
        if e % 512 == 0 { poll(cx)?; }
        let volume = finite(mesh.element_volume(e))?;
        let vertices = tet.map(|v| v as usize);
        let contribution = contractions::source_mass(volume, vertices.map(|v| lambda[v]));
        for (j, v) in vertices.into_iter().enumerate() {
            mass_bar[v] = finite(mass_bar[v] + contribution[j])?;
            lumped[v] = finite(lumped[v] + volume / 4.0)?;
        }
    }
    for v in 0..mass_bar.len() {
        if v % 512 == 0 { poll(cx)?; }
        mass_bar[v] = if lumped[v] > 0.0 { finite(mass_bar[v] / lumped[v])? }
            else { 0.0 };
    }
    let mut result = BTreeMap::new();
    for (e, tet) in mesh.complex().tets.iter().enumerate() {
        if e % 512 == 0 { poll(cx)?; }
        let weight = mesh.element_volume(e) / 4.0;
        let value = result.entry(solved.labels[e]).or_insert(0.0);
        for &v in tet { *value = finite(*value + weight * mass_bar[v as usize])?; }
    }
    Ok(result)
}

fn boundary_pullback(cx: &Cx<'_>, solved: &RungSolved, lambda: &[f64])
    -> Result<BTreeMap<String, [f64; 3]>, SolveRefusal>
{
    let mesh = &solved.mesh;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("missing boundary operator"))?;
    let mut result: BTreeMap<String, [f64; 3]> = data.boundary.region_names().iter()
        .map(|name| (name.clone(), [0.0; 3])).collect();
    for (slot, face) in mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let Some(region) = data.boundary.region_for(slot) else { continue; };
        let vertices = face.vertices.map(|v| v as usize);
        let [a, b, c] = vertices.map(|v| mesh.positions()[v]);
        let u = [b[0]-a[0], b[1]-a[1], b[2]-a[2]];
        let v = [c[0]-a[0], c[1]-a[1], c[2]-a[2]];
        let cross = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]];
        let area = finite(0.5 * fs_math::det::sqrt(cross.iter().map(|x| x*x).sum()))?;
        let (h, reference) = match &data.boundary.conditions()[region] {
            ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(t) } => (*h, *t),
            ThermalBc::Neumann { .. } => (0.0, 0.0),
            ThermalBc::Dirichlet { .. } => continue,
            _ => return Err(bad("native adjoint requires uniform boundary coefficients")),
        };
        let bars = if matches!(&data.boundary.conditions()[region], ThermalBc::Neumann { .. }) {
            // Unused Robin contractions must not manufacture an overflow for a
            // well-defined flux derivative.
            [0.0, 0.0, finite(-(area / 3.0) * vertices.iter().map(|&v| lambda[v]).sum::<f64>())?]
        } else {
            contractions::boundary(area, vertices.map(|v| lambda[v]),
                vertices.map(|v| solved.solution.temperature[v]), h, reference)
        };
        let total = result.get_mut(&data.boundary.region_names()[region]).expect("retained boundary name");
        for i in 0..3 { total[i] = finite(total[i] + bars[i])?; }
    }
    Ok(result)
}

#[cfg(test)]
mod admission_tests {
    use super::*;

    #[test]
    fn natural_primal_and_complete_nominal_adjoint_are_admitted() {
        let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"), "/../../examples/heatsink-fan/heatsink-natural.fsim"
        ))).expect("the native passive-heatsink fixture parses").decoded.spec;
        assert!(!requested(&spec).expect("a primal-only natural solve remains admitted"));
        admit_state_laws(&spec).expect("the native card has a complete differential");
        spec.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
            name: OUTPUT.into(), kind: "report".into(), region: None,
        });
        assert!(requested(&spec).expect("complete natural adjoint"));
        for boundary in &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
            if let fs_project::ThermalBoundaryCondition::NaturalConvection { correlation, .. } = &mut boundary.condition {
                *correlation = "convection.dittus-boelter".into();
            }
        }
        assert!(requested(&spec).is_err(), "an unknown state derivative must never fall back to frozen h");
    }

    #[test]
    fn linear_state_laws_still_admit_nominal_adjoints() {
        let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"), "/../../data/reference-project/cooling-reference.fsim"
        ))).expect("the native reference fixture parses").decoded.spec;
        admit_state_laws(&spec).expect("linear Robin state law");
        spec.outputs.get_or_insert_with(Vec::new).push(fs_project::spec::OutputRequest {
            name: OUTPUT.into(), kind: "report".into(), region: None,
        });
        assert!(requested(&spec).expect("linear adjoint remains admitted"));
    }

    #[test]
    fn contact_extended_request_is_explicit_and_cannot_duplicate_the_same_goal() {
        let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"), "/../../examples/contact-pair/contact-pair.fsim"
        ))).unwrap().decoded.spec;
        assert!(!requested(&spec).unwrap());
        assert!(!contact_controls::requested(&spec));
        let request = fs_project::spec::OutputRequest {
            name: contact_controls::OUTPUT.into(), kind: "report".into(), region: None,
        };
        spec.outputs.get_or_insert_with(Vec::new).push(request.clone());
        assert!(requested(&spec).unwrap());
        assert!(contact_controls::requested(&spec));
        let mut both = spec.clone();
        both.outputs.as_mut().unwrap().push(fs_project::spec::OutputRequest {
            name: OUTPUT.into(), ..request.clone()
        });
        assert!(requested(&both).is_err());
        let mut duplicated = spec.clone();
        duplicated.outputs.as_mut().unwrap().push(request);
        assert!(requested(&duplicated).is_err());
        spec.outputs.as_mut().unwrap().last_mut().unwrap().kind = "scalar".into();
        assert!(requested(&spec).is_err());
    }

    #[test]
    fn prescribed_extended_request_cannot_alias_another_adjoint_report() {
        let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"), "/../../examples/contact-pair/contact-pair.fsim"
        ))).unwrap().decoded.spec;
        let request = fs_project::spec::OutputRequest {
            name: prescribed_controls::OUTPUT.into(), kind: "report".into(), region: None,
        };
        spec.outputs.get_or_insert_with(Vec::new).push(request.clone());
        assert!(requested(&spec).unwrap());
        assert!(prescribed_controls::requested(&spec));
        for other in [OUTPUT, contact_controls::OUTPUT, prescribed_controls::OUTPUT] {
            let mut duplicate = spec.clone();
            duplicate.outputs.as_mut().unwrap().push(fs_project::spec::OutputRequest {
                name: other.into(), ..request.clone()
            });
            assert!(requested(&duplicate).is_err());
        }
    }

    #[test]
    fn combined_request_is_one_goal_and_preserves_duplicate_refusals() {
        let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"), "/../../examples/contact-pair/contact-pair.fsim"
        ))).unwrap().decoded.spec;
        let request = fs_project::spec::OutputRequest {
            name: COMBINED_OUTPUT.into(), kind: "report".into(), region: None,
        };
        spec.outputs.get_or_insert_with(Vec::new).push(request.clone());
        assert!(requested(&spec).unwrap());
        for other in [OUTPUT, contact_controls::OUTPUT, prescribed_controls::OUTPUT, COMBINED_OUTPUT] {
            let mut duplicate = spec.clone();
            duplicate.outputs.as_mut().unwrap().push(fs_project::spec::OutputRequest {
                name: other.into(), ..request.clone()
            });
            assert!(requested(&duplicate).is_err(), "{other}");
        }
        spec.outputs.as_mut().unwrap().last_mut().unwrap().kind = "scalar".into();
        assert!(requested(&spec).is_err());
    }

    #[test]
    fn volume_mean_is_one_explicit_objective_not_an_alias_of_the_maximum() {
        let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"), "/../../examples/contact-pair/contact-pair.fsim"
        ))).unwrap().decoded.spec;
        let request = fs_project::spec::OutputRequest {
            name: volume_mean::OUTPUT.into(), kind: "report".into(), region: None,
        };
        spec.outputs.get_or_insert_with(Vec::new).push(request.clone());
        assert!(requested(&spec).unwrap());
        for other in [OUTPUT, contact_controls::OUTPUT, prescribed_controls::OUTPUT,
            COMBINED_OUTPUT, volume_mean::OUTPUT] {
            let mut duplicate = spec.clone();
            duplicate.outputs.as_mut().unwrap().push(fs_project::spec::OutputRequest {
                name: other.into(), ..request.clone()
            });
            assert!(requested(&duplicate).is_err(), "{other}");
        }
        let mut wrong = spec.clone();
        wrong.outputs.as_mut().unwrap().last_mut().unwrap().kind = "scalar".into();
        assert!(requested(&wrong).is_err());
        wrong = spec.clone();
        wrong.outputs.as_mut().unwrap().last_mut().unwrap().region = Some("hot".into());
        assert!(requested(&wrong).is_err(), "output region syntax cannot silently change this volume goal");
        spec.requirements = None;
        assert!(requested(&spec).is_err());
    }
}
