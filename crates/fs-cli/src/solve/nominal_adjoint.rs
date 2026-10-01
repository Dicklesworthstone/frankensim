//! Opt-in nominal derivatives on the FINAL native conduction mesh and field.
//!
//! A report output named `temperature-max-adjoint` selects the smallest-index
//! hottest vertex of the declared temperature-max region. At a tie this is a
//! selected nodal functional, NOT a unique derivative of a relocating maximum.
//! It is useful as a frozen mean control even when an active branch changes.
//!
//! Reuse production adjoints, contact operators and full affine air feedback.
//! Only contractions of native input laws live here; no perturbed primal or
//! new solver. Nonlinear material and radiation derivatives are NOT frozen.
//! Fan speed, airflow inlet and prescribed-temperature controls remain explicit
//! unsupported rows until their complete physical parameter maps are wired.

use std::collections::BTreeMap;
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc};
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalAnalyzer, RobinFeedbackAnalysisConfig};
use fs_exec::Cx;
use super::{EvidenceWork, ProjectSpec, RungSolved, SolveRefusal, canonical_f64,
    conduction_error, json_string, temperature_maximum_region, trace_qoi_region_vertices};

mod contractions;

const OUTPUT: &str = "temperature-max-adjoint";
const MAX_PARAMETERS: usize = 256;
const SCOPE: &str = "Estimated derivative of a selected hottest nodal temperature on the final accepted native mesh. Fixed geometry, linear conductivity, matching contact and hydraulic operating point; complete affine air-reference feedback is differentiated, not frozen. Not a unique maximum derivative at a tie, a gradient-error enclosure, a continuum/shape derivative, experimental validation, or a parameter-uncertainty bound. Air inlet, fan speed and Dirichlet-temperature derivatives are not supplied.";

fn bad(message: impl Into<String>) -> SolveRefusal {
    conduction_error("cli-solve-nominal-adjoint", message,
        "request a report named temperature-max-adjoint for an admitted linear thermal solve; inspect unsupported parameter rows")
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

pub(super) fn requested(spec: &ProjectSpec) -> Result<bool, SolveRefusal> {
    let rows = spec.outputs.as_deref().unwrap_or(&[]);
    let mut found = false;
    for row in rows.iter().filter(|row| row.name == OUTPUT) {
        if found || row.kind != "report" {
            return Err(bad("temperature-max-adjoint requires exactly one report output"));
        }
        found = true;
    }
    if found && temperature_maximum_region(spec).is_none() {
        return Err(bad("the adjoint needs the existing declared temperature-max requirement and region"));
    }
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
    let region = temperature_maximum_region(spec).ok_or_else(|| bad("missing maximum region"))?;
    let region_id = *ids.get(region).ok_or_else(|| bad("maximum region has no mesh label"))?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("final operator was not retained"))?;
    if data.radiating_boundary.is_some() {
        return Err(bad("radiation needs its complete nonlinear adjoint; a frozen radiative Robin law is not admitted"));
    }
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing conduction setup"))?;
    let n = solved.mesh.vertex_count();
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let count = |divisor| usize::try_from(memory / divisor).unwrap_or(usize::MAX);
    let fan_count = spec.cooling.as_ref().and_then(|c| c.fan_system.as_ref())
        .map_or(0, |system| system.banks.len());
    let parameter_count = spec.power.as_deref().unwrap_or(&[]).len()
        .checked_add(setup.boundaries.len().saturating_mul(2))
        .and_then(|n| n.checked_add(fan_count))
        .ok_or_else(|| bad("adjoint parameter count overflow"))?;
    // Logical extra vectors/records, not a total-allocator/RSS promise. The
    // numerical owners independently enforce their matrix and iteration caps.
    if n > count(128) || solved.mesh.element_count() > count(64)
        || parameter_count > MAX_PARAMETERS {
        return Err(bad("nominal adjoint exceeds declared memory or 256-parameter work envelope"));
    }
    let (vertices, _) = trace_qoi_region_vertices(&solved.labels, &solved.mesh.complex().tets,
        n, region_id, work).map_err(|error| {
            if work.is_requested() { conduction_error("cli-solve-cancelled",
                "nominal region trace interrupted", "resume the accepted pipeline prefix") }
            else { bad(format!("nominal region trace refused: {error:?}")) }
        })?;
    let temperature = &solved.solution.temperature;
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
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &data.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let config = LinearGoalAnalysisConfig {
        residual_limits: fs_solver::goal::GoalResidualLimits { max_rows: count(256), max_nonzeros: count(64) },
        max_stability_iterations: data.linear.max_iterations,
    };
    let (lambda, residual, dual_iterations, stability_iterations, response_iterations, mode) = if data.air_paths.is_empty() {
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
        (gradient.nodal_load, gradient.relative_residual, gradient.iterations,
            None, Some(analyzer.response_iterations()), "linear-solid-full-air-feedback")
    };
    let source_bar = regional_source_pullback(cx, solved, &lambda)?;
    let volumes: BTreeMap<u32, f64> = audited.witness().per_region_auditor.iter()
        .map(|(region, volume)| (region.0, *volume)).collect();
    let mut rows = Vec::new();
    let mut missing = Vec::new();
    for (ordinal, power) in spec.power.as_deref().unwrap_or(&[]).iter().enumerate() {
        poll(cx)?;
        let id = *ids.get(&power.region).ok_or_else(|| bad("power region has no label"))?;
        let volume = *volumes.get(&id).ok_or_else(|| bad("power region has no audited volume"))?;
        if !(volume.is_finite() && volume > 0.0) { return Err(bad("invalid audited source volume")); }
        let value = finite(*source_bar.get(&id).ok_or_else(|| bad("source region has no element contribution"))? * (power.duty / volume))?;
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
            B::FixedTemperature { .. } => missing.push(unsupported("fixed-temperature", &declared.target,
                "prescribed values require their complete lift derivative")),
            B::NaturalConvection { .. } => missing.push(unsupported("natural-convection-ambient", &declared.target,
                "the coefficient's dependence on the wall temperature is not contracted")),
            B::AirflowConvection { branch, .. } => missing.push(unsupported("air-inlet-temperature", branch,
                "full temperature feedback is retained, but inlet and hydraulic parameter contractions are not implemented")),
        }
    }
    if let Some(system) = spec.cooling.as_ref().and_then(|c| c.fan_system.as_ref()) {
        for bank in &system.banks { poll(cx)?; missing.push(unsupported("fan-speed-ratio", &bank.bank_id,
            "requires flow, heat-transfer coefficient and air-capacity derivatives together")); }
    }
    let gap = second.map(|value| number(temperature[selected] - value)).transpose()?
        .unwrap_or_else(|| "null".into());
    poll(cx)?;
    Ok(format!(concat!("{{\"schema\":\"frankensim.cli.nominal-adjoint.v1\",\"output\":\"temperature-max-adjoint\",",
        "\"functional\":\"selected-nodal-temperature\",\"region\":{},\"selected_vertex\":{},",
        "\"value_k\":{},\"tied_maximum_vertices\":{},\"runner_up_gap_k\":{},\"mode\":{},",
        "\"true_relative_residual\":{},\"dual_iterations\":{},\"stability_iterations\":{},\"response_iterations\":{},",
        "\"parameters\":[{}],\"unsupported\":[{}],\"authority\":\"Estimated\",\"scope\":{}}}"),
        json_string(region), selected, number(temperature[selected])?, tied, gap, json_string(mode),
        number(residual)?, dual_iterations,
        stability_iterations.map_or_else(|| "null".into(), |n| n.to_string()),
        response_iterations.map_or_else(|| "null".into(), |n| n.to_string()),
        rows.join(","), missing.join(","), json_string(SCOPE)))
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
