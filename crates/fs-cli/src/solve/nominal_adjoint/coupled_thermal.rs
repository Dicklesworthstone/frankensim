//! One full natural/material/radiation/air adjoint on the unchanged native field.
//! The actual air march supplies references; its analytic wall map supplies
//! feedback. Radiation remains a separate reservoir, never a heat load to air.

use fs_conduction::adjoint::RobinGradient;
use fs_conduction::radiation::AmbientRadiationGradient;
use super::{ConductionProblem, Cx, ProjectSpec, RungSolved, SolveRefusal,
    bad, convective_feedback, lower, poll, row};
use std::collections::BTreeMap;

pub(super) struct Pullback {
    pub gradient: RobinGradient,
    pub natural_ambient: BTreeMap<String, f64>,
    pub air_rows: Vec<String>,
    pub radiation_rows: Vec<String>,
}

pub(super) fn pullback(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, weights: &[f64],
) -> Result<Pullback, SolveRefusal> {
    poll(cx)?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("missing coupled thermal operator"))?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing coupled thermal setup"))?;
    if data.air_paths.is_empty() {
        return Err(bad("coupled thermal adjoints require retained fixed-flow air paths"));
    }
    let n = solved.mesh.vertex_count();
    let temperature = &solved.solution.temperature;
    if temperature.len() != n || weights.len() != n { return Err(bad("coupled thermal field/goal size differs")); }
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let entries = usize::try_from(memory/32).unwrap_or(usize::MAX);
    let bound = convective_feedback::bind(cx, spec, solved, entries)?;
    let names: Vec<_> = bound.names.iter().map(String::as_str).collect();
    let surfaces = setup.radiation.as_ref().map_or(&[][..], |r| r.surfaces.as_slice());
    if surfaces.len() != data.radiation_patches.len() || surfaces.len() > 64 {
        return Err(bad("coupled thermal adjoint did not retain every radiation law"));
    }
    for (surface, patch) in surfaces.iter().zip(&data.radiation_patches) {
        poll(cx)?;
        if surface.target != patch.region()
            || surface.reservoir_temperature.value != patch.ambient_temperature_k()
            || surface.card != patch.emissivity().card_identity().to_hex()
            || surface.query_temperature.value != patch.emissivity().temperature_k() {
            return Err(bad("coupled radiation law differs from its declared physical source"));
        }
    }
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &bound.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let gradient = AmbientRadiationGradient::pullback_with_reference_feedback(cx, problem,
        data.interfaces.as_ref(), data.linear, temperature, &names, &bound.natural.wall_slopes, &bound.feedback,
        &data.radiation_patches, weights, entries).map_err(|e| {
            let mut refusal = lower(e);
            refusal.fix = "tighten the physical solid/air and radiation convergence controls; the unchanged field and complete coupled transpose must satisfy their residual gates".into();
            refusal
        })?;
    let natural_ambient = bound.natural_ambient(&gradient.convection)?;
    let air_rows = bound.air_rows(cx, setup, solved, &gradient.convection)?;
    let mut radiation_rows = Vec::new();
    for (i, surface) in surfaces.iter().enumerate() {
        radiation_rows.push(row("radiation-reservoir-temperature", &surface.name, i, "K",
            gradient.reservoir_temperatures[i])?);
        radiation_rows.push(row("radiation-emissivity", &surface.name, i, "1", gradient.emissivities[i])?);
    }
    poll(cx)?;
    Ok(Pullback { gradient: gradient.convection, natural_ambient, air_rows, radiation_rows })
}
