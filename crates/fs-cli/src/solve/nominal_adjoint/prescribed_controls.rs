//! Native fixed-temperature controls, including the complete eliminated lift.
use std::collections::{BTreeMap, BTreeSet};
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc, ThermalBoundary};
use fs_conduction::adjoint::RobinResponse;
use super::{Cx, ProjectSpec, RungSolved, SolveRefusal, bad, finite, lower, poll, row, unsupported, zeros};

mod law;

pub(super) const OUTPUT: &str = "temperature-max-boundary-adjoint";
pub(super) const SCOPE: &str = "The boundary-extended request additionally differentiates each independently prescribed uniform temperature in kelvin. It includes the direct selected-node objective term and the complete material/contact lift, natural-convection and radiation coefficient/reference feedback, and upstream air-reference feedback. Prescribed nodes remain in every wall mean. Two fixed boundary declarations sharing a vertex cannot be independently varied and are explicitly unsupported rather than assigned duplicate sensitivities. Geometry, contact laws, flow and transport properties remain fixed; these are local discrete derivatives, not continuum or gradient-error bounds.";

pub(super) fn requested(spec: &ProjectSpec) -> bool {
    spec.outputs.as_deref().unwrap_or(&[]).iter().any(|r| r.name == OUTPUT)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn append(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, weights: &[f64], lambda: &[f64],
    rows: &mut Vec<String>, missing: &mut Vec<String>, max_entries: usize,
) -> Result<(), SolveRefusal> {
    poll(cx)?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("prescribed controls require the declared conduction setup"))?;
    if !setup.boundaries.iter().any(|r| matches!(r.condition,
        fs_project::ThermalBoundaryCondition::FixedTemperature {..})) { return Ok(()); }
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("prescribed controls lost the thermal operator"))?;
    let mut vertices = BTreeMap::<usize,BTreeSet<usize>>::new();
    let mut owners = BTreeMap::<usize,BTreeSet<usize>>::new();
    for (slot, face) in solved.mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let Some(region) = data.boundary.region_for(slot) else { continue; };
        if !matches!(data.boundary.conditions()[region],ThermalBc::Dirichlet {..}) { continue; }
        for &v in &face.vertices {
            vertices.entry(region).or_default().insert(v as usize);
            owners.entry(v as usize).or_default().insert(region);
        }
    }
    let bound = law::bind(cx,spec,solved,max_entries)?;
    let names: Vec<_> = bound.names.iter().map(String::as_str).collect();
    let problem = ConductionProblem {mesh:&solved.mesh,boundary:&bound.boundary,
        material:&data.fallback,element_materials:Some(&data.materials),source:&data.source};
    let bars = RobinResponse::prescribed_temperature_pullback(cx,problem,data.interfaces.as_ref(),
        &solved.solution.temperature,weights,lambda,&names,&bound.h_slopes,
        &bound.reference_slopes,&bound.feedback,max_entries).map_err(lower)?;
    let mut seen = BTreeSet::new();
    for (ordinal, declared) in setup.boundaries.iter().enumerate() {
        poll(cx)?;
        let fs_project::ThermalBoundaryCondition::FixedTemperature {temperature}
            = &declared.condition else { continue; };
        let region = data.boundary.region_names().iter().position(|n| n == &declared.target)
            .ok_or_else(|| bad("prescribed control has no native boundary identity"))?;
        if !seen.insert(region) { return Err(bad("prescribed control has duplicate declaration ownership")); }
        let ThermalBc::Dirichlet {temperature:ScalarField::Uniform(actual)}
            = &data.boundary.conditions()[region] else { return Err(bad("prescribed control changed its uniform boundary law")); };
        if *actual != temperature.value { return Err(bad("prescribed control differs from the retained physical value")); }
        let selected = vertices.get(&region).filter(|v| !v.is_empty())
            .ok_or_else(|| bad("prescribed control has no retained trace vertices"))?;
        if selected.iter().any(|v| owners[v].len() != 1) {
            missing.push(unsupported("fixed-temperature",&declared.target,
                "this trace shares prescribed vertices with another fixed-temperature declaration; independent variation would impose conflicting values"));
            continue;
        }
        let mut value = 0.0;
        for (i, &v) in selected.iter().enumerate() {
            if i % 512 == 0 { poll(cx)?; }
            value = finite(value+bars[v])?;
        }
        rows.push(row("fixed-temperature",&declared.target,ordinal,"K",value)?);
    }
    poll(cx)
}
