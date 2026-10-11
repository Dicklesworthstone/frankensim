//! Reverse the actual fine-grid history with one complete endpoint transpose.
use std::collections::BTreeMap;
use std::time::Instant;
use fs_conduction::{ConductionProblem, PowerMap, ScalarField, ThermalBc, ThermalInterfaces};
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_exec::Cx;
use super::{Objective, Selected, Tape, finite, select, storage_bytes};
use super::super::{SolveRefusal, Workload, bad, lower, poll};

mod report;

struct Derivatives {
    selected: Selected,
    initial: Vec<f64>,
    region_names: Vec<String>,
    source_controls: Vec<Vec<f64>>,
    source_multipliers: Vec<f64>,
    capacity: f64,
    ports: Vec<String>,
    references: Vec<f64>,
    coefficients: Vec<f64>,
    primal_iterations: usize,
    dual_iterations: usize,
    worst_dual: f64,
}

fn dot(cx: &Cx<'_>, a: &[f64], b: &[f64], deadline: Option<(Instant, f64)>)
    -> Result<f64, SolveRefusal> {
    if a.len() != b.len() { return Err(bad("transient adjoint contraction shape differs")); }
    let mut total = 0.0;
    for (i, (&a, &b)) in a.iter().zip(b).enumerate() {
        if i % 256 == 0 { poll(cx, deadline)?; }
        total = finite(a.mul_add(b, total))?;
    }
    Ok(total)
}

#[allow(clippy::too_many_arguments)]
pub(in crate::solve::transient) fn finish(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>, ends_s: &[f64], workload: Option<&Workload<'_>>,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>, labels: &[u32],
    region_ids: &BTreeMap<String, u32>, goal_region: u32, tape: &Tape,
    objective: Objective, deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    poll(cx, deadline)?;
    let n = problem.mesh.vertex_count();
    if tape.n != n || tape.states.len() != ends_s.len() + 1
        || labels.len() != problem.mesh.element_count() {
        return Err(bad("transient adjoint mesh, trajectory or time-grid identity differs"));
    }
    let mut selected_vertices = vec![false; n];
    let mut volumes = BTreeMap::<u32, f64>::new();
    for (element, (&id, tet)) in labels.iter().zip(&problem.mesh.complex().tets).enumerate() {
        if element % 256 == 0 { poll(cx, deadline)?; }
        let volume = volumes.entry(id).or_default();
        *volume = finite(*volume + problem.mesh.element_volume(element))?;
        if id == goal_region { for &v in tet { selected_vertices[v as usize] = true; } }
    }
    if volumes.len() > 64 { return Err(bad("transient adjoint exceeds 64 regional source controls")); }
    let vertices: Vec<usize> = selected_vertices.iter().enumerate()
        .filter_map(|(v, &selected)| selected.then_some(v)).collect();
    let selected = select(cx, tape, objective, &vertices, deadline)?;
    let mut bases = Vec::new();
    let mut names = Vec::new();
    for (&id, &volume) in &volumes {
        poll(cx, deadline)?;
        if volume <= 0.0 { return Err(bad("a source-control region has no positive volume")); }
        let mut matches = region_ids.iter().filter(|(_, found)| **found == id);
        let name = matches.next().map(|(name, _)| name.clone())
            .ok_or_else(|| bad("source-control region has no declared name"))?;
        if matches.next().is_some() { return Err(bad("source-control region has ambiguous names")); }
        let densities = volumes.keys().map(|&other|
            (other, if other == id { 1.0 / volume } else { 0.0 })).collect();
        // The same P1 projector as the forward source, including mixed nodes.
        // Baseline power is never a divisor: zero-power regions remain controls.
        let basis = PowerMap::regional_volumetric_source(problem.mesh, labels, &densities).map_err(lower)?;
        let values = match basis {
            ScalarField::Uniform(value) => vec![value; n],
            ScalarField::Nodal(values) => values,
        };
        if values.len() != n { return Err(bad("source-control projector lost its nodal layout")); }
        names.push(name);
        bases.push(values);
    }
    let ports: Vec<&str> = problem.boundary.conditions().iter().enumerate()
        .filter_map(|(i, condition)| matches!(condition, ThermalBc::Robin { .. })
            .then_some(problem.boundary.region_names()[i].as_str())).collect();
    if ports.len() > 64 { return Err(bad("transient adjoint exceeds 64 boundary controls")); }
    let mut carry = vec![0.0; n];
    carry[selected.vertex] = 1.0;
    let mut source_controls = vec![vec![0.0; names.len()]; ends_s.len()];
    let mut source_multipliers = vec![0.0; ends_s.len()];
    let mut capacities = 0.0;
    let mut references = vec![0.0; ports.len()];
    let mut coefficients = vec![0.0; ports.len()];
    let zero_ports = vec![0.0; ports.len()];
    let mut primal_iterations = 0usize;
    let mut dual_iterations = 0usize;
    let mut worst_dual = 0.0_f64;
    for index in (0..selected.state).rev() {
        poll(cx, deadline)?;
        let start = if index == 0 { 0.0 } else { ends_s[index - 1] };
        let dt = ends_s[index] - start;
        let source = workload.map(|history|
            history.source(cx, problem.mesh, labels, start, deadline)).transpose()?;
        let step_problem = ConductionProblem { source: source.as_ref().unwrap_or(problem.source), ..problem };
        let step = engine.linearize_prescribed_step(cx, step_problem, interfaces,
            &tape.states[index], dt, config, nonlinear, &ports).map_err(lower)?;
        let expected = &tape.states[index + 1];
        if step.temperature().len() != expected.len() { return Err(bad("adjoint reconstruction changed field shape")); }
        for (v, (&actual, &expected)) in step.temperature().iter().zip(expected).enumerate() {
            if v % 256 == 0 { poll(cx, deadline)?; }
            if actual.to_bits() != expected.to_bits() {
                return Err(bad(format!("transient adjoint reconstruction diverged at state {}, vertex {v}", index + 1)));
            }
        }
        let gradient = step.pullback(cx, &carry, &zero_ports, &zero_ports).map_err(lower)?;
        if !gradient.relative_residual.is_finite() || gradient.relative_residual >= config.linear.tolerance {
            return Err(bad("native transient transpose failed its recomputed residual gate"));
        }
        let density_bar = step.source_density_pullback(cx, &gradient.nodal_load).map_err(lower)?;
        for (region, basis) in bases.iter().enumerate() {
            source_controls[index][region] = dot(cx, &density_bar, basis, deadline)?;
        }
        source_multipliers[index] = step.source_multiplier_pullback(cx, &gradient.nodal_load).map_err(lower)?;
        capacities = finite(capacities + step.capacity_multiplier_pullback(cx, &gradient.nodal_load).map_err(lower)?)?;
        for i in 0..ports.len() {
            references[i] = finite(references[i] + gradient.references[i])?;
            coefficients[i] = finite(coefficients[i] + gradient.log_htc[i])?;
        }
        carry = step.previous_temperature_pullback(cx, &gradient.nodal_load).map_err(lower)?;
        primal_iterations = primal_iterations.checked_add(step.primal().krylov_iterations)
            .ok_or_else(|| bad("adjoint reconstruction iteration count overflows"))?;
        dual_iterations = dual_iterations.checked_add(gradient.iterations)
            .ok_or_else(|| bad("adjoint iteration count overflows"))?;
        worst_dual = worst_dual.max(gradient.relative_residual);
    }
    let result = Derivatives {
        selected, initial: carry, region_names: names, source_controls, source_multipliers,
        capacity: capacities, ports: ports.iter().map(|s| (*s).to_owned()).collect(),
        references, coefficients, primal_iterations, dual_iterations, worst_dual,
    };
    report::render(cx, &result, ends_s, objective, config, deadline)
}
