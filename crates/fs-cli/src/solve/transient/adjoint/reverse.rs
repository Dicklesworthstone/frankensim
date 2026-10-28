//! Reverse the actual fine-grid history with one complete endpoint transpose.
use std::collections::BTreeMap;
use std::time::Instant;
use fs_conduction::{ConductionProblem, PowerMap, ScalarField, ThermalBc, ThermalInterfaces};
use fs_conduction::adjoint::robin::RobinResponse;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig};
use fs_conduction::radiation::{EnclosureRadiationLinearization, GrayDiffuseEnclosure};
use fs_exec::Cx;
use fs_project::ProjectSpec;
use crate::solve::conjugate::ConjugatePath;
use crate::solve::radiation::LoweredRadiation;
use crate::solve::natural::NaturalLaw;
use super::{Objective, Selected, Tape, finite, select, storage_bytes, flux_controls};
use super::enthalpy::{material_controls, prescribed_controls};
use super::super::{SolveRefusal, Workload, bad, lower, poll};

mod report;
mod radiation;
mod enclosure;
mod natural;
mod air;
mod air_controls;

pub(in crate::solve::transient) struct NaturalCoupling<'a> {
    pub(in crate::solve::transient) laws: &'a [NaturalLaw],
    pub(in crate::solve::transient) pressure_pa: f64,
}

pub(in crate::solve::transient) fn natural_storage_bytes(
    n: usize, ports: usize, natural_ports: usize,
) -> Result<usize, SolveRefusal> {
    natural::storage_bytes(n, ports, natural_ports)
}

pub(in crate::solve::transient) fn radiation_storage_bytes(
    n: usize, ports: usize, patches: usize, is_enclosure: bool,
) -> Result<usize, SolveRefusal> {
    let bytes = if is_enclosure {
        radiation::enclosure_storage_bytes(n, ports, patches)?
    } else {
        radiation::storage_bytes(n, ports, patches)?
    };
    usize::try_from(bytes)
        .map_err(|_| bad("native capacity radiation workspace does not fit this platform"))
}

pub(in crate::solve::transient) fn airflow_storage_bytes(
    spec: &ProjectSpec, n: usize, ports: usize, path: &ConjugatePath,
) -> Result<usize, SolveRefusal> {
    air::storage_bytes(n, ports, path.segments.len())?
        .checked_add(air_controls::storage_bytes(spec, path.segments.len())?)
        .ok_or_else(|| bad("native capacity airflow workspace overflows"))
}

/// Both storage coordinates share the same physical material/contact owner.
/// Capacity has no retained element-load tape; the shared allowance therefore
/// conservatively includes the enthalpy owner's element scratch as well.
pub(in crate::solve::transient) fn material_control_storage_bytes(
    spec: &ProjectSpec, nodes: usize, elements: usize,
) -> Result<usize, SolveRefusal> {
    material_controls::storage_bytes(
        nodes, elements, spec.materials.as_deref().map_or(0, |rows| rows.len()),
        spec.interface_cards.as_deref().map_or(0, |rows| rows.len()),
    )
}

/// The shared prescribed lift assembles the existing sparse spatial Jacobian.
/// Charge its bounded construction scratch and the deduplicated nodal bath map
/// before the forward run; there is no dense nodal matrix or extra linear solve.
pub(in crate::solve::transient) fn prescribed_control_storage_bytes(
    spec: &ProjectSpec, nodes: usize, elements: usize, faces: usize,
) -> Result<usize, SolveRefusal> {
    let count = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .map_or(0, |setup| setup.boundaries.iter().filter(|row| matches!(
            row.condition, fs_project::ThermalBoundaryCondition::FixedTemperature { .. },
        )).count());
    if count == 0 { return Ok(0); }
    nodes.checked_mul(128)
        .and_then(|v| elements.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .and_then(|v| faces.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .and_then(|v| count.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .filter(|&v| isize::try_from(v).is_ok())
        .ok_or_else(|| bad("capacity prescribed-temperature workspace overflows"))
}

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
    max_dual_solves_per_step: usize,
    worst_dual: f64,
    reconstructed_steps: usize,
    revalidated_steps: usize,
    additional_storage_bytes: usize,
    radiation: String,
    convection: String,
    material_controls: String,
    prescribed_temperature: String,
    prescribed_heat_flux: String,
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
    cx: &Cx<'_>, spec: &ProjectSpec, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>, ends_s: &[f64], workload: Option<&Workload<'_>>,
    config: StepConfig, nonlinear: Option<NonlinearStepConfig>, labels: &[u32],
    region_ids: &BTreeMap<String, u32>, goal_region: u32,
    radiation: Option<&LoweredRadiation>, enclosure: Option<&GrayDiffuseEnclosure>,
    natural: Option<&NaturalCoupling<'_>>,
    airflow: Option<&ConjugatePath>, tape: &Tape,
    objective: Objective, deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    poll(cx, deadline)?;
    let n = problem.mesh.vertex_count();
    if tape.n != n || tape.states.len() != ends_s.len() + 1
        || labels.len() != problem.mesh.element_count() {
        return Err(bad("transient adjoint mesh, trajectory or time-grid identity differs"));
    }
    if radiation.is_some_and(LoweredRadiation::is_enclosure) != enclosure.is_some() {
        return Err(bad("capacity adjoint lost its exact bound enclosure model"));
    }
    let coupled = radiation.is_some() || natural.is_some() || airflow.is_some();
    if coupled && tape.physical_targets_j.len() != ends_s.len()
        || !coupled && !tape.physical_targets_j.is_empty()
    {
        return Err(bad("capacity adjoint lost its original physical acceptance targets"));
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
    let radiation_workspace = radiation.map(|lowering|
        radiation_storage_bytes(n, ports.len(), lowering.surface_count(), lowering.is_enclosure()))
        .transpose()?.unwrap_or(0);
    let natural_workspace = natural.map(|coupling|
        natural_storage_bytes(n, ports.len(), coupling.laws.len()))
        .transpose()?.unwrap_or(0);
    let air_workspace = airflow.map(|path|
        airflow_storage_bytes(spec, n, ports.len(), path))
        .transpose()?.unwrap_or(0);
    // This report field is a subset of the natural/air allowances already
    // charged above; it does not add a second charge to the admitted total.
    let natural_control_bytes = natural.map_or(0, |coupling| coupling.laws.len() * 2048);
    let control_workspace_bytes = airflow.map(|path|
        air_controls::storage_bytes(spec, path.segments.len()))
        .transpose()?.unwrap_or(0).checked_add(natural_control_bytes)
        .ok_or_else(|| bad("native capacity cooling-control workspace overflows"))?;
    let material_workspace = if objective.material_controls() {
        material_control_storage_bytes(spec, n, problem.mesh.element_count())?
    } else { 0 };
    let prescribed_workspace = prescribed_control_storage_bytes(
        spec, n, problem.mesh.element_count(), problem.mesh.boundary().len(),
    )?;
    let flux_workspace = flux_controls::storage_bytes(spec, n)?;
    let additional_storage_bytes = storage_bytes(n, ends_s.len(), names.len())?
        .checked_add(radiation_workspace)
        .and_then(|value| value.checked_add(natural_workspace))
        .and_then(|value| value.checked_add(air_workspace))
        .and_then(|value| value.checked_add(material_workspace))
        .and_then(|value| value.checked_add(prescribed_workspace))
        .and_then(|value| value.checked_add(flux_workspace))
        .ok_or_else(|| bad("native capacity adjoint total workspace overflows"))?;
    let patch_count = radiation.map_or(0, LoweredRadiation::surface_count);
    let mut emissivities = vec![0.0; patch_count];
    let mut reservoirs = vec![0.0; patch_count];
    let mut natural_ambients = vec![0.0; natural.map_or(0, |coupling| coupling.laws.len())];
    let mut air_controls = airflow.map(|path| air_controls::Controls::new(cx, spec, path, deadline))
        .transpose()?;
    // Bind even for an initial selected peak so zero controls retain physical
    // identities. Ordinary capacity reports do not acquire this extra work.
    let mut material_controls = if objective.material_controls() {
        Some(material_controls::Controls::new(
            cx, spec, problem, interfaces, labels, region_ids, &tape.states[0], deadline,
        )?)
    } else { None };
    let mut prescribed_controls = prescribed_controls::Controls::new(cx, spec, problem, deadline)?;
    let mut flux_controls = flux_controls::Controls::new(cx, spec, problem, deadline)?;
    let fixed_ports: Vec<usize> = ports.iter().enumerate().filter_map(|(index, &name)| {
        let owned = natural.is_some_and(|coupling|
            coupling.laws.iter().any(|law| law.target == name))
            || airflow.is_some_and(|path| path.segments.iter().any(|segment| segment.target == name));
        (!owned).then_some(index)
    }).collect();
    // Enclosure flux feedback is independent of the convective port layout:
    // vacuum walls need two node-by-patch factors even with no Robin controls.
    let flux_entries = if enclosure.is_some() {
        n.checked_mul(patch_count).and_then(|entries| entries.checked_mul(2))
            .ok_or_else(|| bad("capacity enclosure flux feedback entry count overflows"))?
    } else { 0 };
    let feedback_entry_cap = n.checked_mul(ports.len()).and_then(|v| v.checked_mul(4))
        .and_then(|entries| entries.checked_add(flux_entries))
        .ok_or_else(|| bad("capacity cooling feedback entry count overflows"))?;
    let zero_feedback = if enclosure.is_some() || natural.is_some() || airflow.is_some()
        || prescribed_controls.has_controls() {
        vec![0.0; ports.len() * ports.len()]
    } else {
        Vec::new()
    };
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
    let mut radiosity_iterations = 0usize;
    let mut worst_radiosity_dual = 0.0_f64;
    let max_dual_solves_per_step = enclosure.map_or(1, |model| model.surfaces().len() + 2);
    let prescribed_entries = problem.mesh.element_count().checked_mul(32)
        .and_then(|v| problem.mesh.boundary().len().checked_mul(27).and_then(|b| v.checked_add(b)))
        .and_then(|v| v.checked_add(n))
        .and_then(|v| v.checked_add(ports.len() * ports.len()))
        .and_then(|v| v.checked_add(feedback_entry_cap))
        .ok_or_else(|| bad("capacity prescribed-temperature entry allowance overflows"))?;
    for index in (0..selected.state).rev() {
        poll(cx, deadline)?;
        let start = if index == 0 { 0.0 } else { ends_s[index - 1] };
        let dt = ends_s[index] - start;
        let source = workload.map(|history|
            history.source(cx, problem.mesh, labels, start, deadline)).transpose()?;
        let step_problem = ConductionProblem { source: source.as_ref().unwrap_or(problem.source), ..problem };
        let expected = &tape.states[index + 1];
        let natural_bound = natural.map(|coupling| natural::bind(
            cx, step_problem, expected, coupling.laws, coupling.pressure_pa, &ports, deadline,
        )).transpose()?;
        let natural_problem = ConductionProblem {
            boundary: natural_bound.as_ref().map_or(step_problem.boundary, |bound| &bound.boundary),
            ..step_problem
        };
        let air_bound = airflow.map(|path| air::bind(
            cx, natural_problem, expected, path, &ports, deadline,
        )).transpose()?;
        let convection_problem = ConductionProblem {
            boundary: air_bound.as_ref().map_or(natural_problem.boundary, |bound| &bound.boundary),
            ..natural_problem
        };
        // Revalidate the actual physical endpoint before assembling a tangent.
        // A whole-body heat balance cannot detect cancelling nodal defects.
        // This single audit keeps original convection, gray radiation, and
        // actual transported air heat distinct; no primal solve is repeated.
        let physical = if coupled {
            let audit = if let Some(model) = enclosure {
                let lowering = radiation.ok_or_else(|| bad("capacity enclosure lost its controls"))?;
                engine.evaluate_prescribed_enclosure_endpoint(
                    cx, convection_problem, interfaces, &tape.states[index], expected, dt, config,
                    model, lowering.config,
                ).map_err(lower)?.physical
            } else {
                engine.evaluate_prescribed_endpoint(
                    cx, convection_problem, interfaces, &tape.states[index], expected, dt, config,
                    radiation.map(|lowering| lowering.patches.as_slice()),
                ).map_err(lower)?
            };
            let target = tape.physical_targets_j[index];
            if !target.is_finite() || target < 0.0
                || !audit.residual_norm_j.is_finite() || audit.residual_norm_j > target
                || !audit.energy_residual_j.is_finite()
                || audit.energy_residual_j.abs() > config.energy_tolerance_j
            {
                return Err(bad(format!(
                    "capacity cooling endpoint fails its retained physical target: residual {} J / {target} J, energy {} J / {} J",
                    audit.residual_norm_j, audit.energy_residual_j, config.energy_tolerance_j,
                )));
            }
            if let Some(bound) = &air_bound {
                bound.check_air_energy(&audit, dt, config.energy_tolerance_j)?;
            }
            Some(audit)
        } else {
            None
        };
        let enclosure_bound = enclosure.map(|model| {
            let lowering = radiation.ok_or_else(|| bad("capacity enclosure lost its declaration"))?;
            EnclosureRadiationLinearization::bind(
                cx, convection_problem, config.linear, expected, &ports,
                natural_bound.as_ref().map_or(zero_ports.as_slice(), |bound| &bound.h_mean_slopes),
                air_bound.as_ref().map_or(zero_feedback.as_slice(), |bound| &bound.reference_feedback),
                model, lowering.config.balance_tolerance_w, feedback_entry_cap,
            ).map_err(lower)
        }).transpose()?;
        let bound = if let Some(lowering) = radiation.filter(|lowering| !lowering.is_enclosure()) {
            Some(radiation::bind(
                cx, convection_problem, expected, config, lowering, &ports,
                natural_bound.as_ref().map(|bound| bound.h_mean_slopes.as_slice()),
                air_bound.as_ref().map(|bound| bound.reference_feedback.as_slice()),
                physical.ok_or_else(|| bad("capacity radiation lost its physical endpoint audit"))?,
                deadline,
            )?)
        } else {
            None
        };
        let step = if coupled {
            let boundary = if let Some(bound) = &enclosure_bound {
                bound.boundary()
            } else {
                bound.as_ref().map_or(convection_problem.boundary, |bound| &bound.boundary)
            };
            engine.linearize_prescribed_endpoint(
                cx, ConductionProblem { boundary, ..convection_problem },
                interfaces, &tape.states[index], expected, dt, config,
                tape.physical_targets_j[index], &ports,
            )
        } else {
            engine.linearize_prescribed_step(
                cx, step_problem, interfaces, &tape.states[index], dt, config, nonlinear, &ports,
            )
        }.map_err(lower)?;
        if step.temperature().len() != expected.len() { return Err(bad("adjoint reconstruction changed field shape")); }
        for (v, (&actual, &expected)) in step.temperature().iter().zip(expected).enumerate() {
            if v % 256 == 0 { poll(cx, deadline)?; }
            if actual.to_bits() != expected.to_bits() {
                return Err(bad(format!("transient adjoint reconstruction diverged at state {}, vertex {v}", index + 1)));
            }
        }
        let gradient = if let Some(bound) = &enclosure_bound {
            let combined = step.pullback_with_robin_and_flux_feedback(
                cx, bound.convective_mean_slopes(), bound.local_reference_slopes(),
                bound.reference_feedback(), bound.flux_feedback(), &carry, feedback_entry_cap,
            ).map_err(lower)?;
            let physical = bound.remap_gradient(cx, combined, config.linear).map_err(lower)?;
            let model = enclosure.ok_or_else(|| bad("capacity enclosure gradient lost its model"))?;
            if physical.emissivities.len() != patch_count
                || physical.radiation_regions.len() != patch_count
                || physical.regions.len() != ports.len()
            {
                return Err(bad("capacity enclosure control layout changed during reverse traversal"));
            }
            for (i, surface) in model.surfaces().iter().enumerate() {
                if physical.radiation_regions[i] != surface.name() {
                    return Err(bad("capacity enclosure control order differs from its view factors"));
                }
                emissivities[i] = finite(emissivities[i] + physical.emissivities[i])?;
            }
            radiosity_iterations = radiosity_iterations.checked_add(physical.radiosity_iterations)
                .ok_or_else(|| bad("capacity enclosure reverse work overflows"))?;
            worst_radiosity_dual = worst_radiosity_dual.max(physical.radiosity_relative_residual);
            physical.convection
        } else { match &bound {
            Some(bound) => {
                let combined = step.pullback_with_robin_feedback(
                    cx, &bound.h_mean_slopes, &bound.local_reference_slopes,
                    &bound.reference_feedback, &carry, bound.max_feedback_entries,
                ).map_err(lower)?;
                let physical = bound.remap_gradient(cx, combined, deadline)?;
                if physical.emissivities.len() != patch_count
                    || physical.reservoir_temperatures.len() != patch_count
                {
                    return Err(bad("capacity radiation control layout changed during reverse traversal"));
                }
                for i in 0..patch_count {
                    emissivities[i] = finite(emissivities[i] + physical.emissivities[i])?;
                    reservoirs[i] = finite(reservoirs[i] + physical.reservoir_temperatures[i])?;
                }
                physical.convection
            },
            None if natural_bound.is_some() || air_bound.is_some() => {
                let h_slopes = natural_bound.as_ref()
                    .map_or(zero_ports.as_slice(), |bound| bound.h_mean_slopes.as_slice());
                let reference_feedback = air_bound.as_ref()
                    .map_or(zero_feedback.as_slice(), |bound| bound.reference_feedback.as_slice());
                step.pullback_with_robin_feedback(
                    cx, h_slopes, &zero_ports, reference_feedback, &carry, feedback_entry_cap,
                ).map_err(lower)?
            },
            None => step.pullback(cx, &carry, &zero_ports, &zero_ports).map_err(lower)?,
        } };
        if !gradient.relative_residual.is_finite() || gradient.relative_residual >= config.linear.tolerance {
            return Err(bad("native transient transpose failed its recomputed residual gate"));
        }
        if prescribed_controls.has_controls() {
            let (boundary, h_slopes, reference_slopes, reference_feedback, flux_feedback) =
                if let Some(bound) = &enclosure_bound {
                    (bound.boundary(), bound.convective_mean_slopes(),
                        bound.local_reference_slopes(), bound.reference_feedback(), bound.flux_feedback())
                } else if let Some(bound) = &bound {
                    (&bound.boundary, bound.h_mean_slopes.as_slice(),
                        bound.local_reference_slopes.as_slice(), bound.reference_feedback.as_slice(), None)
                } else {
                    (convection_problem.boundary,
                        natural_bound.as_ref().map_or(zero_ports.as_slice(), |bound| bound.h_mean_slopes.as_slice()),
                        zero_ports.as_slice(),
                        air_bound.as_ref().map_or(zero_feedback.as_slice(), |bound| bound.reference_feedback.as_slice()), None)
                };
            // C is row-sum lumped: C_FD is zero, so the spatial lift is exactly
            // the endpoint lift. The complete physical dual already contains
            // C/dt and every feedback term; never multiply this partial by dt.
            let bars = RobinResponse::prescribed_temperature_pullback_with_flux_feedback(
                cx, ConductionProblem { boundary, ..convection_problem }, interfaces,
                expected, &carry, &gradient.nodal_load, &ports, h_slopes,
                reference_slopes, reference_feedback, flux_feedback, prescribed_entries,
            ).map_err(lower)?;
            prescribed_controls.add(cx, &bars, deadline)?;
        } else {
            // Retain the exact selected-prefix accounting for an empty set.
            prescribed_controls.add(cx, &gradient.nodal_load, deadline)?;
        }
        if let Some(bound) = &natural_bound {
            let ambient = bound.ambient_partials(cx, &gradient, deadline)?;
            if ambient.len() != natural_ambients.len() {
                return Err(bad("capacity natural controls changed during reverse traversal"));
            }
            for (total, value) in natural_ambients.iter_mut().zip(ambient) {
                *total = finite(*total + value)?;
            }
        }
        if let Some(controls) = &mut air_controls {
            controls.add(cx, spec, convection_problem, expected, &ports, &gradient, deadline)?;
        }
        // Fixed-state material partials use the same complete watt-coordinate
        // dual. Neither K(T) nor matching resistance receives another dt.
        if let Some(controls) = &mut material_controls {
            controls.add_capacity(
                cx, convection_problem, interfaces, expected, &gradient.nodal_load, deadline,
            )?;
        }
        let density_bar = step.source_density_pullback(cx, &gradient.nodal_load).map_err(lower)?;
        flux_controls.add(cx, &gradient.nodal_load, deadline)?;
        for (region, basis) in bases.iter().enumerate() {
            source_controls[index][region] = dot(cx, &density_bar, basis, deadline)?;
        }
        source_multipliers[index] = step.source_multiplier_pullback(cx, &gradient.nodal_load).map_err(lower)?;
        capacities = finite(capacities + step.capacity_multiplier_pullback(cx, &gradient.nodal_load).map_err(lower)?)?;
        for i in 0..ports.len() {
            if natural_bound.as_ref().is_some_and(|bound| bound.owns_region(ports[i]))
                || air_bound.as_ref().is_some_and(|bound| bound.owns_region(ports[i])) {
                continue;
            }
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
    dual_iterations = dual_iterations.checked_add(radiosity_iterations)
        .ok_or_else(|| bad("combined capacity enclosure reverse work overflows"))?;
    worst_dual = worst_dual.max(worst_radiosity_dual);
    let reconstructed_steps = if coupled { 0 } else { selected.state };
    let revalidated_steps = if coupled { selected.state } else { 0 };
    let natural_fields = natural.map(|coupling|
        natural::fragment(coupling.laws, &natural_ambients, selected.state))
        .transpose()?.unwrap_or_default();
    let air_fields = air_controls.as_ref().map(|controls| controls.fragment(selected.state))
        .transpose()?.unwrap_or_default();
    let convection = if natural.is_some() || airflow.is_some() {
        let method = if airflow.is_some() {
            "complete-capacity-cooling-tangent"
        } else {
            "complete-natural-capacity-tangent"
        };
        let parameters = if airflow.is_some() {
            "independent-natural-reservoirs-air-inlets-and-absolute-fan-bank-speeds-throughout-trajectory"
        } else {
            "independent-natural-reservoir-temperatures-throughout-trajectory"
        };
        format!(
            ",\"convection\":{{\"method\":\"{method}\",\"control_steps\":{},\"control_workspace_bytes\":{control_workspace_bytes},\"residual_target\":\"original-accepted-forward-physical-tolerance\",\"parameters\":\"{parameters}\"{natural_fields}{air_fields}}}",
            selected.state,
        )
    } else {
        String::new()
    };
    let radiation = radiation.map(|lowering| {
        if let Some(model) = enclosure {
            enclosure::fragment(
                lowering, model, &emissivities, selected.state,
                radiosity_iterations, worst_radiosity_dual,
            )
        } else {
            radiation::fragment(lowering, &emissivities, &reservoirs, selected.state)
        }
    }).transpose()?.unwrap_or_default();
    let material_controls = material_controls.as_ref()
        .map(|controls| controls.fragment(selected.state)).transpose()?.unwrap_or_default();
    let prescribed_temperature = prescribed_controls.fragment_for_capacity(selected.state)?;
    let prescribed_heat_flux = flux_controls.fragment(selected.state)?;
    let result = Derivatives {
        selected, initial: carry, region_names: names, source_controls, source_multipliers,
        capacity: capacities, ports: fixed_ports.iter().map(|&i| ports[i].to_owned()).collect(),
        references: fixed_ports.iter().map(|&i| references[i]).collect(),
        coefficients: fixed_ports.iter().map(|&i| coefficients[i]).collect(),
        primal_iterations, dual_iterations, max_dual_solves_per_step, worst_dual,
        reconstructed_steps, revalidated_steps, additional_storage_bytes, radiation, convection,
        material_controls, prescribed_temperature, prescribed_heat_flux,
    };
    report::render(cx, &result, ends_s, objective, config, deadline)
}

