//! Reverse accepted native enthalpy endpoints without another nonlinear solve.
//! The selected temperature seed is mapped to h once; history carry stays in
//! specific-enthalpy coordinates, including across latent plateaus.

use std::collections::BTreeMap;
use std::time::Instant;

use fs_conduction::transient::enthalpy::adjoint::{
    EnthalpyAdjointError, EnthalpyStepLinearization,
};
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyRadiationStepLinearization,
    HeterogeneousEnthalpyStepLinearization,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
};
use fs_conduction::{
    ConductionError, ConductionProblem, LinearConfig, PowerMap, ScalarField, ThermalInterfaces,
};
use fs_exec::Cx;
use fs_project::{ProjectSpec, TransientEnthalpy};
use fs_solver::NewtonReport;

use super::{Objective, Tape, finite, select};
use super::super::{SolveRefusal, Workload, airflow, bad, buoyancy, json_string,
    lower, number, poll, radiation};

mod convection;
mod prescribed_controls;
pub(super) mod material_controls;
pub(in crate::solve::transient) use material_controls::storage_bytes as material_control_storage_bytes;
pub(in crate::solve::transient) use convection::storage_bytes as convection_storage_bytes;
pub(in crate::solve::transient) use convection::control_storage_bytes;

// Forward native enthalpy also admits at most 24 Arnoldi columns. Keep the
// transpose within the same bounded column scale before allocating its basis.
const REVERSE_RESTART_LIMIT: usize = 24;

// Two accumulated surface controls and bounded JSON rows, independent of the
// number of fine steps. The ordinary radiation allowance covers active core
// parameter scratch, source metadata and face records.
const RADIATION_PARAMETER_BYTES_PER_SURFACE: usize = 2 * std::mem::size_of::<f64>() + 2048;

enum BoundEndpoint<'m> {
    Prescribed(HeterogeneousEnthalpyStepLinearization<'m>),
    Radiative(HeterogeneousEnthalpyRadiationStepLinearization<'m>),
}

impl BoundEndpoint<'_> {
    fn transport(&self) -> &EnthalpyStepLinearization<'_> {
        match self {
            Self::Prescribed(step) => step.transport(),
            Self::Radiative(step) => step.transport(),
        }
    }
}

fn failed(error: EnthalpyAdjointError) -> SolveRefusal {
    match error {
        EnthalpyAdjointError::Enthalpy(EnthalpyError::Conduction(
            error @ ConductionError::Cancelled { .. },
        )) => lower(error),
        other => bad(format!("native enthalpy adjoint refused: {other}")),
    }
}

/// Additional checkpoints, source bases, one revalidation/transpose workspace
/// and derivative output, including capped per-cycle Krylov history/report
/// copies, the prescribed-node owner/trace map and the masked nodal-load dual.
/// Flux/card strings and Newton rows have generous
/// per-element/per-port allowances; this admission estimate is not an RSS proof.
pub(in crate::solve::transient) fn storage_bytes(
    n: usize,
    elements: usize,
    steps: usize,
    regions: usize,
    boundaries: usize,
    max_newton_iterations: usize,
    max_krylov_iterations: usize,
) -> Result<usize, SolveRefusal> {
    let per_endpoint = n.checked_mul(3 * std::mem::size_of::<f64>())
        .and_then(|v| elements.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .and_then(|v| boundaries.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .and_then(|v| max_newton_iterations.checked_mul(256).and_then(|b| v.checked_add(b)))
        .and_then(|v| v.checked_add(std::mem::size_of::<EnthalpyStepSolution>() + 8))
        .ok_or_else(|| bad("enthalpy adjoint checkpoint allowance overflows"))?;
    let bytes = steps.checked_mul(per_endpoint)
        .and_then(|v| n.checked_mul(8 * (1 + 2 * REVERSE_RESTART_LIMIT + 36))
            .and_then(|b| v.checked_add(b)))
        .and_then(|v| v.checked_add(65_536))
        .and_then(|v| max_krylov_iterations.checked_mul(32)
            .and_then(|b| v.checked_add(b)))
        .and_then(|v| elements.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .and_then(|v| v.checked_add(super::storage_bytes(n, steps, regions).ok()?))
        .filter(|&v| isize::try_from(v).is_ok())
        .ok_or_else(|| bad("enthalpy adjoint storage allowance overflows"))?;
    Ok(bytes)
}

fn feedback_entries(n: usize, patches: usize) -> Result<usize, SolveRefusal> {
    if n == 0 || !(1..=64).contains(&patches) {
        return Err(bad("enthalpy radiation adjoints require 1..=64 nonempty patches"));
    }
    n.checked_mul(2).and_then(|v| v.checked_mul(patches))
        .filter(|&v| v.checked_mul(std::mem::size_of::<f64>())
            .is_some_and(|bytes| isize::try_from(bytes).is_ok()))
        .ok_or_else(|| bad("enthalpy radiation adjoint feedback allowance overflows"))
}

/// One active physical radiation revalidation: two retained nodal factors per
/// patch, temporary nodal fields, disjoint boundary faces and patch metadata.
/// Radiation declarations are immutable and need no per-checkpoint copies.
/// Aggregate emissivity/reservoir controls and result rows are charged once.
pub(in crate::solve::transient) fn radiation_storage_bytes(
    n: usize, elements: usize, patches: usize,
) -> Result<usize, SolveRefusal> {
    feedback_entries(n, patches)?.checked_mul(std::mem::size_of::<f64>())
        .and_then(|v| n.checked_mul(96).and_then(|b| v.checked_add(b)))
        .and_then(|v| elements.checked_mul(1024).and_then(|b| v.checked_add(b)))
        .and_then(|v| patches.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .and_then(|v| patches.checked_mul(RADIATION_PARAMETER_BYTES_PER_SURFACE)
            .and_then(|b| v.checked_add(b)))
        .and_then(|v| v.checked_add(2 * std::mem::size_of::<Vec<f64>>()))
        .filter(|&v| isize::try_from(v).is_ok())
        .ok_or_else(|| bad("enthalpy radiation adjoint workspace allowance overflows"))
}

fn clone_values(
    cx: &Cx<'_>, values: &[f64], deadline: Option<(Instant, f64)>,
) -> Result<Vec<f64>, SolveRefusal> {
    let mut copied = Vec::new();
    copied.try_reserve_exact(values.len())
        .map_err(|_| bad("enthalpy adjoint field allocation refused"))?;
    for (index, &value) in values.iter().enumerate() {
        if index % 256 == 0 { poll(cx, deadline)?; }
        copied.push(finite(value)?);
    }
    Ok(copied)
}

fn clone_rows<T: Clone>(
    cx: &Cx<'_>, rows: &[T], deadline: Option<(Instant, f64)>,
) -> Result<Vec<T>, SolveRefusal> {
    let mut copied = Vec::new();
    copied.try_reserve_exact(rows.len())
        .map_err(|_| bad("enthalpy adjoint telemetry allocation refused"))?;
    for row in rows {
        poll(cx, deadline)?;
        copied.push(row.clone());
    }
    Ok(copied)
}

fn zeros(length: usize) -> Result<Vec<f64>, SolveRefusal> {
    let mut values = Vec::new();
    values.try_reserve_exact(length)
        .map_err(|_| bad("enthalpy adjoint reverse allocation refused"))?;
    values.resize(length, 0.0);
    Ok(values)
}

/// Only the accepted fine grid is retained. Temperature checkpoints reuse the
/// existing objective/tie selector; full endpoints are consumed in reverse.
pub(in crate::solve::transient) struct EnthalpyTape {
    temperature: Tape,
    accepted: Vec<EnthalpyStepSolution>,
    initial_h: Vec<f64>,
    physical_tolerances_j: Vec<f64>,
}

impl EnthalpyTape {
    pub(in crate::solve::transient) fn new(
        cx: &Cx<'_>,
        initial_h: &[f64],
        initial_temperature: &[f64],
        steps: usize,
        deadline: Option<(Instant, f64)>,
    ) -> Result<Self, SolveRefusal> {
        poll(cx, deadline)?;
        if initial_h.is_empty() || initial_h.len() != initial_temperature.len() {
            return Err(bad("enthalpy adjoint initial h and temperature shapes differ"));
        }
        let initial_h = clone_values(cx, initial_h, deadline)?;
        let temperature = Tape::new(cx, initial_temperature, steps, deadline)?;
        let mut accepted = Vec::new();
        accepted.try_reserve_exact(steps)
            .map_err(|_| bad("enthalpy adjoint checkpoint allocation refused"))?;
        let mut physical_tolerances_j = Vec::new();
        physical_tolerances_j.try_reserve_exact(steps)
            .map_err(|_| bad("enthalpy adjoint physical-tolerance allocation refused"))?;
        Ok(Self { temperature, accepted, initial_h, physical_tolerances_j })
    }

    pub(in crate::solve::transient) fn record(
        &mut self,
        cx: &Cx<'_>,
        endpoint: &EnthalpyStepSolution,
        physical_tolerance_j: f64,
        deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        let n = self.temperature.n;
        if self.accepted.len().checked_add(1) != Some(self.temperature.states.len())
            || self.temperature.states.len() == self.temperature.limit
            || endpoint.specific_enthalpy_j_kg.len() != n
            || endpoint.temperature.len() != n
            || endpoint.liquid_mass_fraction.len() != n
            || self.physical_tolerances_j.len() != self.accepted.len()
            || !physical_tolerance_j.is_finite() || physical_tolerance_j < 0.0
        {
            return Err(bad("enthalpy adjoint checkpoint shape or count differs"));
        }
        let copied = EnthalpyStepSolution {
            pointwise_radiation: endpoint.pointwise_radiation.clone(),
            specific_enthalpy_j_kg: clone_values(cx, &endpoint.specific_enthalpy_j_kg, deadline)?,
            temperature: clone_values(cx, &endpoint.temperature, deadline)?,
            liquid_mass_fraction: clone_values(cx, &endpoint.liquid_mass_fraction, deadline)?,
            source_w: finite(endpoint.source_w)?,
            neumann_out_w: finite(endpoint.neumann_out_w)?,
            robin_out_w: finite(endpoint.robin_out_w)?,
            dirichlet_in_w: finite(endpoint.dirichlet_in_w)?,
            stored_energy_change_j: finite(endpoint.stored_energy_change_j)?,
            energy_residual_j: finite(endpoint.energy_residual_j)?,
            robin_fluxes: clone_rows(cx, &endpoint.robin_fluxes, deadline)?,
            contact_fluxes: clone_rows(cx, &endpoint.contact_fluxes, deadline)?,
            newton: NewtonReport {
                iterations: endpoint.newton.iterations,
                residual_norm: finite(endpoint.newton.residual_norm)?,
                converged: endpoint.newton.converged,
                history: clone_rows(cx, &endpoint.newton.history, deadline)?,
                diagnosis: endpoint.newton.diagnosis,
            },
        };
        // No fallible operation follows the temperature checkpoint: the endpoint
        // vector reserved its complete capacity in new().
        self.temperature.record(cx, &endpoint.temperature, deadline)?;
        self.accepted.push(copied);
        self.physical_tolerances_j.push(physical_tolerance_j);
        Ok(())
    }
}

fn dot(
    cx: &Cx<'_>, left: &[f64], right: &[f64], deadline: Option<(Instant, f64)>,
) -> Result<f64, SolveRefusal> {
    if left.len() != right.len() {
        return Err(bad("enthalpy adjoint contraction shapes differ"));
    }
    let mut total = 0.0;
    for (index, (&a, &b)) in left.iter().zip(right).enumerate() {
        if index % 256 == 0 { poll(cx, deadline)?; }
        total = finite(a.mul_add(b, total))?;
    }
    Ok(total)
}

fn source_contraction(
    cx: &Cx<'_>, bars: &[f64], source: &ScalarField, deadline: Option<(Instant, f64)>,
) -> Result<f64, SolveRefusal> {
    match source {
        ScalarField::Nodal(values) => dot(cx, bars, values, deadline),
        ScalarField::Uniform(value) => {
            let mut total = 0.0;
            for (index, &bar) in bars.iter().enumerate() {
                if index % 256 == 0 { poll(cx, deadline)?; }
                total = finite(bar.mul_add(*value, total))?;
            }
            Ok(total)
        }
    }
}

fn initial_slope(
    cx: &Cx<'_>,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    vertex: usize,
    h: f64,
    deadline: Option<(Instant, f64)>,
) -> Result<f64, SolveRefusal> {
    poll(cx, deadline)?;
    let curve = engine.phase_curve_at(vertex)
        .ok_or_else(|| bad("initial adjoint vertex has no enthalpy chart"))?;
    let knots = curve.knots();
    let derivative = |value| curve.temperature_derivative_at_specific_enthalpy(value)
        .map_err(|error| bad(format!("initial enthalpy objective derivative refused: {error}")));
    if let Some(index) = knots.iter().position(|knot| knot.specific_enthalpy_j_kg == h) {
        if index == 0 || index + 1 == knots.len() {
            return Err(failed(EnthalpyAdjointError::ChartEndpoint {
                vertex, specific_enthalpy_j_kg: h,
            }));
        }
        if derivative(knots[index - 1].specific_enthalpy_j_kg)? != derivative(h)? {
            return Err(failed(EnthalpyAdjointError::ChartKink {
                vertex, specific_enthalpy_j_kg: h,
            }));
        }
    }
    finite(derivative(h)?)
}

/// Fixed-chart/reference-geometry derivatives of the selected temperature
/// branch. Future fine-interval controls remain exactly zero.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(in crate::solve::transient) fn finish(
    cx: &Cx<'_>,
    spec: &ProjectSpec,
    engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    ends_s: &[f64],
    workload: Option<&Workload<'_>>,
    config: EnthalpyStepConfig,
    linear: LinearConfig,
    labels: &[u32],
    region_ids: &BTreeMap<String, u32>,
    goal_region: u32,
    storage: &TransientEnthalpy,
    radiation: Option<&radiation::LoweredRadiation>,
    airflow: Option<&airflow::ConjugatePath>,
    natural: &[buoyancy::NaturalLaw],
    pressure_pa: f64,
    tape: EnthalpyTape,
    objective: Objective,
    deadline: Option<(Instant, f64)>,
) -> Result<String, SolveRefusal> {
    poll(cx, deadline)?;
    let n = problem.mesh.vertex_count();
    let elements = problem.mesh.element_count();
    if ends_s.is_empty() || tape.temperature.n != n || tape.initial_h.len() != n
        || tape.accepted.len() != ends_s.len()
        || tape.physical_tolerances_j.len() != ends_s.len()
        || tape.temperature.states.len() != ends_s.len().checked_add(1)
            .ok_or_else(|| bad("enthalpy adjoint state count overflows"))?
        || labels.len() != elements
        || engine.vertex_material_ids().len() != n
        || engine.material_count() != storage.materials.len()
        || storage.materials.is_empty() || storage.materials.len() > 64
        || linear.max_iterations == 0 || linear.restart == 0
        || !(linear.tolerance > 0.0 && linear.tolerance < 1.0)
    {
        return Err(bad("enthalpy adjoint mesh, history, material or linear allowance differs"));
    }
    let linear = LinearConfig {
        restart: linear.restart.min(config.newton.linear_restart)
            .min(REVERSE_RESTART_LIMIT).min(n).min(linear.max_iterations),
        ..linear
    };
    if linear.restart == 0 {
        return Err(bad("enthalpy transpose needs a nonzero admitted restart"));
    }
    let iteration_limit = ends_s.len().checked_mul(linear.max_iterations)
        .ok_or_else(|| bad("enthalpy reverse iteration allowance overflows"))?;
    let radiation_workspace = radiation.map(|lowering|
        radiation_storage_bytes(n, elements, lowering.patches.len())).transpose()?.unwrap_or(0);
    let max_feedback_entries = radiation.map(|lowering|
        feedback_entries(n, lowering.patches.len())).transpose()?.unwrap_or(0);
    let convection_ports = convection::ports(natural, airflow)?;
    let convection_workspace = convection_storage_bytes(n, elements, convection_ports)?;
    let control_workspace = control_storage_bytes(spec, convection_ports)?;
    let material_workspace = material_control_storage_bytes(
        n, elements, spec.materials.as_deref().unwrap_or(&[]).len(),
        spec.interface_cards.as_deref().unwrap_or(&[]).len(),
    )?;
    let combined_feedback_entries = max_feedback_entries
        .checked_add(convection::feedback_entries(n, convection_ports)?)
        .ok_or_else(|| bad("combined enthalpy boundary-feedback allowance overflows"))?;
    let estimate = storage_bytes(
        n, elements, ends_s.len(), storage.materials.len(),
        problem.boundary.conditions().len(), config.max_newton_iterations,
        linear.max_iterations,
    )?.checked_add(radiation_workspace)
        .and_then(|v| v.checked_add(convection_workspace))
        .and_then(|v| v.checked_add(control_workspace))
        .and_then(|v| v.checked_add(material_workspace))
        .filter(|&v| isize::try_from(v).is_ok())
        .ok_or_else(|| bad("combined enthalpy adjoint storage allowance overflows"))?;
    let mut start = 0.0;
    for &end in ends_s {
        poll(cx, deadline)?;
        if !(end.is_finite() && end > start) {
            return Err(bad("enthalpy adjoint needs a finite strictly increasing fine grid"));
        }
        start = end;
    }
    for (vertex, (&h, &temperature)) in tape.initial_h.iter()
        .zip(&tape.temperature.states[0]).enumerate()
    {
        if vertex % 256 == 0 { poll(cx, deadline)?; }
        let curve = engine.phase_curve_at(vertex)
            .ok_or_else(|| bad("initial enthalpy adjoint has no assigned chart"))?;
        let knots = curve.knots();
        // Old h enters later endpoint balances affinely, so an internal slope
        // knot is admissible here. Validity endpoints never admit two-sided h
        // controls, including when the selected temperature is in the future.
        if h == knots[0].specific_enthalpy_j_kg
            || h == knots[knots.len() - 1].specific_enthalpy_j_kg
        {
            return Err(failed(EnthalpyAdjointError::ChartEndpoint {
                vertex, specific_enthalpy_j_kg: h,
            }));
        }
        let state = curve.state_at_specific_enthalpy(h)
            .map_err(|error| bad(format!("initial enthalpy adjoint state refused: {error}")))?;
        if state.temperature_k().to_bits() != temperature.to_bits() {
            return Err(bad("initial enthalpy adjoint temperature differs from its chart"));
        }
    }
    let mut selected_vertices = vec![false; n];
    let mut volumes = BTreeMap::<u32, f64>::new();
    for (element, (&id, tet)) in labels.iter().zip(&problem.mesh.complex().tets).enumerate() {
        if element % 256 == 0 { poll(cx, deadline)?; }
        let volume = volumes.entry(id).or_default();
        *volume = finite(*volume + problem.mesh.element_volume(element))?;
        if id == goal_region {
            for &vertex in tet { selected_vertices[vertex as usize] = true; }
        }
    }
    if volumes.len() > 64 {
        return Err(bad("enthalpy adjoint exceeds 64 regional source controls"));
    }
    let vertices: Vec<usize> = selected_vertices.iter().enumerate()
        .filter_map(|(vertex, &selected)| selected.then_some(vertex)).collect();
    let selected = select(cx, &tape.temperature, objective, &vertices, deadline)?;
    let mut names = Vec::new();
    let mut bases = Vec::new();
    for (&id, &volume) in &volumes {
        poll(cx, deadline)?;
        if volume <= 0.0 {
            return Err(bad("enthalpy source-control region has no positive volume"));
        }
        let mut matches = region_ids.iter().filter(|(_, found)| **found == id);
        let name = matches.next().map(|(name, _)| name.clone())
            .ok_or_else(|| bad("enthalpy source-control region has no declared name"))?;
        if matches.next().is_some() {
            return Err(bad("enthalpy source-control region has ambiguous names"));
        }
        let densities = volumes.keys().map(|&other|
            (other, if other == id { 1.0 / volume } else { 0.0 })).collect();
        // Absolute delivered-watt basis, including zero-baseline regions.
        let basis = PowerMap::regional_volumetric_source(problem.mesh, labels, &densities)
            .map_err(lower)?;
        let values = match basis {
            ScalarField::Uniform(value) => {
                let mut values = zeros(n)?;
                values.fill(finite(value)?);
                values
            }
            ScalarField::Nodal(values) => values,
        };
        if values.len() != n {
            return Err(bad("enthalpy source projector lost its nodal layout"));
        }
        names.push(name);
        bases.push(values);
    }
    let mut source_controls = Vec::new();
    source_controls.try_reserve_exact(ends_s.len())
        .map_err(|_| bad("enthalpy source-history allocation refused"))?;
    for _ in ends_s { source_controls.push(zeros(names.len())?); }
    let mut multipliers = zeros(ends_s.len())?;
    let mut densities = zeros(storage.materials.len())?;
    let surface_count = radiation.map_or(0, |lowering| lowering.patches.len());
    let mut emissivity_controls = zeros(surface_count)?;
    let mut reservoir_controls = zeros(surface_count)?;
    let mut carry = zeros(n)?;
    let mut dual_iterations = 0usize;
    let mut worst_dual = 0.0_f64;
    let mut convection_controls = convection::Controls::new(cx, spec, natural, airflow, deadline)?;
    let mut prescribed_controls = prescribed_controls::Controls::new(cx, spec, problem, deadline)?;
    let mut material_controls = material_controls::Controls::new(
        cx, spec, problem, interfaces, labels, region_ids, &tape.temperature.states[0], deadline,
    )?;
    let EnthalpyTape { temperature, mut accepted, initial_h, mut physical_tolerances_j } = tape;
    if selected.state == 0 {
        carry[selected.vertex] = initial_slope(
            cx, engine, selected.vertex, initial_h[selected.vertex], deadline,
        )?;
    }
    while accepted.len() > selected.state {
        poll(cx, deadline)?;
        accepted.pop();
        physical_tolerances_j.pop();
    }
    for index in (0..selected.state).rev() {
        poll(cx, deadline)?;
        let endpoint = accepted.pop()
            .ok_or_else(|| bad("enthalpy reverse lost an accepted endpoint"))?;
        let physical_tolerance_j = physical_tolerances_j.pop()
            .ok_or_else(|| bad("enthalpy reverse lost the original physical residual tolerance"))?;
        let old_h = if index == 0 { &initial_h } else {
            &accepted.last().ok_or_else(|| bad("enthalpy reverse lost its old state"))?
                .specific_enthalpy_j_kg
        };
        for (vertex, &h) in endpoint.specific_enthalpy_j_kg.iter().enumerate() {
            if vertex % 256 == 0 { poll(cx, deadline)?; }
            let state = engine.phase_curve_at(vertex)
                .ok_or_else(|| bad("accepted enthalpy adjoint has no assigned chart"))?
                .state_at_specific_enthalpy(h)
                .map_err(|error| bad(format!("accepted enthalpy adjoint state refused: {error}")))?;
            if state.temperature_k().to_bits() != endpoint.temperature[vertex].to_bits()
                || state.liquid_mass_fraction().to_bits()
                    != endpoint.liquid_mass_fraction[vertex].to_bits()
            {
                return Err(bad("accepted enthalpy adjoint fields differ from their chart"));
            }
        }
        let start = if index == 0 { 0.0 } else { ends_s[index - 1] };
        let dt = ends_s[index] - start;
        let source = workload.map(|history|
            history.source(cx, problem.mesh, labels, start, deadline)).transpose()?;
        let source_problem = ConductionProblem {
            source: source.as_ref().unwrap_or(problem.source), ..problem
        };
        let feedback = convection::bind(
            cx, source_problem, &endpoint.temperature, natural, airflow, pressure_pa, deadline,
        )?;
        let step_problem = ConductionProblem {
            boundary: feedback.as_ref().map_or(problem.boundary, |bound| &bound.boundary),
            ..source_problem
        };
        let mut endpoint_config = config;
        if let Some(bound) = &feedback {
            // The physical target belongs to the accepted forward solve.
            // Reconstructing a target from old h with endpoint-frozen boundary
            // values could loosen that contract. Recheck it explicitly here.
            let physical = engine.evaluate_endpoint(
                cx, step_problem, interfaces, old_h, &endpoint.specific_enthalpy_j_kg,
                dt, config, radiation.map(|lowering| lowering.patches.as_slice()),
            ).map_err(|error| failed(EnthalpyAdjointError::Enthalpy(error)))?;
            if !physical.residual_norm_j.is_finite()
                || physical.residual_norm_j > physical_tolerance_j
                || !physical.energy_residual_j.is_finite()
                || physical.energy_residual_j.abs() > config.energy_tolerance_j
            {
                return Err(bad("coupled enthalpy adjoint endpoint failed its original physical residual or energy gate"));
            }
            bound.check_air_energy(&physical, dt, config.energy_tolerance_j)?;
            // These independent gates have enforced the retained physical
            // target. The shared binder reconstructs its initial residual
            // with endpoint-frozen convection, so that different scale must
            // not introduce a tighter target for this accepted endpoint.
            endpoint_config.newton.absolute_tolerance = config
                .newton.absolute_tolerance.max(physical_tolerance_j);
            poll(cx, deadline)?;
        }
        let mut step = if let Some(lowering) = radiation {
            // The original convection boundary and physical old h are passed
            // unchanged. The shared binder rebuilds radiation at this endpoint
            // and retains its complete mean-dependent feedback in the tangent.
            BoundEndpoint::Radiative(engine.linearize_accepted_parameters_with_ambient_radiation(
                cx, step_problem, interfaces, old_h, dt, endpoint_config,
                &lowering.patches, endpoint, max_feedback_entries,
            ).map_err(failed)?)
        } else {
            BoundEndpoint::Prescribed(engine.linearize_accepted(
                cx, step_problem, interfaces, old_h, dt, endpoint_config, endpoint,
            ).map_err(failed)?)
        };
        if let Some(bound) = &feedback {
            let regions: Vec<&str> = bound.names.iter().map(String::as_str).collect();
            step = match step {
                BoundEndpoint::Prescribed(step) => BoundEndpoint::Prescribed(
                    step.with_robin_feedback(cx, &regions, &bound.references,
                        &bound.wall_slopes, combined_feedback_entries).map_err(failed)?),
                BoundEndpoint::Radiative(step) => BoundEndpoint::Radiative(
                    step.with_robin_feedback(cx, &regions, &bound.references,
                        &bound.wall_slopes, combined_feedback_entries).map_err(failed)?),
            };
        }
        let expected = &temperature.states[index + 1];
        for (vertex, (&actual, &expected)) in step.transport().primal().temperature.iter().zip(expected).enumerate() {
            if vertex % 256 == 0 { poll(cx, deadline)?; }
            if actual.to_bits() != expected.to_bits() {
                return Err(bad(format!(
                    "enthalpy endpoint revalidation changed state {}, vertex {vertex}", index + 1,
                )));
            }
        }
        if index + 1 == selected.state {
            let mut seed = zeros(n)?;
            seed[selected.vertex] = 1.0;
            carry = step.transport().temperature_pullback(cx, &seed).map_err(failed)?;
        }
        // This seed is an h cotangent. Multiplying history carry by dT/dh again
        // would incorrectly erase energy sensitivity on a latent plateau.
        let gradient = match &step {
            BoundEndpoint::Prescribed(step) => step.pullback(cx, &carry, linear).map_err(failed)?,
            BoundEndpoint::Radiative(step) => {
                // The radiative owner contracts its parameters and reference
                // storage with ONE checked physical transpose multiplier.
                let gradient = step.pullback(cx, &carry, linear).map_err(failed)?;
                let lowering = radiation.ok_or_else(|| bad("radiative enthalpy controls lost their source"))?;
                if gradient.radiation_regions.len() != surface_count
                    || gradient.emissivities.len() != surface_count
                    || gradient.reservoir_temperatures.len() != surface_count
                {
                    return Err(bad("radiative enthalpy parameter layout differs"));
                }
                for (surface, patch) in lowering.patches.iter().enumerate() {
                    poll(cx, deadline)?;
                    if gradient.radiation_regions[surface] != patch.region() {
                        return Err(bad("radiative enthalpy control order differs from its bound surfaces"));
                    }
                    emissivity_controls[surface] = finite(
                        emissivity_controls[surface] + gradient.emissivities[surface],
                    )?;
                    reservoir_controls[surface] = finite(
                        reservoir_controls[surface] + gradient.reservoir_temperatures[surface],
                    )?;
                }
                gradient.storage
            }
        };
        if !gradient.transport.relative_residual.is_finite()
            || gradient.transport.relative_residual >= linear.tolerance
            || gradient.material_reference_density.len() != densities.len()
        {
            return Err(bad("enthalpy transpose failed its residual or material-layout gate"));
        }
        // Material and interface parameters contract the SAME checked physical
        // multiplier. The core retains conductivity at this accepted T and
        // phase state; its element bars already include dt exactly once.
        let conductivity = step.transport().conductivity_scale_pullback(
            cx, &gradient.transport.adjoint, elements,
        ).map_err(failed)?;
        let nodal_load_adjoint = step.transport().nodal_load_pullback(
            cx, &gradient.transport.adjoint,
        ).map_err(failed)?;
        material_controls.add(
            cx, step_problem, interfaces, expected, &conductivity,
            &nodal_load_adjoint, deadline,
        )?;
        match (convection_controls.as_mut(), feedback.as_ref()) {
            (Some(controls), Some(bound)) => controls.add(
                cx, spec, step_problem, expected, bound, &nodal_load_adjoint, deadline,
            )?,
            (None, None) => {},
            _ => return Err(bad("enthalpy physical boundary controls lost their bound endpoint")),
        }
        prescribed_controls.add(cx, &gradient.transport.prescribed_temperature, deadline)?;
        for (region, basis) in bases.iter().enumerate() {
            source_controls[index][region] = dot(
                cx, &gradient.transport.source_density, basis, deadline,
            )?;
        }
        multipliers[index] = source_contraction(
            cx, &gradient.transport.source_density, step_problem.source, deadline,
        )?;
        for (total, &value) in densities.iter_mut().zip(&gradient.material_reference_density) {
            *total = finite(*total + value)?;
        }
        dual_iterations = dual_iterations.checked_add(gradient.transport.iterations)
            .ok_or_else(|| bad("enthalpy transpose iteration count overflows"))?;
        worst_dual = worst_dual.max(gradient.transport.relative_residual);
        carry = gradient.transport.previous_specific_enthalpy;
    }
    let mut initial_rows = Vec::new();
    let mut initial_by_material = zeros(storage.materials.len())?;
    let mut initial_uniform = 0.0;
    for (vertex, &value) in carry.iter().enumerate() {
        if vertex % 256 == 0 { poll(cx, deadline)?; }
        initial_rows.push(number(value)?);
        initial_uniform = finite(initial_uniform + value)?;
        let material = engine.vertex_material_ids()[vertex];
        let regional = initial_by_material.get_mut(material)
            .ok_or_else(|| bad("initial enthalpy adjoint material index differs"))?;
        *regional = finite(*regional + value)?;
    }
    let mut regional_storage = Vec::new();
    let mut density_multiplier = 0.0;
    for (material, row) in storage.materials.iter().enumerate() {
        poll(cx, deadline)?;
        density_multiplier = finite(
            row.reference_density.value.mul_add(densities[material], density_multiplier),
        )?;
        regional_storage.push(format!(
            "{{\"region\":{},\"material_index\":{},\"reference_density_kg_m3\":{},\"dtemperature_dinitial_specific_enthalpy_k_kg_per_j\":{},\"dtemperature_dreference_density_k_m3_per_kg\":{}}}",
            json_string(&row.region), material, number(row.reference_density.value)?,
            number(initial_by_material[material])?, number(densities[material])?,
        ));
    }
    let mut totals = zeros(names.len())?;
    let mut intervals = Vec::new();
    let mut start = 0.0;
    for (index, &end) in ends_s.iter().enumerate() {
        poll(cx, deadline)?;
        let mut rows = Vec::new();
        for (region, name) in names.iter().enumerate() {
            let derivative = source_controls[index][region];
            totals[region] = finite(totals[region] + derivative)?;
            rows.push(format!(
                "{{\"region\":{},\"dtemperature_ddelivered_watts_k_per_w\":{}}}",
                json_string(name), number(derivative)?,
            ));
        }
        intervals.push(format!(
            "{{\"step\":{},\"start_s\":{},\"end_s\":{},\"dtemperature_dsource_multiplier_k\":{},\"regional_power\":[{}]}}",
            index + 1, number(start)?, number(end)?, number(multipliers[index])?, rows.join(","),
        ));
        start = end;
    }
    let regional_power = names.iter().zip(&totals).map(|(name, &value)| Ok(format!(
        "{{\"region\":{},\"dtemperature_duniform_added_watts_k_per_w\":{}}}",
        json_string(name), number(value)?,
    ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
    let selected_time = if selected.state == 0 { 0.0 } else { ends_s[selected.state - 1] };
    let radiation_fragment = radiation.map(|lowering| -> Result<String, SolveRefusal> {
        if lowering.surface_labels().len() != surface_count {
            return Err(bad("radiative enthalpy surface identities differ from the bound patches"));
        }
        let mut rows = Vec::new();
        rows.try_reserve_exact(surface_count)
            .map_err(|_| bad("radiative enthalpy control-row allocation refused"))?;
        for (surface, ((name, target), patch)) in lowering.surface_labels()
            .zip(&lowering.patches).enumerate()
        {
            poll(cx, deadline)?;
            if target != patch.region() {
                return Err(bad("radiative enthalpy surface identity order differs"));
            }
            let emissivity = patch.emissivity();
            let directions = if emissivity.value() == 1.0 {
                "nonpositive-at-emissivity-one"
            } else {
                "two-sided-within-zero-exclusive-one-inclusive"
            };
            rows.push(format!(
                "{{\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{},\"reservoir_temperature_k\":{},\"dtemperature_demissivity_k\":{},\"dtemperature_dreservoir_temperature_k_per_k\":{},\"emissivity_feasible_directions\":{}}}",
                json_string(name), json_string(target),
                json_string(&emissivity.card_identity().to_hex()),
                json_string(emissivity.material_state()),
                json_string(&emissivity.receipt().content_hash().to_hex()),
                number(emissivity.value())?, number(emissivity.temperature_k())?,
                number(patch.ambient_temperature_k())?,
                number(emissivity_controls[surface])?, number(reservoir_controls[surface])?,
                json_string(directions),
            ));
        }
        Ok(format!(
            ",\"radiation\":{{\"method\":\"complete-area-mean-gray-enthalpy-tangent\",\"patch_count\":{},\"max_feedback_entries\":{max_feedback_entries},\"additional_workspace_bytes\":{radiation_workspace},\"parameters\":\"absolute-per-surface-emissivity-and-reservoir-temperature-throughout-trajectory\",\"surfaces\":[{}]}}",
            lowering.patches.len(), rows.join(","),
        ))
    }).transpose()?.unwrap_or_default();
    let material_fragment = material_controls.fragment(selected.state)?;
    let prescribed_fragment = prescribed_controls.fragment(selected.state)?;
    let controls_fragment = convection_controls.as_ref()
        .map(|controls| controls.fragment(selected.state)).transpose()?.unwrap_or_default();
    let convection_fragment = if convection_ports == 0 { String::new() } else { format!(
        ",\"convection\":{{\"method\":\"complete-natural-and-streamwise-air-enthalpy-feedback\",\"natural_ports\":{},\"air_ports\":{},\"max_combined_feedback_entries\":{combined_feedback_entries},\"additional_workspace_bytes\":{convection_workspace},\"control_workspace_bytes\":{control_workspace},\"physical_endpoint_rechecks\":{},\"residual_target\":\"original-accepted-forward-physical-tolerance\",\"air_energy\":\"actual-air-enthalpy-gain-with-separate-natural-and-radiation-reservoirs\"{controls_fragment}}}",
        natural.len(), airflow.map_or(0, |path| path.segments.len()), selected.state,
    ) };
    poll(cx, deadline)?;
    Ok(format!(
        "{{\"schema\":\"fs-cli-native-enthalpy-adjoint-v1\",\"output\":{},\"method\":\"discrete-total-enthalpy-backward-euler-reverse\",\"value_k\":{},\"time_s\":{},\"state_index\":{},\"active_vertex\":{},\"active_ties\":{},\"unique_sampled_maximum\":{},\"dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j\":{},\"dtemperature_dinitial_specific_enthalpy_k_kg_per_j\":[{}],\"dtemperature_dreference_density_multiplier_k\":{},\"regional_storage\":[{}],\"regional_power\":[{}],\"intervals\":[{}],\"revalidated_solid_endpoints\":{},\"adjoint_krylov_iterations\":{},\"maximum_dual_relative_residual\":{},\"max_reverse_steps\":{},\"max_krylov_iterations_per_linear_solve\":{},\"max_restart_columns_per_linear_solve\":{},\"max_adjoint_krylov_iterations\":{},\"additional_storage_estimate_bytes\":{}{radiation_fragment}{convection_fragment}{material_fragment}{prescribed_fragment},\"authority\":\"Estimated\",\"scope\":{}}}",
        json_string(objective.name()), number(selected.value)?, number(selected_time)?,
        selected.state, selected.vertex, selected.ties, selected.ties == 1,
        number(initial_uniform)?, initial_rows.join(","), number(density_multiplier)?,
        regional_storage.join(","), regional_power.join(","), intervals.join(","),
        selected.state, dual_iterations, number(worst_dual)?, ends_s.len(),
        linear.max_iterations, linear.restart, iteration_limit, estimate, json_string(SCOPE),
    ))
}

const SCOPE: &str = "Selected final or earliest sampled-peak nodal branch on the accepted fine grid, including the initial state for sampled peaks. Exact ties have no unique maximum derivative. Only the selected prefix is reversed; future interval controls are zero. Each consumed endpoint is independently revalidated through the physical heterogeneous enthalpy balance, energy gate, chart smoothness and complete k(T)/phase-conductivity tangent; no nonlinear endpoint is re-solved or replaced. Initial controls are additive nodal or regional specific-enthalpy offsets, not initial temperatures. Native input edits must preserve initial-temperature/enthalpy consistency or use from-regional-enthalpy. Reference-density controls change frozen storage mass while holding initial enthalpy, charts, conductivity laws and geometry fixed; the multiplier scales all regional reference densities together. Regional power controls add absolute delivered watts through the same P1 source projector, including zero-baseline intervals; sum interval derivatives for one declared workload segment. Prescribed Neumann/Robin data, phase-law parameters and the time grid are fixed. Regional conductivity controls multiply each region's entire effective tensor/curve by one positive scale at one, retaining its declared phase multiplier; contact controls scale every retained face resistance of one named matching interface together at one. Their trajectory derivatives use the same complete physical multiplier and preserve separate region/interface/card identities. These normalized material controls are not absolute tensor-entry or resistance derivatives, new card claims, uncertainty propagation, contact-geometry changes, or derivatives of constitutive knots. Declared ambient radiation uses its complete nonlinear area-mean state feedback in the same h/source/reference-density transpose; the report also differentiates independent absolute surface emissivities and reservoir temperatures held throughout the selected trajectory. Each surface retains its name, target, sourced nominal value, card and usage-receipt identities; sharing a card does not merge controls. Query temperature, selected source claim/support and time grid remain fixed. Emissivity is restricted to (0,1], with only nonpositive feasible increments at one, and reservoirs require positive absolute kelvin. These parameter derivatives are not uncertainty propagation or derivatives of immutable card content. An isothermal plateau has zero direct temperature sensitivity while h history carry remains energetic. Chart/domain endpoints and slope-changing branches refuse classical derivatives. Natural convection retains the complete signed wall and film-density derivative of the in-domain Churchill-Chu vertical-plate law. Quasi-steady air retains the complete upstream-to-downstream reference-temperature feedback of each independent branch at fixed flow and transport properties. Both compose with radiation, phase conductivity and energetic history in the same physical transpose; each coupled endpoint also passes its original forward residual target and actual-air energy gate. Natural-reservoir and independently prescribed air-inlet temperatures are absolute kelvin controls held throughout the selected trajectory. Each natural ambient control includes both its explicit Robin reference and its full film-density/coefficient partial. Independent fan-bank controls differentiate absolute speed ratios through the same retained hydraulic point, capacity rates, smooth in-domain correlation coefficients and downstream air transport; unsupported bank coordinates retain their identities and reasons instead of partial or zero derivatives. Air transport properties, pressure, the hydraulic density envelope, fixed loss geometry and fan-curve data remain fixed. At an initial sampled peak every boundary control has zero influence because no future endpoint contributes. No apparent heat capacity, moving geometry, density-driven expansion, continuous-time maximum, uncertainty or physical-validation certificate is supplied.";
