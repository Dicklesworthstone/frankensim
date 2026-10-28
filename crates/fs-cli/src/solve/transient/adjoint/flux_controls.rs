//! Original prescribed fluxes share the complete accepted thermal multiplier.
//! Both storage coordinates use the conduction owner's exact P1 trace basis.

use std::time::Instant;

use fs_conduction::adjoint::robin::MeanFluxFeedback;
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc};
use fs_exec::Cx;
use fs_project::{ProjectSpec, ThermalBoundaryCondition};

use super::{SolveRefusal, bad, finite, poll};
use super::super::{json_string, lower, number};

struct Control {
    target: String,
    flux_w_m2: f64,
    derivative: f64,
}

/// Retained trace factors, fixed-node indices, temporary zero feedback matrix
/// and bounded output metadata. This is an admission estimate, not peak RSS.
pub(in crate::solve::transient) fn storage_bytes(spec: &ProjectSpec, nodes: usize) -> Result<usize, SolveRefusal> {
    let count = spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
        .map_or(0, |setup| setup.boundaries.iter().filter(|boundary|
            matches!(boundary.condition, ThermalBoundaryCondition::HeatFlux { .. })).count());
    if count == 0 { return Ok(0); }
    if count > 64 { return Err(bad("transient heat-flux adjoints exceed 64 boundary controls")); }
    nodes.checked_mul(count).and_then(|v| v.checked_mul(2))
        .and_then(|v| v.checked_add(nodes))
        .and_then(|v| count.checked_mul(count).and_then(|b| v.checked_add(b)))
        .and_then(|v| v.checked_mul(std::mem::size_of::<f64>()))
        .and_then(|v| count.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .and_then(|v| v.checked_add(512))
        .filter(|&v| isize::try_from(v).is_ok())
        .ok_or_else(|| bad("transient heat-flux control workspace overflows"))
}

pub(super) struct Controls {
    rows: Vec<Control>,
    trace: Option<MeanFluxFeedback>,
    fixed: Vec<usize>,
    steps: usize,
}

impl Controls {
    pub(super) fn new(
        cx: &Cx<'_>, spec: &ProjectSpec, problem: ConductionProblem<'_>,
        deadline: Option<(Instant, f64)>,
    ) -> Result<Self, SolveRefusal> {
        poll(cx, deadline)?;
        let setup = spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
            .ok_or_else(|| bad("transient heat-flux controls lost their declaration"))?;
        storage_bytes(spec, problem.mesh.vertex_count())?;
        let mut rows = Vec::new();
        for boundary in &setup.boundaries {
            poll(cx, deadline)?;
            let ThermalBoundaryCondition::HeatFlux { outward_flux }
                = &boundary.condition else { continue; };
            let region = problem.boundary.region_names().iter()
                .position(|name| name == &boundary.target)
                .ok_or_else(|| bad("prescribed heat-flux control has no physical trace"))?;
            let ThermalBc::Neumann { outward_flux: ScalarField::Uniform(actual) }
                = problem.boundary.conditions()[region] else {
                    return Err(bad("prescribed heat-flux control changed its physical law"));
                };
            if !actual.is_finite() || actual.to_bits() != outward_flux.value.to_bits() {
                return Err(bad("prescribed heat-flux control differs from its declaration"));
            }
            rows.try_reserve(1)
                .map_err(|_| bad("prescribed heat-flux metadata allocation refused"))?;
            rows.push(Control { target: boundary.target.clone(), flux_w_m2: actual, derivative: 0.0 });
        }
        if rows.is_empty() {
            return Ok(Self { rows, trace: None, fixed: Vec::new(), steps: 0 });
        }
        let names: Vec<_> = rows.iter().map(|row| row.target.as_str()).collect();
        let entries = problem.mesh.vertex_count().checked_mul(names.len())
            .and_then(|v| v.checked_mul(2))
            .ok_or_else(|| bad("prescribed heat-flux trace allowance overflows"))?;
        // Zero feedback retains only the existing owner's G_i basis. The
        // physical state feedback is already present in the supplied dual.
        let trace = MeanFluxFeedback::new(
            cx, problem.mesh, problem.boundary, &names,
            &vec![0.0; names.len() * names.len()], entries,
        ).map_err(lower)?;
        let fixed = problem.boundary.dirichlet().iter().map(|&(vertex, _)| vertex).collect();
        poll(cx, deadline)?;
        Ok(Self { rows, trace: Some(trace), fixed, steps: 0 })
    }

    /// `load_weights` already contains the full timestep convention: capacity
    /// supplies its watt-coordinate dual; enthalpy supplies dt times its dual.
    /// Positive outward flux enters the residual positively, hence the minus.
    pub(super) fn add(
        &mut self, cx: &Cx<'_>, load_weights: &[f64], deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        if let Some(trace) = &self.trace {
            if load_weights.len() != trace.vertex_count()
                || self.fixed.iter().any(|&vertex| load_weights.get(vertex) != Some(&0.0))
            {
                return Err(bad("heat-flux controls require the complete masked nodal-load multiplier"));
            }
            let values = trace.flux_pullback(cx, load_weights).map_err(lower)?;
            if values.len() != self.rows.len() {
                return Err(bad("heat-flux controls lost their physical trace ordering"));
            }
            for (row, value) in self.rows.iter_mut().zip(values) {
                row.derivative = finite(row.derivative - value)?;
            }
        }
        self.steps = self.steps.checked_add(1)
            .ok_or_else(|| bad("heat-flux reverse step count overflows"))?;
        poll(cx, deadline)
    }

    pub(super) fn fragment(&self, expected_steps: usize) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps {
            return Err(bad("heat-flux controls lost the selected reverse prefix"));
        }
        if self.rows.is_empty() { return Ok(String::new()); }
        let rows = self.rows.iter().map(|row| Ok(format!(
            "{{\"target\":{},\"flux_w_m2\":{},\"dtemperature_doutward_flux_k_per_w_m2\":{}}}",
            json_string(&row.target), number(row.flux_w_m2)?, number(row.derivative)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        Ok(format!(
            ",\"prescribed_heat_flux\":{{\"control_steps\":{},\"parameter_scope\":\"independent-uniform-outward-heat-fluxes-held-throughout-selected-trajectory\",\"scope\":\"Absolute original Neumann flux densities in W/m2, positive outward, on unchanged named P1 traces. The complete physical thermal multiplier includes storage and all admitted feedback; no additional solve or timestep factor is introduced. Prescribed-node equation weights are zero. Initial state, geometry and time grid stay fixed; an initial sampled peak has zero future flux sensitivity.\",\"boundaries\":[{}]}}",
            self.steps, rows.join(","),
        ))
    }
}
