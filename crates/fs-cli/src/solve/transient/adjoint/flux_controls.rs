//! Prescribed fluxes and surface powers share the accepted thermal multiplier.
//! Both storage coordinates use the conduction owner's exact P1 trace basis.

use std::time::Instant;

use fs_conduction::adjoint::{PrescribedFluxControls, PrescribedFluxGradient};
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc};
use fs_exec::Cx;
use fs_project::{ProjectSpec, ThermalBoundaryCondition, TransientBoundaryTarget};

use super::{SolveRefusal, bad, finite, poll};
use super::super::{json_string, lower, number};

struct Control {
    target: String,
    flux_w_m2: f64,
    derivative: f64,
    surface_area_m2: Option<f64>,
}

struct PowerControl {
    trace: usize,
    source_index: usize,
    watts: f64,
    duty: f64,
}

/// Sparse original-face factors, fixed-node indices and bounded output metadata.
/// This is an admission estimate, not peak RSS; no dense zero-feedback basis.
pub(in crate::solve::transient) fn storage_bytes(spec: &ProjectSpec, nodes: usize,
    boundary_faces: usize) -> Result<usize, SolveRefusal> {
    let fluxes = spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
        .map_or(0, |setup| setup.boundaries.iter().filter(|boundary|
            matches!(boundary.condition, ThermalBoundaryCondition::HeatFlux { .. })).count());
    let heat = crate::solve::surface_heat(spec)?;
    let powers = spec.power.iter().flatten()
        .filter(|row| heat.names.contains(&row.region)).count();
    let count = fluxes.checked_add(heat.watts.len())
        .ok_or_else(|| bad("transient boundary-load control count overflows"))?;
    let metadata = count.checked_add(powers)
        .ok_or_else(|| bad("transient boundary-load metadata count overflows"))?;
    if count == 0 { return Ok(0); }
    if count > 64 { return Err(bad("transient boundary-load adjoints exceed 64 boundary controls")); }
    nodes.checked_mul(std::mem::size_of::<usize>())
        .and_then(|v| boundary_faces.checked_mul(128).and_then(|b| v.checked_add(b)))
        .and_then(|v| metadata.checked_mul(4096).and_then(|b| v.checked_add(b)))
        .and_then(|v| v.checked_add(512))
        .filter(|&v| isize::try_from(v).is_ok())
        .ok_or_else(|| bad("transient heat-flux control workspace overflows"))
}

pub(super) struct Controls {
    rows: Vec<Control>,
    powers: Vec<PowerControl>,
    trace: Option<PrescribedFluxControls>,
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
        storage_bytes(spec, problem.mesh.vertex_count(), problem.mesh.boundary().len())?;
        let mut rows = Vec::new();
        for boundary in &setup.boundaries {
            poll(cx, deadline)?;
            let ThermalBoundaryCondition::HeatFlux { outward_flux }
                = &boundary.condition else { continue; };
            // The authored history replaces this unused static coordinate.
            // Its held flux still enters every original forward/reverse row.
            if setup.transient.as_ref().is_some_and(|time|
                time.boundary_schedules.iter().any(|row| row.target == boundary.target
                    && row.kind == TransientBoundaryTarget::HeatFlux))
            {
                continue;
            }
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
            // An area-mean ambient capacity patch admits explicit absence of
            // nonradiative exchange, not an independently variable Neumann
            // load. Retain that zero-flux identity without inventing a control
            // whose nonzero perturbations the original forward model refuses.
            if actual == 0.0
                && setup.transient.as_ref().is_some_and(|time| time.enthalpy.is_none())
                && setup.radiation.as_ref().is_some_and(|radiation|
                    radiation.spatial_model != Some(fs_project::RadiationSpatialModel::Pointwise)
                        && radiation.surfaces.iter().any(|surface| surface.target == boundary.target))
            {
                continue;
            }
            rows.try_reserve(1)
                .map_err(|_| bad("prescribed heat-flux metadata allocation refused"))?;
            rows.push(Control { target: boundary.target.clone(), flux_w_m2: actual,
                derivative: 0.0, surface_area_m2: None });
        }
        let heat = crate::solve::surface_heat(spec)?;
        for (target, &delivered) in &heat.watts {
            poll(cx, deadline)?;
            if rows.iter().any(|row| row.target == *target) {
                return Err(bad("surface power cannot also own a prescribed flux control"));
            }
            let region = problem.boundary.region_names().iter()
                .position(|name| name == target)
                .ok_or_else(|| bad("surface-power control has no retained boundary owner"))?;
            let ThermalBc::Neumann { outward_flux: ScalarField::Uniform(actual) }
                = problem.boundary.conditions()[region] else {
                    return Err(bad("surface-power control lost its uniform inward Neumann load"));
                };
            let mut area = 0.0;
            for (slot, face) in problem.mesh.boundary().iter().enumerate() {
                if slot % 512 == 0 { poll(cx, deadline)?; }
                if problem.boundary.region_for(slot) == Some(region) {
                    area = finite(area + face.area)?;
                }
            }
            // Use the primal's physical trace and its exact watts/area
            // normalization. A zero baseline never becomes a divisor.
            if area <= 0.0 || !actual.is_finite() || actual > 0.0
                || actual.to_bits() != (-delivered / area).to_bits()
            {
                return Err(bad("surface-power control differs from its declared physical load"));
            }
            rows.try_reserve(1)
                .map_err(|_| bad("surface-power trace allocation refused"))?;
            rows.push(Control { target: target.clone(), flux_w_m2: actual,
                derivative: 0.0, surface_area_m2: Some(area) });
        }
        let mut powers = Vec::new();
        for (source_index, power) in spec.power.iter().flatten().enumerate() {
            poll(cx, deadline)?;
            if !heat.names.contains(&power.region) { continue; }
            let trace = rows.iter().position(|row| row.target == power.region)
                .ok_or_else(|| bad("declared surface power lost its retained trace"))?;
            if !power.watts.value.is_finite() || power.watts.value < 0.0
                || !power.duty.is_finite() || !(0.0..=1.0).contains(&power.duty)
            {
                return Err(bad("surface-power controls require nonnegative watts and duty in [0,1]"));
            }
            powers.try_reserve(1)
                .map_err(|_| bad("surface-power metadata allocation refused"))?;
            powers.push(PowerControl { trace, source_index, watts: power.watts.value,
                duty: power.duty });
        }
        if rows.is_empty() {
            return Ok(Self { rows, powers, trace: None, steps: 0 });
        }
        let names: Vec<_> = rows.iter().map(|row| row.target.as_str()).collect();
        let trace = PrescribedFluxControls::new(cx, problem.mesh, problem.boundary,
            &names, problem.mesh.boundary().len()).map_err(lower)?;
        poll(cx, deadline)?;
        Ok(Self { rows, powers, trace: Some(trace), steps: 0 })
    }

    /// `load_weights` already contains the full timestep convention: capacity
    /// supplies its watt-coordinate dual; enthalpy supplies dt times its dual.
    /// Positive outward flux enters the residual positively, hence the minus.
    pub(super) fn add(
        &mut self, cx: &Cx<'_>, load_weights: &[f64], deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        if let Some(trace) = &self.trace {
            // The native owner checks length, finiteness and zero fixed-node
            // weights, and includes the outward-load minus exactly once.
            let values = trace.pullback(cx, load_weights).map_err(lower)?;
            if values.len() != self.rows.len() {
                return Err(bad("heat-flux controls lost their physical trace ordering"));
            }
            for (row, value) in self.rows.iter_mut().zip(values) {
                let PrescribedFluxGradient::Normal(value) = value
                    else { return Err(bad("project scalar heat flux changed its native coordinate")); };
                row.derivative = finite(row.derivative + value)?;
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
        let rows = self.rows.iter().filter(|row| row.surface_area_m2.is_none()).map(|row| Ok(format!(
            "{{\"target\":{},\"flux_w_m2\":{},\"dtemperature_doutward_flux_k_per_w_m2\":{}}}",
            json_string(&row.target), number(row.flux_w_m2)?, number(row.derivative)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        if rows.is_empty() { return Ok(String::new()); }
        Ok(format!(
            ",\"prescribed_heat_flux\":{{\"control_steps\":{},\"parameter_scope\":\"independent-uniform-outward-heat-fluxes-held-throughout-selected-trajectory\",\"scope\":\"Absolute original Neumann flux densities in W/m2, positive outward, on unchanged named P1 traces. The complete physical thermal multiplier includes storage and all admitted feedback; no additional solve or timestep factor is introduced. Prescribed-node equation weights are zero. Initial state, geometry and time grid stay fixed; an initial sampled peak has zero future flux sensitivity.\",\"boundaries\":[{}]}}",
            self.steps, rows.join(","),
        ))
    }

    pub(super) fn surface_power_fragment(&self, expected_steps: usize) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps {
            return Err(bad("surface-power controls lost the selected reverse prefix"));
        }
        if self.powers.is_empty() { return Ok(String::new()); }
        let rows = self.powers.iter().map(|power| {
            let trace = &self.rows[power.trace];
            let area = trace.surface_area_m2
                .ok_or_else(|| bad("surface-power control lost its physical area"))?;
            let delivered = finite(-trace.derivative / area)?;
            Ok(format!(
                "{{\"target\":{},\"source_index\":{},\"watts\":{},\"duty\":{},\"area_m2\":{},\"dtemperature_ddelivered_watts_k_per_w\":{},\"dtemperature_dwatts_k_per_w\":{},\"dtemperature_dduty_k\":{}}}",
                json_string(&trace.target), power.source_index, number(power.watts)?,
                number(power.duty)?, number(area)?, number(delivered)?,
                number(finite(delivered * power.duty)?)?,
                number(finite(delivered * power.watts)?)?,
            ))
        }).collect::<Result<Vec<_>, SolveRefusal>>()?;
        Ok(format!(
            ",\"surface_power\":{{\"control_steps\":{},\"parameter_scope\":\"independent-declared-watts-and-duty-held-throughout-selected-trajectory\",\"scope\":\"Each original surface power row retains its own watts and duty, including multiple rows on one surface. Delivered watts enter through the original uniform inward flux and exact physical P1 trace area. The complete thermal multiplier already includes storage, admitted feedback and the timestep convention. Prescribed equations have zero load weights; an initial sampled peak has zero future heating sensitivity. Geometry, initial state and time grid stay fixed.\",\"sources\":[{}]}}",
            self.steps, rows.join(","),
        ))
    }
}
