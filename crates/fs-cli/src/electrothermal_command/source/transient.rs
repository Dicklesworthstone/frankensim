//! Time evolution through the same material/electrical/heat owners as the
//! steady command. Cell mass is fixed; storage uses the sourced Cp chart.

use std::path::Path;

use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectrodeDrive, TransientRecord};
use fs_matdb::{PropertyKey, SelectionPolicy};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::Density;

use super::{Problem, ResultField, invalid, poll};
use super::super::{J, array, budget, fields, integer, number, positive, quote, read, scalar, string, values};

const NO_CLAIM: &str = "fixed-mass fixed-geometry single-solid-phase backward-Euler electrothermal evolution; nominal source data with bounded Cp-chart interpolation and optional exterior gray radiation to black isothermal unit-view reservoirs, no propagated uncertainty, phase change, thermal deformation, fluid transport or experimental qualification; temporal accuracy requires time-step refinement";

pub(crate) struct Request {
    dt: f64,
    steps: usize,
    lower_k: f64,
    upper_k: f64,
    interpolation_j_kg: f64,
    absolute_j: f64,
    relative: f64,
    // First step (1-based), scale of all electrode drives, scale of authored power.
    schedule: Vec<(usize, f64, f64)>,
}

pub(crate) fn load(path: &Path, cells: usize) -> Result<Request, String> {
    let bytes = read(path, 256 * 1024)?;
    let root = J::parse(std::str::from_utf8(&bytes).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    fields(&root, &["schema", "time_step_s", "steps", "temperature_min_k", "temperature_max_k",
        "enthalpy_interpolation_error_j_kg", "energy_absolute_tolerance_j", "energy_relative_tolerance", "schedule"])?;
    if string(&root, "schema")? != "frankensim.electrothermal.transient.v1" {
        return Err("unsupported transient request schema".into());
    }
    let dt = positive(&root, "time_step_s")?;
    let steps = budget(&root, "steps", 10_000)?;
    if cells.checked_mul(steps).is_none_or(|n| n > 16_777_216) || !(dt * steps as f64).is_finite() {
        return Err("transient cell-step budget or final-time range exceeded".into());
    }
    let lower_k = positive(&root, "temperature_min_k")?;
    let upper_k = positive(&root, "temperature_max_k")?;
    if lower_k >= upper_k { return Err("storage temperature range must strictly increase".into()); }
    let interpolation_j_kg = positive(&root, "enthalpy_interpolation_error_j_kg")?;
    let absolute_j = positive(&root, "energy_absolute_tolerance_j")?;
    let relative = scalar(&root, "energy_relative_tolerance")?;
    if !(0.0..1.0).contains(&relative) { return Err("energy relative tolerance must be in [0,1)".into()); }
    let mut schedule = Vec::new();
    for row in array(&root, "schedule", 10_000)? {
        let row = row.as_array().filter(|r| r.len() == 3)
            .ok_or("schedule rows must be [first_step, drive_scale, authored_power_scale]")?;
        let first = integer(&row[0])?;
        if first == 0 || first > steps || schedule.last().is_some_and(|&(previous, _, _)| first <= previous) {
            return Err("schedule steps must be strictly increasing inside 1..=steps".into());
        }
        schedule.push((first, number(&row[1])?, number(&row[2])?));
    }
    if schedule.first().map(|r| r.0) != Some(1) {
        return Err("the first schedule row must explicitly start at step 1".into());
    }
    Ok(Request { dt, steps, lower_k, upper_k, interpolation_j_kg, absolute_j, relative, schedule })
}

/// One tentative time-step storage image, never a second material law. The
/// accepted old temperature and mass stay immutable across every nonlinear pass.
pub(super) struct Storage<'a> {
    curves: &'a [ResolvedHeatCapacityCurve],
    indices: &'a [usize],
    masses: &'a [f64],
    pub previous: &'a [f64],
    pub dt: f64,
    absolute_j: f64,
    relative: f64,
}

impl Storage<'_> {
    pub(super) fn capacities(&self, temperature: &[f64], gate: &CancelGate) -> Result<Vec<f64>, ChtError> {
        if temperature.len() != self.previous.len() { return Err(invalid("material.storage", "temperature shape mismatch")); }
        let mut capacities = Vec::with_capacity(temperature.len());
        for (cell, &t) in temperature.iter().enumerate() {
            if cell.is_multiple_of(256) { poll(gate)?; }
            let cp = self.curves[self.indices[cell]].secant_heat_capacity(self.previous[cell], t)
                .map_err(|e| invalid("material.storage", e))?;
            let capacity = self.masses[cell] * cp;
            if !capacity.is_finite() || capacity <= 0.0 {
                return Err(invalid("material.storage", "cell heat capacity is not finite positive J/K"));
            }
            capacities.push(capacity);
        }
        Ok(capacities)
    }

    fn changes(&self, temperature: &[f64], gate: &CancelGate) -> Result<Vec<f64>, ChtError> {
        if temperature.len() != self.previous.len() { return Err(invalid("material.storage", "temperature shape mismatch")); }
        let mut changes = Vec::with_capacity(temperature.len());
        for (cell, &t) in temperature.iter().enumerate() {
            if cell.is_multiple_of(256) { poll(gate)?; }
            let dh = self.curves[self.indices[cell]].specific_enthalpy_change(self.previous[cell], t)
                .map_err(|e| invalid("material.storage", e))?;
            let change = self.masses[cell] * dh;
            if !change.is_finite() || (dh != 0.0 && change == 0.0) {
                return Err(invalid("material.storage", "cell enthalpy increment is unrepresentable"));
            }
            changes.push(change);
        }
        Ok(changes)
    }

    pub(super) fn consistent(&self, temperature: &[f64], capacities: &[f64], gate: &CancelGate) -> Result<bool, ChtError> {
        if capacities.len() != self.previous.len() { return Err(invalid("material.storage", "capacity shape mismatch")); }
        let changes = self.changes(temperature, gate)?;
        let (mut defect, mut scale) = (0.0, 0.0);
        for (cell, &actual) in changes.iter().enumerate() {
            let linearized = capacities[cell] * (temperature[cell] - self.previous[cell]);
            defect += (actual - linearized).abs();
            scale += actual.abs() + linearized.abs();
        }
        let limit = self.absolute_j + self.relative * scale;
        if !defect.is_finite() || !limit.is_finite() { return Err(invalid("material.storage", "storage defect is not finite")); }
        // Sum absolute local defects; opposite errors cannot cancel globally.
        Ok(defect <= limit)
    }
}

pub(crate) struct Trajectory {
    last: ResultField,
    records: Vec<TransientRecord>,
    masses: Vec<f64>,
    mass_states: Vec<String>,
    curves: Vec<ResolvedHeatCapacityCurve>,
    thermal_solves: usize,
}

pub(crate) fn march(problem: &Problem, request: &Request, gate: &CancelGate) -> Result<Trajectory, ChtError> {
    poll(gate)?;
    if let Some(fluid) = &problem.fluid {
        fluid.admit_transient(&problem.materials[fluid.material_index], request.lower_k,
            request.upper_k, &problem.thermal, gate)?;
    }
    let density_key = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    let cp_key = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let density_requirement = ScalarPropertyRequirement::try_with_key(&density_key, ScalarAdmissibility::StrictlyPositive)
        .map_err(|e| invalid("material.storage", e))?;
    let mut curves = Vec::with_capacity(problem.materials.len());
    for (index, source) in problem.materials.iter().enumerate() {
        poll(gate)?;
        let lower = source.point_at(request.lower_k)?;
        let upper = source.point_at(request.upper_k)?;
        // The scene owns the fixed role. This is an explicit single-phase
        // incompressible-fluid assumption, not inferred phase stability.
        let phase = if problem.fluid.as_ref().is_some_and(|f| f.material_index == index) {
            SolidLiquidPhase::Liquid
        } else {
            SolidLiquidPhase::Solid
        };
        curves.push(ResolvedHeatCapacityCurve::try_new(
            source.pack.card(), HeatCapacitySourceQuery {
                heat_capacity: &cp_key, density: &density_key,
                temperature_axis: &source.temperature_axis, lower: &lower, upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            }, phase, 0.0, request.interpolation_j_kg,
        ).map_err(|e| invalid("material.storage", e))?);
    }
    let mut masses = Vec::with_capacity(problem.domain.cell_count());
    let mut mass_states = Vec::with_capacity(problem.domain.cell_count());
    for (cell, &temperature) in problem.initial_temperature.iter().enumerate() {
        if cell.is_multiple_of(256) { poll(gate)?; }
        let m = problem.material_indices[cell];
        let source = &problem.materials[m];
        curves[m].state_at_temperature(temperature).map_err(|e| invalid("material.storage", e))?;
        let state = resolve_material_state_point(
            source.pack.card(), &source.point_at(temperature)?, core::slice::from_ref(&density_requirement),
            MaterialPropertySelection::SingleClaimOnly,
        ).map_err(|e| invalid("material.initial_mass", e))?;
        let mass = state.property(DENSITY_PROPERTY).expect("density was required").value_si() * problem.domain.volume(cell);
        if !mass.is_finite() || mass <= 0.0 { return Err(invalid("material.initial_mass", "cell mass must be finite positive kg")); }
        masses.push(mass);
        mass_states.push(state.identity().to_hex());
    }
    let mut temperature = problem.initial_temperature.clone();
    let mut records = Vec::with_capacity(request.steps);
    let mut last = None;
    let mut thermal_solves = 0usize;
    for step in 1..=request.steps {
        poll(gate)?;
        let index = request.schedule.partition_point(|row| row.0 <= step) - 1;
        let (_, drive_scale, power_scale) = request.schedule[index];
        let mut electric = problem.electric.clone();
        for electrode in &mut electric.electrodes {
            let value = match &mut electrode.drive {
                ElectrodeDrive::Voltage(v) | ElectrodeDrive::Current(v) => v,
            };
            *value *= drive_scale;
            if !value.is_finite() { return Err(invalid("transient.schedule", "electrode drive overflow")); }
        }
        let mut thermal = problem.thermal.clone();
        for power in &mut thermal.power_w {
            *power *= power_scale;
            if !power.is_finite() { return Err(invalid("transient.schedule", "authored heat source overflow")); }
        }
        let storage = Storage {
            curves: &curves, indices: &problem.material_indices, masses: &masses,
            previous: &temperature, dt: request.dt, absolute_j: request.absolute_j, relative: request.relative,
        };
        let result = problem.solve_at(&temperature, &electric, &thermal, Some(&storage), gate)?;
        let next = &result.coupled.thermal;
        let changes = storage.changes(&next.temperature, gate)?;
        let stored: f64 = changes.iter().sum();
        let stored_scale: f64 = changes.iter().map(|e| e.abs()).sum();
        let physical_sink: f64 = thermal.cell_sinks.iter().map(|sink| {
            sink.conductance_w_k * (next.temperature[sink.cell] - sink.temperature_k)
        }).sum();
        // Use electrical heat at the accepted endpoint, not the preceding
        // outer iterate's source, when auditing the complete time step.
        let source_j = request.dt * (thermal.power_w.iter().sum::<f64>() + result.coupled.electric.total_joule_w);
        // Boundary radiation/convection was already present in the heat solve.
        // Audit its independently re-evaluated physical flux, not storage sinks.
        let surface_w: f64 = result.radiation.iter().map(|r| r.outward_w()).sum();
        let radiated_w: f64 = result.radiation.iter().map(|r| r.radiated_w).sum();
        let outflow_j = request.dt * (next.report.balance.boundary_outflow_w + physical_sink + surface_w);
        let closure_j = stored + outflow_j - source_j;
        let tolerance = request.absolute_j + request.relative * (stored_scale + source_j.abs() + outflow_j.abs());
        if [stored, source_j, outflow_j, closure_j, tolerance, result.coupled.electric.delivered_w]
            .iter().any(|v| !v.is_finite()) || closure_j.abs() > tolerance
        {
            return Err(invalid("transient.energy_closure", format!("step {step}: defect {closure_j} J exceeds {tolerance} J")));
        }
        let peak = next.temperature.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let record = TransientRecord {
            time_s: request.dt * step as f64, max_temperature_k: peak,
            max_solid_temperature_k: next.max_where(|cell| !problem.domain.is_fluid(cell))
                .map_or(f64::NAN, |(_, t)| t),
            stored_energy_change_j: stored, source_j, boundary_outflow_j: outflow_j,
            closure_j, radiated_w, iterations: next.report.iterations,
        };
        poll(gate)?;
        temperature.clone_from(&next.temperature);
        thermal_solves += result.thermal_solves;
        records.push(record);
        last = Some(result);
    }
    poll(gate)?;
    Ok(Trajectory { last: last.ok_or_else(|| invalid("transient.steps", "no accepted step"))?,
        records, masses, mass_states, curves, thermal_solves })
}

pub(crate) fn render(problem: &Problem, trajectory: &Trajectory, json: bool) -> String {
    let last = &trajectory.last.coupled;
    let record = trajectory.records.last().expect("a completed trajectory has steps");
    let no_claim = if problem.fluid.is_some() { super::fluid::TRANSIENT_NO_CLAIM } else { NO_CLAIM };
    if !json {
        return format!("Electrothermal transient: {} steps, {} s, {} thermal solves\nPeak temperature: {} K; final Joule power: {} W\nLast stored-energy increment: {} J; energy defect: {} J\nThermal-contact relative defect: {}\n{}\n",
            trajectory.records.len(), record.time_s, trajectory.thermal_solves,
            record.max_temperature_k, last.electric.total_joule_w,
            record.stored_energy_change_j, record.closure_j,
            trajectory.last.thermal_contact_relative_defect, no_claim);
    }
    let records = trajectory.records.iter().map(|r| format!(
        "{{\"time_s\":{},\"max_temperature_k\":{},\"stored_energy_change_j\":{},\"source_j\":{},\"physical_outflow_j\":{},\"closure_j\":{},\"radiated_w\":{}}}",
        r.time_s, r.max_temperature_k, r.stored_energy_change_j, r.source_j, r.boundary_outflow_j, r.closure_j, r.radiated_w,
    )).collect::<Vec<_>>().join(",");
    let materials = problem.materials.iter().zip(&trajectory.curves).map(|(source, chart)| format!(
        "{{\"pack_identity\":{},\"card_identity\":{},\"storage_chart_identity\":{},\"storage_receipts\":{}}}",
        quote(&source.pack.content_hash().to_hex()), quote(&source.pack.card().content_hash().to_hex()),
        quote(&chart.curve().identity().to_hex()), chart.receipts().len(),
    )).collect::<Vec<_>>().join(",");
    let source_states = trajectory.last.sources.iter().enumerate().map(|(cell, state)| format!(
        "{{\"material_index\":{},\"initial_mass_state\":{},\"thermal_state\":{},\"electrical_state\":{}}}",
        problem.material_indices[cell], quote(&trajectory.mass_states[cell]), quote(&state.thermal.identity().to_hex()),
        state.conductor.as_ref().map_or_else(|| "null".into(), |c| quote(&c.material().identity().to_hex())),
    )).collect::<Vec<_>>().join(",");
    format!(concat!(
        "{{\"schema\":\"frankensim.electrothermal.transient.result.v1\",\"authority\":\"estimated-numerical\",",
        "\"time_s\":{},\"steps\":{},\"cells\":{:?},\"cell_mass_kg\":{},\"temperature_k\":{},\"potential_v\":{},\"joule_w\":{},",
        "\"total_joule_w\":{},\"electrode_power_w\":{},\"resolved_conductivity_w_m_k\":{},",
        "\"thermal_solves\":{},\"records\":[{}],\"materials\":[{}],\"source_states\":[{}],",
        "\"electrical_interfaces\":{},\"thermal_interfaces\":{},",
        "\"thermal_contact_relative_defect\":{},\"radiation\":{},\"fluid_reference\":{},\"no_claim\":{}}}\n"),
        record.time_s, trajectory.records.len(), problem.domain.dims(), values(trajectory.masses.iter().copied()),
        values(last.thermal.temperature.iter().copied()), values(last.electric.potential_v.iter().copied()),
        values(last.electric.joule_w.iter().copied()), last.electric.total_joule_w, last.electric.delivered_w,
        values(trajectory.last.sources.iter().map(|s| s.k)), trajectory.thermal_solves, records, materials,
        source_states, super::super::interfaces::render(problem, &trajectory.last.interfaces),
        super::super::interfaces::render_thermal(problem, &trajectory.last),
        trajectory.last.thermal_contact_relative_defect,
        super::super::radiation::render(&trajectory.last.radiation),
        problem.fluid.as_ref().map_or_else(|| "null".into(), super::Fluid::render), quote(no_claim))
}
