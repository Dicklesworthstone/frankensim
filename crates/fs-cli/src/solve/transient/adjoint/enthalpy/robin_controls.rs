//! Physical prescribed convection uses the same masked thermal load multiplier.
//! The consistent triangular mass is retained on nonisothermal traces.

use fs_conduction::{ConductionProblem, ScalarField, ThermalBc};
use fs_exec::Cx;
use fs_project::{ProjectSpec, ThermalBoundaryCondition};
use std::time::Instant;

use super::{SolveRefusal, bad, finite, json_string, number, poll, zeros};

/// Independent original-reference and ln(h) partials, with dt already present
/// in load_weights. Shared by prescribed, natural and transported-air controls.
pub(super) fn partials(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    temperature: &[f64],
    names: &[String],
    load_weights: &[f64],
    deadline: Option<(Instant, f64)>,
) -> Result<(Vec<f64>, Vec<f64>), SolveRefusal> {
    let n = problem.mesh.vertex_count();
    if temperature.len() != n || load_weights.len() != n || names.len() > 64 {
        return Err(bad("enthalpy convection controls lost their complete nodal/port layout"));
    }
    let mut regions = Vec::with_capacity(names.len());
    let mut htc = Vec::with_capacity(names.len());
    let mut references_k = Vec::with_capacity(names.len());
    for name in names {
        poll(cx, deadline)?;
        let region = problem.boundary.region_names().iter().position(|found| found == name)
            .ok_or_else(|| bad("enthalpy convection control has no original trace"))?;
        if regions.contains(&region) {
            return Err(bad("enthalpy convection control repeats a physical trace"));
        }
        let ThermalBc::Robin {
            htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference),
        } = problem.boundary.conditions()[region] else {
            return Err(bad("enthalpy convection controls require original uniform Robin laws"));
        };
        if !(h.is_finite() && h > 0.0 && reference.is_finite() && reference > 0.0) {
            return Err(bad("enthalpy convection control has invalid physical coefficients"));
        }
        regions.push(region);
        htc.push(h);
        references_k.push(reference);
    }
    let mut references = zeros(names.len())?;
    let mut log_htc = zeros(names.len())?;
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % 256 == 0 { poll(cx, deadline)?; }
        let Some(index) = problem.boundary.region_for(slot)
            .and_then(|region| regions.iter().position(|&found| found == region))
        else { continue; };
        for (a, &row) in face.vertices.iter().enumerate() {
            let weight = finite(load_weights[row as usize])?;
            references[index] = finite(
                references[index] + weight * htc[index] * face.area / 3.0,
            )?;
            for (b, &column) in face.vertices.iter().enumerate() {
                let mass = htc[index] * (face.area / 12.0) * if a == b { 2.0 } else { 1.0 };
                log_htc[index] = finite(log_htc[index] + weight * mass
                    * (references_k[index] - finite(temperature[column as usize])?))?;
            }
        }
    }
    poll(cx, deadline)?;
    Ok((references, log_htc))
}

struct Control {
    name: String,
    coefficient: f64,
    reference: f64,
    reference_bar: f64,
    log_coefficient_bar: f64,
}

pub(super) struct Controls {
    rows: Vec<Control>,
    steps: usize,
}

impl Controls {
    pub(super) fn new(
        cx: &Cx<'_>, spec: &ProjectSpec, problem: ConductionProblem<'_>,
        deadline: Option<(Instant, f64)>,
    ) -> Result<Self, SolveRefusal> {
        let setup = spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
            .ok_or_else(|| bad("enthalpy convection controls lost their declaration"))?;
        let mut rows = Vec::new();
        for boundary in &setup.boundaries {
            poll(cx, deadline)?;
            let ThermalBoundaryCondition::Convection { coefficient, reference_temperature }
                = &boundary.condition else { continue; };
            let region = problem.boundary.region_names().iter()
                .position(|name| name == &boundary.target)
                .ok_or_else(|| bad("prescribed convection control has no physical trace"))?;
            let ThermalBc::Robin {
                htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference),
            } = problem.boundary.conditions()[region] else {
                return Err(bad("prescribed convection control changed its physical law"));
            };
            if h != coefficient.value || reference != reference_temperature.value
                || !(h.is_finite() && h > 0.0 && reference.is_finite() && reference > 0.0)
            {
                return Err(bad("prescribed convection control differs from its declaration"));
            }
            rows.push(Control {
                name: boundary.target.clone(), coefficient: h, reference,
                reference_bar: 0.0, log_coefficient_bar: 0.0,
            });
        }
        if rows.len() > 64 { return Err(bad("prescribed convection exceeds 64 controls")); }
        Ok(Self { rows, steps: 0 })
    }

    pub(super) fn add(
        &mut self, cx: &Cx<'_>, problem: ConductionProblem<'_>, temperature: &[f64],
        load_weights: &[f64], deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        let names: Vec<String> = self.rows.iter().map(|row| row.name.clone()).collect();
        let (references, coefficients) = partials(
            cx, problem, temperature, &names, load_weights, deadline,
        )?;
        for ((row, reference), coefficient) in self.rows.iter_mut().zip(references).zip(coefficients) {
            row.reference_bar = finite(row.reference_bar + reference)?;
            row.log_coefficient_bar = finite(row.log_coefficient_bar + coefficient)?;
        }
        self.steps = self.steps.checked_add(1)
            .ok_or_else(|| bad("prescribed convection reverse step count overflows"))?;
        poll(cx, deadline)
    }

    pub(super) fn fragment(&self, expected_steps: usize) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps {
            return Err(bad("prescribed convection controls lost the selected reverse prefix"));
        }
        if self.rows.is_empty() { return Ok(String::new()); }
        let rows = self.rows.iter().map(|row| Ok(format!(
            "{{\"target\":{},\"coefficient_w_m2_k\":{},\"reference_temperature_k\":{},\"dtemperature_dreference_temperature_k_per_k\":{},\"dtemperature_dcoefficient_k_per_w_m2_k\":{},\"dtemperature_dlog_coefficient_k\":{}}}",
            json_string(&row.name), number(row.coefficient)?, number(row.reference)?,
            number(row.reference_bar)?, number(finite(row.log_coefficient_bar / row.coefficient)?)?,
            number(row.log_coefficient_bar)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        Ok(format!(
            ",\"prescribed_convection\":{{\"control_steps\":{},\"parameter_scope\":\"independent-original-reference-and-positive-coefficient-held-throughout-selected-trajectory\",\"boundaries\":[{}]}}",
            self.steps, rows.join(","),
        ))
    }
}
