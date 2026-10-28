//! Bind one unchanged-field convection law for all native coupled adjoints.
//! The production air march supplies physical reference values; its analytic
//! wall map and the existing natural-law owner supply their complete slopes.
use std::collections::BTreeMap;

use fs_airflow::conjugate::goal::{CoupledGoalError, maximum::affine_reference_law};
use fs_conduction::adjoint::RobinGradient;
use fs_conduction::{ScalarField, ThermalBc, ThermalBoundary};
use fs_project::ConductionSetup;

use super::{
    Cx, ProjectSpec, RungSolved, SolveRefusal, air_inlet, bad, conduction_error,
    finite, lower, mixed_feedback, poll, zeros,
};

pub(super) struct Bound {
    pub boundary: ThermalBoundary,
    pub names: Vec<String>,
    pub means: Vec<f64>,
    pub natural: mixed_feedback::Bound,
    pub feedback: Vec<f64>,
    air_names: Vec<String>,
    air_indices: Vec<usize>,
}

impl Bound {
    pub fn natural_ambient(
        &self, gradient: &RobinGradient,
    ) -> Result<BTreeMap<String, f64>, SolveRefusal> {
        let names: Vec<_> = self.names.iter().map(String::as_str).collect();
        self.natural.ambient(&names, gradient)
    }

    pub fn air_rows(
        &self, cx: &Cx<'_>, setup: &ConductionSetup, solved: &RungSolved,
        gradient: &RobinGradient,
    ) -> Result<Vec<String>, SolveRefusal> {
        if gradient.references.len() != self.names.len() {
            return Err(bad("coupled convection pullback changed its original trace order"));
        }
        let data = solved.adjoint_data.as_ref()
            .ok_or_else(|| bad("missing retained air derivative operator"))?;
        if data.air_paths.is_empty() {
            if !super::super::conjugate::airflow_laws(setup)?.is_empty()
                || !self.air_names.is_empty() || !self.air_indices.is_empty()
            {
                return Err(bad("coupled derivative omitted a declared air path"));
            }
            return Ok(Vec::new());
        }
        let names: Vec<_> = self.air_names.iter().map(String::as_str).collect();
        let bars: Vec<_> = self.air_indices.iter().map(|&i| gradient.references[i]).collect();
        air_inlet::rows(cx, setup, &data.air_paths, &names, &bars)
    }
}

pub(super) fn bind(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, max_entries: usize,
) -> Result<Bound, SolveRefusal> {
    poll(cx)?;
    let data = solved.adjoint_data.as_ref()
        .ok_or_else(|| bad("missing coupled convection operator"))?;
    let temperature = &solved.solution.temperature;
    let n = solved.mesh.vertex_count();
    if temperature.len() != n {
        return Err(bad("coupled convection field has the wrong vertex count"));
    }
    let mut names = Vec::new();
    let mut regions = Vec::new();
    for (region, bc) in data.boundary.conditions().iter().enumerate() {
        poll(cx)?;
        match bc {
            ThermalBc::Robin {
                htc: ScalarField::Uniform(_), t_ref: ScalarField::Uniform(_),
            } => {
                names.push(data.boundary.region_names()[region].clone());
                regions.push(region);
            }
            ThermalBc::Robin { .. } =>
                return Err(bad("native coupled derivatives require uniform Robin traces")),
            _ => {},
        }
    }
    let m = names.len();
    if m > 64 || n > max_entries || m*m > max_entries
        || n.checked_mul(4).and_then(|v| v.checked_mul(m))
            .is_none_or(|v| v > max_entries)
    {
        return Err(bad("coupled convection exceeds its declared factor budget or 64 traces"));
    }
    if data.radiating_boundary.as_ref().is_some_and(|b|
        b.region_names() != data.boundary.region_names())
    {
        return Err(bad("radiation changed the original convective region ownership"));
    }
    let mut areas = zeros(m)?;
    let mut means = zeros(m)?;
    for (slot, face) in solved.mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let region = data.boundary.region_for(slot);
        if data.radiating_boundary.as_ref().is_some_and(|b| b.region_for(slot) != region) {
            return Err(bad("radiation changed the original convective face ownership"));
        }
        let Some(i) = region.and_then(|r| regions.iter().position(|&v| v == r))
            else { continue; };
        areas[i] = finite(areas[i] + face.area)?;
        for &v in &face.vertices {
            means[i] = finite(means[i] + face.area/3.0 * finite(temperature[v as usize])?)?;
        }
    }
    for i in 0..m {
        if areas[i] <= 0.0 { return Err(bad("coupled convection has an empty trace")); }
        means[i] = finite(means[i]/areas[i])?;
    }
    let mut replacements = Vec::new();
    let mut feedback = zeros(m*m)?;
    let mut air_names = Vec::new();
    let mut air_indices = Vec::new();
    if !data.air_paths.is_empty() {
        let law = affine_reference_law(cx, &data.air_paths, 64, max_entries)
            .map_err(|e| match e {
                CoupledGoalError::Interrupted => conduction_error("cli-solve-cancelled",
                    "coupled air-reference lowering interrupted",
                    "resume the accepted pipeline prefix"),
                e => bad(format!("complete air-reference law refused: {e}")),
            })?;
        air_names = law.regions().to_vec();
        air_indices = air_names.iter().map(|name| names.iter().position(|n| n == name)
            .ok_or_else(|| bad("an air segment has no original convective trace")))
            .collect::<Result<_, _>>()?;
        let mut start = 0;
        for path in &data.air_paths {
            poll(cx)?;
            let end = start + path.segments().len();
            let selected = air_indices.get(start..end)
                .ok_or_else(|| bad("coupled air law lost stream-wise segment ordering"))?;
            let walls: Vec<_> = selected.iter().map(|&i| means[i]).collect();
            // Match the physical producer's arithmetic, rather than D*T+d.
            let march = path.march(&walls)
                .map_err(|e| bad(format!("retained-field air march refused: {e}")))?;
            if march.segments.len() != selected.len() {
                return Err(bad("retained air march omitted a physical segment"));
            }
            for ((segment, state), &i) in path.segments().iter()
                .zip(&march.segments).zip(selected)
            {
                let ThermalBc::Robin { htc: ScalarField::Uniform(h), .. } =
                    &data.boundary.conditions()[regions[i]]
                    else { return Err(bad("air trace lost its original Robin coefficient")); };
                if segment.region() != names[i] || *h != segment.htc_w_per_m2_k()
                    || (areas[i] - segment.area_m2()).abs()
                        > 128.0*f64::EPSILON*areas[i].max(segment.area_m2())
                {
                    return Err(bad("air and solid must retain the same region, coefficient and wetted area"));
                }
                replacements.push((regions[i], *h, state.reference_temperature_k));
            }
            start = end;
        }
        if start != air_indices.len() {
            return Err(bad("coupled air law omitted retained segments"));
        }
        for (i, &a) in air_indices.iter().enumerate() {
            poll(cx)?;
            for (j, &b) in air_indices.iter().enumerate() {
                feedback[a*m+b] = law.wall_matrix()[i*air_indices.len()+j];
            }
        }
    }
    let refs: Vec<_> = names.iter().map(String::as_str).collect();
    let natural = mixed_feedback::bind(cx, spec, solved, &refs, &regions, &means)?;
    replacements.extend_from_slice(&natural.replacements);
    let boundary = data.boundary.with_uniform_robin_replacements(&replacements).map_err(lower)?;
    poll(cx)?;
    Ok(Bound { boundary, names, means, natural, feedback, air_names, air_indices })
}
