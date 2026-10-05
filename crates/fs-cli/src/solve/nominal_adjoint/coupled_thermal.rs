//! One full material/radiation/air adjoint on the unchanged native field.
//! The actual air march supplies references; its analytic wall map supplies
//! feedback. Radiation remains a separate reservoir, never a heat load to air.

use fs_airflow::conjugate::goal::{CoupledGoalError, maximum::affine_reference_law};
use fs_conduction::{ScalarField, ThermalBc};
use fs_conduction::adjoint::RobinGradient;
use fs_conduction::radiation::AmbientRadiationGradient;
use super::{ConductionProblem, Cx, ProjectSpec, RungSolved, SolveRefusal,
    air_inlet, bad, conduction_error, finite, lower, poll, row, zeros};

pub(super) struct Pullback {
    pub gradient: RobinGradient,
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
    if data.air_paths.is_empty() || setup.boundaries.iter().any(|b| matches!(
        b.condition, fs_project::ThermalBoundaryCondition::NaturalConvection { .. })) {
        return Err(bad("coupled thermal adjoints require fixed-flow air paths and no natural-convection state law"));
    }
    let n = solved.mesh.vertex_count();
    let temperature = &solved.solution.temperature;
    if temperature.len() != n || weights.len() != n { return Err(bad("coupled thermal field/goal size differs")); }
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let entries = usize::try_from(memory/32).unwrap_or(usize::MAX);
    let mut names = Vec::new();
    let mut regions = Vec::new();
    for (region, bc) in data.boundary.conditions().iter().enumerate() {
        poll(cx)?;
        match bc {
            ThermalBc::Robin { htc: ScalarField::Uniform(_), t_ref: ScalarField::Uniform(_) } => {
                names.push(data.boundary.region_names()[region].as_str());
                regions.push(region);
            }
            ThermalBc::Robin { .. } => return Err(bad("coupled thermal adjoints require uniform Robin traces")),
            _ => {},
        }
    }
    let m = names.len();
    if m == 0 || m > 64 || n.checked_mul(4).and_then(|v| v.checked_mul(m)).is_none_or(|v| v>entries)
        || m*m > entries {
        return Err(bad("coupled material/radiation/air adjoint exceeds its declared factor budget or 64 traces"));
    }
    if let Some(combined) = &data.radiating_boundary {
        if combined.region_names() != data.boundary.region_names() {
            return Err(bad("radiation changed the original convective region ownership"));
        }
    }
    let mut areas = zeros(m)?;
    let mut means = zeros(m)?;
    for (slot, face) in solved.mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let region = data.boundary.region_for(slot);
        if data.radiating_boundary.as_ref().is_some_and(|b| b.region_for(slot)!=region) {
            return Err(bad("radiation changed the original convective face ownership"));
        }
        let Some(i) = region.and_then(|r| regions.iter().position(|&v| v==r)) else { continue; };
        areas[i] = finite(areas[i]+face.area)?;
        for &v in &face.vertices {
            means[i] = finite(means[i]+face.area/3.0*finite(temperature[v as usize])?)?;
        }
    }
    for i in 0..m {
        if areas[i] <= 0.0 { return Err(bad("coupled thermal trace has zero area")); }
        means[i] = finite(means[i]/areas[i])?;
    }
    let law = affine_reference_law(cx, &data.air_paths, 64, entries).map_err(|e| match e {
        CoupledGoalError::Interrupted => conduction_error("cli-solve-cancelled",
            "coupled thermal air lowering interrupted", "resume the accepted pipeline prefix"),
        e => bad(format!("complete air-reference law refused: {e}")),
    })?;
    let air_names: Vec<&str> = law.regions().iter().map(String::as_str).collect();
    let indices: Vec<usize> = air_names.iter().map(|name| names.iter().position(|n| n==name)
        .ok_or_else(|| bad("an air segment has no original convective trace"))).collect::<Result<_,_>>()?;
    let mut replacements = Vec::new();
    let mut start = 0;
    for path in &data.air_paths {
        poll(cx)?;
        let end = start+path.segments().len();
        let walls: Vec<_> = indices[start..end].iter().map(|&i| means[i]).collect();
        // March the physical model, not D*T+d rounded in a different order.
        // No solid solve or retained-field replacement is performed here.
        let march = path.march(&walls).map_err(|e| bad(format!("retained-field air march refused: {e}")))?;
        for ((segment, state), &i) in path.segments().iter().zip(&march.segments).zip(&indices[start..end]) {
            let ThermalBc::Robin { htc: ScalarField::Uniform(h), .. } = &data.boundary.conditions()[regions[i]]
                else { return Err(bad("air trace is not an original uniform Robin law")); };
            if *h != segment.htc_w_per_m2_k()
                || (areas[i]-segment.area_m2()).abs() > 128.0*f64::EPSILON*areas[i].max(segment.area_m2()) {
                return Err(bad("air and solid must share the same original coefficient and wetted area"));
            }
            replacements.push((regions[i],*h,state.reference_temperature_k));
        }
        start = end;
    }
    let mut feedback = zeros(m*m)?;
    for (i, &solid_i) in indices.iter().enumerate() {
        poll(cx)?;
        for (j, &solid_j) in indices.iter().enumerate() {
            feedback[solid_i*m+solid_j] = law.wall_matrix()[i*indices.len()+j];
        }
    }
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
    let boundary = data.boundary.with_uniform_robin_replacements(&replacements).map_err(lower)?;
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let gradient = AmbientRadiationGradient::pullback_with_reference_feedback(cx, problem,
        data.interfaces.as_ref(), data.linear, temperature, &names, &vec![0.0;m], &feedback,
        &data.radiation_patches, weights, entries).map_err(|e| {
            let mut refusal = lower(e);
            refusal.fix = "tighten the physical solid/air and radiation convergence controls; the unchanged field and complete coupled transpose must satisfy their residual gates".into();
            refusal
        })?;
    let bars: Vec<_> = indices.iter().map(|&i| gradient.convection.references[i]).collect();
    let air_rows = air_inlet::rows(cx, setup, &data.air_paths, &air_names, &bars)?;
    let mut radiation_rows = Vec::new();
    for (i, surface) in surfaces.iter().enumerate() {
        radiation_rows.push(row("radiation-reservoir-temperature", &surface.name, i, "K",
            gradient.reservoir_temperatures[i])?);
        radiation_rows.push(row("radiation-emissivity", &surface.name, i, "1", gradient.emissivities[i])?);
    }
    poll(cx)?;
    Ok(Pullback { gradient: gradient.convection, air_rows, radiation_rows })
}
