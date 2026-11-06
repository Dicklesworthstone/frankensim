//! Native enriched goals through the primal's complete mean-boundary laws.
//! Each field owns its natural law, physical air march and optional radiation.
//! The core owner supplies the complete transpose, material/contact terms and
//! observed nonlinear remainder. Nothing here certifies continuum error.

use fs_airflow::conjugate::goal::{CoupledGoalError, maximum::affine_reference_law};
use fs_conduction::adjoint::RobinResponse;
use fs_conduction::radiation::{AmbientRadiationGradient, EnclosureRadiationGradient};
use fs_conduction::{ConductionError, ConductionProblem, ScalarField, ThermalBc};
use fs_exec::Cx;

use super::super::{AdaptiveProbe, RungSolved, SolveRefusal, conduction_error, ladder_functional};

const MAX_REGIONS: usize = 64;

fn bad(what: impl Into<String>) -> SolveRefusal {
    conduction_error(
        "cli-solve-conduction-boundary-adaptive",
        what,
        "check the retained natural/radiative laws, original convection/air traces and feedback budget; tighten physical boundary and solid/air convergence tolerances rather than freezing a coefficient",
    )
}

fn interrupted() -> SolveRefusal {
    conduction_error("cli-solve-cancelled", "coupled boundary adaptive comparison cancelled", "resume the retained solve prefix")
}

fn poll(cx: &Cx<'_>, index: usize) -> Result<(), SolveRefusal> {
    if index % 256 == 0 { cx.checkpoint().map_err(|_| interrupted())?; }
    Ok(())
}

fn finite(value: f64) -> Result<f64, SolveRefusal> {
    if value.is_finite() { Ok(value) } else { Err(bad("nonfinite coupled boundary adaptive arithmetic")) }
}

fn enclosure_energy_allowance(fine: &RungSolved, heat_w: f64) -> Result<f64, SolveRefusal> {
    // The fixed-contact primal gates the physical L2 residual, without a
    // separate global-energy stop. Its summed free residual is bounded by
    // sqrt(n_free)*residual_w; add the declared radiosity heat allowance.
    // Live contacts instead retain their original independent energy gate.
    let root = fs_math::next_up((fine.solution.report.free_dofs as f64).sqrt());
    let residual = finite(fs_math::next_up(fine.solution.report.residual_threshold * root))?;
    finite(fs_math::next_up(residual + heat_w))
}

fn zeros<T: Default + Clone>(n: usize) -> Result<Vec<T>, SolveRefusal> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| bad("coupled boundary adaptive scratch allocation refused"))?;
    values.resize(n, T::default());
    Ok(values)
}

fn solid_error(error: ConductionError) -> SolveRefusal {
    match error {
        ConductionError::Cancelled { .. } => interrupted(),
        error => bad(format!("complete boundary-law goal refused: {error}")),
    }
}

fn air_error(error: CoupledGoalError) -> SolveRefusal {
    match error {
        CoupledGoalError::Interrupted
        | CoupledGoalError::Solid(ConductionError::Cancelled { .. })
        | CoupledGoalError::Air(fs_airflow::AirflowError::Cancelled { .. }) => interrupted(),
        error => bad(format!("adaptive air law refused: {error}")),
    }
}

/// Shared binding for radiative or mixed natural/air goals. An incomplete
/// radiation binding refuses instead of falling through to a nonradiating
/// comparison. Natural coefficients and air references are evaluated at BOTH
/// fields; their local and cross-wall derivatives enter one total transpose.
pub(in crate::solve) fn probe(
    cx: &Cx<'_>, coarse: &RungSolved, fine: &RungSolved, approximate: &[f64],
    region: u32, memory_bytes: u64,
) -> Result<AdaptiveProbe, SolveRefusal> {
    poll(cx, 0)?;
    let data = fine.adjoint_data.as_ref().ok_or_else(|| bad("missing coupled boundary goal operator"))?;
    let combined = data.radiating_boundary.as_ref();
    let ambient = !data.radiation_patches.is_empty();
    let enclosed = data.enclosure.is_some() || data.open_enclosure.is_some();
    if combined.is_some() != (ambient || enclosed)
        || (ambient && enclosed)
        || (data.enclosure.is_some() && data.open_enclosure.is_some())
        || enclosed != data.enclosure_heat_tolerance_w.is_some()
    {
        return Err(bad("retain exactly one original radiation model, its accepted boundary, and its declared enclosure heat allowance when applicable"));
    }
    let live_contact = data.temperature_contact.is_some();
    if live_contact && data.interfaces.is_some() {
        return Err(bad("adaptive live contact requires its original source; a separate frozen contact cannot substitute for its coupled tangent"));
    }
    if combined.is_none()
        && !(data.natural_goal.is_some() && !data.air_paths.is_empty())
        && !(live_contact && (data.natural_goal.is_some() || !data.air_paths.is_empty()))
    {
        return Err(bad("a nonradiating coupled boundary comparison needs natural/air feedback, or live contact with either physical convection law"));
    }
    let n = fine.mesh.vertex_count();
    if approximate.len() != n || fine.solution.temperature.len() != n
        || fine.labels.len() != fine.mesh.element_count()
        || coarse.labels.len() != coarse.mesh.element_count()
        || coarse.solution.temperature.len() != coarse.mesh.vertex_count()
        || coarse.mesh.element_count().checked_mul(8) != Some(fine.mesh.element_count())
    {
        return Err(bad("coupled adaptive fields or uniform enrichment lineage differ"));
    }
    for (i, &value) in approximate.iter().chain(&fine.solution.temperature)
        .chain(&coarse.solution.temperature).enumerate()
    { poll(cx, i)?; finite(value)?; }
    if combined.is_some_and(|boundary| boundary.region_names() != data.boundary.region_names()) {
        return Err(bad("radiation changed the original convection region ownership"));
    }
    let mut names = Vec::new();
    let mut regions = Vec::new();
    let mut points = Vec::new();
    let mut other = Vec::new();
    for (i, bc) in data.boundary.conditions().iter().enumerate() {
        poll(cx, i)?;
        match bc {
            ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(r) } => {
                if names.len() == MAX_REGIONS { return Err(bad("coupled adaptive goal exceeds 64 Robin traces")); }
                names.push(data.boundary.region_names()[i].as_str());
                regions.push(i);
                points.push([*h, *r, 0.0, 0.0]);
                other.push([*h, *r]);
            }
            ThermalBc::Robin { .. } => return Err(bad("coupled adaptive goals require uniform original Robin traces")),
            ThermalBc::Neumann { outward_flux: ScalarField::Uniform(flux) }
                if *flux == 0.0 => {
                let name = data.boundary.region_names()[i].as_str();
                if let Some(patch) = data.radiation_patches.iter().find(|patch| patch.region() == name) {
                    if names.len() == MAX_REGIONS { return Err(bad("coupled adaptive goal exceeds 64 thermal traces")); }
                    // This is still the original zero-convection trace. Its
                    // reservoir placeholder is inert; the radiation owner
                    // supplies the actual secant and weighted reference.
                    names.push(name);
                    regions.push(i);
                    points.push([0.0, patch.ambient_temperature_k(), 0.0, 0.0]);
                    other.push([0.0, patch.ambient_temperature_k()]);
                }
            }
            _ => {},
        }
    }
    let m = names.len();
    let entries = usize::try_from(memory_bytes / 32).unwrap_or(usize::MAX);
    // Enclosure patch means and flux factors belong to the radiation owner,
    // independently of convection. A genuine vacuum may have no Robin rows.
    if (m == 0 && !enclosed) || m*m > entries
        || n.checked_mul(4).and_then(|v| v.checked_mul(m)).is_none_or(|v| v > entries)
    { return Err(bad("coupled adaptive goal exceeds the declared feedback-entry budget")); }
    let mut areas = zeros::<f64>(m)?;
    let mut means = zeros::<[f64; 2]>(m)?;
    for (slot, face) in fine.mesh.boundary().iter().enumerate() {
        poll(cx, slot)?;
        let region = data.boundary.region_for(slot);
        if combined.is_some_and(|boundary| boundary.region_for(slot) != region) {
            return Err(bad("radiation changed the original convection face ownership"));
        }
        let Some(i) = region.and_then(|r| regions.iter().position(|&v| v == r)) else { continue; };
        areas[i] = finite(areas[i] + face.area)?;
        // Match energy_balance's original wall observation, including its
        // arithmetic order. A few temperature ulps matter for the contact
        // watt gate on a high-conductance wall.
        let reference_mean = face.vertices.iter()
            .map(|&v| fine.solution.temperature[v as usize]).sum::<f64>() / 3.0;
        let approximate_mean = face.vertices.iter()
            .map(|&v| approximate[v as usize]).sum::<f64>() / 3.0;
        means[i][0] = finite(face.area.mul_add(reference_mean, means[i][0]))?;
        means[i][1] = finite(face.area.mul_add(approximate_mean, means[i][1]))?;
    }
    for i in 0..m {
        if areas[i] <= 0.0 { return Err(bad("coupled adaptive thermal trace has no positive area")); }
        let inverse_area = 1.0 / areas[i];
        means[i][0] = finite(means[i][0] * inverse_area)?;
        means[i][1] = finite(means[i][1] * inverse_area)?;
    }
    if let Some(natural) = &data.natural_goal {
        // h_c, its wall derivative, and the comparison field's h_c all come
        // from the same retained card/pressure/multiplier as the primal.
        // The radiation owner then differentiates BOTH h_total and the
        // weighted reference; adding dh_c to a frozen Robin row is not enough.
        natural.rebind_convection(cx, &names, &means, &mut points, &mut other)?;
    }
    let mut feedback = zeros::<f64>(m*m)?;
    if !data.air_paths.is_empty() {
        let law = affine_reference_law(cx, &data.air_paths, MAX_REGIONS, entries).map_err(air_error)?;
        if data.natural_goal.as_ref().is_some_and(|natural|
            law.regions().iter().any(|name| natural.contains_target(name)))
        {
            return Err(bad("natural reservoirs and forced-air segments must own distinct Robin traces"));
        }
        let indices: Vec<usize> = law.regions().iter().map(|name| names.iter()
            .position(|&selected| selected == name.as_str())
            .ok_or_else(|| bad("an air segment has no original convection trace")))
            .collect::<Result<_, _>>()?;
        let mut start = 0usize;
        for (path_index, path) in data.air_paths.iter().enumerate() {
            poll(cx, path_index)?;
            let end = start.checked_add(path.segments().len())
                .filter(|&end| end <= indices.len())
                .ok_or_else(|| bad("air path and affine law segment ownership differ"))?;
            // March both actual physical states, not the rounded affine map;
            // only the original convective coefficient transfers heat to air.
            for state in 0..2 {
                let walls: Vec<_> = indices[start..end].iter().map(|&i| means[i][state]).collect();
                let march = path.march(&walls).map_err(|error| bad(format!("adaptive air march refused: {error}")))?;
                if march.segments.len() != end-start { return Err(bad("air march lost a selected segment")); }
                for ((segment, result), &i) in path.segments().iter().zip(&march.segments).zip(&indices[start..end]) {
                    if points[i][0] != segment.htc_w_per_m2_k()
                        || (areas[i]-segment.area_m2()).abs() > 128.0*f64::EPSILON*areas[i].max(segment.area_m2())
                    { return Err(bad("air and solid do not share the original convection coefficient and wetted area")); }
                    let reference = finite(result.reference_temperature_k)?;
                    if state == 0 { points[i][1] = reference; } else { other[i][1] = reference; }
                }
            }
            start = end;
        }
        if start != indices.len() || law.wall_matrix().len() != indices.len()*indices.len() {
            return Err(bad("air-reference feedback has inconsistent segment dimensions"));
        }
        for (i, &row) in indices.iter().enumerate() {
            poll(cx, i)?;
            for (j, &column) in indices.iter().enumerate() {
                feedback[row*m+column] = law.wall_matrix()[i*indices.len()+j];
            }
        }
    }
    let mut peak = None;
    for (index, (tet, &label)) in fine.mesh.complex().tets.iter().zip(&fine.labels).enumerate() {
        poll(cx, index)?;
        if label != region { continue; }
        for &v in tet {
            let v = v as usize;
            if peak.is_none_or(|old: usize| fine.solution.temperature[v] > fine.solution.temperature[old]
                || (fine.solution.temperature[v] == fine.solution.temperature[old] && v < old))
            { peak = Some(v); }
        }
    }
    let peak = peak.ok_or_else(|| bad("coupled adaptive goal region has no tetrahedra"))?;
    let mut weights = zeros::<f64>(n)?;
    weights[peak] = 1.0;
    let problem = ConductionProblem { mesh: &fine.mesh, boundary: &data.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let dual = if let Some(enclosure) = &data.enclosure {
        let heat_w = data.enclosure_heat_tolerance_w
            .ok_or_else(|| bad("closed enclosure lost its declared heat allowance"))?;
        if let Some(contact) = &data.temperature_contact {
            contact.compare_goal_with_enclosure(cx, problem,
                &fine.solution.temperature, approximate, data.linear, &names,
                &points, &other, &feedback, enclosure, heat_w, &weights, entries)?
        } else {
            EnclosureRadiationGradient::compare_goal_at(cx, problem, data.interfaces.as_ref(),
                data.linear, &fine.solution.temperature, approximate,
                fine.solution.report.residual_threshold, enclosure_energy_allowance(fine, heat_w)?, &names, &points,
                &other, &feedback, enclosure, heat_w, &weights, entries).map_err(solid_error)?
        }
    } else if let Some((enclosure, bath)) = &data.open_enclosure {
        let heat_w = data.enclosure_heat_tolerance_w
            .ok_or_else(|| bad("open enclosure lost its declared heat allowance"))?;
        if let Some(contact) = &data.temperature_contact {
            contact.compare_goal_with_open_enclosure(cx, problem,
                &fine.solution.temperature, approximate, data.linear, &names,
                &points, &other, &feedback, enclosure, *bath, heat_w, &weights, entries)?
        } else {
            EnclosureRadiationGradient::compare_open_goal_at(cx, problem, data.interfaces.as_ref(),
                data.linear, &fine.solution.temperature, approximate,
                fine.solution.report.residual_threshold, enclosure_energy_allowance(fine, heat_w)?, &names, &points,
                &other, &feedback, enclosure, *bath, heat_w, &weights, entries).map_err(solid_error)?
        }
    } else if let Some(contact) = &data.temperature_contact {
        if combined.is_some() {
            contact.compare_goal_with_ambient_radiation(cx, problem,
                &fine.solution.temperature, approximate, data.linear, &names,
                &points, &other, &feedback, &data.radiation_patches, &weights, entries)?
        } else {
            contact.compare_goal_with_robin_feedback(cx, problem,
                &fine.solution.temperature, approximate, data.linear, &names,
                &points, &other, &feedback, &weights, entries)?
        }
    } else if combined.is_some() {
        AmbientRadiationGradient::compare_goal_at(cx, problem, data.interfaces.as_ref(),
            data.linear, &fine.solution.temperature, approximate, &names, &points, &other,
            &feedback, &data.radiation_patches, &weights, entries).map_err(solid_error)?
    } else {
        RobinResponse::compare_mean_robin_goal_with_reference_feedback_at(
            cx, problem, data.interfaces.as_ref(), data.linear, &fine.solution.temperature,
            approximate, &names, &points, &other, &feedback, &weights, entries).map_err(solid_error)?
    };
    let coarse_max = ladder_functional(&coarse.labels, &coarse.mesh.complex().tets,
        &coarse.solution.temperature, Some(region));
    let remainder = finite(approximate[peak] - coarse_max)?;
    let measured = finite(fine.solution.temperature[peak] - coarse_max)?;
    let estimated = finite(dual.signed_residual_change.abs()
        + dual.linearization_remainder.abs() + remainder.abs())?;
    if dual.nodal_contributions.len() != n { return Err(bad("coupled goal lost nodal contribution ownership")); }
    let mut incidence = zeros::<usize>(n)?;
    for (index, tet) in fine.mesh.complex().tets.iter().enumerate() {
        poll(cx, index)?;
        for &v in tet { incidence[v as usize] += 1; }
    }
    let mut scores = zeros::<f64>(coarse.mesh.element_count())?;
    for (index, tet) in fine.mesh.complex().tets.iter().enumerate() {
        poll(cx, index)?;
        for &v in tet {
            let v = v as usize;
            scores[index/8] = finite(scores[index/8] + dual.nodal_contributions[v].abs()/incidence[v] as f64)?;
        }
    }
    poll(cx, 0)?;
    Ok(AdaptiveProbe {
        signed_linear_change_k: dual.signed_residual_change,
        linearization_remainder_k: dual.linearization_remainder,
        maximum_remainder_k: remainder,
        estimated_change_k: estimated,
        measured_change_k: measured,
        dual_residual: dual.dual_relative_residual,
        dual_iterations: dual.dual_iterations,
        primal_residual: dual.primal_relative_residual,
        uses_nonlinear_jacobian: dual.uses_nonlinear_jacobian,
        // The legacy field specifically describes IQN-ILS. This comparison
        // uses the complete core FGMRES transpose, not that interface solver.
        coupled_adjoint: None,
        scores,
    })
}
