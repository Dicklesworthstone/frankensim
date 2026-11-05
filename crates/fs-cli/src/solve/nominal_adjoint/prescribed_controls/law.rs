//! Re-evaluate existing physical laws for the complete prescribed-value lift.
//! Natural and transported-air feedback share the nominal convection binding;
//! enclosure reflection shares the core's exact unchanged-field linearization.
use super::{
    BTreeSet, ConductionProblem, Cx, ProjectSpec, RungSolved,
    SolveRefusal, ThermalBoundary, bad, lower, poll,
};
use super::super::{convective_feedback, enclosure_feedback};
use fs_conduction::adjoint::MeanFluxFeedback;
use fs_conduction::radiation::AmbientRadiationLinearization;

pub(super) struct Law {
    pub boundary: ThermalBoundary,
    pub names: Vec<String>,
    pub h_slopes: Vec<f64>,
    pub reference_slopes: Vec<f64>,
    pub feedback: Vec<f64>,
    pub flux_feedback: Option<MeanFluxFeedback>,
}

pub(super) fn bind(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, max_entries: usize,
) -> Result<Law, SolveRefusal> {
    poll(cx)?;
    let data = solved.adjoint_data.as_ref()
        .ok_or_else(|| bad("prescribed lift lost its operator"))?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("prescribed lift has no conduction declarations"))?;
    let bound = convective_feedback::bind(cx, spec, solved, max_entries)?;
    if let Some(declaration) = &setup.enclosure_radiation {
        let enclosure = enclosure_feedback::model(cx, spec, solved)?;
        let names: Vec<_> = bound.names.iter().map(String::as_str).collect();
        let problem = ConductionProblem {
            mesh: &solved.mesh, boundary: &bound.boundary, material: &data.fallback,
            element_materials: Some(&data.materials), source: &data.source,
        };
        let radiation = enclosure.bind(
            cx, problem, data.linear, &solved.solution.temperature, &names,
            &bound.natural.wall_slopes, &bound.feedback,
            declaration.heat_tolerance.value, max_entries,
        )?;
        // RobinResponse uses these same complete coefficients for J_FD as
        // the selected-goal solve uses for J_FF. Fixed vertices remain in means.
        return Ok(Law {
            boundary: radiation.boundary().clone(),
            names: radiation.regions().to_vec(),
            h_slopes: radiation.convective_mean_slopes().to_vec(),
            reference_slopes: radiation.local_reference_slopes().to_vec(),
            feedback: radiation.reference_feedback().to_vec(),
            flux_feedback: radiation.flux_feedback().cloned(),
        });
    }
    let surfaces = setup.radiation.as_ref().map_or(&[][..], |r| r.surfaces.as_slice());
    if surfaces.len() != data.radiation_patches.len() || surfaces.len() > 64 {
        return Err(bad("prescribed lift lost original radiation patches"));
    }
    let mut seen = BTreeSet::new();
    for (surface, patch) in surfaces.iter().zip(&data.radiation_patches) {
        poll(cx)?;
        if surface.target != patch.region() || !seen.insert(patch.region())
            || surface.reservoir_temperature.value != patch.ambient_temperature_k()
            || surface.card != patch.emissivity().card_identity().to_hex()
            || surface.query_temperature.value != patch.emissivity().temperature_k()
        {
            return Err(bad("prescribed radiative lift changed its original card or physical law"));
        }
    }
    // Use exactly the same coefficient/reference owner as the selected-goal
    // dual, including genuine radiation-only traces with zero convection.
    let names: Vec<_> = bound.names.iter().map(String::as_str).collect();
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &bound.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let radiation = AmbientRadiationLinearization::bind(cx, problem,
        &solved.solution.temperature, &names, &bound.natural.wall_slopes,
        Some(&bound.feedback), &data.radiation_patches, max_entries).map_err(lower)?;
    poll(cx)?;
    Ok(Law {
        boundary: radiation.boundary().clone(),
        names: radiation.regions().to_vec(),
        h_slopes: radiation.coefficient_mean_slopes().to_vec(),
        reference_slopes: radiation.local_reference_slopes().to_vec(),
        feedback: radiation.reference_feedback()
            .ok_or_else(|| bad("prescribed lift lost original reference feedback"))?.to_vec(),
        flux_feedback: None,
    })
}
