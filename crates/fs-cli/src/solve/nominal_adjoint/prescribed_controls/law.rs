//! Re-evaluate existing physical laws for the complete prescribed-value lift.
//! Natural and transported-air feedback share the nominal convection binding;
//! enclosure reflection shares the core's exact unchanged-field linearization.
use super::{
    BTreeSet, ConductionProblem, Cx, ProjectSpec, RungSolved, ScalarField,
    SolveRefusal, ThermalBc, ThermalBoundary, bad, finite, lower, poll, zeros,
};
use super::super::{convective_feedback, enclosure_feedback};
use fs_conduction::radiation::EnclosureRadiationLinearization;
use fs_conduction::adjoint::MeanFluxFeedback;

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
        let radiation = EnclosureRadiationLinearization::bind(
            cx, problem, data.linear, &solved.solution.temperature, &names,
            &bound.natural.wall_slopes, &bound.feedback, enclosure,
            declaration.heat_tolerance.value, max_entries,
        ).map_err(lower)?;
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
    let m = bound.names.len();
    let mut regions = Vec::with_capacity(m);
    let mut points = Vec::with_capacity(m);
    for (i, name) in bound.names.iter().enumerate() {
        poll(cx)?;
        let region = bound.boundary.region_names().iter().position(|r| r == name)
            .ok_or_else(|| bad("prescribed lift lost an original convection trace"))?;
        let ThermalBc::Robin {
            htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference),
        } = &bound.boundary.conditions()[region] else {
            return Err(bad("prescribed lift requires its original uniform Robin law"));
        };
        regions.push(region);
        points.push([*h, *reference, bound.natural.wall_slopes[i], 0.0]);
    }
    let mut feedback = bound.feedback;
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
        let i = bound.names.iter().position(|name| name == patch.region())
            .ok_or_else(|| bad("prescribed radiative lift has no original convection trace"))?;
        let [hc, tc, hcs, _] = points[i];
        let hr = patch.secant_coefficient_w_m2_k(bound.means[i]).map_err(lower)?;
        let [hrs, _] = patch.secant_partials_w_m2_k2(bound.means[i]).map_err(lower)?;
        let tr = patch.ambient_temperature_k();
        let h = finite(hc + hr)?;
        if h <= 0.0 { return Err(bad("prescribed lift has a nonpositive combined coefficient")); }
        let reference = finite((hc/h)*tc + (hr/h)*tr)?;
        points[i] = [h, reference, finite(hcs + hrs)?,
            finite((hcs/h)*(tc - reference) + (hrs/h)*(tr - reference))?];
        for entry in &mut feedback[i*m..(i+1)*m] {
            *entry = finite((hc/h)*(*entry))?;
        }
    }
    let replacements: Vec<_> = regions.iter().zip(&points)
        .map(|(&r, p)| (r, p[0], p[1])).collect();
    let boundary = bound.boundary.with_uniform_robin_replacements(&replacements).map_err(lower)?;
    let mut h_slopes = zeros(m)?;
    let mut reference_slopes = zeros(m)?;
    for (i, point) in points.iter().enumerate() {
        h_slopes[i] = point[2];
        reference_slopes[i] = point[3];
    }
    poll(cx)?;
    Ok(Law { boundary, names: bound.names, h_slopes, reference_slopes, feedback, flux_feedback: None })
}
