//! A mesh comparison must solve the SAME declared physics on both meshes.
//! Reuse the radiating candidate producer, including its total coupled adjoint,
//! rather than falling through the convection-only Request::evaluate path.
use super::*;

pub(super) struct Solved {
    pub output: String,
    pub objective_k: f64,
    pub source_w: f64,
    pub solid_solves: usize,
    pub adjoint_sweeps: usize,
    pub marking: Option<mark::Marking>,
}

pub(super) fn solve(request: &Request, cx: &Cx<'_>, fraction: Option<f64>) -> Result<Solved> {
    poll(cx)?;
    let flow = request.flow(cx)?;
    let coefficients = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
    let mut evaluated = match &request.radiation {
        Some(policy) => policy.evaluate(request, cx, &flow, &coefficients, fraction.is_some())?,
        None => {
            let value = request.evaluate(cx, &flow, &coefficients, fraction.is_some())?;
            let solid_solves = value.coupled.iterations;
            fan_speed::ThermalEvaluation { value, radiation: None, solid_solves }
        }
    };
    let marking = fraction.map(|fraction| {
        let gradient = evaluated.value.gradient.as_ref()
            .ok_or_else(|| producer("adaptive solve has no total coupled adjoint"))?;
        mark::evaluate(cx, request, &evaluated.value.temperatures, &gradient.nodal_load, fraction)
    }).transpose()?;
    let adjoint_sweeps = evaluated.value.gradient.as_ref().map_or(0, |g| g.iterations);
    // The explicit mesh strategy owns this derivative. Keep ordinary output
    // gradient fields null; neither the temperatures nor mechanism report is
    // recomputed with different physics to render the primal result.
    evaluated.value.gradient = None;
    let output = evaluated.render(request, &flow)?;
    let output = match &request.fan {
        Some(fan) => fan.attach(output, &flow, fan.speed_ratio)?,
        None => output,
    };
    poll(cx)?;
    Ok(Solved {
        output,
        objective_k: evaluated.value.objective,
        source_w: evaluated.value.source_total_w,
        solid_solves: evaluated.solid_solves,
        adjoint_sweeps,
        marking,
    })
}
