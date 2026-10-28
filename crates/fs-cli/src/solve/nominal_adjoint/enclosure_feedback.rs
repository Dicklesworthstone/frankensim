//! Native prescribed-convection enclosure goals through the shared full transpose.
//! Every source and original boundary is retained from the accepted producer.
use fs_conduction::adjoint::RobinGradient;
use fs_conduction::radiation::{ViewFactorEvidence, pullback_gray_diffuse_enclosure};
use fs_conduction::{ScalarField, ThermalBc};
use fs_project::{EnclosureViewFactorEvidence, ThermalBoundaryCondition as B};

use super::{ConductionProblem, Cx, ProjectSpec, RungSolved, SolveRefusal, bad, lower, poll, row};

pub(super) const SCOPE: &str = "Closed-enclosure mode differentiates the complete reflected gray-diffuse radiosity and smooth material/contact Jacobian at the unchanged final field. Every radiating patch retains its original prescribed positive convection and exact view-factor order. The full nonlocal patch tangent and its true transpose use the shared numerical owners; radiative loads are not frozen. Enclosure emissivity rows are local absolute coefficient partials at the selected immutable card query. Enclosure derivatives with natural convection, transported air, prescribed-temperature controls, transient storage, adaptive goals and continuum error bounds are not supplied.";

pub(super) struct Pullback {
    pub gradient: RobinGradient,
    pub radiation_rows: Vec<String>,
}

/// Refuse an incomplete state or prescribed-value derivative before solving.
pub(super) fn admit(spec: &ProjectSpec) -> Result<(), SolveRefusal> {
    let Some(setup) = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref()) else {
        return Ok(());
    };
    if setup.enclosure_radiation.is_none() { return Ok(()); }
    if setup.radiation.is_some() || setup.transient.is_some() {
        return Err(bad("enclosure nominal derivatives require one steady enclosure model"));
    }
    if setup.boundaries.iter().any(|r| matches!(r.condition,
        B::NaturalConvection { .. } | B::AirflowConvection { .. }))
    {
        return Err(bad("native enclosure adjoints require prescribed convection; natural or transported-air feedback cannot be frozen"));
    }
    let prescribed_requested = spec.outputs.iter().flatten().any(|r| matches!(r.name.as_str(),
        "temperature-max-boundary-adjoint" | "temperature-max-contact-boundary-adjoint"
        | "temperature-volume-mean-adjoint"));
    if prescribed_requested && setup.boundaries.iter().any(|r|
        matches!(r.condition, B::FixedTemperature { .. }))
    {
        return Err(bad("native enclosure prescribed-temperature controls require the complete enclosure lift; use a maximum/contact adjoint with fixed prescribed temperatures"));
    }
    Ok(())
}

pub(super) fn pullback(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, weights: &[f64],
) -> Result<Pullback, SolveRefusal> {
    poll(cx)?;
    admit(spec)?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing native enclosure conduction setup"))?;
    let declaration = setup.enclosure_radiation.as_ref()
        .ok_or_else(|| bad("missing native enclosure declaration"))?;
    let data = solved.adjoint_data.as_ref()
        .ok_or_else(|| bad("missing retained enclosure adjoint operator"))?;
    let enclosure = data.enclosure.as_ref()
        .ok_or_else(|| bad("missing retained native enclosure model"))?;
    if !data.air_paths.is_empty() || !data.radiation_patches.is_empty()
        || declaration.surfaces.len() != enclosure.surfaces().len()
        || declaration.view_factors.as_slice() != enclosure.view_factors().factors()
    {
        return Err(bad("retained enclosure law differs from its declared physical model"));
    }
    let same_evidence = match (&declaration.evidence, enclosure.view_factors().evidence()) {
        (EnclosureViewFactorEvidence::Analytic { geometry: a },
            ViewFactorEvidence::Analytic { geometry: b }) => a == b,
        (EnclosureViewFactorEvidence::ExternalQmc { seed: a, samples: b, generator: c },
            ViewFactorEvidence::ExternalQmc { seed: x, samples: y, generator: z }) =>
            a == x && b == y && c == z,
        _ => false,
    };
    if !same_evidence { return Err(bad("enclosure derivative lost its declared view-factor evidence")); }
    let mut names = Vec::with_capacity(declaration.surfaces.len());
    for (declared, surface) in declaration.surfaces.iter().zip(enclosure.surfaces()) {
        poll(cx)?;
        if declared.target != surface.name()
            || declared.card != surface.emissivity().card_identity().to_hex()
            || declared.query_temperature.value != surface.emissivity().temperature_k()
        {
            return Err(bad("enclosure derivative lost the exact emissivity source or target"));
        }
        let original = setup.boundaries.iter().find(|r| r.target == declared.target)
            .ok_or_else(|| bad("enclosure target lost its declared convection"))?;
        let B::Convection { coefficient, reference_temperature } = &original.condition else {
            return Err(bad("native enclosure derivative needs original prescribed convection"));
        };
        let region = data.boundary.region_names().iter().position(|name| name == &declared.target)
            .ok_or_else(|| bad("enclosure target has no original native boundary"))?;
        let ThermalBc::Robin {
            htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference),
        } = &data.boundary.conditions()[region] else {
            return Err(bad("enclosure derivative needs original uniform Robin traces"));
        };
        if *h != coefficient.value || *reference != reference_temperature.value {
            return Err(bad("enclosure derivative cannot substitute a shifted radiative reference for original convection"));
        }
        names.push(declared.target.as_str());
    }
    let problem = ConductionProblem {
        mesh: &solved.mesh, boundary: &data.boundary, material: &data.fallback,
        element_materials: Some(&data.materials), source: &data.source,
    };
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let entries = usize::try_from(memory / 32).unwrap_or(usize::MAX);
    let full = pullback_gray_diffuse_enclosure(
        cx, problem, data.interfaces.as_ref(), data.linear, &solved.solution.temperature,
        &names, enclosure, declaration.heat_tolerance.value, weights, entries,
    ).map_err(|error| {
        let mut refusal = lower(error);
        refusal.fix = "tighten the declared enclosure temperature/heat tolerances and thermal solve tolerance; fresh radiosity, the unchanged physical residual and complete transpose must all converge".into();
        refusal
    })?;
    let mut rows = Vec::with_capacity(declaration.surfaces.len());
    for (i, surface) in declaration.surfaces.iter().enumerate() {
        poll(cx)?;
        rows.push(row("radiation-emissivity", &surface.name, i, "1", full.emissivities[i])?);
    }
    poll(cx)?;
    Ok(Pullback { gradient: full.convection, radiation_rows: rows })
}
