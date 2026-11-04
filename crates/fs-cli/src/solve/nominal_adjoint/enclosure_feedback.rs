//! One native enclosure derivative through the complete shared thermal dual.
//! Fresh original convection, natural slopes and transported-air references
//! are bound by the same owner used by ambient coupled radiation.
use std::collections::BTreeMap;

use fs_conduction::adjoint::{RobinGradient, RobinResponse};
use fs_conduction::radiation::{
    EnclosureRadiationLinearization, GrayDiffuseEnclosure, OpenDiffuseEnclosure,
    RadiationSurface, ViewFactorEvidence, ViewFactorMatrix,
};
use fs_conduction::{ScalarField, ThermalBc};
use fs_project::{EnclosureViewFactorEvidence, ThermalBoundaryCondition as B};

use super::{
    ConductionProblem, Cx, ProjectSpec, RungSolved, SolveRefusal, bad,
    contact_controls, convective_feedback, lower, poll, row,
};

pub(super) const SCOPE: &str = "Enclosure mode differentiates the complete reflected gray-diffuse radiosity, smooth material/contact laws, natural convection and transported-air feedback at the unchanged final field. Uniform prescribed-flux patches, including vacuum, retain their original external heat and need no convection control. Only original convection heats the air; internal enclosure exchange retains its exact view-factor order. An open enclosure additionally has one shared radiation-reservoir-temperature control for its prescribed black environment. This reservoir has no FEM degree of freedom or independent emissivity; its heat stays external. Natural-reservoir and air-inlet controls contract the same complete dual as source, convection, fan, material, contact and emissivity controls. Enclosure emissivity rows are local absolute coefficient partials at the selected immutable card query. This is a local discrete derivative at fixed geometry and view factors, not a continuum or geometric error certificate.";

pub(super) struct Pullback {
    pub gradient: RobinGradient,
    pub natural_ambient: BTreeMap<String, f64>,
    pub air_rows: Vec<String>,
    pub radiation_rows: Vec<String>,
    pub contact_rows: Option<Vec<String>>,
}

pub(super) enum Model<'a> {
    Closed(&'a GrayDiffuseEnclosure),
    Open(&'a OpenDiffuseEnclosure, f64),
}

impl Model<'_> {
    fn surfaces(&self) -> &[RadiationSurface] {
        match self {
            Self::Closed(model) => model.surfaces(),
            Self::Open(model, _) => model.surfaces(),
        }
    }

    fn view_factors(&self) -> &ViewFactorMatrix {
        match self {
            Self::Closed(model) => model.view_factors(),
            Self::Open(model, _) => model.view_factors(),
        }
    }

    fn ambient_temperature(&self) -> Option<f64> {
        match self { Self::Closed(_) => None, Self::Open(_, bath) => Some(*bath) }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn bind(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, linear: fs_conduction::LinearConfig,
        temperature: &[f64], names: &[&str], slopes: &[f64], feedback: &[f64],
        heat_tolerance_w: f64, max_entries: usize,
    ) -> Result<EnclosureRadiationLinearization, SolveRefusal> {
        match self {
            Self::Closed(model) => EnclosureRadiationLinearization::bind(
                cx, problem, linear, temperature, names, slopes, feedback, model,
                heat_tolerance_w, max_entries),
            Self::Open(model, bath) => EnclosureRadiationLinearization::bind_open(
                cx, problem, linear, temperature, names, slopes, feedback, model, *bath,
                heat_tolerance_w, max_entries),
        }.map_err(lower)
    }
}

pub(super) fn admit(spec: &ProjectSpec) -> Result<(), SolveRefusal> {
    let Some(setup) = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref()) else {
        return Ok(());
    };
    if setup.enclosure_radiation.is_none() { return Ok(()); }
    if setup.radiation.is_some() || setup.transient.is_some() {
        return Err(bad("steady enclosure nominal derivatives require one steady enclosure model"));
    }
    Ok(())
}

/// Bind derivative consumers to the exact authored model and physical traces.
/// State-dependent h and air references are re-evaluated by convective_feedback;
/// they must never be read from the radiation-shifted numerical boundary.
pub(super) fn model<'a>(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &'a RungSolved,
) -> Result<Model<'a>, SolveRefusal> {
    poll(cx)?;
    admit(spec)?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing native enclosure conduction setup"))?;
    let declaration = setup.enclosure_radiation.as_ref()
        .ok_or_else(|| bad("missing native enclosure declaration"))?;
    let data = solved.adjoint_data.as_ref()
        .ok_or_else(|| bad("missing retained enclosure adjoint operator"))?;
    let enclosure = match (&data.enclosure, &data.open_enclosure) {
        (Some(model), None) => Model::Closed(model),
        (None, Some((model, bath))) => Model::Open(model, *bath),
        _ => return Err(bad("exactly one retained native enclosure model is required")),
    };
    let p = enclosure.surfaces().len();
    let factors = enclosure.view_factors().factors();
    if !data.radiation_patches.is_empty()
        || declaration.ambient_temperature.map(|bath| bath.value) != enclosure.ambient_temperature()
        || declaration.surfaces.len() != p || declaration.view_factors.len() != p
        || factors.len() != p + usize::from(enclosure.ambient_temperature().is_some())
        || factors.iter().take(p).zip(&declaration.view_factors).any(|(actual, declared)|
            actual.get(..p) != Some(declared.as_slice()))
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
    if !same_evidence {
        return Err(bad("enclosure derivative lost its declared view-factor evidence"));
    }
    for (declared, surface) in declaration.surfaces.iter().zip(enclosure.surfaces()) {
        poll(cx)?;
        if declared.target != surface.name()
            || declared.card != surface.emissivity().card_identity().to_hex()
            || declared.query_temperature.value != surface.emissivity().temperature_k()
        {
            return Err(bad("enclosure derivative lost the exact emissivity source or target"));
        }
        let original = setup.boundaries.iter().find(|r| r.target == declared.target)
            .ok_or_else(|| bad("enclosure target lost its declared boundary law"))?;
        let region = data.boundary.region_names().iter().position(|name| name == &declared.target)
            .ok_or_else(|| bad("enclosure target has no original native boundary"))?;
        match (&original.condition, &data.boundary.conditions()[region]) {
            (B::Convection { coefficient, reference_temperature },
                ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference) }) => {
                if *h != coefficient.value || *reference != reference_temperature.value {
                    return Err(bad("enclosure derivative cannot substitute a shifted radiative reference for original convection"));
                }
            }
            (B::NaturalConvection { ambient_temperature, .. },
                ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference) }) => {
                if !h.is_finite() || *h <= 0.0 || *reference != ambient_temperature.value {
                    return Err(bad("enclosure natural convection lost its original law"));
                }
            }
            (B::AirflowConvection { .. },
                ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference) }) => {
                if !h.is_finite() || *h <= 0.0 || !reference.is_finite() || *reference <= 0.0
                    || !data.air_paths.iter().any(|path|
                        path.segments().iter().any(|segment| segment.region() == declared.target))
                {
                    return Err(bad("enclosure derivative lost its original transported-air path"));
                }
            }
            (B::HeatFlux { outward_flux },
                ThermalBc::Neumann { outward_flux: ScalarField::Uniform(original_flux) }) => {
                if !original_flux.is_finite() || *original_flux != outward_flux.value {
                    return Err(bad("enclosure derivative lost its original prescribed heat flux"));
                }
            }
            _ => return Err(bad("enclosure patch has no matching original uniform convection or prescribed-flux law")),
        }
    }
    poll(cx)?;
    Ok(enclosure)
}

pub(super) fn pullback(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, weights: &[f64],
    contact_requested: bool, max_pairs: usize,
) -> Result<Pullback, SolveRefusal> {
    let enclosure = model(cx, spec, solved)?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("missing native enclosure setup"))?;
    let declaration = setup.enclosure_radiation.as_ref()
        .ok_or_else(|| bad("missing native enclosure declaration"))?;
    let data = solved.adjoint_data.as_ref()
        .ok_or_else(|| bad("missing retained enclosure derivative operator"))?;
    let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
    let entries = usize::try_from(memory / 32).unwrap_or(usize::MAX);
    let bound = convective_feedback::bind(cx, spec, solved, entries)?;
    let names: Vec<_> = bound.names.iter().map(String::as_str).collect();
    let problem = ConductionProblem {
        mesh: &solved.mesh, boundary: &bound.boundary, material: &data.fallback,
        element_materials: Some(&data.materials), source: &data.source,
    };
    let radiation = enclosure.bind(
        cx, problem, data.linear, &solved.solution.temperature, &names,
        &bound.natural.wall_slopes, &bound.feedback, declaration.heat_tolerance.value, entries,
    )?;
    let physical = ConductionProblem { boundary: radiation.boundary(), ..problem };
    let (convection, contact_rows) = if let Some(contacts) = &data.temperature_contact {
        if data.interfaces.is_some() {
            return Err(bad("live enclosure contact cannot also retain a frozen interface operator"));
        }
        // Requery both independent contact temperatures on the freshly bound
        // radiative endpoint. Its complete tangent, including bulk K'(T), is
        // the base for the existing radiosity/air low-rank feedback owner.
        let response = contacts.linearize(cx, physical, &solved.solution.temperature,
            data.linear, &names, max_pairs)?;
        let full = response.pullback_with_robin_and_flux_feedback_and_resistance_controls(
            cx, radiation.convective_mean_slopes(), radiation.local_reference_slopes(),
            radiation.reference_feedback(), radiation.flux_feedback(), weights, entries,
            max_pairs,
        ).map_err(lower)?;
        let rows = if contact_requested {
            contact_controls::rows_from_gradients(cx, spec, full.resistance_multipliers,
                "R''(T_A,T_B) -> s R''(T_A,T_B), evaluated at s=1")?
        } else { Vec::new() };
        (full.transport, Some(rows))
    } else {
        let gradient = RobinResponse::pullback_mean_robin_with_reference_and_flux_feedback(
            cx, physical, data.interfaces.as_ref(), data.linear,
            &solved.solution.temperature, &names,
            radiation.convective_mean_slopes(), radiation.local_reference_slopes(),
            radiation.reference_feedback(), radiation.flux_feedback(), weights, entries,
        ).map_err(|error| {
            let mut refusal = lower(error);
            refusal.fix = "tighten physical solid/air, enclosure temperature/heat and thermal solve tolerances; fresh radiosity, the unchanged physical residual and complete transpose must all converge".into();
            refusal
        })?;
        (gradient, None)
    };
    let full = radiation.remap_gradient(cx, convection, data.linear).map_err(lower)?;
    let natural_ambient = bound.natural_ambient(&full.convection)?;
    let air_rows = bound.air_rows(cx, setup, solved, &full.convection)?;
    let mut rows = Vec::with_capacity(declaration.surfaces.len());
    for (i, surface) in declaration.surfaces.iter().enumerate() {
        poll(cx)?;
        rows.push(row("radiation-emissivity", &surface.name, i, "1", full.emissivities[i])?);
    }
    if let Some(derivative) = full.ambient_temperature {
        rows.push(row("radiation-reservoir-temperature", "enclosure-environment", 0, "K", derivative)?);
    }
    poll(cx)?;
    Ok(Pullback { gradient: full.convection, natural_ambient, air_rows,
        radiation_rows: rows, contact_rows })
}
