//! Native controls through the complete source-backed contact residual.
//! The numerical owner re-queries the original pinned law and revalidates the
//! unchanged accepted field before solving its true nonsymmetric transpose.

use fs_conduction::{ConductionProblem, adjoint::RobinGradient};
use fs_conduction::radiation::AmbientRadiationLinearization;
use std::collections::BTreeMap;
use super::{Cx, ProjectSpec, RungSolved, SolveRefusal, bad, contact_controls,
    convective_feedback, lower, poll, row, zeros};

pub(super) const SCOPE: &str = concat!(
    "Source-backed temperature contact retains the original selected claim and manufactured joint state at each accepted pair of independent P1 face temperatures. ",
    "The complete nonsymmetric Jacobian includes heterogeneous bulk K'(T) and both contact-temperature partials exactly once. ",
    "Prescribed temperature, flux and Robin laws, differentiated Churchill-Chu natural convection and native transported-air paths are admitted. ",
    "Natural coefficients are re-queried at the unchanged physical wall means; signed wall and ambient partials include film-temperature feedback on heated and cooled branches. ",
    "Each retained air path is physically re-marched at those same wall means; its complete ordered wall-to-reference map participates in the same dual. ",
    "Independent branch inlet and admitted fan-speed controls reuse the original transport, convection-card and hydraulic owners. ",
    "The physical endpoint is freshly checked against its original free-node watt residual and heat-balance tolerances before one total adjoint solve. ",
    "Classical derivatives require unique source/material slopes; corners and validity endpoints refuse instead of receiving a frozen or generalized sensitivity. ",
    "Steady area-mean closed/open enclosures add complete reflected-radiosity feedback, including original prescribed-flux patches. ",
    "Ambient radiation adds its exact mean-temperature secant and weighted-reference feedback, including radiation-only traces with no fictitious convection. ",
    "Original convection, emissivity and independent radiative reservoir controls are remapped from that same complete dual. ",
    "Boundary-extended reports contract that same dual with the complete prescribed-temperature lift, including direct objective weights and both contact-temperature columns. ",
    "Selected-node and regional volume-mean goals share this operator; pointwise radiation remains outside this contact adjoint. ",
    "No perturbed primal, solver-iteration derivative, inverse certificate or gradient-error bound is supplied.",
);
pub(super) const CONTROL_SCOPE: &str = "Contact-extended source-backed reports additionally differentiate each named matching interface's whole law R''(T_A,T_B) -> s R''(T_A,T_B) at s=1 through the same complete thermal adjoint. Every actual face resistance and both independent temperature traces participate; different interfaces remain independent even when they share a card. Geometry, pairing, pressure, finish, source selection and all nonthermal context coordinates stay fixed. These are dimensionless whole-law multipliers, not absolute resistance, physical joint-state, material-uncertainty or nonmatching-contact derivatives.";

pub(super) struct Pullback {
    pub(super) gradient: RobinGradient,
    pub(super) contact_rows: Vec<String>,
    pub(super) natural_ambient: BTreeMap<String, f64>,
    pub(super) air_rows: Vec<String>,
    pub(super) radiation_rows: Vec<String>,
}

pub(super) fn pullback(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved,
    weights: &[f64], contact_requested: bool, max_pairs: usize,
) -> Result<Pullback, SolveRefusal> {
    poll(cx)?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("source-backed contact adjoint requires its native conduction setup"))?;
    let data = solved.adjoint_data.as_ref()
        .ok_or_else(|| bad("source-backed contact adjoint lost its accepted operator"))?;
    let contacts = data.temperature_contact.as_ref()
        .ok_or_else(|| bad("source-backed contact adjoint lost its original sources"))?;
    if setup.transient.is_some() || setup.enclosure_radiation.is_some()
        || data.interfaces.is_some()
    {
        return Err(bad("source-backed contact adjoints require complete steady thermal boundaries; a frozen interface or boundary-feedback derivative cannot substitute"));
    }
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &data.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let has_natural = setup.boundaries.iter().any(|row| matches!(row.condition,
        fs_project::ThermalBoundaryCondition::NaturalConvection { .. }));
    let has_air = !data.air_paths.is_empty() || setup.boundaries.iter().any(|row| matches!(row.condition,
        fs_project::ThermalBoundaryCondition::AirflowConvection { .. }));
    let mut radiation_rows = Vec::new();
    let (gradient, resistance_multipliers, natural_ambient, air_rows) = if let Some(declaration) = &setup.radiation {
        if declaration.spatial_model == Some(fs_project::RadiationSpatialModel::Pointwise)
            || declaration.surfaces.is_empty() || declaration.surfaces.len() > 64
            || declaration.surfaces.len() != data.radiation_patches.len()
        {
            return Err(bad("live contact ambient adjoints require every original area-mean radiation patch"));
        }
        for (surface, patch) in declaration.surfaces.iter().zip(&data.radiation_patches) {
            poll(cx)?;
            if surface.target != patch.region()
                || surface.reservoir_temperature.value != patch.ambient_temperature_k()
                || surface.card != patch.emissivity().card_identity().to_hex()
                || surface.query_temperature.value != patch.emissivity().temperature_k()
            {
                return Err(bad("live contact radiation derivative changed its physical source or reservoir"));
            }
        }
        let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
        let entries = usize::try_from(memory / 32).unwrap_or(usize::MAX);
        let bound = convective_feedback::bind(cx, spec, solved, entries)?;
        let names: Vec<_> = bound.names.iter().map(String::as_str).collect();
        let radiation = AmbientRadiationLinearization::bind(cx,
            ConductionProblem { boundary: &bound.boundary, ..problem },
            &solved.solution.temperature, &names, &bound.natural.wall_slopes,
            Some(&bound.feedback), &data.radiation_patches, entries).map_err(lower)?;
        // Requery contacts on the unchanged, freshly bound physical boundary.
        // Its nonsymmetric K'(T) and R'(T) tangent is the base for radiation
        // and upstream-air feedback, with no second or frozen-contact solve.
        let response = contacts.linearize(cx,
            ConductionProblem { boundary: radiation.boundary(), ..problem },
            &solved.solution.temperature, data.linear, &names, max_pairs)?;
        let full = response.pullback_with_robin_feedback_and_resistance_controls(cx,
            radiation.coefficient_mean_slopes(), radiation.local_reference_slopes(),
            radiation.reference_feedback().ok_or_else(|| bad("ambient derivative lost reference feedback"))?,
            weights, entries, max_pairs).map_err(lower)?;
        let physical = radiation.remap_gradient(cx, full.transport, data.linear).map_err(lower)?;
        let ambient = bound.natural_ambient(&physical.convection)?;
        let air = bound.air_rows(cx, setup, solved, &physical.convection)?;
        for (i, surface) in declaration.surfaces.iter().enumerate() {
            radiation_rows.push(row("radiation-reservoir-temperature", &surface.name, i, "K",
                physical.reservoir_temperatures[i])?);
            radiation_rows.push(row("radiation-emissivity", &surface.name, i, "1",
                physical.emissivities[i])?);
        }
        (physical.convection, Some(full.resistance_multipliers), ambient, air)
    } else if has_natural || has_air {
        let memory = spec.budgets.as_ref().map_or(0, |b| b.memory_bytes);
        let entries = usize::try_from(memory / 32).unwrap_or(usize::MAX);
        // This is the existing physical boundary binder used by the coupled
        // thermal adjoint. Requery natural cards and physically re-march air
        // paths at the unchanged field; final frozen rows alone are insufficient.
        let bound = convective_feedback::bind(cx, spec, solved, entries)?;
        let names: Vec<_> = bound.names.iter().map(String::as_str).collect();
        let response = contacts.linearize(cx, ConductionProblem { boundary: &bound.boundary, ..problem },
            &solved.solution.temperature, data.linear, &names, max_pairs)?;
        let full = response.pullback_with_robin_feedback_and_resistance_controls(cx,
            &bound.natural.wall_slopes, &zeros(names.len())?, &bound.feedback,
            weights, entries, max_pairs).map_err(lower)?;
        let ambient = bound.natural_ambient(&full.transport)?;
        // This also rejects a declared air path missing from retained data.
        let air = bound.air_rows(cx, setup, solved, &full.transport)?;
        (full.transport, Some(full.resistance_multipliers), ambient, air)
    } else {
        let response = contacts.linearize(cx, problem, &solved.solution.temperature,
            data.linear, &[], max_pairs)?;
        if contact_requested {
            let full = response.pullback_with_resistance_controls(cx, weights, &[], &[], max_pairs)
                .map_err(lower)?;
            (full.transport, Some(full.resistance_multipliers), BTreeMap::new(), Vec::new())
        } else {
            (response.pullback(cx, weights, &[], &[]).map_err(lower)?, None, BTreeMap::new(), Vec::new())
        }
    };
    let contact_rows = if contact_requested {
        contact_controls::rows_from_gradients(cx, spec,
            resistance_multipliers.ok_or_else(|| bad("complete contact control lost its numerical owner"))?,
            "R''(T_A,T_B) -> s R''(T_A,T_B), evaluated at s=1")?
    } else {
        Vec::new()
    };
    if !gradient.relative_residual.is_finite()
        || gradient.relative_residual > data.linear.tolerance
    {
        return Err(bad("complete source-backed contact dual did not meet its true-residual gate"));
    }
    poll(cx)?;
    Ok(Pullback { gradient, contact_rows, natural_ambient, air_rows, radiation_rows })
}
