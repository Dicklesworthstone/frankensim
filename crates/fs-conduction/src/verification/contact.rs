//! Native matching-interface solves with a broken-H1 continuum mean bound.
//! Bind the original declarations once for both FEM solves. The verifier gets
//! the same face pairs and per-face resistances, not an aggregate K/W value,
//! an inferred temperature continuity condition, or fixed Robin references.
use super::*;
use crate::{InterfaceSurface, ThermalInterfaces};

pub(super) fn admit_contacts(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, surfaces: &[InterfaceSurface], flux: FluxBudget,
) -> Result<(Admitted, ThermalInterfaces)> {
    poll(cx)?;
    // Matching faces have unique ownership. Bound the declarations BEFORE
    // cloning their retained cards; the original binder validates ownership,
    // names, external-BC conflicts and every undeclared coincident pair.
    let cap = problem.mesh.boundary().len() / 2;
    if surfaces.len() > cap { return Err(TetError::Budget.into()); }
    let mut pairs = 0usize;
    for surface in surfaces {
        poll(cx)?;
        pairs = pairs.checked_add(surface.face_pairs().len()).ok_or(TetError::Budget)?;
        if pairs > cap { return Err(TetError::Budget.into()); }
    }
    let mut admitted = admit(cx, problem, flux)?;
    let interfaces = ThermalInterfaces::new(problem.mesh, problem.boundary, surfaces.to_vec())?;
    poll(cx)?;
    for surface in surfaces {
        for (i, pair) in surface.face_pairs().iter().enumerate() {
            poll(cx)?;
            let resistance = if let Some(mapped) = surface.face_resistances() {
                mapped.get(i)
            } else { surface.resistance() }
                .ok_or(TetError::Invalid("missing original per-face resistance"))?
                .value_m2_k_per_w();
            // Admission stores boundary rows in the original mesh-slot order.
            // Keep each declared pair WITH its resistance: independently
            // sorting either list would silently verify a different operator.
            let a = admitted.boundary.get(pair.side_a)
                .ok_or(TetError::Invalid("missing bound side-A face"))?.vertices;
            let b = admitted.boundary.get(pair.side_b)
                .ok_or(TetError::Invalid("missing bound side-B face"))?.vertices;
            admitted.boundary[pair.side_a].condition = BoundaryCondition::Contact { partner: b, resistance };
            admitted.boundary[pair.side_b].condition = BoundaryCondition::Contact { partner: a, resistance };
        }
    }
    poll(cx)?;
    Ok((admitted, interfaces))
}

/// Solve and bound whole-volume mean temperature across declared matching
/// finite-resistance contacts. Each [`InterfaceSurface`] is the original
/// physical declaration, including a mapped resistance for every face pair.
/// The original interface binder checks complete ownership; the original
/// conduction solver supplies both the primal and same-contact unit-source
/// dual. Their reports and card provenance remain in the returned solutions.
///
/// All linear/tensor/source/boundary validity restrictions of
/// [`solve_with_mean_bound`] remain. Temperature traces are independent on
/// either side, and the majorant includes their resistance-weighted jump.
/// The result is conditional on the nominal fixed positive contact values;
/// it is not a contact-law, material-uncertainty, nonlinear, nonmatching,
/// physical-validation or maximum-temperature certificate.
///
/// # Errors
/// Invalid/undeclared contact, original physical/solver failures, verifier
/// refusal, exhausted work, or cancellation return no partial bound. Matching
/// binding retains its existing non-interruptible geometry pass, surrounded
/// by cancellation checks and preceded by mesh/declaration admission.
pub fn solve_with_contact_mean_bound(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, surfaces: &[InterfaceSurface], config: MeanSolveConfig,
) -> Result<MeanTemperatureSolution> {
    let (admitted, interfaces) = admit_contacts(cx, problem, surfaces, config.flux)?;
    let primal = crate::solve_with_interfaces(cx, problem, &interfaces, config.primal)?;
    let result = bound_field(cx, problem, &admitted, &primal.temperature,
        config.dual, config.flux, Some(&interfaces))?;
    Ok(MeanTemperatureSolution { primal, dual: result.dual, bound: result.bound })
}

/// Bound an existing discontinuous-at-contact P1 field without another primal
/// solve. The field need not satisfy equilibrium; its algebraic and contact
/// residuals remain in the majorant and goal correction. Exact prescribed
/// Dirichlet values, complete finite field length and the original contact
/// family are required before the dual solve. No primal report is invented.
///
/// # Errors
/// Same original physics, contact, numerical and resource refusals as
/// [`solve_with_contact_mean_bound`], plus invalid supplied fields. Source,
/// boundary data and all per-face material inputs remain unchanged.
pub fn bound_temperature_mean_with_contacts(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, surfaces: &[InterfaceSurface],
    temperature: &[f64], dual_config: SolveConfig, flux: FluxBudget,
) -> Result<MeanFieldBound> {
    let (admitted, interfaces) = admit_contacts(cx, problem, surfaces, flux)?;
    bound_field(cx, problem, &admitted, temperature, dual_config, flux, Some(&interfaces))
}
