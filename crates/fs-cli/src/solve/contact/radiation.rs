//! Native ambient radiation reuses complete source-contact Newton responses.

use fs_conduction::interface::temperature::TemperatureContacts;
use fs_conduction::{ConductionProblem, SolveConfig};
use fs_exec::Cx;

use super::super::SolveRefusal;
use super::super::radiation::{LoweredRadiation, SolidSolution};
use super::{Prepared, Work, bad, lower, newton, poll};

const SCOPE: &str = "Source-backed matching P1 contact and heterogeneous k(T) use the original complete damped Newton residual and analytic tangent for every ambient-radiation response. The original area-mean gray radiation driver retains its secant iteration, source emissivity, patch temperature and watt gates. The unchanged final field passes a separate freshly queried contact and exact endpoint-radiation nodal residual and heat balance. Newton history describes the last applied radiative secant, not a monolithic radiation Newton solve. Cumulative thermal, nonlinear and Krylov counts include every radiation response and enclosing natural/air callback. Original convective heat alone enters transported air; contact heat remains internal and radiation goes to its declared reservoir. No contact/radiation derivative, arbitrary-law convergence or uniqueness, continuum bound, moving interface or new material authority is inferred.";

#[allow(clippy::too_many_arguments)]
pub(super) fn solve(
    cx: &Cx<'_>,
    prepared: &Prepared,
    problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>,
    initial: &[f64],
    config: &SolveConfig,
    tolerance: f64,
    max_pairs: usize,
    radiation: &LoweredRadiation,
    work: &mut Work,
) -> Result<(SolidSolution, String), SolveRefusal> {
    poll(cx)?;
    if radiation.is_pointwise() || radiation.is_enclosure() || !newton::supported(prepared, config)
    {
        return Err(bad(
            "steady temperature-contact radiation requires the complete Newton owner, scalar/continuous linear source laws and the area-mean ambient radiation model",
        ));
    }
    let accepted = contacts
        .solve_steady_newton_with_ambient_radiation(
            cx,
            problem,
            work.last_temperature.as_deref().unwrap_or(initial),
            newton::controls(config, tolerance, max_pairs, initial.len())?,
            &radiation.patches,
            radiation.config,
        )
        .map_err(lower)?;
    let last_krylov = accepted
        .contact
        .newton
        .history
        .iter()
        .try_fold(0usize, |total, row| {
            total
                .checked_add(row.linear_iterations)
                .ok_or_else(|| bad("counted radiative contact Krylov work overflowed"))
        })?;
    let earlier = |total: usize, last: usize| {
        total
            .checked_sub(last)
            .ok_or_else(|| bad("radiative contact work lost an earlier response"))
    };
    // The shared formatter counts the final complete Newton solve. Account for
    // every earlier radiation response once, including its rejected trials.
    work.record(
        earlier(accepted.radiation.iterations, 1)?,
        earlier(
            accepted.radiation.solid_iterations,
            accepted.contact.newton.iterations,
        )?,
        earlier(accepted.radiation.krylov_iterations, last_krylov)?,
    )?;
    let composition = format!(
        concat!(
            ",\"ambient_radiation\":{{\"method\":\"area-mean-secant-with-complete-contact-newton\",",
            "\"thermal_responses\":{},\"residual_evaluations\":{},\"jacobian_assemblies\":{},",
            "\"physical_endpoint_reassembled\":true}}"
        ),
        accepted.radiation.iterations, accepted.residual_evaluations, accepted.jacobian_assemblies
    );
    let (solid, receipt) = newton::finish(
        cx,
        problem,
        accepted.contact,
        tolerance,
        work,
        SCOPE,
        &composition,
    )?;
    Ok((
        SolidSolution::from_radiative_endpoint(
            solid.conduction,
            accepted.physical_boundary,
            accepted.convective_robin_fluxes,
            accepted.convective_out_w,
            accepted.radiation,
        ),
        receipt,
    ))
}
