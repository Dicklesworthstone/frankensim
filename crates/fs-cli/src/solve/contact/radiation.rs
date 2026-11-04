//! Native ambient radiation reuses the original steady source-contact owners.

use fs_conduction::interface::temperature::TemperatureContacts;
use fs_conduction::{ConductionProblem, SolveConfig};
use fs_exec::Cx;

use super::super::SolveRefusal;
use super::super::radiation::{LoweredRadiation, SolidSolution};
use super::{Prepared, Work, bad, finish_picard, lower, newton, picard_controls, poll};

const SCOPE: &str = "Source-backed matching P1 contact and heterogeneous k(T) use the original complete damped Newton residual and analytic tangent for every ambient-radiation response. The original area-mean gray radiation driver retains its secant iteration, source emissivity, patch temperature and watt gates. The unchanged final field passes a separate freshly queried contact and exact endpoint-radiation nodal residual and heat balance. Newton history describes the last applied radiative secant, not a monolithic radiation Newton solve. Cumulative thermal, nonlinear and Krylov counts include every radiation response and enclosing natural/air callback. Original convective heat alone enters transported air; contact heat remains internal and radiation goes to its declared reservoir. No contact/radiation derivative, arbitrary-law convergence or uniqueness, continuum bound, moving interface or new material authority is inferred.";
const PICARD_SCOPE: &str = "Original source-backed matching P1 contact laws, including temperature/pressure surfaces, use the existing bounded contact Picard iteration inside the original area-mean gray radiation driver. Each contact response retains heterogeneous k(T), its pinned source claims and fixed manufactured context. The unchanged final field passes the same freshly queried contact and exact endpoint-radiation nodal residual and independent heat-balance gates as the Newton composition. Counted spatial responses include all inner contact iterations and remain separate from radiation callbacks. Original convection alone enters transported air; contact transfer is internal and radiation exchanges with its declared reservoir. This nested fixed point supplies no arbitrary-law convergence, uniqueness, source derivative, continuum bound or new material authority.";

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
    if radiation.is_pointwise() || radiation.is_enclosure() {
        return Err(bad(
            "steady temperature-contact radiation requires the area-mean ambient radiation model",
        ));
    }
    if !newton::supported(prepared, config) {
        return solve_picard(cx, problem, contacts, initial, config, tolerance, radiation, work);
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

#[allow(clippy::too_many_arguments)]
fn solve_picard(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, contacts: &TemperatureContacts<'_, '_>,
    initial: &[f64], config: &SolveConfig, tolerance: f64,
    radiation: &LoweredRadiation, work: &mut Work,
) -> Result<(SolidSolution, String), SolveRefusal> {
    let (solver, coupling) = picard_controls(config.clone(), tolerance);
    let accepted = contacts.solve_steady_with_ambient_radiation(cx, problem,
        work.last_temperature.as_deref().unwrap_or(initial), solver, coupling,
        &radiation.patches, radiation.config).map_err(lower)?;
    // The common formatter counts the final Picard response's spatial solves.
    // Count all preceding inner responses exactly once, not merely one per
    // radiative callback; its nonlinear/Krylov counts already span those solves.
    work.record(
        earlier(accepted.thermal_responses, accepted.contact.iterations)?,
        earlier(accepted.radiation.solid_iterations, accepted.contact.nonlinear_iterations)?,
        earlier(accepted.radiation.krylov_iterations, accepted.contact.krylov_iterations)?,
    )?;
    let composition = format!(concat!(
        ",\"ambient_radiation\":{{\"method\":\"area-mean-secant-with-bounded-contact-picard\",",
        "\"radiation_responses\":{},\"thermal_responses\":{},\"physical_endpoint_reassembled\":true}}"),
        accepted.radiation.iterations, accepted.thermal_responses);
    let (solid, receipt) = finish_picard(cx, accepted.contact, tolerance, work,
        PICARD_SCOPE, &composition)?;
    Ok((SolidSolution::from_radiative_endpoint(solid.conduction,
        accepted.physical_boundary, accepted.convective_robin_fluxes,
        accepted.convective_out_w, accepted.radiation)
        .with_radiation_thermal_responses(accepted.thermal_responses), receipt))
}

fn earlier(total: usize, last: usize) -> Result<usize, SolveRefusal> {
    total.checked_sub(last).ok_or_else(|| bad("radiative contact work lost an earlier response"))
}
