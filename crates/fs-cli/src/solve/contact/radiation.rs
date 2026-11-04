//! Native ambient, closed and open enclosure radiation retain live source contacts.

use fs_conduction::interface::temperature::TemperatureContacts;
use fs_conduction::{ConductionProblem, SolveConfig};
use fs_conduction::radiation::{CoupledRadiationConfig, CoupledRadiationReport,
    GrayDiffuseEnclosure, OpenDiffuseEnclosure};
use fs_conduction::interface::temperature::steady::enclosure::ContactEnclosureWork;
use fs_conduction::{ConductionSolution, RobinFlux, ThermalBoundary};
use fs_exec::Cx;

use super::super::SolveRefusal;
use super::super::radiation::{LoweredRadiation, SolidSolution};
use super::{Prepared, Work, bad, finish_picard, lower, newton, picard_controls, poll};

const SCOPE: &str = "Source-backed matching P1 contact and heterogeneous k(T) use the original complete damped Newton residual and analytic tangent for every ambient-radiation response. The original area-mean gray radiation driver retains its secant iteration, source emissivity, patch temperature and watt gates. The unchanged final field passes a separate freshly queried contact and exact endpoint-radiation nodal residual and heat balance. Newton history describes the last applied radiative secant, not a monolithic radiation Newton solve. Cumulative thermal, nonlinear and Krylov counts include every radiation response and enclosing natural/air callback. Original convective heat alone enters transported air; contact heat remains internal and radiation goes to its declared reservoir. No contact/radiation derivative, arbitrary-law convergence or uniqueness, continuum bound, moving interface or new material authority is inferred.";
const PICARD_SCOPE: &str = "Original source-backed matching P1 contact laws, including temperature/pressure surfaces, use the existing bounded contact Picard iteration inside the original area-mean gray radiation driver. Each contact response retains heterogeneous k(T), its pinned source claims and fixed manufactured context. The unchanged final field passes the same freshly queried contact and exact endpoint-radiation nodal residual and independent heat-balance gates as the Newton composition. Counted spatial responses include all inner contact iterations and remain separate from radiation callbacks. Original convection alone enters transported air; contact transfer is internal and radiation exchanges with its declared reservoir. This nested fixed point supplies no arbitrary-law convergence, uniqueness, source derivative, continuum bound or new material authority.";
const ENCLOSURE_SCOPE: &str = "Live matching source-contact laws and heterogeneous k(T) use the original complete Newton or bounded contact Picard response inside the existing closed/open area-mean radiosity iteration. Every response retains its pinned source selection and manufactured context. The unchanged final field requeries contacts and adds fresh physical-surface radiative nodal loads directly to original conduction/convection rows, including imposed reservoir reactions. Full free residual, heat closure, per-patch mismatch and radiosity gates retain their original watt budgets. All inner spatial, nonlinear and Krylov work is counted, including the initial response and discarded trials. Original convection alone heats transported air; contact and closed enclosure transfer remain internal, while open-cavity radiation exchanges with its explicit black reservoir. That reservoir has no solid degree of freedom and its heat never cancels the solid's external radiation. Numerical history describes the final applied load, without a hidden final solve. No contact/radiation derivative, monolithic radiation tangent, arbitrary-law convergence, uniqueness, continuum bound or new material authority is inferred.";

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
    if radiation.is_pointwise() {
        return Err(bad(
            "steady temperature-contact radiation requires area-mean ambient or closed/open enclosure radiation",
        ));
    }
    if radiation.is_enclosure() {
        return solve_enclosure(cx, prepared, problem, contacts, initial, config,
            tolerance, max_pairs, radiation, work);
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

enum BoundEnclosure {
    Closed(GrayDiffuseEnclosure),
    Open(OpenDiffuseEnclosure, f64),
}

impl BoundEnclosure {
    #[allow(clippy::too_many_arguments)]
    fn finish(self, conduction: ConductionSolution, radiation: CoupledRadiationReport,
        applied_boundary: ThermalBoundary, convection: Vec<RobinFlux>, convective_out_w: f64,
        work: ContactEnclosureWork, tolerance: f64,
    ) -> Result<SolidSolution, SolveRefusal> {
        match self {
            Self::Closed(model) => SolidSolution::from_contact_enclosure(conduction, model,
                radiation, applied_boundary, convection, convective_out_w, work, tolerance),
            Self::Open(model, ambient) => SolidSolution::from_contact_open_enclosure(conduction,
                model, ambient, radiation, applied_boundary, convection, convective_out_w, work, tolerance),
        }.map_err(lower)
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_enclosure(
    cx: &Cx<'_>, prepared: &Prepared, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, initial: &[f64], config: &SolveConfig,
    tolerance: f64, max_pairs: usize, radiation: &LoweredRadiation, work: &mut Work,
) -> Result<(SolidSolution, String), SolveRefusal> {
    let model = if let Some(ambient) = radiation.open_ambient_temperature_k() {
        BoundEnclosure::Open(radiation.bind_open_enclosure(cx, problem).map_err(lower)?
            .ok_or_else(|| bad("live contact enclosure lost its open model"))?, ambient)
    } else {
        BoundEnclosure::Closed(radiation.bind_enclosure(cx, problem).map_err(lower)?
            .ok_or_else(|| bad("live contact enclosure lost its closed model"))?)
    };
    let coupling = CoupledRadiationConfig {
        max_iterations: radiation.config.max_iterations,
        surface_temperature_atol_k: radiation.config.temperature_tolerance_k,
        surface_temperature_rtol: 0.0,
        relaxation: radiation.config.relaxation,
    };
    let initial = work.last_temperature.as_deref().unwrap_or(initial);
    if newton::supported(prepared, config) {
        let solver = newton::controls(config, tolerance, max_pairs, initial.len())?;
        let accepted = match &model {
            BoundEnclosure::Closed(model) => contacts.solve_steady_newton_with_gray_diffuse_enclosure(
                cx, problem, initial, solver, model, coupling),
            BoundEnclosure::Open(model, ambient) => contacts.solve_steady_newton_with_open_gray_diffuse_enclosure(
                cx, problem, initial, solver, model, *ambient, coupling),
        }.map_err(lower)?;
        let last_krylov = accepted.contact.newton.history.iter().try_fold(0usize, |n, row|
            n.checked_add(row.linear_iterations).ok_or_else(|| bad("enclosure Krylov work overflowed")))?;
        work.record(earlier(accepted.work.thermal_responses, 1)?,
            earlier(accepted.work.nonlinear_iterations, accepted.contact.newton.iterations)?,
            earlier(accepted.work.krylov_iterations, last_krylov)?)?;
        let composition = format!(concat!(
            ",\"enclosure_radiation\":{{\"method\":\"area-mean-radiosity-with-complete-contact-newton\",",
            "\"radiation_responses\":{},\"thermal_responses\":{},\"residual_evaluations\":{},",
            "\"jacobian_assemblies\":{},\"physical_endpoint_reassembled\":true}}"),
            accepted.radiation.iterations + 1, accepted.work.thermal_responses,
            accepted.work.residual_evaluations, accepted.work.jacobian_assemblies);
        let (solid, receipt) = newton::finish(cx, problem, accepted.contact, tolerance,
            work, ENCLOSURE_SCOPE, &composition)?;
        let solid = model.finish(solid.conduction,
            accepted.radiation, accepted.applied_boundary, accepted.convective_robin_fluxes,
            accepted.convective_out_w, accepted.work, radiation.config.balance_tolerance_w)?;
        Ok((solid, receipt))
    } else {
        let (solver, contacts_config) = picard_controls(config.clone(), tolerance);
        let accepted = match &model {
            BoundEnclosure::Closed(model) => contacts.solve_steady_with_gray_diffuse_enclosure(
                cx, problem, initial, solver, contacts_config, model, coupling),
            BoundEnclosure::Open(model, ambient) => contacts.solve_steady_with_open_gray_diffuse_enclosure(
                cx, problem, initial, solver, contacts_config, model, *ambient, coupling),
        }.map_err(lower)?;
        work.record(earlier(accepted.work.thermal_responses, accepted.contact.iterations)?,
            earlier(accepted.work.nonlinear_iterations, accepted.contact.nonlinear_iterations)?,
            earlier(accepted.work.krylov_iterations, accepted.contact.krylov_iterations)?)?;
        let composition = format!(concat!(
            ",\"enclosure_radiation\":{{\"method\":\"area-mean-radiosity-with-bounded-contact-picard\",",
            "\"radiation_responses\":{},\"thermal_responses\":{},\"physical_endpoint_reassembled\":true}}"),
            accepted.radiation.iterations + 1, accepted.work.thermal_responses);
        let (solid, receipt) = finish_picard(cx, accepted.contact, tolerance, work,
            ENCLOSURE_SCOPE, &composition)?;
        let solid = model.finish(solid.conduction,
            accepted.radiation, accepted.applied_boundary, accepted.convective_robin_fluxes,
            accepted.convective_out_w, accepted.work, radiation.config.balance_tolerance_w)?;
        Ok((solid, receipt))
    }
}
