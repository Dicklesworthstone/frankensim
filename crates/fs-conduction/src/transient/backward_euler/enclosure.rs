//! Capacity storage with implicit reflected radiation and immutable temperature history.
//! The existing linear and full k(T) Newton owners remain the only solid solvers.
//! Shared enclosure traces separate internal exchange from original external heat.

use fs_exec::Cx;
use super::{
    BackwardEuler, EndpointEvaluation, NonlinearStepConfig, NonlinearStepSolution,
    StepConfig, StepSolution, finite, invalid, poll,
};
use crate::{
    AmbientRadiationConfig, ConductionError, ConductionProblem, RobinFlux, ScalarField,
    ThermalBc, ThermalBoundary, ThermalInterfaces,
    radiation::{
        EnclosureExchangeReport, GrayDiffuseEnclosure, RadiosityReport,
        enclosure_exchange::{self as exchange, Patch, PatchLaw},
    },
};

/// A supplied temperature-capacity candidate audited under fresh enclosure exchange.
/// No Newton or Krylov solve, timestep advance, or history mutation is performed.
#[derive(Debug, Clone)]
pub struct EnclosureEndpointEvaluation {
    /// Complete physical residual, stored energy, original convection/Neumann
    /// exchanges and prescribed reaction. Radiation is internal signed exchange.
    pub physical: EndpointEvaluation,
    /// Fresh radiosity on the supplied candidate's exact patch means.
    pub radiosity: RadiosityReport,
}

/// One accepted capacity step with all original external heat kept separate.
#[derive(Debug, Clone)]
pub struct EnclosureStepSolution {
    /// Final inner solve, whose Robin/Neumann heat includes applied radiation.
    pub conduction: StepSolution,
    /// Exact final inner nonlinear solve for explicitly admitted k(T).
    pub nonlinear: Option<NonlinearStepSolution>,
    /// Rejected Newton line-search trials across every inner enclosure solve.
    pub nonlinear_backtracks: usize,
    /// Frozen boundary actually used by the final inner solid solve.
    pub combined_boundary: ThermalBoundary,
    /// Original convection only, in boundary declaration order.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Original convection leaving the solid, W; this alone heats coolant.
    pub convective_out_w: f64,
    /// Authored outward Neumann heat, W, excluding internal radiation.
    pub prescribed_neumann_out_w: f64,
    /// Fresh reflected radiation and actual cumulative solid/Krylov work.
    pub radiation: EnclosureExchangeReport,
    /// Full free capacity residual on the fresh physical boundary, J.
    pub physical_residual_norm_j: f64,
    /// Original fixed predictor-based physical residual target, J.
    pub physical_residual_tolerance_j: f64,
    /// Actual storage minus dt times complete physical net input, J.
    pub physical_energy_residual_j: f64,
    /// Fresh prescribed reaction including fixed-node capacity, W inward.
    pub physical_dirichlet_in_w: f64,
}

/// One complete thermal response to the shared enclosure's applied boundary.
/// Retained state belongs to the numerical response, not the fresh audit.
pub(crate) struct EnclosureStepResponse<S> {
    pub conduction: StepSolution,
    pub nonlinear: Option<NonlinearStepSolution>,
    pub retained: S,
    pub thermal_responses: usize,
    pub nonlinear_updates: usize,
    pub nonlinear_backtracks: usize,
    pub krylov_iterations: usize,
}

pub(crate) struct EnclosureStepOutcome<S, A> {
    pub response: EnclosureStepResponse<S>,
    pub physical: EndpointEvaluation,
    pub physical_state: A,
    pub combined_boundary: ThermalBoundary,
    pub radiation: EnclosureExchangeReport,
    pub physical_residual_tolerance_j: f64,
    pub thermal_responses: usize,
    pub nonlinear_updates: usize,
    pub nonlinear_backtracks: usize,
    pub krylov_iterations: usize,
}

impl BackwardEuler<'_> {
    /// Re-evaluate a candidate under original convection and reflected radiation.
    ///
    /// Arbitrary off-equilibrium candidates retain their actual residual; the
    /// outer owner applies its unchanged physical acceptance targets. The
    /// existing endpoint assembler handles C(T-T_old), full k(T), contacts and
    /// changing prescribed temperatures. Fresh radiation enters that same bath
    /// reaction; returned convection and Neumann heat exclude internal exchange.
    ///
    /// # Errors
    /// Refuses invalid step/fields, mismatched traces, invalid physical laws,
    /// failed radiosity heat gates, nonfinite arithmetic or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_prescribed_enclosure_endpoint(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], candidate: &[f64],
        dt_s: f64, config: StepConfig, enclosure: &GrayDiffuseEnclosure,
        radiation_config: AmbientRadiationConfig,
    ) -> Result<EnclosureEndpointEvaluation, ConductionError> {
        self.admit_endpoint_step(cx, problem, old, dt_s, config)?;
        radiation_config.validate()?;
        let patches = exchange::bind(cx, problem, enclosure, candidate)?;
        let means = exchange::means(cx, problem, enclosure, candidate)?;
        let fresh = exchange::radiosity(cx, enclosure, &means, radiation_config)?;
        self.enclosure_endpoint_from_radiosity(
            cx, problem, interfaces, old, candidate, dt_s, config, &patches, fresh,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn enclosure_endpoint_from_radiosity(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], candidate: &[f64],
        dt_s: f64, config: StepConfig, patches: &[Patch], fresh: RadiosityReport,
    ) -> Result<EnclosureEndpointEvaluation, ConductionError> {
        let (endpoint, ()) = endpoint_from_radiosity(
            cx, problem, candidate, dt_s, patches, fresh,
            &mut |combined, field| self.evaluate_prescribed_endpoint(
                cx, combined, interfaces, old, field, dt_s, config, None,
            ).map(|physical| (physical, ())),
        )?;
        Ok(endpoint)
    }

    /// Advance genuine volumetric-capacity storage with closed-enclosure exchange.
    ///
    /// Every inner solve receives the SAME old temperature field. Only this
    /// method's complete accepted endpoint may become the next physical history.
    /// Supply the explicit nonlinear policy to retain the full k(T) Newton
    /// tangent; the linear entry point still refuses temperature-dependent k.
    /// Exact named Robin or Neumann traces share one reflected-radiation owner
    /// with enthalpy storage, without introducing a density or phase chart.
    ///
    /// Fixed physical residual targets use the old free temperatures and new
    /// prescribed temperatures as predictor. Inner policies only tighten.
    /// Surface change, applied/fresh patch heat, complete nodal residual and
    /// independent energy balance must all satisfy their original allowances.
    ///
    /// # Errors
    /// Invalid input, unmatched traces, material/reflection failures, exhausted
    /// work, cancellation or incomplete physical balance publish no endpoint.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_prescribed_with_enclosure_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        step_config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
        enclosure: &GrayDiffuseEnclosure, radiation_config: AmbientRadiationConfig,
    ) -> Result<EnclosureStepSolution, ConductionError> {
        let accepted = self.advance_enclosure_with_responses(
            cx, problem, old, dt_s, step_config, nonlinear, enclosure, radiation_config,
            None, 1,
            |combined, inner_step, policy, _| {
                let (conduction, nonlinear) = match policy {
                    Some(policy) => {
                        let solved = self.advance_nonlinear_prescribed(
                            cx, combined, interfaces, old, dt_s, inner_step, policy,
                        )?;
                        (solved.step.clone(), Some(solved))
                    }
                    None => (self.advance_prescribed(
                        cx, combined, interfaces, old, dt_s, inner_step,
                    )?, None),
                };
                Ok(EnclosureStepResponse {
                    nonlinear_updates: nonlinear.as_ref().map_or(0, |step| step.nonlinear_iterations),
                    nonlinear_backtracks: nonlinear.as_ref().map_or(0, |step| step.backtracks),
                    krylov_iterations: conduction.krylov_iterations,
                    conduction, nonlinear, retained: (), thermal_responses: 1,
                })
            },
            |combined, field| self.evaluate_prescribed_endpoint(
                cx, combined, interfaces, old, field, dt_s, step_config, None,
            ).map(|physical| (physical, ())),
        )?;
        let physical = accepted.physical;
        Ok(EnclosureStepSolution {
            conduction: accepted.response.conduction,
            nonlinear: accepted.response.nonlinear,
            nonlinear_backtracks: accepted.nonlinear_backtracks,
            combined_boundary: accepted.combined_boundary,
            convective_robin_fluxes: physical.convective_robin_fluxes,
            convective_out_w: physical.convective_out_w,
            prescribed_neumann_out_w: physical.neumann_out_w,
            radiation: accepted.radiation,
            physical_residual_norm_j: physical.residual_norm_j,
            physical_residual_tolerance_j: accepted.physical_residual_tolerance_j,
            physical_energy_residual_j: physical.energy_residual_j,
            physical_dirichlet_in_w: physical.dirichlet_in_w,
        })
    }

    /// The one capacity enclosure iterator. Response/evaluation hooks select
    /// the complete thermal/contact owner; radiosity, shifted boundary laws,
    /// physical heat separation and convergence remain owned here.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(crate) fn advance_enclosure_with_responses<S, A>(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>, old: &[f64], dt_s: f64,
        step_config: StepConfig, nonlinear: Option<NonlinearStepConfig>,
        enclosure: &GrayDiffuseEnclosure, radiation_config: AmbientRadiationConfig,
        absolute_target_j: Option<f64>, responses_per_trial: usize,
        mut respond: impl FnMut(ConductionProblem<'_>, StepConfig, Option<NonlinearStepConfig>, f64)
            -> Result<EnclosureStepResponse<S>, ConductionError>,
        mut evaluate: impl FnMut(ConductionProblem<'_>, &[f64])
            -> Result<(EndpointEvaluation, A), ConductionError>,
    ) -> Result<EnclosureStepOutcome<S, A>, ConductionError> {
        let dofs = self.admit_endpoint_step(cx, problem, old, dt_s, step_config)?;
        radiation_config.validate()?;
        if responses_per_trial == 0 || absolute_target_j.is_some_and(|value| !value.is_finite() || value <= 0.0) {
            return Err(invalid("enclosure response budget and explicit residual target must be positive"));
        }
        let max_responses = radiation_config.max_iterations.checked_mul(responses_per_trial)
            .ok_or_else(|| invalid("enclosure capacity response allowance overflows"))?;
        let max_krylov = max_responses.checked_mul(step_config.linear.max_iterations)
            .ok_or_else(|| invalid("enclosure capacity Krylov work budget overflow"))?;
        let max_updates = max_responses.checked_mul(nonlinear.map_or(0, |policy| policy.max_iterations))
            .ok_or_else(|| invalid("enclosure capacity Newton work budget overflow"))?;
        let max_backtracks = max_updates.checked_mul(nonlinear.map_or(0, |policy| policy.line_search.max_backtracks))
            .ok_or_else(|| invalid("enclosure capacity backtrack work budget overflow"))?;
        if let Some(policy) = nonlinear {
            policy.validate()?;
            policy.line_search.max_backtracks.checked_add(1)
                .and_then(|n| n.checked_mul(max_updates))
                .ok_or_else(|| invalid("enclosure capacity Newton trial budget overflow"))?;
        }
        let mut predictor = old.to_vec();
        for (index, &vertex) in dofs.fixed().iter().enumerate() {
            if index % 512 == 0 { poll(cx, index)?; }
            predictor[vertex] = dofs.prescribed()[vertex];
        }
        let patches = exchange::bind(cx, problem, enclosure, &predictor)?;
        let prepared = enclosure.prepare_radiosity(cx)?;
        let initial_means = exchange::means(cx, problem, enclosure, &predictor)?;
        let initial_radiosity =
            exchange::prepared_radiosity(cx, &prepared, &initial_means, radiation_config)?;
        let (initial, _) = endpoint_from_radiosity(
            cx, problem, &predictor, dt_s, &patches, initial_radiosity, &mut evaluate,
        )?;
        let initial_norm_j = initial.physical.residual_norm_j;
        let target = match absolute_target_j {
            Some(target) => target,
            None => match nonlinear {
                Some(policy) => finite(policy.residual_atol_j + policy.residual_rtol * initial_norm_j)?,
                None if dofs.n() == 0 => 0.0,
                None => finite(step_config.linear.tolerance * initial_norm_j
                    + 0.01 * step_config.energy_tolerance_j / (dofs.n() as f64).sqrt())?,
            },
        };
        if absolute_target_j.is_some() && target * 0.125 <= 0.0 {
            return Err(invalid("enclosure contact inner residual allowance is not representable"));
        }
        let mut driving = initial.radiosity.surface_temperatures_k;
        let mut thermal_responses = 0usize;
        let mut solid_iterations = 0usize;
        let mut krylov_iterations = 0usize;
        let mut nonlinear_backtracks = 0usize;
        let mut last_change = 0.0_f64;
        let mut last_mismatch = 0.0_f64;
        for iteration in 0..radiation_config.max_iterations {
            poll(cx, iteration)?;
            let applied = exchange::prepared_radiosity(cx, &prepared, &driving, radiation_config)?;
            let boundary = exchange::shifted_boundary(problem, &patches, &applied)?;
            let combined = ConductionProblem { boundary: &boundary, ..problem };
            let (initial, _) = evaluate(combined, &predictor)?;
            let inner_norm_j = initial.residual_norm_j;
            let mut inner_step = step_config;
            let inner_nonlinear = match nonlinear {
                Some(mut policy) => {
                    policy.residual_atol_j = policy.residual_atol_j.min(target * 0.125);
                    policy.residual_rtol =
                        tighter_relative(policy.residual_rtol, target * 0.125, inner_norm_j)?;
                    Some(policy)
                }
                None => {
                    inner_step.linear.tolerance =
                        tighter_relative(inner_step.linear.tolerance, target * 0.25, inner_norm_j)?;
                    None
                }
            };
            let response = respond(combined, inner_step, inner_nonlinear, target * 0.25)?;
            if response.thermal_responses == 0 || response.thermal_responses > responses_per_trial {
                return Err(invalid("enclosure thermal response count exceeds its declared allowance"));
            }
            count(&mut thermal_responses, response.thermal_responses, max_responses)?;
            count(&mut solid_iterations, response.nonlinear_updates, max_updates)?;
            count(&mut nonlinear_backtracks, response.nonlinear_backtracks, max_backtracks)?;
            count(&mut krylov_iterations, response.krylov_iterations, max_krylov)?;
            let conduction = &response.conduction;
            let actual = exchange::means(cx, problem, enclosure, &conduction.temperature)?;
            let fresh = exchange::prepared_radiosity(cx, &prepared, &actual, radiation_config)?;
            let (audited, physical_state) = endpoint_from_radiosity(
                cx, problem, &conduction.temperature, dt_s, &patches, fresh, &mut evaluate,
            )?;
            let physical = audited.physical;
            let fresh = audited.radiosity;
            let mut heat_ok = true;
            last_mismatch = 0.0;
            for (index, patch) in patches.iter().enumerate() {
                poll(cx, index)?;
                let actual_applied = match patch.law {
                    PatchLaw::Convection { .. } => {
                        let name = &problem.boundary.region_names()[patch.region];
                        let full = conduction.robin_fluxes.iter().find(|row| row.region == *name)
                            .ok_or_else(|| invalid("enclosure inner solve omitted its Robin patch"))?;
                        let original = physical.convective_robin_fluxes.iter()
                            .find(|row| row.region == *name)
                            .ok_or_else(|| invalid("enclosure audit omitted original convection"))?;
                        finite(full.heat_rate_w - original.heat_rate_w)?
                    }
                    PatchLaw::PrescribedFlux { outward_flux } => {
                        let ThermalBc::Neumann { outward_flux: ScalarField::Uniform(combined) }
                            = boundary.conditions()[patch.region]
                        else { return Err(invalid("enclosure changed Neumann patch ownership")); };
                        finite((combined - outward_flux) * enclosure.surfaces()[index].area_m2())?
                    }
                };
                if finite(actual_applied - applied.net_outward_heat_w[index])?.abs()
                    > exchange::heat_tolerance(radiation_config,
                        actual_applied, applied.net_outward_heat_w[index])? {
                    return Err(invalid("enclosure capacity boundary cannot separate authored heat and radiation"));
                }
                let mismatch = finite(actual_applied - fresh.net_outward_heat_w[index])?.abs();
                last_mismatch = last_mismatch.max(mismatch);
                heat_ok &= mismatch <= exchange::heat_tolerance(radiation_config,
                    actual_applied, fresh.net_outward_heat_w[index])?;
            }
            last_change = actual.iter().zip(&driving).try_fold(0.0_f64, |largest, (&a, &b)| {
                Ok::<_, ConductionError>(largest.max(finite(a - b)?.abs()))
            })?;
            if heat_ok && last_change <= radiation_config.temperature_tolerance_k
                && physical.residual_norm_j <= target
                && physical.energy_residual_j.abs() <= step_config.energy_tolerance_j {
                poll(cx, iteration)?;
                return Ok(EnclosureStepOutcome {
                    response, physical, physical_state, nonlinear_backtracks,
                    combined_boundary: boundary,
                    radiation: EnclosureExchangeReport {
                        radiosity: fresh, applied_heat_w: applied.net_outward_heat_w,
                        max_heat_mismatch_w: last_mismatch, iterations: iteration + 1,
                        solid_iterations, krylov_iterations,
                    },
                    physical_residual_tolerance_j: target,
                    thermal_responses, nonlinear_updates: solid_iterations, krylov_iterations,
                });
            }
            for (value, next) in driving.iter_mut().zip(actual) {
                *value = finite((1.0 - radiation_config.relaxation) * *value
                    + radiation_config.relaxation * next)?;
            }
        }
        Err(ConductionError::Radiation {
            surface: "capacity-enclosure".into(),
            what: format!("{} enclosure responses exhausted; final temperature change {last_change} K and patch mismatch {last_mismatch} W; physical residual or energy may remain unresolved",
                radiation_config.max_iterations),
            fix: "provide adequate explicit enclosure and solid budgets or reduce the timestep".into(),
        })
    }
}

/// One physical heat decomposition for both fixed and freshly queried contacts.
/// The callback assembles capacity on the fresh reflected boundary; only this
/// owner removes internal radiation from authored convection and Neumann heat.
#[allow(clippy::too_many_arguments)]
fn endpoint_from_radiosity<A>(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, candidate: &[f64], dt_s: f64,
    patches: &[Patch], fresh: RadiosityReport,
    evaluate: &mut impl FnMut(ConductionProblem<'_>, &[f64])
        -> Result<(EndpointEvaluation, A), ConductionError>,
) -> Result<(EnclosureEndpointEvaluation, A), ConductionError> {
    let boundary = exchange::shifted_boundary(problem, patches, &fresh)?;
    let (mut physical, state) = evaluate(
        ConductionProblem { boundary: &boundary, ..problem }, candidate,
    )?;
    let (rows, convective_out_w) =
        exchange::convection(cx, problem, &physical.convective_robin_fluxes)?;
    physical.convective_robin_fluxes = rows;
    physical.convective_out_w = convective_out_w;
    physical.neumann_out_w = exchange::prescribed_neumann_heat(cx, problem)?;
    physical.radiation_out_w = fresh.enclosure_energy_closure_w;
    physical.energy_residual_j = finite(physical.stored_energy_change_j - dt_s * (
        physical.source_w + physical.dirichlet_in_w - physical.neumann_out_w
            - physical.convective_out_w - physical.radiation_out_w
    ))?;
    poll(cx, candidate.len())?;
    Ok((EnclosureEndpointEvaluation { physical, radiosity: fresh }, state))
}

fn count(total: &mut usize, added: usize, maximum: usize) -> Result<(), ConductionError> {
    *total = total.checked_add(added).filter(|&value| value <= maximum)
        .ok_or_else(|| invalid("enclosure capacity work exceeds its combined allowance"))?;
    Ok(())
}

fn tighter_relative(requested: f64, allowance_j: f64, initial_norm_j: f64)
    -> Result<f64, ConductionError> {
    let tolerance = if initial_norm_j == 0.0 { requested }
        else { requested.min(allowance_j / initial_norm_j) };
    if tolerance.is_finite() && tolerance > 0.0 { Ok(tolerance) }
    else { Err(invalid("enclosure capacity residual tolerance is not representable")) }
}
