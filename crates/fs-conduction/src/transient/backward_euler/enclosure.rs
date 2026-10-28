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
        let boundary = exchange::shifted_boundary(problem, patches, &fresh)?;
        let mut physical = self.evaluate_prescribed_endpoint(cx,
            ConductionProblem { boundary: &boundary, ..problem }, interfaces,
            old, candidate, dt_s, config, None)?;
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
        Ok(EnclosureEndpointEvaluation { physical, radiosity: fresh })
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
        let dofs = self.admit_endpoint_step(cx, problem, old, dt_s, step_config)?;
        radiation_config.validate()?;
        radiation_config.max_iterations.checked_mul(step_config.linear.max_iterations)
            .ok_or_else(|| invalid("enclosure capacity Krylov work budget overflow"))?;
        if let Some(policy) = nonlinear {
            policy.validate()?;
            policy.line_search.max_backtracks.checked_add(1)
                .and_then(|n| n.checked_mul(policy.max_iterations))
                .and_then(|n| n.checked_mul(radiation_config.max_iterations))
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
        let initial = self.enclosure_endpoint_from_radiosity(cx, problem, interfaces,
            old, &predictor, dt_s, step_config, &patches, initial_radiosity)?;
        let initial_norm_j = initial.physical.residual_norm_j;
        let target = match nonlinear {
            Some(policy) => finite(policy.residual_atol_j + policy.residual_rtol * initial_norm_j)?,
            None if dofs.n() == 0 => 0.0,
            None => finite(step_config.linear.tolerance * initial_norm_j
                + 0.01 * step_config.energy_tolerance_j / (dofs.n() as f64).sqrt())?,
        };
        let mut driving = initial.radiosity.surface_temperatures_k;
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
            let inner_norm_j = self.evaluate_endpoint(
                cx, combined, interfaces, old, &predictor, dt_s, &dofs,
            )?.norm;
            let (conduction, final_nonlinear) = match nonlinear {
                Some(mut policy) => {
                    policy.residual_atol_j = policy.residual_atol_j.min(target * 0.125);
                    policy.residual_rtol =
                        tighter_relative(policy.residual_rtol, target * 0.125, inner_norm_j)?;
                    let solved = self.advance_nonlinear_prescribed(
                        cx, combined, interfaces, old, dt_s, step_config, policy,
                    )?;
                    nonlinear_backtracks = nonlinear_backtracks.checked_add(solved.backtracks)
                        .ok_or_else(|| invalid("enclosure capacity backtrack count overflow"))?;
                    solid_iterations = solid_iterations.checked_add(solved.nonlinear_iterations)
                        .ok_or_else(|| invalid("enclosure capacity Newton count overflow"))?;
                    (solved.step.clone(), Some(solved))
                }
                None => {
                    let mut config = step_config;
                    config.linear.tolerance =
                        tighter_relative(config.linear.tolerance, target * 0.25, inner_norm_j)?;
                    (self.advance_prescribed(cx, combined, interfaces, old, dt_s, config)?, None)
                }
            };
            krylov_iterations = krylov_iterations.checked_add(conduction.krylov_iterations)
                .ok_or_else(|| invalid("enclosure capacity Krylov count overflow"))?;
            let actual = exchange::means(cx, problem, enclosure, &conduction.temperature)?;
            let fresh = exchange::prepared_radiosity(cx, &prepared, &actual, radiation_config)?;
            let audited = self.enclosure_endpoint_from_radiosity(cx, problem, interfaces,
                old, &conduction.temperature, dt_s, step_config, &patches, fresh)?;
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
                return Ok(EnclosureStepSolution {
                    conduction, nonlinear: final_nonlinear, nonlinear_backtracks,
                    combined_boundary: boundary,
                    convective_robin_fluxes: physical.convective_robin_fluxes,
                    convective_out_w: physical.convective_out_w,
                    prescribed_neumann_out_w: physical.neumann_out_w,
                    radiation: EnclosureExchangeReport {
                        radiosity: fresh, applied_heat_w: applied.net_outward_heat_w,
                        max_heat_mismatch_w: last_mismatch, iterations: iteration + 1,
                        solid_iterations, krylov_iterations,
                    },
                    physical_residual_norm_j: physical.residual_norm_j,
                    physical_residual_tolerance_j: target,
                    physical_energy_residual_j: physical.energy_residual_j,
                    physical_dirichlet_in_w: physical.dirichlet_in_w,
                });
            }
            for (value, next) in driving.iter_mut().zip(actual) {
                *value = finite((1.0 - radiation_config.relaxation) * *value
                    + radiation_config.relaxation * next)?;
            }
        }
        Err(ConductionError::Radiation {
            surface: "capacity-enclosure".into(),
            what: format!("{} solid solves exhausted; final temperature change {last_change} K and patch mismatch {last_mismatch} W; physical residual or energy may remain unresolved",
                radiation_config.max_iterations),
            fix: "provide adequate explicit enclosure and solid budgets or reduce the timestep".into(),
        })
    }
}

fn tighter_relative(requested: f64, allowance_j: f64, initial_norm_j: f64)
    -> Result<f64, ConductionError> {
    let tolerance = if initial_norm_j == 0.0 { requested }
        else { requested.min(allowance_j / initial_norm_j) };
    if tolerance.is_finite() && tolerance > 0.0 { Ok(tolerance) }
    else { Err(invalid("enclosure capacity residual tolerance is not representable")) }
}
