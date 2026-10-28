//! Closed gray-diffuse exchange on existing natural enthalpy boundaries.
//!
//! Each patch emits at its P1 area mean and applies a uniform outward flux.
//! The shifted Robin reference `T_air-q_rad/htc` is algebraic, never an air
//! temperature. On a prescribed-flux or insulated trace, radiation adds to the
//! original Neumann load without introducing convection. Every trial solves
//! the same immutable specific-enthalpy
//! history. Acceptance checks the complete nonlinear residual in joules and
//! the independent energy balance; a closed enclosure transfers heat internally.

/// Physical accepted-endpoint tangent, including reflected radiation.
pub mod adjoint;

use fs_exec::Cx;
use fs_solver::{NewtonKrylovState, norm2};

use super::{
    EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution, StepContext,
    finite, poll,
};
use crate::{
    AmbientRadiationConfig, ConductionError, ConductionProblem, RobinFlux, ScalarField, ThermalBc,
    ThermalBoundary, ThermalInterfaces,
    radiation::{GrayDiffuseEnclosure, PreparedRadiosity, RadiosityReport},
};

pub use crate::radiation::EnclosureExchangeReport as EnthalpyEnclosureReport;
use crate::radiation::enclosure_exchange::{self as exchange, Patch, PatchLaw};

/// One accepted physical step with separately identified convection and exchange.
#[derive(Debug, Clone)]
pub struct EnthalpyEnclosureStepSolution {
    /// Accepted h/phase state; radiating Robin and Neumann rows include the
    /// applied enclosure flux. Use the separated accounts below externally.
    pub conduction: EnthalpyStepSolution,
    /// Exact shifted boundary used by the accepted inner solve.
    pub combined_boundary: ThermalBoundary,
    /// Original convection only, in the original boundary declaration order.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Original convection heat, W outward; this alone is passed to air.
    pub convective_out_w: f64,
    /// Original prescribed outward flux integrated over all Neumann faces, W.
    /// Excludes internal enclosure exchange, including on vacuum boundaries.
    pub prescribed_neumann_out_w: f64,
    /// Reflected radiation and all actual inner work.
    pub radiation: EnthalpyEnclosureReport,
    /// Complete nonlinear residual at the accepted h field, J.
    pub physical_residual_norm_j: f64,
    /// Original physical initial-residual target, J.
    pub physical_residual_tolerance_j: f64,
    /// Storage minus net input using fresh enclosure exchange, J.
    pub physical_energy_residual_j: f64,
    /// Inward bath reaction using fresh enclosure exchange, W, distinct from
    /// the applied-boundary reaction in `conduction`.
    pub physical_dirichlet_in_w: f64,
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Advance closed enclosure exchange on existing natural boundary traces.
    /// Each named radiation surface must exactly match one positive uniform
    /// Robin or uniform Neumann region. Zero Neumann flux declares vacuum;
    /// nonzero authored flux remains additive and is never replaced by exchange.
    /// The existing mesh, source, contact, phase and nonlinear conductivity
    /// owners remain authoritative. At most 64 patches are admitted.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_with_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        enclosure: &GrayDiffuseEnclosure,
        radiation_config: AmbientRadiationConfig,
    ) -> Result<EnthalpyEnclosureStepSolution, EnthalpyError> {
        self.admit_step(cx, problem, old_h, dt_s, step_config)?;
        radiation_config.validate()?;
        radiation_config
            .max_iterations
            .checked_mul(step_config.max_newton_iterations)
            .and_then(|v| v.checked_mul(step_config.newton.linear_restart))
            .and_then(|v| v.checked_mul(step_config.newton.max_linear_cycles))
            .ok_or(EnthalpyError::InvalidInput(
                "enclosure Newton work overflow",
            ))?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_h,
            dt: dt_s,
            pointwise_radiation: None,
        };
        let old_temperature = context.temperatures(old_h)?;
        let patches = bind(cx, problem, enclosure, &old_temperature)?;
        let prescribed_neumann_out_w = prescribed_neumann_heat(cx, problem)?;
        let prepared = enclosure.prepare_radiosity(cx)?;
        let target = initial_target(&context, step_config, &prepared, &patches, radiation_config)?;
        let mut driving = means(cx, problem, enclosure, &old_temperature)?;
        let mut solid_iterations = 0_usize;
        let mut krylov_iterations = 0_usize;
        let mut last_change = 0.0_f64;
        let mut last_mismatch = 0.0_f64;
        for iteration in 0..radiation_config.max_iterations {
            poll(cx, iteration)?;
            let applied = prepared_radiosity(cx, &prepared, &driving, radiation_config)?;
            let boundary = shifted_boundary(problem, &patches, &applied)?;
            let trial_context = StepContext {
                problem: ConductionProblem {
                    boundary: &boundary,
                    ..problem
                },
                ..context
            };
            let mut initial = vec![0.0; old_h.len()];
            trial_context.residual(&trial_context.initial_values()?, &mut initial)?;
            let norm = finite(norm2(&initial))?;
            let mut inner = step_config;
            inner.newton.absolute_tolerance = target * 0.25;
            if target != 0.0 || norm != 0.0 {
                inner.newton.relative_tolerance = inner
                    .newton
                    .relative_tolerance
                    .min((target * 0.25) / norm.max(f64::MIN_POSITIVE));
            }
            if inner.newton.relative_tolerance <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "enclosure residual tolerance is not representable",
                ));
            }
            let solid = self.advance(cx, trial_context.problem, interfaces, old_h, dt_s, inner)?;
            solid_iterations = solid_iterations
                .checked_add(solid.newton.iterations)
                .ok_or(EnthalpyError::InvalidInput(
                    "enclosure Newton work overflow",
                ))?;
            for row in &solid.newton.history {
                krylov_iterations = krylov_iterations.checked_add(row.linear_iterations).ok_or(
                    EnthalpyError::InvalidInput("enclosure Krylov work overflow"),
                )?;
            }
            let actual = means(cx, problem, enclosure, &solid.temperature)?;
            let fresh = prepared_radiosity(cx, &prepared, &actual, radiation_config)?;
            let (convection, convective_out) = convection(cx, problem, &solid)?;
            let mut heat_ok = true;
            last_mismatch = 0.0;
            for (index, patch) in patches.iter().enumerate() {
                poll(cx, index)?;
                let actual_applied = match patch.law {
                    PatchLaw::Convection { .. } => {
                        let name = &problem.boundary.region_names()[patch.region];
                        let full = solid.robin_fluxes.iter().find(|row| row.region == *name)
                            .ok_or(EnthalpyError::InvalidInput(
                                "enclosure inner solve omitted a convective patch",
                            ))?;
                        let air = convection.iter().find(|row| row.region == *name)
                            .ok_or(EnthalpyError::InvalidInput(
                                "enclosure convection omitted a patch",
                            ))?;
                        finite(full.heat_rate_w - air.heat_rate_w)?
                    }
                    PatchLaw::PrescribedFlux { outward_flux } => {
                        let ThermalBc::Neumann { outward_flux: ScalarField::Uniform(combined) }
                            = boundary.conditions()[patch.region]
                        else {
                            return Err(EnthalpyError::InvalidInput(
                                "enclosure changed prescribed-flux boundary ownership",
                            ));
                        };
                        finite((combined - outward_flux) * enclosure.surfaces()[index].area_m2())?
                    }
                };
                if finite(actual_applied - applied.net_outward_heat_w[index])?.abs()
                    > heat_tolerance(
                        radiation_config,
                        actual_applied,
                        applied.net_outward_heat_w[index],
                    )?
                {
                    return Err(EnthalpyError::InvalidInput(
                        "enclosure boundary heat does not separate into authored load and radiation",
                    ));
                }
                let mismatch = finite(actual_applied - fresh.net_outward_heat_w[index])?.abs();
                last_mismatch = last_mismatch.max(mismatch);
                heat_ok &= mismatch
                    <= heat_tolerance(
                        radiation_config,
                        actual_applied,
                        fresh.net_outward_heat_w[index],
                    )?;
            }
            last_change = actual
                .iter()
                .zip(&driving)
                .try_fold(0.0_f64, |m, (&a, &b)| {
                    Ok::<_, EnthalpyError>(m.max(finite(a - b)?.abs()))
                })?;
            let physical_boundary = shifted_boundary(problem, &patches, &fresh)?;
            let physical = StepContext {
                problem: ConductionProblem {
                    boundary: &physical_boundary,
                    ..problem
                },
                ..context
            };
            let mut residual = vec![0.0; old_h.len()];
            physical.residual(&solid.specific_enthalpy_j_kg, &mut residual)?;
            let residual_norm = finite(norm2(&residual))?;
            let physical_dirichlet_in_w =
                physical.physical_prescribed_heat_input_w(&solid.specific_enthalpy_j_kg)?;
            let energy = finite(
                solid.stored_energy_change_j
                    + dt_s
                        * (prescribed_neumann_out_w + convective_out + fresh.enclosure_energy_closure_w
                            - solid.source_w - physical_dirichlet_in_w),
            )?;
            if heat_ok
                && last_change <= radiation_config.temperature_tolerance_k
                && residual_norm <= target
                && energy.abs() <= step_config.energy_tolerance_j
            {
                poll(cx, iteration)?;
                return Ok(EnthalpyEnclosureStepSolution {
                    conduction: solid,
                    combined_boundary: boundary,
                    convective_robin_fluxes: convection,
                    convective_out_w: convective_out,
                    prescribed_neumann_out_w,
                    radiation: EnthalpyEnclosureReport {
                        radiosity: fresh,
                        applied_heat_w: applied.net_outward_heat_w,
                        max_heat_mismatch_w: last_mismatch,
                        iterations: iteration + 1,
                        solid_iterations,
                        krylov_iterations,
                    },
                    physical_residual_norm_j: residual_norm,
                    physical_residual_tolerance_j: target,
                    physical_energy_residual_j: energy,
                    physical_dirichlet_in_w,
                });
            }
            for (old, new) in driving.iter_mut().zip(actual) {
                *old = finite(
                    (1.0 - radiation_config.relaxation) * *old + radiation_config.relaxation * new,
                )?;
            }
        }
        Err(ConductionError::Radiation {surface:"enthalpy-enclosure".into(),
            what:format!("{} solid solves exhausted; final temperature change {last_change} K and patch mismatch {last_mismatch} W; physical residual or energy may also remain unresolved",radiation_config.max_iterations),
            fix:"provide adequate explicit enclosure/Newton budgets or reduce the timestep".into()}.into())
    }
}

fn bind(cx: &Cx<'_>, problem: ConductionProblem<'_>, enclosure: &GrayDiffuseEnclosure,
    temperature: &[f64]) -> Result<Vec<Patch>, EnthalpyError> {
    Ok(exchange::bind(cx, problem, enclosure, temperature)?)
}

fn means(cx: &Cx<'_>, problem: ConductionProblem<'_>, enclosure: &GrayDiffuseEnclosure,
    temperature: &[f64]) -> Result<Vec<f64>, EnthalpyError> {
    Ok(exchange::means(cx, problem, enclosure, temperature)?)
}

fn heat_tolerance(config: AmbientRadiationConfig, a: f64, b: f64) -> Result<f64, EnthalpyError> {
    Ok(exchange::heat_tolerance(config, a, b)?)
}

fn radiosity(cx: &Cx<'_>, enclosure: &GrayDiffuseEnclosure, temperature: &[f64],
    config: AmbientRadiationConfig) -> Result<RadiosityReport, EnthalpyError> {
    Ok(exchange::radiosity(cx, enclosure, temperature, config)?)
}

fn prepared_radiosity(cx: &Cx<'_>, prepared: &PreparedRadiosity<'_>, temperature: &[f64],
    config: AmbientRadiationConfig) -> Result<RadiosityReport, EnthalpyError> {
    Ok(exchange::prepared_radiosity(cx, prepared, temperature, config)?)
}

fn shifted_boundary(problem: ConductionProblem<'_>, patches: &[Patch],
    radiosity: &RadiosityReport) -> Result<ThermalBoundary, EnthalpyError> {
    Ok(exchange::shifted_boundary(problem, patches, radiosity)?)
}

fn prescribed_neumann_heat(cx: &Cx<'_>, problem: ConductionProblem<'_>)
    -> Result<f64, EnthalpyError> {
    Ok(exchange::prescribed_neumann_heat(cx, problem)?)
}

fn convection(cx: &Cx<'_>, problem: ConductionProblem<'_>, solid: &EnthalpyStepSolution)
    -> Result<(Vec<RobinFlux>, f64), EnthalpyError> {
    Ok(exchange::convection(cx, problem, &solid.robin_fluxes)?)
}

fn initial_target(
    context: &StepContext<'_, '_, '_, '_>,
    config: EnthalpyStepConfig,
    prepared: &PreparedRadiosity<'_>,
    patches: &[Patch],
    radiation_config: AmbientRadiationConfig,
) -> Result<f64, EnthalpyError> {
    let initial_h = context.initial_values()?;
    let temperature = context.temperatures(&initial_h)?;
    let means = means(context.cx, context.problem, prepared.enclosure(), &temperature)?;
    let initial = prepared_radiosity(context.cx, prepared, &means, radiation_config)?;
    let boundary = shifted_boundary(context.problem, patches, &initial)?;
    let physical = StepContext {
        problem: ConductionProblem {
            boundary: &boundary,
            ..context.problem
        },
        storage: context.storage,
        cx: context.cx,
        interfaces: context.interfaces,
        old: context.old,
        dt: context.dt,
        pointwise_radiation: None,
    };
    let stage = physical.stage(&initial_h)?;
    let initial = NewtonKrylovState::new(&stage, initial_h, config.newton)
        .map_err(|error| stage.take_failure().unwrap_or(EnthalpyError::Newton(error)))?;
    finite(
        config
            .newton
            .absolute_tolerance
            .max(config.newton.relative_tolerance * initial.residual_norm().max(f64::MIN_POSITIVE)),
    )
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Re-evaluate a supplied enclosure endpoint under the original boundary laws.
    ///
    /// No solve is performed and neither candidate nor old specific enthalpy is
    /// modified. Fresh reflected radiation enters the same physical residual and
    /// prescribed-temperature reaction as the primal. Returned convection and
    /// Neumann heat retain their authored values, excluding internal radiation.
    /// The caller applies its original fixed residual and energy allowances.
    ///
    /// # Errors
    /// Refuses invalid storage/candidate/boundary data, failed radiosity gates,
    /// mismatched traces, nonfinite arithmetic or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_enclosure(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_h: &[f64], candidate_h: &[f64],
        dt_s: f64, config: EnthalpyStepConfig, enclosure: &GrayDiffuseEnclosure,
        radiation_config: AmbientRadiationConfig,
    ) -> Result<super::EnthalpyEndpointEvaluation, EnthalpyError> {
        self.admit_step(cx, problem, old_h, dt_s, config)?;
        radiation_config.validate()?;
        if candidate_h.len() != self.masses.len() {
            return Err(EnthalpyError::InvalidInput(
                "candidate enclosure enthalpy must match the reference mass mesh",
            ));
        }
        let context = StepContext {
            storage: self, cx, problem, interfaces, old: old_h, dt: dt_s,
            pointwise_radiation: None,
        };
        context.check_prescribed_values(candidate_h)?;
        let temperature = context.temperatures(candidate_h)?;
        let patches = bind(cx, problem, enclosure, &temperature)?;
        let fresh = radiosity(cx, enclosure,
            &means(cx, problem, enclosure, &temperature)?, radiation_config)?;
        let boundary = shifted_boundary(problem, &patches, &fresh)?;
        let mut physical = self.evaluate_endpoint(cx,
            ConductionProblem { boundary: &boundary, ..problem }, interfaces,
            old_h, candidate_h, dt_s, config, None)?;
        let (rows, convection_w) =
            exchange::convection(cx, problem, &physical.convective_robin_fluxes)?;
        physical.convective_robin_fluxes = rows;
        physical.convective_out_w = convection_w;
        physical.neumann_out_w = prescribed_neumann_heat(cx, problem)?;
        physical.radiation_out_w = fresh.enclosure_energy_closure_w;
        physical.energy_residual_j = finite(physical.stored_energy_change_j - dt_s * (
            physical.source_w + physical.dirichlet_in_w - physical.neumann_out_w
                - physical.convective_out_w - physical.radiation_out_w
        ))?;
        poll(cx, candidate_h.len())?;
        Ok(physical)
    }
}
