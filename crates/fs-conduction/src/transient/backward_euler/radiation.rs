//! Immutable-history temperature-capacity steps with ambient gray radiation.
//!
//! The existing area-mean patch driver supplies frozen Robin secants. Every
//! trial solves the SAME physical history, including nonlinear conductivity
//! and prescribed endpoint temperatures. Acceptance reassembles the secants
//! at the returned field and checks both the free residual in joules and the
//! storage balance with actual nonlinear radiative heat.

use std::collections::BTreeSet;

use super::{
    BackwardEuler, NonlinearStepConfig, NonlinearStepSolution, StepConfig,
    StepSolution, energy_balance, finite, invalid, poll, sum,
};
use crate::{
    AmbientRadiationConfig, AmbientRadiationPatch, AmbientRadiationReport,
    ConductionError, ConductionProblem, RobinFlux, ScalarField, ThermalBc,
    ThermalBoundary, ThermalInterfaces,
    radiation::{AmbientRadiationTrial, solve_ambient_radiation_with},
};
use fs_exec::Cx;

/// One accepted physical step with the original convection separated from radiation.
#[derive(Debug, Clone)]
pub struct RadiationStepSolution {
    /// Final inner endpoint with its actual frozen-secant Robin heat/reaction.
    /// This field never substitutes nonlinear radiation into an inner report.
    pub conduction: StepSolution,
    /// Final accepted inner k(T) solve, when that explicit policy was supplied.
    /// Its work counts cover this inner solve; cumulative counts are below.
    pub nonlinear: Option<NonlinearStepSolution>,
    /// Rejected Newton trials across ALL successful radiative inner solves.
    pub nonlinear_backtracks: usize,
    /// Exact frozen-secant boundary used by the final inner solve.
    pub combined_boundary: ThermalBoundary,
    /// Original convection only, in boundary declaration order.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Sum of original convective heat, watts outward.
    pub convective_out_w: f64,
    /// Patch physics and all inner Newton/Krylov work across the outer solve.
    pub radiation: AmbientRadiationReport,
    /// Free endpoint balance reassembled at its own radiative coefficients, J.
    pub physical_residual_norm_j: f64,
    /// Fixed target derived from the complete physical predictor residual, J.
    pub physical_residual_tolerance_j: f64,
    /// Stored energy minus step time times actual net input, J.
    pub physical_energy_residual_j: f64,
    /// Reassembled prescribed reaction including fixed-node capacity, W inward.
    /// It can differ from the frozen inner reaction on shared boundary vertices.
    pub physical_dirichlet_in_w: f64,
}

struct Trial {
    conduction: StepSolution,
    nonlinear: Option<NonlinearStepSolution>,
    physical_residual_norm_j: f64,
    physical_energy_residual_j: f64,
    physical_dirichlet_in_w: f64,
}

impl BackwardEuler<'_> {
    /// Advance fixed-capacity conduction with implicit ambient gray radiation.
    ///
    /// The full old field remains immutable physical history on every outer
    /// trial. Endpoint Dirichlet values may differ from history, exactly as in
    /// advance_prescribed. Supply an explicit nonlinear policy for k(T);
    /// otherwise the ordinary constant-conductivity refusal remains in force.
    ///
    /// Patches bind distinct uniform Robin regions and retain their sourced
    /// emissivity validity. This is the existing area-mean T^4 patch model,
    /// applied through the pointwise Robin trace, with a black ambient reservoir.
    ///
    /// The nonlinear physical target is atol + rtol * ||R_predictor||. The
    /// linear target is linear_tolerance * ||R_predictor|| plus 1% of the caller's
    /// energy budget divided by sqrt(free DOFs), a rounding floor for the
    /// published absolute temperatures. The predictor has old free temperatures
    /// and new prescribed temperatures, while its capacity history is unchanged.
    /// Inner policies only tighten. Raw patch-temperature and watt convergence,
    /// the complete free residual, and the original energy gate must all pass.
    ///
    /// # Errors
    ///
    /// Invalid input, patch/material validity, exhausted work, cancellation or
    /// incomplete physical balance return no endpoint. Radiative trials do not
    /// advance time, widen an energy budget, or add latent storage.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_prescribed_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old: &[f64],
        dt_s: f64,
        step_config: StepConfig,
        nonlinear: Option<NonlinearStepConfig>,
        patches: &[AmbientRadiationPatch],
        radiation_config: AmbientRadiationConfig,
    ) -> Result<RadiationStepSolution, ConductionError> {
        let dofs = self.admit_endpoint_step(cx, problem, old, dt_s, step_config)?;
        radiation_config.validate()?;
        radiation_config
            .max_iterations
            .checked_mul(step_config.linear.max_iterations)
            .ok_or_else(|| invalid("radiative Krylov work budget overflow"))?;
        if let Some(policy) = nonlinear {
            policy.validate()?;
            let trials = policy
                .line_search
                .max_backtracks
                .checked_add(1)
                .and_then(|v| v.checked_mul(policy.max_iterations))
                .and_then(|v| v.checked_mul(radiation_config.max_iterations));
            if trials.is_none() {
                return Err(invalid("radiative Newton trial budget overflow"));
            }
        }
        let mut predictor = old.to_vec();
        for (slot, &vertex) in dofs.fixed().iter().enumerate() {
            if slot % 512 == 0 {
                poll(cx, slot)?;
            }
            predictor[vertex] = dofs.prescribed()[vertex];
        }
        let initial_norm_j = {
            let initial_boundary = boundary_at_temperature(cx, problem, patches, &predictor)?;
            self.evaluate_endpoint(
                cx,
                ConductionProblem {
                    boundary: &initial_boundary,
                    ..problem
                },
                interfaces,
                old,
                &predictor,
                dt_s,
                &dofs,
            )?.norm
        };
        let target = match nonlinear {
            Some(policy) => finite(
                policy.residual_atol_j + policy.residual_rtol * initial_norm_j,
            )?,
            None => finite(
                step_config.linear.tolerance * initial_norm_j
                    + (step_config.energy_tolerance_j * 0.01) / (dofs.n() as f64).sqrt(),
            )?,
        };
        let mut nonlinear_backtracks = 0_usize;
        let result = solve_ambient_radiation_with::<_, ConductionError>(
            cx,
            problem,
            patches,
            radiation_config,
            |combined_problem| {
                let inner_initial_norm_j = self.evaluate_endpoint(
                    cx, combined_problem, interfaces, old, &predictor, dt_s, &dofs,
                )?.norm;
                let (conduction, nonlinear) = match nonlinear {
                    Some(mut policy) => {
                        // Two additive Newton tolerances together use at most
                        // one quarter of the fixed physical endpoint target.
                        policy.residual_atol_j = policy.residual_atol_j.min(target * 0.125);
                        policy.residual_rtol = tighter_relative(
                            policy.residual_rtol, target * 0.125, inner_initial_norm_j,
                        )?;
                        let solved = self.advance_nonlinear_prescribed(
                            cx, combined_problem, interfaces, old, dt_s, step_config, policy,
                        )?;
                        nonlinear_backtracks = nonlinear_backtracks
                            .checked_add(solved.backtracks)
                            .ok_or_else(|| invalid("radiative backtrack count overflow"))?;
                        (solved.step.clone(), Some(solved))
                    }
                    None => {
                        let mut inner_config = step_config;
                        inner_config.linear.tolerance = tighter_relative(
                            step_config.linear.tolerance, target * 0.25, inner_initial_norm_j,
                        )?;
                        (
                            self.advance_prescribed(
                                cx, combined_problem, interfaces, old, dt_s, inner_config,
                            )?,
                            None,
                        )
                    }
                };
                Ok(AmbientRadiationTrial {
                    robin_fluxes: conduction.robin_fluxes.clone(),
                    robin_out_w: conduction.robin_out_w,
                    solid_iterations: nonlinear
                        .as_ref()
                        .map_or(0, |solution| solution.nonlinear_iterations),
                    krylov_iterations: conduction.krylov_iterations,
                    solution: Trial {
                        conduction,
                        nonlinear,
                        physical_residual_norm_j: 0.0,
                        physical_energy_residual_j: 0.0,
                        physical_dirichlet_in_w: 0.0,
                    },
                })
            },
            |trial, rows, convective_out_w| {
                let physical_boundary = boundary_at_temperature(
                    cx, problem, patches, &trial.conduction.temperature,
                )?;
                let physical_problem = ConductionProblem {
                    boundary: &physical_boundary,
                    ..problem
                };
                let physical = self.evaluate_endpoint(
                    cx, physical_problem, interfaces, old,
                    &trial.conduction.temperature, dt_s, &dofs,
                )?;
                trial.physical_residual_norm_j = physical.norm;
                let (energy, _) = energy_balance(
                    self.mesh, &physical_boundary, problem.source, &physical.system,
                    &dofs, &trial.conduction.temperature,
                );
                let delta = trial
                    .conduction
                    .temperature
                    .iter()
                    .zip(old)
                    .map(|(t, previous)| finite(t - previous))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut storage = vec![0.0; old.len()];
                self.capacity.spmv(&delta, &mut storage);
                trial.physical_dirichlet_in_w = finite(
                    energy.dirichlet_in_w
                        + sum(dofs.fixed().iter().map(|&v| storage[v]))? / dt_s,
                )?;
                let radiation_out_w = sum(rows.iter().map(|row| row.nonlinear_heat_w))?;
                let net_input_w = finite(
                    energy.source_w + trial.physical_dirichlet_in_w
                        - energy.neumann_out_w - convective_out_w - radiation_out_w,
                )?;
                trial.physical_energy_residual_j = finite(
                    trial.conduction.stored_energy_change_j - dt_s * net_input_w,
                )?;
                poll(cx, old.len())?;
                Ok(trial.physical_residual_norm_j <= target
                    && trial.physical_energy_residual_j.abs() <= step_config.energy_tolerance_j)
            },
        )?;
        Ok(RadiationStepSolution {
            conduction: result.solution.conduction,
            nonlinear: result.solution.nonlinear,
            nonlinear_backtracks,
            combined_boundary: result.combined_boundary,
            convective_robin_fluxes: result.convective_robin_fluxes,
            convective_out_w: result.convective_out_w,
            radiation: result.radiation,
            physical_residual_norm_j: result.solution.physical_residual_norm_j,
            physical_residual_tolerance_j: target,
            physical_energy_residual_j: result.solution.physical_energy_residual_j,
            physical_dirichlet_in_w: result.solution.physical_dirichlet_in_w,
        })
    }
}

fn tighter_relative(
    requested: f64,
    absolute_allowance_j: f64,
    initial_norm_j: f64,
) -> Result<f64, ConductionError> {
    let tolerance = if initial_norm_j == 0.0 {
        requested
    } else {
        requested.min(absolute_allowance_j / initial_norm_j)
    };
    if tolerance > 0.0 && tolerance.is_finite() {
        Ok(tolerance)
    } else {
        Err(invalid("radiative residual tolerance is not representable"))
    }
}

/// Evaluate the same mean-temperature law on a complete physical endpoint.
/// The existing Robin assembler still owns its pointwise trace and outward sign.
pub(super) fn boundary_at_temperature(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    patches: &[AmbientRadiationPatch],
    temperature: &[f64],
) -> Result<ThermalBoundary, ConductionError> {
    if patches.is_empty() || patches.len() > problem.boundary.region_names().len() {
        return Err(invalid("nonempty distinct radiation regions required"));
    }
    let mut seen = BTreeSet::new();
    let mut replacements = Vec::with_capacity(patches.len());
    for patch in patches {
        poll(cx, replacements.len())?;
        if !seen.insert(patch.region()) {
            return Err(invalid("duplicate transient radiation region"));
        }
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == patch.region())
            .ok_or_else(|| invalid("unknown transient radiation region"))?;
        let ThermalBc::Robin {
            htc: ScalarField::Uniform(h),
            t_ref: ScalarField::Uniform(reference),
        } = problem.boundary.conditions()[region]
        else {
            return Err(invalid("radiation requires a uniform Robin region"));
        };
        let mut area = 0.0;
        let mut integral = 0.0;
        for (slot, face) in problem.mesh.boundary().iter().enumerate() {
            if slot % crate::assemble::ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            if problem.boundary.region_for(slot) == Some(region) {
                area = finite(area + face.area)?;
                let mean = face
                    .vertices
                    .iter()
                    .map(|&v| temperature[v as usize] / 3.0)
                    .sum::<f64>();
                integral = finite(integral + face.area * mean)?;
            }
        }
        if area <= 0.0 {
            return Err(invalid("radiation region has no exterior area"));
        }
        let h_rad = patch.secant_coefficient_w_m2_k(finite(integral / area)?)?;
        let combined = finite(h + h_rad)?;
        let reference = finite(
            (h / combined) * reference + (h_rad / combined) * patch.ambient_temperature_k(),
        )?;
        replacements.push((region, combined, reference));
    }
    Ok(problem.boundary.with_uniform_robin_replacements(&replacements)?)
}
