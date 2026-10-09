//! Spatial enthalpy stepping with implicit ambient gray radiation.
//!
//! The shared ambient-patch driver solves the same backward-Euler history on
//! every radiative trial. Its area-mean law is epsilon sigma A (T_mean^4-Ta^4),
//! applied through the existing pointwise Robin trace. This is not quadrature
//! of T(x)^4, enclosure exchange, or an enthalpy-radiation adjoint.
//!
//! Acceptance additionally reassembles the boundary at the returned endpoint,
//! checks the full nodal residual in joules, and closes stored energy against
//! nonlinear radiative heat. Small relaxed updates alone never accept a step.

use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_solver::{NewtonKrylovState, norm2};

use super::{
    EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution, StepContext,
    finite, poll,
};
use crate::{
    AmbientRadiationConfig, AmbientRadiationPatch, AmbientRadiationReport, ConductionProblem,
    RobinFlux, ScalarField, ThermalBc, ThermalBoundary, ThermalInterfaces,
    radiation::{AmbientRadiationTrial, solve_ambient_radiation_with},
};

/// One accepted physical time step; all radiative trials share immutable history.
#[derive(Debug, Clone)]
pub struct EnthalpyRadiationStepSolution {
    /// Actual enthalpy/phase step. Its Robin rows contain applied radiation and
    /// convection together, so only the separated convection rows feed air.
    pub conduction: EnthalpyStepSolution,
    /// Exact frozen-secant boundary used by the final inner enthalpy solve.
    pub combined_boundary: ThermalBoundary,
    /// Original convection only, in boundary declaration order.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Sum of original convection heat, watts outward.
    pub convective_out_w: f64,
    /// Nonlinear radiation and cumulative work across all radiative trials.
    pub radiation: AmbientRadiationReport,
    /// Euclidean nodal balance at the returned field with its own radiative
    /// coefficients, joules; independent of the frozen inner Newton report.
    pub physical_residual_norm_j: f64,
    /// Caller Newton target from the complete radiative residual at history.
    pub physical_residual_tolerance_j: f64,
    /// Stored energy minus dt times net input using nonlinear radiative heat,
    /// joules; must also satisfy the caller's original absolute energy budget.
    pub physical_energy_residual_j: f64,
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Advance natural-boundary enthalpy transport with implicit ambient
    /// radiation, preserving heterogeneous conductivity and finite contact.
    ///
    /// Patches bind distinct uniform Robin regions. Initial and actual surface
    /// area-mean temperatures must lie in each emissivity claim's support. Each inner
    /// Newton solve starts from the SAME immutable old enthalpy; radiative
    /// iterations do not advance time or accumulate latent heat repeatedly.
    ///
    /// The coupled residual target is `max(atol, rtol*||R(h_old)||)` in joules
    /// with radiation evaluated at history. Inner residual budgets tighten to
    /// one quarter of that target. The public watt and raw-temperature gates
    /// remain independent, and the accepted nonlinear heat must close the
    /// original energy budget. Failure returns no endpoint.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        patches: &[AmbientRadiationPatch],
        radiation_config: AmbientRadiationConfig,
    ) -> Result<EnthalpyRadiationStepSolution, EnthalpyError> {
        self.admit_step(cx, problem, old_h, dt_s, step_config)?;
        radiation_config.validate()?;
        radiation_config
            .max_iterations
            .checked_mul(step_config.max_newton_iterations)
            .and_then(|v| v.checked_mul(step_config.newton.linear_restart))
            .and_then(|v| v.checked_mul(step_config.newton.max_linear_cycles))
            .ok_or(EnthalpyError::InvalidInput(
                "radiative Newton work budget overflow",
            ))?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_h,
            dt: dt_s,
        };
        let old_temperature = context.temperatures(old_h)?;
        let initial_boundary = boundary_at_temperature(cx, problem, patches, &old_temperature)?;
        let initial_context = StepContext {
            problem: ConductionProblem {
                boundary: &initial_boundary,
                ..problem
            },
            ..context
        };
        // The shared constructor validates all original Newton controls and
        // measures the complete physical initial residual without advancing it.
        let stage = initial_context.stage(old_h)?;
        let initial = NewtonKrylovState::new(&stage, old_h.to_vec(), step_config.newton)
            .map_err(|error| stage.take_failure().unwrap_or(EnthalpyError::Newton(error)))?;
        let target = finite(step_config.newton.absolute_tolerance.max(
            step_config.newton.relative_tolerance * initial.residual_norm().max(f64::MIN_POSITIVE),
        ))?;
        drop(initial);
        drop(stage);
        let result = solve_ambient_radiation_with(
            cx,
            problem,
            patches,
            radiation_config,
            |combined_problem| {
                let context = StepContext {
                    storage: self,
                    cx,
                    problem: combined_problem,
                    interfaces,
                    old: old_h,
                    dt: dt_s,
                };
                let mut initial_residual = vec![0.0; old_h.len()];
                context.residual(old_h, &mut initial_residual)?;
                let initial_norm = finite(norm2(&initial_residual))?;
                let mut inner_config = step_config;
                inner_config.newton.absolute_tolerance = target * 0.25;
                if target != 0.0 || initial_norm != 0.0 {
                    inner_config.newton.relative_tolerance = step_config
                        .newton
                        .relative_tolerance
                        .min((target * 0.25) / initial_norm.max(f64::MIN_POSITIVE));
                }
                if inner_config.newton.relative_tolerance <= 0.0 {
                    return Err(EnthalpyError::InvalidInput(
                        "radiative residual tolerance is not representable",
                    ));
                }
                let conduction =
                    self.advance(cx, combined_problem, interfaces, old_h, dt_s, inner_config)?;
                let krylov_iterations = conduction
                    .newton
                    .history
                    .iter()
                    .try_fold(0_usize, |sum, row| sum.checked_add(row.linear_iterations))
                    .ok_or(EnthalpyError::InvalidInput(
                        "radiative Krylov work overflow",
                    ))?;
                Ok(AmbientRadiationTrial {
                    robin_fluxes: conduction.robin_fluxes.clone(),
                    robin_out_w: conduction.robin_out_w,
                    solid_iterations: conduction.newton.iterations,
                    krylov_iterations,
                    solution: (conduction, 0.0, 0.0),
                })
            },
            |trial, rows, convective_out_w| {
                let physical_boundary =
                    boundary_at_temperature(cx, problem, patches, &trial.0.temperature)?;
                let physical = StepContext {
                    storage: self,
                    cx,
                    problem: ConductionProblem {
                        boundary: &physical_boundary,
                        ..problem
                    },
                    interfaces,
                    old: old_h,
                    dt: dt_s,
                };
                let mut residual = vec![0.0; old_h.len()];
                physical.residual(&trial.0.specific_enthalpy_j_kg, &mut residual)?;
                trial.1 = finite(norm2(&residual))?;
                let radiation_out_w = rows
                    .iter()
                    .try_fold(0.0, |sum, row| finite(sum + row.nonlinear_heat_w))?;
                let net_out_w = finite(
                    trial.0.neumann_out_w + convective_out_w + radiation_out_w - trial.0.source_w,
                )?;
                trial.2 = finite(trial.0.stored_energy_change_j + dt_s * net_out_w)?;
                poll(cx, old_h.len())?;
                Ok(trial.1 <= target && trial.2.abs() <= step_config.energy_tolerance_j)
            },
        )?;
        Ok(EnthalpyRadiationStepSolution {
            conduction: result.solution.0,
            combined_boundary: result.combined_boundary,
            convective_robin_fluxes: result.convective_robin_fluxes,
            convective_out_w: result.convective_out_w,
            radiation: result.radiation,
            physical_residual_norm_j: result.solution.1,
            physical_residual_tolerance_j: target,
            physical_energy_residual_j: result.solution.2,
        })
    }
}

/// Evaluate the nonlinear boundary on a complete physical field. The actual
/// integration and outward sign stay owned by the existing Robin assembler.
fn boundary_at_temperature(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    patches: &[AmbientRadiationPatch],
    temperature: &[f64],
) -> Result<ThermalBoundary, EnthalpyError> {
    if patches.is_empty() || patches.len() > problem.boundary.region_names().len() {
        return Err(EnthalpyError::InvalidInput(
            "nonempty distinct radiation regions required",
        ));
    }
    let mut seen = BTreeSet::new();
    let mut replacements = Vec::with_capacity(patches.len());
    for patch in patches {
        poll(cx, replacements.len())?;
        if !seen.insert(patch.region()) {
            return Err(EnthalpyError::InvalidInput(
                "duplicate enthalpy radiation region",
            ));
        }
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == patch.region())
            .ok_or(EnthalpyError::InvalidInput(
                "unknown enthalpy radiation region",
            ))?;
        let ThermalBc::Robin {
            htc: ScalarField::Uniform(h),
            t_ref: ScalarField::Uniform(reference),
        } = problem.boundary.conditions()[region]
        else {
            return Err(EnthalpyError::InvalidInput(
                "radiation requires a uniform Robin region",
            ));
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
            return Err(EnthalpyError::InvalidInput(
                "radiation region has no exterior area",
            ));
        }
        let h_rad = patch.secant_coefficient_w_m2_k(finite(integral / area)?)?;
        let combined = finite(h + h_rad)?;
        let reference = finite(
            (h / combined) * reference + (h_rad / combined) * patch.ambient_temperature_k(),
        )?;
        replacements.push((region, combined, reference));
    }
    Ok(problem
        .boundary
        .with_uniform_robin_replacements(&replacements)?)
}
