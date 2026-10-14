//! Temperature-dependent contact with primary enthalpy and immutable history.

use fs_exec::Cx;
use crate::ConductionProblem;
use crate::transient::enthalpy::{
    EnthalpyEndpointEvaluation, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
};
use crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyBackwardEuler;
use super::{ContactIterationConfig, TemperatureContactEvaluation, TemperatureContacts, invalid, poll};

/// One accepted enthalpy field and the original and re-evaluated contact laws.
#[derive(Debug)]
pub struct ContactEnthalpyStep {
    /// Last counted nonlinear response. Advance history with its specific h,
    /// never by attempting to invert its possibly isothermal temperature field.
    pub conduction: EnthalpyStepSolution,
    /// Source states and operator actually used by the final inner response.
    pub applied_contacts: TemperatureContactEvaluation,
    /// Fresh source states at the unchanged, chart-resolved endpoint temperatures.
    pub physical_contacts: TemperatureContactEvaluation,
    /// Full nodal enthalpy residual and independent physical energy audit.
    pub physical: EnthalpyEndpointEvaluation,
    /// All contact trials, including discarded nonlinear responses.
    pub iterations: usize,
    /// Sum of actual Newton iterations over every contact trial.
    pub nonlinear_iterations: usize,
    /// Sum of actual Krylov iterations over every contact trial.
    pub krylov_iterations: usize,
}

impl TemperatureContacts<'_, '_> {
    /// Advance heterogeneous enthalpy with source-backed R(T_interface).
    ///
    /// Every contact trial calls the existing spatial enthalpy solver from the
    /// SAME old h. That owner retains heterogeneous reference masses, distinct
    /// charts, nonlinear conductivity and optional phase-conductivity feedback.
    /// Only the contact-query temperature guess is relaxed. Neither h nor the
    /// final field is mixed or reconstructed from temperature.
    ///
    /// Requery every face at the candidate's actual temperature and reassemble
    /// the complete nodal joule residual before acceptance. A latent plateau or
    /// a closed total-energy sum cannot conceal a local enthalpy imbalance.
    /// The original contact residual target and energy tolerance never widen.
    /// Inner Newton tolerances must be tight enough to satisfy those gates.
    ///
    /// Prescribed Neumann/Robin data and fixed contact geometry are supported.
    /// Dirichlet enthalpy, radiation/air feedback and nonlinear-contact adjoints
    /// are not supplied by this entry point. No state is returned on failure.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_enthalpy(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        problem: ConductionProblem<'_>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        coupling: ContactIterationConfig,
    ) -> Result<ContactEnthalpyStep, EnthalpyError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        step_config.newton.validate().map_err(EnthalpyError::Newton)?;
        let per_response = step_config.max_newton_iterations
            .checked_mul(step_config.newton.linear_restart)
            .and_then(|n| n.checked_mul(step_config.newton.max_linear_cycles))
            .ok_or(EnthalpyError::InvalidInput("contact enthalpy inner work allowance overflows"))?;
        let max_krylov = coupling.admit(per_response)?;
        let max_newton = coupling.max_iterations.checked_mul(step_config.max_newton_iterations)
            .ok_or(EnthalpyError::InvalidInput("contact enthalpy Newton allowance overflows"))?;
        if old_h.len() != problem.mesh.vertex_count()
            || engine.vertex_material_ids().len() != old_h.len()
        {
            return Err(EnthalpyError::InvalidInput("contact enthalpy history and reference mesh differ"));
        }
        let mut guess = Vec::new();
        guess.try_reserve_exact(old_h.len())
            .map_err(|_| EnthalpyError::InvalidInput("contact enthalpy query field allocation refused"))?;
        for (vertex, &h) in old_h.iter().enumerate() {
            if vertex % 256 == 0 { poll(cx, vertex)?; }
            let curve = engine.phase_curve_at(vertex)
                .ok_or(EnthalpyError::InvalidInput("contact enthalpy vertex has no phase chart"))?;
            let state = curve.state_at_specific_enthalpy(h)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            guess.push(state.temperature_k());
        }
        let mut newton = 0usize;
        let mut krylov = 0usize;
        let mut last_residual = f64::INFINITY;
        for iteration in 1..=coupling.max_iterations {
            poll(cx, iteration)?;
            let applied = self.evaluate(cx, problem.boundary, &guess)?;
            let solved = engine.advance(
                cx, problem, Some(applied.frozen_interfaces()), old_h, dt_s, step_config,
            )?;
            newton = newton.checked_add(solved.newton.iterations)
                .filter(|&n| n <= max_newton)
                .ok_or(EnthalpyError::InvalidInput("contact enthalpy Newton work exceeds allowance"))?;
            for row in &solved.newton.history {
                poll(cx, krylov)?;
                krylov = krylov.checked_add(row.linear_iterations)
                    .filter(|&n| n <= max_krylov)
                    .ok_or(EnthalpyError::InvalidInput("contact enthalpy Krylov work exceeds allowance"))?;
            }
            let actual = self.evaluate(cx, problem.boundary, &solved.temperature)?;
            let physical = engine.evaluate_endpoint(
                cx, problem, Some(actual.frozen_interfaces()), old_h,
                &solved.specific_enthalpy_j_kg, dt_s, step_config, None,
            )?;
            poll(cx, iteration)?;
            last_residual = physical.residual_norm_j;
            if last_residual.is_finite() && last_residual <= coupling.residual_tolerance_j
                && physical.energy_residual_j.is_finite()
                && physical.energy_residual_j.abs() <= step_config.energy_tolerance_j
            {
                return Ok(ContactEnthalpyStep {
                    conduction: solved, applied_contacts: applied, physical_contacts: actual,
                    physical, iterations: iteration, nonlinear_iterations: newton,
                    krylov_iterations: krylov,
                });
            }
            coupling.relax(cx, &mut guess, &solved.temperature)?;
        }
        Err(invalid("<contact-enthalpy-step>", format!(
            "temperature-dependent contact exhausted {} enthalpy responses; last physical residual {last_residual:e} J, target {:e} J",
            coupling.max_iterations, coupling.residual_tolerance_j,
        )).into())
    }
}
