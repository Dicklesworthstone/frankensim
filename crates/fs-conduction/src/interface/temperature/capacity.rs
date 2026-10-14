//! Backward Euler with freshly queried R(T_interface), not a frozen contact law.

use fs_exec::Cx;
use crate::{ConductionError, ConductionProblem};
use crate::transient::backward_euler::{
    BackwardEuler, EndpointEvaluation, NonlinearStepConfig, StepConfig, StepSolution,
};
use super::{ContactIterationConfig, TemperatureContactEvaluation, TemperatureContacts, invalid, poll};

/// One counted field, its applied contact operator, and a fresh physical audit.
#[derive(Debug)]
pub struct ContactCapacityStep {
    /// Last counted inner solve; its contact fluxes describe `applied_contacts`.
    pub conduction: StepSolution,
    /// Contact source states actually used by that inner thermal response.
    pub applied_contacts: TemperatureContactEvaluation,
    /// Source states queried at the unchanged returned temperature field.
    pub physical_contacts: TemperatureContactEvaluation,
    /// Full nodal residual and energy under the actual endpoint contact law.
    pub physical: EndpointEvaluation,
    /// All thermal responses, including discarded contact trials.
    pub iterations: usize,
    /// Cumulative inner Krylov work, not just the last response's work.
    pub krylov_iterations: usize,
    /// Cumulative Newton updates across all contact trials.
    pub nonlinear_updates: usize,
    /// Cumulative line-search backtracks across all contact trials.
    pub nonlinear_backtracks: usize,
}

impl TemperatureContacts<'_, '_> {
    /// Advance source-backed temperature-dependent matching contacts and k(T).
    ///
    /// Every inner linear/Newton solve receives the same immutable `old` field.
    /// Contact is evaluated at each face pair's two-trace mean. The candidate
    /// passes only when its COMPLETE nodal balance, reassembled with freshly
    /// queried endpoint resistances, meets the original absolute joule target
    /// and the existing independent energy gate. Opposing local heat errors
    /// cannot cancel into a false success. No final replacement solve is hidden.
    ///
    /// Prescribed Neumann/Robin and endpoint Dirichlet data are supported.
    /// Supply `nonlinear` for temperature-dependent bulk conductivity. This
    /// entry point supplies no radiation/airflow feedback or contact-law adjoint.
    /// Inner tolerances must be tight enough to meet the physical outer target;
    /// failure never relaxes it, extrapolates source data or publishes history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_capacity(
        &self,
        cx: &Cx<'_>,
        engine: &BackwardEuler<'_>,
        problem: ConductionProblem<'_>,
        old: &[f64],
        dt_s: f64,
        step_config: StepConfig,
        nonlinear: Option<NonlinearStepConfig>,
        coupling: ContactIterationConfig,
    ) -> Result<ContactCapacityStep, ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        let max_krylov = coupling.admit(step_config.linear.max_iterations)?;
        let mut guess = old.to_vec();
        let mut krylov = 0usize;
        let mut updates = 0usize;
        let mut backtracks = 0usize;
        let mut last_residual = f64::INFINITY;
        for iteration in 1..=coupling.max_iterations {
            poll(cx, iteration)?;
            let applied = self.evaluate(cx, problem.boundary, &guess)?;
            let solved = match nonlinear {
                Some(policy) => {
                    let response = engine.advance_nonlinear_prescribed(
                        cx, problem, Some(applied.frozen_interfaces()), old, dt_s, step_config, policy,
                    )?;
                    updates = updates.checked_add(response.nonlinear_iterations)
                        .ok_or_else(|| invalid("<contact-step>", "contact Newton counter overflows"))?;
                    backtracks = backtracks.checked_add(response.backtracks)
                        .ok_or_else(|| invalid("<contact-step>", "contact backtrack counter overflows"))?;
                    response.step
                },
                None => engine.advance_prescribed(
                    cx, problem, Some(applied.frozen_interfaces()), old, dt_s, step_config,
                )?,
            };
            krylov = krylov.checked_add(solved.krylov_iterations)
                .filter(|&n| n <= max_krylov)
                .ok_or_else(|| invalid("<contact-step>", "contact thermal work exceeds its combined allowance"))?;
            let actual = self.evaluate(cx, problem.boundary, &solved.temperature)?;
            let physical = engine.evaluate_prescribed_endpoint(
                cx, problem, Some(actual.frozen_interfaces()), old, &solved.temperature,
                dt_s, step_config, None,
            )?;
            poll(cx, iteration)?;
            last_residual = physical.residual_norm_j;
            if last_residual.is_finite() && last_residual <= coupling.residual_tolerance_j
                && physical.energy_residual_j.is_finite()
                && physical.energy_residual_j.abs() <= step_config.energy_tolerance_j
            {
                return Ok(ContactCapacityStep {
                    conduction: solved, applied_contacts: applied, physical_contacts: actual,
                    physical, iterations: iteration, krylov_iterations: krylov,
                    nonlinear_updates: updates, nonlinear_backtracks: backtracks,
                });
            }
            coupling.relax(cx, &mut guess, &solved.temperature)?;
        }
        Err(invalid("<contact-step>", format!(
            "temperature-dependent contact exhausted {} thermal responses; last physical residual {last_residual:e} J, target {:e} J",
            coupling.max_iterations, coupling.residual_tolerance_j,
        )))
    }
}
