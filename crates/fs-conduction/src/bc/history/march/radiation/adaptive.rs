//! Adaptive intervals retaining the original radiative material and heat reports.

use super::{RadiativeThermalHistoryRun, ThermalHistoryMethod, ConductionError,
    ConductionProblem, RadiationStepSolution, StepConfig, patches_before};
use super::super::events::{self, adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep,
    admit_proposal, select_step, temperature_error}};
use fs_exec::Cx;

impl RadiativeThermalHistoryRun<'_, '_> {
    /// Compare one full radiative endpoint with two half endpoints, never omitting
    /// the retained surfaces or changing their sourced emissivities. Every source,
    /// boundary and reservoir knot caps the proposed interval. Both accepted half
    /// reports retain separate convection/radiation and physical residual/energy
    /// checks; their stored-energy changes and integrated supplies may be summed.
    ///
    /// The pair consumes two physical-step slots and publishes only its second
    /// endpoint after all checks. Local temperature error is estimated, not bounded
    /// globally. The original radiation/Newton policies and refusal behavior remain.
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        step: StepConfig, adaptive: AdaptiveHistoryConfig)
        -> Result<AdaptiveHistoryStep<RadiationStepSolution>, ConductionError>
    {
        let start = self.run.time_s;
        let cap = admit_proposal(cx, self.run.accepted_steps, self.run.max_steps, proposed_step_s, adaptive)?;
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.patches.iter().map(|p| &p.ambient_temperature), start, target_s, cap)?;
        let nonlinear = match self.run.method {
            ThermalHistoryMethod::Linear => None,
            ThermalHistoryMethod::Nonlinear(policy) => Some(policy),
        };
        let solve = |old: &[f64], from: f64, to: f64| -> Result<RadiationStepSolution, ConductionError> {
            let boundary = events::boundary_before(cx, self.run.problem.boundary, to)?;
            let source = events::source_before(&self.run.problem.source, to, old.len())?;
            let patches = patches_before(cx, &self.patches, to)?;
            self.run.stepper.advance_prescribed_with_ambient_radiation(cx, ConductionProblem {
                mesh: self.run.problem.boundary.mesh(), boundary: &boundary, source: &source,
                material: self.run.problem.material, element_materials: self.run.problem.element_materials },
                self.run.problem.interfaces, old, to - from, step, nonlinear, &patches, self.config)
        };
        let selected = select_step(cx, start, end, adaptive, |e| e, |middle, end| {
            let coarse = solve(&self.run.temperature, start, end)?;
            let first = solve(&self.run.temperature, start, middle)?;
            let second = solve(&first.conduction.temperature, middle, end)?;
            let error = temperature_error(cx, &self.run.temperature, &coarse.conduction.temperature,
                &second.conduction.temperature, adaptive)?;
            Ok((first, second, error))
        })?;
        self.run.publish(cx, selected.end_time_s, &selected.second_half.conduction)?;
        self.run.accepted_steps += 1;
        Ok(selected)
    }
}
