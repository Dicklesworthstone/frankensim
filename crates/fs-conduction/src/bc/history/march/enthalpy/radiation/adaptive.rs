//! Latent-heat-aware adaptation of the retained radiative enthalpy trajectory.

use super::{RadiativeEnthalpyHistoryRun, HeterogeneousEnthalpyError, EnthalpyStepConfig,
    EnthalpyRadiationStepSolution, ConductionProblem, conduction, poll};
use super::super::adaptive::{PhaseErrorTolerance, phase_error};
use super::super::super::{events::{self, adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep,
    admit_proposal, select_step}}, radiation::patches_before};
use fs_exec::Cx;

impl RadiativeEnthalpyHistoryRun<'_, '_, '_> {
    /// Adapt against local temperature, specific-enthalpy and liquid-fraction
    /// discrepancies while retaining every radiative surface and reservoir knot.
    /// Temperature alone cannot admit a latent plateau: all caloric gates pass.
    ///
    /// Each radiative/Newton candidate uses its own immutable input history. The
    /// accepted second half starts from the first half's actual enthalpy, not its
    /// temperature and not an extrapolated phase. Both complete physical reports
    /// are returned. Two lifetime slots are reserved; the pair publishes atomically.
    /// The same estimated-local-error and material-corner limits as the nonradiative
    /// adaptive enthalpy method apply. Solver refusals never widen a physical gate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        step: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
        -> Result<AdaptiveHistoryStep<EnthalpyRadiationStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps).map_err(conduction)?;
        phase.validate()?;
        let start = self.run.time_s;
        let cap = admit_proposal(cx, self.run.accepted_steps, self.run.max_steps, proposed_step_s, adaptive).map_err(conduction)?;
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.patches.iter().map(|p| &p.ambient_temperature), start, target_s, cap).map_err(conduction)?;
        let solve = |old: &[f64], from: f64, to: f64| -> Result<EnthalpyRadiationStepSolution, HeterogeneousEnthalpyError> {
            let boundary = events::boundary_before(cx, self.run.problem.boundary, to).map_err(conduction)?;
            self.run.storage.admit_prescribed_boundary(cx, &boundary)?;
            let source = events::source_before(&self.run.problem.source, to, old.len()).map_err(conduction)?;
            let patches = patches_before(cx, &self.patches, to).map_err(conduction)?;
            Ok(self.run.storage.advance_with_ambient_radiation(cx, ConductionProblem {
                mesh: self.run.problem.boundary.mesh(), boundary: &boundary, source: &source,
                material: self.run.problem.material, element_materials: self.run.problem.element_materials },
                self.run.problem.interfaces, old, to - from, step, &patches, self.config)?)
        };
        let selected = select_step(cx, start, end, adaptive, conduction, |middle, end| {
            let coarse = solve(&self.run.specific_h, start, end)?;
            let first = solve(&self.run.specific_h, start, middle)?;
            let second = solve(&first.conduction.specific_enthalpy_j_kg, middle, end)?;
            let error = phase_error(cx, &self.run, &coarse.conduction, &second.conduction, adaptive, phase)?;
            Ok((first, second, error))
        })?;
        self.run.publish(cx, selected.end_time_s, &selected.second_half.conduction)?;
        self.run.accepted_steps += 1;
        Ok(selected)
    }
}
