//! Time adaptation that cannot confuse a latent plateau with an unchanged state.

use super::{EnthalpyHistoryRun, EnthalpyStepConfig, EnthalpyStepSolution,
    HeterogeneousEnthalpyError, conduction, invalid, poll};
use super::super::events::{self, adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep,
    admit_proposal, component_error, select_step, temperature_error}};
use crate::ConductionProblem;
use fs_exec::Cx;

/// Additional error coordinates for a primary-enthalpy trajectory.
/// All three gates (temperature, enthalpy, liquid fraction) must pass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseErrorTolerance {
    /// Positive absolute local specific-enthalpy allowance, J/kg.
    pub absolute_enthalpy_j_kg: f64,
    /// Relative allowance against the larger coarse/fine enthalpy CHANGE from
    /// history, not the arbitrary enthalpy reference offset. Must be in [0,1).
    pub relative_enthalpy: f64,
    /// Positive absolute allowance on liquid mass fraction, at most one.
    pub absolute_liquid_fraction: f64,
}
impl PhaseErrorTolerance {
    /// Validate every caloric tolerance before any candidate solve.
    pub fn validate(self) -> Result<(), HeterogeneousEnthalpyError> {
        if !self.absolute_enthalpy_j_kg.is_finite() || self.absolute_enthalpy_j_kg <= 0.0
            || !self.relative_enthalpy.is_finite() || !(0.0..1.0).contains(&self.relative_enthalpy)
            || !self.absolute_liquid_fraction.is_finite()
            || self.absolute_liquid_fraction <= 0.0 || self.absolute_liquid_fraction > 1.0
        { return Err(invalid("phase time-error control requires positive finite J/kg and fraction allowances and relative tolerance in [0,1)")); }
        Ok(())
    }
}

impl EnthalpyHistoryRun<'_, '_, '_> {
    /// Advance two original implicit half steps, accepting only a small estimated
    /// discrepancy in temperature, specific enthalpy AND liquid fraction.
    ///
    /// Each candidate uses the same retained phase charts, reference masses,
    /// k(T), phase multipliers and contacts. The event controller limits the
    /// interval first; held loads retain their exact interval-end left limits.
    /// No apparent capacity, extrapolated caloric state or phase reset is used.
    /// The accepted pair consumes two lifetime slots and publishes atomically.
    ///
    /// This is a local estimator, not a phase-transition-time certificate or a
    /// global error bound. Kinks in constitutive charts need not have smooth
    /// first-order asymptotics. Original Newton/material/energy refusals propagate
    /// unchanged; they never cause an untested field to be accepted.
    ///
    /// # Errors
    /// Invalid policy/time, cancellation, insufficient accepted slots, original
    /// solver refusal, or unmet temporal accuracy at the floor/attempt limit.
    /// Every failure preserves enthalpy, temperature, liquid fraction and time.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        step: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
        -> Result<AdaptiveHistoryStep<EnthalpyStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.accepted_steps).map_err(conduction)?;
        phase.validate()?;
        let start = self.time_s;
        let cap = admit_proposal(cx, self.accepted_steps, self.max_steps, proposed_step_s, adaptive).map_err(conduction)?;
        let end = events::next_endpoint(cx, self.problem.boundary, &self.problem.source,
            std::iter::empty(), start, target_s, cap).map_err(conduction)?;
        let solve = |old: &[f64], from: f64, to: f64| -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError> {
            let boundary = events::boundary_before(cx, self.problem.boundary, to).map_err(conduction)?;
            self.storage.admit_prescribed_boundary(cx, &boundary)?;
            let source = events::source_before(&self.problem.source, to, old.len()).map_err(conduction)?;
            Ok(self.storage.advance(cx, ConductionProblem { mesh: self.problem.boundary.mesh(),
                boundary: &boundary, material: self.problem.material, element_materials: self.problem.element_materials,
                source: &source }, self.problem.interfaces, old, to - from, step)?)
        };
        let selected = select_step(cx, start, end, adaptive, conduction, |middle, end| {
            let coarse = solve(&self.specific_h, start, end)?;
            let first = solve(&self.specific_h, start, middle)?;
            let second = solve(&first.specific_enthalpy_j_kg, middle, end)?;
            let error = phase_error(cx, self, &coarse, &second, adaptive, phase)?;
            Ok((first, second, error))
        })?;
        self.publish(cx, selected.end_time_s, &selected.second_half)?;
        self.accepted_steps += 1; // Two reserved physical half steps, one atomic publication.
        Ok(selected)
    }
}

pub(super) fn phase_error(cx: &Cx<'_>, run: &EnthalpyHistoryRun<'_, '_, '_>,
    coarse: &EnthalpyStepSolution, fine: &EnthalpyStepSolution,
    adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
    -> Result<f64, HeterogeneousEnthalpyError>
{
    let n = run.specific_h.len();
    if [coarse.specific_enthalpy_j_kg.len(), fine.specific_enthalpy_j_kg.len(),
        coarse.liquid_mass_fraction.len(), fine.liquid_mass_fraction.len()].iter().any(|&len| len != n)
    { return Err(invalid("adaptive caloric fields must preserve the original nodal shape")); }
    let mut error = temperature_error(cx, &run.temperature, &coarse.temperature,
        &fine.temperature, adaptive).map_err(conduction)?;
    for i in 0..n {
        if i % 512 == 0 { poll(cx, i).map_err(conduction)?; }
        error = error.max(component_error(run.specific_h[i], coarse.specific_enthalpy_j_kg[i],
            fine.specific_enthalpy_j_kg[i], phase.absolute_enthalpy_j_kg, phase.relative_enthalpy).map_err(conduction)?);
        error = error.max(component_error(run.liquid_mass_fraction[i], coarse.liquid_mass_fraction[i],
            fine.liquid_mass_fraction[i], phase.absolute_liquid_fraction, 0.0).map_err(conduction)?);
    }
    Ok(error)
}
