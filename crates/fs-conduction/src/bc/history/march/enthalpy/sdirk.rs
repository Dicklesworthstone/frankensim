//! Scheduled second-order stages on the original primary-enthalpy history.
//!
//! This owner changes only time quadrature. Storage, source k, material contact
//! and reflected radiation retain their original spatial residuals and source
//! identities. Stage carry is internal numerical history, never external heat.

use super::{EnthalpyHistoryRun, EnthalpyStepConfig, HeterogeneousEnthalpyError,
    conduction, invalid, poll};
use super::adaptive::{PhaseErrorTolerance, PhaseFields, phase_field_error};
use super::radiation::enclosure::MaterialEnclosureHistoryBinding;
use super::super::events::{self, HistoryEventStep,
    adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep, admit_proposal, select_step}};
use crate::interface::temperature::TemperatureContacts;
use crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use crate::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use crate::transient::enthalpy::sdirk::{
    EnthalpySdirk2Input, EnthalpySdirk2Stage, EnthalpySdirk2Step, SDIRK2_GAMMA,
};
use crate::{MaterialEnclosureRadiationPatch, MaterialPointwiseEnclosureRadiation};
use fs_exec::Cx;
use fs_scenario::signal::TimeSignal;

/// Persistent SDIRK2 evolution of an already admitted physical history.
///
/// Construct with `into_sdirk2` on an original `EnthalpyHistoryRun`, coupled
/// Newton `ContactEnthalpyHistoryRun`, or `MaterialEnclosureEnthalpyHistoryRun`.
/// Conversion preserves accepted fields, time, source declarations and budget.
/// Only the stiffly accurate second stage is published. Internal stages do not
/// consume accepted-step slots or appear as independently accepted timesteps.
pub struct Sdirk2EnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    run: EnthalpyHistoryRun<'h, 'm, 'c>,
    contacts: Option<(TemperatureContacts<'m, 'i>, ContactEnthalpyNewtonConfig)>,
    radiation: Option<MaterialEnclosureHistoryBinding<'h, 'c>>,
}

impl<'h, 'm, 'c, 'i> Sdirk2EnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    pub(super) fn from_parts(run: EnthalpyHistoryRun<'h, 'm, 'c>,
        contacts: Option<(TemperatureContacts<'m, 'i>, ContactEnthalpyNewtonConfig)>,
        radiation: Option<MaterialEnclosureHistoryBinding<'h, 'c>>) -> Self
    { Self { run, contacts, radiation } }

    /// Read-only accepted physical state; it cannot advance without its sources.
    #[must_use]
    pub fn conduction(&self) -> &EnthalpyHistoryRun<'h, 'm, 'c> { &self.run }

    /// Original material contact set, if present.
    #[must_use]
    pub fn contacts(&self) -> Option<&TemperatureContacts<'m, 'i>> {
        self.contacts.as_ref().map(|(contacts, _)| contacts)
    }

    /// Pinned material enclosure sources in their original surface order.
    #[must_use]
    pub fn material_enclosure_patches(&self) -> Option<&[MaterialEnclosureRadiationPatch<'c>]> {
        self.radiation.as_ref().map(MaterialEnclosureHistoryBinding::patches)
    }

    /// Sole external black reservoir; absent for a closed or missing enclosure.
    #[must_use]
    pub fn material_enclosure_ambient_temperature(&self) -> Option<&TimeSignal> {
        self.radiation.as_ref().and_then(MaterialEnclosureHistoryBinding::ambient_temperature)
    }

    /// Extend the original physical-interval budget without resetting history.
    pub fn extend_step_budget(&mut self, maximum: usize) -> Result<(), HeterogeneousEnthalpyError> {
        self.run.extend_step_budget(maximum)
    }

    /// Advance one complete SDIRK2 interval, sampling authored schedules at
    /// `start + gamma*dt` and the requested endpoint. Both samples use right
    /// limits. Caller-supplied intervals are not automatically event-split.
    /// Source, nonlinear, heat-balance or cancellation refusal publishes nothing.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time_s: f64, config: EnthalpyStepConfig)
        -> Result<EnthalpySdirk2Step, HeterogeneousEnthalpyError>
    {
        self.admit_one(cx)?;
        let solved = self.solve(cx, &self.run.specific_h, self.run.time_s, time_s, config, false)?;
        self.publish(cx, time_s, &solved.second_stage)?;
        Ok(solved)
    }

    /// End at the next source, boundary or open-bath knot. The final stage uses
    /// the exact interval-end left limit, preserving short held heat pulses.
    pub fn advance_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: EnthalpyStepConfig)
        -> Result<HistoryEventStep<EnthalpySdirk2Step>, HeterogeneousEnthalpyError>
    {
        self.admit_one(cx)?;
        let start_time_s = self.run.time_s;
        let end_time_s = self.next_endpoint(cx, target_s, max_step_s)?;
        let report = self.solve(cx, &self.run.specific_h, start_time_s, end_time_s, config, true)?;
        self.publish(cx, end_time_s, &report.second_stage)?;
        Ok(HistoryEventStep { start_time_s, end_time_s, report })
    }

    /// Compare one SDIRK2 interval with two actual SDIRK2 half intervals using
    /// the same retained phase-aware controller. Only the accepted pair is
    /// published, consuming two physical-interval slots and retaining each
    /// interval's separately integrated source/boundary/radiation accounts.
    ///
    /// The full coarse/fine discrepancy is used, with the original conservative
    /// unequal-half inflation. It is deliberately not divided by the smooth
    /// second-order factor three: phase and source knots may reduce regularity.
    /// This is local control, not a global or phase-transition-time certificate.
    /// `endpoint_solves` counts SDIRK2 interval solves; each has two implicit
    /// stage solves. Original stage-domain and solver refusals propagate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        config: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
        -> Result<AdaptiveHistoryStep<EnthalpySdirk2Step>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps).map_err(conduction)?;
        phase.validate()?;
        let start = self.run.time_s;
        let cap = admit_proposal(cx, self.run.accepted_steps, self.run.max_steps,
            proposed_step_s, adaptive).map_err(conduction)?;
        let end = self.next_endpoint(cx, target_s, cap)?;
        let selected = select_step(cx, start, end, adaptive, conduction, |middle, end| {
            let coarse = self.solve(cx, &self.run.specific_h, start, end, config, true)?;
            let first = self.solve(cx, &self.run.specific_h, start, middle, config, true)?;
            let second = self.solve(cx, &first.second_stage.specific_enthalpy_j_kg,
                middle, end, config, true)?;
            let error = phase_field_error(cx, &self.run,
                fields(&coarse.second_stage), fields(&second.second_stage), adaptive, phase)?;
            Ok((first, second, error))
        })?;
        self.publish(cx, selected.end_time_s, &selected.second_half.second_stage)?;
        self.run.accepted_steps += 1;
        Ok(selected)
    }

    fn admit_one(&self, cx: &Cx<'_>) -> Result<(), HeterogeneousEnthalpyError> {
        poll(cx, self.run.accepted_steps).map_err(conduction)?;
        if self.run.accepted_steps >= self.run.max_steps {
            return Err(invalid("SDIRK2 enthalpy history step budget exhausted"));
        }
        Ok(())
    }

    fn next_endpoint(&self, cx: &Cx<'_>, target: f64, cap: f64)
        -> Result<f64, HeterogeneousEnthalpyError>
    {
        events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.material_enclosure_ambient_temperature().into_iter(), self.run.time_s,
            target, cap).map_err(conduction)
    }

    #[allow(clippy::too_many_arguments)]
    fn solve(&self, cx: &Cx<'_>, old: &[f64], from: f64, to: f64,
        config: EnthalpyStepConfig, before: bool) -> Result<EnthalpySdirk2Step, HeterogeneousEnthalpyError>
    {
        let dt = to - from;
        let stage_time = SDIRK2_GAMMA.mul_add(dt, from);
        if ![from, to, dt, stage_time].iter().all(|x| x.is_finite()) || dt <= 0.0
            || !(from < stage_time && stage_time < to)
        {
            return Err(invalid("SDIRK2 history needs a finite positive interval and representable interior stage time"));
        }
        let history = self.run.problem.boundary;
        let first_boundary = history.at(cx, stage_time).map_err(conduction)?;
        let second_boundary = if before { events::boundary_before(cx, history, to) }
            else { history.at(cx, to) }.map_err(conduction)?;
        let first_source = self.run.problem.source.at(stage_time, old.len()).map_err(conduction)?;
        let second_source = if before { events::source_before(&self.run.problem.source, to, old.len()) }
            else { self.run.problem.source.at(to, old.len()) }.map_err(conduction)?;
        let first_problem = self.run.conduction_problem(&first_boundary, &first_source)?;
        let second_problem = self.run.conduction_problem(&second_boundary, &second_source)?;
        let contacts = self.contacts.as_ref().map(|(owner, config)| (owner, *config));
        let inputs = |radiation| EnthalpySdirk2Input { problem: first_problem,
            interfaces: self.run.problem.interfaces, radiation, contacts };
        let endpoint = |radiation| EnthalpySdirk2Input { problem: second_problem,
            interfaces: self.run.problem.interfaces, radiation, contacts };
        match &self.radiation {
            None => Ok(self.run.storage.advance_sdirk2(cx,
                [inputs(EnthalpyRadiation::None), endpoint(EnthalpyRadiation::None)], old, dt, config)?),
            Some(binding) => {
                let first_base = binding.sample(cx, history.mesh(), &first_boundary, stage_time, false)?;
                let second_base = binding.sample(cx, history.mesh(), &second_boundary, to, before)?;
                let first = MaterialPointwiseEnclosureRadiation::new(cx, &first_base, binding.patches()).map_err(conduction)?;
                let second = MaterialPointwiseEnclosureRadiation::new(cx, &second_base, binding.patches()).map_err(conduction)?;
                Ok(self.run.storage.advance_sdirk2(cx,
                    [inputs(EnthalpyRadiation::MaterialEnclosure(&first)),
                     endpoint(EnthalpyRadiation::MaterialEnclosure(&second))], old, dt, config)?)
            }
        }
    }

    fn publish(&mut self, cx: &Cx<'_>, time_s: f64, stage: &EnthalpySdirk2Stage)
        -> Result<(), HeterogeneousEnthalpyError>
    {
        self.run.publish_fields(cx, time_s, &stage.specific_enthalpy_j_kg,
            &stage.temperature, &stage.liquid_mass_fraction)
    }
}

fn fields(stage: &EnthalpySdirk2Stage) -> PhaseFields<'_> {
    PhaseFields { h: &stage.specific_enthalpy_j_kg, temperature: &stage.temperature,
        fraction: &stage.liquid_mass_fraction }
}
