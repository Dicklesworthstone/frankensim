//! Scheduled material reflection through the original primary-enthalpy owner.
//!
//! Closed enclosure exchange remains internal; an open enclosure owns one
//! explicitly scheduled black reservoir. Original traces, view factors and
//! source claims survive every candidate. Only accepted h/T/f advances time.

use super::super::{
    EnthalpyHistoryRun, EnthalpyStepConfig, EnthalpyStepSolution,
    HeterogeneousEnthalpyError, conduction, invalid, poll,
    adaptive::{PhaseErrorTolerance, phase_error},
};
use super::super::super::events::{
    self, HistoryEventStep,
    adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep, admit_proposal, select_step},
};
use crate::{
    ConductionError, MaterialEnclosureRadiationPatch, MaterialPointwiseEnclosureRadiation,
    PointwiseEnclosureRadiation, SURFACE_EMISSIVITY_PROPERTY, ScalarField,
    SurfaceEmissivity, TEMPERATURE_DIMS, ThermalBoundary,
};
use fs_exec::Cx;
use fs_material::state_point::MaterialPropertySelection;
use fs_scenario::signal::TimeSignal;

/// Persistent material h/T/f with gray reflected exchange and optional bath.
///
/// The supplied pointwise owner fixes the original mesh, named traces, view
/// factors and radiation tolerances. Its immutable geometry and optical factors
/// are shared across interval candidates. Material emissivity and the reflected
/// operator are refreshed by the original material owner at each Newton state.
/// Pressure, surface finish and other source conditions stay explicitly fixed.
pub struct MaterialEnclosureEnthalpyHistoryRun<'h, 'm, 'c> {
    run: EnthalpyHistoryRun<'h, 'm, 'c>,
    base: &'h PointwiseEnclosureRadiation<'h>,
    patches: Vec<MaterialEnclosureRadiationPatch<'c>>,
    reference_emissivities: Vec<SurfaceEmissivity>,
    ambient_temperature: Option<TimeSignal>,
}

impl<'h, 'm, 'c> MaterialEnclosureEnthalpyHistoryRun<'h, 'm, 'c> {
    /// Attach a closed or open material enclosure before the first endpoint.
    ///
    /// An open owner requires exactly one absolute-temperature schedule; a
    /// closed owner requires `None`. Source selections are pinned here. Every
    /// trace must use the exact original caloric object assigned to its storage
    /// nodes. The preparation snapshot's emissivity never replaces that law.
    ///
    /// # Errors
    /// Advanced history, wrong mesh/partition/chart, invalid source or schedule,
    /// unsupported spectral data, mismatched reservoir kind or cancellation.
    /// Construction does not advance any physical state.
    pub fn new(
        cx: &Cx<'_>, run: EnthalpyHistoryRun<'h, 'm, 'c>,
        base: &'h PointwiseEnclosureRadiation<'h>,
        mut patches: Vec<MaterialEnclosureRadiationPatch<'c>>,
        ambient_temperature: Option<TimeSignal>,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0).map_err(conduction)?;
        if run.accepted_steps() != 0 {
            return Err(invalid("material enclosure must be attached before the first enthalpy step"));
        }
        if base.is_closed() != ambient_temperature.is_none() {
            return Err(invalid("only an open material enclosure requires one black-reservoir history"));
        }
        let mut reference_emissivities = Vec::with_capacity(patches.len());
        for (index, patch) in patches.iter_mut().enumerate() {
            poll(cx, index).map_err(conduction)?;
            let (reference, _) = SurfaceEmissivity::from_condition(
                patch.region, patch.material, patch.selection.clone(),
            ).map_err(conduction)?;
            patch.selection = MaterialPropertySelection::PinnedByProperty(vec![(
                SURFACE_EMISSIVITY_PROPERTY.into(), reference.receipt().selected,
            )]);
            reference_emissivities.push(reference);
        }
        let prepared = Self { run, base, patches, reference_emissivities, ambient_temperature };
        let boundary = prepared.run.problem.boundary.at(cx, prepared.run.time_s())
            .map_err(conduction)?;
        let base = prepared.sample(cx, &boundary, prepared.run.time_s(), false)?;
        let radiation = MaterialPointwiseEnclosureRadiation::new(cx, &base, &prepared.patches)
            .map_err(conduction)?;
        if !radiation.matches_curves(cx, |vertex| prepared.run.storage.phase_curve_at(vertex)
            .expect("admitted history vertex has an original caloric chart"))
            .map_err(conduction)?
        {
            return Err(invalid("material enclosure history requires the exact storage caloric objects"));
        }
        radiation.load_and_report(cx, &prepared.run.specific_h, &prepared.run.temperature)
            .map_err(conduction)?;
        poll(cx, prepared.patches.len()).map_err(conduction)?;
        Ok(prepared)
    }

    /// Accepted enthalpy, temperature, fractions, reference storage and clock.
    #[must_use]
    pub fn conduction(&self) -> &EnthalpyHistoryRun<'h, 'm, 'c> { &self.run }

    /// Original material declarations with pinned source claims in surface order.
    #[must_use]
    pub fn patches(&self) -> &[MaterialEnclosureRadiationPatch<'c>] { &self.patches }

    /// Preparation snapshots only; accepted reports carry actual nodal states.
    #[must_use]
    pub fn reference_emissivities(&self) -> &[SurfaceEmissivity] {
        &self.reference_emissivities
    }

    /// The sole black-bath history, absent for a closed enclosure.
    #[must_use]
    pub fn ambient_temperature(&self) -> Option<&TimeSignal> {
        self.ambient_temperature.as_ref()
    }

    /// Resume with more accepted endpoints without replacing any physical owner.
    pub fn extend_step_budget(&mut self, maximum: usize) -> Result<(), HeterogeneousEnthalpyError> {
        self.run.extend_step_budget(maximum)
    }

    /// Advance a supplied endpoint with right-continuous schedule sampling.
    ///
    /// # Errors
    /// Original source, time, budget, nonlinear, energy or cancellation refusal.
    /// A failed proposal leaves every accepted field and the clock unchanged.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time_s: f64, config: EnthalpyStepConfig)
        -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError>
    {
        let (dt, boundary, source) = self.run.prepare_endpoint(cx, time_s)?;
        let report = self.solve(cx, &self.run.specific_h, &boundary, &source,
            dt, time_s, false, config)?;
        self.run.publish(cx, time_s, &report)?;
        Ok(report)
    }

    /// Clip to source, boundary and bath events, using interval-end left limits.
    /// Held radiation pulses therefore retain heat delivered before a switch.
    ///
    /// # Errors
    /// Invalid target/cap or an original physical/source/budget refusal. No
    /// accepted state changes when any part of the interval is refused.
    pub fn advance_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: EnthalpyStepConfig)
        -> Result<HistoryEventStep<EnthalpyStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        if self.run.accepted_steps() >= self.run.max_steps() {
            return Err(invalid("enthalpy history step budget exhausted"));
        }
        let start_time_s = self.run.time_s();
        let end_time_s = events::next_endpoint(cx, self.run.problem.boundary,
            &self.run.problem.source, self.ambient_temperature.iter(),
            start_time_s, target_s, max_step_s).map_err(conduction)?;
        let report = self.solve_before(cx, &self.run.specific_h, start_time_s, end_time_s, config)?;
        self.run.publish(cx, end_time_s, &report)?;
        Ok(HistoryEventStep { start_time_s, end_time_s, report })
    }

    /// Control local time error in h, T and liquid fraction, including latent
    /// evolution at constant temperature. Accepted half-step heat reports stay
    /// separate; the two endpoints publish atomically through the original run.
    ///
    /// # Errors
    /// Invalid policy, insufficient endpoint allowance, original physical/source
    /// refusal, cancellation or an unmet error target at the floor/attempt cap.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        step: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
        -> Result<AdaptiveHistoryStep<EnthalpyStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        phase.validate()?;
        let start = self.run.time_s();
        let cap = admit_proposal(cx, self.run.accepted_steps(), self.run.max_steps(),
            proposed_step_s, adaptive).map_err(conduction)?;
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.ambient_temperature.iter(), start, target_s, cap).map_err(conduction)?;
        let selected = select_step(cx, start, end, adaptive, conduction, |middle, end| {
            let coarse = self.solve_before(cx, &self.run.specific_h, start, end, step)?;
            let first = self.solve_before(cx, &self.run.specific_h, start, middle, step)?;
            let second = self.solve_before(cx, &first.specific_enthalpy_j_kg, middle, end, step)?;
            let error = phase_error(cx, &self.run, &coarse, &second, adaptive, phase)?;
            Ok((first, second, error))
        })?;
        self.run.publish(cx, selected.end_time_s, &selected.second_half)?;
        self.run.accepted_steps += 1;
        Ok(selected)
    }

    fn sample<'s>(&'s self, cx: &Cx<'_>, boundary: &'s ThermalBoundary, time_s: f64, before: bool)
        -> Result<PointwiseEnclosureRadiation<'s>, HeterogeneousEnthalpyError>
    {
        let base = self.base.at_history_boundary(cx, self.run.problem.boundary.mesh(), boundary)
            .map_err(conduction)?;
        let Some(signal) = &self.ambient_temperature else { return Ok(base); };
        let value = if before { events::value_before(signal, time_s).map_err(conduction)? }
        else { signal.eval(time_s).map_err(|error| conduction(ConductionError::Radiation {
            surface: "open-enclosure-bath".into(),
            what: format!("radiative reservoir history: {error}"),
            fix: "supply a valid absolute-temperature history".into(),
        }))? };
        if value.dims != TEMPERATURE_DIMS {
            return Err(conduction(ConductionError::Dimensions {
                context: "material enclosure black-reservoir history".into(),
                expected: TEMPERATURE_DIMS.0, found: value.dims.0,
            }));
        }
        base.with_ambient_temperature_k(value.value).map_err(conduction)
    }

    fn solve_before(&self, cx: &Cx<'_>, old: &[f64], from: f64, to: f64,
        config: EnthalpyStepConfig) -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError>
    {
        let boundary = events::boundary_before(cx, self.run.problem.boundary, to).map_err(conduction)?;
        self.run.storage.admit_prescribed_boundary(cx, &boundary)?;
        let source = events::source_before(&self.run.problem.source, to, old.len()).map_err(conduction)?;
        self.solve(cx, old, &boundary, &source, to - from, to, true, config)
    }

    #[allow(clippy::too_many_arguments)]
    fn solve(&self, cx: &Cx<'_>, old: &[f64], boundary: &ThermalBoundary,
        source: &ScalarField, dt_s: f64, time_s: f64, before: bool, config: EnthalpyStepConfig)
        -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError>
    {
        let base = self.sample(cx, boundary, time_s, before)?;
        let radiation = MaterialPointwiseEnclosureRadiation::new(cx, &base, &self.patches)
            .map_err(conduction)?;
        Ok(self.run.storage.advance_with_material_enclosure_radiation(cx,
            self.run.conduction_problem(boundary, source)?, self.run.problem.interfaces,
            old, dt_s, config, &radiation)?)
    }
}
