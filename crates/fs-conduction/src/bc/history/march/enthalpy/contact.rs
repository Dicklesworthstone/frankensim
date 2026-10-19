//! Scheduled phase-changing bodies coupled by source-backed contact resistance.
//!
//! Every accepted endpoint closes R(T_interface), k(T, phase), enthalpy storage
//! and optional gray radiation through their existing owners. Contact snapshots
//! are inner trials, not frozen replacements for the material law. The contact
//! set owns ALL fixed and variable interfaces; a second set is never discarded.

use super::{EnthalpyHistoryRun, conduction, invalid, natural_only, poll};
use super::adaptive::{PhaseErrorTolerance, phase_error};
use super::super::{events, radiation::{RadiationRegionHistory, admit_histories, patches_before}};
use events::adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep, admit_proposal, select_step};
use crate::{AmbientRadiationConfig, ConductionProblem};
use crate::interface::temperature::{ContactIterationConfig, TemperatureContacts};
use crate::interface::temperature::enthalpy::ContactEnthalpyStep;
use crate::transient::enthalpy::EnthalpyStepConfig;
use crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyError;
use fs_exec::Cx;

/// One phase-history run and its complete, immutable material-contact model.
/// No mutable accessor can advance the inner run while omitting its contacts.
/// Only temperature coordinates of the contact law evolve; geometry, pressure,
/// selected source claims and reference mass remain fixed.
pub struct ContactEnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    run: EnthalpyHistoryRun<'h, 'm, 'c>,
    contacts: TemperatureContacts<'m, 'i>,
    coupling: ContactIterationConfig,
    radiation: Option<(Vec<RadiationRegionHistory>, AmbientRadiationConfig)>,
}

impl<'h, 'm, 'c, 'i> ContactEnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    /// Attach the full fixed/variable contact set before the first physical step.
    /// Declare any constant interfaces in `TemperatureContacts::new`, not in the
    /// base problem as well. Initial contact source queries are checked without
    /// advancing thermal history. The production solver validates complete work
    /// products against each later endpoint's actual Newton/radiation budgets.
    ///
    /// # Errors
    /// Advanced run, duplicate interface ownership, foreign mesh, invalid scalar
    /// contact policy, source-domain refusal or cancellation. No step is taken.
    pub fn new(cx: &Cx<'_>, run: EnthalpyHistoryRun<'h, 'm, 'c>,
        contacts: TemperatureContacts<'m, 'i>, coupling: ContactIterationConfig)
        -> Result<Self, HeterogeneousEnthalpyError>
    {
        poll(cx, 0).map_err(conduction)?;
        if run.accepted_steps != 0 || run.problem.interfaces.is_some() {
            return Err(invalid("attach the complete contact set to an unadvanced run with no separately declared interfaces"));
        }
        if coupling.max_iterations == 0 || !coupling.relaxation.is_finite()
            || coupling.relaxation <= 0.0 || coupling.relaxation > 1.0
            || !coupling.residual_tolerance_j.is_finite() || coupling.residual_tolerance_j <= 0.0
        {
            return Err(invalid("contact history needs positive iteration/residual budgets and relaxation in (0,1]"));
        }
        contacts.check_mesh(run.problem.boundary.mesh()).map_err(conduction)?;
        let boundary = run.problem.boundary.at(cx, run.time_s).map_err(conduction)?;
        contacts.evaluate(cx, &boundary, &run.temperature).map_err(conduction)?;
        poll(cx, 0).map_err(conduction)?;
        Ok(Self { run, contacts, coupling, radiation: None })
    }

    /// Retain sourced radiating surfaces and independent reservoir histories.
    /// Radiation cannot be attached twice or after a physical step. The original
    /// patch admission and surface-validity checks apply; no emissivity is inferred.
    pub fn with_radiation(mut self, cx: &Cx<'_>, patches: Vec<RadiationRegionHistory>,
        config: AmbientRadiationConfig) -> Result<Self, HeterogeneousEnthalpyError>
    {
        poll(cx, 0).map_err(conduction)?;
        if self.run.accepted_steps != 0 || self.radiation.is_some() {
            return Err(invalid("attach radiation once, before the first material-contact history step"));
        }
        config.validate().map_err(conduction)?;
        admit_histories(cx, self.run.problem.boundary, self.run.time_s, &patches).map_err(conduction)?;
        self.radiation = Some((patches, config));
        Ok(self)
    }

    /// Accepted caloric fields, reference materials, time and physical-step count.
    #[must_use]
    pub fn conduction(&self) -> &EnthalpyHistoryRun<'h, 'm, 'c> { &self.run }

    /// Original contact geometry and source laws, never a mutable frozen snapshot.
    #[must_use]
    pub fn contacts(&self) -> &TemperatureContacts<'m, 'i> { &self.contacts }

    /// Original outer contact policy, unchanged during retries and continuation.
    #[must_use]
    pub const fn coupling(&self) -> ContactIterationConfig { self.coupling }

    /// Original radiative declarations, absent when radiation was not selected.
    #[must_use]
    pub fn radiation(&self) -> Option<(&[RadiationRegionHistory], AmbientRadiationConfig)> {
        self.radiation.as_ref().map(|(patches, config)| (patches.as_slice(), *config))
    }

    /// Extend the accepted-step allowance without changing contacts or history.
    pub fn extend_step_budget(&mut self, total: usize) -> Result<(), HeterogeneousEnthalpyError> {
        self.run.extend_step_budget(total)
    }

    /// Solve one explicit right-continuous endpoint through the original contact
    /// solver. Reports include both applied and re-evaluated physical contacts,
    /// their source receipts, latent state, local residual and global joule audit.
    /// Internal contact heat is never counted as an external source.
    ///
    /// # Errors
    /// Budget/time/source/phase/contact/radiation/solver or cancellation refusal.
    /// Failure publishes no temperature, enthalpy, liquid fraction or physical time.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time_s: f64, step: EnthalpyStepConfig)
        -> Result<ContactEnthalpyStep, HeterogeneousEnthalpyError>
    {
        self.admit_one(cx)?;
        let dt = time_s - self.run.time_s;
        if !time_s.is_finite() || !dt.is_finite() || dt <= 0.0 {
            return Err(invalid("contact history endpoint requires a later finite time and positive finite interval"));
        }
        let solved = self.solve(cx, &self.run.specific_h, self.run.time_s, time_s, step, false)?;
        self.run.publish(cx, time_s, &solved.conduction)?;
        Ok(solved)
    }

    /// Advance one interval capped by every source, boundary and reservoir knot.
    /// Use exact pre-switch loads at the interval end, preserving short pulses.
    /// The full contact/radiation/Newton solve shares one immutable old enthalpy.
    pub fn advance_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, step: EnthalpyStepConfig)
        -> Result<events::HistoryEventStep<ContactEnthalpyStep>, HeterogeneousEnthalpyError>
    {
        self.admit_one(cx)?;
        let start_time_s = self.run.time_s;
        let end_time_s = self.next_endpoint(cx, target_s, max_step_s)?;
        let report = self.solve(cx, &self.run.specific_h, start_time_s, end_time_s, step, true)?;
        self.run.publish(cx, end_time_s, &report.conduction)?;
        Ok(events::HistoryEventStep { start_time_s, end_time_s, report })
    }

    /// Adapt using temperature, specific-enthalpy AND liquid-fraction errors.
    /// Each of the three endpoint trials closes its own material-contact and
    /// radiative feedback. Publish only the real two-half-step trajectory after
    /// all local-error and physical gates pass. Two lifetime slots are reserved.
    /// Both accepted reports retain complete contact and physical heat accounts.
    /// The estimator is local, not a global bound or a phase-corner certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64,
        proposed_step_s: f64, step: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig,
        phase: PhaseErrorTolerance)
        -> Result<AdaptiveHistoryStep<ContactEnthalpyStep>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps).map_err(conduction)?;
        phase.validate()?;
        let cap = admit_proposal(cx, self.run.accepted_steps, self.run.max_steps,
            proposed_step_s, adaptive).map_err(conduction)?;
        let start = self.run.time_s;
        let end = self.next_endpoint(cx, target_s, cap)?;
        let selected = select_step(cx, start, end, adaptive, conduction, |middle, end| {
            let coarse = self.solve(cx, &self.run.specific_h, start, end, step, true)?;
            let first = self.solve(cx, &self.run.specific_h, start, middle, step, true)?;
            let second = self.solve(cx, &first.conduction.specific_enthalpy_j_kg, middle, end, step, true)?;
            let error = phase_error(cx, &self.run, &coarse.conduction, &second.conduction, adaptive, phase)?;
            Ok((first, second, error))
        })?;
        self.run.publish(cx, selected.end_time_s, &selected.second_half.conduction)?;
        self.run.accepted_steps += 1;
        Ok(selected)
    }

    fn admit_one(&self, cx: &Cx<'_>) -> Result<(), HeterogeneousEnthalpyError> {
        poll(cx, self.run.accepted_steps).map_err(conduction)?;
        if self.run.accepted_steps >= self.run.max_steps {
            return Err(invalid("material-contact history step budget exhausted"));
        }
        Ok(())
    }

    fn next_endpoint(&self, cx: &Cx<'_>, target: f64, cap: f64)
        -> Result<f64, HeterogeneousEnthalpyError>
    {
        let reservoirs = self.radiation.iter().flat_map(|(patches, _)|
            patches.iter().map(|patch| &patch.ambient_temperature));
        events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            reservoirs, self.run.time_s, target, cap).map_err(conduction)
    }

    fn solve(&self, cx: &Cx<'_>, old: &[f64], from: f64, to: f64,
        step: EnthalpyStepConfig, before: bool) -> Result<ContactEnthalpyStep, HeterogeneousEnthalpyError>
    {
        let history = self.run.problem.boundary;
        let boundary = if before { events::boundary_before(cx, history, to) }
            else { history.at(cx, to) }.map_err(conduction)?;
        natural_only(&boundary)?;
        let source = if before { events::source_before(&self.run.problem.source, to, old.len()) }
            else { self.run.problem.source.at(to, old.len()) }.map_err(conduction)?;
        let problem = ConductionProblem { mesh: history.mesh(), boundary: &boundary, source: &source,
            material: self.run.problem.material, element_materials: self.run.problem.element_materials };
        match &self.radiation {
            None => Ok(self.contacts.advance_enthalpy(cx, &self.run.storage, problem,
                old, to - from, step, self.coupling)?),
            Some((histories, config)) => {
                let patches = if before { patches_before(cx, histories, to).map_err(conduction)? }
                    else {
                        let mut patches = Vec::with_capacity(histories.len());
                        for (index, history) in histories.iter().enumerate() {
                            poll(cx, index).map_err(conduction)?;
                            patches.push(history.at(to).map_err(conduction)?);
                        }
                        patches
                    };
                Ok(self.contacts.advance_enthalpy_with_ambient_radiation(cx, &self.run.storage,
                    problem, old, to - from, step, self.coupling, &patches, *config)?)
            }
        }
    }
}
