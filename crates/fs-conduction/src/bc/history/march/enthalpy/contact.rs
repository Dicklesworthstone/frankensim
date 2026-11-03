//! Scheduled phase-changing bodies coupled by source-backed contact resistance.
//!
//! Every accepted endpoint closes R(T_interface), k(T, phase), enthalpy storage
//! and optional gray radiation through their existing owners. Contact snapshots
//! are inner trials, not frozen replacements for the material law. The contact
//! set owns ALL fixed and variable interfaces; a second set is never discarded.

mod enclosure;

use super::{EnthalpyHistoryRun, conduction, invalid, natural_only, poll};
use super::adaptive::{PhaseErrorTolerance, phase_error};
use super::radiation::enclosure::MaterialEnclosureHistoryBinding;
use super::super::{events, radiation::{RadiationRegionHistory, admit_histories, patches_before}};
use events::adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep, admit_proposal, select_step};
use crate::AmbientRadiationConfig;
use crate::interface::temperature::{ContactIterationConfig, TemperatureContacts};
use crate::interface::temperature::enthalpy::{ContactEnthalpyNewtonConfig, ContactEnthalpyStep};
use crate::transient::enthalpy::{EnthalpyError, EnthalpyStepConfig};
use crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyError;
use fs_exec::Cx;

/// One phase-history run and its complete, immutable material-contact model.
/// No mutable accessor can advance the inner run while omitting its contacts.
/// The coupled Newton method evolves each body's own temperature and liquid
/// fraction through its original caloric chart. Geometry, pressure, selected
/// source claims and reference mass remain fixed. The original outer-iteration
/// method remains available for explicitly temperature-only contact laws.
pub struct ContactEnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    run: EnthalpyHistoryRun<'h, 'm, 'c>,
    contacts: TemperatureContacts<'m, 'i>,
    coupling: ContactHistoryCoupling,
    radiation: Option<(Vec<RadiationRegionHistory>, AmbientRadiationConfig)>,
    material_enclosure: Option<MaterialEnclosureHistoryBinding<'h, 'c>>,
}

#[derive(Debug, Clone, Copy)]
enum ContactHistoryCoupling {
    Legacy(ContactIterationConfig),
    Newton(ContactEnthalpyNewtonConfig),
}

impl<'h, 'm, 'c, 'i> ContactEnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    /// Retain the original coupled contact/enclosure model and accepted state
    /// while switching to second-order SDIRK2 time integration.
    ///
    /// # Errors
    /// The outer temperature-only iteration and mean-patch radiation cannot
    /// supply the original full-h stage residual required by this method.
    pub fn into_sdirk2(self)
        -> Result<super::sdirk::Sdirk2EnthalpyHistoryRun<'h, 'm, 'c, 'i>, HeterogeneousEnthalpyError>
    {
        let config = self.newton_config()
            .ok_or_else(|| invalid("SDIRK2 contact history requires the retained coupled enthalpy Newton policy"))?;
        if self.radiation.is_some() {
            return Err(invalid("SDIRK2 contact history requires an explicit pointwise radiation owner"));
        }
        Ok(super::sdirk::Sdirk2EnthalpyHistoryRun::from_parts(
            self.run, Some((self.contacts, config)), self.material_enclosure,
        ))
    }

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
        Ok(Self { run, contacts, coupling: ContactHistoryCoupling::Legacy(coupling),
            radiation: None, material_enclosure: None })
    }

    /// Attach original phase-aware contacts to the existing full enthalpy Newton
    /// solver before the first physical step. Every trial uses each body's
    /// independent P1 T(h) and f(h) trace and the complete source contact tangent.
    /// The latent fraction remains active when dT/dh is zero.
    ///
    /// The complete coupled policy is retained for continuation. Every advance
    /// must pass the same `config.step`; mismatched controls refuse explicitly.
    /// Matching/nonmatching geometry, source claims, reference masses and
    /// source-owned k(h_h) remain with their original physical owners.
    ///
    /// # Errors
    /// Advanced history, duplicate interfaces, invalid Newton work/physical
    /// gates, foreign mesh or caloric objects, original source refusal or
    /// cancellation. Construction performs no thermal step.
    pub fn new_newton(cx: &Cx<'_>, run: EnthalpyHistoryRun<'h, 'm, 'c>,
        contacts: TemperatureContacts<'m, 'i>, config: ContactEnthalpyNewtonConfig)
        -> Result<Self, HeterogeneousEnthalpyError>
    {
        poll(cx, 0).map_err(conduction)?;
        if run.accepted_steps != 0 || run.problem.interfaces.is_some() {
            return Err(invalid("attach the complete contact set to an unadvanced run with no separately declared interfaces"));
        }
        config.validate()?;
        config.step.newton.validate().map_err(EnthalpyError::Newton)?;
        if config.step.max_newton_iterations == 0 || !config.step.energy_tolerance_j.is_finite()
            || config.step.energy_tolerance_j <= 0.0
        {
            return Err(invalid("contact history requires positive nonlinear work and finite positive energy tolerance"));
        }
        contacts.check_mesh(run.problem.boundary.mesh()).map_err(conduction)?;
        let boundary = run.problem.boundary.at(cx, run.time_s).map_err(conduction)?;
        contacts.evaluate_enthalpy(cx, &boundary, &run.storage, &run.specific_h)?;
        poll(cx, 0).map_err(conduction)?;
        Ok(Self { run, contacts, coupling: ContactHistoryCoupling::Newton(config),
            radiation: None, material_enclosure: None })
    }

    /// Retain sourced radiating surfaces and independent reservoir histories.
    /// Radiation cannot be attached twice or after a physical step. The original
    /// patch admission and surface-validity checks apply; no emissivity is inferred.
    pub fn with_radiation(mut self, cx: &Cx<'_>, patches: Vec<RadiationRegionHistory>,
        config: AmbientRadiationConfig) -> Result<Self, HeterogeneousEnthalpyError>
    {
        poll(cx, 0).map_err(conduction)?;
        if self.run.accepted_steps != 0 || self.radiation.is_some() || self.material_enclosure.is_some() {
            return Err(invalid("attach radiation once, before the first material-contact history step"));
        }
        if matches!(self.coupling, ContactHistoryCoupling::Newton(_)) {
            return Err(invalid("area-mean radiation requires the iterated contact history method; coupled Newton needs an explicit pointwise radiation owner"));
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

    /// Original outer contact policy; absent for the coupled enthalpy Newton method.
    #[must_use]
    pub const fn coupling(&self) -> Option<ContactIterationConfig> {
        match self.coupling { ContactHistoryCoupling::Legacy(config) => Some(config),
            ContactHistoryCoupling::Newton(_) => None }
    }

    /// Complete retained coupled Newton policy, absent for outer iteration.
    #[must_use]
    pub const fn newton_config(&self) -> Option<ContactEnthalpyNewtonConfig> {
        match self.coupling { ContactHistoryCoupling::Legacy(_) => None,
            ContactHistoryCoupling::Newton(config) => Some(config) }
    }

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
            patches.iter().map(|patch| &patch.ambient_temperature))
            .chain(self.material_enclosure.iter().filter_map(
                MaterialEnclosureHistoryBinding::ambient_temperature));
        events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            reservoirs, self.run.time_s, target, cap).map_err(conduction)
    }

    fn solve(&self, cx: &Cx<'_>, old: &[f64], from: f64, to: f64,
        step: EnthalpyStepConfig, before: bool) -> Result<ContactEnthalpyStep, HeterogeneousEnthalpyError>
    {
        let history = self.run.problem.boundary;
        let boundary = if before { events::boundary_before(cx, history, to) }
            else { history.at(cx, to) }.map_err(conduction)?;
        match self.coupling {
            ContactHistoryCoupling::Legacy(_) => natural_only(&boundary)?,
            ContactHistoryCoupling::Newton(_) => self.run.storage.admit_prescribed_boundary(cx, &boundary)?,
        }
        let source = if before { events::source_before(&self.run.problem.source, to, old.len()) }
            else { self.run.problem.source.at(to, old.len()) }.map_err(conduction)?;
        let problem = self.run.conduction_problem(&boundary, &source)?;
        if let Some(radiation) = &self.material_enclosure {
            return self.solve_material_enclosure(cx, radiation, problem, old, to - from, to, step, before);
        }
        let coupling = match self.coupling {
            ContactHistoryCoupling::Newton(_) => return Ok(self.contacts.advance_enthalpy_newton(
                cx, &self.run.storage, problem, old, to - from, self.newton_config_for(step)?)?),
            ContactHistoryCoupling::Legacy(coupling) => coupling,
        };
        match &self.radiation {
            None => Ok(self.contacts.advance_enthalpy(cx, &self.run.storage, problem,
                old, to - from, step, coupling)?),
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
                    problem, old, to - from, step, coupling, &patches, *config)?)
            }
        }
    }

    fn newton_config_for(&self, step: EnthalpyStepConfig)
        -> Result<ContactEnthalpyNewtonConfig, HeterogeneousEnthalpyError>
    {
        let config = self.newton_config().ok_or_else(||
            invalid("this contact history requires the coupled enthalpy Newton method"))?;
        if config.step != step {
            return Err(invalid("advance controls must match the retained contact enthalpy Newton policy"));
        }
        Ok(config)
    }
}
