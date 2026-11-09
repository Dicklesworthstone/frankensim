//! Scheduled second-order stages on the original primary-enthalpy history.
//!
//! This owner changes only time quadrature. Storage, source k, material contact
//! and reflected radiation retain their original spatial residuals and source
//! identities. Stage carry is internal numerical history, never external heat.
//! Optional moving absorbed-power sources retain their original normalized P1
//! projection, sampled separately at each physical stage time.

/// Fixed-grid native stage adjoints and authored moving-source controls.
pub mod sensitivity;

use super::{EnthalpyHistoryRun, EnthalpyStepConfig, HeterogeneousEnthalpyError,
    conduction, invalid, poll};
use super::adaptive::{
    EntropyAdaptiveHistoryError, EntropyAdaptiveHistoryStep, PhaseErrorTolerance, PhaseFields,
    accepted_entropy, admit_entropy_tolerance, entropy_conduction, phase_field_error,
    recoverable_history, recoverable_entropy_history,
};
use crate::transient::enthalpy::entropy::EnthalpyEntropyBalance;
use super::radiation::enclosure::MaterialEnclosureHistoryBinding;
use super::radiation::material::{MaterialAmbientHistoryBinding, MaterialRadiationRegionHistory};
use super::scanning::radiation::sample as pointwise_ambient_samples;
use super::super::events::{self, HistoryEventStep,
    adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep, admit_proposal, select_step_with_recovery}};
use super::super::moving::{MovingHeatSource, assembly, motion::ScanMotionLimit};
use super::super::radiation::RadiationRegionHistory;
use crate::interface::temperature::TemperatureContacts;
use crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use crate::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use crate::transient::enthalpy::sdirk::{
    EnthalpySdirk2Input, EnthalpySdirk2Stage, EnthalpySdirk2Step, SDIRK2_GAMMA,
};
use crate::{MaterialEnclosureRadiationPatch, MaterialPointwiseAmbientRadiation,
    MaterialPointwiseEnclosureRadiation, PointwiseAmbientRadiation, ScalarField};
use fs_exec::Cx;
use fs_scenario::signal::TimeSignal;
use std::sync::Arc;

/// Persistent SDIRK2 evolution of an already admitted physical history.
///
/// Construct with `into_sdirk2` on an original `EnthalpyHistoryRun`, coupled
/// Newton `ContactEnthalpyHistoryRun`, `MaterialRadiativeEnthalpyHistoryRun`,
/// `MaterialEnclosureEnthalpyHistoryRun`, or an admitted pointwise-radiative
/// `ScanningEnthalpyRun`.
/// Conversion preserves accepted fields, time, source declarations and budget.
/// Only the stiffly accurate second stage is published. Internal stages do not
/// consume accepted-step slots or appear as independently accepted timesteps.
pub struct Sdirk2EnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    run: EnthalpyHistoryRun<'h, 'm, 'c>,
    contacts: Option<(TemperatureContacts<'m, 'i>, ContactEnthalpyNewtonConfig)>,
    radiation: Option<MaterialRadiationHistory<'h, 'c>>,
    moving_sources: Arc<[MovingHeatSource]>,
    motion_limit: Option<ScanMotionLimit>,
}

enum MaterialRadiationHistory<'h, 'c> {
    PointwiseAmbient(Vec<RadiationRegionHistory>),
    Ambient(MaterialAmbientHistoryBinding<'c>),
    Enclosure(MaterialEnclosureHistoryBinding<'h, 'c>),
}

impl<'h, 'm, 'c, 'i> Sdirk2EnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    pub(super) fn from_parts(run: EnthalpyHistoryRun<'h, 'm, 'c>,
        contacts: Option<(TemperatureContacts<'m, 'i>, ContactEnthalpyNewtonConfig)>,
        radiation: Option<MaterialEnclosureHistoryBinding<'h, 'c>>) -> Self
    { Self { run, contacts, radiation: radiation.map(MaterialRadiationHistory::Enclosure),
        moving_sources: Arc::from([]), motion_limit: None } }

    // Transfer an already admitted scan, including an evolved accepted state.
    // The public attachment API correctly refuses an already evolved history.
    pub(super) fn from_scanning(
        run: EnthalpyHistoryRun<'h, 'm, 'c>,
        contacts: Option<(TemperatureContacts<'m, 'i>, ContactEnthalpyNewtonConfig)>,
        moving_sources: Vec<MovingHeatSource>, motion_limit: Option<ScanMotionLimit>,
        pointwise_ambient: Option<Vec<RadiationRegionHistory>>,
    ) -> Self {
        Self { run, contacts,
            radiation: pointwise_ambient.map(MaterialRadiationHistory::PointwiseAmbient),
            moving_sources: moving_sources.into(), motion_limit }
    }

    pub(super) fn from_material_ambient(
        run: EnthalpyHistoryRun<'h, 'm, 'c>, radiation: MaterialAmbientHistoryBinding<'c>,
    ) -> Self {
        Self { run, contacts: None, radiation: Some(MaterialRadiationHistory::Ambient(radiation)),
            moving_sources: Arc::from([]), motion_limit: None }
    }

    /// Attach independently normalized moving absorbed-power sources before the
    /// first accepted interval. Original background generation, phase materials,
    /// contact laws and radiation remain active. Every SDIRK stage samples its
    /// own physical position and power; no endpoint field is reused across stages.
    /// Native stage/interval source accounts include the complete summed field.
    ///
    /// # Errors
    /// Advanced history, duplicate attachment, source/work-domain refusal or
    /// cancellation. Attachment takes no physical step and changes no caloric state.
    pub fn with_moving_sources(
        mut self, cx: &Cx<'_>, sources: Vec<MovingHeatSource>, max_sources: usize,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0).map_err(conduction)?;
        if self.run.accepted_steps != 0 || !self.moving_sources.is_empty() {
            return Err(invalid("attach moving sources once, before the first SDIRK2 history interval"));
        }
        let mesh = self.run.problem.boundary.mesh();
        assembly::admit(cx, mesh, &sources, max_sources).map_err(conduction)?;
        let background = self.run.problem.source.at(self.run.time_s, mesh.vertex_count())
            .map_err(conduction)?;
        assembly::project(cx, mesh, &background, &sources, self.run.time_s, false, |_| {})
            .map_err(conduction)?;
        poll(cx, sources.len()).map_err(conduction)?;
        self.moving_sources = sources.into();
        Ok(self)
    }

    /// Original source declarations, immutable across stage trials and resumes.
    #[must_use]
    pub fn moving_sources(&self) -> &[MovingHeatSource] { &self.moving_sources }

    /// Resolve moving footprints independently of the thermal error estimate.
    /// Select before the first interval. All sources participate, including
    /// zero-power travel. Event/adaptive intervals shorten; direct endpoints
    /// refuse if they exceed the motion bound or cross an unresolved path knot.
    pub fn set_motion_limit(&mut self, cx: &Cx<'_>, policy: ScanMotionLimit)
        -> Result<(), HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps).map_err(conduction)?;
        if self.run.accepted_steps != 0 || self.moving_sources.is_empty() {
            return Err(invalid("SDIRK2 motion policy requires attached sources and an unadvanced history"));
        }
        assembly::admit_motion(cx, &self.moving_sources, policy).map_err(conduction)?;
        self.motion_limit = Some(policy);
        Ok(())
    }

    /// Selected footprint-relative motion resolution, when explicitly supplied.
    #[must_use]
    pub const fn motion_limit(&self) -> Option<ScanMotionLimit> { self.motion_limit }

    /// Retain fixed source-backed emissivities and independent black-reservoir
    /// schedules before the first interval. Each physical stage integrates the
    /// original local P1 T^4 law on its actual boundary, adding to its existing
    /// natural heat flux. Reservoir events use the same end-left sampling as
    /// sources; stage radiation remains separate from internal contact heat.
    ///
    /// # Errors
    /// Advanced history, another radiation owner, invalid patch/bath/trace or
    /// initial surface-material support, and cancellation. No step is taken.
    pub fn with_pointwise_ambient(
        mut self, cx: &Cx<'_>, patches: Vec<RadiationRegionHistory>,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0).map_err(conduction)?;
        if self.run.accepted_steps != 0 || self.radiation.is_some() {
            return Err(invalid("attach exclusive pointwise ambient radiation before the first SDIRK2 interval"));
        }
        super::scanning::radiation::admit_pointwise_ambient(cx, &self.run, &patches)?;
        self.radiation = Some(MaterialRadiationHistory::PointwiseAmbient(patches));
        Ok(self)
    }

    /// Original fixed emissivity declarations and black-reservoir schedules.
    #[must_use]
    pub fn pointwise_ambient_patches(&self) -> Option<&[RadiationRegionHistory]> {
        match &self.radiation {
            Some(MaterialRadiationHistory::PointwiseAmbient(patches)) => Some(patches),
            _ => None,
        }
    }

    /// Attach independent material black reservoirs to an unadvanced original
    /// history, including a retained native phase-aware contact history.
    ///
    /// Every stage evaluates emissivity from its own original h/T/f and samples
    /// each bath at that stage's actual time. Bath knots join source/boundary
    /// events. Radiation is counted once, separately from internal contact,
    /// source heat and convection; existing Newton/energy/entropy gates apply.
    ///
    /// # Errors
    /// Advanced history, another radiation owner, invalid source/bath/trace,
    /// nonidentical caloric charts, or cancellation. No physical step is taken.
    pub fn with_material_ambient(
        mut self, cx: &Cx<'_>, patches: Vec<MaterialRadiationRegionHistory<'c>>,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0).map_err(conduction)?;
        if self.radiation.is_some() {
            return Err(invalid("SDIRK2 history requires exclusive ambient or enclosure radiation ownership"));
        }
        let radiation = MaterialAmbientHistoryBinding::new(cx, &self.run, patches)?;
        self.radiation = Some(MaterialRadiationHistory::Ambient(radiation));
        Ok(self)
    }

    /// Original material declarations and per-surface black-reservoir histories.
    #[must_use]
    pub fn material_ambient_patches(&self) -> Option<&[MaterialRadiationRegionHistory<'c>]> {
        match &self.radiation {
            Some(MaterialRadiationHistory::Ambient(binding)) => Some(binding.patches()),
            _ => None,
        }
    }

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
        match &self.radiation {
            Some(MaterialRadiationHistory::Enclosure(binding)) => Some(binding.patches()),
            _ => None,
        }
    }

    /// Sole external black reservoir; absent for a closed or missing enclosure.
    #[must_use]
    pub fn material_enclosure_ambient_temperature(&self) -> Option<&TimeSignal> {
        match &self.radiation {
            Some(MaterialRadiationHistory::Enclosure(binding)) => binding.ambient_temperature(),
            _ => None,
        }
    }

    /// Extend the original physical-interval budget without resetting history.
    pub fn extend_step_budget(&mut self, maximum: usize) -> Result<(), HeterogeneousEnthalpyError> {
        self.run.extend_step_budget(maximum)
    }

    /// Advance one complete SDIRK2 interval, sampling authored schedules at
    /// `start + gamma*dt` and the requested endpoint. Both samples use right
    /// limits. Caller-supplied intervals are not automatically event-split.
    /// An attached motion policy refuses intervals exceeding source resolution.
    /// Source, nonlinear, heat-balance or cancellation refusal publishes nothing.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time_s: f64, config: EnthalpyStepConfig)
        -> Result<EnthalpySdirk2Step, HeterogeneousEnthalpyError>
    {
        self.admit_one(cx)?;
        if self.motion_endpoint(cx, self.run.time_s, time_s)? < time_s {
            return Err(invalid("explicit SDIRK2 interval exceeds the moving-source resolution; shorten it or use event/adaptive stepping"));
        }
        let solved = self.solve(cx, &self.run.specific_h, self.run.time_s, time_s, config, false)?;
        self.publish(cx, time_s, &solved.second_stage)?;
        Ok(solved)
    }

    /// End at the next source, scan-position, boundary or radiative-bath knot,
    /// additionally capped by any attached footprint-relative motion policy.
    /// The final stage uses
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
    /// A moving-source cap can shorten a trial below the nominal minimum, like
    /// an event; the same phase-error acceptance still applies to that interval.
    /// This is local control, not a global or phase-transition-time certificate.
    /// `endpoint_solves` counts SDIRK2 interval solves; each has two implicit
    /// stage solves. Native nonlinear nonconvergence and residual/energy
    /// refusals bisect within the same bounded controller. Stage-domain,
    /// material, budget and cancellation refusals propagate immediately.
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
        let selected = select_step_with_recovery(cx, start, end, adaptive, conduction,
            recoverable_history, |middle, end, solves| {
            *solves += 1;
            let coarse = self.solve(cx, &self.run.specific_h, start, end, config, true)?;
            *solves += 1;
            let first = self.solve(cx, &self.run.specific_h, start, middle, config, true)?;
            *solves += 1;
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

    /// Adapt with native entropy accounts on the actual two SDIRK2 fine
    /// intervals. The energy-only coarse comparison is never entropy-admitted.
    ///
    /// Each half retains its authentic two-stage quadrature and numerical
    /// defect. The accepted pair also passes one whole-window entropy allowance,
    /// using exact original-to-final storage. No discarded solve, internal carry
    /// or stage is counted as an independently accepted physical interval.
    ///
    /// Schedule sampling, material/contact identity and unsupported radiation
    /// refusals are shared with the native solver. Recoverable physical
    /// failures bisect as on the energy-only path; entropy refusals propagate.
    /// Any final refusal preserves all accepted fields, time and counters.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step_with_entropy(
        &mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        config: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig,
        phase: PhaseErrorTolerance, entropy_tolerance_j_k: f64,
    ) -> Result<EntropyAdaptiveHistoryStep<EnthalpySdirk2Step>, EntropyAdaptiveHistoryError> {
        poll(cx, self.run.accepted_steps).map_err(entropy_conduction)?;
        phase.validate()?;
        admit_entropy_tolerance(entropy_tolerance_j_k)?;
        let start = self.run.time_s;
        let cap = admit_proposal(cx, self.run.accepted_steps, self.run.max_steps,
            proposed_step_s, adaptive).map_err(entropy_conduction)?;
        let end = self.next_endpoint(cx, target_s, cap)?;
        let selected = select_step_with_recovery(cx, start, end, adaptive, entropy_conduction,
            recoverable_entropy_history, |middle, end, solves| -> Result<_, EntropyAdaptiveHistoryError> {
                *solves += 1;
                let coarse = self.solve(cx, &self.run.specific_h, start, end, config, true)?;
                *solves += 1;
                let first = self.solve_with_entropy(cx, &self.run.specific_h, start, middle,
                    config, true, entropy_tolerance_j_k)?;
                *solves += 1;
                let second = self.solve_with_entropy(cx, &first.0.second_stage.specific_enthalpy_j_kg,
                    middle, end, config, true, entropy_tolerance_j_k)?;
                let error = phase_field_error(cx, &self.run,
                    fields(&coarse.second_stage), fields(&second.0.second_stage), adaptive, phase)?;
                Ok((first, second, error))
            })?;
        let entropy = accepted_entropy(cx, &self.run,
            &selected.second_half.0.second_stage.specific_enthalpy_j_kg,
            [selected.first_half.1, selected.second_half.1], entropy_tolerance_j_k)?;
        let selected = EntropyAdaptiveHistoryStep::from_halves(selected, entropy);
        self.publish(cx, selected.physical.end_time_s, &selected.physical.second_half.second_stage)?;
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
        let reservoirs = self.material_enclosure_ambient_temperature().into_iter()
            .chain(self.material_ambient_patches().into_iter().flatten()
                .map(|patch| &patch.ambient_temperature))
            .chain(self.pointwise_ambient_patches().into_iter().flatten()
                .map(|patch| &patch.ambient_temperature))
            .chain(self.moving_sources.iter().flat_map(MovingHeatSource::signals));
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            reservoirs, self.run.time_s, target, cap).map_err(conduction)?;
        self.motion_endpoint(cx, self.run.time_s, end)
    }

    fn motion_endpoint(&self, cx: &Cx<'_>, start: f64, end: f64)
        -> Result<f64, HeterogeneousEnthalpyError>
    {
        assembly::motion_endpoint(cx, &self.moving_sources, self.motion_limit, start, end)
            .map_err(conduction)
    }

    fn source_at(&self, cx: &Cx<'_>, time_s: f64, before: bool)
        -> Result<ScalarField, HeterogeneousEnthalpyError>
    {
        let mesh = self.run.problem.boundary.mesh();
        let background = if before {
            events::source_before(&self.run.problem.source, time_s, mesh.vertex_count())
        } else { self.run.problem.source.at(time_s, mesh.vertex_count()) }.map_err(conduction)?;
        if self.moving_sources.is_empty() { return Ok(background); }
        assembly::project(cx, mesh, &background, &self.moving_sources, time_s, before, |_| {})
            .map_err(conduction)
    }

    #[allow(clippy::too_many_arguments)]
    fn solve(&self, cx: &Cx<'_>, old: &[f64], from: f64, to: f64,
        config: EnthalpyStepConfig, before: bool) -> Result<EnthalpySdirk2Step, HeterogeneousEnthalpyError>
    {
        let solved = self.with_stage_inputs(cx, from, to, before,
            |inputs, dt| self.run.storage.advance_sdirk2(cx, inputs, old, dt, config))?;
        Ok(solved?)
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_with_entropy(
        &self, cx: &Cx<'_>, old: &[f64], from: f64, to: f64,
        config: EnthalpyStepConfig, before: bool, entropy_tolerance_j_k: f64,
    ) -> Result<(EnthalpySdirk2Step, EnthalpyEntropyBalance), EntropyAdaptiveHistoryError> {
        let solved = self.with_stage_inputs(cx, from, to, before,
            |inputs, dt| self.run.storage.advance_sdirk2_with_entropy(
                cx, inputs, old, dt, config, entropy_tolerance_j_k))?;
        Ok(solved?)
    }

    /// One schedule sampler supplies either native finish mode. Every local
    /// boundary/radiation owner stays alive throughout its actual stage solves.
    fn with_stage_inputs<R>(
        &self, cx: &Cx<'_>, from: f64, to: f64, before: bool,
        solve: impl for<'a> FnOnce([EnthalpySdirk2Input<'a>; 2], f64) -> R,
    ) -> Result<R, HeterogeneousEnthalpyError> {
        let (dt, stage_time) = stage_clock(from, to)?;
        let history = self.run.problem.boundary;
        let first_boundary = history.at(cx, stage_time).map_err(conduction)?;
        let second_boundary = if before { events::boundary_before(cx, history, to) }
            else { history.at(cx, to) }.map_err(conduction)?;
        let first_source = self.source_at(cx, stage_time, false)?;
        let second_source = self.source_at(cx, to, before)?;
        let first_problem = self.run.conduction_problem(&first_boundary, &first_source)?;
        let second_problem = self.run.conduction_problem(&second_boundary, &second_source)?;
        let contacts = self.contacts.as_ref().map(|(owner, config)| (owner, *config));
        let inputs = |radiation| EnthalpySdirk2Input { problem: first_problem,
            interfaces: self.run.problem.interfaces, radiation, contacts };
        let endpoint = |radiation| EnthalpySdirk2Input { problem: second_problem,
            interfaces: self.run.problem.interfaces, radiation, contacts };
        match &self.radiation {
            None => Ok(solve(
                [inputs(EnthalpyRadiation::None), endpoint(EnthalpyRadiation::None)], dt)),
            Some(MaterialRadiationHistory::PointwiseAmbient(patches)) => {
                let first_patches = pointwise_ambient_samples(cx, patches, stage_time, false)
                    .map_err(conduction)?;
                let second_patches = pointwise_ambient_samples(cx, patches, to, before)
                    .map_err(conduction)?;
                let first = PointwiseAmbientRadiation::new(
                    cx, history.mesh(), &first_boundary, &first_patches).map_err(conduction)?;
                let second = PointwiseAmbientRadiation::new(
                    cx, history.mesh(), &second_boundary, &second_patches).map_err(conduction)?;
                Ok(solve(
                    [inputs(EnthalpyRadiation::PointwiseAmbient(&first)),
                     endpoint(EnthalpyRadiation::PointwiseAmbient(&second))], dt))
            }
            Some(MaterialRadiationHistory::Ambient(binding)) => {
                let first_patches = binding.sample(cx, stage_time, false).map_err(conduction)?;
                let second_patches = binding.sample(cx, to, before).map_err(conduction)?;
                let first = MaterialPointwiseAmbientRadiation::new(
                    cx, history.mesh(), &first_boundary, &first_patches).map_err(conduction)?;
                let second = MaterialPointwiseAmbientRadiation::new(
                    cx, history.mesh(), &second_boundary, &second_patches).map_err(conduction)?;
                Ok(solve(
                    [inputs(EnthalpyRadiation::MaterialAmbient(&first)),
                     endpoint(EnthalpyRadiation::MaterialAmbient(&second))], dt))
            }
            Some(MaterialRadiationHistory::Enclosure(binding)) => {
                let first_base = binding.sample(cx, history.mesh(), &first_boundary, stage_time, false)?;
                let second_base = binding.sample(cx, history.mesh(), &second_boundary, to, before)?;
                let first = MaterialPointwiseEnclosureRadiation::new(cx, &first_base, binding.patches()).map_err(conduction)?;
                let second = MaterialPointwiseEnclosureRadiation::new(cx, &second_base, binding.patches()).map_err(conduction)?;
                Ok(solve(
                    [inputs(EnthalpyRadiation::MaterialEnclosure(&first)),
                     endpoint(EnthalpyRadiation::MaterialEnclosure(&second))], dt))
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

fn stage_clock(from: f64, to: f64) -> Result<(f64, f64), HeterogeneousEnthalpyError> {
    let dt = to - from;
    let stage_time = SDIRK2_GAMMA.mul_add(dt, from);
    if ![from, to, dt, stage_time].iter().all(|x| x.is_finite()) || dt <= 0.0
        || !(from < stage_time && stage_time < to)
    {
        return Err(invalid("SDIRK2 history needs a finite positive interval and representable interior stage time"));
    }
    Ok((dt, stage_time))
}

fn fields(stage: &EnthalpySdirk2Stage) -> PhaseFields<'_> {
    PhaseFields { h: &stage.specific_enthalpy_j_kg, temperature: &stage.temperature,
        fraction: &stage.liquid_mass_fraction }
}
