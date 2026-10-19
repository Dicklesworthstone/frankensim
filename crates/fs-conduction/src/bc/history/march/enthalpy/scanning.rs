//! Moving absorbed-power footprints through the original enthalpy/radiation solves.
//!
//! This owns the scan declarations for the entire trajectory. Background source
//! histories, distinct phase charts, reference masses, k(T), phase conductivity
//! and fixed matching contacts stay with the original run. There is no alternate
//! heat equation, phase update or implicit reset of history. A source is a
//! discretely normalized deposition profile, not an optical or melt-flow model.

/// Implicit endpoint derivatives with respect to moving-source controls.
pub mod sensitivity;

use super::{EnthalpyHistoryRun, EnthalpyStepConfig, EnthalpyStepSolution,
    HeterogeneousEnthalpyError, conduction, invalid, natural_only, poll};
use super::adaptive::{PhaseErrorTolerance, phase_error};
use super::super::{events, moving::MovingHeatSource};
use super::super::moving::motion::ScanMotionLimit;
use super::super::events::adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep,
    admit_proposal, select_step};
use super::super::radiation::{RadiationRegionHistory, admit_histories, patches_before};
use crate::{AmbientRadiationConfig, ConductionProblem, ScalarField};
use crate::transient::enthalpy::radiation::EnthalpyRadiationStepSolution;
use fs_exec::Cx;
use fs_scenario::signal::TimeSignal;

/// Optional sourced ambient radiation, retained for every scan trial and resume.
pub struct ScanRadiation {
    /// Original source-backed emissivities and independent Kelvin reservoirs.
    pub patches: Vec<RadiationRegionHistory>,
    /// Existing radiative iteration controls, never relaxed by the scan driver.
    pub config: AmbientRadiationConfig,
}

/// One moving-source projection, in the original source declaration order.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanDeposition {
    /// Source center in world metres at this endpoint.
    pub center_m: [f64; 3],
    /// Declared absorbed watts for this source alone.
    pub absorbed_power_w: f64,
    /// Its independently integrated P1 field power, before superposition.
    pub delivered_power_w: f64,
    /// Its resolved positively weighted nodes; zero for an off source.
    pub active_vertices: usize,
}

/// The original physical response; radiative reports are not flattened.
#[derive(Debug, Clone)]
pub enum ScanTransport {
    /// Natural-boundary nonlinear conduction and latent storage.
    Conduction(EnthalpyStepSolution),
    /// The same storage with separate convection and physical radiation audits.
    Radiation(Box<EnthalpyRadiationStepSolution>),
}

/// Accepted endpoint plus source observations. The solver's source_w includes
/// the background AND the summed moving-source fields it actually received.
#[derive(Debug, Clone)]
pub struct ScanStep {
    /// Original solver result, including nonlinear/radiative evidence.
    pub transport: ScanTransport,
    /// Evaluated positions and source powers, without retaining redundant fields.
    pub deposition: Vec<ScanDeposition>,
}
impl ScanStep {
    /// Common physical enthalpy, temperature, phase and energy response.
    #[must_use]
    pub fn solution(&self) -> &EnthalpyStepSolution {
        match &self.transport {
            ScanTransport::Conduction(step) => step,
            ScanTransport::Radiation(step) => &step.conduction,
        }
    }
}

/// A persistent scan-driven phase-change run. No mutable escape can advance the
/// base without its attached scan/radiation. Geometry and reference mass stay fixed.
pub struct ScanningEnthalpyRun<'h, 'm, 'c> {
    run: EnthalpyHistoryRun<'h, 'm, 'c>,
    sources: Vec<MovingHeatSource>,
    radiation: Option<ScanRadiation>,
    motion_limit: Option<ScanMotionLimit>,
}
impl<'h, 'm, 'c> ScanningEnthalpyRun<'h, 'm, 'c> {
    /// Attach one or more independently moving/controlled sources before the
    /// first accepted step. Overlapping source fields add; background generation
    /// is retained. Sampling obeys each source's original spatial/power budget.
    /// The explicit source cap bounds repeated projection work.
    pub fn new(cx: &Cx<'_>, run: EnthalpyHistoryRun<'h, 'm, 'c>,
        sources: Vec<MovingHeatSource>, max_sources: usize, radiation: Option<ScanRadiation>)
        -> Result<Self, HeterogeneousEnthalpyError>
    {
        poll(cx, 0).map_err(conduction)?;
        if run.accepted_steps() != 0 || sources.is_empty() || sources.len() > max_sources {
            return Err(invalid("scanning requires an unadvanced run and a nonempty source set within its explicit cap"));
        }
        sources.len().checked_mul(run.problem.boundary.mesh().vertex_count())
            .and_then(|v| v.checked_add(sources.len().checked_mul(run.problem.boundary.mesh().element_count())?))
            .ok_or_else(|| invalid("combined moving-source work count overflows"))?;
        if let Some(radiation) = &radiation {
            radiation.config.validate().map_err(conduction)?;
            admit_histories(cx, run.problem.boundary, run.time_s(), &radiation.patches).map_err(conduction)?;
        }
        let model = Self { run, sources, radiation, motion_limit: None };
        model.source_at(cx, model.run.time_s(), false)?;
        poll(cx, 0).map_err(conduction)?;
        Ok(model)
    }

    /// Accepted state, immutable phase charts, fixed interfaces and original clock.
    #[must_use]
    pub fn conduction(&self) -> &EnthalpyHistoryRun<'h, 'm, 'c> { &self.run }

    /// Original source declarations, unchanged by all trial solves.
    #[must_use]
    pub fn sources(&self) -> &[MovingHeatSource] { &self.sources }

    /// Attach a footprint-relative motion limit before the first accepted step.
    /// This independently bounds source travel even when thermal step doubling
    /// sees identical endpoint loads. All moving sources participate, including
    /// off-source travel. Refusal changes neither policy nor caloric state.
    ///
    /// Event and adaptive steps shorten automatically. Explicit `advance_to`
    /// instead refuses an interval exceeding the bound or crossing a path knot;
    /// it never silently changes the requested endpoint or its right-side rule.
    /// The limit cannot be changed or removed after stepping begins.
    pub fn set_motion_limit(&mut self, cx: &Cx<'_>, policy: ScanMotionLimit)
        -> Result<(), HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        if self.run.accepted_steps() != 0 {
            return Err(invalid("scan motion policy must be selected before the first accepted step"));
        }
        policy.validate().map_err(conduction)?;
        for source in &self.sources { source.admit_motion_limit(cx, policy).map_err(conduction)?; }
        poll(cx, self.sources.len()).map_err(conduction)?;
        self.motion_limit = Some(policy);
        Ok(())
    }

    /// Selected geometric motion resolution; None preserves caller-owned steps.
    #[must_use]
    pub const fn motion_limit(&self) -> Option<ScanMotionLimit> { self.motion_limit }

    /// Increase the physical step allowance without changing scan or phase history.
    pub fn extend_step_budget(&mut self, maximum: usize) -> Result<(), HeterogeneousEnthalpyError> {
        self.run.extend_step_budget(maximum)
    }

    /// Sample all data at an explicitly requested right-continuous endpoint.
    /// Original material/solver/radiation refusals publish no state or time.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time: f64, config: EnthalpyStepConfig)
        -> Result<ScanStep, HeterogeneousEnthalpyError>
    {
        self.admit(cx, time)?;
        if self.motion_endpoint(cx, self.run.time_s(), time)? < time {
            return Err(invalid("explicit scan interval exceeds the motion limit; shorten it or use event/adaptive stepping"));
        }
        let step = self.solve(cx, &self.run.specific_h, self.run.time_s(), time, false, config)?;
        self.run.publish(cx, time, step.solution())?;
        Ok(step)
    }

    /// Advance one interval capped by all background, boundary, scan-position,
    /// scan-power and reservoir knots. Held data use exact end-left values.
    /// A selected motion limit also caps travel relative to every footprint.
    /// Without it, the caller's maximum time step must resolve continuous motion.
    pub fn advance_event_step(&mut self, cx: &Cx<'_>, target: f64,
        max_step: f64, config: EnthalpyStepConfig)
        -> Result<events::HistoryEventStep<ScanStep>, HeterogeneousEnthalpyError>
    {
        self.admit(cx, target)?;
        let start = self.run.time_s();
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.signals(), start, target, max_step).map_err(conduction)?;
        let end = self.motion_endpoint(cx, start, end)?;
        let report = self.solve(cx, &self.run.specific_h, start, end, true, config)?;
        self.run.publish(cx, end, report.solution())?;
        Ok(events::HistoryEventStep { start_time_s: start, end_time_s: end, report })
    }

    /// Use the original event-limited step-doubling controller and all three
    /// phase-error coordinates. Each coarse/half solve resamples the moving
    /// field at its own endpoint. Publish only the actual two-half-step path,
    /// atomically, reserving two lifetime slots. This is a local estimate, not a
    /// solution-error bound. A selected motion limit restricts the initial trial
    /// independently of its endpoint discrepancy. Like a nearby event, this may
    /// fall below the adaptive minimum; the original error test still must pass.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target: f64, proposed: f64,
        step: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
        -> Result<AdaptiveHistoryStep<ScanStep>, HeterogeneousEnthalpyError>
    {
        self.admit(cx, target)?;
        phase.validate()?;
        let cap = admit_proposal(cx, self.run.accepted_steps(), self.run.max_steps(),
            proposed, adaptive).map_err(conduction)?;
        let start = self.run.time_s();
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.signals(), start, target, cap).map_err(conduction)?;
        let end = self.motion_endpoint(cx, start, end)?;
        let selected = select_step(cx, start, end, adaptive, conduction, |mid, end| {
            let coarse = self.solve(cx, &self.run.specific_h, start, end, true, step)?;
            let first = self.solve(cx, &self.run.specific_h, start, mid, true, step)?;
            let second = self.solve(cx, &first.solution().specific_enthalpy_j_kg, mid, end, true, step)?;
            let error = phase_error(cx, &self.run, coarse.solution(), second.solution(), adaptive, phase)?;
            Ok((first, second, error))
        })?;
        self.run.publish(cx, selected.end_time_s, selected.second_half.solution())?;
        self.run.accepted_steps += 1;
        Ok(selected)
    }

    fn motion_endpoint(&self, cx: &Cx<'_>, start: f64, mut end: f64)
        -> Result<f64, HeterogeneousEnthalpyError>
    {
        if let Some(policy) = self.motion_limit {
            for source in &self.sources {
                end = source.motion_limited_endpoint(cx, start, end, policy).map_err(conduction)?;
            }
        }
        Ok(end)
    }

    fn admit(&self, cx: &Cx<'_>, target: f64) -> Result<(), HeterogeneousEnthalpyError> {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        if self.run.accepted_steps() >= self.run.max_steps()
            || !target.is_finite() || !(target - self.run.time_s()).is_finite()
            || target <= self.run.time_s()
        { return Err(invalid("scan endpoint requires remaining physical steps and a finite positive time interval")); }
        Ok(())
    }

    fn signals(&self) -> impl Iterator<Item = &TimeSignal> {
        self.sources.iter().flat_map(|source| source.signals()).chain(
            self.radiation.iter().flat_map(|r| r.patches.iter().map(|p| &p.ambient_temperature)))
    }

    fn source_at(&self, cx: &Cx<'_>, time: f64, left: bool)
        -> Result<(ScalarField, Vec<ScanDeposition>), HeterogeneousEnthalpyError>
    {
        let mesh = self.run.problem.boundary.mesh();
        let background = if left {
            events::source_before(&self.run.problem.source, time, mesh.vertex_count())
        } else { self.run.problem.source.at(time, mesh.vertex_count()) }.map_err(conduction)?;
        let mut field = Vec::new();
        field.try_reserve_exact(mesh.vertex_count()).map_err(|_| invalid("scan field allocation refused"))?;
        for i in 0..mesh.vertex_count() {
            if i % 256 == 0 { poll(cx, i).map_err(conduction)?; }
            field.push(background.at(i));
        }
        let mut deposition = Vec::with_capacity(self.sources.len());
        for source in &self.sources {
            poll(cx, deposition.len()).map_err(conduction)?;
            let sample = if left { source.before(cx, mesh, time) } else { source.at(cx, mesh, time) }
                .map_err(conduction)?;
            for (i, value) in field.iter_mut().enumerate() {
                if i % 256 == 0 { poll(cx, i).map_err(conduction)?; }
                *value += sample.source.at(i);
                crate::require_finite("combined moving heat source", *value).map_err(conduction)?;
            }
            deposition.push(ScanDeposition { center_m: sample.center_m,
                absorbed_power_w: sample.absorbed_power_w, delivered_power_w: sample.delivered_power_w,
                active_vertices: sample.active_vertices });
        }
        poll(cx, field.len()).map_err(conduction)?;
        Ok((ScalarField::Nodal(field), deposition))
    }

    fn solve(&self, cx: &Cx<'_>, old: &[f64], start: f64, end: f64,
        left: bool, config: EnthalpyStepConfig) -> Result<ScanStep, HeterogeneousEnthalpyError>
    {
        let boundary = if left { events::boundary_before(cx, self.run.problem.boundary, end) }
            else { self.run.problem.boundary.at(cx, end) }.map_err(conduction)?;
        natural_only(&boundary)?;
        let (source, deposition) = self.source_at(cx, end, left)?;
        let problem = ConductionProblem { mesh: self.run.problem.boundary.mesh(), boundary: &boundary,
            material: self.run.problem.material, element_materials: self.run.problem.element_materials, source: &source };
        let transport = match &self.radiation {
            None => ScanTransport::Conduction(self.run.storage.advance(cx, problem,
                self.run.problem.interfaces, old, end - start, config)?),
            Some(radiation) => {
                let patches = if left { patches_before(cx, &radiation.patches, end) } else {
                    radiation.patches.iter().map(|p| p.at(end)).collect()
                }.map_err(conduction)?;
                ScanTransport::Radiation(Box::new(self.run.storage.advance_with_ambient_radiation(cx,
                    problem, self.run.problem.interfaces, old, end - start, config, &patches, radiation.config)?))
            }
        };
        Ok(ScanStep { transport, deposition })
    }
}
