//! Native SDIRK2 stage adjoints under the original source and boundary histories.
//!
//! The storage owner solves the physical interval once and retains both stage
//! tangents, including phase contacts and pointwise material radiation. Moving
//! sources use those same stage source-density cotangents. Source and prescribed
//! bath schedules are contracted at the original tableau times and event side;
//! their contributions already contain every tableau factor and are summed once.
//! Parameter-aware pullbacks additionally retain the original per-stage contact
//! and radiation controls. Geometry, caloric charts and the selected time grid
//! stay fixed. Neither an adaptive decision nor the event that selected an
//! endpoint is differentiated. Source declarations are shared immutably with
//! the original run.

use super::{EnthalpyStepConfig, Sdirk2EnthalpyHistoryRun, poll, stage_clock};
use super::super::scanning::sensitivity::{ScanSensitivityError, fraction};
use super::super::scanning::sensitivity::objective::{ScanPrescribedTemperatureGradient, baths};
use crate::bc::history::ThermalBoundaryHistory;
use crate::bc::history::march::moving::MovingHeatSource;
use crate::bc::history::march::moving::sensitivity::{MovingHeatGradient, MovingHeatLinearization};
use crate::bc::history::march::moving::sensitivity::schedule::MovingHeatScheduleGradient;
use crate::transient::enthalpy::adjoint::EnthalpyAdjointError;
use crate::transient::enthalpy::adjoint::sdirk::{EnthalpySdirk2Gradient, EnthalpySdirk2Linearization,
    EnthalpySdirk2ParameterBudget, EnthalpySdirk2ParameterGradient, EnthalpySdirk2StageParameters};
use crate::transient::enthalpy::sdirk::EnthalpySdirk2Step;
use crate::LinearConfig;
use fs_exec::Cx;
use std::sync::Arc;

/// Additional retention allowances for one actual two-stage interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sdirk2SensitivityBudget {
    /// Two stages times sources times mesh vertices. Each source binding retains
    /// up to three scalar arrays per vertex; this is a shape count, not an RSS cap.
    pub max_stage_source_vertices: usize,
    /// Sum of the original radiation-adjoint workspace allowances at both stages.
    /// Zero is valid when no radiation is attached.
    pub max_radiation_entries: usize,
}

/// Complete native history/source pullback; arrays use first, second stage order.
#[derive(Debug)]
pub struct Sdirk2HistoryGradient {
    /// Original previous-h seed, physical source/bath cotangents and checked
    /// transpose work. Carry `previous_specific_enthalpy` into the prior interval.
    pub transport: EnthalpySdirk2Gradient,
    /// Independent sampled center, width and absorbed-power controls at each
    /// stage, in the original source declaration order.
    pub sources: [Vec<MovingHeatGradient>; 2],
}

/// One interval's authored controls, summed over its two actual physical stages.
#[derive(Debug)]
pub struct Sdirk2HistoryScheduleGradient {
    /// Complete native and sampled physical gradient; no transpose is discarded.
    pub physical: Sdirk2HistoryGradient,
    /// Shared source values/coefficients and widths, in declaration order.
    pub sources: Vec<MovingHeatScheduleGradient>,
    /// Original authored Dirichlet rows, retaining boundary ownership/precedence.
    /// Static partition values remain in the native nodal gradient, not new rows.
    pub prescribed_temperatures: Vec<ScanPrescribedTemperatureGradient>,
}

/// The same history/source transpose with original contact and radiation controls.
#[derive(Debug)]
pub struct Sdirk2HistoryParameterGradient {
    pub physical: Sdirk2HistoryGradient,
    /// Original independent stage laws and sampled radiative reservoir controls.
    /// These are already weighted physical contributions; add matching controls
    /// only when their original declarations are shared between stages.
    pub stages: [EnthalpySdirk2StageParameters; 2],
}

/// Shared source/bath schedules and independently retained physical stage controls.
#[derive(Debug)]
pub struct Sdirk2HistoryParameterScheduleGradient {
    pub physical: Sdirk2HistoryScheduleGradient,
    /// Retains material/radiation identity without inferring an authored radiative
    /// schedule from the stage values. No second tableau factor is needed.
    pub stages: [EnthalpySdirk2StageParameters; 2],
}

/// Immutable actual two-stage interval, with independently bound source controls.
///
/// The boundary borrow is the original external history, not a borrow of the
/// mutable run. Retaining this binding therefore permits subsequent advances.
/// Source schedules use a shared immutable allocation rather than one table copy
/// per retained interval. No trajectory tape is allocated by this object.
#[derive(Debug)]
pub struct Sdirk2HistoryLinearization<'h, 'm> {
    transport: EnthalpySdirk2Linearization<'m>,
    sources: [Vec<MovingHeatLinearization<'m>>; 2],
    declarations: Arc<[MovingHeatSource]>,
    boundary: &'h ThermalBoundaryHistory<'m>,
    fractions: [Vec<fraction::FractionDerivative>; 2],
    start_time_s: f64,
    stage_times_s: [f64; 2],
    end_before: bool,
}

impl Sdirk2HistoryLinearization<'_, '_> {
    /// Original interval report with both physical stages and integrated heat.
    #[must_use]
    pub fn primal(&self) -> &EnthalpySdirk2Step { self.transport.primal() }

    /// Complete native tangent, including contact and material-radiation feedback.
    #[must_use]
    pub const fn transport(&self) -> &EnthalpySdirk2Linearization<'_> { &self.transport }

    #[must_use]
    pub const fn start_time_s(&self) -> f64 { self.start_time_s }
    #[must_use]
    pub const fn end_time_s(&self) -> f64 { self.stage_times_s[1] }
    #[must_use]
    pub const fn stage_times_s(&self) -> [f64; 2] { self.stage_times_s }
    /// Whether final-stage data use the exact interval-end left limit.
    /// The interior first stage always uses its right-continuous value.
    #[must_use]
    pub const fn end_before(&self) -> bool { self.end_before }

    /// Endpoint enthalpy objective or reverse-history carry.
    pub fn pullback(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<Sdirk2HistoryGradient, ScanSensitivityError>
    {
        self.with_sources(cx, self.transport.pullback(cx, seed, linear)?)
    }

    /// Independent h objectives at both authentic stages. The native eliminated
    /// block transpose includes the complete interstage storage carry.
    pub fn pullback_stages(&self, cx: &Cx<'_>, seeds: [&[f64]; 2], linear: LinearConfig)
        -> Result<Sdirk2HistoryGradient, ScanSensitivityError>
    {
        self.with_sources(cx, self.transport.pullback_stages(cx, seeds, linear)?)
    }

    /// Endpoint objective with native contact/radiation controls. Uses the same
    /// two thermal transposes as the source-only entry, without another solve.
    pub fn pullback_with_parameters(&self, cx: &Cx<'_>, seed: &[f64],
        linear: LinearConfig, budget: EnthalpySdirk2ParameterBudget)
        -> Result<Sdirk2HistoryParameterGradient, ScanSensitivityError>
    {
        self.with_parameters(cx, self.transport.pullback_with_parameters(cx, seed, linear, budget)?)
    }

    /// Both stage objectives through the same native two-stage block transpose.
    pub fn pullback_stages_with_parameters(&self, cx: &Cx<'_>, seeds: [&[f64]; 2],
        linear: LinearConfig, budget: EnthalpySdirk2ParameterBudget)
        -> Result<Sdirk2HistoryParameterGradient, ScanSensitivityError>
    {
        self.with_parameters(cx, self.transport.pullback_stages_with_parameters(cx, seeds, linear, budget)?)
    }

    /// Convert independent temperature observations through each original T(h).
    pub fn temperature_pullback(&self, cx: &Cx<'_>, seeds: [&[f64]; 2])
        -> Result<[Vec<f64>; 2], ScanSensitivityError>
    { Ok(self.transport.temperature_pullback(cx, seeds)?) }

    /// Convert only an accepted endpoint temperature observation to h coordinates.
    pub fn endpoint_temperature_pullback(&self, cx: &Cx<'_>, seed: &[f64])
        -> Result<Vec<f64>, ScanSensitivityError>
    {
        let zero = zeros(self.primal().first_stage.specific_enthalpy_j_kg.len())?;
        let [_, endpoint] = self.temperature_pullback(cx, [&zero, seed])?;
        Ok(endpoint)
    }

    /// Convert stage liquid-mass-fraction observations through each original
    /// caloric chart. A seeded fraction-only corner refuses without inventing a
    /// classical derivative; an unobserved corner does not disable h/T objectives.
    pub fn liquid_mass_fraction_pullback(&self, cx: &Cx<'_>, seeds: [&[f64]; 2])
        -> Result<[Vec<f64>; 2], ScanSensitivityError>
    {
        Ok([fraction::pullback(cx, &self.fractions[0],
            &self.primal().first_stage.specific_enthalpy_j_kg, seeds[0])?,
            self.endpoint_liquid_mass_fraction_pullback(cx, seeds[1])?])
    }

    /// Convert an accepted endpoint liquid-fraction observation to h coordinates.
    pub fn endpoint_liquid_mass_fraction_pullback(&self, cx: &Cx<'_>, seed: &[f64])
        -> Result<Vec<f64>, ScanSensitivityError>
    {
        fraction::pullback(cx, &self.fractions[1],
            &self.primal().second_stage.specific_enthalpy_j_kg, seed)
    }

    /// Full endpoint-temperature objective through the same native two-stage dual.
    pub fn pullback_temperature(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<Sdirk2HistoryGradient, ScanSensitivityError>
    {
        self.pullback(cx, &self.endpoint_temperature_pullback(cx, seed)?, linear)
    }

    /// Full endpoint liquid-fraction objective, including latent contact feedback.
    pub fn pullback_liquid_mass_fraction(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<Sdirk2HistoryGradient, ScanSensitivityError>
    {
        self.pullback(cx, &self.endpoint_liquid_mass_fraction_pullback(cx, seed)?, linear)
    }

    /// Endpoint objective contracted into this binding's ORIGINAL schedules.
    /// `max_controls` counts shared authored source/bath coordinates once, not
    /// once per stage. Unsupported chirp control maps refuse before transposes;
    /// instantaneous source pullbacks remain available for those trajectories.
    pub fn pullback_schedule(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig,
        max_controls: usize) -> Result<Sdirk2HistoryScheduleGradient, ScanSensitivityError>
    {
        self.admit_schedules(cx, max_controls)?;
        self.with_schedules(cx, self.pullback(cx, seed, linear)?, max_controls)
    }

    /// Both stage objectives contracted at their actual sampled times/sides.
    /// The complete native transpose already includes every tableau coefficient;
    /// summing these controls adds no second time or integration weight.
    pub fn pullback_schedule_stages(&self, cx: &Cx<'_>, seeds: [&[f64]; 2],
        linear: LinearConfig, max_controls: usize)
        -> Result<Sdirk2HistoryScheduleGradient, ScanSensitivityError>
    {
        self.admit_schedules(cx, max_controls)?;
        self.with_schedules(cx, self.pullback_stages(cx, seeds, linear)?, max_controls)
    }

    /// Endpoint source/bath schedules with original contact/radiation controls.
    /// All outputs derive from the same admitted native thermal transpose.
    pub fn pullback_schedule_with_parameters(&self, cx: &Cx<'_>, seed: &[f64],
        linear: LinearConfig, budget: EnthalpySdirk2ParameterBudget, max_controls: usize)
        -> Result<Sdirk2HistoryParameterScheduleGradient, ScanSensitivityError>
    {
        self.admit_schedules(cx, max_controls)?;
        let Sdirk2HistoryParameterGradient { physical, stages } =
            self.pullback_with_parameters(cx, seed, linear, budget)?;
        Ok(Sdirk2HistoryParameterScheduleGradient {
            physical: self.with_schedules(cx, physical, max_controls)?, stages })
    }

    /// Both stage observations, shared authored schedules and independent original
    /// stage laws, retaining every physical control without another transpose.
    pub fn pullback_schedule_stages_with_parameters(&self, cx: &Cx<'_>, seeds: [&[f64]; 2],
        linear: LinearConfig, budget: EnthalpySdirk2ParameterBudget, max_controls: usize)
        -> Result<Sdirk2HistoryParameterScheduleGradient, ScanSensitivityError>
    {
        self.admit_schedules(cx, max_controls)?;
        let Sdirk2HistoryParameterGradient { physical, stages } =
            self.pullback_stages_with_parameters(cx, seeds, linear, budget)?;
        Ok(Sdirk2HistoryParameterScheduleGradient {
            physical: self.with_schedules(cx, physical, max_controls)?, stages })
    }

    fn with_parameters(&self, cx: &Cx<'_>, gradient: EnthalpySdirk2ParameterGradient)
        -> Result<Sdirk2HistoryParameterGradient, ScanSensitivityError>
    {
        Ok(Sdirk2HistoryParameterGradient {
            physical: self.with_sources(cx, gradient.transport)?, stages: gradient.stages })
    }

    fn with_sources(&self, cx: &Cx<'_>, transport: EnthalpySdirk2Gradient)
        -> Result<Sdirk2HistoryGradient, ScanSensitivityError>
    {
        let mut sources = [Vec::new(), Vec::new()];
        for stage in 0..2 {
            sources[stage].try_reserve_exact(self.sources[stage].len())
                .map_err(|_| invalid("SDIRK2 source gradient allocation refused"))?;
            for (i, source) in self.sources[stage].iter().enumerate() {
                poll(cx, i)?;
                sources[stage].push(source.pullback(cx, &transport.source_density[stage])?);
            }
        }
        poll(cx, self.declarations.len())?;
        Ok(Sdirk2HistoryGradient { transport, sources })
    }

    fn admit_schedules(&self, cx: &Cx<'_>, maximum: usize) -> Result<(), ScanSensitivityError> {
        let mut count = baths::control_count(cx, self.boundary)?;
        for (i, source) in self.declarations.iter().enumerate() {
            poll(cx, i)?;
            count = count.checked_add(source.schedule_control_count()?)
                .ok_or_else(|| invalid("SDIRK2 authored-control count overflows"))?;
        }
        if count > maximum { return Err(invalid("SDIRK2 schedules exceed their shared control allowance")); }
        Ok(())
    }

    fn with_schedules(&self, cx: &Cx<'_>, physical: Sdirk2HistoryGradient, max_controls: usize)
        -> Result<Sdirk2HistoryScheduleGradient, ScanSensitivityError>
    {
        let mut sources = Vec::new();
        sources.try_reserve_exact(self.declarations.len())
            .map_err(|_| invalid("SDIRK2 schedule gradient allocation refused"))?;
        let mut prescribed_temperatures = Vec::new();
        for stage in 0..2 {
            let time = self.stage_times_s[stage];
            let before = stage == 1 && self.end_before;
            for (i, (source, &gradient)) in self.declarations.iter()
                .zip(&physical.sources[stage]).enumerate()
            {
                poll(cx, i)?;
                let contribution = if before {
                    source.pullback_schedule_before(cx, time, gradient, max_controls)
                } else { source.pullback_schedule_at(cx, time, gradient, max_controls) }?;
                if stage == 0 { sources.push(contribution); }
                else { accumulate_source(cx, &mut sources[i], contribution)?; }
            }
            baths::accumulate(cx, &mut prescribed_temperatures, baths::pullback(cx,
                self.boundary, time, before, &physical.transport.prescribed_temperature[stage])?)?;
        }
        poll(cx, sources.len())?;
        Ok(Sdirk2HistoryScheduleGradient { physical, sources, prescribed_temperatures })
    }
}

impl<'h, 'm, 'c, 'i> Sdirk2EnthalpyHistoryRun<'h, 'm, 'c, 'i> {
    /// Solve and bind the actual native two-stage interval at a fixed explicit
    /// endpoint. The original primal runs once; no field/time/counter is published
    /// until both stage tangents and source bindings have passed admission.
    /// Original derivative corners may refuse when the primal alone is valid.
    pub fn advance_linearized_to(&mut self, cx: &Cx<'_>, time_s: f64,
        config: EnthalpyStepConfig, budget: Sdirk2SensitivityBudget)
        -> Result<Sdirk2HistoryLinearization<'h, 'm>, ScanSensitivityError>
    {
        self.admit_one(cx)?;
        if self.motion_endpoint(cx, self.run.time_s, time_s)? < time_s {
            return Err(invalid("linearized SDIRK2 endpoint exceeds the selected source motion limit"));
        }
        let binding = self.linearize_interval(cx, time_s, false, config, budget)?;
        self.publish(cx, time_s, &binding.primal().second_stage)?;
        Ok(binding)
    }

    /// Bind one original event/motion-limited interval, holding its selected times
    /// and sampling sides fixed. The first stage remains right-continuous and
    /// the second uses exact pre-switch data. No derivative of event selection,
    /// motion caps or adaptive time control is claimed.
    pub fn advance_linearized_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: EnthalpyStepConfig, budget: Sdirk2SensitivityBudget)
        -> Result<Sdirk2HistoryLinearization<'h, 'm>, ScanSensitivityError>
    {
        self.admit_one(cx)?;
        let end = self.next_endpoint(cx, target_s, max_step_s)?;
        let binding = self.linearize_interval(cx, end, true, config, budget)?;
        self.publish(cx, end, &binding.primal().second_stage)?;
        Ok(binding)
    }

    fn linearize_interval(&self, cx: &Cx<'_>, end: f64, before: bool,
        config: EnthalpyStepConfig, budget: Sdirk2SensitivityBudget)
        -> Result<Sdirk2HistoryLinearization<'h, 'm>, ScanSensitivityError>
    {
        poll(cx, 0)?;
        let start = self.run.time_s;
        let (_, first_time) = stage_clock(start, end)?;
        let mesh = self.run.problem.boundary.mesh();
        if self.moving_sources.len().checked_mul(mesh.vertex_count())
            .and_then(|n| n.checked_mul(2)).is_none_or(|n| n > budget.max_stage_source_vertices)
        { return Err(invalid("SDIRK2 source tangents exceed their combined stage/source/vertex allowance")); }
        let times = [first_time, end];
        let mut sources = [Vec::new(), Vec::new()];
        for stage in 0..2 {
            sources[stage].try_reserve_exact(self.moving_sources.len())
                .map_err(|_| invalid("SDIRK2 source binding allocation refused"))?;
            for (i, source) in self.moving_sources.iter().enumerate() {
                poll(cx, i)?;
                sources[stage].push(if stage == 1 && before {
                    source.linearize_before(cx, mesh, times[stage])?
                } else { source.linearize_at(cx, mesh, times[stage])? });
            }
        }
        let transport = self.with_stage_inputs(cx, start, end, before, |inputs, dt|
            self.run.storage.linearize_sdirk2(cx, inputs, &self.run.specific_h,
                dt, config, budget.max_radiation_entries))??;
        let fractions = [fraction::bind(cx, &self.run.storage,
            &transport.primal().first_stage.specific_enthalpy_j_kg)?,
            fraction::bind(cx, &self.run.storage,
                &transport.primal().second_stage.specific_enthalpy_j_kg)?];
        poll(cx, self.moving_sources.len())?;
        Ok(Sdirk2HistoryLinearization { transport, sources,
            declarations: Arc::clone(&self.moving_sources), boundary: self.run.problem.boundary,
            fractions, start_time_s: start, stage_times_s: times, end_before: before })
    }
}

fn accumulate_source(cx: &Cx<'_>, target: &mut MovingHeatScheduleGradient,
    contribution: MovingHeatScheduleGradient) -> Result<(), ScanSensitivityError>
{
    for axis in 0..3 {
        add(cx, &mut target.center_values_m[axis], &contribution.center_values_m[axis])?;
        target.half_width_m[axis] = finite(target.half_width_m[axis] + contribution.half_width_m[axis])?;
    }
    add(cx, &mut target.absorbed_power_values_w, &contribution.absorbed_power_values_w)
}

fn add(cx: &Cx<'_>, target: &mut [f64], contribution: &[f64]) -> Result<(), ScanSensitivityError> {
    if target.len() != contribution.len() { return Err(invalid("SDIRK2 authored source layout changed")); }
    for (i, (a, b)) in target.iter_mut().zip(contribution).enumerate() {
        if i % 256 == 0 { poll(cx, i)?; }
        *a = finite(*a + *b)?;
    }
    Ok(())
}
fn zeros(count: usize) -> Result<Vec<f64>, ScanSensitivityError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| invalid("SDIRK2 objective seed allocation refused"))?;
    values.resize(count, 0.0);
    Ok(values)
}
fn finite(value: f64) -> Result<f64, ScanSensitivityError> {
    if value.is_finite() { Ok(value) } else { Err(EnthalpyAdjointError::NonFiniteArithmetic.into()) }
}
fn invalid(what: &'static str) -> ScanSensitivityError { EnthalpyAdjointError::InvalidInput(what).into() }
