//! Finite-cell thermal feedback on the original rotor Maxwell memories.
//!
//! The existing thermal owner selects every interval after conduction, event
//! and user step limits. Its accepted start temperatures set the original
//! calibrated branch clocks; the native rotor interval then integrates motion,
//! source work and each branch's heat. That heat enters the same caloric cells
//! as interval-average watts. The mechanical endpoint is accepted only after
//! the thermal step and endpoint calibration admission both succeed.
//!
//! This is first-order lagged-temperature coupling, requiring independent
//! max_dt_s refinement. Exact held-coefficient mechanics and energy closure do
//! not make continuously varying temperature feedback exact. Cell uniformity,
//! original caloric laws, conductances and calibration support remain caller
//! obligations. C/K/spin, branch stiffness/strain and geometry stay fixed; no
//! Kelvin scaling, thermoelastic work or passivity is inferred. Only Maxwell
//! branch heat is transferred; signed C and circulatory-K work remain work.

use crate::rotor_kinetics::{RotorTemperatureError, RotorTemperatureKinetics};
use fs_material::visco::temperature::TemperatureShift;
use fs_math::dd::Dd;
use fs_mbd::rotordynamics::periodic::power::RotorCyclePowerError;
use fs_mbd::rotordynamics::periodic::{
    RotorFloquetBudget, RotorIntervalForcing, RotorMaxwellIntervalPower, RotorMaxwellState,
    RotorPeriodicSegment, prepare_maxwell_power_interval,
};
use fs_mbd::rotordynamics::relaxation::{MAX_ROTOR_MAXWELL_BRANCHES, MaxwellDampingBranch};
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig,
};
use std::collections::BTreeMap;

/// Original held mechanics, accepted history and explicit physical heat sites.
/// Continuations must retain these original storage laws and branch identities.
#[derive(Clone, Copy, Debug)]
pub struct RotorThermalInput<'a> {
    /// Must equal the original thermal timeline's accepted time.
    pub time_s: f64,
    pub state: &'a RotorMaxwellState,
    pub model: &'a RotorBearingModel,
    pub spin_rad_s: f64,
    pub branches: &'a [MaxwellDampingBranch],
    /// Each reference time belongs to this calibration's reference temperature.
    pub temperature_shifts: &'a [TemperatureShift],
    /// One cell per original branch: both its temperature source and receiver
    /// of ALL its heat. Repeated cells are additive. No temperature averaging,
    /// partition fractions or unmodeled heat sink is introduced.
    pub thermal_cells: &'a [CellId],
    /// One source row for this held model. Frequencies are physical rad/s and
    /// phase uses global time, independently of the selected thermal timestep.
    pub forcing: RotorIntervalForcing<'a>,
}

/// Independent native accuracy, thermal stepping and complete-call limits.
#[derive(Clone, Copy, Debug)]
pub struct RotorThermalConfig {
    /// Per-interval native dense work/payload admission.
    pub rotor: RotorFloquetBudget,
    pub relative_generator_residual_tolerance: f64,
    pub relative_map_residual_tolerance: f64,
    pub thermal: StepConfig,
    /// Whole-call interval ceiling, also tightened by the thermal owner's
    /// budget. Together with rotor.max_multiply_terms this bounds total work.
    pub maximum_intervals: usize,
    /// Gates summed absolute mechanical and heat-handoff defects, and the
    /// combined balance, independently [J]. Numerical defects never become heat.
    pub energy_tolerance_j: f64,
}

/// Native thermal accounting and the original mechanical ports, kept distinct.
#[derive(Clone, Debug)]
pub struct RotorThermalReport {
    pub thermal: EnergyAudit,
    pub applied_work_j: f64,
    pub velocity_force_work_j: f64,
    pub nonconservative_stiffness_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub relaxed_elastic_energy_change_j: f64,
    pub branch_elastic_energy_change_j: Vec<f64>,
    /// Independently integrated physical heat in original branch order.
    pub branch_heat_j: Vec<f64>,
    /// Integral of the generated interval source rates, already included once
    /// in thermal.external_heat_j. This is an INTERNAL mechanical/thermal transfer.
    pub deposited_branch_heat_j: f64,
    /// Raw native external heat less the represented branch-source integral.
    pub external_thermal_heat_j: f64,
    pub mechanical_balance_residual_j: f64,
    pub absolute_mechanical_residual_j: f64,
    pub heat_handoff_residual_j: f64,
    pub absolute_heat_handoff_residual_j: f64,
    /// Delta mechanical + delta thermal storage + C/K_skew work - applied work
    /// - incoming thermal inventory - external thermal heat. Computed through
    /// native residuals and handoff defects, never by subtracting large absolute
    /// caloric reference energies. This is conservation admission, not passivity.
    pub combined_energy_residual_j: f64,
}

/// No mechanical or thermal half-step is published on any refusal.
#[must_use]
pub struct RotorThermalCandidate<'thermal> {
    state: RotorMaxwellState,
    thermal: DepositionSimulation<'thermal>,
    report: RotorThermalReport,
}
impl<'thermal> RotorThermalCandidate<'thermal> {
    pub fn state(&self) -> &RotorMaxwellState {
        &self.state
    }
    pub fn thermal(&self) -> &DepositionSimulation<'thermal> {
        &self.thermal
    }
    pub fn report(&self) -> &RotorThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        RotorMaxwellState,
        DepositionSimulation<'thermal>,
        RotorThermalReport,
    ) {
        (self.state, self.thermal, self.report)
    }
}

#[derive(Debug)]
pub enum RotorThermalError {
    Invalid(&'static str),
    Rotor(RotorError),
    Power(RotorCyclePowerError),
    Temperature(RotorTemperatureError),
    Thermal(DepositionError),
    Energy { residual_j: f64, allowance_j: f64 },
}
impl core::fmt::Display for RotorThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rotor thermal coupling refused: {self:?}")
    }
}
impl std::error::Error for RotorThermalError {}
impl From<RotorError> for RotorThermalError {
    fn from(error: RotorError) -> Self {
        Self::Rotor(error)
    }
}
impl From<RotorCyclePowerError> for RotorThermalError {
    fn from(error: RotorCyclePowerError) -> Self {
        Self::Power(error)
    }
}
impl From<RotorTemperatureError> for RotorThermalError {
    fn from(error: RotorTemperatureError) -> Self {
        Self::Temperature(error)
    }
}
impl From<DepositionError> for RotorThermalError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// Advance full q/v/stress and the ORIGINAL finite thermal timeline together.
///
/// Same-time native births finish before receiver admission. Thermal events and
/// conduction select exact intervals; every original stress is carried forward
/// without reset. Endpoint temperatures and effective times/rates must still
/// be within calibration support before any candidate is returned. Authored
/// thermal sources and incoming birth energy retain their native accounting.
/// Neither input is mutated, including on late cancellation or work exhaustion.
pub fn advance_candidate<'thermal>(
    input: RotorThermalInput<'_>,
    thermal: &DepositionSimulation<'thermal>,
    target_s: f64,
    config: RotorThermalConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RotorThermalCandidate<'thermal>, RotorThermalError> {
    let mut check = || {
        if cancelled() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    if input.time_s.to_bits() != thermal.time_s().to_bits()
        || input.branches.len() > MAX_ROTOR_MAXWELL_BRANCHES
        || input.temperature_shifts.len() != input.branches.len()
        || input.thermal_cells.len() != input.branches.len()
        || config.maximum_intervals == 0
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j < 0.0
        || !relative_tolerance(config.relative_generator_residual_tolerance)
        || !relative_tolerance(config.relative_map_residual_tolerance)
        || input.state.displacement.len() != input.model.dofs()
        || input.state.velocity.len() != input.model.dofs()
        || input.state.branch_stress.len() != input.branches.len()
        || input
            .state
            .displacement
            .iter()
            .chain(&input.state.velocity)
            .chain(&input.state.branch_stress)
            .any(|x| !x.is_finite())
    {
        return Err(RotorThermalError::Invalid(
            "complete bounded original rotor, thermal bindings and controls required",
        ));
    }
    native(&mut check, |cancel| {
        input
            .model
            .at_speed_with_maxwell_memory(input.spin_rad_s, input.branches, cancel)
    })?;
    native(&mut check, |cancel| {
        input.forcing.validate(input.model.dofs(), 1, cancel)
    })?;
    let observer = Coupled {
        input,
        config,
        time_s: input.time_s,
        target_s,
        state: input.state.clone(),
        pending: None,
        terms: [Dd::ZERO; 5],
        branch_storage: vec![Dd::ZERO; input.branches.len()],
        branch_heat: vec![Dd::ZERO; input.branches.len()],
        deposited: Dd::ZERO,
        mechanical_residual: Dd::ZERO,
        absolute_mechanical: Dd::ZERO,
        absolute_handoff: Dd::ZERO,
    };
    let (thermal, coupled, audit) = thermal.advance_observed_limited(
        &observer,
        target_s,
        config.thermal,
        config.maximum_intervals,
        &mut check,
    )?;
    if coupled.pending.is_some() || coupled.time_s.to_bits() != thermal.time_s().to_bits() {
        return Err(RotorThermalError::Invalid(
            "unaccepted rotor thermal endpoint",
        ));
    }
    // Also admits zero-duration calls and temperatures after final-time events.
    coupled.kinetics(&thermal, &mut check)?;
    let heat = coupled
        .branch_heat
        .iter()
        .copied()
        .fold(Dd::ZERO, |sum, x| sum + x);
    let handoff = coupled.deposited - heat;
    let combined = coupled.mechanical_residual + Dd::from_f64(audit.residual_j) + handoff;
    for residual in [
        coupled.mechanical_residual,
        coupled.absolute_mechanical,
        handoff,
        coupled.absolute_handoff,
        combined,
    ] {
        admit(residual, config.energy_tolerance_j)?;
    }
    let report = RotorThermalReport {
        thermal: audit,
        applied_work_j: value(coupled.terms[0])?,
        velocity_force_work_j: value(coupled.terms[1])?,
        nonconservative_stiffness_work_j: value(coupled.terms[2])?,
        kinetic_energy_change_j: value(coupled.terms[3])?,
        relaxed_elastic_energy_change_j: value(coupled.terms[4])?,
        branch_elastic_energy_change_j: coupled
            .branch_storage
            .iter()
            .copied()
            .map(value)
            .collect::<Result<_, _>>()?,
        branch_heat_j: coupled
            .branch_heat
            .iter()
            .copied()
            .map(value)
            .collect::<Result<_, _>>()?,
        deposited_branch_heat_j: value(coupled.deposited)?,
        external_thermal_heat_j: value(Dd::from_f64(audit.external_heat_j) - coupled.deposited)?,
        mechanical_balance_residual_j: value(coupled.mechanical_residual)?,
        absolute_mechanical_residual_j: value(coupled.absolute_mechanical)?,
        heat_handoff_residual_j: value(handoff)?,
        absolute_heat_handoff_residual_j: value(coupled.absolute_handoff)?,
        combined_energy_residual_j: value(combined)?,
    };
    check()?;
    Ok(RotorThermalCandidate {
        state: coupled.state,
        thermal,
        report,
    })
}

#[derive(Clone)]
struct Pending {
    start_s: f64,
    dt_s: f64,
    end_s: f64,
    state: RotorMaxwellState,
    power: RotorMaxwellIntervalPower,
    deposited: Dd,
    absolute_handoff: Dd,
}

#[derive(Clone)]
struct Coupled<'a> {
    input: RotorThermalInput<'a>,
    config: RotorThermalConfig,
    time_s: f64,
    target_s: f64,
    state: RotorMaxwellState,
    pending: Option<Pending>,
    // Applied, signed C, signed K_skew, kinetic change, relaxed elastic change.
    terms: [Dd; 5],
    branch_storage: Vec<Dd>,
    branch_heat: Vec<Dd>,
    deposited: Dd,
    mechanical_residual: Dd,
    absolute_mechanical: Dd,
    absolute_handoff: Dd,
}
impl<'a> Coupled<'a> {
    fn kinetics(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<RotorTemperatureKinetics<'a>, RotorThermalError> {
        check()?;
        let mut temperatures = Vec::with_capacity(self.input.branches.len());
        for &cell in self.input.thermal_cells {
            check()?;
            temperatures.push(thermal.temperature_k(cell)?);
        }
        native(check, |cancel| {
            RotorTemperatureKinetics::try_new(
                self.input.branches,
                self.input.temperature_shifts,
                &[temperatures],
                cancel,
            )
        })
    }
}

impl DepositionObserver for Coupled<'_> {
    type Error = RotorThermalError;
    fn requires_step_interval(&self) -> bool {
        true
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        check()?;
        let stop = thermal
            .next_event_time_s()
            .map_or(self.target_s, |event| event.min(self.target_s));
        let remaining = stop - thermal.time_s();
        let half = 0.5 * remaining;
        if proposed_s < remaining && half < proposed_s {
            // Balance the final two intervals before native event/end clipping.
            // Repeated floating addition must not create an artificial terminal
            // sliver whose state increment cannot be resolved. Both intervals
            // remain below the actual CFL/user bound; no time or heat is dropped.
            let middle = thermal.time_s() + half;
            if !half.is_finite() || half <= 0.0 || middle <= thermal.time_s() || middle >= stop {
                return Err(DepositionError::Arithmetic.into());
            }
            Ok(half)
        } else {
            Ok(proposed_s)
        }
    }

    fn visit_step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        check()?;
        if self.pending.is_some()
            || self.time_s.to_bits() != start_s.to_bits()
            || thermal.time_s().to_bits() != start_s.to_bits()
        {
            return Err(RotorThermalError::Invalid(
                "rotor interval already pending or stale",
            ));
        }
        let kinetics = self.kinetics(thermal, check)?;
        let segments = [RotorPeriodicSegment {
            model: self.input.model,
            duration_seconds: dt_s,
        }];
        let plan = native(check, |cancel| {
            prepare_maxwell_power_interval(
                self.input.spin_rad_s,
                start_s,
                &segments,
                self.input.branches,
                Some(kinetics.as_kinetics()),
                self.input.forcing,
                self.config.relative_generator_residual_tolerance,
                self.config.rotor,
                cancel,
            )
        })?;
        let (trajectory, power) = native(check, |cancel| {
            plan.advance(
                &self.state,
                self.config.relative_map_residual_tolerance,
                self.config.energy_tolerance_j,
                cancel,
            )
        })?;
        if !trajectory
            .boundary_times_seconds
            .last()
            .is_some_and(|&time| time.to_bits() == end_s.to_bits())
        {
            return Err(RotorThermalError::Invalid(
                "native rotor endpoint differs from the selected thermal endpoint",
            ));
        }
        if power.branch_heat_j.len() != self.input.branches.len()
            || power.branch_elastic_energy_change_j.len() != self.input.branches.len()
        {
            return Err(RotorThermalError::Invalid(
                "native rotor account lost an original branch",
            ));
        }
        let mut deposited = Dd::ZERO;
        let mut absolute_handoff = Dd::ZERO;
        for (&cell, &heat) in self.input.thermal_cells.iter().zip(&power.branch_heat_j) {
            check()?;
            let watts = heat / dt_s;
            if !heat.is_finite()
                || heat < 0.0
                || !watts.is_finite()
                || (heat != 0.0 && watts == 0.0)
            {
                return Err(DepositionError::Arithmetic.into());
            }
            // Stream every original contribution. The thermal owner performs
            // compensated repeated-cell summation, including authored sources.
            emit(cell, watts)?;
            let delivered = Dd::from_f64(watts) * Dd::from_f64(dt_s);
            deposited = deposited + delivered;
            absolute_handoff =
                absolute_handoff + Dd::from_f64(value(delivered - Dd::from_f64(heat))?.abs());
        }
        admit(
            self.absolute_handoff + absolute_handoff,
            self.config.energy_tolerance_j,
        )?;
        let state =
            trajectory
                .boundary_states
                .last()
                .cloned()
                .ok_or(RotorThermalError::Invalid(
                    "native rotor interval lacks endpoint",
                ))?;
        self.pending = Some(Pending {
            start_s,
            dt_s,
            end_s,
            state,
            power,
            deposited,
            absolute_handoff,
        });
        Ok(())
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if self.pending.is_some() || self.time_s.to_bits() != thermal.time_s().to_bits() {
            return Err(RotorThermalError::Invalid(
                "event outside accepted rotor thermal clock",
            ));
        }
        // All same-time births finish before any branch needs its receiver.
        Ok(())
    }

    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        check()?;
        let pending = self.pending.take().ok_or(RotorThermalError::Invalid(
            "thermal step lacks rotor interval",
        ))?;
        if pending.start_s.to_bits() != start_s.to_bits()
            || pending.dt_s.to_bits() != dt_s.to_bits()
            || pending.end_s.to_bits() != thermal.time_s().to_bits()
        {
            return Err(RotorThermalError::Invalid(
                "thermal step changed the rotor interval",
            ));
        }
        self.kinetics(thermal, check)?;
        let power = pending.power;
        for (sum, x) in self.terms.iter_mut().zip([
            power.applied_work_j,
            power.velocity_force_work_j,
            power.nonconservative_stiffness_work_j,
            power.kinetic_energy_change_j,
            power.relaxed_elastic_energy_change_j,
        ]) {
            *sum = *sum + Dd::from_f64(x);
        }
        for ((storage, heat), (&delta, &loss)) in self
            .branch_storage
            .iter_mut()
            .zip(&mut self.branch_heat)
            .zip(
                power
                    .branch_elastic_energy_change_j
                    .iter()
                    .zip(&power.branch_heat_j),
            )
        {
            check()?;
            *storage = *storage + Dd::from_f64(delta);
            *heat = *heat + Dd::from_f64(loss);
        }
        self.deposited = self.deposited + pending.deposited;
        self.mechanical_residual =
            self.mechanical_residual + Dd::from_f64(power.balance_residual_j);
        self.absolute_mechanical =
            self.absolute_mechanical + Dd::from_f64(power.absolute_segment_residual_j);
        self.absolute_handoff = self.absolute_handoff + pending.absolute_handoff;
        admit(self.absolute_mechanical, self.config.energy_tolerance_j)?;
        admit(self.absolute_handoff, self.config.energy_tolerance_j)?;
        self.state = pending.state;
        self.time_s = pending.end_s;
        check()?;
        Ok(())
    }
}

fn relative_tolerance(x: f64) -> bool {
    x.is_finite() && x > 0.0 && x < 1.0
}
fn value(x: Dd) -> Result<f64, RotorThermalError> {
    let x = x.hi + x.lo;
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn admit(residual: Dd, allowance_j: f64) -> Result<(), RotorThermalError> {
    let residual_j = value(residual)?;
    if residual_j.abs() <= allowance_j {
        Ok(())
    } else {
        Err(RotorThermalError::Energy {
            residual_j,
            allowance_j,
        })
    }
}

// Native boolean cancellation retains the original thermal checkpoint error.
fn native<T, E: Into<RotorThermalError>>(
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, E>,
) -> Result<T, RotorThermalError> {
    let mut refusal = None;
    let result = operation(&mut || match check() {
        Ok(()) => false,
        Err(error) => {
            refusal = Some(error);
            true
        }
    });
    if let Some(error) = refusal {
        Err(error.into())
    } else {
        result.map_err(Into::into)
    }
}
