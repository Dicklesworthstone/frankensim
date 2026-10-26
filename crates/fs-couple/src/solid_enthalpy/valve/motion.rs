//! Finite spatial heat and the existing source-bound modal trajectory.
//!
//! One interval advances pressure-driven motion with the accepted coefficients,
//! then advances the original thermal owner and transfers that same physical
//! displacement/momentum to the new thermal equilibrium and tangent modes.
//! This first-order split preserves the native state-transfer error/work
//! account; it does not convert parameter work or projection loss into heat.

use std::convert::Infallible;

use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_modal::transfer::{ModalTransfer, ModalTransferControls};
use fs_solid::reduce::thermal::ThermalPressureLoadedValve;
use fs_solid::reduce::thermal::motion::{
    ThermalValveModalMotion, ThermalValveMotionError, ThermalValveMotionFrame,
};
use fs_thermal::enthalpy::{DepositionSimulation, EnergyAudit, StepConfig};

use super::{HeatedValveError, HeatedValveReduction};
use crate::solid_enthalpy::SolidEnthalpy;

/// A failure leaves both accepted owners and their sample clock unchanged.
#[derive(Debug)]
pub enum HeatedValveMotionError<E = Infallible> {
    /// Original thermal/material/equilibrium/reduction refusal.
    Valve(HeatedValveError),
    /// Original modal transition or physical state-transfer refusal.
    Motion(ThermalValveMotionError),
    /// Caller-supplied observation cannot be formed on the new source basis.
    Observation(E),
    /// The absolute thermal and relative modal clocks cannot share an interval.
    Clock,
    /// Cancellation before joint publication.
    Cancelled,
}
impl<E> From<HeatedValveError> for HeatedValveMotionError<E> {
    fn from(error: HeatedValveError) -> Self {
        Self::Valve(error)
    }
}
impl<E> From<ThermalValveMotionError> for HeatedValveMotionError<E> {
    fn from(error: ThermalValveMotionError) -> Self {
        Self::Motion(error)
    }
}
impl<E: core::fmt::Display> core::fmt::Display for HeatedValveMotionError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Valve(error) => write!(f, "heated valve motion: {error}"),
            Self::Motion(error) => write!(f, "heated valve motion: {error}"),
            Self::Observation(error) => write!(f, "heated valve observation: {error}"),
            Self::Clock => f.write_str("heated valve thermal/modal intervals differ"),
            Self::Cancelled => f.write_str("heated valve motion cancelled"),
        }
    }
}
impl<E: core::fmt::Debug + core::fmt::Display> std::error::Error for HeatedValveMotionError<E> {}

fn poll<E>(cx: &Cx<'_>) -> Result<(), HeatedValveMotionError<E>> {
    cx.checkpoint()
        .map_err(|_| HeatedValveMotionError::Cancelled)
}

/// One finite thermal inventory and its continuing physical modal state.
///
/// The source mapping, reference geometry, supports, pressure footprint and
/// sample rate remain fixed. Every source update rebuilds the compatible
/// equilibrium/tangent and uses the existing mass-inner-product transfer,
/// including the changed equilibrium offset. Mode order is never state identity.
pub struct HeatedValveMotion<'s, 'm, 't, 'r> {
    source: &'s SolidEnthalpy<'m>,
    thermal: DepositionSimulation<'t>,
    reduction: HeatedValveReduction<'r>,
    motion: ThermalValveModalMotion,
    initial_time_s: f64,
}

/// Separate original physical accounts at one accepted absolute time.
#[derive(Clone, Debug)]
pub struct HeatedValveMotionStep {
    /// Absolute time shared with the accepted thermal inventory [s].
    pub time_s: f64,
    /// Native held-pressure transition before the endpoint source update.
    /// Its local time is elapsed from `initial_time_s`, not a second clock.
    pub pressure_step: ThermalValveMotionFrame,
    /// Original caloric balance. Parameter work is not debited as heat.
    pub thermal: EnergyAudit,
    /// Signed incremental-storage change and discarded-state/error accounts.
    /// This is not the complete prestressed body's thermodynamic energy.
    pub transfer: ModalTransfer,
    /// Observation evaluated on the accepted NEW source basis and state [Pa].
    /// The supplied narrow-band transfer does not imply radiation back-reaction.
    pub endpoint_observer_pressure_pa: f64,
}

impl<'s, 'm, 't, 'r> HeatedValveMotion<'s, 'm, 't, 'r> {
    /// Bind both owners to the same finite thermal state. Initial coordinates
    /// belong to the reduction made here. The observation callback sees that
    /// exact reduction and must return one explicit transfer per retained mode.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new<E>(
        source: &'s SolidEnthalpy<'m>,
        thermal: DepositionSimulation<'t>,
        reduction: HeatedValveReduction<'r>,
        sample_rate_hz: u32,
        initial: &[ModalAcousticState],
        observe: impl FnOnce(&ThermalPressureLoadedValve) -> Result<Vec<C64>, E>,
        budget: ModalAcousticTimeBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, HeatedValveMotionError<E>> {
        poll(cx)?;
        let state = source
            .resolve(&thermal, cx)
            .map_err(HeatedValveError::from)?;
        let valve = state.reduce_pressure_valve(&reduction, cx)?;
        let transfers = observe(&valve).map_err(HeatedValveMotionError::Observation)?;
        let motion =
            ThermalValveModalMotion::new(valve, sample_rate_hz, initial, &transfers, budget, cx)?;
        let initial_time_s = thermal.time_s();
        if !initial_time_s.is_finite() {
            return Err(HeatedValveMotionError::Clock);
        }
        poll(cx)?;
        Ok(Self {
            source,
            thermal,
            reduction,
            motion,
            initial_time_s,
        })
    }

    /// Accepted finite heat inventory, with its original source laws.
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }

    /// Actual accepted native motion and source reduction. Its local clock is
    /// elapsed time from `initial_time_s`; `time_s` supplies absolute time.
    pub const fn motion(&self) -> &ThermalValveModalMotion {
        &self.motion
    }

    /// Absolute origin of the native modal sample clock [s].
    pub const fn initial_time_s(&self) -> f64 {
        self.initial_time_s
    }

    /// Absolute accepted physical time [s].
    pub fn time_s(&self) -> f64 {
        self.thermal.time_s()
    }

    /// Advance one pressure/heat interval and publish both owners together.
    ///
    /// The original thermal observer admits every substep and event. The new
    /// observation is formed only after thermal/material/reduction acceptance,
    /// and is evaluated after physical state transfer. Any late refusal leaves
    /// heat, vibration, basis, observation and sample count unchanged.
    ///
    /// This is prescribed one-way thermal loading with explicit parameter work.
    /// Viscous loss and basis projection loss are not sent to the caloric owner;
    /// reciprocal thermoelastic heat, pressure feedback and moving footprints
    /// remain outside this model. Refine physical time and modal budget.
    pub fn step<E>(
        &mut self,
        pressure_pa: f64,
        thermal_config: StepConfig,
        controls: ModalTransferControls,
        observe: impl FnOnce(&ThermalPressureLoadedValve) -> Result<Vec<C64>, E>,
        cx: &Cx<'_>,
    ) -> Result<HeatedValveMotionStep, HeatedValveMotionError<E>> {
        poll(cx)?;
        let mut motion = self.motion.clone();
        let pressure_step = motion.step_pressure(pressure_pa, cx)?;
        let target_s = self.initial_time_s + pressure_step.time_s;
        let duration_s = target_s - self.thermal.time_s();
        let modal_dt_s = motion.model().sample_period_s();
        let absolute_clock_roundoff =
            16.0 * f64::EPSILON * (target_s.abs() + self.thermal.time_s().abs() + modal_dt_s);
        // A coarse absolute epoch must not authorize a physically different
        // interval. Retain at least half the significand's relative precision
        // in the actual sample duration, independently of the clock origin.
        let clock_tolerance = absolute_clock_roundoff.min(f64::EPSILON.sqrt() * modal_dt_s);
        if !target_s.is_finite()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || (duration_s - modal_dt_s).abs() > clock_tolerance
        {
            return Err(HeatedValveMotionError::Clock);
        }
        let candidate = self.source.advance_pressure_valve_to(
            &self.thermal,
            target_s,
            thermal_config,
            &self.reduction,
            cx,
        )?;
        let transfers = observe(&candidate.valve).map_err(HeatedValveMotionError::Observation)?;
        let (motion, transfer) =
            motion.rebind_candidate(candidate.valve, &transfers, controls, cx)?;
        let endpoint_observer_pressure_pa = motion
            .model()
            .observer_pressure_with_transfers(&transfers)
            .map_err(ThermalValveMotionError::from)?;
        poll(cx)?;
        self.thermal = candidate.thermal;
        self.motion = motion;
        Ok(HeatedValveMotionStep {
            time_s: target_s,
            pressure_step,
            thermal: candidate.thermal_audit,
            transfer,
            endpoint_observer_pressure_pa,
        })
    }
}
