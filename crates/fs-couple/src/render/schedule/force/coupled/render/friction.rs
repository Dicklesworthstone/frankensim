//! Shared block rendering and typed physical controls for implicit friction.
//! The force/reaction solve remains in the existing modal/friction owners;
//! event timing remains in ScheduledRenderer. No second interpolation loop.

use super::super::contact::friction::{FrictionDrive, FrictionModalSystem};
use crate::render::{ControlDelta, RenderContext, RenderError, RenderVoice};
use crate::render::schedule::{GestureCompileError, ScheduledControl, ScheduledRenderer};
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureValue};

/// One real two-way network, with held external modal forces and surface inputs.
/// A callback refusal poisons this wrapper; the lower physical owner still
/// retains its last complete sample transaction. No allocation qualification
/// or pressure-observer accuracy beyond the supplied modal images is implied.
pub struct FrictionModalVoice {
    system: FrictionModalSystem,
    held_force: Vec<f64>,
    drive: FrictionDrive,
    poisoned: bool,
}
impl FrictionModalVoice {
    /// Retain every component, connection, and physical clock without a restart.
    pub fn new(system: FrictionModalSystem, held_force: Vec<f64>, drive: FrictionDrive)
        -> Result<Self, RenderError>
    {
        if held_force.len() != system.mode_count() || held_force.iter().any(|x| !x.is_finite()) {
            return Err(invalid("friction voice requires finite external forces for every mode"));
        }
        let voice = Self { system, held_force, drive, poisoned: false };
        voice.validate_drive(drive)?;
        Ok(voice)
    }
    /// Read-only accepted mechanics and complete work/loss diagnostics.
    #[must_use]
    pub const fn system(&self) -> &FrictionModalSystem { &self.system }
    /// Actual shared mechanical sample period [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { self.system.sample_period_s() }
    /// Validate a future input without claiming its as-yet-unsolved force is safe.
    pub fn validate_drive(&self, drive: FrictionDrive) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if !drive.speed_m_s.is_finite() || !drive.normal_force_n.is_finite() || drive.normal_force_n < 0.0 {
            return Err(invalid("friction control requires finite speed and nonnegative normal force"));
        }
        Ok(())
    }
    /// External force is separate from the internal friction reaction.
    pub fn validate_force(&self, mode: usize, force: f64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if mode >= self.held_force.len() || !force.is_finite() {
            return Err(invalid("friction external-mode control is outside its finite basis"));
        }
        Ok(())
    }
    pub(in crate::render) fn set_drive_admitted(&mut self, drive: FrictionDrive) { self.drive = drive; }
    pub(in crate::render) fn set_force_admitted(&mut self, mode: usize, force: f64) { self.held_force[mode] = force; }
    /// Render original modal pressure with the jointly solved physical reaction.
    /// On failure discard the callback, including any completed earlier samples.
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if out.is_empty() { return Err(RenderError::EmptyBlock); }
        let count = u64::try_from(out.len()).map_err(|_| invalid("friction callback exceeds u64"))?;
        self.system.samples_rendered().checked_add(count)
            .ok_or_else(|| invalid("friction callback would overflow the physical clock"))?;
        for sample in out {
            match self.system.step(&self.held_force, self.drive) {
                Ok(frame) => *sample = frame.observer_pressure_pa,
                Err(error) => { self.poisoned = true; return Err(RenderError::Friction(error)); }
            }
        }
        Ok(())
    }
}
fn invalid(what: &'static str) -> RenderError { RenderError::Control { what } }
fn source_invalid(what: &'static str) -> GestureCompileError { GestureCompileError::Invalid { what } }

/// Explicit physical binding and offline work/storage limits.
#[derive(Clone, Copy, Debug)]
pub struct FrictionGestureConfig {
    /// Must match the actual network's rate, never retunes or resamples it.
    pub sample_rate_hz: u32,
    /// Half-open compilation horizon [0,samples), relative to the binding.
    /// Rendering later holds the last input; author a zero-load event to release.
    pub samples: u64,
    /// Largest host callback [samples].
    pub max_block: usize,
    /// Upper bound on stateless source track/event visits during compilation.
    pub max_work: u64,
    /// Maximum emitted drive assignments, including the initial input.
    pub max_events: usize,
    /// Caller-declared bow station represented by the fixed attachment column.
    /// All source values must keep it unchanged; no moving-basis fiction.
    pub station_fraction: f64,
}

impl ScheduledRenderer {
    /// Bind one explicit BowStroke track to the network's physical friction port.
    ///
    /// The caller supplies the mechanically valid attachment shapes for the
    /// declared station. This adapter cannot prove they came from that geometry.
    /// Other source tracks are not executed. The selected track's speed/load
    /// ramps use GestureSchedule::sample_value, including interrupted ramps.
    /// Tick k applies at ceil(k*audio/control), and is held to the next tick.
    /// Existing vibration and the physical clock are preserved; only the new
    /// outer performance timeline begins at zero. No pitch or force is inferred
    /// from a note, and no normal collision is invented from the prescribed load.
    ///
    /// After the compile horizon the last input remains held, as with the shared
    /// pressure/force compilers. Use EnsembleRenderer for a finite output window.
    pub fn from_friction_gesture(
        system: FrictionModalSystem, held_force: Vec<f64>, schedule: &GestureSchedule,
        track_id: &str, config: FrictionGestureConfig,
    ) -> Result<Self, GestureCompileError> {
        let rate = schedule.control_rate_hz;
        if rate == 0 || config.sample_rate_hz == 0 || rate > config.sample_rate_hz || config.max_block == 0
            || system.sample_period_s().to_bits() != f64::from(config.sample_rate_hz).recip().to_bits()
            || !config.station_fraction.is_finite() || !(0.0 < config.station_fraction && config.station_fraction < 1.0) {
            return Err(source_invalid("friction gesture requires matching clocks, positive capacity and an interior station"));
        }
        system.samples_rendered().checked_add(config.samples)
            .ok_or_else(|| source_invalid("friction gesture horizon overflows the physical clock"))?;
        let track = schedule.tracks().iter().find(|t| t.id == track_id)
            .ok_or_else(|| GestureCompileError::Track { track: track_id.into(), what: "no source track with this id" })?;
        if !matches!(track.target, GestureTarget::BowStroke { .. }) {
            return Err(source_invalid("friction gesture requires an explicit bow-stroke track"));
        }
        // Validate every authored station, including future commands, before
        // lowering any values. A fixed port never silently ignores a moving bow.
        for value in std::iter::once(&track.initial).chain(track.events.iter().map(|e| &e.value)) {
            let GestureValue::Bow { station, .. } = value else {
                return Err(source_invalid("bow track contains a non-bow value"));
            };
            if station.to_bits() != config.station_fraction.to_bits() {
                return Err(source_invalid("moving bow stations require a new physical attachment model"));
            }
        }
        let audio = u128::from(config.sample_rate_hz);
        let control = u128::from(rate);
        let ticks = if config.samples == 0 { 0 } else { u128::from(config.samples-1)*control/audio + 1 };
        let visits = schedule.tracks().len() as u128 + track.events.len() as u128 + 1;
        let required = ticks.checked_mul(visits).ok_or_else(|| source_invalid("friction gesture work overflows"))?;
        if required > u128::from(config.max_work) {
            return Err(GestureCompileError::WorkBudget { required, allowed: config.max_work });
        }
        let voice = FrictionModalVoice::new(system, held_force, FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 })
            .map_err(GestureCompileError::Render)?;
        let mut events = Vec::new();
        let mut previous = None;
        for tick in 0..ticks as u64 {
            let GestureValue::Bow { velocity_m_per_s, normal_force_n, .. } =
                schedule.sample_value(track_id, tick).map_err(GestureCompileError::Gesture)?
            else { return Err(source_invalid("bow track sampled a non-bow value")); };
            let drive = FrictionDrive { speed_m_s: velocity_m_per_s, normal_force_n };
            voice.validate_drive(drive).map_err(GestureCompileError::Render)?;
            let bits = [velocity_m_per_s.to_bits(), normal_force_n.to_bits()];
            if previous == Some(bits) { continue; }
            if events.len() == config.max_events {
                return Err(source_invalid("friction gesture exceeds its emitted-control budget"));
            }
            events.try_reserve(1).map_err(|_| source_invalid("cannot reserve friction gesture controls"))?;
            events.push(ScheduledControl {
                sample: (u128::from(tick)*audio).div_ceil(control) as u64,
                delta: ControlDelta::SetFrictionDrive { voice: 0, speed_m_s: velocity_m_per_s, normal_force_n },
            });
            previous = Some(bits);
        }
        let context = RenderContext::new(vec![RenderVoice::FrictionModal(Box::new(voice))], config.max_block);
        Self::new(context, events, config.max_events).map_err(GestureCompileError::Render)
    }
}
