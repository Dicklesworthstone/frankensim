//! Authored pressure phrases on retained coupled valve/tube mechanics.
//!
//! This adapter owns the existing system, not a replacement voice. Plate shape,
//! spatial closure, material memory, wall/load state and traveling waves remain
//! in their physical owners. Only pressure assignments are scheduled. The output
//! is either a named internal pressure or an explicitly selected one-way baffled
//! outlet receiver. The latter observes actual terminal FLOW, not bore pressure.

/// Supplied mesh, closure, material memory and canonical pressure sources.
pub mod file;

use super::dynamic::DynamicAperture;
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use super::network::ApertureNetwork;
use super::tube::{ApertureTube, TubeDrive};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::{ControlDelta, RenderError};
use crate::render::schedule::{
    GestureCompileError, PressureGestureBinding, ScheduledControl, compile_pressure_gestures,
};
use fs_scenario::gesture::GestureSchedule;

/// Existing reciprocal physical systems; both retain their accepted state.
pub enum CoupledAperture {
    /// One geometry-bound uniform tube and its termination.
    Tube(ApertureTube),
    /// Explicit connected sections, junctions and any already-attached loads.
    Network(ApertureNetwork),
}
impl CoupledAperture {
    /// Immutable access to the actual valve, including plate/contact/material data.
    #[must_use]
    pub fn aperture(&self) -> &DynamicAperture {
        match self { Self::Tube(m) => m.aperture(), Self::Network(m) => m.aperture() }
    }
}

/// Which physical pressure is sampled once per complete mechanical step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ApertureObservation {
    /// Junction pressure used in the implicit pressure/structure solve.
    Inlet,
    /// Actual terminal pressure after propagation through a uniform tube.
    TubeTerminal,
    /// Actual accepted node pressure in a connected tube network.
    NetworkNode(usize),
    /// One-way exterior receiver of a uniform tube's outward terminal flow.
    /// The circular mouth lies in an infinite rigid baffle. The authored
    /// terminal reflection still drives the mechanics; no additional radiation
    /// reaction is silently applied or claimed to match that terminal load.
    TubeBaffled(CircularOutletReceiver),
    /// Exterior observation of flow into a degree-one passive network terminal.
    /// The adjacent section supplies the actual outlet radius. Its existing load
    /// determines reaction; observation does not add a second load or assume that
    /// an arbitrary supplied impedance represents this baffled mouth.
    NetworkBaffled {
        /// Degree-one passive terminal, not a junction, shunt, series node or inlet.
        node: usize,
        /// Position/quadrature relative to that outlet's outward-facing baffle.
        receiver: CircularOutletReceiver,
    },
}

/// Finite window and bounded offline pressure compilation.
#[derive(Clone, Copy, Debug)]
pub struct AperturePerformanceConfig {
    /// Must exactly match dt = 1 / sample_rate_hz in the supplied physical system.
    pub sample_rate_hz: u32,
    /// Positive half-open source window [0, samples); must fit its step budget.
    pub samples: u64,
    /// Largest source callback; independent of control timing.
    pub max_block: usize,
    /// Worst-case visits in the existing canonical gesture sampler.
    pub max_compile_work: u64,
    /// Bounds sampled control ticks AND stored assignments before compilation.
    pub max_controls: usize,
}

/// Pressure playback preserving one complete physical state and source clock.
///
/// Starts with an unadvanced system, which may have nonzero physical initial
/// motion/contact/material energy. The supplied track owns the pressure history;
/// no fixture attack envelope, reset at release, note-to-pressure calibration or
/// prescribed body-flow source is added. Continue by retaining this whole object.
/// Scheduling itself allocates nothing inside callbacks; physical-owner
/// allocations are unchanged. A failed physical callback poisons this adapter,
/// as required by PressureRenderer, even though the owner rolls back its failed
/// individual step. Completed earlier samples remain inspectable.
pub struct AperturePerformance {
    system: CoupledAperture,
    observation: ApertureObservation,
    receiver: Option<BaffledPressure>,
    config: AperturePerformanceConfig,
    schedule: GestureSchedule,
    controls: Vec<ScheduledControl>,
    next: usize,
    held_pressure_pa: f64,
    poisoned: bool,
}
fn receiver_error(error: String) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(error))
}
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }
fn invalid(what: &'static str) -> GestureCompileError { GestureCompileError::Invalid { what } }

impl AperturePerformance {
    /// Compile and admit the complete pressure phrase without advancing physics.
    /// No source command may begin after the last observed control tick. Ramps
    /// extending beyond the window retain only their actually observed portion.
    ///
    /// # Errors
    /// Clock/horizon/observer mismatch, unsupported or multiple tracks, exhausted
    /// compile/storage budgets, or a source gesture refusal. No runtime is returned
    /// with an incomplete schedule or a substituted pressure observation.
    pub fn new(system: CoupledAperture, observation: ApertureObservation,
        schedule: GestureSchedule, config: AperturePerformanceConfig)
        -> Result<Self, GestureCompileError>
    {
        let a = system.aperture();
        if config.sample_rate_hz == 0 || config.samples == 0 || config.max_block == 0
            || a.accepted_steps() != 0 || config.samples > a.spec().max_steps
            || a.spec().time_step_s.to_bits() != (1.0 / f64::from(config.sample_rate_hz)).to_bits()
        { return Err(invalid("aperture performance requires a matching clock, positive capacity/window and an unadvanced complete physical system")); }
        let receiver = match (&system, observation) {
            (_, ApertureObservation::Inlet) | (CoupledAperture::Tube(_), ApertureObservation::TubeTerminal) => None,
            (CoupledAperture::Network(m), ApertureObservation::NetworkNode(n)) if n < m.spec().nodes.len() => None,
            (CoupledAperture::Tube(m), ApertureObservation::TubeBaffled(location)) => {
                // Radius and gas are the actual physical tube's, never an
                // independently tuned observer area or invented gain. Waveguide
                // construction is quiescent, so the initial flow history is zero.
                Some(BaffledPressure::circular_outlet(m.spec().radius_m,
                    config.sample_rate_hz, location, RayleighMedium {
                        density: a.spec().density_kg_m3, sound_speed: m.spec().sound_speed_m_s,
                    }).map_err(|e| GestureCompileError::Render(receiver_error(e)))?)
            }
            (CoupledAperture::Network(m), ApertureObservation::NetworkBaffled {node,receiver}) => {
                use super::network::NetworkNode;
                if !matches!(m.spec().nodes.get(node), Some(NetworkNode::Termination {..}
                    | NetworkNode::Impedance {..} | NetworkNode::Relaxation {..})) {
                    return Err(invalid("baffled network observation requires a degree-one passive terminal"));
                }
                let section=m.spec().sections.iter().find(|s|s.nodes.contains(&node))
                    .ok_or_else(||invalid("baffled network terminal has no physical section"))?;
                Some(BaffledPressure::circular_outlet(section.radius_m,config.sample_rate_hz,receiver,
                    RayleighMedium {density:a.spec().density_kg_m3,sound_speed:m.spec().sound_speed_m_s})
                    .map_err(|e|GestureCompileError::Render(receiver_error(e)))?)
            }
            _ => return Err(invalid("pressure observation is not a node or terminal of the supplied system")),
        };
        let [track] = schedule.tracks() else {
            return Err(invalid("one explicitly typed blowing-pressure track is required"));
        };
        let rate = schedule.control_rate_hz;
        if rate == 0 || rate > config.sample_rate_hz {
            return Err(invalid("control rate must be positive and no faster than mechanics"));
        }
        let last_tick = u128::from(config.samples - 1) * u128::from(rate)
            / u128::from(config.sample_rate_hz);
        if last_tick + 1 > config.max_controls as u128 {
            return Err(invalid("pressure sampling exceeds the declared control-tick/storage budget"));
        }
        if track.events.iter().any(|e| e.time_s > last_tick as f64 / f64::from(rate)) {
            return Err(invalid("pressure command begins after the last observed control tick"));
        }
        let controls = compile_pressure_gestures(&schedule, &[PressureGestureBinding {
            track: track.id.clone(), voice: 0,
        }], config.sample_rate_hz, config.samples, config.max_compile_work)?;
        Ok(Self { system, observation, receiver, config, schedule, controls, next: 0,
            held_pressure_pa: 0.0, poisoned: false })
    }

    /// Read-only physical system; callers cannot step it behind the sample clock.
    #[must_use]
    pub const fn system(&self) -> &CoupledAperture { &self.system }
    /// Original source gestures, unchanged by audio resampling or block partition.
    #[must_use]
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    /// Exact finite clock/compilation choices.
    #[must_use]
    pub const fn config(&self) -> &AperturePerformanceConfig { &self.config }
    /// Named physical location of the output trace.
    #[must_use]
    pub const fn observation(&self) -> ApertureObservation { self.observation }
    /// Immutable receiver geometry/delay information, absent for internal traces.
    /// The physical propagation delay remains inside this source; downstream
    /// rate-conversion or ensemble alignment must not remove it.
    #[must_use]
    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> { self.receiver.as_ref() }
    /// Controls committed with successful mechanical samples, never ahead of them.
    #[must_use]
    pub fn applied_controls(&self) -> &[ScheduledControl] { &self.controls[..self.next] }
    /// Remaining controls on the original mechanical clock.
    #[must_use]
    pub fn pending_controls(&self) -> &[ScheduledControl] { &self.controls[self.next..] }
    /// Finite remaining source samples, not resampled output frames.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.config.samples - self.samples_rendered() }

    fn step_pressure(&mut self, pressure_pa: f64) -> Result<f64, RenderError> {
        let drive = TubeDrive { upstream_pressure_pa: pressure_pa, body_flow_m3_s: 0.0 };
        match &mut self.system {
            CoupledAperture::Tube(m) => {
                let f = m.step(drive).map_err(RenderError::Voice)?;
                Ok(match self.observation {
                    ApertureObservation::Inlet => f.aperture.bore_pressure_pa,
                    ApertureObservation::TubeTerminal => f.waveguide.terminal_pressure_pa,
                    ApertureObservation::TubeBaffled(_) => self.receiver.as_mut()
                        .expect("baffled observation constructed its complete receiver")
                        .step(&[f.waveguide.terminal_flow_m3_s]).map_err(receiver_error)?,
                    _ => unreachable!("observer admitted against the exclusively owned system"),
                })
            }
            CoupledAperture::Network(m) => {
                let f = m.step(drive).map_err(RenderError::Voice)?;
                Ok(match self.observation {
                    ApertureObservation::Inlet => f.aperture.bore_pressure_pa,
                    ApertureObservation::NetworkNode(n) => m.node_frame(n)
                        .expect("admitted node remains in the owned topology").pressure_pa,
                    ApertureObservation::NetworkBaffled {node,..} => {
                        let flow=m.node_frame(node).expect("admitted terminal remains in owned topology")
                            .net_flow_into_node_m3_s;
                        self.receiver.as_mut().expect("admitted network outlet receiver")
                            .step(&[flow]).map_err(receiver_error)?
                    }
                    _ => unreachable!("observer admitted against the exclusively owned system"),
                })
            }
        }
    }
}
impl PressureRenderer for AperturePerformance {
    fn samples_rendered(&self) -> u64 { self.system.aperture().accepted_steps() }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz { return Err(sizing("pressure output must use the physical aperture clock")); }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() { return Err(sizing("pressure request exceeds the complete aperture performance window")); }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block { return Err(sizing("aperture callback exceeds its construction-time capacity")); }
        for slot in output {
            let mut next = self.next;
            let mut pressure = self.held_pressure_pa;
            while let Some(event) = self.controls.get(next) {
                if event.sample != self.samples_rendered() { break; }
                if let ControlDelta::SetBlowingPressure { pressure_pa, .. } = event.delta {
                    pressure = pressure_pa;
                } else { unreachable!("existing compiler emits only pressure assignments"); }
                next += 1;
            }
            match self.step_pressure(pressure) {
                Ok(value) => { *slot = value; self.next = next; self.held_pressure_pa = pressure; }
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}
