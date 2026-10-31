//! Continuing observation of resolved flow, complete fields and finite ducts.
//!
//! Every source control belongs to the original physical interval. The shared
//! NSF observation clock projects actual interval means when the output clock
//! differs; the existing DG exterior and Rayleigh receiver supply microphone
//! pressure. No physical state is retimed or reconstructed from audio samples.

use super::{NsfThermoviscousDg, NsfThermoviscousDgFrame};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::performance::ApertureObservation;
use crate::nsf_dg::{NsfDgDrive, NsfDgError, same_clock};
use crate::nsf_dg::performance::{clock::IntervalClock, source_intervals_for};
use crate::pcm_wav::baffled::BaffledPressure;
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_lbm::lattice3::Lattice3;
use std::sync::Arc;

fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }
fn physical(error: NsfDgError) -> RenderError {
    RenderError::Voice(match error {
        NsfDgError::Cancelled => AcousticRealizeError::Cancelled,
        other => AcousticRealizeError::Nonlinear(other.to_string()),
    })
}

/// A finite stream observing one attached duct while all participants continue.
/// A failed step leaves every physical owner unchanged. Observation failure
/// retains the accepted physical step; either failure poisons playback and
/// preserves its completed output prefix. Candidate allocation remains owned
/// by the original solvers; no allocation-free or real-time claim is implied.
pub struct NsfThermoviscousDgPerformance<'a, 'c, L: Lattice3<27>> {
    system: NsfThermoviscousDg<'a, 'c, L>,
    observation_port: usize,
    observation: ApertureObservation,
    controls: Vec<NsfDgDrive>,
    receiver: Option<BaffledPressure>,
    observer_gas: Option<GasPropagationReport>,
    terminal_flows: Vec<f64>,
    sample_rate_hz: u32,
    maximum_block: usize,
    gate: Arc<CancelGate>,
    last: Option<NsfThermoviscousDgFrame>,
    clock: Option<IntervalClock>,
    projected_observer: bool,
    output_samples: u64,
    output_cursor: u64,
    poisoned: bool,
}

impl<'a, 'c, L: Lattice3<27>> NsfThermoviscousDgPerformance<'a, 'c, L> {
    /// One physical interval per output sample. `observation_port` is the FE
    /// identity of an attached DG, not the external NSF connection or an ordinal.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_drives(system: NsfThermoviscousDg<'a, 'c, L>, observation_port: usize,
        observation: ApertureObservation, sample_rate_hz: u32, controls: Vec<NsfDgDrive>,
        maximum_block: usize, gate: Arc<CancelGate>) -> Result<Self, RenderError>
    {
        if !same_clock(system.load().duration_s(), 1.0 / f64::from(sample_rate_hz)) {
            return Err(sizing("resolved complete-field playback requires its original matching clock"));
        }
        let output_samples = controls.len() as u64;
        Self::build(system, observation_port, observation, sample_rate_hz, controls,
            output_samples, maximum_block, gate, None, false)
    }

    /// Project the original physical interval means onto an independent output
    /// clock using the existing conservative overlap owner. Controls must cover
    /// exactly `source_intervals_for` original intervals. This piecewise-constant
    /// projection is not a band-limited anti-alias filter; source and observer
    /// clocks need refinement, and shared output decimation remains available.
    #[allow(clippy::too_many_arguments)]
    pub fn new_resampled_with_drives(system: NsfThermoviscousDg<'a, 'c, L>, observation_port: usize,
        observation: ApertureObservation, output_rate_hz: u32, controls: Vec<NsfDgDrive>,
        output_samples: u64, maximum_block: usize, gate: Arc<CancelGate>) -> Result<Self, RenderError>
    {
        let dt = system.load().duration_s();
        if controls.len() != source_intervals_for(dt, output_rate_hz, output_samples)? {
            return Err(sizing("resolved field controls must cover the original physical output window"));
        }
        if same_clock(dt, 1.0 / f64::from(output_rate_hz)) {
            if dt.to_bits() == (1.0 / f64::from(output_rate_hz)).to_bits() {
                return Self::new_with_drives(system, observation_port, observation,
                    output_rate_hz, controls, maximum_block, gate);
            }
            return Self::build(system, observation_port, observation, output_rate_hz,
                controls, output_samples, maximum_block, gate, None, true);
        }
        let clock = IntervalClock::new(dt, output_rate_hz, output_samples).map_err(sizing)?;
        Self::build(system, observation_port, observation, output_rate_hz, controls,
            output_samples, maximum_block, gate, Some(clock), true)
    }

    #[allow(clippy::too_many_arguments)]
    fn build(system: NsfThermoviscousDg<'a, 'c, L>, observation_port: usize,
        observation: ApertureObservation, sample_rate_hz: u32, controls: Vec<NsfDgDrive>,
        output_samples: u64, maximum_block: usize, gate: Arc<CancelGate>,
        clock: Option<IntervalClock>, projected_observer: bool) -> Result<Self, RenderError>
    {
        if gate.is_requested() { return Err(RenderError::Voice(AcousticRealizeError::Cancelled)); }
        if system.accepted_steps() != 0 || sample_rate_hz == 0
            || !(1..=65_536).contains(&maximum_block) || controls.is_empty()
            || controls.len() as u64 > system.limits().maximum_steps
            || controls.len() as u64 > system.load().system().limits().maximum_steps
            || controls.len() != source_intervals_for(system.load().duration_s(), sample_rate_hz, output_samples)?
        {
            return Err(sizing("resolved field playback requires cold owners and a bounded complete control history"));
        }
        for &drive in &controls {
            if gate.is_requested() { return Err(RenderError::Voice(AcousticRealizeError::Cancelled)); }
            system.validate_drive(drive).map_err(physical)?;
        }
        let load = system.load().system().load_at_port(observation_port)
            .ok_or_else(|| sizing("resolved field observation must name an attached original duct"))?;
        let receiver = match observation {
            ApertureObservation::Inlet => None,
            ApertureObservation::TubeBaffled(location) => Some(if projected_observer {
                load.baffled_receiver_projected(sample_rate_hz, location)
            } else { load.baffled_receiver(sample_rate_hz, location) }.map_err(RenderError::Voice)?),
            _ => return Err(sizing("resolved field observation requires a duct inlet or its computed microphone")),
        };
        let terminal_flows = if receiver.is_some() { vec![0.0; load.terminal_ports().len()] }
            else { Vec::new() };
        if gate.is_requested() { return Err(RenderError::Voice(AcousticRealizeError::Cancelled)); }
        Ok(Self { system, observation_port, observation, controls, receiver,
            observer_gas: None, terminal_flows, sample_rate_hz, maximum_block, gate,
            last: None, clock, projected_observer, output_samples, output_cursor: 0, poisoned: false })
    }

    /// Retain the selected duct's original exterior gas propagation on every
    /// microphone path, including attenuation, phase and causal travel time.
    pub fn with_observer_gas(mut self, propagation: &GasPropagation) -> Result<Self, RenderError> {
        if self.poisoned || self.system.accepted_steps() != 0 || self.observer_gas.is_some() {
            return Err(sizing("resolved field atmosphere requires one cold microphone attachment"));
        }
        let ApertureObservation::TubeBaffled(location) = self.observation else {
            return Err(sizing("resolved field atmosphere requires an exterior microphone"));
        };
        let load = self.system.load().system().load_at_port(self.observation_port)
            .ok_or_else(|| sizing("resolved field microphone lost its original duct"))?;
        let (receiver, report) = if self.projected_observer {
            load.baffled_receiver_projected_with_gas(self.sample_rate_hz, location, propagation, &self.gate)
        } else { load.baffled_receiver_with_gas(self.sample_rate_hz, location, propagation, &self.gate) }
            .map_err(RenderError::Voice)?;
        if self.gate.is_requested() { return Err(RenderError::Voice(AcousticRealizeError::Cancelled)); }
        self.receiver = Some(receiver);
        self.observer_gas = Some(report);
        Ok(self)
    }

    pub const fn system(&self) -> &NsfThermoviscousDg<'a, 'c, L> { &self.system }
    pub fn drive_history(&self) -> &[NsfDgDrive] { &self.controls }
    pub const fn last_frame(&self) -> Option<&NsfThermoviscousDgFrame> { self.last.as_ref() }
    pub const fn observation_port(&self) -> usize { self.observation_port }
    pub const fn observation(&self) -> ApertureObservation { self.observation }
    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> { self.receiver.as_ref() }
    pub fn observer_gas(&self) -> Option<&GasPropagationReport> { self.observer_gas.as_ref() }
    pub fn cancellation_gate(&self) -> &CancelGate { &self.gate }
    pub fn remaining_samples(&self) -> u64 { self.output_samples - self.output_cursor }

    fn advance_source(&mut self) -> Result<(), RenderError> {
        let drive = *self.controls.get(self.system.accepted_steps() as usize)
            .ok_or_else(|| sizing("resolved complete-field source history is exhausted"))?;
        self.last = Some(self.system.step_with_drive(drive, &self.gate).map_err(physical)?);
        Ok(())
    }

    fn accumulate(&mut self, weight: f64) -> Result<f64, RenderError> {
        let system = self.system.load().system();
        let frame = system.last_frame().ok_or_else(|| sizing("resolved acoustic field has no accepted interval"))?;
        let pressure = if self.observation_port == system.port() { frame.mean_pressure_pa }
            else { frame.additional_ducts.iter().find(|duct| duct.port == self.observation_port)
                .ok_or_else(|| sizing("resolved acoustic field omitted its observed duct"))?.mean_pressure_pa };
        if self.receiver.is_some() {
            let terminals = system.load_at_port(self.observation_port)
                .ok_or_else(|| sizing("resolved field microphone lost its original duct"))?.last_terminal_frames();
            if terminals.is_empty() || terminals.len() != self.terminal_flows.len() {
                return Err(sizing("resolved field microphone requires every original terminal flow"));
            }
            for (flow, terminal) in self.terminal_flows.iter_mut().zip(terminals) {
                *flow += weight * terminal.outward_flow_m3_s;
            }
        }
        Ok(weight * pressure)
    }

    fn advance(&mut self) -> Result<f64, RenderError> {
        if self.gate.is_requested() { return Err(RenderError::Voice(AcousticRealizeError::Cancelled)); }
        self.terminal_flows.fill(0.0);
        let pressure = if let Some(clock) = self.clock {
            let (mut cursor, end) = clock.output_bounds(self.output_cursor).map_err(sizing)?;
            let mut pressure = 0.0;
            while cursor < end {
                while self.last.is_none() || clock.source_end(self.system.accepted_steps()) <= cursor {
                    self.advance_source()?;
                }
                let (weight, next) = clock.overlap(cursor, end, self.system.accepted_steps()).map_err(sizing)?;
                pressure += self.accumulate(weight)?;
                cursor = next;
            }
            pressure
        } else {
            self.advance_source()?;
            self.accumulate(1.0)?
        };
        let pressure = if let Some(receiver) = &mut self.receiver {
            receiver.step_trace(&self.terminal_flows)
                .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))?
        } else { pressure };
        if !pressure.is_finite() { return Err(sizing("resolved field observation produced nonfinite pressure")); }
        self.output_cursor += 1;
        Ok(pressure)
    }
}

impl<L: Lattice3<27>> PressureRenderer for NsfThermoviscousDgPerformance<'_, '_, L> {
    fn samples_rendered(&self) -> u64 { self.output_cursor }
    fn max_block_len(&self) -> usize { self.maximum_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.sample_rate_hz { return Err(sizing("resolved field pressure uses its admitted observer clock")); }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() { return Err(sizing("resolved field output exceeds its finite history")); }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.maximum_block { return Err(sizing("resolved field callback exceeds its declared capacity")); }
        for slot in output {
            match self.advance() {
                Ok(value) => *slot = value,
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}
