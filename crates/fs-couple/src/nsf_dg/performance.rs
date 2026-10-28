//! Finite physical-pressure playback of the original resolved zone and DG bore.
//! The caller supplies incoming characteristic samples on the OTHER resolved
//! face, in m/s. This is physical boundary data, not a note, envelope or pressure
//! calibration. The coupled face is still solved against the live DG history.
//! Initial resolved motion and thermochemical material state are retained.

use super::{NsfDg, NsfDgFrame, same_clock};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::performance::ApertureObservation;
use crate::pcm_wav::baffled::BaffledPressure;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_lbm::lattice3::Lattice3;
use std::sync::Arc;

fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }

/// A finite source on the existing pressure, decimation, ensemble and PCM path.
/// Failed individual physical steps publish neither owner. A failed callback
/// poisons playback, retaining any already accepted prefix and observer state.
pub struct NsfDgPerformance<L: Lattice3<27>> {
    system: NsfDg<L>,
    incoming: Vec<f64>,
    observation: ApertureObservation,
    receiver: Option<BaffledPressure>,
    terminal_flows: Vec<f64>,
    sample_rate_hz: u32,
    maximum_block: usize,
    gate: Arc<CancelGate>,
    last: Option<NsfDgFrame>,
    poisoned: bool,
}

impl<L: Lattice3<27>> NsfDgPerformance<L> {
    /// Bind one held incoming-characteristic value per coupled interval.
    /// The sample rate matches the existing physical clock; explicit output
    /// decimation belongs to the shared PressureRenderer adapter.
    pub fn new(
        system: NsfDg<L>, observation: ApertureObservation, sample_rate_hz: u32,
        incoming_characteristics_m_s: Vec<f64>, maximum_block: usize,
        gate: Arc<CancelGate>,
    ) -> Result<Self, RenderError> {
        if gate.is_requested() || system.accepted_steps() != 0
            || sample_rate_hz == 0 || !(1..=65_536).contains(&maximum_block)
            || incoming_characteristics_m_s.is_empty()
            || incoming_characteristics_m_s.len() as u64 > system.limits().maximum_steps
            || !incoming_characteristics_m_s.iter().all(|x| x.is_finite())
            || !same_clock(system.load().propagator().dt(), 1.0 / f64::from(sample_rate_hz))
        {
            return Err(sizing("resolved performance requires a cold matching clock and finite boundary samples/capacity"));
        }
        let receiver = match observation {
            ApertureObservation::Inlet => None,
            ApertureObservation::TubeBaffled(location) => Some(
                system.load().baffled_receiver(sample_rate_hz, location)
                    .map_err(RenderError::Voice)?),
            _ => return Err(sizing("resolved pressure observation must name its DG inlet or computed baffled outlet")),
        };
        let terminal_flows = if receiver.is_some() {
            vec![0.0; system.load().terminal_ports().len()]
        } else { Vec::new() };
        Ok(Self {
            system, incoming: incoming_characteristics_m_s, observation, receiver,
            terminal_flows, sample_rate_hz, maximum_block, gate, last: None, poisoned: false,
        })
    }
    /// Original resolved source, material state and acoustic load.
    pub const fn system(&self) -> &NsfDg<L> { &self.system }
    /// Supplied physical incoming characteristic on the noncoupled face [m/s].
    pub fn incoming_characteristics_m_s(&self) -> &[f64] { &self.incoming }
    /// Last accepted physical interval, even if its later observation failed.
    pub const fn last_frame(&self) -> Option<&NsfDgFrame> { self.last.as_ref() }
    /// Physical location of the emitted pressure.
    pub const fn observation(&self) -> ApertureObservation { self.observation }
    /// Original reaction-bound retarded receiver and its explicit delay.
    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> { self.receiver.as_ref() }
    /// Shared cancellation owner, also polled during each nonlinear interval.
    pub fn cancellation_gate(&self) -> &CancelGate { &self.gate }
    /// Remaining physical source samples.
    pub fn remaining_samples(&self) -> u64 {
        self.incoming.len() as u64 - self.system.accepted_steps()
    }

    fn advance(&mut self) -> Result<f64, RenderError> {
        let incoming = self.incoming[self.system.accepted_steps() as usize];
        let frame = self.system.step(incoming, &self.gate)
            .map_err(|e| RenderError::Voice(AcousticRealizeError::Nonlinear(e.to_string())))?;
        // Accepted source samples are consumed before their fallible observer.
        // The cursor derives from the sole accepted physical clock.
        self.last = Some(frame);
        let frame = self.last.as_ref().expect("accepted original resolved interval");
        let pressure = match self.observation {
            ApertureObservation::Inlet => frame.inlet_pressure_pa,
            ApertureObservation::TubeBaffled(_) => {
                if frame.terminals.len() != self.terminal_flows.len() || frame.terminals.is_empty() {
                    return Err(sizing("resolved DG microphone requires every accepted terminal flow"));
                }
                for (flow, terminal) in self.terminal_flows.iter_mut().zip(&frame.terminals) {
                    *flow = terminal.outward_flow_m3_s;
                }
                self.receiver.as_mut().expect("admitted original DG receiver")
                    .step_trace(&self.terminal_flows)
                    .map_err(|e| RenderError::Voice(AcousticRealizeError::Nonlinear(e)))?
            }
            _ => unreachable!("observation admitted at construction"),
        };
        if !pressure.is_finite() {
            return Err(sizing("resolved observer produced nonfinite physical pressure"));
        }
        Ok(pressure)
    }
}

impl<L: Lattice3<27>> PressureRenderer for NsfDgPerformance<L> {
    fn samples_rendered(&self) -> u64 { self.system.accepted_steps() }
    fn max_block_len(&self) -> usize { self.maximum_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.sample_rate_hz {
            return Err(sizing("resolved pressure uses its admitted physical clock"));
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() {
            return Err(sizing("resolved pressure exceeds its finite source window"));
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.maximum_block {
            return Err(sizing("resolved pressure callback exceeds its declared capacity"));
        }
        for slot in output {
            match self.advance() {
                Ok(value) => *slot = value,
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}
