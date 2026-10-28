//! Scheduled physical pressure from the shared DG jet composition.
//!
//! The default observation is the original interval-mean inlet pressure.
//! An explicit baffled microphone instead consumes every accepted exterior
//! terminal flows through the original retarded Rayleigh pressure owner, using
//! the same bound apertures and material gas as the computed reaction source.
//! This keeps the existing flue control compiler, accepted source histories,
//! finite pressure ceiling and common PressureRenderer/PCM interface.
//! An optional physical rectangular mouth shares that exterior and microphone.
//! Its compact self reaction can replace the legacy mouth load; it does not
//! include mutual mouth/outlet loading or resolve the surrounding jet flow.

use super::super::network::FlueResonator;
use super::super::network::performance::controls::{
    self, FlueControlLimits, ScheduledFlueControl,
};
use super::super::{FlueError, FlueVoice};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver};
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::GestureCompileError;
use fs_scenario::gesture::GestureSchedule;

/// Explicit source grammar using the same material-gas and DG input owners.
pub mod file;
mod mouth;

pub use crate::flue_loop::network::performance::FlueMouthRadiation;

fn invalid(what: &'static str) -> GestureCompileError {
    GestureCompileError::Render(RenderError::Voice(
        AcousticRealizeError::InvalidDescription { what },
    ))
}
fn source_error(error: FlueError) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string()))
}
fn sizing(what: &'static str) -> RenderError {
    RenderError::Voice(AcousticRealizeError::InvalidDescription { what })
}

/// The single admitted mechanical clock and its finite pressure window.
#[derive(Clone, Copy, Debug)]
pub struct FlueDgConfig {
    pub sample_rate_hz: u32,
    pub samples: u64,
    pub max_block: usize,
    pub max_controls: usize,
    pub max_compile_work: u64,
    pub maximum_abs_pressure_pa: f64,
}

struct BaffledOutletObservation {
    config: CircularOutletReceiver,
    receiver: BaffledPressure,
    terminal_flows_m3_s: Vec<f64>,
}

pub struct FlueDgPerformance {
    voice: FlueVoice,
    schedule: GestureSchedule,
    config: FlueDgConfig,
    controls: Vec<ScheduledFlueControl>,
    next: usize,
    launched: Vec<f64>,
    outlet: Option<BaffledOutletObservation>,
    mouth: Option<mouth::MouthPressure>,
    observer_gas: Option<GasPropagationReport>,
    poisoned: bool,
}

impl FlueDgPerformance {
    /// Compile original pressure/angle gestures for a cold, source-bound DG
    /// voice, observing its accepted mean inlet pressure.
    pub fn new(
        voice: FlueVoice,
        schedule: GestureSchedule,
        config: FlueDgConfig,
    ) -> Result<Self, GestureCompileError> {
        Self::build(voice, schedule, config, None)
    }

    /// Observe the held outward flows of every computed baffled DG terminal.
    ///
    /// The original radiation record supplies all apertures and exterior gas;
    /// neither is inferred from the inlet or independently re-baked. One accepted
    /// source sample enters the existing differentiating Rayleigh observer once,
    /// in the shared terminal owner's exact channel order.
    /// Both clocks are the declared source rate. Receiver band, quadrature and
    /// causal-history limits belong to that owner, including its 8..192 kHz
    /// same-clock range. Numerical filter delay is exposed by baffled_receiver().
    ///
    /// These are coplanar circular openings in an infinite rigid baffle, not free-field
    /// instrument-body/mouth radiation or a new qualification of the jet card.
    pub fn with_baffled_outlet(
        voice: FlueVoice,
        schedule: GestureSchedule,
        config: FlueDgConfig,
        receiver: CircularOutletReceiver,
    ) -> Result<Self, GestureCompileError> {
        Self::build(voice, schedule, config, Some(receiver))
    }

    fn build(
        voice: FlueVoice,
        schedule: GestureSchedule,
        config: FlueDgConfig,
        receiver: Option<CircularOutletReceiver>,
    ) -> Result<Self, GestureCompileError> {
        if config.sample_rate_hz == 0
            || config.samples == 0
            || config.samples > (1_u64 << 53)
            || !(1..=65_536).contains(&config.max_block)
            || voice.samples_rendered() != 0
            || voice.dg().is_none()
            || voice.sample_period_s().to_bits()
                != (1.0 / f64::from(config.sample_rate_hz)).to_bits()
            || !config.maximum_abs_pressure_pa.is_finite()
            || config.maximum_abs_pressure_pa <= 0.0
        {
            return Err(invalid("DG flue pressure requires a fresh source, matching clock and finite limits"));
        }
        voice.validate_blowing_pressure(0.0)
            .map_err(|error| GestureCompileError::Render(source_error(error)))?;
        let outlet = if let Some(receiver) = receiver {
            let load = voice.dg().expect("admitted original DG voice").load();
            let pressure = load.baffled_receiver(config.sample_rate_hz, receiver)
                .map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
            let mut terminal_flows_m3_s = Vec::new();
            terminal_flows_m3_s.try_reserve_exact(load.terminal_ports().len())
                .map_err(|_| invalid("DG flue terminal flow trace exceeds available allocation"))?;
            terminal_flows_m3_s.resize(load.terminal_ports().len(), 0.0);
            Some(BaffledOutletObservation { config: receiver, receiver: pressure,
                terminal_flows_m3_s })
        } else { None };
        let controls = controls::compile(
            &voice,
            &schedule,
            FlueControlLimits {
                max_controls: config.max_controls,
                max_compile_work: config.max_compile_work,
            },
            config.sample_rate_hz,
            config.samples,
        )?;
        let mut launched = Vec::new();
        launched.try_reserve_exact(config.max_block)
            .map_err(|_| invalid("DG flue callback exceeds available allocation"))?;
        launched.resize(config.max_block, 0.0);
        Ok(Self { voice, schedule, config, controls, next: 0, launched, outlet,
            mouth: None, observer_gas: None, poisoned: false })
    }

    pub const fn voice(&self) -> &FlueVoice { &self.voice }
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    pub const fn config(&self) -> &FlueDgConfig { &self.config }
    pub fn compiled_controls(&self) -> usize { self.controls.len() }

    /// The original exterior observer, including explicit numerical and travel delays.
    /// Absent for the compatible internal mean-inlet observation.
    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> {
        self.outlet.as_ref().map(|outlet| &outlet.receiver)
    }

    /// Original receiver position, surface resolution and declared physical band.
    pub fn baffled_receiver_config(&self) -> Option<CircularOutletReceiver> {
        self.outlet.as_ref().map(|outlet| outlet.config)
    }

    fn advance_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        let end = self.voice.samples_rendered() + output.len() as u64;
        let mut written = 0;
        while self.voice.samples_rendered() < end {
            let now = self.voice.samples_rendered();
            while let Some(control) = self.controls.get(self.next) {
                if control.sample != now { break; }
                self.voice.apply(control.delta).map_err(source_error)?;
                self.next += 1;
            }
            let until = self.controls.get(self.next).map_or(end, |control| control.sample.min(end));
            let count = (until - now) as usize;
            let ceiling = self.config.maximum_abs_pressure_pa;
            let outlet = &mut self.outlet;
            let mouth = &mut self.mouth;
            let mut destinations = output[written..written + count].iter_mut();
            self.voice.step_block_observed(&mut self.launched[..count], |line| {
                let FlueResonator::Dg(duct) = line else {
                    return Err(FlueError::Invalid { what: "DG pressure lost its original bore" });
                };
                let mut pressure = if let Some(observer) = outlet.as_mut() {
                    let frames = duct.load().last_terminal_frames();
                    if frames.is_empty() || frames.len() != observer.terminal_flows_m3_s.len() {
                        return Err(FlueError::Invalid {
                            what: "baffled DG pressure requires every accepted terminal flow in source order",
                        });
                    }
                    for (flow, frame) in observer.terminal_flows_m3_s.iter_mut().zip(frames) {
                        *flow = frame.outward_flow_m3_s;
                    }
                    observer.receiver.step_trace(&observer.terminal_flows_m3_s)
                        .map_err(|error| FlueError::Physics(
                            AcousticRealizeError::Nonlinear(error),
                        ))?
                } else {
                    duct.last_frame().ok_or(FlueError::Invalid {
                        what: "DG pressure has no accepted original port frame",
                    })?.inlet_pressure_pa
                };
                if let Some(mouth) = mouth.as_mut() {
                    let inward_flow = duct.last_frame().ok_or(FlueError::Invalid {
                        what: "DG mouth pressure has no accepted original inlet flow",
                    })?.inlet_flow_m3_s;
                    pressure += mouth.push(inward_flow)?;
                }
                if !pressure.is_finite() || pressure.abs() > ceiling {
                    return Err(FlueError::Invalid {
                        what: "DG observed pressure exceeds the declared finite ceiling",
                    });
                }
                *destinations.next().expect("one pressure per accepted physical sample") = pressure;
                Ok(())
            }).map_err(source_error)?;
            written += count;
        }
        Ok(())
    }
}

impl PressureRenderer for FlueDgPerformance {
    fn samples_rendered(&self) -> u64 { self.voice.samples_rendered() }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz {
            return Err(sizing("DG flue pressure uses its admitted physical source clock"));
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.config.samples - self.voice.samples_rendered() {
            return Err(sizing("DG flue pressure exceeds its declared window"));
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block {
            return Err(sizing("DG flue callback exceeds its declared capacity"));
        }
        let result = self.advance_block(output);
        if result.is_err() { self.poisoned = true; }
        result
    }
}

#[cfg(test)]
mod tests;

