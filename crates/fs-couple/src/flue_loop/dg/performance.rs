//! Scheduled internal mean-inlet pressure from the shared DG jet composition.
//!
//! This keeps the existing flue control compiler, original accepted source
//! history and common PressureRenderer/PCM interface. It adds no exterior
//! radiation or source physics. A caller may use the shared decimating observer
//! explicitly when the mechanical rate exceeds the output rate.

use super::super::network::FlueResonator;
use super::super::network::performance::controls::{
    self, FlueControlLimits, ScheduledFlueControl,
};
use super::super::{FlueError, FlueVoice};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::GestureCompileError;
use fs_scenario::gesture::GestureSchedule;

/// Explicit source grammar using the same material-gas and DG input owners.
pub mod file;

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

pub struct FlueDgPerformance {
    voice: FlueVoice,
    schedule: GestureSchedule,
    config: FlueDgConfig,
    controls: Vec<ScheduledFlueControl>,
    next: usize,
    launched: Vec<f64>,
    poisoned: bool,
}

impl FlueDgPerformance {
    /// Compile original pressure/angle gestures for a cold, source-bound DG
    /// voice. Only mean inlet pressure is observed; no terminal-flow surrogate
    /// or acoustic normalization is inferred.
    pub fn new(
        voice: FlueVoice,
        schedule: GestureSchedule,
        config: FlueDgConfig,
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
        Ok(Self { voice, schedule, config, controls, next: 0, launched, poisoned: false })
    }

    pub const fn voice(&self) -> &FlueVoice { &self.voice }
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    pub const fn config(&self) -> &FlueDgConfig { &self.config }
    pub fn compiled_controls(&self) -> usize { self.controls.len() }

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
            let mut destinations = output[written..written + count].iter_mut();
            self.voice.step_block_observed(&mut self.launched[..count], |line| {
                let FlueResonator::Dg(duct) = line else {
                    return Err(FlueError::Invalid { what: "DG pressure lost its original bore" });
                };
                let pressure = duct.last_frame().ok_or(FlueError::Invalid {
                    what: "DG pressure has no accepted original port frame",
                })?.inlet_pressure_pa;
                if !pressure.is_finite() || pressure.abs() > ceiling {
                    return Err(FlueError::Invalid {
                        what: "DG mean inlet pressure exceeds the declared finite ceiling",
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
