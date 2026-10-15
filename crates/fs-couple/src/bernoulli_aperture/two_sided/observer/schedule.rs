//! Canonical pressure phrases on the original two-sided valve and receiver.
//! Reuses the scalar/paired performance compiler and its sample-boundary rules.
//! Pressure is an ideal SERIES supply, not a prescribed pressure on the reed
//! face: both acoustic networks continue to determine that actual face pressure.
use super::TwoSidedPressure;
use crate::bernoulli_aperture::performance::{
    AperturePerformanceConfig, compile_schedule, pressure_candidate,
};
use crate::bernoulli_aperture::two_sided::TwoSidedDrive;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::{RenderError, schedule::{GestureCompileError, ScheduledControl}};
use fs_scenario::gesture::GestureSchedule;

/// Finite canonical source program retaining all valve, duct and receiver state.
/// Every control is consumed at its actual mechanical sample, even inside a host
/// callback or a decimation group. A zero-pressure event releases the source;
/// it does not reset blades, material/contact histories or either acoustic side.
/// Physical owners keep their original allocations and approximation limits.
pub struct ScheduledTwoSidedPressure {
    source: TwoSidedPressure,
    schedule: GestureSchedule,
    config: AperturePerformanceConfig,
    controls: Vec<ScheduledControl>,
    next: usize,
    held_pressure_pa: f64,
    poisoned: bool,
}
fn invalid(what: &'static str) -> GestureCompileError { GestureCompileError::Invalid { what } }
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }

impl ScheduledTwoSidedPressure {
    /// Compile the existing typed blowing-pressure track without stepping physics.
    /// Starts at sample zero; an already held nonzero source pressure refuses
    /// rather than silently being replaced by a different initial condition.
    /// Previously declared independent blade forces are retained. No receiver,
    /// gas, geometry, modal state or material history is rebuilt on attachment.
    ///
    /// # Errors
    /// Invalid rate/window/capacity, late attachment, compiler work/storage
    /// exhaustion, unsupported source gestures or a prior physical refusal.
    pub fn new(source: TwoSidedPressure, schedule: GestureSchedule,
        config: AperturePerformanceConfig) -> Result<Self, GestureCompileError>
    {
        if config.samples == 0 || config.max_block == 0
            || config.max_block > source.max_block_len() || source.samples_rendered() != 0
            || source.held.source_pressure_pa != 0.0 {
            return Err(invalid("two-sided phrase needs a finite positive window, admitted capacity and unadvanced zero-pressure source"));
        }
        source.validate_sample_rate(config.sample_rate_hz).map_err(GestureCompileError::Render)?;
        source.validate_sample_count(config.samples).map_err(GestureCompileError::Render)?;
        let controls = compile_schedule(&schedule, &config)?;
        Ok(Self { source, schedule, config, controls, next: 0, held_pressure_pa: 0.0, poisoned: false })
    }
    /// Original loaded mechanics and pressure receiver. Mutable access would
    /// permit stepping around the source schedule and is intentionally absent.
    #[must_use]
    pub const fn source(&self) -> &TwoSidedPressure { &self.source }
    /// Original canonical source gestures, not a second envelope implementation.
    #[must_use]
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    /// Existing source-clock and compilation choices.
    #[must_use]
    pub const fn config(&self) -> &AperturePerformanceConfig { &self.config }
    /// Source assignments committed with completely observed physical samples.
    #[must_use]
    pub fn applied_controls(&self) -> &[ScheduledControl] { &self.controls[..self.next] }
    /// Source assignments still on the original mechanical clock.
    #[must_use]
    pub fn pending_controls(&self) -> &[ScheduledControl] { &self.controls[self.next..] }
    /// Remaining finite source window, independent of downstream decimation.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.config.samples - self.source.samples_rendered() }
    /// Set independent opening-positive mechanical forces [N] between callbacks.
    /// The pressure phrase is unchanged. Zero releases the actuator, not motion.
    /// The original source refuses nonexistent/nonfinite force ports; there is
    /// no new contact or lip law and no duplicated source of mechanical work.
    pub fn set_blade_forces(&mut self, forces_n: [f64; 2]) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        self.source.set_drive(TwoSidedDrive {
            source_pressure_pa: self.held_pressure_pa, blade_forces_n: forces_n,
        })
    }
}
impl PressureRenderer for ScheduledTwoSidedPressure {
    fn samples_rendered(&self) -> u64 { self.source.samples_rendered() }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        self.source.validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        self.source.validate_sample_count(samples)?;
        if samples > self.remaining_samples() {
            return Err(sizing("two-sided pressure request exceeds the canonical phrase window"));
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        let count = u64::try_from(output.len()).map_err(|_| sizing("phrase sample count overflow"))?;
        self.validate_sample_count(count)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block { return Err(sizing("two-sided phrase callback exceeds its admitted capacity")); }
        for slot in output {
            let (next, pressure) = pressure_candidate(
                &self.controls, self.next, self.held_pressure_pa, self.samples_rendered());
            let drive = TwoSidedDrive {
                source_pressure_pa: pressure, blade_forces_n: self.source.held.blade_forces_n,
            };
            // There is one physical/receiver loop, already owned by the source.
            // This adapter only supplies the next canonical control before it.
            let result = self.source.set_drive(drive)
                .and_then(|()| self.source.block(core::slice::from_mut(slot)));
            match result {
                Ok(()) => { self.next = next; self.held_pressure_pa = pressure; }
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}
