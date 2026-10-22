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
use crate::bernoulli_aperture::performance::force::ForceSchedule;
use crate::bernoulli_aperture::performance::boundaries::BoundarySchedule;
use crate::bernoulli_aperture::two_sided::ValveSide;

mod plate_forces;
/// Mechanical-sample amplitudes of source-bound spatial pressure traces.
pub mod spatial_pressures;
/// Shared physical boundary programs and typed source/key gestures.
pub mod boundaries;

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
    plate_forces: [Option<ForceSchedule>; 2],
    spatial_pressures: [Option<spatial_pressures::PressureSchedule>; 2],
    boundaries: [Option<BoundarySchedule>; 2],
    gesture_work: u64,
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
        // Public window/callbacks count OUTPUT samples. The original compiler
        // sees the actual mechanical clock, including controls within an audio
        // frame. Work and stored-control budgets are not renewed per substep.
        let ratio = source.mechanics_per_output();
        let mechanical_config = AperturePerformanceConfig {
            sample_rate_hz: source.mechanics_rate_hz(),
            samples: config.samples.checked_mul(ratio as u64)
                .ok_or_else(|| invalid("two-sided mechanical phrase window overflow"))?,
            max_block: config.max_block.checked_mul(ratio)
                .ok_or_else(|| invalid("two-sided mechanical callback size overflow"))?,
            ..config
        };
        let controls = compile_schedule(&schedule, &mechanical_config)?;
        Ok(Self { source, schedule, config, controls, next: 0, held_pressure_pa: 0.0,
            plate_forces: [None, None], spatial_pressures: [None, None],
            boundaries: [None, None], gesture_work: 0, poisoned: false })
    }
    /// Original loaded mechanics and pressure receiver. Mutable access would
    /// permit stepping around the source schedule and is intentionally absent.
    #[must_use]
    pub const fn source(&self) -> &TwoSidedPressure { &self.source }
    /// Original canonical source gestures, not a second envelope implementation.
    #[must_use]
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    /// Output-clock window and the original whole-phrase compilation budgets.
    #[must_use]
    pub const fn config(&self) -> &AperturePerformanceConfig { &self.config }
    /// Source assignments committed with completely observed physical samples.
    #[must_use]
    pub fn applied_controls(&self) -> &[ScheduledControl] { &self.controls[..self.next] }
    /// Source assignments still on the original mechanical clock.
    #[must_use]
    pub fn pending_controls(&self) -> &[ScheduledControl] { &self.controls[self.next..] }
    /// Remaining complete output samples, including any declared receiver tail.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.config.samples - self.source.samples_rendered() }
    /// Set independent opening-positive mechanical forces [N] between callbacks.
    /// The pressure phrase is unchanged. Zero releases the actuator, not motion.
    /// The original source refuses nonexistent/nonfinite force ports; there is
    /// no new contact or lip law and no duplicated source of mechanical work.
    /// A bound source-mesh force program owns its inputs for the whole phrase;
    /// manual force changes then refuse instead of silently overriding the program.
    pub fn set_blade_forces(&mut self, forces_n: [f64; 2]) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if self.plate_forces.iter().any(Option::is_some) {
            return Err(RenderError::Control { what: "source-mesh force programs cannot be overridden by manual blade forces" });
        }
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
        if samples != 0 {
            for side in [ValveSide::Upstream, ValveSide::Downstream] {
                if self.source.system.boundary_apertures(side).is_some()
                    && self.boundaries[side.index()].is_none() {
                    return Err(RenderError::Control {
                        what: "each keyed acoustic side requires its explicit boundary program",
                    });
                }
            }
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        let count = u64::try_from(output.len()).map_err(|_| sizing("phrase sample count overflow"))?;
        self.validate_sample_count(count)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block { return Err(sizing("two-sided phrase callback exceeds its admitted capacity")); }
        for slot in output {
            for _ in 0..self.source.mechanics_per_output() {
                let sample = self.source.mechanical_samples_rendered();
                let (next, pressure) = pressure_candidate(&self.controls, self.next,
                    self.held_pressure_pa, sample);
                let forces = std::array::from_fn(|i| self.plate_forces[i].as_ref()
                    .map_or(self.source.held.blade_forces_n[i], |f| f.candidate(sample)));
                let drive = TwoSidedDrive {
                    source_pressure_pa: pressure, blade_forces_n: forces,
                };
                let spatial = std::array::from_fn(|i| self.spatial_pressures[i].as_ref()
                    .map_or(self.source.held_spatial_pressures_pa[i], |p| p.candidate(sample)));
                for boundaries in self.boundaries.iter_mut().flatten() { boundaries.prepare(sample); }
                let boundary_drives = self.boundaries.each_ref().map(|side|
                    side.as_ref().map_or(&[][..], BoundarySchedule::candidate));
                // Same physical/receiver step as held playback. Partial traces
                // survive every control boundary instead of endpoint sampling.
                let result = self.source.set_spatial_pressures(spatial)
                    .and_then(|()| self.source.set_drive(drive))
                    .and_then(|()| self.source.advance_mechanical_with_boundaries(boundary_drives));
                match result {
                    Ok(output) => {
                        self.next = next;
                        self.held_pressure_pa = pressure;
                        for force in self.plate_forces.iter_mut().flatten() { force.accept(sample); }
                        for pressure in self.spatial_pressures.iter_mut().flatten() { pressure.accept(sample); }
                        for boundaries in self.boundaries.iter_mut().flatten() { boundaries.accept(); }
                        if let Some(p) = output { *slot = p; }
                    }
                    Err(error) => { self.poisoned = true; return Err(error); }
                }
            }
        }
        Ok(())
    }
}
