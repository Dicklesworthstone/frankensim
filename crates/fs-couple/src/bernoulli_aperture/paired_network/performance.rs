//! Scheduled paired-blade mechanics through the common physical-pressure sink.
//!
//! Owns the original PairedApertureNetwork, including both material histories,
//! contact, section gases and reactive loads. The existing canonical pressure
//! compiler sets the upstream supply. Independent blade-force assignments are
//! physical newtons on the source clock, never pitch, gain or state resets.
//! PressureRenderer supplies the existing streaming, decimation and ensemble
//! seam; it does not replace the mechanical, material or radiation owners.

use super::{PairedApertureNetwork, PairedNetworkFrame};
use crate::bernoulli_aperture::performance::{
    ApertureObservation, AperturePerformanceConfig, compile_schedule, pressure_candidate,
};
use crate::bernoulli_aperture::performance::force::ForceSchedule;

mod plate_forces;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::pcm_wav::baffled::BaffledPressure;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::{
    GestureCompileError, ScheduledControl,
};
use fs_scenario::gesture::GestureSchedule;

/// One held-force change at a physical sample boundary. Forces act independently
/// on each opening-positive blade coordinate, not on the total slit opening.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BladeForceControl {
    /// Zero-based mechanical sample; this force applies to that sample's solve.
    pub sample: u64,
    /// Original source blade, 0 or 1.
    pub blade: usize,
    /// Signed force [N]. Zero releases the actuator, not stored vibration.
    pub force_n: f64,
}

/// A finite performance retaining the complete paired physical system.
///
/// Source pressure and blade forces are applied within callbacks at their exact
/// admitted sample boundaries. Different block sizes cannot smear those events.
/// Observations are explicitly internal pressure or one existing baffled outlet
/// receiver of actual volume flow. There is no normalization or pressure gain.
///
/// This adapter's event loop allocates nothing per callback. Mechanical and
/// material owners keep their existing allocation behavior. A failed physics or
/// receiver step poisons the adapter: an accepted prefix must not be replayed.
/// Cancellation belongs to the existing stream writer's callback boundaries.
pub struct PairedAperturePerformance {
    system: PairedApertureNetwork,
    observation: ApertureObservation,
    receiver: Option<BaffledPressure>,
    config: AperturePerformanceConfig,
    schedule: GestureSchedule,
    controls: Vec<ScheduledControl>,
    next: usize,
    held_pressure_pa: f64,
    forces: Option<Vec<BladeForceControl>>,
    // Source-mesh force programs are an alternative to already-generalized
    // blade forces, never a second hidden supply or duplicated force law.
    plate_forces: [Option<ForceSchedule>; 2],
    next_force: usize,
    held_forces_n: [f64; 2],
    last: Option<PairedNetworkFrame>,
    poisoned: bool,
}
fn invalid(what: &'static str) -> GestureCompileError { GestureCompileError::Invalid { what } }
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }
fn receiver_error(error: String) -> RenderError {
    RenderError::Voice(crate::acoustic_realize::AcousticRealizeError::Nonlinear(error))
}

impl PairedAperturePerformance {
    /// Compile one canonical blowing-pressure track against an unadvanced,
    /// complete paired network. Nonzero initial blade/material energy is retained.
    /// The existing compiler owns interpolation, typed pressure admission and
    /// control-to-mechanics clock conversion. No replacement pressure envelope.
    ///
    /// # Errors
    /// Wrong clock/window, unsupported observation, missing/multiple tracks,
    /// late commands, compiler work/storage budget or existing physical refusal.
    pub fn new(system: PairedApertureNetwork, observation: ApertureObservation,
        schedule: GestureSchedule, config: AperturePerformanceConfig)
        -> Result<Self, GestureCompileError>
    {
        let a = system.aperture();
        if config.sample_rate_hz == 0 || config.samples == 0 || config.max_block == 0
            || a.accepted_steps() != 0 || config.samples > a.spec().max_steps
            || a.spec().time_step_s.to_bits() != (1.0 / f64::from(config.sample_rate_hz)).to_bits()
        { return Err(invalid("paired performance needs an unadvanced complete system and matching finite source clock")); }
        let controls = compile_schedule(&schedule, &config)?;
        let receiver = match observation {
            ApertureObservation::Inlet => None,
            ApertureObservation::NetworkNode(n) if system.node_frame(n).is_some() => None,
            ApertureObservation::NetworkBaffled { node, receiver } => Some(
                system.baffled_receiver(node, config.sample_rate_hz, receiver)
                    .map_err(|e| GestureCompileError::Render(RenderError::Voice(e)))?),
            _ => return Err(invalid("paired observation must name this network's pressure or an admitted outlet receiver")),
        };
        Ok(Self { system, observation, receiver, config, schedule, controls, next: 0,
            held_pressure_pa: 0.0, forces: None, plate_forces: [None, None],
            next_force: 0, held_forces_n: [0.0; 2],
            last: None, poisoned: false })
    }

    /// Bind sample-addressed, independently held opening forces before playback.
    /// Strict ordering is (sample, blade); duplicate assignments refuse rather
    /// than silently replacing a load. Forces start at zero. Pressure and force
    /// assignments together must fit the original control-storage allowance.
    /// This supplies actuator work, not an inferred lip/contact constitutive law.
    /// These are already-generalized forces in each retained opening coordinate.
    /// For physical node/patch loads use `with_plate_forces` instead; the two
    /// representations cannot be attached together and accidentally double a load.
    ///
    /// # Errors
    /// Late/repeated attachment, unordered/duplicate/out-of-window controls,
    /// invalid blade, nonfinite force or combined storage-budget exhaustion.
    pub fn with_blade_forces(mut self, controls: Vec<BladeForceControl>)
        -> Result<Self, GestureCompileError>
    {
        if self.poisoned || self.system.aperture().accepted_steps() != 0 || self.forces.is_some()
            || self.plate_forces.iter().any(Option::is_some) {
            return Err(invalid("paired force program requires one unadvanced attachment"));
        }
        if self.controls.len().checked_add(controls.len())
            .is_none_or(|n| n > self.config.max_controls) {
            return Err(invalid("paired pressure and force controls exceed the combined storage budget"));
        }
        let mut previous = None;
        for control in &controls {
            let key = (control.sample, control.blade);
            if control.sample >= self.config.samples || control.blade >= 2
                || !control.force_n.is_finite() || previous.is_some_and(|p| key <= p) {
                return Err(invalid("paired force assignments require finite newtons and unique ordered in-window blade samples"));
            }
            previous = Some(key);
        }
        self.forces = Some(controls);
        Ok(self)
    }

    /// Complete read-only mechanics, material sources/histories and acoustic graph.
    #[must_use]
    pub const fn system(&self) -> &PairedApertureNetwork { &self.system }
    /// Original source gesture program, not a resampled or normalized substitute.
    #[must_use]
    pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    /// Fixed source clock, horizon and work/storage limits.
    #[must_use]
    pub const fn config(&self) -> &AperturePerformanceConfig { &self.config }
    /// Physical location and type of the emitted pressure.
    #[must_use]
    pub const fn observation(&self) -> ApertureObservation { self.observation }
    /// Original outlet receiver and its retained propagation-delay metadata.
    #[must_use]
    pub fn baffled_receiver(&self) -> Option<&BaffledPressure> { self.receiver.as_ref() }
    /// Last completely observed mechanical sample, including both blade balances.
    /// On a failed callback consult the poisoned status before using its prefix.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&PairedNetworkFrame> { self.last.as_ref() }
    /// Pressure assignments committed with successfully observed samples.
    #[must_use]
    pub fn applied_controls(&self) -> &[ScheduledControl] { &self.controls[..self.next] }
    /// Unconsumed canonical pressure assignments on the mechanical clock.
    #[must_use]
    pub fn pending_controls(&self) -> &[ScheduledControl] { &self.controls[self.next..] }
    /// Blade assignments committed with successfully observed samples.
    #[must_use]
    pub fn applied_force_controls(&self) -> &[BladeForceControl] {
        &self.forces.as_deref().unwrap_or(&[])[..self.next_force]
    }
    /// Unconsumed physical blade assignments, not interpolated audio gains.
    #[must_use]
    pub fn pending_force_controls(&self) -> &[BladeForceControl] {
        &self.forces.as_deref().unwrap_or(&[])[self.next_force..]
    }
    /// Remaining mechanical samples, unaffected by observation resampling.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 {
        self.config.samples - self.system.aperture().accepted_steps()
    }

    fn step_pressure(&mut self, pressure: f64, forces: [f64; 2])
        -> Result<(f64, PairedNetworkFrame), RenderError>
    {
        let frame = self.system.step_with_forces(TubeDrive {
            upstream_pressure_pa: pressure, body_flow_m3_s: 0.0,
        }, forces).map_err(RenderError::Voice)?;
        let pressure = match self.observation {
            ApertureObservation::Inlet => frame.mechanics.aperture.bore_pressure_pa,
            ApertureObservation::NetworkNode(node) => self.system.node_frame(node)
                .expect("owned admitted node").pressure_pa,
            ApertureObservation::NetworkBaffled { node, .. } => self.receiver.as_mut()
                .expect("owned admitted receiver").step(&[
                    self.system.node_frame(node).expect("owned admitted outlet").net_flow_into_node_m3_s,
                ]).map_err(receiver_error)?,
            _ => unreachable!("construction refused non-network observation"),
        };
        if !pressure.is_finite() {
            return Err(receiver_error("paired observer returned nonfinite physical pressure".into()));
        }
        Ok((pressure, frame))
    }
}

impl PressureRenderer for PairedAperturePerformance {
    fn samples_rendered(&self) -> u64 { self.system.aperture().accepted_steps() }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.config.sample_rate_hz { return Err(sizing("paired output must use its actual mechanical sample clock")); }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() { return Err(sizing("paired pressure request exceeds the finite source horizon")); }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block { return Err(sizing("paired callback exceeds its admitted capacity")); }
        for slot in output {
            let sample = self.samples_rendered();
            let (next, pressure) = pressure_candidate(
                &self.controls, self.next, self.held_pressure_pa, sample);
            let (mut next_force, mut forces) = (self.next_force, self.held_forces_n);
            if let Some(program) = &self.forces {
                while let Some(event) = program.get(next_force) {
                    if event.sample != sample { break; }
                    forces[event.blade] = event.force_n;
                    next_force += 1;
                }
            } else {
                forces=std::array::from_fn(|i|self.plate_forces[i].as_ref()
                    .map_or(0.0,|f|f.candidate(sample)));
            }
            match self.step_pressure(pressure, forces) {
                Ok((value, frame)) => {
                    *slot = value; self.last = Some(frame);
                    self.next = next; self.held_pressure_pa = pressure;
                    self.next_force = next_force; self.held_forces_n = forces;
                    for program in self.plate_forces.iter_mut().flatten() {
                        program.accept(sample);
                    }
                }
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod plate_tests;
