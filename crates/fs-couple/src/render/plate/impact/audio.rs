//! Physical far-field pressure from the existing nonlinear impact mechanics.
//!
//! Offline velocity-to-far-field samples come from an acoustic solver, not
//! from an instrument name or an arbitrary modal gain. Reuse the existing SH
//! change of basis, vector fit, withheld complex-response gate and Tustin bank.
//! Divide by -i*omega to drive that bank with actual velocity increments. This
//! avoids a fitted DC velocity gain radiating steady translation forever.
//!
//! A finite source can have an ADVANCE in its origin-referenced far field.
//! Add a/c to the sampled filter and propagate the remaining (r-a)/c in a
//! delay line: F exp(i*k*a) exp(i*k*(r-a))/r = F exp(i*k*r)/r. Omitting either
//! factor corrupts phase/arrival time. a must enclose the acoustic source about
//! its declared origin. It is caller-provided geometry, not fitted latency.
//!
//! The causal backward velocity difference adds a half-mechanical-step delay
//! and sinc attenuation relative to an exact derivative. Propagation rounds UP
//! to whole mechanical samples within an explicit tolerance. Both errors are
//! inspectable and converge with sample rate. Neither is time compensated.
//! This is stationary-reference, ONE-WAY far-field observation: no acoustic
//! reaction, near field, room, air absorption or moving-surface claim. Original
//! mechanics remain unchanged and are still an allocating reference solver.

use super::{ImpactFrame, ImpactSystem};
use crate::broadband_radiation::{
    BroadbandRadiationArtifact, BroadbandRadiationControls, BroadbandRadiationError,
    BroadbandRadiationRuntime, SampledRadiationData, build_broadband_radiation_artifact,
    evaluate_real_tesseral,
};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

/// A velocity-response sample bank converted to causal acceleration response.
/// Kept immutable so a caller cannot remove the source-time shift after fitting.
pub struct ImpactRadiation {
    bank: BroadbandRadiationArtifact,
    source_radius_m: f64,
    sound_speed_m_s: f64,
}
impl ImpactRadiation {
    /// Fit supplied F(omega)/velocity data under exp(-i*omega*t), including its
    /// ORIGINAL spatial phase about the source origin. Input IDs must describe
    /// the same velocity coordinates later supplied by VelocityProjection.
    /// The enclosing radius comes from that actual geometry. Both training and
    /// withheld data undergo the identical exact phase/unit conversion; the
    /// existing fitter still checks withheld complex data, not just magnitudes.
    pub fn from_velocity_samples(
        samples: &SampledRadiationData,
        controls: BroadbandRadiationControls,
        source_radius_m: f64,
        sound_speed_m_s: f64,
    ) -> Result<Self, BroadbandRadiationError> {
        if !source_radius_m.is_finite() || source_radius_m < 0.0
            || !sound_speed_m_s.is_finite() || sound_speed_m_s <= 0.0
            || !(source_radius_m / sound_speed_m_s).is_finite()
        {
            return Err(BroadbandRadiationError::InvalidInput("invalid source radius or acoustic speed"));
        }
        let mut converted = samples.clone();
        let factor = |omega: f64| -> Result<C64, BroadbandRadiationError> {
            let phase = omega * (source_radius_m / sound_speed_m_s);
            if !omega.is_finite() || omega <= 0.0 || !phase.is_finite() || !omega.recip().is_finite() {
                return Err(BroadbandRadiationError::InvalidInput("radiation frequencies must have a finite positive reciprocal"));
            }
            // 1/(-i*omega) = +i/omega under the acoustic owner's convention.
            Ok(C64::new(-det::sin(phase) / omega, det::cos(phase) / omega))
        };
        for row in &mut converted.training {
            let scale = factor(row.omega_rad_s)?;
            for value in row.coefficients_by_input.iter_mut().flatten() { *value = *value * scale; }
        }
        for row in &mut converted.held_out {
            let scale = factor(row.omega_rad_s)?;
            for value in row.far_field_by_input.iter_mut().flatten() { *value = *value * scale; }
        }
        Ok(Self { bank: build_broadband_radiation_artifact(&converted, controls)?,
            source_radius_m, sound_speed_m_s })
    }
    /// Inspect the actual transformed acceleration-response fit and its limits.
    #[must_use]
    pub const fn acceleration_bank(&self) -> &BroadbandRadiationArtifact { &self.bank }
    /// Geometric time shift inserted into the fitted filter [s].
    #[must_use]
    pub fn source_time_shift_s(&self) -> f64 { self.source_radius_m / self.sound_speed_m_s }
}

/// One radiation input as a linear projection of the accepted mechanical velocity.
#[derive(Clone, Debug)]
pub struct VelocityProjection {
    /// Must exactly match the same-index acoustic input ID; never inferred.
    pub input_id: String,
    /// Mechanical-mode order, excluding internal felt recovery coordinates.
    /// Identity selects generalized velocity; shape weights give surface velocity.
    /// The acoustic sample producer must use the same units/normalization.
    pub weights: Vec<f64>,
}
/// Stationary listener relative to the acoustic sample producer's origin/frame.
#[derive(Clone, Copy, Debug)]
pub struct ImpactListener {
    /// Position [m]. Must lie outside the source's enclosing sphere. This check
    /// is NOT a far-field accuracy certificate; the caller must justify range.
    pub position_m: [f64; 3],
    /// Maximum added propagation rounding delay [s]; zero requires exact transit.
    pub maximum_delay_error_s: f64,
    /// Maximum delay-line payload in samples, admitted before allocation.
    pub maximum_delay_samples: usize,
    /// Physical pressure ceiling [Pa], a refusal threshold, never a gain control.
    pub maximum_abs_pressure_pa: f64,
}
/// Explicit observation timing on the mechanical clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpactObservationTiming {
    /// Source radius/speed already included in the fitted bank [s].
    pub source_time_shift_s: f64,
    /// Remaining geometrical propagation time before rounding [s].
    pub remaining_travel_time_s: f64,
    /// Delay samples actually retained after the bank.
    pub propagation_samples: usize,
    /// Additional nonnegative delay due to integer propagation [s].
    pub propagation_rounding_s: f64,
    /// Backward velocity difference's half-step phase delay [s].
    pub velocity_difference_delay_s: f64,
}

/// Nonlinear mechanics -> causal radiation -> geometric travel -> pressure [Pa].
/// Implements the existing stream surface; no second WAV encoder or scheduler.
/// Attach at sample zero with explicitly quiescent prior acoustics. Retaining
/// this object retains mechanical, felt, radiation and propagation histories.
/// A failed callback poisons this host because a physical prefix may have run.
/// Callers must discard the failed callback's output, as PressureRenderer requires.
pub struct ImpactPressureRenderer<'a> {
    mechanics: ImpactSystem,
    radiation: BroadbandRadiationRuntime<'a>,
    projections: Vec<VelocityProjection>,
    forces: Vec<f64>,
    previous_velocity: Vec<f64>,
    velocity: Vec<f64>,
    acceleration: Vec<f64>,
    directional_weights: Vec<f64>,
    travel: Vec<f64>,
    head: usize,
    timing: ImpactObservationTiming,
    sample_rate_hz: u32,
    max_block: usize,
    pressure_limit: f64,
    completed: u64,
    poisoned: bool,
    gate: CancelGate,
    last_frame: Option<ImpactFrame>,
}
fn input(what: &'static str) -> RenderError { RenderError::Control { what } }
fn owner(error: impl core::fmt::Display) -> RenderError {
    RenderError::Voice(crate::acoustic_realize::AcousticRealizeError::Nonlinear(error.to_string()))
}
impl<'a> ImpactPressureRenderer<'a> {
    /// Bind the exact mechanical and radiation clocks, coordinate IDs and fixed
    /// listener. Forces are held until explicitly replaced between callbacks.
    /// Startup radiating velocity must be zero: otherwise its acoustic history
    /// is missing. A moving nonradiating striker remains fully legal.
    pub fn new(
        mechanics: ImpactSystem, radiation: &'a ImpactRadiation,
        projections: Vec<VelocityProjection>, forces: Vec<f64>,
        listener: ImpactListener, sample_rate_hz: u32, max_block: usize,
    ) -> Result<Self, RenderError> {
        let dt = 1.0 / f64::from(sample_rate_hz);
        if sample_rate_hz == 0 || max_block == 0 || mechanics.sample != 0
            || dt.to_bits() != mechanics.config.dt_s.to_bits()
            || dt.to_bits() != radiation.bank.sample_interval_s.to_bits()
        { return Err(input("impact/radiation clock mismatch, nonzero start or invalid block capacity")); }
        if projections.len() != radiation.bank.inputs.len()
            || projections.iter().zip(&radiation.bank.inputs).any(|(p, a)|
                p.input_id != a.id || p.weights.len() != mechanics.modes || p.weights.iter().any(|v| !v.is_finite()))
            || forces.len() != mechanics.modes
            || forces.iter().any(|f| !f.is_finite() || f.abs() > mechanics.config.maximum_generalized_force)
        { return Err(input("impact radiation input identities, projections or forces disagree")); }
        let range = det::hypot(det::hypot(listener.position_m[0], listener.position_m[1]), listener.position_m[2]);
        if !range.is_finite() || range <= radiation.source_radius_m
            || !listener.maximum_delay_error_s.is_finite() || listener.maximum_delay_error_s < 0.0
            || !listener.maximum_abs_pressure_pa.is_finite() || listener.maximum_abs_pressure_pa <= 0.0
        { return Err(input("listener needs finite exterior range, nonnegative delay tolerance and positive pressure limit")); }
        let time = (range - radiation.source_radius_m) / radiation.sound_speed_m_s;
        let count = (time / dt).ceil();
        if !time.is_finite() || !count.is_finite() || count < 0.0
            || count >= usize::MAX as f64 || count > listener.maximum_delay_samples as f64
        { return Err(input("listener propagation exceeds representable delay budget")); }
        let delay = count as usize;
        let rounding = count * dt - time;
        if rounding < 0.0 || rounding > listener.maximum_delay_error_s {
            return Err(input("listener propagation rounding exceeds its declared time tolerance"));
        }
        let mut previous_velocity = Vec::with_capacity(projections.len());
        for p in &projections {
            let v: f64 = p.weights.iter().enumerate().map(|(i,b)| b * mechanics.x[2*i+1]).sum();
            if !v.is_finite() || v != 0.0 { return Err(input("nonzero initial radiating velocity needs an acoustic-state initialization")); }
            previous_velocity.push(v);
        }
        let mut travel = Vec::new();
        travel.try_reserve_exact(delay).map_err(|_| input("pressure propagation allocation failed"))?;
        travel.resize(delay, 0.0);
        // Cold-only SH evaluation. No spherical harmonics or allocation in the
        // observation loop. Keep the bank's declared real-channel order.
        let n = radiation.bank.channels.len();
        let mut unit = vec![C64::ZERO; n];
        let mut directional_weights = Vec::with_capacity(n);
        for i in 0..n {
            unit[i] = C64::ONE;
            let value = evaluate_real_tesseral(radiation.bank.l_max, &unit, listener.position_m).map_err(owner)?.re / range;
            if !value.is_finite() { return Err(input("directional pressure weight overflow")); }
            directional_weights.push(value); unit[i] = C64::ZERO;
        }
        let runtime = radiation.bank.try_runtime().map_err(owner)?;
        let timing = ImpactObservationTiming { source_time_shift_s: radiation.source_time_shift_s(),
            remaining_travel_time_s: time, propagation_samples: delay, propagation_rounding_s: rounding,
            velocity_difference_delay_s: 0.5 * dt };
        Ok(Self { mechanics, radiation: runtime, velocity: vec![0.0; projections.len()],
            acceleration: vec![0.0; projections.len()], projections, forces, previous_velocity,
            directional_weights, travel, head: 0, timing, sample_rate_hz, max_block,
            pressure_limit: listener.maximum_abs_pressure_pa, completed: 0, poisoned: false,
            gate: CancelGate::new_clock_free(), last_frame: None })
    }
    /// Accepted physical mechanics. No mutation bypasses acoustic histories.
    #[must_use]
    pub const fn mechanics(&self) -> &ImpactSystem { &self.mechanics }
    /// The exact travel and observation approximation actually used.
    #[must_use]
    pub const fn timing(&self) -> ImpactObservationTiming { self.timing }
    /// Last accepted mechanical sample; not a fluid-energy or radiation audit.
    #[must_use]
    pub const fn last_mechanical_frame(&self) -> Option<&ImpactFrame> { self.last_frame.as_ref() }
    /// Atomically replace held physical forces without erasing vibration/history.
    pub fn set_forces(&mut self, forces: &[f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if forces.len() != self.forces.len() || forces.iter().any(|f|
            !f.is_finite() || f.abs() > self.mechanics.config.maximum_generalized_force)
        { return Err(input("held impact force shape or physical limit failed")); }
        self.forces.copy_from_slice(forces); Ok(())
    }
    fn sample(&mut self) -> Result<f64, RenderError> {
        let frame = self.mechanics.step(&self.forces, &self.gate).map_err(owner)?;
        for (k, p) in self.projections.iter().enumerate() {
            self.velocity[k] = p.weights.iter().enumerate().map(|(i,b)| b * self.mechanics.x[2*i+1]).sum();
            self.acceleration[k] = (self.velocity[k] - self.previous_velocity[k]) / self.mechanics.config.dt_s;
            if !self.acceleration[k].is_finite() { return Err(input("radiating acceleration left finite set")); }
        }
        let coefficients = self.radiation.step(&self.acceleration).map_err(owner)?;
        let emitted: f64 = coefficients.iter().zip(&self.directional_weights).map(|(a,b)|a*b).sum();
        if !emitted.is_finite() || emitted.abs() > self.pressure_limit {
            return Err(input("physical pressure exceeded its declared ceiling"));
        }
        let pressure = if self.travel.is_empty() { emitted } else {
            let old = self.travel[self.head]; self.travel[self.head] = emitted;
            self.head = (self.head + 1) % self.travel.len(); old
        };
        self.previous_velocity.copy_from_slice(&self.velocity);
        self.last_frame = Some(frame); self.completed += 1; Ok(pressure)
    }
}
impl PressureRenderer for ImpactPressureRenderer<'_> {
    fn samples_rendered(&self) -> u64 { self.completed }
    fn max_block_len(&self) -> usize { self.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if rate != self.sample_rate_hz { return Err(input("sink and physical pressure clocks differ")); } Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if self.completed.checked_add(samples).is_none_or(|n| n > self.mechanics.config.max_steps) {
            return Err(input("pressure request exceeds remaining mechanical horizon"));
        } Ok(())
    }
    fn block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if out.is_empty() { return Err(RenderError::EmptyBlock); }
        if out.len() > self.max_block { return Err(input("pressure callback exceeds capacity")); }
        self.validate_sample_count(out.len() as u64)?;
        for value in out { match self.sample() {
            Ok(p) => *value = p,
            Err(e) => { self.poisoned = true; return Err(e); }
        }} Ok(())
    }
}

#[cfg(test)]
mod tests;
