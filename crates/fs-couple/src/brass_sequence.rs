//! Sample-indexed performances on the original multimodal lip/cup/valve voice.
//!
//! Events carry the existing [`BrassControl`] values, not notes or oscillator
//! frequencies. Callback boundaries do not move an event. The output is the
//! actual MOUTHPIECE pressure in Pa; it is not an exterior microphone. Optional
//! outgoing-wave traces expose the same incident pressure used by downstream
//! acoustic observers, without synthesizing a radiation law here.

use crate::acoustic_realize::AcousticRealizeError;
use crate::brass_loop::{BrassControl, BrassVoice, LipIsland};
use crate::mm_line::{MmBankLoad, MmLineConfig};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use fs_duct::Duct;
use fs_material::gas::GasState;

/// Every original mouth mode observed from the same loaded bore and valve history.
pub mod mouth;

/// Maximum number of authored control assignments in a sequence.
pub const MAX_BRASS_EVENTS: usize = 16_384;
/// Maximum number of independently realized valve paths.
pub const MAX_BRASS_COMBOS: usize = 16;
const MAX_BLOCK: usize = 65_536;
const MAX_SECONDS: u64 = 600;

/// An assignment applied immediately BEFORE its zero-based source sample.
#[derive(Clone, Copy, Debug)]
pub struct BrassEvent {
    /// Source-clock sample index, not a callback or output-decimator index.
    pub sample: u64,
    /// Original physical control and valve-history lift.
    pub control: BrassControl,
}

/// Explicit finite horizon, callback capacity and initial physical controls.
#[derive(Clone, Copy, Debug)]
pub struct BrassSequenceConfig {
    /// Total source-clock samples, at most 600 seconds.
    pub samples: u64,
    /// Largest admitted callback, in 1..=65536 samples.
    pub max_block_len: usize,
    /// Initial nonnegative mouth pressure [Pa].
    pub initial_blowing_pressure_pa: f64,
    /// Initial positive dimensionless lip pre-stress multiplier.
    pub initial_lip_tension: f64,
    /// Initial original valve-path index.
    pub initial_combo: usize,
}

/// A continuing physical voice with immutable, pre-admitted control events.
///
/// All paths share the original sample clock, medium and input impedance.
/// Each may carry its own referenced strict table through [`MmBankLoad`].
/// No callback pre-renders future samples and no valve event resets mechanics.
/// The voice's existing per-chunk diagnostics remain callback-dependent;
/// pressure, outgoing waves and event placement do not. The stream writer owns
/// cancellation boundaries. This wrapper adds no hard-real-time claim.
pub struct BrassSequence {
    voice: BrassVoice,
    config: BrassSequenceConfig,
    rate_hz: u32,
    events: Vec<BrassEvent>,
    next_event: usize,
    completed: u64,
    poisoned: bool,
}

fn invalid(what: &'static str) -> RenderError {
    RenderError::Control { what }
}
fn physical(error: impl core::fmt::Display) -> RenderError {
    RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string()))
}

fn validate_control(control: BrassControl, lip: &LipIsland, combos: usize)
    -> Result<u8, RenderError>
{
    match control {
        BrassControl::SetBlowingPressure(p) if p.is_finite() && p >= 0.0 => Ok(1),
        BrassControl::SetLipTension(t) if t.is_finite() && t > 0.0 => {
            let stiffness = lip.stiffness_n_m * t * t;
            let damping = lip.damping_n_s_m * t;
            if !stiffness.is_finite() || stiffness <= 0.0 || !damping.is_finite() {
                return Err(invalid("scheduled tension has unrepresentable lip coefficients"));
            }
            Ok(2)
        }
        BrassControl::SetValve { combo, .. } if combo < combos => Ok(4),
        _ => Err(invalid("invalid scheduled pressure, tension or valve identity")),
    }
}

fn validate_events(config: BrassSequenceConfig, rate_hz: u32, combos: usize,
    lip: &LipIsland, events: &[BrassEvent]) -> Result<(), RenderError>
{
    if rate_hz == 0 || config.samples == 0
        || config.samples > MAX_SECONDS * u64::from(rate_hz)
        || !(1..=MAX_BLOCK).contains(&config.max_block_len)
        || !(1..=MAX_BRASS_COMBOS).contains(&combos) || events.len() > MAX_BRASS_EVENTS
    {
        return Err(invalid("brass sequence exceeds its clock, horizon, callback, path or event limits"));
    }
    if ![lip.mass_kg, lip.stiffness_n_m, lip.width_m, lip.rest_gap_m, lip.face_area_m2]
        .iter().all(|x| x.is_finite() && *x > 0.0)
        || !lip.damping_n_s_m.is_finite() || lip.damping_n_s_m < 0.0
    {
        return Err(invalid("brass sequence requires finite physical lip parameters"));
    }
    for control in [BrassControl::SetBlowingPressure(config.initial_blowing_pressure_pa),
        BrassControl::SetLipTension(config.initial_lip_tension),
        BrassControl::SetValve { combo: config.initial_combo, fade_samples: 0 }]
    {
        validate_control(control, lip, combos)?;
    }
    let mut previous = None;
    let mut assigned = 0_u8;
    let mut active = config.initial_combo;
    let mut fade_end = 0;
    for event in events {
        if event.sample >= config.samples || previous.is_some_and(|sample| event.sample < sample) {
            return Err(invalid("brass events must be ordered and inside the observed horizon"));
        }
        if previous != Some(event.sample) { assigned = 0; }
        let kind = validate_control(event.control, lip, combos)?;
        if assigned & kind != 0 {
            return Err(invalid("a brass control is assigned twice at the same sample"));
        }
        assigned |= kind;
        previous = Some(event.sample);
        if let BrassControl::SetValve { combo, fade_samples } = event.control {
            // A no-op assignment retains the original bank's no-op semantics.
            if combo != active {
                let fade = u64::try_from(fade_samples)
                    .map_err(|_| invalid("valve fade is outside the source-clock range"))?;
                if event.sample < fade_end || fade > config.samples - event.sample {
                    return Err(invalid("valve fades must not overlap or extend past the sequence"));
                }
                fade_end = event.sample + fade;
                active = combo;
            }
        }
    }
    Ok(())
}

impl BrassSequence {
    /// Admit the complete schedule before realizing any acoustic line. The
    /// original bank validates all load references and builds one shared ladder.
    /// Simultaneous DIFFERENT controls follow authored order; duplicate controls
    /// and overlapping valve fades refuse rather than silently rescheduling.
    ///
    /// # Errors
    /// Schedule/size errors or the original acoustic/lip owner's refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn new<'a>(ducts: &[Duct], labels: &[&str], gas: &GasState,
        loads: impl Into<MmBankLoad<'a>>, line: &MmLineConfig, lip: LipIsland,
        cup_volume_m3: f64, config: BrassSequenceConfig, events: Vec<BrassEvent>)
        -> Result<Self, RenderError>
    {
        validate_events(config, line.sample_rate_hz, ducts.len(), &lip, &events)?;
        let mut voice = BrassVoice::new(ducts, labels, gas, loads, line, lip, cup_volume_m3)
            .map_err(physical)?;
        voice.apply(BrassControl::SetBlowingPressure(config.initial_blowing_pressure_pa))
            .map_err(physical)?;
        voice.apply(BrassControl::SetLipTension(config.initial_lip_tension)).map_err(physical)?;
        voice.apply(BrassControl::SetValve { combo: config.initial_combo, fade_samples: 0 })
            .map_err(physical)?;
        Ok(Self { voice, config, rate_hz: line.sample_rate_hz, events,
            next_event: 0, completed: 0, poisoned: false })
    }

    /// Original voice, immutable so external writes cannot bypass the schedule.
    #[must_use]
    pub const fn voice(&self) -> &BrassVoice { &self.voice }
    /// Original source-clock assignments in their authored order.
    #[must_use]
    pub fn events(&self) -> &[BrassEvent] { &self.events }
    /// Unrendered source samples in the admitted finite horizon.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.config.samples - self.completed }

    /// Render mouthpiece pressure and the SAME sample's outgoing inlet wave.
    /// Neither stream is a radiated microphone signal. Length/size refusals
    /// occur before advancement; a physical failure poisons continuation because
    /// a prefix of the original physical state may already have advanced.
    ///
    /// # Errors
    /// Mismatched buffers, exhausted horizon, invalid callback size or a
    /// refusal from the original physical voice.
    pub fn block_outgoing(&mut self, output: &mut [f64], outgoing: &mut [f64])
        -> Result<(), RenderError>
    {
        if output.len() != outgoing.len() {
            return Err(invalid("outgoing trace must match the mouthpiece callback length"));
        }
        self.render(output, Some(outgoing))
    }

    fn render(&mut self, output: &mut [f64], outgoing: Option<&mut [f64]>)
        -> Result<(), RenderError>
    {
        self.validate_sample_count(output.len() as u64)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block_len {
            return Err(RenderError::Sizing { what: "brass sequence callback exceeds its capacity" });
        }
        if let Err(error) = self.advance(output, outgoing) {
            self.poisoned = true;
            return Err(error);
        }
        self.completed += output.len() as u64;
        Ok(())
    }

    fn advance(&mut self, output: &mut [f64], mut outgoing: Option<&mut [f64]>)
        -> Result<(), RenderError>
    {
        let mut done = 0;
        while done < output.len() {
            let sample = self.completed + done as u64;
            while let Some(event) = self.events.get(self.next_event) {
                if event.sample != sample { break; }
                self.voice.apply(event.control).map_err(physical)?;
                self.next_event += 1;
            }
            let next = self.events.get(self.next_event)
                .map_or(self.config.samples, |event| event.sample);
            let count = usize::try_from(next - sample).unwrap_or(usize::MAX)
                .min(output.len() - done);
            let end = done + count;
            if let Some(trace) = outgoing.as_deref_mut() {
                self.voice.step_block_outgoing(&mut output[done..end], &mut trace[done..end])
                    .map_err(physical)?;
                if trace[done..end].iter().any(|x| !x.is_finite()) {
                    return Err(invalid("brass outgoing pressure left the finite set"));
                }
            } else {
                self.voice.step_block(&mut output[done..end]).map_err(physical)?;
            }
            if output[done..end].iter().any(|x| !x.is_finite()) {
                return Err(invalid("brass mouthpiece pressure left the finite set"));
            }
            done = end;
        }
        Ok(())
    }
}

impl PressureRenderer for BrassSequence {
    fn samples_rendered(&self) -> u64 { self.completed }
    fn max_block_len(&self) -> usize { self.config.max_block_len }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if rate != self.rate_hz {
            return Err(RenderError::Sizing { what: "brass sequence output clock differs" });
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() {
            return Err(RenderError::Sizing { what: "brass sequence request exceeds its horizon" });
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.render(output, None)
    }
}

#[cfg(test)]
mod tests;
