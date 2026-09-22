//! Observe independent physical systems on one latency-aligned pressure clock.
//!
//! Each source keeps its own mechanics, sample-timed inputs and existing causal
//! decimator. A pure integer delay aligns the faster observation paths to the
//! largest filter delay. Nothing resamples forces or physical state, invents a
//! gain, or interconnects the mechanics. Sources must describe the SAME observer
//! and physical time origin; type/clock admission cannot establish that fact.
//!
//! A new mix starts at output sample zero with zero pressure prehistory. Retain
//! the whole object to continue: reconstructing delay lines over running parts
//! would lose their histories and is refused. The common causal latency remains
//! in output. The finite horizon includes startup latency; no synthetic tail is
//! flushed and delayed pressure beyond that horizon remains unobserved.

use super::{DecimatedRenderer, DecimationInfo, PressureRenderer};
use crate::render::{RenderError, schedule::ScheduledRenderer};

/// Finite output and storage limits, independent of each source's own budgets.
#[derive(Clone, Copy, Debug)]
pub struct PressureEnsembleConfig {
    /// Common OUTPUT pressure sample rate [Hz], not the rate of every solver.
    pub sample_rate_hz: u32,
    /// Additional output samples from the common time origin.
    pub samples: u64,
    /// Largest output callback. Every observed part must admit this size.
    pub max_block: usize,
    /// Maximum number of independent observed parts.
    pub max_parts: usize,
}

/// Exact observation policy for a part, in caller-supplied summation order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlignedPartInfo {
    /// The existing decimator's clock, filter profile and original latency.
    pub observation: DecimationInfo,
    /// Extra pure delay on the OUTPUT clock; never applied to the mechanics.
    pub alignment_delay_samples: usize,
}

struct DelayLine {
    history: Vec<f64>,
    head: usize,
}
impl DelayLine {
    fn new(samples: usize) -> Result<Self, RenderError> {
        let mut history = Vec::new();
        history.try_reserve_exact(samples).map_err(|_| sizing("cannot reserve pressure alignment history"))?;
        history.resize(samples, 0.0);
        Ok(Self { history, head: 0 })
    }
    fn push(&mut self, value: f64) -> f64 {
        if self.history.is_empty() { return value; }
        let delayed = self.history[self.head];
        self.history[self.head] = value;
        self.head = (self.head + 1) % self.history.len();
        delayed
    }
}

/// Mixed-rate mechanics observed through one finite physical-pressure stream.
///
/// Use `S = Box<dyn PressureRenderer>` for different physical producer types.
/// Allocation is confined to construction in this layer; source steppers retain
/// their disclosed allocation behavior. A failed callback permanently poisons
/// the mix, since another part or a filter history may already have advanced.
/// Read-only source access supports diagnosis, not partial-callback recovery.
pub struct PressureEnsemble<S = ScheduledRenderer> {
    parts: Vec<DecimatedRenderer<S>>,
    info: Vec<AlignedPartInfo>,
    delays: Vec<DelayLine>,
    scratch: Vec<f64>,
    config: PressureEnsembleConfig,
    common_delay: usize,
    completed: u64,
    poisoned: bool,
}
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }

impl<S: PressureRenderer> PressureEnsemble<S> {
    /// Admit all source clocks and the COMPLETE window before rendering.
    ///
    /// Integer rate conversion is chosen explicitly when constructing each
    /// `DecimatedRenderer`. The current shared profiles have integral group
    /// delays of at most 80 output samples. A different/fractional profile is
    /// refused rather than rounded or compensated by a second resampler.
    pub fn new(parts: Vec<DecimatedRenderer<S>>, config: PressureEnsembleConfig)
        -> Result<Self, RenderError>
    {
        if parts.is_empty() || parts.len() > config.max_parts
            || config.max_block == 0 || config.sample_rate_hz == 0
        {
            return Err(sizing("pressure ensemble needs parts, a positive clock and explicit capacities"));
        }
        let mut common_delay = 0_usize;
        for part in &parts {
            part.validate_sample_rate(config.sample_rate_hz)?;
            part.validate_sample_count(config.samples)?;
            if part.samples_rendered() != 0 || config.max_block > part.max_block_len() {
                return Err(sizing("pressure alignment needs sample-zero parts with sufficient callback capacity"));
            }
            let delay = part.info().delay_output_samples;
            if !delay.is_finite() || !(0.0..=80.0).contains(&delay) || delay.fract() != 0.0 {
                return Err(sizing("pressure alignment requires an admitted integral decimator latency"));
            }
            common_delay = common_delay.max(delay as usize);
        }
        let mut info = Vec::new();
        let mut delays = Vec::new();
        let mut scratch = Vec::new();
        info.try_reserve_exact(parts.len()).map_err(|_| sizing("cannot reserve pressure part descriptions"))?;
        delays.try_reserve_exact(parts.len()).map_err(|_| sizing("cannot reserve pressure alignment lines"))?;
        scratch.try_reserve_exact(config.max_block).map_err(|_| sizing("cannot reserve pressure mixing scratch"))?;
        scratch.resize(config.max_block, 0.0);
        for part in &parts {
            let observation = part.info();
            let alignment_delay_samples = common_delay - observation.delay_output_samples as usize;
            info.push(AlignedPartInfo { observation, alignment_delay_samples });
            delays.push(DelayLine::new(alignment_delay_samples)?);
        }
        Ok(Self { parts, info, delays, scratch, config, common_delay, completed: 0, poisoned: false })
    }

    /// Exact source/observation clocks and added delays in summation order.
    #[must_use]
    pub fn part_info(&self) -> &[AlignedPartInfo] { &self.info }
    /// Inspect retained solvers and decimator histories without advancing them.
    #[must_use]
    pub fn parts(&self) -> &[DecimatedRenderer<S>] { &self.parts }
    /// Common causal delay, retained in the output rather than removed.
    #[must_use]
    pub const fn delay_output_samples(&self) -> usize { self.common_delay }
    /// Completed pressure samples on the common output clock.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 { self.completed }
    /// Samples left in the admitted observation window.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 { self.config.samples - self.completed }
    /// Largest admitted output callback, in samples.
    #[must_use]
    pub const fn max_block_len(&self) -> usize { self.config.max_block }

    /// Render one pressure block using the original per-part steppers.
    /// Shape/window refusals leave state and output untouched. Physical or sum
    /// overflow failures poison the mix; the caller must discard that callback.
    pub fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.max_block_len() {
            return Err(sizing("pressure block exceeds ensemble callback capacity"));
        }
        let samples = u64::try_from(output.len()).map_err(|_| sizing("pressure block length exceeds u64"))?;
        self.validate_sample_count(samples)?;
        output.fill(0.0);
        for (part, delay) in self.parts.iter_mut().zip(&mut self.delays) {
            let scratch = &mut self.scratch[..output.len()];
            if let Err(error) = part.block(scratch) {
                self.poisoned = true;
                return Err(error);
            }
            for (sum, &sample) in output.iter_mut().zip(scratch.iter()) {
                if !sample.is_finite() {
                    self.poisoned = true;
                    return Err(sizing("a pressure part produced a non-finite observation"));
                }
                *sum += delay.push(sample);
                if !sum.is_finite() {
                    self.poisoned = true;
                    return Err(sizing("physical pressure superposition overflowed"));
                }
            }
        }
        self.completed += samples;
        Ok(())
    }
}

impl<S: PressureRenderer> PressureRenderer for PressureEnsemble<S> {
    fn samples_rendered(&self) -> u64 { Self::samples_rendered(self) }
    fn max_block_len(&self) -> usize { Self::max_block_len(self) }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if rate != self.config.sample_rate_hz {
            return Err(sizing("sink rate must match the pressure ensemble output clock"));
        }
        for part in &self.parts { part.validate_sample_rate(rate)?; }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() {
            return Err(sizing("pressure request exceeds the complete ensemble window"));
        }
        for part in &self.parts { part.validate_sample_count(samples)?; }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> { Self::block(self, output) }
}
