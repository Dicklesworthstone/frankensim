//! Observe independent physical systems on one latency-aligned pressure clock.
//!
//! Each source keeps its own mechanics, sample-timed inputs and existing causal
//! decimator. Pure integer alignment retains the original exact-sample path.
//! Explicitly band-limited alignment also admits fractional receiver latency
//! on the output clock, using the shared causal fractional-delay owner and a
//! common retained reconstruction delay. Physical travel times stay in each
//! receiver. Nothing resamples forces or physical state, invents a
//! gain, or interconnects the mechanics. Sources must describe the SAME observer
//! and physical time origin; type/clock admission cannot establish that fact.
//! Three or more signed pressures accumulate in double-double precision and
//! round once after every source has contributed, retaining weak fields at
//! destructive nulls. One/two-source mixes need only the original binary64 sum.
//!
//! A new mix starts at output sample zero with zero pressure prehistory. Retain
//! the whole object to continue: reconstructing delay lines over running parts
//! would lose their histories and is refused. The common causal latency remains
//! in output. The finite horizon includes startup latency; no synthetic tail is
//! flushed and delayed pressure beyond that horizon remains unobserved.

use super::{DecimatedRenderer, DecimationInfo, PressureRenderer};
use crate::pcm_wav::{baffled::path::HALF, delay::FractionalDelay};
use crate::render::{RenderError, schedule::ScheduledRenderer};
use fs_math::dd::Dd;

/// Bounded per-part alignment storage, in final output samples. Source receiver
/// latency plus the outer decimator and any reconstruction must fit this domain.
/// 8192 is the n/4 pre-arrival latency of a 32768-point baffled reed outlet
/// grid, the largest a site card converges on: at most 64 KiB of history per part.
pub const MAX_TOTAL_FILTER_DELAY_OUTPUT_SAMPLES: usize = 8192;

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
    /// Whole part of the already-applied source receiver latency, on the final
    /// output clock. Add `source_filter_delay_fraction_samples` for exact timing.
    pub source_filter_delay_output_samples: usize,
    /// Fractional part of that existing receiver latency, in [0, 1).
    pub source_filter_delay_fraction_samples: f64,
    /// Whole source receiver plus outer decimator latency, before alignment.
    /// The source fraction above remains additional; no fraction is rounded up.
    pub total_filter_delay_output_samples: usize,
    /// Whole part of extra padding on the OUTPUT clock, never the mechanics.
    pub alignment_delay_samples: usize,
    /// Fractional padding on the output clock, in [0, 1).
    pub alignment_delay_fraction_samples: f64,
    /// Additional common reconstruction latency: 32 output samples when any
    /// part needs fractional alignment, otherwise zero. Never compensated.
    pub alignment_filter_delay_samples: usize,
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

enum Alignment {
    Integer(Vec<DelayLine>),
    Bandlimited(Vec<FractionalDelay>),
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
    delays: Alignment,
    scratch: Vec<f64>,
    sums: Vec<Dd>,
    config: PressureEnsembleConfig,
    common_delay: usize,
    band_limit_hz: Option<f64>,
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
        Self::build(parts, config, None, None)
    }

    /// Align existing receiver filters as well as the outer decimators.
    /// Each supplied delay is ALREADY present in that source and is measured
    /// in final output samples, not mechanical samples or seconds. It excludes
    /// physical source-to-listener propagation. No filter is rebuilt or removed.
    ///
    /// There must be one finite, nonnegative, integral delay per part. Each
    /// combined latency is bounded by MAX_TOTAL_FILTER_DELAY_OUTPUT_SAMPLES;
    /// fractional output delays refuse instead of being rounded or resampled.
    pub fn new_with_source_delays(
        parts: Vec<DecimatedRenderer<S>>,
        config: PressureEnsembleConfig,
        source_delay_output_samples: &[f64],
    ) -> Result<Self, RenderError> {
        Self::build(parts, config, Some(source_delay_output_samples), None)
    }

    /// Align existing receiver latencies that become fractional after rate
    /// conversion, retaining their physical phases and source clocks.
    ///
    /// The declared pressure band must already be resolved on the output clock
    /// and lie at or below 0.45 of that rate. This declaration is not a filter or
    /// proof that the source contains no higher frequencies. The existing rate
    /// conversion must suppress those before this fractional reconstruction.
    /// Delays describe numerical filters ONLY, excluding physical travel.
    ///
    /// A fractional request pads each original latency to the next whole sample
    /// at/after their maximum, then retains the same 32-sample reconstruction
    /// latency on every part. All-integral requests use the original integer
    /// path bit for bit. Complete history and final latency remain bounded by
    /// MAX_TOTAL_FILTER_DELAY_OUTPUT_SAMPLES; no future sample or tail flush.
    pub fn new_with_bandlimited_source_delays(
        parts: Vec<DecimatedRenderer<S>>,
        config: PressureEnsembleConfig,
        source_delay_output_samples: &[f64],
        maximum_pressure_frequency_hz: f64,
    ) -> Result<Self, RenderError> {
        Self::build(parts, config, Some(source_delay_output_samples),
            Some(maximum_pressure_frequency_hz))
    }

    fn build(
        parts: Vec<DecimatedRenderer<S>>,
        config: PressureEnsembleConfig,
        source_delays: Option<&[f64]>,
        band_limit_hz: Option<f64>,
    ) -> Result<Self, RenderError> {
        if parts.is_empty() || parts.len() > config.max_parts
            || config.max_block == 0 || config.sample_rate_hz == 0
        {
            return Err(sizing("pressure ensemble needs parts, a positive clock and explicit capacities"));
        }
        if source_delays.is_some_and(|delays| delays.len() != parts.len()) {
            return Err(sizing("pressure alignment needs one declared source filter delay per part"));
        }
        if band_limit_hz.is_some_and(|frequency| !frequency.is_finite() || frequency <= 0.0
            || frequency > 0.45 * f64::from(config.sample_rate_hz))
        {
            return Err(sizing("fractional pressure alignment needs a resolved positive band at or below 0.45 of the output rate"));
        }
        let mut common_delay = 0_usize;
        let mut fractional = false;
        for (index, part) in parts.iter().enumerate() {
            part.validate_sample_rate(config.sample_rate_hz)?;
            part.validate_sample_count(config.samples)?;
            if part.samples_rendered() != 0 || config.max_block > part.max_block_len() {
                return Err(sizing("pressure alignment needs sample-zero parts with sufficient callback capacity"));
            }
            let delay = part.info().delay_output_samples;
            if !delay.is_finite() || !(0.0..=80.0).contains(&delay) || delay.fract() != 0.0 {
                return Err(sizing("pressure alignment requires an admitted integral decimator latency"));
            }
            let source_delay = source_delays.map_or(0.0, |delays| delays[index]);
            if !source_delay.is_finite() || source_delay < 0.0
                || band_limit_hz.is_none() && source_delay.fract() != 0.0
                || delay + source_delay > MAX_TOTAL_FILTER_DELAY_OUTPUT_SAMPLES as f64
            {
                return Err(sizing("pressure alignment needs finite nonnegative source delays within its admitted output-clock policy"));
            }
            common_delay = common_delay.max((delay + source_delay).ceil() as usize);
            fractional |= source_delay.fract() != 0.0;
        }
        let reconstruction = if fractional { HALF } else { 0 };
        if common_delay + reconstruction > MAX_TOTAL_FILTER_DELAY_OUTPUT_SAMPLES {
            return Err(sizing("fractional pressure reconstruction exceeds the total latency budget"));
        }
        let mut info = Vec::new();
        let mut delays = Vec::new();
        let mut fractional_delays = Vec::new();
        let mut scratch = Vec::new();
        let mut sums = Vec::new();
        info.try_reserve_exact(parts.len()).map_err(|_| sizing("cannot reserve pressure part descriptions"))?;
        if fractional {
            fractional_delays.try_reserve_exact(parts.len()).map_err(|_| sizing("cannot reserve fractional pressure alignment lines"))?;
        } else {
            delays.try_reserve_exact(parts.len()).map_err(|_| sizing("cannot reserve pressure alignment lines"))?;
        }
        scratch.try_reserve_exact(config.max_block).map_err(|_| sizing("cannot reserve pressure mixing scratch"))?;
        scratch.resize(config.max_block, 0.0);
        let coherent_frames = if parts.len() > 2 { config.max_block } else { 0 };
        sums.try_reserve_exact(coherent_frames).map_err(|_| sizing("cannot reserve coherent pressure sums"))?;
        sums.resize(coherent_frames, Dd::ZERO);
        for (index, part) in parts.iter().enumerate() {
            let observation = part.info();
            let source_delay = source_delays.map_or(0.0, |delays| delays[index]);
            let source_filter_delay_output_samples = source_delay.floor() as usize;
            let source_filter_delay_fraction_samples = source_delay - source_filter_delay_output_samples as f64;
            let total_filter_delay_output_samples = source_filter_delay_output_samples
                + observation.delay_output_samples as usize;
            let padding = common_delay as f64 - (source_delay + observation.delay_output_samples);
            let alignment_delay_samples = padding.floor() as usize;
            let alignment_delay_fraction_samples = padding - alignment_delay_samples as f64;
            info.push(AlignedPartInfo { observation, source_filter_delay_output_samples,
                source_filter_delay_fraction_samples, total_filter_delay_output_samples,
                alignment_delay_samples, alignment_delay_fraction_samples,
                alignment_filter_delay_samples: reconstruction });
            if fractional {
                fractional_delays.push(FractionalDelay::new(padding,
                    MAX_TOTAL_FILTER_DELAY_OUTPUT_SAMPLES).map_err(sizing)?);
            } else {
                delays.push(DelayLine::new(alignment_delay_samples)?);
            }
        }
        let delays = if fractional { Alignment::Bandlimited(fractional_delays) }
            else { Alignment::Integer(delays) };
        Ok(Self { parts, info, delays, scratch, sums, config,
            common_delay: common_delay + reconstruction, band_limit_hz, completed: 0, poisoned: false })
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
    /// Explicitly admitted pressure band for fractional alignment; the original
    /// integer-only constructors impose no new bandwidth restriction.
    #[must_use]
    pub const fn alignment_band_limit_hz(&self) -> Option<f64> { self.band_limit_hz }
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
        if !self.sums.is_empty() { self.sums[..output.len()].fill(Dd::ZERO); }
        for (index, part) in self.parts.iter_mut().enumerate() {
            let scratch = &mut self.scratch[..output.len()];
            if let Err(error) = part.block(scratch) {
                self.poisoned = true;
                return Err(error);
            }
            let result = match &mut self.delays {
                Alignment::Integer(lines) => accumulate(scratch, output, &mut self.sums,
                    |sample| lines[index].push(sample)),
                Alignment::Bandlimited(lines) => accumulate(scratch, output, &mut self.sums,
                    |sample| lines[index].step(sample)),
            };
            if let Err(error) = result {
                self.poisoned = true;
                return Err(error);
            }
        }
        for (sample, sum) in output.iter_mut().zip(&self.sums) {
            *sample = sum.to_f64();
        }
        self.completed += samples;
        Ok(())
    }
}

// Select the prepared delay outside the sample loop. This retains the original
// integer fast path and compensated multi-source sum without a virtual callback.
fn accumulate(samples: &[f64], output: &mut [f64], sums: &mut [Dd],
    mut delay: impl FnMut(f64) -> f64) -> Result<(), RenderError>
{
    if sums.is_empty() {
        for (sum, &sample) in output.iter_mut().zip(samples) {
            if !sample.is_finite() { return Err(sizing("a pressure part produced a non-finite observation")); }
            *sum += delay(sample);
            if !sum.is_finite() { return Err(sizing("physical pressure superposition overflowed")); }
        }
    } else {
        for (sum, &sample) in sums[..output.len()].iter_mut().zip(samples) {
            if !sample.is_finite() { return Err(sizing("a pressure part produced a non-finite observation")); }
            *sum = *sum + Dd::from_f64(delay(sample));
            if !sum.hi.is_finite() || !sum.lo.is_finite() {
                return Err(sizing("physical pressure superposition overflowed"));
            }
        }
    }
    Ok(())
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

#[cfg(test)]
#[path = "ensemble/fractional_tests.rs"]
mod fractional_tests;
