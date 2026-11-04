//! Simultaneous original physical receivers and authored directional readouts.
//!
//! One retained realization owns the shared real/proper source draws. Every
//! output channel uses that realization and the same absolute sample clock.

use crate::pcm_wav::{MAX_WAV_CHANNELS, observation::PressureFrameRenderer};
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::observation::{
    HarmonicSampleBudget, HarmonicSampleClock, PreparedHarmonicField,
    array::{HarmonicRealizationBudget, PreparedHarmonicArrayRealization, PreparedHarmonicReadout},
};
use fs_math::dd::Dd;

/// An original physical receiver, or an already authored pressure-equivalent
/// microphone/beamformer. Readouts retain their actual p/u projection, signed
/// weights and causal delays; those operations are never applied a second time.
#[derive(Clone, Copy, Debug)]
pub enum HarmonicOutputSource<'a> {
    Receiver(usize),
    Readout(&'a PreparedHarmonicReadout<'a>),
}

/// Stable identity and physical source of one output channel, in file order.
/// IDs are unique, nonempty, at most 128 UTF-8 bytes, with no control characters.
/// Retain this mapping beside exported WAV bytes; no speaker layout is inferred.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicOutputChannel<'a> {
    pub id: &'a str,
    pub source: HarmonicOutputSource<'a>,
}

enum ChannelField<'a, const D: usize> {
    Receiver(&'a PreparedHarmonicField<D>),
    Readout(PreparedHarmonicField<0>),
}

impl<const D: usize> ChannelField<'_, D> {
    fn frequency_count(&self) -> usize {
        match self {
            Self::Receiver(field) => field.frequency_count(),
            Self::Readout(field) => field.frequency_count(),
        }
    }

    fn validate_range(
        &self,
        clock: HarmonicSampleClock,
        first: u64,
        count: u64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        match self {
            Self::Receiver(field) => {
                field.validate_sample_range_cancellable(clock, first, count, cancelled)
            }
            Self::Readout(field) => {
                field.validate_sample_range_cancellable(clock, first, count, cancelled)
            }
        }
    }

    fn pressure(&self, time: Dd, cancelled: &mut impl FnMut() -> bool) -> Result<f64, AeroacError> {
        match self {
            Self::Receiver(field) => field.pressure_at_clock_dd_cancellable(time, cancelled),
            Self::Readout(field) => field.pressure_at_clock_dd_cancellable(time, cancelled),
        }
    }
}

/// Bounded simultaneous pressure output of one coherent array realization.
///
/// Receivers are borrowed; readouts are realized once with the SAME named
/// coefficient draws. All channels retain their relative signs, phases and Pa
/// amplitudes. Samples are frame-interleaved and share one absolute integer/DD
/// clock. They are settled harmonic continuation, not a transient onset or a
/// newly inferred stochastic process. No source randomness, normalization,
/// pressure-to-velocity approximation or propagation is introduced.
///
/// Construction retains bounded channel fields and one callback buffer. Every
/// callback allocates nothing and publishes only after all channels succeed.
/// Cancellation/refusal leaves the caller's buffer and completed frame count
/// unchanged, so this stateless physical observer can safely retry.
pub struct HarmonicPressureFrameRenderer<'a, const D: usize> {
    realization: &'a PreparedHarmonicArrayRealization<'a, D>,
    channels: Vec<HarmonicOutputChannel<'a>>,
    fields: Vec<ChannelField<'a, D>>,
    clock: HarmonicSampleClock,
    budget: HarmonicSampleBudget,
    scratch: Vec<f64>,
    completed: u64,
}

impl<'a, const D: usize> HarmonicPressureFrameRenderer<'a, D> {
    /// Bind named channels to one original realization and sample clock.
    ///
    /// The realization budget bounds selected channels, shared factors, total
    /// retained receiver/readout lines, and the SUM of squared readout line
    /// counts before coalescing. Borrowed receiver fields need no coalescing.
    /// The sample budget counts SCALAR output values across all channels and
    /// must contain an integer number of frames. Tone work is frame count times
    /// the sum of distinct physical frequencies in the selected channels.
    pub fn new(
        realization: &'a PreparedHarmonicArrayRealization<'a, D>,
        channels: &[HarmonicOutputChannel<'a>],
        clock: HarmonicSampleClock,
        preparation: HarmonicRealizationBudget,
        sampling: HarmonicSampleBudget,
    ) -> Result<Self, RenderError> {
        Self::new_cancellable(realization, channels, clock, preparation, sampling, || {
            false
        })
    }

    /// Whole-request admission and bounded cancellation before publication.
    pub fn new_cancellable(
        realization: &'a PreparedHarmonicArrayRealization<'a, D>,
        channels: &[HarmonicOutputChannel<'a>],
        clock: HarmonicSampleClock,
        preparation: HarmonicRealizationBudget,
        sampling: HarmonicSampleBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        poll(&mut cancelled)?;
        let count = channels.len();
        if count == 0
            || count > usize::from(MAX_WAV_CHANNELS)
            || count > preparation.max_channels
            || realization.factor_draws().len() > preparation.max_factors
            || sampling.max_samples == 0
            || sampling.max_samples % count != 0
        {
            return Err(sizing(
                "harmonic output channel, factor or scalar callback budget",
            ));
        }
        let mut total_lines = 0usize;
        let mut pairs = 0usize;
        for (index, channel) in channels.iter().enumerate() {
            poll(&mut cancelled)?;
            if channel.id.trim().is_empty()
                || channel.id.len() > 128
                || channel.id.chars().any(char::is_control)
                || channels[..index].iter().any(|other| other.id == channel.id)
            {
                return Err(sizing(
                    "harmonic output channel IDs must be unique bounded nonempty text",
                ));
            }
            let lines = match channel.source {
                HarmonicOutputSource::Receiver(receiver) => realization
                    .receivers()
                    .get(receiver)
                    .ok_or_else(|| sizing("harmonic output receiver index"))?
                    .field
                    .frequency_count(),
                HarmonicOutputSource::Readout(readout) => {
                    if realization.source_clock() != readout.source_clock()
                        || realization.factor_coordinates() != readout.factor_coordinates()
                    {
                        return Err(sizing(
                            "harmonic output readout clock or factor coordinates differ",
                        ));
                    }
                    let lines = readout.ensemble().retained_line_count();
                    pairs = lines
                        .checked_mul(lines)
                        .and_then(|n| pairs.checked_add(n))
                        .ok_or_else(|| sizing("harmonic output realization work overflow"))?;
                    lines
                }
            };
            total_lines = total_lines
                .checked_add(lines)
                .ok_or_else(|| sizing("harmonic output line count overflow"))?;
        }
        if total_lines > preparation.max_total_lines || pairs > preparation.max_line_pairs {
            return Err(sizing(
                "harmonic output exceeds whole-channel realization budget",
            ));
        }
        // Preflight every selected readout before allocating or realizing any
        // channel. Each readout uses the original physical coefficient owner.
        let mut fields = Vec::new();
        fields
            .try_reserve_exact(count)
            .map_err(|_| sizing("harmonic output channel allocation failed"))?;
        for channel in channels {
            poll(&mut cancelled)?;
            fields.push(match channel.source {
                HarmonicOutputSource::Receiver(index) => {
                    ChannelField::Receiver(&realization.receivers()[index].field)
                }
                HarmonicOutputSource::Readout(readout) => ChannelField::Readout(
                    realization
                        .realize_readout(readout, preparation, &mut cancelled)
                        .map_err(RenderError::Acoustic)?,
                ),
            });
        }
        let mut tones = 0usize;
        for field in &fields {
            tones = tones
                .checked_add(field.frequency_count())
                .ok_or_else(|| sizing("harmonic output tone count overflow"))?;
            field
                .validate_range(clock, 0, 0, &mut cancelled)
                .map_err(RenderError::Acoustic)?;
        }
        if (sampling.max_samples / count)
            .checked_mul(tones)
            .is_none_or(|work| work > sampling.max_tone_samples)
        {
            return Err(sizing("harmonic output callback exceeds tone work budget"));
        }
        let mut retained_channels = Vec::new();
        retained_channels
            .try_reserve_exact(count)
            .map_err(|_| sizing("harmonic output identity allocation failed"))?;
        retained_channels.extend_from_slice(channels);
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(sampling.max_samples)
            .map_err(|_| sizing("harmonic output callback allocation failed"))?;
        scratch.resize(sampling.max_samples, 0.0);
        poll(&mut cancelled)?;
        Ok(Self {
            realization,
            channels: retained_channels,
            fields,
            clock,
            budget: sampling,
            scratch,
            completed: 0,
        })
    }

    pub const fn realization(&self) -> &PreparedHarmonicArrayRealization<'a, D> {
        self.realization
    }

    pub fn channels(&self) -> &[HarmonicOutputChannel<'a>] {
        &self.channels
    }

    pub const fn clock(&self) -> HarmonicSampleClock {
        self.clock
    }

    /// Atomically publish one complete frame block. The original field owner
    /// polls through frequencies, including long per-channel reductions.
    pub fn frame_block_cancellable(
        &mut self,
        output: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        poll(&mut cancelled)?;
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        let channels = self.channels.len();
        if output.len() % channels != 0 || output.len() > self.budget.max_samples {
            return Err(sizing(
                "harmonic output needs complete frames within callback capacity",
            ));
        }
        let count = u64::try_from(output.len() / channels)
            .map_err(|_| sizing("harmonic output frame count exceeds u64"))?;
        for field in &self.fields {
            field
                .validate_range(self.clock, self.completed, count, &mut cancelled)
                .map_err(RenderError::Acoustic)?;
        }
        for (offset, frame) in self.scratch[..output.len()]
            .chunks_exact_mut(channels)
            .enumerate()
        {
            let time = self
                .clock
                .time_at(self.completed + offset as u64)
                .map_err(RenderError::Acoustic)?;
            for (sample, field) in frame.iter_mut().zip(&self.fields) {
                *sample = field
                    .pressure(time, &mut cancelled)
                    .map_err(RenderError::Acoustic)?;
            }
        }
        poll(&mut cancelled)?;
        output.copy_from_slice(&self.scratch[..output.len()]);
        self.completed += count;
        Ok(())
    }
}

impl<const D: usize> PressureFrameRenderer for HarmonicPressureFrameRenderer<'_, D> {
    fn channel_count(&self) -> usize {
        self.channels.len()
    }
    fn frames_rendered(&self) -> u64 {
        self.completed
    }
    fn max_frame_block_len(&self) -> usize {
        self.budget.max_samples / self.channels.len()
    }

    fn validate_frame_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate != self.clock.sample_rate_hz() {
            return Err(sizing("harmonic output and WAV frame clocks differ"));
        }
        Ok(())
    }

    fn validate_frame_count(&self, frames: u64) -> Result<(), RenderError> {
        for field in &self.fields {
            field
                .validate_range(self.clock, self.completed, frames, &mut || false)
                .map_err(RenderError::Acoustic)?;
        }
        Ok(())
    }

    fn frame_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.frame_block_cancellable(output, || false)
    }
}

fn sizing(what: &'static str) -> RenderError {
    RenderError::Sizing { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RenderError> {
    if cancelled() {
        Err(RenderError::Acoustic(AeroacError::Cancelled {
            what: "harmonic pressure frame callback",
        }))
    } else {
        Ok(())
    }
}
