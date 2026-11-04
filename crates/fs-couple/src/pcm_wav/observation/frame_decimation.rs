//! Causal rate conversion for simultaneous physical pressure channels.
//!
//! All channels use the original Decimator's independently retained delay lines
//! and identical coefficients. Only pressure is filtered; each source's motion,
//! event times, constitutive laws, spatial phases and absolute Pa scale are intact.

use super::{DecimationInfo, PressureFrameRenderer, sizing};
use crate::pcm_wav::{MAX_WAV_CHANNELS, decimate::Decimator};
use crate::render::RenderError;

/// One native multichannel pressure stream on a lower integer output clock.
///
/// The admitted producer supplies complete frame-interleaved values. The SAME
/// original causal filter owns every channel, with no mixing, normalization,
/// inferred speaker layout, independent resampling phase, or delay compensation.
/// A single input group is retained irrespective of callback size or duration.
/// After construction, this adapter allocates nothing; source callback allocation
/// guarantees remain with the source. Filter size is bounded by 1..=16 rate ratio
/// and the existing decimator's 1..=128 channels.
///
/// Start at source frame zero with zero filter history, and retain this owner
/// for continuation. Source pressure must already resolve its actual Doppler and
/// modulation bandwidth. Decimation cannot repair aliasing at the source clock.
/// The original filter's passband/transition and startup delay remain explicit.
///
/// As with the mono DecimatedRenderer, a failed source/filter callback can have
/// advanced an internal prefix. It poisons the adapter: discard that callback's
/// output and do not resume by resetting a clock or discarding filter history.
/// Input-shape/rate/horizon admission failures occur before source advancement
/// and do not poison. The existing PCM writer polls cancellation between complete
/// callbacks and drains an accepted callback before stopping.
pub struct DecimatedFrameRenderer<S> {
    source: S,
    filter: Decimator,
    input: Vec<f64>,
    info: DecimationInfo,
    channels: usize,
    max_block_frames: usize,
    completed: u64,
    poisoned: bool,
}

impl<S: PressureFrameRenderer> DecimatedFrameRenderer<S> {
    /// Observe existing channels without rebuilding or retuning their mechanics.
    /// `max_block_frames` counts output FRAMES, not interleaved scalar values.
    /// Ratio one delegates the source's original callback without arithmetic.
    pub fn new(
        source: S,
        source_sample_rate_hz: u32,
        output_sample_rate_hz: u32,
        max_block_frames: usize,
    ) -> Result<Self, RenderError> {
        source.validate_frame_rate(source_sample_rate_hz)?;
        if output_sample_rate_hz == 0 || source_sample_rate_hz % output_sample_rate_hz != 0 {
            return Err(sizing("frame decimation requires a positive integer source/output rate ratio"));
        }
        let ratio = (source_sample_rate_hz / output_sample_rate_hz) as usize;
        let channels = source.channel_count();
        if !(1..=16).contains(&ratio)
            || channels == 0 || channels > usize::from(MAX_WAV_CHANNELS) || channels > 128
            || max_block_frames == 0 || source.max_frame_block_len() == 0
            || max_block_frames.checked_mul(channels).is_none()
        {
            return Err(sizing("frame decimation requires admitted channel, rate and callback counts"));
        }
        if source.frames_rendered() != 0 {
            return Err(sizing("frame decimation must start at source frame zero with complete filter history"));
        }
        if ratio == 1 && max_block_frames > source.max_frame_block_len() {
            return Err(sizing("bypass frame capacity must fit the source callback"));
        }
        let filter = Decimator::new(ratio, channels).map_err(|what| RenderError::Control { what })?;
        let info = DecimationInfo {
            mechanics_sample_rate_hz: source_sample_rate_hz,
            output_sample_rate_hz,
            ratio,
            delay_output_samples: filter.delay_output_frames(),
            first_output_source_index: ratio - 1,
            filter_profile: if ratio == 1 { "identity" } else { "blackman-harris-integer-decimator-v1" },
        };
        Ok(Self { source, filter, input: vec![0.0; ratio * channels], info,
            channels, max_block_frames, completed: 0, poisoned: false })
    }

    /// Original geometry/source state, without mutable access that bypasses the filter.
    pub const fn source(&self) -> &S { &self.source }
    /// Common physical delay and group sampling phase of ALL output channels.
    pub const fn info(&self) -> DecimationInfo { self.info }
    pub const fn channel_count(&self) -> usize { self.channels }
    pub const fn frames_rendered(&self) -> u64 { self.completed }
    pub const fn max_frame_block_len(&self) -> usize { self.max_block_frames }

    fn admit_source(&self) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        self.source.validate_frame_rate(self.info.mechanics_sample_rate_hz)?;
        if self.source.channel_count() != self.channels || self.source.max_frame_block_len() == 0 {
            return Err(sizing("original decimated source changed its channel layout or capacity"));
        }
        Ok(())
    }

    pub fn validate_frame_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.admit_source()?;
        if rate != self.info.output_sample_rate_hz {
            return Err(sizing("sink rate must match the declared decimated frame clock"));
        }
        Ok(())
    }

    /// Admit the entire expanded source-frame horizon before any callback.
    pub fn validate_frame_count(&self, frames: u64) -> Result<(), RenderError> {
        self.admit_source()?;
        self.completed.checked_add(frames)
            .ok_or_else(|| sizing("decimated frame clock overflows"))?;
        let source_frames = frames.checked_mul(self.info.ratio as u64)
            .ok_or_else(|| sizing("expanded source frame count overflows"))?;
        self.source.validate_frame_count(source_frames)
    }

    /// Convert a remaining source horizon with no dropping, padding or synthetic tail.
    pub fn output_frames_for(&self, source_frames: u64) -> Result<u64, RenderError> {
        let ratio = self.info.ratio as u64;
        if source_frames % ratio != 0 {
            return Err(sizing("source horizon must contain complete decimation frame groups"));
        }
        let frames = source_frames / ratio;
        self.validate_frame_count(frames)?;
        Ok(frames)
    }

    pub fn frame_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() % self.channels != 0 || output.len() / self.channels > self.max_block_frames {
            return Err(sizing("decimated output must contain complete frames within callback capacity"));
        }
        let count = u64::try_from(output.len() / self.channels)
            .map_err(|_| sizing("decimated output frame count exceeds u64"))?;
        self.validate_frame_count(count)?;
        if self.info.ratio == 1 {
            if output.len() / self.channels > self.source.max_frame_block_len() {
                return Err(sizing("bypass frame capacity no longer fits the source callback"));
            }
            if let Err(error) = self.source.frame_block(output) {
                self.poisoned = true;
                return Err(error);
            }
        } else {
            // Cap before multiplication, so even a source with an enormous
            // declared capacity cannot overflow this fixed small group shape.
            let source_values = self.source.max_frame_block_len().min(self.info.ratio) * self.channels;
            for frame in output.chunks_exact_mut(self.channels) {
                for chunk in self.input.chunks_mut(source_values) {
                    if let Err(error) = self.source.frame_block(chunk) {
                        self.poisoned = true;
                        return Err(error);
                    }
                }
                match self.filter.preview(&self.input) {
                    Ok(candidate) => frame.copy_from_slice(candidate),
                    Err(what) => {
                        self.poisoned = true;
                        return Err(RenderError::Control { what });
                    }
                }
                self.filter.commit();
            }
        }
        self.completed += count;
        Ok(())
    }
}

impl<S: PressureFrameRenderer> PressureFrameRenderer for DecimatedFrameRenderer<S> {
    fn channel_count(&self) -> usize { self.channels }
    fn frames_rendered(&self) -> u64 { self.completed }
    fn max_frame_block_len(&self) -> usize { self.max_block_frames }
    fn validate_frame_rate(&self, rate: u32) -> Result<(), RenderError> {
        Self::validate_frame_rate(self, rate)
    }
    fn validate_frame_count(&self, frames: u64) -> Result<(), RenderError> {
        Self::validate_frame_count(self, frames)
    }
    fn frame_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        Self::frame_block(self, output)
    }
}
