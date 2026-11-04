//! Independent physical microphones, observed simultaneously without panning.
//!
//! Each channel retains its original accepted receiver trajectory, source paths,
//! Doppler/clearance admission and Pa calibration. This owner composes their
//! existing transactional pressure evaluators; it does not interpolate one
//! microphone into another, add a velocity-dependent gain, or re-run mechanics.

use super::{CompactMotion, CubicInterval, MovingCompactRenderer, PreparedTrack, check_cancel, sizing};
use crate::pcm_wav::{MAX_WAV_CHANNELS, observation::PressureFrameRenderer};
use crate::render::RenderError;

/// Aggregate admission of already prepared channels. Limits apply to the full
/// original horizon, including an aligned prefix already observed before joining.
/// Original per-channel source and geometry budgets remain in effect as well.
#[derive(Clone, Copy, Debug)]
pub struct MovingPressureFrameBudget {
    /// Frames times channels, not frames alone.
    pub maximum_output_values: u64,
    /// Frames times the sum of all direct/reflected source paths in all channels.
    pub maximum_source_receiver_pairs: u64,
    /// Sum of the original lifetime root-iteration bounds. Analytic roots add zero.
    pub maximum_retarded_iterations: u64,
    /// This owner plus original track/coefficient/scratch payloads. Borrowed
    /// endpoints, paths and IDs, allocator metadata and rounding are excluded.
    pub maximum_payload_bytes: usize,
}

/// Frame-interleaved pressure from independently prepared physical receivers.
///
/// Takes ownership of the original mono observers and REUSES their interpolants
/// and callback scratch. Construction and callbacks allocate no new buffers.
/// The channel count is fixed; IDs preserve file order without inferring speaker
/// positions. Use shared immutable source histories to observe one coherent scene.
///
/// Channels must have exactly the same arrival-time origin, rate, horizon and
/// next sample index. Bandwidth/clearance/root tolerances may differ because
/// they were admitted for each actual microphone separately. A partially read
/// aligned array may be joined; no samples are replayed or silently discarded.
///
/// Every channel stages its complete block before ANY output or counter is
/// published. A failure/cancellation in the last channel leaves all accepted
/// clocks and the caller's output unchanged. Retrying overwrites all scratch.
/// No mutable child access is exposed while clocks are coupled.
pub struct MovingPressureFrameRenderer<'id, 'history, M: CompactMotion, const CHANNELS: usize> {
    channel_ids: [&'id str; CHANNELS],
    channels: [MovingCompactRenderer<'history, M>; CHANNELS],
    max_block: usize,
    payload_bytes: usize,
    source_receiver_pairs: u64,
    retarded_iterations: u64,
}

impl<'id, 'history, M: CompactMotion, const CHANNELS: usize>
    MovingPressureFrameRenderer<'id, 'history, M, CHANNELS>
{
    /// Bind existing physical observers without cloning or recompiling any path.
    /// Failed admission consumes the supplied owners but advances no mechanics
    /// or observer clock. Retained IDs are borrowed, bounded, unique UTF-8 text.
    pub fn new(
        channel_ids: [&'id str; CHANNELS],
        channels: [MovingCompactRenderer<'history, M>; CHANNELS],
        budget: MovingPressureFrameBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        if CHANNELS == 0 || CHANNELS > usize::from(MAX_WAV_CHANNELS) {
            return Err(sizing("moving pressure frame channel count is outside WAV limits"));
        }
        let first = &channels[0];
        let clock = first.config;
        let completed = first.completed;
        let values = clock.samples.checked_mul(CHANNELS as u64)
            .ok_or_else(|| sizing("moving pressure scalar output count overflows"))?;
        let mut max_block = usize::MAX;
        let mut payload_bytes = size_of::<Self>();
        let mut path_count = 0_u64;
        let mut retarded_iterations = 0_u64;
        for (index, channel) in channels.iter().enumerate() {
            check_cancel(cancelled)?;
            let id = channel_ids[index];
            if id.trim().is_empty() || id.len() > 128 || id.chars().any(char::is_control)
                || channel_ids[..index].contains(&id)
            {
                return Err(sizing("moving pressure channel IDs must be unique bounded nonempty text"));
            }
            if channel.config.first_receiver_time_s.to_bits() != clock.first_receiver_time_s.to_bits()
                || channel.config.sample_rate_hz != clock.sample_rate_hz
                || channel.config.samples != clock.samples || channel.completed != completed
            {
                return Err(sizing("moving pressure channels need identical physical clocks and cursors"));
            }
            max_block = max_block.min(channel.config.max_block);
            // Logical allocation sizes agree with the original mono admission.
            // Self already includes every child owner and its Vec headers.
            for (count, width) in [
                (channel.tracks.len(), size_of::<PreparedTrack<'history, M>>()),
                (channel.intervals.len(), size_of::<CubicInterval>()),
                (channel.scratch.len(), size_of::<f64>()),
            ] {
                payload_bytes = count.checked_mul(width)
                    .and_then(|bytes| payload_bytes.checked_add(bytes))
                    .ok_or_else(|| sizing("moving pressure frame payload overflows"))?;
            }
            path_count = u64::try_from(channel.source_path_count).ok()
                .and_then(|count| path_count.checked_add(count))
                .ok_or_else(|| sizing("moving pressure frame path count overflows"))?;
            retarded_iterations = retarded_iterations.checked_add(channel.maximum_retarded_iterations)
                .ok_or_else(|| sizing("moving pressure frame root work overflows"))?;
        }
        let source_receiver_pairs = path_count.checked_mul(clock.samples)
            .ok_or_else(|| sizing("moving pressure frame lifetime work overflows"))?;
        if values > budget.maximum_output_values
            || source_receiver_pairs > budget.maximum_source_receiver_pairs
            || retarded_iterations > budget.maximum_retarded_iterations
            || payload_bytes > budget.maximum_payload_bytes
        {
            return Err(sizing("moving pressure channels exceed aggregate work or payload budget"));
        }
        check_cancel(cancelled)?;
        Ok(Self { channel_ids, channels, max_block, payload_bytes,
            source_receiver_pairs, retarded_iterations })
    }

    /// Named physical receiver order in every interleaved output frame.
    pub const fn channel_ids(&self) -> &[&'id str; CHANNELS] { &self.channel_ids }
    /// Read original paths, observer identities, budgets and synchronized cursors.
    pub const fn channels(&self) -> &[MovingCompactRenderer<'history, M>; CHANNELS] { &self.channels }
    pub const fn channel_count(&self) -> usize { CHANNELS }
    pub fn frames_rendered(&self) -> u64 { self.channels[0].completed }
    pub fn remaining_frames(&self) -> u64 { self.channels[0].remaining_samples() }
    pub const fn max_frame_block_len(&self) -> usize { self.max_block }
    pub const fn retained_payload_bytes(&self) -> usize { self.payload_bytes }
    pub const fn maximum_source_receiver_pairs(&self) -> u64 { self.source_receiver_pairs }
    pub const fn maximum_retarded_iterations(&self) -> u64 { self.retarded_iterations }
    /// Original integer-indexed clock, not an accumulated or resampled timestamp.
    pub fn receiver_time_at(&self, frame: u64) -> Result<f64, RenderError> {
        self.channels[0].receiver_time_at(frame)
    }

    /// Release the original observers at their common accepted cursor. No history
    /// or filter reset is performed; each can continue on its own native clock.
    pub fn into_channels(self) -> [MovingCompactRenderer<'history, M>; CHANNELS] { self.channels }

    pub fn frame_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.frame_block_controlled(output, &mut || false)
    }

    /// Stage all channels, then publish complete frames and advance every clock.
    /// Final publication contains only bounded copies/counter additions; no
    /// cancellation or fallible work can split a coherent frame transaction.
    pub fn frame_block_controlled(
        &mut self,
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        check_cancel(cancelled)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() % CHANNELS != 0 || output.len() / CHANNELS > self.max_block {
            return Err(sizing("moving pressure output needs complete frames within callback capacity"));
        }
        let count = output.len() / CHANNELS;
        let frames = u64::try_from(count)
            .map_err(|_| sizing("moving pressure frame count exceeds u64"))?;
        self.validate_frame_count(frames)?;
        for channel in &mut self.channels {
            channel.preview_block(count, cancelled)?;
        }
        check_cancel(cancelled)?;
        for (index, frame) in output.chunks_exact_mut(CHANNELS).enumerate() {
            for (value, channel) in frame.iter_mut().zip(&self.channels) {
                *value = channel.scratch[index];
            }
        }
        for channel in &mut self.channels { channel.completed += frames; }
        Ok(())
    }
}

impl<M: CompactMotion, const CHANNELS: usize> PressureFrameRenderer
    for MovingPressureFrameRenderer<'_, '_, M, CHANNELS>
{
    fn channel_count(&self) -> usize { CHANNELS }
    fn frames_rendered(&self) -> u64 { Self::frames_rendered(self) }
    fn max_frame_block_len(&self) -> usize { self.max_block }
    fn validate_frame_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate != self.channels[0].config.sample_rate_hz {
            return Err(sizing("sink frame rate must match the original moving receiver clocks"));
        }
        Ok(())
    }
    fn validate_frame_count(&self, frames: u64) -> Result<(), RenderError> {
        if frames > self.remaining_frames() {
            return Err(sizing("moving pressure frames exceed retained physical history coverage"));
        }
        Ok(())
    }
    fn frame_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        Self::frame_block(self, output)
    }
}
