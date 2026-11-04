//! Bounded-memory PCM16 file export for scheduled physical performances.
//!
//! Blocks use the existing encoder's exact pressure scaling, quantization and
//! clip counting. Only their payloads are appended; one RIFF header is finalized
//! at the end. Memory is O(max_block), not O(performance duration). This is an
//! OFFLINE I/O path: block encoding allocates and file writes can block. It is
//! not a real-time device callback or a filesystem-durability guarantee.

use std::io::{Seek, SeekFrom, Write};
use fs_exec::CancelGate;
use crate::render::RenderError;
use crate::render::schedule::ScheduledRenderer;
use super::{MAX_WAV_CHANNELS, WavError, encode_pcm16_wav_interleaved, pcm16_header_len};
use super::observation::{PressureFrameRenderer, PressureRenderer};

/// Largest mono PCM16 history fitting this RIFF/WAVE format (not RF64).
pub const MAX_PCM16_WAV_SAMPLES: u64 = (u32::MAX as u64 - 36) / 2;

/// Input, physics or output failure. I/O failures poison the stream because a
/// write may have partially changed the sink; it must not then be finalized.
#[derive(Debug)]
pub enum WavStreamError {
    /// A size, clock or stream-position requirement was not met.
    Invalid { what: &'static str },
    /// A prior partial I/O failure prevents further writes/finalization.
    Poisoned,
    /// The sole PCM encoder refused the whole input block before output changed.
    Encode(WavError),
    /// The sink failed; its bytes may be incomplete.
    Io(std::io::Error),
    /// Physical rendering failed. The failed callback was not written.
    Render(RenderError),
}
impl core::fmt::Display for WavStreamError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid { what } => write!(f, "WAV stream: {what}"),
            Self::Poisoned => write!(f, "WAV stream has incomplete I/O; cannot continue or finalize"),
            Self::Encode(e) => write!(f, "WAV block: {e}"),
            Self::Io(e) => write!(f, "WAV I/O: {e}"),
            Self::Render(e) => write!(f, "WAV render: {e}"),
        }
    }
}
impl std::error::Error for WavStreamError {}

/// Statistics of a successfully finalized stream, with no implicit normalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WavStreamSummary {
    /// Fixed sample rate of the file [Hz].
    pub sample_rate_hz: u32,
    /// Number of interleaved physical output channels.
    pub channels: u16,
    /// Complete simultaneous frames across all writes.
    pub frames: u64,
    /// Scalar channel samples across all writes: frames times channels.
    pub samples: u64,
    /// Clip count using exactly the existing encoder's threshold rules.
    pub clipped_samples: u64,
    /// Authored pascals mapped to positive full-scale PCM.
    pub full_scale_pa: f64,
}

/// Incremental interleaved PCM16 output, owning an append-positioned seekable sink.
/// No mutable sink access is exposed while header/payload counts are live.
pub struct Pcm16WavStream<W> {
    output: W,
    header: [u8; 68],
    header_len: usize,
    start: u64,
    sample_rate_hz: u32,
    channels: u16,
    full_scale_pa: f64,
    max_block: usize,
    samples: u64,
    clips: u64,
    poisoned: bool,
}

impl<W: Write + Seek> Pcm16WavStream<W> {
    /// Begin a WAV at the sink's current END, refusing to overwrite existing
    /// bytes. A prefix is allowed (e.g. an enclosing container), but a cursor
    /// positioned inside existing data is not. `finish` must succeed before
    /// this new WAV is called complete; the provisional RIFF length is zero.
    pub fn new(output: W, sample_rate_hz: u32, full_scale_pa: f64, max_block: usize)
        -> Result<Self, WavStreamError>
    {
        Self::new_interleaved(output, sample_rate_hz, 1, full_scale_pa, max_block)
    }

    /// Begin 1..=64 channels on one frame clock and one authored Pa/PCM scale.
    /// `max_block_frames` counts complete simultaneous frames, not scalars.
    /// Mono/stereo headers retain their original bytes. Larger arrays use the
    /// original encoder's extensible direct-output layout and supplied order.
    pub fn new_interleaved(mut output: W, sample_rate_hz: u32, channels: u16,
        full_scale_pa: f64, max_block_frames: usize) -> Result<Self, WavStreamError>
    {
        if !(1..=MAX_WAV_CHANNELS).contains(&channels) {
            return Err(invalid("WAV requires 1..=64 channels"));
        }
        let header_len = pcm16_header_len(channels);
        let max_block = max_block_frames.checked_mul(usize::from(channels))
            .ok_or_else(|| invalid("WAV callback sample count overflow"))?;
        if max_block_frames == 0
            || max_block.checked_mul(2).and_then(|n| n.checked_add(header_len)).is_none() {
            return Err(invalid("max_block_frames must be positive with representable encoded size"));
        }
        let zero_frame = [0.0; MAX_WAV_CHANNELS as usize];
        // Reuse the sole encoder's configuration admission and canonical header.
        let (encoded, _) = encode_pcm16_wav_interleaved(&zero_frame[..usize::from(channels)],
            sample_rate_hz, channels, full_scale_pa).map_err(WavStreamError::Encode)?;
        let mut header = [0_u8; 68];
        header[..header_len].copy_from_slice(&encoded[..header_len]);
        header[4..8].fill(0);
        header[header_len-4..header_len].fill(0);
        let start = output.stream_position().map_err(WavStreamError::Io)?;
        let end = output.seek(SeekFrom::End(0)).map_err(WavStreamError::Io)?;
        output.seek(SeekFrom::Start(start)).map_err(WavStreamError::Io)?;
        if start != end || start.checked_add(header_len as u64).is_none() {
            return Err(invalid("WAV output must begin at a representable append position"));
        }
        output.write_all(&header[..header_len]).map_err(WavStreamError::Io)?;
        Ok(Self { output, header, header_len, start, sample_rate_hz, channels, full_scale_pa, max_block,
            samples: 0, clips: 0, poisoned: false })
    }

    /// Scalar channel samples currently appended (not yet a finalized WAV).
    #[must_use]
    pub const fn samples_written(&self) -> u64 { self.samples }

    /// Complete simultaneous frames currently appended.
    #[must_use]
    pub const fn frames_written(&self) -> u64 { self.samples / self.channels as u64 }

    /// Fixed number of output channels, in the supplied interleave order.
    #[must_use]
    pub const fn channels(&self) -> u16 { self.channels }

    /// Inspect the sink without moving its cursor or writing through this API.
    #[must_use]
    pub const fn sink(&self) -> &W { &self.output }

    /// Append a complete finite pressure block, preserving the existing encoder's
    /// bits. Empty blocks are no-ops. Invalid input or length is refused before
    /// touching the sink or counters; I/O errors can be partial and poison it.
    pub fn write_block(&mut self, pressure_pa: &[f64]) -> Result<(), WavStreamError> {
        if pressure_pa.len() % usize::from(self.channels) != 0 {
            return Err(invalid("WAV block must contain complete interleaved frames"));
        }
        let count = u64::try_from(pressure_pa.len()).map_err(|_| invalid("block length exceeds u64"))?;
        let total = self.validate_append(count)?;
        if pressure_pa.len() > self.max_block { return Err(invalid("block exceeds max_block")); }
        if pressure_pa.is_empty() { return Ok(()); }
        let (encoded, clips) = encode_pcm16_wav_interleaved(pressure_pa, self.sample_rate_hz,
            self.channels, self.full_scale_pa).map_err(WavStreamError::Encode)?;
        if let Err(e) = self.output.write_all(&encoded[self.header_len..]) {
            self.poisoned = true;
            return Err(WavStreamError::Io(e));
        }
        self.samples = total;
        self.clips += clips as u64;
        Ok(())
    }

    /// Finalize lengths for the actually written prefix, flush, and return the
    /// owned sink and statistics. Zero samples are legal after cancellation.
    /// This does not synchronize a file to durable storage or claim that an
    /// interrupted performance reached its requested duration.
    pub fn finish(mut self) -> Result<(W, WavStreamSummary), WavStreamError> {
        self.validate_append(0)?;
        let data_bytes = self.samples * 2;
        self.header[4..8].copy_from_slice(&((self.header_len as u64 - 8 + data_bytes) as u32).to_le_bytes());
        self.header[self.header_len-4..self.header_len].copy_from_slice(&(data_bytes as u32).to_le_bytes());
        let end = self.start + self.header_len as u64 + data_bytes;
        self.output.seek(SeekFrom::Start(self.start)).map_err(WavStreamError::Io)?;
        self.output.write_all(&self.header[..self.header_len]).map_err(WavStreamError::Io)?;
        self.output.seek(SeekFrom::Start(end)).map_err(WavStreamError::Io)?;
        self.output.flush().map_err(WavStreamError::Io)?;
        let frames = self.frames_written();
        Ok((self.output, WavStreamSummary { sample_rate_hz: self.sample_rate_hz,
            channels: self.channels, frames, samples: self.samples, clipped_samples: self.clips,
            full_scale_pa: self.full_scale_pa }))
    }

    fn validate_append(&self, samples: u64) -> Result<u64, WavStreamError> {
        if self.poisoned { return Err(WavStreamError::Poisoned); }
        let total = self.samples.checked_add(samples).ok_or_else(|| invalid("sample count overflow"))?;
        let max_samples = (u64::from(u32::MAX) - (self.header_len as u64 - 8)) / 2;
        if total > max_samples { return Err(invalid("performance exceeds PCM16 RIFF size; use separate files")); }
        self.start.checked_add(self.header_len as u64).and_then(|s| s.checked_add(2*total))
            .ok_or_else(|| invalid("WAV sink offset overflow"))?;
        Ok(total)
    }
}

/// Outcome of THIS render call. A cancelled prefix can continue with the same
/// renderer and stream and a fresh/unrequested cancellation gate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScheduledWavProgress {
    /// All requested samples were rendered and appended.
    Completed { samples: u64 },
    /// Cancellation was observed before the next callback.
    Cancelled { samples: u64 },
}

/// Complete simultaneous frames written by THIS call, independent of the
/// number of physical channels. A cancelled prefix can continue on the same
/// renderer and stream with an unrequested gate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WavFrameProgress {
    Completed { frames: u64 },
    Cancelled { frames: u64 },
}

/// Render directly to a WAV stream using one caller-owned pressure block.
///
/// Refuse incompatible rates, timelines, whole-request RIFF overflow and invalid
/// scratch sizes BEFORE any voice or output advances. Cancellation is observed
/// only between host callbacks: an in-flight callback drains to the file. Its
/// internal control splits do not create extra cancellation boundaries. The
/// final short callback is retained, not zero-padded or dropped.
///
/// A physics failure leaves the previous complete file prefix available for
/// explicit finalization, but the renderer cannot resume. An I/O failure also
/// poisons the stream because its last payload may be incomplete.
pub fn render_scheduled_pcm16<W: Write + Seek>(
    renderer: &mut ScheduledRenderer,
    stream: &mut Pcm16WavStream<W>,
    gate: &CancelGate,
    scratch: &mut [f64],
    samples: u64,
) -> Result<ScheduledWavProgress, WavStreamError> {
    render_pressure_pcm16(renderer, stream, gate, scratch, samples)
}

/// Stream any admitted pressure producer on its output clock. In particular,
/// a DecimatedRenderer drains all high-rate mechanics and filter stages before
/// a callback is written. Continuation retains BOTH the producer and the sink;
/// no lost filter history, synthetic tail, or partial output interval is inferred.
/// The PCM writer, prefix/error rules and cancellation loop are shared with
/// render_scheduled_pcm16, not an alternative encoder.
pub fn render_pressure_pcm16<W: Write + Seek>(
    renderer: &mut impl PressureRenderer,
    stream: &mut Pcm16WavStream<W>,
    gate: &CancelGate,
    scratch: &mut [f64],
    samples: u64,
) -> Result<ScheduledWavProgress, WavStreamError> {
    if stream.channels != 1 {
        return Err(invalid("mono pressure rendering requires a mono WAV stream"));
    }
    Ok(match render_pressure_frames_pcm16(renderer, stream, gate, scratch, samples)? {
        WavFrameProgress::Completed { frames } => ScheduledWavProgress::Completed { samples: frames },
        WavFrameProgress::Cancelled { frames } => ScheduledWavProgress::Cancelled { samples: frames },
    })
}

/// Stream coherent simultaneous physical pressure frames through the original
/// PCM16 writer. Channel count, rate, completed clock, complete-frame scratch,
/// whole-request RIFF size and producer admission are checked before advancing
/// either owner. All channels share one Pa/PCM scale. Cancellation is observed
/// between complete callbacks; each in-flight callback drains before stopping.
/// The final short callback keeps every channel without padding or downmixing.
pub fn render_pressure_frames_pcm16<W: Write + Seek>(
    renderer: &mut (impl PressureFrameRenderer + ?Sized),
    stream: &mut Pcm16WavStream<W>,
    gate: &CancelGate,
    scratch: &mut [f64],
    frames: u64,
) -> Result<WavFrameProgress, WavStreamError> {
    let channels = usize::from(stream.channels);
    if renderer.channel_count() != channels {
        return Err(invalid("pressure producer and WAV channel counts differ"));
    }
    renderer.validate_frame_rate(stream.sample_rate_hz).map_err(WavStreamError::Render)?;
    let samples = frames.checked_mul(u64::from(stream.channels))
        .ok_or_else(|| invalid("WAV request sample count overflow"))?;
    stream.validate_append(samples)?;
    if renderer.frames_rendered() != stream.frames_written() {
        return Err(invalid("renderer and WAV stream must have the same completed frame clock"));
    }
    if scratch.is_empty() || scratch.len() % channels != 0
        || scratch.len() / channels > renderer.max_frame_block_len()
        || scratch.len() > stream.max_block {
        return Err(invalid("pressure scratch must fit both capacities and contain complete nonempty frames"));
    }
    let block = u64::try_from(scratch.len() / channels).map_err(|_| invalid("scratch frame count exceeds u64"))?;
    renderer.validate_frame_count(frames).map_err(|error| match error {
        RenderError::Sizing { what } => invalid(what),
        other => WavStreamError::Render(other),
    })?;
    let mut written = 0;
    while written < frames {
        if gate.is_requested() { return Ok(WavFrameProgress::Cancelled { frames: written }); }
        let count = (frames-written).min(block) as usize;
        let len = count * channels; // bounded by the admitted scratch length
        renderer.frame_block(&mut scratch[..len]).map_err(WavStreamError::Render)?;
        stream.write_block(&scratch[..len])?;
        written += count as u64;
    }
    Ok(WavFrameProgress::Completed { frames: written })
}

fn invalid(what: &'static str) -> WavStreamError { WavStreamError::Invalid { what } }
