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
use super::{WavError, encode_pcm16_wav};

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
    /// Complete encoded samples across all writes.
    pub samples: u64,
    /// Clip count using exactly the existing encoder's threshold rules.
    pub clipped_samples: u64,
    /// Authored pascals mapped to positive full-scale PCM.
    pub full_scale_pa: f64,
}

/// Incremental mono PCM16 output, owning an append-positioned seekable sink.
/// No mutable sink access is exposed while header/payload counts are live.
pub struct Pcm16WavStream<W> {
    output: W,
    header: [u8; 44],
    start: u64,
    sample_rate_hz: u32,
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
    pub fn new(mut output: W, sample_rate_hz: u32, full_scale_pa: f64, max_block: usize)
        -> Result<Self, WavStreamError>
    {
        if max_block == 0 || max_block.checked_mul(2).and_then(|n| n.checked_add(44)).is_none() {
            return Err(invalid("max_block must be positive with representable encoded size"));
        }
        // Reuse the sole encoder's configuration admission and canonical header.
        let (encoded, _) = encode_pcm16_wav(&[0.0], sample_rate_hz, full_scale_pa)
            .map_err(WavStreamError::Encode)?;
        let mut header = [0_u8; 44];
        header.copy_from_slice(&encoded[..44]);
        header[4..8].fill(0);
        header[40..44].fill(0);
        let start = output.stream_position().map_err(WavStreamError::Io)?;
        let end = output.seek(SeekFrom::End(0)).map_err(WavStreamError::Io)?;
        output.seek(SeekFrom::Start(start)).map_err(WavStreamError::Io)?;
        if start != end || start.checked_add(44).is_none() {
            return Err(invalid("WAV output must begin at a representable append position"));
        }
        output.write_all(&header).map_err(WavStreamError::Io)?;
        Ok(Self { output, header, start, sample_rate_hz, full_scale_pa, max_block,
            samples: 0, clips: 0, poisoned: false })
    }

    /// Complete samples currently appended (not yet a finalized WAV).
    #[must_use]
    pub const fn samples_written(&self) -> u64 { self.samples }

    /// Inspect the sink without moving its cursor or writing through this API.
    #[must_use]
    pub const fn sink(&self) -> &W { &self.output }

    /// Append a complete finite pressure block, preserving the existing encoder's
    /// bits. Empty blocks are no-ops. Invalid input or length is refused before
    /// touching the sink or counters; I/O errors can be partial and poison it.
    pub fn write_block(&mut self, pressure_pa: &[f64]) -> Result<(), WavStreamError> {
        let count = u64::try_from(pressure_pa.len()).map_err(|_| invalid("block length exceeds u64"))?;
        let total = self.validate_append(count)?;
        if pressure_pa.len() > self.max_block { return Err(invalid("block exceeds max_block")); }
        if pressure_pa.is_empty() { return Ok(()); }
        let (encoded, clips) = encode_pcm16_wav(pressure_pa, self.sample_rate_hz, self.full_scale_pa)
            .map_err(WavStreamError::Encode)?;
        if let Err(e) = self.output.write_all(&encoded[44..]) {
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
        self.header[4..8].copy_from_slice(&((36 + data_bytes) as u32).to_le_bytes());
        self.header[40..44].copy_from_slice(&(data_bytes as u32).to_le_bytes());
        let end = self.start + 44 + data_bytes;
        self.output.seek(SeekFrom::Start(self.start)).map_err(WavStreamError::Io)?;
        self.output.write_all(&self.header).map_err(WavStreamError::Io)?;
        self.output.seek(SeekFrom::Start(end)).map_err(WavStreamError::Io)?;
        self.output.flush().map_err(WavStreamError::Io)?;
        Ok((self.output, WavStreamSummary { sample_rate_hz: self.sample_rate_hz,
            samples: self.samples, clipped_samples: self.clips, full_scale_pa: self.full_scale_pa }))
    }

    fn validate_append(&self, samples: u64) -> Result<u64, WavStreamError> {
        if self.poisoned { return Err(WavStreamError::Poisoned); }
        let total = self.samples.checked_add(samples).ok_or_else(|| invalid("sample count overflow"))?;
        if total > MAX_PCM16_WAV_SAMPLES { return Err(invalid("performance exceeds PCM16 RIFF size; use separate files")); }
        self.start.checked_add(44).and_then(|s| s.checked_add(2*total))
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
    renderer.validate_sample_rate(stream.sample_rate_hz).map_err(WavStreamError::Render)?;
    stream.validate_append(samples)?;
    if renderer.samples_rendered() != stream.samples_written() {
        return Err(invalid("renderer and WAV stream must have the same completed sample clock"));
    }
    if scratch.is_empty() || scratch.len() > renderer.context().max_block_len()
        || scratch.len() > stream.max_block {
        return Err(invalid("pressure scratch must fit both callback capacities and be nonempty"));
    }
    let block = u64::try_from(scratch.len()).map_err(|_| invalid("scratch length exceeds u64"))?;
    if renderer.samples_rendered().checked_add(samples).is_none()
        || renderer.context().blocks_rendered().checked_add(samples).is_none() {
        return Err(invalid("requested render could overflow a renderer clock"));
    }
    let mut written = 0;
    while written < samples {
        if gate.is_requested() { return Ok(ScheduledWavProgress::Cancelled { samples: written }); }
        let len = (samples-written).min(block) as usize;
        renderer.block(&mut scratch[..len]).map_err(WavStreamError::Render)?;
        stream.write_block(&scratch[..len])?;
        written += len as u64;
    }
    Ok(ScheduledWavProgress::Completed { samples: written })
}

fn invalid(what: &'static str) -> WavStreamError { WavStreamError::Invalid { what } }
