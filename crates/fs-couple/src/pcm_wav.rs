//! Deterministic little-endian PCM16 WAV encoding of a pressure history.
//!
//! Samples are physical pascals mapped through a declared full-scale.
//! The writer never peak-normalizes: that would hide a material or
//! temperature change.

/// Typed WAV refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WavError {
    /// Empty history, zero rate, or non-positive full-scale.
    InvalidInput {
        /// Which check failed.
        what: &'static str,
    },
}

impl core::fmt::Display for WavError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput { what } => write!(f, "FS-COUPLE-WAV: {what}"),
        }
    }
}

impl std::error::Error for WavError {}

/// Bounded physical output channels in one WAV. Channels above stereo use
/// WAVEFORMATEXTENSIBLE with a zero speaker mask (ordered direct outputs).
pub const MAX_WAV_CHANNELS: u16 = 64;

/// Encode mono PCM16 WAV bytes.
///
/// `full_scale_pa` maps to `i16::MAX`. Values outside ±full-scale saturate
/// after a finite check — saturation is reported by counting clips, not
/// by rewriting physics.
///
/// # Errors
/// [`WavError::InvalidInput`] for empty samples, a zero rate, a
/// non-positive full-scale, or a non-finite sample.
pub fn encode_pcm16_wav(
    pressure_pa: &[f64],
    sample_rate_hz: u32,
    full_scale_pa: f64,
) -> Result<(Vec<u8>, usize), WavError> {
    encode_pcm16_wav_interleaved(pressure_pa, sample_rate_hz, 1, full_scale_pa)
}

/// Encode frame-interleaved pressure without changing relative channel gain.
///
/// Stereo input is `[left_0, right_0, left_1, right_1, ...]`; the sample rate
/// is the frame rate of EACH channel, not the number of scalar values per second.
/// All channels use the same explicit pressure full-scale and the original
/// PCM quantizer. The returned clip count counts scalar channel samples.
/// There is no downmix, panning, normalization, or implicit channel padding.
///
/// # Errors
/// Refuses channel counts outside 1..=64, empty/incomplete frames, invalid
/// pressure/rate/full-scale, and RIFF or byte-rate overflow before allocation.
/// Mono/stereo retain their original 44-byte headers and exact payloads.
/// Larger arrays use WAVEFORMATEXTENSIBLE direct outputs, in the supplied
/// channel order, without assigning physical microphones to speaker positions.
/// Retain the application's channel identities alongside the file.
pub fn encode_pcm16_wav_interleaved(
    pressure_pa: &[f64],
    sample_rate_hz: u32,
    channels: u16,
    full_scale_pa: f64,
) -> Result<(Vec<u8>, usize), WavError> {
    if !(1..=MAX_WAV_CHANNELS).contains(&channels)
        || pressure_pa.len() % usize::from(channels) != 0 {
        return Err(WavError::InvalidInput {
            what: "PCM16 requires complete frames of 1..=64 channels",
        });
    }
    if pressure_pa.is_empty() {
        return Err(WavError::InvalidInput {
            what: "pressure history is empty",
        });
    }
    if sample_rate_hz == 0 {
        return Err(WavError::InvalidInput {
            what: "sample rate must be positive",
        });
    }
    if !(full_scale_pa > 0.0 && full_scale_pa.is_finite()) {
        return Err(WavError::InvalidInput {
            what: "full-scale pressure must be positive and finite",
        });
    }
    let n = pressure_pa.len();
    let data_bytes = n.checked_mul(2).ok_or(WavError::InvalidInput {
        what: "WAV data length overflow",
    })?;
    let data_length = u32::try_from(data_bytes).map_err(|_| WavError::InvalidInput {
        what: "WAV payload exceeds u32",
    })?;
    let header_len = pcm16_header_len(channels);
    let riff_length = data_length.checked_add((header_len - 8) as u32).ok_or(WavError::InvalidInput {
        what: "WAV payload exceeds u32",
    })?;
    let block_align = channels * 2;
    let byte_rate = sample_rate_hz.checked_mul(u32::from(block_align))
        .ok_or(WavError::InvalidInput { what: "WAV byte rate overflows u32" })?;
    let capacity = data_bytes.checked_add(header_len).ok_or(WavError::InvalidInput {
        what: "WAV data length overflow",
    })?;
    if pressure_pa.iter().any(|p| !p.is_finite()) {
        return Err(WavError::InvalidInput { what: "pressure sample is not finite" });
    }
    let mut out = Vec::with_capacity(capacity);
    out.extend_from_slice(b"RIFF");
    write_u32_le(&mut out, riff_length);
    out.extend_from_slice(b"WAVE");
    out.extend_from_slice(b"fmt ");
    write_u32_le(&mut out, if channels <= 2 { 16 } else { 40 });
    write_u16_le(&mut out, if channels <= 2 { 1 } else { 0xfffe });
    write_u16_le(&mut out, channels);
    write_u32_le(&mut out, sample_rate_hz);
    write_u32_le(&mut out, byte_rate);
    write_u16_le(&mut out, block_align);
    write_u16_le(&mut out, 16);
    if channels > 2 {
        write_extensible_fields(&mut out, channels, 16, 1);
    }
    out.extend_from_slice(b"data");
    write_u32_le(&mut out, data_length);
    let mut clips = 0usize;
    for &p in pressure_pa {
        let scaled = p / full_scale_pa * f64::from(i16::MAX);
        let quantized = if scaled >= f64::from(i16::MAX) {
            clips += 1;
            i16::MAX
        } else if scaled <= f64::from(i16::MIN) {
            clips += 1;
            i16::MIN
        } else {
            #[allow(clippy::cast_possible_truncation)]
            {
                scaled.round() as i16
            }
        };
        write_i16_le(&mut out, quantized);
    }
    Ok((out, clips))
}

fn write_u16_le(out: &mut Vec<u8>, value: u16) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn write_u32_le(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn write_i16_le(out: &mut Vec<u8>, value: i16) {
    out.extend_from_slice(&value.to_le_bytes());
}

const fn pcm16_header_len(channels: u16) -> usize {
    if channels <= 2 { 44 } else { 68 }
}

// Shared by PCM16, PCM24 and float32. A zero mask is the documented direct-out
// layout: retain channel order without inventing a speaker configuration.
// https://learn.microsoft.com/en-us/windows-hardware/drivers/audio/dsspeaker-directout-speaker-configuration
fn write_extensible_fields(out: &mut Vec<u8>, channels: u16, bits: u16, subformat: u32) {
    write_u16_le(out, 22);
    write_u16_le(out, bits);
    write_u32_le(out, match channels { 1 => 4, 2 => 3, _ => 0 });
    write_u32_le(out, subformat);
    out.extend_from_slice(&[0, 0, 0x10, 0, 0x80, 0, 0, 0xaa, 0, 0x38, 0x9b, 0x71]);
}

/// Incremental file export using the same PCM encoder and scheduled physics.
pub mod stream;

/// Shared causal integer-rate pressure filter (also used by the piano example).
pub mod decimate;
/// Output-clock adapters over the existing scheduled mechanics.
pub mod observation;

#[cfg(test)]
mod stereo_tests;

/// Shared retarded-time baffled-surface receiver, before PCM conversion.
#[path = "baffled_pressure.rs"]
pub mod baffled;

/// Explicit PCM24/float32 export of physical pressure, preserving quiet tails.
pub mod high_resolution;
