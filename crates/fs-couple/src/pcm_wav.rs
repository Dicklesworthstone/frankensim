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
    let mut out = Vec::with_capacity(44 + data_bytes);
    out.extend_from_slice(b"RIFF");
    write_u32_le(
        &mut out,
        u32::try_from(36 + data_bytes).map_err(|_| WavError::InvalidInput {
            what: "WAV payload exceeds u32",
        })?,
    );
    out.extend_from_slice(b"WAVE");
    out.extend_from_slice(b"fmt ");
    write_u32_le(&mut out, 16);
    write_u16_le(&mut out, 1);
    write_u16_le(&mut out, 1);
    write_u32_le(&mut out, sample_rate_hz);
    write_u32_le(&mut out, sample_rate_hz * 2);
    write_u16_le(&mut out, 2);
    write_u16_le(&mut out, 16);
    out.extend_from_slice(b"data");
    write_u32_le(
        &mut out,
        u32::try_from(data_bytes).map_err(|_| WavError::InvalidInput {
            what: "WAV payload exceeds u32",
        })?,
    );
    let mut clips = 0usize;
    for &p in pressure_pa {
        if !p.is_finite() {
            return Err(WavError::InvalidInput {
                what: "pressure sample is not finite",
            });
        }
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
