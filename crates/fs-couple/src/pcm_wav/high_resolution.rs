//! Explicit high-resolution pressure export; no normalization or invented tails.
//! WAVEFORMATEXTENSIBLE preserves sample precision and mono/stereo channel order.
//! https://learn.microsoft.com/en-us/windows/win32/api/mmreg/ns-mmreg-waveformatextensible
use super::WavError;

/// Output storage, independent of the mechanics and physical pressure scale.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Encoding {
    /// Signed little-endian packed PCM, 24 valid bits in 24 container bits.
    Pcm24,
    /// IEEE binary32, including finite values outside nominal digital full scale.
    Float32,
}

fn invalid(what: &'static str) -> WavError { WavError::InvalidInput { what } }
fn u16_le(out: &mut Vec<u8>, x: u16) { out.extend_from_slice(&x.to_le_bytes()); }
fn u32_le(out: &mut Vec<u8>, x: u32) { out.extend_from_slice(&x.to_le_bytes()); }

/// Encode complete mono or interleaved left/right pressure frames.
///
/// The SAME explicit `full_scale_pa` maps to +1 in either channel. PCM24 maps
/// this to 8388607; -full scale maps to -8388607, keeping symmetric pressure
/// scaling. Out-of-range PCM refuses instead of saturating. Float32 retains
/// representable excursions beyond +/-1, but a player may clip them: attenuate
/// explicitly before fixed-point playback. Neither format changes the source
/// clock, peak-normalizes, adds dither, or reconstructs previously quantized data.
/// Floating export rounds f64 to f32 and is not a bit-exact f64 archive.
///
/// # Errors
/// Invalid/empty/incomplete frames, nonfinite pressure, invalid clock/scale,
/// RIFF/byte-rate/allocation-size overflow, PCM clipping, or f32 range overflow.
/// All samples and lengths are admitted before output allocation.
pub fn encode_pressure_wav(pressure_pa: &[f64], sample_rate_hz: u32, channels: u16,
    full_scale_pa: f64, encoding: Encoding) -> Result<Vec<u8>, WavError> {
    if !matches!(channels, 1 | 2) || pressure_pa.is_empty()
        || pressure_pa.len() % usize::from(channels) != 0 {
        return Err(invalid("high-resolution WAV needs complete nonempty mono/stereo frames"));
    }
    if sample_rate_hz == 0 || !full_scale_pa.is_finite() || full_scale_pa <= 0.0 {
        return Err(invalid("WAV needs a positive clock and finite positive pressure full scale"));
    }
    let floating = encoding == Encoding::Float32;
    let bytes: u16 = if floating { 4 } else { 3 };
    let header: usize = if floating { 80 } else { 68 };
    let data = pressure_pa.len().checked_mul(usize::from(bytes))
        .ok_or_else(|| invalid("WAV data size overflow"))?;
    let padded = data.checked_add(data % 2).ok_or_else(|| invalid("WAV pad overflow"))?;
    let length = u32::try_from(data).map_err(|_| invalid("WAV data exceeds RIFF capacity"))?;
    let total = header.checked_add(padded).ok_or_else(|| invalid("WAV size overflow"))?;
    let riff = u32::try_from(total - 8).map_err(|_| invalid("WAV file exceeds RIFF capacity"))?;
    let align = channels * bytes;
    let byte_rate = sample_rate_hz.checked_mul(u32::from(align))
        .ok_or_else(|| invalid("WAV byte rate overflow"))?;
    let frames = u32::try_from(pressure_pa.len() / usize::from(channels))
        .map_err(|_| invalid("WAV frame count overflow"))?;
    for &p in pressure_pa {
        let normalized = p / full_scale_pa;
        if !p.is_finite() || !normalized.is_finite() {
            return Err(invalid("pressure or normalized sample is nonfinite"));
        }
        if floating {
            if normalized.abs() > f64::from(f32::MAX) {
                return Err(invalid("pressure exceeds finite float32 WAV range"));
            }
        } else if normalized.abs() > 1.0 {
            return Err(invalid("PCM24 would clip; increase the explicit pressure full scale"));
        }
    }
    let mut out = Vec::with_capacity(total);
    out.extend_from_slice(b"RIFF"); u32_le(&mut out, riff); out.extend_from_slice(b"WAVEfmt ");
    u32_le(&mut out, 40); u16_le(&mut out, 0xfffe); u16_le(&mut out, channels);
    u32_le(&mut out, sample_rate_hz); u32_le(&mut out, byte_rate); u16_le(&mut out, align);
    u16_le(&mut out, bytes * 8); u16_le(&mut out, 22); u16_le(&mut out, bytes * 8);
    u32_le(&mut out, if channels == 1 { 4 } else { 3 }); // center or left/right
    // PCM/IEEE_FLOAT subtype GUID, serialized by its little-endian fields.
    u32_le(&mut out, if floating { 3 } else { 1 });
    out.extend_from_slice(&[0, 0, 0x10, 0, 0x80, 0, 0, 0xaa, 0, 0x38, 0x9b, 0x71]);
    if floating {
        out.extend_from_slice(b"fact"); u32_le(&mut out, 4); u32_le(&mut out, frames);
    }
    out.extend_from_slice(b"data"); u32_le(&mut out, length);
    for &p in pressure_pa {
        let normalized = p / full_scale_pa;
        if floating {
            out.extend_from_slice(&(normalized as f32).to_le_bytes());
        } else {
            let sample = (normalized * 8_388_607.0).round() as i32;
            out.extend_from_slice(&sample.to_le_bytes()[..3]);
        }
    }
    if data % 2 != 0 { out.push(0); }
    debug_assert_eq!(out.len(), total);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn read_u32(b: &[u8], at: usize) -> u32 { u32::from_le_bytes(b[at..at+4].try_into().unwrap()) }
    fn pcm(b: &[u8], i: usize) -> i32 {
        let p = &b[68 + i*3..68 + i*3 + 3];
        i32::from_le_bytes([p[0], p[1], p[2], if p[2] & 128 != 0 {255} else {0}])
    }
    #[test]
    fn packed_24_bit_extensible_layout_preserves_signed_samples_and_odd_padding() {
        let p = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let b = encode_pressure_wav(&p, 48_000, 1, 2.0, Encoding::Pcm24).unwrap();
        assert_eq!(&b[..4], b"RIFF"); assert_eq!(read_u32(&b, 4) as usize + 8, b.len());
        assert_eq!(read_u32(&b, 16), 40); assert_eq!(&b[20..24], &[254,255,1,0]);
        assert_eq!(read_u32(&b, 28), 144_000); assert_eq!(&b[32..40], &[3,0,24,0,22,0,24,0]);
        assert_eq!(read_u32(&b, 40), 4); assert_eq!(read_u32(&b, 44), 1);
        assert_eq!(&b[60..64], b"data"); assert_eq!(read_u32(&b, 64), 15);
        assert_eq!(b.len(), 84); assert_eq!(*b.last().unwrap(), 0);
        assert_eq!((0..5).map(|i| pcm(&b, i)).collect::<Vec<_>>(),
            [-8_388_607, -4_194_304, 0, 4_194_304, 8_388_607]);
    }
    #[test]
    fn quiet_pressure_survives_where_the_existing_pcm16_quantizer_returns_zero() {
        let p = [1e-6, -1e-6];
        let (old, clips) = super::super::encode_pcm16_wav(&p, 48_000, 1.0).unwrap();
        assert_eq!(clips, 0); assert_eq!(&old[44..], &[0; 4]);
        let b = encode_pressure_wav(&p, 48_000, 1, 1.0, Encoding::Pcm24).unwrap();
        assert_eq!((pcm(&b, 0), pcm(&b, 1)), (8, -8));
        for (i, p) in p.iter().enumerate() {
            assert!((f64::from(pcm(&b, i)) / 8_388_607.0 - p).abs() <= 0.5/8_388_607.0);
        }
    }
    #[test]
    fn float_headroom_and_quiet_tails_are_not_clipped_or_peak_normalized() {
        let p = [25.0, -25.0, 2e-25, -2e-25];
        let b = encode_pressure_wav(&p, 96_000, 2, 2.0, Encoding::Float32).unwrap();
        assert_eq!(read_u32(&b, 40), 3); assert_eq!(read_u32(&b, 44), 3);
        assert_eq!(&b[60..64], b"fact"); assert_eq!(read_u32(&b, 68), 2);
        assert_eq!(&b[72..76], b"data"); assert_eq!(read_u32(&b, 76), 16);
        for (i, p) in p.iter().enumerate() {
            let got = f32::from_le_bytes(b[80+i*4..84+i*4].try_into().unwrap());
            assert_eq!(got, (p / 2.0) as f32);
        }
        assert_eq!(read_u32(&b, 4) as usize + 8, b.len());
    }
    #[test]
    fn stereo_is_interleaved_with_one_pressure_scale_not_independent_channel_gains() {
        let b = encode_pressure_wav(&[1.0, -0.25, 0.5, -1.0], 48_000, 2, 2.0, Encoding::Pcm24).unwrap();
        assert_eq!(read_u32(&b, 28), 288_000); assert_eq!(&b[32..34], &[6, 0]);
        assert_eq!(read_u32(&b, 40), 3);
        assert_eq!((0..4).map(|i| pcm(&b,i)).collect::<Vec<_>>(),
            [4_194_304, -1_048_576, 2_097_152, -4_194_304]);
    }
    #[test]
    fn clipping_and_nonrepresentable_samples_refuse_instead_of_corrupting_audio() {
        assert!(encode_pressure_wav(&[1.01], 48_000, 1, 1.0, Encoding::Pcm24).is_err());
        assert!(encode_pressure_wav(&[-1.01], 48_000, 1, 1.0, Encoding::Pcm24).is_err());
        for format in [Encoding::Pcm24, Encoding::Float32] {
            for p in [f64::NAN, f64::INFINITY, f64::MAX] {
                assert!(encode_pressure_wav(&[p], 48_000, 1, 1.0, format).is_err());
            }
            assert!(encode_pressure_wav(&[1.0], 48_000, 1, f64::MIN_POSITIVE, format).is_err());
        }
    }
    #[test]
    fn invalid_clocks_channels_frames_and_scales_refuse() {
        for format in [Encoding::Pcm24, Encoding::Float32] {
            for channels in [0, 2, 3] {
                assert!(encode_pressure_wav(&[0.0], 48_000, channels, 1.0, format).is_err());
            }
            assert!(encode_pressure_wav(&[], 48_000, 1, 1.0, format).is_err());
            for rate in [0, u32::MAX] {
                assert!(encode_pressure_wav(&[0.0], rate, 1, 1.0, format).is_err());
            }
            for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                assert!(encode_pressure_wav(&[0.0], 48_000, 1, scale, format).is_err());
            }
        }
    }
}
