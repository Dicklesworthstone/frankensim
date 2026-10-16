//! Encode the computed pressure directly, before any intermediate quantization.
//! One explicit pressure scale applies to both channels and every format.
use fs_couple::pcm_wav::{encode_pcm16_wav_interleaved, high_resolution::{self, Encoding}};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Format {
    #[default]
    Pcm16,
    Pcm24,
    Float32,
}

impl Format {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "pcm16" => Ok(Self::Pcm16),
            "pcm24" => Ok(Self::Pcm24),
            "float32" => Ok(Self::Float32),
            _ => Err("--wav-format requires pcm16, pcm24 or float32".into()),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Pcm16 => "PCM16",
            Self::Pcm24 => "PCM24",
            Self::Float32 => "Float32",
        }
    }

    /// This changes storage precision only. Fixed-point formats refuse clipping;
    /// float retains finite excursions beyond +/-1 for subsequent attenuation.
    /// Neither format normalizes peaks or changes the physical pressure history.
    pub fn encode(self, pressure: &[f64], rate: u32, channels: u16, full_scale_pa: f64)
        -> Result<Vec<u8>, String> {
        match self {
            Self::Pcm16 => {
                let (wav, clips) = encode_pcm16_wav_interleaved(pressure, rate, channels, full_scale_pa)
                    .map_err(|e| e.to_string())?;
                if clips != 0 {
                    return Err(format!("PCM16 refuses {clips} clipped samples at {full_scale_pa} Pa full scale; increase the explicit pressure full scale; no WAV published"));
                }
                Ok(wav)
            }
            Self::Pcm24 | Self::Float32 => high_resolution::encode_pressure_wav(
                pressure, rate, channels, full_scale_pa,
                if self == Self::Pcm24 { Encoding::Pcm24 } else { Encoding::Float32 },
            ).map_err(|e| e.to_string()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_formats_select_the_existing_encoders_without_changing_legacy_pcm() {
        let pressure = [0.0, -0.5, 0.25, 1.0];
        assert_eq!(Format::default(), Format::Pcm16);
        assert_eq!(Format::parse("pcm16").unwrap().encode(&pressure, 48_000, 2, 2.0).unwrap(),
            encode_pcm16_wav_interleaved(&pressure, 48_000, 2, 2.0).unwrap().0);
        for (name, encoding) in [("pcm24", Encoding::Pcm24), ("float32", Encoding::Float32)] {
            assert_eq!(Format::parse(name).unwrap().encode(&pressure, 48_000, 2, 2.0).unwrap(),
                high_resolution::encode_pressure_wav(&pressure, 48_000, 2, 2.0, encoding).unwrap());
        }
        for invalid in ["", "24", "float64", "PCM24"] {
            assert!(Format::parse(invalid).is_err());
        }
    }

    #[test]
    fn only_float_export_accepts_representable_pressure_above_digital_full_scale() {
        let pressure = [3.0, -3.0];
        for fixed in [Format::Pcm16, Format::Pcm24] {
            assert!(fixed.encode(&pressure, 48_000, 2, 2.0).is_err());
        }
        let wav = Format::Float32.encode(&pressure, 48_000, 2, 2.0).unwrap();
        assert_eq!(f32::from_le_bytes(wav[80..84].try_into().unwrap()), 1.5);
        assert_eq!(f32::from_le_bytes(wav[84..88].try_into().unwrap()), -1.5);
        for format in [Format::Pcm16, Format::Pcm24, Format::Float32] {
            assert!(format.encode(&[f64::NAN], 48_000, 1, 2.0).is_err());
            assert!(format.encode(&pressure, 48_000, 2, 0.0).is_err());
        }
    }
}
