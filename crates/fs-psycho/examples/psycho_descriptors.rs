//! Grader bridge example for `fs-psycho` (`M7.5-W1`, `pexfp.14`).
//!
//! Reads one mono PCM16 or IEEE float32 WAV file with an explicit calibration
//! (`--pa-per-full-scale <Pa>`, `--calibration-source <sidecar|assumed>`) and
//! explicit sound field (`--field <free|diffuse>`), evaluates `fs-psycho`
//! descriptors, and prints a single canonical JSON line (`fs-psycho-descriptors/1`)
//! using Rust's round-trip-exact `f64` formatting.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::ExitCode;

use fs_psycho::fluctuation::{FS_FRAME, fluctuation_strength_frame};
use fs_psycho::roughness::{DW_BLOCK, roughness_dw_block};
use fs_psycho::signal::{loudness_stationary_from_pcm, loudness_time_varying};
use fs_psycho::tonality::{prominence_ratio_ecma, tone_to_noise_ecma};
use fs_psycho::{PsychoError, SoundField, log_attack_time, sharpness_din};

/// Maximum admitted WAV byte size (64 MiB).
pub const MAX_WAV_BYTES: usize = 64 * 1024 * 1024;

/// Calibration provenance source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationSource {
    /// Calibrated from a simulation render sidecar (`full_scale_pa`).
    Sidecar,
    /// Assumed calibration level (e.g. uncalibrated reference recordings).
    Assumed,
}

impl CalibrationSource {
    fn parse(raw: &str) -> Result<Self, BridgeError> {
        match raw {
            "sidecar" => Ok(Self::Sidecar),
            "assumed" => Ok(Self::Assumed),
            other => Err(BridgeError::InvalidCalibrationSource(other.to_owned())),
        }
    }

    const fn as_str(self) -> &'static str {
        match self {
            Self::Sidecar => "sidecar",
            Self::Assumed => "assumed",
        }
    }

    const fn level_unverified(self) -> bool {
        matches!(self, Self::Assumed)
    }
}

/// CLI and WAV-bridge refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum BridgeError {
    /// Missing or invalid `--pa-per-full-scale` (`PsychoError::UncalibratedAbsolute`).
    MissingOrInvalidFullScale(String),
    /// Missing `--calibration-source`.
    MissingCalibrationSource,
    /// Unrecognized `--calibration-source`.
    InvalidCalibrationSource(String),
    /// Missing `--field`.
    MissingField,
    /// Unrecognized `--field`.
    InvalidField(String),
    /// Missing `--wav`.
    MissingWavPath,
    /// Unknown CLI flag.
    UnknownArgument(String),
    /// File I/O failure.
    Io(String),
    /// WAV file exceeds `MAX_WAV_BYTES`.
    WavOversize {
        /// Observed byte length.
        bytes: usize,
        /// Maximum allowed byte length.
        max_bytes: usize,
    },
    /// WAV file is not mono.
    Multichannel {
        /// Observed channel count.
        channels: u16,
    },
    /// Malformed or unsupported WAV container/encoding.
    WavFormat(&'static str),
}

impl core::fmt::Display for BridgeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingOrInvalidFullScale(msg) => {
                write!(f, "{}: {msg}", PsychoError::UncalibratedAbsolute)
            }
            Self::MissingCalibrationSource => {
                write!(
                    f,
                    "missing required --calibration-source <sidecar|assumed>"
                )
            }
            Self::InvalidCalibrationSource(v) => {
                write!(
                    f,
                    "invalid --calibration-source '{v}': expected 'sidecar' or 'assumed'"
                )
            }
            Self::MissingField => write!(f, "missing required --field <free|diffuse>"),
            Self::InvalidField(v) => {
                write!(f, "invalid --field '{v}': expected 'free' or 'diffuse'")
            }
            Self::MissingWavPath => write!(f, "missing required --wav <path>"),
            Self::UnknownArgument(arg) => write!(f, "unknown argument '{arg}'"),
            Self::Io(e) => write!(f, "io error: {e}"),
            Self::WavOversize { bytes, max_bytes } => {
                write!(f, "wav size {bytes} bytes exceeds maximum {max_bytes} bytes")
            }
            Self::Multichannel { channels } => {
                write!(f, "mono wav required (got {channels} channels)")
            }
            Self::WavFormat(what) => write!(f, "invalid wav format: {what}"),
        }
    }
}

/// Parsed CLI configuration.
#[derive(Debug, Clone)]
pub struct BridgeArgs {
    /// Path to the mono WAV file.
    pub wav_path: PathBuf,
    /// Full-scale sound pressure [Pa].
    pub pa_per_full_scale: f64,
    /// Calibration provenance (`sidecar` or `assumed`).
    pub calibration_source: CalibrationSource,
    /// Sound field (`Free` or `Diffuse`).
    pub field: SoundField,
    /// Stationary loudness `time_skip` [s] (default `0.0`).
    pub time_skip_s: f64,
}

/// Parse CLI arguments without defaults for calibration or sound field.
pub fn parse_args<I, S>(iter: I) -> Result<BridgeArgs, BridgeError>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let args: Vec<String> = iter.into_iter().map(Into::into).collect();
    let mut wav_path: Option<PathBuf> = None;
    let mut pa_per_full_scale: Option<f64> = None;
    let mut calibration_source: Option<CalibrationSource> = None;
    let mut field: Option<SoundField> = None;
    let mut time_skip_s: f64 = 0.0;

    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        match arg.as_str() {
            "--wav" => {
                i += 1;
                let val = args.get(i).ok_or(BridgeError::MissingWavPath)?;
                wav_path = Some(PathBuf::from(val));
            }
            "--pa-per-full-scale" => {
                i += 1;
                let val = args.get(i).ok_or_else(|| {
                    BridgeError::MissingOrInvalidFullScale("missing value".to_owned())
                })?;
                let parsed: f64 = val.parse().map_err(|_| {
                    BridgeError::MissingOrInvalidFullScale(format!("invalid float '{val}'"))
                })?;
                if !parsed.is_finite() || parsed <= 0.0 {
                    return Err(BridgeError::MissingOrInvalidFullScale(format!(
                        "pa_per_full_scale must be finite and > 0 (got {val})"
                    )));
                }
                pa_per_full_scale = Some(parsed);
            }
            "--calibration-source" => {
                i += 1;
                let val = args.get(i).ok_or(BridgeError::MissingCalibrationSource)?;
                calibration_source = Some(CalibrationSource::parse(val)?);
            }
            "--field" => {
                i += 1;
                let val = args.get(i).ok_or(BridgeError::MissingField)?;
                field = Some(match val.as_str() {
                    "free" => SoundField::Free,
                    "diffuse" => SoundField::Diffuse,
                    other => return Err(BridgeError::InvalidField(other.to_owned())),
                });
            }
            "--time-skip" => {
                i += 1;
                let val = args
                    .get(i)
                    .ok_or(BridgeError::UnknownArgument("--time-skip".to_owned()))?;
                let parsed: f64 = val
                    .parse()
                    .map_err(|_| BridgeError::UnknownArgument(val.clone()))?;
                time_skip_s = parsed;
            }
            other if !other.starts_with('-') && wav_path.is_none() => {
                wav_path = Some(PathBuf::from(other));
            }
            other => return Err(BridgeError::UnknownArgument(other.to_owned())),
        }
        i += 1;
    }

    let pa_per_full_scale = pa_per_full_scale.ok_or_else(|| {
        BridgeError::MissingOrInvalidFullScale("missing --pa-per-full-scale".to_owned())
    })?;
    let calibration_source = calibration_source.ok_or(BridgeError::MissingCalibrationSource)?;
    let field = field.ok_or(BridgeError::MissingField)?;
    let wav_path = wav_path.ok_or(BridgeError::MissingWavPath)?;

    Ok(BridgeArgs {
        wav_path,
        pa_per_full_scale,
        calibration_source,
        field,
        time_skip_s,
    })
}

/// Decoded mono WAV buffer.
#[derive(Debug, Clone)]
pub struct DecodedWav {
    /// Sample rate [Hz].
    pub sample_rate: u32,
    /// Normalized digital samples widened to `f64` (prior to `pa_per_full_scale` scaling).
    pub samples_unit: Vec<f64>,
    /// SHA-256 hex digest of the raw WAV bytes.
    pub wav_sha256: String,
}

/// Parse a bounded mono PCM16 (format 1) or IEEE float32 (format 3) RIFF WAVE byte slice.
pub fn parse_mono_wav(bytes: &[u8]) -> Result<DecodedWav, BridgeError> {
    if bytes.is_empty() {
        return Err(BridgeError::WavFormat("empty wav buffer"));
    }
    if bytes.len() > MAX_WAV_BYTES {
        return Err(BridgeError::WavOversize {
            bytes: bytes.len(),
            max_bytes: MAX_WAV_BYTES,
        });
    }
    if bytes.len() < 44 || &bytes[0..4] != b"RIFF" || &bytes[8..12] != b"WAVE" {
        return Err(BridgeError::WavFormat("missing RIFF/WAVE header"));
    }

    let mut pos = 12usize;
    let mut fmt_info: Option<(u16, u16, u32, u16)> = None;
    let mut data_slice: Option<&[u8]> = None;

    while pos + 8 <= bytes.len() {
        let chunk_id = &bytes[pos..pos + 4];
        let chunk_len = u32::from_le_bytes([
            bytes[pos + 4],
            bytes[pos + 5],
            bytes[pos + 6],
            bytes[pos + 7],
        ]) as usize;
        pos += 8;
        if pos + chunk_len > bytes.len() {
            return Err(BridgeError::WavFormat("truncated wav chunk"));
        }
        let chunk = &bytes[pos..pos + chunk_len];
        if chunk_id == b"fmt " {
            if chunk.len() < 16 {
                return Err(BridgeError::WavFormat("fmt chunk shorter than 16 bytes"));
            }
            let format_tag = u16::from_le_bytes([chunk[0], chunk[1]]);
            let channels = u16::from_le_bytes([chunk[2], chunk[3]]);
            let sample_rate = u32::from_le_bytes([chunk[4], chunk[5], chunk[6], chunk[7]]);
            let bits_per_sample = u16::from_le_bytes([chunk[14], chunk[15]]);
            if channels != 1 {
                return Err(BridgeError::Multichannel { channels });
            }
            if sample_rate == 0 {
                return Err(BridgeError::WavFormat("sample_rate must be positive"));
            }
            fmt_info = Some((format_tag, channels, sample_rate, bits_per_sample));
        } else if chunk_id == b"data" {
            data_slice = Some(chunk);
        }
        // RIFF chunks are word-aligned (padded to even byte boundary).
        pos += chunk_len + (chunk_len & 1);
    }

    let (format_tag, _channels, sample_rate, bits_per_sample) =
        fmt_info.ok_or(BridgeError::WavFormat("missing fmt chunk"))?;
    let raw_data = data_slice.ok_or(BridgeError::WavFormat("missing data chunk"))?;

    let samples_unit: Vec<f64> = match (format_tag, bits_per_sample) {
        (1, 16) => {
            if raw_data.len() % 2 != 0 {
                return Err(BridgeError::WavFormat("unaligned pcm16 data chunk"));
            }
            raw_data
                .as_chunks::<2>()
                .0
                .iter()
                .map(|b| f64::from(i16::from_le_bytes(*b)) / 32768.0)
                .collect()
        }
        (3, 32) => {
            if raw_data.len() % 4 != 0 {
                return Err(BridgeError::WavFormat("unaligned float32 data chunk"));
            }
            raw_data
                .as_chunks::<4>()
                .0
                .iter()
                .map(|b| f64::from(f32::from_bits(u32::from_le_bytes(*b))))
                .collect()
        }
        _ => {
            return Err(BridgeError::WavFormat(
                "unsupported wav encoding (only mono PCM16 format=1 and float32 format=3 are supported)",
            ));
        }
    };

    if samples_unit.is_empty() {
        return Err(BridgeError::WavFormat("data chunk contains zero samples"));
    }

    Ok(DecodedWav {
        sample_rate,
        samples_unit,
        wav_sha256: sha256_hex(bytes),
    })
}

/// Evaluated psychoacoustic descriptor report (`fs-psycho-descriptors/1`).
#[derive(Debug, Clone, PartialEq)]
pub struct PsychoDescriptorReport {
    /// SHA-256 hex digest of the WAV file.
    pub wav_sha256: String,
    /// Sample rate [Hz].
    pub sample_rate: u32,
    /// Sound field (`Free` or `Diffuse`).
    pub field: SoundField,
    /// Full-scale sound pressure [Pa].
    pub pa_per_full_scale: f64,
    /// Calibration provenance (`Sidecar` or `Assumed`).
    pub calibration_source: CalibrationSource,
    /// Whether absolute level is unverified (`true` for `Assumed`).
    pub level_unverified: bool,
    /// Stationary Zwicker loudness [sone] (ISO 532-1).
    pub loudness_sone: Option<f64>,
    /// Time-varying loudness N5 percentile [sone] (ISO 532-1).
    pub n5_sone: Option<f64>,
    /// Time-varying loudness N50 median percentile [sone] (ISO 532-1).
    pub n50_sone: Option<f64>,
    /// DIN 45692 sharpness [acum].
    pub sharpness_acum: Option<f64>,
    /// Daniel-Weber roughness [asper] averaged over consecutive `DW_BLOCK` blocks.
    pub roughness_asper: Option<f64>,
    /// ECMA-74 tone-to-noise ratio total over prominent tones [dB].
    pub tnr_db: Option<f64>,
    /// ECMA-74 prominence ratio total over prominent tones [dB].
    pub pr_db: Option<f64>,
    /// Power-of-two sample count analyzed by ECMA-74 tonality (`None` if refused).
    pub tonality_block: Option<usize>,
    /// Fluctuation strength [vacil] over the first `FS_FRAME` (131,072) samples.
    pub fluctuation_vacil: Option<f64>,
    /// Timbre-toolbox log-attack-time [log10 s].
    pub log_attack_time: Option<f64>,
    /// Per-metric typed refusal messages (`<metric>` -> `<PsychoError Display>`).
    pub refusals: BTreeMap<&'static str, String>,
}

/// ISO 532-1 `f_calc_percentile` estimator for percentile `p_percent` (e.g. 50.0 for N50).
fn percentile_iso(series: &[f64], p_percent: f64) -> f64 {
    let n = series.len();
    let mut sorted = series.to_vec();
    sorted.sort_by(f64::total_cmp);
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    let np = (((1.0 - p_percent / 100.0) * n as f64) as usize)
        .max(1)
        .min(n - 1);
    f64::midpoint(sorted[np - 1], sorted[np])
}

/// Evaluate all `fs-psycho` descriptors on `pcm_pa` (sound pressure in Pa),
/// recording any individual metric's [`PsychoError`] in `refusals`.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn evaluate_from_pa(
    pcm_pa: &[f64],
    sample_rate: u32,
    field: SoundField,
    pa_per_full_scale: f64,
    calibration_source: CalibrationSource,
    time_skip_s: f64,
    wav_sha256: String,
) -> PsychoDescriptorReport {
    let sr = f64::from(sample_rate);
    let mut refusals: BTreeMap<&'static str, String> = BTreeMap::new();

    let (loudness_sone, sharpness_acum) =
        match loudness_stationary_from_pcm(pcm_pa, sr, time_skip_s, field) {
            Ok(stat) => {
                let sones = Some(stat.sones);
                let sharp = match sharpness_din(&stat.specific) {
                    Ok(s) => Some(s),
                    Err(e) => {
                        refusals.insert("sharpness_acum", e.to_string());
                        None
                    }
                };
                (sones, sharp)
            }
            Err(e) => {
                let msg = e.to_string();
                refusals.insert("loudness_sone", msg.clone());
                refusals.insert("sharpness_acum", msg);
                (None, None)
            }
        };

    let (n5_sone, n50_sone) = match loudness_time_varying(pcm_pa, sr, field) {
        Ok(tv) => (Some(tv.n5), Some(percentile_iso(&tv.loudness, 50.0))),
        Err(e) => {
            let msg = e.to_string();
            refusals.insert("n5_sone", msg.clone());
            refusals.insert("n50_sone", msg);
            (None, None)
        }
    };

    let n_blocks = pcm_pa.len() / DW_BLOCK;
    let roughness_asper = if n_blocks == 0 {
        let err = PsychoError::DegenerateSignal {
            what: "no whole roughness block in the signal",
        };
        refusals.insert("roughness_asper", err.to_string());
        None
    } else {
        let mut r_sum = 0.0f64;
        let mut failed: Option<PsychoError> = None;
        for b in 0..n_blocks {
            let start = b * DW_BLOCK;
            match roughness_dw_block(&pcm_pa[start..start + DW_BLOCK], sr) {
                Ok(r) => r_sum += r,
                Err(e) => {
                    failed = Some(e);
                    break;
                }
            }
        }
        if let Some(e) = failed {
            refusals.insert("roughness_asper", e.to_string());
            None
        } else {
            #[allow(clippy::cast_precision_loss)]
            Some(r_sum / n_blocks as f64)
        }
    };

    let tonality_slice = if pcm_pa.len() < 4096 || pcm_pa.len().is_power_of_two() {
        pcm_pa
    } else {
        let pow2 = 1usize << (usize::BITS - 1 - pcm_pa.len().leading_zeros());
        &pcm_pa[..pow2]
    };

    let mut tonality_block = None;
    let tnr_db = match tone_to_noise_ecma(tonality_slice, sr) {
        Ok(t) => {
            tonality_block = Some(tonality_slice.len());
            Some(t.total_db)
        }
        Err(e) => {
            refusals.insert("tnr_db", e.to_string());
            None
        }
    };
    let pr_db = match prominence_ratio_ecma(tonality_slice, sr) {
        Ok(t) => {
            tonality_block = Some(tonality_slice.len());
            Some(t.total_db)
        }
        Err(e) => {
            refusals.insert("pr_db", e.to_string());
            None
        }
    };

    let fluctuation_vacil = match fluctuation_strength_frame(pcm_pa, sr) {
        Ok(fs) => {
            let _ = FS_FRAME;
            Some(fs)
        }
        Err(e) => {
            refusals.insert("fluctuation_vacil", e.to_string());
            None
        }
    };

    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    let lat_window = ((0.01 * sr).round() as usize).max(1);
    let log_attack_time_val = match log_attack_time(pcm_pa, sr, lat_window) {
        Ok(lat) => Some(lat),
        Err(e) => {
            refusals.insert("log_attack_time", e.to_string());
            None
        }
    };

    PsychoDescriptorReport {
        wav_sha256,
        sample_rate,
        field,
        pa_per_full_scale,
        calibration_source,
        level_unverified: calibration_source.level_unverified(),
        loudness_sone,
        n5_sone,
        n50_sone,
        sharpness_acum,
        roughness_asper,
        tnr_db,
        pr_db,
        tonality_block,
        fluctuation_vacil,
        log_attack_time: log_attack_time_val,
        refusals,
    }
}

impl PsychoDescriptorReport {
    /// Serialize as a single canonical JSON line using Rust's round-trip-exact `f64` formatting.
    #[must_use]
    pub fn to_json_line(&self) -> String {
        fn opt_f64(v: Option<f64>) -> String {
            match v {
                Some(x) if x.is_finite() => format!("{x}"),
                _ => "null".to_owned(),
            }
        }
        fn opt_usize(v: Option<usize>) -> String {
            match v {
                Some(x) => format!("{x}"),
                None => "null".to_owned(),
            }
        }
        let field_str = match self.field {
            SoundField::Free => "free",
            SoundField::Diffuse => "diffuse",
        };
        let mut refusals_entries: Vec<String> = Vec::with_capacity(self.refusals.len());
        for (k, v) in &self.refusals {
            refusals_entries.push(format!("\"{}\":\"{}\"", json_escape(k), json_escape(v)));
        }
        let refusals_json = format!("{{{}}}", refusals_entries.join(","));
        format!(
            "{{\"schema\":\"fs-psycho-descriptors/1\",\"wav_sha256\":\"{}\",\"sample_rate\":{},\"field\":\"{}\",\"calibration\":{{\"pa_per_full_scale\":{},\"source\":\"{}\",\"level_unverified\":{}}},\"loudness_sone\":{},\"n5_sone\":{},\"n50_sone\":{},\"sharpness_acum\":{},\"roughness_asper\":{},\"tnr_db\":{},\"pr_db\":{},\"tonality_block\":{},\"fluctuation_vacil\":{},\"log_attack_time\":{},\"refusals\":{}}}",
            json_escape(&self.wav_sha256),
            self.sample_rate,
            field_str,
            self.pa_per_full_scale,
            self.calibration_source.as_str(),
            self.level_unverified,
            opt_f64(self.loudness_sone),
            opt_f64(self.n5_sone),
            opt_f64(self.n50_sone),
            opt_f64(self.sharpness_acum),
            opt_f64(self.roughness_asper),
            opt_f64(self.tnr_db),
            opt_f64(self.pr_db),
            opt_usize(self.tonality_block),
            opt_f64(self.fluctuation_vacil),
            opt_f64(self.log_attack_time),
            refusals_json,
        )
    }
}

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c => out.push(c),
        }
    }
    out
}

/// Pure-Rust FIPS 180-4 SHA-256 hex digest.
fn sha256_hex(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let bit_len = (data.len() as u64) * 8;
    let mut msg = data.to_vec();
    msg.push(0x80);
    while (msg.len() % 64) != 56 {
        msg.push(0x00);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in msg.as_chunks::<64>().0 {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                chunk[4 * i],
                chunk[4 * i + 1],
                chunk[4 * i + 2],
                chunk[4 * i + 3],
            ]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let temp1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }

    let mut out = String::with_capacity(64);
    for word in h {
        use core::fmt::Write as _;
        let _ = write!(out, "{word:08x}");
    }
    out
}

/// Run the full WAV -> `PsychoDescriptorReport` pipeline from raw WAV bytes.
pub fn evaluate_wav_bytes(
    wav_bytes: &[u8],
    pa_per_full_scale: f64,
    calibration_source: CalibrationSource,
    field: SoundField,
    time_skip_s: f64,
) -> Result<PsychoDescriptorReport, BridgeError> {
    if !pa_per_full_scale.is_finite() || pa_per_full_scale <= 0.0 {
        return Err(BridgeError::MissingOrInvalidFullScale(format!(
            "pa_per_full_scale must be finite and > 0 (got {pa_per_full_scale})"
        )));
    }
    let decoded = parse_mono_wav(wav_bytes)?;
    let pcm_pa: Vec<f64> = decoded
        .samples_unit
        .iter()
        .map(|&s| s * pa_per_full_scale)
        .collect();
    Ok(evaluate_from_pa(
        &pcm_pa,
        decoded.sample_rate,
        field,
        pa_per_full_scale,
        calibration_source,
        time_skip_s,
        decoded.wav_sha256,
    ))
}

fn main() -> ExitCode {
    let args = match parse_args(std::env::args().skip(1)) {
        Ok(a) => a,
        Err(e) => {
            println!(
                "{{\"schema\":\"fs-psycho-descriptors/1\",\"error\":\"{}\"}}",
                json_escape(&e.to_string())
            );
            return ExitCode::from(2);
        }
    };
    let bytes = match std::fs::read(&args.wav_path) {
        Ok(b) => b,
        Err(e) => {
            let err = BridgeError::Io(e.to_string());
            println!(
                "{{\"schema\":\"fs-psycho-descriptors/1\",\"error\":\"{}\"}}",
                json_escape(&err.to_string())
            );
            return ExitCode::from(2);
        }
    };
    match evaluate_wav_bytes(
        &bytes,
        args.pa_per_full_scale,
        args.calibration_source,
        args.field,
        args.time_skip_s,
    ) {
        Ok(report) => {
            println!("{}", report.to_json_line());
            ExitCode::SUCCESS
        }
        Err(e) => {
            println!(
                "{{\"schema\":\"fs-psycho-descriptors/1\",\"error\":\"{}\"}}",
                json_escape(&e.to_string())
            );
            ExitCode::from(2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::det;

    fn encode_f32_wav(sample_rate: u32, channels: u16, samples: &[f32]) -> Vec<u8> {
        let data_bytes = samples.len() * 4;
        let byte_rate = sample_rate * u32::from(channels) * 4;
        let block_align = channels * 4;
        let riff_size = (36 + data_bytes) as u32;
        let mut buf = Vec::with_capacity(44 + data_bytes);
        buf.extend_from_slice(b"RIFF");
        buf.extend_from_slice(&riff_size.to_le_bytes());
        buf.extend_from_slice(b"WAVE");
        buf.extend_from_slice(b"fmt ");
        buf.extend_from_slice(&16u32.to_le_bytes());
        buf.extend_from_slice(&3u16.to_le_bytes()); // IEEE float
        buf.extend_from_slice(&channels.to_le_bytes());
        buf.extend_from_slice(&sample_rate.to_le_bytes());
        buf.extend_from_slice(&byte_rate.to_le_bytes());
        buf.extend_from_slice(&block_align.to_le_bytes());
        buf.extend_from_slice(&32u16.to_le_bytes());
        buf.extend_from_slice(b"data");
        buf.extend_from_slice(&(data_bytes as u32).to_le_bytes());
        for &s in samples {
            buf.extend_from_slice(&s.to_bits().to_le_bytes());
        }
        buf
    }

    fn extract_json_f64(line: &str, key: &str) -> Option<f64> {
        let pat = format!("\"{key}\":");
        let start = line.find(&pat)? + pat.len();
        let rest = &line[start..];
        let end = rest.find([',', '}'])?;
        let token = rest[..end].trim();
        if token == "null" {
            None
        } else {
            Some(token.parse::<f64>().expect("valid f64 token"))
        }
    }

    /// G0: JSON values produced by the bridge parse back bit-identical to direct
    /// `fs-psycho` library calls on the same float32-quantized samples.
    #[test]
    fn bridge_is_bitwise_equal_to_direct_call() {
        let sr = 48_000u32;
        let sr_f64 = f64::from(sr);
        let n = FS_FRAME; // 131_072 samples: exercises all metrics including fluctuation_vacil
        let pa_per_full_scale = 2.0;
        let a = (det::sqrt(2.0) * 2.0e-5 * det::pow(10.0, 60.0 / 20.0)) / pa_per_full_scale;

        let samples_f32: Vec<f32> = (0..n)
            .map(|i| {
                let t = i as f64 / sr_f64;
                let env = (t / 0.05).min(1.0);
                let carrier = det::sin(2.0 * core::f64::consts::PI * 1000.0 * t)
                    + 0.25 * det::sin(2.0 * core::f64::consts::PI * 3000.0 * t);
                let am = 1.0 + 0.5 * det::sin(2.0 * core::f64::consts::PI * 70.0 * t);
                #[allow(clippy::cast_possible_truncation)]
                let s = (env * a * am * carrier) as f32;
                s
            })
            .collect();

        let wav = encode_f32_wav(sr, 1, &samples_f32);
        let report = evaluate_wav_bytes(
            &wav,
            pa_per_full_scale,
            CalibrationSource::Sidecar,
            SoundField::Free,
            0.0,
        )
        .expect("bridge evaluation");
        assert!(
            report.refusals.is_empty(),
            "unexpected refusals: {:?}",
            report.refusals
        );
        let json = report.to_json_line();

        // Direct calls on the exact f32-quantized Pa samples.
        let pcm_pa: Vec<f64> = samples_f32
            .iter()
            .map(|&s| f64::from(s) * pa_per_full_scale)
            .collect();
        let stat =
            loudness_stationary_from_pcm(&pcm_pa, sr_f64, 0.0, SoundField::Free).expect("stat");
        let sharp = sharpness_din(&stat.specific).expect("sharpness");
        let tv = loudness_time_varying(&pcm_pa, sr_f64, SoundField::Free).expect("tv");
        let n50 = percentile_iso(&tv.loudness, 50.0);
        let n_blocks = pcm_pa.len() / DW_BLOCK;
        let mut r_sum = 0.0;
        for b in 0..n_blocks {
            r_sum += roughness_dw_block(&pcm_pa[b * DW_BLOCK..(b + 1) * DW_BLOCK], sr_f64)
                .expect("roughness");
        }
        let rough = r_sum / n_blocks as f64;
        let tnr = tone_to_noise_ecma(&pcm_pa, sr_f64).expect("tnr").total_db;
        let pr = prominence_ratio_ecma(&pcm_pa, sr_f64).expect("pr").total_db;
        let fluc = fluctuation_strength_frame(&pcm_pa, sr_f64).expect("fluctuation");
        let lat = log_attack_time(&pcm_pa, sr_f64, 480).expect("lat");

        let checks = [
            ("loudness_sone", stat.sones),
            ("n5_sone", tv.n5),
            ("n50_sone", n50),
            ("sharpness_acum", sharp),
            ("roughness_asper", rough),
            ("tnr_db", tnr),
            ("pr_db", pr),
            ("fluctuation_vacil", fluc),
            ("log_attack_time", lat),
        ];
        for (key, direct) in checks {
            let parsed = extract_json_f64(&json, key)
                .unwrap_or_else(|| panic!("missing numeric key {key} in {json}"));
            assert_eq!(
                parsed.to_bits(),
                direct.to_bits(),
                "{key} parsed {parsed:.17} != direct {direct:.17}"
            );
        }
    }

    /// G0: Missing or invalid `--pa-per-full-scale`, `--field`, or `--calibration-source` refuses.
    #[test]
    fn refuses_without_full_scale() {
        let err = parse_args(["--wav", "x.wav", "--calibration-source", "sidecar", "--field", "free"])
            .expect_err("missing --pa-per-full-scale must refuse");
        assert!(matches!(err, BridgeError::MissingOrInvalidFullScale(_)));
        assert!(
            err.to_string()
                .contains(&PsychoError::UncalibratedAbsolute.to_string())
        );

        for bad in ["0", "-1.0", "NaN", "inf"] {
            let err = parse_args([
                "--wav",
                "x.wav",
                "--pa-per-full-scale",
                bad,
                "--calibration-source",
                "sidecar",
                "--field",
                "free",
            ])
            .expect_err("invalid --pa-per-full-scale must refuse");
            assert!(matches!(err, BridgeError::MissingOrInvalidFullScale(_)));
        }

        assert!(matches!(
            parse_args(["--wav", "x.wav", "--pa-per-full-scale", "1.0", "--field", "free"]),
            Err(BridgeError::MissingCalibrationSource)
        ));
        assert!(matches!(
            parse_args([
                "--wav",
                "x.wav",
                "--pa-per-full-scale",
                "1.0",
                "--calibration-source",
                "sidecar"
            ]),
            Err(BridgeError::MissingField)
        ));
    }

    /// G0: `DW_BLOCK - 1` (8191) refuses roughness, while `DW_BLOCK` (8192) passes (cap / cap+1).
    #[test]
    fn refuses_short_block_for_roughness() {
        let sr = 48_000u32;
        let samples_short: Vec<f32> = (0..(DW_BLOCK - 1))
            .map(|i| {
                #[allow(clippy::cast_possible_truncation)]
                let s = (0.05 * det::sin(2.0 * core::f64::consts::PI * 1000.0 * i as f64 / 48_000.0))
                    as f32;
                s
            })
            .collect();
        let wav_short = encode_f32_wav(sr, 1, &samples_short);
        let rep_short = evaluate_wav_bytes(
            &wav_short,
            1.0,
            CalibrationSource::Assumed,
            SoundField::Free,
            0.0,
        )
        .expect("short wav evaluates with per-metric refusal");
        assert!(rep_short.level_unverified);
        assert!(rep_short.roughness_asper.is_none());
        assert!(
            rep_short.refusals.contains_key("roughness_asper"),
            "DW_BLOCK - 1 must record roughness_asper refusal: {:?}",
            rep_short.refusals
        );

        let mut samples_ok = samples_short;
        samples_ok.push(0.01);
        assert_eq!(samples_ok.len(), DW_BLOCK);
        let wav_ok = encode_f32_wav(sr, 1, &samples_ok);
        let rep_ok = evaluate_wav_bytes(
            &wav_ok,
            1.0,
            CalibrationSource::Assumed,
            SoundField::Free,
            0.0,
        )
        .expect("DW_BLOCK wav evaluates");
        assert!(rep_ok.roughness_asper.is_some());
        assert!(!rep_ok.refusals.contains_key("roughness_asper"));
    }

    /// G0: WAV parser rejects multichannel and oversized inputs.
    #[test]
    fn wav_parser_rejects_multichannel_and_oversize() {
        let stereo = encode_f32_wav(48_000, 2, &[0.1, 0.1, -0.1, -0.1]);
        assert_eq!(
            parse_mono_wav(&stereo).expect_err("stereo must refuse"),
            BridgeError::Multichannel { channels: 2 }
        );

        let oversize = vec![0u8; MAX_WAV_BYTES + 1];
        assert_eq!(
            parse_mono_wav(&oversize).expect_err("oversize must refuse"),
            BridgeError::WavOversize {
                bytes: MAX_WAV_BYTES + 1,
                max_bytes: MAX_WAV_BYTES
            }
        );
    }
}
