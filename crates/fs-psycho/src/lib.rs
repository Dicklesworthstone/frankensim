//! # fs-psycho — standards-anchored psychoacoustic metrics
//!
//! Perceptual objective axes for the simulate-vs-listen loop:
//! stationary Zwicker loudness (ISO 532-1, third-octave-level method,
//! data tables transcribed MECHANICALLY from the standard's own free
//! reference implementation — see [`tables`]), DIN 45692 sharpness
//! over the specific-loudness pattern, and log-attack-time.
//!
//! CALIBRATION CONTRACT (explicit in types): metrics carrying
//! absolute meaning (loudness in sones) require absolute band levels
//! in dB SPL; the PCM-to-SPL bridge [`spl_from_pcm_rms`] REFUSES
//! without a [`Calibration`] — uncalibrated digital amplitude cannot
//! make an absolute claim. Level-relative metrics (sharpness in acum
//! is a ratio over the loudness pattern; log-attack-time is a time
//! measurement) work uncalibrated and their signatures say so.
//!
//! THESE METRICS ARE NEVER A SUBSTITUTE FOR HUMAN LISTENING. They are
//! objective columns for Pareto fits and listening-evidence tables;
//! a human ear adjudicates (the program's listening law — this
//! statement is load-bearing and pinned by a test).
//!
//! Honest scope (stated; the bead stays open until the rest lands):
//! stationary loudness from BAND LEVELS only; Daniel-Weber roughness
//! per analysis block ([`roughness`]); the PCM third-octave
//! filterbank, time-varying loudness, fluctuation strength, and
//! tonality are not yet implemented — no placeholder claims.

pub mod dw_tables;
pub mod roughness;
pub mod tables;

use fs_math::det;
use tables::{A0, DCB, DDF, DLL, LTQ, RAP, RNS, USL, ZUP};

/// Number of third-octave input bands (25 Hz .. 12.5 kHz).
pub const N_THIRD_OCTAVE_BANDS: usize = 28;
/// Specific-loudness grid: 0.1-Bark steps to 24 Bark.
pub const N_BARK_STEPS: usize = 240;

/// Sound field of the level measurement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SoundField {
    /// Frontal free field.
    Free,
    /// Diffuse field.
    Diffuse,
}

/// Digital-to-absolute calibration: the SPL a full-scale (|x| = 1)
/// sine produces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Calibration {
    /// dB SPL at digital full scale.
    pub db_spl_at_full_scale: f64,
}

/// Typed refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PsychoError {
    /// Band-level vector has the wrong length.
    Shape {
        /// What disagreed.
        what: &'static str,
    },
    /// Non-finite input.
    NonFinite {
        /// Where.
        what: &'static str,
    },
    /// An ABSOLUTE metric was requested without calibration — refused
    /// by name (the honest-scope law: uncalibrated digital amplitude
    /// carries no SPL meaning).
    UncalibratedAbsolute,
    /// Degenerate PCM (empty, or silent where a level is required).
    DegenerateSignal {
        /// Why.
        what: &'static str,
    },
}

impl core::fmt::Display for PsychoError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PsychoError::Shape { what } => write!(f, "shape mismatch: {what}"),
            PsychoError::NonFinite { what } => write!(f, "non-finite input: {what}"),
            PsychoError::UncalibratedAbsolute => write!(
                f,
                "absolute-level metric requested without calibration: refusing"
            ),
            PsychoError::DegenerateSignal { what } => write!(f, "degenerate signal: {what}"),
        }
    }
}

impl std::error::Error for PsychoError {}

/// A stationary loudness result. (Sones only: the ISO reference
/// implementation outputs sones; a phon conversion would need its own
/// verified source and is deliberately absent rather than recalled
/// from memory.)
#[derive(Debug, Clone)]
pub struct Loudness {
    /// Total loudness [sone].
    pub sones: f64,
    /// Specific loudness on the 0.1-Bark grid [sone/Bark].
    pub specific: Vec<f64>,
}

/// Stationary Zwicker loudness from 28 third-octave band levels
/// (25 Hz .. 12.5 kHz, dB SPL) — the ISO 532-1 clause-5 method,
/// ported statement-for-statement from the standard's reference
/// implementation with the tables in [`tables`].
///
/// # Errors
/// [`PsychoError::Shape`] / [`PsychoError::NonFinite`].
#[allow(clippy::too_many_lines)] // one clause-5 pipeline, told in the reference's order
pub fn loudness_stationary(
    third_octave_levels_db: &[f64],
    field: SoundField,
) -> Result<Loudness, PsychoError> {
    if third_octave_levels_db.len() != N_THIRD_OCTAVE_BANDS {
        return Err(PsychoError::Shape {
            what: "28 third-octave levels required (25 Hz .. 12.5 kHz)",
        });
    }
    if third_octave_levels_db.iter().any(|v| !v.is_finite()) {
        return Err(PsychoError::NonFinite {
            what: "third-octave level",
        });
    }
    // --- Low-frequency correction and intensities (bands 0..11). ---
    let mut intens = [0.0f64; 11];
    for i in 0..11 {
        let level = third_octave_levels_db[i];
        let mut range = 0usize;
        while level > RAP[range] - DLL[range * 11 + i] && range < 7 {
            range += 1;
        }
        let corr = level + DLL[range * 11 + i];
        intens[i] = det::pow(10.0, corr / 10.0);
    }
    // --- First three critical-band levels (LCB). ---
    let sum1: f64 = intens[0..6].iter().sum();
    let sum2: f64 = intens[6..9].iter().sum();
    let sum3: f64 = intens[9..11].iter().sum();
    let lcb = [sum1, sum2, sum3].map(|s| {
        if s > 0.0 {
            10.0 * det::ln(s) / det::ln(10.0)
        } else {
            0.0
        }
    });
    // --- Core loudness per critical band (21 entries, last zero). ---
    let mut core = [0.0f64; 21];
    for idx in 0..20 {
        let mut le = if idx < 3 {
            lcb[idx]
        } else {
            third_octave_levels_db[idx + 8]
        };
        le -= A0[idx];
        if field == SoundField::Diffuse {
            le += DDF[idx];
        }
        if le > LTQ[idx] {
            le -= DCB[idx];
            let s = 0.25;
            let mp1 = 0.0635 * det::pow(10.0, 0.025 * LTQ[idx]);
            let mp2 = det::pow(1.0 - s + s * det::pow(10.0, 0.1 * (le - LTQ[idx])), 0.25) - 1.0;
            core[idx] = (mp1 * mp2).max(0.0);
        }
    }
    // --- Lowest-band threshold correction. ---
    let corr_cl = 0.4 + 0.32 * det::pow(core[0], 0.2);
    if corr_cl < 1.0 {
        core[0] *= corr_cl;
    }
    // --- Slopes: specific loudness pattern + total loudness. ---
    let mut specific = vec![0.0f64; N_BARK_STEPS];
    let mut loud = 0.0f64;
    let mut n1 = 0.0f64;
    let mut z = 0.1f64;
    let mut z1 = 0.0f64;
    let mut idx_rns = 0usize;
    let mut idx_ns = 0usize;
    for idx_cl in 0..21 {
        let core_l = core[idx_cl];
        let zup = ZUP[idx_cl] + 0.0001;
        let idx_cbn = idx_cl.saturating_sub(1).min(7);
        let mut n2;
        loop {
            let next_band;
            if n1 > core_l {
                let usl = USL[idx_rns * 8 + idx_cbn];
                n2 = RNS[idx_rns].max(core_l);
                let mut dz = (n1 - n2) / usl;
                let mut z2 = z1 + dz;
                if z2 > zup {
                    next_band = true;
                    z2 = zup;
                    dz = z2 - z1;
                    n2 = n1 - dz * usl;
                } else {
                    next_band = false;
                }
                loud += dz * (n1 + n2) / 2.0;
                let mut zk = z;
                while zk <= z2 {
                    if idx_ns < N_BARK_STEPS {
                        specific[idx_ns] = n1 - (zk - z1) * usl;
                        idx_ns += 1;
                    }
                    zk += 0.1;
                }
                z = zk;
                z1 = z2;
                n1 = n2;
            } else {
                if n1 < core_l {
                    idx_rns = 0;
                    while idx_rns < 18 && RNS[idx_rns] >= core_l {
                        idx_rns += 1;
                    }
                }
                next_band = true;
                let z2 = zup;
                n2 = core_l;
                loud += n2 * (z2 - z1);
                let mut zk = z;
                while zk <= z2 {
                    if idx_ns < N_BARK_STEPS {
                        specific[idx_ns] = n2;
                        idx_ns += 1;
                    }
                    zk += 0.1;
                }
                z = zk;
                z1 = z2;
                n1 = n2;
            }
            while n2 <= RNS[idx_rns] && idx_rns < 17 {
                idx_rns += 1;
            }
            idx_rns = idx_rns.min(17);
            if next_band {
                break;
            }
        }
        if loud < 0.0 {
            loud = 0.0;
        }
    }
    Ok(Loudness {
        sones: loud,
        specific,
    })
}

/// DIN 45692 sharpness [acum] from a specific-loudness pattern on the
/// 0.1-Bark grid: `S = 0.11 * sum(N' g(z) z dz) / N` with
/// `g(z) = 1` for `z <= 15.8` and
/// `g(z) = 0.15 exp(0.42 (z - 15.8)) + 0.85` above (the standard's
/// weighting as published; cross-checked against the Apache-2.0
/// MoSQITo reference implementation). Level-RELATIVE: no calibration
/// needed — sharpness is a ratio over the pattern.
///
/// # Errors
/// [`PsychoError::Shape`] on a wrong grid;
/// [`PsychoError::DegenerateSignal`] on zero total loudness (a
/// sharpness of silence is undefined — refused, not fabricated).
pub fn sharpness_din(specific_loudness: &[f64]) -> Result<f64, PsychoError> {
    if specific_loudness.len() != N_BARK_STEPS {
        return Err(PsychoError::Shape {
            what: "specific loudness must be on the 240-step 0.1-Bark grid",
        });
    }
    let dz = 0.1;
    let mut total = 0.0f64;
    let mut weighted = 0.0f64;
    for (i, &ns) in specific_loudness.iter().enumerate() {
        let z = (i as f64 + 1.0) * dz;
        let g = if z <= 15.8 {
            1.0
        } else {
            0.15 * det::exp(0.42 * (z - 15.8)) + 0.85
        };
        total += ns * dz;
        weighted += ns * g * z * dz;
    }
    if total <= 0.0 {
        return Err(PsychoError::DegenerateSignal {
            what: "zero loudness has no sharpness",
        });
    }
    Ok(0.11 * weighted / total)
}

/// Log-attack-time [log10 seconds] per the timbre-toolbox definition:
/// `log10(t_90 - t_10)` where `t_10`/`t_90` are the first crossings
/// of 10%/90% of the peak of the amplitude envelope (here: rectified
/// signal smoothed by a moving average of `env_window` samples).
/// Level-relative: no calibration needed.
///
/// # Errors
/// [`PsychoError::DegenerateSignal`] on empty/silent input or a
/// window longer than the signal.
pub fn log_attack_time(
    pcm: &[f64],
    sample_rate: f64,
    env_window: usize,
) -> Result<f64, PsychoError> {
    if pcm.is_empty() || env_window == 0 || env_window > pcm.len() {
        return Err(PsychoError::DegenerateSignal {
            what: "empty signal or bad envelope window",
        });
    }
    if sample_rate.is_nan() || sample_rate <= 0.0 {
        return Err(PsychoError::NonFinite {
            what: "sample rate",
        });
    }
    // Moving-average envelope of |x| (deterministic, no allocation
    // tricks).
    let mut env = vec![0.0f64; pcm.len() - env_window + 1];
    let mut acc: f64 = pcm[..env_window].iter().map(|v| v.abs()).sum();
    env[0] = acc / env_window as f64;
    for i in 1..env.len() {
        acc += pcm[i + env_window - 1].abs() - pcm[i - 1].abs();
        env[i] = acc / env_window as f64;
    }
    let peak = env.iter().fold(0.0f64, |a, &v| a.max(v));
    if peak <= 0.0 {
        return Err(PsychoError::DegenerateSignal {
            what: "silent signal has no attack",
        });
    }
    let t_of =
        |threshold: f64| -> Option<usize> { env.iter().position(|&v| v >= threshold * peak) };
    let (Some(i10), Some(i90)) = (t_of(0.1), t_of(0.9)) else {
        return Err(PsychoError::DegenerateSignal {
            what: "envelope never crosses attack thresholds",
        });
    };
    let dt = ((i90.max(i10 + 1)) - i10) as f64 / sample_rate;
    Ok(det::ln(dt) / det::ln(10.0))
}

/// The only PCM-to-ABSOLUTE bridge in v1: RMS level in dB SPL, which
/// REQUIRES calibration (the typed uncalibrated-absolute refusal the
/// honest-scope law demands).
///
/// # Errors
/// [`PsychoError::UncalibratedAbsolute`] when `calibration` is
/// `None`; degenerate-signal refusals.
pub fn spl_from_pcm_rms(pcm: &[f64], calibration: Option<Calibration>) -> Result<f64, PsychoError> {
    let Some(cal) = calibration else {
        return Err(PsychoError::UncalibratedAbsolute);
    };
    if pcm.is_empty() {
        return Err(PsychoError::DegenerateSignal { what: "empty pcm" });
    }
    let ms: f64 = pcm.iter().map(|v| v * v).sum::<f64>() / pcm.len() as f64;
    if ms <= 0.0 {
        return Err(PsychoError::DegenerateSignal { what: "silent pcm" });
    }
    // Full-scale SINE reference: rms_fs = 1/sqrt(2).
    let rms = det::sqrt(ms);
    Ok(cal.db_spl_at_full_scale + 20.0 * det::ln(rms * det::sqrt(2.0)) / det::ln(10.0))
}

/// The listening law, pinned as data so a test can assert it never
/// disappears from the crate.
pub const LISTENING_LAW: &str = "psychoacoustic metrics are never a substitute for human listening";
