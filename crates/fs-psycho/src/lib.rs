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
//! stationary loudness from band levels or from 48 kHz PCM through
//! the reference filterbank ([`signal`]); TIME-VARYING loudness with
//! Nmax/N5 ([`signal::loudness_time_varying`]); the verified phon
//! conversion ([`signal::phon_from_sone`]); Daniel-Weber roughness
//! per analysis block ([`roughness`]); ECMA-74 discrete-tone
//! tonality ([`tonality`]); fluctuation strength from the published
//! Osses/Garcia/Kohlrausch model equations ([`fluctuation`] — no
//! open reference code is license-compatible, so the equations were
//! implemented directly and validated against the published
//! Fastl-Zwicker values).

pub mod dw_tables;
pub mod filter_tables;
pub mod fluctuation;
pub mod receipt;
pub mod roughness;
pub mod signal;
pub mod tables;
pub mod tonality;

use fs_math::det;
use tables::{A0, DCB, DDF, DLL, LTQ, RAP, RNS, USL, ZUP};

/// Number of third-octave input bands (25 Hz .. 12.5 kHz).
pub const N_THIRD_OCTAVE_BANDS: usize = 28;
/// Specific-loudness grid: 0.1-Bark steps to 24 Bark.
pub const N_BARK_STEPS: usize = 240;
/// Conventional reference pressure for airborne SPL, in pascals (20 µPa).
/// Other media or measurement conventions can supply their own reference to
/// [`spl_from_pressure_rms`] and [`Calibration::from_pressure_scale`].
pub const AIR_REFERENCE_PRESSURE_PA: f64 = 20.0e-6;

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
    /// dB SPL of a digital sine with PEAK amplitude 1 (RMS = 1/sqrt(2)).
    /// This is not the level of a constant or square wave with magnitude 1.
    pub db_spl_at_full_scale: f64,
}

impl Calibration {
    /// Construct a full-scale-sine calibration from the pressure represented
    /// by digital sample 1 and an explicit RMS reference pressure, both in Pa.
    ///
    /// `pa_per_full_scale` is a linear sample scale, not an RMS sine level:
    /// `p[j] = pa_per_full_scale * pcm[j]`. Thus the calibration includes the
    /// sine's 1/sqrt(2) RMS factor exactly once. It does not certify a sensor or
    /// physical model; callers supply the provenance of their pressure scale.
    /// Use [`AIR_REFERENCE_PRESSURE_PA`] for conventional airborne dB SPL.
    pub fn from_pressure_scale(
        pa_per_full_scale: f64,
        reference_pressure_pa: f64,
    ) -> Result<Self, PsychoError> {
        admit_positive_pressure(pa_per_full_scale)?;
        admit_positive_pressure(reference_pressure_pa)?;
        // Separate logarithms avoid overflowing the ratio or losing a tiny
        // representable pressure while converting peak to RMS.
        Ok(Self {
            db_spl_at_full_scale: 20.0 / core::f64::consts::LN_10
                * (det::ln(pa_per_full_scale)
                    - det::ln(reference_pressure_pa)
                    - 0.5 * core::f64::consts::LN_2),
        })
    }
}

/// Typed refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PsychoError {
    /// An explicit physical reference or scale is outside its domain.
    InvalidParameter {
        /// The violated condition.
        what: &'static str,
    },
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
    /// A sampling rate the method's coefficient tables do not cover —
    /// refused rather than resampled or extrapolated (the ISO
    /// reference ships 48 kHz filter tables only).
    UnsupportedRate {
        /// What was required.
        what: &'static str,
    },
}

impl core::fmt::Display for PsychoError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PsychoError::InvalidParameter { what } => write!(f, "invalid parameter: {what}"),
            PsychoError::Shape { what } => write!(f, "shape mismatch: {what}"),
            PsychoError::NonFinite { what } => write!(f, "non-finite input: {what}"),
            PsychoError::UncalibratedAbsolute => write!(
                f,
                "absolute-level metric requested without calibration: refusing"
            ),
            PsychoError::DegenerateSignal { what } => write!(f, "degenerate signal: {what}"),
            PsychoError::UnsupportedRate { what } => {
                write!(f, "unsupported sampling rate: {what}")
            }
        }
    }
}

impl std::error::Error for PsychoError {}

/// A stationary loudness result in sones. (The phon conversion lives
/// in [`signal::phon_from_sone`], ported from the reference's own
/// `f_sone_to_phon` — it was absent until that verified source was
/// found; nothing here was recalled from memory.)
#[derive(Debug, Clone)]
pub struct Loudness {
    /// Total loudness [sone].
    pub sones: f64,
    /// Specific loudness on the 0.1-Bark grid [sone/Bark].
    pub specific: Vec<f64>,
}

/// numpy-style clamped linear interpolation over a table (used by
/// the fluctuation chain; `roughness` keeps its own byte-identical
/// private copy, transcribed with its reference).
pub(crate) fn tables_interp(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[xs.len() - 1] {
        return ys[ys.len() - 1];
    }
    let mut i = 0;
    while xs[i + 1] < x {
        i += 1;
    }
    let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    ys[i] + t * (ys[i + 1] - ys[i])
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
    let core = core_loudness_21(third_octave_levels_db, field);
    let (sones, specific) = slopes_240(&core);
    Ok(Loudness { sones, specific })
}

/// Core loudness per critical band (21 entries, last zero) from one
/// frame of 28 third-octave levels — the reference's
/// `f_corr_third_octave_intensities` + `f_calc_lcbs` +
/// `f_calc_core_loudness` + `f_corr_loudness` chain, shared between
/// the stationary and time-varying methods.
pub(crate) fn core_loudness_21(third_octave_levels_db: &[f64], field: SoundField) -> [f64; 21] {
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
            // C literal .0635f: float32-rounded, like the tables.
            let mp1 = f64::from(0.0635f32) * det::pow(10.0, 0.025 * LTQ[idx]);
            let mp2 = det::pow(1.0 - s + s * det::pow(10.0, 0.1 * (le - LTQ[idx])), 0.25) - 1.0;
            core[idx] = (mp1 * mp2).max(0.0);
        }
    }
    // --- Lowest-band threshold correction. ---
    // C literals 0.4f/0.32f: float32-rounded, like the tables.
    let corr_cl = f64::from(0.4f32) + f64::from(0.32f32) * det::pow(core[0], 0.2);
    if corr_cl < 1.0 {
        core[0] *= corr_cl;
    }
    core
}

/// Specific-loudness pattern + total loudness from one frame of core
/// loudness — the reference's `f_calc_slopes`.
#[allow(clippy::too_many_lines)] // one reference routine, told in its order
pub(crate) fn slopes_240(core: &[f64; 21]) -> (f64, Vec<f64>) {
    // --- Slopes: specific loudness pattern + total loudness. ---
    let mut specific = vec![0.0f64; N_BARK_STEPS];
    let mut loud = 0.0f64;
    // The reference's Bark-walk constants are float LITERALS widened
    // to double (0.1f, .0001f) — replicated exactly; the walk's
    // comparisons are sensitive to them.
    let step = f64::from(0.1f32);
    let mut n1 = 0.0f64;
    let mut z = step;
    let mut z1 = 0.0f64;
    let mut idx_rns = 0usize;
    let mut idx_ns = 0usize;
    for (idx_cl, &core_l) in core.iter().enumerate() {
        let zup = ZUP[idx_cl] + f64::from(0.0001f32);
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
                    zk += step;
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
                    zk += step;
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
    (loud, specific)
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
    if !sample_rate.is_finite() || sample_rate <= 0.0 {
        return Err(PsychoError::NonFinite {
            what: "sample rate",
        });
    }
    let scale = finite_nonzero_peak(pcm)?;
    // Moving-average envelope of |x| (deterministic, no allocation
    // tricks). The descriptor is level-relative: normalize BEFORE summing,
    // so changing physical units or recording gain cannot overflow its
    // envelope or underflow its threshold crossings.
    let mut env = vec![0.0f64; pcm.len() - env_window + 1];
    let mut acc: f64 = pcm[..env_window].iter().map(|v| v.abs() / scale).sum();
    env[0] = acc / env_window as f64;
    for i in 1..env.len() {
        acc += pcm[i + env_window - 1].abs() / scale - pcm[i - 1].abs() / scale;
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
    let rise_samples = ((i90.max(i10 + 1)) - i10) as f64;
    Ok((det::ln(rise_samples) - det::ln(sample_rate)) / core::f64::consts::LN_10)
}

/// The only PCM-to-ABSOLUTE bridge in v1: RMS level in dB SPL, which
/// REQUIRES calibration (the typed uncalibrated-absolute refusal the
/// honest-scope law demands).
/// The supplied samples are evaluated together, retaining coherent source
/// interference and the actual waveform's crest factor. No peak normalization
/// or mean subtraction is performed. A constant offset contributes to RMS;
/// callers must provide the intended acoustic signal rather than ambient
/// absolute pressure or an unwanted sensor offset.
///
/// # Errors
/// [`PsychoError::UncalibratedAbsolute`] when `calibration` is
/// `None`; degenerate-signal refusals.
pub fn spl_from_pcm_rms(pcm: &[f64], calibration: Option<Calibration>) -> Result<f64, PsychoError> {
    let Some(cal) = calibration else {
        return Err(PsychoError::UncalibratedAbsolute);
    };
    if !cal.db_spl_at_full_scale.is_finite() {
        return Err(PsychoError::NonFinite {
            what: "full-scale-sine calibration",
        });
    }
    // Full-scale SINE reference: rms_fs = 1/sqrt(2). Work in logarithms
    // because the level is representable even when x² or physical RMS is not.
    Ok(cal.db_spl_at_full_scale
        + 20.0 / core::f64::consts::LN_10 * (log_rms(pcm)? + 0.5 * core::f64::consts::LN_2))
}

/// Unweighted equivalent SPL of acoustic pressure samples in Pa over the
/// supplied record: `20 log10(p_rms / reference_pressure_pa)`.
///
/// The reference is explicit; [`AIR_REFERENCE_PRESSURE_PA`] is conventional
/// for airborne sound. Input is acoustic pressure relative to ambient, not
/// thermodynamic absolute pressure. The caller owns that physical meaning and
/// any microphone calibration. There is no weighting filter, DC removal,
/// inferred calibration, or peak normalization. Silence and empty records
/// refuse, as in [`spl_from_pcm_rms`].
///
/// Definition: NPL, "Sound Pressure Level and Intensity Level",
/// <https://resource.npl.co.uk/acoustics/techguides/concepts/spl.html>.
pub fn spl_from_pressure_rms(
    pressure_pa: &[f64],
    reference_pressure_pa: f64,
) -> Result<f64, PsychoError> {
    admit_positive_pressure(reference_pressure_pa)?;
    Ok(20.0 / core::f64::consts::LN_10 * (log_rms(pressure_pa)? - det::ln(reference_pressure_pa)))
}

fn admit_positive_pressure(pressure_pa: f64) -> Result<(), PsychoError> {
    if !(pressure_pa.is_finite() && pressure_pa > 0.0) {
        return Err(PsychoError::InvalidParameter {
            what: "pressure scale and reference must be positive and finite",
        });
    }
    Ok(())
}

fn finite_nonzero_peak(samples: &[f64]) -> Result<f64, PsychoError> {
    if samples.is_empty() {
        return Err(PsychoError::DegenerateSignal { what: "empty pcm" });
    }
    let mut peak = 0.0_f64;
    for &sample in samples {
        if !sample.is_finite() {
            return Err(PsychoError::NonFinite { what: "pcm sample" });
        }
        peak = peak.max(sample.abs());
    }
    if peak == 0.0 {
        return Err(PsychoError::DegenerateSignal { what: "silent pcm" });
    }
    Ok(peak)
}

// Scaled sum of squares; O(N) work, constant storage, no intermediate x² in
// physical units. At least one normalized square is 1, including subnormal
// nonzero input, so finite signal amplitude cannot be misreported as silence.
fn log_rms(samples: &[f64]) -> Result<f64, PsychoError> {
    let scale = finite_nonzero_peak(samples)?;
    let sum_squares: f64 = samples.iter().map(|&x| (x / scale).powi(2)).sum();
    Ok(det::ln(scale) + 0.5 * det::ln(sum_squares / samples.len() as f64))
}

/// The listening law, pinned as data so a test can assert it never
/// disappears from the crate.
pub const LISTENING_LAW: &str = "psychoacoustic metrics are never a substitute for human listening";
