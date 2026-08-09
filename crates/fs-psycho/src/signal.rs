//! ISO 532-1 loudness from a PCM signal: the third-octave filterbank
//! path (clause 5 via filtered levels) and the TIME-VARYING method
//! (clause 6), ported statement-for-statement from the standard's
//! free reference implementation (Annex A.4 electronic insert,
//! `ISO_532-1.c` / `ISO_532-1_helper.c`), with the 48 kHz biquad
//! tables in [`crate::filter_tables`] transcribed mechanically.
//!
//! Scope facts, stated:
//! - Inputs are sound pressure in Pa (the caller's calibration turns
//!   digital samples into Pa; the uncalibrated-absolute refusal
//!   applies upstream at [`crate::spl_from_pcm_rms`]).
//! - The reference ships filter coefficients for EXACTLY 48 kHz;
//!   other rates REFUSE (`UnsupportedRate`) — no resampling, no
//!   coefficient redesign, no silent degradation.
//! - The time-varying loudness series is sampled at 2 kHz
//!   ([`SR_LEVEL`]), like the reference (decimation by 24).
//! - `phon_from_sone` is the reference's own `f_sone_to_phon`
//!   (previously deliberately absent while unverified). Its two
//!   branches meet at 1 sone with a small step (40.007 -> 40.000
//!   phon), a reference behavior kept for faithfulness.

use crate::filter_tables::{FILTER_DIFF, FILTER_GAIN, FILTER_REF};
use crate::{
    Loudness, N_THIRD_OCTAVE_BANDS, PsychoError, SoundField, core_loudness_21, slopes_240,
};
use fs_math::det;

/// The only sampling rate the reference filter tables cover [Hz].
pub const SIGNAL_SAMPLE_RATE: f64 = 48_000.0;
/// Sampling rate of the decimated level/loudness time series [Hz].
pub const SR_LEVEL: f64 = 2_000.0;
/// Decimation factor from signal rate to level rate.
const DEC_FACTOR: usize = 24;
/// Reference intensity (2e-5 Pa)^2.
const I_REF: f64 = 4e-10;
/// Additive floor before the level log (reference TINY_VALUE).
const TINY_VALUE: f64 = 1e-12;
/// Inner iterations of the nonlinear-decay and temporal-weighting
/// low-passes (virtual upsampling).
const NL_ITER: usize = 24;
const LP_ITER: usize = 24;
/// Time constants of the nonlinear temporal decay [s].
const T_SHORT: f64 = 0.005;
const T_LONG: f64 = 0.015;
const T_VAR: f64 = 0.075;

/// A time-varying loudness result (clause-6 method).
#[derive(Debug, Clone)]
pub struct TimeVaryingLoudness {
    /// Loudness time series [sone], sampled at [`SR_LEVEL`].
    pub loudness: Vec<f64>,
    /// Maximum of the series [sone].
    pub n_max: f64,
    /// N5: the loudness exceeded 5% of the time [sone] — the
    /// standard's single-number summary for time-varying sounds.
    pub n5: f64,
}

/// One biquad stage in the reference's direct-form-II arrangement
/// (`f_filter_2ndOrder`): coefficients enter as `REF - DIFF`, the
/// input is pre-scaled by the stage gain.
fn filter_2nd_order(input: &[f64], output: &mut [f64], coeffs: &[f64; 6], gain: f64) {
    let mut wn1 = 0.0f64;
    let mut wn2 = 0.0f64;
    for (x, y) in input.iter().zip(output.iter_mut()) {
        let wn0 = x * gain - coeffs[4] * wn1 - coeffs[5] * wn2;
        *y = coeffs[0] * wn0 + coeffs[1] * wn1 + coeffs[2] * wn2;
        wn2 = wn1;
        wn1 = wn0;
    }
}

/// First-order low-pass, in place (`f_lowpass` with aliased buffers —
/// the reference calls it input==output).
fn lowpass_in_place(x: &mut [f64], tau: f64, sample_rate: f64) {
    let a1 = det::exp(-1.0 / (sample_rate * tau));
    let b0 = 1.0 - a1;
    let mut y1 = 0.0f64;
    for v in x.iter_mut() {
        y1 = b0 * *v + a1 * y1;
        *v = y1;
    }
}

/// First-order low-pass with `LP_ITER` linear-interpolation
/// sub-steps per sample (`f_lowpass_intp`): each output sample is the
/// filter state after the FIRST sub-step; the remaining 23 sub-steps
/// advance the state toward the next input sample.
fn lowpass_interpolated(input: &[f64], output: &mut [f64], tau: f64, sample_rate: f64) {
    let a1 = det::exp(-1.0 / (sample_rate * LP_ITER as f64 * tau));
    let b0 = 1.0 - a1;
    let mut y1 = 0.0f64;
    let n = input.len();
    for t in 0..n {
        let x0 = input[t];
        y1 = b0 * x0 + a1 * y1;
        output[t] = y1;
        if t < n - 1 {
            let xd = (input[t + 1] - x0) / LP_ITER as f64;
            let mut x = x0;
            for _ in 1..LP_ITER {
                x += xd;
                y1 = b0 * x + a1 * y1;
            }
        }
    }
}

/// State + coefficients of the nonlinear temporal-decay element
/// (`NlLpData` / `f_init_nl_lp`), built for the virtually upsampled
/// rate `SR_LEVEL * NL_ITER`.
struct NlLp {
    b: [f64; 6],
    uo_last: f64,
    u2_last: f64,
}

impl NlLp {
    fn new(sample_rate: f64) -> Self {
        let delta_t = 1.0 / sample_rate;
        let p = (T_VAR + T_LONG) / (T_VAR * T_SHORT);
        let q = 1.0 / (T_SHORT * T_VAR);
        let lambda1 = -p / 2.0 + det::sqrt(p * p / 4.0 - q);
        let lambda2 = -p / 2.0 - det::sqrt(p * p / 4.0 - q);
        let den = T_VAR * (lambda1 - lambda2);
        let e1 = det::exp(lambda1 * delta_t);
        let e2 = det::exp(lambda2 * delta_t);
        NlLp {
            b: [
                (e1 - e2) / den,
                ((T_VAR * lambda2 + 1.0) * e1 - (T_VAR * lambda1 + 1.0) * e2) / den,
                ((T_VAR * lambda1 + 1.0) * e1 - (T_VAR * lambda2 + 1.0) * e2) / den,
                (T_VAR * lambda1 + 1.0) * (T_VAR * lambda2 + 1.0) * (e1 - e2) / den,
                det::exp(-delta_t / T_LONG),
                det::exp(-delta_t / T_VAR),
            ],
            uo_last: 0.0,
            u2_last: 0.0,
        }
    }

    /// One step of the nonlinear element (`f_nl_lp`), case structure
    /// exactly as the reference (including its 1e-5 equality band).
    fn step(&mut self, ui: f64) -> f64 {
        let (uo, u2);
        if ui < self.uo_last {
            if self.uo_last > self.u2_last {
                // Case 1.1: two-capacitor discharge.
                let mut u2c = self.uo_last * self.b[0] - self.u2_last * self.b[1];
                let mut uoc = self.uo_last * self.b[2] - self.u2_last * self.b[3];
                if uoc < ui {
                    uoc = ui;
                }
                if u2c > uoc {
                    u2c = uoc;
                }
                uo = uoc;
                u2 = u2c;
            } else {
                // Case 1.2: single-capacitor discharge.
                let mut uoc = self.uo_last * self.b[4];
                if uoc < ui {
                    uoc = ui;
                }
                uo = uoc;
                u2 = uoc;
            }
        } else if (ui - self.uo_last).abs() < 1e-5 {
            // Case 2: input steady.
            uo = ui;
            if uo > self.u2_last {
                u2 = (self.u2_last - ui) * self.b[5] + ui;
            } else {
                u2 = ui;
            }
        } else {
            // Case 3: input rising.
            uo = ui;
            u2 = (self.u2_last - ui) * self.b[5] + ui;
        }
        self.uo_last = uo;
        self.u2_last = u2;
        uo
    }
}

/// Nonlinear temporal decay over one core-loudness channel's time
/// series (`f_nl`): NL_ITER linear-interpolation sub-steps per
/// sample; only the first sub-step's output is stored.
fn nl_channel(x: &mut [f64]) {
    let mut nl = NlLp::new(SR_LEVEL * NL_ITER as f64);
    let n = x.len();
    for t in 0..n - 1 {
        let delta = (x[t + 1] - x[t]) / NL_ITER as f64;
        let mut ui = x[t];
        x[t] = nl.step(ui);
        ui += delta;
        for _ in 1..NL_ITER {
            let _ = nl.step(ui);
            ui += delta;
        }
    }
    x[n - 1] = nl.step(x[n - 1]);
}

/// Duration-dependent temporal weighting of the loudness series
/// (`f_temporal_weight_loudness`): 0.47 x (3.5 ms interpolated
/// low-pass) + 0.53 x (70 ms interpolated low-pass).
fn temporal_weight(loudness: &mut [f64]) {
    let mut l1 = vec![0.0f64; loudness.len()];
    let mut l2 = vec![0.0f64; loudness.len()];
    lowpass_interpolated(loudness, &mut l1, 3.5e-3, SR_LEVEL);
    lowpass_interpolated(loudness, &mut l2, 70e-3, SR_LEVEL);
    for ((v, a), b) in loudness.iter_mut().zip(&l1).zip(&l2) {
        *v = 0.47 * a + 0.53 * b;
    }
}

/// The standard's percentile estimator (`f_calc_percentile`, P = 5):
/// sort ascending, `Np = floor(0.95 n)`, average of the two samples
/// around the cut. Requires at least two frames (the reference
/// indexes out of bounds below that; the public API's two-frame
/// refusal guards it here).
fn percentile_n5(series: &[f64]) -> f64 {
    let n = series.len();
    debug_assert!(
        n >= 2,
        "percentile needs >= 2 frames (guarded by the public API)"
    );
    let mut sorted = series.to_vec();
    sorted.sort_by(f64::total_cmp);
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    let np = (((1.0 - 5.0 / 100.0) * n as f64) as usize).max(1);
    f64::midpoint(sorted[np - 1], sorted[np])
}

/// Validate the PCM + rate common to both signal paths.
fn validate_signal(pcm_pa: &[f64], sample_rate: f64) -> Result<(), PsychoError> {
    if sample_rate != SIGNAL_SAMPLE_RATE {
        return Err(PsychoError::UnsupportedRate {
            what: "the ISO 532-1 reference filter tables are 48 kHz only",
        });
    }
    if pcm_pa.iter().any(|v| !v.is_finite()) {
        return Err(PsychoError::NonFinite { what: "pcm sample" });
    }
    Ok(())
}

/// Run the 28-band filterbank; for each band return the full-rate
/// filtered signal handed to the caller's squaring/smoothing stage.
/// (`f_calc_third_octave_levels`, filtering stages only.)
fn filter_band(pcm_pa: &[f64], band: usize, output: &mut [f64], scratch: &mut [f64]) {
    // Stage 1 reads the input; stages 2 and 3 read the previous
    // stage's output (double-buffered here; the reference reuses one
    // buffer since each stage reads sequentially).
    let mut coeffs = [0.0f64; 6];
    for stage in 0..3 {
        for (c, (r, d)) in coeffs
            .iter_mut()
            .zip(FILTER_REF[stage].iter().zip(&FILTER_DIFF[band][stage]))
        {
            *c = r - d;
        }
        let gain = FILTER_GAIN[band][stage];
        if stage == 0 {
            filter_2nd_order(pcm_pa, output, &coeffs, gain);
        } else {
            scratch.copy_from_slice(output);
            filter_2nd_order(scratch, output, &coeffs, gain);
        }
    }
}

/// Band center frequency exactly as the reference computes it.
#[allow(clippy::cast_precision_loss)]
fn center_frequency(band: usize) -> f64 {
    det::pow(10.0, (band as f64 - 16.0) / 10.0) * 1000.0
}

/// Stationary loudness from PCM (clause 5 via the 48 kHz filterbank):
/// mean-square band levels after `time_skip` seconds, then the
/// third-octave-level method. This is the reference's
/// stationary-from-signal path (`f_loudness_from_signal`,
/// `LoudnessMethodStationary`).
///
/// # Errors
/// [`PsychoError::UnsupportedRate`] off 48 kHz;
/// [`PsychoError::NonFinite`] on bad samples;
/// [`PsychoError::DegenerateSignal`] when `time_skip` consumes the
/// whole signal or `time_skip` is outside the reference's accepted
/// 0..=1 s.
pub fn loudness_stationary_from_pcm(
    pcm_pa: &[f64],
    sample_rate: f64,
    time_skip: f64,
    field: SoundField,
) -> Result<Loudness, PsychoError> {
    validate_signal(pcm_pa, sample_rate)?;
    if !(0.0..=1.0).contains(&time_skip) {
        return Err(PsychoError::DegenerateSignal {
            what: "time_skip must be within 0..=1 s (reference bound)",
        });
    }
    let n = pcm_pa.len();
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    let num_skip = (time_skip * sample_rate).floor() as usize;
    if num_skip >= n {
        return Err(PsychoError::DegenerateSignal {
            what: "time_skip consumes the whole signal",
        });
    }
    let mut levels = [0.0f64; N_THIRD_OCTAVE_BANDS];
    let mut output = vec![0.0f64; n];
    let mut scratch = vec![0.0f64; n];
    for (band, level) in levels.iter_mut().enumerate() {
        filter_band(pcm_pa, band, &mut output, &mut scratch);
        let ms: f64 = output[num_skip..].iter().map(|v| v * v).sum::<f64>() / (n - num_skip) as f64;
        *level = 10.0 * det::ln((ms + TINY_VALUE) / I_REF) / det::ln(10.0);
    }
    crate::loudness_stationary(&levels, field)
}

/// Time-varying loudness from PCM (the clause-6 method): filterbank,
/// squaring + three frequency-dependent smoothing low-passes,
/// decimation to 2 kHz levels, per-frame core loudness, nonlinear
/// temporal decay, slopes, and duration-dependent temporal weighting;
/// summarized by `Nmax` and the standard's `N5` percentile.
///
/// # Errors
/// [`PsychoError::UnsupportedRate`] off 48 kHz;
/// [`PsychoError::NonFinite`] on bad samples;
/// [`PsychoError::DegenerateSignal`] when the signal is shorter than
/// two 0.5 ms level frames.
pub fn loudness_time_varying(
    pcm_pa: &[f64],
    sample_rate: f64,
    field: SoundField,
) -> Result<TimeVaryingLoudness, PsychoError> {
    validate_signal(pcm_pa, sample_rate)?;
    let n = pcm_pa.len();
    let num_frames = n / DEC_FACTOR;
    if num_frames < 2 {
        return Err(PsychoError::DegenerateSignal {
            what: "time-varying loudness needs at least two 0.5 ms frames",
        });
    }
    // --- Filterbank + squaring + smoothing + decimated levels. ---
    let mut levels = vec![[0.0f64; N_THIRD_OCTAVE_BANDS]; num_frames];
    let mut output = vec![0.0f64; n];
    let mut scratch = vec![0.0f64; n];
    for band in 0..N_THIRD_OCTAVE_BANDS {
        filter_band(pcm_pa, band, &mut output, &mut scratch);
        // Frequency-dependent smoothing time constant, capped at the
        // 1 kHz value above 1 kHz.
        let fc = center_frequency(band);
        let tau = if fc <= 1000.0 {
            2.0 / (3.0 * fc)
        } else {
            2.0 / (3.0 * 1000.0)
        };
        for v in &mut output {
            *v *= *v;
        }
        for _ in 0..3 {
            lowpass_in_place(&mut output, tau, sample_rate);
        }
        for (frame, lv) in levels.iter_mut().enumerate() {
            let ms = output[frame * DEC_FACTOR];
            lv[band] = 10.0 * det::ln((ms + TINY_VALUE) / I_REF) / det::ln(10.0);
        }
    }
    // Finite samples can still OVERFLOW the squaring stage (executed:
    // 1e200 Pa read Ok(n_max = inf, n5 = NaN) before this gate) — the
    // stationary path refuses through loudness_stationary's level
    // check; this is the time-varying twin of that refusal.
    if levels.iter().flatten().any(|v| !v.is_finite()) {
        return Err(PsychoError::NonFinite {
            what: "third-octave level (squaring overflow)",
        });
    }
    // --- Core loudness per frame, then nonlinear decay per channel. ---
    let mut core_t: Vec<[f64; 21]> = levels
        .iter()
        .map(|lv| core_loudness_21(lv, field))
        .collect();
    let mut channel = vec![0.0f64; num_frames];
    for ch in 0..21 {
        for (t, c) in core_t.iter().enumerate() {
            channel[t] = c[ch];
        }
        nl_channel(&mut channel);
        for (t, c) in core_t.iter_mut().enumerate() {
            c[ch] = channel[t];
        }
    }
    // --- Slopes per frame, then temporal weighting. ---
    let mut loudness: Vec<f64> = core_t.iter().map(|c| slopes_240(c).0).collect();
    temporal_weight(&mut loudness);
    let n_max = loudness.iter().copied().fold(0.0f64, f64::max);
    let n5 = percentile_n5(&loudness);
    Ok(TimeVaryingLoudness {
        loudness,
        n_max,
        n5,
    })
}

/// Loudness level [phon] from loudness [sone] — the reference's own
/// `f_sone_to_phon` (verified source; previously deliberately
/// absent). Below 1 sone: `40 (N + 0.0005)^0.35`, floored at 3 phon;
/// at and above 1 sone: `10 log2(N) + 40`.
///
/// # Errors
/// [`PsychoError::NonFinite`] on a non-finite input;
/// [`PsychoError::DegenerateSignal`] on a negative loudness.
pub fn phon_from_sone(sones: f64) -> Result<f64, PsychoError> {
    if !sones.is_finite() {
        return Err(PsychoError::NonFinite { what: "loudness" });
    }
    if sones < 0.0 {
        return Err(PsychoError::DegenerateSignal {
            what: "loudness cannot be negative",
        });
    }
    if sones < 1.0 {
        Ok((40.0 * det::pow(sones + 0.0005, 0.35)).max(3.0))
    } else {
        Ok(10.0 * det::ln(sones) / det::ln(2.0) + 40.0)
    }
}

/// The batch metric set for Pareto/listening-evidence consumers: one
/// call, one PCM block in Pa at 48 kHz, every metric this crate can
/// currently claim. Aggregation-exact BY CONTRACT: each field equals
/// the corresponding standalone call on the same input (pinned by
/// test — a wiring mistake like swapped fields is the failure mode
/// this contract exists to catch).
#[derive(Debug, Clone)]
pub struct ParetoMetrics {
    /// Stationary loudness [sone] via the PCM filterbank path
    /// (`time_skip` as passed).
    pub sones_stationary: f64,
    /// Loudness level [phon] of the stationary loudness.
    pub phon_stationary: f64,
    /// Time-varying N5 [sone] (the standard's percentile summary).
    pub n5: f64,
    /// Time-varying Nmax [sone].
    pub n_max: f64,
    /// DIN 45692 sharpness [acum] over the stationary
    /// specific-loudness pattern.
    pub sharpness_acum: f64,
    /// Daniel-Weber roughness [asper], mean over the consecutive
    /// whole [`crate::roughness::DW_BLOCK`] blocks in the signal.
    pub roughness_asper_mean: f64,
    /// Number of whole roughness blocks averaged.
    pub roughness_blocks: usize,
    /// Timbre-toolbox log-attack-time [log10 s] with the given
    /// envelope window.
    pub log_attack_time: f64,
}

/// Compute the full [`ParetoMetrics`] set from one calibrated PCM
/// block (Pa, 48 kHz). Refusals from any component metric propagate
/// unchanged — a batch call never papers over a member's typed
/// refusal with a partial result.
///
/// # Errors
/// Any error of the component metrics
/// ([`loudness_stationary_from_pcm`], [`loudness_time_varying`],
/// [`crate::sharpness_din`], [`crate::roughness::roughness_dw_block`],
/// [`crate::log_attack_time`], [`phon_from_sone`]).
pub fn pareto_metrics(
    pcm_pa: &[f64],
    sample_rate: f64,
    time_skip: f64,
    field: SoundField,
    lat_env_window: usize,
) -> Result<ParetoMetrics, PsychoError> {
    let stationary = loudness_stationary_from_pcm(pcm_pa, sample_rate, time_skip, field)?;
    let tv = loudness_time_varying(pcm_pa, sample_rate, field)?;
    let sharpness = crate::sharpness_din(&stationary.specific)?;
    let n_blocks = pcm_pa.len() / crate::roughness::DW_BLOCK;
    if n_blocks == 0 {
        return Err(PsychoError::DegenerateSignal {
            what: "no whole roughness block in the signal",
        });
    }
    let mut r_sum = 0.0f64;
    for b in 0..n_blocks {
        let start = b * crate::roughness::DW_BLOCK;
        r_sum += crate::roughness::roughness_dw_block(
            &pcm_pa[start..start + crate::roughness::DW_BLOCK],
            sample_rate,
        )?;
    }
    #[allow(clippy::cast_precision_loss)]
    let roughness_mean = r_sum / n_blocks as f64;
    let lat = crate::log_attack_time(pcm_pa, sample_rate, lat_env_window)?;
    Ok(ParetoMetrics {
        sones_stationary: stationary.sones,
        phon_stationary: phon_from_sone(stationary.sones)?,
        n5: tv.n5,
        n_max: tv.n_max,
        sharpness_acum: sharpness,
        roughness_asper_mean: roughness_mean,
        roughness_blocks: n_blocks,
        log_attack_time: lat,
    })
}
