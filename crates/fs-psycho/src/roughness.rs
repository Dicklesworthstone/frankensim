//! Daniel-Weber roughness — the modulation-analysis metric whose
//! band-pass character peaks near 70 Hz modulation.
//!
//! Ported from the Apache-2.0 MoSQITo reference implementation
//! (github.com/Eomys/MoSQITo, commit d990c33; algorithm from Daniel &
//! Weber 1997), with its data in [`crate::dw_tables`]. The port is
//! numerically EXACT against the reference at the same block length
//! (the AM validation matches to 10+ digits at every modulation
//! frequency). Stated deviation: the block is 8192 samples (fs-fft
//! is power-of-two; the reference default is 0.2 s = 9600 at 48 kHz;
//! the reference itself was re-run at 8192 to establish the pinned
//! values). Level-RELATIVE in shape but the
//! asper scale is anchored at 1 asper for a 1 kHz 60 dB tone with
//! 100% 70 Hz AM — inputs are treated as sound pressure in Pa (the
//! caller's calibration turns digital samples into Pa; the
//! uncalibrated-absolute refusal applies upstream).

use crate::PsychoError;
use crate::dw_tables::{
    A0_BARK, A0_DB, BARK_FREQS, GZI_Y, H2_X, H2_Y, H5_X, H5_Y, H16_X, H16_Y, H21_X, H21_Y, H42_X,
    H42_Y, LTQ_R_BARK, LTQ_R_DB,
};
use fs_fft::{C64, Fft};
use fs_math::det;

/// Real scaling for fs-fft's `C64` (its own `scale` is private).
fn cscale(c: C64, k: f64) -> C64 {
    C64::new(c.re * k, c.im * k)
}

/// Analysis block length (samples): power-of-two for fs-fft; 170.7 ms
/// at 48 kHz vs the reference's 200 ms (stated deviation).
pub const DW_BLOCK: usize = 8192;

/// Linear interpolation over a table (numpy-interp semantics: clamps
/// at the ends).
fn interp(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
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

fn freq_to_bark(f: f64) -> f64 {
    let barks: Vec<f64> = (0..BARK_FREQS.len()).map(|i| i as f64 * 0.5).collect();
    interp(f, &BARK_FREQS, &barks)
}

fn bark_to_freq(z: f64) -> f64 {
    let barks: Vec<f64> = (0..BARK_FREQS.len()).map(|i| i as f64 * 0.5).collect();
    interp(z, &barks, &BARK_FREQS)
}

/// The H modulation-bandpass weight for channel `ch` (0-based) at
/// modulation frequency `fm` — the reference's channel-group mapping:
/// channels 0..4 use the H2 anchor, 4..15 H5, 15..20 H16, 20..41 H21,
/// 41..47 H42 (piecewise-constant across groups, as implemented by
/// the reference).
fn h_weight(ch: usize, fm: f64) -> f64 {
    let (xs, ys): (&[f64], &[f64]) = match ch {
        0..=3 => (&H2_X, &H2_Y),
        4..=14 => (&H5_X, &H5_Y),
        15..=19 => (&H16_X, &H16_Y),
        20..=40 => (&H21_X, &H21_Y),
        _ => (&H42_X, &H42_Y),
    };
    if fm <= 0.0 || fm >= xs[xs.len() - 1] {
        return 0.0;
    }
    interp(fm, xs, ys)
}

/// Pearson correlation of two equal-length slices.
fn corrcoef(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        num += (x - ma) * (y - mb);
        da += (x - ma) * (x - ma);
        db += (y - mb) * (y - mb);
    }
    let denom = det::sqrt(da * db);
    if denom > 0.0 { num / denom } else { 0.0 }
}

/// Daniel-Weber roughness of one PCM block (sound pressure in Pa),
/// [asper]. Uses the first [`DW_BLOCK`] samples; refuses shorter
/// input.
///
/// # Errors
/// [`PsychoError::DegenerateSignal`] on short input;
/// [`PsychoError::NonFinite`] on bad samples/rate.
#[allow(clippy::too_many_lines)] // one D&W pipeline, told in the reference's stage order
pub fn roughness_dw_block(pcm_pa: &[f64], sample_rate: f64) -> Result<f64, PsychoError> {
    if pcm_pa.len() < DW_BLOCK {
        return Err(PsychoError::DegenerateSignal {
            what: "roughness needs a full analysis block",
        });
    }
    if pcm_pa.iter().take(DW_BLOCK).any(|v| !v.is_finite()) {
        return Err(PsychoError::NonFinite { what: "pcm sample" });
    }
    if sample_rate.is_nan() || sample_rate <= 0.0 {
        return Err(PsychoError::NonFinite {
            what: "sample rate",
        });
    }
    let n = DW_BLOCK;
    let fs = sample_rate;
    // Blackman window, amplitude-normalized by its sum (reference
    // comp_spectrum convention), then one-sided doubling.
    let mut win = vec![0.0f64; n];
    let mut wsum = 0.0;
    for (i, w) in win.iter_mut().enumerate() {
        let t = i as f64 / (n as f64 - 1.0);
        *w = 0.42 - 0.5 * det::cos(2.0 * core::f64::consts::PI * t)
            + 0.08 * det::cos(4.0 * core::f64::consts::PI * t);
        wsum += *w;
    }
    let fft = Fft::new(n);
    let mut scratch = vec![C64::new(0.0, 0.0); n];
    let mut buf: Vec<C64> = (0..n)
        .map(|i| C64::new(pcm_pa[i] * win[i] / wsum, 0.0))
        .collect();
    fft.forward(&mut buf, &mut scratch);
    // One-sided spectrum exactly as the reference slices it: bins
    // 0..n/2 (DC included) scaled by the reference's empirical 1.42
    // factor, with the frequency axis labeled (k+1)*fs/n — the
    // reference's own one-bin offset, kept for faithfulness (a 2.0
    // doubling instead of 1.42 read 2.59 asper at the calibration
    // anchor — executed).
    let half = n / 2;
    let spec: Vec<C64> = (0..half).map(|k| cscale(buf[k], 1.42)).collect();
    let freq_axis: Vec<f64> = (0..half)
        .map(|k| (k as f64 + 1.0) * fs / n as f64)
        .collect();
    let bark_axis: Vec<f64> = freq_axis.iter().map(|&f| freq_to_bark(f)).collect();
    // a0 outer/middle-ear weighting (applied as amplitude factor).
    let spec_a0: Vec<C64> = spec
        .iter()
        .zip(&bark_axis)
        .map(|(s, &z)| cscale(*s, det::pow(10.0, interp(z, &A0_BARK, &A0_DB) / 20.0)))
        .collect();
    let module: Vec<f64> = spec_a0.iter().map(|s| det::sqrt(s.norm_sq())).collect();
    let spec_db: Vec<f64> = module
        .iter()
        .map(|&m| {
            if m > 0.0 {
                20.0 * det::ln(m / 2.0e-5) / det::ln(10.0)
            } else {
                -400.0
            }
        })
        .collect();
    // Audible components above the roughness-reference threshold.
    let threshold: Vec<f64> = bark_axis
        .iter()
        .map(|&z| interp(z, &LTQ_R_BARK, &LTQ_R_DB))
        .collect();
    let audible: Vec<usize> = (0..half).filter(|&k| spec_db[k] > threshold[k]).collect();
    if audible.is_empty() {
        return Ok(0.0);
    }
    // Terhardt slopes.
    let s1 = -27.0f64;
    let s2: Vec<f64> = audible
        .iter()
        .map(|&k| (-24.0 - 230.0 / freq_axis[k] + 0.2 * spec_db[k]).min(0.0))
        .collect();
    let n_channel = 47usize;
    // Channel center Bark values 0.5..23.5 and their bin positions.
    let zb: Vec<f64> = (1..=n_channel)
        .map(|i| bark_to_freq(i as f64 / 2.0) * n as f64 / fs)
        .collect();
    let min_excit_db: Vec<f64> = zb
        .iter()
        .map(|&b| {
            // interp over the positive-bin index axis (1..=half).
            let axis: f64 = b;
            let k = axis.clamp(1.0, half as f64);
            let lo = (k.floor() as usize - 1).min(half - 1);
            let hi = (lo + 1).min(half - 1);
            let t = k - k.floor();
            threshold[lo] * (1.0 - t) + threshold[hi] * t
        })
        .collect();
    let ch_low: Vec<isize> = audible
        .iter()
        .map(|&k| (2.0 * bark_axis[k]).floor() as isize - 1)
        .collect();
    let ch_high: Vec<isize> = audible
        .iter()
        .map(|&k| (2.0 * bark_axis[k]).ceil() as isize - 1)
        .collect();
    // Excitation slope amplitudes per (component, channel).
    let n_aud = audible.len();
    let mut slopes = vec![0.0f64; n_aud * n_channel];
    for (a, &k) in audible.iter().enumerate() {
        let lev = spec_db[k];
        let b = bark_axis[k];
        let low_end = (ch_low[a] + 1).max(0) as usize;
        for j in 0..low_end.min(n_channel) {
            let sl = s1 * (b - f64::midpoint(j as f64, 1.0)) + lev;
            if sl > min_excit_db[j] {
                slopes[a * n_channel + j] = 2.0e-5 * det::pow(10.0, sl / 20.0);
            }
        }
        let high_start = ch_high[a].max(0) as usize;
        for j in high_start..n_channel {
            let sl = s2[a] * (f64::midpoint(j as f64, 1.0) - b) + lev;
            if sl > min_excit_db[j] {
                slopes[a * n_channel + j] = 2.0e-5 * det::pow(10.0, sl / 20.0);
            }
        }
    }
    // Per-channel: excitation reconstruction, envelope, modulation
    // depth via the H bandpass.
    let mut h_bp = vec![vec![0.0f64; n]; n_channel];
    let mut mod_depth = vec![0.0f64; n_channel];
    for ch in 0..n_channel {
        let mut exc = vec![C64::new(0.0, 0.0); n];
        for (a, &k) in audible.iter().enumerate() {
            let ampl = if ch_low[a] == ch.cast_signed() || ch_high[a] == ch.cast_signed() {
                1.0
            } else if ch_high[a] > ch.cast_signed() {
                if ch + 1 < n_channel && module[k] > 0.0 {
                    slopes[a * n_channel + ch + 1] / module[k]
                } else {
                    0.0
                }
            } else if ch >= 1 && module[k] > 0.0 {
                slopes[a * n_channel + ch - 1] / module[k]
            } else {
                0.0
            };
            // Reconstruct at the reference's own bin index.
            exc[k] = cscale(spec_a0[k], ampl);
        }
        let mut tmp = exc.clone();
        fft.inverse(&mut tmp, &mut scratch);
        let temporal: Vec<f64> = tmp.iter().map(|c| (n as f64 * c.re).abs()).collect();
        let h0 = temporal.iter().sum::<f64>() / n as f64;
        // Envelope spectrum, H-weighted, back to time.
        let mut env: Vec<C64> = temporal.iter().map(|&v| C64::new(v - h0, 0.0)).collect();
        fft.forward(&mut env, &mut scratch);
        // The reference's H array covers ONLY positive-frequency bins
        // from 2 upward (bins 0..2 and the entire negative half are
        // ZERO); the factor 2 in hBP compensates the dropped
        // conjugate half. Mirroring the weight over both halves AND
        // doubling gave 2x envelope amplitude ~ 4x roughness
        // (executed: R(20 Hz) read 0.85 vs the reference's 0.197).
        for (idx, e) in env.iter_mut().enumerate() {
            let w = if (2..=n / 2).contains(&idx) {
                h_weight(ch, idx as f64 * fs / n as f64)
            } else {
                0.0
            };
            *e = cscale(*e, w);
        }
        fft.inverse(&mut env, &mut scratch);
        let hbp: Vec<f64> = env.iter().map(|c| 2.0 * c.re).collect();
        let rms = det::sqrt(hbp.iter().map(|v| v * v).sum::<f64>() / n as f64);
        mod_depth[ch] = if h0 > 0.0 { (rms / h0).min(1.0) } else { 0.0 };
        h_bp[ch] = hbp;
    }
    // Cross-correlation between channels i and i+2, then the Aures
    // gzi weighting and the published 0.25 calibration.
    let mut ki = vec![0.0f64; n_channel];
    for i in 0..n_channel - 2 {
        let any_a = h_bp[i].iter().any(|&v| v != 0.0);
        let any_b = h_bp[i + 2].iter().any(|&v| v != 0.0);
        if any_a && any_b {
            ki[i] = corrcoef(&h_bp[i], &h_bp[i + 2]);
        }
    }
    let gzi_at = |i: usize| -> f64 {
        let z = f64::midpoint(i as f64, 1.0);
        let xs: Vec<f64> = (0..GZI_Y.len()).map(|k| k as f64).collect();
        interp(z, &xs, &GZI_Y)
    };
    let mut r_total = 0.0f64;
    for i in 0..n_channel {
        let r_spec = match i {
            0 | 1 => gzi_at(i) * det::powi(mod_depth[i] * ki[i], 2),
            45 => gzi_at(i) * det::powi(mod_depth[i] * ki[43], 2),
            46 => gzi_at(i) * det::powi(mod_depth[i] * ki[44], 2),
            _ => gzi_at(i) * det::powi(mod_depth[i] * ki[i] * ki[i - 2], 2),
        };
        r_total += r_spec;
    }
    Ok(0.25 * r_total)
}
